// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* HypreDoFLinearSystem.cc                                     (C) 2022-2025 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "DoFLinearSystem.h"

#include <arcane/accelerator/RunCommandLoop.h>
#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/PlatformUtils.h>
#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/MemoryUtils.h>
#include <arcane/utils/MemoryView.h>
#include <arcane/utils/ITraceMng.h>
#include <arcane/utils/NumArray.h>
#include <arcane/utils/CommandLineArguments.h>

#include <arcane/core/ServiceFactory.h>
#include <arcane/core/VariableTypes.h>
#include <arcane/core/BasicService.h>
#include <arcane/core/VariableUtils.h>
#include <arcane/core/IParallelMng.h>
#include <arcane/core/IItemFamily.h>
#include <arcane/core/ItemPrinter.h>
#include <arcane/core/Timer.h>
#include <arcane/aleph/AlephTypesSolver.h> // TO MODIFY

#include <arcane/accelerator/VariableViews.h>
#include <arcane/accelerator/core/Runner.h>
#include <arcane/accelerator/core/Memory.h>

#include "IDoFLinearSystemFactory.h"
#include "internal/CsrDoFLinearSystemImpl.h"

#include <petsc.h>
#include <petscpctypes.h>
#include <krylov.h>
#include "PETScDoFLinearSystemFactory_axl.h"

namespace Arcane::FemUtils
{

using namespace Arcane;

class PETScDoFLinearSystemImpl
: public CsrDoFLinearSystemImpl
{
 public:

 PETScDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
  : CsrDoFLinearSystemImpl(dof_family, solver_name)
  , m_dof_matrix_numbering(VariableBuildInfo(dof_family, solver_name + "MatrixNumbering"))
  {
    info() << "[PETSc-Info] Creating PETScDoFLinearSystemImpl()";
  }

  ~PETScDoFLinearSystemImpl() override
  {
    info() << "[PETSc-Info] Calling PETScDoFLinearSystemImpl destructor";
    PetscFinalize();
  }

 public:

  void build() {
    // _computeMatrixNumeration();
    PetscFunctionBeginUser;
  }

  void solve() override;

  void setSolverCommandLineArguments(const CommandLineArguments& args) override
  {
    PetscInitialize(args.commandLineArgc(), args.commandLineArgv(), nullptr, nullptr);
    info() << "[PETSc-Info] initialize command lines arguments";
    auto argv = *args.commandLineArgv();
    auto o = info() << "[PETSc-Info] ./" << argv[0];

    for (int i = 1; i < *args.commandLineArgc(); i++)
      o << ' ' << argv[i];
  }

  void setMaxIter(Int32 v) { m_max_iter = v; }
  void setVerbosityLevel(Int32 v) { m_verbosity = v; }
  void setAmgCoarsener(Int32 v) { m_amg_coarsener = v; }
  void setAmgInterpType(Int32 v) { m_amg_interp_type = v; }
  void setAmgSmoother(Int32 v) { m_amg_smoother = v; }
  void setKrylovDim(Int32 v) { m_krylov_dim = v; }

  void setRelTolerance(Real v) { m_rtol = v; }
  void setAbsTolerance(Real v) { m_atol = v; }
  void setAmgThreshold(Real v) { m_amg_threshold = v; }

  CaseOptionsPETScDoFLinearSystemFactory *options;

 private:

  KSP m_petsc_solver_context;
  PC m_petsc_preconditioner_context;
  Vec m_petsc_solution_vector;
  Vec m_petsc_rhs_vector;
  Mat m_petsc_matrix;

 private:

  VariableDoFInt32 m_dof_matrix_numbering;

  NumArray<Int32, MDDim1> m_parallel_columns_index;
  NumArray<Int32, MDDim1> m_parallel_rows_index;
  //! Work array to store values of solution vector in parallel
  NumArray<Real, MDDim1> m_result_work_values;

  Int32 m_first_own_row;
  Int32 m_nb_own_row;
  Int32 m_max_iter;
  Int32 m_verbosity;
  Int32 m_amg_coarsener;
  Int32 m_amg_interp_type;
  Int32 m_amg_smoother;
  Int32 m_krylov_dim;

  Real m_amg_threshold;
  Real m_rtol;
  Real m_atol;

 private:

  void _computeMatrixNumeration();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Compute global numeration of the matrix.
 *
 * Each rank owns consecutive rows of the matrix in increasing order.
 */
void PETScDoFLinearSystemImpl::
_computeMatrixNumeration()
{
  IItemFamily* dof_family = dofFamily();
  IParallelMng* pm = dof_family->parallelMng();
  const bool is_parallel = pm->isParallel();
  const Int32 nb_rank = pm->commSize();
  const Int32 my_rank = pm->commRank();

  DoFGroup all_dofs = dof_family->allItems();
  DoFGroup own_dofs = all_dofs.own();
  const Int32 nb_own_row = own_dofs.size();

  Int32 own_first_index = 0;

  if (is_parallel) {
    // TODO: utiliser un Scan lorsque ce sera disponible dans Arcane
    UniqueArray<Int32> parallel_rows_index(nb_rank, 0);
    pm->allGather(ConstArrayView<Int32>(1, &nb_own_row), parallel_rows_index);
    info() << "ALL_NB_ROW = " << parallel_rows_index;
    for (Int32 i = 0; i < my_rank; ++i)
      own_first_index += parallel_rows_index[i];
  }

  info() << "OwnFirstIndex=" << own_first_index << " NbOwnRow=" << nb_own_row;

  m_first_own_row = own_first_index;
  m_nb_own_row = nb_own_row;

  // TODO: Faire avec API accelerateur
  ENUMERATE_DOF (idof, own_dofs) {
    DoF dof = *idof;
    m_dof_matrix_numbering[idof] = own_first_index + idof.index();
    //info() << "Numbering dof_uid=" << dof.uniqueId() << " M=" << m_dof_matrix_numbering[idof];
  }
  info() << " nb_own_row=" << nb_own_row << " nb_item=" << dof_family->nbItem();
  m_dof_matrix_numbering.synchronize();

  m_parallel_rows_index.resize(nb_own_row);
  m_result_work_values.resize(nb_own_row);
}

void PETScDoFLinearSystemImpl::
solve()
{
  info() << "[PETSc-Info] Calling PETSc solver";
  _applyRowElimination();
  _applyForcedValuesToLhs();

  PetscBool isInitialized;
  PetscInitialized(&isInitialized);

  if (!isInitialized) // no command line arguments were given
  {
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);
    KSPSetTolerances(m_petsc_solver_context, m_rtol, m_atol, PETSC_DEFAULT, m_max_iter);
  }
  else
    KSPSetFromOptions(m_petsc_solver_context);

 // TODO
}

#include "PETScDoFLinearSystemFactory_axl.h"

class PETScDoFLinearSystemFactoryService
: public ArcanePETScDoFLinearSystemFactoryObject
{
 public:

  explicit PETScDoFLinearSystemFactoryService(const ServiceBuildInfo& sbi)
  : ArcanePETScDoFLinearSystemFactoryObject(sbi)
  {
    info() << "[PETSc-Info] Create PETScDoF";
  }
;
  IDoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name) override
  {
    auto* x = new PETScDoFLinearSystemImpl(dof_family, solver_name);

    x->options = options();

    x->build();
    x->setRelTolerance(options()->rtol());
    x->setAbsTolerance(options()->atol());
    x->setAmgThreshold(options()->amgThreshold());
    x->setMaxIter(options()->maxIter());
    x->setAmgCoarsener(options()->amgCoarsener());
    x->setAmgInterpType(options()->amgInterpType());
    x->setAmgSmoother(options()->amgSmoother());
    x->setKrylovDim(options()->krylovDim());
    x->setVerbosityLevel(options()->verbosity());
    // x->setSolver(options()->solver());
    // x->setPreconditioner(options()->preconditioner());
    return x;
  }
};


ARCANE_REGISTER_SERVICE_PETSCDOFLINEARSYSTEMFACTORY(PETScLinearSystem,
                                                    PETScDoFLinearSystemFactoryService);

} // namespace Arcane::FemUtils

