// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* PETScDoFLinearSystem.cc                                     (C) 2022-2025 */
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
  void setRelTolerance(Real v) { m_rtol = v; }
  void setAbsTolerance(Real v) { m_atol = v; }

  void setSolver(String v) { m_solver_method = std::string{v.localstr()}; }
  void setPreconditioner(String v) { m_preconditionner_method = std::string{v.localstr()}; }

  CaseOptionsPETScDoFLinearSystemFactory *options;

 private:

  KSP m_petsc_solver_context;
  Vec m_petsc_solution_vector;
  Vec m_petsc_rhs_vector;
  Mat m_petsc_matrix;

 private:

  VariableDoFInt32 m_dof_matrix_numbering;

  NumArray<Int32, MDDim1> m_parallel_columns_index;
  NumArray<Int32, MDDim1> m_parallel_rows_index;
  //! Work array to store values of solution vector in parallel
  NumArray<Real, MDDim1> m_result_work_values;

  Int32 m_nb_total_row;
  Int32 m_nb_own_row;
  Int32 m_max_iter;

  Real m_rtol;
  Real m_atol;

  std::string m_solver_method; // cannot use String type because we need this to be mutable
  std::string m_preconditionner_method;

 private:

  void _computeMatrixNumeration();
  void _handleParameters();
};

void PETScDoFLinearSystemImpl::_handleParameters()
{
  PetscBool is_initialized;
  PetscInitialized(&is_initialized);

  if (!is_initialized) // no command line arguments were given
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);

#define MAX_STRING_LENGTH 256

#define X_OPTION(type, petsc_string, variable) \
  PetscOptionsGet##type(nullptr, nullptr, petsc_string, &variable, &set); \
  if (!set) \
    PetscOptionsSetValue(nullptr, petsc_string, std::to_string(variable).c_str());

#define X_OPTION_STRING(petsc_string, variable) \
  variable.resize(MAX_STRING_LENGTH); \
  PetscOptionsGetString(nullptr, nullptr, petsc_string, variable.data(), MAX_STRING_LENGTH, &set); \
  if (!set) \
    PetscOptionsSetValue(nullptr, petsc_string, variable.c_str());

  // we do variable.resize(MAX_STRING_LENGTH) to make sure there is enough memory to store the option name.
  // may otherwise result in a segfault

  PetscBool set;
  X_OPTION(Real, "-ksp_rtol", m_rtol);
  X_OPTION(Real, "-ksp_atol", m_atol);
  X_OPTION(Int, "-ksp_max_it", m_max_iter);
  X_OPTION_STRING("-ksp_type", m_solver_method);
  X_OPTION_STRING("-pc_type", m_preconditionner_method);
}

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
  m_nb_own_row = own_dofs.size();
  m_nb_total_row = all_dofs.size();
  Int32 own_first_index = 0;

  if (is_parallel) {
    // TODO: utiliser un Scan lorsque ce sera disponible dans Arcane
    UniqueArray<Int32> parallel_rows_index(nb_rank, 0);
    pm->allGather(ConstArrayView<Int32>(1, &m_nb_own_row), parallel_rows_index);
    info() << "ALL_NB_ROW = " << parallel_rows_index;
    m_nb_total_row = 0;
    // TODO optimize partial and total sum
    for (Int32 v : parallel_rows_index) m_nb_total_row += v;
    for (Int32 i = 0; i < my_rank; ++i)
      own_first_index += parallel_rows_index[i];
  }

  ENUMERATE_DOF (idof, own_dofs) {
    DoF dof = *idof;
    m_dof_matrix_numbering[idof] = own_first_index + idof.index();
    info() << "Numbering dof_uid=" << dof.uniqueId() << " M=" << m_dof_matrix_numbering[idof];
  }
  m_dof_matrix_numbering.synchronize();

  info() << "my rank: " << my_rank;
  info() << "Total " << m_nb_total_row << " local: " << m_nb_own_row;
  pm->barrier();

  m_parallel_rows_index.resize(m_nb_own_row);
  m_result_work_values.resize(m_nb_own_row);

}
void PETScDoFLinearSystemImpl::
solve()
{
  info() << "[PETSc-Info] Calling PETSc solver";
  _handleParameters();
  _computeMatrixNumeration();

  IItemFamily* dof_family = dofFamily();
  IParallelMng* pm = dof_family->parallelMng();
  Parallel::Communicator arcane_comm = pm->communicator();
  bool is_parallel = pm->isParallel();
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  if (arcane_comm.isValid())
    mpi_comm = static_cast<MPI_Comm>(arcane_comm);

  Span<const Int32> rows_index_span = m_dof_matrix_numbering.asArray();
  // for (auto a : rows_index_span) info() << a;
  info() << "size: " << rows_index_span.size() << "local rows: " << m_nb_own_row;
  CSRFormatView csr_view = this->getCSRValues();

  PetscInt local_rows = m_nb_own_row;          // rows this rank owns
  PetscInt global_rows = m_nb_total_row; // total rows across all ranks

  // PetscSplitOwnership(mpi_comm, &local_rows, &global_rows);

  MatCreate(mpi_comm, &m_petsc_matrix);
  MatSetSizes(m_petsc_matrix, local_rows, local_rows, global_rows, global_rows);
  MatSetFromOptions(m_petsc_matrix);
  MatSetUp(m_petsc_matrix);

  // m_csr_view.columns() use matrix coordinates local to sub-domain
  // We need to translate them to global matrix coordinates
  Span<const Int32> columns_index_span = csr_view.columns();

  info() << "COLUMNS=" << csr_view.columns();

  if (is_parallel)
  {
    // TODO: Faire sur accélérateur et ne faire qu'une fois si la structure
    // ne change pas.
    Int64 nb_column = columns_index_span.size();
    m_parallel_columns_index.resize(nb_column);
    for (Int64 i = 0; i < nb_column; ++i)
    {
      DoFLocalId lid(columns_index_span[i]);
      info() << "I=" << i << " index=" << columns_index_span[i] << " lid: " << lid;
      // Si lid correspond à une entité nulle, alors la valeur de la matrice
      // ne sera pas utilisée.
      if (!lid.isNull())
        m_parallel_columns_index[i] = m_dof_matrix_numbering[lid];
      else
        m_parallel_columns_index[i] = 0;
    }
    columns_index_span = m_parallel_columns_index.to1DSpan();
  }

  Span<const Real> matrix_values = csr_view.values();


  const PetscInt* rows_index_data = rows_index_span.data();
  const PetscInt* columns_index_data = columns_index_span.data();
  const PetscReal* matrix_values_data = matrix_values.data();

  pm->barrier();

  // assemble matrix row by row
  ENUMERATE_ (DoF, idof, dof_family->allItems())
  {
    DoF dof = *idof;

    if (!dof.isOwn()) // for obscure reasons, allItems().own() doesnt work :'(
      continue;

    PetscInt nb_col = csr_view.rowsNbColumn()[idof.index()];
    PetscInt row_csr_index = csr_view.rows()[idof.index()];
    PetscInt global_row = rows_index_span[idof.index()];
    // info() << "global_row: " << global_row << " nb_col: " << nb_col << " row_csr_index: " << row_csr_index;
    const PetscInt* cols = &columns_index_data[row_csr_index];
    const PetscScalar* vals = &matrix_values_data[row_csr_index];
    MatSetValues(m_petsc_matrix, 1, &global_row, nb_col, cols, vals, INSERT_VALUES);
  }

  pm->barrier();

  MatAssemblyBegin(m_petsc_matrix, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m_petsc_matrix, MAT_FINAL_ASSEMBLY);

  if (is_parallel)
  {
    // Fill 'm_parallel_rows_index' with only rows we owns
    // NOTE: This is only needed if matrix structure has changed.
    Int32 index = 0;
    ENUMERATE_ (DoF, idof, dof_family->allItems())
    {
      DoF dof = *idof;
      if (!dof.isOwn())
        continue;
      m_parallel_rows_index[index] = rows_index_span[idof.index()];
      ++index;
    }
  }

  info() << "[PETSc-Info] Created Matrix";


  const Real* rhs_data = rhsVariable().asArray().data();
  const Real* result_data = solutionVariable().asArray().data();
  VariableDoFReal& dof_variable = this->solutionVariable();

  VecCreateMPI(mpi_comm, local_rows, global_rows, &m_petsc_rhs_vector);
  VecSetSizes(m_petsc_rhs_vector, local_rows, global_rows);
  VecSetFromOptions(m_petsc_rhs_vector);
  VecSetValues(m_petsc_rhs_vector, local_rows, rows_index_data, rhs_data, INSERT_VALUES);
  VecAssemblyBegin(m_petsc_rhs_vector);
  VecAssemblyEnd(m_petsc_rhs_vector);

  VecCreateMPI(mpi_comm, local_rows, global_rows, &m_petsc_solution_vector);
  VecSetSizes(m_petsc_solution_vector, local_rows, global_rows);
  VecSetFromOptions(m_petsc_solution_vector);
  VecSetValues(m_petsc_solution_vector, local_rows, rows_index_data, result_data, INSERT_VALUES);
  VecAssemblyBegin(m_petsc_solution_vector);
  VecAssemblyEnd(m_petsc_solution_vector);

  info() << "[PETSc-Info] Created vectors";

  KSPCreate(mpi_comm, &m_petsc_solver_context);
  KSPSetOperators(m_petsc_solver_context, m_petsc_matrix, m_petsc_matrix);
  KSPSetFromOptions(m_petsc_solver_context);
  KSPSolve(m_petsc_solver_context, m_petsc_rhs_vector, m_petsc_solution_vector);

  info() << "[PETSc-Info] Solved linear system";

  if (is_parallel) {
    Int32 nb_wanted_row = m_parallel_rows_index.extent0();
    PetscInt* idx = m_parallel_rows_index.to1DSpan().data();

    VecGetValues(m_petsc_solution_vector, nb_wanted_row, idx, m_result_work_values.to1DSpan().data());

    ENUMERATE_ (DoF, idof, dof_family->allItems().own()) {
      info() << "index: " << idof.index();
      info() << "value: " << m_result_work_values[idof.index()];
      dof_variable[idof] = m_result_work_values[idof.index()];
    }
    info() << "ayaaaaa";
  }
  else {
    PetscScalar* vals = nullptr;
    VecGetArray(m_petsc_solution_vector, &vals);

    // Copy directly into Arcane DoF variable
    dof_variable.asArray().copy(Span<const Real>(vals, m_nb_own_row));

    VecRestoreArray(m_petsc_solution_vector, &vals);
  }

  info() << "[PETSc-Info] Wrote solution in solution_variable";
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
    x->setMaxIter(options()->maxIter());
    x->setSolver(options()->solver());
    x->setPreconditioner(options()->preconditioner());
    return x;
  }
};


ARCANE_REGISTER_SERVICE_PETSCDOFLINEARSYSTEMFACTORY(PETScLinearSystem,
                                                    PETScDoFLinearSystemFactoryService);

} // namespace Arcane::FemUtils

