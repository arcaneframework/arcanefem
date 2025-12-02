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

#include <mpi.h>
#include <arccore/trace/TraceAccessor.h>

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
#include <arcane/core/IParallelMng.h>
#include <arcane/core/IItemFamily.h>
#include <arcane/core/ItemPrinter.h>
#include <arcane/core/Timer.h>

#include <arcane/accelerator/VariableViews.h>
#include <arcane/accelerator/core/Runner.h>
#include <arcane/accelerator/core/Memory.h>

#include "IDoFLinearSystemFactory.h"
#include "internal/CsrDoFLinearSystemImpl.h"

#include <alien/arcane_tools/accessors/ItemVectorAccessor.h>
#include <alien/core/block/VBlock.h>

#include <alien/arcane_tools/IIndexManager.h>
#include <alien/arcane_tools/indexManager/BasicIndexManager.h>
#include <alien/arcane_tools/indexManager/SimpleAbstractFamily.h>
#include <alien/arcane_tools/distribution/DistributionFabric.h>
#include <alien/arcane_tools/indexSet/IndexSetFabric.h>
#include <alien/arcane_tools/data/Space.h>

#include <alien/kernels/simple_csr/algebra/SimpleCSRLinearAlgebra.h>
#include <alien/kernels/simple_csr/algebra/SimpleCSRInternalLinearAlgebra.h>

#include <alien/ref/AlienRefSemantic.h>

#include <alien/kernels/redistributor/Redistributor.h>
#include <alien/ref/data/scalar/RedistributedVector.h>
#include <alien/ref/data/scalar/RedistributedMatrix.h>
#include <alien/ref/import_export/MatrixMarketSystemWriter.h>

#include <alien/expression/solver/SolverStater.h>
#include <alien/AlienLegacyConfig.h>

#include <AlienDoFLinearSystemFactory_axl.h>

#ifdef ALIEN_USE_PETSC
#include <alien/kernels/petsc/io/AsciiDumper.h>
#include <alien/kernels/petsc/algebra/PETScLinearAlgebra.h>
#endif

#ifdef ALIEN_USE_HYPRE
#include <alien/kernels/hypre/HypreBackEnd.h>
#include <alien/kernels/hypre/data_structure/HypreMatrix.h>
#include <alien/kernels/hypre/data_structure/HypreVector.h>
#include <alien/kernels/hypre/algebra/HypreLinearAlgebra.h>
#endif

namespace Arcane::FemUtils
{

using namespace Arcane;
using namespace Alien;

class AlienDoFLinearSystemImpl
: public CsrDoFLinearSystemImpl
{
 public:

  AlienDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
  : CsrDoFLinearSystemImpl(dof_family, solver_name)
  , m_dof_matrix_numbering(VariableBuildInfo(dof_family, solver_name + "MatrixNumbering"))
  {
    info() << "[Alien-Info] Creating AlienDoFLinearSystemImpl()";
  }

  ~AlienDoFLinearSystemImpl() override
  {
    info() << "[Alien-Info] Calling AlienDoFLinearSystemImpl destructor";
  }

 public:

  void solve() override;

  void setSolverCommandLineArguments(const CommandLineArguments& args) override
  {
    info() << "[Alien-Info] initialize command lines arguments";
    auto argv = *args.commandLineArgv();
    auto o = info() << "[Alien-Info] ./" << argv[0];

    for (int i = 1; i < *args.commandLineArgc(); i++)
      o << ' ' << argv[i];
  }

  void setMaxIter(Int32 v) { m_ksp_max_it = v; }
  void setRelTolerance(Real v) { m_ksp_rtol = v; }
  void setAbsTolerance(Real v) { m_ksp_atol = v; }

  void setSolver(Alien::ILinearSolver* s) { m_solver_backend = s; }
  void setPreconditioner(String v) { m_pc_type = std::string{ v.localstr() }; }

  CaseOptionsAlienDoFLinearSystemFactory* options;

 private:

  Alien::ILinearSolver* m_solver_backend;
  VariableDoFInt32 m_dof_matrix_numbering;

  NumArray<Int32, MDDim1> m_parallel_columns_index;
  NumArray<Int32, MDDim1> m_parallel_rows_index;
  NumArray<Real, MDDim1> m_rhs_work_values;
  //! Work array to store values of solution vector in parallel
  NumArray<Real, MDDim1> m_result_work_values;

  Int32 m_nb_total_row;
  Int32 m_nb_own_row;
  Int32 m_first_row;
  Int32 m_ksp_max_it;

  Real m_ksp_rtol;
  Real m_ksp_atol;

  std::string m_ksp_type; // cannot use String type because we need this to be mutable
  std::string m_pc_type;
  std::string m_mat_type;
  std::string m_vec_type;

 private:

  void _computeMatrixNumeration();
};

void AlienDoFLinearSystemImpl::
_computeMatrixNumeration()
{
  // TODO use Index manager structure
  IItemFamily* dof_family = dofFamily();
  IParallelMng* pm = dof_family->parallelMng();
  const bool is_parallel = pm->isParallel();
  const Int32 nb_rank = pm->commSize();
  const Int32 my_rank = pm->commRank();

  DoFGroup all_dofs = dof_family->allItems();
  DoFGroup own_dofs = all_dofs.own();
  m_nb_own_row = own_dofs.size();
  m_nb_total_row = all_dofs.size();
  m_first_row = 0;

  if (is_parallel) {
    // TODO: utiliser un Scan lorsque ce sera disponible dans Arcane
    UniqueArray<Int32> parallel_rows_index(nb_rank, 0);
    pm->allGather(ConstArrayView<Int32>(1, &m_nb_own_row), parallel_rows_index);
    // info() << "ALL_NB_ROW = " << parallel_rows_index;
    m_nb_total_row = 0;
    // TODO optimize partial and total sum
    for (Int32 v : parallel_rows_index)
      m_nb_total_row += v;
    for (Int32 i = 0; i < my_rank; ++i)
      m_first_row += parallel_rows_index[i];
  }

  // TODO api accelerator
  ENUMERATE_DOF (idof, own_dofs) {
    DoF dof = *idof;
    m_dof_matrix_numbering[idof] = m_first_row + idof.index();
    // info() << "Numbering dof_uid=" << dof.uniqueId() << " M=" << m_dof_matrix_numbering[idof];
  }

  m_dof_matrix_numbering.synchronize();
  pm->barrier();

  // info() << "my rank: " << my_rank;
  // info() << "Total " << m_nb_total_row << " local: " << m_nb_own_row;

  m_parallel_rows_index.resize(m_nb_own_row);
  m_result_work_values.resize(m_nb_own_row);
  m_rhs_work_values.resize(m_nb_own_row);
}

void AlienDoFLinearSystemImpl::
solve()
{
  info() << "[Alien-Info] Calling Alien solver";

  // _computeMatrixNumeration();

  IItemFamily* dof_family = dofFamily();
  IParallelMng* pm = dof_family->parallelMng();
  Runner runner = this->runner();
  // MPI_Comm mpi_comm = MPI_COMM_WORLD;
  CSRFormatView csr_view = this->getCSRValues();

  auto areaU = dof_family->allItems();

  Alien::ArcaneTools::BasicIndexManager index_manager(pm);
  index_manager.setTraceMng(traceMng());

  auto indexSetU = index_manager.buildScalarIndexSet("U", areaU);
  index_manager.prepare();
  UniqueArray<Arccore::Integer> allUIndex = index_manager.getIndexes(indexSetU);

  Alien::ArcaneTools::Space space(&index_manager, "TestSpace");
  auto mdist = Alien::ArcaneTools::createMatrixDistribution(space);
  auto vdist = Alien::ArcaneTools::createVectorDistribution(space);

  Alien::Vector vectorB(vdist);
  Alien::Vector vectorX(vdist);

  Alien::Matrix matrixA(mdist); // local matrix for exact measure without side effect
  // (however, you can reuse a matrix with several
  // builder)

  {
    Alien::MatrixProfiler profiler(matrixA);
    ///////////////////////////////////////////////////////////////////////////
    //
    // DEFINE PROFILE
    //
    for (int i = 0; i < csr_view.rows().size() ; i++) {
      for (CsrRowColumnIndex csr_index : csr_view.rowRange(i)) {
        Int32 column_index = csr_view.column(csr_index);
        profiler.addMatrixEntry(i, column_index);
      }
    }
  }
  {
    Alien::ProfiledMatrixBuilder builder(
    matrixA, Alien::ProfiledMatrixOptions::eResetValues);

    for (int i = 0; i < csr_view.rows().size() ; i++) {
      for (CsrRowColumnIndex csr_index : csr_view.rowRange(i)) {
        Int32 column_index = csr_view.column(csr_index);
        // info() << "row: " << i << " col: " << column_index << " val: " << csr_view.value(csr_index);
        builder(i, column_index) += csr_view.value(csr_index);
      }
    }
    builder.finalize();
  }

  VariableDoFReal& rhs_variable = this->rhsVariable();
  VariableDoFReal& dof_variable = this->solutionVariable();
  auto rhs_data = rhs_variable.asArray();
  auto result_data = dof_variable.asArray();

  {
    Alien::VectorWriter writer(vectorX);

    ENUMERATE_DOF (idof, areaU.own()) {
      const Integer iIndex = allUIndex[idof->localId()];
      // info() << "local " << idof->localId( )<< " iIndex " << iIndex << " res val: " << result_data[idof->localId()];
      writer[iIndex] = result_data[idof->localId()];
    }
  }

  {
    Alien::LocalVectorReader reader(vectorX);
    ENUMERATE_DOF(idof, areaU.own()) {
      const Integer iIndex = allUIndex[idof->localId()];
      // info() << "res val:" <<  reader[iIndex];
    }
  }

  {
    Alien::VectorWriter writer(vectorB);

    ENUMERATE_DOF (idof, areaU.own()) {
      const Integer iIndex = allUIndex[idof->localId()];
      // info() << "local " << idof->localId( )<< " iIndex " << iIndex << " rhs val: " << rhs_data[idof->localId()];
      writer[iIndex] = rhs_data[idof->localId()];
    }
  }

  m_solver_backend->solve(matrixA, vectorB, vectorX);

  Alien::SolverStatus status = m_solver_backend->getStatus();

  if (status.succeeded) {
    info()<<"RESOLUTION SUCCEED";

    Alien::VectorReader reader(vectorX);
    // TODO understand this
    ENUMERATE_DOF(idof, areaU.own()) {
      const Integer iIndex = allUIndex[idof->localId()];
      // info() << "val:" <<  reader[iIndex];
      dof_variable[idof] = reader[iIndex];
    }
  }
  else
    info()<<"SOLVER FAILED";
  m_solver_backend->getSolverStat().print(Universe().traceMng(), status, "Linear Solver : ");
}

class AlienDoFLinearSystemFactoryService
: public ArcaneAlienDoFLinearSystemFactoryObject
{
 public:

  explicit AlienDoFLinearSystemFactoryService(const ServiceBuildInfo& sbi)
  : ArcaneAlienDoFLinearSystemFactoryObject(sbi)
  {
    info() << "[Alien-Info] Create AlienDoF";
  };
  IDoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name) override
  {
    auto* x = new AlienDoFLinearSystemImpl(dof_family, solver_name);
    x->options = options();

    if (options()->linearSolver.size() > 0) {
      Alien::ILinearSolver* solver_backend = options()->linearSolver[0];
      x->setSolver(solver_backend);
      solver_backend->init();
    }

    return x;
  }
};

ARCANE_REGISTER_SERVICE_ALIENDOFLINEARSYSTEMFACTORY(AlienLinearSystem,
                                                    AlienDoFLinearSystemFactoryService);

} // namespace Arcane::FemUtils