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

  void setMaxIter(Int32 v) { m_ksp_max_it = v; }
  void setRelTolerance(Real v) { m_ksp_rtol = v; }
  void setAbsTolerance(Real v) { m_ksp_atol = v; }

  void setSolver(String v) { m_ksp_type = std::string{v.localstr()}; }
  void setPreconditioner(String v) { m_pc_type = std::string{v.localstr()}; }

  CaseOptionsPETScDoFLinearSystemFactory *options;

 private:

  KSP m_petsc_solver_context;
  Vec m_petsc_solution_vector;
  Vec m_petsc_rhs_vector;
  Mat m_petsc_matrix;
  ISLocalToGlobalMapping m_petsc_map;

 private:

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

  void _computeMatrixNumeration(MPI_Comm mpi_comm);
  void _handleParameters(IParallelMng* pm);
};

void PETScDoFLinearSystemImpl::_handleParameters(IParallelMng* pm)
{
  PetscBool is_initialized;
  PetscInitialized(&is_initialized);

  if (!is_initialized) // no petsc_flags were given
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);

#define MAX_STRING_LENGTH 256

#define PETSC_OPTION(type, petsc_string, variable) \
  PetscOptionsGet##type(nullptr, nullptr, petsc_string, &variable, &set); \
  if (!set) \
    PetscOptionsSetValue(nullptr, petsc_string, std::to_string(variable).c_str());

#define PETSC_OPTION_STRING(petsc_string, variable) \
  variable.resize(MAX_STRING_LENGTH); \
  PetscOptionsGetString(nullptr, nullptr, petsc_string, variable.data(), MAX_STRING_LENGTH, &set); \
  if (!set) \
    PetscOptionsSetValue(nullptr, petsc_string, variable.c_str());                                 \
  variable.resize(std::strlen(variable.data()));
  // we do variable.resize(MAX_STRING_LENGTH) to make sure there is enough memory to store the option name.
  // may otherwise result in a segfault
  // after this, we resize to the right size to avoid bugs with info()

  Runner runner = this->runner();

  if (pm->isParallel())
  {
    m_mat_type = "mpiaij";
    m_vec_type = "mpi";
  }
  else
  {
    m_mat_type = "seqaij";
    m_vec_type = "seq";
  }

  if (runner.isInitialized())
  {
    switch (runner.executionPolicy()) {
    case Accelerator::eExecutionPolicy::CUDA:
      m_vec_type += "cuda";
      m_mat_type += "cusparse";
      break;
    case Accelerator::eExecutionPolicy::HIP:
      m_vec_type += "hip";
      m_mat_type += "hipsparse";
      break;
    default:
      break;
    }
  }

  PetscBool set;
  PETSC_OPTION(Real, "-ksp_rtol", m_ksp_rtol);
  PETSC_OPTION(Real, "-ksp_atol", m_ksp_atol);
  PETSC_OPTION(Int, "-ksp_max_it", m_ksp_max_it);
  PETSC_OPTION_STRING("-ksp_type", m_ksp_type);
  PETSC_OPTION_STRING("-pc_type", m_pc_type);
  PETSC_OPTION_STRING("-mat_type", m_mat_type);
  PETSC_OPTION_STRING("-vec_type", m_vec_type);

  info() << "[PETSc-Info] Using " << m_mat_type << " matrix type";
  info() << "[PETSc-Info] Using " << m_vec_type << " vector type";
}

void PETScDoFLinearSystemImpl::
_computeMatrixNumeration(MPI_Comm mpi_comm)
{
  // TODO use ISLocalToGlobalMappingCreate PETSc struct
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
    for (Int32 v : parallel_rows_index) m_nb_total_row += v;
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

  NumArray<PetscInt, MDDim1> indices{all_dofs.size()};
  // info() << "nb total row " << m_nb_total_row;

  ENUMERATE_DOF (idof, all_dofs)
  {
    indices[idof.index()] = m_dof_matrix_numbering[idof];
    // info() << "local index: " << idof.index() << " global index: " << indices[idof.index()];
  }

  PetscCallAbort(mpi_comm, ISLocalToGlobalMappingCreate(mpi_comm, 1, all_dofs.size(), indices._internalData(), PETSC_COPY_VALUES, &m_petsc_map));

  // info() << "my rank: " << my_rank;
  // info() << "Total " << m_nb_total_row << " local: " << m_nb_own_row;
  pm->barrier();

  m_parallel_rows_index.resize(m_nb_own_row);
  m_result_work_values.resize(m_nb_own_row);
  m_rhs_work_values.resize(m_nb_own_row);
}

void PETScDoFLinearSystemImpl::
solve()
{
  info() << "[PETSc-Info] Calling PETSc solver";

  IItemFamily* dof_family = dofFamily();
  IParallelMng* pm = dof_family->parallelMng();
  Runner runner = this->runner();
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  _handleParameters(pm);
  _computeMatrixNumeration(mpi_comm);

  Parallel::Communicator arcane_comm = pm->communicator();
  ITimeStats* tstat = pm->timeStats();
  bool is_parallel = pm->isParallel();

  if (arcane_comm.isValid())
    mpi_comm = static_cast<MPI_Comm>(arcane_comm);

  CSRFormatView csr_view = this->getCSRValues();

  PetscInt local_rows = m_nb_own_row;          // rows this rank owns
  PetscInt global_rows = m_nb_total_row; // total rows across all ranks
  Span<const Int32> columns_index_span = csr_view.columns();

  // info() << "COLUMNS=" << csr_view.columns();

  if (is_parallel)
  {
    // TODO: Faire sur accélérateur et ne faire qu'une fois si la structure
    // ne change pas.
    Int64 nb_column = columns_index_span.size();
    m_parallel_columns_index.resize(nb_column);
    for (Int64 i = 0; i < nb_column; ++i)
    {
      DoFLocalId lid(columns_index_span[i]);
      // info() << "I=" << i << " index=" << columns_index_span[i] << " lid: " << lid;
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

  const PetscInt* columns_index_data = columns_index_span.data();
  const PetscReal* matrix_values_data = matrix_values.data();
  Span<const Int32> rows_index_span = m_dof_matrix_numbering.asArray();

  pm->barrier();
  bool is_use_device = false;

  if (runner.isInitialized())
  {
    is_use_device = isAcceleratorPolicy(runner.executionPolicy());
    info() << "[PETSc-Info] Runner for PETSc=" << runner.executionPolicy() << " wanted_is_device=" << is_use_device;
  }

  Real c1 = platform::getRealTime();

  PetscCallAbort(mpi_comm, MatCreate(mpi_comm, &m_petsc_matrix));
  PetscCallAbort(mpi_comm, MatSetSizes(m_petsc_matrix, local_rows, local_rows, global_rows, global_rows));
  PetscCallAbort(mpi_comm, MatSetFromOptions(m_petsc_matrix));

  // m_csr_view.columns() use matrix coordinates local to sub-domain
  // We need to translate them to global matrix coordinates

  PetscCallAbort(mpi_comm, MatSetLocalToGlobalMapping(m_petsc_matrix, m_petsc_map, m_petsc_map));

  if (true)
  {
    // TODO optimize for mutli processes gpu
    // info() << "nb cols: " << csr_view.nbColumn() << "nb rows: " << csr_view.nbRow() << "nb vals: " << csr_view.nbValue();
    std::vector<PetscInt> coo_rows;
    auto csr_rows = csr_view.rows();

    for (int i = 0; i < csr_rows.size() - 1; i++)
    {
      // info() << "csr_rows[ " << i << " ]: " << csr_rows[i];
      int diff = csr_rows[i + 1] - csr_rows[i];

      for (int j = 0; j < diff; j++)
        coo_rows.push_back(i);
    }

    for (int j = 0; j < csr_view.nbValue() - csr_rows[csr_rows.size() - 1]; j++)
      coo_rows.push_back(csr_rows.size() - 1);

    std::vector<PetscInt> coo_cols;
    coo_cols.assign(csr_view.columns().begin(), csr_view.columns().end()); // copy columns array

    PetscCallAbort(mpi_comm, MatSetPreallocationCOOLocal(m_petsc_matrix, csr_view.nbValue(),coo_rows.data(), coo_cols.data()));
    PetscCallAbort(mpi_comm, MatSetValuesCOO(m_petsc_matrix, csr_view.values().data(), INSERT_VALUES));
  }
  else {
    // assemble matrix row by row
    ENUMERATE_ (DoF, idof, dof_family->allItems())
    {
      DoF dof = *idof;

      if (!dof.isOwn()) // for obscure reasons, allItems().own() doesnt work :'(
        continue;

      // info() << idof.index();

      PetscInt nb_col = csr_view.rowsNbColumn()[idof.index()];
      PetscInt row_csr_index = csr_view.rows()[idof.index()];
      PetscInt global_row = rows_index_span[idof.index()];
      // info() << "global_row: " << global_row << " nb_col: " << nb_col << " row_csr_index: " << row_csr_index;
      const PetscInt* cols = &columns_index_data[row_csr_index];
      const PetscScalar* vals = &matrix_values_data[row_csr_index];
      PetscCallAbort(mpi_comm, MatSetValues(m_petsc_matrix, 1, &global_row, nb_col, cols, vals, INSERT_VALUES));
    }
  }

  pm->barrier();

  PetscCallAbort(mpi_comm, MatAssemblyBegin(m_petsc_matrix, MAT_FINAL_ASSEMBLY));
  PetscCallAbort(mpi_comm, MatAssemblyEnd(m_petsc_matrix, MAT_FINAL_ASSEMBLY));

  Real b2 = platform::getRealTime();

  info() << "[PETSc-Timer] Time to create matrix = " << (b2 - c1);

  VariableDoFReal& rhs_variable = this->rhsVariable();
  VariableDoFReal& dof_variable = this->solutionVariable();
  const Real* rhs_data = rhs_variable.asArray().data();
  const Real* result_data = dof_variable.asArray().data();
  const Int32* rows_index_data = rows_index_span.data();

  Real b1 = platform::getRealTime();

  PetscCallAbort(mpi_comm, VecCreateMPI(mpi_comm, local_rows, global_rows, &m_petsc_rhs_vector));
  PetscCallAbort(mpi_comm, VecSetFromOptions(m_petsc_rhs_vector));
  PetscCallAbort(mpi_comm, VecCreateMPI(mpi_comm, local_rows, global_rows, &m_petsc_solution_vector));
  PetscCallAbort(mpi_comm, VecSetFromOptions(m_petsc_solution_vector));

  PetscCallAbort(mpi_comm, VecSetValues(m_petsc_rhs_vector, dof_variable.asArray().size(), rows_index_data, rhs_data, ADD_VALUES));
  PetscCallAbort(mpi_comm, VecSetValues(m_petsc_solution_vector, dof_variable.asArray().size(), rows_index_data, result_data, ADD_VALUES));

  pm->barrier();

  PetscCallAbort(mpi_comm, VecAssemblyBegin(m_petsc_rhs_vector));
  PetscCallAbort(mpi_comm, VecAssemblyEnd(m_petsc_rhs_vector));
  PetscCallAbort(mpi_comm, VecAssemblyBegin(m_petsc_solution_vector));
  PetscCallAbort(mpi_comm, VecAssemblyEnd(m_petsc_solution_vector));

  Real a1 = platform::getRealTime();
  info() << "[PETSc-Timer] Time to create vectors = " << (a1 - b1);

  PetscCallAbort(mpi_comm, KSPCreate(mpi_comm, &m_petsc_solver_context));
  PetscCallAbort(mpi_comm, KSPSetOperators(m_petsc_solver_context, m_petsc_matrix, m_petsc_matrix));
  PetscCallAbort(mpi_comm, KSPSetFromOptions(m_petsc_solver_context));
  PetscCallAbort(mpi_comm, KSPSolve(m_petsc_solver_context, m_petsc_rhs_vector, m_petsc_solution_vector));
  Real a2 = platform::getRealTime();
  info() << "[PETSc-Timer] Time to solve = " << (a2 - a1);

  PetscInt iteration_idx;
  PetscCallAbort(mpi_comm, KSPGetIterationNumber(m_petsc_solver_context, &iteration_idx));

  info() << "[PETSc-Info] Used " << m_pc_type << " preconditionner. Converged in " << iteration_idx + 1 << " iterations";


  if (is_parallel) {
    // Fill 'm_parallel_rows_index' with only rows we owns
    // NOTE: This is only needed if matrix structure has changed.
    Int32 index = 0;
    ENUMERATE_ (DoF, idof, dof_family->allItems()) {
      DoF dof = *idof;
      if (!dof.isOwn())
        continue;
      m_parallel_rows_index[index] = rows_index_span[idof.index()];
      ++index;
    }
  }

  if (is_parallel) {
    Int32 nb_wanted_row = m_parallel_rows_index.extent0();

    PetscCallAbort(mpi_comm, VecGetValues(m_petsc_solution_vector, nb_wanted_row, m_parallel_rows_index.to1DSpan().data(), m_result_work_values.to1DSpan().data()));
    // for (int i = 0; i < nb_wanted_row; i++)
    //   info() << "rows: " << m_parallel_rows_index[i];

    ENUMERATE_ (DoF, idof, dof_family->allItems().own()) {
      Int32 global_idx = rows_index_span[idof.index()];
      // info() << "u[" << global_idx << "] = " << m_result_work_values[idof.index()];
      dof_variable[idof] = m_result_work_values[idof.index()];
    }
  }
  else {
    // TODO check architecture
    const PetscScalar* vals = nullptr;
    PetscCallAbort(mpi_comm, VecGetArrayRead(m_petsc_solution_vector, &vals));

    // Copy directly into Arcane DoF variable
    dof_variable.asArray().copy(Span<const Real>(vals, m_nb_own_row));
    PetscCallAbort(mpi_comm, VecRestoreArrayRead(m_petsc_solution_vector, &vals));
  }

  info() << "[PETSc-Info] Wrote solution in solution_variable";

  PetscCallAbort(mpi_comm, VecDestroy(&m_petsc_solution_vector));
  PetscCallAbort(mpi_comm, VecDestroy(&m_petsc_rhs_vector));
  PetscCallAbort(mpi_comm, MatDestroy(&m_petsc_matrix));
  PetscCallAbort(mpi_comm, KSPDestroy(&m_petsc_solver_context));
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
    x->setPreconditioner(options()->pcType());
    return x;
  }
};


ARCANE_REGISTER_SERVICE_PETSCDOFLINEARSYSTEMFACTORY(PETScLinearSystem,
                                                    PETScDoFLinearSystemFactoryService);

} // namespace Arcane::FemUtils

