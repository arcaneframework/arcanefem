// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* AlephDoFLinearSystem.cc                                     (C) 2000-2026 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "DoFLinearSystem.h"

#include <arcane/utils/FatalErrorException.h>

#include <arcane/core/VariableTypes.h>
#include <arcane/core/IItemFamily.h>

#include <arcane/accelerator/core/Runner.h>

#include <arcane/aleph/AlephTypesSolver.h>
#include <arcane/aleph/Aleph.h>

#ifdef FEMUTILS_HAS_PETSC
#include <petsclog.h>
#endif

#include "FemUtils.h"
#include "internal/DoKDoFLinearSystemImpl.h"
#include "IDoFLinearSystemFactory.h"
#include "CsrFormatMatrixView.h"

namespace Arcane::FemUtils
{
enum class eSolverBackend
{
  Hypre = 2,
  Trilinos = 3,
  Petsc = 5,
};
}

#include "AlephDoFLinearSystemFactory_axl.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AlephDoFLinearSystemImpl
: public DoKDoFLinearSystemImpl
{
  using RowColumn = DoKMatrix::RowColumn;

 public:

  // TODO: do not use subDomain() but we need to modify aleph before
  AlephDoFLinearSystemImpl(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name)
  : DoKDoFLinearSystemImpl(dof_family, solver_name)
  , m_sub_domain(sd)
  , m_dof_matrix_indexes(VariableBuildInfo(dof_family, solver_name + "DoFMatrixIndexes"))
  {
    info() << "Creating AlephDoFLinearSystemImpl()";
  }

  ~AlephDoFLinearSystemImpl() override
  {
#ifdef FEMUTILS_HAS_PETSC
    // Aleph initializes PETSc, but currently does not call PetscFinalize().
    // PETSc normally prints -log_view during PetscFinalize(), so explicitly
    // process the log-view option while Aleph's PETSc communicator is valid.
    if (m_solver_backend == eSolverBackend::Petsc) {
      PetscBool is_initialized = PETSC_FALSE;
      PetscBool is_finalized = PETSC_FALSE;

      // Fetch PETSc initialization and finalization status
      PetscCallAbort(PETSC_COMM_WORLD, PetscInitialized(&is_initialized));
      PetscCallAbort(PETSC_COMM_WORLD, PetscFinalized(&is_finalized));
      if (is_initialized && !is_finalized)
        PetscCallAbort(PETSC_COMM_WORLD, PetscLogViewFromOptions());
    }
#endif
    delete m_aleph_params;
    if (m_need_destroy_matrix_and_vector) {
      delete m_aleph_matrix;
      delete m_aleph_rhs_vector;
      delete m_aleph_solution_vector;
    }
    // Doit être fait dans Arcane
    // delete m_aleph_kernel->factory();
    delete m_aleph_kernel;
  }

 public:

  void build()
  {
    _computeMatrixInfo();
    m_aleph_params = _createAlephParam();
    DoKDoFLinearSystemImpl::clearValues();
  }

  AlephParams* params() const { return m_aleph_params; }

  void setSolverBackend(eSolverBackend v) { m_solver_backend = v; }

 private:

  void _computeMatrixInfo();

 public:

  void applyMatrixTransformation() override;
  void solve() override;

  void setSolverCommandLineArguments(const CommandLineArguments& args) override
  {
    m_aleph_kernel->solverInitializeArgs().setCommandLineArguments(args);
  }

  void clearValues() override
  {
    info() << "[Aleph] Clear values of current solver";
    DoKDoFLinearSystemImpl::clearValues();
    _computeMatrixInfo();
  }

 private:

  ISubDomain* m_sub_domain = nullptr;
  VariableDoFInt32 m_dof_matrix_indexes;
  AlephKernel* m_aleph_kernel = nullptr;
  AlephMatrix* m_aleph_matrix = nullptr;
  AlephVector* m_aleph_rhs_vector = nullptr;
  AlephVector* m_aleph_solution_vector = nullptr;
  AlephParams* m_aleph_params = nullptr;
  eSolverBackend m_solver_backend = eSolverBackend::Hypre;

  //! True to print matrix values during filling
  bool m_do_print_filling = true;

  //! True is we need to manually destroy the matrix/vector
  bool m_need_destroy_matrix_and_vector = true;

 private:

  AlephParams* _createAlephParam() const;
  void _applyMatrixTransformationAndFillAlephMatrix();
  void _fillRHSVector();
  void _fillSolutionVector();
  void _setMatrixValue(DoF row, DoF column, Real value)
  {
    if (m_do_print_filling)
      info() << "SET MATRIX VALUE (" << std::setw(4) << row.localId()
             << "," << std::setw(4) << column.localId() << ")"
             << " v=" << std::setw(25) << value;
    VariableDoFReal& solution_variable = solutionVariable();
    m_aleph_matrix->setValue(solution_variable, row, solution_variable, column, value);
  }
  void _createRHSAndSolutionVector();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AlephDoFLinearSystemFactoryService
: public ArcaneAlephDoFLinearSystemFactoryObject
{
 public:

  explicit AlephDoFLinearSystemFactoryService(const ServiceBuildInfo& sbi)
  : ArcaneAlephDoFLinearSystemFactoryObject(sbi)
  {
  }

  IDoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name) override
  {
    auto* x = new AlephDoFLinearSystemImpl(sd, dof_family, solver_name);
    x->setSolverBackend(options()->solverBackend());

    x->build();

    auto* p = x->params();
    p->setEpsilon(options()->epsilon());
    p->setPrecond(options()->preconditioner());
    p->setMethod(options()->solverMethod());

    return x;
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

extern "C++" IDoFLinearSystemImpl*
createAlephDoFLinearSystemImpl(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name)
{
  auto* x = new AlephDoFLinearSystemImpl(sd, dof_family, solver_name);
  x->build();
  return x;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
_applyMatrixTransformationAndFillAlephMatrix()
{
  fillRowColumnEliminationInfos();
  // We provide two ways to fill the Aleph matrix.
  // The first one (currently the default) fill the matrix using the DoK Matrix.
  // The second one converts the DoKMatrix to a CSR Matrix then fill the Aleph Matrix.
  // The second one will be used if we want to directly call other linear solver
  // like PETSc or Hypre when using a DoK Matrix.
  bool do_with_csr = false;
  if (do_with_csr) {
    convertToCSRMatrix();
    CsrFormatMatrixView csr_view = getCsrFormatMatrixView();
    Int32 nb_row = csr_view.nbRow();

    IItemFamily* dof_family = dofFamily();
    DoFInfoListView item_list_view(dof_family);

    // Fill the Aleph Matrix
    for (Int32 row_id = 0; row_id < nb_row; ++row_id) {
      for (CsrRowColumnIndex rc : csr_view.rowRange(row_id)){
        Int32 column_id = csr_view.column(rc);
        Real value = csr_view.value(rc);
        //info() << "ROW_ID=" << row_id << " column=" << column_id << " value=" << value;
        _setMatrixValue(item_list_view[row_id], item_list_view[column_id], value);
      }
    }
  }
  else {
    auto set_matrix_value = [&](DoF row, DoF column, Real value) {
      _setMatrixValue(row, column, value);
    };
    visitDoKMatrix(set_matrix_value);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
_createRHSAndSolutionVector()
{
  // We need to call createSolverVector() two times.
  // The first time returns the RHS and the second the solution vector
  m_aleph_rhs_vector = m_aleph_kernel->createSolverVector();
  m_aleph_solution_vector = m_aleph_kernel->createSolverVector();

  m_aleph_rhs_vector->create();
  m_aleph_solution_vector->create();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
_fillRHSVector()
{
  ARCANE_CHECK_POINTER(m_aleph_rhs_vector);

  // For the LinearSystem class we need an array
  // with only the values for the ownNodes().
  // The values of 'rhs_values' should not be updated after
  // this call.
  UniqueArray<Real> rhs_values_for_linear_system;
  VariableDoFReal& rhs_values(rhsVariable());
  IItemFamily* dof_family = dofFamily();
  ENUMERATE_ (DoF, idof, dof_family->allItems().own()) {
    Real v = rhs_values[idof];
    if (m_do_print_filling)
      info() << "SET VECTOR VALUE (" << std::setw(4) << idof.itemLocalId() << ") = " << v;
    rhs_values_for_linear_system.add(rhs_values[idof]);
  }

  m_aleph_rhs_vector->setLocalComponents(rhs_values_for_linear_system.view());
  m_aleph_rhs_vector->assemble();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
_fillSolutionVector()
{
  ARCANE_CHECK_POINTER(m_aleph_solution_vector);

  UniqueArray<Real> solution_values_for_linear_system;
  VariableDoFReal& solution_values(solutionVariable());
  IItemFamily* dof_family = dofFamily();
  ENUMERATE_ (DoF, idof, dof_family->allItems().own()) {
    Real v = solution_values[idof];
    if (m_do_print_filling)
      info() << "SET SOLUTION VECTOR VALUE (" << std::setw(4) << idof.itemLocalId() << ") = " << v;
    solution_values_for_linear_system.add(solution_values[idof]);
  }

  m_aleph_solution_vector->setLocalComponents(solution_values_for_linear_system.view());
  m_aleph_solution_vector->assemble();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

AlephParams* AlephDoFLinearSystemImpl::
_createAlephParam() const
{
  auto* p = new AlephParams(traceMng(),
                            1.0e-15, // m_param_epsilon epsilon de convergence
                            2000, // m_param_max_iteration nb max iterations
                            TypesSolver::AMG, // m_param_preconditioner_method préconditionnement: DIAGONAL, AMG, IC
                            TypesSolver::PCG, // m_param_solver_method méthode de résolution
                            -1, // m_param_gamma
                            -1.0, // m_param_alpha
                            false, // m_param_xo_user par défaut Xo n'est pas égal à 0
                            false, // m_param_check_real_residue
                            false, // m_param_print_real_residue
                            // Default: false
                            false, // m_param_debug_info
                            1.e-40, // m_param_min_rhs_norm
                            false, // m_param_convergence_analyse
                            true, // m_param_stop_error_strategy
                            true, // m_param_write_matrix_to_file_error_strategy
                            "SolveErrorAlephMatrix.dbg", // m_param_write_matrix_name_error_strategy
                            false, // m_param_listing_output
                            0., // m_param_threshold
                            false, // m_param_print_cpu_time_resolution
                            0, // m_param_amg_coarsening_method: par défault celui de Sloop,
                            100, // m_param_output_level
                            1, // m_param_amg_cycle: 1-cycle amg en V, 2= cycle amg en W, 3=cycle en Full Multigrid V
                            1, // m_param_amg_solver_iterations
                            1, // m_param_amg_smoother_iterations
                            TypesSolver::SymHybGSJ_smoother, // m_param_amg_smootherOption
                            TypesSolver::ParallelRugeStuben, // m_param_amg_coarseningOption
                            TypesSolver::CG_coarse_solver, // m_param_amg_coarseSolverOption
                            // Default: false
                            true, // m_param_keep_solver_structure
                            false, // m_param_sequential_solver
                            TypesSolver::RB); // m_param_criteria_stop
  return p;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
_computeMatrixInfo()
{
  int solver_backend = static_cast<int>(m_solver_backend);
  info() << "[AlephFem] COMPUTE_MATRIX_INFO solver_backend=" << solver_backend;
  IParallelMng* pm = m_sub_domain->parallelMng();
  // Aleph solver:
  // Hypre = 2
  // Trilinos = 3
  // Cuda = 4 (not available)
  // Petsc = 5
  // We need to compile Arcane with the needed library and link
  // the code with the associated aleph library (see CMakeLists.txt)
  // TODO: Linear algebra backend should be accessed from arc file.
  if (!m_aleph_kernel) {
    info() << "Creating Aleph Kernel";
    // We can use less than the number of MPI ranks
    // but for the moment we use all the available cores.
    Int32 nb_core = pm->commSize();
    m_aleph_kernel = new AlephKernel(m_sub_domain, solver_backend, nb_core);
  }
  else {
    //
    m_need_destroy_matrix_and_vector = false;
  }
  IItemFamily* dof_family = dofFamily();
  VariableDoFReal& solution_variable(solutionVariable());
  DoFGroup own_dofs = dof_family->allItems().own();
  m_dof_matrix_indexes.fill(-1);
  AlephIndexing* indexing = m_aleph_kernel->indexing();
  ENUMERATE_ (DoF, idof, own_dofs) {
    DoF dof = *idof;
    Integer row = indexing->get(solution_variable, dof);
    m_dof_matrix_indexes[dof] = row;
  }

  // Do not print information about setting matrix if matrix is too big
  if (own_dofs.size() > 200) {
    m_do_print_filling = false;
    setPrintFilling(!m_do_print_filling);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
applyMatrixTransformation()
{
  info() << "[AlephFem] Assemble matrix ptr=" << m_aleph_matrix;
  // Check this is called only one time
  if (m_aleph_matrix)
    ARCANE_FATAL("applyMatrixTransformation() has already been called");
  m_aleph_matrix = m_aleph_kernel->createSolverMatrix();
  m_aleph_matrix->create();

  // Matrix transformation
  _applyMatrixTransformationAndFillAlephMatrix();
  m_aleph_matrix->assemble();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
solve()
{
  // The matrix always need to be created if this is not explicitly done
  if (!m_aleph_matrix)
    applyMatrixTransformation();

  _createRHSAndSolutionVector();
  _fillRHSVector();
  _fillSolutionVector();

  info() << "Calling AlephDoFLinearSystemImpl::solve()";
  UniqueArray<Real> aleph_result;

  IItemFamily* dof_family = dofFamily();
  DoFGroup own_dofs = dof_family->allItems().own();
  const Int32 nb_dof = own_dofs.size();

  Int32 nb_iteration = 0;
  Real residual_norm = 0.0;
  info() << "[AlephFem] BEGIN SOLVING WITH ALEPH solver_backend=" << static_cast<int>(m_solver_backend);

  // Post the solver. The call is asynchronous, and we wait for the result
  // when calling syncSolver().
  // The values nb_iteration and residual_norm are not used in this case.
  // We get them during the call to syncSolver()
  m_aleph_matrix->solve(m_aleph_solution_vector, m_aleph_rhs_vector,
                        nb_iteration, &residual_norm,
                        m_aleph_params, true);

  // Reset matrix and vectors because there can no longer be used
  // They will be re-created when needed if we call solve() again.
  // NOTE: it is the aleph library we do not need to call delete() because
  m_aleph_rhs_vector = nullptr;
  m_aleph_solution_vector = nullptr;
  m_aleph_matrix = nullptr;

  info() << "[AlephFem] END SOLVING WITH ALEPH r=" << residual_norm
         << " nb_iter=" << nb_iteration;

  // Wait for the solver to finish and get solution vector
  auto* solution_vector = m_aleph_kernel->syncSolver(0, nb_iteration, &residual_norm);

  solution_vector->getLocalComponents(aleph_result);

  const bool do_verbose = (nb_dof < 200);
  Int32 index = 0;

  VariableDoFReal& solution_variable(this->solutionVariable());
  ENUMERATE_ (DoF, idof, dofFamily()->allItems().own()) {
    DoF dof = *idof;

    solution_variable[dof] = aleph_result[m_aleph_kernel->indexing()->get(solution_variable, dof)];
    if (do_verbose)
      info() << "Node uid=" << dof.uniqueId() << " V=" << aleph_result[index] << " T=" << solution_variable[dof];
    ++index;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ALEPHDOFLINEARSYSTEMFACTORY(AlephLinearSystem,
                                                    AlephDoFLinearSystemFactoryService);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
