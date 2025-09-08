// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* AlephDoFLinearSystem.cc                                     (C) 2022-2025 */
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

#include "FemUtils.h"
#include "internal/DoFLinearSystemImplBase.h"
#include "IDoFLinearSystemFactory.h"
#include "arcane_version.h"

namespace Arcane::FemUtils
{
enum class eSolverBackend
{
  Hypre = 2,
  Trilinos = 3,
  PETSc = 5,
};
}

#include "AlephDoFLinearSystemFactory_axl.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class AlephDoFLinearSystemImpl
: public DoFLinearSystemImplBase
{

  /*!
   * \brief Map to store values by Row/Column.
   *
   * This map has to be sorted if we want to reuse the internal structure because
   * the matrix filled has to be in the same order when we reuse it.
   */
  using RowColumnMap = OrderedRowColumnMap;
  using RowColumn = RowColumnMap::RowColumn;

 public:

  // TODO: do not use subDomain() but we need to modify aleph before
  AlephDoFLinearSystemImpl(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name)
  : DoFLinearSystemImplBase(dof_family, solver_name)
  , m_sub_domain(sd)
  , m_dof_matrix_indexes(VariableBuildInfo(dof_family, solver_name + "DoFMatrixIndexes"))
  {
    info() << "Creating AlephDoFLinearSystemImpl()";
  }

  ~AlephDoFLinearSystemImpl() override
  {
    delete m_aleph_params;
    if (m_need_destroy_matrix_and_vector){
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
    DoFLinearSystemImplBase::clearValues();
  }

  AlephParams* params() const { return m_aleph_params; }

  void setSolverBackend(eSolverBackend v) { m_solver_backend = v; }

 private:

  void _computeMatrixInfo();

 public:

  void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (column.isNull())
      ARCANE_FATAL("Column is null");
    if (value == 0.0)
      return;
    RowColumn rc{ row.localId(), column.localId() };
    m_values_map.addValue(rc, value);
  }

  void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (column.isNull())
      ARCANE_FATAL("Column is null");
    m_forced_set_values_map[{ row.localId(), column.localId() }] = value;
  }

  void eliminateRow(DoFLocalId row, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    getEliminationInfo()[row] = ELIMINATE_ROW;
    getEliminationValue()[row] = value;
    info() << "EliminateRow row=" << row.localId() << " v=" << value;
  }

  void eliminateRowColumn(DoFLocalId row, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    getEliminationInfo()[row] = ELIMINATE_ROW_COLUMN;
    getEliminationValue()[row] = value;
    info() << "EliminateRowColumn row=" << row.localId() << " v=" << value;
  }

  void applyMatrixTransformation() override;
  void applyRHSTransformation() override;
  void solve() override;

  void setSolverCommandLineArguments(const CommandLineArguments& args) override
  {
#if ARCANE_VERSION >= 31002
    m_aleph_kernel->solverInitializeArgs().setCommandLineArguments(args);
#else
    pwarning() << "Call to setSolverCommandLineArguments() is not used because version of Arcane is too old (3.10.2+ required)";
#endif
  }

  void clearValues() override
  {
    info() << "[Aleph] Clear values of current solver";
    DoFLinearSystemImplBase::clearValues();
    m_values_map.clear();
    m_forced_set_values_map.clear();
    _computeMatrixInfo();
  }

  void setCSRValues(const CSRFormatView& csr_view) override
  {
    ARCANE_THROW(NotImplementedException,"");
  }
  CSRFormatView& getCSRValues() override
  {
    ARCANE_THROW(NotImplementedException, "");
  }
  bool hasSetCSRValues() const override { return false; }

 private:

  ISubDomain* m_sub_domain = nullptr;
  VariableDoFInt32 m_dof_matrix_indexes;
  AlephKernel* m_aleph_kernel = nullptr;
  AlephMatrix* m_aleph_matrix = nullptr;
  AlephVector* m_aleph_rhs_vector = nullptr;
  AlephVector* m_aleph_solution_vector = nullptr;
  AlephParams* m_aleph_params = nullptr;
  eSolverBackend m_solver_backend = eSolverBackend::Hypre;
  //! List of (i,j) values added to the matrix
  RowColumnMap m_values_map;
  //! List of (i,j) whose value is fixed. This will override added values in m_values_map.
  RowColumnMap m_forced_set_values_map;

  //! True to print matrix values during filling
  bool m_do_print_filling = true;

  //! True is we need to manually destroy the matrix/vector
  bool m_need_destroy_matrix_and_vector = true;

  UniqueArray<Real> m_vector_zero;

 private:

  AlephParams* _createAlephParam() const;
  void _applyMatrixTransformationAndFillAlephMatrix();
  void _fillRHSVector();
  void _applyRHSTransformation();
  void _setMatrixValue(DoF row, DoF column, Real value)
  {
    if (m_do_print_filling)
      info() << "SET MATRIX VALUE (" << std::setw(4) << row.localId()
             << "," << std::setw(4) << column.localId() << ")"
             << " v=" << std::setw(25) << value;
    VariableDoFReal& solution_variable = solutionVariable();
    m_aleph_matrix->setValue(solution_variable, row, solution_variable, column, value);
  }
  void _fillRowColumnEliminationInfos();
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

void AlephDoFLinearSystemImpl::
_fillRowColumnEliminationInfos()
{
  OrderedRowColumnMap& rc_elimination_map = _rowColumnEliminationMap();
  rc_elimination_map.clear();
  DoFInfoListView item_list_view(dofFamily());

  auto& dof_elimination_info = getEliminationInfo();
  auto& dof_elimination_value = getEliminationValue();

  for (const auto& rc_value : m_values_map) {
    RowColumn rc = rc_value.first;
    Real value = rc_value.second;
    DoF dof_row = item_list_view[rc.row_id];
    DoF dof_column = item_list_view[rc.column_id];
    Byte row_elimination_info = dof_elimination_info[dof_row];
    Byte column_elimination_info = dof_elimination_info[dof_column];
    if (row_elimination_info == ELIMINATE_ROW_COLUMN || column_elimination_info == ELIMINATE_ROW_COLUMN)
      rc_elimination_map[{ rc.row_id, rc.column_id }] = value;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
_applyMatrixTransformationAndFillAlephMatrix()
{
  _fillRowColumnEliminationInfos();
  OrderedRowColumnMap& rc_elimination_map = _rowColumnEliminationMap();

  IItemFamily* dof_family = dofFamily();
  DoFInfoListView item_list_view(dof_family);

  auto& dof_elimination_info = getEliminationInfo();
  auto& dof_elimination_value = getEliminationValue();

  // Fill the matrix from the values of \a m_values_map
  // Skip (row,column) values which are part of an elimination.
  for (const auto& rc_value : m_values_map) {
    RowColumn rc = rc_value.first;
    Real value = rc_value.second;
    DoF dof_row = item_list_view[rc.row_id];
    DoF dof_column = item_list_view[rc.column_id];

    Byte row_elimination_info = dof_elimination_info[dof_row];

    if (row_elimination_info == ELIMINATE_ROW)
      // Will be computed after this loop
      continue;
    if (rc_elimination_map.contains(rc))
      continue;

    // Check if value is forced for current RowColumn
    auto x = m_forced_set_values_map.find(rc);
    if (x != m_forced_set_values_map.end()) {
      info(4) << "FORCED VALUE R=" << rc.row_id << " C=" << rc.column_id
              << " old=" << value << " new=" << x->second;
      value = x->second;
    }

    _setMatrixValue(dof_row, dof_column, value);
  }

  const bool do_print_filling = m_do_print_filling;

  // Apply Row or Row+Column elimination on Matrix
  // Phase 2: set the diagonal value for elimination row to 1.0
  ENUMERATE_ (DoF, idof, dof_family->allItems()) {
    DoF dof = *idof;
    if (!dof.isOwn())
      continue;
    Byte elimination_info = dof_elimination_info[dof];
    if (elimination_info == ELIMINATE_ROW || elimination_info == ELIMINATE_ROW_COLUMN) {
      Real elimination_value = dof_elimination_value[dof];
      if (do_print_filling)
        info() << "EliminateMatrix info=" << static_cast<int>(elimination_info) << " row="
               << std::setw(4) << dof.localId() << " value=" << elimination_value;
      _setMatrixValue(dof, dof, 1.0);
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
_applyRHSTransformation()
{
  const bool do_print_filling = m_do_print_filling;

  // Apply Row+Column elimination
  // Phase 1:
  // - subtract values of the RHS vector if Row+Column elimination
  _applyRowColumnEliminationToRHS(do_print_filling);

  IItemFamily* dof_family = dofFamily();

  auto& dof_elimination_info = getEliminationInfo();
  auto& dof_elimination_value = getEliminationValue();
  auto& rhs_variable = rhsVariable();

  // Apply Row or Row+Column elimination on RHS
  ENUMERATE_ (DoF, idof, dof_family->allItems()) {
    DoF dof = *idof;
    if (!dof.isOwn())
      continue;
    Byte elimination_info = dof_elimination_info[dof];
    if (elimination_info == ELIMINATE_ROW || elimination_info == ELIMINATE_ROW_COLUMN) {
      Real elimination_value = dof_elimination_value[dof];
      rhs_variable[dof] = elimination_value;
      if (do_print_filling)
        info() << "EliminateRHS info=" << static_cast<int>(elimination_info) << " row="
               << std::setw(4) << dof.localId() << " value=" << elimination_value;
    }
  }
}

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
  // PETSc = 5
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
  if (own_dofs.size() > 200)
    m_do_print_filling = false;
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
applyRHSTransformation()
{
  // RHS Transformation
  _applyRHSTransformation();
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

  info() << "Calling AlephDoFLinearSystemImpl::solve()";
  UniqueArray<Real> aleph_result;

  IItemFamily* dof_family = dofFamily();
  DoFGroup own_dofs = dof_family->allItems().own();
  const Int32 nb_dof = own_dofs.size();
  m_vector_zero.resize(nb_dof);
  m_vector_zero.fill(0.0);

  m_aleph_solution_vector->setLocalComponents(m_vector_zero);
  m_aleph_solution_vector->assemble();

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
