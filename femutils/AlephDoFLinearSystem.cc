// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* AlephDoFLinearSystem.cc                                     (C) 2022-2023 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "DoFLinearSystem.h"

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/NumArray.h>

#include <arcane/VariableTypes.h>
#include <arcane/IItemFamily.h>

#include <arcane/aleph/AlephTypesSolver.h>
#include <arcane/aleph/Aleph.h>
#include <arccore/base/NotImplementedException.h>

#include "FemUtils.h"
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

#include <map>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AlephDoFLinearSystemImpl
: public TraceAccessor
, public DoFLinearSystemImpl
{
  static constexpr Byte ELIMINATE_NONE = 0;
  static constexpr Byte ELIMINATE_ROW = 1;
  static constexpr Byte ELIMINATE_ROW_COLUMN = 2;

  struct RowColumn
  {
    Int32 row_id = 0;
    Int32 column_id = 0;
    friend bool operator==(RowColumn rc1, RowColumn rc2)
    {
      if (rc1.row_id != rc2.row_id)
        return false;
      return rc1.column_id == rc2.column_id;
    }
    friend bool operator<(RowColumn rc1, RowColumn rc2)
    {
      if (rc1.row_id == rc2.row_id)
        return rc1.column_id < rc2.column_id;
      return rc1.row_id < rc2.row_id;
    }
  };

  struct RowColumnHash
  {
    size_t operator()(const RowColumn& s) const noexcept
    {
      std::size_t h1 = std::hash<Int32>{}(s.row_id);
      std::size_t h2 = std::hash<Int32>{}(s.column_id);

      return h1 ^ (h2 << 1);
    }
  };

  /*!
   * \brief Map to store values by Row/Column.
   *
   * This map has to be sorted if we want to reuse the internal structure because
   * the matrix filled has to be in the same order when we reuse it.
   */
  using RowColumnMap = std::map<RowColumn, Real>;

 public:

  // TODO: do not use subDomain() but we need to modify aleph before
  AlephDoFLinearSystemImpl(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name)
  : TraceAccessor(sd->traceMng())
  , m_sub_domain(sd)
  , m_dof_family(dof_family)
  , m_rhs_variable(VariableBuildInfo(dof_family, solver_name + "RHSVariable"))
  , m_dof_variable(VariableBuildInfo(dof_family, solver_name + "SolutionVariable"))
  , m_dof_matrix_indexes(VariableBuildInfo(m_dof_family, solver_name + "DoFMatrixIndexes"))
  , m_dof_elimination_info(VariableBuildInfo(m_dof_family, solver_name + "DoFEliminationInfo"))
  , m_dof_elimination_value(VariableBuildInfo(m_dof_family, solver_name + "DoFEliminationValue"))
  {
    info() << "Creating AlephDoFLinearSystemImpl()";
  }

  ~AlephDoFLinearSystemImpl()
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
    m_dof_elimination_info.fill(ELIMINATE_NONE);
    m_dof_elimination_info.fill(0.0);
  }

  AlephParams* params() const { return m_aleph_params; }

  void setSolverBackend(eSolverBackend v) { m_solver_backend = v; }

 private:

  void _computeMatrixInfo()
  {
    int solver_backend = (int)m_solver_backend;
    info() << "[AlephFem] COMPUTE_MATRIX_INFO solver_backend=" << solver_backend;
    // Aleph solver:
    // Hypre = 2
    // Trilinos = 3
    // Cuda = 4 (not available)
    // PETSc = 5
    // We need to compile Arcane with the needed library and link
    // the code with the associated aleph library (see CMakeLists.txt)
    // TODO: Linear algebra backend should be accessed from arc file.
    if (!m_aleph_kernel){
      info() << "Creating Aleph Kernel";
      // We can use less than the number of MPI ranks
      // but for the moment we use all the available cores.
      Int32 nb_core = m_sub_domain->parallelMng()->commSize();
      m_aleph_kernel = new AlephKernel(m_sub_domain, solver_backend, nb_core);
    }
    else{
      //
      m_need_destroy_matrix_and_vector = false;
    }

    DoFGroup own_dofs = m_dof_family->allItems().own();
    //Int32 nb_node = own_nodes.size();
    //Int32 total_nb_node = m_sub_domain->parallelMng()->reduce(Parallel::ReduceSum, nb_node);
    m_dof_matrix_indexes.fill(-1);
    AlephIndexing* indexing = m_aleph_kernel->indexing();
    ENUMERATE_ (DoF, idof, own_dofs) {
      DoF dof = *idof;
      Integer row = indexing->get(m_dof_variable, dof);
      //info() << "ROW=" << row;
      m_dof_matrix_indexes[dof] = row;
    }
    // Do not print informations about setting matrix if matrix is too big
    if (own_dofs.size()>200)
      m_do_print_filling = false;

    m_aleph_matrix = m_aleph_kernel->createSolverMatrix();
    m_aleph_rhs_vector = m_aleph_kernel->createSolverVector();
    m_aleph_solution_vector = m_aleph_kernel->createSolverVector();
    m_aleph_matrix->create();
    m_aleph_rhs_vector->create();
    m_aleph_solution_vector->create();
  }

 public:

  void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (column.isNull())
      ARCANE_FATAL("Column is null");
    if (value == 0.0)
      return;
    if (m_use_value_map) {
      RowColumn rc{ row.localId(), column.localId() };
      auto x = m_values_map.find(rc);
      if (x == m_values_map.end())
        m_values_map.insert(std::make_pair(rc, value));
      else
        x->second += value;
    }
    else {
      ItemInfoListView item_list_view(m_dof_family);
      info() << "AlephAdd R=" << row.localId() << " C=" << column.localId() << " V=" << value;
      m_aleph_matrix->addValue(m_dof_variable, item_list_view[row], m_dof_variable, item_list_view[column], value);
    }
  }

  void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (column.isNull())
      ARCANE_FATAL("Column is null");
    if (!m_use_value_map)
      ARCANE_FATAL("matrixSetValue() is only allowed if 'm_use_value_map' is true");
    m_forced_set_values_map[{ row.localId(), column.localId() }] = value;
  }

  void eliminateRow(DoFLocalId row, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (!m_use_value_map)
      ARCANE_FATAL("matrixEliminateRow() is only allowed if 'm_use_value_map' is true");
    m_dof_elimination_info[row] = ELIMINATE_ROW;
    m_dof_elimination_value[row] = value;
    info() << "EliminateRow row=" << row.localId() << " v=" << value;
  }

  void eliminateRowColumn(DoFLocalId row, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (!m_use_value_map)
      ARCANE_FATAL("matrixEliminateRowColumn() is only allowed if 'm_use_value_map' is true");
    m_dof_elimination_info[row] = ELIMINATE_ROW_COLUMN;
    m_dof_elimination_value[row] = value;
    info() << "EliminateRowColumn row=" << row.localId() << " v=" << value;
  }

  void setEliminationArrays(VariableDoFByte& dof_elimination_info, VariableDoFReal& dof_elimination_value) override
  {
    ARCANE_THROW(NotImplementedException, "");
  };

  void solve() override
  {
    UniqueArray<Real> aleph_result;

    // _fillMatrix() may change the values of RHS vector
    // with row or row-column elimination so we has to fill the RHS vector
    // before the matrix.
    _fillMatrix();
    _fillRHSVector();

    info() << "[AlephFem] Assemble matrix ptr=" << m_aleph_matrix;
    m_aleph_matrix->assemble();
    m_aleph_rhs_vector->assemble();
    auto* aleph_solution_vector = m_aleph_solution_vector;
    DoFGroup own_dofs = m_dof_family->allItems().own();
    const Int32 nb_dof = own_dofs.size();
    m_vector_zero.resize(nb_dof);
    m_vector_zero.fill(0.0);

    aleph_solution_vector->setLocalComponents(m_vector_zero);
    aleph_solution_vector->assemble();

    Int32 nb_iteration = 0;
    Real residual_norm = 0.0;
    info() << "[AlephFem] BEGIN SOLVING WITH ALEPH solver_backend=" << (int)m_solver_backend;

    m_aleph_matrix->solve(aleph_solution_vector,
                          m_aleph_rhs_vector,
                          nb_iteration,
                          &residual_norm,
                          m_aleph_params,
                          false);
    info() << "[AlephFem] END SOLVING WITH ALEPH r=" << residual_norm
           << " nb_iter=" << nb_iteration;
    auto* rhs_vector = m_aleph_kernel->createSolverVector();
    auto* solution_vector = m_aleph_kernel->createSolverVector();

    UniqueArray<Real> rhs_results;
    rhs_vector->getLocalComponents(rhs_results);
    solution_vector->getLocalComponents(aleph_result);

    bool do_verbose = nb_dof < 200;
    do_verbose = false;
    Int32 index = 0;
    ENUMERATE_ (DoF, idof, m_dof_family->allItems().own()) {
      DoF dof = *idof;

      m_dof_variable[dof] = aleph_result[m_aleph_kernel->indexing()->get(m_dof_variable, dof)];
      if (do_verbose)
        info() << "Node uid=" << dof.uniqueId() << " V=" << aleph_result[index] << " T=" << m_dof_variable[dof]
               << " RHS=" << rhs_results[index];

      ++index;
    }
  }

  VariableDoFReal& solutionVariable() override
  {
    return m_dof_variable;
  }

  VariableDoFReal& rhsVariable() override
  {
    return m_rhs_variable;
  }

  void setSolverCommandLineArguments(const CommandLineArguments& args) override
  {
#if ARCANE_VERSION >= 31002
    m_aleph_kernel->solverInitializeArgs().setCommandLineArguments(args);
#else
    pwarning() << "Call to setSolverCommandLineArguments() is not used because version of Arcane is too old (3.10.2+ required)";
#endif
  }

  void clearValues()
  {
    info() << "[Aleph] Clear values of current solver";
    m_dof_elimination_info.fill(ELIMINATE_NONE);
    m_dof_elimination_info.fill(0.0);
    m_values_map.clear();
    m_forced_set_values_map.clear();
    _computeMatrixInfo();
  }

  void setCSRValues(const CSRFormatView& csr_view) override
  {
    ARCANE_THROW(NotImplementedException,"");
  }

  bool hasSetCSRValues() const { return false; }
  void setRunner(Runner* r) override { m_runner = r; }
  Runner* runner() const { return m_runner; }

 private:

  AlephParams* _createAlephParam()
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

 private:

  ISubDomain* m_sub_domain = nullptr;
  IItemFamily* m_dof_family = nullptr;
  VariableDoFReal m_rhs_variable;
  VariableDoFReal m_dof_variable;
  VariableDoFInt32 m_dof_matrix_indexes;
  VariableDoFByte m_dof_elimination_info;
  VariableDoFReal m_dof_elimination_value;
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
  /*!
   * \brief True is we use 'm_values_map' to mix add and set.
   *
   * You may set the value to 'false' if you want the old behavior when
   * there is only matrixAddValue() calls.
   */
  bool m_use_value_map = true;

  //! True to print matrix values during filling
  bool m_do_print_filling = true;

  //! True is we need to manually destroy the matrix/vector
  bool m_need_destroy_matrix_and_vector = true;

  UniqueArray<Real> m_vector_zero;

  Runner* m_runner = nullptr;

 private:

  void _fillMatrix();
  void _fillRHSVector();
  void _setMatrixValue(DoF row, DoF column, Real value)
  {
    if (m_do_print_filling)
      info() << "SET MATRIX VALUE (" << std::setw(4) << row.localId()
             << "," << std::setw(4) << column.localId() << ")"
             << " v=" << std::setw(25) << value;
    m_aleph_matrix->setValue(m_dof_variable, row, m_dof_variable, column, value);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AlephDoFLinearSystemFactoryService
: public ArcaneAlephDoFLinearSystemFactoryObject
{
 public:

  AlephDoFLinearSystemFactoryService(const ServiceBuildInfo& sbi)
  : ArcaneAlephDoFLinearSystemFactoryObject(sbi)
  {
  }

  DoFLinearSystemImpl*
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

extern "C++" DoFLinearSystemImpl*
createAlephDoFLinearSystemImpl(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name)
{
  auto* x = new AlephDoFLinearSystemImpl(sd, dof_family, solver_name);
  x->build();
  return x;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
_fillMatrix()
{
  // Fill the matrix
  if (!m_use_value_map)
    return;

  RowColumnMap row_column_elimination_map;

  DoFInfoListView item_list_view(m_dof_family);
  for (const auto& rc_value : m_values_map) {
    RowColumn rc = rc_value.first;
    Real value = rc_value.second;
    DoF dof_row = item_list_view[rc.row_id];
    DoF dof_column = item_list_view[rc.column_id];

    Byte row_elimination_info = m_dof_elimination_info[dof_row];
    Byte column_elimination_info = m_dof_elimination_info[dof_column];

    if (row_elimination_info == ELIMINATE_ROW_COLUMN || column_elimination_info == ELIMINATE_ROW_COLUMN) {
      row_column_elimination_map[{ rc.row_id, rc.column_id }] = value;
      continue;
    }

    if (row_elimination_info == ELIMINATE_ROW)
      // Will be computed after this loop
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

  // Apply Row+Column elimination
  // Phase 1:
  // - substract values of the RHS vector if Row+Column elimination
  for (const auto& rc_value : row_column_elimination_map) {
    RowColumn rc = rc_value.first;
    Real matrix_value = rc_value.second;
    DoF dof_row = item_list_view[rc.row_id];
    DoF dof_column = item_list_view[rc.column_id];
    if (dof_row == dof_column)
      continue;
    if (!dof_column.isOwn())
      continue;
    Byte row_elimination_info = m_dof_elimination_info[dof_row];
    Real elimination_value = m_dof_elimination_value[dof_row];
    // Substract the value of RHS vector for current column.
    if (row_elimination_info == ELIMINATE_ROW_COLUMN) {
      Real v = m_rhs_variable[dof_column];
      m_rhs_variable[dof_column] = v - matrix_value * elimination_value;
      if (m_do_print_filling)
        info() << "EliminateRowColumn (" << std::setw(4) << rc.row_id
               << "," << std::setw(4) << rc.column_id << ")"
               << " elimination_value=" << std::setw(25) << elimination_value
               << "  old_rhs=" << std::setw(25) << v
               << "  new_rhs=" << std::setw(25) << m_rhs_variable[dof_column];
    }
  }

  // Apply Row or Row+Column elimination
  // Phase 2: set the value of the RHS if Row elimination
  // Phase 2: fill the diagonal with 1.0
  ENUMERATE_ (DoF, idof, m_dof_family->allItems()) {
    DoF dof = *idof;
    if (!dof.isOwn())
      continue;
    Byte elimination_info = m_dof_elimination_info[dof];
    if (elimination_info == ELIMINATE_ROW || elimination_info == ELIMINATE_ROW_COLUMN) {
      Real elimination_value = m_dof_elimination_value[dof];
      m_rhs_variable[dof] = elimination_value;
      info() << "Eliminate info=" << (int)elimination_info << " row="
             << std::setw(4) << dof.localId() << " value=" << elimination_value;
      _setMatrixValue(dof, dof, 1.0);
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlephDoFLinearSystemImpl::
_fillRHSVector()
{
  // For the LinearSystem class we need an array
  // with only the values for the ownNodes().
  // The values of 'rhs_values' should not be updated after
  // this call.
  UniqueArray<Real> rhs_values_for_linear_system;
  VariableDoFReal& rhs_values(rhsVariable());
  ENUMERATE_ (DoF, idof, m_dof_family->allItems().own()) {
    Real v = rhs_values[idof];
    if (m_do_print_filling)
      info() << "SET VECTOR VALUE (" << std::setw(4) << idof.itemLocalId() << ") = " << v;
    rhs_values_for_linear_system.add(rhs_values[idof]);
  }
  m_aleph_rhs_vector->setLocalComponents(rhs_values_for_linear_system.view());
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
