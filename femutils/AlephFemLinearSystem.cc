// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* NodeLinearSystem.cc                                          (C) 2022-2022 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "NodeLinearSystem.h"

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/NumArray.h>

#include <arcane/VariableTypes.h>
#include <arcane/IItemFamily.h>

#include <arcane/aleph/AlephTypesSolver.h>
#include <arcane/aleph/Aleph.h>

#include "FemUtils.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AlephFemLinearSystemImpl
: public TraceAccessor
, public NodeLinearSystemImpl
{
 public:

  // TODO: do not use subDomain() but we need to modify aleph before
  AlephFemLinearSystemImpl(ISubDomain* sd, const Arcane::VariableNodeReal& node_variable)
  : TraceAccessor(sd->traceMng())
  , m_sub_domain(sd)
  , m_node_family(node_variable.variable()->itemFamily())
  , m_node_variable(node_variable)
  , m_node_matrix_indexes(VariableBuildInfo(sd->defaultMeshHandle(), "NodeMatrixIndexes"))
  {}

 public:

  void build() { _computeMatrixInfo(); }

 private:

  void _computeMatrixInfo()
  {
    info() << "[AlephFem] COMPUTE_MATRIX_INFO\n";
    m_aleph_kernel = new AlephKernel(m_sub_domain, 2, 1);
    NodeGroup own_nodes = m_node_family->allItems().own();
    //Int32 nb_node = own_nodes.size();
    //Int32 total_nb_node = m_sub_domain->parallelMng()->reduce(Parallel::ReduceSum, nb_node);
    m_node_matrix_indexes.fill(-1);
    ENUMERATE_ (Node, inode, own_nodes) {
      Node node = *inode;
      Integer row = m_aleph_kernel->indexing()->get(m_node_variable, node);
      //info() << "ROW=" << row;
      m_node_matrix_indexes[node] = row;
    }
    //m_aleph_kernel->initialize(total_nb_node, nb_node);

#if 0
  const Integer row_offset = m_aleph_kernel->topology()->part()[m_aleph_kernel->rank()];
  
  ENUMERATE_ (Node, inode, ownNodes()) {
    Node node = *inode;
    m_node_matrix_indexes[node] += row_offset;
    info() << "ROW2=" << m_node_matrix_indexes[node];
  }
#endif

    m_aleph_matrix = m_aleph_kernel->createSolverMatrix();
    m_aleph_rhs_vector = m_aleph_kernel->createSolverVector();
    m_aleph_solution_vector = m_aleph_kernel->createSolverVector();
    m_aleph_matrix->create();
    m_aleph_rhs_vector->create();
    m_aleph_solution_vector->create();
  }

  void matrixAddValue(NodeLocalId row, NodeLocalId column, Real value) override
  {
    ItemInfoListView item_list_view(m_node_family);
    m_aleph_matrix->addValue(m_node_variable, item_list_view[row], m_node_variable, item_list_view[column], value);
  }

  void solve() override
  {
    UniqueArray<Real> aleph_result;

    m_aleph_matrix->assemble();
    m_aleph_rhs_vector->assemble();
    auto* aleph_solution_vector = m_aleph_solution_vector;
    NodeGroup own_nodes = m_node_family->allItems().own();
    UniqueArray<Real> vector_zero(own_nodes.size());
    vector_zero.fill(0.0);
    //aleph_solution_vector->create();
    aleph_solution_vector->setLocalComponents(vector_zero);
    aleph_solution_vector->assemble();

    auto* aleph_params = _createAlephParam();

    Int32 nb_iteration = 0;
    Real residual_norm = 0.0;
    info() << "[AlephFem] BEGIN SOLVING WITH ALEPH";

    m_aleph_matrix->solve(aleph_solution_vector,
                          m_aleph_rhs_vector,
                          nb_iteration,
                          &residual_norm,
                          aleph_params,
                          false);
    info() << "[AlephFem] END SOLVING WITH ALEPH r=" << residual_norm
           << " nb_iter=" << nb_iteration;
    auto* rhs_vector = m_aleph_kernel->createSolverVector();
    auto* solution_vector = m_aleph_kernel->createSolverVector();

    UniqueArray<Real> rhs_results;
    rhs_vector->getLocalComponents(rhs_results);
    solution_vector->getLocalComponents(aleph_result);

    bool do_verbose = own_nodes.size() < 200;
    Int32 index = 0;
    ENUMERATE_ (Node, inode, own_nodes) {
      Node node = *inode;

      m_node_variable[node] = aleph_result[m_aleph_kernel->indexing()->get(m_node_variable, node)];
      if (do_verbose)
        info() << "Node uid=" << node.uniqueId() << " V=" << aleph_result[index] << " T=" << m_node_variable[node]
               << " RHS=" << rhs_results[index];

      ++index;
    }
  }

  void setRHSValues(Span<const Real> values) override
  {
    m_aleph_rhs_vector->setLocalComponents(ConstArrayView(values.size(), values.data()));
  }

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
                              false, // m_param_write_matrix_to_file_error_strategy
                              "SolveErrorAlephMatrix.dbg", // m_param_write_matrix_name_error_strategy
                              false, // m_param_listing_output
                              0., // m_param_threshold
                              false, // m_param_print_cpu_time_resolution
                              0, // m_param_amg_coarsening_method: par défault celui de Sloop,
                              0, // m_param_output_level
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
  IItemFamily* m_node_family = nullptr;
  VariableNodeReal m_node_variable;
  VariableNodeInt32 m_node_matrix_indexes;
  AlephKernel* m_aleph_kernel = nullptr;
  AlephMatrix* m_aleph_matrix = nullptr;
  AlephVector* m_aleph_rhs_vector = nullptr;
  AlephVector* m_aleph_solution_vector = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

extern "C++" NodeLinearSystemImpl*
createAlephFemLinearSystemImpl(ISubDomain* sd, const Arcane::VariableNodeReal& node_variable)
{
  auto* x = new AlephFemLinearSystemImpl(sd, node_variable);
  x->build();
  return x;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
