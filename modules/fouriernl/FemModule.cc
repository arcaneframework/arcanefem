// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2026 */
/*                                                                           */
/* Simple module to solve Nonlinear Fourier's equation using FEM.            */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include "ConductivityCoefficient.h"
#include "ElementMatrix.h"
#include "ElementMatrixHexQuad.h"
#include <arcane/IParallelMng.h>

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModuleFourierNL at the start of the simulation.
 *
 *  - initializes degrees of freedom (DoFs) on nodes.
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
startInit()
{
  info() << "[ArcaneFem-Info] Started module startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();
  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->hasSolutionComparisonFile();
  m_petsc_flags = options()->petscFlags();
  m_hex_quad_mesh = options()->hexQuadMesh();

  m_perform_fixed_point_iters = options()->performFpIters();
  if(!m_perform_fixed_point_iters) {
    m_max_fp_iters = 1;
  } else {
    m_max_fp_iters = options()->maxFpIters();
  }
  m_fp_tol = options()->fpTol();

  m_qdot = options()->qdot();


  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs the main computation for the FemModuleFourierNL.
 *
 * - Stops the time loop after 1 iteration since the equation is steady state.
 * - Resets, configures, and initializes the linear system.
 * - Executes the stationary solve.
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
compute()
{
  info() << "[ArcaneFem-Info] Started module compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(),  acceleratorMng()->defaultRunner(),  m_dofs_on_nodes.dofFamily(), "Solver");

  if (m_petsc_flags != NULL){
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") {
    bool use_csr_in_linear_system =
    options()->linearSystem.serviceName() == "HypreLinearSystem" ||
    options()->linearSystem.serviceName() == "AlienLinearSystem" ||
    options()->linearSystem.serviceName() == "PETScLinearSystem";
    if (m_matrix_format == "BSR")
      m_bsr_format.initialize(mesh(), 1, use_csr_in_linear_system, 0);
    else
      m_bsr_format.initialize(mesh(), 1, use_csr_in_linear_system, 1);
    m_bsr_format.computeSparsity();
  }

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows via the following steps:
 *   1. _getMaterialParameters()     Retrieves material parameters
 *   2. _assembleBilinearOperator()  Assembles the FEM  matrix 𝐀
 *   3. _assembleLinearOperator()    Assembles the FEM RHS vector 𝐛
 *   4. _solve()                     Solves for solution vector 𝐮 = 𝐀⁻¹𝐛
 *   5. _updateVariables()           Updates FEM variables 𝐮 = 𝐱
 *   6. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
_doStationarySolve()
{
  info() << "[ArcaneFem-Info] Started module _doStationarySolve()";
  _updatePreviousIterationVariables();
  while (m_fp_iter < m_max_fp_iters) {
    if (m_assemble_linear_system) {

      if (m_linear_system.isInitialized() && m_fp_iter != 0) {
        m_linear_system.clearValues();

        if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR")
          m_bsr_format.resetMatrixValues();

        _assembleBilinearOperator();

        /* TODO : We should ideally not update the matrix row/columns concerning Dirichlet if Dirichlet BC are fixed for all iterations */
        /* need to create a separate function _assembleLinearOperatorLHSOnly. */
        _assembleLinearOperator(); // We use _assembleLinearOperator now for simplicity
      }
      else {
        _assembleBilinearOperator();
        _assembleLinearOperator();
      }
    }
    if (m_solve_linear_system) {
      _solve();
      _updateVariables();
    }

    ++m_fp_iter;
    _checkConvergence();

    if (m_converged) {
      info() << "[ArcaneFem-Info] Fixed-point iterations converged after " << m_fp_iter << " iterations";
      break;
    }
    else {
      _updatePreviousIterationVariables(); // copy u into uk for next iteration convergence check
      // m_uk.copy(m_u);
      _updateSolutionFromVariables(); // copy u into u_dof to update initial guess for linear solve TODO See how to use swap instead of deep copy
    }
  }
  if (m_fp_iter == m_max_fp_iters && !m_converged) {
    info() << "[ArcaneFem-Info] Fixed-point iterations did not converge after maximum (" << m_max_fp_iters << ") iterations";
    ARCANE_FATAL("Fixed-point iterations diverged after max iters");
  }
  if (m_cross_validation) {
    _validateResults();
  }
}

void FemModuleFourierNL::
_assembleLinearOperator()
{
  if (options()->linearSystem.serviceName() == "HypreLinearSystem" ||
      options()->linearSystem.serviceName() == "PETScLinearSystem")
    _assembleLinearOperatorGpu();
  else
    _assembleLinearOperatorCpu();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief FEM linear operator for the current simulation step.      
 *
 * This method constructs the linear  system by  assembling the LHS matrix 
 * and  RHS vector, applying various boundary conditions and source terms.
 *
 * Steps involved:
 *  1. The RHS vector is initialized to zero before applying any conditions.
 *  2. If a constant source term is specified (`qdot`), apply it to the RHS.
 *  3. If Neumann BC are specified applied to the RHS.
 *  4. If Dirichlet BC are specified apply to the LHS & RHS. 
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
_assembleLinearOperatorCpu()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR")
    m_bsr_format.toLinearSystem(m_linear_system);

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (options()->qdot.isPresent()) {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh)
        ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhsQuad4(m_qdot, mesh(), node_dof, m_node_coord, rhs_values);
      else
        ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhs(m_qdot, mesh(), node_dof, m_node_coord, rhs_values);
    }
    else {
      if (m_hex_quad_mesh)
        ArcaneFemFunctions::BoundaryConditions3D::applyConstantSourceToRhsHexa8(m_qdot, mesh(), node_dof, m_node_coord, rhs_values);
      else
        ArcaneFemFunctions::BoundaryConditions3D::applyConstantSourceToRhs(m_qdot, mesh(), node_dof, m_node_coord, rhs_values);
    }
  }

  // Helper lambda to apply boundary conditions
  auto applyBoundaryConditions = [&](auto BCFunctions) {
    BC::IArcaneFemBC* bc = options()->boundaryConditions();
    if (bc) {

      // Neumann
      for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions()) {
        if (mesh()->dimension() == 2) {
          if (m_hex_quad_mesh)
            ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhsQuad4(bs, node_dof, m_node_coord, rhs_values);
          else
            ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);
        }
        if (mesh()->dimension() == 3) {
          if (m_hex_quad_mesh)
            ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhsHexa8(bs, node_dof, m_node_coord, rhs_values);
          else
            ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);
        }
      }

      // Dirichlet
      for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
        ArcaneFemFunctions::BoundaryConditions::applyDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);
    }
  };

   if (mesh()->dimension() == 3)
     applyBoundaryConditions(ArcaneFemFunctions::BoundaryConditions3D());
   else
     applyBoundaryConditions(ArcaneFemFunctions::BoundaryConditions2D());

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief FEM linear operator for the current simulation step.
 * GPU compatible. Currently working with HypreDoFLinearSystem.
 *
 * This method constructs the linear  system by  assembling the LHS matrix
 * and  RHS vector, applying various boundary conditions and source terms.
 *
 * Steps involved:
 *  1. The RHS vector is initialized to zero before applying any conditions.
 *  2. If Neumann BC are specified applied to the RHS.
 *  3. If Dirichlet BC/Point are specified apply to the LHS & RHS.
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::_assembleLinearOperatorGpu()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperatorGpu()";
  Real elapsedTime = platform::getRealTime();

  m_bsr_format.toLinearSystem(m_linear_system);

  auto& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  // Because of Dirichlet (Penalty) implementation in Hypre.
  m_linear_system.getForcedInfo().fill(false);

  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  auto mesh_ptr = mesh();

  auto applyBoundaryConditions = [&](auto BCFunctions) {
    if (options()->qdot.isPresent())
      BCFunctions.applyConstantSourceToRhs(m_qdot, m_dofs_on_nodes, m_node_coord, rhs_values, mesh_ptr, queue);

    BC::IArcaneFemBC* bc = options()->boundaryConditions();
    if (bc) {
      for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions())
        BCFunctions.applyNeumannToRhs(bs, m_dofs_on_nodes, m_node_coord, rhs_values, mesh_ptr, queue);

      for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
        FemUtils::Gpu::BoundaryConditions::applyDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);

    }
  };

  if (mesh()->dimension() == 3)
    applyBoundaryConditions(FemUtils::Gpu::BoundaryConditions3D());
  else
    applyBoundaryConditions(FemUtils::Gpu::BoundaryConditions2D());

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"rhs-vector-assembly-gpu", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly given as mesh type.
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
_assembleBilinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto m_queue = subDomain()->acceleratorMng()->defaultQueue();
    auto command = makeCommand(m_queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    auto in_node_uk    = ax::viewIn(command, m_uk);

    if (mesh()->dimension() == 2)
      if (m_matrix_format == "BSR")
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord, in_node_uk); });
      else
        m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, in_node_uk, node_lid); });
    else
      if (m_matrix_format == "BSR")
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_node_uk); });
      else
        m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_node_uk, node_lid); });
  }

  if (m_matrix_format == "DOK") {
    if (mesh()->dimension() == 3)
      if(m_hex_quad_mesh)
        _assembleBilinear<8>([this](const Cell& cell) { return _computeElementMatrixHexa8(cell); });
      else
        _assembleBilinear<4>([this](const Cell& cell) { return _computeElementMatrixTetra4(cell); });

    if (mesh()->dimension() == 2)
      if(m_hex_quad_mesh)
        _assembleBilinear<4>([this](const Cell& cell) { return _computeElementMatrixQuad4(cell); });
      else
        _assembleBilinear<3>([this](const Cell& cell) { return _computeElementMatrixTria3(cell); });
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "lhs-matrix-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system.
 *
 * The method performs the following steps:
 *   1. Computes element matrix using provided `compute_element_matrix` function.
 *   2. Assembles global matrix by adding contributions from each cell's element
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

template <int N>
void FemModuleFourierNL::
_assembleBilinear(const std::function<RealMatrix<N, N>(const Cell&)>& compute_element_matrix)
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    auto K_e = compute_element_matrix(cell); // element matrix based on the provided function
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        Real v = K_e(n1_index, n2_index);
        if (node1.isOwn()) {
          m_linear_system.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Solves the linear system and updates the solution vector.
 *
 * This method performs the following actions:
 *   1. Solves the linear system to compute the solution.
 *   2. Copies the computed solution from the DoF to the node values.
 *   3. Synchronizes the updated node values.
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
_solve()
{
  info() << "[ArcaneFem-Info] Started module _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.applyLinearSystemTransformationAndSolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "solve-linear-system", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Check for convergence and Update the FEM variables.
 *
 * This method performs the following actions:
 *   1. Updates the FEM solutions from the solutions of linear solver on DOF
 *   2. Evaluates the error w.r.t the guess FEM variables using max norm
 *   3. Checks for convergence
 *   4. Updates guess or ends fixed point iterations
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
_checkConvergence()
{
  info() << "[ArcaneFem-Info] Started module _checkConvergence()";
  Real elapsedTime = platform::getRealTime();

  m_u.synchronize();
  m_uk.synchronize();

  Real max_error = 0.0;
  // Real l1_error = 0.0;
  {
    ENUMERATE_ (Node, inode, ownNodes()) {
      const Real error = abs(m_u[inode] - m_uk[inode]);

      max_error = math::max(error, max_error);
      // l1_error  += error;
    }
  }
  IParallelMng* pm = defaultMesh()->parallelMng();
  max_error = pm->reduce(Parallel::ReduceMax, max_error);
  // l1_error  = pm->reduce(Parallel::ReduceSum, l1_error);

  // if ( max_error < m_fp_tol || l1_error < m_fp_tol){
  if ( max_error < m_fp_tol){
  // if ( l1_error < m_fp_tol ){
    m_converged = true;
  } else {
    m_converged = false;
  }

  // info() << "[ArcaneFem-FP-iters] At fixed-point iteration "<< m_fp_iter <<": linf(max)-error = " << max_error << " and l1-error = " << l1_error;
  info() << "[ArcaneFem-Info] At fixed-point iteration "<< m_fp_iter <<": linf(max)-error = " << max_error;
  // info() << "[ArcaneFem-FP-iters] At fixed-point iteration "<< m_fp_iter <<": l1-error = " << l1_error;

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "check-convergence", elapsedTime);
}


/*---------------------------------------------------------------------------*/
/**
 * @brief Update the FEM variables.
 *
 * This method performs the following actions:
 *   1. Fetches values of solution from solved linear system to FEM variables,
 *      i.e., it copies RHS DOF to u.
 *   2. Performs synchronize of FEM variables across subdomains.
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
_updateVariables(bool verbose)
{
  info() << "[ArcaneFem-Info] Started module _updateVariables()";
  Real elapsedTime = platform::getRealTime();

  { // copies solution (and optionally exact solution) to FEM output
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_u[node_dof.dofId(node, 0)];
      Real m = options()->expNlin;
      m_u[node] = v;
      if (verbose) {
        info() << "u[" << node.uniqueId() << "] = " << m_u[node];
      }
      m_u_exact[node] = math::pow((math::pow(2.0, m + 1) - 1) * m_node_coord[node].x + 1, 1 / (m + 1)) - 1.0;
    }

  }

  m_u.synchronize();
  m_u_exact.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "update-variables", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Update the previous fixed-point iteration FEM variables.
 *
 * This method performs the following actions:
 *   1. Fetches values of linear solution vector to the previous
 *      fixed-point iteration FEM variable, i.e., it copies RHS DOF to uk.
 *   2. Performs synchronize of the previous fixed-point iteration FEM
 *      variable across subdomains.
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
_updatePreviousIterationVariables(bool verbose)
{
  info() << "[ArcaneFem-Info] Started module _updatePreviousIterationVariables()";
  Real elapsedTime = platform::getRealTime();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      m_uk[node] = dof_u[node_dof.dofId(node, 0)];
      if (verbose) {
        info() << "uk[" << node.uniqueId() << "] = " << m_uk[node];
      }
    }
  }

  m_uk.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"update-previous-iteration-variables", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Reinitialize the solution vector of the linear solve with the FEM variables.
 *
 * This method performs the following actions:
 *   1. Performs synchronize of FEM variables across subdomains.
 *   2. Fetches the FEM variables to the solution vector of the
 *      linear solver for next fixed point iteration.
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
_updateSolutionFromVariables()
{
  info() << "[ArcaneFem-Info] Started module _updateSolutionFromVariables()";
  Real elapsedTime = platform::getRealTime();

  m_u.synchronize();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      dof_u[node_dof.dofId(node, 0)] = m_u[node];
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "_update-solution-from-variables", elapsedTime);
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/**
 * @brief Validates and prints the results of the FEM computation.
 *
 * This method performs the following actions:
 *   1. If number of nodes < 200, prints the computed values for each node.
 *   2. Retrieves the filename for the result file from options.
 *   3. If a filename is provided, checks the computed results against result file.
 *
 * @note The result comparison uses a tolerance of 1.0e-4.
 */
/*---------------------------------------------------------------------------*/

void FemModuleFourierNL::
_validateResults()
{
  info() << "[ArcaneFem-Info] Started module _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200)
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "u[" << node.uniqueId() << "] = " << m_u[node];
    }

  String filename = options()->solutionComparisonFile();

  checkNodeResultFile(traceMng(), filename, m_u, 1.0e-4);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModuleFourierNL);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
