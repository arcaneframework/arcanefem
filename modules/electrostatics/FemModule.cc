// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Simple module to solve electrostatics using FEM.                          */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModule at the start of the simulation.
 *
 *  - initializes degrees of freedom (DoFs) on nodes.
 *  - builds support for manufactured test case (optional).
 */
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "[ArcaneFem-Info] Started module compute()";
  Real elapsedTime = platform::getRealTime();

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs the main computation for the FemModule.
 *
 * - Stops the time loop after 1 iteration since the equation is steady state.
 * - Resets, configures, and initializes the linear system.
 * - Executes the stationary solve.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "[ArcaneFem-Info] Started module compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  if (m_petsc_flags != NULL){
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows via the following steps:
 *   1. _assembleBilinearOperator()  Assembles the FEM  matrix 𝐀
 *   2. _assembleLinearOperator()    Assembles the FEM RHS vector 𝐛
 *   3. _solve()                     Solves for solution vector 𝐮 = 𝐀⁻¹𝐛
 *   4. _updateVariables()           Updates FEM variables 𝐮 = 𝐱
 *   5. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  _getMaterialParameters();
  if(m_assemble_linear_system){
    _assembleBilinearOperator();
    _assembleLinearOperator();
  }
  if(m_solve_linear_system){
    _solve();
    _updateVariables();
  }
  if(m_cross_validation){
    _validateResults();
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Retrieves and sets the material parameters for the simulation.
 *
 * This method initializes material properties:
 *   - charge density (`rho`)
 *   - free space permittivity (`epsilon`)
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  rho = options()->rho();
  epsilon = options()->epsilon();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"get-material-param", elapsedTime);
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
 *  2. If a constant source term is specified (`rho`), apply it to the RHS.
 *  3. If Neumann BC are specified applied to the RHS.
 *  4. If Dirichlet BC are specified apply to the LHS & RHS. 
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (options()->rho.isPresent()) {
    Real qdot = -rho / epsilon;
    ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhs(qdot, mesh(), node_dof, m_node_coord, rhs_values);
  }

  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions()){
    ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);
  }

  for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
    ArcaneFemFunctions::BoundaryConditions2D::applyDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"rhs-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (P1 FE).
 *
 * This function calculates the integral of the expression:
 * (u.dx * v.dx + u.dy * v.dy)
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the gradients of the shape functions.
 * 3. Return (u.dx * v.dx + u.dy * v.dy);
 */
/*---------------------------------------------------------------------------*/

RealMatrix<3, 3> FemModule::
_computeElementMatrixTria3(Cell cell)
{
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  Real3 dxU = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyU = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

  return area * (dxU ^ dxU) + area * (dyU ^ dyU);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly given as mesh type.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (options()->meshType == "QUAD4")
    ARCANE_FATAL("Non supported meshType");
  else if (options()->meshType == "TRIA3")
    _assembleBilinear<3>([this](const Cell& cell) {
      return _computeElementMatrixTria3(cell);
    });
  else
    ARCANE_FATAL("Non supported meshType");

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"lhs-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system.
 *
 * The method performs the following steps:
 *   1. For each cell, retrieves the cell-specific constant `lambda`.
 *   2. Computes element matrix using provided `compute_element_matrix` function.
 *   3. Assembles global matrix by adding contributions from each cell's element
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

template <int N>
void FemModule::
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

void FemModule::
_solve()
{
  info() << "[ArcaneFem-Module] _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.solve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"solve", elapsedTime);
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

void FemModule::
_updateVariables()
{
  info() << "[ArcaneFem-Module] _updateVariables()";
  Real elapsedTime = platform::getRealTime();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    // Copy RHS DoF to Node u
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_u[node_dof.dofId(node, 0)];
      m_phi[node] = v;
    }
  }

  m_phi.synchronize();

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    m_E[cell] = ArcaneFemFunctions::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_phi);
  }

  m_E.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"update-variables", elapsedTime);
}

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

void FemModule::
_validateResults()
{
  info() << "[ArcaneFem-Module] _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200)
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "phi[" << node.localId() << "][" << node.uniqueId() << "] = " << m_phi[node];
    }

  String filename = options()->resultFile();
  info() << "ValidateResultFile filename=" << filename;

  if (!filename.empty())
    checkNodeResultFile(traceMng(), filename, m_phi, 1.0e-4);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "cross-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/