// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Simple FEM solver to solve acoustics.                                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModule at the start of the simulation.
 *
 * This method initializes degrees of freedom (DoFs) on nodes.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "[ArcaneFem-Info] Started module startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs the main computation for the FemModule.
 *
 * This method:
 *   1. Stops the time loop after 1 iteration since the equation is steady state.
 *   2. Resets, configures, and initializes the linear system.
 *   3. Executes the stationary solve.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "[ArcaneFem-Info] Started module compute()";
  Real elapsedTime = platform::getRealTime();

  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows a sequence of steps to solve FEM system:
 *   1. _getMaterialParameters()          Retrieves material parameters via
 *   2. _assembleBilinearOperator()       Assembles the FEM  matrix A
 *   3. _assembleLinearOperator()         Assembles the FEM RHS vector b
 *   4.  _solve()                         Solves for solution vector u = A^-1*b
 *   5. _updateVariables()                Updates FEM variables u = x
 *   6. _validateResults()                Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  _getMaterialParameters();
  _assembleBilinearOperatorTria3();
  _assembleLinearOperator();
  _solve();
  _updateVariables();
  _validateResults();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  m_kc2 = options()->kc2();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] get-material-parms", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the FEM linear operator.
 *
 * This method assembles the RHS 'b' of the linear system by:
 *   1. Initializing the RHS values to zero.
 *   2. Iterating over Neumann boundary conditions to apply it to RHS.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable for RHS values
  rhs_values.fill(0.0);

  const auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  BC::IArcaneFemBC* bc = options()->boundaryConditions();
  if(bc){
    for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions()){
      ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (P1 FE).
 *
 * This function calculates the integral of the expression:
 * -(u.dx * v.dx + u.dy * v.dy) + kc2 * u * v
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the integral U*V term.
 * 3. Compute the gradients of the shape functions.
 * 4. Return -(u.dx * v.dx + u.dy * v.dy) + kc2 * u * v ;
 */
/*---------------------------------------------------------------------------*/

FixedMatrix<3, 3> FemModule::
_computeElementMatrixTria3(Cell cell)
{

  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  Real3x3 UV = ArcaneFemFunctions::FeOperation2D::computeUVTria3(cell, m_node_coord);

  Real3 dxU = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyU = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

  return -area * (dxU ^ dxU) - area * (dyU ^ dyU) + m_kc2 * area * UV;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for triangular elements.
 *
 * This method computes and assembles the global  matrix A for the FEM linear
 * system. For each cell (triangle), it calculates the local element matrix &
 * adds the contributions to the global matrix based on the nodes of the cell.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTria3()
{
  info() << "[ArcaneFem-Info] Started module _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    auto K_e = _computeElementMatrixTria3(cell); // element matrix
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

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] lhs-matrix-assembly", elapsedTime);
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
  info() << "[ArcaneFem-Info] Started module _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.solve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] solve-linear-system", elapsedTime);
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
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      m_u[node] = dof_u[node_dof.dofId(node, 0)];
    }
  }

  m_u.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] update-variables", elapsedTime);
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
  info() << "[ArcaneFem-Info] Started module _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200) {
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "u[" << node.localId() << "][" << node.uniqueId() << "] = " << m_u[node];
    }
  }

  String filename = options()->resultFile();
  const double epsilon = 1.0e-4;
  const double min_value_to_test = 1.0e-16;

  info() << "[ArcaneFem-Info] Validating results filename=" << filename << " epsilon =" << epsilon;

  if (!filename.empty())
    checkNodeResultFile(traceMng(), filename, m_u, epsilon, min_value_to_test);

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Function to prints the execution time `value` of phase `label`
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_printArcaneFemTime(const String label, const Real value)
{
  info() << std::left << std::setw(40) << label << " = " << value;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
