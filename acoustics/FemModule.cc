// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2024 */
/*                                                                           */
/* Simple module to test simple FEM mechanism.                               */
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
  info() << "Module Fem INIT";

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();
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
  info() << "Module Fem COMPUTE";

  // step 1
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  // step 2
  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();

  // step 3
  _doStationarySolve();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows a sequence of steps to solve FEM system:
 *   1. _getMaterialParameters()          Retrieves material parameters via
 *   2. _assembleBilinearOperatorTRIA3()  Assembles the FEM  matrix A
 *   3. _assembleLinearOperator()         Assembles the FEM RHS vector b
 *   4.  _solve()                         Solves for solution vector u = A^-1*b
 *   5. _validateResults()                Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  _getMaterialParameters();
  _assembleBilinearOperatorTRIA3();
  _assembleLinearOperator();
  _solve();
  _validateResults();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "Get material parameters...";
  m_kc2 = options()->kc2();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the FEM linear operator.
 *
 * This method assembles the RHS 'b' of the linear system by:
 *   1. Initializing the RHS values to zero.
 *   2. Iterating over Neumann boundary conditions to apply constant flux terms.
 *      2.1 The integral of the flux term is computed over the boundary faces.
 *      2.2 Handles both scalar and vector flux components.
 *   3. For each boundary face, the edge length and normal are computed.
 *   4. The contribution to the RHS is calculated and accumulated for each node.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "Assembly of FEM linear operator ";

  // step 1
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable for RHS values
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // setp 2
  for (const auto& bs : options()->neumannBoundaryCondition()) {
    FaceGroup group = bs->surface();

    Real value = 0.0;
    Real valueX = 0.0;
    Real valueY = 0.0;
    bool hasValue = bs->value.isPresent();
    bool hasValueX = bs->valueX.isPresent();
    bool hasValueY = bs->valueY.isPresent();

    if (hasValue) {
      value = bs->value();
    }
    else {
      if (hasValueX)
        valueX = bs->valueX();
      if (hasValueY)
        valueY = bs->valueY();
    }

    // setp 2.1
    ENUMERATE_ (Face, iface, group) {
      Face face = *iface;

      // step 3.
      Real length = ArcaneFemFunctions::computeEdgeLength2(face, m_node_coord);
      Real2 normal = ArcaneFemFunctions::computeEdgeNormal2(face, m_node_coord);

      // step 4
      for (Node node : iface->nodes()) {
        if (!node.isOwn())
          continue;
        Real rhs_value = 0.0;

        // step 2.2
        if (hasValue) {
          rhs_value = value * length / 2.0;
        }
        else {
          rhs_value = (normal.x * valueX + normal.y * valueY) * length / 2.0;
        }

        rhs_values[node_dof.dofId(node, 0)] += rhs_value;
      }
    }
  }
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
_computeElementMatrixTRIA3(Cell cell)
{
  // step 1
  Real area = ArcaneFemFunctions::computeAreaTriangle3(cell, m_node_coord);

  // step 2
  Real3x3 UV = ArcaneFemFunctions::computeUVTria3(cell, m_node_coord);

  // step 3
  Real3 dxU = ArcaneFemFunctions::computeGradientXTria3(cell, m_node_coord);
  Real3 dyU = ArcaneFemFunctions::computeGradientYTria3(cell, m_node_coord);

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
_assembleBilinearOperatorTRIA3()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    auto K_e = _computeElementMatrixTRIA3(cell); // element matrix
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
  // step 1
  m_linear_system.solve();

  // step 2
  VariableDoFReal& dof_u(m_linear_system.solutionVariable());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Node, inode, ownNodes()) {
    Node node = *inode;
    m_u[node] = dof_u[node_dof.dofId(node, 0)];
  }

  // step 3
  m_u.synchronize();
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
  // setp 1
  if (allNodes().size() < 200) {
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "u[" << node.localId() << "][" << node.uniqueId() << "] = " << m_u[node];
    }
  }

  // setp 2
  String filename = options()->resultFile();
  info() << "ValidateResultFile filename=" << filename;

  // setp 3
  if (!filename.empty())
    checkNodeResultFile(traceMng(), filename, m_u, 1.0e-4);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
