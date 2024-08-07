﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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

#include <arcane/utils/NumArray.h>
#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/StringList.h>

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>

#include "IDoFLinearSystemFactory.h"
#include "Fem_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/**
 * @brief A module for finite element method.
 *
 * This class handles the initialization and computation for finite element
 * method (FEM) simulations, providing methods to  set  up and solve linear
 * systems, assemble FEM operators, and perform result checks.
 */
/*---------------------------------------------------------------------------*/

class FemModule
: public ArcaneFemObject
{
 public:

  explicit FemModule(const ModuleBuildInfo& mbi)
  : ArcaneFemObject(mbi)
  , m_dofs_on_nodes(mbi.subDomain()->traceMng())
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
  }

  void compute() override;   //! Method called at each iteration 
  void startInit() override; //! Method called at the beginning of the simulation
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 private:

  Real m_kc2;

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnNodes m_dofs_on_nodes;

  void _doStationarySolve();
  void _getMaterialParameters();
  void _assembleBilinearOperatorTRIA3();
  void _solve();
  void _assembleLinearOperator();
  void _validateResults();

  FixedMatrix<3, 3> _computeElementMatrixTRIA3(Cell cell);

  Real _computeAreaTriangle3(Cell cell);
  Real _computeEdgeLength2(Face face);
  Real2 _computeEdgeNormal2(Face face);
};

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
      Real length = _computeEdgeLength2(face);
      Real2 normal = _computeEdgeNormal2(face);

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
 * @brief Computes the area of a triangle defined by three nodes.
 *
 * This method calculates the area using the determinant formula for a triangle.
 * The area is computed as half the value of the determinant of the matrix
 * formed by the coordinates of the triangle's vertices.
 */
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeAreaTriangle3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  return 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the length of the edge defined by a given face.
 * 
 * This method calculates Euclidean distance between the two nodes of the face.
 */
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeEdgeLength2(Face face)
{
  Real3 m0 = m_node_coord[face.nodeId(0)];
  Real3 m1 = m_node_coord[face.nodeId(1)];

  Real dx = m1.x - m0.x;
  Real dy = m1.y - m0.y;

  return math::sqrt(dx * dx + dy * dy);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the normalized edge normal for a given face.
 *
 * This method calculates normal vector to the edge defined by nodes of the face,
 * normalizes it, and ensures the correct orientation.
 */
/*---------------------------------------------------------------------------*/

Real2 FemModule::
_computeEdgeNormal2(Face face)
{
  Real3 m0 = m_node_coord[face.nodeId(0)];
  Real3 m1 = m_node_coord[face.nodeId(1)];

  if (!face.isSubDomainBoundaryOutside())
    std::swap(m0, m1);

  Real dx = m1.x - m0.x;
  Real dy = m1.y - m0.y;
  Real norm_N = math::sqrt(dx * dx + dy * dy);

  return { dy / norm_N, -dx / norm_N };
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (P1 FE).
 *
 * This function calculates the integral of the expression:
 * -(u.dx * v.dx + u.dy * v.dy) + kc2 * u * v
 *
 * Steps involved:
 * 1. Get the coordinates of the triangle vertices.
 * 2. Calculate the area of the triangle.
 * 3. Compute the gradients of the shape functions.
 * 4. Normalize the gradients by the area.
 * 5. Compute the element matrix:
 *    5.1 first compute -(u.dx * v.dx + u.dy * v.dy) terms
 *    5.2 then (kc2 * u * v) term adjusted
 */
/*---------------------------------------------------------------------------*/

FixedMatrix<3, 3> FemModule::
_computeElementMatrixTRIA3(Cell cell)
{
  // step 1
  Real3 vertex0 = m_node_coord[cell.nodeId(0)];
  Real3 vertex1 = m_node_coord[cell.nodeId(1)];
  Real3 vertex2 = m_node_coord[cell.nodeId(2)];

  // step 2
  Real area = _computeAreaTriangle3(cell);

  // step 3
  Real2 gradN0(vertex1.y - vertex2.y, vertex2.x - vertex1.x);
  Real2 gradN1(vertex2.y - vertex0.y, vertex0.x - vertex2.x);
  Real2 gradN2(vertex0.y - vertex1.y, vertex1.x - vertex0.x);

  // step 4
  Real A2 = 2.0 * area;
  gradN0 /= A2;
  gradN1 /= A2;
  gradN2 /= A2;

  // step 5
  FixedMatrix<3, 3> elementMatrix;

  // step 5.1
  elementMatrix(0, 0) = -area * Arcane::math::dot(gradN0, gradN0);
  elementMatrix(0, 1) = -area * Arcane::math::dot(gradN0, gradN1);
  elementMatrix(0, 2) = -area * Arcane::math::dot(gradN0, gradN2);

  elementMatrix(1, 0) = -area * Arcane::math::dot(gradN1, gradN0);
  elementMatrix(1, 1) = -area * Arcane::math::dot(gradN1, gradN1);
  elementMatrix(1, 2) = -area * Arcane::math::dot(gradN1, gradN2);

  elementMatrix(2, 0) = -area * Arcane::math::dot(gradN2, gradN0);
  elementMatrix(2, 1) = -area * Arcane::math::dot(gradN2, gradN1);
  elementMatrix(2, 2) = -area * Arcane::math::dot(gradN2, gradN2);

  // step 5 .2
  Real uvTerm = m_kc2 * area / 12.0;
  elementMatrix(0, 0) += uvTerm * 2.;
  elementMatrix(0, 1) += uvTerm;
  elementMatrix(0, 2) += uvTerm;

  elementMatrix(1, 0) += uvTerm;
  elementMatrix(1, 1) += uvTerm * 2.;
  elementMatrix(1, 2) += uvTerm;

  elementMatrix(2, 0) += uvTerm;
  elementMatrix(2, 1) += uvTerm;
  elementMatrix(2, 2) += uvTerm * 2.;

  return elementMatrix;
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
