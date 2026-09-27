// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BoundaryConditions.cc                                       (C) 2000-2026 */
/*                                                                           */
/* Helper functions for handling boundary conditions.                        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "ArcaneFemFunctions.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Dirichlet boundary conditions to RHS and LHS.
 *
 * Updates the LHS matrix and RHS vector to enforce Dirichlet conditions.
 *
 * - For LHS matrix `𝐀`, the diagonal term for the Dirichlet DOF is set to `𝑃`.
 * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
 *
 * @param bs Boundary condition values.
 * @param node_dof DOF connectivity view.
 * @param node_coord Node coordinates.
 * @param m_linear_system Linear system for LHS.
 * @param rhs_values RHS values to update.
 */
void ArcaneFemFunctions::BoundaryConditions::
applyDirichletToLhsAndRhs(BC::IDirichletBoundaryCondition* bs,
                          const IndexedNodeDoFConnectivityView& node_dof,
                          DoFLinearSystem& m_linear_system,
                          VariableDoFReal& rhs_values)
{
  FaceGroup face_group = bs->getSurface();
  NodeGroup node_group = face_group.nodeGroup();
  const StringConstArrayView u_dirichlet_string = bs->getValue();
  for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
    if (u_dirichlet_string[dof_index] != "NULL") {
      Real value = std::stod(u_dirichlet_string[dof_index].localstr());
      if (bs->getEnforceDirichletMethod() == "Penalty") {
        Real penalty = bs->getPenalty();
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, node_dof, m_linear_system, rhs_values, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowElimination") {
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowColumnElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
      }
      else {
        ARCANE_FATAL("Unknown Dirichlet method");
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Point Dirichlet boundary conditions to RHS and LHS.
 *
 * Updates the LHS matrix and RHS vector to enforce the Dirichlet.
 *
 * - For LHS matrix `𝐀`, the diagonal term for the Dirichlet DOF is set to `𝑃`.
 * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
 *
 * @param bs Boundary condition values.
 * @param node_dof DOF connectivity view.
 * @param node_coord Node coordinates.
 * @param m_linear_system Linear system for LHS.
 * @param rhs_values RHS values to update.
 */
/*---------------------------------------------------------------------------*/
void ArcaneFemFunctions::BoundaryConditions::
applyPointDirichletToLhsAndRhs(BC::IDirichletPointCondition* bs,
                               const IndexedNodeDoFConnectivityView& node_dof,
                               DoFLinearSystem& m_linear_system,
                               VariableDoFReal& rhs_values)
{
  NodeGroup node_group = bs->getNode();
  const StringConstArrayView u_dirichlet_string = bs->getValue();
  for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
    if (u_dirichlet_string[dof_index] != "NULL") {
      Real value = std::stod(u_dirichlet_string[dof_index].localstr());
      if (bs->getEnforceDirichletMethod() == "Penalty") {
        Real penalty = bs->getPenalty();
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, node_dof, m_linear_system, rhs_values, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowElimination") {
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowColumnElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
      }
      else {
        ARCANE_FATAL("Unknown Dirichlet method");
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Dirichlet boundary conditions to RHS.
 *
 * Updates the RHS vector to enforce Dirichlet conditions.
 *
 * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
 *
 * @param bs Boundary condition values.
 * @param node_dof DOF connectivity view.
 * @param node_coord Node coordinates.
 * @param rhs_values RHS values to update.
 */
/*---------------------------------------------------------------------------*/
void ArcaneFemFunctions::BoundaryConditions::
applyDirichletToRhs(BC::IDirichletBoundaryCondition* bs,
                    const IndexedNodeDoFConnectivityView& node_dof,
                    VariableDoFReal& rhs_values)
{
  FaceGroup face_group = bs->getSurface();
  NodeGroup node_group = face_group.nodeGroup();
  const StringConstArrayView u_dirichlet_string = bs->getValue();
  for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
    if (u_dirichlet_string[dof_index] != "NULL") {
      Real value = std::stod(u_dirichlet_string[dof_index].localstr());
      if (bs->getEnforceDirichletMethod() == "Penalty") {
        Real penalty = bs->getPenalty();
        value = value * penalty;
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowElimination") {
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
      }
      else {
        ARCANE_FATAL("Unknown Dirichlet method");
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Point Dirichlet boundary conditions to RHS.
 *
 * Updates the RHS vector to enforce the Dirichlet.
 *
 * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
 *
 * @param bs Boundary condition values.
 * @param node_dof DOF connectivity view.
 * @param node_coord Node coordinates.
 * @param rhs_values RHS values to update.
 */
void ArcaneFemFunctions::BoundaryConditions::
applyPointDirichletToRhs(BC::IDirichletPointCondition* bs,
                         const IndexedNodeDoFConnectivityView& node_dof,
                         VariableDoFReal& rhs_values)
{
  NodeGroup node_group = bs->getNode();
  const StringConstArrayView u_dirichlet_string = bs->getValue();
  for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
    if (u_dirichlet_string[dof_index] != "NULL") {
      Real value = std::stod(u_dirichlet_string[dof_index].localstr());
      if (bs->getEnforceDirichletMethod() == "Penalty") {
        Real penalty = bs->getPenalty();
        value = value * penalty;
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowElimination") {
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
      }
      else {
        ARCANE_FATAL("Unknown Dirichlet method");
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Neumann conditions to the right-hand side (RHS) values.
 *
 * This method updates the RHS values of the finite element method equations
 * based on the provided Neumann boundary condition. The boundary condition
 * can specify a value or its components along the x and y directions.
 *
 * @param bs The Neumann boundary condition values.
 * @param node_dof Connectivity view for degrees of freedom at nodes.
 * @param node_coord Coordinates of the nodes in the mesh.
 * @param rhs_values The right-hand side values to be updated.
 */
/*---------------------------------------------------------------------------*/
void ArcaneFemFunctions::BoundaryConditions::
applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, IMesh* mesh,
                  const IndexedNodeDoFConnectivityView& node_dof,
                  const VariableNodeReal3& node_coord,
                  VariableDoFReal& rhs_values)
{
  ARCANE_CHECK_PTR(bs);
  ARCANE_CHECK_PTR(mesh);

  // Get mesh type via the number of nodes of fist cell works only for unifrom mesh
  Int32 nb_nodes = 0;
  {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();
    CellLocalId first_cell_lid(0);
    nb_nodes = cell_node_cv.nbNode(first_cell_lid);
  }

  if (mesh->dimension() == 2 && nb_nodes == 3) { // Triangular mesh
    ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhsTria3(bs, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 4) { // Quadrilateral mesh
    ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhsQuad4(bs, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 8) { // Quadrilateral mesh Quad8
    ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhsQuad8(bs, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 9) { // Quadrilateral mesh Quad9
    ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhsQuad9(bs, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 4) { // Tetrahedral mesh
    ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhsTetra4(bs, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 8) { // Hexahedral mesh
    ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhsHexa8(bs, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 20) { // Hexa20 mesh (Quad8 faces)
    ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhsHexa20(bs, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 27) { // Hexa27 mesh (Quad9 faces)
    ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhsHexa27(bs, node_dof, node_coord, rhs_values);
  }
  else {
    ARCANE_FATAL("Unsupported cell type in applyConstantNeumannToRhs()");
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector.
 *
 * This method adds a constant source term `qdot` to the RHS vector for each
 * node in the mesh. The contribution to each node is weighted by the area of
 * the cell and evenly distributed among the number of nodes of the cell.
 *
 * @param qdot The constant source term.
 * @param mesh The mesh containing all cells.
 * @param node_dof DOF connectivity view.
 * @param node_coord The coordinates of the nodes.
 * @param rhs_values The RHS values to update.
 */
void ArcaneFemFunctions::BoundaryConditions::
applyConstantSourceToRhs(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
{
  ARCANE_CHECK_PTR(mesh);

  // Get mesh type via the number of nodes of fist cell works only for unifrom mesh
  Int32 nb_nodes = 0;
  {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();
    CellLocalId first_cell_lid(0);
    nb_nodes = cell_node_cv.nbNode(first_cell_lid);
  }

  if (mesh->dimension() == 2 && nb_nodes == 3) { // Triangular mesh
    ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhsTria3(qdot, mesh, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 4) { // Quadrilateral mesh Quad4
    ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhsQuad4(qdot, mesh, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 8) { // Quadrilateral mesh Quad8
    ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhsQuad8(qdot, mesh, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 9) { // Quadrilateral mesh Quad9
    ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhsQuad9(qdot, mesh, node_dof, node_coord, rhs_values);
  }

  else if (mesh->dimension() == 3 && nb_nodes == 4) { // Tetrahedral mesh
    ArcaneFemFunctions::BoundaryConditions3D::applyConstantSourceToRhsTetra4(qdot, mesh, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 8) { // Hexahedral mesh
    ArcaneFemFunctions::BoundaryConditions3D::applyConstantSourceToRhsHexa8(qdot, mesh, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 20) { // Hexahedral mesh Hexa20
    ArcaneFemFunctions::BoundaryConditions3D::applyConstantSourceToRhsHexa20(qdot, mesh, node_dof, node_coord, rhs_values);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 27) { // Hexahedral mesh Hexa27
    ArcaneFemFunctions::BoundaryConditions3D::applyConstantSourceToRhsHexa27(qdot, mesh, node_dof, node_coord, rhs_values);
  }
  else {
    ARCANE_FATAL("Unsupported cell type in applyConstantSourceToRhs()");
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
