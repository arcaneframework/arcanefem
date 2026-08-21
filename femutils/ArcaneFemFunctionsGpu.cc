// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArcaneFemFunctionsGpu.cc                                    (C) 2000-2026 */
/*                                                                           */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "ArcaneFemFunctionsGpu.h"
#include "FemUtilsGlobal.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils::Gpu
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BoundaryConditionsHelpers::
applyDirichletToNodeGroupRhsOnly(const Int32 dof_index, Real value, RunQueue* queue, IMesh* mesh,
                                 DoFLinearSystem& linear_system,
                                 const FemDoFsOnNodes& dofs_on_nodes, const NodeGroup& node_group)
{
  ARCANE_CHECK_PTR(queue);
  ARCANE_CHECK_PTR(mesh);

  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());

  auto command = makeCommand(queue);
  auto in_out_forced_info = viewInOut(command, linear_system.getForcedInfo());
  auto in_out_rhs_variable = viewInOut(command, linear_system.rhsVariable());

  command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, node_group)
  {
    if (nodes_infos.isOwn(node_lid)) {
      DoFLocalId dof_id = node_dof.dofId(node_lid, dof_index);
      in_out_forced_info[dof_id] = true;
      in_out_rhs_variable[dof_id] = value;
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BoundaryConditionsHelpers::
applyDirichletToNodeGroupViaPenalty(const Int32 dof_index, Real value, Real penalty, RunQueue* queue, IMesh* mesh,
                                    DoFLinearSystem& linear_system, const FemDoFsOnNodes& dofs_on_nodes, const NodeGroup& node_group)
{
  ARCANE_CHECK_PTR(queue);
  ARCANE_CHECK_PTR(mesh);

  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());

  auto command = makeCommand(queue);
  auto in_out_forced_info = viewInOut(command, linear_system.getForcedInfo());
  auto in_out_forced_value = viewInOut(command, linear_system.getForcedValue());
  auto in_out_rhs_variable = viewInOut(command, linear_system.rhsVariable());

  command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, node_group)
  {
    if (nodes_infos.isOwn(node_lid)) {
      DoFLocalId dof_id = node_dof.dofId(node_lid, dof_index);
      in_out_forced_info[dof_id] = true;
      in_out_forced_value[dof_id] = penalty;
      in_out_rhs_variable[dof_id] = penalty * value;
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BoundaryConditionsHelpers::
applyDirichletToNodeGroupViaRowOrRowColumnElimination(Byte elimination_type, const Int32 dof_index, Real value, RunQueue* queue,
                                                      DoFLinearSystem& linear_system, const FemDoFsOnNodes& dofs_on_nodes,
                                                      const NodeGroup& node_group)
{
  ARCANE_CHECK_PTR(queue);

  DoFLinearSystemRowEliminationHelper elimination_helper(linear_system.rowEliminationHelper());
  NodeInfoListView nodes_infos(node_group.itemFamily());
  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());

  auto command = makeCommand(queue);
  auto in_out_elimination_info = viewInOut(command, elimination_helper.getEliminationInfo());
  auto in_out_elimination_value = viewInOut(command, elimination_helper.getEliminationValue());
  command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, node_group)
  {
    if (nodes_infos.isOwn(node_lid)) {
      DoFLocalId dof_id = node_dof.dofId(node_lid, dof_index);
      in_out_elimination_info[dof_id] = elimination_type;
      in_out_elimination_value[dof_id] = value;
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Dirichlet boundary conditions to RHS and LHS.
 *
 * - For LHS matrix `𝐀`, the diagonal term for the Dirichlet DOF is set to `𝑃`.
 * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
 */
void BoundaryConditions::
applyDirichletToLhsAndRhs(BC::IDirichletBoundaryCondition* bs,
                          const FemDoFsOnNodes& dofs_on_nodes, DoFLinearSystem& linear_system,
                          IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);

  FaceGroup face_group = bs->getSurface();
  NodeGroup node_group = face_group.nodeGroup();

  const StringConstArrayView u_dirichlet_string = bs->getValue();

  for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
    if (u_dirichlet_string[dof_index] != "NULL") {

      Real value = std::stod(u_dirichlet_string[dof_index].localstr());

      if (bs->getEnforceDirichletMethod() == "Penalty") {
        Real penalty = bs->getPenalty();
        BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, queue, mesh, linear_system, dofs_on_nodes, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowElimination") {
        BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowOrRowColumnElimination(ELIMINATE_ROW, dof_index, value, queue, linear_system, dofs_on_nodes, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
        BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowOrRowColumnElimination(ELIMINATE_ROW_COLUMN, dof_index, value, queue, linear_system, dofs_on_nodes, node_group);
      }
      else {
        ARCANE_FATAL("Unknown method to enforce Dirichlet BC: '{0}'", bs->getEnforceDirichletMethod());
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Point Dirichlet boundary conditions to RHS and LHS.
 *
 * - For LHS matrix `𝐀`, the diagonal term for the Dirichlet DOF is set to `𝑃`.
 * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
 */
void BoundaryConditions::
applyPointDirichletToLhsAndRhs(BC::IDirichletPointCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                               DoFLinearSystem& linear_system, IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);
  NodeGroup node_group = bs->getNode();

  const StringConstArrayView u_dirichlet_str = bs->getValue();

  for (Int32 dof_index = 0; dof_index < u_dirichlet_str.size(); ++dof_index) {
    if (u_dirichlet_str[dof_index] != "NULL") {
      Real value = std::stod(u_dirichlet_str[dof_index].localstr());

      if (bs->getEnforceDirichletMethod() == "Penalty") {
        Real penalty = bs->getPenalty();
        BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, queue, mesh, linear_system, dofs_on_nodes, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowElimination") {
        BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowOrRowColumnElimination(ELIMINATE_ROW, dof_index, value, queue, linear_system, dofs_on_nodes, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
        BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowOrRowColumnElimination(ELIMINATE_ROW_COLUMN, dof_index, value, queue, linear_system, dofs_on_nodes, node_group);
      }
      else {
        ARCANE_FATAL("Unknown method to enforce Dirichlet BC: '{0}'", bs->getEnforceDirichletMethod());
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Dirichlet boundary conditions to RHS.
 *
 * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
 */
void BoundaryConditions::
applyDirichletToRhs(BC::IDirichletBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                    DoFLinearSystem& linear_system, IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);

  FaceGroup face_group = bs->getSurface();
  NodeGroup node_group = face_group.nodeGroup();

  const StringConstArrayView u_dirichlet_string = bs->getValue();

  for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
    if (u_dirichlet_string[dof_index] != "NULL") {

      Real value = std::stod(u_dirichlet_string[dof_index].localstr());

      if (bs->getEnforceDirichletMethod() == "Penalty") {
        Real penalty = bs->getPenalty();
        value = value * penalty;
        BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, queue, mesh, linear_system, dofs_on_nodes, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowElimination") {
        BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, queue, mesh, linear_system, dofs_on_nodes, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
        ARCANE_THROW(Arccore::NotImplementedException, "RowColumnElimination is not supported.");
      }
      else {
        ARCANE_FATAL("Unknown method to enforce Dirichlet BC: '{0}'", bs->getEnforceDirichletMethod());
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Point Dirichlet boundary conditions to RHS.
 *
 * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
 */
void BoundaryConditions::
applyPointDirichletToRhs(BC::IDirichletPointCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                         DoFLinearSystem& linear_system, IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);
  NodeGroup node_group = bs->getNode();

  const StringConstArrayView u_dirichlet_str = bs->getValue();

  for (Int32 dof_index = 0; dof_index < u_dirichlet_str.size(); ++dof_index) {
    if (u_dirichlet_str[dof_index] != "NULL") {
      Real value = std::stod(u_dirichlet_str[dof_index].localstr());

      if (bs->getEnforceDirichletMethod() == "Penalty") {
        Real penalty = bs->getPenalty();
        value = value * penalty;
        BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, queue, mesh, linear_system, dofs_on_nodes, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowElimination") {
        BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, queue, mesh, linear_system, dofs_on_nodes, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
        ARCANE_THROW(Arccore::NotImplementedException, "RowColumnElimination is not supported.");
      }
      else {
        ARCANE_FATAL("Unknown method to enforce Dirichlet BC: '{0}'", bs->getEnforceDirichletMethod());
      }
    }
  }
}

/*---------------------------------------------------------------------------*/

void BoundaryConditions::
applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                  const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                  IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);
  ARCANE_CHECK_PTR(mesh);
  ARCANE_CHECK_PTR(queue);

  // Get mesh type via the number of nodes of fist cell works only for unifrom mesh
  Int32 nb_nodes = 0;
  {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();
    CellLocalId first_cell_lid(0);
    nb_nodes = cell_node_cv.nbNode(first_cell_lid);
  }

  if (mesh->dimension() == 2 && nb_nodes == 3) { // Triangular mesh
    BoundaryConditions2D::applyNeumannToRhsTria3(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 4) { // Quad mesh
    BoundaryConditions2D::applyNeumannToRhsQuad4(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 8) { // Quadratic serendipity quad mesh
    BoundaryConditions2D::applyNeumannToRhsQuad8(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 9) { // Quadratic Lagrange quad mesh
    BoundaryConditions2D::applyNeumannToRhsQuad9(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 4) { // Tetra mesh
    BoundaryConditions3D::applyNeumannToRhsTetra4(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 8) { // Hexa mesh
    BoundaryConditions3D::applyNeumannToRhsHexa8(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 20) { // Quadratic serendipity hexa mesh
    BoundaryConditions3D::applyNeumannToRhsHexa20(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 27) { // Quadratic Lagrange hexa mesh
    BoundaryConditions3D::applyNeumannToRhsHexa27(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else {
    ARCANE_FATAL("Unknown mesh type only works for uniform TRIA3, QUAD4, QUAD8, QUAD9, TETRA4, HEXA8, HEXA20, HEXA27");
  }
}
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector.
 *
 * This method adds a constant source term `qdot` to the RHS vector for each
 * node in the mesh. The contribution to each node is weighted by the area of
 * the cell and evenly distributed among the number of nodes of the cell.
 *
 */
void BoundaryConditions::
applyConstantSourceToRhs(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
                         const VariableNodeReal3& node_coord,
                         VariableDoFReal& rhs_variable_na,
                         IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(mesh);
  ARCANE_CHECK_PTR(queue);

  // Get mesh type via the number of nodes of fist cell works only for unifrom mesh
  Int32 nb_nodes = 0;
  {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();
    CellLocalId first_cell_lid(0);
    nb_nodes = cell_node_cv.nbNode(first_cell_lid);
  }

  if (mesh->dimension() == 2 && nb_nodes == 3) { // Triangular mesh
    BoundaryConditions2D::applyConstantSourceToRhsTria3(qdot, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 4) { // Quad mesh
    BoundaryConditions2D::applyConstantSourceToRhsQuad4(qdot, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 8) { // Quadratic serendipity quad mesh
    BoundaryConditions2D::applyConstantSourceToRhsQuad8(qdot, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 2 && nb_nodes == 9) { // Quadratic Lagrange quad mesh
    BoundaryConditions2D::applyConstantSourceToRhsQuad9(qdot, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 4) { // Tetra mesh
    BoundaryConditions3D::applyConstantSourceToRhsTetra4(qdot, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 8) { // Hexa mesh
    BoundaryConditions3D::applyConstantSourceToRhsHexa8(qdot, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 20) { // Quadratic serendipity hexa mesh
    BoundaryConditions3D::applyConstantSourceToRhsHexa20(qdot, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else if (mesh->dimension() == 3 && nb_nodes == 27) { // Quadratic Lagrange hexa mesh
    BoundaryConditions3D::applyConstantSourceToRhsHexa27(qdot, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
  }
  else {
    ARCANE_FATAL("Unknown mesh type only works for uniform TRIA3, QUAD4, QUAD8, QUAD9, TETRA4, HEXA8, HEXA20, HEXA27");
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector.
 *
 * This method adds a constant source term `qdot` to the RHS vector for each
 * node in the mesh. The contribution to each node is weighted by the area of
 * the cell and evenly distributed among the number of nodes of the cell.
 *
 */
void BoundaryConditions2D::
applyConstantSourceToRhsTria3(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
                         const VariableNodeReal3& node_coord,
                         VariableDoFReal& rhs_variable_na,
                         IMesh* mesh, RunQueue* queue)
{
  auto func = [] ARCCORE_HOST_DEVICE(CellLocalId cell_lid,
                                     const IndexedCellNodeConnectivityView& cn_cv,
                                     const Accelerator::VariableNodeReal3InView& in_node_coord) {
    return MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);
  };
  BoundaryConditionsHelpers::applyConstantSourceToRhsBase(func, qdot, queue, rhs_variable_na, mesh, dofs_on_nodes, node_coord);
}

// TODO Atomic-free version of applyConstantSourceToRhsQuad4 for better performance
void BoundaryConditions2D::
applyConstantSourceToRhsQuad4(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
                              const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                              IMesh* mesh, RunQueue* queue)
{
  UnstructuredMeshConnectivityView connectivity_view;
  connectivity_view.setMesh(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
  auto cn_cv = connectivity_view.cellNode();

  auto command = makeCommand(queue);
  auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
  auto in_node_coord = viewIn(command, node_coord);

  command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh->allCells())
  {
    // Gauss points and weights (2x2 Gauss Quadrature)
    const Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 };
    const Real weights[2] = { 1.0, 1.0 };

    // Fetch cell nodes into a local array for quick access inside the integration loops
    NodeLocalId cell_nodes[4];
    Int32 index = 0;
    for (NodeLocalId node_lid : cn_cv.nodes(cell_lid)) {
      if (index < 4) {
        cell_nodes[index++] = node_lid;
      }
    }

    for (Int32 ixi = 0; ixi < 2; ++ixi) {
      for (Int32 ieta = 0; ieta < 2; ++ieta) {

        Real xi = gp[ixi];
        Real eta = gp[ieta];
        Real weight = weights[ixi] * weights[ieta];

        // Shape functions 𝐍 for Quad4
        RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);

        // Shape function derivatives
        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad4(xi, eta);
        // Jacobian calculation
        Real J00 = 0.0, J01 = 0.0, J10 = 0.0, J11 = 0.0;
        for (Int8 a = 0; a < 4; ++a) {
          Real3 coord = in_node_coord[cell_nodes[a]];
          J00 += reference_gradients.dN_dxi[a] * coord.x;
          J01 += reference_gradients.dN_dxi[a] * coord.y;
          J10 += reference_gradients.dN_deta[a] * coord.x;
          J11 += reference_gradients.dN_deta[a] * coord.y;
        }

        // Determinant of the Jacobian
        Real detJ = J00 * J11 - J01 * J10;
        Real integration_weight = weight * detJ;

        // Assemble RHS via Atomic Operations
        for (Int32 i = 0; i < 4; ++i) {
          NodeLocalId node_lid = cell_nodes[i];
          if (nodes_infos.isOwn(node_lid)) {
            Real rhs_value = N[i] * qdot * integration_weight;
            Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
          }
        }
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS quadratic quad elements.
 *
 * This method adds a constant source term `qdot` to the RHS vector for each
 * node in the mesh. The contribution to each node is weighted by the area of
 * the cell and evenly distributed among the number of nodes of the cell.
 *
 */
void _applyConstantSourceToRhsQuadraticQuad(Real qdot, Int32 nb_cell_node,
                                            const FemDoFsOnNodes& dofs_on_nodes,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_variable_na,
                                            IMesh* mesh, RunQueue* queue)
{
  UnstructuredMeshConnectivityView connectivity_view(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof = dofs_on_nodes.nodeDoFConnectivityView();
  auto cn_cv = connectivity_view.cellNode();

  if (nb_cell_node == 8) {
    auto command = makeCommand(queue);
    auto in_out_rhs = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);

    command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh->allCells())
    {
      // Gauss points and weights (3x3 Gauss Quadrature)
      constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
      constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

      // Fetch cell nodes into a local array for quick access inside the integration loops
      NodeLocalId cell_nodes[8];
      for (Int32 i = 0; i < 8; ++i)
        cell_nodes[i] = cn_cv.nodeId(cell_lid, i);

      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          const Real xi = gp[ixi];
          const Real eta = gp[ieta];

          // Shape functions 𝐍 for Quad8
          RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad8(xi, eta);

          // Determinant of the Jacobian
          const Real det_j = FeOperation2D::computeGradientsAndJacobianQuad8Gpu(cell_lid, cn_cv, in_node_coord, xi, eta).det_j;
          const Real integration_weight = weights[ixi] * weights[ieta] * det_j;

          // Assemble RHS via Atomic Operations
          for (Int32 i = 0; i < 8; ++i) {
            const NodeLocalId node_lid = cell_nodes[i];
            if (nodes_infos.isOwn(node_lid)) {
              const Real rhs_value = N[i] * qdot * integration_weight;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        }
      }
    };
  }
  else if (nb_cell_node == 9) {
    auto command = makeCommand(queue);
    auto in_out_rhs = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);

    command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh->allCells())
    {
      // Gauss points and weights (3x3 Gauss Quadrature)
      constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
      constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

      // Fetch cell nodes into a local array for quick access inside the integration loops
      NodeLocalId cell_nodes[9];
      for (Int32 i = 0; i < 9; ++i)
        cell_nodes[i] = cn_cv.nodeId(cell_lid, i);

      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          const Real xi = gp[ixi];
          const Real eta = gp[ieta];

          // Shape functions 𝐍 for Quad9
          RealVector<9> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad9(xi, eta);

          // Determinant of the Jacobian
          const Real det_j = FeOperation2D::computeGradientsAndJacobianQuad9Gpu(cell_lid, cn_cv, in_node_coord, xi, eta).det_j;
          const Real integration_weight = weights[ixi] * weights[ieta] * det_j;

          // Assemble RHS via Atomic Operations
          for (Int32 i = 0; i < 9; ++i) {
            const NodeLocalId node_lid = cell_nodes[i];
            if (nodes_infos.isOwn(node_lid)) {
              const Real rhs_value = N[i] * qdot * integration_weight;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        }
      }
    };
  }
  else
    ARCANE_FATAL("Unsupported quadratic quadrilateral with '{0}' nodes", nb_cell_node);
}

void _applyNeumannToRhsLine3(BC::INeumannBoundaryCondition* bs,
                             const FemDoFsOnNodes& dofs_on_nodes,
                             const VariableNodeReal3& node_coord,
                             VariableDoFReal& rhs_variable_na,
                             IMesh* mesh, RunQueue* queue)
{
  const FaceGroup group = bs->getSurface();
  const StringConstArrayView neumann_str = bs->getValue();
  const bool scalar_neumann = neumann_str.size() == 1 && neumann_str[0] != "NULL";
  const Real scalar_value = scalar_neumann ? std::stod(neumann_str[0].localstr()) : 0.0;
  const Real value_x = !scalar_neumann && neumann_str.size() > 0 && neumann_str[0] != "NULL"
  ? std::stod(neumann_str[0].localstr())
  : 0.0;
  const Real value_y = !scalar_neumann && neumann_str.size() > 1 && neumann_str[1] != "NULL"
  ? std::stod(neumann_str[1].localstr())
  : 0.0;

  UnstructuredMeshConnectivityView connectivity_view(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  FaceInfoListView faces_infos(mesh->faceFamily());
  auto node_dof = dofs_on_nodes.nodeDoFConnectivityView();
  auto fn_cv = connectivity_view.faceNode();

  auto command = makeCommand(queue);
  auto in_out_rhs = viewInOut(command, rhs_variable_na);
  auto in_node_coord = viewIn(command, node_coord);

  command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
  {
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    NodeLocalId face_nodes[3];
    Real3 coords[3];
    for (Int32 i = 0; i < 3; ++i) {
      face_nodes[i] = fn_cv.nodeId(face_lid, i);
      coords[i] = in_node_coord[face_nodes[i]];
    }
    const Real orientation = faces_infos.isSubDomainBoundaryOutside(face_lid) ? 1.0 : -1.0;

    for (Int32 igauss = 0; igauss < 3; ++igauss) {
      const Real xi = gp[igauss];
      RealVector<3> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsLine3(xi);

      const Real dN[3] = { xi - 0.5, xi + 0.5, -2.0 * xi };

      Real dx_dxi = 0.0;
      Real dy_dxi = 0.0;
      for (Int32 i = 0; i < 3; ++i) {
        dx_dxi += dN[i] * coords[i].x;
        dy_dxi += dN[i] * coords[i].y;
      }

      const Real jacobian = math::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);
      const Real normal_x = orientation * dy_dxi / jacobian;
      const Real normal_y = orientation * -dx_dxi / jacobian;
      const Real flux = scalar_neumann ? scalar_value : normal_x * value_x + normal_y * value_y;
      const Real integration_weight = weights[igauss] * jacobian;

      for (Int32 i = 0; i < 3; ++i) {
        const NodeLocalId node_lid = face_nodes[i];
        if (nodes_infos.isOwn(node_lid)) {
          const Real rhs_value = flux * N[i] * integration_weight;
          Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs[node_dof.dofId(node_lid, 0)], rhs_value);
        }
      }
    }
  };
}

void BoundaryConditions2D::
applyConstantSourceToRhsQuad8(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
                              const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                              IMesh* mesh, RunQueue* queue)
{
  _applyConstantSourceToRhsQuadraticQuad(qdot, 8, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
}

void BoundaryConditions2D::
applyConstantSourceToRhsQuad9(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
                              const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                              IMesh* mesh, RunQueue* queue)
{
  _applyConstantSourceToRhsQuadraticQuad(qdot, 9, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
}

void BoundaryConditions2D::
applyNeumannToRhsQuad8(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                       const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                       IMesh* mesh, RunQueue* queue)
{
  _applyNeumannToRhsLine3(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
}

void BoundaryConditions2D::
applyNeumannToRhsQuad9(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                       const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                       IMesh* mesh, RunQueue* queue)
{
  _applyNeumannToRhsLine3(bs, dofs_on_nodes, node_coord, rhs_variable_na, mesh, queue);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Neumann conditions to the right-hand side (RHS) values.
 *
 * This method updates the RHS values of the finite element method equations
 * based on the provided Neumann boundary condition. The boundary condition
 * can specify a value or its components along the x and y directions.
 */
void BoundaryConditions2D::
applyNeumannToRhsTria3(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                  const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                  IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);
  ARCANE_CHECK_PTR(queue);
  ARCANE_CHECK_PTR(mesh);

  FaceGroup group = bs->getSurface();
  bool scalarNeumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalarNeumann = true;
  }

  UnstructuredMeshConnectivityView connectivity_view;
  connectivity_view.setMesh(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
  auto fn_cv = connectivity_view.faceNode();

  if (scalarNeumann) {
    Real value = std::stod(neumann_str[0].localstr());

    auto command = makeCommand(queue);
    auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);
    command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
    {
      Real length = MeshOperation::computeLengthFace(face_lid, fn_cv, in_node_coord);
      for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
        if (nodes_infos.isOwn(node_lid)) {
          Real rhs_value = value * length / 2.0;
          Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
        }
      }
    };
  }
  else {
    Real value_x = neumann_str[0] != "NULL" ? std::stod(neumann_str[0].localstr()) : 0.0;
    Real value_y = neumann_str[1] != "NULL" ? std::stod(neumann_str[1].localstr()) : 0.0;

    auto command = makeCommand(queue);
    auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);
    FaceInfoListView faces_infos(mesh->faceFamily());
    command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
    {
      Real length = MeshOperation::computeLengthFace(face_lid, fn_cv, in_node_coord);
      Real2 normal = MeshOperation::computeNormalFace(face_lid, fn_cv, in_node_coord, faces_infos);
      for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
        if (nodes_infos.isOwn(node_lid)) {
          Real rhs_value = (normal.x * value_x + normal.y * value_y) * length / 2.0;
          Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
        }
      }
    };
  }
}

void BoundaryConditions2D::
applyNeumannToRhsQuad4(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                       const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                       IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);
  ARCANE_CHECK_PTR(queue);
  ARCANE_CHECK_PTR(mesh);

  FaceGroup group = bs->getSurface();
  bool scalarNeumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalarNeumann = true;
  }

  UnstructuredMeshConnectivityView connectivity_view;
  connectivity_view.setMesh(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  FaceInfoListView faces_infos(mesh->faceFamily());
  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
  auto fn_cv = connectivity_view.faceNode();

  if (scalarNeumann) {
    Real value = std::stod(neumann_str[0].localstr());

    auto command = makeCommand(queue);
    auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);

    command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
    {
      const Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 };
      const Real weights[2] = { 1.0, 1.0 };

      Real length = MeshOperation::computeLengthFace(face_lid, fn_cv, in_node_coord);

      // Cache face nodes locally
      NodeLocalId face_nodes[2];
      Int32 index = 0;
      for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
        if (index < 2) {
          face_nodes[index++] = node_lid;
        }
      }

      for (Int32 i = 0; i < 2; ++i) {
        Real xi = gp[i];
        Real weight = weights[i];

        RealVector<2> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsLine2(xi);

        Real integration_weight = weight * length * 0.5;

        for (Int32 j = 0; j < 2; ++j) {
          NodeLocalId node_lid = face_nodes[j];
          if (nodes_infos.isOwn(node_lid)) {
            Real rhs_value = value * N[j] * integration_weight;
            Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
          }
        }
      }
    };
  }
  else {
    Real valueX = neumann_str[0] != "NULL" ? std::stod(neumann_str[0].localstr()) : 0.0;
    Real valueY = neumann_str[1] != "NULL" ? std::stod(neumann_str[1].localstr()) : 0.0;

    auto command = makeCommand(queue);
    auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);

    command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
    {
      const Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 };
      const Real weights[2] = { 1.0, 1.0 };

      Real length = MeshOperation::computeLengthFace(face_lid, fn_cv, in_node_coord);
      Real2 normal = MeshOperation::computeNormalFace(face_lid, fn_cv, in_node_coord, faces_infos);

      // Cache face nodes locally
      NodeLocalId face_nodes[2];
      Int32 index = 0;
      for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
        if (index < 2) {
          face_nodes[index++] = node_lid;
        }
      }

      for (Int32 i = 0; i < 2; ++i) {
        Real xi = gp[i];
        Real weight = weights[i];

        RealVector<2> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsLine2(xi);

        Real integration_weight = weight * length * 0.5;

        for (Int32 j = 0; j < 2; ++j) {
          NodeLocalId node_lid = face_nodes[j];
          if (nodes_infos.isOwn(node_lid)) {
            Real rhs_value = (normal.x * valueX + normal.y * valueY) * N[j] * integration_weight;
            Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
          }
        }
      }
    };
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
 */
/*---------------------------------------------------------------------------*/

void BoundaryConditions3D::
applyConstantSourceToRhsTetra4(Real qdot, const FemDoFsOnNodes& dofs_on_nodes, const VariableNodeReal3& node_coord,
                         VariableDoFReal& rhs_variable_na, IMesh* mesh, RunQueue* queue)
{
  auto func = [] ARCCORE_HOST_DEVICE(CellLocalId cell_lid,
                                     const IndexedCellNodeConnectivityView& cn_cv,
                                     const Accelerator::VariableNodeReal3InView& in_node_coord) {
    return MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);
  };
  BoundaryConditionsHelpers::applyConstantSourceToRhsBase(func, qdot, queue, rhs_variable_na, mesh, dofs_on_nodes, node_coord);
}

void BoundaryConditions3D::
applyConstantSourceToRhsHexa8(Real qdot, const FemDoFsOnNodes& dofs_on_nodes, const VariableNodeReal3& node_coord,
                              VariableDoFReal& rhs_variable_na, IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(mesh);
  ARCANE_CHECK_PTR(queue);

  UnstructuredMeshConnectivityView connectivity_view;
  connectivity_view.setMesh(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
  auto cn_cv = connectivity_view.cellNode();

  auto command = makeCommand(queue);
  auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
  auto in_node_coord = viewIn(command, node_coord);

  command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh->allCells())
  {

    // Gauss quadrature for Hexa8
    // Using 2x2x2 Gauss points for integration
    constexpr Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 }; // {-1/sqrt(3) 1/sqrt(3)}
    constexpr Real weights[2] = { 1.0, 1.0 };

    // Fetch cell nodes into a local array for quick access inside the integration loops
    NodeLocalId cell_nodes[8];
    Int32 index = 0;
    for (NodeLocalId node_lid : cn_cv.nodes(cell_lid)) {
      if (index < 8) {
        cell_nodes[index++] = node_lid;
      }
    }

    for (Int32 ixi = 0; ixi < 2; ++ixi) {
      for (Int32 ieta = 0; ieta < 2; ++ieta) {
        for (Int32 izeta = 0; izeta < 2; ++izeta) {

          // Gauss point coordinates in reference space
          Real xi = gp[ixi]; // ξ coordinate
          Real eta = gp[ieta]; // η coordinate
          Real zeta = gp[izeta]; // ζ coordinate
          Real weight = weights[ixi] * weights[ieta] * weights[izeta];

          // Shape functions 𝐍 for Hexa8
          RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa8(xi, eta, zeta);

          // Shape function derivatives in reference space
          //  ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ  ∂𝑁₅/∂ξ  ∂𝑁₆/∂ξ  ∂𝑁₇/∂ξ  ∂𝑁₈/∂ξ ]
          //  ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η  ∂𝑁₅/∂η  ∂𝑁₆/∂η  ∂𝑁₇/∂η  ∂𝑁₈/∂η ]
          //  ∂𝐍/∂ζ = [ ∂𝑁₁/∂ζ  ∂𝑁₂/∂ζ  ∂𝑁₃/∂ζ  ∂𝑁₄/∂ζ  ∂𝑁₅/∂ζ  ∂𝑁₆/∂ζ  ∂𝑁₇/∂ζ  ∂𝑁₈/∂ζ ]
          const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsHexa8(xi, eta, zeta);
          // Jacobian for 3D (using your working stiffness matrix approach)
          Real3x3 J;
          for (Int8 a = 0; a < 8; ++a) {
            const Real3& n = in_node_coord[cell_nodes[a]];
            J[0][0] += reference_gradients.dN_dxi[a] * n.x; // ∂𝑥/∂ξ
            J[0][1] += reference_gradients.dN_dxi[a] * n.y; // ∂𝑦/∂ξ
            J[0][2] += reference_gradients.dN_dxi[a] * n.z; // ∂𝑧/∂ξ
            J[1][0] += reference_gradients.dN_deta[a] * n.x; // ∂𝑥/∂η
            J[1][1] += reference_gradients.dN_deta[a] * n.y; // ∂𝑦/∂η
            J[1][2] += reference_gradients.dN_deta[a] * n.z; // ∂𝑧/∂η
            J[2][0] += reference_gradients.dN_dzeta[a] * n.x; // ∂𝑥/∂ζ
            J[2][1] += reference_gradients.dN_dzeta[a] * n.y; // ∂𝑦/∂ζ
            J[2][2] += reference_gradients.dN_dzeta[a] * n.z; // ∂𝑧/∂ζ
          }

          // Compute determinant of Jacobian
          Real detJ = math::matrixDeterminant(J);

          // Compute integration weight
          Real integration_weight = weight * detJ;

          // Assemble RHS
          for (Int32 i = 0; i < 8; ++i) {
            NodeLocalId node_lid = cell_nodes[i];
            if (nodes_infos.isOwn(node_lid)) {
              Real rhs_value = N[i] * qdot * integration_weight;

              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        }
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BoundaryConditions3D::
applyConstantSourceToRhsHexa20(Real qdot, const FemDoFsOnNodes& dofs_on_nodes, const VariableNodeReal3& node_coord,
                               VariableDoFReal& rhs_variable_na, IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(mesh);
  ARCANE_CHECK_PTR(queue);

  UnstructuredMeshConnectivityView connectivity_view(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof = dofs_on_nodes.nodeDoFConnectivityView();
  auto cn_cv = connectivity_view.cellNode();

  auto command = makeCommand(queue);
  auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
  auto in_node_coord = viewIn(command, node_coord);

  command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh->allCells())
  {
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    NodeLocalId cell_nodes[20];
    for (Int32 i = 0; i < 20; ++i)
      cell_nodes[i] = cn_cv.nodeId(cell_lid, i);

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {
        for (Int32 izeta = 0; izeta < 3; ++izeta) {
          const Real xi = gp[ixi];
          const Real eta = gp[ieta];
          const Real zeta = gp[izeta];
          const Real weight = weights[ixi] * weights[ieta] * weights[izeta];

          // Shape functions 𝐍 for Hexa20
          const RealVector<20> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa20(xi, eta, zeta);
          const auto gp_info = FeOperation3D::computeGradientsAndJacobianHexa20Gpu(cell_lid, cn_cv, in_node_coord, xi, eta, zeta);
          const Real integration_weight = weight * gp_info.det_j;

          for (Int32 i = 0; i < 20; ++i) {
            const NodeLocalId node_lid = cell_nodes[i];
            if (nodes_infos.isOwn(node_lid)) {
              const Real rhs_value = N(i) * qdot * integration_weight;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(
              in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        }
      }
    }
  };
}

/*---------------------------------------------------------------------------*/

void BoundaryConditions3D::
applyConstantSourceToRhsHexa27(Real qdot, const FemDoFsOnNodes& dofs_on_nodes, const VariableNodeReal3& node_coord,
                               VariableDoFReal& rhs_variable_na, IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(mesh);
  ARCANE_CHECK_PTR(queue);

  UnstructuredMeshConnectivityView connectivity_view(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof = dofs_on_nodes.nodeDoFConnectivityView();
  auto cn_cv = connectivity_view.cellNode();

  auto command = makeCommand(queue);
  auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
  auto in_node_coord = viewIn(command, node_coord);

  command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh->allCells())
  {
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    NodeLocalId cell_nodes[27];
    for (Int32 i = 0; i < 27; ++i)
      cell_nodes[i] = cn_cv.nodeId(cell_lid, i);

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {
        for (Int32 izeta = 0; izeta < 3; ++izeta) {
          const Real xi = gp[ixi];
          const Real eta = gp[ieta];
          const Real zeta = gp[izeta];
          const Real weight = weights[ixi] * weights[ieta] * weights[izeta];

          // Shape functions 𝐍 for Hexa27
          const RealVector<27> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa27(xi, eta, zeta);
          const auto gp_info = FeOperation3D::computeGradientsAndJacobianHexa27Gpu(cell_lid, cn_cv, in_node_coord, xi, eta, zeta);
          const Real integration_weight = weight * gp_info.det_j;

          for (Int32 i = 0; i < 27; ++i) {
            const NodeLocalId node_lid = cell_nodes[i];
            if (nodes_infos.isOwn(node_lid)) {
              const Real rhs_value = N(i) * qdot * integration_weight;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(
              in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        }
      }
    }
  };
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
 */
void BoundaryConditions3D::
applyNeumannToRhsTetra4(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                  const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                  IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);
  ARCANE_CHECK_PTR(queue);
  ARCANE_CHECK_PTR(mesh);

  FaceGroup group = bs->getSurface();

  bool scalarNeumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalarNeumann = true;
  }

  UnstructuredMeshConnectivityView connectivity_view;
  connectivity_view.setMesh(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
  auto fn_cv = connectivity_view.faceNode();

  if (scalarNeumann) {
    Real value = std::stod(neumann_str[0].localstr());

    auto command = makeCommand(queue);
    auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);
    command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
    {
      Real area = MeshOperation::computeAreaTria(face_lid, fn_cv, in_node_coord);
      for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
        if (nodes_infos.isOwn(node_lid)) {
          Real rhs_value = value * area / 3.0;
          Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
        }
      }
    };
  }
  else {
    Real value_x = neumann_str[0] != "NULL" ? std::stod(neumann_str[0].localstr()) : 0.0;
    Real value_y = neumann_str[1] != "NULL" ? std::stod(neumann_str[1].localstr()) : 0.0;
    Real value_z = neumann_str[2] != "NULL" ? std::stod(neumann_str[2].localstr()) : 0.0;

    auto command = makeCommand(queue);
    auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);
    FaceInfoListView faces_infos(mesh->faceFamily());
    command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
    {
      Real area = MeshOperation::computeAreaTria(face_lid, fn_cv, in_node_coord);
      Real3 normal = MeshOperation::computeNormalTriangle(face_lid, fn_cv, in_node_coord, faces_infos);
      for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
        if (nodes_infos.isOwn(node_lid)) {
          Real rhs_value = (normal.x * value_x + normal.y * value_y + normal.z * value_z) * area / 3.0;
          Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
        }
      }
    };
  }
}

void BoundaryConditions3D::
applyNeumannToRhsHexa8(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                    const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                    IMesh* mesh, Accelerator::RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);
  ARCANE_CHECK_PTR(queue);
  ARCANE_CHECK_PTR(mesh);

  FaceGroup group = bs->getSurface();
  bool scalarNeumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalarNeumann = true;
  }

  UnstructuredMeshConnectivityView connectivity_view;
  connectivity_view.setMesh(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
  auto fn_cv = connectivity_view.faceNode();

  if (scalarNeumann) {
    Real value = std::stod(neumann_str[0].localstr());

    auto command = makeCommand(queue);
    auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);

    command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
    {
      const Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 };
      const Real w = 1.0;

      // Cache the 4 face nodes locally
      NodeLocalId face_nodes[4];
      Int32 index = 0;
      for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
        if (index < 4)
          face_nodes[index++] = node_lid;
      }

      // Cache node coordinates locally
      Real3 coords[4];
      for (Int32 i = 0; i < 4; ++i)
        coords[i] = in_node_coord[face_nodes[i]];

      // 2x2 Gauss quadrature over quad face
      for (Int32 ixi = 0; ixi < 2; ++ixi) {
        for (Int32 ieta = 0; ieta < 2; ++ieta) {
          Real xi = gp[ixi];
          Real eta = gp[ieta];

          // Shape functions 𝐍 for Quad4
          RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);

          // Shape function derivatives w.r.t. natural coordinates
          const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad4(xi, eta);
          // Tangent vectors: t1 = ∂r/∂ξ,  t2 = ∂r/∂η
          Real3 t1(0.0, 0.0, 0.0);
          Real3 t2(0.0, 0.0, 0.0);
          for (Int32 i = 0; i < 4; ++i) {
            t1.x += reference_gradients.dN_dxi[i] * coords[i].x;
            t1.y += reference_gradients.dN_dxi[i] * coords[i].y;
            t1.z += reference_gradients.dN_dxi[i] * coords[i].z;

            t2.x += reference_gradients.dN_deta[i] * coords[i].x;
            t2.y += reference_gradients.dN_deta[i] * coords[i].y;
            t2.z += reference_gradients.dN_deta[i] * coords[i].z;
          }

          // Normal vector via cross product t1 × t2
          Real3 normal;
          normal.x = t1.y * t2.z - t1.z * t2.y;
          normal.y = t1.z * t2.x - t1.x * t2.z;
          normal.z = t1.x * t2.y - t1.y * t2.x;

          // Surface Jacobian = |t1 × t2|
          Real detJ = math::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

          // Integration weight (Gauss weights are both 1.0, so w*w = 1.0)
          Real integration_weight = w * w * detJ;

          for (Int32 j = 0; j < 4; ++j) {
            NodeLocalId node_lid = face_nodes[j];
            if (nodes_infos.isOwn(node_lid)) {
              Real rhs_value = value * N[j] * integration_weight;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(
              in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        }
      }
    };
  }
  else {
    Real valueX = neumann_str[0] != "NULL" ? std::stod(neumann_str[0].localstr()) : 0.0;
    Real valueY = neumann_str[1] != "NULL" ? std::stod(neumann_str[1].localstr()) : 0.0;
    Real valueZ = neumann_str[2] != "NULL" ? std::stod(neumann_str[2].localstr()) : 0.0;

    auto command = makeCommand(queue);
    auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
    auto in_node_coord = viewIn(command, node_coord);

    command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
    {
      const Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 };
      const Real w = 1.0;

      // Cache the 4 face nodes locally
      NodeLocalId face_nodes[4];
      Int32 index = 0;
      for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
        if (index < 4)
          face_nodes[index++] = node_lid;
      }

      // Cache node coordinates locally
      Real3 coords[4];
      for (Int32 i = 0; i < 4; ++i)
        coords[i] = in_node_coord[face_nodes[i]];

      // 2x2 Gauss quadrature over quad face
      for (Int32 ixi = 0; ixi < 2; ++ixi) {
        for (Int32 ieta = 0; ieta < 2; ++ieta) {
          Real xi = gp[ixi];
          Real eta = gp[ieta];

          // Shape functions 𝐍 for Quad4
          RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);

          // Shape function derivatives w.r.t. natural coordinates
          const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad4(xi, eta);
          // Tangent vectors: t1 = ∂r/∂ξ,  t2 = ∂r/∂η
          Real3 t1(0.0, 0.0, 0.0);
          Real3 t2(0.0, 0.0, 0.0);
          for (Int32 i = 0; i < 4; ++i) {
            t1.x += reference_gradients.dN_dxi[i] * coords[i].x;
            t1.y += reference_gradients.dN_dxi[i] * coords[i].y;
            t1.z += reference_gradients.dN_dxi[i] * coords[i].z;

            t2.x += reference_gradients.dN_deta[i] * coords[i].x;
            t2.y += reference_gradients.dN_deta[i] * coords[i].y;
            t2.z += reference_gradients.dN_deta[i] * coords[i].z;
          }

          // Normal vector via cross product t1 × t2
          Real3 normal;
          normal.x = t1.y * t2.z - t1.z * t2.y;
          normal.y = t1.z * t2.x - t1.x * t2.z;
          normal.z = t1.x * t2.y - t1.y * t2.x;

          // Surface Jacobian = |t1 × t2|
          Real detJ = math::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

          // Unit normal
          Real3 unit_normal;
          unit_normal.x = normal.x / detJ;
          unit_normal.y = normal.y / detJ;
          unit_normal.z = normal.z / detJ;

          // Integration weight
          Real integration_weight = w * w * detJ;

          for (Int32 j = 0; j < 4; ++j) {
            NodeLocalId node_lid = face_nodes[j];
            if (nodes_infos.isOwn(node_lid)) {
              Real rhs_value = (unit_normal.x * valueX + unit_normal.y * valueY + unit_normal.z * valueZ) * N[j] * integration_weight;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(
              in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        }
      }
    };
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BoundaryConditions3D::
applyNeumannToRhsHexa20(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                        const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                        IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);
  ARCANE_CHECK_PTR(mesh);
  ARCANE_CHECK_PTR(queue);

  const FaceGroup group = bs->getSurface();
  const StringConstArrayView neumann_str = bs->getValue();
  const bool scalar_neumann = neumann_str.size() == 1 && neumann_str[0] != "NULL";
  const Real value = scalar_neumann ? std::stod(neumann_str[0].localstr()) : 0.0;
  const Real value_x = !scalar_neumann && neumann_str.size() > 0 && neumann_str[0] != "NULL"
  ? std::stod(neumann_str[0].localstr())
  : 0.0;
  const Real value_y = !scalar_neumann && neumann_str.size() > 1 && neumann_str[1] != "NULL"
  ? std::stod(neumann_str[1].localstr())
  : 0.0;
  const Real value_z = !scalar_neumann && neumann_str.size() > 2 && neumann_str[2] != "NULL"
  ? std::stod(neumann_str[2].localstr())
  : 0.0;

  UnstructuredMeshConnectivityView connectivity_view(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof = dofs_on_nodes.nodeDoFConnectivityView();
  auto fn_cv = connectivity_view.faceNode();

  auto command = makeCommand(queue);
  auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
  auto in_node_coord = viewIn(command, node_coord);

  command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
  {
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    NodeLocalId face_nodes[8];
    Real3 coords[8];
    for (Int32 i = 0; i < 8; ++i) {
      face_nodes[i] = fn_cv.nodeId(face_lid, i);
      coords[i] = in_node_coord[face_nodes[i]];
    }

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real weight = weights[ixi] * weights[ieta];

        RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad8(xi, eta);

        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad8(xi, eta);
        Real3 tangent_xi(0.0, 0.0, 0.0);
        Real3 tangent_eta(0.0, 0.0, 0.0);
        for (Int32 i = 0; i < 8; ++i) {
          tangent_xi += reference_gradients.dN_dxi[i] * coords[i];
          tangent_eta += reference_gradients.dN_deta[i] * coords[i];
        }

        const Real3 normal = math::cross(tangent_xi, tangent_eta);
        const Real surface_jacobian = math::sqrt(math::dot(normal, normal));
        const Real flux = scalar_neumann ? value
                                          : (normal.x * value_x + normal.y * value_y + normal.z * value_z) / surface_jacobian;
        const Real integration_weight = weight * surface_jacobian;

        for (Int32 i = 0; i < 8; ++i) {
          const NodeLocalId node_lid = face_nodes[i];
          if (nodes_infos.isOwn(node_lid)) {
            const Real rhs_value = flux * N[i] * integration_weight;
            Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(
            in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
          }
        }
      }
    }
  };
}

/*---------------------------------------------------------------------------*/

void BoundaryConditions3D::
applyNeumannToRhsHexa27(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                        const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                        IMesh* mesh, RunQueue* queue)
{
  ARCANE_CHECK_PTR(bs);
  ARCANE_CHECK_PTR(mesh);
  ARCANE_CHECK_PTR(queue);

  const FaceGroup group = bs->getSurface();
  const StringConstArrayView neumann_str = bs->getValue();
  const bool scalar_neumann = neumann_str.size() == 1 && neumann_str[0] != "NULL";
  const Real value = scalar_neumann ? std::stod(neumann_str[0].localstr()) : 0.0;
  const Real value_x = !scalar_neumann && neumann_str.size() > 0 && neumann_str[0] != "NULL"
  ? std::stod(neumann_str[0].localstr())
  : 0.0;
  const Real value_y = !scalar_neumann && neumann_str.size() > 1 && neumann_str[1] != "NULL"
  ? std::stod(neumann_str[1].localstr())
  : 0.0;
  const Real value_z = !scalar_neumann && neumann_str.size() > 2 && neumann_str[2] != "NULL"
  ? std::stod(neumann_str[2].localstr())
  : 0.0;

  UnstructuredMeshConnectivityView connectivity_view(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());
  auto node_dof = dofs_on_nodes.nodeDoFConnectivityView();
  auto fn_cv = connectivity_view.faceNode();

  auto command = makeCommand(queue);
  auto in_out_rhs_variable_na = viewInOut(command, rhs_variable_na);
  auto in_node_coord = viewIn(command, node_coord);

  command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
  {
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    NodeLocalId face_nodes[9];
    Real3 coords[9];
    for (Int32 i = 0; i < 9; ++i) {
      face_nodes[i] = fn_cv.nodeId(face_lid, i);
      coords[i] = in_node_coord[face_nodes[i]];
    }

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real weight = weights[ixi] * weights[ieta];

        RealVector<9> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad9(xi, eta);

        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad9(xi, eta);
        Real3 tangent_xi(0.0, 0.0, 0.0);
        Real3 tangent_eta(0.0, 0.0, 0.0);
        for (Int32 i = 0; i < 9; ++i) {
          tangent_xi += reference_gradients.dN_dxi[i] * coords[i];
          tangent_eta += reference_gradients.dN_deta[i] * coords[i];
        }

        const Real3 normal = math::cross(tangent_xi, tangent_eta);
        const Real surface_jacobian = math::sqrt(math::dot(normal, normal));
        const Real flux = scalar_neumann ? value
                                          : (normal.x * value_x + normal.y * value_y + normal.z * value_z) / surface_jacobian;
        const Real integration_weight = weight * surface_jacobian;

        for (Int32 i = 0; i < 9; ++i) {
          const NodeLocalId node_lid = face_nodes[i];
          if (nodes_infos.isOwn(node_lid)) {
            const Real rhs_value = flux * N[i] * integration_weight;
            Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(
            in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
          }
        }
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}; // namespace Arcane::FemUtils::Gpu

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
