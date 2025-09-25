// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArcaneFemFunctionsGpu.cc                                    (C) 2022-2025 */
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

  auto command = Accelerator::makeCommand(queue);
  auto in_out_forced_info = Accelerator::viewInOut(command, linear_system.getForcedInfo());
  auto in_out_rhs_variable = Accelerator::viewInOut(command, linear_system.rhsVariable());

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

  auto command = Accelerator::makeCommand(queue);
  auto in_out_forced_info = Accelerator::viewInOut(command, linear_system.getForcedInfo());
  auto in_out_forced_value = Accelerator::viewInOut(command, linear_system.getForcedValue());
  auto in_out_rhs_variable = Accelerator::viewInOut(command, linear_system.rhsVariable());

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

  auto command = Accelerator::makeCommand(queue);
  auto in_out_elimination_info = Accelerator::viewInOut(command, elimination_helper.getEliminationInfo());
  auto in_out_elimination_value = Accelerator::viewInOut(command, elimination_helper.getEliminationValue());
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
applyConstantSourceToRhs(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
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
applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
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

    auto command = Accelerator::makeCommand(queue);
    auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
    auto in_node_coord = Accelerator::viewIn(command, node_coord);
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

    auto command = Accelerator::makeCommand(queue);
    auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
    auto in_node_coord = Accelerator::viewIn(command, node_coord);
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
applyConstantSourceToRhs(Real qdot, const FemDoFsOnNodes& dofs_on_nodes, const VariableNodeReal3& node_coord,
                         VariableDoFReal& rhs_variable_na, IMesh* mesh, RunQueue* queue)
{
  auto func = [] ARCCORE_HOST_DEVICE(CellLocalId cell_lid,
                                     const IndexedCellNodeConnectivityView& cn_cv,
                                     const Accelerator::VariableNodeReal3InView& in_node_coord) {
    return MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);
  };
  BoundaryConditionsHelpers::applyConstantSourceToRhsBase(func, qdot, queue, rhs_variable_na, mesh, dofs_on_nodes, node_coord);
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
applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
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

    auto command = Accelerator::makeCommand(queue);
    auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
    auto in_node_coord = Accelerator::viewIn(command, node_coord);
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

    auto command = Accelerator::makeCommand(queue);
    auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
    auto in_node_coord = Accelerator::viewIn(command, node_coord);
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}; // namespace Arcane::FemUtils::Gpu

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
