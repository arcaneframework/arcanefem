// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* InternalBodyForce.h                                 (C) 2000-2026         */
/*                                                                           */
/* Contains functions to compute and assemble source term contribution to RHS*/
/* corresponding to the internal force term                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies nonlinear internal body force term to RHS vector of
 * the linear system for Von Mises plasticity law.
 * 
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::
_applyInternalBodyForce(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  auto use_gpu = options()->linearSystem.serviceName() == "HypreLinearSystem" ||
    options()->linearSystem.serviceName() == "PetscLinearSystem";

  if (use_gpu && m_use_gpu_functions) {
    auto queue = subDomain()->acceleratorMng()->defaultQueue();
    auto mesh_ptr = mesh();
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {
        if (m_nodes_per_cell == 4)
          _applyInternalBodyForceQuad4Cpu(rhs_values, node_dof);  // TODO: Implement GPU version for Quad4
        else if (m_nodes_per_cell == 8)
          _applyInternalBodyForceQuad8Cpu(rhs_values, node_dof);  // TODO: Implement GPU version for Quad8
        else
          _applyInternalBodyForceQuad9Cpu(rhs_values, node_dof);  // TODO: Implement GPU version for Quad9
      }
      else {
        _applyInternalBodyForceTria3Gpu(rhs_values, m_dofs_on_nodes, m_node_coord, mesh_ptr, queue);
      }
    }
    else {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Not IMPLEMENTED");
      }
      else {
        ARCANE_FATAL("Not IMPLEMENTED");
      }
    }
  }
  else {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {
        if (m_nodes_per_cell == 4)
          _applyInternalBodyForceQuad4Cpu(rhs_values, node_dof);
        else if (m_nodes_per_cell == 8)
          _applyInternalBodyForceQuad8Cpu(rhs_values, node_dof);
        else
          _applyInternalBodyForceQuad9Cpu(rhs_values, node_dof);
      }
      else {
        _applyInternalBodyForceTria3Cpu(rhs_values, node_dof);
      }
    }
    else {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Not IMPLEMENTED");
      }
      else {
        ARCANE_FATAL("Not IMPLEMENTED");
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies internal body force for Quad4 elements on CPU.
 *
 * This function computes the contribution of internal body forces to the RHS
 * vector for Quad4 elements. It iterates over all cells, evaluates the shape
 * functions at Gauss points, computes the Jacobian determinant, and updates
 * the RHS vector accordingly.
 *
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 *                 degrees of freedom (DoFs).
 */
/*---------------------------------------------------------------------------*/

inline void FemModuleElastoplasticity::
_applyInternalBodyForceQuad4Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 };

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Int8 iGP = 0;
    for (Int8 ixi = 0; ixi < 2; ++ixi) {
      for (Int8 ieta = 0; ieta < 2; ++ieta) {
        const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, gp[ixi], gp[ieta]);
        const Real sigma_xx = m_sigma_gp(cell, iGP, 0);
        const Real sigma_yy = m_sigma_gp(cell, iGP, 1);
        const Real sigma_xy = m_sigma_gp(cell, iGP, 2);

        for (Int32 i = 0; i < 4; ++i) {
          const Node node = cell.node(i);
          if (!node.isOwn())
            continue;
          rhs_values[node_dof.dofId(node, 0)] -= gp_info.det_j * (sigma_xx * gp_info.dN_dx(i) + M_SQRT1_2 * sigma_xy * gp_info.dN_dy(i));
          rhs_values[node_dof.dofId(node, 1)] -= gp_info.det_j * (sigma_yy * gp_info.dN_dy(i) + M_SQRT1_2 * sigma_xy * gp_info.dN_dx(i));
        }
        ++iGP;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies internal body force for Quad8 elements on CPU.
 *
 * This function computes the contribution of internal body forces to the RHS
 * vector for Quad8 elements. It iterates over all cells, evaluates the shape
 * functions at Gauss points, computes the Jacobian determinant, and updates
 * the RHS vector accordingly.
 *
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 *                 degrees of freedom (DoFs).
 */
/*---------------------------------------------------------------------------*/

inline void FemModuleElastoplasticity::
_applyInternalBodyForceQuad8Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Int8 iGP = 0;
    for (Int8 ixi = 0; ixi < 3; ++ixi) {
      for (Int8 ieta = 0; ieta < 3; ++ieta) {
        const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad8(cell, m_node_coord, gp[ixi], gp[ieta]);
        const Real integration_weight = gp_info.det_j * weights[ixi] * weights[ieta];
        const Real sigma_xx = m_sigma_gp(cell, iGP, 0);
        const Real sigma_yy = m_sigma_gp(cell, iGP, 1);
        const Real sigma_xy = m_sigma_gp(cell, iGP, 2);

        for (Int32 i = 0; i < 8; ++i) {
          const Node node = cell.node(i);
          if (!node.isOwn())
            continue;
          rhs_values[node_dof.dofId(node, 0)] -= integration_weight * (sigma_xx * gp_info.dN_dx(i) + M_SQRT1_2 * sigma_xy * gp_info.dN_dy(i));
          rhs_values[node_dof.dofId(node, 1)] -= integration_weight * (sigma_yy * gp_info.dN_dy(i) + M_SQRT1_2 * sigma_xy * gp_info.dN_dx(i));
        }
        ++iGP;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies internal body force for Quad9 elements on CPU.
 *
 * This function computes the contribution of internal body forces to the RHS
 * vector for Quad9 elements. It iterates over all cells, evaluates the shape
 * functions at Gauss points, computes the Jacobian determinant, and updates
 * the RHS vector accordingly.
 *
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 *                 degrees of freedom (DoFs).
 */
/*---------------------------------------------------------------------------*/

inline void FemModuleElastoplasticity::
_applyInternalBodyForceQuad9Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Int8 iGP = 0;
    for (Int8 ixi = 0; ixi < 3; ++ixi) {
      for (Int8 ieta = 0; ieta < 3; ++ieta) {
        const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad9(cell, m_node_coord, gp[ixi], gp[ieta]);
        const Real integration_weight = gp_info.det_j * weights[ixi] * weights[ieta];
        const Real sigma_xx = m_sigma_gp(cell, iGP, 0);
        const Real sigma_yy = m_sigma_gp(cell, iGP, 1);
        const Real sigma_xy = m_sigma_gp(cell, iGP, 2);

        for (Int32 i = 0; i < 9; ++i) {
          const Node node = cell.node(i);
          if (!node.isOwn())
            continue;
          rhs_values[node_dof.dofId(node, 0)] -= integration_weight * (sigma_xx * gp_info.dN_dx(i) + M_SQRT1_2 * sigma_xy * gp_info.dN_dy(i));
          rhs_values[node_dof.dofId(node, 1)] -= integration_weight * (sigma_yy * gp_info.dN_dy(i) + M_SQRT1_2 * sigma_xy * gp_info.dN_dx(i));
        }
        ++iGP;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
* @brief Evaluates and assembles the rhs vector components corresponding
*        to integration of the internal force for a triangular element (ℙ1 FE)
*        using CPUs.
*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline RealVector<6>
computeInternalBodyForceTria3Base(Real3 dxu,
                                          Real3 dyu,
                                          Real area,
                                          RealVector<3> sigma_2d)
{

  RealVector<6> epsxx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  RealVector<6> epsyy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };
  RealVector<6> epsxy = { dyu[0], dxu[0], dyu[1], dxu[1], dyu[2], dxu[2] };
  epsxy = 0.70710678118654746172 * epsxy;

  RealVector<6> rhs = - area * (sigma_2d[0] * epsxx + sigma_2d[1] * epsyy + sigma_2d[2] * epsxy);

  return rhs;
}

inline void FemModuleElastoplasticity::
_applyInternalBodyForceTria3Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  info() << "[ArcaneFem-Info] Started module  _applyInternalBodyForceTria3Cpu()";

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
    Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
    Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

    Int8 iGP = 0; // for tria P1 elements nGP=1
    Real sigma_xx = m_sigma_gp(cell , iGP, 0);
    Real sigma_yy = m_sigma_gp(cell , iGP, 1);
    Real sigma_xy = m_sigma_gp(cell , iGP, 2);

    RealVector<6> rhs = computeInternalBodyForceTria3Base(dxu, dyu, area, { sigma_xx, sigma_yy, sigma_xy });

    rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
    rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
    rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(2);
    rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(3);
    rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(4);
    rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(5);
  }
}

/*---------------------------------------------------------------------------*/
/**
* @brief Evaluates and assembles the rhs vector components corresponding
*        to integration of the internal force for a triangular element (ℙ1 FE)
*        using GPUs.
*
*/
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::
_applyInternalBodyForceTria3Gpu(VariableDoFReal& rhs_values,
                                        const FemDoFsOnNodes& dofs_on_nodes,
                                        const VariableNodeReal3& node_coord,
                                        IMesh* mesh, RunQueue* queue)
{
  info() << "[ArcaneFem-Info] Started module  _applyInternalBodyForceTria3Gpu()";
  ARCANE_CHECK_PTR(queue);
  ARCANE_CHECK_PTR(mesh);

  UnstructuredMeshConnectivityView connectivity_view;
  connectivity_view.setMesh(mesh);
  NodeInfoListView nodes_infos(mesh->nodeFamily());

  auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
  auto cn_cv = connectivity_view.cellNode();

  auto command = Accelerator::makeCommand(queue);

  auto in_out_rhs_values = Accelerator::viewInOut(command, rhs_values);
  auto in_node_coord = Accelerator::viewIn(command, node_coord);

  auto in_sigma_gp = Accelerator::viewIn(command, m_sigma_gp);

  command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh->allCells())
  {
    Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);
    Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
    Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);

    Int8 iGP = 0; // for tria P1 elements nGP=1
    Real sigma_xx = in_sigma_gp(cell_lid , iGP, 0);
    Real sigma_yy = in_sigma_gp(cell_lid , iGP, 1);
    Real sigma_xy = in_sigma_gp(cell_lid , iGP, 2);

    RealVector<6> rhs = computeInternalBodyForceTria3Base(dxu, dyu, area, { sigma_xx, sigma_yy, sigma_xy });

    NodeLocalId cell_nodes[3];
    Int32 index = 0;
    for (NodeLocalId node_lid : cn_cv.nodes(cell_lid)) {
      if (index < 3) {
        cell_nodes[index++] = node_lid;
      }
    }

    for (Int8 i = 0; i < 3; ++i) {
      NodeLocalId node_lid = cell_nodes[i];
      if (nodes_infos.isOwn(node_lid)) {
        Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_values[node_dof.dofId(node_lid, 0)], rhs(2*i));
        Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_values[node_dof.dofId(node_lid, 1)], rhs(2*i + 1));
      }
    }
  };

}
