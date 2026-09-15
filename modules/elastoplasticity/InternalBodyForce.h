// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* InternalBodyForce.h                                 (C) 2000-2026 */
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
        ARCANE_FATAL("Not IMPLEMENTED");
      } else {
        _applyInternalBodyForceTria3Gpu(rhs_values, m_dofs_on_nodes, m_node_coord, mesh_ptr, queue);
      }
    } else {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Not IMPLEMENTED");
      } else {
        ARCANE_FATAL("Not IMPLEMENTED");
      }
    }
  } else {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Not IMPLEMENTED");
      } else {
        _applyInternalBodyForceTria3Cpu(rhs_values, node_dof);
      }
    } else {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Not IMPLEMENTED");
      } else {
        ARCANE_FATAL("Not IMPLEMENTED");
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝑈,𝐯) = ∫∫ σ(𝑈):ε(𝐯)dΩ        with  𝑈 = (𝑈𝑥,𝑈𝑦) and 𝐯 = (𝑣𝑥,𝑣𝑦)
 *   σ(𝑈) is known stress tensor    with  σᵢⱼ = Cᵢⱼₖₗεₖₗ
 *   ε(𝐯) is strain tensor          with  εᵢⱼ = 0.5 (∂𝑣ᵢ/∂xⱼ + ∂𝑣ⱼ/∂xᵢ)
 *
 *   the linear integral expands to
 *
 *      a(𝑈,𝐯) = ∫∫ [σ_𝑥𝑥 ε_𝑥𝑥 + σ_𝑦𝑦 ε_𝑦𝑦 + 2 σ_𝑥𝑦 ε_𝑥𝑦]dΩ
 *
 *   this further expands to
 *
 *      a(𝐮,𝐯) =   ∫∫ C11 ∂𝑈𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + C12 ∂𝑈𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + C13 (∂𝑈𝑦/∂𝑥 + ∂𝑈𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ C12 ∂𝑈𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + C22 ∂𝑈𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + C23 (∂𝑈𝑦/∂𝑥 + ∂𝑈𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ C13 ∂𝑈𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C23 ∂𝑈𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C33 (∂𝑈𝑦/∂𝑥 + ∂𝑈𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
 *
 *
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