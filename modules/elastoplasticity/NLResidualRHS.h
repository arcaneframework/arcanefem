// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ResidualRHS.h                                                (C) 2000-2026 */
/*                                                                           */
/* Contains functions to compute and assemble source term contribution to RHS*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies nonlinear residual term to RHS vector of the linear system.
 * 
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 */
/*---------------------------------------------------------------------------*/

inline void FemModuleElastoplasticity::
_applyResidualRHS(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  if (mesh()->dimension() == 2)
    if (m_hex_quad_mesh)
      _applyResidualRHSQuad4(rhs_values, node_dof);
    else
      _applyResidualRHSTria3(rhs_values, node_dof);

  if (mesh()->dimension() == 3)
    if (m_hex_quad_mesh)
      _applyResidualRHSHexa8(rhs_values, node_dof);
    else
      _applyResidualRHSTetra4(rhs_values, node_dof);
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

inline void FemModuleElastoplasticity::
_applyResidualRHSTria3(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
    Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
    Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

    //----------------------------------------------------------------------
    //  ∫∫∫ (σ(𝑈)ε(𝐯) = ∫∫∫ (ε(𝑈):𝐶) ε(𝐯)
    //----------------------------------------------------------------------

    RealVector<6> epsxx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
    RealVector<6> epsyy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };
    RealVector<6> epsxy = { dyu[0], dxu[0], dyu[1], dxu[1], dyu[2], dxu[2] };

    Real3x3 grad_U = ArcaneFemFunctions::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_U);

    Real epsxx_U = grad_U(0, 0);
    Real epsyy_U = grad_U(1, 1);
    Real epsxy_U = grad_U(0, 1) + grad_U(1, 0);

    Real sigmaxx_U = m_C_2d_cell(cell, 0, 0) * epsxx_U
                    + m_C_2d_cell(cell, 0, 1) * epsyy_U
                    + m_C_2d_cell(cell, 0, 2) * epsxy_U;
    Real sigmayy_U = m_C_2d_cell(cell, 1, 0) * epsxx_U
                    + m_C_2d_cell(cell, 1, 1) * epsyy_U
                    + m_C_2d_cell(cell, 1, 2) * epsxy_U;
    Real sigmaxy_U = m_C_2d_cell(cell, 2, 0) * epsxx_U
                    + m_C_2d_cell(cell, 2, 1) * epsyy_U
                    + m_C_2d_cell(cell, 2, 2) * epsxy_U;

    RealVector<6> rhs = - area * (sigmaxx_U * epsxx + sigmayy_U * epsyy + sigmaxy_U * epsxy);

    //----------------------------------------------------------------------
    //  ∫∫∫ (ε(𝑈):𝐶)ε(𝐯) = ∫∫∫ (𝑈 ε(𝘶):𝐶) ε(𝐯)
    //----------------------------------------------------------------------
    // RealVector<6> Uk =  { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y,
    //                      m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y,
    //                      m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y };
    // RealVector<6> rhs = - area *
    //                       ( Uk * ((m_C_2d_cell(cell, 0, 0) * epsxx
    //                               + m_C_2d_cell(cell, 0, 1) * epsyy
    //                               + m_C_2d_cell(cell, 0, 2) * epsxy) ^ epsxx)
    //                       + Uk * ((m_C_2d_cell(cell, 1, 0) * epsxx
    //                               + m_C_2d_cell(cell, 1, 1) * epsyy
    //                               + m_C_2d_cell(cell, 1, 2) * epsxy) ^ epsyy)
    //                       + Uk * ((m_C_2d_cell(cell, 2, 0) * epsxx
    //                               + m_C_2d_cell(cell, 2, 1) * epsyy
    //                               + m_C_2d_cell(cell, 2, 2) * epsxy) ^ epsxy));  // To verify mathematically

    rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
    rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
    rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(2);
    rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(3);
    rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(4);
    rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(5);
  }
}

void FemModuleElastoplasticity::
_applyResidualRHSQuad4(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    // Initialize RHS contributions (2 dof/node for 4 quad nodes)
    // Real rhs_x_contributions[4] = { 0., 0., 0., 0. };
    // Real rhs_y_contributions[4] = { 0., 0., 0., 0. };

    // 2x2 Gauss integration for quadrilateral element
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 };
    constexpr Real w = 1.0;

    for (Int8 ixi = 0; ixi < 2; ++ixi) {
      for (Int8 ieta = 0; ieta < 2; ++ieta) {

        // Get the coordinates of the Gauss point
        Real xi = gp[ixi]; // Get the ξ
        Real eta = gp[ieta]; // Get the η
        Real weight = w * w; // Weight

        // Shape functions  𝐍 for Quad4
        RealVector<4> N = ArcaneFemFunctions::FeOperation2D::computeShapeFunctionsQuad4(xi, eta);

        // compute the det(Jacobian)
        auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, xi, eta);
        const RealVector<4>& dxu = gp_info.dN_dx;
        const RealVector<4>& dyu = gp_info.dN_dy;
        const Real detJ = gp_info.det_j;

        // compute integration weight
        Real integration_weight = weight * detJ;

        RealVector<8> epsxx = { dxu(0), 0., dxu(1), 0., dxu(2), 0., dxu(3), 0. };
        RealVector<8> epsyy = { 0., dyu(0), 0., dyu(1), 0., dyu(2), 0., dyu(3) };
        RealVector<8> epsxy = { dyu(0), dxu(0), dyu(1), dxu(1), dyu(2), dxu(2), dyu(3), dxu(3) };

        Real3x3 grad_U = ArcaneFemFunctions::FeOperation2D::computeGradientQuad4(cell, m_node_coord, m_U, xi, eta);

        Real epsxx_U = grad_U(0, 0);
        Real epsyy_U = grad_U(1, 1);
        Real epsxy_U = grad_U(0, 1) + grad_U(1, 0);

        Real sigmaxx_U =  m_C_2d_cell(cell, 0, 0) * epsxx_U
                        + m_C_2d_cell(cell, 0, 1) * epsyy_U
                        + m_C_2d_cell(cell, 0, 2) * epsxy_U;
        Real sigmayy_U =  m_C_2d_cell(cell, 1, 0) * epsxx_U
                        + m_C_2d_cell(cell, 1 ,1) * epsyy_U
                        + m_C_2d_cell(cell, 1, 2) * epsxy_U;
        Real sigmaxy_U =  m_C_2d_cell(cell, 2, 0) * epsxx_U
                        + m_C_2d_cell(cell, 2, 1) * epsyy_U
                        + m_C_2d_cell(cell, 2, 2) * epsxy_U;

        RealVector<8> rhs = - integration_weight * ( sigmaxx_U * epsxx + sigmayy_U * epsyy + sigmaxy_U * epsxy);

        // RealVector<8> Uk = {m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y,
        //                     m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y,
        //                     m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y,
        //                     m_U[cell.nodeId(3)].x, m_U[cell.nodeId(3)].y};
        //
        // RealVector<8> rhs_1 = - integration_weight *
        //               ( Uk * ((m_C_2d_cell(cell, 0, 0) * epsxx
        //                        + m_C_2d_cell(cell, 0, 1) * epsyy
        //                        + m_C_2d_cell(cell, 0, 2) * epsxy) ^ epsxx)
        //                + Uk * ((m_C_2d_cell(cell, 1, 0) * epsxx
        //                        + m_C_2d_cell(cell, 1, 1) * epsyy
        //                        + m_C_2d_cell(cell, 1, 2) * epsxy) ^ epsyy)
        //                + Uk * ((m_C_2d_cell(cell, 2, 0) * epsxx
        //                        + m_C_2d_cell(cell, 2, 1) * epsyy
        //                        + m_C_2d_cell(cell, 2, 2) * epsxy) ^ epsxy)); // verified
        //
        // Real diff = 0.0;
        // Real maxx = 0.0;
        //
        // for (Int8 i = 0; i < 8; ++i) {
        //   Real err = math::abs(rhs[i] - rhs_1[i]);
        //   diff += err;
        //   maxx = math::max(maxx, err);
        // }
        // diff /= 8.0;
        // info() << "Let us check: dif = " << diff << " max = " << maxx;

        rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
        rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
        rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(2);
        rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(3);
        rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(4);
        rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(5);
        rhs_values[node_dof.dofId(cell.nodeId(3), 0)] += rhs(6);
        rhs_values[node_dof.dofId(cell.nodeId(3), 1)] += rhs(7);

      }
    }

    // Add contributions to global RHS
    // for (Int8 a = 0; a < 4; ++a) {
    //   Node node = cell.node(a);
    //   if (node.isOwn()) {
    //     rhs_values[node_dof.dofId(node, 0)] += rhs_x_contributions[a];
    //     rhs_values[node_dof.dofId(node, 1)] += rhs_y_contributions[a];
    //   }
    // }
  }
}

inline void FemModuleElastoplasticity::
_applyResidualRHSTetra4(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);
    Real4 dxu = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
    Real4 dyu = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
    Real4 dzu = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);

    //----------------------------------------------------------------------
    //  ∫∫∫ (σ(𝑈)ε(𝐯) = ∫∫∫ (ε(𝑈):𝐶) ε(𝐯)
    //----------------------------------------------------------------------

    RealVector<12> epsxx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
    RealVector<12> epsyy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
    RealVector<12> epszz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

    RealVector<12> epsyz = { 0., dzu[0], dyu[0],    0., dzu[1], dyu[1],    0., dzu[2], dyu[2],    0., dzu[3], dyu[3] };
    RealVector<12> epszx = { dzu[0], 0., dxu[0],    dzu[1], 0., dxu[1],    dzu[2], 0., dxu[2],    dzu[3], 0., dxu[3] };
    RealVector<12> epsxy = { dyu[0], dxu[0], 0.,    dyu[1], dxu[1], 0.,    dyu[2], dxu[2], 0.,    dyu[3], dxu[3], 0. };


    Real3x3 grad_U = ArcaneFemFunctions::FeOperation3D::computeGradientTetra4(cell, m_node_coord, m_U);

    Real epsxx_U = grad_U(0, 0);
    Real epsyy_U = grad_U(1, 1);
    Real epszz_U = grad_U(2, 2);
    Real epsyz_U = grad_U(1, 2) + grad_U(2, 1);
    Real epszx_U = grad_U(0, 2) + grad_U(2, 0);
    Real epsxy_U = grad_U(0, 1) + grad_U(1, 0);

    Real sigmaxx_U =  m_C_3d_cell(cell, 0, 0) * epsxx_U
                    + m_C_3d_cell(cell, 0, 1) * epsyy_U
                    + m_C_3d_cell(cell, 0, 2) * epszz_U
                    + m_C_3d_cell(cell, 0, 3) * epsyz_U
                    + m_C_3d_cell(cell, 0, 4) * epszx_U
                    + m_C_3d_cell(cell, 0, 5) * epsxy_U;
    Real sigmayy_U =  m_C_3d_cell(cell, 1, 0) * epsxx_U
                    + m_C_3d_cell(cell, 1, 1) * epsyy_U
                    + m_C_3d_cell(cell, 1, 2) * epszz_U
                    + m_C_3d_cell(cell, 1, 3) * epsyz_U
                    + m_C_3d_cell(cell, 1, 4) * epszx_U
                    + m_C_3d_cell(cell, 1, 5) * epsxy_U;
    Real sigmazz_U =  m_C_3d_cell(cell, 2, 0) * epsxx_U
                    + m_C_3d_cell(cell, 2, 1) * epsyy_U
                    + m_C_3d_cell(cell, 2, 2) * epszz_U
                    + m_C_3d_cell(cell, 2, 3) * epsyz_U
                    + m_C_3d_cell(cell, 2, 4) * epszx_U
                    + m_C_3d_cell(cell, 2, 5) * epsxy_U;
    Real sigmayz_U =  m_C_3d_cell(cell, 3, 0) * epsxx_U
                    + m_C_3d_cell(cell, 3, 1) * epsyy_U
                    + m_C_3d_cell(cell, 3, 2) * epszz_U
                    + m_C_3d_cell(cell, 3, 3) * epsyz_U
                    + m_C_3d_cell(cell, 3, 4) * epszx_U
                    + m_C_3d_cell(cell, 3, 5) * epsxy_U;
    Real sigmazx_U =  m_C_3d_cell(cell, 4, 0) * epsxx_U
                    + m_C_3d_cell(cell, 4, 1) * epsyy_U
                    + m_C_3d_cell(cell, 4, 2) * epszz_U
                    + m_C_3d_cell(cell, 4, 3) * epsyz_U
                    + m_C_3d_cell(cell, 4, 4) * epszx_U
                    + m_C_3d_cell(cell, 4, 5) * epsxy_U;
    Real sigmaxy_U =  m_C_3d_cell(cell, 5, 0) * epsxx_U
                    + m_C_3d_cell(cell, 5, 1) * epsyy_U
                    + m_C_3d_cell(cell, 5, 2) * epszz_U
                    + m_C_3d_cell(cell, 5, 3) * epsyz_U
                    + m_C_3d_cell(cell, 5, 4) * epszx_U
                    + m_C_3d_cell(cell, 5, 5) * epsxy_U;


    RealVector<12> rhs_1 = - volume * ( sigmaxx_U * epsxx + sigmayy_U * epsyy + sigmazz_U * epszz
                                      + sigmayz_U * epsyz + sigmazx_U * epszx + sigmaxy_U * epsxy);

    //----------------------------------------------------------------------
    //  ∫∫∫ (ε(𝑈):𝐶)ε(𝐯) = ∫∫∫ (𝑈 ε(𝘶):𝐶) ε(𝐯)
    //----------------------------------------------------------------------
    RealVector<12> Uk = { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y, m_U[cell.nodeId(0)].z,
                          m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y, m_U[cell.nodeId(1)].z,
                          m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y, m_U[cell.nodeId(2)].z,
                          m_U[cell.nodeId(3)].x, m_U[cell.nodeId(3)].y, m_U[cell.nodeId(3)].z};

    RealVector<12> rhs = - volume *
                           (Uk * ((m_C_3d_cell(cell, 0, 0) * epsxx
                                  + m_C_3d_cell(cell, 0, 1) * epsyy
                                  + m_C_3d_cell(cell, 0, 2) * epszz
                                  + m_C_3d_cell(cell, 0, 3) * epsyz
                                  + m_C_3d_cell(cell, 0, 4) * epszx
                                  + m_C_3d_cell(cell, 0, 5) * epsxy) ^ epsxx)
                          + Uk * ((m_C_3d_cell(cell, 0, 1) * epsxx
                                  + m_C_3d_cell(cell, 1, 1) * epsyy
                                  + m_C_3d_cell(cell, 1, 2) * epszz
                                  + m_C_3d_cell(cell, 1, 3) * epsyz
                                  + m_C_3d_cell(cell, 1, 4) * epszx
                                  + m_C_3d_cell(cell, 1, 5) * epsxy) ^ epsyy)
                          + Uk * ((m_C_3d_cell(cell, 0, 2) * epsxx
                                  + m_C_3d_cell(cell, 1, 2) * epsyy
                                  + m_C_3d_cell(cell, 2, 2) * epszz
                                  + m_C_3d_cell(cell, 2, 3) * epsyz
                                  + m_C_3d_cell(cell, 2, 4) * epszx
                                  + m_C_3d_cell(cell, 2, 5) * epsxy) ^ epszz)
                          + Uk * ((m_C_3d_cell(cell, 0, 3) * epsxx
                                  + m_C_3d_cell(cell, 1, 3) * epsyy
                                  + m_C_3d_cell(cell, 2, 3) * epszz
                                  + m_C_3d_cell(cell, 3, 3) * epsyz
                                  + m_C_3d_cell(cell, 3, 4) * epszx
                                  + m_C_3d_cell(cell, 3, 5) * epsxy) ^ epsyz)
                          + Uk * ((m_C_3d_cell(cell, 0, 4) * epsxx
                                  + m_C_3d_cell(cell, 1, 4) * epsyy
                                  + m_C_3d_cell(cell, 2, 4) * epszz
                                  + m_C_3d_cell(cell, 3, 4) * epsyz
                                  + m_C_3d_cell(cell, 4, 4) * epszx
                                  + m_C_3d_cell(cell, 4, 5) * epsxy) ^ epszx)
                          + Uk * ((m_C_3d_cell(cell, 0, 5) * epsxx
                                  + m_C_3d_cell(cell, 1, 5) * epsyy
                                  + m_C_3d_cell(cell, 2, 5) * epszz
                                  + m_C_3d_cell(cell, 3, 5) * epsyz
                                  + m_C_3d_cell(cell, 4, 5) * epszx
                                  + m_C_3d_cell(cell, 5, 5) * epsxy) ^ epsxy));

    Real diff = 0.0;
    Real maxx = 0.0;

    for (Int8 i = 0; i < 12; ++i) {
      Real err = math::abs(rhs[i] - rhs_1[i]);
      diff += err;
      maxx = math::max(maxx, err);
    }
    diff /= 12.0;
    info() << "Let us check: diff = " << diff << " maxx = " << maxx;

    rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
    rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
    rhs_values[node_dof.dofId(cell.nodeId(0), 2)] += rhs(2);
    rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(3);
    rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(4);
    rhs_values[node_dof.dofId(cell.nodeId(1), 2)] += rhs(5);
    rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(6);
    rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(7);
    rhs_values[node_dof.dofId(cell.nodeId(2), 2)] += rhs(8);
    rhs_values[node_dof.dofId(cell.nodeId(3), 0)] += rhs(9);
    rhs_values[node_dof.dofId(cell.nodeId(3), 1)] += rhs(10);
    rhs_values[node_dof.dofId(cell.nodeId(3), 2)] += rhs(11);
  }
}

void FemModuleElastoplasticity::
_applyResidualRHSHexa8(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    // Initialize RHS contributions (2 dof/node for 8 hexa nodes)
    // Real rhs_x_contributions[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };
    // Real rhs_y_contributions[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };
    // Real rhs_z_contributions[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };

    // 2x2 Gauss integration for quadrilateral element
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 };
    constexpr Real w = 1.0;

    for (Int8 ixi = 0; ixi < 2; ++ixi) {
      for (Int8 ieta = 0; ieta < 2; ++ieta) {
        for (Int8 izeta = 0; izeta < 2; ++izeta) {

          // Get the coordinates of the Gauss point
          Real xi = gp[ixi]; // Get the ξ
          Real eta = gp[ieta]; // Get the η
          Real zeta = gp[izeta]; // Get the ζ
          Real weight = w * w * w; // Weight

          // Shape functions  𝐍 for Hexa8
          RealVector<8> N = ArcaneFemFunctions::FeOperation3D::computeShapeFunctionsHexa8(xi, eta, zeta);

          // compute the det(Jacobian)
          auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa8(cell, m_node_coord, xi, eta, zeta);
          const RealVector<8>& dxu = gp_info.dN_dx;
          const RealVector<8>& dyu = gp_info.dN_dy;
          const RealVector<8>& dzu = gp_info.dN_dz;
          const Real detJ = gp_info.det_j;

          // compute integration weight
          Real integration_weight = weight * detJ;

          RealVector<24> epsxx = { dxu(0), 0., 0.,    dxu(1), 0., 0.,    dxu(2), 0., 0.,    dxu(3), 0., 0.,
                                   dxu(4), 0., 0.,    dxu(5), 0., 0.,    dxu(6), 0., 0.,    dxu(7), 0., 0. };

          RealVector<24> epsyy = { 0., dyu(0), 0.,    0., dyu(1), 0.,    0., dyu(2), 0.,    0., dyu(3), 0.,
                                   0., dyu(4), 0.,    0., dyu(5), 0.,    0., dyu(6), 0.,    0., dyu(7), 0. };

          RealVector<24> epszz = { 0., 0., dzu(0),    0., 0., dzu(1),    0., 0., dzu(2),    0., 0., dzu(3),
                                   0., 0., dzu(4),    0., 0., dzu(5),    0., 0., dzu(6),    0., 0., dzu(7) };

          RealVector<24> epsyz = { 0., dzu(0), dyu(0),    0., dzu(1), dyu(1),
                                   0., dzu(2), dyu(2),    0., dzu(3), dyu(3),
                                   0., dzu(4), dyu(4),    0., dzu(5), dyu(5),
                                   0., dzu(6), dyu(6),    0., dzu(7), dyu(7) };

          RealVector<24> epszx = { dzu(0), 0., dxu(0),    dzu(1), 0., dxu(1),
                                   dzu(2), 0., dxu(2),    dzu(3), 0., dxu(3),
                                   dzu(4), 0., dxu(4),    dzu(5), 0., dxu(5),
                                   dzu(6), 0., dxu(6),    dzu(7), 0., dxu(7) };

          RealVector<24> epsxy = { dyu(0), dxu(0), 0.,    dyu(1), dxu(1), 0.,
                                   dyu(2), dxu(2), 0.,    dyu(3), dxu(3), 0.,
                                   dyu(4), dxu(4), 0.,    dyu(5), dxu(5), 0.,
                                   dyu(6), dxu(6), 0.,    dyu(7), dxu(7), 0. };

          RealVector<24> Uk = { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y, m_U[cell.nodeId(0)].z,
                                m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y, m_U[cell.nodeId(1)].z,
                                m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y, m_U[cell.nodeId(2)].z,
                                m_U[cell.nodeId(3)].x, m_U[cell.nodeId(3)].y, m_U[cell.nodeId(3)].z,
                                m_U[cell.nodeId(4)].x, m_U[cell.nodeId(4)].y, m_U[cell.nodeId(4)].z,
                                m_U[cell.nodeId(5)].x, m_U[cell.nodeId(5)].y, m_U[cell.nodeId(5)].z,
                                m_U[cell.nodeId(6)].x, m_U[cell.nodeId(6)].y, m_U[cell.nodeId(6)].z,
                                m_U[cell.nodeId(7)].x, m_U[cell.nodeId(7)].y, m_U[cell.nodeId(7)].z };

          RealVector<24> rhs = - integration_weight *
                        ( Uk * ((  m_C_3d_cell(cell, 0, 0) * epsxx
                                 + m_C_3d_cell(cell, 0, 1) * epsyy
                                 + m_C_3d_cell(cell, 0, 2) * epszz
                                 + m_C_3d_cell(cell, 0, 3) * epsyz
                                 + m_C_3d_cell(cell, 0, 4) * epszx
                                 + m_C_3d_cell(cell, 0, 5) * epsxy) ^ epsxx)
                         + Uk * (( m_C_3d_cell(cell, 1, 0) * epsxx
                                 + m_C_3d_cell(cell, 1, 1) * epsyy
                                 + m_C_3d_cell(cell, 1, 2) * epszz
                                 + m_C_3d_cell(cell, 1, 3) * epsyz
                                 + m_C_3d_cell(cell, 1, 4) * epszx
                                 + m_C_3d_cell(cell, 1, 5) * epsxy) ^ epsyy)
                         + Uk * (( m_C_3d_cell(cell, 2, 0) * epsxx
                                 + m_C_3d_cell(cell, 2, 1) * epsyy
                                 + m_C_3d_cell(cell, 2, 2) * epszz
                                 + m_C_3d_cell(cell, 2, 3) * epsyz
                                 + m_C_3d_cell(cell, 2, 4) * epszx
                                 + m_C_3d_cell(cell, 2, 5) * epsxy) ^ epszz)
                         + Uk * (( m_C_3d_cell(cell, 3, 0) * epsxx
                                 + m_C_3d_cell(cell, 3, 1) * epsyy
                                 + m_C_3d_cell(cell, 3, 2) * epszz
                                 + m_C_3d_cell(cell, 3, 3) * epsyz
                                 + m_C_3d_cell(cell, 3, 4) * epszx
                                 + m_C_3d_cell(cell, 3, 5) * epsxy) ^ epsyz)
                         + Uk * (( m_C_3d_cell(cell, 4, 0) * epsxx
                                 + m_C_3d_cell(cell, 4, 1) * epsyy
                                 + m_C_3d_cell(cell, 4, 2) * epszz
                                 + m_C_3d_cell(cell, 4, 3) * epsyz
                                 + m_C_3d_cell(cell, 4, 4) * epszx
                                 + m_C_3d_cell(cell, 4, 5) * epsxy) ^ epszx)
                         + Uk * (( m_C_3d_cell(cell, 5, 0) * epsxx
                                 + m_C_3d_cell(cell, 5, 1) * epsyy
                                 + m_C_3d_cell(cell, 5, 2) * epszz
                                 + m_C_3d_cell(cell, 5, 3) * epsyz
                                 + m_C_3d_cell(cell, 5, 4) * epszx
                                 + m_C_3d_cell(cell, 5, 5) * epsxy) ^ epsxy)
                         ); // verified

          rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
          rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
          rhs_values[node_dof.dofId(cell.nodeId(0), 2)] += rhs(2);
          rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(3);
          rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(4);
          rhs_values[node_dof.dofId(cell.nodeId(1), 2)] += rhs(5);
          rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(6);
          rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(7);
          rhs_values[node_dof.dofId(cell.nodeId(2), 2)] += rhs(8);
          rhs_values[node_dof.dofId(cell.nodeId(3), 0)] += rhs(9);
          rhs_values[node_dof.dofId(cell.nodeId(3), 1)] += rhs(10);
          rhs_values[node_dof.dofId(cell.nodeId(3), 2)] += rhs(11);
          rhs_values[node_dof.dofId(cell.nodeId(4), 0)] += rhs(12);
          rhs_values[node_dof.dofId(cell.nodeId(4), 1)] += rhs(13);
          rhs_values[node_dof.dofId(cell.nodeId(4), 2)] += rhs(14);
          rhs_values[node_dof.dofId(cell.nodeId(5), 0)] += rhs(15);
          rhs_values[node_dof.dofId(cell.nodeId(5), 1)] += rhs(16);
          rhs_values[node_dof.dofId(cell.nodeId(5), 2)] += rhs(17);
          rhs_values[node_dof.dofId(cell.nodeId(6), 0)] += rhs(18);
          rhs_values[node_dof.dofId(cell.nodeId(6), 1)] += rhs(19);
          rhs_values[node_dof.dofId(cell.nodeId(6), 2)] += rhs(20);
          rhs_values[node_dof.dofId(cell.nodeId(7), 0)] += rhs(21);
          rhs_values[node_dof.dofId(cell.nodeId(7), 1)] += rhs(22);
          rhs_values[node_dof.dofId(cell.nodeId(7), 2)] += rhs(23);

        }
      }
    }

    // Add contributions to global RHS
    // for (Int8 a = 0; a < 8; ++a) {
    //   Node node = cell.node(a);
    //   if (node.isOwn()) {
    //     rhs_values[node_dof.dofId(node, 0)] += rhs_x_contributions[a];
    //     rhs_values[node_dof.dofId(node, 1)] += rhs_y_contributions[a];
    //     rhs_values[node_dof.dofId(node, 2)] += rhs_z_contributions[a];
    //   }
    // }
  }
}