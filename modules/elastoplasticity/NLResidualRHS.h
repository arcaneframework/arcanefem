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

    RealVector<6> epsxx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
    RealVector<6> epsyy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };
    RealVector<6> epsxy = { dyu[0], dxu[0], dyu[1], dxu[1], dyu[2], dxu[2] };

    Real3x3 grad_U = ArcaneFemFunctions::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_U);

    Real epsxx_U = grad_U(0, 0);
    Real epsyy_U = grad_U(1, 1);
    Real epsxy_U = grad_U(0, 1) + grad_U(1, 0);

    // Real sigmaxx_U = m_C_2d_cell[cell](0, 0) * epsxx_U + m_C_2d_cell[cell](0, 1) * epsyy_U + m_C_2d_cell[cell](0, 2) * epsxy_U;
    // Real sigmayy_U = m_C_2d_cell[cell](1, 0) * espxx_U + m_C_2d_cell[cell](1, 1) * epsyy_U + m_C_2d_cell[cell](1, 2) * epsxy_U;
    // Real sigmaxy_U = m_C_2d_cell[cell](2, 0) * espxx_U + m_C_2d_cell[cell](2, 1) * epsyy_U + m_C_2d_cell[cell](2, 2) * epsxy_U;

    // test
    Real sigmaxx_U = m_C_tang_2d(0, 0) * epsxx_U + m_C_tang_2d(0, 1) * epsyy_U + m_C_tang_2d(0, 2) * epsxy_U;
    Real sigmayy_U = m_C_tang_2d(1, 0) * epsxx_U + m_C_tang_2d(1, 1) * epsyy_U + m_C_tang_2d(1, 2) * epsxy_U;
    Real sigmaxy_U = m_C_tang_2d(2, 0) * epsxx_U + m_C_tang_2d(2, 1) * epsyy_U + m_C_tang_2d(2, 2) * epsxy_U;

    //----------------------------------------------------------------------
    //  ∫∫∫ (c₅)(∇𝐮ₙ.∇𝐯) + ∫∫∫ (c₆)(ε(𝐮ₙ):ε(𝐯))
    //----------------------------------------------------------------------
    // RealVector<6> rhs = Uk * ((m_C_2d_cell[cell](0, 0) * epsxx
    //                         + m_C_2d_cell[cell](0, 1) * epsyy
    //                         + m_C_2d_cell[cell](0, 2) * epsxy) ^ epsxx)
    //                   + Uk * ((m_C_2d_cell[cell](1, 0) * epsxx
    //                         + m_C_2d_cell[cell](1, 1) * epsyy
    //                         + m_C_2d_cell[cell](1, 2) * epsxy) ^ epsyy)
    //                   + Uk * ((m_C_2d_cell[cell](2, 0) * epsxx
    //                         + m_C_2d_cell[cell](2, 1) * epsyy
    //                         + m_C_2d_cell[cell](2, 2) * epsxy) ^ epsxy);  // To verify mathematically

    RealVector<6> rhs = - area * (sigmaxx_U * epsxx + sigmayy_U * epsyy + sigmaxy_U * epsxy);

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
    Real rhs_x_contributions[4] = { 0., 0., 0., 0. };
    Real rhs_y_contributions[4] = { 0., 0., 0., 0. };

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
        const Real detJ = gp_info.det_j;

        // compute integration weight
        Real integration_weight = weight * detJ;

        // Interpolate fields (𝐮ₙ,𝐮ᵗₙ,𝐮ᵗᵗ) at the quadrature point: (.)_gp = ∑ 𝑁ᵢ * (.)

        // Source force term 𝐟

        // ∫∫ (𝐟.𝐯) + ∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯)
      }
    }

    // Add contributions to global RHS
    for (Int8 a = 0; a < 4; ++a) {
      Node node = cell.node(a);
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += rhs_x_contributions[a];
        rhs_values[node_dof.dofId(node, 1)] += rhs_y_contributions[a];
      }
    }
  }
}

inline void FemModuleElastoplasticity::
_applyResidualRHSTetra4(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  RealVector<12> Uy = { 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0. };
  RealVector<12> Ux = { 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0. };
  RealVector<12> Uz = { 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1. };
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);
    Real4 dxu = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
    Real4 dyu = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
    Real4 dzu = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);
    RealVector<12> dxUx = { dxu[0], 0., 0., dxu[1], 0., 0., dxu[2], 0., 0., dxu[3], 0., 0. };
    RealVector<12> dyUx = { dyu[0], 0., 0., dyu[1], 0., 0., dyu[2], 0., 0., dyu[3], 0., 0. };
    RealVector<12> dzUx = { dzu[0], 0., 0., dzu[1], 0., 0., dzu[2], 0., 0., dzu[3], 0., 0. };
    RealVector<12> dxUy = { 0., dxu[0], 0., 0., dxu[1], 0., 0., dxu[2], 0., 0., dxu[3], 0. };
    RealVector<12> dyUy = { 0., dyu[0], 0., 0., dyu[1], 0., 0., dyu[2], 0., 0., dyu[3], 0. };
    RealVector<12> dzUy = { 0., dzu[0], 0., 0., dzu[1], 0., 0., dzu[2], 0., 0., dzu[3], 0. };
    RealVector<12> dxUz = { 0., 0., dxu[0], 0., 0., dxu[1], 0., 0., dxu[2], 0., 0., dxu[3] };
    RealVector<12> dyUz = { 0., 0., dyu[0], 0., 0., dyu[1], 0., 0., dyu[2], 0., 0., dyu[3] };
    RealVector<12> dzUz = { 0., 0., dzu[0], 0., 0., dzu[1], 0., 0., dzu[2], 0., 0., dzu[3] };
    //----------------------------------------------------------------------
    //  ∫∫∫ (c₅)(∇𝐮ₙ.∇𝐯) + ∫∫∫ (c₆)(ε(𝐮ₙ):ε(𝐯)) +
    //----------------------------------------------------------------------
    RealVector<12> rhs = {0., 0., 0., 0.,
                          0., 0., 0., 0.,
                          0., 0., 0., 0.};
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
    Real rhs_x_contributions[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };
    Real rhs_y_contributions[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };
    Real rhs_z_contributions[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };

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
          const Real detJ = gp_info.det_j;

          // compute integration weight
          Real integration_weight = weight * detJ;

          // Interpolate fields (𝐮ₙ,𝐮ᵗₙ,𝐮ᵗᵗ) at the quadrature point: (.)_gp = ∑ 𝑁ᵢ * (.)

          // Source force term 𝐟
          Real3 f_gp(f[0], f[1], f[2]);

          // ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯)

        }
      }
    }
    // Add contributions to global RHS
    for (Int8 a = 0; a < 8; ++a) {
      Node node = cell.node(a);
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += rhs_x_contributions[a];
        rhs_values[node_dof.dofId(node, 1)] += rhs_y_contributions[a];
        rhs_values[node_dof.dofId(node, 2)] += rhs_z_contributions[a];
      }
    }
  }
}