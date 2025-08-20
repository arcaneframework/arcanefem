// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrixHexQuad.h                                      (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Acoustics      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD4, ℙ1 FE).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫ (∂𝑢/∂𝑥 ∂𝑣/∂𝑥  + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦)dΩ  + ∫∫ 𝑘²𝑢𝑣 dΩ
 *
 * Steps involved:
 * 1. Define Gauss points (2x2) and weights.
 * 2. Loop over Gauss points to compute the gradients in physical space
 *    and the determinant of the Jacobian, via computeGradientsAndJacobianQuad4.
 * 3. Compute the integration weight.
 * 4. Assemble the element matrix using the computed gradients.
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<4, 4> FemModule::_computeElementMatrixQuad4(Cell cell)
{
  // Gauss points and weights for 2x2 quadrature
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // [-1/sqrt(3) , 1/sqrt(3)]
  constexpr Real w = 1.0;

  // Initialize the element matrix
  RealMatrix<4, 4> ae;
  ae.fill(0.0);

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {

      // Get the coordinates of the Gauss point in natural coordinates (ξ,η,ζ)
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian
      const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, xi, eta);
      const RealVector<4>& dxU = gp_info.dN_dx;
      const RealVector<4>& dyU = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      // Shape functions at the Gauss point
      RealVector<4> N = ArcaneFemFunctions::FeOperation2D::computeShapeFunctionsQuad4(xi, eta);

      // Integration weight
      const Real integration_weight = detJ * w * w;

      // stiffness matrix assembly
      ae += (dxU ^ dxU) * integration_weight + (dyU ^ dyU) * integration_weight +(N ^ N) * m_kc2 * integration_weight;
    }
  }
  return ae;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a hexahedral element (HEXA8, ℙ1 FE).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫∫ (∂𝑢/∂𝑥 ∂𝑣/∂𝑥 + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦 + ∂𝑢/∂𝑧 ∂𝑣/∂𝑧)dΩ  + ∫∫∫ 𝑘²𝑢𝑣 dΩ
 *
 * Steps involved:
 * 1. Define Gauss points (2x2x2) and weights.
 * 2. Loop over Gauss points to compute the gradients in physical space (𝑥,𝑦,𝑧)
 *    and the determinant of the Jacobian via computeGradientsAndJacobianHexa8.
 * 3. Compute the integration weight.
 * 4. Assemble the element matrix using the computed gradients.
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<8, 8> FemModule::_computeElementMatrixHexa8(Cell cell)
{
  // 2x2x2 Gauss points and weights for [-1,1]^3
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
  constexpr Real w = 1.0;

  // Initialize the element matrix
  RealMatrix<8, 8> ae;
  ae.fill(0.0);

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      for (Int8 izeta = 0; izeta < 2; ++izeta) {

        // Get the coordinates of Gauss points in natural coordinates (ξ,η,ζ)
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        // Get shape function gradients w.r.t (𝑥,𝑦,𝑧) and determinant of Jacobian
        const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa8(cell, m_node_coord, xi, eta, zeta);
        const RealVector<8>& dxU = gp_info.dN_dx;
        const RealVector<8>& dyU = gp_info.dN_dy;
        const RealVector<8>& dzU = gp_info.dN_dz;
        const Real detJ = gp_info.det_j;

        // Shape functions at the Gauss point
        RealVector<8> N = ArcaneFemFunctions::FeOperation3D::computeShapeFunctionsHexa8(xi, eta, zeta);

        // Integration weight
        const Real integration_weight = detJ * w * w;

        // Assemble element matrix (variational form)
        ae += (dxU ^ dxU) * integration_weight + (dyU ^ dyU) * integration_weight + (dzU ^ dzU) * integration_weight +
        (N ^ N) * m_kc2 * integration_weight;
      }
    }
  }
  return ae;
}