// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrixHexQuad.h                                      (C) 2022-2026 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Nonlinear      */
/* Fourier                                                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD4, ℙ1 FE).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫ λ(𝑢)(∂𝑢/∂𝑥 ∂𝑣/∂𝑥  + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦)dΩ
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

RealMatrix<4, 4> FemModuleFourierNL::_computeElementMatrixQuad4(Cell cell)
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

      const Real lambda_0 = m_node_lambda[cell.nodeId(0)];
      const Real lambda_1 = m_node_lambda[cell.nodeId(1)];
      const Real lambda_2 = m_node_lambda[cell.nodeId(2)];
      const Real lambda_3 = m_node_lambda[cell.nodeId(3)];

      const Real lambda_dxU_0 = lambda_0 * dxU[0];
      const Real lambda_dxU_1 = lambda_1 * dxU[1];
      const Real lambda_dxU_2 = lambda_2 * dxU[2];
      const Real lambda_dxU_3 = lambda_3 * dxU[3];

      const Real lambda_dyU_0 = lambda_0 * dyU[0];
      const Real lambda_dyU_1 = lambda_1 * dyU[1];
      const Real lambda_dyU_2 = lambda_2 * dyU[2];
      const Real lambda_dyU_3 = lambda_3 * dyU[3];

      const RealVector<4>& lambda_dxU {lambda_dxU_0, lambda_dxU_1, lambda_dxU_2, lambda_dxU_3};
      const RealVector<4>& lambda_dyU {lambda_dyU_0, lambda_dyU_1, lambda_dyU_2, lambda_dyU_3};

      // Integration weight
      const Real lambda_cell = (lambda_0 + lambda_1 + lambda_2 + lambda_3) / 4.0;
      const Real integration_weight = detJ * w * w;

      // stiffness matrix assembly
      // ae += (dxU ^ dxU) * integration_weight * lambda + (dyU ^ dyU) * integration_weight * lambda;
      ae += (dxU ^ dxU) * integration_weight * lambda_cell + (dyU ^ dyU) * integration_weight * lambda_cell;
      // ae += (lambda_dxU ^ dxU) * integration_weight + (lambda_dyU ^ dyU) * integration_weight;
    }
  }
  return ae;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a hexahedral element (HEXA8, ℙ1 FE).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫∫ λ(∂𝑢/∂𝑥 ∂𝑣/∂𝑥 + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦 + ∂𝑢/∂𝑧 ∂𝑣/∂𝑧)dΩ
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

RealMatrix<8, 8> FemModuleFourierNL::_computeElementMatrixHexa8(Cell cell)
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

        // Integration weight
        const Real integration_weight = detJ * w * w;

        // Assemble element matrix (variational form)
        ae += (dxU ^ dxU) * integration_weight * lambda + (dyU ^ dyU) * integration_weight * lambda + (dzU ^ dzU) * integration_weight * lambda;
      }
    }
  }
  return ae;
}