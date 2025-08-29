// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrixHexQuad.h                                      (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Elastodynamics */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD4, ℙ1 FE).
 *
 * Theory:
 *   a(𝐮,𝐯) = ∫∫ [(∂²𝐮/∂𝑡²).(𝐯)]dΩ + ∫∫ [σ(𝐮):ε(𝐯)]dΩ
 *
 *   with  trial func 𝐮 = (𝑢𝑥,𝑢𝑦) and test func 𝐯 = (𝑣𝑥,𝑣𝑦),
 *   σ(𝐮) is stress tensor with     σᵢⱼ = λδᵢⱼεₖₖ + 2μεᵢⱼ
 *   ε(𝐯) is strain tensor with     εᵢⱼ = 0.5 (∂𝑣ᵢ/∂xⱼ + ∂𝑣ⱼ/∂xᵢ)
 *
 *   the bilinear integral after Newmark-Beta and damping terms expands to:
 *
 *      a(𝐮,𝐯) =   ∫∫ (c₀)(𝐮.𝐯)
 *               + ∫∫ (c₁)(∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦)
 *               + ∫∫ (c₁+2c₂)(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦)
 *               + ∫∫ (c₂)(∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
 *
 *   with c₀ = (ρ)/(β δ𝑡²) + (ηₘ ρ γ)/(β δ𝑡)
 *        c₁ = λ + (λ ηₖ γ)/(β δ𝑡)
 *        c₂ = 2μ + (2μ ηₖ γ)/(β δ𝑡)
 *
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

RealMatrix<8, 8> FemModule::_computeElementMatrixQuad4(Cell cell)
{
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 };
  constexpr Real w = 1.0;

  RealMatrix<8, 8> ae; ae.fill(0.0);

  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, xi, eta);
      const RealVector<4>& dN_dx = gp_info.dN_dx;
      const RealVector<4>& dN_dy = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      const Real weight = detJ * w * w;

      RealVector<8> dxUx = {dN_dx(0),0, dN_dx(1),0, dN_dx(2),0, dN_dx(3),0};
      RealVector<8> dyUx = {dN_dy(0),0, dN_dy(1),0, dN_dy(2),0, dN_dy(3),0};
      RealVector<8> dxUy = {0,dN_dx(0), 0,dN_dx(1), 0,dN_dx(2), 0,dN_dx(3)};
      RealVector<8> dyUy = {0,dN_dy(0), 0,dN_dy(1), 0,dN_dy(2), 0,dN_dy(3)};

      RealVector<4> N = ArcaneFemFunctions::FeOperation2D::computeShapeFunctionsQuad4(xi, eta);
      RealVector<8> Nx = {N(0),0, N(1),0, N(2),0, N(3),0};
      RealVector<8> Ny = {0,N(0), 0,N(1), 0,N(2), 0,N(3)};

      ae += c0 * ((Nx ^ Nx) + (Ny ^ Ny)) * weight
            + c1 * ((dyUy ^ dxUx) + (dxUx ^ dyUy)) * weight
            + (2*c2 + c1) * ((dxUx ^ dxUx) + (dyUy ^ dyUy)) * weight
            + c2 *  ((dxUy + dyUx) ^ (dyUx + dxUy))  * weight;
    }
  }

  return ae;
}