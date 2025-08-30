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

  RealMatrix<8, 8> ae;
  ae.fill(0.0);

  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, xi, eta);
      const RealVector<4>& dN_dx = gp_info.dN_dx;
      const RealVector<4>& dN_dy = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      const Real weight = detJ * w * w;

      RealVector<8> dxUx = { dN_dx(0), 0, dN_dx(1), 0, dN_dx(2), 0, dN_dx(3), 0 };
      RealVector<8> dyUx = { dN_dy(0), 0, dN_dy(1), 0, dN_dy(2), 0, dN_dy(3), 0 };
      RealVector<8> dxUy = { 0, dN_dx(0), 0, dN_dx(1), 0, dN_dx(2), 0, dN_dx(3) };
      RealVector<8> dyUy = { 0, dN_dy(0), 0, dN_dy(1), 0, dN_dy(2), 0, dN_dy(3) };

      RealVector<4> N = ArcaneFemFunctions::FeOperation2D::computeShapeFunctionsQuad4(xi, eta);
      RealVector<8> Nx = { N(0), 0, N(1), 0, N(2), 0, N(3), 0 };
      RealVector<8> Ny = { 0, N(0), 0, N(1), 0, N(2), 0, N(3) };

      ae += c0 * ((Nx ^ Nx) + (Ny ^ Ny)) * weight + c1 * ((dyUy ^ dxUx) + (dxUx ^ dyUy)) * weight + (2 * c2 + c1) * ((dxUx ^ dxUx) + (dyUy ^ dyUy)) * weight + c2 * ((dxUy + dyUx) ^ (dyUx + dxUy)) * weight;
    }
  }

  return ae;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a hexahedral element (HEXA8, ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) = ∫∫∫ [(∂²𝐮/∂𝑡²).(𝐯)] dΩ + ∫∫∫ [σ(𝐮):ε(𝐯)] dΩ
 *
 *   with trial function 𝐮 = (𝑢𝑥, 𝑢𝑦, 𝑢𝑧) and test function 𝐯 = (𝑣𝑥, 𝑣𝑦, 𝑣𝑧),
 *   σ(𝐮) is the stress tensor, given by     σᵢⱼ = λδᵢⱼεₖₖ + 2μεᵢⱼ
 *   ε(𝐯) is the strain tensor, defined as    εᵢⱼ = 0.5 (∂𝑣ᵢ/∂xⱼ + ∂𝑣ⱼ/∂xᵢ)
 *
 *   The bilinear integral after applying the Newmark-Beta scheme and damping terms expands to:
 *
 *      a(𝐮,𝐯) =   ∫∫∫ (c₀)(𝐮 ⋅ 𝐯) dΩ
 *               + ∫∫∫ (c₁) (∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧 +
 *                           ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 +
 *                           ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 +
 *                           ∂𝑢𝑦/∂𝑦 ∂𝑣𝑧/∂𝑧 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 )
 *               + ∫∫∫ (c₂)(2(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧) +
 *                           (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) +
 *                           (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)(∂𝑣𝑦/∂𝑧 + ∂𝑣𝑧/∂𝑦) +
 *                           (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)(∂𝑣𝑧/∂𝑥 + ∂𝑣𝑥/∂𝑧) )
 *
 *   with c₀ = (ρ)/(β δ𝑡²) + (ηₘ ρ γ)/(β δ𝑡)
 *        c₁ = λ + (λ ηₖ γ)/(β δ𝑡)
 *        c₂ = 2μ + (2μ ηₖ γ)/(β δ𝑡)
 *
 * Steps involved:
 * 1. Define Gauss points (2x2) and weights.
 * 2. Loop over Gauss points to compute the gradients in physical space
 *    and the determinant of the Jacobian, via computeGradientsAndJacobianHexa8
 * 3. Compute the integration weight.
 * 4. Assemble the element matrix using the computed gradients.
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<24, 24> FemModule::_computeElementMatrixHexa8(Cell cell)
{
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 };
  constexpr Real w = 1.0;

  RealMatrix<24, 24> ae;
  ae.fill(0.0);

  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      for (Int8 izeta = 0; izeta < 2; ++izeta) {

        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa8(cell, m_node_coord, xi, eta, zeta);
        const RealVector<8>& dxu = gp_info.dN_dx;
        const RealVector<8>& dyu = gp_info.dN_dy;
        const RealVector<8>& dzu = gp_info.dN_dz;
        const Real detJ = gp_info.det_j;

        const Real weight = detJ * w * w * w;

        // X-displacement gradients
        RealVector<24> dxUx = { dxu(0), 0., 0., dxu(1), 0., 0., dxu(2), 0., 0., dxu(3), 0., 0.,
                                dxu(4), 0., 0., dxu(5), 0., 0., dxu(6), 0., 0., dxu(7), 0., 0. };
        RealVector<24> dyUx = { dyu(0), 0., 0., dyu(1), 0., 0., dyu(2), 0., 0., dyu(3), 0., 0.,
                                dyu(4), 0., 0., dyu(5), 0., 0., dyu(6), 0., 0., dyu(7), 0., 0. };
        RealVector<24> dzUx = { dzu(0), 0., 0., dzu(1), 0., 0., dzu(2), 0., 0., dzu(3), 0., 0.,
                                dzu(4), 0., 0., dzu(5), 0., 0., dzu(6), 0., 0., dzu(7), 0., 0. };

        // Y-displacement gradients
        RealVector<24> dxUy = { 0., dxu(0), 0., 0., dxu(1), 0., 0., dxu(2), 0., 0., dxu(3), 0.,
                                0., dxu(4), 0., 0., dxu(5), 0., 0., dxu(6), 0., 0., dxu(7), 0. };
        RealVector<24> dyUy = { 0., dyu(0), 0., 0., dyu(1), 0., 0., dyu(2), 0., 0., dyu(3), 0.,
                                0., dyu(4), 0., 0., dyu(5), 0., 0., dyu(6), 0., 0., dyu(7), 0. };
        RealVector<24> dzUy = { 0., dzu(0), 0., 0., dzu(1), 0., 0., dzu(2), 0., 0., dzu(3), 0.,
                                0., dzu(4), 0., 0., dzu(5), 0., 0., dzu(6), 0., 0., dzu(7), 0. };

        // Z-displacement gradients
        RealVector<24> dxUz = { 0., 0., dxu(0), 0., 0., dxu(1), 0., 0., dxu(2), 0., 0., dxu(3),
                                0., 0., dxu(4), 0., 0., dxu(5), 0., 0., dxu(6), 0., 0., dxu(7) };
        RealVector<24> dyUz = { 0., 0., dyu(0), 0., 0., dyu(1), 0., 0., dyu(2), 0., 0., dyu(3),
                                0., 0., dyu(4), 0., 0., dyu(5), 0., 0., dyu(6), 0., 0., dyu(7) };
        RealVector<24> dzUz = { 0., 0., dzu(0), 0., 0., dzu(1), 0., 0., dzu(2), 0., 0., dzu(3),
                                0., 0., dzu(4), 0., 0., dzu(5), 0., 0., dzu(6), 0., 0., dzu(7) };

        RealVector<8> N = ArcaneFemFunctions::FeOperation3D::computeShapeFunctionsHexa8(xi, eta, zeta);
        RealVector<24> Nx = { N(0), 0, 0, N(1), 0, 0, N(2), 0, 0, N(3), 0, 0,
                              N(4), 0, 0, N(5), 0, 0, N(6), 0, 0, N(7), 0, 0 };
        RealVector<24> Ny = { 0, N(0), 0, 0, N(1), 0, 0, N(2), 0, 0, N(3), 0,
                              0, N(4), 0, 0, N(5), 0, 0, N(6), 0, 0, N(7), 0 };
        RealVector<24> Nz = { 0, 0, N(0), 0, 0, N(1), 0, 0, N(2), 0, 0, N(3),
                              0, 0, N(4), 0, 0, N(5), 0, 0, N(6), 0, 0, N(7) };

        ae += c0 * ((Nx ^ Nx) + (Ny ^ Ny) + (Nz ^ Nz)) * weight + (c1) * ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) + (dyUy ^ dxUx) + (dxUx ^ dyUy) + (dzUz ^ dxUx) + (dxUx ^ dzUz) + (dyUy ^ dzUz) + (dzUz ^ dyUy)) * weight + (c2) * (2. * ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz)) + (((dxUy + dyUx) ^ (dyUx + dxUy)) + ((dzUy + dyUz) ^ (dyUz + dzUy)) + ((dxUz + dzUx) ^ (dzUx + dxUz)))) * weight;
      }
    }
  }

  return ae;
}