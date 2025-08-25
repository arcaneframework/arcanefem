// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrixHexQuad.h                                      (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Elasticity     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD4, ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) = ∫∫ σ(𝐮):ε(𝐯)dΩ     with  𝐮 = (𝑢𝑥,𝑢𝑦) and 𝐯 = (𝑣𝑥,𝑣𝑦)
 *   σ(𝐮) is stress tensor       with  σᵢⱼ = λδᵢⱼεₖₖ + 2μεᵢⱼ
 *   ε(𝐯) is strain tensor       with  εᵢⱼ = 0.5 (∂𝑣ᵢ/∂xⱼ + ∂𝑣ⱼ/∂xᵢ)
 *
 *   the bilinear integral expands to
 *
 *      a(𝐮,𝐯) = ∫∫ [σ_𝑥𝑥ε_𝑥𝑥 + σ_𝑦𝑦ε_𝑦𝑦 + 2σ_𝑥𝑦ε_𝑥𝑦]dΩ
 *
 *   this further expands to
 *
 *      a(𝐮,𝐯) =   ∫∫ (λ+2μ)(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦)
 *               + ∫∫ λ(∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦)
 *               + ∫∫ μ(∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
 *
 *   - The first term is "normal strain energy"
 *   - The second term is "compressibility effect"
 *   - The third term is "shear energy"
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

ARCCORE_HOST_DEVICE RealMatrix<8, 8> computeElementMatrixQuad4Base(
const RealVector<4>& dxu, const RealVector<4>& dyu, Real integration_weight, Real lambda, Real mu)
{
  // Create displacement gradient vectors for x and y components
  // For quad4: 4 nodes × 2 DOF = 8 total DOF
  // Pattern: [u1x, u1y, u2x, u2y, u3x, u3y, u4x, u4y]
  RealVector<8> dxUx = { dxu(0), 0., dxu(1), 0., dxu(2), 0., dxu(3), 0. };
  RealVector<8> dyUx = { dyu(0), 0., dyu(1), 0., dyu(2), 0., dyu(3), 0. };
  RealVector<8> dxUy = { 0., dxu(0), 0., dxu(1), 0., dxu(2), 0., dxu(3) };
  RealVector<8> dyUy = { 0., dyu(0), 0., dyu(1), 0., dyu(2), 0., dyu(3) };

  // ∫∫ (λ+2μ)(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦)
  RealMatrix<8, 8> normal_strain_energy = (lambda + 2 * mu) * ((dxUx ^ dxUx) + (dyUy ^ dyUy)) * integration_weight;

  // ∫∫ λ(∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦)
  RealMatrix<8, 8> compressibility_effect = (lambda) * ((dyUy ^ dxUx) + (dxUx ^ dyUy)) * integration_weight;

  // ∫∫ μ(∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
  RealMatrix<8, 8> shear_energy = (mu) * ((dxUy + dyUx) ^ (dyUx + dxUy)) * integration_weight;

  return (normal_strain_energy + compressibility_effect + shear_energy);
}

RealMatrix<8, 8> FemModule::_computeElementMatrixQuad4(Cell cell)
{
  // Gauss points and weights for 2x2 quadrature
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // [-1/sqrt(3) , 1/sqrt(3)]
  constexpr Real w = 1.0;

  // Initialize the element matrix
  RealMatrix<8, 8> ae;
  ae.fill(0.0);

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      // Get the coordinates of the Gauss point in natural coordinates (ξ,η)
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian
      const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, xi, eta);
      const RealVector<4>& dxU = gp_info.dN_dx;
      const RealVector<4>& dyU = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      // Integration weight
      const Real integration_weight = detJ * w * w;

      // Add contribution from this Gauss point
      ae += computeElementMatrixQuad4Base(dxU, dyU, integration_weight, lambda, mu);
    }
  }

  return ae;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (HEXA8, ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) = ∫∫∫ [σ(𝐮):ε(𝐯)dΩ    with  𝐮 = (𝑢𝑥, 𝑢𝑦, 𝑢𝑧) and 𝐯 = (𝑣𝑥, 𝑣𝑦, 𝑣𝑧)
 *
 * where:
 *
 *   σ(𝐮) is stress tensor       with  σᵢⱼ = λδᵢⱼεₖₖ + 2μεᵢⱼ
 *   ε(𝐯) is strain tensor       with  εᵢⱼ = 0.5 (∂𝑣ᵢ/∂xⱼ + ∂𝑣ⱼ/∂xᵢ)
 *
 *   the bilinear integral expands to
 *
 *      a(𝐮,𝐯) = ∫∫∫ [σ_𝑥𝑥ε_𝑥𝑥 + σ_𝑦𝑦ε_𝑦𝑦 + σ_𝑧𝑧ε_𝑧𝑧 + 2σ_𝑥𝑦ε_𝑥𝑦 + 2σ_𝑦𝑧ε_𝑦𝑧 + 2σ_𝑧𝑥ε_𝑧𝑥]dΩ
 *
 *   this further expands to
 *
 *      a(𝐮,𝐯) =   ∫∫∫ (λ+2μ)(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧)
 *               + ∫∫∫ λ(∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑧/∂𝑧 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 )
 *               + ∫∫∫ μ(∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + μ(∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)(∂𝑣𝑦/∂𝑧 + ∂𝑣𝑧/∂𝑦) + μ(∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)(∂𝑣𝑧/∂𝑥 + ∂𝑣𝑥/∂𝑧)
 *
 *   - The first term is "normal strain energy"
 *   - The second term is "compressibility effect"
 *   - The third term is "shear energy"
 *
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<24, 24> computeElementMatrixHexa8Base(
const RealVector<8>& dxu, const RealVector<8>& dyu, const RealVector<8>& dzu,
Real integration_weight, Real lambda, Real mu)
{

  // Create displacement gradient vectors for x, y, z components
  // For hex8: 8 nodes × 3 DOF = 24 total DOF
  // Pattern: [u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, ..., u8x, u8y, u8z]

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

  // Normal strain energy terms: ∫∫∫ (λ+2μ)(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧)
  RealMatrix<24, 24> normal_strain_energy = (lambda + 2 * mu) *
  ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz)) * integration_weight;

  // Compressibility effects: ∫∫∫ λ(∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑧/∂𝑧 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 )
  RealMatrix<24, 24> compressibility_effect = lambda *
  ((dxUx ^ (dyUy + dzUz)) + (dyUy ^ (dxUx + dzUz)) + (dzUz ^ (dxUx + dyUy))) * integration_weight;

  // Shear energy terms: ∫∫∫ μ(∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 
  //                         μ(∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)(∂𝑣𝑦/∂𝑧 + ∂𝑣𝑧/∂𝑦) + 
  //                         μ(∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)(∂𝑣𝑧/∂𝑥 + ∂𝑣𝑥/∂𝑧)
  RealMatrix<24, 24> shear_energy = mu *
  (((dyUx + dxUy) ^ (dyUx + dxUy)) +
   ((dzUx + dxUz) ^ (dzUx + dxUz)) +
   ((dzUy + dyUz) ^ (dzUy + dyUz))) * integration_weight;

  return (normal_strain_energy + compressibility_effect + shear_energy);
}

RealMatrix<24, 24> FemModule::_computeElementMatrixHexa8(Cell cell)
{
  // Gauss points and weights for 2x2x2 quadrature
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // [-1/sqrt(3), 1/sqrt(3)]
  constexpr Real w = 1.0;

  // Initialize the element matrix
  RealMatrix<24, 24> ae;
  ae.fill(0.0);

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      for (Int8 izeta = 0; izeta < 2; ++izeta) {
        // Get the coordinates of the Gauss point in natural coordinates (ξ,η,ζ)
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        // Get shape function gradients w.r.t (x,y,z) and determinant of Jacobian
        const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa8(
        cell, m_node_coord, xi, eta, zeta);
        const RealVector<8>& dxU = gp_info.dN_dx;
        const RealVector<8>& dyU = gp_info.dN_dy;
        const RealVector<8>& dzU = gp_info.dN_dz;
        const Real detJ = gp_info.det_j;

        // Integration weight
        const Real integration_weight = detJ * w * w * w;

        // Add contribution from this Gauss point
        ae += computeElementMatrixHexa8Base(dxU, dyU, dzU, integration_weight, lambda, mu);
      }
    }
  }

  return ae;
}