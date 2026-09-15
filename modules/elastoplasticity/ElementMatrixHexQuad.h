// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrixHexQuad.h                                      (C) 2000-2026 */
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
 *      a(𝐮,𝐯) =   ∫∫ C_tang11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + C_tang12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + C_tang13 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ C_tang12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + C_tang22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + C_tang23 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ C_tang13 ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C_tang23 ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C_tang33 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<8, 8> computeElementMatrixQuad4Base(
const RealVector<4>& dxu, const RealVector<4>& dyu, Real integration_weight, RealMatrix<3, 3> C_tang)
{
  RealVector<8> epsxx = { dxu(0), 0., dxu(1), 0., dxu(2), 0., dxu(3), 0. };
  RealVector<8> epsyy = { 0., dyu(0), 0., dyu(1), 0., dyu(2), 0., dyu(3) };
  RealVector<8> epsxy = { dyu(0), dxu(0), dyu(1), dxu(1), dyu(2), dxu(2), dyu(3), dxu(3) };

  // Kelvin notation: scale shear strains by 1/sqrt(2)
  epsxy = M_SQRT1_2 * epsxy;

  // ∫∫ C_tang11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + C_tang12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + C_tang13 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<8, 8> sigmaXepsxx = (C_tang(0, 0) * epsxx + C_tang(0, 1) * epsyy + C_tang(0, 2) * epsxy) ^ epsxx;

  // ∫∫ C_tang12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + C_tang22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + C_tang23 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<8, 8> sigmaXepsyy = (C_tang(0, 1) * epsxx + C_tang(1, 1) * epsyy + C_tang(1, 2) * epsxy) ^ epsyy;

  // ∫∫ C_tang13 ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C_tang23 ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C_tang33 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
  RealMatrix<8, 8> sigmaXepsxy = (C_tang(0, 2) * epsxx + C_tang(1, 2) * epsyy + C_tang(2, 2) * epsxy) ^ epsxy;

  return integration_weight * (sigmaXepsxx + sigmaXepsyy + sigmaXepsxy);
}

RealMatrix<8, 8> FemModuleElastoplasticity::_computeElementMatrixQuad4(Cell cell)
{
  // Gauss points and weights for 2x2 quadrature
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // [-1/sqrt(3) , 1/sqrt(3)]
  constexpr Real w = 1.0;

  // Initialize the element matrix
  RealMatrix<8, 8> ae;
  ae.fill(0.0);

  // Loop over Gauss points
  Int8 iGP = 0; // TODO verify order with ixi, ieta, izeta.
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
      if (m_gp_material_tensor_strategy == "local") {
        ae += computeElementMatrixQuad4Base(dxU, dyU, integration_weight, m_C_tang_2d);
      } else {
        RealMatrix<3, 3> C_tang_2d;
        for (Int32 ix = 0; ix < 3; ++ix) {
          for (Int32 iy = 0; iy < 3; ++iy) {
            C_tang_2d(ix, iy) = m_C_tang_gp(cell, iGP, ix, iy);
          }
        }
        iGP++;
        ae += computeElementMatrixQuad4Base(dxU, dyU, integration_weight, C_tang_2d);
      }
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
 *      a(𝐮,𝐯) =   ∫∫∫ C_tang11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + C_tang12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + C_tang13 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + C_tang14 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑥/∂𝑥 + C_tang15 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑥/∂𝑥 + C_tang16 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫∫ C_tang12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + C_tang22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + C_tang23 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 + C_tang24 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑦/∂𝑦 + C_tang25 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑦/∂𝑦 + C_tang26 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑦/∂𝑦
 *               + ∫∫∫ C_tang13 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 + C_tang23 ∂𝑢𝑧/∂𝑦 ∂𝑣𝑧/∂𝑧 + C_tang33 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧 + C_tang34 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑧/∂𝑧 + C_tang35 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑧/∂𝑧 + C_tang36 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑧/∂𝑧
 *               + ∫∫∫ C_tang14 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang24 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang34 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang44 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang45 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang46 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)
 *               + ∫∫∫ C_tang15 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang25 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang35 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang45 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang55 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang56 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)
 *               + ∫∫∫ C_tang16 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang26 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang36 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang46 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang56 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang66 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)
 *
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<24, 24> computeElementMatrixHexa8Base(
const RealVector<8>& dxu, const RealVector<8>& dyu, const RealVector<8>& dzu,
Real integration_weight, RealMatrix<6, 6> C_tang)
{
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

  // Kelvin notation: scale shear strains by 1/sqrt(2)
  epsyz = M_SQRT1_2 * epsyz;
  epszx = M_SQRT1_2 * epszx;
  epsxy = M_SQRT1_2 * epsxy;

  // ∫∫∫ C_tang11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + C_tang12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + C_tang13 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + C_tang14 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑥/∂𝑥 + C_tang15 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑥/∂𝑥 + C_tang16 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<24, 24> sigmaXepsxx = (C_tang(0, 0) * epsxx + C_tang(0, 1) * epsyy + C_tang(0, 2) * epszz + C_tang(0, 3) * epsyz + C_tang(0, 4) * epszx + C_tang(0, 5) * epsxy) ^ epsxx;
  // ∫∫∫ C_tang12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + C_tang22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + C_tang23 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 + C_tang24 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑦/∂𝑦 + C_tang25 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑦/∂𝑦 + C_tang26 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑦/∂𝑦
  RealMatrix<24, 24> sigmaXepsyy = (C_tang(0, 1) * epsxx + C_tang(1, 1) * epsyy + C_tang(1, 2) * epszz + C_tang(1, 3) * epsyz + C_tang(1, 4) * epszx + C_tang(1, 5) * epsxy) ^ epsyy;
  // ∫∫∫ C_tang13 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 + C_tang23 ∂𝑢𝑧/∂𝑦 ∂𝑣𝑧/∂𝑧 + C_tang33 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧 + C_tang34 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑧/∂𝑧 + C_tang35 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑧/∂𝑧 + C_tang36 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑧/∂𝑧
  RealMatrix<24, 24> sigmaXepszz = (C_tang(0, 2) * epsxx + C_tang(1, 2) * epsyy + C_tang(2, 2) * epszz + C_tang(2, 3) * epsyz + C_tang(2, 4) * epszx + C_tang(2, 5) * epsxy) ^ epszz;
  // ∫∫∫ C_tang14 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang24 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang34 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang44 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang45 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang46 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)
  RealMatrix<24, 24> sigmaXepsyz = (C_tang(0, 3) * epsxx + C_tang(1, 3) * epsyy + C_tang(2, 3) * epszz + C_tang(3, 3) * epsyz + C_tang(3, 4) * epszx + C_tang(3, 5) * epsxy) ^ epsyz;
  // ∫∫∫ C_tang15 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang25 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang35 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang45 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang55 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang56 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)
  RealMatrix<24, 24> sigmaXepszx = (C_tang(0, 4) * epsxx + C_tang(1, 4) * epsyy + C_tang(2, 4) * epszz + C_tang(3, 4) * epsyz + C_tang(4, 4) * epszx + C_tang(4, 5) * epsxy) ^ epszx;
  // ∫∫∫ C_tang16 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang26 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang36 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang46 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang56 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang66 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)
  RealMatrix<24, 24> sigmaXepsxy = (C_tang(0, 5) * epsxx + C_tang(1, 5) * epsyy + C_tang(2, 5) * epszz + C_tang(3, 5) * epsyz + C_tang(4, 5) * epszx + C_tang(5, 5) * epsxy) ^ epsxy;

  return integration_weight * ( sigmaXepsxx + sigmaXepsyy + sigmaXepszz + sigmaXepsyz + sigmaXepszx + sigmaXepsxy);
}

RealMatrix<24, 24> FemModuleElastoplasticity::_computeElementMatrixHexa8(Cell cell)
{
  // Gauss points and weights for 2x2x2 quadrature
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // [-1/sqrt(3), 1/sqrt(3)]
  constexpr Real w = 1.0;

  // Initialize the element matrix
  RealMatrix<24, 24> ae;
  ae.fill(0.0);

  // Loop over Gauss points
  Int8 iGP = 0; // TODO verify order with ixi, ieta, izeta.
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
        if (m_gp_material_tensor_strategy == "local") {
          ae += computeElementMatrixHexa8Base(dxU, dyU, dzU, integration_weight, m_C_tang_3d);
        } else {
          RealMatrix<6, 6> C_tang_3d;
          for (Int32 ix = 0; ix < 6; ++ix) {
            for (Int32 iy = 0; iy < 6; ++iy) {
              C_tang_3d(ix, iy) = m_C_tang_gp(cell, iGP, ix, iy);
            }
          }
          iGP++;
          ae += computeElementMatrixHexa8Base(dxU, dyU, dzU, integration_weight,  C_tang_3d);
        }
      }
    }
  }

  return ae;
}