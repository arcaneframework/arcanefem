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
#ifndef ARCANFEM_ELASTOPLATICITY2_ELEMENTMATRIXHEXQUAD_H
#define ARCANFEM_ELASTOPLATICITY2_ELEMENTMATRIXHEXQUAD_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "femutils/FemUtils.h"
#include "femutils/ArcaneFemFunctions.h"
#include "femutils/ArcaneFemFunctionsGpu.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::ArcaneFem
{

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
 *      a(𝐮,𝐯) =   ∫∫ 𝑪ᵗ₀₀ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ 𝑪ᵗ₀₁ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ 𝑪ᵗ₀₂ ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₁₂ ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₂₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

inline ARCCORE_HOST_DEVICE RealMatrix<8, 8>
computeElementMatrixQuad4Base(const RealVector<4>& dxu, const RealVector<4>& dyu,
                              Real integration_weight, RealMatrix<3, 3> C_tang)
{
  RealVector<8> epsxx = { dxu(0), 0., dxu(1), 0., dxu(2), 0., dxu(3), 0. };
  RealVector<8> epsyy = { 0., dyu(0), 0., dyu(1), 0., dyu(2), 0., dyu(3) };
  RealVector<8> epsxy = { dyu(0), dxu(0), dyu(1), dxu(1), dyu(2), dxu(2), dyu(3), dxu(3) };

  // Kelvin notation: scale shear strains by 1/sqrt(2)
  epsxy = M_SQRT1_2 * epsxy;

  // ∫∫ 𝑪ᵗ₀₀ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<8, 8> sigmaXepsxx = (C_tang(0, 0) * epsxx + C_tang(0, 1) * epsyy + C_tang(0, 2) * epsxy) ^ epsxx;

  // ∫∫ 𝑪ᵗ₀₁ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<8, 8> sigmaXepsyy = (C_tang(0, 1) * epsxx + C_tang(1, 1) * epsyy + C_tang(1, 2) * epsxy) ^ epsyy;

  // ∫∫ 𝑪ᵗ₀₂ ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₁₂ ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₂₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
  RealMatrix<8, 8> sigmaXepsxy = (C_tang(0, 2) * epsxx + C_tang(1, 2) * epsyy + C_tang(2, 2) * epsxy) ^ epsxy;

  return integration_weight * (sigmaXepsxx + sigmaXepsyy + sigmaXepsxy);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline Real3x3
computeDisplacementGradientQuad4(Cell cell, const VariableNodeReal3& node_coord,
                                 const VariableNodeReal3& displacement, Real xi, Real eta)
{
  const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, node_coord, xi, eta);
  Real3x3 gradient{};
  for (Int32 i = 0; i < 4; ++i) {
    const Real3 u = displacement[cell.nodeId(i)];
    gradient(0, 0) += u.x * gp_info.dN_dx(i);
    gradient(0, 1) += u.x * gp_info.dN_dy(i);
    gradient(1, 0) += u.y * gp_info.dN_dx(i);
    gradient(1, 1) += u.y * gp_info.dN_dy(i);
  }
  return gradient;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline RealMatrix<8, 8> Elastoplasticity2Module::
_computeElementMatrixQuad4(Cell cell)
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

      RealMatrix<3, 3> C_tang_2d;
      for (Int8 ix = 0; ix < 3; ++ix) {
        for (Int8 iy = 0; iy < 3; ++iy) {
          C_tang_2d(ix, iy) = m_C_tang_gp(cell, iGP, ix, iy);
        }
      }
      iGP++;
      ae += computeElementMatrixQuad4Base(dxU, dyU, integration_weight, C_tang_2d);
    }
  }
  return ae;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD8, ℙ2 FE).
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
 *      a(𝐮,𝐯) =   ∫∫ 𝑪ᵗ₀₀ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ 𝑪ᵗ₀₁ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ 𝑪ᵗ₀₂ ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₁₂ ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₂₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂x)
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
inline ARCCORE_HOST_DEVICE RealMatrix<16, 16>
computeElementMatrixQuad8Base(const RealVector<8>& dxu, const RealVector<8>& dyu,
                              Real integration_weight, const RealMatrix<3, 3>& C_tang)
{
  RealVector<16> epsxx = { dxu(0), 0., dxu(1), 0., dxu(2), 0., dxu(3), 0.,
                           dxu(4), 0., dxu(5), 0., dxu(6), 0., dxu(7), 0. };
  RealVector<16> epsyy = { 0., dyu(0), 0., dyu(1), 0., dyu(2), 0., dyu(3),
                           0., dyu(4), 0., dyu(5), 0., dyu(6), 0., dyu(7) };
  RealVector<16> epsxy = { dyu(0), dxu(0), dyu(1), dxu(1), dyu(2), dxu(2), dyu(3), dxu(3),
                           dyu(4), dxu(4), dyu(5), dxu(5), dyu(6), dxu(6), dyu(7), dxu(7) };

  // Kelvin notation: scale shear strains by 1/sqrt(2)
  epsxy = M_SQRT1_2 * epsxy;

  // ∫∫ 𝑪ᵗ₀₀ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  const RealMatrix<16, 16> sigma_epsxx = (C_tang(0, 0) * epsxx + C_tang(0, 1) * epsyy + C_tang(0, 2) * epsxy) ^ epsxx;
  // ∫∫ 𝑪ᵗ₀₁ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  const RealMatrix<16, 16> sigma_epsyy = (C_tang(1, 0) * epsxx + C_tang(1, 1) * epsyy + C_tang(1, 2) * epsxy) ^ epsyy;
  // ∫∫ 𝑪ᵗ₀₂ ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₁₂ ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₂₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂x)
  const RealMatrix<16, 16> sigma_epsxy = (C_tang(2, 0) * epsxx + C_tang(2, 1) * epsyy + C_tang(2, 2) * epsxy) ^ epsxy;

  return integration_weight * (sigma_epsxx + sigma_epsyy + sigma_epsxy);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline Real3x3
computeDisplacementGradientQuad8(Cell cell, const VariableNodeReal3& node_coord,
                                 const VariableNodeReal3& displacement, Real xi, Real eta)
{
  const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad8(cell, node_coord, xi, eta);
  Real3x3 gradient{};
  for (Int32 i = 0; i < 8; ++i) {
    const Real3 u = displacement[cell.nodeId(i)];
    gradient(0, 0) += u.x * gp_info.dN_dx(i);
    gradient(0, 1) += u.x * gp_info.dN_dy(i);
    gradient(1, 0) += u.y * gp_info.dN_dx(i);
    gradient(1, 1) += u.y * gp_info.dN_dy(i);
  }
  return gradient;
}

inline RealMatrix<16, 16> Elastoplasticity2Module::
_computeElementMatrixQuad8(Cell cell)
{
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
  RealMatrix<16, 16> ae;
  ae.fill(0.0);
  Int8 iGP = 0;
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {
      const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad8(cell, m_node_coord, gp[ixi], gp[ieta]);
      RealMatrix<3, 3> C_tang;
      for (Int8 i = 0; i < 3; ++i)
        for (Int8 j = 0; j < 3; ++j)
          C_tang(i, j) = m_C_tang_gp(cell, iGP, i, j);
      ae += computeElementMatrixQuad8Base(gp_info.dN_dx, gp_info.dN_dy, gp_info.det_j * weights[ixi] * weights[ieta], C_tang);
      ++iGP;
    }
  }
  return ae;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD9, ℙ2 FE).
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
 *      a(𝐮,𝐯) =   ∫∫ 𝑪ᵗ₀₀ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ 𝑪ᵗ₀₁ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ 𝑪ᵗ₀₂ ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₁₂ ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂x) + 𝑪ᵗ₂₂ (∂𝑢𝑦/∂x + ∂𝑢x/∂y)(∂𝑣x/∂y + ∂𝑣y/∂x)
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

inline ARCCORE_HOST_DEVICE RealMatrix<18, 18>
computeElementMatrixQuad9Base(const RealVector<9>& dxu, const RealVector<9>& dyu, Real integration_weight, const RealMatrix<3, 3>& C_tang)
{
  RealVector<18> epsxx = { dxu(0), 0., dxu(1), 0., dxu(2), 0., dxu(3), 0., dxu(4), 0.,
                           dxu(5), 0., dxu(6), 0., dxu(7), 0., dxu(8), 0. };

  RealVector<18> epsyy = { 0., dyu(0), 0., dyu(1), 0., dyu(2), 0., dyu(3), 0., dyu(4),
                           0., dyu(5), 0., dyu(6), 0., dyu(7), 0., dyu(8) };

  RealVector<18> epsxy = { dyu(0), dxu(0), dyu(1), dxu(1), dyu(2), dxu(2), dyu(3), dxu(3), dyu(4), dxu(4),
                           dyu(5), dxu(5), dyu(6), dxu(6), dyu(7), dxu(7), dyu(8), dxu(8) };

  // Kelvin notation: scale shear strains by 1/sqrt(2)
  epsxy = M_SQRT1_2 * epsxy;

  // ∫∫ 𝑪ᵗ₀₀ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + 𝑪ᵗ₀₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  const RealMatrix<18, 18> sigma_epsxx = (C_tang(0, 0) * epsxx + C_tang(0, 1) * epsyy + C_tang(0, 2) * epsxy) ^ epsxx;
  // ∫∫ 𝑪ᵗ₀₁ ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₁ ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + 𝑪ᵗ₁₂ (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  const RealMatrix<18, 18> sigma_epsyy = (C_tang(1, 0) * epsxx + C_tang(1, 1) * epsyy + C_tang(1, 2) * epsxy) ^ epsyy;
  // ∫∫ 𝑪ᵗ₀₂ ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + 𝑪ᵗ₁₂ ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂x) + 𝑪ᵗ₂₂ (∂𝑢y/∂x + ∂𝑢x/∂y)(∂𝑣x/∂y + ∂𝑣y/∂x)
  const RealMatrix<18, 18> sigma_epsxy = (C_tang(2, 0) * epsxx + C_tang(2, 1) * epsyy + C_tang(2, 2) * epsxy) ^ epsxy;

  return integration_weight * (sigma_epsxx + sigma_epsyy + sigma_epsxy);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline Real3x3
computeDisplacementGradientQuad9(Cell cell, const VariableNodeReal3& node_coord,
                                 const VariableNodeReal3& displacement, Real xi, Real eta)
{
  const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad9(cell, node_coord, xi, eta);
  Real3x3 gradient{};
  for (Int32 i = 0; i < 9; ++i) {
    const Real3 u = displacement[cell.nodeId(i)];
    gradient(0, 0) += u.x * gp_info.dN_dx(i);
    gradient(0, 1) += u.x * gp_info.dN_dy(i);
    gradient(1, 0) += u.y * gp_info.dN_dx(i);
    gradient(1, 1) += u.y * gp_info.dN_dy(i);
  }
  return gradient;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

RealMatrix<18, 18> Elastoplasticity2Module::
_computeElementMatrixQuad9(Cell cell)
{
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
  RealMatrix<18, 18> ae;
  ae.fill(0.0);
  Int8 iGP = 0;
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {
      const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad9(cell, m_node_coord, gp[ixi], gp[ieta]);
      RealMatrix<3, 3> C_tang;
      for (Int8 i = 0; i < 3; ++i)
        for (Int8 j = 0; j < 3; ++j)
          C_tang(i, j) = m_C_tang_gp(cell, iGP, i, j);
      ae += computeElementMatrixQuad9Base(gp_info.dN_dx, gp_info.dN_dy, gp_info.det_j * weights[ixi] * weights[ieta], C_tang);
      ++iGP;
    }
  }
  return ae;
}

/*---------------------------------------------------------------------------*/
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

inline ARCCORE_HOST_DEVICE RealMatrix<24, 24>
computeElementMatrixHexa8Base(const RealVector<8>& dxu, const RealVector<8>& dyu, const RealVector<8>& dzu,
                              Real integration_weight, RealMatrix<6, 6> C_tang)
{
  RealVector<24> epsxx = { dxu(0), 0., 0., dxu(1), 0., 0., dxu(2), 0., 0., dxu(3), 0., 0.,
                           dxu(4), 0., 0., dxu(5), 0., 0., dxu(6), 0., 0., dxu(7), 0., 0. };

  RealVector<24> epsyy = { 0., dyu(0), 0., 0., dyu(1), 0., 0., dyu(2), 0., 0., dyu(3), 0.,
                           0., dyu(4), 0., 0., dyu(5), 0., 0., dyu(6), 0., 0., dyu(7), 0. };

  RealVector<24> epszz = { 0., 0., dzu(0), 0., 0., dzu(1), 0., 0., dzu(2), 0., 0., dzu(3),
                           0., 0., dzu(4), 0., 0., dzu(5), 0., 0., dzu(6), 0., 0., dzu(7) };

  RealVector<24> epsyz = { 0., dzu(0), dyu(0), 0., dzu(1), dyu(1),
                           0., dzu(2), dyu(2), 0., dzu(3), dyu(3),
                           0., dzu(4), dyu(4), 0., dzu(5), dyu(5),
                           0., dzu(6), dyu(6), 0., dzu(7), dyu(7) };

  RealVector<24> epszx = { dzu(0), 0., dxu(0), dzu(1), 0., dxu(1),
                           dzu(2), 0., dxu(2), dzu(3), 0., dxu(3),
                           dzu(4), 0., dxu(4), dzu(5), 0., dxu(5),
                           dzu(6), 0., dxu(6), dzu(7), 0., dxu(7) };

  RealVector<24> epsxy = { dyu(0), dxu(0), 0., dyu(1), dxu(1), 0.,
                           dyu(2), dxu(2), 0., dyu(3), dxu(3), 0.,
                           dyu(4), dxu(4), 0., dyu(5), dxu(5), 0.,
                           dyu(6), dxu(6), 0., dyu(7), dxu(7), 0. };

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

  return integration_weight * (sigmaXepsxx + sigmaXepsyy + sigmaXepszz + sigmaXepsyz + sigmaXepszx + sigmaXepsxy);
}

RealMatrix<24, 24> Elastoplasticity2Module::_computeElementMatrixHexa8(Cell cell)
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
        RealMatrix<6, 6> C_tang_3d;
        for (Int32 ix = 0; ix < 6; ++ix) {
          for (Int32 iy = 0; iy < 6; ++iy) {
            C_tang_3d(ix, iy) = m_C_tang_gp(cell, iGP, ix, iy);
          }
        }
        iGP++;
        ae += computeElementMatrixHexa8Base(dxU, dyU, dzU, integration_weight, C_tang_3d);
      }
    }
  }

  return ae;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::ArcaneFem

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
