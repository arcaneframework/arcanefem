// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrix.h                                             (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Soildynamics   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (ℙ1 FE).
 *
 * Theory:
 *
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
 */
/*---------------------------------------------------------------------------*/

FixedMatrix<6, 6> FemModule::
_computeElementMatrixTria3(Cell cell)
{
  Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  FixedMatrix<1, 6> Uy = {0., 1., 0., 1., 0., 1.};
  FixedMatrix<1, 6> Ux = {1., 0., 1., 0., 1., 0.};
  FixedMatrix<1, 6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  FixedMatrix<1, 6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
  FixedMatrix<1, 6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
  FixedMatrix<1, 6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

  FixedMatrix<6, 6> int_Omega_i = (c0 / 12.) * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy)) * area +
                                  (c1) * ((dyUy ^ dxUx) + (dxUx ^ dyUy)) * area +
                                  (2*c2 + c1) * ((dxUx ^ dxUx) + (dyUy ^ dyUy)) * area +
                                  (c2) * ((dxUy + dyUx) ^ (dyUx + dxUy)) * area;

  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a edge element (ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) =  ∫ (c₇)(cₚ 𝑁𝑥² + cₛ 𝑁𝑦²)(𝑢𝑥 𝑣𝑥) +
 *             ∫ (c₇)(cₚ 𝑁𝑦² + cₛ 𝑁𝑥²)(𝑢𝑦 𝑣𝑦) +
 *             ∫ (c₇)(𝑁𝑦 𝑁𝑥 (cₚ - cₛ))(𝑢𝑥 𝑣𝑦) +
 *             ∫ (c₇)(𝑁𝑥 𝑁𝑦 (cₚ - cₛ))(𝑢𝑦 𝑣𝑥) ;
 */
/*---------------------------------------------------------------------------*/

FixedMatrix<4, 4> FemModule::
_computeElementMatrixEdge2(Face face)
{
  Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);
  Real2 N   = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, m_node_coord);

  FixedMatrix<1, 4> Uy = {0., 1., 0., 1.};
  FixedMatrix<1, 4> Ux = {1., 0., 1., 0.};

  FixedMatrix<4, 4> int_Omega_i = (c7*(N.x*N.x*cp + N.y*N.y*cs)) * (massMatrix(Ux,Ux)) * length/6. +
                                  (c7*(N.y*N.y*cp + N.x*N.x*cs)) * (massMatrix(Uy,Uy)) * length/6. +
                                  (c7*(N.x*N.y*(cp - cs))) * (massMatrix(Ux,Uy)) * length/6. +
                                  (c7*(N.x*N.y*(cp - cs))) * (massMatrix(Uy,Ux)) * length/6. ;
  return int_Omega_i;
}