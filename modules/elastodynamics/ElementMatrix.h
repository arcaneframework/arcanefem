﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrix.h                                             (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Elastodynamics */
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

RealMatrix<6, 6> FemModule::
_computeElementMatrixTria3(Cell cell)
{
  Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  RealVector<6> Uy = {0., 1., 0., 1., 0., 1.};
  RealVector<6> Ux = {1., 0., 1., 0., 1., 0.};

  RealVector<6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  RealVector<6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
  RealVector<6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
  RealVector<6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

  RealMatrix<6, 6> int_Omega_i = (c0 / 12.) * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy)) * area +
                                  (c1) * ((dyUy ^ dxUx) + (dxUx ^ dyUy)) * area +
                                  (2*c2 + c1) * ((dxUx ^ dxUx) + (dyUy ^ dyUy)) * area +
                                  (c2) * ((dxUy + dyUx) ^ (dyUx + dxUy)) * area;

  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a tetrahedral element (ℙ1 FE).
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
 */
/*---------------------------------------------------------------------------*/

RealMatrix<12, 12> FemModule::_computeElementMatrixTetra4(Cell cell)
{
  Real4 dxu = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
  Real4 dyu = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
  Real4 dzu = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);

  Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);

  RealVector<12> Uy = { 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0. };
  RealVector<12> Ux = { 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0. };
  RealVector<12> Uz = { 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1. };

  RealVector<12> dxUx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
  RealVector<12> dyUx = { dyu[0], 0., 0.,    dyu[1], 0., 0.,    dyu[2], 0., 0.,    dyu[3], 0., 0. };
  RealVector<12> dzUx = { dzu[0], 0., 0.,    dzu[1], 0., 0.,    dzu[2], 0., 0.,    dzu[3], 0., 0. };

  RealVector<12> dxUy = { 0., dxu[0], 0.,    0., dxu[1], 0.,    0., dxu[2], 0.,    0., dxu[3], 0. };
  RealVector<12> dyUy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
  RealVector<12> dzUy = { 0., dzu[0], 0.,    0., dzu[1], 0.,    0., dzu[2], 0.,    0., dzu[3], 0. };

  RealVector<12> dxUz = { 0., 0., dxu[0],    0., 0., dxu[1],    0., 0., dxu[2],    0., 0., dxu[3] };
  RealVector<12> dyUz = { 0., 0., dyu[0],    0., 0., dyu[1],    0., 0., dyu[2],    0., 0., dyu[3] };
  RealVector<12> dzUz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

  RealMatrix<12, 12> int_Omega_i = (c0 / 20.) * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy) + massMatrix(Uz,Uz)) * volume +
                                    (c1)*((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) +
                                          (dyUy ^ dxUx) + (dxUx ^ dyUy) +
                                          (dzUz ^ dxUx) + (dxUx ^ dzUz) +
                                          (dyUy ^ dzUz) + (dzUz ^ dyUy) ) * volume +
                                    (c2)*(2.*((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) ) +
                                          ( ((dxUy + dyUx) ^ (dyUx + dxUy)) +
                                            ((dzUy + dyUz) ^ (dyUz + dzUy)) +
                                            ((dxUz + dzUx) ^ (dzUx + dxUz)) ) )*volume;

  return int_Omega_i;
}