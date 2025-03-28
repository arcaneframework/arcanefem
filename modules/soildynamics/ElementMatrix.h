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
 * @brief Computes 2D problem element matrix for a triangular element (ℙ1 FE)
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
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<6, 6> computeElementMatrixTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real c0, Real c1, Real c2)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

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
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<2, 6> computeElementVectorTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real c0, Real c1, Real c2, Int32 node_lid)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  RealVector<6> Uy = {0., 1., 0., 1., 0., 1.};
  RealVector<6> Ux = {1., 0., 1., 0., 1., 0.};

  RealVector<6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  RealVector<6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
  RealVector<6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
  RealVector<6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

  RealMatrix<6, 6> massMat = (c0 / 12.) * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy)) * area;
  RealVector<6> massVect_x = {massMat(node_lid*2,0) , massMat(node_lid*2,1) , massMat(node_lid*2,2) ,
                              massMat(node_lid*2,3) , massMat(node_lid*2,4) , massMat(node_lid*2,5) };

  RealVector<6> massVect_y = {massMat(node_lid*2+1,0) , massMat(node_lid*2+1,1) , massMat(node_lid*2+1,2) ,
                              massMat(node_lid*2+1,3) , massMat(node_lid*2+1,4) , massMat(node_lid*2+1,5)};

  RealVector <6> result_x = (c1) * ((dyUy(node_lid*2) * dxUx) + (dxUx(node_lid*2) * dyUy)) * area +
                            (2*c2 + c1) * ((dxUx(node_lid*2) * dxUx) + (dyUy(node_lid*2) * dyUy)) * area +
                            (c2) * ((dxUy(node_lid*2) + dyUx(node_lid*2)) * (dyUx + dxUy)) * area;
  RealVector <6> result_y = (c1) * ((dyUy(node_lid*2+1) * dxUx) + (dxUx(node_lid*2+1) * dyUy)) * area + 
                            (2*c2 + c1) * ((dxUx(node_lid*2+1) * dxUx) + (dyUy(node_lid*2+1) * dyUy)) * area +
                            (c2) * ((dxUy(node_lid*2+1) + dyUx(node_lid*2+1)) * (dyUx + dxUy)) * area;

  result_x = result_x + massVect_x;
  result_y = result_y + massVect_y;

  RealMatrix<2, 6> result = {
    { result_x(0), result_x(1), result_x(2), result_x(3), result_x(4), result_x(5) },
    { result_y(0), result_y(1), result_y(2), result_y(3), result_y(4), result_y(5) }
  };

  return result;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes 2D paraxial element matrix for a edge element (ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) =  ∫ (cₚ 𝑁𝑥² + cₛ 𝑁𝑦²)(𝑢𝑥 𝑣𝑥) +
 *             ∫ (cₚ 𝑁𝑦² + cₛ 𝑁𝑥²)(𝑢𝑦 𝑣𝑦) +
 *             ∫ (𝑁𝑦 𝑁𝑥 (cₚ - cₛ))(𝑢𝑥 𝑣𝑦) +
 *             ∫ (𝑁𝑥 𝑁𝑦 (cₚ - cₛ))(𝑢𝑦 𝑣𝑥) ;
 *
 *   with  trial func 𝐮 = (𝑢𝑥,𝑢𝑦) and test func 𝐯 = (𝑣𝑥,𝑣𝑦)
 */
/*---------------------------------------------------------------------------*/

RealMatrix<4, 4> FemModule::
_computeParaxialElementMatrixEdge2(Face face)
{
  Real2 N   = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, m_node_coord);

  RealVector<4> Uy = {0., 1., 0., 1.};
  RealVector<4> Ux = {1., 0., 1., 0.};

  RealMatrix<4, 4> int_Omega_i = (((N.x*N.x*cp + N.y*N.y*cs)) * (massMatrix(Ux,Ux)) +
                                  ((N.y*N.y*cp + N.x*N.x*cs)) * (massMatrix(Uy,Uy)) +
                                  ((N.x*N.y*(cp - cs))) * (massMatrix(Ux,Uy)) +
                                  ((N.x*N.y*(cp - cs))) * (massMatrix(Uy,Ux)) )/6. ;
  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes 3D problem element matrix for a tetrahedral element (ℙ1 FE).
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

ARCCORE_HOST_DEVICE RealMatrix<12, 12> computeElementMatrixTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real c0, Real c1, Real c2)
{
  Real4 dxu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);

  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

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

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes 3D paraxial element matrix for a triangular element (ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) =  ∫∫ (cₚ 𝑁𝑥² + cₛ (1 - 𝑁𝑥²))(𝑢𝑥 𝑣𝑥) +
 *             ∫∫ (cₚ 𝑁𝑦² + cₛ (1 - 𝑁𝑦²))(𝑢𝑦 𝑣𝑦) +
 *             ∫∫ (cₚ 𝑁𝑧² + cₛ (1 - 𝑁𝑧²))(𝑢𝑧 𝑣𝑧) +
 *             ∫∫ (𝑁𝑥 𝑁𝑦 (cₚ - cₛ))(𝑢𝑥 𝑣𝑦) + ∫∫ (c₇)(𝑁𝑥 𝑁𝑧 (cₚ - cₛ))(𝑢𝑥 𝑣𝑧) +
 *             ∫∫ (𝑁𝑦 𝑁𝑥 (cₚ - cₛ))(𝑢𝑦 𝑣𝑥) + ∫∫ (c₇)(𝑁𝑦 𝑁𝑧 (cₚ - cₛ))(𝑢𝑦 𝑣𝑧) +
 *             ∫∫ (𝑁𝑧 𝑁𝑥 (cₚ - cₛ))(𝑢𝑧 𝑣𝑥) + ∫∫ (c₇)(𝑁𝑧 𝑁𝑦 (cₚ - cₛ))(𝑢𝑧 𝑣𝑦) ;
 *
 *   with trial function 𝐮 = (𝑢𝑥, 𝑢𝑦, 𝑢𝑧) and test function 𝐯 = (𝑣𝑥, 𝑣𝑦, 𝑣𝑧)
 */
/*---------------------------------------------------------------------------*/

RealMatrix<9, 9> FemModule::
_computeParaxialElementMatrixTria3(Face face)
{
  Real3 N = ArcaneFemFunctions::MeshOperation::computeNormalTriangle(face, m_node_coord);

  RealVector<9> Ux = {1., 0., 0., 1., 0., 0., 1., 0., 0.};
  RealVector<9> Uy = {0., 1., 0., 0., 1., 0., 0., 1., 0.};
  RealVector<9> Uz = {0., 0., 1., 0., 0., 1., 0., 0., 1.};

  RealMatrix<9, 9> int_Omega_i =  (((N.x*N.x*cp + (1.-N.x*N.x)*cs)) * (massMatrix(Ux,Ux)) +
                                   ((N.y*N.y*cp + (1.-N.y*N.y)*cs)) * (massMatrix(Uy,Uy)) +
                                   ((N.z*N.z*cp + (1.-N.z*N.z)*cs)) * (massMatrix(Uz,Uz)) +
                                   ((N.x*N.y*(cp - cs))) * (massMatrix(Ux,Uy))  +
                                   ((N.x*N.z*(cp - cs))) * (massMatrix(Ux,Uz))  +
                                   ((N.y*N.x*(cp - cs))) * (massMatrix(Uy,Ux))  +
                                   ((N.y*N.z*(cp - cs))) * (massMatrix(Uy,Uz))  +
                                   ((N.z*N.x*(cp - cs))) * (massMatrix(Uz,Ux))  +
                                   ((N.z*N.y*(cp - cs))) * (massMatrix(Uz,Uy)) ) / 12.
                                  ;

  return int_Omega_i;
}