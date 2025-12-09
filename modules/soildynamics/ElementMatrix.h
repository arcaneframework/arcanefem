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

RealMatrix<6, 6> FemModuleSoildynamics::
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

RealMatrix<12, 12> FemModuleSoildynamics::_computeElementMatrixTetra4(Cell cell)
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
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<3, 12> computeElementVectorTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord,  Real c0, Real c1, Real c2, Int32 node_lid)
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

  RealMatrix<12, 12> massMat = (c0 / 20.) * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy) + massMatrix(Uz,Uz)) * volume;
  RealVector<12> massVect_x = {massMat(node_lid*3,0) , massMat(node_lid*3,1) , massMat(node_lid*3,2) ,
                               massMat(node_lid*3,3) , massMat(node_lid*3,4) , massMat(node_lid*3,5) ,
                               massMat(node_lid*3,6) , massMat(node_lid*3,7) , massMat(node_lid*3,8) ,
                               massMat(node_lid*3,9) , massMat(node_lid*3,10), massMat(node_lid*3,11) };
  RealVector<12> massVect_y = {massMat(node_lid*3+1,0) , massMat(node_lid*3+1,1) , massMat(node_lid*3+1,2) ,
                               massMat(node_lid*3+1,3) , massMat(node_lid*3+1,4) , massMat(node_lid*3+1,5) ,
                               massMat(node_lid*3+1,6) , massMat(node_lid*3+1,7) , massMat(node_lid*3+1,8) ,
                               massMat(node_lid*3+1,9) , massMat(node_lid*3+1,10), massMat(node_lid*3+1,11) };
  RealVector<12> massVect_z =  {massMat(node_lid*3+2,0) , massMat(node_lid*3+2,1) , massMat(node_lid*3+2,2) ,
                                massMat(node_lid*3+2,3) , massMat(node_lid*3+2,4) , massMat(node_lid*3+2,5) ,
                                massMat(node_lid*3+2,6) , massMat(node_lid*3+2,7) , massMat(node_lid*3+2,8) ,
                                massMat(node_lid*3+2,9) , massMat(node_lid*3+2,10), massMat(node_lid*3+2,11) };

  RealVector<12> result_x = (c1)*((dxUx(node_lid*3) * dxUx) + (dyUy(node_lid*3) * dyUy) + (dzUz (node_lid*3) * dzUz) +
                                  (dyUy(node_lid*3) * dxUx) + (dxUx(node_lid*3) * dyUy) +
                                  (dzUz(node_lid*3) * dxUx) + (dxUx(node_lid*3) * dzUz) +
                                  (dyUy(node_lid*3) * dzUz) + (dzUz(node_lid*3) * dyUy) ) * volume +
                            (c2)*(2.*((dxUx(node_lid*3) * dxUx) + (dyUy(node_lid*3) * dyUy) + (dzUz(node_lid*3) * dzUz) ) +
                                 (((dxUy(node_lid*3) + dyUx(node_lid*3)) * (dyUx + dxUy)) +
                                  ((dzUy(node_lid*3) + dyUz(node_lid*3)) * (dyUz + dzUy)) +
                                  ((dxUz(node_lid*3) + dzUx(node_lid*3)) * (dzUx + dxUz)) ) )*volume;

  RealVector<12> result_y = (c1)*((dxUx(node_lid*3+1) * dxUx) + (dyUy(node_lid*3+1) * dyUy) + (dzUz (node_lid*3+1) * dzUz) +
                                (dyUy(node_lid*3+1) * dxUx) + (dxUx(node_lid*3+1) * dyUy) +
                                (dzUz(node_lid*3+1) * dxUx) + (dxUx(node_lid*3+1) * dzUz) +
                                (dyUy(node_lid*3+1) * dzUz) + (dzUz(node_lid*3+1) * dyUy) ) * volume +
                          (c2)*(2.*((dxUx(node_lid*3+1) * dxUx) + (dyUy(node_lid*3+1) * dyUy) + (dzUz(node_lid*3+1) * dzUz) ) +
                                (((dxUy(node_lid*3+1) + dyUx(node_lid*3+1)) * (dyUx + dxUy)) +
                                 ((dzUy(node_lid*3+1) + dyUz(node_lid*3+1)) * (dyUz + dzUy)) +
                                 ((dxUz(node_lid*3+1) + dzUx(node_lid*3+1)) * (dzUx + dxUz)) ) )*volume;

  RealVector<12> result_z = (c1)*((dxUx(node_lid*3+2) * dxUx) + (dyUy(node_lid*3+2) * dyUy) + (dzUz (node_lid*3+2) * dzUz) +
                                (dyUy(node_lid*3+2) * dxUx) + (dxUx(node_lid*3+2) * dyUy) +
                                (dzUz(node_lid*3+2) * dxUx) + (dxUx(node_lid*3+2) * dzUz) +
                                (dyUy(node_lid*3+2) * dzUz) + (dzUz(node_lid*3+2) * dyUy) ) * volume +
                          (c2)*(2.*((dxUx(node_lid*3+2) * dxUx) + (dyUy(node_lid*3+2) * dyUy) + (dzUz(node_lid*3+2) * dzUz) ) +
                                (((dxUy(node_lid*3+2) + dyUx(node_lid*3+2)) * (dyUx + dxUy)) +
                                 ((dzUy(node_lid*3+2) + dyUz(node_lid*3+2)) * (dyUz + dzUy)) +
                                 ((dxUz(node_lid*3+2) + dzUx(node_lid*3+2)) * (dzUx + dxUz)) ) )*volume;

  result_x = result_x + massVect_x;
  result_y = result_y + massVect_y;
  result_z = result_z + massVect_z;

  RealMatrix<3, 12> result = {
    { result_x(0), result_x(1), result_x(2), result_x(3), result_x(4), result_x(5),
      result_x(6), result_x(7), result_x(8), result_x(9), result_x(10), result_x(11) },

    { result_y(0), result_y(1), result_y(2), result_y(3), result_y(4), result_y(5),
      result_y(6), result_y(7), result_y(8), result_y(9), result_y(10), result_y(11) },

    { result_z(0), result_z(1), result_z(2), result_z(3), result_z(4), result_z(5),
      result_z(6), result_z(7), result_z(8), result_z(9), result_z(10), result_z(11) }
  };

  return result;
}