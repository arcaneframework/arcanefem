// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrix.h                                             (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Elasticity     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) = ∫∫ [σ(𝐮):ε(𝐯)dΩ    with  𝐮 = (𝑢𝑥,𝑢𝑦) and 𝐯 = (𝑣𝑥,𝑣𝑦)
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
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<6, 6> computeElementMatrixTria3Base(Real3 dxu, Real3 dyu, Real area, Real lambda, Real mu)
{
  FixedMatrix<1, 6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  FixedMatrix<1, 6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
  FixedMatrix<1, 6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
  FixedMatrix<1, 6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

  // ∫∫ (λ+2μ)(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦)
  FixedMatrix<6, 6> normal_strain_energy = (lambda + 2 * mu) * ((dxUx ^ dxUx) + (dyUy ^ dyUy)) * area;

  // ∫∫ λ(∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦)
  FixedMatrix<6, 6> compressibility_effect = (lambda) * ((dyUy ^ dxUx) + (dxUx ^ dyUy)) * area;

  // ∫∫ μ(∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
  FixedMatrix<6, 6> shear_energy = (mu) * ((dxUy + dyUx) ^ (dyUx + dxUy)) * area;

  return (normal_strain_energy + compressibility_effect + shear_energy);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<6, 6> computeElementMatrixTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  return computeElementMatrixTria3Base(dxu, dyu, area, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<6, 6> FemModule::_computeElementMatrixTria3(Cell cell)
{
  Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  return computeElementMatrixTria3Base(dxu, dyu, area, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<2, 6> computeElementVectorTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu, Int32 node_lid)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  FixedMatrix<1, 6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  FixedMatrix<1, 6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
  FixedMatrix<1, 6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
  FixedMatrix<1, 6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

  // ∫∫ (λ+2μ)(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦)
  FixedMatrix<1, 6> normal_strain_energy_x = (lambda + 2 * mu) * ((dxUx(0,node_lid*2) * dxUx) + (dyUy(0,node_lid*2) * dyUy)) * area;
  FixedMatrix<1, 6> normal_strain_energy_y = (lambda + 2 * mu) * ((dxUx(0,node_lid*2+1) * dxUx) + (dyUy(0,node_lid*2+1) * dyUy)) * area;

  // ∫∫ λ(∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦)
  FixedMatrix<1, 6> compressibility_effect_x = (lambda) * ((dyUy(0,node_lid*2) * dxUx) + (dxUx(0,node_lid*2) * dyUy)) * area;
  FixedMatrix<1, 6> compressibility_effect_y = (lambda) * ((dyUy(0,node_lid*2+1) * dxUx) + (dxUx(0,node_lid*2+1) * dyUy)) * area;

  // ∫∫ μ(∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
  FixedMatrix<1, 6> shear_energy_x = (mu) * ((dxUy(0,node_lid*2) + dyUx(0,node_lid*2)) * (dyUx + dxUy)) * area;
  FixedMatrix<1, 6> shear_energy_y = (mu) * ((dxUy(0,node_lid*2+1) + dyUx(0,node_lid*2+1)) * (dyUx + dxUy)) * area;

  FixedMatrix <1, 6> result_x = normal_strain_energy_x + compressibility_effect_x + shear_energy_x;
  FixedMatrix <1, 6> result_y = normal_strain_energy_y + compressibility_effect_y + shear_energy_y;

  FixedMatrix <2, 6> result;
  result(0,0)=result_x(0,0); result(0,1)=result_x(0,1); result(0,2)=result_x(0,2); result(0,3)=result_x(0,3); result(0,4)=result_x(0,4); result(0,5)=result_x(0,5);
  result(1,0)=result_y(0,0); result(1,1)=result_y(0,1); result(1,2)=result_y(0,2); result(1,3)=result_y(0,3); result(1,4)=result_y(0,4); result(1,5)=result_y(0,5);

  return result;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a tetrahedral element (ℙ1 FE).
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

ARCCORE_HOST_DEVICE FixedMatrix<12, 12> computeElementMatrixTetra4Base(Real4 dxu, Real4 dyu, Real4 dzu, Real volume, Real lambda, Real mu)
{
  FixedMatrix<1, 12> dxUx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
  FixedMatrix<1, 12> dyUx = { dyu[0], 0., 0.,    dyu[1], 0., 0.,    dyu[2], 0., 0.,    dyu[3], 0., 0. };
  FixedMatrix<1, 12> dzUx = { dzu[0], 0., 0.,    dzu[1], 0., 0.,    dzu[2], 0., 0.,    dzu[3], 0., 0. };

  FixedMatrix<1, 12> dxUy = { 0., dxu[0], 0.,    0., dxu[1], 0.,    0., dxu[2], 0.,    0., dxu[3], 0. };
  FixedMatrix<1, 12> dyUy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
  FixedMatrix<1, 12> dzUy = { 0., dzu[0], 0.,    0., dzu[1], 0.,    0., dzu[2], 0.,    0., dzu[3], 0. };

  FixedMatrix<1, 12> dxUz = { 0., 0., dxu[0],    0., 0., dxu[1],    0., 0., dxu[2],    0., 0., dxu[3] };
  FixedMatrix<1, 12> dyUz = { 0., 0., dyu[0],    0., 0., dyu[1],    0., 0., dyu[2],    0., 0., dyu[3] };
  FixedMatrix<1, 12> dzUz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

  // ∫∫∫ (λ+2μ)(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧)
  FixedMatrix<12, 12> normal_strain_energy = (lambda + 2 * mu) * ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) ) * volume;

  // ∫∫∫ λ(∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦
  //     + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧
  //     + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑧/∂𝑧 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦)
  FixedMatrix<12, 12> compressibility_effect = (lambda) * ((dyUy ^ dxUx) + (dxUx ^ dyUy) +
                                                           (dzUz ^ dxUx) + (dxUx ^ dzUz) +
                                                           (dyUy ^ dzUz) + (dzUz ^ dyUy) ) * volume;

  // ∫∫∫ μ(∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) +
  //     μ(∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)(∂𝑣𝑦/∂𝑧 + ∂𝑣𝑧/∂𝑦) +
  //     μ(∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)(∂𝑣𝑧/∂𝑥 + ∂𝑣𝑥/∂𝑧)
  FixedMatrix<12, 12> shear_energy = (mu) * ( ((dxUy + dyUx) ^ (dyUx + dxUy)) +
                                              ((dzUy + dyUz) ^ (dyUz + dzUy)) +
                                              ((dxUz + dzUx) ^ (dzUx + dxUz)) ) * volume;

  return (normal_strain_energy + compressibility_effect + shear_energy);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<12, 12> FemModule::_computeElementMatrixTetra4(Cell cell)
{
  Real4 dxu = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
  Real4 dyu = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
  Real4 dzu = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);

  Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);

  return computeElementMatrixTetra4Base(dxu, dyu, dzu, volume, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<12, 12> computeElementMatrixTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu)
{
  Real4 dxu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);

  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  return computeElementMatrixTetra4Base(dxu, dyu, dzu, volume, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<3, 12> computeElementVectorTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu, Int32 node_lid)
{
  Real4 dxu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);

  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  FixedMatrix<1, 12> dxUx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
  FixedMatrix<1, 12> dyUx = { dyu[0], 0., 0.,    dyu[1], 0., 0.,    dyu[2], 0., 0.,    dyu[3], 0., 0. };
  FixedMatrix<1, 12> dzUx = { dzu[0], 0., 0.,    dzu[1], 0., 0.,    dzu[2], 0., 0.,    dzu[3], 0., 0. };

  FixedMatrix<1, 12> dxUy = { 0., dxu[0], 0.,    0., dxu[1], 0.,    0., dxu[2], 0.,    0., dxu[3], 0. };
  FixedMatrix<1, 12> dyUy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
  FixedMatrix<1, 12> dzUy = { 0., dzu[0], 0.,    0., dzu[1], 0.,    0., dzu[2], 0.,    0., dzu[3], 0. };

  FixedMatrix<1, 12> dxUz = { 0., 0., dxu[0],    0., 0., dxu[1],    0., 0., dxu[2],    0., 0., dxu[3] };
  FixedMatrix<1, 12> dyUz = { 0., 0., dyu[0],    0., 0., dyu[1],    0., 0., dyu[2],    0., 0., dyu[3] };
  FixedMatrix<1, 12> dzUz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

 // ∫∫∫ (λ+2μ)(∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 +
 //            ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 +
 //            ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧 )
 FixedMatrix<1, 12> normal_strain_energy_x = (lambda + 2 * mu) * ( (dxUx(0,node_lid*3) * dxUx) +
                                                                   (dyUy(0,node_lid*3) * dyUy) +
                                                                   (dzUz(0,node_lid*3) * dzUz) ) * volume;

 FixedMatrix<1, 12> normal_strain_energy_y = (lambda + 2 * mu) * ( (dxUx(0,node_lid*3+1) * dxUx) +
                                                                   (dyUy(0,node_lid*3+1) * dyUy) +
                                                                   (dzUz(0,node_lid*3+1) * dzUz) ) * volume;

 FixedMatrix<1, 12> normal_strain_energy_z = (lambda + 2 * mu) * ( (dxUx(0,node_lid*3+2) * dxUx) +
                                                                   (dyUy(0,node_lid*3+2) * dyUy) +
                                                                   (dzUz(0,node_lid*3+2) * dzUz) ) * volume;

  // ∫∫∫ λ(∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦
  //     + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧
  //     + ∂𝑢𝑦/∂𝑦 ∂𝑣𝑧/∂𝑧 + ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦)
  FixedMatrix<1, 12> compressibility_effect_x = (lambda) * ((dyUy(0,node_lid*3) * dxUx) + (dxUx(0,node_lid*3) * dyUy) +
                                                            (dzUz(0,node_lid*3) * dxUx) + (dxUx(0,node_lid*3) * dzUz) +
                                                            (dyUy(0,node_lid*3) * dzUz) + (dzUz(0,node_lid*3) * dyUy) ) * volume;

  FixedMatrix<1, 12> compressibility_effect_y = (lambda) * ((dyUy(0,node_lid*3+1) * dxUx) + (dxUx(0,node_lid*3+1) * dyUy) +
                                                            (dzUz(0,node_lid*3+1) * dxUx) + (dxUx(0,node_lid*3+1) * dzUz) +
                                                            (dyUy(0,node_lid*3+1) * dzUz) + (dzUz(0,node_lid*3+1) * dyUy) ) * volume;

  FixedMatrix<1, 12> compressibility_effect_z = (lambda) * ((dyUy(0,node_lid*3+2) * dxUx) + (dxUx(0,node_lid*3+2) * dyUy) +
                                                            (dzUz(0,node_lid*3+2) * dxUx) + (dxUx(0,node_lid*3+2) * dzUz) +
                                                            (dyUy(0,node_lid*3+2) * dzUz) + (dzUz(0,node_lid*3+2) * dyUy) ) * volume;


  // ∫∫∫ μ(∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) +
  //     μ(∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)(∂𝑣𝑦/∂𝑧 + ∂𝑣𝑧/∂𝑦) +
  //     μ(∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)(∂𝑣𝑧/∂𝑥 + ∂𝑣𝑥/∂𝑧)
  FixedMatrix<1, 12> shear_energy_x = (mu) * ( ((dxUy(0,node_lid*3) + dyUx(0,node_lid*3)) * (dyUx + dxUy)) +
                                               ((dzUy(0,node_lid*3) + dyUz(0,node_lid*3)) * (dyUz + dzUy)) +
                                               ((dxUz(0,node_lid*3) + dzUx(0,node_lid*3)) * (dzUx + dxUz)) ) * volume;

  FixedMatrix<1, 12> shear_energy_y = (mu) * ( ((dxUy(0,node_lid*3+1) + dyUx(0,node_lid*3+1)) * (dyUx + dxUy)) +
                                               ((dzUy(0,node_lid*3+1) + dyUz(0,node_lid*3+1)) * (dyUz + dzUy)) +
                                               ((dxUz(0,node_lid*3+1) + dzUx(0,node_lid*3+1)) * (dzUx + dxUz)) ) * volume;



  FixedMatrix<1, 12> shear_energy_z = (mu) * ( ((dxUy(0,node_lid*3+2) + dyUx(0,node_lid*3+2)) * (dyUx + dxUy)) +
                                               ((dzUy(0,node_lid*3+2) + dyUz(0,node_lid*3+2)) * (dyUz + dzUy)) +
                                               ((dxUz(0,node_lid*3+2) + dzUx(0,node_lid*3+2)) * (dzUx + dxUz)) ) * volume;


  FixedMatrix <1, 12> result_x = normal_strain_energy_x + compressibility_effect_x + shear_energy_x;
  FixedMatrix <1, 12> result_y = normal_strain_energy_y + compressibility_effect_y + shear_energy_y;
  FixedMatrix <1, 12> result_z = normal_strain_energy_z + compressibility_effect_z + shear_energy_z;

  FixedMatrix <3, 12> result;
  result(0,0)=result_x(0,0); result(0,1)=result_x(0,1); result(0,2) =result_x(0,2);  result(0,3)=result_x(0,3);
  result(0,4)=result_x(0,4); result(0,5)=result_x(0,5); result(0,6) =result_x(0,6);  result(0,7)=result_x(0,7);
  result(0,8)=result_x(0,8); result(0,9)=result_x(0,9); result(0,10)=result_x(0,10); result(0,11)=result_x(0,11);

  result(1,0)=result_y(0,0); result(1,1)=result_y(0,1); result(1,2) =result_y(0,2);  result(1,3) =result_y(0,3);
  result(1,4)=result_y(0,4); result(1,5)=result_y(0,5); result(1,6) =result_y(0,6);  result(1,7) =result_y(0,7);
  result(1,8)=result_y(0,8); result(1,9)=result_y(0,9); result(1,10)=result_y(0,10); result(1,11)=result_y(0,11);

  result(2,0)=result_z(0,0); result(2,1)=result_z(0,1); result(2,2) =result_z(0,2);  result(2,3) =result_z(0,3);
  result(2,4)=result_z(0,4); result(2,5)=result_z(0,5); result(2,6) =result_z(0,6);  result(2,7) =result_z(0,7);
  result(2,8)=result_z(0,8); result(2,9)=result_z(0,9); result(2,10)=result_z(0,10); result(2,11)=result_z(0,11);

  return result;
}