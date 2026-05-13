// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrix.h                                             (C) 2000-2026 */
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
 *   a(𝐮,𝐯) = ∫∫ σ(𝐮):ε(𝐯)dΩ     with  𝐮 = (𝑢𝑥,𝑢𝑦) and 𝐯 = (𝑣𝑥,𝑣𝑦)
 *   σ(𝐮) is stress tensor       with  σᵢⱼ = Mtᵢⱼₖₗεₖₗ
 *   ε(𝐯) is strain tensor       with  εᵢⱼ = 0.5 (∂𝑣ᵢ/∂xⱼ + ∂𝑣ⱼ/∂xᵢ)
 *
 *   the bilinear integral expands to
 *
 *      a(𝐮,𝐯) = ∫∫ [σ_𝑥𝑥ε_𝑥𝑥 + σ_𝑦𝑦ε_𝑦𝑦 + 2σ_𝑥𝑦ε_𝑥𝑦]dΩ
 *
 *   this further expands to
 *
 *      a(𝐮,𝐯) =   ∫∫ Mt11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + Mt12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + Mt13 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ Mt12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + Mt22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + Mt23 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫ Mt13 ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + Mt23 ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + Mt33 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
 *
 *
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<6, 6> computeElementMatrixTria3Base(Real3 dxu, Real3 dyu, Real area, Real lambda, Real mu)
{

  // RealVector<6> Mt_vec = {lambda + 2. * mu, lambda + 2. * mu, mu, lambda, 0.};

  RealMatrix<3, 3> Mt = {
    {lambda + 2. * mu, lambda,           0.},
          {lambda,           lambda + 2. * mu, 0.},
          {0.,               0.,               mu}
  };

  RealVector<6> epsxx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  RealVector<6> epsyy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };
  RealVector<6> epsxy = { dyu[0], dxu[0], dyu[1], dxu[1], dyu[2], dxu[2] };

  // ∫∫ Mt11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + Mt12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + Mt13 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<6, 6> sigmaXepsxx = (Mt(0, 0) * epsxx + Mt(0, 1) * epsyy + Mt(0, 2) * epsxy) ^ epsxx;

  // ∫∫ Mt12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + Mt22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + Mt23 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<6, 6> sigmaXepsyy = (Mt(0, 1) * epsxx + Mt(1, 1) * epsyy + Mt(1, 2) * epsxy) ^ epsyy;

  // ∫∫ Mt13 ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + Mt23 ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + Mt33 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
  RealMatrix<6, 6> sigmaXepsxy = (Mt(0, 2) * epsxx + Mt(1, 2) * epsyy + Mt(2, 2) * epsxy) ^ epsxy;

  return area * (sigmaXepsxx + sigmaXepsyy + sigmaXepsxy);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<6, 6> computeElementMatrixTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  return computeElementMatrixTria3Base(dxu, dyu, area, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

RealMatrix<6, 6> FemModuleElastoplasticity::_computeElementMatrixTria3(Cell cell)
{
  Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  return computeElementMatrixTria3Base(dxu, dyu, area, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<2, 6> computeElementVectorTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu, Int32 node_lid)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  RealMatrix<3, 3> Mt = {
    {lambda + 2. * mu, lambda,           0.},
          {lambda,           lambda + 2. * mu, 0.},
          {0.,               0.,               mu}
  };

  RealVector<6> epsxx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  RealVector<6> epsyy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };
  RealVector<6> epsxy = { dyu[0], dxu[0], dyu[1], dxu[1], dyu[2], dxu[2] };

  // ∫∫ Mt11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + Mt12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + Mt13 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealVector<6> sigmaXepsxx_x = (Mt(0, 0) * epsxx(node_lid*2) + Mt(0, 1) * epsyy(node_lid*2) + Mt(0, 2) * epsxy(node_lid*2)) * epsxx;
  RealVector<6> sigmaXepsxx_y = (Mt(0, 0) * epsxx(node_lid*2+1) + Mt(0, 1) * epsyy(node_lid*2+1) + Mt(0, 2) * epsxy(node_lid*2+1)) * epsxx;

  // ∫∫ Mt12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + Mt22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + Mt23 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealVector<6> sigmaXepsyy_x = (Mt(0, 1) * epsxx(node_lid*2) + Mt(1, 1) * epsyy(node_lid*2) + Mt(1, 2) * epsxy(node_lid*2)) * epsyy;
  RealVector<6> sigmaXepsyy_y = (Mt(0, 1) * epsxx(node_lid*2+1) + Mt(1, 1) * epsyy(node_lid*2+1) + Mt(1, 2) * epsxy(node_lid*2+1)) * epsyy;

  // ∫∫ Mt13 ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + Mt23 ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + Mt33 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
  RealVector<6> sigmaXepsxy_x = (Mt(0, 2) * epsxx(node_lid*2) + Mt(1, 2) * epsyy(node_lid*2) + Mt(2, 2) * epsxy(node_lid*2)) * epsxy;
  RealVector<6> sigmaXepsxy_y = (Mt(0, 2) * epsxx(node_lid*2+1) + Mt(1, 2) * epsyy(node_lid*2+1) + Mt(2, 2) * epsxy(node_lid*2+1)) * epsxy;

  RealVector <6> result_x = (sigmaXepsxx_x + sigmaXepsyy_x + sigmaXepsxy_x) * area;
  RealVector <6> result_y = (sigmaXepsxx_y + sigmaXepsyy_y + sigmaXepsxy_y) * area;

  RealMatrix<2, 6> result = {
    { result_x(0), result_x(1), result_x(2), result_x(3), result_x(4), result_x(5) },
    { result_y(0), result_y(1), result_y(2), result_y(3), result_y(4), result_y(5) }
  };

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
 *   σ(𝐮) is stress tensor       with  σᵢⱼ = Mtᵢⱼₖₗεₖₗ
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
 *      a(𝐮,𝐯) =   ∫∫∫ Mt11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + Mt12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + Mt13 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + Mt14 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑥/∂𝑥 + Mt15 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑥/∂𝑥 + Mt16 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
 *               + ∫∫∫ Mt12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + Mt22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + Mt23 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 + Mt24 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑦/∂𝑦 + Mt25 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑦/∂𝑦 + Mt26 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑦/∂𝑦
 *               + ∫∫∫ Mt13 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 + Mt23 ∂𝑢𝑧/∂𝑦 ∂𝑣𝑧/∂𝑧 + Mt33 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧 + Mt34 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑧/∂𝑧 + Mt35 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑧/∂𝑧 + Mt36 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑧/∂𝑧
 *               + ∫∫∫ Mt14 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt24 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt34 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt44 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt45 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt46 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)
 *               + ∫∫∫ Mt15 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt25 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt35 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt45 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt55 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt56 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)
 *               + ∫∫∫ Mt16 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt26 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt36 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt46 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt56 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt66 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)
 *
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<12, 12> computeElementMatrixTetra4Base(Real4 dxu, Real4 dyu, Real4 dzu, Real volume, Real lambda, Real mu)
{
  RealMatrix<6, 6> Mt = {
    {lambda + 2.*mu, lambda,         lambda,         0.,  0.,  0.},
          {lambda,         lambda + 2.*mu, lambda,         0.,  0.,  0.},
          {lambda,         lambda,         lambda + 2.*mu, 0.,  0.,  0.},
          {0.,             0.,             0.,             mu,  0.,  0.},
          {0.,             0.,             0.,             0.,  mu,  0.},
          {0.,             0.,             0.,             0.,  0.,  mu}
  };

  RealVector<12> epsxx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
  RealVector<12> epsyy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
  RealVector<12> epszz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

  RealVector<12> epsyz = { 0., dzu[0], dyu[0],    0., dzu[1], dyu[1],    0., dzu[2], dyu[2],    0., dzu[3], dyu[3] };
  RealVector<12> epszx = { dzu[0], 0., dxu[0],    dzu[1], 0., dxu[1],    dzu[2], 0., dxu[2],    dzu[3], 0., dxu[3] };
  RealVector<12> epsxy = { dyu[0], dxu[0], 0.,    dyu[1], dxu[1], 0.,    dyu[2], dxu[2], 0.,    dyu[3], dxu[3], 0. };

  // ∫∫∫ Mt11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + Mt12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + Mt13 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + Mt14 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑥/∂𝑥 + Mt15 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑥/∂𝑥 + Mt16 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<12, 12> sigmaXepsxx = (Mt(0, 0) * epsxx + Mt(0, 1) * epsyy + Mt(0, 2) * epszz + Mt(0, 3) * epsyz + Mt(0, 4) * epszx + Mt(0, 5) * epsxy) ^ epsxx;
  // ∫∫∫ Mt12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + Mt22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + Mt23 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 + Mt24 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑦/∂𝑦 + Mt25 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑦/∂𝑦 + Mt26 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑦/∂𝑦
  RealMatrix<12, 12> sigmaXepsyy = (Mt(0, 1) * epsxx + Mt(1, 1) * epsyy + Mt(1, 2) * epszz + Mt(1, 3) * epsyz + Mt(1, 4) * epszx + Mt(1, 5) * epsxy) ^ epsyy;
  // ∫∫∫ Mt13 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 + Mt23 ∂𝑢𝑧/∂𝑦 ∂𝑣𝑧/∂𝑧 + Mt33 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧 + Mt34 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑧/∂𝑧 + Mt35 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑧/∂𝑧 + Mt36 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑧/∂𝑧
  RealMatrix<12, 12> sigmaXepszz = (Mt(0, 2) * epsxx + Mt(1, 2) * epsyy + Mt(2, 2) * epszz + Mt(2, 3) * epsyz + Mt(2, 4) * epszx + Mt(2, 5) * epsxy) ^ epszz;
  // ∫∫∫ Mt14 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt24 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt34 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt44 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt45 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt46 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)
  RealMatrix<12, 12> sigmaXepsyz = (Mt(0, 3) * epsxx + Mt(1, 3) * epsyy + Mt(2, 3) * epszz + Mt(3, 3) * epsyz + Mt(3, 4) * epszx + Mt(3, 5) * epsxy) ^ epsyz;
  // ∫∫∫ Mt15 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt25 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt35 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt45 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt55 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt56 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)
  RealMatrix<12, 12> sigmaXepszx = (Mt(0, 4) * epsxx + Mt(1, 4) * epsyy + Mt(2, 4) * epszz + Mt(3, 4) * epsyz + Mt(4, 4) * epszx + Mt(4, 5) * epsxy) ^ epszx;
  // ∫∫∫ Mt16 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt26 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt36 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt46 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt56 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt66 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)
  RealMatrix<12, 12> sigmaXepsxy = (Mt(0, 5) * epsxx + Mt(1, 5) * epsyy + Mt(2, 5) * epszz + Mt(3, 5) * epsyz + Mt(4, 5) * epszx + Mt(5, 5) * epsxy) ^ epsxy;

  return volume * ( sigmaXepsxx + sigmaXepsyy + sigmaXepszz + sigmaXepsyz + sigmaXepszx + sigmaXepsxy);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

RealMatrix<12, 12> FemModuleElastoplasticity::_computeElementMatrixTetra4(Cell cell)
{
  Real4 dxu = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
  Real4 dyu = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
  Real4 dzu = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);

  Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);

  return computeElementMatrixTetra4Base(dxu, dyu, dzu, volume, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<12, 12> computeElementMatrixTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu)
{
  Real4 dxu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);

  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  return computeElementMatrixTetra4Base(dxu, dyu, dzu, volume, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<3, 12> computeElementVectorTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu, Int32 node_lid)
{
  Real4 dxu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);

  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  Int32 idx_x = node_lid * 3;
  Int32 idx_y = idx_x + 1;
  Int32 idx_z = idx_x + 2;

  RealMatrix<6, 6> Mt = {
    {lambda + 2.*mu, lambda,         lambda,         0.,  0.,  0.},
          {lambda,         lambda + 2.*mu, lambda,         0.,  0.,  0.},
          {lambda,         lambda,         lambda + 2.*mu, 0.,  0.,  0.},
          {0.,             0.,             0.,             mu,  0.,  0.},
          {0.,             0.,             0.,             0.,  mu,  0.},
          {0.,             0.,             0.,             0.,  0.,  mu}
  };

  RealVector<12> epsxx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
  RealVector<12> epsyy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
  RealVector<12> epszz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

  RealVector<12> epsyz = { 0., dzu[0], dyu[0],    0., dzu[1], dyu[1],    0., dzu[2], dyu[2],    0., dzu[3], dyu[3] };
  RealVector<12> epszx = { dzu[0], 0., dxu[0],    dzu[1], 0., dxu[1],    dzu[2], 0., dxu[2],    dzu[3], 0., dxu[3] };
  RealVector<12> epsxy = { dyu[0], dxu[0], 0.,    dyu[1], dxu[1], 0.,    dyu[2], dxu[2], 0.,    dyu[3], dxu[3], 0. };

  // ∫∫∫ Mt11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + Mt12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + Mt13 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + Mt14 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑥/∂𝑥 + Mt15 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑥/∂𝑥 + Mt16 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealVector<12> sigmaXepsxx_x = (Mt(0, 0) * epsxx(idx_x) + Mt(0, 1) * epsyy(idx_x) + Mt(0, 2) * epszz(idx_x) + Mt(0, 3) * epsyz(idx_x) + Mt(0, 4) * epszx(idx_x) + Mt(0, 5) * epsxy(idx_x)) * epsxx;
  RealVector<12> sigmaXepsxx_y = (Mt(0, 0) * epsxx(idx_y) + Mt(0, 1) * epsyy(idx_y) + Mt(0, 2) * epszz(idx_y) + Mt(0, 3) * epsyz(idx_y) + Mt(0, 4) * epszx(idx_y) + Mt(0, 5) * epsxy(idx_y)) * epsxx;
  RealVector<12> sigmaXepsxx_z = (Mt(0, 0) * epsxx(idx_z) + Mt(0, 1) * epsyy(idx_z) + Mt(0, 2) * epszz(idx_z) + Mt(0, 3) * epsyz(idx_z) + Mt(0, 4) * epszx(idx_z) + Mt(0, 5) * epsxy(idx_z)) * epsxx;

  // ∫∫∫ Mt12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + Mt22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + Mt23 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 + Mt24 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑦/∂𝑦 + Mt25 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑦/∂𝑦 + Mt26 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑦/∂𝑦
  RealVector<12> sigmaXepsyy_x = (Mt(0, 1) * epsxx(idx_x) + Mt(1, 1) * epsyy(idx_x) + Mt(1, 2) * epszz(idx_x) + Mt(1, 3) * epsyz(idx_x) + Mt(1, 4) * epszx(idx_x) + Mt(1, 5) * epsxy(idx_x)) * epsyy;
  RealVector<12> sigmaXepsyy_y = (Mt(0, 1) * epsxx(idx_y) + Mt(1, 1) * epsyy(idx_y) + Mt(1, 2) * epszz(idx_y) + Mt(1, 3) * epsyz(idx_y) + Mt(1, 4) * epszx(idx_y) + Mt(1, 5) * epsxy(idx_y)) * epsyy;
  RealVector<12> sigmaXepsyy_z = (Mt(0, 1) * epsxx(idx_z) + Mt(1, 1) * epsyy(idx_z) + Mt(1, 2) * epszz(idx_z) + Mt(1, 3) * epsyz(idx_z) + Mt(1, 4) * epszx(idx_z) + Mt(1, 5) * epsxy(idx_z)) * epsyy;

  // ∫∫∫ Mt13 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 + Mt23 ∂𝑢𝑧/∂𝑦 ∂𝑣𝑧/∂𝑧 + Mt33 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧 + Mt34 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑧/∂𝑧 + Mt35 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑧/∂𝑧 + Mt36 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑧/∂𝑧
  RealVector<12> sigmaXepszz_x = (Mt(0, 2) * epsxx(idx_x) + Mt(1, 2) * epsyy(idx_x) + Mt(2, 2) * epszz(idx_x) + Mt(2, 3) * epsyz(idx_x) + Mt(2, 4) * epszx(idx_x) + Mt(2, 5) * epsxy(idx_x)) * epszz;
  RealVector<12> sigmaXepszz_y = (Mt(0, 2) * epsxx(idx_y) + Mt(1, 2) * epsyy(idx_y) + Mt(2, 2) * epszz(idx_y) + Mt(2, 3) * epsyz(idx_y) + Mt(2, 4) * epszx(idx_y) + Mt(2, 5) * epsxy(idx_y)) * epszz;
  RealVector<12> sigmaXepszz_z = (Mt(0, 2) * epsxx(idx_z) + Mt(1, 2) * epsyy(idx_z) + Mt(2, 2) * epszz(idx_z) + Mt(2, 3) * epsyz(idx_z) + Mt(2, 4) * epszx(idx_z) + Mt(2, 5) * epsxy(idx_z)) * epszz;

  // ∫∫∫ Mt14 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt24 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt34 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt44 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt45 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + Mt46 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)
  RealVector<12> sigmaXepsyz_x = (Mt(0, 3) * epsxx(idx_x) + Mt(1, 3) * epsyy(idx_x) + Mt(2, 3) * epszz(idx_x) + Mt(3, 3) * epsyz(idx_x) + Mt(3, 4) * epszx(idx_x) + Mt(3, 5) * epsxy(idx_x)) * epsyz;
  RealVector<12> sigmaXepsyz_y = (Mt(0, 3) * epsxx(idx_y) + Mt(1, 3) * epsyy(idx_y) + Mt(2, 3) * epszz(idx_y) + Mt(3, 3) * epsyz(idx_y) + Mt(3, 4) * epszx(idx_y) + Mt(3, 5) * epsxy(idx_y)) * epsyz;
  RealVector<12> sigmaXepsyz_z = (Mt(0, 3) * epsxx(idx_z) + Mt(1, 3) * epsyy(idx_z) + Mt(2, 3) * epszz(idx_z) + Mt(3, 3) * epsyz(idx_z) + Mt(3, 4) * epszx(idx_z) + Mt(3, 5) * epsxy(idx_z)) * epsyz;

  // ∫∫∫ Mt15 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt25 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt35 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt45 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt55 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + Mt56 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)
  RealVector<12> sigmaXepszx_x = (Mt(0, 4) * epsxx(idx_x) + Mt(1, 4) * epsyy(idx_x) + Mt(2, 4) * epszz(idx_x) + Mt(3, 4) * epsyz(idx_x) + Mt(4, 4) * epszx(idx_x) + Mt(4, 5) * epsxy(idx_x)) * epszx;
  RealVector<12> sigmaXepszx_y = (Mt(0, 4) * epsxx(idx_y) + Mt(1, 4) * epsyy(idx_y) + Mt(2, 4) * epszz(idx_y) + Mt(3, 4) * epsyz(idx_y) + Mt(4, 4) * epszx(idx_y) + Mt(4, 5) * epsxy(idx_y)) * epszx;
  RealVector<12> sigmaXepszx_z = (Mt(0, 4) * epsxx(idx_z) + Mt(1, 4) * epsyy(idx_z) + Mt(2, 4) * epszz(idx_z) + Mt(3, 4) * epsyz(idx_z) + Mt(4, 4) * epszx(idx_z) + Mt(4, 5) * epsxy(idx_z)) * epszx;

  // ∫∫∫ Mt16 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt26 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt36 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt46 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt56 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + Mt66 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)
  RealVector<12> sigmaXepsxy_x = (Mt(0, 5) * epsxx(idx_x) + Mt(1, 5) * epsyy(idx_x) + Mt(2, 5) * epszz(idx_x) + Mt(3, 5) * epsyz(idx_x) + Mt(4, 5) * epszx(idx_x) + Mt(5, 5) * epsxy(idx_x)) * epsxy;
  RealVector<12> sigmaXepsxy_y = (Mt(0, 5) * epsxx(idx_y) + Mt(1, 5) * epsyy(idx_y) + Mt(2, 5) * epszz(idx_y) + Mt(3, 5) * epsyz(idx_y) + Mt(4, 5) * epszx(idx_y) + Mt(5, 5) * epsxy(idx_y)) * epsxy;
  RealVector<12> sigmaXepsxy_z = (Mt(0, 5) * epsxx(idx_z) + Mt(1, 5) * epsyy(idx_z) + Mt(2, 5) * epszz(idx_z) + Mt(3, 5) * epsyz(idx_z) + Mt(4, 5) * epszx(idx_z) + Mt(5, 5) * epsxy(idx_z)) * epsxy;

  RealVector<12> result_x = volume * ( sigmaXepsxx_x + sigmaXepsyy_x + sigmaXepszz_x + sigmaXepsyz_x + sigmaXepszx_x + sigmaXepsxy_x );
  RealVector<12> result_y = volume * ( sigmaXepsxx_y + sigmaXepsyy_y + sigmaXepszz_y + sigmaXepsyz_y + sigmaXepszx_y + sigmaXepsxy_y );
  RealVector<12> result_z = volume * ( sigmaXepsxx_z + sigmaXepsyy_z + sigmaXepszz_z + sigmaXepsyz_z + sigmaXepszx_z + sigmaXepsxy_z );

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