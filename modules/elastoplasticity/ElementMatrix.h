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
 *   σ(𝐮) is stress tensor       with  σᵢⱼ = C_tangᵢⱼₖₗεₖₗ
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
 *
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<6, 6> computeElementMatrixTria3Base(Real3 dxu, Real3 dyu, Real area, RealMatrix<3, 3> C_tang)
{
  RealVector<6> epsxx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  RealVector<6> epsyy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };
  RealVector<6> epsxy = { dyu[0], dxu[0], dyu[1], dxu[1], dyu[2], dxu[2] };

  // ∫∫ C_tang11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + C_tang12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + C_tang13 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<6, 6> sigmaXepsxx = (C_tang(0, 0) * epsxx + C_tang(0, 1) * epsyy + C_tang(0, 2) * epsxy) ^ epsxx;

  // ∫∫ C_tang12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + C_tang22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + C_tang23 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<6, 6> sigmaXepsyy = (C_tang(0, 1) * epsxx + C_tang(1, 1) * epsyy + C_tang(1, 2) * epsxy) ^ epsyy;

  // ∫∫ C_tang13 ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C_tang23 ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C_tang33 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
  RealMatrix<6, 6> sigmaXepsxy = (C_tang(0, 2) * epsxx + C_tang(1, 2) * epsyy + C_tang(2, 2) * epsxy) ^ epsxy;

  return area * (sigmaXepsxx + sigmaXepsyy + sigmaXepsxy);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<6, 6> computeHookeElementMatrixTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, RealVector<2> params)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  Real lambda = params[0];
  Real mu = params[1];
  RealMatrix<3, 3> C_tang{0.};

  C_tang(0, 0) = lambda + 2. * mu;
  C_tang(0, 1) = lambda;
  C_tang(0, 2) = 0.;
  C_tang(1, 0) = lambda;
  C_tang(1, 1) = lambda + 2. * mu;
  C_tang(1, 2) = 0.;
  C_tang(2, 0) = 0.;
  C_tang(2, 1) = 0.;
  C_tang(2, 2) = mu;

  return computeElementMatrixTria3Base(dxu, dyu, area, C_tang);
}

ARCCORE_HOST_DEVICE RealMatrix<6, 6> computeElementMatrixTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1, -1>>& in_C_tang)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  RealMatrix<3, 3> C_tang;

  // a flattened sequence is faster than a for loop on GPUs
  C_tang(0, 0) = in_C_tang(cell_lid, 0, 0);
  C_tang(0, 1) = in_C_tang(cell_lid, 0, 1);
  C_tang(0, 2) = in_C_tang(cell_lid, 0, 2);
  C_tang(1, 0) = in_C_tang(cell_lid, 1, 0);
  C_tang(1, 1) = in_C_tang(cell_lid, 1, 1);
  C_tang(1, 2) = in_C_tang(cell_lid, 1, 2);
  C_tang(2, 0) = in_C_tang(cell_lid, 2, 0);
  C_tang(2, 1) = in_C_tang(cell_lid, 2, 1);
  C_tang(2, 2) = in_C_tang(cell_lid, 2, 2);

  return computeElementMatrixTria3Base(dxu, dyu, area, C_tang);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

RealMatrix<6, 6> FemModuleElastoplasticity::_computeElementMatrixTria3(Cell cell)
{
  Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
  if (m_gp_material_tensor_strategy == "local") {
    return computeElementMatrixTria3Base(dxu, dyu, area, m_C_tang_2d);
  } else {
    RealMatrix<3, 3> C_tang_2d;
    for (Int32 ix = 0; ix < 3; ++ix) {
      for (Int32 iy = 0; iy < 3; ++iy) {
        C_tang_2d(ix, iy) = m_C_tang_2d_cell(cell, ix, iy);
      }
    }
    return computeElementMatrixTria3Base(dxu, dyu, area, C_tang_2d);
  }

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<2, 6> computeElementVectorTria3GpuBase(Real3 dxu, Real3 dyu, Real area, RealMatrix<3, 3> C_tang, Int32 node_lid)
{

  RealVector<6> epsxx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  RealVector<6> epsyy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };
  RealVector<6> epsxy = { dyu[0], dxu[0], dyu[1], dxu[1], dyu[2], dxu[2] };

  // ∫∫ C_tang11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + C_tang12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + C_tang13 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealVector<6> sigmaXepsxx_x = (C_tang(0, 0) * epsxx(node_lid*2) + C_tang(0, 1) * epsyy(node_lid*2) + C_tang(0, 2) * epsxy(node_lid*2)) * epsxx;
  RealVector<6> sigmaXepsxx_y = (C_tang(0, 0) * epsxx(node_lid*2+1) + C_tang(0, 1) * epsyy(node_lid*2+1) + C_tang(0, 2) * epsxy(node_lid*2+1)) * epsxx;

  // ∫∫ C_tang12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + C_tang22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + C_tang23 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealVector<6> sigmaXepsyy_x = (C_tang(0, 1) * epsxx(node_lid*2) + C_tang(1, 1) * epsyy(node_lid*2) + C_tang(1, 2) * epsxy(node_lid*2)) * epsyy;
  RealVector<6> sigmaXepsyy_y = (C_tang(0, 1) * epsxx(node_lid*2+1) + C_tang(1, 1) * epsyy(node_lid*2+1) + C_tang(1, 2) * epsxy(node_lid*2+1)) * epsyy;

  // ∫∫ C_tang13 ∂𝑢𝑥/∂𝑥 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C_tang23 ∂𝑢𝑦/∂𝑦 (∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥) + C_tang33 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)(∂𝑣𝑥/∂𝑦 + ∂𝑣𝑦/∂𝑥)
  RealVector<6> sigmaXepsxy_x = (C_tang(0, 2) * epsxx(node_lid*2) + C_tang(1, 2) * epsyy(node_lid*2) + C_tang(2, 2) * epsxy(node_lid*2)) * epsxy;
  RealVector<6> sigmaXepsxy_y = (C_tang(0, 2) * epsxx(node_lid*2+1) + C_tang(1, 2) * epsyy(node_lid*2+1) + C_tang(2, 2) * epsxy(node_lid*2+1)) * epsxy;

  RealVector <6> result_x = (sigmaXepsxx_x + sigmaXepsyy_x + sigmaXepsxy_x) * area;
  RealVector <6> result_y = (sigmaXepsxx_y + sigmaXepsyy_y + sigmaXepsxy_y) * area;

  RealMatrix<2, 6> result = {
    { result_x(0), result_x(1), result_x(2), result_x(3), result_x(4), result_x(5) },
    { result_y(0), result_y(1), result_y(2), result_y(3), result_y(4), result_y(5) }
  };

  return result;
}

ARCCORE_HOST_DEVICE RealMatrix<2, 6> computeHookeElementVectorTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, RealVector<2> params, Int32 node_lid)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  Real lambda = params[0];
  Real mu = params[1];
  RealMatrix<3, 3> C_tang;
  C_tang(0, 0) = lambda + 2. * mu;
  C_tang(0, 1) = lambda;
  C_tang(0, 2) = 0.;
  C_tang(1, 0) = lambda;
  C_tang(1, 1) = lambda + 2. * mu;
  C_tang(1, 2) = 0.;
  C_tang(2, 0) = 0.;
  C_tang(2, 1) = 0.;
  C_tang(2, 2) = mu;


  return computeElementVectorTria3GpuBase(dxu, dyu, area, C_tang, node_lid);
}


ARCCORE_HOST_DEVICE RealMatrix<2, 6> computeElementVectorTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1, -1>>& in_C_tang, Int32 node_lid)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  RealMatrix<3, 3> C_tang;

  // a flattened sequence is faster than a for loop on GPUs
  C_tang(0, 0) = in_C_tang(cell_lid, 0, 0);
  C_tang(0, 1) = in_C_tang(cell_lid, 0, 1);
  C_tang(0, 2) = in_C_tang(cell_lid, 0, 2);
  C_tang(1, 0) = in_C_tang(cell_lid, 1, 0);
  C_tang(1, 1) = in_C_tang(cell_lid, 1, 1);
  C_tang(1, 2) = in_C_tang(cell_lid, 1, 2);
  C_tang(2, 0) = in_C_tang(cell_lid, 2, 0);
  C_tang(2, 1) = in_C_tang(cell_lid, 2, 1);
  C_tang(2, 2) = in_C_tang(cell_lid, 2, 2);

  return computeElementVectorTria3GpuBase(dxu, dyu, area, C_tang, node_lid);
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
 *   σ(𝐮) is stress tensor       with  σᵢⱼ = C_tangᵢⱼₖₗεₖₗ
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

ARCCORE_HOST_DEVICE RealMatrix<12, 12> computeElementMatrixTetra4Base(Real4 dxu, Real4 dyu, Real4 dzu, Real volume, RealMatrix<6, 6> C_tang)
{
  RealVector<12> epsxx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
  RealVector<12> epsyy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
  RealVector<12> epszz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

  RealVector<12> epsyz = { 0., dzu[0], dyu[0],    0., dzu[1], dyu[1],    0., dzu[2], dyu[2],    0., dzu[3], dyu[3] };
  RealVector<12> epszx = { dzu[0], 0., dxu[0],    dzu[1], 0., dxu[1],    dzu[2], 0., dxu[2],    dzu[3], 0., dxu[3] };
  RealVector<12> epsxy = { dyu[0], dxu[0], 0.,    dyu[1], dxu[1], 0.,    dyu[2], dxu[2], 0.,    dyu[3], dxu[3], 0. };

  // ∫∫∫ C_tang11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + C_tang12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + C_tang13 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + C_tang14 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑥/∂𝑥 + C_tang15 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑥/∂𝑥 + C_tang16 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealMatrix<12, 12> sigmaXepsxx = (C_tang(0, 0) * epsxx + C_tang(0, 1) * epsyy + C_tang(0, 2) * epszz + C_tang(0, 3) * epsyz + C_tang(0, 4) * epszx + C_tang(0, 5) * epsxy) ^ epsxx;
  // ∫∫∫ C_tang12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + C_tang22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + C_tang23 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 + C_tang24 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑦/∂𝑦 + C_tang25 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑦/∂𝑦 + C_tang26 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑦/∂𝑦
  RealMatrix<12, 12> sigmaXepsyy = (C_tang(0, 1) * epsxx + C_tang(1, 1) * epsyy + C_tang(1, 2) * epszz + C_tang(1, 3) * epsyz + C_tang(1, 4) * epszx + C_tang(1, 5) * epsxy) ^ epsyy;
  // ∫∫∫ C_tang13 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 + C_tang23 ∂𝑢𝑧/∂𝑦 ∂𝑣𝑧/∂𝑧 + C_tang33 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧 + C_tang34 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑧/∂𝑧 + C_tang35 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑧/∂𝑧 + C_tang36 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑧/∂𝑧
  RealMatrix<12, 12> sigmaXepszz = (C_tang(0, 2) * epsxx + C_tang(1, 2) * epsyy + C_tang(2, 2) * epszz + C_tang(2, 3) * epsyz + C_tang(2, 4) * epszx + C_tang(2, 5) * epsxy) ^ epszz;
  // ∫∫∫ C_tang14 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang24 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang34 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang44 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang45 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang46 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)
  RealMatrix<12, 12> sigmaXepsyz = (C_tang(0, 3) * epsxx + C_tang(1, 3) * epsyy + C_tang(2, 3) * epszz + C_tang(3, 3) * epsyz + C_tang(3, 4) * epszx + C_tang(3, 5) * epsxy) ^ epsyz;
  // ∫∫∫ C_tang15 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang25 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang35 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang45 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang55 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang56 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)
  RealMatrix<12, 12> sigmaXepszx = (C_tang(0, 4) * epsxx + C_tang(1, 4) * epsyy + C_tang(2, 4) * epszz + C_tang(3, 4) * epsyz + C_tang(4, 4) * epszx + C_tang(4, 5) * epsxy) ^ epszx;
  // ∫∫∫ C_tang16 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang26 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang36 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang46 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang56 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang66 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)
  RealMatrix<12, 12> sigmaXepsxy = (C_tang(0, 5) * epsxx + C_tang(1, 5) * epsyy + C_tang(2, 5) * epszz + C_tang(3, 5) * epsyz + C_tang(4, 5) * epszx + C_tang(5, 5) * epsxy) ^ epsxy;

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

  if (m_gp_material_tensor_strategy == "local") {
    return computeElementMatrixTetra4Base(dxu, dyu, dzu, volume, m_C_tang_3d);
  } else {
    RealMatrix<6, 6> C_tang_3d;
    for (Int32 ix = 0; ix < 6; ++ix) {
      for (Int32 iy = 0; iy < 6; ++iy) {
        C_tang_3d(ix, iy) = m_C_tang_3d_cell(cell, ix, iy);
      }
    }
    return computeElementMatrixTetra4Base(dxu, dyu, dzu, volume, C_tang_3d);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<12, 12> computeElementMatrixTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, RealMatrix<6, 6> C_tang)
{
  Real4 dxu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);

  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  return computeElementMatrixTetra4Base(dxu, dyu, dzu, volume, C_tang);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<3, 12> computeElementVectorTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, RealMatrix<6, 6> C_tang, Int32 node_lid)
{
  Real4 dxu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzu = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);

  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  Int32 idx_x = node_lid * 3;
  Int32 idx_y = idx_x + 1;
  Int32 idx_z = idx_x + 2;

  RealVector<12> epsxx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
  RealVector<12> epsyy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
  RealVector<12> epszz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

  RealVector<12> epsyz = { 0., dzu[0], dyu[0],    0., dzu[1], dyu[1],    0., dzu[2], dyu[2],    0., dzu[3], dyu[3] };
  RealVector<12> epszx = { dzu[0], 0., dxu[0],    dzu[1], 0., dxu[1],    dzu[2], 0., dxu[2],    dzu[3], 0., dxu[3] };
  RealVector<12> epsxy = { dyu[0], dxu[0], 0.,    dyu[1], dxu[1], 0.,    dyu[2], dxu[2], 0.,    dyu[3], dxu[3], 0. };

  // ∫∫∫ C_tang11 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑥/∂𝑥 + C_tang12 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑥/∂𝑥 + C_tang13 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑥/∂𝑥 + C_tang14 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑥/∂𝑥 + C_tang15 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑥/∂𝑥 + C_tang16 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑥/∂𝑥
  RealVector<12> sigmaXepsxx_x = (C_tang(0, 0) * epsxx(idx_x) + C_tang(0, 1) * epsyy(idx_x) + C_tang(0, 2) * epszz(idx_x) + C_tang(0, 3) * epsyz(idx_x) + C_tang(0, 4) * epszx(idx_x) + C_tang(0, 5) * epsxy(idx_x)) * epsxx;
  RealVector<12> sigmaXepsxx_y = (C_tang(0, 0) * epsxx(idx_y) + C_tang(0, 1) * epsyy(idx_y) + C_tang(0, 2) * epszz(idx_y) + C_tang(0, 3) * epsyz(idx_y) + C_tang(0, 4) * epszx(idx_y) + C_tang(0, 5) * epsxy(idx_y)) * epsxx;
  RealVector<12> sigmaXepsxx_z = (C_tang(0, 0) * epsxx(idx_z) + C_tang(0, 1) * epsyy(idx_z) + C_tang(0, 2) * epszz(idx_z) + C_tang(0, 3) * epsyz(idx_z) + C_tang(0, 4) * epszx(idx_z) + C_tang(0, 5) * epsxy(idx_z)) * epsxx;

  // ∫∫∫ C_tang12 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑦 + C_tang22 ∂𝑢𝑦/∂𝑦 ∂𝑣𝑦/∂𝑦 + C_tang23 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑦/∂𝑦 + C_tang24 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑦/∂𝑦 + C_tang25 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑦/∂𝑦 + C_tang26 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑦/∂𝑦
  RealVector<12> sigmaXepsyy_x = (C_tang(0, 1) * epsxx(idx_x) + C_tang(1, 1) * epsyy(idx_x) + C_tang(1, 2) * epszz(idx_x) + C_tang(1, 3) * epsyz(idx_x) + C_tang(1, 4) * epszx(idx_x) + C_tang(1, 5) * epsxy(idx_x)) * epsyy;
  RealVector<12> sigmaXepsyy_y = (C_tang(0, 1) * epsxx(idx_y) + C_tang(1, 1) * epsyy(idx_y) + C_tang(1, 2) * epszz(idx_y) + C_tang(1, 3) * epsyz(idx_y) + C_tang(1, 4) * epszx(idx_y) + C_tang(1, 5) * epsxy(idx_y)) * epsyy;
  RealVector<12> sigmaXepsyy_z = (C_tang(0, 1) * epsxx(idx_z) + C_tang(1, 1) * epsyy(idx_z) + C_tang(1, 2) * epszz(idx_z) + C_tang(1, 3) * epsyz(idx_z) + C_tang(1, 4) * epszx(idx_z) + C_tang(1, 5) * epsxy(idx_z)) * epsyy;

  // ∫∫∫ C_tang13 ∂𝑢𝑥/∂𝑥 ∂𝑣𝑧/∂𝑧 + C_tang23 ∂𝑢𝑧/∂𝑦 ∂𝑣𝑧/∂𝑧 + C_tang33 ∂𝑢𝑧/∂𝑧 ∂𝑣𝑧/∂𝑧 + C_tang34 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) ∂𝑣𝑧/∂𝑧 + C_tang35 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) ∂𝑣𝑧/∂𝑧 + C_tang36 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) ∂𝑣𝑧/∂𝑧
  RealVector<12> sigmaXepszz_x = (C_tang(0, 2) * epsxx(idx_x) + C_tang(1, 2) * epsyy(idx_x) + C_tang(2, 2) * epszz(idx_x) + C_tang(2, 3) * epsyz(idx_x) + C_tang(2, 4) * epszx(idx_x) + C_tang(2, 5) * epsxy(idx_x)) * epszz;
  RealVector<12> sigmaXepszz_y = (C_tang(0, 2) * epsxx(idx_y) + C_tang(1, 2) * epsyy(idx_y) + C_tang(2, 2) * epszz(idx_y) + C_tang(2, 3) * epsyz(idx_y) + C_tang(2, 4) * epszx(idx_y) + C_tang(2, 5) * epsxy(idx_y)) * epszz;
  RealVector<12> sigmaXepszz_z = (C_tang(0, 2) * epsxx(idx_z) + C_tang(1, 2) * epsyy(idx_z) + C_tang(2, 2) * epszz(idx_z) + C_tang(2, 3) * epsyz(idx_z) + C_tang(2, 4) * epszx(idx_z) + C_tang(2, 5) * epsxy(idx_z)) * epszz;

  // ∫∫∫ C_tang14 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang24 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang34 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang44 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang45 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) + C_tang46 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧)
  RealVector<12> sigmaXepsyz_x = (C_tang(0, 3) * epsxx(idx_x) + C_tang(1, 3) * epsyy(idx_x) + C_tang(2, 3) * epszz(idx_x) + C_tang(3, 3) * epsyz(idx_x) + C_tang(3, 4) * epszx(idx_x) + C_tang(3, 5) * epsxy(idx_x)) * epsyz;
  RealVector<12> sigmaXepsyz_y = (C_tang(0, 3) * epsxx(idx_y) + C_tang(1, 3) * epsyy(idx_y) + C_tang(2, 3) * epszz(idx_y) + C_tang(3, 3) * epsyz(idx_y) + C_tang(3, 4) * epszx(idx_y) + C_tang(3, 5) * epsxy(idx_y)) * epsyz;
  RealVector<12> sigmaXepsyz_z = (C_tang(0, 3) * epsxx(idx_z) + C_tang(1, 3) * epsyy(idx_z) + C_tang(2, 3) * epszz(idx_z) + C_tang(3, 3) * epsyz(idx_z) + C_tang(3, 4) * epszx(idx_z) + C_tang(3, 5) * epsxy(idx_z)) * epsyz;

  // ∫∫∫ C_tang15 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang25 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang35 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang45 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang55 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) + C_tang56 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥)
  RealVector<12> sigmaXepszx_x = (C_tang(0, 4) * epsxx(idx_x) + C_tang(1, 4) * epsyy(idx_x) + C_tang(2, 4) * epszz(idx_x) + C_tang(3, 4) * epsyz(idx_x) + C_tang(4, 4) * epszx(idx_x) + C_tang(4, 5) * epsxy(idx_x)) * epszx;
  RealVector<12> sigmaXepszx_y = (C_tang(0, 4) * epsxx(idx_y) + C_tang(1, 4) * epsyy(idx_y) + C_tang(2, 4) * epszz(idx_y) + C_tang(3, 4) * epsyz(idx_y) + C_tang(4, 4) * epszx(idx_y) + C_tang(4, 5) * epsxy(idx_y)) * epszx;
  RealVector<12> sigmaXepszx_z = (C_tang(0, 4) * epsxx(idx_z) + C_tang(1, 4) * epsyy(idx_z) + C_tang(2, 4) * epszz(idx_z) + C_tang(3, 4) * epsyz(idx_z) + C_tang(4, 4) * epszx(idx_z) + C_tang(4, 5) * epsxy(idx_z)) * epszx;

  // ∫∫∫ C_tang16 ∂𝑢𝑥/∂𝑥 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang26 ∂𝑢𝑧/∂𝑦 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang36 ∂𝑢𝑧/∂𝑧 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang46 (∂𝑢𝑧/∂𝑦 + ∂𝑢𝑦/∂𝑧) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang56 (∂𝑢𝑥/∂𝑧 + ∂𝑢𝑧/∂𝑥) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) + C_tang66 (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦) (∂𝑢𝑦/∂𝑥 + ∂𝑢𝑥/∂𝑦)
  RealVector<12> sigmaXepsxy_x = (C_tang(0, 5) * epsxx(idx_x) + C_tang(1, 5) * epsyy(idx_x) + C_tang(2, 5) * epszz(idx_x) + C_tang(3, 5) * epsyz(idx_x) + C_tang(4, 5) * epszx(idx_x) + C_tang(5, 5) * epsxy(idx_x)) * epsxy;
  RealVector<12> sigmaXepsxy_y = (C_tang(0, 5) * epsxx(idx_y) + C_tang(1, 5) * epsyy(idx_y) + C_tang(2, 5) * epszz(idx_y) + C_tang(3, 5) * epsyz(idx_y) + C_tang(4, 5) * epszx(idx_y) + C_tang(5, 5) * epsxy(idx_y)) * epsxy;
  RealVector<12> sigmaXepsxy_z = (C_tang(0, 5) * epsxx(idx_z) + C_tang(1, 5) * epsyy(idx_z) + C_tang(2, 5) * epszz(idx_z) + C_tang(3, 5) * epsyz(idx_z) + C_tang(4, 5) * epszx(idx_z) + C_tang(5, 5) * epsxy(idx_z)) * epsxy;

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