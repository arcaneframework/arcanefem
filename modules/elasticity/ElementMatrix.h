// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrix.h                                             (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matric for Elasticity       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (P1 FE).
 *
 * Theory:
 *
 *   a(U,V) = ∫∫ [σ(U):ε(V)dΩ    with  U = (u𝑥,u𝑦) and V = (v𝑥,v𝑦)
 *   σ(U) is stress tensor       with  σᵢⱼ = λδᵢⱼεₖₖ + 2μεᵢⱼ
 *   ε(V) is strain tensor       with  εᵢⱼ = 0.5 (∂vᵢ/∂xⱼ + ∂vⱼ/∂xᵢ)
 *
 *   the bilinear integral expands to
 *
 *      a(U,V) = ∫∫ [σ_𝑥𝑥ε_𝑥𝑥 + σ_𝑦𝑦ε_𝑦𝑦 + 2σ_𝑥𝑦ε_𝑥𝑦]dΩ
 *
 *   this further expands to
 *
 *      a(U,V) =   ∫∫ (λ+2μ)(∂u𝑥/∂𝑥 ∂v𝑥/∂𝑥 + ∂u𝑦/∂𝑦 ∂v𝑦/∂𝑦)
 *               + ∫∫ λ(∂u𝑦/∂𝑦 ∂v𝑥/∂𝑥 + ∂u𝑥/∂𝑥 ∂v𝑦/∂𝑦)
 *               + ∫∫ μ(∂u𝑦/∂𝑥 + ∂u𝑥/∂𝑦)(∂v𝑥/∂𝑦 + ∂v𝑦/∂𝑥)
 *
 *   - The first term is "normal strain energy"
 *   - The second term is "compressibility effect"
 *   - The third term is "shear energy"
 *
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<6, 6> computeElementMatrixTRIA3Base(Real3 dxu, Real3 dyu, Real area, Real lambda, Real mu)
{
  FixedMatrix<1, 6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  FixedMatrix<1, 6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
  FixedMatrix<1, 6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
  FixedMatrix<1, 6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

  // ∫∫ (λ+2μ)(∂u𝑥/∂𝑥 ∂v𝑥/∂𝑥 + ∂u𝑦/∂𝑦 ∂v𝑦/∂𝑦)
  FixedMatrix<6, 6> normal_strain_energy = (lambda + 2 * mu) * ((dxUx ^ dxUx) + (dyUy ^ dyUy)) * area;

  // ∫∫ λ(∂u𝑦/∂𝑦 ∂v𝑥/∂𝑥 + ∂u𝑥/∂𝑥 ∂v𝑦/∂𝑦)
  FixedMatrix<6, 6> compressibility_effect = (lambda) * ((dyUy ^ dxUx) + (dxUx ^ dyUy)) * area;

  // ∫∫ μ(∂u𝑦/∂𝑥 + ∂u𝑥/∂𝑦)(∂v𝑥/∂𝑦 + ∂v𝑦/∂𝑥)
  FixedMatrix<6, 6> shear_energy = (mu) * ((dxUy + dyUx) ^ (dyUx + dxUy)) * area;

  return (normal_strain_energy + compressibility_effect + shear_energy);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<6, 6> computeElementMatrixTRIA3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  return computeElementMatrixTRIA3Base(dxu, dyu, area, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<6, 6> FemModule::_computeElementMatrixTRIA3(Cell cell)
{
  Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  return computeElementMatrixTRIA3Base(dxu, dyu, area, lambda, mu);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<2, 6> computeElementVectorTRIA3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu, Int32 node_lid)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  FixedMatrix<1, 6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  FixedMatrix<1, 6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
  FixedMatrix<1, 6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
  FixedMatrix<1, 6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

  // ∫∫ (λ+2μ)(∂u𝑥/∂𝑥 ∂v𝑥/∂𝑥 + ∂u𝑦/∂𝑦 ∂v𝑦/∂𝑦)
  FixedMatrix<1, 6> normal_strain_energy_x = (lambda + 2 * mu) * ((dxUx(0,node_lid*2) * dxUx) + (dyUy(0,node_lid*2) * dyUy)) * area;
  FixedMatrix<1, 6> normal_strain_energy_y = (lambda + 2 * mu) * ((dxUx(0,node_lid*2+1) * dxUx) + (dyUy(0,node_lid*2+1) * dyUy)) * area;

  // ∫∫ λ(∂u𝑦/∂𝑦 ∂v𝑥/∂𝑥 + ∂u𝑥/∂𝑥 ∂v𝑦/∂𝑦)
  FixedMatrix<1, 6> compressibility_effect_x = (lambda) * ((dyUy(0,node_lid*2) * dxUx) + (dxUx(0,node_lid*2) * dyUy)) * area;
  FixedMatrix<1, 6> compressibility_effect_y = (lambda) * ((dyUy(0,node_lid*2+1) * dxUx) + (dxUx(0,node_lid*2+1) * dyUy)) * area;

  // ∫∫ μ(∂u𝑦/∂𝑥 + ∂u𝑥/∂𝑦)(∂v𝑥/∂𝑦 + ∂v𝑦/∂𝑥)
  FixedMatrix<1, 6> shear_energy_x = (mu) * ((dxUy(0,node_lid*2) + dyUx(0,node_lid*2)) * (dyUx + dxUy)) * area;
  FixedMatrix<1, 6> shear_energy_y = (mu) * ((dxUy(0,node_lid*2+1) + dyUx(0,node_lid*2+1)) * (dyUx + dxUy)) * area;

  FixedMatrix <1, 6> result_x = normal_strain_energy_x + compressibility_effect_x + shear_energy_x;
  FixedMatrix <1, 6> result_y = normal_strain_energy_y + compressibility_effect_y + shear_energy_y;

  FixedMatrix <2, 6> result;
  result(0,0)=result_x(0,0); result(0,1)=result_x(0,1); result(0,2)=result_x(0,2); result(0,3)=result_x(0,3); result(0,4)=result_x(0,4); result(0,5)=result_x(0,5);
  result(1,0)=result_y(0,0); result(1,1)=result_y(0,1); result(1,2)=result_y(0,2); result(1,3)=result_y(0,3); result(1,4)=result_y(0,4); result(1,5)=result_y(0,5);

  return result;
}