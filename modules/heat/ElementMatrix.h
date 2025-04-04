// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrix.h                                             (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Heat transfer   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (ℙ1 FE).
 *
 * This function calculates the integral of the expression:
 *       a(𝑢,𝑣) = ∫∫ λ(∂𝑢/∂𝑥 ∂𝑣/∂𝑥  + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦)dΩ + ∫∫ (𝑢𝑣/δ𝑡)dΩ
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the gradients of the shape functions.
 * 3. Return a(𝑢,𝑣);
 */
/*---------------------------------------------------------------------------*/

RealMatrix<3, 3> FemModule::
_computeElementMatrixTria3(Cell cell)
{
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  RealVector<3> U = { 1., 1., 1. };
  Real3 dxU = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyU = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real lambda = m_cell_lambda[cell];

  return lambda * (area * (dxU ^ dxU) + area * (dyU ^ dyU)) + (1 / 12.) * massMatrix(U, U) * area / dt;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<3, 3> _computeElementMatrixTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord, const ax::VariableCellRealInView& in_cell_lambda, Real in_dt)
{
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  RealVector<3> U = { 1., 1., 1. };
  Real3 dxU = FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyU = FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real lambda = in_cell_lambda[cell_lid];

  return lambda * (area * (dxU ^ dxU) + area * (dyU ^ dyU)) + (1 / 12.) * massMatrix(U, U) * area / in_dt;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<1, 3> _computeElementVectorTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord, const ax::VariableCellRealInView& in_cell_lambda, Real in_dt, Int32 node_lid)
{
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  RealVector<3> U = { 1., 1., 1. };
  Real3 dxU = FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyU = FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real lambda = in_cell_lambda[cell_lid];

  RealMatrix<3, 3> massMat = (1 / 12.) * (massMatrix(U, U)) * area / in_dt;
  Real3 massVect = {massMat(node_lid,0) , massMat(node_lid,1) , massMat(node_lid,2) };

  Real3 node_vector_integral = lambda * area * dxU[node_lid] * dxU + lambda * area * dyU[node_lid] * dyU + massVect;
  return { node_vector_integral[0], node_vector_integral[1], node_vector_integral[2] };
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a tetrahedral element (ℙ1 FE).
 *
 * This function calculates the integral of the expression:
 *       a(𝑢,𝑣) = ∫∫∫ λ(∂𝑢/∂𝑥 ∂𝑣/∂𝑥  + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦 + ∂𝑢/∂𝑧 ∂𝑣/∂𝑧)dΩ + ∫∫∫ (𝑢𝑣/δ𝑡)dΩ
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the gradients of the shape functions.
 * 3. Return a(𝑢,𝑣);
 */
/*---------------------------------------------------------------------------*/

RealMatrix<4, 4> FemModule::
_computeElementMatrixTetra4(Cell cell)
{
  Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);

  RealVector<4> U = { 1., 1., 1., 1. };
  Real4 dxU = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
  Real4 dyU = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
  Real4 dzU = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);

  return lambda * (volume * (dxU ^ dxU) + volume * (dyU ^ dyU) + volume * (dzU ^ dzU))
         + (1 / 20.) * massMatrix(U, U) * volume / dt;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<4, 4> _computeElementMatrixTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord, const ax::VariableCellRealInView& in_cell_lambda, Real in_dt)
{
  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  RealVector<4> U = { 1., 1., 1., 1. };
  Real4 dxU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);
  Real lambda = in_cell_lambda[cell_lid];

  return  lambda * (volume * (dxU ^ dxU) + volume * (dyU ^ dyU) + volume * (dzU ^ dzU))
         + (1 / 20.) * massMatrix(U, U) * volume / in_dt;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE RealMatrix<1, 4> _computeElementVectorTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord, const ax::VariableCellRealInView& in_cell_lambda, Real in_dt, Int32 node_lid)
{
  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  RealVector<4> U = { 1., 1., 1., 1. };
  Real4 dxU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);
  Real lambda = in_cell_lambda[cell_lid];
  RealMatrix<4, 4> massMat = (1 / 20.) * (massMatrix(U, U)) * volume / in_dt;
  Real4 massVect = {massMat(node_lid,0) , massMat(node_lid,1) , massMat(node_lid,2), massMat(node_lid,3) };

  Real4 node_vector_integral = lambda * volume * dxU[node_lid] * dxU + lambda * volume * dyU[node_lid] * dyU + lambda * volume * dzU[node_lid] * dzU + massVect;

  return { node_vector_integral[0], node_vector_integral[1], node_vector_integral[2], node_vector_integral[3] };
}
