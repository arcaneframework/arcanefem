// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrix.h                                             (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Fourier        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (ℙ1 FE).
 *
 * This function calculates the expression:
 *       a(𝑢,𝑣) = ∫∫ λ(∂𝑢/∂𝑥 ∂𝑣/∂𝑥  + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦)dΩ
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the gradients of the shape functions.
 * 3. Return a(𝑢,𝑣);
 */
/*---------------------------------------------------------------------------*/

RealMatrix<3, 3> FemModule::_computeElementMatrixTria3(Cell cell)
{
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  Real3 dxU = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyU = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

  return area * lambda * (dxU ^ dxU) + area * lambda * (dyU ^ dyU);
}

ARCCORE_HOST_DEVICE RealMatrix<3, 3> computeElementMatrixTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord, const ax::VariableCellRealInView& in_cell_lambda)
{
  Real area = FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  Real3 dxU = FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyU = FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real in_lambda = in_cell_lambda[cell_lid];

  return area * in_lambda * (dxU ^ dxU) + area * in_lambda * (dyU ^ dyU);
}

ARCCORE_HOST_DEVICE RealMatrix<1, 3> computeElementVectorTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord, const ax::VariableCellRealInView& in_cell_lambda, Int32 node_lid)
{
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  Real3 dxU = FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyU = FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real in_lambda = in_cell_lambda[cell_lid];

  Real3 node_vector_integral = area * in_lambda * dxU[node_lid] * dxU + area * in_lambda * dyU[node_lid] * dyU;
  return { node_vector_integral[0], node_vector_integral[1], node_vector_integral[2] };
}
/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a tetrahedral element (ℙ1 FE).
 *
 * This function calculates the integral of the expression:
 *       a(𝑢,𝑣) = ∫∫∫ λ(∂𝑢/∂𝑥 ∂𝑣/∂𝑥  + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦 + ∂𝑢/∂𝑧 ∂𝑣/∂𝑧)dΩ
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the gradients of the shape functions.
 * 3. Return a(𝑢,𝑣);
 */
/*---------------------------------------------------------------------------*/

RealMatrix<4, 4> FemModule::_computeElementMatrixTetra4(Cell cell)
{
  Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);

  Real4 dxU = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
  Real4 dyU = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
  Real4 dzU = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);

  return volume * lambda* (dxU ^ dxU) + volume * lambda * (dyU ^ dyU) + volume * lambda * (dzU ^ dzU);
}

ARCCORE_HOST_DEVICE RealMatrix<4, 4> computeElementMatrixTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord, const ax::VariableCellRealInView& in_cell_lambda)
{
  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  Real4 dxU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);
  Real in_lambda = in_cell_lambda[cell_lid];

  return volume * in_lambda * (dxU ^ dxU) + volume * in_lambda * (dyU ^ dyU) + volume * in_lambda  * (dzU ^ dzU);
}

ARCCORE_HOST_DEVICE RealMatrix<1, 4> computeElementVectorTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord, const ax::VariableCellRealInView& in_cell_lambda, Int32 node_lid)
{
  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  Real4 dxU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);
  Real in_lambda = in_cell_lambda[cell_lid];

  Real4 node_vector_integral = volume * in_lambda * dxU[node_lid] * dxU + volume * in_lambda * dyU[node_lid] * dyU + volume * in_lambda * dzU[node_lid] * dzU;

  return { node_vector_integral[0], node_vector_integral[1], node_vector_integral[2], node_vector_integral[3] };
}