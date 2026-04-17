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

RealMatrix<3, 3> FemModuleFourierNL::_computeElementMatrixTria3(Cell cell)
{
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  Real3 dxU = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyU = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

  Real lambda_cell = (m_node_lambda[cell.nodeId(0)] + m_node_lambda[cell.nodeId(1)] + m_node_lambda[cell.nodeId(2)]) / 3.;
  info() << "lambda_cell[" << cell.uniqueId() << "] = "<< lambda_cell ;
  return (area * lambda_cell *  (dxU ^ dxU) + area  * lambda_cell * (dyU ^ dyU));

  // Real uk_cell = (m_uk[cell.nodeId(0)] + m_uk[cell.nodeId(1)] + m_uk[cell.nodeId(2)]) / 3.;
  // Real lambda_cell = math::pow(1 + uk_cell, options()->expNlin);

  // Real3 lambda_this_node{ lambda_cell, lambda_cell, lambda_cell};


  // Real3 lambda_this_node{ m_node_lambda[cell.nodeId(0)], m_node_lambda[cell.nodeId(1)], m_node_lambda[cell.nodeId(2)]};

  // info() << "-------";
  // info() << "Cell: " << cell.uniqueId() ;
  // info() << "lambda0: "<< lambda_this_node[0] << " lambda1: " << lambda_this_node[1] << " lambda2: " << lambda_this_node[2];
  // info() << "dxu0: "<< dxU[0] << " dxu1: " << dxU[1] << " dxu2: " << dxU[2];
  // info() << "dyu0: "<< dyU[0] << " dxy1: " << dyU[1] << " dyu2: " << dyU[2];

  // Real3 lmbda_dxU{lambda_this_node[0] * dxU[0], lambda_this_node[1] * dxU[1], lambda_this_node[2] * dxU[2] };
  // Real3 lmbda_dyU{lambda_this_node[0] * dyU[0], lambda_this_node[1] * dyU[1], lambda_this_node[2] * dyU[2] };

  // info() << "lmbda_dxu0: "<< lmbda_dxU[0] << " lmbda_dxu1: " << lmbda_dxU[1] << " lmbda_dxu2: " << lmbda_dxU[2];
  // info() << "lmbda_dyu0: "<< lmbda_dyU[0] << " lmbda_dyu1: " << lmbda_dyU[1] << " lmbda_dyu2: " << lmbda_dyU[2];
  // info() << "--------------";

  // Real uk_cell = (m_uk[cell.nodeId(0)] + m_uk[cell.nodeId(1)] + m_uk[cell.nodeId(2)]) / 3.;
  // Real lambda_cell = math::pow(1 + uk_cell, options()->expNlin);


  // return (area *  (lmbda_dxU ^ dxU) + area  * (lmbda_dyU ^ dyU));

  // RealMatrix<3, 3> eMat;
  // eMat.fill(0.0);
  //
  // eMat(0, 0) = m_node_lambda[cell.nodeId(0)] * dxU[0] * dxU[0]
  //                 + m_node_lambda[cell.nodeId(0)] * dyU[0] * dyU[0];
  // eMat(0, 1) = m_node_lambda[cell.nodeId(0)] * dxU[0] * dxU[1]
  //                 + m_node_lambda[cell.nodeId(0)] * dyU[0] * dyU[1];
  // eMat(0, 2) = m_node_lambda[cell.nodeId(0)] * dxU[0] * dxU[2]
  //                 + m_node_lambda[cell.nodeId(0)] * dyU[0] * dyU[2];
  //
  // eMat(1, 0) = m_node_lambda[cell.nodeId(1)] * dxU[1] * dxU[0]
  //                 + m_node_lambda[cell.nodeId(1)] * dyU[1] * dyU[0];
  // eMat(1, 1) = m_node_lambda[cell.nodeId(1)] * dxU[1] * dxU[1]
  //                 + m_node_lambda[cell.nodeId(1)] * dyU[1] * dyU[1];
  // eMat(1, 2) = m_node_lambda[cell.nodeId(1)] * dxU[1] * dxU[2]
  //                 + m_node_lambda[cell.nodeId(1)] * dyU[1] * dyU[2];
  //
  // eMat(2, 0) = m_node_lambda[cell.nodeId(2)] * dxU[2] * dxU[0]
  //                 + m_node_lambda[cell.nodeId(2)] * dyU[2] * dyU[0];
  // eMat(2, 1) = m_node_lambda[cell.nodeId(2)] * dxU[2] * dxU[1]
  //                 + m_node_lambda[cell.nodeId(2)] * dyU[2] * dyU[1];
  // eMat(2, 2 ) = m_node_lambda[cell.nodeId(2)] * dxU[2] * dxU[2]
  //                 + m_node_lambda[cell.nodeId(2)] * dyU[2] * dyU[2];
  //
  // return eMat;

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

RealMatrix<4, 4> FemModuleFourierNL::_computeElementMatrixTetra4(Cell cell)
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