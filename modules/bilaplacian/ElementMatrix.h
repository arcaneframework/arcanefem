// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrix.h                                             (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for bilaplacian    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (ℙ1 FE).
 *
 * This function calculates the integral of the expression:
 *       a(𝐮,𝐯) =  ∫∫ (∂𝑢𝑥/∂𝑥 ∂𝑣𝑦/∂𝑥 + ∂𝑢𝑦/∂𝑥 ∂𝑣𝑦/∂𝑦)dΩ
 *               + ∫∫ (∂𝑢𝑦/∂𝑥 ∂𝑣𝑥/∂𝑥 + ∂𝑣𝑦/∂𝑦 ∂𝑣𝑥/∂𝑦)dΩ
 *               + ∫∫ (𝑢𝑦 𝑣𝑦)dΩ
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the gradients of the shape functions.
 * 3. Return a(𝑢,𝑣);
 */
/*---------------------------------------------------------------------------*/

RealMatrix<6, 6> FemModule::
_computeElementMatrixTria3(Cell cell)
{
  Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  RealVector<6> Uy = { 0., 1., 0., 1., 0., 1. };
  RealVector<6> Ux = { 1., 0., 1., 0., 1., 0. };

  RealVector<6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
  RealVector<6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
  RealVector<6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
  RealVector<6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

  RealMatrix<6, 6> int_Omega_i = ((dxUx ^ dxUy) + (dyUx ^ dyUy)) * area + ((dxUy ^ dxUx) + (dyUy ^ dyUx)) * area + massMatrix(Uy, Uy) * area;

  return int_Omega_i;
}