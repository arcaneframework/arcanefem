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
_computeElementMatrixTRIA3(Cell cell)
{
  // Get coordinates of the triangle element  TRI3
  //------------------------------------------------
  //                  0 o
  //                   . .
  //                  .   .
  //                 .     .
  //              1 o . . . o 2
  //------------------------------------------------
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real area =  ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  RealMatrix<1, 3> b_matrix;
  RealMatrix<3, 1> bT_matrix;
  RealMatrix<3, 3> int_Omega_i;

  for (Int32 i = 0; i<3; i++)
    for (Int32 j = 0; j<3; j++)
      int_Omega_i(i,j) = 0.;

// -----------------------------------------------------------------------------
//  lambda*(dx(u)dx(v) + dy(u)dy(v)) + uv/dt
//------------------------------------------------------------------------------


  // dx(u)dx(v) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = dPhi1.x/area;
  b_matrix(0, 2) = dPhi2.x/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = dPhi1.x;
  bT_matrix(2, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.5f);

  RealMatrix<3, 3> int_dxUdxV = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxUdxV);


  // dy(u)dy(v) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = dPhi1.y/area;
  b_matrix(0, 2) = dPhi2.y/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = dPhi1.y;
  bT_matrix(2, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5f);

  RealMatrix<3, 3> int_dyUdyV = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyUdyV);

  int_Omega_i.multInPlace(lambda);

  // uv //
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 1.;
  b_matrix(0, 2) = 1.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = area/3.;
  bT_matrix(1, 0) = area/3.;
  bT_matrix(2, 0) = area/3.;

  bT_matrix.multInPlace(0.5f);

  RealMatrix<3, 3> int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<3; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(1./dt);

  int_Omega_i = matrixAddition( int_Omega_i, int_UV);

  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

RealMatrix<2, 2> FemModule::
_computeElementMatrixEDGE2(Face face)
{
  // Get coordinates of the triangle element  EDGE2
  //------------------------------------------------
  //                   1         0
  //                   o . . . . o
  //
  //------------------------------------------------
  Real3 m0 = m_node_coord[face.nodeId(0)];
  Real3 m1 = m_node_coord[face.nodeId(1)];

  Real area =  ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);

  Real2 dPhi0(m0.y - m1.y, m1.x - m0.x);
  Real2 dPhi1(m1.y - m0.y, m0.x - m1.x);

  RealMatrix<1, 2> b_matrix;
  RealMatrix<2, 1> bT_matrix;
  RealMatrix<2, 2> int_DOmega_i;

  for (Int32 i = 0; i<2; i++)
    for (Int32 j = 0; j<2; j++)
      int_DOmega_i(i,j) = 0.;

  // uv //
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 1.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = area/3;
  bT_matrix(1, 0) = area/3;

  bT_matrix.multInPlace(1.);

  RealMatrix<2, 2> int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<2; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(h);
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  return int_DOmega_i;
}