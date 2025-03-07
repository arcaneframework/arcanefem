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
/*---------------------------------------------------------------------------*/

FixedMatrix<6, 6> FemModule::
_computeElementMatrixTRIA3(Cell cell)
{
  // Get coordiantes of the triangle element  TRI3
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

  FixedMatrix<1, 6> b_matrix;
  FixedMatrix<6, 1> bT_matrix;
  FixedMatrix<6, 6> int_Omega_i;

  for (Int32 i = 0; i<6; i++)
    for (Int32 j = 0; j<6; j++)
      int_Omega_i(i,j) = 0.;

// -----------------------------------------------------------------------------
//  c1( dx(du1)dx(v1) + dy(du2)dx(v1) + dx(du1)dy(v2) + dy(du2)dy(v2) )
//------------------------------------------------------------------------------


  // dx(u1)dx(v1) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_dxU1dxV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxU1dxV1);

  // dy(u2)dx(v1) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_dyU1dyV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyU1dyV1);

  // dx(u1)dy(v2) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_dxU2dxV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxU2dxV1);

  // dy(u2)dy(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_dyU2dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyU2dyV1);

  // lambda * (.....)
  int_Omega_i.multInPlace(c1);


// -----------------------------------------------------------------------------
//  c2( dx(u1)dx(v1) + dy(u2)dy(v2) + 0.5*(   dy(u1)dy(v1) + dx(u2)dy(v1)
//                                            + dy(u1)dx(v2) + dx(u2)dx(v2) )
//      )
//------------------------------------------------------------------------------

  // mu*dx(u1)dx(v1) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudxU1dxV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_mudxU1dxV1.multInPlace(c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU1dxV1);


  // mu*dy(u2)dy(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudyU2dyV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudyU2dyV2.multInPlace(c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU2dyV2);

  // 0.5*0.5*mu*dy(u1)dy(v1) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.y/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.y/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.y;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.y;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudyU1dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudyU1dyV1.multInPlace(0.5*c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU1dyV1);

  // 0.5*mu*dx(u2)dy(v1) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.x/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.x/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.x/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.y;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.y;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudxU2dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudxU2dyV1.multInPlace(0.5*c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU2dyV1);

  // 0.5*mu*dy(u1)dx(v2) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.y/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.y/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.x;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.x;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudyU1dxV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudyU1dxV2.multInPlace(0.5*c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU1dxV2);

  // 0.5*mu*dx(u2)dx(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.x/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.x/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.x/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.x;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.x;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudxU2dxV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudxU2dxV2.multInPlace(0.5*c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU2dxV2);

// -----------------------------------------------------------------------------
//   c0(du1v1 + du2v2)
// -----------------------------------------------------------------------------
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = 1.;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = 1.;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = area/3.;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = area/3.;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = area/3.;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_U2V2   = matrixMultiplication(bT_matrix, b_matrix);
  int_U2V2.multInPlace(c0);

  for (Int32 i = 0; i<6; i++)
    int_U2V2(i,i) *= 2.;

  int_Omega_i = matrixAddition( int_Omega_i, int_U2V2);



  // du2v2 //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = 1.;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = 1.;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = 1.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = area/3.;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = area/3.;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = area/3.;

  bT_matrix.multInPlace(0.5);

  int_U2V2   = matrixMultiplication(bT_matrix, b_matrix);
  int_U2V2.multInPlace(c0);

  for (Int32 i = 0; i<6; i++)
    int_U2V2(i,i) *= 2.;


  int_Omega_i = matrixAddition( int_Omega_i, int_U2V2);

  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<4, 4> FemModule::
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

  Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);
  Real2 N   = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, m_node_coord);

  FixedMatrix<1, 4> b_matrix;
  FixedMatrix<4, 1> bT_matrix;
  FixedMatrix<4, 4> int_DOmega_i;

  for (Int32 i = 0; i<4; i++)
    for (Int32 j = 0; j<4; j++)
      int_DOmega_i(i,j) = 0.;

  // c7*u1v1 //
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = 1.;
  b_matrix(0, 3) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = length/3;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = length/3;
  bT_matrix(3, 0) = 0.;

  bT_matrix.multInPlace(1.);

  FixedMatrix<4, 4> int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<4; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(c7*(N.x*N.x*cp + N.y*N.y*cs));
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  // c7*u2v2 //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = 1.;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = 1.;

  b_matrix.multInPlace(0.5f);


  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = length/3;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = length/3;

  bT_matrix.multInPlace(1.);

  int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<4; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(c7*(N.y*N.y*cp + N.x*N.x*cs));
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  // c7*u1v2 //
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = 1.;
  b_matrix(0, 3) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = length/3;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = length/3;

  bT_matrix.multInPlace(1.);

  int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  int_UV.multInPlace(c7*(N.x*N.y*(cp - cs)) );
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  // c7*u2v1 //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = 1.;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = 1.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = length/3;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = length/3;
  bT_matrix(3, 0) = 0.;

  bT_matrix.multInPlace(1.);

  int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  int_UV.multInPlace( c7*(N.x*N.y*(cp - cs)) );
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  return int_DOmega_i;
}