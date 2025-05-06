// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* utilFEM.h                                                   (C) 2022-2024 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "utilFEM.h"

using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
Real real3x3Trace(const Real3x3& mat) {
  return mat[0][0] + mat[1][1] + mat[2][2];
}
/*---------------------------------------------------------------------------*/
Real3 real3x3GetSupOutdiagonal(const Real3x3& mat) {
  return {mat[0][1], mat[0][2], mat[1][2]};
}
/*---------------------------------------------------------------------------*/
Real3 real3x3GetLowOutdiagonal(const Real3x3& mat){
  return {mat[1][0], mat[2][0], mat[2][1]};
}
/*---------------------------------------------------------------------------*/
Real3x3 diagonalReal3x3(const Real3x3& mat){
  Real3x3 newmat;
  newmat[0][0] = mat[0][0];
  newmat[1][1] = mat[1][1];
  newmat[2][2] = mat[2][2];
  return newmat;
}
/*---------------------------------------------------------------------------*/
Real3x3 outdiagonalReal3x3(const Real3x3& mat){
  return (mat - diagonalReal3x3(mat));
}

bool	real3x3IsSym(const Real3x3& mat)
{
  Real3 matsup = real3x3GetSupOutdiagonal(mat);
  Real3 matlow = real3x3GetLowOutdiagonal(mat);
  return (matsup == matlow);
}