// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* GaussQuadrature.cpp                                         (C) 2022-2024 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "arcane/MathUtils.h"
#include <arcane/IParallelMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/geometry/IGeometry.h>
#include <arcane/VariableTypes.h>

#include "GaussQuadrature.h"

using namespace Arcane;
using namespace Arcane::FemUtils;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real getRefPosition(const Integer& indx,const Integer& ordre){
  Real x = xgauss1; // default is order 1

  switch(ordre){
  case 2: x = xgauss2[indx]; break;
  case 3: x = xgauss3[indx]; break;
  case 4: x = xgauss4[indx]; break;
  case 5: x = xgauss5[indx]; break;
  case 6: x = xgauss6[indx]; break;
  case 7: x = xgauss7[indx]; break;
  case 8: x = xgauss8[indx]; break;
  case 9: x = xgauss9[indx]; break;
  default: break;
  }
  return x;
}

Real getWeight(const Integer& indx,const Integer& ordre){
  Real w = wgauss1; // default is order 1

  switch(ordre){
  case 2: w = wgauss2[indx]; break;
  case 3: w = wgauss3[indx]; break;
  case 4: w = wgauss4[indx]; break;
  case 5: w = wgauss5[indx]; break;
  case 6: w = wgauss6[indx]; break;
  case 7: w = wgauss7[indx]; break;
  case 8: w = wgauss8[indx]; break;
  case 9: w = wgauss9[indx]; break;
  default: break;
  }
  return w;
}

/*---------------------------------------------------------------------------*/

Real3 LineRefPosition(const Integer3& indices,const Integer3& ordre){
  return {getRefPosition(indices[0],ordre[0]),0.,0.};
}

Real LineWeight(const Integer3& indices,const Integer3& ordre){
  return getWeight(indices[0],ordre[0]);
}

/*---------------------------------------------------------------------------*/

Real3 TriRefPosition(const Integer3& indices,const Integer3& ordre){
  Integer o = ordre[0]-1;
  Integer i = indices[0];
  return {xg1[o][i],xg2[o][i],0.};
}

Real TriWeight(const Integer3& indices,const Integer3& ordre){
  return wg[ordre[0]-1][indices[0]];
}

/*---------------------------------------------------------------------------*/

Real3 QuadRefPosition(const Integer3& indices,const Integer3& ordre){
  Real3 pos;
  for (Integer i = 0; i < 2; i++) pos[i] = getRefPosition(indices[i],ordre[i]);
  return pos;
}

Real QuadWeight(const Integer3& indices,const Integer3& ordre){
  Real w = 1.;
  for (Integer i = 0; i < 2; i++) w *= getWeight(indices[i],ordre[i]);
  return w;
}

/*---------------------------------------------------------------------------*/

Real3 HexaRefPosition(const Integer3& indices,const Integer3& ordre){
  Real3 pos;
  for (Integer i = 0; i < 3; i++) pos[i] = getRefPosition(indices[i],ordre[i]);
  return pos;
}

Real HexaWeight(const Integer3& indices,const Integer3& ordre){
  Real w = 1.;
  for (Integer i = 0; i < 3; i++) w *= getWeight(indices[i],ordre[i]);
  return w;
}

/*---------------------------------------------------------------------------*/

Real3 TetraRefPosition(const Integer3& indices,const Integer3& /*ordre*/){
  Integer i = indices[0];
  return {xit[i],yit[i],zit[i]};
}

Real TetraWeight(const Integer3& indices,const Integer3& /*ordre*/){
  return wgtetra;
}

/*---------------------------------------------------------------------------*/

Real3 PentaRefPosition(const Integer3& indices,const Integer3& ordre){

  // Same as TriRefPosition on reference coordinate plane (r,s)
  // and LineRefPosition along reference coordinate t (vertical)
  auto pos = TriRefPosition(indices,ordre);
  pos.z = getRefPosition(indices[2],ordre[2]);

  return pos;
}

Real PentaWeight(const Integer3& indices,const Integer3& ordre){

  // Same as TriWeight on reference coordinate plane (r,s)
  // and LineWeight with ordre[2] to account for reference coordinate t (vertical)
  Real wgpenta = TriWeight(indices,ordre)*getWeight(indices[2],ordre[2]);
  return wgpenta;
}

/*---------------------------------------------------------------------------*/

Integer getNbGaussPointsfromOrder(const Int16& cell_type, const Integer& ninteg){
  Integer nbgauss{0};
  auto ninteg2{ ninteg * ninteg };
  auto ninteg3{ ninteg2 * ninteg };

  if (ninteg <= 1)
    nbgauss = 1;
  else if (cell_type == IT_Line2 || cell_type == IT_Line3)
    nbgauss = ninteg;
  else if (cell_type == IT_Quad4 || cell_type == IT_Quad8)
    nbgauss = ninteg2;
  else if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20)
    nbgauss = ninteg3;

  else if (ninteg == 2) {
    switch (cell_type) {
    default: break;

    case IT_Triangle3:
    case IT_Triangle6: nbgauss = 3;
      break;

    case IT_Tetraedron4:
    case IT_Tetraedron10: nbgauss = 4;
      break;

    case IT_Pentaedron6: nbgauss = 6;
      break;

    case IT_Pyramid5: nbgauss = 5;
      break;
    }
  }
  else if (ninteg == 3) {
    switch (cell_type) {
    default: break;

    case IT_Triangle3:
    case IT_Triangle6: nbgauss = 4;
      break;

    case IT_Tetraedron4:
    case IT_Tetraedron10: nbgauss = 5;
      break;

    case IT_Pentaedron6: nbgauss = 8;
      break;

    case IT_Pyramid5: nbgauss = 6;
      break;
    }
  }
  else if (ninteg >= 4) {
    switch (cell_type) {
    default: break;

    case IT_Triangle3:
    case IT_Triangle6: nbgauss = 7;
      break;

    case IT_Tetraedron4:
    case IT_Tetraedron10: nbgauss = 15;
      break;

    case IT_Pentaedron6: nbgauss = 21;
      break;

    case IT_Pyramid5: nbgauss = 27;
      break;
    }
  }
  return nbgauss;
}

Real3 getGaussRefPosition(const ItemWithNodes& cell, const Integer& ninteg, const Int32& rank)
{
  auto cell_type = cell.type();
  Integer nint{ninteg};

  if (nint < 1)
    nint = 1;

  if (cell_type == IT_Line2 || cell_type == IT_Line3)
    return LineRefPosition({rank,-1,-1},{nint,0,0});

  if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
    auto in{0};
    for (Int32 i1 = 0; i1 < nint; ++i1) {
      for (Int32 i2 = 0; i2 < nint; ++i2) {

        if (rank == in)
          return QuadRefPosition({i1,i2,-1},{nint,nint,0});

        ++in;
      }
    }
  }

  if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
    auto in{0};
    for (Int32 i1 = 0; i1 < nint; ++i1) {
      for (Int32 i2 = 0; i2 < nint; ++i2) {
        for (Int32 i3 = 0; i3 < nint; ++i3) {

          if (rank == in)
            return HexaRefPosition({ i1, i2, i3 }, { nint, nint, nint });

          ++in;
        }
      }
    }
  }

  if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
    auto o {3};
    if (nint <= 3)
      o = nint-1;

    return {xg1[o][rank],xg2[o][rank],0.};
  }

  if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10) {
    auto o {3};
    if (nint <= 3)
      o = nint-1;

    return {xtet[o][rank],ytet[o][rank],ztet[o][rank]};
  }

  if (cell_type == IT_Pyramid5) {
    auto o {1};
    if (nint <= 2)
      o = 0;

    return {xpyr[o][rank],ypyr[o][rank],zpyr[o][rank]};
  }

  if (cell_type == IT_Pentaedron6) {
    auto o {1};
    if (nint <= 2)
      o = 0;

    return {xpent[o][rank],ypent[o][rank],zpent[o][rank]};
  }
  return {};
}

Real getGaussWeight(const ItemWithNodes& cell, const Integer& ninteg, const Int32& rank)
{
  auto cell_type = cell.type();
  Integer nint{ninteg};

  if (nint < 1)
    nint = 1;

  if (cell_type == IT_Line2 || cell_type == IT_Line3)
    return LineWeight({rank,-1,-1},{nint,0,0});

  if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
    auto in{0};
    for (Int32 i1 = 0; i1 < nint; ++i1) {
      for (Int32 i2 = 0; i2 < nint; ++i2) {

        if (rank == in)
          return QuadWeight({i1,i2,-1},{nint,nint,0});

        ++in;
      }
    }
  }

  if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
    auto in{0};
    for (Int32 i1 = 0; i1 < nint; ++i1) {
      for (Int32 i2 = 0; i2 < nint; ++i2) {
        for (Int32 i3 = 0; i3 < nint; ++i3) {

          if (rank == in)
            return HexaWeight({ i1, i2, i3 }, { nint, nint, nint });

          ++in;
        }
      }
    }
  }

  if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
    auto o {3};
    if (nint <= 3)
      o = nint-1;

    return wg[o][rank];
  }

  if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10) {
    if (nint == 1)
      return wgtet1;
    if (nint == 2)
      return wgtet2;
    if (nint == 3) {
      auto i = npwgtet3[rank];
      return wgtet3[i];
    }
    // nint >= 4
    auto i = npwgtet4[rank];
    return wgtet4[i];
  }

  if (cell_type == IT_Pyramid5) {
    if (nint <= 2)
      return wgpyr2;

    // nint >= 3
    auto i = npwgpyr3[rank];
    return wgpyr3[i];
  }

  if (cell_type == IT_Pentaedron6) {
    if (nint <= 2)
      return wgpent2;

    // nint >= 3
    auto i = npwgpent3[rank];
    return wgpent3[i];
  }

  return 1.;
}

/*---------------------------------------------------------------------------*/
