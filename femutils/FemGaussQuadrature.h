// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemGaussQuadrature.h                                        (C) 2000-2026 */
/*                                                                           */
/* Provides methods for Gauss quadrature.                                    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMTUTILS_FEMGAUSSQUADRATURE_H
#define ARCANEFEM_FEMTUTILS_FEMGAUSSQUADRATURE_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemUtils.h"
#include "Integer3std.h"
#include "GaussQuadrature.h"

#include <arcane/core/Item.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace ArcaneFemFunctions
{

using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods for Gauss quadrature.
 *
 * This class includes static methods for computing Gauss-Legendre integration
 * depending on finite element types (coming from PASSMO).
 */
class FemGaussQuadrature
{
 public:

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the number of Gauss Points for a given finite element type,
   * depending on the integration order chosen by user (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  static inline Integer getNbGaussPointsfromOrder(Int16 cell_type, Integer ninteg)
  {
    Integer nbgauss{ 0 };
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
      default:
        break;

      case IT_Triangle3:
      case IT_Triangle6:
        nbgauss = 3;
        break;

      case IT_Tetraedron4:
      case IT_Tetraedron10:
        nbgauss = 4;
        break;

      case IT_Pentaedron6:
        nbgauss = 6;
        break;

      case IT_Pyramid5:
        nbgauss = 5;
        break;
      }
    }
    else if (ninteg == 3) {
      switch (cell_type) {
      default:
        break;

      case IT_Triangle3:
      case IT_Triangle6:
        nbgauss = 4;
        break;

      case IT_Tetraedron4:
      case IT_Tetraedron10:
        nbgauss = 5;
        break;

      case IT_Pentaedron6:
        nbgauss = 8;
        break;

      case IT_Pyramid5:
        nbgauss = 6;
        break;
      }
    }
    else if (ninteg >= 4) {
      switch (cell_type) {
      default:
        break;

      case IT_Triangle3:
      case IT_Triangle6:
        nbgauss = 7;
        break;

      case IT_Tetraedron4:
      case IT_Tetraedron10:
        nbgauss = 15;
        break;

      case IT_Pentaedron6:
        nbgauss = 21;
        break;

      case IT_Pyramid5:
        nbgauss = 27;
        break;
      }
    }
    return nbgauss;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the position of a Gauss Point in the reference (local) element
   * for a given finite element, depending on the rank of the point in the
   * element loop (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 getGaussRefPosition(ItemWithNodes cell, Integer ninteg, Integer rank)
  {
    auto cell_type = cell.type();
    Integer nint{ ninteg };

    if (nint < 1)
      nint = 1;

    if (cell_type == IT_Line2 || cell_type == IT_Line3)
      return lineRefPosition({ rank, -1, -1 }, { nint, 0, 0 });

    if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
      auto in{ 0 };
      for (Int32 i1 = 0; i1 < nint; ++i1) {
        for (Int32 i2 = 0; i2 < nint; ++i2) {

          if (rank == in)
            return quadRefPosition({ i1, i2, -1 }, { nint, nint, 0 });

          ++in;
        }
      }
    }

    if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
      auto in{ 0 };
      for (Int32 i1 = 0; i1 < nint; ++i1) {
        for (Int32 i2 = 0; i2 < nint; ++i2) {
          for (Int32 i3 = 0; i3 < nint; ++i3) {

            if (rank == in)
              return hexaRefPosition({ i1, i2, i3 }, { nint, nint, nint });

            ++in;
          }
        }
      }
    }

    if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
      auto o{ 3 };
      if (nint <= 3)
        o = nint - 1;

      return { xg1[o][rank], xg2[o][rank], 0. };
    }

    if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10) {
      auto o{ 3 };
      if (nint <= 3)
        o = nint - 1;

      return { xtet[o][rank], ytet[o][rank], ztet[o][rank] };
    }

    if (cell_type == IT_Pyramid5) {
      auto o{ 1 };
      if (nint <= 2)
        o = 0;

      return { xpyr[o][rank], ypyr[o][rank], zpyr[o][rank] };
    }

    if (cell_type == IT_Pentaedron6) {
      auto o{ 1 };
      if (nint <= 2)
        o = 0;

      return { xpent[o][rank], ypent[o][rank], zpent[o][rank] };
    }
    return {};
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for a given
   * finite element, depending on the rank of the point in the element loop
   * (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  static inline Real getGaussWeight(ItemWithNodes cell, Integer ninteg, Integer rank)
  {
    auto cell_type = cell.type();
    Integer nint{ ninteg };

    if (nint < 1)
      nint = 1;

    if (cell_type == IT_Line2 || cell_type == IT_Line3)
      return lineWeight({ rank, -1, -1 }, { nint, 0, 0 });

    if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
      auto in{ 0 };
      for (Int32 i1 = 0; i1 < nint; ++i1) {
        for (Int32 i2 = 0; i2 < nint; ++i2) {

          if (rank == in)
            return quadWeight({ i1, i2, -1 }, { nint, nint, 0 });

          ++in;
        }
      }
    }

    if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
      auto in{ 0 };
      for (Int32 i1 = 0; i1 < nint; ++i1) {
        for (Int32 i2 = 0; i2 < nint; ++i2) {
          for (Int32 i3 = 0; i3 < nint; ++i3) {

            if (rank == in)
              return hexaWeight({ i1, i2, i3 }, { nint, nint, nint });

            ++in;
          }
        }
      }
    }

    if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
      auto o{ 3 };
      if (nint <= 3)
        o = nint - 1;

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
  /**
   * @brief Provides the position of a Gauss Point in the reference (local) element
   * coordinates, depending on the index of the Gauss Point and integration order
   * chosen by user (coming fro PASSMO).
   * This method is generic (not depending on the finite element type) and is called
   * by specialized methods (depending on FE type)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real getRefPosition(Integer indx, Integer ordre)
  {
    Real x = xgauss1; // default is order 1

    switch (ordre) {
    case 2:
      x = xgauss2[indx];
      break;
    case 3:
      x = xgauss3[indx];
      break;
    case 4:
      x = xgauss4[indx];
      break;
    case 5:
      x = xgauss5[indx];
      break;
    case 6:
      x = xgauss6[indx];
      break;
    case 7:
      x = xgauss7[indx];
      break;
    case 8:
      x = xgauss8[indx];
      break;
    case 9:
      x = xgauss9[indx];
      break;
    default:
      break;
    }
    return x;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight depending on the index of the Gauss Point
   * and integration order chosen by user (coming fro PASSMO).
   * This method is generic (not depending on the finite element type) and is called
   * by specialized methods (depending on FE type)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real getWeight(Integer indx, Integer ordre)
  {
    Real w = wgauss1; // default is order 1

    switch (ordre) {
    case 2:
      w = wgauss2[indx];
      break;
    case 3:
      w = wgauss3[indx];
      break;
    case 4:
      w = wgauss4[indx];
      break;
    case 5:
      w = wgauss5[indx];
      break;
    case 6:
      w = wgauss6[indx];
      break;
    case 7:
      w = wgauss7[indx];
      break;
    case 8:
      w = wgauss8[indx];
      break;
    case 9:
      w = wgauss9[indx];
      break;
    default:
      break;
    }
    return w;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for edge elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/

  static inline Real3 lineRefPosition(Integer3 indices, Integer3 ordre)
  {
    return { getRefPosition(indices[0], ordre[0]), 0., 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for triangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 triRefPosition(Integer3 indices, Integer3 ordre)
  {
    Integer o = ordre[0] - 1;
    Integer i = indices[0];
    return { xg1[o][i], xg2[o][i], 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for quadrangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 quadRefPosition(Integer3 indices, Integer3 ordre)
  {
    Real3 pos;
    for (Integer i = 0; i < 2; i++)
      pos[i] = getRefPosition(indices[i], ordre[i]);
    return pos;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for hexaedron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 hexaRefPosition(Integer3 indices, Integer3 ordre)
  {
    Real3 pos;
    for (Integer i = 0; i < 3; i++)
      pos[i] = getRefPosition(indices[i], ordre[i]);
    return pos;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for tetraedron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  [[maybe_unused]] static inline Real3 tetraRefPosition(Integer3 indices, Integer3 /*ordre*/)
  {
    Integer i = indices[0];
    return { xit[i], yit[i], zit[i] };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for wedge (pentaedron)
   * elements (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  [[maybe_unused]] static inline Real3 pentaRefPosition(Integer3 indices, Integer3 ordre)
  {
    // Same as TriRefPosition on reference coordinate plane (r,s)
    // and LineRefPosition along reference coordinate t (vertical)
    auto pos = triRefPosition(indices, ordre);
    pos.z = getRefPosition(indices[2], ordre[2]);

    return pos;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for edge elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real lineWeight(Integer3 indices, Integer3 ordre)
  {
    return getWeight(indices[0], ordre[0]);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for triangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real triWeight(Integer3 indices, Integer3 ordre)
  {
    return wg[ordre[0] - 1][indices[0]];
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for quadrangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real quadWeight(Integer3 indices, Integer3 ordre)
  {
    Real w = 1.;
    for (Integer i = 0; i < 2; i++)
      w *= getWeight(indices[i], ordre[i]);
    return w;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for hexadron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real hexaWeight(Integer3 indices, Integer3 ordre)
  {
    Real w = 1.;
    for (Integer i = 0; i < 3; i++)
      w *= getWeight(indices[i], ordre[i]);
    return w;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for tetraedron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  [[maybe_unused]] static inline Real tetraWeight(Integer3 /*indices*/, Integer3 /*ordre*/)
  {
    return wgtetra;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for wedge (pentaedron)
   * elements (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  [[maybe_unused]] static inline Real pentaWeight(Integer3 indices, Integer3 ordre)
  {

    // Same as TriWeight on reference coordinate plane (r,s)
    // and LineWeight with ordre[2] to account for reference coordinate t (vertical)
    Real wgpenta = triWeight(indices, ordre) * getWeight(indices[2], ordre[2]);
    return wgpenta;
  }

  /*---------------------------------------------------------------------------*/
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace ArcaneFemFunctions

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
