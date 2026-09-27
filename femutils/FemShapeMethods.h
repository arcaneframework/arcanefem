// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemGaussMethods.h                                           (C) 2000-2026 */
/*                                                                           */
/* Provides methods for various FEM-related operations.                      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMTUTILS_FEMGAUSSMETHODS_H
#define ARCANEFEM_FEMTUTILS_FEMGAUSSMETHODS_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemUtils.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace ArcaneFemFunctions
{

using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods for various FEM-related operations.
 *
 * This class includes static methods for computing shape functions and their
 * derivatives depending on finite element types. These methods are used within
 * the Dispatcher mechanism available in Arcane through the class
 * CellFEMDispatcher (coming from PASSMO).
 */
/*---------------------------------------------------------------------------*/
class FemShapeMethods
{
 public:

  /*---------------------------------------------------------------------------*/
  /**
     * @brief Provides methods for reference linear (P1) edge finite-element
     * The "Line2" reference element is assumed as follows:
     *  0           1
     *  o-----------o---> x
     * -1           1
     * direct local numbering : 0->1
     */
  /*---------------------------------------------------------------------------*/
  static inline Real line2ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 2);
#endif

    Real r = ref_coord[0];
    if (inod == 1)
      return (0.5 * (1 + r));
    return (0.5 * (1 - r));
  }

  static inline Real3 line2ShapeFuncDeriv(Integer inod, Real3)
  {
    if (inod == 1)
      return { 0.5, 0., 0. };
    return { -0.5, 0., 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
     * @brief Provides methods for reference quadratic (P2) edge finite-element
     * The "Line3" reference element is assumed as follows:
     *  0     2      1
     *  o-----o------o---> x
     * -1     0      1
     * direct local numbering : 0->1->2
     */
  /*---------------------------------------------------------------------------*/
  static inline Real line3ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 3);
#endif

    Real ri = ref_coord[0];
    if (inod == 0)
      ri *= -1;

    if (inod < 2)
      return 0.5 * ri * (1 + ri); // nodes 0 or 1
    return (1 - ri * ri); // middle node
  }

  static inline Real3 line3ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    Real ri = ref_coord[0];
    if (!inod)
      return { -0.5 + ri, 0., 0. };
    if (inod == 1)
      return { 0.5 + ri, 0., 0. };
    return { -2. * ri, 0., 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
     * @brief Provides methods for reference linear (P1) triangle finite-element
     * The "Tri3" reference element is assumed as follows:
     *
     *   ^ s
     *   |
     *  2 (1,0)
     *   o
     *   . .
     *   .   .
     *   .     .
     *   .       .
     *   .         .
     *   .           .
     *   o-------------o---------> r
     *  0 (0,0)         1 (1,0)
     * direct local numbering : 0->1->2
     */
  /*---------------------------------------------------------------------------*/
  static inline Real tri3ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 3);
#endif
    Real r = ref_coord[0];
    Real s = ref_coord[1];
    if (!inod)
      return (1 - r - s);
    if (inod == 1)
      return r;
    return s;
  }

  static inline Real3 tri3ShapeFuncDeriv(Integer inod, Real3)
  {
    if (!inod)
      return { -1., -1., 0. };
    if (inod == 1)
      return { 1., 0., 0. };
    return { 0., 1., 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
     * @brief Provides methods for reference quadratic (P2) triangle finite-element
     * The "Tri6" reference element is assumed as follows:
     *   ^ s
     *   |
     *  2 (1,0)
     *   o
     *   .  .
     *   .    .
     *   .      .
     *   o 6      o 5(0.5;0.5)
     *   .(0;0.5)   .
     *   .            .
     *   .              .
     *   o-------o-------o---------> r
     * 0(0,0)  4(0.5;0)  1(1,0)
     * direct local numbering : 0->1->2->3->4->5
     */
  /*---------------------------------------------------------------------------*/
  static inline Real tri6ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 6);
#endif
    auto wi = 0., ri = ref_coord[0], si = ref_coord[1];
    auto ri2 = 2. * ri - 1.;
    auto si2 = 2. * si - 1.;
    auto ti = 1. - ri - si, ti2 = 2. * ti - 1.;

    switch (inod) {
    default:
      break;
    case 0:
      wi = ti * ti2;
      break;
    case 1:
      wi = ri * ri2;
      break;
    case 2:
      wi = si * si2;
      break;
    case 3:
      wi = 4. * ri * ti;
      break;
    case 4:
      wi = 4. * ri * si;
      break;
    case 5:
      wi = 4. * si * ti;
      break;
    }
    return wi;
  }

  static inline Real3 tri6ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    auto ri = ref_coord[0], si = ref_coord[1];
    auto ti = 1. - ri - si;

    if (!inod) {
      auto wi = -3. + 4. * (ri + si);
      return { wi, wi, 0. };
    }
    if (inod == 1)
      return { -1. + 4. * ri, 0., 0. };
    if (inod == 2)
      return { 0., -1. + 4. * si, 0. };

    if (inod == 3)
      return { 4. * (ti - ri), -4. * ri, 0. };
    if (inod == 4)
      return { 4. * si, 4. * ri, 0. };
    return { -4. * si, 4. * (ti - si), 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) quadrangle finite-element
   * The "Quad4" reference element is assumed as follows:
   *         ^y
   *          |
   *  1 o-----1-----o 0
   *    |     |     |
   *    |     |     |
   *    |     |     |
   *   -1 ----|---- 1 ---> x
   *    |     |     |
   *    |     |     |
   *    |     |     |
   *  2 o--- -1 ----o 3
   *
   * direct local numbering : 0->1->2->3
   */
  /*---------------------------------------------------------------------------*/
  static inline Real quad4ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 4);
#endif

    auto r{ ref_coord[0] }, s{ ref_coord[1] };
    auto ri{ 1. }, si{ 1. };

    switch (inod) {
    default:
      break; // default is first node (index 0)
    case 2:
      si = -1;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 3:
      si = -1;
      break;
    }
    return ((1 + ri * r) * (1 + si * s) / 4.);
  }

  static inline Real3 quad4ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 4);
#endif

    auto r{ ref_coord[0] }, s{ ref_coord[1] };
    auto ri{ 1. }, si{ 1. }; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

    switch (inod) {
    default:
      break; // default is first node (index 0)
    case 2:
      si = -1;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 3:
      si = -1;
      break;
    }
    return { 0.25 * ri * (1 + si * s), 0.25 * si * (1 + ri * r), 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference quadratic (P2) quadrangle finite-element
   * The "Quad8" reference element is assumed as follows:
   *         ^y
   *          |
   *          4
   *  1 o-----o-----o 0
   *    |     |     |
   *    |     |     |
   *    |     |     |
   *  5 o ----|---- o 7 ---> x
   *    |     |     |
   *    |     |     |
   *    |     |     |
   *  2 o-----o-----o 3
   *          6
   * Normalized coordinates (x, y) vary between -1/+1
   * Nodes 4, 6 are on line (x = 0)
   * direct local numbering :  0->1->2->...->5->6->7
   */
  /*---------------------------------------------------------------------------*/
  static inline Real quad8ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 8);
#endif
    Real tol{ 1.0e-15 };

    auto r{ ref_coord[0] }, s{ ref_coord[1] };
    auto ri{ 1. }, si{ 1. };

    switch (inod) {
    default:
      break; // default is first node (index 0)
    case 2:
      si = -1;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 3:
      si = -1;
      break;

    case 6:
      si = -1;
      [[fallthrough]];
    case 4:
      ri = 0;
      break;

    case 5:
      ri = -1;
      [[fallthrough]];
    case 7:
      si = 0;
      break;
    }

    auto r0{ r * ri }, s0{ s * si };
    Real Phi{ 0. };
    auto t0{ r0 + s0 - 1. };

    if (inod < 4) // Corner nodes
      Phi = (1 + r0) * (1 + s0) * t0 / 4.;

    else { // Middle nodes
      if (fabs(ri) < tol)
        Phi = (1 - r * r) * (1 + s0) / 2.;
      else if (fabs(si) < tol)
        Phi = (1 - s * s) * (1 + r0) / 2.;
    }
    return Phi;
  }

  static inline Real3 quad8ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    Real tol{ 1.0e-15 };

    auto r{ ref_coord[0] }, s{ ref_coord[1] };
    auto ri{ 1. }, si{ 1. };

    switch (inod) {
    default:
      break; // default is first node (index 0)
    case 2:
      si = -1;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 3:
      si = -1;
      break;

    case 6:
      si = -1;
      [[fallthrough]];
    case 4:
      ri = 0;
      break;

    case 5:
      ri = -1;
      [[fallthrough]];
    case 7:
      si = 0;
      break;
    }

    auto r0{ r * ri }, s0{ s * si };
    Real3 dPhi;
    auto t0{ r0 + s0 - 1. };

    if (inod < 4) { // Corner nodes
      dPhi.x = ri * (1 + s0) * (t0 + 1. + r0) / 4.;
      dPhi.y = si * (1 + r0) * (t0 + 1. + s0) / 4.;
    }
    else { // Middle nodes
      if (fabs(ri) < tol) {
        dPhi.x = -r * (1 + s0);
        dPhi.y = si * (1 - r * r) / 2.;
      }
      else if (fabs(si) < tol) {
        dPhi.x = -s * (1 + r0);
        dPhi.y = ri * (1 - s * s) / 2.;
      }
    }
    dPhi.z = 0.;
    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) hexaedron finite-element
   * The "Hexa8" reference element is assumed as follows:
   *     (-1, 1,1)
   *         1-------------0 (1,1,1)
   *        /|            /|
   *       / |           / |
   *     /   |          /  |
   *    2----|---------3   |   z   y
   *  (-1,-1,1)        |   |   | /
   *    |    |         |   |   |/--->x
   *    |    |         |   |
   *    |    |         |   |
   *    |    5---------|---4 (1,1,-1)
   *    |  /           |  /
   *    | /            | /
   *    |/             |/
   *    6--------------7 (1,-1,-1)
   * (-1,-1,-1)
   * Normalized coordinates (x, y, z) vary between -1/+1
   * direct local numbering : 0->1->2->3->4->5->6->7
   */
  /*---------------------------------------------------------------------------*/
  static inline Real hexa8ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 8);
#endif
    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

    switch (inod) {
    default:
      break;
    case 3:
    case 2:
      ri = -1;
      break;
    case 0:
    case 1:
      ri = -1;
      si = -1;
      break;
    case 4:
    case 5:
      si = -1;

      break;
    }
    if (inod == 1 || inod == 2 || inod == 5 || inod == 6)
      ti = -1;

    auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
    auto Phi = (1 + r0) * (1 + s0) * (1 + t0) / 8.;

    return Phi;
  }

  static inline Real3 hexa8ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {

    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

    switch (inod) {
    default:
      break;
    case 3:
    case 2:
      ri = -1;
      break;
    case 0:
    case 1:
      ri = -1;
      si = -1;
      break;
    case 4:
    case 5:
      si = -1;
      break;
    }
    if (inod == 1 || inod == 2 || inod == 5 || inod == 6)
      ti = -1;

    auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
    Real3 dPhi;
    dPhi.x = ri * (1 + s0) * (1 + t0) / 8.;
    dPhi.y = si * (1 + r0) * (1 + t0) / 8.;
    dPhi.z = ti * (1 + r0) * (1 + s0) / 8.;
    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference quadratic (P2) hexaedron finite-element
   * The "Hexa20" reference element is assumed as follows:
   *     (-1, 1,1)
   *         1------8------0 (1,1,1)
   *        /|            /|
   *      9  |          11 |
   *     /   |          /  |
   *    2----|--10-----3   |   z   y
   *  (-1,-1,1)        |   |   | /
   *    |   17         |  16   |/--->x
   *    |    |         |   |
   *   18    |        19   |
   *    |    5----12---|---4 (1,1,-1)
   *    |  /           |  /
   *    | 13           | 15
   *    |/             |/
   *    6-----14-------7 (1,-1,-1)
   * (-1,-1,-1)
   * Normalized coordinates (x, y, z) vary between -1/+1
   * direct local numbering : 0->1->2->3->...->18->19
   */
  /*---------------------------------------------------------------------------*/

  static inline Real hexa20ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 20);
#endif
    Real tol{ 1.0e-15 };

    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 0 = (1,1,1)

    switch (inod) {
    default:
      break;

    case 5:
      ti = -1.;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 6:
      ti = -1.;
      [[fallthrough]];
    case 2:
      ri = -1;
      si = -1;
      break;

    case 7:
      ti = -1.;
      [[fallthrough]];
    case 3:
      si = -1;
      break;

    case 4:
      ti = -1.;
      break;

    case 9:
      ri = -1.;
      [[fallthrough]];
    case 11:
      si = 0.;
      break;

    case 10:
      si = -1.;
      [[fallthrough]];
    case 8:
      ri = 0.;
      break;

    case 14:
      si = -1.;
      [[fallthrough]];
    case 12:
      ri = 0.;
      ti = -1.;
      break;

    case 17:
      ri = -1.;
      [[fallthrough]];
    case 16:
      ti = 0.;
      break;

    case 18:
      ri = -1.;
      [[fallthrough]];
    case 19:
      si = -1.;
      ti = 0.;
      break;
    }

    auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
    Real Phi{ 0. };
    auto t{ r0 + s0 + t0 - 2. };

    if (inod < 8) // Corner nodes
      Phi = (1 + r0) * (1 + s0) * (1 + t0) * t / 8.;

    else { // Middle nodes
      if (math::abs(ri) < tol)
        Phi = (1 - x * x) * (1 + s0) * (1 + t0) / 4.;
      else if (math::abs(si) < tol)
        Phi = (1 - y * y) * (1 + r0) * (1 + t0) / 4.;
      else if (math::abs(ti) < tol)
        Phi = (1 - z * z) * (1 + r0) * (1 + s0) / 4.;
    }
    return Phi;
  }

  static inline Real3 hexa20ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    Real tol{ 1.0e-15 };

    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 0 = (1,1,1)

    switch (inod) {
    default:
      break;

    case 5:
      ti = -1.;
    case 1:
      ri = -1;
      break;

    case 6:
      ti = -1.;
    case 2:
      ri = -1;
      si = -1;
      break;

    case 7:
      ti = -1.;
    case 3:
      si = -1;
      break;

    case 4:
      ti = -1.;
      break;

    case 9:
      ri = -1.;
    case 11:
      si = 0.;
      break;

    case 10:
      si = -1.;
    case 8:
      ri = 0.;
      break;

    case 14:
      si = -1.;
    case 12:
      ri = 0.;
      ti = -1.;
      break;

    case 17:
      ri = -1.;
      [[fallthrough]];
    case 16:
      ti = 0.;
      break;

    case 18:
      ri = -1.;
      [[fallthrough]];
    case 19:
      si = -1.;
      ti = 0.;
      break;
    }

    auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
    auto t{ r0 + s0 + t0 - 2. };
    Real3 dPhi;

    if (inod < 8) { // Corner nodes
      dPhi = hexa8ShapeFuncDeriv(inod, ref_coord);
      dPhi.x *= (t + 1. + r0);
      dPhi.y *= (t + 1. + s0);
      dPhi.z *= (t + 1. + t0);
    }
    else { // Middle nodes
      auto x2{ x * x }, y2{ y * y }, z2{ z * z };
      if (math::abs(ri) < tol) {
        dPhi.x = -x * (1 + s0) * (1 + t0) / 2.;
        dPhi.y = si * (1 - x2) * (1 + t0) / 4.;
        dPhi.z = ti * (1 - x2) * (1 + s0) / 4.;
      }
      else if (math::abs(si) < tol) {
        dPhi.x = ri * (1 - y2) * (1 + t0) / 4.;
        dPhi.y = -y * (1 + r0) * (1 + t0) / 2.;
        dPhi.z = ti * (1 - y2) * (1 + r0) / 4.;
      }
      else if (math::abs(ti) < tol) {
        dPhi.x = ri * (1 - z2) * (1 + s0) / 4.;
        dPhi.y = si * (1 - z2) * (1 + r0) / 4.;
        dPhi.z = -z * (1 + r0) * (1 + s0) / 2.;
      }
    }
    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) tetrahedral finite-element
   * The "Tetra4" reference element is assumed as follows:
   *
   *    (0,0,1)                     3
   *       .                        *.*
   *       .                        * . *
   *       .                        *  .  *
   *       .                        *   .   *
   *       Z   (0,1,0)              *    .    *
   *       .    .                   *     2     *
   *       .   .                    *   .    .    *
   *       .  Y                     *  .        .   *
   *       . .                      * .            .  *
   *       ..           (1,0,0)     *.                . *
   *       --------X------>         0********************1
   *
   * direct local numbering : 0->1->2->3
   */
  /*---------------------------------------------------------------------------*/

  static inline Real tetra4ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 4);
#endif

    auto ri = ref_coord[0], si = ref_coord[1], ti = ref_coord[2]; // default is first node (index 3)

    switch (inod) {
    default:
      break;
    case 1:
      return ri;
    case 2:
      return si;
    case 0:
      return (1. - ri - si - ti);
    }
    return ti;
  }

  static inline Real3 tetra4ShapeFuncDeriv(Integer inod, Real3 /*ref_coord*/)
  {

    if (inod == 3)
      return { 0., 0., 1. };
    if (inod == 1)
      return { 1., 0., 0. };
    if (inod == 2)
      return { 0., 1., 0. };
    return { -1., -1., -1. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference quadratic (P2) tetrahedral finite-element
   * The "Tetra10" reference element is assumed as follows:
   *
   *    (0,0,1)                     x 3
   *       .                        *.*
   *       .                        * . *
   *       .                        *  .  *
   *       .                        *  9   *
   *       Z   (0,1,0)              *    .    *
   *       .    .                   *     x     8
   *       .   .                    7   . 2  .    *
   *       .  Y                     *  6       5    *
   *       . .                      * .            .  *
   *       ..           (1,0,0)     *.                . *
   *       --------X------>       0 x ******* 4 ******** x 1
   *
   * direct local numbering : 0->1->2->3->...->8->9
   */
  /*---------------------------------------------------------------------------*/

  static inline Real tetra10ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 10);
#endif

    auto x = ref_coord[0], y = ref_coord[1], z = ref_coord[2],
         t = 1. - x - y - z,
         wi{ 0. };

    switch (inod) {
    default:
      break;

    // Corner nodes
    case 0:
      wi = t * (2 * t - 1.);
      break; //=(1. - 2*x - 2*y - 2*z) * t
    case 1:
      wi = x * (2 * x - 1.);
      break; //=(1. - 2*t - 2*y - 2*z)*x
    case 2:
      wi = y * (2 * y - 1.);
      break; //=(1. - 2*x - 2*t - 2*z)*y
    case 3:
      wi = z * (2 * z - 1.);
      break; //=(1. - 2*t - 2*x - 2*y)*z

    // Middle nodes
    case 4:
      wi = 4 * x * t;
      break;
    case 5:
      wi = 4 * x * y;
      break;
    case 6:
      wi = 4 * y * t;
      break;
    case 7:
      wi = 4 * z * t;
      break;
    case 8:
      wi = 4 * z * x;
      break;
    case 9:
      wi = 4 * z * y;
      break;
    }
    return wi;
  }

  static inline Real3 tetra10ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] },
    t{ 1. - x - y - z },
    x4{ 4 * x },
    y4{ 4 * y },
    z4{ 4 * z },
    t4{ 4 * t };

    // Corner nodes
    /*
   if (inod == 3) return {0.,0.,1. + 2*t - 2*x - 2*y + 2*z};
   if (inod == 1) return {1. - 2*t - 2*y - 2*z + 2*x,0.,0.};
   if (inod == 2) return {0.,1. - 2*x - 2*t - 2*z + 2*y,0.};
   if (!inod) return {-1. - 2*t + 2*x + 2*y + 2*z,-1. - 2*t+ 2*x + 2*y + 2*z,-1. - 2*t + 2*x + 2*y + 2*z};
*/
    if (!inod)
      return { 1. - t4, 1. - t4, 1. - t4 };
    if (inod == 1)
      return { x4 - 1., 0., 0. };
    if (inod == 2)
      return { 0., y4 - 1., 0. };
    if (inod == 3)
      return { 0., 0., z4 - 1. };

    // Middle nodes
    if (inod == 4)
      return { t4 - x4, -x4, -x4 };
    if (inod == 5)
      return { y4, x4, 0. };
    if (inod == 6)
      return { -y4, t4 - y4, -y4 };
    if (inod == 8)
      return { z4, 0., x4 };
    if (inod == 9)
      return { 0., z4, y4 };
    return { -z4, -z4, t4 - z4 }; //inod == 7
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) pentaedron finite-element
   * The "Penta6" reference element is assumed as follows:
   *
   *                     5 (0,1,1)
   *                   . |  .
   *                  .  |     .
   *                 .   Z        .
   *                .    |           .
   *               .     |             .
   *       (0,0,1) 3 ------------------ 4 (1,0,1)
   *               |     |              |
   *               |     |              |
   *               |     |              |
   *               |     |              |
   *               |     2 (0,1,-1)     |
   *               |   .    .           |
   *               |  Y        .        |
   *               | .            .     |
   *               |.                .  |
   *      (0,0,-1) 0 -------- X ------- 1 (1,0,-1)
   *
   * direct local numbering : 0->1->2->3->4->5->6
   */
  /*---------------------------------------------------------------------------*/

  static inline Real penta6ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 6);
#endif
    auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
    auto r0{ 1. }, s0{ 1. }, ti{ -1. };
    auto rs{ 1. - r - s };

    if (inod >= 3)
      ti = 1.;
    auto t0{ 1 + ti * t };

    switch (inod) {
    default:
      break; // Node 0
    case 4:
    case 1:
      r0 = r;
      rs = 1.;
      break;
    case 5:
    case 2:
      s0 = s;
      rs = 1.;
      break;
    }

    return 0.5 * r0 * s0 * rs * t0;
  }

  static inline Real3 penta6ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {

#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 6);
#endif
    auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. };
    auto r0{ 1. }, s0{ 1. }, ti{ -1. };
    auto rs{ 1. - r - s };

    if (inod >= 3)
      ti = 1.;
    auto t0{ 1 + ti * t };

    switch (inod) {
    default:
      break;
    case 3:
    case 0:
      ri = -1.;
      si = -1.;
      break;
    case 4:
    case 1:
      r0 = r;
      si = 0.;
      rs = 1.;
      break;
    case 5:
    case 2:
      s0 = s;
      rs = 1.;
      break;
    }

    Real3 dPhi;
    dPhi.x = 0.5 * ri * t0;
    dPhi.y = 0.5 * si * t0;
    dPhi.z = 0.5 * ti * rs * r0 * s0;
    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) pyramid finite-element
   * The "Pyramid5" reference element is assumed as follows:
   *
   *                               ^
   *                               |
   *                               Z
   *                               |
   *                               4 (0,0,1)
   *                              *
   *                             * **
   *                            ** *  *
   *                           * * |*   *
   *                          * *  | *    *
   *                         * *   |  *     *
   *                        * *    |   *      *      .Y
   *                       * *     |    *       *  .
   *            (-1,0,0)  * 2 -----|-----*------ 1 (0,1,0)
   *                     * .       |      *  .  .
   *                    * .        |     . *   .
   *                   *.             .     * .
   *                  *                  X . *
   *        (0,-1,0) 3 --------------------- 0 (1,0,0)
   *
   * direct local numbering : 0->1->2->3->4
   */
  /*---------------------------------------------------------------------------*/

  static inline Real pyramid5ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 5);
#endif
    Real tol{ 1.0e-15 };

    auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
    auto r1{ -1. }, s1{ 1. }, r2{ -1. }, s2{ -1. };

    if (inod == 4)
      return t;
    auto ti{ t - 1. };
    auto t0{ 0. };

    if (math::abs(ti) < tol)
      ti = 0.;
    else
      t0 = -1. / ti / 4.;

    switch (inod) {
    case 1:
      s1 = -1.;
      r2 = 1.;
      break;
    case 2:
      r1 = 1.;
      r2 = 1.;
      break;
    case 3:
      r1 = 1.;
      s2 = 1.;
      break;
    default:
      break; // default is for node 0
    }

    return (r1 * r + s1 * s + ti) * (r2 * r + s2 * s + ti) * t0;
  }

  static inline Real3 pyramid5ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {

#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 5);
#endif
    Real tol{ 1.0e-15 };

    auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
    auto r1{ -1. }, s1{ 1. }, r2{ -1. }, s2{ -1. };

    auto ti{ t - 1. };
    auto t0{ 0. };

    if (math::abs(ti) < tol)
      ti = 0.;
    else
      t0 = -1. / ti / 4.;

    switch (inod) {
    case 1:
      s1 = -1.;
      r2 = 1.;
      break;
    case 2:
      r1 = 1.;
      r2 = 1.;
      break;
    case 3:
      r1 = 1.;
      s2 = 1.;
      break;
    default:
      break; // default is for node 0
    }

    if (inod == 4)
      return { 0., 0., 1. };

    Real3 dPhi;
    auto r12{ r1 + r2 }, rr{ 2. * r1 * r2 }, s12{ s1 + s2 }, ss{ 2. * s1 * s2 }, rs{ r1 * s2 + r2 * s1 }, t02{ 4. * t0 * t0 };

    dPhi.x = t0 * (rr * r + rs * s + r12 * ti);
    dPhi.y = t0 * (rs * r + ss * s + s12 * ti);

    if (math::abs(ti) < tol)
      dPhi.z = 0.;
    else
      dPhi.z = t0 * (r12 * r + s12 * s + 2. * ti) + t02 * (r1 * r + s1 * s + ti) * (r2 * r + s2 * s + ti);

    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace ArcaneFemFunctions

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
