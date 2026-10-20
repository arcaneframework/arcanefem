// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMapping.h                                            (C) 2022-2026 */
/*                                                                           */
/* Finite element mapping for various element types.                         */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_ELEMENTMAPPING_H
#define ARCANEFEM_FEMUTILS_ELEMENTMAPPING_H

#include "ShapeFunctions.h"

namespace Arcane::FemUtils
{

//! Element tag used for Quad4 mapping.
struct Quad4 {};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

//! Template for element mappings for each supported element type.
template <typename Element> class ElementMapping;

/*---------------------------------------------------------------------------*/
/**
 * @brief Quadrilateral element mapping for 2D finite elements (Quad4).
 *
 * This class maps natural coordinates (ξ,η) to physical coordinates (x,y)
 * and computes the Jacobian matrix and gradients of shape functions.
 *
 * The mapping is defined as:
 *   p(ξ, η) = center + A*ξ + B*η + C*ξ*η
 *
 * where A, B, C are computed from the physical coordinates of the four nodes.
 */
/*---------------------------------------------------------------------------*/
template <> class ElementMapping<Quad4>
{
 public:
  ARCCORE_HOST_DEVICE ElementMapping(const Real3& p0, const Real3& p1,
                                    const Real3& p2, const Real3& p3)
  : m_ax(0.25 * ((p1.x - p0.x) + (p2.x - p3.x)))
  , m_ay(0.25 * ((p1.y - p0.y) + (p2.y - p3.y)))
  , m_bx(0.25 * ((p3.x - p0.x) + (p2.x - p1.x)))
  , m_by(0.25 * ((p3.y - p0.y) + (p2.y - p1.y)))
  , m_cx(0.25 * ((p0.x - p1.x) + (p2.x - p3.x)))
  , m_cy(0.25 * ((p0.y - p1.y) + (p2.y - p3.y)))
  {}

  struct JacobianInfo
  {
    Real2x2 j;
    Real det_j;
  };

  ARCCORE_HOST_DEVICE JacobianInfo jacobianAndDeterminant(Real xi, Real eta) const
  {
    // p(xi, eta) = center + A*xi + B*eta + C*xi*eta.
    // Coordinate differences avoid cancellation from a large origin.
    Real2x2 j;
    j[0][0] = m_ax + eta * m_cx;
    j[0][1] = m_ay + eta * m_cy;
    j[1][0] = m_bx + xi * m_cx;
    j[1][1] = m_by + xi * m_cy;
    return { j, j[0][0] * j[1][1] - j[0][1] * j[1][0] };
  }

  struct Gradients
  {
    RealVector<4> dN_dx;
    RealVector<4> dN_dy;
  };

  // The caller can validate det_j before inversion (the CPU path does so).
  ARCCORE_HOST_DEVICE static Gradients gradients(const JacobianInfo& jacobian_info,
                                                Real xi, Real eta)
  {
    const auto ref = ShapeFunctions::computeReferenceGradientsQuad4(xi, eta);
    const Real2x2& j = jacobian_info.j;
    const Real inv_det = 1.0 / jacobian_info.det_j;
    const Real inv00 = j[1][1] * inv_det;
    const Real inv01 = -j[0][1] * inv_det;
    const Real inv10 = -j[1][0] * inv_det;
    const Real inv11 = j[0][0] * inv_det;
    Gradients result;
    for (int a = 0; a < 4; ++a) {
      result.dN_dx[a] = inv00 * ref.dN_dxi[a] + inv01 * ref.dN_deta[a];
      result.dN_dy[a] = inv10 * ref.dN_dxi[a] + inv11 * ref.dN_deta[a];
    }
    return result;
  }

 private:
  Real m_ax, m_ay;
  Real m_bx, m_by;
  Real m_cx, m_cy;
};

} // namespace Arcane::FemUtils
#endif
