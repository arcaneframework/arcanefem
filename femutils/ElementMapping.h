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
//! Element tag used for Hexa8 mapping.
struct Hexa8 {};

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

/*---------------------------------------------------------------------------*/
/**
 * @brief Trilinear mapping for a three-dimensional Hexa8 element.
 *
 * The seven nonconstant coefficients of
 * p(xi,eta,zeta) = center + A*xi + B*eta + C*zeta + D*xi*eta
 *                   + E*xi*zeta + F*eta*zeta + G*xi*eta*zeta
 * are cached from node-coordinate differences. The center is unnecessary for
 * Jacobians and gradients.
 */
/*---------------------------------------------------------------------------*/
template <> class ElementMapping<Hexa8>
{
 public:
  ARCCORE_HOST_DEVICE ElementMapping(const Real3& p0, const Real3& p1,
                                    const Real3& p2, const Real3& p3,
                                    const Real3& p4, const Real3& p5,
                                    const Real3& p6, const Real3& p7)
  : m_ax(0.125 * ((p1.x - p0.x) + (p2.x - p3.x) + (p5.x - p4.x) + (p6.x - p7.x)))
  , m_ay(0.125 * ((p1.y - p0.y) + (p2.y - p3.y) + (p5.y - p4.y) + (p6.y - p7.y)))
  , m_az(0.125 * ((p1.z - p0.z) + (p2.z - p3.z) + (p5.z - p4.z) + (p6.z - p7.z)))
  , m_bx(0.125 * ((p3.x - p0.x) + (p2.x - p1.x) + (p7.x - p4.x) + (p6.x - p5.x)))
  , m_by(0.125 * ((p3.y - p0.y) + (p2.y - p1.y) + (p7.y - p4.y) + (p6.y - p5.y)))
  , m_bz(0.125 * ((p3.z - p0.z) + (p2.z - p1.z) + (p7.z - p4.z) + (p6.z - p5.z)))
  , m_cx(0.125 * ((p4.x - p0.x) + (p5.x - p1.x) + (p6.x - p2.x) + (p7.x - p3.x)))
  , m_cy(0.125 * ((p4.y - p0.y) + (p5.y - p1.y) + (p6.y - p2.y) + (p7.y - p3.y)))
  , m_cz(0.125 * ((p4.z - p0.z) + (p5.z - p1.z) + (p6.z - p2.z) + (p7.z - p3.z)))
  , m_dx(0.125 * ((p0.x - p1.x) + (p2.x - p3.x) + (p4.x - p5.x) + (p6.x - p7.x)))
  , m_dy(0.125 * ((p0.y - p1.y) + (p2.y - p3.y) + (p4.y - p5.y) + (p6.y - p7.y)))
  , m_dz(0.125 * ((p0.z - p1.z) + (p2.z - p3.z) + (p4.z - p5.z) + (p6.z - p7.z)))
  , m_ex(0.125 * ((p0.x - p1.x) + (p3.x - p2.x) + (p5.x - p4.x) + (p6.x - p7.x)))
  , m_ey(0.125 * ((p0.y - p1.y) + (p3.y - p2.y) + (p5.y - p4.y) + (p6.y - p7.y)))
  , m_ez(0.125 * ((p0.z - p1.z) + (p3.z - p2.z) + (p5.z - p4.z) + (p6.z - p7.z)))
  , m_fx(0.125 * ((p0.x - p3.x) + (p1.x - p2.x) + (p6.x - p5.x) + (p7.x - p4.x)))
  , m_fy(0.125 * ((p0.y - p3.y) + (p1.y - p2.y) + (p6.y - p5.y) + (p7.y - p4.y)))
  , m_fz(0.125 * ((p0.z - p3.z) + (p1.z - p2.z) + (p6.z - p5.z) + (p7.z - p4.z)))
  , m_gx(0.125 * ((p1.x - p0.x) + (p3.x - p2.x) + (p4.x - p5.x) + (p6.x - p7.x)))
  , m_gy(0.125 * ((p1.y - p0.y) + (p3.y - p2.y) + (p4.y - p5.y) + (p6.y - p7.y)))
  , m_gz(0.125 * ((p1.z - p0.z) + (p3.z - p2.z) + (p4.z - p5.z) + (p6.z - p7.z)))
  {}

  struct JacobianInfo
  {
    Real3x3 j;
    Real det_j;
  };

  ARCCORE_HOST_DEVICE JacobianInfo jacobianAndDeterminant(Real xi, Real eta, Real zeta) const
  {
    const Real xi_eta = xi * eta;
    const Real xi_zeta = xi * zeta;
    const Real eta_zeta = eta * zeta;
    Real3x3 j;
    j[0][0] = m_ax + eta * m_dx + zeta * m_ex + eta_zeta * m_gx;
    j[0][1] = m_ay + eta * m_dy + zeta * m_ey + eta_zeta * m_gy;
    j[0][2] = m_az + eta * m_dz + zeta * m_ez + eta_zeta * m_gz;
    j[1][0] = m_bx + xi * m_dx + zeta * m_fx + xi_zeta * m_gx;
    j[1][1] = m_by + xi * m_dy + zeta * m_fy + xi_zeta * m_gy;
    j[1][2] = m_bz + xi * m_dz + zeta * m_fz + xi_zeta * m_gz;
    j[2][0] = m_cx + xi * m_ex + eta * m_fx + xi_eta * m_gx;
    j[2][1] = m_cy + xi * m_ey + eta * m_fy + xi_eta * m_gy;
    j[2][2] = m_cz + xi * m_ez + eta * m_fz + xi_eta * m_gz;
    return { j, math::matrixDeterminant(j) };
  }

  struct Gradients
  {
    RealVector<8> dN_dx;
    RealVector<8> dN_dy;
    RealVector<8> dN_dz;
  };

  // The caller can validate det_j before inversion (the CPU path does so).
  ARCCORE_HOST_DEVICE static Gradients gradients(const JacobianInfo& jacobian_info,
                                                Real xi, Real eta, Real zeta)
  {
    const auto ref = ShapeFunctions::computeReferenceGradientsHexa8(xi, eta, zeta);
    const Real3x3 inv_j = math::inverseMatrix(jacobian_info.j, jacobian_info.det_j);
    Gradients result;
    for (int a = 0; a < 8; ++a) {
      result.dN_dx[a] = inv_j[0][0] * ref.dN_dxi[a] + inv_j[0][1] * ref.dN_deta[a] + inv_j[0][2] * ref.dN_dzeta[a];
      result.dN_dy[a] = inv_j[1][0] * ref.dN_dxi[a] + inv_j[1][1] * ref.dN_deta[a] + inv_j[1][2] * ref.dN_dzeta[a];
      result.dN_dz[a] = inv_j[2][0] * ref.dN_dxi[a] + inv_j[2][1] * ref.dN_deta[a] + inv_j[2][2] * ref.dN_dzeta[a];
    }
    return result;
  }

 private:
  Real m_ax, m_ay, m_az;
  Real m_bx, m_by, m_bz;
  Real m_cx, m_cy, m_cz;
  Real m_dx, m_dy, m_dz;
  Real m_ex, m_ey, m_ez;
  Real m_fx, m_fy, m_fz;
  Real m_gx, m_gy, m_gz;
};

} // namespace Arcane::FemUtils
#endif
