// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FeOperation.h                                               (C) 2000-2026 */
/*                                                                           */
/* Various 2D and 3D finite element operations on mesh items.                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMTUTILS_FEOPERATION_H
#define ARCANEFEM_FEMTUTILS_FEOPERATION_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/VariableTypes.h>

#include "FemUtils.h"
#include "ShapeFunctions.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace ArcaneFemFunctions
{

using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
  * @brief Provides methods for finite element operations in 2D.
  *
  * This class includes static methods for calculating gradients of basis
  * functions and integrals for P1 triangles in 2D finite element analysis.
  *
  * Reference Tri3 and Quad4 element in Arcane
  *
  *               0 o                    1 o . . . . o 0
  *                . .                     .         .
  *               .   .                    .         .
  *              .     .                   .         .
  *           1 o  . .  o 2              2 o . . . . o 3
  */
/*---------------------------------------------------------------------------*/
class FeOperation2D
{
 public:

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the gradients of given scalar function U for P1 triangles.
   *
   * This method calculates gradient operator ∇ Ui for a given P1 cell
   * with i = 1,..,3 for the three scalar values of Ui hence at cell nodes.
   * The output is ∇ Ui is P0 (piece-wise constant) hence Real3 value
   * per cell
   *
   *         ∇ Ui = [ ∂U/∂𝑥   ∂U/∂𝑦   ∂U/∂𝑧 ]
   *
   * @param cell The Tria3 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param u The variable Ui defined at the nodes.
   *
   * @return A Real3 vector of the gradient ∇Ui = {∂U/∂𝑥, ∂U/∂𝑦, 0} at the cell center.
   *
   * @note we can adapt the same for 3D by filling the third component
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3
  computeGradientTria3(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal& u)
  {
    Real3 n0 = node_coord[cell.nodeId(0)];
    Real3 n1 = node_coord[cell.nodeId(1)];
    Real3 n2 = node_coord[cell.nodeId(2)];

    Real u0 = u[cell.nodeId(0)];
    Real u1 = u[cell.nodeId(1)];
    Real u2 = u[cell.nodeId(2)];

    Real A2 = ((n1.x - n0.x) * (n2.y - n0.y) - (n2.x - n0.x) * (n1.y - n0.y));

    return { (u0 * (n1.y - n2.y) + u1 * (n2.y - n0.y) + u2 * (n0.y - n1.y)) / A2, (u0 * (n2.x - n1.x) + u1 * (n0.x - n2.x) + u2 * (n1.x - n0.x)) / A2, 0 };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the gradients of given vector function U for P1 triangles.
   *
   * This method calculates gradient operator ∇ Ui for a given P1 cell
   * with i = 1,..,3 for the three vector values of Ui hence at cell nodes.
   * The output is ∇ Ui is P0 (piece-wise constant) hence Real3x3 value
   * per cell
   *
   *         ∇ Ui = [ ∂U/∂𝑥   ∂U/∂𝑦   ∂U/∂𝑧 ]
   *
   * @param cell The Tria3 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param u The variable Ui defined at the nodes.
   *
   * @return A Real3x3 vector of the gradient ∇Ui = {∂U/∂𝑥, ∂U/∂𝑦, 0} at the cell center.
   *
   * @note we can adapt the same for 3D by filling the third component
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3x3
  computeGradientTria3(Cell cell, const VariableNodeReal3& node_coord,
                       const VariableNodeReal3& u)
  {
    Real3 n0 = node_coord[cell.nodeId(0)];
    Real3 n1 = node_coord[cell.nodeId(1)];
    Real3 n2 = node_coord[cell.nodeId(2)];

    Real3 u0 = u[cell.nodeId(0)];
    Real3 u1 = u[cell.nodeId(1)];
    Real3 u2 = u[cell.nodeId(2)];

    Real A2 = ((n1.x - n0.x) * (n2.y - n0.y) - (n2.x - n0.x) * (n1.y - n0.y));

    Real3 d_ux = { (u0.x * (n1.y - n2.y) + u1.x * (n2.y - n0.y) + u2.x * (n0.y - n1.y)) / A2, (u0.x * (n2.x - n1.x) + u1.x * (n0.x - n2.x) + u2.x * (n1.x - n0.x)) / A2, 0 };
    Real3 d_uy = { (u0.y * (n1.y - n2.y) + u1.y * (n2.y - n0.y) + u2.y * (n0.y - n1.y)) / A2, (u0.y * (n2.x - n1.x) + u1.y * (n0.x - n2.x) + u2.y * (n1.x - n0.x)) / A2, 0 };
    Real3 d_uz = { 0., 0., 0. };
    Real3x3 grad = { d_ux, d_uy, d_uz };

    return grad;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the 𝑥 gradients of basis functions 𝐍 for ℙ1 triangles.
   *
   * This method calculates gradient operator ∂/∂𝑥 of 𝑁ᵢ for a given ℙ1
   * cell with i = 1,..,3 for the three shape function  𝑁ᵢ  hence output
   * is a vector of size 3
   *
   *         ∂𝐍/∂𝑥 = [ ∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ∂𝑁₃/∂𝑥 ]
   *
   *         ∂𝐍/∂𝑥 = 1/(2𝐴) [ 𝑦₂-𝑦₃  𝑦₃-𝑦₁  𝑦₁-𝑦₂ ]
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 computeGradientXTria3(Cell cell, const VariableNodeReal3& node_coord)
  {
    Real3 vertex0 = node_coord[cell.nodeId(0)];
    Real3 vertex1 = node_coord[cell.nodeId(1)];
    Real3 vertex2 = node_coord[cell.nodeId(2)];

    auto A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

    return { (vertex1.y - vertex2.y) / A2, (vertex2.y - vertex0.y) / A2, (vertex0.y - vertex1.y) / A2 };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the 𝑦 gradients of basis functions 𝐍 for ℙ1 triangles.
   *
   * This method calculates gradient operator ∂/∂𝑦 of 𝑁ᵢ for a given ℙ1
   * cell with i = 1,..,3 for the three shape function  𝑁ᵢ  hence output
   * is a vector of size 3
   *
   *         ∂𝐍/∂𝑦 = [ ∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ∂𝑁₃/∂𝑦 ]
   *
   *         ∂𝐍/∂𝑦 = 1/(2𝐴) [ 𝑥₃−𝑥₂  𝑥₁−𝑥₃  𝑥₂−𝑥₁ ]
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 computeGradientYTria3(Cell cell, VariableNodeReal3 node_coord)
  {
    Real3 vertex0 = node_coord[cell.nodeId(0)];
    Real3 vertex1 = node_coord[cell.nodeId(1)];
    Real3 vertex2 = node_coord[cell.nodeId(2)];

    auto A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

    return { (vertex2.x - vertex1.x) / A2, (vertex0.x - vertex2.x) / A2, (vertex1.x - vertex0.x) / A2 };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Holds information for a quadrilateral element at a single Gauss point.
   *
   * This includes the gradients of the shape functions in the physical space (x, y)
   * and the determinant of the Jacobian matrix.
   */
  /*---------------------------------------------------------------------------*/
  template <Int32 n> struct QuadGaussPointInfo
  {
    RealVector<n> dN_dx; // Derivatives of shape functions in x {∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ...  ∂𝑁ₙ/∂𝑥 }
    RealVector<n> dN_dy; // Derivatives of shape functions in y {∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ...  ∂𝑁ₙ/∂𝑦}
    Real det_j; // Determinant of the Jacobian matrix at the Gauss point.
  };

  using Quad4GaussPointInfo = QuadGaussPointInfo<4>;
  using Quad8GaussPointInfo = QuadGaussPointInfo<8>;
  using Quad9GaussPointInfo = QuadGaussPointInfo<9>;

  template <Int32 N> static inline QuadGaussPointInfo<N>
  _computeQuadGradientsAndJacobian(Cell cell,
                                   const VariableNodeReal3& node_coord,
                                   const Arcane::FemUtils::ShapeFunctions::ReferenceGradients2D<N>& reference_gradients)
  {
    // Jacobian calculation 𝑱
    //    𝑱 = [ 𝒋₀₀  𝒋₀₁ ] = [ ∂𝑥/∂ξ  ∂𝑦/∂ξ ]
    //        [ 𝒋₁₀  𝒋₁₁ ]   [ ∂𝑥/∂η  ∂𝑦/∂η ]
    Real2x2 J;
    J[0][0] = 0.0;
    J[0][1] = 0.0;
    J[1][0] = 0.0;
    J[1][1] = 0.0;

    for (Int8 a = 0; a < N; ++a) {
      const auto& coord = node_coord[cell.nodeId(a)];
      J[0][0] += reference_gradients.dN_dxi[a] * coord.x;
      J[0][1] += reference_gradients.dN_dxi[a] * coord.y;
      J[1][0] += reference_gradients.dN_deta[a] * coord.x;
      J[1][1] += reference_gradients.dN_deta[a] * coord.y;
    }

    const Real detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    if (detJ <= 0.0) {
      ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
    }

    const Real invJ00 = J[1][1] / detJ;
    const Real invJ01 = -J[0][1] / detJ;
    const Real invJ10 = -J[1][0] / detJ;
    const Real invJ11 = J[0][0] / detJ;

    RealVector<N> dN_dx_result;
    RealVector<N> dN_dy_result;

    for (Int8 a = 0; a < N; ++a) {
      dN_dx_result(a) = invJ00 * reference_gradients.dN_dxi[a] + invJ01 * reference_gradients.dN_deta[a];
      dN_dy_result(a) = invJ10 * reference_gradients.dN_dxi[a] + invJ11 * reference_gradients.dN_deta[a];
    }

    return { dN_dx_result, dN_dy_result, detJ };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes shape function gradients and the Jacobian determinant for a Quad8 element.
   * https://www.meil.pw.edu.pl/content/download/50088/264514/file/FEM_1_9_8node_2D.pdf
   *
   *   3 ---- 6 ---- 2
   *   |             |
   *   7             5
   *   |             |
   *   0 ---- 4 ---- 1
   *
   *   Reference coordinates:
   *  0: (-1,-1), 1: (1,-1), 2: (1,1), 3: (-1,1)
   *  4: (0,-1),  5: (1,0),  6: (0,1), 7: (-1,0)
   *
   * @param cell The Quad8 (serendipity) cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param xi The ξ coordinate of the evaluation point (-1 to 1).
   * @param eta The η coordinate of the evaluation point (-1 to 1).
   * @return A Quad8GaussPointInfo struct containing {∂𝐍/∂𝑥, ∂𝐍/∂𝑦, det(𝑱)}.
   */
  /*---------------------------------------------------------------------------*/

  static inline Quad8GaussPointInfo
  computeGradientsAndJacobianQuad8(Cell cell, const VariableNodeReal3& node_coord, Real xi, Real eta)
  {
    const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad8(xi, eta);
    return _computeQuadGradientsAndJacobian<8>(cell, node_coord, reference_gradients);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes shape function gradients and the Jacobian determinant for a Quad9 element.
   * https://www.meil.pw.edu.pl/content/download/50088/264514/file/FEM_1_9_8node_2D.pdf
   *
   *   3 ---- 6 ---- 2
   *   |             |
   *   7      8      5
   *   |             |
   *   0 ---- 4 ---- 1
   *
   *   Reference coordinates:
   *  0: (-1,-1), 1: (1,-1), 2: (1,1), 3: (-1,1)
   *  4: (0,-1),  5: (1,0),  6: (0,1), 7: (-1,0), 8: (0,0)
   *
   * @param cell The Quad9 (Lagrange) cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param xi The ξ coordinate of the evaluation point (-1 to 1).
   * @param eta The η coordinate of the evaluation point (-1 to 1).
   * @return A Quad9GaussPointInfo struct containing {∂𝐍/∂𝑥, ∂𝐍/∂𝑦, det(𝑱)}.
   */
  /*---------------------------------------------------------------------------*/

  static inline Quad9GaussPointInfo
  computeGradientsAndJacobianQuad9(Cell cell, const VariableNodeReal3& node_coord, Real xi, Real eta)
  {
    const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad9(xi, eta);
    return _computeQuadGradientsAndJacobian<9>(cell, node_coord, reference_gradients);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes shape function gradients and the Jacobian determinant for a Quad4 element.
   *
   * @param cell The Quad4 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param xi The ξ coordinate of the evaluation point (-1 to 1).
   * @param eta The η coordinate of the evaluation point (-1 to 1).
   * @return A Quad4GaussPointInfo struct containing {∂𝐍/∂𝑥, ∂𝐍/∂𝑦, det(𝑱)}.
   */
  /*---------------------------------------------------------------------------*/
  static inline Quad4GaussPointInfo
  computeGradientsAndJacobianQuad4(Cell cell, const VariableNodeReal3& node_coord, Real xi, Real eta)
  {
    const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad4(xi, eta);
    return _computeQuadGradientsAndJacobian<4>(cell, node_coord, reference_gradients);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the gradient of a scalar field 'u' for a Quad4 element.
   *
   * This method calculates gradient operator ∇ Ui for a given P1 cell
   * with i = 1,..,4 for the four scalar values of Ui hence at cell nodes.
   * The output is ∇ Ui is P0 (piece-wise constant) hence Real3 value
   * per cell
   *
   *         ∇ Ui = [ ∂U/∂𝑥   ∂U/∂𝑦   ∂U/∂𝑧 ]
   *
   * @param cell The Quad4 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param u The variable Ui defined at the nodes.
   *
   * @return A Real3 vector of the gradient ∇Ui = {∂U/∂𝑥, ∂U/∂𝑦, 0} at the cell center.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3
  computeGradientQuad4(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal& u /*, Real xi, Real eta*/)
  {
    // get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian at (ξ,η) = (0,0)
    const auto gp_util = computeGradientsAndJacobianQuad4(cell, node_coord, 0.0, 0.0);
    const RealVector<4>& dN_dx = gp_util.dN_dx;
    const RealVector<4>& dN_dy = gp_util.dN_dy;

    // get the nodal values of the variable 𝑢ᵢ ∀ 𝑖= 1,……,4 for the cell.
    const Real u_nodes[4] = {
      u[cell.nodeId(0)],
      u[cell.nodeId(1)],
      u[cell.nodeId(2)],
      u[cell.nodeId(3)]
    };

    // Compute the gradient components using shape function
    //    ∂𝑢/∂𝑥 = Σ (∂𝑁ᵢ/∂𝑥 * uᵢ) ∀ 𝑖= 1,……,4
    //    ∂𝑢/∂𝑦 = Σ (∂𝑁ᵢ/∂𝑦 * uᵢ) ∀ 𝑖= 1,……,4
    Real grad_x = 0.0;
    Real grad_y = 0.0;
    for (Int8 a = 0; a < 4; ++a) {
      grad_x += dN_dx(a) * u_nodes[a];
      grad_y += dN_dy(a) * u_nodes[a];
    }

    return { grad_x, grad_y, 0.0 };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the gradient of a vector field 'u' for a Quad4 element.
   *
   * This method calculates gradient operator ∇ Ui for a given P1 cell
   * with i = 1,..,4 for the four vector values of Ui hence at cell nodes.
   * The output is ∇ Ui is P0 (piece-wise constant) hence Real3x3 value
   * per cell
   *
   *         ∇ Ui = [ ∂U/∂𝑥   ∂U/∂𝑦   ∂U/∂𝑧 ]
   *
   * @param cell The Quad4 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param u The variable Ui defined at the nodes.
   * @param xi The first coordinate of quadrature point.
   * @param eta The second coordinate of quadrature point.
   *
   * @return A Real3x3 vector of the gradient ∇Ui = {∂U/∂𝑥, ∂U/∂𝑦, 0} at the cell center.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3x3
  computeGradientQuad4(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal3& u, const Real xi, const Real eta)
  {
    // get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian at (ξ,η)
    const auto gp_util = computeGradientsAndJacobianQuad4(cell, node_coord, xi, eta);
    const RealVector<4>& dN_dx = gp_util.dN_dx;
    const RealVector<4>& dN_dy = gp_util.dN_dy;

    // get the nodal values of the variable 𝑢ᵢ ∀ 𝑖= 1,……,4 for the cell.
    const Real3 u_nodes[4] = {
      u[cell.nodeId(0)],
      u[cell.nodeId(1)],
      u[cell.nodeId(2)],
      u[cell.nodeId(3)]
    };

    // Compute the gradient components using shape function for each component (ₖ) of the vector field
    //    ∂𝑢ₖ/∂𝑥 = Σ (∂𝑁ᵢ/∂𝑥 * uᵢ) ∀ 𝑖= 1,……,4
    //    ∂𝑢ₖ/∂𝑦 = Σ (∂𝑁ᵢ/∂𝑦 * uᵢ) ∀ 𝑖= 1,……,4
    Real3 d_ux = { 0., 0., 0. };
    Real3 d_uy = { 0., 0., 0. };
    Real3 d_uz = { 0., 0., 0. };

    for (Int8 a = 0; a < 4; ++a) {
      d_ux[0] += dN_dx(a) * u_nodes[a].x;
      d_ux[1] += dN_dy(a) * u_nodes[a].x;
      d_uy[0] += dN_dx(a) * u_nodes[a].y;
      d_uy[1] += dN_dy(a) * u_nodes[a].y;
    }

    return { d_ux, d_uy, d_uz };
  }
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods for finite element operations in 3D.
 *
 * This class includes static methods for calculating gradients of basis
 * functions and integrals for P1 triangles in 3D finite element analysis.
 *
 * Reference Tetra4 element in Arcane
 *
 *               3 o
 *                /|  \
 *               / |   \
 *              /  |    \
 *             /   o 2   \
 *            / .    .                          \
 *         0 o-----------o 1
 */
/*---------------------------------------------------------------------------*/
class FeOperation3D
{
 public:

  /*-------------------------------------------------------------------------*/
  /**
   * @brief Computes the X gradients of basis functions 𝐍 for ℙ1 Tetrahedron.
   *
   * This method calculates gradient operator ∂/∂𝑥 of 𝑁ᵢ for a given ℙ1
   * cell with i = 1,..,4 for the four shape function  𝑁ᵢ  hence output
   * is a vector of size 4
   *
   *         ∂𝐍/∂𝑥 = [ ∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ∂𝑁₃/∂𝑥  ∂𝑁₄/∂𝑥 ]
   *
   *         ∂𝐍/∂𝑥 = 1/(6𝑉) [ dx0  dx1  dx2  dx3 ]
   *
   * where:
   *    dx0 = (n1.y * (n3.z - n2.z) + n2.y * (n1.z - n3.z) + n3.y * (n2.z - n1.z)),
   *    dx1 = (n0.y * (n2.z - n3.z) + n2.y * (n3.z - n0.z) + n3.y * (n0.z - n2.z)),
   *    dx2 = (n0.y * (n3.z - n1.z) + n1.y * (n0.z - n3.z) + n3.y * (n1.z - n0.z)),
   *    dx3 = (n0.y * (n1.z - n2.z) + n1.y * (n2.z - n0.z) + n2.y * (n0.z - n1.z)).
   *
   */
  /*-------------------------------------------------------------------------*/

  static inline Real4 computeGradientXTetra4(Cell cell, const VariableNodeReal3& node_coord)
  {
    Real3 n0 = node_coord[cell.nodeId(0)];
    Real3 n1 = node_coord[cell.nodeId(1)];
    Real3 n2 = node_coord[cell.nodeId(2)];
    Real3 n3 = node_coord[cell.nodeId(3)];

    Real3 v0 = n1 - n0;
    Real3 v1 = n2 - n0;
    Real3 v2 = n3 - n0;

    // 6 x Volume of tetrahedron
    Real V6 = math::abs(math::dot(v0, math::cross(v1, v2)));

    Real4 dx{};

    dx[0] = (n1.y * (n3.z - n2.z) + n2.y * (n1.z - n3.z) + n3.y * (n2.z - n1.z)) / V6;
    dx[1] = (n0.y * (n2.z - n3.z) + n2.y * (n3.z - n0.z) + n3.y * (n0.z - n2.z)) / V6;
    dx[2] = (n0.y * (n3.z - n1.z) + n1.y * (n0.z - n3.z) + n3.y * (n1.z - n0.z)) / V6;
    dx[3] = (n0.y * (n1.z - n2.z) + n1.y * (n2.z - n0.z) + n2.y * (n0.z - n1.z)) / V6;

    return dx;
  }

  /*-------------------------------------------------------------------------*/
  /**
   * @brief Computes the Y gradients of basis functions 𝐍 for ℙ1 Tetrahedron.
   *
   * This method calculates gradient operator ∂/∂𝑦 of 𝑁ᵢ for a given ℙ1
   * cell with i = 1,..,4 for the four shape functions 𝑁ᵢ, hence the output
   * is a vector of size 4.
   *
   *         ∂𝐍/∂𝑦 = [ ∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ∂𝑁₃/∂𝑦  ∂𝑁₄/∂𝑦 ]
   *
   *         ∂𝐍/∂𝑦 = 1/(6𝑉) [ dy0  dy1  dy2  dy3 ]
   *
   * where:
   *    dy0 = (n1.z * (n3.x - n2.x) + n2.z * (n1.x - n3.x) + n3.z * (n2.x - n1.x)),
   *    dy1 = (n0.z * (n2.x - n3.x) + n2.z * (n3.x - n0.x) + n3.z * (n0.x - n2.x)),
   *    dy2 = (n0.z * (n3.x - n1.x) + n1.z * (n0.x - n3.x) + n3.z * (n1.x - n0.x)),
   *    dy3 = (n0.z * (n1.x - n2.x) + n1.z * (n2.x - n0.x) + n2.z * (n0.x - n1.x)).
   *
   */
  /*-------------------------------------------------------------------------*/

  static inline Real4 computeGradientYTetra4(Cell cell, const VariableNodeReal3& node_coord)
  {
    Real3 n0 = node_coord[cell.nodeId(0)];
    Real3 n1 = node_coord[cell.nodeId(1)];
    Real3 n2 = node_coord[cell.nodeId(2)];
    Real3 n3 = node_coord[cell.nodeId(3)];

    Real3 v0 = n1 - n0;
    Real3 v1 = n2 - n0;
    Real3 v2 = n3 - n0;

    // 6 x Volume of tetrahedron
    Real V6 = math::abs(math::dot(v0, math::cross(v1, v2)));

    Real4 dy{};

    dy[0] = (n1.z * (n3.x - n2.x) + n2.z * (n1.x - n3.x) + n3.z * (n2.x - n1.x)) / V6;
    dy[1] = (n0.z * (n2.x - n3.x) + n2.z * (n3.x - n0.x) + n3.z * (n0.x - n2.x)) / V6;
    dy[2] = (n0.z * (n3.x - n1.x) + n1.z * (n0.x - n3.x) + n3.z * (n1.x - n0.x)) / V6;
    dy[3] = (n0.z * (n1.x - n2.x) + n1.z * (n2.x - n0.x) + n2.z * (n0.x - n1.x)) / V6;

    return dy;
  }

  /*-------------------------------------------------------------------------*/
  /**
   * @brief Computes the 𝑧 gradients of basis functions 𝐍 for ℙ1 Tetrahedron.
   *
   * This method calculates gradient operator ∂/∂𝑧 of 𝑁ᵢ for a given ℙ1
   * cell with i = 1,..,4 for the four shape functions 𝑁ᵢ, hence the output
   * is a vector of size 4.
   *
   *         ∂𝐍/∂𝑧 = [ ∂𝑁₁/∂𝑧  ∂𝑁₂/∂𝑧  ∂𝑁₃/∂𝑧  ∂𝑁₄/∂𝑧 ]
   *
   *         ∂𝐍/∂𝑧 = 1/(6𝑉) [ dz0  dz1  dz2  dz3 ]
   *
   * where:
   *    dz0 = (n1.x * (n3.y - n2.y) + n2.x * (n1.y - n3.y) + n3.x * (n2.y - n1.y)),
   *    dz1 = (n0.x * (n2.y - n3.y) + n2.x * (n3.y - n0.y) + n3.x * (n0.y - n2.y)),
   *    dz2 = (n0.x * (n3.y - n1.y) + n1.x * (n0.y - n3.y) + n3.x * (n1.y - n0.y)),
   *    dz3 = (n0.x * (n1.y - n2.y) + n1.x * (n2.y - n0.y) + n2.x * (n0.y - n1.y)).
   *
   */
  /*-------------------------------------------------------------------------*/

  static inline Real4 computeGradientZTetra4(Cell cell, const VariableNodeReal3& node_coord)
  {
    Real3 n0 = node_coord[cell.nodeId(0)];
    Real3 n1 = node_coord[cell.nodeId(1)];
    Real3 n2 = node_coord[cell.nodeId(2)];
    Real3 n3 = node_coord[cell.nodeId(3)];

    auto v0 = n1 - n0;
    auto v1 = n2 - n0;
    auto v2 = n3 - n0;

    // 6 x Volume of tetrahedron
    Real V6 = math::abs(math::dot(v0, math::cross(v1, v2)));

    Real4 dz{};

    dz[0] = (n1.x * (n3.y - n2.y) + n2.x * (n1.y - n3.y) + n3.x * (n2.y - n1.y)) / V6;
    dz[1] = (n0.x * (n2.y - n3.y) + n2.x * (n3.y - n0.y) + n3.x * (n0.y - n2.y)) / V6;
    dz[2] = (n0.x * (n3.y - n1.y) + n1.x * (n0.y - n3.y) + n3.x * (n1.y - n0.y)) / V6;
    dz[3] = (n0.x * (n1.y - n2.y) + n1.x * (n2.y - n0.y) + n2.x * (n0.y - n1.y)) / V6;

    return dz;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the gradients of given scalar function U for P1 tetrahedrals.
   *
   * This method calculates gradient operator ∇ Ui for a given P1 cell
   * with i = 1,..,4 for the four scalar values of Ui hence at cell nodes.
   * The output is ∇ Ui is P0 (piece-wise constant) hence Real3 value
   * per cell
   *
   *         ∇ Ui = [ ∂U/∂𝑥   ∂U/∂𝑦   ∂U/∂𝑧 ]
   *
   * @param cell The Tetra4 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param u The variable Ui defined at the nodes.
   *
   * @return A Real3 vector of the gradient ∇Ui = {∂U/∂𝑥, ∂U/∂𝑦, ∂U/∂𝑧}
   * at the cell center.
   *
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 computeGradientTetra4(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal& u)
  {
    Real3 m0 = node_coord[cell.nodeId(0)];
    Real3 m1 = node_coord[cell.nodeId(1)];
    Real3 m2 = node_coord[cell.nodeId(2)];
    Real3 m3 = node_coord[cell.nodeId(3)];

    Real f0 = u[cell.nodeId(0)];
    Real f1 = u[cell.nodeId(1)];
    Real f2 = u[cell.nodeId(2)];
    Real f3 = u[cell.nodeId(3)];

    Real3 v0 = m1 - m0;
    Real3 v1 = m2 - m0;
    Real3 v2 = m3 - m0;

    // 6 x Volume of tetrahedron
    Real V6 = math::abs(math::dot(v0, math::cross(v1, v2)));

    // Compute gradient components
    Real3 grad;
    grad.x = (f0 * (m1.y * m2.z + m2.y * m3.z + m3.y * m1.z - m3.y * m2.z - m2.y * m1.z - m1.y * m3.z) - f1 * (m0.y * m2.z + m2.y * m3.z + m3.y * m0.z - m3.y * m2.z - m2.y * m0.z - m0.y * m3.z) + f2 * (m0.y * m1.z + m1.y * m3.z + m3.y * m0.z - m3.y * m1.z - m1.y * m0.z - m0.y * m3.z) - f3 * (m0.y * m1.z + m1.y * m2.z + m2.y * m0.z - m2.y * m1.z - m1.y * m0.z - m0.y * m2.z)) / V6;
    grad.y = (f0 * (m1.z * m2.x + m2.z * m3.x + m3.z * m1.x - m3.z * m2.x - m2.z * m1.x - m1.z * m3.x) - f1 * (m0.z * m2.x + m2.z * m3.x + m3.z * m0.x - m3.z * m2.x - m2.z * m0.x - m0.z * m3.x) + f2 * (m0.z * m1.x + m1.z * m3.x + m3.z * m0.x - m3.z * m1.x - m1.z * m0.x - m0.z * m3.x) - f3 * (m0.z * m1.x + m1.z * m2.x + m2.z * m0.x - m2.z * m1.x - m1.z * m0.x - m0.z * m2.x)) / V6;
    grad.z = (f0 * (m1.x * m2.y + m2.x * m3.y + m3.x * m1.y - m3.x * m2.y - m2.x * m1.y - m1.x * m3.y) - f1 * (m0.x * m2.y + m2.x * m3.y + m3.x * m0.y - m3.x * m2.y - m2.x * m0.y - m0.x * m3.y) + f2 * (m0.x * m1.y + m1.x * m3.y + m3.x * m0.y - m3.x * m1.y - m1.x * m0.y - m0.x * m3.y) - f3 * (m0.x * m1.y + m1.x * m2.y + m2.x * m0.y - m2.x * m1.y - m1.x * m0.y - m0.x * m2.y)) / V6;
    return grad;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the gradients of given vector function U for P1 tetrahedrals.
   *
   * This method calculates gradient operator ∇ Ui for a given P1 cell
   * with i = 1,..,4 for the four vector values of Ui hence at cell nodes.
   * The output is ∇ Ui is P0 (piece-wise constant) hence Real3x3 value
   * per cell
   *
   *         ∇ Ui = [ ∂U/∂𝑥   ∂U/∂𝑦   ∂U/∂𝑧 ]
   *
   * @param cell The Tetra4 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param u The variable Ui defined at the nodes.
   *
   * @return A Real3x3 vector of the gradient ∇Ui = {∂U/∂𝑥, ∂U/∂𝑦, ∂U/∂𝑧}
   * at the cell center.
   *
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3x3 computeGradientTetra4(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal3& u)
  {
    Real3 m0 = node_coord[cell.nodeId(0)];
    Real3 m1 = node_coord[cell.nodeId(1)];
    Real3 m2 = node_coord[cell.nodeId(2)];
    Real3 m3 = node_coord[cell.nodeId(3)];

    Real3 f0 = u[cell.nodeId(0)];
    Real3 f1 = u[cell.nodeId(1)];
    Real3 f2 = u[cell.nodeId(2)];
    Real3 f3 = u[cell.nodeId(3)];

    Real3 v0 = m1 - m0;
    Real3 v1 = m2 - m0;
    Real3 v2 = m3 - m0;

    // 6 x Volume of tetrahedron
    Real V6 = math::abs(math::dot(v0, math::cross(v1, v2)));

    // Compute gradient components

    Real3 d_ux = { (f0.x * (m1.y * m2.z + m2.y * m3.z + m3.y * m1.z - m3.y * m2.z - m2.y * m1.z - m1.y * m3.z) - f1.x * (m0.y * m2.z + m2.y * m3.z + m3.y * m0.z - m3.y * m2.z - m2.y * m0.z - m0.y * m3.z) + f2.x * (m0.y * m1.z + m1.y * m3.z + m3.y * m0.z - m3.y * m1.z - m1.y * m0.z - m0.y * m3.z) - f3.x * (m0.y * m1.z + m1.y * m2.z + m2.y * m0.z - m2.y * m1.z - m1.y * m0.z - m0.y * m2.z)) / V6,
                   (f0.x * (m1.z * m2.x + m2.z * m3.x + m3.z * m1.x - m3.z * m2.x - m2.z * m1.x - m1.z * m3.x) - f1.x * (m0.z * m2.x + m2.z * m3.x + m3.z * m0.x - m3.z * m2.x - m2.z * m0.x - m0.z * m3.x) + f2.x * (m0.z * m1.x + m1.z * m3.x + m3.z * m0.x - m3.z * m1.x - m1.z * m0.x - m0.z * m3.x) - f3.x * (m0.z * m1.x + m1.z * m2.x + m2.z * m0.x - m2.z * m1.x - m1.z * m0.x - m0.z * m2.x)) / V6,
                   (f0.x * (m1.x * m2.y + m2.x * m3.y + m3.x * m1.y - m3.x * m2.y - m2.x * m1.y - m1.x * m3.y) - f1.x * (m0.x * m2.y + m2.x * m3.y + m3.x * m0.y - m3.x * m2.y - m2.x * m0.y - m0.x * m3.y) + f2.x * (m0.x * m1.y + m1.x * m3.y + m3.x * m0.y - m3.x * m1.y - m1.x * m0.y - m0.x * m3.y) - f3.x * (m0.x * m1.y + m1.x * m2.y + m2.x * m0.y - m2.x * m1.y - m1.x * m0.y - m0.x * m2.y)) / V6 };
    Real3 d_uy = { (f0.y * (m1.y * m2.z + m2.y * m3.z + m3.y * m1.z - m3.y * m2.z - m2.y * m1.z - m1.y * m3.z) - f1.y * (m0.y * m2.z + m2.y * m3.z + m3.y * m0.z - m3.y * m2.z - m2.y * m0.z - m0.y * m3.z) + f2.y * (m0.y * m1.z + m1.y * m3.z + m3.y * m0.z - m3.y * m1.z - m1.y * m0.z - m0.y * m3.z) - f3.y * (m0.y * m1.z + m1.y * m2.z + m2.y * m0.z - m2.y * m1.z - m1.y * m0.z - m0.y * m2.z)) / V6,
                   (f0.y * (m1.z * m2.x + m2.z * m3.x + m3.z * m1.x - m3.z * m2.x - m2.z * m1.x - m1.z * m3.x) - f1.y * (m0.z * m2.x + m2.z * m3.x + m3.z * m0.x - m3.z * m2.x - m2.z * m0.x - m0.z * m3.x) + f2.y * (m0.z * m1.x + m1.z * m3.x + m3.z * m0.x - m3.z * m1.x - m1.z * m0.x - m0.z * m3.x) - f3.y * (m0.z * m1.x + m1.z * m2.x + m2.z * m0.x - m2.z * m1.x - m1.z * m0.x - m0.z * m2.x)) / V6,
                   (f0.y * (m1.x * m2.y + m2.x * m3.y + m3.x * m1.y - m3.x * m2.y - m2.x * m1.y - m1.x * m3.y) - f1.y * (m0.x * m2.y + m2.x * m3.y + m3.x * m0.y - m3.x * m2.y - m2.x * m0.y - m0.x * m3.y) + f2.y * (m0.x * m1.y + m1.x * m3.y + m3.x * m0.y - m3.x * m1.y - m1.x * m0.y - m0.x * m3.y) - f3.y * (m0.x * m1.y + m1.x * m2.y + m2.x * m0.y - m2.x * m1.y - m1.x * m0.y - m0.x * m2.y)) / V6 };
    Real3 d_uz = { (f0.z * (m1.y * m2.z + m2.y * m3.z + m3.y * m1.z - m3.y * m2.z - m2.y * m1.z - m1.y * m3.z) - f1.z * (m0.y * m2.z + m2.y * m3.z + m3.y * m0.z - m3.y * m2.z - m2.y * m0.z - m0.y * m3.z) + f2.z * (m0.y * m1.z + m1.y * m3.z + m3.y * m0.z - m3.y * m1.z - m1.y * m0.z - m0.y * m3.z) - f3.z * (m0.y * m1.z + m1.y * m2.z + m2.y * m0.z - m2.y * m1.z - m1.y * m0.z - m0.y * m2.z)) / V6,
                   (f0.z * (m1.z * m2.x + m2.z * m3.x + m3.z * m1.x - m3.z * m2.x - m2.z * m1.x - m1.z * m3.x) - f1.z * (m0.z * m2.x + m2.z * m3.x + m3.z * m0.x - m3.z * m2.x - m2.z * m0.x - m0.z * m3.x) + f2.z * (m0.z * m1.x + m1.z * m3.x + m3.z * m0.x - m3.z * m1.x - m1.z * m0.x - m0.z * m3.x) - f3.z * (m0.z * m1.x + m1.z * m2.x + m2.z * m0.x - m2.z * m1.x - m1.z * m0.x - m0.z * m2.x)) / V6,
                   (f0.z * (m1.x * m2.y + m2.x * m3.y + m3.x * m1.y - m3.x * m2.y - m2.x * m1.y - m1.x * m3.y) - f1.z * (m0.x * m2.y + m2.x * m3.y + m3.x * m0.y - m3.x * m2.y - m2.x * m0.y - m0.x * m3.y) + f2.z * (m0.x * m1.y + m1.x * m3.y + m3.x * m0.y - m3.x * m1.y - m1.x * m0.y - m0.x * m3.y) - f3.z * (m0.x * m1.y + m1.x * m2.y + m2.x * m0.y - m2.x * m1.y - m1.x * m0.y - m0.x * m2.y)) / V6 };

    Real3x3 grad = { d_ux, d_uy, d_uz };
    return grad;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Holds information for a hexahedral element at a single Gauss point.
   *
   * This includes the gradients of the shape functions in the physical space (𝑥,𝑦,𝑧)
   * and the determinant of the Jacobian matrix.
   */
  /*---------------------------------------------------------------------------*/
  template <Int32 N> struct HexaGaussPointInfo
  {
    RealVector<N> dN_dx; // Derivatives of shape functions in x.
    RealVector<N> dN_dy; // Derivatives of shape functions in y.
    RealVector<N> dN_dz; // Derivatives of shape functions in z.
    Real det_j; // Determinant of the Jacobian matrix at the Gauss point.
  };

  using Hexa8GaussPointInfo = HexaGaussPointInfo<8>;
  using Hexa20GaussPointInfo = HexaGaussPointInfo<20>;
  using Hexa27GaussPointInfo = HexaGaussPointInfo<27>;

  template <Int32 N> static inline HexaGaussPointInfo<N>
  _computeHexaGradientsAndJacobian(Cell cell,
                                   const VariableNodeReal3& node_coord,
                                   const Arcane::FemUtils::ShapeFunctions::ReferenceGradients3D<N>& reference_gradients)
  {
    Real3x3 J;
    for (Int8 a = 0; a < N; ++a) {
      const Real3& n_coord = node_coord[cell.nodeId(a)];
      J[0][0] += reference_gradients.dN_dxi[a] * n_coord.x;
      J[0][1] += reference_gradients.dN_dxi[a] * n_coord.y;
      J[0][2] += reference_gradients.dN_dxi[a] * n_coord.z;
      J[1][0] += reference_gradients.dN_deta[a] * n_coord.x;
      J[1][1] += reference_gradients.dN_deta[a] * n_coord.y;
      J[1][2] += reference_gradients.dN_deta[a] * n_coord.z;
      J[2][0] += reference_gradients.dN_dzeta[a] * n_coord.x;
      J[2][1] += reference_gradients.dN_dzeta[a] * n_coord.y;
      J[2][2] += reference_gradients.dN_dzeta[a] * n_coord.z;
    }

    const Real detJ = math::matrixDeterminant(J);
    if (detJ <= 0.0) {
      ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
    }
    const Real3x3 invJ = math::inverseMatrix(J, detJ);

    RealVector<N> dN_dx_result;
    RealVector<N> dN_dy_result;
    RealVector<N> dN_dz_result;
    for (Int8 a = 0; a < N; ++a) {
      dN_dx_result(a) = invJ[0][0] * reference_gradients.dN_dxi[a] + invJ[0][1] * reference_gradients.dN_deta[a] + invJ[0][2] * reference_gradients.dN_dzeta[a];
      dN_dy_result(a) = invJ[1][0] * reference_gradients.dN_dxi[a] + invJ[1][1] * reference_gradients.dN_deta[a] + invJ[1][2] * reference_gradients.dN_dzeta[a];
      dN_dz_result(a) = invJ[2][0] * reference_gradients.dN_dxi[a] + invJ[2][1] * reference_gradients.dN_deta[a] + invJ[2][2] * reference_gradients.dN_dzeta[a];
    }

    return { dN_dx_result, dN_dy_result, dN_dz_result, detJ };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes shape function gradients and the Jacobian determinant for a Hexa8 element.
   *
   * @param cell The Hexa8 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param xi The ξ coordinate of the evaluation point (-1 to 1).
   * @param eta The η coordinate of the evaluation point (-1 to 1).
   * @param zeta The ζ coordinate of the evaluation point (-1 to 1).
   * @return A Hexa8GaussPointInfo struct containing the gradients and Jacobian determinant.
   */
  /*---------------------------------------------------------------------------*/

  static inline Hexa8GaussPointInfo
  computeGradientsAndJacobianHexa8(Cell cell, const VariableNodeReal3& node_coord, Real xi, Real eta, Real zeta)
  {
    const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsHexa8(xi, eta, zeta);
    return _computeHexaGradientsAndJacobian<8>(cell, node_coord, reference_gradients);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the gradient of a scalar field 'u' for a Hexa8 element.
   *
   * This method calculates gradient operator ∇ Ui for a given P1 cell
   * with i = 1,..,8 for the eight vector values of Ui hence at cell nodes.
   * The output is ∇ Ui is P0 (piece-wise constant) hence Real3 value
   * per cell
   *
   *         ∇ Ui = [ ∂U/∂𝑥   ∂U/∂𝑦   ∂U/∂𝑧 ]
   *
   * @param cell The Hexa8 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param u The variable Ui defined at the nodes.
   *
   * @return A Real3 vector of the gradient ∇Ui = {∂U/∂𝑥, ∂U/∂𝑦, ∂U/∂𝑧}
   * at the cell center.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3
  computeGradientHexa8(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal& u /*, Real xi, Real eta*/)
  {
    // get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian at (ξ,η,ζ) = (0,0,0)
    const auto gp_util = computeGradientsAndJacobianHexa8(cell, node_coord, 0.0, 0.0, 0.0);
    const RealVector<8>& dN_dx = gp_util.dN_dx;
    const RealVector<8>& dN_dy = gp_util.dN_dy;
    const RealVector<8>& dN_dz = gp_util.dN_dz;

    // get the nodal values of the variable 𝑢ᵢ ∀ 𝑖= 1,……,8 for the cell.
    const Real u_nodes[8] = {
      u[cell.nodeId(0)],
      u[cell.nodeId(1)],
      u[cell.nodeId(2)],
      u[cell.nodeId(3)],
      u[cell.nodeId(4)],
      u[cell.nodeId(5)],
      u[cell.nodeId(6)],
      u[cell.nodeId(7)]
    };

    // Compute the gradient components using shape function
    //    ∂𝑢/∂𝑥 = Σ (∂𝑁ᵢ/∂𝑥 * uᵢ) ∀ 𝑖= 1,……,8
    //    ∂𝑢/∂𝑦 = Σ (∂𝑁ᵢ/∂𝑦 * uᵢ) ∀ 𝑖= 1,……,8
    //    ∂𝑢/∂𝑧 = Σ (∂𝑁ᵢ/∂𝑧 * uᵢ) ∀ 𝑖= 1,……,8
    Real grad_x = 0.0;
    Real grad_y = 0.0;
    Real grad_z = 0.0;
    for (Int8 a = 0; a < 8; ++a) {
      grad_x += dN_dx(a) * u_nodes[a];
      grad_y += dN_dy(a) * u_nodes[a];
      grad_z += dN_dz(a) * u_nodes[a];
    }

    return { grad_x, grad_y, grad_z };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the gradient of a vector field 'u' for a Hexa8 element.
   *
   * This method calculates gradient operator ∇ Ui for a given P1 cell
   * with i = 1,..,8 for the eight vector values of Ui hence at cell nodes.
   * The output is ∇ Ui is P0 (piece-wise constant) hence Real3x3 value
   * per cell
   *
   *         ∇ Ui = [ ∂U/∂𝑥   ∂U/∂𝑦   ∂U/∂𝑧 ]
   *
   * @param cell The Hexa8 cell entity.
   * @param node_coord The coordinates of the mesh nodes.
   * @param u The variable Ui defined at the nodes.
   * @param xi The first coordinate of quadrature point.
   * @param eta The second coordinate of quadrature point.
   * @param zeta The third coordinate of quadrature point.
   *
   * @return A Real3x3 vector of the gradient ∇Ui = {∂U/∂𝑥, ∂U/∂𝑦, ∂U/∂𝑧}
   * at the cell center.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3x3
  computeGradientHexa8(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal3& u, const Real xi, const Real eta, const Real zeta)
  {
    // get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian at (ξ,η,ζ)
    const auto gp_util = computeGradientsAndJacobianHexa8(cell, node_coord, xi, eta, zeta);
    const RealVector<8>& dN_dx = gp_util.dN_dx;
    const RealVector<8>& dN_dy = gp_util.dN_dy;
    const RealVector<8>& dN_dz = gp_util.dN_dz;

    // get the nodal values of the variable 𝑢ᵢ ∀ 𝑖= 1,……,8 for the cell.
    const Real3 u_nodes[8] = {
      u[cell.nodeId(0)],
      u[cell.nodeId(1)],
      u[cell.nodeId(2)],
      u[cell.nodeId(3)],
      u[cell.nodeId(4)],
      u[cell.nodeId(5)],
      u[cell.nodeId(6)],
      u[cell.nodeId(7)]
    };

    // Compute the gradient components using shape function for each component (ₖ) of the vector field
    //    ∂𝑢ₖ/∂𝑥 = Σ (∂𝑁ᵢ/∂𝑥 * uᵢ) ∀ 𝑖= 1,……,8
    //    ∂𝑢ₖ/∂𝑦 = Σ (∂𝑁ᵢ/∂𝑦 * uᵢ) ∀ 𝑖= 1,……,8
    //    ∂𝑢ₖ/∂𝑧 = Σ (∂𝑁ᵢ/∂𝑧 * uᵢ) ∀ 𝑖= 1,……,8

    Real3 d_ux = { 0., 0., 0. };
    Real3 d_uy = { 0., 0., 0. };
    Real3 d_uz = { 0., 0., 0. };

    for (Int8 a = 0; a < 8; ++a) {
      d_ux[0] += dN_dx(a) * u_nodes[a].x;
      d_ux[1] += dN_dy(a) * u_nodes[a].x;
      d_ux[2] += dN_dz(a) * u_nodes[a].x;
      d_uy[0] += dN_dx(a) * u_nodes[a].y;
      d_uy[1] += dN_dy(a) * u_nodes[a].y;
      d_uy[2] += dN_dz(a) * u_nodes[a].y;
      d_uz[0] += dN_dx(a) * u_nodes[a].z;
      d_uz[1] += dN_dy(a) * u_nodes[a].z;
      d_uz[2] += dN_dz(a) * u_nodes[a].z;
    }

    return { d_ux, d_uy, d_uz };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes shape function gradients and the Jacobian determinant for a Hexa20 element.
   *
   * Node ordering follows ItemTypeMng.cc (identical to the VTK convention):
   *   0-7 corners, 8-11 bottom edges, 12-15 top edges, 16-19 vertical edges.
   *
   * @return A Hexa20GaussPointInfo struct containing {∂𝐍/∂𝑥, ∂𝐍/∂𝑦, ∂𝐍/∂𝑧, det(𝑱)}.
   */
  /*---------------------------------------------------------------------------*/
  static inline Hexa20GaussPointInfo
  computeGradientsAndJacobianHexa20(Cell cell, const VariableNodeReal3& node_coord, Real xi, Real eta, Real zeta)
  {
    const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsHexa20(xi, eta, zeta);
    return _computeHexaGradientsAndJacobian<20>(cell, node_coord, reference_gradients);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes shape function gradients and the Jacobian determinant for a Hexa27 element.
   *
   * Node ordering follows ItemTypeMng.cc (identical to the VTK convention):
   *   0-7 corners, 8-19 edges, 20-25 face centers, 26 body center.
   *
   * @return A Hexa27GaussPointInfo struct containing {∂𝐍/∂𝑥, ∂𝐍/∂𝑦, ∂𝐍/∂𝑧, det(𝑱)}.
   */
  /*---------------------------------------------------------------------------*/
  static inline Hexa27GaussPointInfo
  computeGradientsAndJacobianHexa27(Cell cell, const VariableNodeReal3& node_coord, Real xi, Real eta, Real zeta)
  {
    const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsHexa27(xi, eta, zeta);
    return _computeHexaGradientsAndJacobian<27>(cell, node_coord, reference_gradients);
  }
};

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

} // namespace ArcaneFemFunctions

#endif
