// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8 -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ShapeFunctions.h                                           (C) 2000-2026 */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_SHAPEFUNCTIONS_H
#define ARCANEFEM_FEMUTILS_SHAPEFUNCTIONS_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemUtils.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils::ShapeFunctions
{

/*---------------------------------------------------------------------------*/
/**
 * @brief Shape-function derivatives in a two-dimensional reference element.
 *
 * Each vector contains one derivative per element node. For example,
 * `ReferenceGradients2D<9>` contains 9 derivatives with respect to ξ and
 * 9 derivatives with respect to η.
 */
/*---------------------------------------------------------------------------*/

template <Int32 N> struct ReferenceGradients2D
{
  RealVector<N> dN_dxi;
  RealVector<N> dN_deta;
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Shape-function derivatives in a three-dimensional reference element.
 *
 * Each vector contains one derivative per element node. For example,
 * `ReferenceGradients3D<27>` contains three vectors of 27 derivatives: one
 * for each of the ξ, η and ζ reference-coordinate directions.
 */
/*---------------------------------------------------------------------------*/

template <Int32 N> struct ReferenceGradients3D
{
  RealVector<N> dN_dxi;
  RealVector<N> dN_deta;
  RealVector<N> dN_dzeta;
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the shape functions for a Line2 element at a given point (ξ).
 *
 * @param xi The ξ coordinate of the evaluation point (-1 to 1).
 * @return A RealVector<2> containing the shape functions {𝑁₁, 𝑁₂}.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline RealVector<2>
computeShapeFunctionsLine2(Real xi)
{
  RealVector<2> N;
  N(0) = 0.5 * (1.0 - xi);
  N(1) = 0.5 * (1.0 + xi);
  return N;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the shape functions for a Line3 element at a given point (ξ).
 *
 * @param xi The ξ coordinate of the evaluation point (-1 to 1).
 * @return A RealVector<3> containing the shape functions {𝑁₁, 𝑁₂, 𝑁₃}.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline RealVector<3>
computeShapeFunctionsLine3(Real xi)
{
  RealVector<3> N;
  N(0) = 0.5 * xi * (xi - 1.0);
  N(1) = 0.5 * xi * (xi + 1.0);
  N(2) = 1.0 - xi * xi;
  return N;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the shape functions for a Quad4 element at a given point (ξ, η).
 *
 * @param xi The ξ coordinate of the evaluation point (-1 to 1).
 * @param eta The η coordinate of the evaluation point (-1 to 1).
 * @return A RealVector<4> containing the shape functions {𝑁₁, 𝑁₂, 𝑁₃, 𝑁₄}.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline RealVector<4>
computeShapeFunctionsQuad4(Real xi, Real eta)
{
  RealVector<4> N;
  N(0) = 0.25 * (1.0 - xi) * (1.0 - eta);
  N(1) = 0.25 * (1.0 + xi) * (1.0 - eta);
  N(2) = 0.25 * (1.0 + xi) * (1.0 + eta);
  N(3) = 0.25 * (1.0 - xi) * (1.0 + eta);
  return N;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the reference derivatives for a Quad4 element.
 *
 * @return Four derivatives with respect to ξ and four with respect to η.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline ReferenceGradients2D<4>
computeReferenceGradientsQuad4(Real xi, Real eta)
{
  ReferenceGradients2D<4> gradients;
  gradients.dN_dxi = { -0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta) };
  gradients.dN_deta = { -0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi) };
  return gradients;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the shape functions for a Quad8 element at a given point (ξ, η).
 *
 * @param xi The ξ coordinate of the evaluation point (-1 to 1).
 * @param eta The η coordinate of the evaluation point (-1 to 1).
 * @return A RealVector<8> containing the shape functions {𝑁₁, ..., 𝑁₈}.
 *
 * Shape functions 𝐍 for Quad8 (serendipity)
 * 𝐍 = [𝑁₁  𝑁₂  𝑁₃  𝑁₄  𝑁₅  𝑁₆  𝑁₇  𝑁₈]
 * 𝑁₁ = 1/4 * (1-ξ)(1-η)(-ξ-η-1)     𝑁₅ = 1/2 * (1-ξ²)(1-η)
 * 𝑁₂ = 1/4 * (1+ξ)(1-η)( ξ-η-1)     𝑁₆ = 1/2 * (1+ξ)(1-η²)
 * 𝑁₃ = 1/4 * (1+ξ)(1+η)( ξ+η-1)     𝑁₇ = 1/2 * (1-ξ²)(1+η)
 * 𝑁₄ = 1/4 * (1-ξ)(1+η)(-ξ+η-1)     𝑁₈ = 1/2 * (1-ξ)(1-η²)
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline RealVector<8>
computeShapeFunctionsQuad8(Real xi, Real eta)
{
  RealVector<8> N;
  N[0] = 0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0);
  N[1] = 0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0);
  N[2] = 0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0);
  N[3] = 0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0);
  N[4] = 0.5 * (1.0 - xi * xi) * (1.0 - eta);
  N[5] = 0.5 * (1.0 + xi) * (1.0 - eta * eta);
  N[6] = 0.5 * (1.0 - xi * xi) * (1.0 + eta);
  N[7] = 0.5 * (1.0 - xi) * (1.0 - eta * eta);
  return N;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the reference derivatives for a Quad8 element.
 *
 * @return Eight derivatives with respect to ξ and eight with respect to η.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline ReferenceGradients2D<8>
computeReferenceGradientsQuad8(Real xi, Real eta)
{
  ReferenceGradients2D<8> gradients;
  gradients.dN_dxi = {
    0.25 * (1.0 - eta) * (2.0 * xi + eta),
    0.25 * (1.0 - eta) * (2.0 * xi - eta),
    0.25 * (1.0 + eta) * (2.0 * xi + eta),
    0.25 * (1.0 + eta) * (2.0 * xi - eta),
    -xi * (1.0 - eta),
    0.5 * (1.0 - eta * eta),
    -xi * (1.0 + eta),
    -0.5 * (1.0 - eta * eta)
  };
  gradients.dN_deta = {
    0.25 * (1.0 - xi) * (xi + 2.0 * eta),
    0.25 * (1.0 + xi) * (-xi + 2.0 * eta),
    0.25 * (1.0 + xi) * (xi + 2.0 * eta),
    0.25 * (1.0 - xi) * (-xi + 2.0 * eta),
    -0.5 * (1.0 - xi * xi),
    -eta * (1.0 + xi),
    0.5 * (1.0 - xi * xi),
    -eta * (1.0 - xi)
  };
  return gradients;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the shape functions for a Quad9 element at a given point (ξ, η).
 * @param xi The ξ coordinate of the evaluation point (-1 to 1).
 * @param eta The η coordinate of the evaluation point (-1 to 1).
 * @return A RealVector<9> containing {𝑁₁, ..., 𝑁₉}.
 *
 * Shape functions 𝐍 for Quad9 (Lagrange)
 *   𝐍 = [𝑁₁  𝑁₂  𝑁₃  𝑁₄  𝑁₅  𝑁₆  𝑁₇  𝑁₈  𝑁₉]
 *   𝑁₁ = 1/4 * ξ(ξ-1)η(η-1)     𝑁₅ = 1/2 * (1-ξ²)η(η-1)
 *   𝑁₂ = 1/4 * ξ(ξ+1)η(η-1)     𝑁₆ = 1/2 * ξ(ξ+1)(1-η²)
 *   𝑁₃ = 1/4 * ξ(ξ+1)η(η+1)     𝑁₇ = 1/2 * (1-ξ²)η(η+1)
 *   𝑁₄ = 1/4 * ξ(ξ-1)η(η+1)     𝑁₈ = 1/2 * ξ(ξ-1)(1-η²)
 *   𝑁₉ = (1-ξ²)(1-η²)
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline RealVector<9>
computeShapeFunctionsQuad9(Real xi, Real eta)
{
  RealVector<9> N;
  N[0] = 0.25 * xi * (xi - 1.0) * eta * (eta - 1.0);
  N[1] = 0.25 * xi * (xi + 1.0) * eta * (eta - 1.0);
  N[2] = 0.25 * xi * (xi + 1.0) * eta * (eta + 1.0);
  N[3] = 0.25 * xi * (xi - 1.0) * eta * (eta + 1.0);
  N[4] = 0.5 * (1.0 - xi * xi) * eta * (eta - 1.0);
  N[5] = 0.5 * xi * (xi + 1.0) * (1.0 - eta * eta);
  N[6] = 0.5 * (1.0 - xi * xi) * eta * (eta + 1.0);
  N[7] = 0.5 * xi * (xi - 1.0) * (1.0 - eta * eta);
  N[8] = (1.0 - xi * xi) * (1.0 - eta * eta);
  return N;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the reference derivatives for a Quad9 element.
 *
 * @return Nine derivatives with respect to ξ and nine with respect to η.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline ReferenceGradients2D<9>
computeReferenceGradientsQuad9(Real xi, Real eta)
{
  ReferenceGradients2D<9> gradients;
  gradients.dN_dxi = {
    0.25 * (2.0 * xi - 1.0) * eta * (eta - 1.0),
    0.25 * (2.0 * xi + 1.0) * eta * (eta - 1.0),
    0.25 * (2.0 * xi + 1.0) * eta * (eta + 1.0),
    0.25 * (2.0 * xi - 1.0) * eta * (eta + 1.0),
    -xi * eta * (eta - 1.0),
    0.5 * (2.0 * xi + 1.0) * (1.0 - eta * eta),
    -xi * eta * (eta + 1.0),
    0.5 * (2.0 * xi - 1.0) * (1.0 - eta * eta),
    -2.0 * xi * (1.0 - eta * eta)
  };
  gradients.dN_deta = {
    0.25 * xi * (xi - 1.0) * (2.0 * eta - 1.0),
    0.25 * xi * (xi + 1.0) * (2.0 * eta - 1.0),
    0.25 * xi * (xi + 1.0) * (2.0 * eta + 1.0),
    0.25 * xi * (xi - 1.0) * (2.0 * eta + 1.0),
    0.5 * (1.0 - xi * xi) * (2.0 * eta - 1.0),
    -eta * xi * (xi + 1.0),
    0.5 * (1.0 - xi * xi) * (2.0 * eta + 1.0),
    -eta * xi * (xi - 1.0),
    -2.0 * eta * (1.0 - xi * xi)
  };
  return gradients;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the shape functions for a Hexa8 element at a given point (ξ,η,ζ).
 *
 * @param xi The ξ coordinate of the evaluation point (-1 to 1).
 * @param eta The η coordinate of the evaluation point (-1 to 1).
 * @param zeta The ζ coordinate of the evaluation point (-1 to 1).
 * @return A RealVector<8> containing the shape functions {𝑁₁, ..., 𝑁₈}.
 * 
 * Shape functions 𝐍 for Hexa8
 *   𝐍 = [𝑁₁  𝑁₂  𝑁₃  𝑁₄  𝑁₅  𝑁₆  𝑁₇  𝑁₈]
 *   𝑁₁ = 1/8 * (1 - ξ) * (1 - η) * (1 - ζ)
 *   𝑁₂ = 1/8 * (1 + ξ) * (1 - η) * (1 - ζ)
 *   𝑁₃ = 1/8 * (1 + ξ) * (1 + η) * (1 - ζ)
 *   𝑁₄ = 1/8 * (1 - ξ) * (1 + η) * (1 - ζ)
 *   𝑁₅ = 1/8 * (1 - ξ) * (1 - η) * (1 + ζ)
 *   𝑁₆ = 1/8 * (1 + ξ) * (1 - η) * (1 + ζ)
 *   𝑁₇ = 1/8 * (1 + ξ) * (1 + η) * (1 + ζ)
 *   𝑁₈ = 1/8 * (1 - ξ) * (1 + η) * (1 + ζ)
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline RealVector<8>
computeShapeFunctionsHexa8(Real xi, Real eta, Real zeta)
{
  RealVector<8> N;
  N(0) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta);
  N(1) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta);
  N(2) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta);
  N(3) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta);
  N(4) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta);
  N(5) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta);
  N(6) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta);
  N(7) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta);
  return N;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the reference derivatives for a Hexa8 element.
 *
 * @return Eight derivatives in each of the ξ, η and ζ directions.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline ReferenceGradients3D<8>
computeReferenceGradientsHexa8(Real xi, Real eta, Real zeta)
{
  ReferenceGradients3D<8> gradients;
  gradients.dN_dxi = {
    -0.125 * (1.0 - eta) * (1.0 - zeta),
    0.125 * (1.0 - eta) * (1.0 - zeta),
    0.125 * (1.0 + eta) * (1.0 - zeta),
    -0.125 * (1.0 + eta) * (1.0 - zeta),
    -0.125 * (1.0 - eta) * (1.0 + zeta),
    0.125 * (1.0 - eta) * (1.0 + zeta),
    0.125 * (1.0 + eta) * (1.0 + zeta),
    -0.125 * (1.0 + eta) * (1.0 + zeta)
  };
  gradients.dN_deta = {
    -0.125 * (1.0 - xi) * (1.0 - zeta),
    -0.125 * (1.0 + xi) * (1.0 - zeta),
    0.125 * (1.0 + xi) * (1.0 - zeta),
    0.125 * (1.0 - xi) * (1.0 - zeta),
    -0.125 * (1.0 - xi) * (1.0 + zeta),
    -0.125 * (1.0 + xi) * (1.0 + zeta),
    0.125 * (1.0 + xi) * (1.0 + zeta),
    0.125 * (1.0 - xi) * (1.0 + zeta)
  };
  gradients.dN_dzeta = {
    -0.125 * (1.0 - xi) * (1.0 - eta),
    -0.125 * (1.0 + xi) * (1.0 - eta),
    -0.125 * (1.0 + xi) * (1.0 + eta),
    -0.125 * (1.0 - xi) * (1.0 + eta),
    0.125 * (1.0 - xi) * (1.0 - eta),
    0.125 * (1.0 + xi) * (1.0 - eta),
    0.125 * (1.0 + xi) * (1.0 + eta),
    0.125 * (1.0 - xi) * (1.0 + eta)
  };
  return gradients;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the shape functions for a Hexa20 element at a given point (ξ,η,ζ).
 *
 * @param xi The ξ coordinate of the evaluation point (-1 to 1).
 * @param eta The η coordinate of the evaluation point (-1 to 1).
 * @param zeta The ζ coordinate of the evaluation point (-1 to 1).
 * @return A RealVector<20> containing the shape functions {𝑁₁, ..., 𝑁₂₀}.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline RealVector<20>
computeShapeFunctionsHexa20(Real xi, Real eta, Real zeta)
{
  RealVector<20> N;
  N[0] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta) * (-xi - eta - zeta - 2.0);
  N[1] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta) * (xi - eta - zeta - 2.0);
  N[2] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta) * (xi + eta - zeta - 2.0);
  N[3] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta) * (-xi + eta - zeta - 2.0);
  N[4] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta) * (-xi - eta + zeta - 2.0);
  N[5] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta) * (xi - eta + zeta - 2.0);
  N[6] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta) * (xi + eta + zeta - 2.0);
  N[7] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta) * (-xi + eta + zeta - 2.0);
  N[8] = 0.25 * (1.0 - xi * xi) * (1.0 - eta) * (1.0 - zeta);
  N[9] = 0.25 * (1.0 + xi) * (1.0 - eta * eta) * (1.0 - zeta);
  N[10] = 0.25 * (1.0 - xi * xi) * (1.0 + eta) * (1.0 - zeta);
  N[11] = 0.25 * (1.0 - xi) * (1.0 - eta * eta) * (1.0 - zeta);
  N[12] = 0.25 * (1.0 - xi * xi) * (1.0 - eta) * (1.0 + zeta);
  N[13] = 0.25 * (1.0 + xi) * (1.0 - eta * eta) * (1.0 + zeta);
  N[14] = 0.25 * (1.0 - xi * xi) * (1.0 + eta) * (1.0 + zeta);
  N[15] = 0.25 * (1.0 - xi) * (1.0 - eta * eta) * (1.0 + zeta);
  N[16] = 0.25 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta * zeta);
  N[17] = 0.25 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta * zeta);
  N[18] = 0.25 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta * zeta);
  N[19] = 0.25 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta * zeta);
  return N;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the reference derivatives for a Hexa20 element.
 *
 * @return Twenty derivatives in each of the ξ, η and ζ directions.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline ReferenceGradients3D<20>
computeReferenceGradientsHexa20(Real xi, Real eta, Real zeta)
{
  ReferenceGradients3D<20> gradients;
  gradients.dN_dxi = {
    0.125 * (1 - eta) * (1 - zeta) * (2 * xi + eta + zeta + 1),
    0.125 * (1 - eta) * (1 - zeta) * (2 * xi - eta - zeta - 1),
    0.125 * (1 + eta) * (1 - zeta) * (2 * xi + eta - zeta - 1),
    0.125 * (1 + eta) * (1 - zeta) * (2 * xi - eta + zeta + 1),
    0.125 * (1 - eta) * (1 + zeta) * (2 * xi + eta - zeta + 1),
    0.125 * (1 - eta) * (1 + zeta) * (2 * xi - eta + zeta - 1),
    0.125 * (1 + eta) * (1 + zeta) * (2 * xi + eta + zeta - 1),
    0.125 * (1 + eta) * (1 + zeta) * (2 * xi - eta - zeta + 1),
    -0.5 * xi * (1 - eta) * (1 - zeta),
    0.25 * (1 - eta * eta) * (1 - zeta),
    -0.5 * xi * (1 + eta) * (1 - zeta),
    -0.25 * (1 - eta * eta) * (1 - zeta),
    -0.5 * xi * (1 - eta) * (1 + zeta),
    0.25 * (1 - eta * eta) * (1 + zeta),
    -0.5 * xi * (1 + eta) * (1 + zeta),
    -0.25 * (1 - eta * eta) * (1 + zeta),
    -0.25 * (1 - eta) * (1 - zeta * zeta),
    0.25 * (1 - eta) * (1 - zeta * zeta),
    0.25 * (1 + eta) * (1 - zeta * zeta),
    -0.25 * (1 + eta) * (1 - zeta * zeta)
  };
  gradients.dN_deta = {
    0.125 * (1 - xi) * (1 - zeta) * (xi + 2 * eta + zeta + 1),
    0.125 * (1 + xi) * (1 - zeta) * (-xi + 2 * eta + zeta + 1),
    0.125 * (1 + xi) * (1 - zeta) * (xi + 2 * eta - zeta - 1),
    0.125 * (1 - xi) * (1 - zeta) * (-xi + 2 * eta - zeta - 1),
    0.125 * (1 - xi) * (1 + zeta) * (xi + 2 * eta - zeta + 1),
    0.125 * (1 + xi) * (1 + zeta) * (-xi + 2 * eta - zeta + 1),
    0.125 * (1 + xi) * (1 + zeta) * (xi + 2 * eta + zeta - 1),
    0.125 * (1 - xi) * (1 + zeta) * (-xi + 2 * eta + zeta - 1),
    -0.25 * (1 - xi * xi) * (1 - zeta),
    -0.5 * eta * (1 + xi) * (1 - zeta),
    0.25 * (1 - xi * xi) * (1 - zeta),
    -0.5 * eta * (1 - xi) * (1 - zeta),
    -0.25 * (1 - xi * xi) * (1 + zeta),
    -0.5 * eta * (1 + xi) * (1 + zeta),
    0.25 * (1 - xi * xi) * (1 + zeta),
    -0.5 * eta * (1 - xi) * (1 + zeta),
    -0.25 * (1 - xi) * (1 - zeta * zeta),
    -0.25 * (1 + xi) * (1 - zeta * zeta),
    0.25 * (1 + xi) * (1 - zeta * zeta),
    0.25 * (1 - xi) * (1 - zeta * zeta)
  };
  gradients.dN_dzeta = {
    0.125 * (1 - xi) * (1 - eta) * (xi + eta + 2 * zeta + 1),
    0.125 * (1 + xi) * (1 - eta) * (-xi + eta + 2 * zeta + 1),
    0.125 * (1 + xi) * (1 + eta) * (-xi - eta + 2 * zeta + 1),
    0.125 * (1 - xi) * (1 + eta) * (xi - eta + 2 * zeta + 1),
    0.125 * (1 - xi) * (1 - eta) * (-xi - eta + 2 * zeta - 1),
    0.125 * (1 + xi) * (1 - eta) * (xi - eta + 2 * zeta - 1),
    0.125 * (1 + xi) * (1 + eta) * (xi + eta + 2 * zeta - 1),
    0.125 * (1 - xi) * (1 + eta) * (-xi + eta + 2 * zeta - 1),
    -0.25 * (1 - xi * xi) * (1 - eta),
    -0.25 * (1 + xi) * (1 - eta * eta),
    -0.25 * (1 - xi * xi) * (1 + eta),
    -0.25 * (1 - xi) * (1 - eta * eta),
    0.25 * (1 - xi * xi) * (1 - eta),
    0.25 * (1 + xi) * (1 - eta * eta),
    0.25 * (1 - xi * xi) * (1 + eta),
    0.25 * (1 - xi) * (1 - eta * eta),
    -0.5 * zeta * (1 - xi) * (1 - eta),
    -0.5 * zeta * (1 + xi) * (1 - eta),
    -0.5 * zeta * (1 + xi) * (1 + eta),
    -0.5 * zeta * (1 - xi) * (1 + eta)
  };
  return gradients;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the shape functions for a Hexa27 element at a given point (ξ,η,ζ).
 *
 * @param xi The ξ coordinate of the evaluation point (-1 to 1).
 * @param eta The η coordinate of the evaluation point (-1 to 1).
 * @param zeta The ζ coordinate of the evaluation point (-1 to 1).
 * @return A RealVector<27> containing the shape functions {𝑁₁, ..., 𝑁₂₇}.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline RealVector<27>
computeShapeFunctionsHexa27(Real xi, Real eta, Real zeta)
{
  RealVector<27> N;
  N[0] = 0.125 * xi * (xi - 1.0) * eta * (eta - 1.0) * zeta * (zeta - 1.0);
  N[1] = 0.125 * xi * (xi + 1.0) * eta * (eta - 1.0) * zeta * (zeta - 1.0);
  N[2] = 0.125 * xi * (xi + 1.0) * eta * (eta + 1.0) * zeta * (zeta - 1.0);
  N[3] = 0.125 * xi * (xi - 1.0) * eta * (eta + 1.0) * zeta * (zeta - 1.0);
  N[4] = 0.125 * xi * (xi - 1.0) * eta * (eta - 1.0) * zeta * (zeta + 1.0);
  N[5] = 0.125 * xi * (xi + 1.0) * eta * (eta - 1.0) * zeta * (zeta + 1.0);
  N[6] = 0.125 * xi * (xi + 1.0) * eta * (eta + 1.0) * zeta * (zeta + 1.0);
  N[7] = 0.125 * xi * (xi - 1.0) * eta * (eta + 1.0) * zeta * (zeta + 1.0);
  N[8] = 0.25 * (1.0 - xi * xi) * eta * (eta - 1.0) * zeta * (zeta - 1.0);
  N[9] = 0.25 * xi * (xi + 1.0) * (1.0 - eta * eta) * zeta * (zeta - 1.0);
  N[10] = 0.25 * (1.0 - xi * xi) * eta * (eta + 1.0) * zeta * (zeta - 1.0);
  N[11] = 0.25 * xi * (xi - 1.0) * (1.0 - eta * eta) * zeta * (zeta - 1.0);
  N[12] = 0.25 * (1.0 - xi * xi) * eta * (eta - 1.0) * zeta * (zeta + 1.0);
  N[13] = 0.25 * xi * (xi + 1.0) * (1.0 - eta * eta) * zeta * (zeta + 1.0);
  N[14] = 0.25 * (1.0 - xi * xi) * eta * (eta + 1.0) * zeta * (zeta + 1.0);
  N[15] = 0.25 * xi * (xi - 1.0) * (1.0 - eta * eta) * zeta * (zeta + 1.0);
  N[16] = 0.25 * xi * (xi - 1.0) * eta * (eta - 1.0) * (1.0 - zeta * zeta);
  N[17] = 0.25 * xi * (xi + 1.0) * eta * (eta - 1.0) * (1.0 - zeta * zeta);
  N[18] = 0.25 * xi * (xi + 1.0) * eta * (eta + 1.0) * (1.0 - zeta * zeta);
  N[19] = 0.25 * xi * (xi - 1.0) * eta * (eta + 1.0) * (1.0 - zeta * zeta);
  N[20] = 0.5 * xi * (xi - 1.0) * (1.0 - eta * eta) * (1.0 - zeta * zeta);
  N[21] = 0.5 * xi * (xi + 1.0) * (1.0 - eta * eta) * (1.0 - zeta * zeta);
  N[22] = 0.5 * (1.0 - xi * xi) * eta * (eta - 1.0) * (1.0 - zeta * zeta);
  N[23] = 0.5 * (1.0 - xi * xi) * eta * (eta + 1.0) * (1.0 - zeta * zeta);
  N[24] = 0.5 * (1.0 - xi * xi) * (1.0 - eta * eta) * zeta * (zeta - 1.0);
  N[25] = 0.5 * (1.0 - xi * xi) * (1.0 - eta * eta) * zeta * (zeta + 1.0);
  N[26] = (1.0 - xi * xi) * (1.0 - eta * eta) * (1.0 - zeta * zeta);
  return N;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Computes the reference derivatives for a Hexa27 element.
 *
 * The derivatives are listed explicitly in Arcane's Hexa27 node ordering.
 *
 * @return Twenty-seven derivatives in each of the ξ, η and ζ directions.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline ReferenceGradients3D<27>
computeReferenceGradientsHexa27(Real xi, Real eta, Real zeta)
{
  ReferenceGradients3D<27> gradients;
  gradients.dN_dxi = {
    0.125 * (2 * xi - 1) * eta * (eta - 1) * zeta * (zeta - 1),
    0.125 * (2 * xi + 1) * eta * (eta - 1) * zeta * (zeta - 1),
    0.125 * (2 * xi + 1) * eta * (eta + 1) * zeta * (zeta - 1),
    0.125 * (2 * xi - 1) * eta * (eta + 1) * zeta * (zeta - 1),
    0.125 * (2 * xi - 1) * eta * (eta - 1) * zeta * (zeta + 1),
    0.125 * (2 * xi + 1) * eta * (eta - 1) * zeta * (zeta + 1),
    0.125 * (2 * xi + 1) * eta * (eta + 1) * zeta * (zeta + 1),
    0.125 * (2 * xi - 1) * eta * (eta + 1) * zeta * (zeta + 1),
    0.25 * (-2 * xi) * eta * (eta - 1) * zeta * (zeta - 1),
    0.25 * (2 * xi + 1) * (1 - eta * eta) * zeta * (zeta - 1),
    0.25 * (-2 * xi) * eta * (eta + 1) * zeta * (zeta - 1),
    0.25 * (2 * xi - 1) * (1 - eta * eta) * zeta * (zeta - 1),
    0.25 * (-2 * xi) * eta * (eta - 1) * zeta * (zeta + 1),
    0.25 * (2 * xi + 1) * (1 - eta * eta) * zeta * (zeta + 1),
    0.25 * (-2 * xi) * eta * (eta + 1) * zeta * (zeta + 1),
    0.25 * (2 * xi - 1) * (1 - eta * eta) * zeta * (zeta + 1),
    0.25 * (2 * xi - 1) * eta * (eta - 1) * (1 - zeta * zeta),
    0.25 * (2 * xi + 1) * eta * (eta - 1) * (1 - zeta * zeta),
    0.25 * (2 * xi + 1) * eta * (eta + 1) * (1 - zeta * zeta),
    0.25 * (2 * xi - 1) * eta * (eta + 1) * (1 - zeta * zeta),
    0.5 * (2 * xi - 1) * (1 - eta * eta) * (1 - zeta * zeta),
    0.5 * (2 * xi + 1) * (1 - eta * eta) * (1 - zeta * zeta),
    0.5 * (-2 * xi) * eta * (eta - 1) * (1 - zeta * zeta),
    0.5 * (-2 * xi) * eta * (eta + 1) * (1 - zeta * zeta),
    0.5 * (-2 * xi) * (1 - eta * eta) * zeta * (zeta - 1),
    0.5 * (-2 * xi) * (1 - eta * eta) * zeta * (zeta + 1),
    (-2 * xi) * (1 - eta * eta) * (1 - zeta * zeta)
  };

  gradients.dN_deta = {
    0.125 * xi * (xi - 1) * (2 * eta - 1) * zeta * (zeta - 1),
    0.125 * xi * (xi + 1) * (2 * eta - 1) * zeta * (zeta - 1),
    0.125 * xi * (xi + 1) * (2 * eta + 1) * zeta * (zeta - 1),
    0.125 * xi * (xi - 1) * (2 * eta + 1) * zeta * (zeta - 1),
    0.125 * xi * (xi - 1) * (2 * eta - 1) * zeta * (zeta + 1),
    0.125 * xi * (xi + 1) * (2 * eta - 1) * zeta * (zeta + 1),
    0.125 * xi * (xi + 1) * (2 * eta + 1) * zeta * (zeta + 1),
    0.125 * xi * (xi - 1) * (2 * eta + 1) * zeta * (zeta + 1),
    0.25 * (1 - xi * xi) * (2 * eta - 1) * zeta * (zeta - 1),
    0.25 * xi * (xi + 1) * (-2 * eta) * zeta * (zeta - 1),
    0.25 * (1 - xi * xi) * (2 * eta + 1) * zeta * (zeta - 1),
    0.25 * xi * (xi - 1) * (-2 * eta) * zeta * (zeta - 1),
    0.25 * (1 - xi * xi) * (2 * eta - 1) * zeta * (zeta + 1),
    0.25 * xi * (xi + 1) * (-2 * eta) * zeta * (zeta + 1),
    0.25 * (1 - xi * xi) * (2 * eta + 1) * zeta * (zeta + 1),
    0.25 * xi * (xi - 1) * (-2 * eta) * zeta * (zeta + 1),
    0.25 * xi * (xi - 1) * (2 * eta - 1) * (1 - zeta * zeta),
    0.25 * xi * (xi + 1) * (2 * eta - 1) * (1 - zeta * zeta),
    0.25 * xi * (xi + 1) * (2 * eta + 1) * (1 - zeta * zeta),
    0.25 * xi * (xi - 1) * (2 * eta + 1) * (1 - zeta * zeta),
    0.5 * xi * (xi - 1) * (-2 * eta) * (1 - zeta * zeta),
    0.5 * xi * (xi + 1) * (-2 * eta) * (1 - zeta * zeta),
    0.5 * (1 - xi * xi) * (2 * eta - 1) * (1 - zeta * zeta),
    0.5 * (1 - xi * xi) * (2 * eta + 1) * (1 - zeta * zeta),
    0.5 * (1 - xi * xi) * (-2 * eta) * zeta * (zeta - 1),
    0.5 * (1 - xi * xi) * (-2 * eta) * zeta * (zeta + 1),
    (1 - xi * xi) * (-2 * eta) * (1 - zeta * zeta)
  };

  gradients.dN_dzeta = {
    0.125 * xi * (xi - 1) * eta * (eta - 1) * (2 * zeta - 1),
    0.125 * xi * (xi + 1) * eta * (eta - 1) * (2 * zeta - 1),
    0.125 * xi * (xi + 1) * eta * (eta + 1) * (2 * zeta - 1),
    0.125 * xi * (xi - 1) * eta * (eta + 1) * (2 * zeta - 1),
    0.125 * xi * (xi - 1) * eta * (eta - 1) * (2 * zeta + 1),
    0.125 * xi * (xi + 1) * eta * (eta - 1) * (2 * zeta + 1),
    0.125 * xi * (xi + 1) * eta * (eta + 1) * (2 * zeta + 1),
    0.125 * xi * (xi - 1) * eta * (eta + 1) * (2 * zeta + 1),
    0.25 * (1 - xi * xi) * eta * (eta - 1) * (2 * zeta - 1),
    0.25 * xi * (xi + 1) * (1 - eta * eta) * (2 * zeta - 1),
    0.25 * (1 - xi * xi) * eta * (eta + 1) * (2 * zeta - 1),
    0.25 * xi * (xi - 1) * (1 - eta * eta) * (2 * zeta - 1),
    0.25 * (1 - xi * xi) * eta * (eta - 1) * (2 * zeta + 1),
    0.25 * xi * (xi + 1) * (1 - eta * eta) * (2 * zeta + 1),
    0.25 * (1 - xi * xi) * eta * (eta + 1) * (2 * zeta + 1),
    0.25 * xi * (xi - 1) * (1 - eta * eta) * (2 * zeta + 1),
    0.25 * xi * (xi - 1) * eta * (eta - 1) * (-2 * zeta),
    0.25 * xi * (xi + 1) * eta * (eta - 1) * (-2 * zeta),
    0.25 * xi * (xi + 1) * eta * (eta + 1) * (-2 * zeta),
    0.25 * xi * (xi - 1) * eta * (eta + 1) * (-2 * zeta),
    0.5 * xi * (xi - 1) * (1 - eta * eta) * (-2 * zeta),
    0.5 * xi * (xi + 1) * (1 - eta * eta) * (-2 * zeta),
    0.5 * (1 - xi * xi) * eta * (eta - 1) * (-2 * zeta),
    0.5 * (1 - xi * xi) * eta * (eta + 1) * (-2 * zeta),
    0.5 * (1 - xi * xi) * (1 - eta * eta) * (2 * zeta - 1),
    0.5 * (1 - xi * xi) * (1 - eta * eta) * (2 * zeta + 1),
    (1 - xi * xi) * (1 - eta * eta) * (-2 * zeta)
  };

  return gradients;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils::ShapeFunctions

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
