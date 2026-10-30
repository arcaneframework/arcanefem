// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrixHexQuad.h                                      (C) 2000-2026 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Poisson        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD4, Q1 FE).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫ (∂𝑢/∂𝑥 ∂𝑣/∂𝑥  + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦)dΩ
 *
 * Steps involved:
 * 1. Define Gauss points (2x2) and weights.
 * 2. Loop over Gauss points to compute the gradients in physical space
 *    and the determinant of the Jacobian, via computeGradientsAndJacobianQuad4.
 * 3. Compute the integration weight.
 * 4. Assemble the element matrix using the computed gradients.
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

#include "ElementMapping.h"

RealMatrix<4, 4> FemModulePoisson::_computeElementMatrixQuad4(Cell cell)
{
  // Gauss points and weights for 2x2 quadrature
  constexpr Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 };
  constexpr Real w = 1.0;

  // Element mapping for the quadrilateral element using the coordinates of its nodes
  const FemUtils::ElementMapping<FemUtils::Quad4> mapping(
  m_node_coord[cell.nodeId(0)], m_node_coord[cell.nodeId(1)],
  m_node_coord[cell.nodeId(2)], m_node_coord[cell.nodeId(3)]);

  // Initialize the element matrix (default constructor sets all elements to zero)
  RealMatrix<4, 4> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {

      // Get the coordinates of the Gauss point in natural coordinates (ξ,η,ζ)
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Calculate the Jacobian and its determinant
      const auto jacobian_info = mapping.jacobianAndDeterminant(xi, eta);
      const Real detJ = jacobian_info.det_j;
      if (detJ <= 0.0)
        ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);

      // Compute the gradients of the shape functions in physical space (𝑥,𝑦)
      const auto gp_gradients = mapping.gradients(jacobian_info, xi, eta);
      const RealVector<4>& dxU = gp_gradients.dN_dx;
      const RealVector<4>& dyU = gp_gradients.dN_dy;

      // Integration weight
      const Real integration_weight = detJ * w * w;

      // stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ)   with  𝑎 ∈ ℝ⁴ˣ⁴
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j]);
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<4, 4>
_computeElementMatrixQuad4Gpu(CellLocalId cell_lid,
                              const IndexedCellNodeConnectivityView& cn_cv,
                              const ax::VariableNodeReal3InView& in_node_coord)
{
  // Gauss points and weights for 2x2 quadrature
  constexpr Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 };
  constexpr Real w = 1.0;

  // Element mapping for the quadrilateral element using the coordinates of its nodes
  const FemUtils::ElementMapping<FemUtils::Quad4> mapping(
  in_node_coord[cn_cv.nodeId(cell_lid, 0)], in_node_coord[cn_cv.nodeId(cell_lid, 1)],
  in_node_coord[cn_cv.nodeId(cell_lid, 2)], in_node_coord[cn_cv.nodeId(cell_lid, 3)]);

  // Initialize the element matrix (default constructor sets all elements to zero)
  RealMatrix<4, 4> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {

      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Calculate the Jacobian and its determinant
      const auto jacobian_info = mapping.jacobianAndDeterminant(xi, eta);
      const Real detJ = jacobian_info.det_j;

      // Compute the gradients of the shape functions in physical space (𝑥,𝑦)
      const auto gp_gradients = mapping.gradients(jacobian_info, xi, eta);
      const RealVector<4>& dxU = gp_gradients.dN_dx;
      const RealVector<4>& dyU = gp_gradients.dN_dy;

      // Integration weight
      const Real integration_weight = detJ * w * w;

      // stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ)   with  𝑎 ∈ ℝ⁴ˣ⁴
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j]);
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<1, 4>
_computeElementVectorQuad4Gpu(CellLocalId cell_lid,
                              const IndexedCellNodeConnectivityView& cn_cv,
                              const ax::VariableNodeReal3InView& in_node_coord, Int32 node_lid)
{
  // Gauss points and weights for 2x2 quadrature
  constexpr Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 };
  constexpr Real w = 1.0;

  // Element mapping for the quadrilateral element using the coordinates of its nodes
  const FemUtils::ElementMapping<FemUtils::Quad4> mapping(
  in_node_coord[cn_cv.nodeId(cell_lid, 0)], in_node_coord[cn_cv.nodeId(cell_lid, 1)],
  in_node_coord[cn_cv.nodeId(cell_lid, 2)], in_node_coord[cn_cv.nodeId(cell_lid, 3)]);

  // Initialize the element vector
  RealVector<4> a_i;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {

      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Calculate the Jacobian and its determinant
      const auto jacobian_info = mapping.jacobianAndDeterminant(xi, eta);
      const Real detJ = jacobian_info.det_j;

      // Compute the gradients of the shape functions in physical space (𝑥,𝑦)
      const auto gp_gradients = mapping.gradients(jacobian_info, xi, eta);
      const RealVector<4>& dxU = gp_gradients.dN_dx;
      const RealVector<4>& dyU = gp_gradients.dN_dy;

      // Integration weight
      const Real integration_weight = detJ * w * w;

      // Stiffness vector 𝑎(𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ), with fixed 𝑖 and 𝑎 ∈ ℝ⁴
      for (int j = 0; j < 4; ++j)
        a_i[j] += integration_weight * (dxU[node_lid] * dxU[j] + dyU[node_lid] * dyU[j]);
    }
  }
  return a_i;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a hexahedral element (HEXA8, ℚ1 FE).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫∫ (∂𝑢/∂𝑥 ∂𝑣/∂𝑥 + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦 + ∂𝑢/∂𝑧 ∂𝑣/∂𝑧)dΩ
 *
 * Steps involved:
 * 1. Define Gauss points (2x2x2) and weights.
 * 2. Loop over Gauss points to compute the gradients in physical space (𝑥,𝑦,𝑧)
 *    and the determinant of the Jacobian via ElementMapping<Hexa8>.
 * 3. Compute the integration weight.
 * 4. Assemble the element matrix using the computed gradients.
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<8, 8> FemModulePoisson::_computeElementMatrixHexa8(Cell cell)
{
  // 2x2x2 Gauss points and weights for [-1,1]^3
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
  constexpr Real w = 1.0;

  const FemUtils::ElementMapping<FemUtils::Hexa8> mapping(
  m_node_coord[cell.nodeId(0)], m_node_coord[cell.nodeId(1)],
  m_node_coord[cell.nodeId(2)], m_node_coord[cell.nodeId(3)],
  m_node_coord[cell.nodeId(4)], m_node_coord[cell.nodeId(5)],
  m_node_coord[cell.nodeId(6)], m_node_coord[cell.nodeId(7)]);

  // Initialize the element matrix
  RealMatrix<8, 8> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      for (Int8 izeta = 0; izeta < 2; ++izeta) {

        // Get the coordinates of Gauss points in natural coordinates (ξ,η,ζ)
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        // Get shape function gradients w.r.t (𝑥,𝑦,𝑧) and determinant of Jacobian
        const auto jacobian_info = mapping.jacobianAndDeterminant(xi, eta, zeta);
        const Real detJ = jacobian_info.det_j;
        if (detJ <= 0.0)
          ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
        const auto gp_gradients = mapping.gradients(jacobian_info, xi, eta, zeta);
        const RealVector<8>& dxU = gp_gradients.dN_dx;
        const RealVector<8>& dyU = gp_gradients.dN_dy;
        const RealVector<8>& dzU = gp_gradients.dN_dz;

        // Integration weight
        const Real integration_weight = detJ * w * w;

        // Stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ + ∂𝓏𝑁ᵢ ∂𝓏𝑁ⱼ), with 𝑎 ∈ ℝ⁸ˣ⁸
        for (Int8 i = 0; i < 8; ++i)
          for (Int8 j = 0; j < 8; ++j)
            a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j] + dzU[i] * dzU[j]);
      }
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<8, 8>
_computeElementMatrixHexa8Gpu(CellLocalId cell_lid,
                              const IndexedCellNodeConnectivityView& cn_cv,
                              const ax::VariableNodeReal3InView& in_node_coord)
{
  // 2x2x2 Gauss points and weights for [-1,1]^3
  constexpr Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 }; // -1/sqrt(3), 1/sqrt(3)
  constexpr Real w = 1.0;

  const FemUtils::ElementMapping<FemUtils::Hexa8> mapping(
  in_node_coord[cn_cv.nodeId(cell_lid, 0)], in_node_coord[cn_cv.nodeId(cell_lid, 1)],
  in_node_coord[cn_cv.nodeId(cell_lid, 2)], in_node_coord[cn_cv.nodeId(cell_lid, 3)],
  in_node_coord[cn_cv.nodeId(cell_lid, 4)], in_node_coord[cn_cv.nodeId(cell_lid, 5)],
  in_node_coord[cn_cv.nodeId(cell_lid, 6)], in_node_coord[cn_cv.nodeId(cell_lid, 7)]);

  // Initialize the element matrix
  RealMatrix<8, 8> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      for (Int8 izeta = 0; izeta < 2; ++izeta) {

        // Get the coordinates of Gauss points in natural coordinates (ξ,η,ζ)
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        // Get shape function gradients w.r.t (𝑥,𝑦,𝑧) and determinant of Jacobian
        const auto jacobian_info = mapping.jacobianAndDeterminant(xi, eta, zeta);
        const Real detJ = jacobian_info.det_j;
        const auto gp_gradients = mapping.gradients(jacobian_info, xi, eta, zeta);
        const RealVector<8>& dxU = gp_gradients.dN_dx;
        const RealVector<8>& dyU = gp_gradients.dN_dy;
        const RealVector<8>& dzU = gp_gradients.dN_dz;

        // Integration weight
        const Real integration_weight = detJ * w * w;

        // Stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ + ∂𝓏𝑁ᵢ ∂𝓏𝑁ⱼ), with 𝑎 ∈ ℝ⁸ˣ⁸
        for (Int8 i = 0; i < 8; ++i)
          for (Int8 j = 0; j < 8; ++j)
            a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j] + dzU[i] * dzU[j]);
      }
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<1, 8>
_computeElementVectorHexa8Gpu(CellLocalId cell_lid,
                              const IndexedCellNodeConnectivityView& cn_cv,
                              const ax::VariableNodeReal3InView& in_node_coord,
                              Int32 node_lid)
{
  constexpr Real gp[2] = { -0.57735026918962576451, 0.57735026918962576451 };
  constexpr Real w = 1.0;

  const FemUtils::ElementMapping<FemUtils::Hexa8> mapping(
  in_node_coord[cn_cv.nodeId(cell_lid, 0)], in_node_coord[cn_cv.nodeId(cell_lid, 1)],
  in_node_coord[cn_cv.nodeId(cell_lid, 2)], in_node_coord[cn_cv.nodeId(cell_lid, 3)],
  in_node_coord[cn_cv.nodeId(cell_lid, 4)], in_node_coord[cn_cv.nodeId(cell_lid, 5)],
  in_node_coord[cn_cv.nodeId(cell_lid, 6)], in_node_coord[cn_cv.nodeId(cell_lid, 7)]);

  RealVector<8> a_i;

  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      for (Int8 izeta = 0; izeta < 2; ++izeta) {
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];
        const auto jacobian_info = mapping.jacobianAndDeterminant(xi, eta, zeta);
        const Real detJ = jacobian_info.det_j;
        const auto gp_gradients = mapping.gradients(jacobian_info, xi, eta, zeta);
        const RealVector<8>& dxU = gp_gradients.dN_dx;
        const RealVector<8>& dyU = gp_gradients.dN_dy;
        const RealVector<8>& dzU = gp_gradients.dN_dz;

        // Integration weight
        const Real integration_weight = detJ * w * w;

        // Stiffness vector 𝑎(𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ + ∂𝓏𝑁ᵢ ∂𝓏𝑁ⱼ), with fixed 𝑖 and 𝑎 ∈ ℝ⁸
        for (Int8 j = 0; j < 8; ++j)
          a_i[j] +=  integration_weight * (dxU[node_lid] * dxU[j] + dyU[node_lid] * dyU[j] + dzU[node_lid] * dzU[j]);
      }
    }
  }

  return a_i;
}
/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a hexahedral element (HEXA20, ℚ2 FE, serendipity).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫∫ (∂𝑢/∂𝑥 ∂𝑣/∂𝑥 + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦 + ∂𝑢/∂𝑧 ∂𝑣/∂𝑧)dΩ
 *
 * Steps involved:
 * 1. Define Gauss points (3x3x3) and weights.
 * 2. Loop over Gauss points to compute the gradients in physical space (𝑥,𝑦,𝑧)
 *    and the determinant of the Jacobian via computeGradientsAndJacobianHexa20.
 * 3. Compute the integration weight.
 * 4. Assemble the element matrix using the computed gradients.
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<20, 20> FemModulePoisson::_computeElementMatrixHexa20(Cell cell)
{
  // Gauss points and weights for 3x3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(3/5) , 0 , sqrt(3/5)]
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element matrix
  RealMatrix<20, 20> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {
      for (Int8 izeta = 0; izeta < 3; ++izeta) {

        // Get the coordinates of Gauss points in natural coordinates (ξ,η,ζ)
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        // Get shape function gradients w.r.t (𝑥,𝑦,𝑧) and determinant of Jacobian
        const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa20(cell, m_node_coord, xi, eta, zeta);
        const RealVector<20>& dxU = gp_info.dN_dx;
        const RealVector<20>& dyU = gp_info.dN_dy;
        const RealVector<20>& dzU = gp_info.dN_dz;
        const Real detJ = gp_info.det_j;

        // Integration weight
        const Real integration_weight = detJ * weight[ixi] * weight[ieta] * weight[izeta];

        // Assemble element matrix (variational form)
        for (Int8 i = 0; i < 20; ++i)
          for (Int8 j = 0; j < 20; ++j)
            a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j] + dzU[i] * dzU[j]);
      }
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<20, 20>
_computeElementMatrixHexa20Gpu(CellLocalId cell_lid,
                               const IndexedCellNodeConnectivityView& cn_cv,
                               const ax::VariableNodeReal3InView& in_node_coord)
{
  // Gauss points and weights for 3x3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element matrix
  RealMatrix<20, 20> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {
      for (Int8 izeta = 0; izeta < 3; ++izeta) {

        // Get the coordinates of Gauss points in natural coordinates (ξ,η,ζ)
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        // Get shape function gradients w.r.t (𝑥,𝑦,𝑧) and determinant of Jacobian
        const auto gp_info = FemUtils::Gpu::FeOperation3D::computeGradientsAndJacobianHexa20Gpu(cell_lid, cn_cv, in_node_coord, xi, eta, zeta);
        const RealVector<20>& dxU = gp_info.dN_dx;
        const RealVector<20>& dyU = gp_info.dN_dy;
        const RealVector<20>& dzU = gp_info.dN_dz;
        const Real detJ = gp_info.det_j;

        // Integration weight
        const Real integration_weight = detJ * weight[ixi] * weight[ieta] * weight[izeta];

        // Assemble element matrix (variational form)
        for (Int8 i = 0; i < 20; ++i) {
          for (Int8 j = 0; j < 20; ++j) {
            a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j] + dzU[i] * dzU[j]);
          }
        }
      }
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<1, 20>
_computeElementVectorHexa20Gpu(CellLocalId cell_lid,
                               const IndexedCellNodeConnectivityView& cn_cv,
                               const ax::VariableNodeReal3InView& in_node_coord,
                               Int32 node_lid)
{
  RealVector<20> a_i;

  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {
      for (Int8 izeta = 0; izeta < 3; ++izeta) {
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        const auto gp_info = FemUtils::Gpu::FeOperation3D::computeGradientsAndJacobianHexa20Gpu(cell_lid, cn_cv, in_node_coord, xi, eta, zeta);
        const RealVector<20>& dxU = gp_info.dN_dx;
        const RealVector<20>& dyU = gp_info.dN_dy;
        const RealVector<20>& dzU = gp_info.dN_dz;
        const Real detJ = gp_info.det_j;

        // Integration weight
        const Real integration_weight = detJ * weight[ixi] * weight[ieta] * weight[izeta];

        // Stiffness vector 𝑎(𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ + ∂𝓏𝑁ᵢ ∂𝓏𝑁ⱼ), with fixed 𝑖 and 𝑎 ∈ ℝ²⁰
        for (Int32 j = 0; j < 20; ++j)
          a_i[j] += integration_weight * (dxU[node_lid] * dxU[j] + dyU[node_lid] * dyU[j] + dzU[node_lid] * dzU[j]);
      }
    }
  }

  return a_i;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a hexahedral element (HEXA27, ℚ2 FE).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫∫ (∂𝑢/∂𝑥 ∂𝑣/∂𝑥 + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦 + ∂𝑢/∂𝑧 ∂𝑣/∂𝑧)dΩ
 *
 * Steps involved:
 * 1. Define Gauss points (3x3x3) and weights.
 * 2. Loop over Gauss points to compute the gradients in physical space (𝑥,𝑦,𝑧)
 *    and the determinant of the Jacobian via computeGradientsAndJacobianHexa27.
 * 3. Compute the integration weight.
 * 4. Assemble the element matrix using the computed gradients.
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<27, 27> FemModulePoisson::_computeElementMatrixHexa27(Cell cell)
{
  // Gauss points and weights for 3x3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(3/5) , 0 , sqrt(3/5)]
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element matrix
  RealMatrix<27, 27> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {
      for (Int8 izeta = 0; izeta < 3; ++izeta) {

        // Get the coordinates of Gauss points in natural coordinates (ξ,η,ζ)
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        // Get shape function gradients w.r.t (𝑥,𝑦,𝑧) and determinant of Jacobian
        const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa27(cell, m_node_coord, xi, eta, zeta);
        const RealVector<27>& dxU = gp_info.dN_dx;
        const RealVector<27>& dyU = gp_info.dN_dy;
        const RealVector<27>& dzU = gp_info.dN_dz;
        const Real detJ = gp_info.det_j;

        // Integration weight
        const Real integration_weight = detJ * weight[ixi] * weight[ieta] * weight[izeta];

        // Stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ + ∂𝓏𝑁ᵢ ∂𝓏𝑁ⱼ), with 𝑎 ∈ ℝ²⁷ˣ²⁷
        for (Int8 i = 0; i < 27; ++i)
          for (Int8 j = 0; j < 27; ++j)
            a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j] + dzU[i] * dzU[j]);
      }
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<27, 27>
_computeElementMatrixHexa27Gpu(CellLocalId cell_lid,
                               const IndexedCellNodeConnectivityView& cn_cv,
                               const ax::VariableNodeReal3InView& in_node_coord)
{
  // Gauss points and weights for 3x3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element matrix
  RealMatrix<27, 27> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {
      for (Int8 izeta = 0; izeta < 3; ++izeta) {

        // Get the coordinates of Gauss points in natural coordinates (ξ,η,ζ)
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        // Get shape function gradients w.r.t (𝑥,𝑦,𝑧) and determinant of Jacobian
        const auto gp_info = FemUtils::Gpu::FeOperation3D::computeGradientsAndJacobianHexa27Gpu(cell_lid, cn_cv, in_node_coord, xi, eta, zeta);
        const RealVector<27>& dxU = gp_info.dN_dx;
        const RealVector<27>& dyU = gp_info.dN_dy;
        const RealVector<27>& dzU = gp_info.dN_dz;
        const Real detJ = gp_info.det_j;

        // Integration weight
        const Real integration_weight = detJ * weight[ixi] * weight[ieta] * weight[izeta];

        // Stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ + ∂𝓏𝑁ᵢ ∂𝓏𝑁ⱼ), with 𝑎 ∈ ℝ²⁷ˣ²⁷
        for (Int8 i = 0; i < 27; ++i)
          for (Int8 j = 0; j < 27; ++j)
            a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j] + dzU[i] * dzU[j]);
      }
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<1, 27>
_computeElementVectorHexa27Gpu(CellLocalId cell_lid,
                               const IndexedCellNodeConnectivityView& cn_cv,
                               const ax::VariableNodeReal3InView& in_node_coord,
                               Int32 node_lid)
{
  RealMatrix<1, 27> a_i;

  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {
      for (Int8 izeta = 0; izeta < 3; ++izeta) {
        const Real xi = gp[ixi];
        const Real eta = gp[ieta];
        const Real zeta = gp[izeta];

        const auto gp_info = FemUtils::Gpu::FeOperation3D::computeGradientsAndJacobianHexa27Gpu(cell_lid, cn_cv, in_node_coord, xi, eta, zeta);
        const RealVector<27>& dxU = gp_info.dN_dx;
        const RealVector<27>& dyU = gp_info.dN_dy;
        const RealVector<27>& dzU = gp_info.dN_dz;
        const Real detJ = gp_info.det_j;

        // Integration weight
        const Real integration_weight = detJ * weight[ixi] * weight[ieta] * weight[izeta];

        // Stiffness vector 𝑎(𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ + ∂𝓏𝑁ᵢ ∂𝓏𝑁ⱼ), with fixed 𝑖 and 𝑎 ∈ ℝ²⁷
        for (Int32 j = 0; j < 27; ++j)
          a_i(0, j) += integration_weight * (dxU(node_lid) * dxU(j) + dyU(node_lid) * dyU(j) + dzU(node_lid) * dzU(j));
      }
    }
  }

  return a_i;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD8, Q2 FE).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫ (∂𝑢/∂𝑥 ∂𝑣/∂𝑥  + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦)dΩ
 *
 * Steps involved:
 * 1. Define Gauss points (3x3) and weights.
 * 2. Loop over Gauss points to compute the gradients in physical space
 *    and the determinant of the Jacobian, via computeGradientsAndJacobianQuad8.
 * 3. Compute the integration weight.
 * 4. Assemble the element matrix using the computed gradients.
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<8, 8> FemModulePoisson::_computeElementMatrixQuad8(Cell cell)
{
  // Gauss points and weights for 3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(3/5) , 0 , sqrt(3/5)]
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element matrix
  RealMatrix<8, 8> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {

      // Get the coordinates of the Gauss point in natural coordinates (ξ,η,ζ)
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian
      const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad8(cell, m_node_coord, xi, eta);
      const RealVector<8>& dxU = gp_info.dN_dx;
      const RealVector<8>& dyU = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      // Integration weight
      const Real integration_weight = detJ * weight[ixi] * weight[ieta];

      // Stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ), with 𝑎 ∈ ℝ⁸ˣ⁸
      for (Int8 i = 0; i < 8; ++i)
        for (Int8 j = 0; j < 8; ++j)
          a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j]);
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<8, 8>
_computeElementMatrixQuad8Gpu(CellLocalId cell_lid,
                              const IndexedCellNodeConnectivityView& cn_cv,
                              const ax::VariableNodeReal3InView& in_node_coord)
{

  // Gauss points and weights for 3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element matrix
  RealMatrix<8, 8> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {

      // Get the coordinates of the Gauss point in natural coordinates (ξ,η,ζ)
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian
      const auto gp_info = FemUtils::Gpu::FeOperation2D::computeGradientsAndJacobianQuad8Gpu(cell_lid, cn_cv, in_node_coord, xi, eta);
      const RealVector<8>& dxU = gp_info.dN_dx;
      const RealVector<8>& dyU = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      // Integration weight
      const Real integration_weight = detJ * weight[ixi] * weight[ieta];

      // Stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ), with 𝑎 ∈ ℝ⁸ˣ⁸
      for (Int8 i = 0; i < 8; ++i)
        for (Int8 j = 0; j < 8; ++j)
          a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j]);
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<1, 8>
_computeElementVectorQuad8Gpu(CellLocalId cell_lid,
                              const IndexedCellNodeConnectivityView& cn_cv,
                              const ax::VariableNodeReal3InView& in_node_coord,
                              Int32 node_lid)
{

  // Gauss points and weights for 3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element vector
  RealVector<8> a_i;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {

      // Get the coordinates of the Gauss point in natural coordinates (ξ,η,ζ)
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian
      const auto gp_info = FemUtils::Gpu::FeOperation2D::computeGradientsAndJacobianQuad8Gpu(cell_lid, cn_cv, in_node_coord, xi, eta);
      const RealVector<8>& dxU = gp_info.dN_dx;
      const RealVector<8>& dyU = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      // Integration weight
      const Real integration_weight = detJ * weight[ixi] * weight[ieta];

      // Stiffness vector 𝑎(𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ), with fixed 𝑖 and 𝑎 ∈ ℝ⁸
      for (Int8 j = 0; j < 8; ++j)
        a_i[j] += integration_weight * (dxU[node_lid] * dxU[j] + dyU[node_lid] * dyU[j]);
    }
  }
  return a_i;
}

/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD9, Q2 FE).
 *
 * This function calculates the integral of:
 *       𝑎(𝑢,𝑣) = ∫∫ (∂𝑢/∂𝑥 ∂𝑣/∂𝑥  + ∂𝑢/∂𝑦 ∂𝑣/∂𝑦)dΩ
 *
 * Steps involved:
 * 1. Define Gauss points (3x3) and weights.
 * 2. Loop over Gauss points to compute the gradients in physical space
 *    and the determinant of the Jacobian, via computeGradientsAndJacobianQuad8.
 * 3. Compute the integration weight.
 * 4. Assemble the element matrix using the computed gradients.
 *
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<9, 9> FemModulePoisson::
_computeElementMatrixQuad9(Cell cell)
{
  // Gauss points and weights for 3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(3/5) , 0 , sqrt(3/5)]
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element matrix
  RealMatrix<9, 9> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {

      // Get the coordinates of the Gauss point in natural coordinates (ξ,η,ζ)
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian
      const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad9(cell, m_node_coord, xi, eta);
      const RealVector<9>& dxU = gp_info.dN_dx;
      const RealVector<9>& dyU = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      // Integration weight
      const Real integration_weight = detJ * weight[ixi] * weight[ieta];

      // Stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ), with 𝑎 ∈ ℝ⁹ˣ⁹
      for (Int8 i = 0; i < 9; ++i)
        for (Int8 j = 0; j < 9; ++j)
          a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j]);
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<9, 9>
_computeElementMatrixQuad9Gpu(CellLocalId cell_lid,
                              const IndexedCellNodeConnectivityView& cn_cv,
                              const ax::VariableNodeReal3InView& in_node_coord)
{

  // Gauss points and weights for 3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(3/5) , 0 , sqrt(3/5)]
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element matrix
  RealMatrix<9, 9> a;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {

      // Get the coordinates of the Gauss point in natural coordinates (ξ,η,ζ)
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian
      const auto gp_info = FemUtils::Gpu::FeOperation2D::computeGradientsAndJacobianQuad9Gpu(cell_lid, cn_cv, in_node_coord, xi, eta);
      const RealVector<9>& dxU = gp_info.dN_dx;
      const RealVector<9>& dyU = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      // Integration weight
      const Real integration_weight = detJ * weight[ixi] * weight[ieta];

      // Stiffness matrix 𝑎(𝑖,𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ), with 𝑎 ∈ ℝ⁹ˣ⁹
      for (Int8 i = 0; i < 9; ++i)
        for (Int8 j = 0; j < 9; ++j)
          a(i, j) += integration_weight * (dxU[i] * dxU[j] + dyU[i] * dyU[j]);
    }
  }
  return a;
}

ARCCORE_HOST_DEVICE RealMatrix<1, 9>
_computeElementVectorQuad9Gpu(CellLocalId cell_lid,
                              const IndexedCellNodeConnectivityView& cn_cv,
                              const ax::VariableNodeReal3InView& in_node_coord,
                              Int32 node_lid)
{

  // Gauss points and weights for 3x3 quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  // Initialize the element vector
  RealVector<9> a_i;

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 3; ++ixi) {
    for (Int8 ieta = 0; ieta < 3; ++ieta) {

      // Get the coordinates of the Gauss point in natural coordinates (ξ,η,ζ)
      const Real xi = gp[ixi];
      const Real eta = gp[ieta];

      // Get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian
      const auto gp_info = FemUtils::Gpu::FeOperation2D::computeGradientsAndJacobianQuad9Gpu(cell_lid, cn_cv, in_node_coord, xi, eta);
      const RealVector<9>& dxU = gp_info.dN_dx;
      const RealVector<9>& dyU = gp_info.dN_dy;
      const Real detJ = gp_info.det_j;

      // Integration weight
      const Real integration_weight = detJ * weight[ixi] * weight[ieta];

      // Stiffness vector 𝑎(𝑗) += 𝑤 (∂ₓ𝑁ᵢ ∂ₓ𝑁ⱼ + ∂ᵧ𝑁ᵢ ∂ᵧ𝑁ⱼ), with fixed 𝑖 and 𝑎 ∈ ℝ⁹
      for (Int8 j = 0; j < 9; ++j)
        a_i[j] += integration_weight * (dxU[node_lid] * dxU[j] + dyU[node_lid] * dyU[j]);
    }
  }
  return a_i;
}
