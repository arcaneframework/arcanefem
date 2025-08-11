// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrixHexQuad.h                                      (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute the FEM element matrices for Poisson        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a quadrilateral element (QUAD4, ℙ1 FE).
 *
 * This function calculates the integral of:
 *       a(u,v) = ∫∫ (∂u/∂x ∂v/∂x  + ∂u/∂y ∂v/∂y)dΩ
 * 
 * steps involved:
 * 1. Define Gauss points (2x2) and weights.
 * 2. Loop over Gauss points to compute the shape function derivatives.
 * 3. Compute the Jacobian matrix and its determinant.
 * 4. Compute the inverse of the Jacobian.
 * 5. Compute the gradients in physical space.
 * 6. Assemble the element matrix using the computed gradients.
 * 
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<4, 4> FemModule::_computeElementMatrixQuad4(Cell cell)
{
  // Gauss points and weights for 2x2 quadrature
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // [-1/sqrt(3) , 1/sqrt(3)]
  constexpr Real w = 1.0;

  // Initialize the element matrix
  RealMatrix<4, 4> ae;
  ae.fill(0.0);

  // Loop over Gauss points
  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {

      // Get the coordinates of the Gauss point
      Real xi = gp[ixi]; // Get the ξ coordinate of the Gauss point
      Real eta = gp[ieta]; // Get the ∂η coordinate of the Gauss point

      // Shape function derivatives ∂𝐍/∂ξ and ∂𝐍/∂η
      //     ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ ]
      //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η ]

      Real dN_dxi[4] = { -0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta) };
      Real dN_deta[4] = { -0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi) };

      // Jacobian calculation 𝑱
      //    𝑱 = [ 𝒋₀₀  𝒋₀₁ ] = [ ∂x/∂ξ  ∂y/∂ξ ]
      //        [ 𝒋₁₀  𝒋₁₁ ]   [ ∂x/∂η  ∂y/∂η ]
      //
      // The Jacobian is computed as follows:
      //   𝒋₀₀ = ∑ (∂𝑁ᵢ/∂ξ * xᵢ) ∀ 𝑖= 𝟏,……,𝟒
      //   𝒋₀₁ = ∑ (∂𝑁ᵢ/∂ξ * yᵢ) ∀ 𝑖= 𝟏,……,𝟒
      //   𝒋₁₀ = ∑ (∂𝑁ᵢ/∂η * yᵢ) ∀ 𝑖= 𝟏,……,𝟒
      //   𝒋₁₁ = ∑ (∂𝑁ᵢ/∂η * xᵢ) ∀ 𝑖= 𝟏,……,𝟒

      Real J00 = 0, J01 = 0, J10 = 0, J11 = 0;
      for (Int8 a = 0; a < 4; ++a) {
        J00 += dN_dxi[a] * m_node_coord[cell.nodeId(a)].x;
        J01 += dN_dxi[a] * m_node_coord[cell.nodeId(a)].y;
        J10 += dN_deta[a] * m_node_coord[cell.nodeId(a)].x;
        J11 += dN_deta[a] * m_node_coord[cell.nodeId(a)].y;
      }

      // Determinant of the Jacobian
      Real detJ = J00 * J11 - J01 * J10;

      if (detJ <= 0.0) {
        ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
      }

      // Inverse of the Jacobian
      //    𝑱⁻¹ = [ invJ00 invJ01 ]
      //          [ invJ10 invJ11 ]
      Real invJ00 = J11 / detJ, invJ01 = -J01 / detJ;
      Real invJ10 = -J10 / detJ, invJ11 = J00 / detJ;

      // Gradients in physical space at this Gauss point
      Real4 dxU, dyU;
      for (Int8 a = 0; a < 4; ++a) {
        dxU[a] = invJ00 * dN_dxi[a] + invJ01 * dN_deta[a];
        dyU[a] = invJ10 * dN_dxi[a] + invJ11 * dN_deta[a];
      }

      // stiffness matrix assembly
      ae += (dxU ^ dxU) * detJ * w * w + (dyU ^ dyU) * detJ * w * w;
    }
  }
  return ae;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a hexahedral element (HEXA8, ℙ1 FE).
 * This function calculates the integral of:
 *       a(u,v) = ∫∫∫ (∂u/∂x ∂v/∂x + ∂u/∂y ∂v/∂y + ∂u/∂z ∂v/∂z)dΩ
 * steps involved:
 * 1. Define Gauss points (2x2x2) and weights.
 * 2. Loop over Gauss points to compute the shape function derivatives.
 * 3. Compute the Jacobian matrix and its determinant.
 * 4. Compute the inverse of the Jacobian.
 * 5. Compute the gradients in physical space.
 * 6. Assemble the element matrix using the computed gradients.
 * @param cell The cell for which the element matrix is computed.
 * @return The computed element matrix.
 */
/*---------------------------------------------------------------------------*/

RealMatrix<8, 8> FemModule::_computeElementMatrixHexa8(Cell cell)
{
  // 2x2x2 Gauss points and weights for [-1,1]^3
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
  constexpr Real w = 1.0;

  RealMatrix<8, 8> ae;
  ae.fill(0.0); // Zero initialize

  for (Int8 ixi = 0; ixi < 2; ++ixi) {
    for (Int8 ieta = 0; ieta < 2; ++ieta) {
      for (Int8 izeta = 0; izeta < 2; ++izeta) {
        Real xi = gp[ixi];
        Real eta = gp[ieta];
        Real zeta = gp[izeta];

        // Shape function derivatives in reference space
        Real dN_dxi[8], dN_deta[8], dN_dzeta[8];
        dN_dxi[0] = -0.125 * (1 - eta) * (1 - zeta);
        dN_dxi[1] = 0.125 * (1 - eta) * (1 - zeta);
        dN_dxi[2] = 0.125 * (1 + eta) * (1 - zeta);
        dN_dxi[3] = -0.125 * (1 + eta) * (1 - zeta);
        dN_dxi[4] = -0.125 * (1 - eta) * (1 + zeta);
        dN_dxi[5] = 0.125 * (1 - eta) * (1 + zeta);
        dN_dxi[6] = 0.125 * (1 + eta) * (1 + zeta);
        dN_dxi[7] = -0.125 * (1 + eta) * (1 + zeta);

        dN_deta[0] = -0.125 * (1 - xi) * (1 - zeta);
        dN_deta[1] = -0.125 * (1 + xi) * (1 - zeta);
        dN_deta[2] = 0.125 * (1 + xi) * (1 - zeta);
        dN_deta[3] = 0.125 * (1 - xi) * (1 - zeta);
        dN_deta[4] = -0.125 * (1 - xi) * (1 + zeta);
        dN_deta[5] = -0.125 * (1 + xi) * (1 + zeta);
        dN_deta[6] = 0.125 * (1 + xi) * (1 + zeta);
        dN_deta[7] = 0.125 * (1 - xi) * (1 + zeta);

        dN_dzeta[0] = -0.125 * (1 - xi) * (1 - eta);
        dN_dzeta[1] = -0.125 * (1 + xi) * (1 - eta);
        dN_dzeta[2] = -0.125 * (1 + xi) * (1 + eta);
        dN_dzeta[3] = -0.125 * (1 - xi) * (1 + eta);
        dN_dzeta[4] = 0.125 * (1 - xi) * (1 - eta);
        dN_dzeta[5] = 0.125 * (1 + xi) * (1 - eta);
        dN_dzeta[6] = 0.125 * (1 + xi) * (1 + eta);
        dN_dzeta[7] = 0.125 * (1 - xi) * (1 + eta);

        // Jacobian matrix (default-initialized to zero see Real3x3.h)
        //    𝑱 = [ 𝒋₀₀  𝒋₀₁  𝒋₀₂ ]
        //        [ 𝒋₁₀  𝒋₁₁  𝒋₁₂ ]
        //        [ 𝒋₂₀  𝒋₂₁  𝒋₂₂ ]
        //
        // The Jacobian is computed as follows:
        //   𝒋₀₀ = ∑ (∂x/∂ξ * xᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₀₁ = ∑ (∂x/∂ξ * yᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₀₂ = ∑ (∂x/∂ξ * zᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₁₀ = ∑ (∂y/∂η * xᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₁₁ = ∑ (∂y/∂η * yᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₁₂ = ∑ (∂y/∂η * zᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₂₀ = ∑ (∂z/∂ζ * xᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₂₁ = ∑ (∂z/∂ζ * yᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₂₂ = ∑ (∂z/∂ζ * zᵢ) ∀ 𝑖= 𝟏,……,𝟖
        Real3x3 J;
        for (Int8 a = 0; a < 8; ++a) {
          const Real3& n = m_node_coord[cell.nodeId(a)];
          J[0][0] += dN_dxi[a] * n.x; // ∂x/∂ξ
          J[0][1] += dN_dxi[a] * n.y; // ∂y/∂ξ
          J[0][2] += dN_dxi[a] * n.z; // ∂z/∂ξ
          J[1][0] += dN_deta[a] * n.x; // ∂x/∂η
          J[1][1] += dN_deta[a] * n.y; // ∂y/∂η
          J[1][2] += dN_deta[a] * n.z; // ∂z/∂η
          J[2][0] += dN_dzeta[a] * n.x; // ∂x/∂ζ
          J[2][1] += dN_dzeta[a] * n.y; // ∂y/∂ζ
          J[2][2] += dN_dzeta[a] * n.z; // ∂z/∂ζ
        }

        // Determinant of the Jacobian
        Real detJ = math::matrixDeterminant(J);

        // Inverse Jacobian
        Real3x3 invJ;
        invJ = math::inverseMatrix(J, detJ);

        if (detJ <= 0.0) {
          ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
        }

        // Gradients in physical space
        Real dxU[8], dyU[8], dzU[8];
        for (Int8 a = 0; a < 8; ++a) {
          dxU[a] = invJ[0][0] * dN_dxi[a] + invJ[0][1] * dN_deta[a] + invJ[0][2] * dN_dzeta[a];
          dyU[a] = invJ[1][0] * dN_dxi[a] + invJ[1][1] * dN_deta[a] + invJ[1][2] * dN_dzeta[a];
          dzU[a] = invJ[2][0] * dN_dxi[a] + invJ[2][1] * dN_deta[a] + invJ[2][2] * dN_dzeta[a];
        }

        // Assemble element matrix (variational form)
        for (Int8 i = 0; i < 8; ++i) {
          for (Int8 j = 0; j < 8; ++j) {
            ae(i, j) += (dxU[i] * dxU[j] + dyU[i] * dyU[j] + dzU[i] * dzU[j]) * detJ * w * w * w;
          }
        }
      }
    }
  }
  return ae;
}