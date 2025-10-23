// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* SourceTerm.h                                                (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute and assemble source term contribution to RHS*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies source term to RHS vector of the linear system.
 * 
 * This function computes the contribution of the source term to the 
 * right-hand side (RHS) vector of the linear system based on the mesh
 * dimension. 
 * 
 * source term = ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯)
 * 
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_applySourceTerm(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  if (mesh()->dimension() == 2)
    if (m_hex_quad_mesh)
      _applySourceTermQuad4(rhs_values, node_dof);
    else
      _applySourceTermTria3(rhs_values, node_dof);

  if (mesh()->dimension() == 3)
    if (m_hex_quad_mesh)
      _applySourceTermHexa8(rhs_values, node_dof);
    else
      _applySourceTermTetra4(rhs_values, node_dof);
}


void FemModule::
_applySourceTermTria3(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{

  // Pre-compute basis vectors and mass matrix outside the loop
  RealVector<6> Uy = { 0., 1., 0., 1., 0., 1. };
  RealVector<6> Ux = { 1., 0., 1., 0., 1., 0. };
  RealVector<6> F = { f[0], f[1], f[0], f[1], f[0], f[1] };

  RealMatrix<6, 6> massMatrixXplusY = massMatrix(Ux, Ux) + massMatrix(Uy, Uy);

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

    RealVector<6> Un = { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y,
                         m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y,
                         m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y };

    RealVector<6> Vn = { m_V[cell.nodeId(0)].x, m_V[cell.nodeId(0)].y,
                         m_V[cell.nodeId(1)].x, m_V[cell.nodeId(1)].y,
                         m_V[cell.nodeId(2)].x, m_V[cell.nodeId(2)].y };

    RealVector<6> An = { m_A[cell.nodeId(0)].x, m_A[cell.nodeId(0)].y,
                         m_A[cell.nodeId(1)].x, m_A[cell.nodeId(1)].y,
                         m_A[cell.nodeId(2)].x, m_A[cell.nodeId(2)].y };

    //  ∫∫ (𝐟.𝐯) + ∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯)
    RealVector<6> rhs = ( F * (1/3.)
                            + Un * (massMatrixXplusY)*(c0*1/12.)
                            + Vn * (massMatrixXplusY)*(c3*1/12.)
                            + An * (massMatrixXplusY)*(c4*1/12.)
                            ) * area;

    rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
    rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
    rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(2);
    rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(3);
    rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(4);
    rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(5);
  }
}


void FemModule::
_applySourceTermQuad4(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    // Initialize RHS contributions (2 dof/node for 4 quad nodes)
    Real rhs_x_contributions[4] = { 0., 0., 0., 0. };
    Real rhs_y_contributions[4] = { 0., 0., 0., 0. };

    // 2x2 Gauss integration for quadrilateral element
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 };
    constexpr Real w = 1.0;

    for (Int8 ixi = 0; ixi < 2; ++ixi) {
      for (Int8 ieta = 0; ieta < 2; ++ieta) {

        // Get the coordinates of the Gauss point
        Real xi = gp[ixi]; // Get the ξ
        Real eta = gp[ieta]; // Get the η
        Real weight = w * w; // Weight

        // Shape functions  𝐍 for Quad4
        RealVector<4> N = ArcaneFemFunctions::FeOperation2D::computeShapeFunctionsQuad4(xi, eta);

        // compute the det(Jacobian)
        auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, xi, eta);
        const Real detJ = gp_info.det_j;

        // compute integration weight
        Real integration_weight = weight * detJ;

        // Interpolate fields (𝐮ₙ,𝐮ᵗₙ,𝐮ᵗᵗ) at the quadrature point: (.)_gp = ∑ 𝑁ᵢ * (.)
        Real2 u_gp = { 0, 0 };
        Real2 v_gp = { 0, 0 };
        Real2 a_gp = { 0, 0 };
        for (Int8 a = 0; a < 4; ++a) {
          u_gp.x += N[a] * m_U[cell.nodeId(a)].x;
          u_gp.y += N[a] * m_U[cell.nodeId(a)].y;

          v_gp.x += N[a] * m_V[cell.nodeId(a)].x;
          v_gp.y += N[a] * m_V[cell.nodeId(a)].y;

          a_gp.x += N[a] * m_A[cell.nodeId(a)].x;
          a_gp.y += N[a] * m_A[cell.nodeId(a)].y;
        }

        // Source force term 𝐟
        Real2 f_gp(f[0], f[1]);

        // ∫∫ (𝐟.𝐯) + ∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯)
        for (Int8 a = 0; a < 4; ++a) {
          rhs_x_contributions[a] += (f_gp.x + c0 * u_gp.x + c3 * v_gp.x + c4 * a_gp.x) * N[a] * integration_weight;
          rhs_y_contributions[a] += (f_gp.y + c0 * u_gp.y + c3 * v_gp.y + c4 * a_gp.y) * N[a] * integration_weight;
        }
      }
    }

    // Add contributions to global RHS
    for (Int8 a = 0; a < 4; ++a) {
      Node node = cell.node(a);
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += rhs_x_contributions[a];
        rhs_values[node_dof.dofId(node, 1)] += rhs_y_contributions[a];
      }
    }
  }
}

void FemModule::
_applySourceTermTetra4(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{

  // Pre-compute basis vectors and mass matrix outside the loop
  RealVector<12> Uy = { 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0. };
  RealVector<12> Ux = { 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0. };
  RealVector<12> Uz = { 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1. };
  RealVector<12> F = { f[0], f[1], f[2], f[0], f[1], f[2], f[0], f[1], f[2], f[0], f[1], f[2] };

  RealMatrix<12, 12> massMatrixXplusYplusZ = massMatrix(Ux, Ux) + massMatrix(Uy, Uy) + massMatrix(Uz, Uz);

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);

    RealVector<12> Un = { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y, m_U[cell.nodeId(0)].z,
                          m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y, m_U[cell.nodeId(1)].z,
                          m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y, m_U[cell.nodeId(2)].z,
                          m_U[cell.nodeId(3)].x, m_U[cell.nodeId(3)].y, m_U[cell.nodeId(3)].z };

    RealVector<12> Vn = { m_V[cell.nodeId(0)].x, m_V[cell.nodeId(0)].y, m_V[cell.nodeId(0)].z,
                          m_V[cell.nodeId(1)].x, m_V[cell.nodeId(1)].y, m_V[cell.nodeId(1)].z,
                          m_V[cell.nodeId(2)].x, m_V[cell.nodeId(2)].y, m_V[cell.nodeId(2)].z,
                          m_V[cell.nodeId(3)].x, m_V[cell.nodeId(3)].y, m_V[cell.nodeId(3)].z };

    RealVector<12> An = { m_A[cell.nodeId(0)].x, m_A[cell.nodeId(0)].y, m_A[cell.nodeId(0)].z,
                          m_A[cell.nodeId(1)].x, m_A[cell.nodeId(1)].y, m_A[cell.nodeId(1)].z,
                          m_A[cell.nodeId(2)].x, m_A[cell.nodeId(2)].y, m_A[cell.nodeId(2)].z,
                          m_A[cell.nodeId(3)].x, m_A[cell.nodeId(3)].y, m_A[cell.nodeId(3)].z };

    //  ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯) +
    RealVector<12> rhs = ( F * (1/4.)
                            + Un * (massMatrixXplusYplusZ)*(c0*1/20.)
                            + Vn * (massMatrixXplusYplusZ)*(c3*1/20.)
                            + An * (massMatrixXplusYplusZ)*(c4*1/20.)
                            ) * volume;

    rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
    rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
    rhs_values[node_dof.dofId(cell.nodeId(0), 2)] += rhs(2);
    rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(3);
    rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(4);
    rhs_values[node_dof.dofId(cell.nodeId(1), 2)] += rhs(5);
    rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(6);
    rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(7);
    rhs_values[node_dof.dofId(cell.nodeId(2), 2)] += rhs(8);
    rhs_values[node_dof.dofId(cell.nodeId(3), 0)] += rhs(9);
    rhs_values[node_dof.dofId(cell.nodeId(3), 1)] += rhs(10);
    rhs_values[node_dof.dofId(cell.nodeId(3), 2)] += rhs(11);
  }
}

void FemModule::
_applySourceTermHexa8(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    // Initialize RHS contributions (2 dof/node for 8 hexa nodes)
    Real rhs_x_contributions[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };
    Real rhs_y_contributions[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };
    Real rhs_z_contributions[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };

    // 2x2 Gauss integration for quadrilateral element
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 };
    constexpr Real w = 1.0;

    for (Int8 ixi = 0; ixi < 2; ++ixi) {
      for (Int8 ieta = 0; ieta < 2; ++ieta) {
        for (Int8 izeta = 0; izeta < 2; ++izeta) {

          // Get the coordinates of the Gauss point
          Real xi = gp[ixi]; // Get the ξ
          Real eta = gp[ieta]; // Get the η
          Real zeta = gp[izeta]; // Get the ζ
          Real weight = w * w * w; // Weight

          // Shape functions  𝐍 for Hexa8
          RealVector<8> N = ArcaneFemFunctions::FeOperation3D::computeShapeFunctionsHexa8(xi, eta, zeta);

          // compute the det(Jacobian)
          auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa8(cell, m_node_coord, xi, eta, zeta);
          const Real detJ = gp_info.det_j;

          // compute integration weight
          Real integration_weight = weight * detJ;

          // Interpolate fields (𝐮ₙ,𝐮ᵗₙ,𝐮ᵗᵗ) at the quadrature point: (.)_gp = ∑ 𝑁ᵢ * (.)
          Real3 u_gp = { 0, 0, 0 };
          Real3 v_gp = { 0, 0, 0 };
          Real3 a_gp = { 0, 0, 0 };
          for (Int8 a = 0; a < 8; ++a) {
            u_gp.x += N[a] * m_U[cell.nodeId(a)].x;
            u_gp.y += N[a] * m_U[cell.nodeId(a)].y;
            u_gp.z += N[a] * m_U[cell.nodeId(a)].z;

            v_gp.x += N[a] * m_V[cell.nodeId(a)].x;
            v_gp.y += N[a] * m_V[cell.nodeId(a)].y;
            v_gp.z += N[a] * m_V[cell.nodeId(a)].z;

            a_gp.x += N[a] * m_A[cell.nodeId(a)].x;
            a_gp.y += N[a] * m_A[cell.nodeId(a)].y;
            a_gp.z += N[a] * m_A[cell.nodeId(a)].z;
          }

          // Source force term 𝐟
          Real3 f_gp(f[0], f[1], f[2]);

          // ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯)
          for (Int8 a = 0; a < 8; ++a) {
            rhs_x_contributions[a] += (f_gp.x + c0 * u_gp.x + c3 * v_gp.x + c4 * a_gp.x) * N[a] * integration_weight;
            rhs_y_contributions[a] += (f_gp.y + c0 * u_gp.y + c3 * v_gp.y + c4 * a_gp.y) * N[a] * integration_weight;
            rhs_z_contributions[a] += (f_gp.z + c0 * u_gp.z + c3 * v_gp.z + c4 * a_gp.z) * N[a] * integration_weight;
          }
        }
      }
    }
    // Add contributions to global RHS
    for (Int8 a = 0; a < 8; ++a) {
      Node node = cell.node(a);
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += rhs_x_contributions[a];
        rhs_values[node_dof.dofId(node, 1)] += rhs_y_contributions[a];
        rhs_values[node_dof.dofId(node, 2)] += rhs_z_contributions[a];
      }
    }
  }
}