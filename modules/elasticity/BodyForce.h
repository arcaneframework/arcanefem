// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BodyForce.h                                                 (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute and assemble body force contribution to RHS */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies body force to the RHS vector of the linear system.
 * 
 * This function computes the contribution of body forces to the RHS vector 
 * of the linear system. It iterates over all cells in the mesh, calculates 
 * the appropriate force  contributions  based on the element type and mesh 
 * dimension, and updates the RHS vector accordingly.
 * 
 * body force ∫∫∫ (𝐟.𝐯)  with 𝐟 = (𝑓𝑥, 𝑓𝑦, 𝑓𝑧) = (f[0], f[1], f[2])
 * 
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding 
 *                 degrees of freedom (DoFs).
 * 
/*---------------------------------------------------------------------------*/

inline void FemModule::
_applyBodyForce(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{

  // get bodyforce vector
  bool applyBodyForce = false;
  const UniqueArray<String> f_string = options()->f();
  info() << "[ArcaneFem-Info] Applying Bodyforce " << f_string;
  for (Int32 i = 0; i < f_string.size(); ++i) {
    f[i] = 0.0;
    if (f_string[i] != "NULL") {
      applyBodyForce = true;
      f[i] = std::stod(f_string[i].localstr());
    }
  }

  // no bodyforce to apply hence return
  if (!applyBodyForce)
    return;

  // apply bodyforce based on dimension and mesh type
  if (mesh()->dimension() == 2) {
    if (m_hex_quad_mesh) {
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;

        constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; //(ξ,η)
        constexpr Real weights[2] = { 1.0, 1.0 };

        for (Int32 ixi = 0; ixi < 2; ++ixi) {
          for (Int32 ieta = 0; ieta < 2; ++ieta) {

            // set the coordinates of the Gauss point
            Real xi = gp[ixi]; // Get the ξ coordinate of the Gauss point
            Real eta = gp[ieta]; // Get the η coordinate of the Gauss point

            // set the weight of the Gauss point
            Real weight = weights[ixi] * weights[ieta];

            // get shape functions 𝐍 for Quad4: 𝐍(ξ,η) = [𝑁₁  𝑁₂  𝑁₃  𝑁₄]
            RealVector<4> N = ArcaneFemFunctions::FeOperation2D::computeShapeFunctionsQuad4(xi, eta);

            // get determinant of Jacobian
            const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, xi, eta);
            const Real detJ = gp_info.det_j;

            // compute integration weight
            Real integration_weight = weight * detJ;

            // Assemble RHS
            for (Int32 i = 0; i < 4; ++i) {
              Node node = cell.node(i);
              if (node.isOwn()) {
                rhs_values[node_dof.dofId(node, 0)] += N[i] * f[0] * integration_weight;
                rhs_values[node_dof.dofId(node, 1)] += N[i] * f[1] * integration_weight;
              }
            }
          }
        }
      }
    }
    else {
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
        for (Node node : cell.nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += f[0] * area / 3;
            rhs_values[node_dof.dofId(node, 1)] += f[1] * area / 3;
          }
        }
      }
    }
  }
  if (mesh()->dimension() == 3) {
    if (m_hex_quad_mesh) {
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;

        constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // [-1/sqrt(3), 1/sqrt(3)]
        constexpr Real weights[2] = { 1.0, 1.0 };

        for (Int32 ixi = 0; ixi < 2; ++ixi) {
          for (Int32 ieta = 0; ieta < 2; ++ieta) {
            for (Int32 izeta = 0; izeta < 2; ++izeta) {

              // set the coordinates of the Gauss point
              Real xi = gp[ixi];
              Real eta = gp[ieta];
              Real zeta = gp[izeta];

              // set the weight of the Gauss point
              Real weight = weights[ixi] * weights[ieta] * weights[izeta];

              // Get shape functions for Hexa8: 𝐍(ξ,η,ζ) = [𝑁₁  𝑁₂  𝑁₃  𝑁₄  𝑁₅  𝑁₆  𝑁₇  𝑁₈]
              RealVector<8> N = ArcaneFemFunctions::FeOperation3D::computeShapeFunctionsHexa8(xi, eta, zeta);

              // get determinant of Jacobian
              const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa8(cell, m_node_coord, xi, eta, zeta);
              const Real detJ = gp_info.det_j;

              // compute integration weight
              Real integration_weight = detJ * weight;

              // Assemble RHS
              for (Int32 i = 0; i < 8; ++i) {
                Node node = cell.node(i);
                if (node.isOwn()) {
                  rhs_values[node_dof.dofId(node, 0)] += N[i] * f[0] * integration_weight;
                  rhs_values[node_dof.dofId(node, 1)] += N[i] * f[1] * integration_weight;
                  rhs_values[node_dof.dofId(node, 2)] += N[i] * f[2] * integration_weight;
                }
              }
            }
          }
        }
      }
    }
    else {
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);
        for (Node node : cell.nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += f[0] * volume / 4;
            rhs_values[node_dof.dofId(node, 1)] += f[1] * volume / 4;
            rhs_values[node_dof.dofId(node, 2)] += f[2] * volume / 4;
          }
        }
      }
    }
  }
}