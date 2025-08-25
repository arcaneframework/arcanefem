// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Traction.h                                                  (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute and assemble traction contribution to RHS   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies traction to the RHS vector of the linear system.
 * 
 * This function  computes  the  contribution of traction to the RHS vector 
 * of the linear system. It iterates over all cells in the mesh, calculates 
 * the appropriate force  contributions  based on the element type and mesh 
 * dimension, and updates the RHS vector accordingly.
 * 
 * traction term ∫∫ (𝐭.𝐯)  with 𝐭 = (𝑡𝑥, 𝑡𝑦, 𝑡𝑧) = (t[0], t[1], t[2])
 * 
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding 
 *                 degrees of freedom (DoFs).
 * 
/*---------------------------------------------------------------------------*/

inline void FemModule::
_applyTractionToRHS(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  // loop over all traction boundries
  for (const auto& bs : options()->tractionBoundaryCondition()) {

    // mesh boundary group on which traction is applied
    FaceGroup group = bs->surface();

    // get traction force vector
    bool applyTraction = false;
    const UniqueArray<String> t_string = bs->t();
    for (Int32 i = 0; i < t_string.size(); ++i) {
      t[i] = 0.0;
      if (t_string[i] != "NULL") {
        applyTraction = true;
        t[i] = std::stod(t_string[i].localstr());
      }
    }

    // no traction to apply hence return
    if (!applyTraction)
      return;

    // print traction info
    info() << "[ArcaneFem-Info] Applying Traction " << t_string;
    info() << "[ArcaneFem-Info] Traction surface '" << bs->surface().name() << "'";

    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {

        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;

          // 2-point Gauss integration for line element
          constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
          constexpr Real weights[2] = { 1.0, 1.0 };

          Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);
          Real2 normal = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, m_node_coord);

          Node node0 = face.node(0);
          Node node1 = face.node(1);

          for (Int32 i = 0; i < 2; ++i) {
            Real xi = gp[i];
            Real weight = weights[i];

            // Linear shape functions for Line2
            Real N[2];
            N[0] = 0.5 * (1 - xi);
            N[1] = 0.5 * (1 + xi);

            // Integration weight: weight * jacobian (length/2 for reference element [-1,1])
            Real integration_weight = weight * length * 0.5;

            // Apply to both nodes
            Node nodes[2] = { node0, node1 };
            for (Int32 j = 0; j < 2; ++j) {
              Node node = nodes[j];
              if (!node.isOwn())
                continue;

              rhs_values[node_dof.dofId(node, 0)] += t[0] * N[j] * integration_weight;
              rhs_values[node_dof.dofId(node, 1)] += t[1] * N[j] * integration_weight;
            }
          }
        }
      }
      else {
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);
          for (Node node : iface->nodes()) {
            if (node.isOwn()) {
              rhs_values[node_dof.dofId(node, 0)] += t[0] * length / 2.;
              rhs_values[node_dof.dofId(node, 1)] += t[1] * length / 2.;
            }
          }
        }
      }
    }

    if (mesh()->dimension() == 3)
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, m_node_coord);
        for (Node node : iface->nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += t[0] * area / 3.;
            rhs_values[node_dof.dofId(node, 1)] += t[1] * area / 3.;
            rhs_values[node_dof.dofId(node, 2)] += t[2] * area / 3.;
          }
        }
      }
  }
}