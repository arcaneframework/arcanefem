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

    if (mesh()->dimension() == 3) {
      if (m_hex_quad_mesh) {
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;

          // 2x2 Gauss integration for quadrilateral face
          constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
          constexpr Real w = 1.0;

          // Get face nodes (assuming quad4 face)
          Node node0 = face.node(0);
          Node node1 = face.node(1);
          Node node2 = face.node(2);
          Node node3 = face.node(3);
          Node nodes[4] = { node0, node1, node2, node3 };

          // Get node coordinates
          Real3 coords[4];
          for (Int32 i = 0; i < 4; ++i) {
            coords[i] = m_node_coord[nodes[i]];
          }

          // Loop through 2x2 Gauss points
          for (Int32 ixi = 0; ixi < 2; ++ixi) {
            for (Int32 ieta = 0; ieta < 2; ++ieta) {
              Real xi = gp[ixi];
              Real eta = gp[ieta];

              // Quad4 shape functions
              Real N[4];
              N[0] = 0.25 * (1 - xi) * (1 - eta);
              N[1] = 0.25 * (1 + xi) * (1 - eta);
              N[2] = 0.25 * (1 + xi) * (1 + eta);
              N[3] = 0.25 * (1 - xi) * (1 + eta);

              // Shape function derivatives w.r.t. natural coordinates
              Real dN_dxi[4], dN_deta[4];
              dN_dxi[0] = -0.25 * (1 - eta);
              dN_dxi[1] = 0.25 * (1 - eta);
              dN_dxi[2] = 0.25 * (1 + eta);
              dN_dxi[3] = -0.25 * (1 + eta);

              dN_deta[0] = -0.25 * (1 - xi);
              dN_deta[1] = -0.25 * (1 + xi);
              dN_deta[2] = 0.25 * (1 + xi);
              dN_deta[3] = 0.25 * (1 - xi);

              // Compute tangent vectors
              Real3 t1(0.0, 0.0, 0.0); // ∂r/∂ξ
              Real3 t2(0.0, 0.0, 0.0); // ∂r/∂η

              for (Int32 i = 0; i < 4; ++i) {
                t1.x += dN_dxi[i] * coords[i].x;
                t1.y += dN_dxi[i] * coords[i].y;
                t1.z += dN_dxi[i] * coords[i].z;

                t2.x += dN_deta[i] * coords[i].x;
                t2.y += dN_deta[i] * coords[i].y;
                t2.z += dN_deta[i] * coords[i].z;
              }

              // Normal vector (cross product of tangent vectors)
              Real3 normal;
              normal.x = t1.y * t2.z - t1.z * t2.y;
              normal.y = t1.z * t2.x - t1.x * t2.z;
              normal.z = t1.x * t2.y - t1.y * t2.x;

              // Jacobian (magnitude of normal vector for surface integration)
              Real detJ = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

              // Integration weight
              Real integration_weight = w * w * detJ;

              // Apply to all four nodes of the face
              for (Int32 j = 0; j < 4; ++j) {
                Node node = nodes[j];
                if (!node.isOwn())
                  continue;

                rhs_values[node_dof.dofId(node, 0)] += t[0] * N[j] * integration_weight;
                rhs_values[node_dof.dofId(node, 1)] += t[1] * N[j] * integration_weight;
                rhs_values[node_dof.dofId(node, 2)] += t[2] * N[j] * integration_weight;
              }
            }
          }
        }
      }
      else {

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
  }
}