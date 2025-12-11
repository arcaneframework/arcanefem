// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Paraxial.h                                                   (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute & assemble paraxial contribution to LHS/RHS */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes 2D paraxial element matrix for a edge element (ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) =  ∫ (cₚ 𝑁𝑥² + cₛ 𝑁𝑦²)(𝑢𝑥 𝑣𝑥) +
 *             ∫ (cₚ 𝑁𝑦² + cₛ 𝑁𝑥²)(𝑢𝑦 𝑣𝑦) +
 *             ∫ (𝑁𝑦 𝑁𝑥 (cₚ - cₛ))(𝑢𝑥 𝑣𝑦) +
 *             ∫ (𝑁𝑥 𝑁𝑦 (cₚ - cₛ))(𝑢𝑦 𝑣𝑥) ;
 *
 *   with  trial func 𝐮 = (𝑢𝑥,𝑢𝑦) and test func 𝐯 = (𝑣𝑥,𝑣𝑦)
 */
/*---------------------------------------------------------------------------*/

RealMatrix<4, 4> FemModuleSoildynamics::
_computeParaxialElementMatrixEdge2(Face face)
{
  Real2 N = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, m_node_coord);

  RealVector<4> Uy = { 0., 1., 0., 1. };
  RealVector<4> Ux = { 1., 0., 1., 0. };

  RealMatrix<4, 4> int_Omega_i = (((N.x * N.x * cp + N.y * N.y * cs)) * (massMatrix(Ux, Ux)) +
                                  ((N.y * N.y * cp + N.x * N.x * cs)) * (massMatrix(Uy, Uy)) +
                                  ((N.x * N.y * (cp - cs))) * (massMatrix(Ux, Uy)) +
                                  ((N.x * N.y * (cp - cs))) * (massMatrix(Uy, Ux))) / 6.;
  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes 3D paraxial element matrix for a triangular element (ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) =  ∫∫ (cₚ 𝑁𝑥² + cₛ (1 - 𝑁𝑥²))(𝑢𝑥 𝑣𝑥) +
 *             ∫∫ (cₚ 𝑁𝑦² + cₛ (1 - 𝑁𝑦²))(𝑢𝑦 𝑣𝑦) +
 *             ∫∫ (cₚ 𝑁𝑧² + cₛ (1 - 𝑁𝑧²))(𝑢𝑧 𝑣𝑧) +
 *             ∫∫ (𝑁𝑥 𝑁𝑦 (cₚ - cₛ))(𝑢𝑥 𝑣𝑦) + ∫∫ (c₇)(𝑁𝑥 𝑁𝑧 (cₚ - cₛ))(𝑢𝑥 𝑣𝑧) +
 *             ∫∫ (𝑁𝑦 𝑁𝑥 (cₚ - cₛ))(𝑢𝑦 𝑣𝑥) + ∫∫ (c₇)(𝑁𝑦 𝑁𝑧 (cₚ - cₛ))(𝑢𝑦 𝑣𝑧) +
 *             ∫∫ (𝑁𝑧 𝑁𝑥 (cₚ - cₛ))(𝑢𝑧 𝑣𝑥) + ∫∫ (c₇)(𝑁𝑧 𝑁𝑦 (cₚ - cₛ))(𝑢𝑧 𝑣𝑦) ;
 *
 *   with trial function 𝐮 = (𝑢𝑥, 𝑢𝑦, 𝑢𝑧) and test function 𝐯 = (𝑣𝑥, 𝑣𝑦, 𝑣𝑧)
 */
/*---------------------------------------------------------------------------*/

RealMatrix<9, 9> FemModuleSoildynamics::
_computeParaxialElementMatrixTria3(Face face)
{
  Real3 N = ArcaneFemFunctions::MeshOperation::computeNormalTriangle(face, m_node_coord);

  RealVector<9> Ux = { 1., 0., 0., 1., 0., 0., 1., 0., 0. };
  RealVector<9> Uy = { 0., 1., 0., 0., 1., 0., 0., 1., 0. };
  RealVector<9> Uz = { 0., 0., 1., 0., 0., 1., 0., 0., 1. };

  RealMatrix<9, 9> int_Omega_i = (((N.x * N.x * cp + (1. - N.x * N.x) * cs)) * (massMatrix(Ux, Ux)) +
                                  ((N.y * N.y * cp + (1. - N.y * N.y) * cs)) * (massMatrix(Uy, Uy)) +
                                  ((N.z * N.z * cp + (1. - N.z * N.z) * cs)) * (massMatrix(Uz, Uz)) +
                                  ((N.x * N.y * (cp - cs))) * (massMatrix(Ux, Uy)) +
                                  ((N.x * N.z * (cp - cs))) * (massMatrix(Ux, Uz)) +
                                  ((N.y * N.x * (cp - cs))) * (massMatrix(Uy, Ux)) +
                                  ((N.y * N.z * (cp - cs))) * (massMatrix(Uy, Uz)) +
                                  ((N.z * N.x * (cp - cs))) * (massMatrix(Uz, Ux)) +
                                  ((N.z * N.y * (cp - cs))) * (massMatrix(Uz, Uy))) / 12.;
  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes 3D paraxial element matrix for a quadrangular element (ℙ1 FE).
 *
 * Theory:
 *
 *   a(𝐮,𝐯) =  ∫∫ (cₚ 𝑁𝑥² + cₛ (1 - 𝑁𝑥²))(𝑢𝑥 𝑣𝑥) +
 *             ∫∫ (cₚ 𝑁𝑦² + cₛ (1 - 𝑁𝑦²))(𝑢𝑦 𝑣𝑦) +
 *             ∫∫ (cₚ 𝑁𝑧² + cₛ (1 - 𝑁𝑧²))(𝑢𝑧 𝑣𝑧) +
 *             ∫∫ (𝑁𝑥 𝑁𝑦 (cₚ - cₛ))(𝑢𝑥 𝑣𝑦) + ∫∫ (c₇)(𝑁𝑥 𝑁𝑧 (cₚ - cₛ))(𝑢𝑥 𝑣𝑧) +
 *             ∫∫ (𝑁𝑦 𝑁𝑥 (cₚ - cₛ))(𝑢𝑦 𝑣𝑥) + ∫∫ (c₇)(𝑁𝑦 𝑁𝑧 (cₚ - cₛ))(𝑢𝑦 𝑣𝑧) +
 *             ∫∫ (𝑁𝑧 𝑁𝑥 (cₚ - cₛ))(𝑢𝑧 𝑣𝑥) + ∫∫ (c₇)(𝑁𝑧 𝑁𝑦 (cₚ - cₛ))(𝑢𝑧 𝑣𝑦) ;
 *
 *   with trial function 𝐮 = (𝑢𝑥, 𝑢𝑦, 𝑢𝑧) and test function 𝐯 = (𝑣𝑥, 𝑣𝑦, 𝑣𝑧)
 */
/*---------------------------------------------------------------------------*/

RealMatrix<12, 12> FemModuleSoildynamics::
_computeParaxialElementMatrixQuad4(Face face)
{
  Real3 N = ArcaneFemFunctions::MeshOperation::computeNormalQuad(face, m_node_coord);

  RealVector<12> Ux = { 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0. };
  RealVector<12> Uy = { 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0. };
  RealVector<12> Uz = { 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1. };

  RealMatrix<12, 12> int_Omega_i = (((N.x * N.x * cp + (1. - N.x * N.x) * cs)) * (massMatrix(Ux, Ux)) +
                                    ((N.y * N.y * cp + (1. - N.y * N.y) * cs)) * (massMatrix(Uy, Uy)) +
                                    ((N.z * N.z * cp + (1. - N.z * N.z) * cs)) * (massMatrix(Uz, Uz)) +
                                    ((N.x * N.y * (cp - cs))) * (massMatrix(Ux, Uy)) +
                                    ((N.x * N.z * (cp - cs))) * (massMatrix(Ux, Uz)) +
                                    ((N.y * N.x * (cp - cs))) * (massMatrix(Uy, Ux)) +
                                    ((N.y * N.z * (cp - cs))) * (massMatrix(Uy, Uz)) +
                                    ((N.z * N.x * (cp - cs))) * (massMatrix(Uz, Ux)) +
                                    ((N.z * N.y * (cp - cs))) * (massMatrix(Uz, Uy))) / 20.;
  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies paraxial boundary conditions to the RHS and LHS
 * of the linear system.
 * 
 * @param rhs_values The RHS values to be updated.
 * @param node_dof The node DoF connectivity view.
 * @param bsr_matrix The BSR matrix to be updated (optional).
 * 
/*---------------------------------------------------------------------------*/

void FemModuleSoildynamics::
_applyParaxial(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof, BSRMatrix* bsr_matrix)
{
  for (const auto& bs : options()->paraxialBoundaryCondition()) {
    FaceGroup group = bs->surface();

    info() << "[ArcaneFem-Info] Applying constant Paraxial boundary conditions for surface " << group.name();

    if (mesh()->dimension() == 2) {
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);

        RealVector<4> Un = { m_U[face.nodeId(0)].x, m_U[face.nodeId(0)].y, m_U[face.nodeId(1)].x, m_U[face.nodeId(1)].y };
        RealVector<4> Vn = { m_V[face.nodeId(0)].x, m_V[face.nodeId(0)].y, m_V[face.nodeId(1)].x, m_V[face.nodeId(1)].y };
        RealVector<4> An = { m_A[face.nodeId(0)].x, m_A[face.nodeId(0)].y, m_A[face.nodeId(1)].x, m_A[face.nodeId(1)].y };

        RealMatrix<4, 4> ParaxialElementMatrixEdge2 = _computeParaxialElementMatrixEdge2(face);
        RealVector<4> rhs = length * (c7 * Un * ParaxialElementMatrixEdge2 - c8 * Vn * ParaxialElementMatrixEdge2 + c9 * An * ParaxialElementMatrixEdge2);

        rhs_values[node_dof.dofId(face.nodeId(0), 0)] += rhs(0);
        rhs_values[node_dof.dofId(face.nodeId(0), 1)] += rhs(1);
        rhs_values[node_dof.dofId(face.nodeId(1), 0)] += rhs(2);
        rhs_values[node_dof.dofId(face.nodeId(1), 1)] += rhs(3);

        if (t <= dt) {
          ParaxialElementMatrixEdge2 = c7 * length * ParaxialElementMatrixEdge2;
          Int32 n1_index = 0;
          for (Node node1 : face.nodes()) {
            Int32 n2_index = 0;
            for (Node node2 : face.nodes()) {
              Real v1 = ParaxialElementMatrixEdge2(2 * n1_index, 2 * n2_index);
              Real v2 = ParaxialElementMatrixEdge2(2 * n1_index, 2 * n2_index + 1);
              Real v3 = ParaxialElementMatrixEdge2(2 * n1_index + 1, 2 * n2_index);
              Real v4 = ParaxialElementMatrixEdge2(2 * n1_index + 1, 2 * n2_index + 1);
              if (node1.isOwn()) {
                DoFLocalId node1_dof1 = node_dof.dofId(node1, 0);
                DoFLocalId node1_dof2 = node_dof.dofId(node1, 1);
                DoFLocalId node2_dof1 = node_dof.dofId(node2, 0);
                DoFLocalId node2_dof2 = node_dof.dofId(node2, 1);

                if (bsr_matrix) {
                  bsr_matrix->addValue(node1_dof1, node2_dof1, v1);
                  bsr_matrix->addValue(node1_dof1, node2_dof2, v2);
                  bsr_matrix->addValue(node1_dof2, node2_dof1, v3);
                  bsr_matrix->addValue(node1_dof2, node2_dof2, v4);
                }
                else {
                  m_linear_system.matrixAddValue(node1_dof1, node2_dof1, v1);
                  m_linear_system.matrixAddValue(node1_dof1, node2_dof2, v2);
                  m_linear_system.matrixAddValue(node1_dof2, node2_dof1, v3);
                  m_linear_system.matrixAddValue(node1_dof2, node2_dof2, v4);
                }
              }
              ++n2_index;
            }
            ++n1_index;
          }
        }
      }
    }

    if (mesh()->dimension() == 3) {
      if (m_hex_quad_mesh)
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real area = ArcaneFemFunctions::MeshOperation::computeAreaQuad4(face, m_node_coord);

          RealVector<12> Un = {
            m_U[face.nodeId(0)].x,
            m_U[face.nodeId(0)].y,
            m_U[face.nodeId(0)].z,
            m_U[face.nodeId(1)].x,
            m_U[face.nodeId(1)].y,
            m_U[face.nodeId(1)].z,
            m_U[face.nodeId(2)].x,
            m_U[face.nodeId(2)].y,
            m_U[face.nodeId(2)].z,
            m_U[face.nodeId(3)].x,
            m_U[face.nodeId(3)].y,
            m_U[face.nodeId(3)].z
          };

          RealVector<12> Vn = {
            m_V[face.nodeId(0)].x,
            m_V[face.nodeId(0)].y,
            m_V[face.nodeId(0)].z,
            m_V[face.nodeId(1)].x,
            m_V[face.nodeId(1)].y,
            m_V[face.nodeId(1)].z,
            m_V[face.nodeId(2)].x,
            m_V[face.nodeId(2)].y,
            m_V[face.nodeId(2)].z,
            m_V[face.nodeId(3)].x,
            m_V[face.nodeId(3)].y,
            m_V[face.nodeId(3)].z
          };

          RealVector<12> An = {
            m_A[face.nodeId(0)].x,
            m_A[face.nodeId(0)].y,
            m_A[face.nodeId(0)].z,
            m_A[face.nodeId(1)].x,
            m_A[face.nodeId(1)].y,
            m_A[face.nodeId(1)].z,
            m_A[face.nodeId(2)].x,
            m_A[face.nodeId(2)].y,
            m_A[face.nodeId(2)].z,
            m_A[face.nodeId(3)].x,
            m_A[face.nodeId(3)].y,
            m_A[face.nodeId(3)].z
          };

          RealMatrix<12, 12> ParaxialElementMatrix = _computeParaxialElementMatrixQuad4(face);
          RealVector<12> rhs = area * (c7 * Un * ParaxialElementMatrix - c8 * Vn * ParaxialElementMatrix + c9 * An * ParaxialElementMatrix);

          rhs_values[node_dof.dofId(face.nodeId(0), 0)] += rhs(0);
          rhs_values[node_dof.dofId(face.nodeId(0), 1)] += rhs(1);
          rhs_values[node_dof.dofId(face.nodeId(0), 2)] += rhs(2);
          rhs_values[node_dof.dofId(face.nodeId(1), 0)] += rhs(3);
          rhs_values[node_dof.dofId(face.nodeId(1), 1)] += rhs(4);
          rhs_values[node_dof.dofId(face.nodeId(1), 2)] += rhs(5);
          rhs_values[node_dof.dofId(face.nodeId(2), 0)] += rhs(6);
          rhs_values[node_dof.dofId(face.nodeId(2), 1)] += rhs(7);
          rhs_values[node_dof.dofId(face.nodeId(2), 2)] += rhs(8);
          rhs_values[node_dof.dofId(face.nodeId(3), 0)] += rhs(9);
          rhs_values[node_dof.dofId(face.nodeId(3), 1)] += rhs(10);
          rhs_values[node_dof.dofId(face.nodeId(3), 2)] += rhs(11);

          if (t <= dt) {
            ParaxialElementMatrix = area * c7 * ParaxialElementMatrix;
            Int32 n1_index = 0;
            for (Node node1 : face.nodes()) {
              Int32 n2_index = 0;
              for (Node node2 : face.nodes()) {
                Real v1 = ParaxialElementMatrix(3 * n1_index, 3 * n2_index);
                Real v2 = ParaxialElementMatrix(3 * n1_index, 3 * n2_index + 1);
                Real v3 = ParaxialElementMatrix(3 * n1_index, 3 * n2_index + 2);

                Real v4 = ParaxialElementMatrix(3 * n1_index + 1, 3 * n2_index);
                Real v5 = ParaxialElementMatrix(3 * n1_index + 1, 3 * n2_index + 1);
                Real v6 = ParaxialElementMatrix(3 * n1_index + 1, 3 * n2_index + 2);

                Real v7 = ParaxialElementMatrix(3 * n1_index + 2, 3 * n2_index);
                Real v8 = ParaxialElementMatrix(3 * n1_index + 2, 3 * n2_index + 1);
                Real v9 = ParaxialElementMatrix(3 * n1_index + 2, 3 * n2_index + 2);
                if (node1.isOwn()) {
                  DoFLocalId node1_dof1 = node_dof.dofId(node1, 0);
                  DoFLocalId node1_dof2 = node_dof.dofId(node1, 1);
                  DoFLocalId node1_dof3 = node_dof.dofId(node1, 2);

                  DoFLocalId node2_dof1 = node_dof.dofId(node2, 0);
                  DoFLocalId node2_dof2 = node_dof.dofId(node2, 1);
                  DoFLocalId node2_dof3 = node_dof.dofId(node2, 2);

                  if (bsr_matrix) {
                    bsr_matrix->addValue(node1_dof1, node2_dof1, v1);
                    bsr_matrix->addValue(node1_dof1, node2_dof2, v2);
                    bsr_matrix->addValue(node1_dof1, node2_dof3, v3);

                    bsr_matrix->addValue(node1_dof2, node2_dof1, v4);
                    bsr_matrix->addValue(node1_dof2, node2_dof2, v5);
                    bsr_matrix->addValue(node1_dof2, node2_dof3, v6);

                    bsr_matrix->addValue(node1_dof3, node2_dof1, v7);
                    bsr_matrix->addValue(node1_dof3, node2_dof2, v8);
                    bsr_matrix->addValue(node1_dof3, node2_dof3, v9);
                  }
                  else {
                    m_linear_system.matrixAddValue(node1_dof1, node2_dof1, v1);
                    m_linear_system.matrixAddValue(node1_dof1, node2_dof2, v2);
                    m_linear_system.matrixAddValue(node1_dof1, node2_dof3, v3);

                    m_linear_system.matrixAddValue(node1_dof2, node2_dof1, v4);
                    m_linear_system.matrixAddValue(node1_dof2, node2_dof2, v5);
                    m_linear_system.matrixAddValue(node1_dof2, node2_dof3, v6);

                    m_linear_system.matrixAddValue(node1_dof3, node2_dof1, v7);
                    m_linear_system.matrixAddValue(node1_dof3, node2_dof2, v8);
                    m_linear_system.matrixAddValue(node1_dof3, node2_dof3, v9);
                  }
                }
                ++n2_index;
              }
              ++n1_index;
            }
          }
        }
      else
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, m_node_coord);

          RealVector<9> Un = {
            m_U[face.nodeId(0)].x,
            m_U[face.nodeId(0)].y,
            m_U[face.nodeId(0)].z,
            m_U[face.nodeId(1)].x,
            m_U[face.nodeId(1)].y,
            m_U[face.nodeId(1)].z,
            m_U[face.nodeId(2)].x,
            m_U[face.nodeId(2)].y,
            m_U[face.nodeId(2)].z,
          };

          RealVector<9> Vn = {
            m_V[face.nodeId(0)].x,
            m_V[face.nodeId(0)].y,
            m_V[face.nodeId(0)].z,
            m_V[face.nodeId(1)].x,
            m_V[face.nodeId(1)].y,
            m_V[face.nodeId(1)].z,
            m_V[face.nodeId(2)].x,
            m_V[face.nodeId(2)].y,
            m_V[face.nodeId(2)].z,
          };

          RealVector<9> An = {
            m_A[face.nodeId(0)].x,
            m_A[face.nodeId(0)].y,
            m_A[face.nodeId(0)].z,
            m_A[face.nodeId(1)].x,
            m_A[face.nodeId(1)].y,
            m_A[face.nodeId(1)].z,
            m_A[face.nodeId(2)].x,
            m_A[face.nodeId(2)].y,
            m_A[face.nodeId(2)].z,
          };

          RealMatrix<9, 9> ParaxialElementMatrixTria3 = _computeParaxialElementMatrixTria3(face);
          RealVector<9> rhs = area * (c7 * Un * ParaxialElementMatrixTria3 - c8 * Vn * ParaxialElementMatrixTria3 + c9 * An * ParaxialElementMatrixTria3);

          rhs_values[node_dof.dofId(face.nodeId(0), 0)] += rhs(0);
          rhs_values[node_dof.dofId(face.nodeId(0), 1)] += rhs(1);
          rhs_values[node_dof.dofId(face.nodeId(0), 2)] += rhs(2);
          rhs_values[node_dof.dofId(face.nodeId(1), 0)] += rhs(3);
          rhs_values[node_dof.dofId(face.nodeId(1), 1)] += rhs(4);
          rhs_values[node_dof.dofId(face.nodeId(1), 2)] += rhs(5);
          rhs_values[node_dof.dofId(face.nodeId(2), 0)] += rhs(6);
          rhs_values[node_dof.dofId(face.nodeId(2), 1)] += rhs(7);
          rhs_values[node_dof.dofId(face.nodeId(2), 2)] += rhs(8);

          if (t <= dt) {
            ParaxialElementMatrixTria3 = area * c7 * ParaxialElementMatrixTria3;
            Int32 n1_index = 0;
            for (Node node1 : face.nodes()) {
              Int32 n2_index = 0;
              for (Node node2 : face.nodes()) {
                Real v1 = ParaxialElementMatrixTria3(3 * n1_index, 3 * n2_index);
                Real v2 = ParaxialElementMatrixTria3(3 * n1_index, 3 * n2_index + 1);
                Real v3 = ParaxialElementMatrixTria3(3 * n1_index, 3 * n2_index + 2);

                Real v4 = ParaxialElementMatrixTria3(3 * n1_index + 1, 3 * n2_index);
                Real v5 = ParaxialElementMatrixTria3(3 * n1_index + 1, 3 * n2_index + 1);
                Real v6 = ParaxialElementMatrixTria3(3 * n1_index + 1, 3 * n2_index + 2);

                Real v7 = ParaxialElementMatrixTria3(3 * n1_index + 2, 3 * n2_index);
                Real v8 = ParaxialElementMatrixTria3(3 * n1_index + 2, 3 * n2_index + 1);
                Real v9 = ParaxialElementMatrixTria3(3 * n1_index + 2, 3 * n2_index + 2);
                if (node1.isOwn()) {
                  DoFLocalId node1_dof1 = node_dof.dofId(node1, 0);
                  DoFLocalId node1_dof2 = node_dof.dofId(node1, 1);
                  DoFLocalId node1_dof3 = node_dof.dofId(node1, 2);

                  DoFLocalId node2_dof1 = node_dof.dofId(node2, 0);
                  DoFLocalId node2_dof2 = node_dof.dofId(node2, 1);
                  DoFLocalId node2_dof3 = node_dof.dofId(node2, 2);

                  if (bsr_matrix) {
                    bsr_matrix->addValue(node1_dof1, node2_dof1, v1);
                    bsr_matrix->addValue(node1_dof1, node2_dof2, v2);
                    bsr_matrix->addValue(node1_dof1, node2_dof3, v3);

                    bsr_matrix->addValue(node1_dof2, node2_dof1, v4);
                    bsr_matrix->addValue(node1_dof2, node2_dof2, v5);
                    bsr_matrix->addValue(node1_dof2, node2_dof3, v6);

                    bsr_matrix->addValue(node1_dof3, node2_dof1, v7);
                    bsr_matrix->addValue(node1_dof3, node2_dof2, v8);
                    bsr_matrix->addValue(node1_dof3, node2_dof3, v9);
                  }
                  else {
                    m_linear_system.matrixAddValue(node1_dof1, node2_dof1, v1);
                    m_linear_system.matrixAddValue(node1_dof1, node2_dof2, v2);
                    m_linear_system.matrixAddValue(node1_dof1, node2_dof3, v3);

                    m_linear_system.matrixAddValue(node1_dof2, node2_dof1, v4);
                    m_linear_system.matrixAddValue(node1_dof2, node2_dof2, v5);
                    m_linear_system.matrixAddValue(node1_dof2, node2_dof3, v6);

                    m_linear_system.matrixAddValue(node1_dof3, node2_dof1, v7);
                    m_linear_system.matrixAddValue(node1_dof3, node2_dof2, v8);
                    m_linear_system.matrixAddValue(node1_dof3, node2_dof3, v9);
                  }
                }
                ++n2_index;
              }
              ++n1_index;
            }
          }
        }
    }
  }
}