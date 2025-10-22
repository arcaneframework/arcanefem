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
 * @brief Applies paraxial boundary conditions to the RHS and LHS
 * of the linear system.
 * 
 * @param rhs_values The RHS values to be updated.
 * @param node_dof The node DoF connectivity view.
 * @param bsr_matrix The BSR matrix to be updated (optional).
 * 
/*---------------------------------------------------------------------------*/

void FemModule::
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