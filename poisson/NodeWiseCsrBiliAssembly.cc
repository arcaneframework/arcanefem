// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* NWCSRiliAssembly.hxx                                      (C) 2022-2024   */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/* which avoid to add in the global matrix by iterating through the node.    */
/* It supports GPU Parallelization                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Builds the node-wise CSR matrix for the finite element method.
 *
 * This function initializes and fills the CSR (Compressed Sparse Row) matrix
 * based on the node-wise connectivity of the mesh. It supports different mesh
 * types, specifically "TRIA3" and "TETRA4".
 *
 * For "TRIA3" mesh type, the function iterates over all nodes and sets the
 * coordinates in the CSR matrix based on the node's faces.
 *
 * For "TETRA4" mesh type, the function iterates over all nodes and sets the
 * coordinates in the CSR matrix based on the node's edges.
 *
 * The function also initializes the memory space required for the CSR matrix
 * based on the number of nodes and edges/faces in the mesh.
 */
/*---------------------------------------------------------------------------*/

void FemModule::_buildMatrixNodeWiseCsr()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // Compute the number of nnz and initialize the memory space
  Integer nbnde = nbNode();
  Int64 nedge = options()->meshType == "TETRA4" ? nbEdge() : nbFace();
  Int32 nnz = nedge * 2 + nbnde;
  m_csr_matrix.initialize(m_dof_family, nnz, nbnde);

  /*removing the neoighbouring currently as it is useless
  // Creating a connectivity from node to their neighbouring nodes
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  cn = idx_cn->connectivity();
  */

  if (options()->meshType == "TRIA3") {
    ENUMERATE_NODE (inode, allNodes()) {

      //Since we compute the neighbouring connectivity here, we also fill the csr matrix

      Node node = *inode;

      m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

      for (Face face : node.faces()) {
        if (face.nodeId(0) == node.localId()) {
          //    cn->addConnectedItem(node, face.node(0));
          m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
        }
        else {
          //  cn->addConnectedItem(node, face.node(1));
          m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
        }
      }
    }
  }
  else if (options()->meshType == "TETRA4") {
    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;
      Int32 node_dof_id = node_dof.dofId(node, 0);
      ItemLocalIdT<DoF> diagonal_entry(node_dof_id);

      m_csr_matrix.setCoordinates(diagonal_entry, diagonal_entry);

      for (Edge edge : node.edges()) {
        if (edge.nodeId(0) == node.localId())
          m_csr_matrix.setCoordinates(diagonal_entry, node_dof.dofId(edge.nodeId(1), 0));
        else
          m_csr_matrix.setCoordinates(diagonal_entry, node_dof.dofId(edge.nodeId(0), 0));
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_assembleNodeWiseCsrBilinearOperatorTria3()
{
  Timer::Action timer_nwcsr_bili(m_time_stats, "AssembleNodeWiseCsrBilinearOperatorTria3");
  {
    Timer::Action timer_nwcsr_build(m_time_stats, "NodeWiseCsrBuildMatrix");
    // Build the csr matrix
    _buildMatrixNodeWiseCsr();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();

  // Boucle sur les noeuds déportée sur accélérateur
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
  Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
  auto in_out_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);

  auto in_node_coord = ax::viewIn(command, m_node_coord);

  UnstructuredMeshConnectivityView m_connectivity_view;
  m_connectivity_view.setMesh(this->mesh());
  auto ncc = m_connectivity_view.nodeCell();
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
  Arcane::ItemGenericInfoListView cells_infos(this->mesh()->cellFamily());

  Timer::Action timer_nwcsr_add_compute(m_time_stats, "NodeWiseCsrAddAndCompute");
  command << RUNCOMMAND_ENUMERATE(Node, inode, allNodes())
  {
    Int32 inode_index = 0;
    for (auto cell : ncc.cells(inode)) {

      // How can I know the right index ?
      // By checking in the global id ?
      // Working currently, but maybe only because p = 1 ?
      if (inode == cnc.nodeId(cell, 1)) {
        inode_index = 1;
      }
      else if (inode == cnc.nodeId(cell, 2)) {
        inode_index = 2;
      }
      else {
        inode_index = 0;
      }
      Real3 m0 = in_node_coord[cnc.nodeId(cell, 0)];
      Real3 m1 = in_node_coord[cnc.nodeId(cell, 1)];
      Real3 m2 = in_node_coord[cnc.nodeId(cell, 2)];

      Real area = 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y)); // calculate area

      Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
      Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
      Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

      Real b_matrix[3][2] = { 0 };
      Real mul = (1.0 / (2.0 * area));
      b_matrix[0][0] = dPhi0.x * mul;
      b_matrix[0][1] = dPhi0.y * mul;

      b_matrix[1][0] = dPhi1.x * mul;
      b_matrix[1][1] = dPhi1.y * mul;

      b_matrix[2][0] = dPhi2.x * mul;
      b_matrix[2][1] = dPhi2.y * mul;

      Int32 i = 0;
      for (NodeLocalId node2 : cnc.nodes(cell)) {
        Real x = 0.0;
        for (Int32 k = 0; k < 2; k++) {
          x += b_matrix[inode_index][k] * b_matrix[i][k];
        }
        if (nodes_infos.isOwn(inode)) {

          Int32 row = node_dof.dofId(inode, 0).localId();
          Int32 col = node_dof.dofId(node2, 0).localId();
          Int32 begin = in_row_csr[row];
          Int32 end;
          if (row == row_csr_size - 1) {
            end = col_csr_size;
          }
          else {
            end = in_row_csr[row + 1];
          }
          while (begin < end) {
            if (in_col_csr[begin] == col) {
              in_out_val_csr[begin] += x * area;
              break;
            }
            begin++;
          }
        }
        i++;
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the node-wise CSR bilinear operator for TETRA4 mesh type.
 *
 * This function builds the CSR matrix and computes the bilinear operator
 * for a TETRA4 mesh type. It initializes the CSR matrix, sets up the 
 * necessary views, and performs the computation on the accelerator.
 */
void FemModule::_assembleNodeWiseCsrBilinearOperatorTetra4()
{
  Timer::Action timer_nwcsr_bili(m_time_stats, "AssembleNodeWiseCsrBilinearOperatorTetra4");
  {
    Timer::Action timer_nwcsr_build(m_time_stats, "NodeWiseCsrBuildMatrix");
    _buildMatrixNodeWiseCsr();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
  Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
  auto in_out_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);
  auto in_node_coord = ax::viewIn(command, m_node_coord);

  UnstructuredMeshConnectivityView m_connectivity_view;
  m_connectivity_view.setMesh(this->mesh());
  auto ncc = m_connectivity_view.nodeCell();
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
  Arcane::ItemGenericInfoListView cells_infos(this->mesh()->cellFamily());

  Timer::Action timer_nwcsr_add_compute(m_time_stats, "NodeWiseCsrAddAndCompute");
  command << RUNCOMMAND_ENUMERATE(Node, inode, allNodes())
  {
    for (auto cell : ncc.cells(inode)) {
      Int32 inode_index = 0;
      for (Int32 i = 0; i < 4; ++i) {
        if (inode == cnc.nodeId(cell, i)) {
          inode_index = i;
          break;
        }
      }

      Real3 m0 = in_node_coord[cnc.nodeId(cell, 0)];
      Real3 m1 = in_node_coord[cnc.nodeId(cell, 1)];
      Real3 m2 = in_node_coord[cnc.nodeId(cell, 2)];
      Real3 m3 = in_node_coord[cnc.nodeId(cell, 3)];

      Real3 v0 = m1 - m0;
      Real3 v1 = m2 - m0;
      Real3 v2 = m3 - m0;
      Real volume = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2))) / 6.0;

      Real3 dPhi[4] = {
        Arcane::math::cross(m2 - m1, m1 - m3),
        Arcane::math::cross(m3 - m0, m0 - m2),
        Arcane::math::cross(m1 - m0, m0 - m3),
        Arcane::math::cross(m0 - m1, m1 - m2)
      };

      Real b_matrix[4][3] = { 0 };
      Real mul = (1.0 / (6.0 * volume));
      for (Int32 i = 0; i < 4; ++i) {
        b_matrix[i][0] = dPhi[i].x * mul;
        b_matrix[i][1] = dPhi[i].y * mul;
        b_matrix[i][2] = dPhi[i].z * mul;
      }

      for (Int32 i = 0; i < 4; ++i) {
        Real x = 0.0;
        for (Int32 k = 0; k < 3; ++k) {
          x += b_matrix[inode_index][k] * b_matrix[i][k];
        }
        if (nodes_infos.isOwn(inode)) {
          Int32 row = node_dof.dofId(inode, 0).localId();
          Int32 col = node_dof.dofId(cnc.nodeId(cell, i), 0).localId();
          Int32 begin = in_row_csr[row];
          Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];
          while (begin < end) {
            if (in_col_csr[begin] == col) {
              in_out_val_csr[begin] += x * volume;
              break;
            }
            ++begin;
          }
        }
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/