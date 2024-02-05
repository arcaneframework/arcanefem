// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BlcsrBiliAssembly.hxx                                     (C) 2022-2023   */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/* which avoid to add in the global matrix by iterating through the node.    */
/* It supports GPU Parallelization                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

void FemModule::_buildMatrixNodeWiseCsr()
{

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // Compute the number of nnz and initialize the memory space
  Int32 nnz = nbFace() * 2 + nbNode();
  m_csr_matrix.initialize(m_dof_family, nnz, nbNode());

  /*removing the neoighbouring currently as it is useless
  // Creating a connectivity from node to their neighbouring nodes
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  cn = idx_cn->connectivity();
  */
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_assembleNodeWiseCsrBilinearOperatorTria3()
{
  Timer::Action timer_blcsr_bili(m_time_stats, "AssembleNodeWiseCsrBilinearOperatorTria3");
  {
    Timer::Action timer_blcsr_build(m_time_stats, "NodeWiseCsrBuildMatrix");
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

  Timer::Action timer_blcsr_add_compute(m_time_stats, "NodeWiseCsrAddAndCompute");
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
