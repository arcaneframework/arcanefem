// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrBiliAssembly.hxx                                         (C) 2022-2023 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/**
 * @brief Initialization of the csr matrix. It only works for p=1 since there is
 * one node per Edge.
 * 
 * 
 */
void FemModule::
_buildMatrixCsr()
{
  //Initialization of the coo matrix;
  //This formula only works in p=1

  /*
  //Create a connection between nodes through the faces
  //Useless here because we only need this information once
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  auto idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  auto* cn = idx_cn->connectivity();
  ENUMERATE_NODE (node, allNodes()) {
  }
  */

  Int32 nnz;
  if (options()->meshType == "TETRA4")
    nnz = nbFace() * 2 + nbNode() + 2;   // For testing only 
  else if (options()->meshType == "TRIA3")
    nnz = nbFace() * 2 + nbNode();
  else
    ARCANE_THROW(NotImplementedException, "");

  m_csr_matrix.initialize(m_dof_family, nnz, nbNode());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  //We iterate through the node, and we do not sort anymore : we assume the nodes ID are sorted, and we will iterate throught the column to avoid making < and > comparison
  if (options()->meshType == "TRIA3"){
    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;

      //info() << "DEBUG Add:   (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(node, 0) << " )";
      m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

      for (Face face : node.faces()) {
        if (face.nodeId(0) == node.localId()){
          //info() << "DEBUG Add: 0 (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(1), 0) << " )";
          m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
        }else{
          //info() << "DEBUG Add:   (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(0), 0) << " )";
          m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
        }
      }
    }
  }
  else if (options()->meshType == "TETRA"){
    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;

      m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

      for (Face face : node.faces()) {
        if (face.nodeId(0) == node.localId()){
          info() << "DEBUG Add: 0 (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(1), 0) << " )";
          info() << "DEBUG Add:   (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(2), 0) << " )";
        //m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
        //m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(2), 0));
        }
        else if (face.nodeId(1) == node.localId()){
          info() << "DEBUG Add: 1 (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(0), 0) << " )";
          info() << "DEBUG Add:   (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(2), 0) << " )";
        //m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
        //m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(2), 0));
        }
        else if (face.nodeId(2) == node.localId()){
          info() << "DEBUG Add: 2 (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(0), 0) << " )";
          info() << "DEBUG Add:   (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(1), 0) << " )";
        //m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
        //m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
        }
      }
    }
  }
  else{
    ARCANE_THROW(NotImplementedException, "");
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrBilinearOperatorTRIA3()
{

  Timer::Action timer_csr_bili(m_time_stats, "AssembleCsrBilinearOperatorTria3");
  {
    Timer::Action timer_csr_build(m_time_stats, "CsrBuildMatrix");
    // Build the csr matrix
    _buildMatrixCsr();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    FixedMatrix<3, 3> K_e;
    {
      //Timer::Action timer_csr_compute_add(m_time_stats, "CsrComputeElementMatrixTria3");
      K_e = _computeElementMatrixTRIA3(cell); // element stifness matrix
    }
    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]

    //Timer::Action timer_action(m_time_stats, "CsrAddToGlobalMatrix");
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e(n1_index, n2_index);
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (node1.isOwn()) {
          m_csr_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrBilinearOperatorTETRA4()
{

  Timer::Action timer_csr_bili(m_time_stats, "AssembleCsrBilinearOperatorTetra4");
  {
    _buildMatrixCsr();
  }
/*
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    FixedMatrix<4, 4> K_e;
    {
      K_e = _computeElementMatrixTETRA4(cell); // element stifness matrix
    }

    //Timer::Action timer_action(m_time_stats, "CsrAddToGlobalMatrix");
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e(n1_index, n2_index);
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (node1.isOwn()) {
          m_csr_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }
  }
*/
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
