﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrBiliAssembly.hxx                                         (C) 2022-2025 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/* which handle the parallelization on CPU's                                 */
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
  //ARCANE_FATAL("BUILD_CSR");
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
  IMesh* mesh = defaultMesh();
  Int64 nedge = 0;
  Int64 nbnde = nbNode();
  Int32 mesh_dim = mesh->dimension();
  if (mesh_dim == 3) {
    // m_nb_edge should have been computed in startInit()
    nedge = m_nb_edge;
  }
  else if (mesh_dim == 2)
    nedge = nbFace();
  else
    ARCANE_THROW(NotSupportedException, "Only mesh of dimension 2 or 3 are supported");

  Int32 nnz = nedge * 2 + nbnde;

  m_csr_matrix.initialize(m_dof_family, nnz, nbNode(), m_queue);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  //We iterate through the node, and we do not sort anymore : we assume the nodes ID are sorted, and we will iterate throught the column to avoid making < and > comparison
  if (mesh_dim == 2) {
    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;

      //info() << "DEBUG Add:   (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(node, 0) << " )";
      m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

      for (Face face : node.faces()) {
        if (face.nodeId(0) == node.localId()) {
          //info() << "DEBUG Add: 0 (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(1), 0) << " )";
          m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
        }
        else {
          //info() << "DEBUG Add:   (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(face.nodeId(0), 0) << " )";
          m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
        }
      }
    }
  }
  else if (mesh_dim == 3) {
    // Connectivity should have been build during init (in startInit() entry point)
    auto* connectivity_ptr = m_node_node_via_edge_connectivity.get();
    ARCANE_CHECK_POINTER(connectivity_ptr);
    IndexedNodeNodeConnectivityView nn_cv = connectivity_ptr->view();

    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;
      DoFLocalId dof = node_dof.dofId(node, 0);
      //info() << "DEBUG Add:   (" << node_dof.dofId(node, 0) << ", " << node_dof.dofId(node, 0) << " )";
      m_csr_matrix.setCoordinates(dof, node_dof.dofId(node, 0));
      for (NodeLocalId other_node : nn_cv.nodeIds(node))
        m_csr_matrix.setCoordinates(dof, node_dof.dofId(other_node, 0));
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrBilinearOperatorTRIA3()
{

  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_Csr");
  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _buildMatrixCsr();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    RealMatrix<3, 3> K_e;
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
  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_Csr");
  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _buildMatrixCsr();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    RealMatrix<4, 4> K_e;
    {
      K_e = _computeElementMatrixTETRA4(cell); // element stiffness matrix
    }

    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        Real v = K_e(n1_index, n2_index);
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
