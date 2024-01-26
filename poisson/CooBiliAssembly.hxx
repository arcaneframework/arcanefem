// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CooBiliAssembly.hxx                                         (C) 2022-2023 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the coo data structure       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Initialization of the coo matrix. It only works for p=1 since there is
 * one node per Edge.
 * 
 * 
 */
void FemModule::
_buildMatrix()
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

  Int32 nnz = nbFace() * 2 + nbNode();
  m_coo_matrix.initialize(m_dof_family, nnz);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  //We iterate through the node, and we do not sort anymore : we assume the nodes ID are sorted, and we will iterate throught the column to avoid making < and > comparison
  ENUMERATE_NODE (inode, allNodes()) {
    Node node = *inode;

    m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

    for (Face face : node.faces()) {
      if (face.nodeId(0) == node.localId())
        m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
      else
        m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCooBilinearOperatorTRIA3()
{

  Timer::Action timer_coo_bili(m_time_stats, "AssembleCooBilinearOperatorTria3");

  {
    Timer::Action timer_coo_build(m_time_stats, "CooBuildMatrix");
    // Build the coo matrix
    _buildMatrix();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    FixedMatrix<3, 3> K_e;
    {
      //Timer::Action timer_coo_compute(m_time_stats, "CooComputeElementMatrix");
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

    //Timer::Action timer_coo_add(m_time_stats, "CooAddToGlobalMatrix");
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e(n1_index, n2_index);
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (node1.isOwn()) {
          m_coo_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
