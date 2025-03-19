// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CooSortBiliAssembly.hxx                                     (C) 2022-2025 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the coo data structure       */
/* followed by a sort                                                        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Builds the COO Sort matrix by first filling the diagonal and then the off-diagonal entries.
 *
 * This method initializes the diagonal elements of the matrix and then proceeds to fill the remaining
 * entries based on the connectivity of the mesh entities (nodes, edges, faces). After filling the matrix,
 * a sorting step is performed to optimize further operations.
 * 
 * @note This implementation is specific to 2D and 3D meshes.
 *
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_buildMatrixCooSort()
{
  Int8 mesh_dim = mesh()->dimension();
  Int64 nbEdge = mesh_dim == 3 ? m_nb_edge : nbFace();
  Int32 nnz = nbEdge * 2 + nbNode();
  m_coo_matrix.initialize(m_dof_family, nnz, m_queue);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (mesh_dim == 2) {
    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;
      DoFLocalId dof = node_dof.dofId(node, 0);
      m_coo_matrix.setCoordinates(dof, dof);

      for (Face face : node.faces()) {
        Node other_node = (face.nodeId(0) == node.localId()) ? face.node(1) : face.node(0);
        m_coo_matrix.setCoordinates(dof, node_dof.dofId(other_node, 0));
      }
    }
  }
  else {
    auto* connectivity_ptr = m_node_node_via_edge_connectivity.get();
    ARCANE_CHECK_POINTER(connectivity_ptr);
    IndexedNodeNodeConnectivityView nn_cv = connectivity_ptr->view();

    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;
      DoFLocalId dof = node_dof.dofId(node, 0);
      m_coo_matrix.setCoordinates(dof, dof);

      for (NodeLocalId other_node : nn_cv.nodeIds(node))
        m_coo_matrix.setCoordinates(dof, node_dof.dofId(other_node, 0));
    }
  }

  // Sort both row and column arrays of the matrix
  m_coo_matrix.sort();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system with
 * the COO Sorted sparse matrix format for TRIA3 elements.
 *
 * The method performs the following steps:
 *   1. Builds the COO Sorted matrix using _buildMatrixCooSort.
 *   2. For each cell, computes element matrix using _computeElementMatrixTRIA3.
 *   3. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 * 
 * @note This implementation is the same as the one in CooBiliAssembly.cc, but
 * uses the sorted COO matrix.
 * 
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCooSortBilinearOperatorTRIA3()
{
  info() << "Assembling COO Sort Bilinear Operator TRIA3";

  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_CooSort");

  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _buildMatrixCooSort();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    RealMatrix<3, 3> K_e;
    {
      K_e = _computeElementMatrixTRIA3(cell);
    }

    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        Real v = K_e(n1_index, n2_index);

        if (node1.isOwn())
          m_coo_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);

        ++n2_index;
      }
      ++n1_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system with
 * the COO Sorted sparse matrix format for TETRA4 elements.
 *
 * The method performs the following steps:
 *   1. Builds the COO Sorted matrix using _buildMatrixCooSort.
 *   2. For each cell, computes element matrix using _computeElementMatrixTETRA4.
 *   3. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 * 
 * @note This implementation is the same as the one in CooBiliAssembly.cc, but
 * uses the sorted COO matrix.
 * 
 */
/*---------------------------------------------------------------------------*/

void FemModule::_assembleCooSortBilinearOperatorTETRA4()
{
  info() << "Assembling COO Sort Bilinear Operator TETRA4";

  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_CooSort");

  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _buildMatrixCooSort();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    RealMatrix<4, 4> K_e;
    {
      K_e = _computeElementMatrixTETRA4(cell);
    }

    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        Real v = K_e(n1_index, n2_index);

        if (node1.isOwn())
          m_coo_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);

        ++n2_index;
      }
      ++n1_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
