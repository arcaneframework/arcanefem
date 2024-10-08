// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CooBiliAssembly.hxx                                         (C) 2022-2024 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the coo data structure       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Builds the coordinate (COO) matrix for the finite element method (FEM) module.
 *
 * This function initializes and populates the COO matrix based on the mesh's nodes and their connectivity.
 * It assumes a polynomial degree (p) of 1. The matrix is constructed by iterating over all nodes in the mesh
 * and setting the coordinates for the degrees of freedom (DOFs) associated with each node and its connected
 * edges or faces, depending on the mesh dimension.
 *
 * @note This implementation is specific to 2D and 3D meshes.
 */
/*---------------------------------------------------------------------------*/

void FemModule::_buildMatrixCoo()
{
  Int8 mesh_dim = mesh()->dimension();
  Int64 nbEdge = mesh_dim == 3 ? m_nb_edge : nbFace();
  Int32 nnz = nbEdge * 2 + nbNode();
  m_coo_matrix.initialize(m_dof_family, nnz);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_NODE (inode, allNodes()) {
    Node node = *inode;

    m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

    if (mesh_dim == 2) {
      for (Face face : node.faces()) {
        Node other_node = (face.nodeId(0) == node.localId()) ? face.node(1) : face.node(0);
        m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(other_node, 0));
      }
    }
    else {
      for (Edge edge : node.edges()) {
        Node other_node = (edge.nodeId(0) == node.localId()) ? edge.node(1) : edge.node(0);
        m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(other_node, 0));
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system with
 * the COO sparse matrix format for TRIA3 elements.
 *
 * The method performs the following steps:
 *   1. Builds the COO matrix using _buildMatrixCoo.
 *   2. For each cell, computes element matrix using _computeElementMatrixTRIA3.
 *   3. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

void FemModule::_assembleCooBilinearOperatorTRIA3()
{
  info() << "Assembling COO Bilinear Operator TRIA3";

  Timer::Action timer_coo_bili(m_time_stats, "AssembleCooBilinearOperatorTria3");

  {
    Timer::Action timer_coo_build(m_time_stats, "BuildMatrixCoo");
    _buildMatrixCoo();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    FixedMatrix<3, 3> K_e;
    {
      Timer::Action timer_coo_compute(m_time_stats, "CooComputeElementMatrixTria3");
      K_e = _computeElementMatrixTRIA3(cell);
    }

    Timer::Action timer_coo_add(m_time_stats, "CooAddToGlobalMatrix");
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
 * the COO sparse matrix format.
 *
 * The method performs the following steps:
 *   1. Builds the COO matrix using _buildMatrixCoo.
 *   2. For each cell, computes element matrix using _computeElementMatrixTETRA4.
 *   3. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

void FemModule::_assembleCooBilinearOperatorTETRA4()
{
  info() << "Assembling COO Bilinear Operator TETRA4";

  Timer::Action timer_coo_bili(m_time_stats, "AssembleCooBilinearOperatorTetra4");

  {
    Timer::Action timer_coo_build(m_time_stats, "BuildMatrixCoo");
    _buildMatrixCoo();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    FixedMatrix<4, 4> K_e;
    {
      Timer::Action timer_coo_compute(m_time_stats, "CooComputeElementMatrixTetra4");
      K_e = _computeElementMatrixTETRA4(cell);
    }

    Timer::Action timer_coo_add(m_time_stats, "CooAddToGlobalMatrix");
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