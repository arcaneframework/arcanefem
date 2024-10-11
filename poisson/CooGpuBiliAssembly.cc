
// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CooGpuiBiliAssembly.hxx                                     (C) 2022-2023 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the coo data structure       */
/* which handle the parallelization on GPU using Arcane accelerator API and  */
/* an atomic operation for adding the value into the global matrix           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * 
 * !!!! Whereas this function is called GPU, it is not actually executed on the GPU.
 * It is only used in the GPU version of the bilinear assembly to build the matrix,
 * this function is in fact the same as in CooBiliAssembly, it is used here only to
 * reiterate the old terminology and structure of the code !!!!
 * 
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

void FemModule::_buildMatrixCooGPU()
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
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCooGPUBilinearOperatorTRIA3()
{
  info() << "Assembling COO GPU Bilinear Operator for TRIA3 elements";

  Timer::Action timer_coo_gpu_bili(m_time_stats, "AssembleCooGpuBilinearOperatorTria3");

  {
    Timer::Action timer_coo_gpu_build(m_time_stats, "BuildMatrixCooGpu");
    _buildMatrixCooGPU();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  Int32 row_length = m_coo_matrix.m_matrix_row.totalNbElement();

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto in_row_coo = ax::viewIn(command, m_coo_matrix.m_matrix_row);
  auto in_col_coo = ax::viewIn(command, m_coo_matrix.m_matrix_column);
  auto in_out_val_coo = ax::viewInOut(command, m_coo_matrix.m_matrix_value);
  UnstructuredMeshConnectivityView m_connectivity_view;
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_connectivity_view.setMesh(this->mesh());
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
  Arcane::ItemGenericInfoListView cells_infos(this->mesh()->cellFamily());

  command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
  {
    Real K_e[9] = { 0 };
    {
        _computeElementMatrixTRIA3GPU(icell, cnc, in_node_coord, K_e);
    }

    Int32 n1_index = 0;
    for (NodeLocalId node1 : cnc.nodes(icell)) {
      Int32 n2_index = 0;
      for (NodeLocalId node2 : cnc.nodes(icell)) {
        Real v = K_e[n1_index * 3 + n2_index];

        if (nodes_infos.isOwn(node1)) {
          Int32 row_index = node_dof.dofId(node1, 0);
          Int32 col_index = node_dof.dofId(node2, 0);
          Int32 value_index;

          // Find the index of the value in the coo matrix
          for (value_index = 0; value_index < row_length; value_index++) {
            if (in_row_coo(value_index) == row_index && in_col_coo(value_index) == col_index) {
              ax::doAtomic<ax::eAtomicOperation::Add>(in_out_val_coo(value_index), v);
              break;
            }
          }
        }
        ++n2_index;
      }
      ++n1_index;
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCooGPUBilinearOperatorTETRA4() {
  info() << "Assembling COO GPU Bilinear Operator for TETRA4 elements";

  Timer::Action timer_coo_gpu_bili(m_time_stats, "AssembleCooGpuBilinearOperatorTetra4");

  {
    Timer::Action timer_coo_gpu_build(m_time_stats, "BuildMatrixCooGpu");
    _buildMatrixCooGPU();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  Int32 row_length = m_coo_matrix.m_matrix_row.totalNbElement();

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto in_row_coo = ax::viewIn(command, m_coo_matrix.m_matrix_row);
  auto in_col_coo = ax::viewIn(command, m_coo_matrix.m_matrix_column);
  auto in_out_val_coo = ax::viewInOut(command, m_coo_matrix.m_matrix_value);
  UnstructuredMeshConnectivityView m_connectivity_view;
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_connectivity_view.setMesh(this->mesh());
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
  Arcane::ItemGenericInfoListView cells_infos(this->mesh()->cellFamily());

  command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
  {
    Real K_e[16] = { 0 };
    {
        _computeElementMatrixTETRA4GPU(icell, cnc, in_node_coord, K_e);
    }

    Int32 n1_index = 0;
    for (NodeLocalId node1 : cnc.nodes(icell)) {
      Int32 n2_index = 0;
      for (NodeLocalId node2 : cnc.nodes(icell)) {
        if (nodes_infos.isOwn(node1)) {
          Real v = K_e[n1_index * 4 + n2_index];
          Int32 row_index = node_dof.dofId(node1, 0);
          Int32 col_index = node_dof.dofId(node2, 0);
          Int32 value_index;

          // Find the index of the value in the coo matrix
          for (value_index = 0; value_index < row_length; value_index++) {
            if (in_row_coo(value_index) == row_index && in_col_coo(value_index) == col_index) {
              ax::doAtomic<ax::eAtomicOperation::Add>(in_out_val_coo(value_index), v);
              break;
            }
          }
        }
        ++n2_index;
      }
      ++n1_index;
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/