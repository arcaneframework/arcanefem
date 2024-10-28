
// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CooGpuiBiliAssembly.hxx                                     (C) 2022-2024 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the coo data structure       */
/* which handle the parallelization on GPU using Arcane accelerator API and  */
/* an atomic operation for adding the value into the global matrix           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE Int32 findIndexBinarySearch(Int32 row, Int32 col,
                                                const auto& in_row_coo,
                                                const auto& in_col_coo,
                                                Int32 row_length);

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

ARCCORE_HOST_DEVICE void _fillMatrix(Int32 id, Int64 nb_edge,
                                     auto nodes,
                                     auto inout_m_matrix_row,
                                     auto inout_m_matrix_column)
{
  Int32 node1_idx = nodes[0];
  Int32 node2_idx = nodes[1];

  // Upper triangular part
  inout_m_matrix_row[id] = node1_idx;
  inout_m_matrix_column[id] = node2_idx;

  // Matrix is symmetrical in Poisson
  // Lower triangular part
  inout_m_matrix_row[nb_edge + id] = node2_idx;
  inout_m_matrix_column[nb_edge + id] = node1_idx;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_fillDiagonal(Int64 nb_edge, NodeGroup nodes)
{
  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto inout_m_matrix_row = viewInOut(command, m_coo_matrix.m_matrix_row);
  auto inout_m_matrix_column = viewInOut(command, m_coo_matrix.m_matrix_column);

  Int32 offset = 2 * nb_edge;
  command << RUNCOMMAND_ENUMERATE(Node, node_id, nodes)
  {
    inout_m_matrix_row(offset + node_id) = node_id;
    inout_m_matrix_column(offset + node_id) = node_id;
  };
}

void FemModule::_buildMatrixCooGPU()
{
  Int8 mesh_dim = mesh()->dimension();

  if (mesh_dim != 2 && mesh_dim != 3)
    ARCANE_THROW(NotSupportedException, "Only mesh of dimension 2 or 3 are supported");

  Int64 nb_edge = (mesh_dim == 2) ? nbFace() : m_nb_edge;
  Int32 nb_non_zero = nb_edge * 2 + nbNode();
  m_coo_matrix.initialize(m_dof_family, nb_non_zero);

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  UnstructuredMeshConnectivityView m_connectivity_view;
  m_connectivity_view.setMesh(mesh());

  auto inout_m_matrix_row = viewInOut(command, m_coo_matrix.m_matrix_row);
  auto inout_m_matrix_column = viewInOut(command, m_coo_matrix.m_matrix_column);

  if (mesh_dim == 2) {

    if (m_use_coo_sort_gpu) {
      info()
      << "_buildMatrixCooSortGPU for 2D mesh with face-node connectivity";
      NumArray<uint, MDDim1> neighbors(nbNode());
      SmallSpan<uint> in_data = neighbors.to1DSmallSpan();

      UnstructuredMeshConnectivityView connectivity_view(mesh());
      auto node_face_connectivity_view = connectivity_view.nodeFace();

      {
        auto command = makeCommand(queue);
        command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
        {
          in_data[node_id] = node_face_connectivity_view.nbFace(node_id) + 1;
        };
      }
      queue->barrier();

      NumArray<uint, MDDim1> copy_output_data(nbNode());
      SmallSpan<uint> out_data = copy_output_data.to1DSmallSpan();
      Accelerator::Scanner<uint> scanner;
      scanner.exclusiveSum(queue, in_data, out_data);

      auto face_node_connectivity_view = connectivity_view.faceNode();

      {
        auto command = makeCommand(queue);
        // Fill the neighbors relation (including node with itself) into the
        // matrix
        command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
        {
          Int32 offset = out_data[node_id];

          for (auto edge_id : node_face_connectivity_view.faceIds(node_id)) {
            auto nodes = face_node_connectivity_view.nodes(edge_id);
            inout_m_matrix_row[offset] = node_id;
            inout_m_matrix_column[offset] =
            nodes[0] == node_id ? nodes[1] : nodes[0];
            ++offset;
          }

          inout_m_matrix_row[offset] = node_id;
          inout_m_matrix_column[offset] = node_id;
        };
      }
    }
    else {
      info() << "_buildMatrixCooGPU for 2D mesh with edge-node connectivity";

      IndexedFaceNodeConnectivityView face_node_connectivity_view =
      m_connectivity_view.faceNode();

      command << RUNCOMMAND_ENUMERATE(Face, face_id, allFaces())
      {
        auto nodes = face_node_connectivity_view.nodes(face_id);
        _fillMatrix(face_id, nb_edge, nodes, inout_m_matrix_row,
                    inout_m_matrix_column);
      };

      _fillDiagonal(nb_edge, allNodes());
    }
  }
  else if (options()
           ->createEdges()) { // 3D mesh without node-node connectivity

    if (m_use_coo_sort_gpu) {
      info()
      << "_buildMatrixCooSortGPU for 3D mesh with edge-node connectivity";
      NumArray<uint, MDDim1> neighbors(nbNode());
      SmallSpan<uint> in_data = neighbors.to1DSmallSpan();

      UnstructuredMeshConnectivityView connectivity_view(mesh());
      auto node_edge_connectivity_view = connectivity_view.nodeEdge();

      {
        auto command = makeCommand(queue);
        command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
        {
          in_data[node_id] = node_edge_connectivity_view.nbEdge(node_id) + 1;
        };
      }
      queue->barrier();

      NumArray<uint, MDDim1> copy_output_data(nbNode());
      SmallSpan<uint> out_data = copy_output_data.to1DSmallSpan();
      Accelerator::Scanner<uint> scanner;
      scanner.exclusiveSum(queue, in_data, out_data);

      auto edge_node_connectivity_view = connectivity_view.edgeNode();

      {
        auto command = makeCommand(queue);
        // Fill the neighbors relation (including node with itself) into the
        // matrix
        command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
        {
          Int32 offset = out_data[node_id];

          for (auto edge_id : node_edge_connectivity_view.edgeIds(node_id)) {
            auto nodes = edge_node_connectivity_view.nodes(edge_id);
            inout_m_matrix_row[offset] = node_id;
            inout_m_matrix_column[offset] =
            nodes[0] == node_id ? nodes[1] : nodes[0];
            ++offset;
          }

          inout_m_matrix_row[offset] = node_id;
          inout_m_matrix_column[offset] = node_id;
        };
      }
    }
    else {
      info() << "_buildMatrixCooGPU for 3D mesh with edge-node connectivity";
      IndexedEdgeNodeConnectivityView edge_node_connectivity_view =
      m_connectivity_view.edgeNode();

      command << RUNCOMMAND_ENUMERATE(Edge, edge_id, allEdges())
      {
        auto nodes = edge_node_connectivity_view.nodes(edge_id);
        _fillMatrix(edge_id, nb_edge, nodes, inout_m_matrix_row,
                    inout_m_matrix_column);
      };

      _fillDiagonal(nb_edge, allNodes());
    }
  }
  else { // 3D mesh with node-node connectivity

    auto* connectivity_ptr = m_node_node_via_edge_connectivity.get();
    ARCANE_CHECK_POINTER(connectivity_ptr);
    IndexedNodeNodeConnectivityView nn_cv = connectivity_ptr->view();

    // We allow the use of the accelerated algorithm only if the user provides
    // an accelerator runtime. On Cpu, the not-accelerated version is faster.

    if (queue->isAcceleratorPolicy()) {
      info() << "Using accelerated version of _buildMatrixCooGPU for 3D mesh "
                "with node-node connectivity";

      // This array will contain the number of neighbors of each node
      // (type uint is enough for counting neighbors).
      NumArray<uint, MDDim1> nb_neighbor_arr;
      nb_neighbor_arr.resize(nbNode());
      auto inout_nb_neighbor_arr = viewInOut(command, nb_neighbor_arr);

      command << RUNCOMMAND_ENUMERATE(Node, node_idx, allNodes())
      {
        inout_nb_neighbor_arr[node_idx] = nn_cv.nbNode(node_idx) + 1;
        // We add 1 to count the node's relation with itself.
      };

      // Do the exclusive cumulative sum of the neighbors array
      SmallSpan<uint> input = nb_neighbor_arr.to1DSmallSpan();
      NumArray<uint, MDDim1> tmp_output;
      tmp_output.resize(nbNode());
      SmallSpan<uint> output = tmp_output.to1DSmallSpan();
      Accelerator::Scanner<uint> scanner;
      scanner.exclusiveSum(queue, input, output);

      // Fill the neighbors relation (including node with itself) into the matrix
      command << RUNCOMMAND_ENUMERATE(Node, node_idx, allNodes())
      {
        Int32 offset = output[node_idx];

        for (NodeLocalId other_node_idx : nn_cv.nodeIds(node_idx)) {
          inout_m_matrix_row[offset] = node_idx;
          inout_m_matrix_column[offset] = other_node_idx;
          ++offset;
        }

        inout_m_matrix_row[offset] = node_idx;
        inout_m_matrix_column[offset] = node_idx;
      };
    }
    else {
      info() << "Using unaccelerated version of _buildMatrixCooGPU for 3D mesh "
                "with node-node connectivity";

      auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

      ENUMERATE_NODE (inode, allNodes()) {
        Node node = *inode;
        DoFLocalId dof = node_dof.dofId(node, 0);

        m_coo_matrix.setCoordinates(dof, node_dof.dofId(node, 0));

        for (NodeLocalId other_node : nn_cv.nodeIds(node))
          m_coo_matrix.setCoordinates(dof, node_dof.dofId(other_node, 0));
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

  Timer::Action timer_coo_gpu_bili(
  m_time_stats, m_use_coo_gpu ? "AssembleCooGpuBilinearOperatorTria3" : "AssembleCooSortGpuBilinearOperatorTria3");

  {
    Timer::Action timer_coo_gpu_build(m_time_stats,
                                      m_use_coo_gpu ? "BuildMatrixCooGpu"
                                                    : "BuildMatrixCooSortGpu");
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

  bool is_row_array_sorted = m_use_coo_sort_gpu || (mesh()->dimension() == 3 && !options()->createEdges());

  Timer::Action timer_coo_gpu_compute_add(m_time_stats, "CooGpuComputeAndAdd");

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

          // Find the index of the value in the coo matrix
          if (is_row_array_sorted) {
            auto value_index = findIndexBinarySearch(
            row_index, col_index, in_row_coo, in_col_coo, row_length);

            ax::doAtomic<ax::eAtomicOperation::Add>(in_out_val_coo(value_index),
                                                    v);
          }
          else {
            for (auto value_index = 0; value_index < row_length;
                 value_index++) {
              if (in_row_coo(value_index) == row_index &&
                  in_col_coo(value_index) == col_index) {
                ax::doAtomic<ax::eAtomicOperation::Add>(
                in_out_val_coo(value_index), v);
                break;
              }
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
_assembleCooGPUBilinearOperatorTETRA4()
{
  info() << "Assembling COO GPU Bilinear Operator for TETRA4 elements";

  Timer::Action timer_coo_gpu_bili(
  m_time_stats, m_use_coo_gpu ? "AssembleCooGpuBilinearOperatorTetra4" : "AssembleCooSortGpuBilinearOperatorTetra4");

  {
    Timer::Action timer_coo_gpu_build(m_time_stats,
                                      m_use_coo_gpu ? "BuildMatrixCooGpu"
                                                    : "BuildMatrixCooSortGpu");
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

  bool is_row_array_sorted = m_use_coo_sort_gpu || (mesh()->dimension() == 3 && !options()->createEdges());

  Timer::Action timer_coo_gpu_compute_add(m_time_stats, "CooGpuComputeAndAdd");

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

          // Find the index of the value in the coo matrix
          if (is_row_array_sorted) {
            auto value_index = findIndexBinarySearch(
            row_index, col_index, in_row_coo, in_col_coo, row_length);

            ax::doAtomic<ax::eAtomicOperation::Add>(in_out_val_coo(value_index),
                                                    v);
          }
          else {
            for (auto value_index = 0; value_index < row_length;
                 value_index++) {
              if (in_row_coo(value_index) == row_index &&
                  in_col_coo(value_index) == col_index) {
                ax::doAtomic<ax::eAtomicOperation::Add>(
                in_out_val_coo(value_index), v);
                break;
              }
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

ARCCORE_HOST_DEVICE Int32 findIndexBinarySearch(Int32 row, Int32 col,
                                                const auto& in_row_coo,
                                                const auto& in_col_coo,
                                                Int32 row_length)
{
  Int32 begin = 0;
  Int32 end = row_length - 1;
  Int32 i = -1; // Default value in case the index is not found

  // Binary search to find the first occurrence of `row`
  while (begin <= end) {
    Int32 mid = begin + (end - begin) / 2;

    if (row == in_row_coo(mid)) {
      // Move back to the first occurrence of the row
      while (mid > 0 && in_row_coo(mid - 1) == row) {
        mid--;
      }
      i = mid;
      break;
    }

    if (row > in_row_coo(mid)) {
      begin = mid + 1;
    }
    else {
      end = mid - 1;
    }
  }

  // If no occurrence of the row was found, return -1
  if (i == -1) {
    return -1;
  }

  // Search for the matching column in the found row's block
  while (i < row_length && in_row_coo(i) == row) {
    if (in_col_coo(i) == col) {
      return i;
    }
    i++;
  }

  // If no matching column is found, return -1
  return -1;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
