// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrGpuiBiliAssembly.hxx                                     (C) 2022-2023 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/* which handle the parallelization on Nvidia GPU                            */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include "arcane/core/IndexedItemConnectivityView.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_buildMatrixCsrGPUNodeNodeConnectivity()
{
  info() << "_buildMatrixCsrGPU: 3D mesh with node-node connectivity";

  Timer::Action timer_gpu_build_nncv(m_time_stats, "BuildMatrixCsrGPUNodeNodeConnectivity");

  Accelerator::RunQueue* queue = acceleratorMng()->defaultQueue();

  auto* connectivity_ptr = m_node_node_via_edge_connectivity.get();
  ARCANE_CHECK_POINTER(connectivity_ptr);
  IndexedNodeNodeConnectivityView nn_cv = connectivity_ptr->view();

  Int32 nb_node = nbNode();
  NumArray<uint, MDDim1> neighbors(nb_node + 1);
  neighbors[0] = 0;
  NumArray<uint, MDDim1> offsets(nb_node + 1);
  SmallSpan<uint> in_data = neighbors.to1DSmallSpan();
  SmallSpan<uint> out_data = offsets.to1DSmallSpan();

  {
    auto command = makeCommand(queue);
    command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
    {
      in_data[node_id + 1] = nn_cv.nbNode(node_id) + 1;
    };
  }
  queue->barrier();

  Accelerator::Scanner<uint> scanner;
  scanner.inclusiveSum(queue, in_data, out_data);

  {
    auto command = makeCommand(queue);

    auto out_m_matrix_row = viewInOut(command, m_csr_matrix.m_matrix_row);
    auto out_m_matrix_column = viewInOut(command, m_csr_matrix.m_matrix_column);

    command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
    {
      auto offset = out_data[node_id];
      out_m_matrix_row[node_id] = offset;
      for (auto neighbor_idx : nn_cv.nodeIds(node_id)) {
        out_m_matrix_column[offset] = neighbor_idx;
        ++offset;
      }
      out_m_matrix_column[offset] = node_id;
    };
  }
  queue->barrier();
}

void FemModule::_buildMatrixCsrGPUNodeEdgeConnectivity()
{
  info() << "_buildMatrixCsrGPU: 3D mesh with node-edge connectivity";

  Timer::Action timer_gpu_build_necv(m_time_stats, "BuildMatrixCsrGPUNodeEdgeConnectivity");

  Accelerator::RunQueue* queue = acceleratorMng()->defaultQueue();

  Int32 nb_node = nbNode();
  NumArray<uint, MDDim1> neighbors(nb_node + 1);
  neighbors[0] = 0;
  NumArray<uint, MDDim1> offsets(nb_node + 1);
  SmallSpan<uint> in_data = neighbors.to1DSmallSpan();
  SmallSpan<uint> out_data = offsets.to1DSmallSpan();

  UnstructuredMeshConnectivityView m_connectivity_view;
  m_connectivity_view.setMesh(mesh());

  ItemGenericInfoListView nodes_infos(mesh()->nodeFamily());
  IndexedNodeEdgeConnectivityView node_edge_connectivity_view = m_connectivity_view.nodeEdge();

  {
    auto command = makeCommand(queue);
    command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
    {
      in_data[node_id + 1] = node_edge_connectivity_view.nbEdge(node_id) + 1;
    };
  }
  queue->barrier();

  Accelerator::Scanner<uint> scanner;
  scanner.inclusiveSum(queue, in_data, out_data);

  {
    auto command = makeCommand(queue);

    auto out_m_matrix_column = viewInOut(command, m_csr_matrix.m_matrix_column);
    auto out_m_matrix_row = viewInOut(command, m_csr_matrix.m_matrix_row);

    IndexedEdgeNodeConnectivityView edge_node_connectivity_view = m_connectivity_view.edgeNode();

    command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
    {
      // Insertion sort: array is always sorted and small (max number is the neighbor of a node) so its linear in time
      auto offset = out_data[node_id];
      out_m_matrix_row[node_id] = offset;

      out_m_matrix_column[offset] = node_id;

      auto size = offset + 1;
      for (auto edge_id : node_edge_connectivity_view.edgeIds(node_id))
        for (auto neighbor_id : edge_node_connectivity_view.nodeIds(edge_id))
          if (neighbor_id != node_id) {
            auto position = size;
            for (auto i = offset; i < size; ++i)
              if (out_m_matrix_column[i] > neighbor_id || out_m_matrix_column[i] == -1) {
                position = i;
                break;
              }

            for (auto j = size - 1; j >= position; --j)
              out_m_matrix_column[j + 1] = out_m_matrix_column[j];

            out_m_matrix_column[position] = neighbor_id;
            size++;
          }

      /* auto offset = out_data[node_id];
      out_m_matrix_row[node_id] = offset;
      for (auto edge_id : node_edge_connectivity_view.edgeIds(node_id))
        for (auto neighbor_id : edge_node_connectivity_view.nodeIds(edge_id))
          if (neighbor_id != node_id) {
            out_m_matrix_column[offset] = neighbor_id;
            ++offset;
          }
      out_m_matrix_column[offset] = node_id; */
    };
  }
  queue->barrier();
}

void FemModule::
_buildMatrixCsrGPU()
{
  Int8 mesh_dim = mesh()->dimension();

  Int32 nb_node = nbNode();
  Int32 nb_non_zero = m_nb_edge * 2 + nb_node;
  m_csr_matrix.initialize(m_dof_family, nb_non_zero, nb_node);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (mesh_dim == 2) {
    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;
      Int32 node_dof_id = node_dof.dofId(node, 0);
      ItemLocalIdT<DoF> diagonal_entry(node_dof_id);

      m_csr_matrix.setCoordinates(diagonal_entry, diagonal_entry);

      for (Face face : node.faces()) {
        if (face.nodeId(0) == node.localId())
          m_csr_matrix.setCoordinates(diagonal_entry, node_dof.dofId(face.nodeId(1), 0));
        else
          m_csr_matrix.setCoordinates(diagonal_entry, node_dof.dofId(face.nodeId(0), 0));
      }
    }
  }
  else if (options()->createEdges()) {
    _buildMatrixCsrGPUNodeEdgeConnectivity();
  }
  else {
    _buildMatrixCsrGPUNodeNodeConnectivity();
  }
}

/*---------------------------------------------------------------------------*/
/**
* @brief Assembles the bilinear operator matrix for the FEM linear system with
 * the CSR sparse matrix format for TRIA3 elements.
 *
 * The method performs the following steps:
 *   1. Builds the CSR matrix using _buildMatrixCsrGPU.
 *   2. For each cell, computes element matrix using _computeElementMatrixTRIA3GPU.
 *   3. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrGPUBilinearOperatorTRIA3()
{

  Timer::Action timer_gpu_bili(m_time_stats, "AssembleCsrGpuBilinearOperatorTria3");

  {
    Timer::Action timer_gpu_build(m_time_stats, "CsrGpuBuildMatrix");
    // Build the csr matrix
    _buildMatrixCsrGPU();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();
  // Boucle sur les mailles déportée sur accélérateur
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
  Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
  auto in_out_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);
  UnstructuredMeshConnectivityView m_connectivity_view;
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_connectivity_view.setMesh(this->mesh());
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());

  Timer::Action timer_add_compute(m_time_stats, "CsrGpuAddComputeLoop");

  command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
  {

    Real K_e[9] = { 0 };

    _computeElementMatrixTRIA3GPU(icell, cnc, in_node_coord, K_e); // element stifness matrix
    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
    Int32 n1_index = 0;
    for (NodeLocalId node1 : cnc.nodes(icell)) {
      Int32 n2_index = 0;
      for (NodeLocalId node2 : cnc.nodes(icell)) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        double v = K_e[n1_index * 3 + n2_index];
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (nodes_infos.isOwn(node1)) {

          Int32 row = node_dof.dofId(node1, 0).localId();
          Int32 col = node_dof.dofId(node2, 0).localId();
          Int32 begin = in_row_csr[row];
          Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];

          while (begin < end) {
            if (in_col_csr[begin] == col) {
              // t is necessary to get the right type for the atomicAdd (but that means that we have more operations ?)
              // The Macro is there to avoid compilation error if not in c++ 20
              ax::doAtomic<ax::eAtomicOperation::Add>(in_out_val_csr(begin), v);
              break;
            }
            begin++;
          }
        }
        ++n2_index;
      }
      ++n1_index;
    }
  };
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system with
 * the CSR sparse matrix format for TETRA4 elements.
 *
 * The method performs the following steps:
 *   1. Builds the CSR matrix using _buildMatrixCsrGPU.
 *   2. For each cell, computes element matrix using _computeElementMatrixTETRA4GPU.
 *   3. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrGPUBilinearOperatorTETRA4()
{
  Timer::Action timer_gpu_bili(m_time_stats, "AssembleCsrGpuBilinearOperatorTetra4");

  {
    Timer::Action timer_gpu_build(m_time_stats, "CsrGpuBuildMatrix");
    _buildMatrixCsrGPU();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
  Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
  auto in_out_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);
  UnstructuredMeshConnectivityView m_connectivity_view;
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_connectivity_view.setMesh(this->mesh());
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());

  info() << "Start Add Compute Long";

  Timer::Action timer_add_compute(m_time_stats, "CsrGpuAddComputeLoop");

  command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
  {

    Real K_e[16] = { 0 };
    _computeElementMatrixTETRA4GPU(icell, cnc, in_node_coord, K_e);

    Int32 node1_idx_in_cell = 0;
    for (NodeLocalId node1_id : cnc.nodes(icell)) {

      Int32 node2_idx_in_cell = 0;
      for (NodeLocalId node2_id : cnc.nodes(icell)) {

        if (nodes_infos.isOwn(node1_id)) {
          double v = K_e[node1_idx_in_cell * 4 + node2_idx_in_cell];

          Int32 row = node1_id.localId();
          Int32 col = node2_id.localId();
          Int32 begin = in_row_csr[row];
          Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];

          // Binary search
          bool found = false;
          Int32 position = -1;
          while (!found && begin <= end) {
            Int32 mid = begin + (end - begin) / 2;
            if (col == in_col_csr(mid)) {
              found = true;
              position = mid;
            }
            if (col > in_col_csr(mid)) {
              begin = mid + 1;
            }
            if (col < in_col_csr(mid)) {
              end = mid - 1;
            }
          }

          if (position != -1)
            ax::doAtomic<ax::eAtomicOperation::Add>(in_out_val_csr(position), v);

          // while (begin < end) {
          //   if (in_col_csr[begin] == col) {
          //     ax::doAtomic<ax::eAtomicOperation::Add>(in_out_val_csr(begin), v);
          //     break;
          //   }
          //   begin++;
          // }
        }
        ++node2_idx_in_cell;
      }
      ++node1_idx_in_cell;
    }
  };
  queue->barrier();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
