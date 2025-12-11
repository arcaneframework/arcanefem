// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrGpuBiliAssembly.hxx                                      (C) 2022-2025 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/* which handle the parallelization on GPU's or CPU's                        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

#include <arcane/accelerator/Scan.h>
#include <arcane/accelerator/Sort.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static UInt64 pack(Int32 n0, Int32 n1)
{
  Int32 min = n0 > n1 ? n1 : n0;
  Int32 max = n0 > n1 ? n0 : n1;
  return ((UInt64)min << 32) | (UInt64)max;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static void unpack(UInt64 packed_edge, Int32& n0, Int32& n1)
{
  n0 = (Int32)(packed_edge >> 32);
  n1 = (Int32)(packed_edge & 0xFFFFFFFF);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModuleTestlab::_computeSortedEdges(Int8 edges_per_element, Int64 nb_edge_total, SmallSpan<UInt64>& sorted_edges_ss)
{
  auto mem_ressource = m_queue.memoryRessource();
  NumArray<UInt64, MDDim1> edges(mem_ressource);
  edges.resize(nb_edge_total);

  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cell_node_cv = m_connectivity_view.cellNode();

  {
    auto command = makeCommand(m_queue);
    auto inout_edges = viewInOut(command, edges);

    if (mesh()->dimension() == 2) {
      command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh()->allCells())
      {
        auto n0 = cell_node_cv.nodeId(cell_lid, 0);
        auto n1 = cell_node_cv.nodeId(cell_lid, 1);
        auto n2 = cell_node_cv.nodeId(cell_lid, 2);

        auto start = cell_lid * edges_per_element;
        inout_edges[start] = pack(n0, n1);
        inout_edges[start + 1] = pack(n0, n2);
        inout_edges[start + 2] = pack(n1, n2);
      };
    }
    else {
      command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh()->allCells())
      {
        auto n0 = cell_node_cv.nodeId(cell_lid, 0);
        auto n1 = cell_node_cv.nodeId(cell_lid, 1);
        auto n2 = cell_node_cv.nodeId(cell_lid, 2);
        auto n3 = cell_node_cv.nodeId(cell_lid, 3);

        auto start = cell_lid * edges_per_element;
        inout_edges[start] = pack(n0, n1);
        inout_edges[start + 1] = pack(n0, n2);
        inout_edges[start + 2] = pack(n0, n3);
        inout_edges[start + 3] = pack(n1, n2);
        inout_edges[start + 4] = pack(n1, n3);
        inout_edges[start + 5] = pack(n2, n3);
      };
    }
  }
  m_queue.barrier();

  auto sorter = Accelerator::GenericSorter(m_queue);
  SmallSpan<const UInt64> edges_ss = edges.to1DSmallSpan();
  sorter.apply(edges_ss, sorted_edges_ss);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModuleTestlab::_computeNeighbors(Int8 edges_per_element, Int64 nb_edge_total, NumArray<Int32, MDDim1>& neighbors, SmallSpan<UInt64>& sorted_edges_ss)
{
  auto command = makeCommand(m_queue);
  auto inout_neighbors = viewInOut(command, neighbors);

  command << RUNCOMMAND_LOOP1(iter, nb_edge_total)
  {
    auto [thread_id] = iter();
    auto cur_edge = sorted_edges_ss[thread_id];
    if (thread_id == (nb_edge_total - 1) || cur_edge != sorted_edges_ss[thread_id + 1]) {
      Int32 n0, n1 = 0;
      unpack(cur_edge, n0, n1);
      Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(inout_neighbors[n0], 1);
      Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(inout_neighbors[n1], 1);
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModuleTestlab::_computeRowIndex(Int8 edges_per_element, Int64 nb_edge_total, SmallSpan<UInt64>& sorted_edges_ss) {
    auto mem_ressource = m_queue.memoryRessource();
    NumArray<Int32, MDDim1> neighbors(mem_ressource);
    neighbors.resize(mesh()->nbNode());
    neighbors.fill(1, &m_queue);
    _computeNeighbors(edges_per_element, nb_edge_total, neighbors, sorted_edges_ss);

    Accelerator::Scanner<Int32> scanner;
    SmallSpan<Int32> neighbors_ss = neighbors.to1DSmallSpan();
    scanner.exclusiveSum(&m_queue, neighbors_ss, m_csr_matrix.m_matrix_row.to1DSmallSpan());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static void registerEdgeInColumns(Int32 src, Int32 dst, Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> offsets, Accelerator::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> row_index, Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> columns)
{
  Int32 start = row_index[src];
  Int32 offset = Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(offsets[src], 1);
  columns[start + offset] = dst;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModuleTestlab::_computeColumns(Int8 edges_per_element, Int64 nb_edge_total, SmallSpan<uint64_t>& sorted_edges_ss) {
    auto nb_node = mesh()->nbNode();

    {
      auto command = makeCommand(m_queue);
      auto in_row_index = viewIn(command, m_csr_matrix.m_matrix_row);
      auto inout_columns = viewInOut(command, m_csr_matrix.m_matrix_column);

      command << RUNCOMMAND_LOOP1(iter, nb_node)
      {
        auto [n] = iter();
        auto start = in_row_index[n];
        inout_columns[start] = n;
      };
    }
    m_queue.barrier();

    auto mem_ressource = m_queue.memoryRessource();
    NumArray<Int32, MDDim1> offsets(mem_ressource);
    offsets.resize(nb_node);
    offsets.fill(1, &m_queue);

    {
      auto command = makeCommand(m_queue);
      auto inout_columns = viewInOut(command, m_csr_matrix.m_matrix_column);
      auto in_row_index = viewIn(command, m_csr_matrix.m_matrix_row);
      auto inout_offsets = viewInOut(command, offsets);

      command << RUNCOMMAND_LOOP1(iter, nb_edge_total)
      {
        auto [thread_id] = iter();
        auto cur_edge = sorted_edges_ss[thread_id];
        if (thread_id == (nb_edge_total - 1) || cur_edge != sorted_edges_ss[thread_id + 1]) {
          Int32 n0, n1 = 0;
          unpack(cur_edge, n0, n1);
          registerEdgeInColumns(n0, n1, inout_offsets, in_row_index, inout_columns);
          registerEdgeInColumns(n1, n0, inout_offsets, in_row_index, inout_columns);
        }
      };
    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModuleTestlab::_computeSparsity()
{

  Int8 mesh_dim = mesh()->dimension();

  Int32 nb_node = nbNode();
  Int32 nb_non_zero = nb_node + 2 * (mesh_dim == 2 ? nbFace() : m_nb_edge);
  m_csr_matrix.initialize(m_dof_family, nb_non_zero, nb_node, m_queue);

  auto edges_per_element = mesh_dim == 2 ? 3 : 6;
  auto nb_edge_total = mesh()->nbCell() * edges_per_element;

  auto mem_ressource = m_queue.memoryRessource();
  NumArray<UInt64, MDDim1> sorted_edges(mem_ressource);
  sorted_edges.resize(nb_edge_total);
  auto sorted_edges_ss = sorted_edges.to1DSmallSpan();

  _computeSortedEdges(edges_per_element, nb_edge_total, sorted_edges_ss);
  _computeRowIndex(edges_per_element, nb_edge_total, sorted_edges_ss);
  _computeColumns(edges_per_element, nb_edge_total, sorted_edges_ss);
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

void FemModuleTestlab::
_assembleCsrGPUBilinearOperatorTRIA3()
{

  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_Csr_Gpu");

  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    // Build the csr matrix
    _computeSparsity();
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

  Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");

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

void FemModuleTestlab::
_assembleCsrGPUBilinearOperatorTETRA4()
{
  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_Csr_Gpu");

  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _computeSparsity();
  }

  {
    auto command = makeCommand(m_queue);

    Int32 row_csr_size = m_csr_matrix.m_matrix_row.extent0();
    Int32 col_csr_size = m_csr_matrix.m_matrix_column.extent0();

    auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
    auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
    auto inout_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);

    auto in_node_coord = ax::viewIn(command, m_node_coord);

    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cell_node_connectivity_view = m_connectivity_view.cellNode();

    ItemGenericInfoListView nodes_infos(mesh()->nodeFamily());

    Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");
    ax::ProfileRegion ps_region(m_queue,"AddAndComputeBilinearTetra4",0x00FF7F);
    command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
    {

      Real K_e[16]{};
      _computeElementMatrixTETRA4GPU(icell, cell_node_connectivity_view, in_node_coord, K_e);

      Int32 node1_idx_in_cell = 0;
      for (NodeLocalId node1_id : cell_node_connectivity_view.nodes(icell)) {

        Int32 node2_idx_in_cell = 0;
        for (NodeLocalId node2_id : cell_node_connectivity_view.nodes(icell)) {

          if (nodes_infos.isOwn(node1_id)) {
            double v = K_e[node1_idx_in_cell * 4 + node2_idx_in_cell];

            Int32 row = node1_id.localId();
            Int32 col = node2_id.localId();
            Int32 begin = in_row_csr[row];

            Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];

            while (begin < end) {
              if (in_col_csr[begin] == col) {
                ax::doAtomic<ax::eAtomicOperation::Add>(inout_val_csr(begin), v);
                break;
              }
              begin++;
            }
          }
          ++node2_idx_in_cell;
        }
        ++node1_idx_in_cell;
      }
    };
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
