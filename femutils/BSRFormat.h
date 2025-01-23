// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BSRFormat.h                                                 (C) 2022-2025 */
/*                                                                           */
/* BSRFormat class definition.                                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef BSRFORMAT_H
#define BSRFORMAT_H

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <ios>
#include <iomanip>

#include <arccore/trace/TraceAccessor.h>

#include <arccore/base/NotImplementedException.h>
#include <arccore/base/ArgumentException.h>
#include <arccore/base/ArccoreGlobal.h>

#include <arcane/core/IIndexedIncrementalItemConnectivityMng.h>
#include <arcane/core/IIndexedIncrementalItemConnectivity.h>
#include <arcane/core/UnstructuredMeshConnectivity.h>
#include <arcane/core/IIncrementalItemConnectivity.h>
#include <arcane/core/IndexedItemConnectivityView.h>
#include <arcane/core/VariableTypedef.h>
#include <arcane/core/ItemEnumerator.h>
#include <arcane/core/ItemTypes.h>
#include <arcane/core/MeshUtils.h>
#include <arcane/core/DataView.h>
#include <arcane/core/IMesh.h>
#include <arcane/core/Item.h>

#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/ArrayLayout.h>
#include <arcane/utils/UtilsTypes.h>
#include <arcane/utils/NumArray.h>
#include <arcane/utils/MDDim.h>

#include <arcane/accelerator/RunCommandEnumerate.h>
#include <arcane/accelerator/RunCommandLoop.h>
#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/accelerator/NumArrayViews.h>
#include <arcane/accelerator/GenericSorter.h>
#include <arcane/accelerator/Atomic.h>
#include <arcane/accelerator/Scan.h>
#include <arcane/accelerator/Sort.h>

#include "DoFLinearSystem.h"
#include "CsrFormatMatrix.h"
#include "FemDoFsOnNodes.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/**
 * @brief A class representing a Block Sparse Row (BSR) matrix.
 *
 * The `BSRMatrix` class stores sparse matrices in the Block Sparse Row (BSR)
 * format, suitable for block-structured computations.
 *
 * The BSR format groups non-zero elements into fixed-size blocks, enabling
 * efficient storage and computation for multi-DoFs meshes.
 * Below is an example for a 6x6 matrix with a block size of 2x2 (2 DoFs per node).
 *
 *     Original Matrix (6x6, Block Size = 2x2):
 *
 *     [ 1  2  0  0  0  0 ]
 *     [ 3  4  0  0  0  0 ]
 *     [ 0  0  5  6  7  8 ]
 *     [ 0  0  9 10 11 12 ]
 *     [ 0  0  0  0 13 14 ]
 *     [ 0  0  0  0 15 16 ]
 *
 *     BSR Representation:
 *       - Values:
 *          if "Ordered per block":
 *            [ 1, 2, 3, 4, 5, 6, 9, 10, 7, 8, 11, 12, 13, 14, 15, 16 ]
 *          if "Ordered per row" (same as CSR):
 *            [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
 *       - Columns:   [ 0, 1, 2, 2 ]  // Column indices of each block
 *       - Row Index: [ 0, 1, 3 ]     // Starting indices of blocks per row
 *
 * @note The class supports flexible memory management for computations on
 * different hardware backends (e.g., CPU, GPU). Having the value array
 * stored in per row (which the way CSR format stores matrices) allow us to
 * pass it directly to HYPRE (which support CSR but not BSR), whithout the need
 * for a translation step on it.
 */
/*---------------------------------------------------------------------------*/

template <int BLOCK_SIZE>
class BSRMatrix : public TraceAccessor
{
 public:

  BSRMatrix(ITraceMng* tm, const eMemoryRessource& mem_ressource, RunQueue& queue)
  : TraceAccessor(tm)
  , m_values(mem_ressource)
  , m_columns(mem_ressource)
  , m_row_index(mem_ressource)
  , m_queue(queue) {};

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void initialize(Int32 nb_non_zero_value, Int32 nb_col, Int32 nb_row, bool order_values_per_block)
  {
    if (BLOCK_SIZE <= 0 || nb_non_zero_value <= 0 || nb_row <= 0)
      ARCANE_THROW(ArgumentException, "BSRMatrix(initialize): arguments should be positive and not null (block_size={0}, nb_non_zero_value={1} and nb_row={2})", BLOCK_SIZE, nb_non_zero_value, nb_row);

    if (BLOCK_SIZE > nb_row)
      ARCANE_THROW(ArgumentException, "BSRMatrix(initialize): block_size should be less than nb_row");

    info() << "BSRMatrix(initialize): Initialize BSRMatrix with block_size=" << BLOCK_SIZE << ", nb_non_zero_value=" << nb_non_zero_value << ", nb_col=" << nb_col << ", nb_row=" << nb_row << ", order_values_per_block=" << std::boolalpha << order_values_per_block;

    m_order_values_per_block = order_values_per_block;

    m_nb_non_zero_value = nb_non_zero_value;
    m_nb_col = nb_col;
    m_nb_row = nb_row;

    m_values.resize(nb_non_zero_value);
    m_values.fill(0, &m_queue);
    m_columns.resize(nb_col);
    m_row_index.resize(nb_row);
    m_nb_nz_per_row.resize(nb_row);
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  Int32 findValueIndex(DoFLocalId row, DoFLocalId col) const
  {
    auto block_row = row / BLOCK_SIZE;
    auto block_col = col / BLOCK_SIZE;

    auto block_start = m_row_index[block_row];
    auto block_end = (block_row == m_nb_row - 1) ? m_nb_col : m_row_index[block_row + 1];

    auto row_offset = row % BLOCK_SIZE;
    auto col_offset = col % BLOCK_SIZE;

    constexpr int BLOCK_SIZE_SQ = BLOCK_SIZE * BLOCK_SIZE;
    auto block_start_in_value = block_start * BLOCK_SIZE_SQ;
    auto col_index = 0;
    while (block_start < block_end) {
      if (m_columns[block_start] == block_col) {
        if (!m_order_values_per_block)
          return block_start_in_value + (BLOCK_SIZE * col_index) + (row_offset * BLOCK_SIZE * m_nb_nz_per_row[block_row]) + col_offset;
        else
          return (block_start * BLOCK_SIZE_SQ) + ((row_offset * BLOCK_SIZE) + col_offset);
      }
      ++block_start;
      ++col_index;
    }

    throw std::runtime_error("BSRMatrix(findValueIndex): Value not found");
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  Real getValue(DoFLocalId row, DoFLocalId col) const
  {
    auto value_idx = findValueIndex(row, col);
    return m_values[value_idx];
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void setValue(DoFLocalId row, DoFLocalId col, Real value)
  {
    auto value_idx = findValueIndex(row, col);
    m_values[value_idx] = value;
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void toCsr(CsrFormat* csr_matrix)
  {
    info() << "BSRMatrix(toCsr): Convert matrix to CSR";

    auto startTime = platform::getRealTime();

    auto nb_block_rows = m_row_index.extent0();
    auto nb_rows = nb_block_rows * BLOCK_SIZE;
    auto total_non_zero_elements = m_nb_non_zero_value;

    csr_matrix->initialize(nullptr, total_non_zero_elements, nb_rows, m_queue);

    if constexpr (BLOCK_SIZE == 1) {
      csr_matrix->m_matrix_row = m_row_index;
      csr_matrix->m_matrix_column = m_columns;
      csr_matrix->m_matrix_rows_nb_column = m_nb_nz_per_row;
    }
    else {
      // Translate `row_index`
      csr_matrix->m_matrix_row[0] = 0;
      auto offset = 1;
      for (auto i = 0; i < m_row_index.extent0() - 1; ++i) {
        for (auto j = 0; j < BLOCK_SIZE; ++j) {
          auto start = m_row_index[i];
          auto end = i == m_row_index.extent0() - 1 ? m_nb_col : m_row_index[i + 1];
          auto nombre_de_dof_dans_la_rangee = (end - start) * BLOCK_SIZE;
          csr_matrix->m_matrix_row[offset] = csr_matrix->m_matrix_row[offset - 1] + nombre_de_dof_dans_la_rangee;
          offset++;
        }
      }

      // Translate `columns`
      offset = 0;
      for (auto i = 0; i < m_row_index.extent0(); ++i) {
        for (auto j = 0; j < BLOCK_SIZE; ++j) {
          auto start = m_row_index[i];
          auto end = i == m_row_index.extent0() - 1 ? m_nb_col : m_row_index[i + 1];
          for (auto block_index = start; block_index < end; ++block_index) {
            for (auto k = 0; k < BLOCK_SIZE; ++k) {
              csr_matrix->m_matrix_column[offset] = m_columns[block_index] * BLOCK_SIZE + k;
              offset++;
            }
          }
        }
      }

      // Translate `nb_nz_per_row`
      offset = 0;
      for (auto i = 0; i < m_nb_nz_per_row.extent0(); ++i) {
        for (auto j = 0; j < BLOCK_SIZE; ++j) {
          csr_matrix->m_matrix_rows_nb_column[offset++] = m_nb_nz_per_row[i] * BLOCK_SIZE;
        }
      }

      // Fill last entry of CSR matrix row_index array
      auto row_index_size = csr_matrix->m_matrix_row.extent0();
      csr_matrix->m_matrix_row[row_index_size - 1] = csr_matrix->m_matrix_row[row_index_size - 2] + csr_matrix->m_matrix_rows_nb_column[csr_matrix->m_matrix_rows_nb_column.extent0() - 1];
    }

    // NOTE: If we don't want to keep bsr matrix coefficients we could move the data instead of copying it.
    csr_matrix->m_matrix_value = m_values;
    info() << "[ArcaneFem-Timer] Time to translate BSR to CSR = " << (platform::getRealTime() - startTime);
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void toLinearSystem(DoFLinearSystem& linear_system)
  {
    info() << "BSRMatrix(toLinearSystem): Translate matrix to linear system using `matrixAddValue`";

    for (auto row = 0; row < m_nb_row; ++row) {
      auto row_start = m_row_index[row];
      auto row_end = (row + 1 < m_nb_row) ? m_row_index[row + 1] : m_nb_col;
      for (auto block_idx = row_start; block_idx < row_end; ++block_idx) {
        auto col = m_columns[block_idx];
        for (auto i = 0; i < BLOCK_SIZE; ++i) {
          for (auto j = 0; j < BLOCK_SIZE; ++j) {
            auto global_row = (row * BLOCK_SIZE) + i;
            auto global_col = (col * BLOCK_SIZE) + j;
            auto value = getValue(DoFLocalId(global_row), DoFLocalId(global_col));
            linear_system.matrixAddValue(DoFLocalId(global_row), DoFLocalId(global_col), value);
          }
        }
      }
    }
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void dump(std::string filename)
  {
    info() << "BSRMatrix(dump): Dump matrix in " << std::quoted(filename);
    ofstream file(filename);

    file << "size :" << nbNz() << "\n";
    for (auto i = 0; i < nbRow(); ++i) {
      file << m_row_index(i) << " ";
      for (Int32 j = m_row_index(i) + 1; (i + 1 < m_row_index.dim1Size() && j < m_row_index(i + 1)) || (i + 1 == m_row_index.dim1Size() && j < m_columns.dim1Size()); j++)
        file << "  ";
    }
    file << "\n";

    for (auto i = 0; i < m_columns.extent0(); ++i)
      file << m_columns(i) << " ";
    file << "\n";

    for (auto i = 0; i < nbNz(); ++i)
      file << m_values(i) << " ";
    file << "\n";

    file.close();
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  bool orderValuePerBlock() { return m_order_values_per_block; }
  Int32 nbNz() { return m_nb_non_zero_value; };
  Int32 nbCol() { return m_nb_col; };
  Int32 nbRow() { return m_nb_row; };

  NumArray<Real, MDDim1>& values() { return m_values; }
  NumArray<Int32, MDDim1>& columns() { return m_columns; }
  NumArray<Int32, MDDim1>& rowIndex() { return m_row_index; }
  NumArray<Int32, MDDim1>& nbNzPerRow() { return m_nb_nz_per_row; }

 private:

  bool m_order_values_per_block = true;

  Int32 m_nb_non_zero_value;
  Int32 m_nb_col;
  Int32 m_nb_row;

  NumArray<Real, MDDim1> m_values;
  NumArray<Int32, MDDim1> m_columns;
  NumArray<Int32, MDDim1> m_row_index;
  NumArray<Int32, MDDim1> m_nb_nz_per_row;

  RunQueue& m_queue;
};

/*---------------------------------------------------------------------------*/
/**
 * @brief A class for assembling Block Sparse Row (BSR) matrices from mesh data.
 *
 * The `BSRFormat` class is designed to construct and manage BSR matrices using 
 * mesh connectivity and degree-of-freedom (DoF) information.
 *
 * In 2D, the sparsity is computed based on node-face connectivity, while in 3D, 
 * node-node connectivity is used.
 *
 * @note This class uses Arcane's accelerator api and will use GPU is possible.
 * It uses a `BSRMatrix` under the hood for representation and operations.
 */
/*---------------------------------------------------------------------------*/

template <int NB_DOF>
class BSRFormat : public TraceAccessor
{
 public:

  BSRFormat(ITraceMng* tm, RunQueue& queue, const FemDoFsOnNodes& dofs_on_nodes)
  : TraceAccessor(tm)
  , m_dofs_on_nodes(dofs_on_nodes)
  , m_queue(queue)
  , m_bsr_matrix(tm, queue.memoryRessource(), queue)
  , m_csr_matrix(tm) {};

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

 private:

  Int64 computeNbEdge(IMesh* mesh)
  {
    Int64 nb_edge = 0;
    if (mesh->dimension() == 2)
      nb_edge = mesh->nbFace();
    else {
      auto nn_via_edge_cv = MeshUtils::computeNodeNodeViaEdgeConnectivity(mesh, "NodeNodeViaEdge");
      IndexedNodeNodeConnectivityView nn_cv = nn_via_edge_cv->view();
      ENUMERATE_NODE (inode, mesh->allNodes()) {
        Node node = *inode;
        nb_edge += nn_cv.nbNode(node);
      }
      nb_edge /= 2;
    }
    return nb_edge;
  }

 public:

  void initialize(IMesh* mesh, bool does_linear_system_use_csr, bool use_atomic_free = false)
  {
    ARCANE_CHECK_POINTER(mesh);

    if (mesh->dimension() != 2 && mesh->dimension() != 3)
      ARCANE_THROW(NotImplementedException, "BSRFormat(initialize): Only supports 2D and 3D");

    auto startTime = platform::getRealTime();

    m_mesh = mesh;
    Int64 nb_edge = computeNbEdge(mesh);
    Int32 nb_node = m_mesh->nbNode();
    Int32 nb_col = 2 * nb_edge + nb_node;
    Int32 nb_non_zero_value = (NB_DOF * NB_DOF) * (2 * nb_edge + nb_node);

    m_use_csr_in_linear_system = does_linear_system_use_csr;
    bool order_values_per_block = !does_linear_system_use_csr;
    m_bsr_matrix.initialize(nb_non_zero_value, nb_col, nb_node, order_values_per_block);
    m_use_atomic_free = use_atomic_free;
    info() << "[ArcaneFem-Timer] Time to initialize BSR format = " << (platform::getRealTime() - startTime);
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void toLinearSystem(DoFLinearSystem& linear_system)
  {
    auto startTime = platform::getRealTime();
    if (m_use_csr_in_linear_system) {
      if (!linear_system.hasSetCSRValues())
        ARCANE_THROW(ArgumentException, "BSRFormat(toLinearSystem): Linear system was set to use CSR but is incompatible");

      m_bsr_matrix.toCsr(&m_csr_matrix);

      info() << "BSRFormat(toLinearSystem): Set CSR values into linear system";
      CSRFormatView csr_view(m_csr_matrix.m_matrix_row.to1DSpan(), m_csr_matrix.m_matrix_rows_nb_column.to1DSpan(), m_csr_matrix.m_matrix_column.to1DSpan(), m_csr_matrix.m_matrix_value.to1DSpan());
      linear_system.setCSRValues(csr_view);
    }
    else
      m_bsr_matrix.toLinearSystem(linear_system);
    info() << "[ArcaneFem-Timer] Time to translate to linear system the BSR matrix = " << (platform::getRealTime() - startTime);
  };

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void computeNzPerRowArray()
  {
    info() << "BSRFormat(computeNzPerRowArray): Compute nb_nz_per_row BSR matrix array";

    auto startTime = platform::getRealTime();
    auto command = makeCommand(m_queue);
    NumArray<Int32, MDDim1>& nb_nz_per_row = m_bsr_matrix.nbNzPerRow();
    auto inout_nb_nz_per_row = viewInOut(command, nb_nz_per_row);
    auto inout_row_index = viewInOut(command, m_bsr_matrix.rowIndex());

    auto nb_row = m_bsr_matrix.nbRow();
    // Copy `row_index` in `nz_per_row`
    // NOTE: Does Arcane provide a copy-on-gpu primitive ?
    {
      auto command = makeCommand(m_queue);
      command << RUNCOMMAND_LOOP1(iter, nb_row)
      {
        auto [i] = iter();
        inout_nb_nz_per_row[i] = inout_row_index[i];
      };
    }
    m_queue.barrier();

    {
      auto nb_col = m_bsr_matrix.nbCol();
      auto command = makeCommand(m_queue);
      command << RUNCOMMAND_LOOP1(iter, nb_row)
      {
        auto [i] = iter();
        auto x = inout_nb_nz_per_row[i];
        if (i == nb_row - 1)
          inout_nb_nz_per_row[i] = nb_col - inout_row_index[i];
        else
          inout_nb_nz_per_row[i] = inout_row_index[i + 1] - x;
      };
    }
    m_queue.barrier();

    info() << "[ArcaneFem-Timer] Time to compute nb_nz_per_row = " << (platform::getRealTime() - startTime);
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void computeNeighborsAtomicFree(SmallSpan<Int32>& neighbors_ss)
  {
    auto command = makeCommand(m_queue);
    if (m_mesh->dimension() == 2) {
      UnstructuredMeshConnectivityView connectivity_view(m_mesh);
      auto node_face_cv = connectivity_view.nodeFace();
      command << RUNCOMMAND_ENUMERATE(Node, node_id, m_mesh->allNodes())
      {
        neighbors_ss[node_id] = node_face_cv.nbFace(node_id) + 1;
      };
    }
    else {
      auto connectivity_mng = m_mesh->indexedConnectivityMng();
      auto connectivity_ptr = connectivity_mng->findOrCreateConnectivity(m_mesh->nodeFamily(), m_mesh->nodeFamily(), "NodeNodeViaEdge");
      IndexedNodeNodeConnectivityView node_node_cv = connectivity_ptr->view();

      command << RUNCOMMAND_ENUMERATE(Node, node_id, m_mesh->allNodes())
      {
        neighbors_ss[node_id] = node_node_cv.nbNode(node_id) + 1;
      };
    }
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void computeRowIndexAtomicFree()
  {
    auto mem_ressource = m_queue.memoryRessource();
    NumArray<Int32, MDDim1> neighbors(mem_ressource);
    neighbors.resize(m_mesh->nbNode());
    SmallSpan<Int32> neighbors_ss = neighbors.to1DSmallSpan();
    computeNeighborsAtomicFree(neighbors_ss);

    Accelerator::Scanner<Int32> scanner;
    scanner.exclusiveSum(&m_queue, neighbors_ss, m_bsr_matrix.rowIndex().to1DSmallSpan());
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void computeColumnsAtomicFree()
  {
    auto command = makeCommand(m_queue);
    auto inout_columns = viewInOut(command, m_bsr_matrix.columns());
    auto row_index = m_bsr_matrix.rowIndex().to1DSmallSpan();

    if (m_mesh->dimension() == 2) {
      UnstructuredMeshConnectivityView connectivity_view(m_mesh);
      auto node_face_cv = connectivity_view.nodeFace();
      auto face_node_cv = connectivity_view.faceNode();

      command << RUNCOMMAND_ENUMERATE(Node, node_id, m_mesh->allNodes())
      {
        auto offset = row_index[node_id];

        for (auto face_lid : node_face_cv.faceIds(node_id)) {
          auto nodes = face_node_cv.nodes(face_lid);
          inout_columns[offset] = nodes[0] == node_id ? nodes[1] : nodes[0];
          ++offset;
        }

        inout_columns[offset] = node_id;
      };
    }
    else {
      auto connectivity_mng = m_mesh->indexedConnectivityMng();
      auto connectivity_ptr = connectivity_mng->findOrCreateConnectivity(m_mesh->nodeFamily(), m_mesh->nodeFamily(), "NodeNodeViaEdge");
      IndexedNodeNodeConnectivityView node_node_cv = connectivity_ptr->view();

      command << RUNCOMMAND_ENUMERATE(Node, node_id, m_mesh->allNodes())
      {
        auto offset = row_index[node_id];

        for (auto neighbor_idx : node_node_cv.nodeIds(node_id)) {
          inout_columns[offset] = neighbor_idx;
          ++offset;
        }

        inout_columns[offset] = node_id;
      };
    }
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void computeSparsityAtomicFree()
  {
    info() << "BSRFormat(computeSparsityAtomicFree): Compute sparsity of BSR matrix using Arcane connectivities (e.g node-face for 2D and node-node for 3D)";
    auto startTime = platform::getRealTime();

    computeRowIndexAtomicFree();
    computeColumnsAtomicFree();

    info() << "[ArcaneFem-Timer] Time to compute the sparsity of BSR matrix (with Arcane connectivities implementation) = " << (platform::getRealTime() - startTime);

    if (m_use_csr_in_linear_system)
      computeNzPerRowArray();
  }

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

  void computeSortedEdges(Int8 edges_per_element, Int64 nb_edge_total, SmallSpan<UInt64>& sorted_edges_ss)
  {
    auto mem_ressource = m_queue.memoryRessource();
    NumArray<UInt64, MDDim1> edges(mem_ressource);
    edges.resize(nb_edge_total);

    UnstructuredMeshConnectivityView m_connectivity_view(m_mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();

    {
      auto command = makeCommand(m_queue);
      auto inout_edges = viewInOut(command, edges);

      if (m_mesh->dimension() == 2) {
        command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, m_mesh->allCells())
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
        command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, m_mesh->allCells())
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

  void computeNeighbors(Int8 edges_per_element, Int64 nb_edge_total, NumArray<Int32, MDDim1>& neighbors, SmallSpan<UInt64>& sorted_edges_ss)
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

  void computeRowIndex(Int8 edges_per_element, Int64 nb_edge_total, SmallSpan<UInt64>& sorted_edges_ss)
  {
    auto mem_ressource = m_queue.memoryRessource();
    NumArray<Int32, MDDim1> neighbors(mem_ressource);
    neighbors.resize(m_mesh->nbNode());
    neighbors.fill(1, &m_queue);
    computeNeighbors(edges_per_element, nb_edge_total, neighbors, sorted_edges_ss);

    Accelerator::Scanner<Int32> scanner;
    SmallSpan<Int32> neighbors_ss = neighbors.to1DSmallSpan();
    scanner.exclusiveSum(&m_queue, neighbors_ss, m_bsr_matrix.rowIndex().to1DSmallSpan());
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

  void computeColumns(Int8 edges_per_element, Int64 nb_edge_total, SmallSpan<uint64_t>& sorted_edges_ss)
  {
    auto nb_node = m_mesh->nbNode();

    {
      auto command = makeCommand(m_queue);
      auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
      auto inout_columns = viewInOut(command, m_bsr_matrix.columns());

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
      auto inout_columns = viewInOut(command, m_bsr_matrix.columns());
      auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
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
  };

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void computeSparsityAtomic()
  {
    info() << "BSRFormat(computeSparsityAtomic): Computing sparsity of BSR matrix without Arcane connectivities (e.g with atomics)...";
    auto startTime = platform::getRealTime();

    auto edges_per_element = m_mesh->dimension() == 2 ? 3 : 6;
    auto nb_edge_total = m_mesh->nbCell() * edges_per_element;

    auto mem_ressource = m_queue.memoryRessource();
    NumArray<UInt64, MDDim1> sorted_edges(mem_ressource);
    sorted_edges.resize(nb_edge_total);
    auto sorted_edges_ss = sorted_edges.to1DSmallSpan();

    computeSortedEdges(edges_per_element, nb_edge_total, sorted_edges_ss);
    computeRowIndex(edges_per_element, nb_edge_total, sorted_edges_ss);
    computeColumns(edges_per_element, nb_edge_total, sorted_edges_ss);

    info() << "[ArcaneFem-Timer] Time to compute the sparsity of BSR matrix (without Arcane connectivities implementation, e.g with atomics) = " << (platform::getRealTime() - startTime);

    if (m_use_csr_in_linear_system)
      computeNzPerRowArray();
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  void computeSparsity()
  {
    if (m_use_atomic_free)
      computeSparsityAtomicFree();
    else
      computeSparsityAtomic();
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  template <class Function>
  void assembleBilinearOrderedPerBlock(Function compute_element_matrix)
  {
    UnstructuredMeshConnectivityView m_connectivity_view(m_mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();

    ItemGenericInfoListView nodes_infos(m_mesh->nodeFamily());

    constexpr int NB_DOF_SQ = NB_DOF * NB_DOF;
    auto matrix_nb_row = m_bsr_matrix.nbRow();
    auto matrix_nb_column = m_bsr_matrix.nbCol();
    auto matrix_nb_nz = m_bsr_matrix.nbNz();

    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    auto command = makeCommand(m_queue);
    auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
    auto in_columns = viewIn(command, m_bsr_matrix.columns());
    auto inout_values = viewInOut(command, m_bsr_matrix.values());
    auto in_nz_per_row = viewIn(command, m_bsr_matrix.nbNzPerRow());

    command << RUNCOMMAND_ENUMERATE(Cell, cell, m_mesh->allCells())
    {
      auto element_matrix = compute_element_matrix(cell);

      auto cur_row_node_idx = 0;
      for (NodeLocalId row_node_lid : cell_node_cv.nodes(cell)) {
        auto cur_col_node_idx = 0;
        for (NodeLocalId col_node_lid : cell_node_cv.nodes(cell)) {
          if (nodes_infos.isOwn(row_node_lid)) {
            auto begin = in_row_index[row_node_lid];
            auto end = (row_node_lid == matrix_nb_row - 1) ? matrix_nb_column : in_row_index[row_node_lid + 1];
            while (begin < end) {
              if (in_columns[begin] == col_node_lid) {
                auto block_start = begin * NB_DOF_SQ;
                for (auto i = 0; i < NB_DOF; ++i) {
                  for (auto j = 0; j < NB_DOF; ++j) {
                    double value = element_matrix(NB_DOF * cur_row_node_idx + i, NB_DOF * cur_col_node_idx + j);
                    Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(inout_values[block_start + (i * NB_DOF + j)], value);
                  }
                }
                break;
              }
              ++begin;
            }
          }
          ++cur_col_node_idx;
        }
        ++cur_row_node_idx;
      }
    };
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  template <class Function> void assembleBilinearOrderedPerRow(Function compute_element_matrix)
  {
    UnstructuredMeshConnectivityView m_connectivity_view(m_mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();

    ItemGenericInfoListView nodes_infos(m_mesh->nodeFamily());

    constexpr int NB_DOF_SQ = NB_DOF * NB_DOF;
    auto matrix_nb_row = m_bsr_matrix.nbRow();
    auto matrix_nb_column = m_bsr_matrix.nbCol();
    auto matrix_nb_nz = m_bsr_matrix.nbNz();

    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    auto command = makeCommand(m_queue);
    auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
    auto in_columns = viewIn(command, m_bsr_matrix.columns());
    auto inout_values = viewInOut(command, m_bsr_matrix.values());
    auto in_nz_per_row = viewIn(command, m_bsr_matrix.nbNzPerRow());

    command << RUNCOMMAND_ENUMERATE(Cell, cell, m_mesh->allCells())
    {
      auto element_matrix = compute_element_matrix(cell);

      auto cur_row_node_idx = 0;
      for (NodeLocalId row_node_lid : cell_node_cv.nodes(cell)) {
        auto cur_col_node_idx = 0;
        for (NodeLocalId col_node_lid : cell_node_cv.nodes(cell)) {
          if (nodes_infos.isOwn(row_node_lid)) {

            auto g_block_start = in_row_index[row_node_lid] * NB_DOF_SQ;
            auto begin = in_row_index[row_node_lid];
            auto end = (row_node_lid == matrix_nb_row - 1) ? matrix_nb_column : in_row_index[row_node_lid + 1];

            auto x = 0;
            while (begin < end) {
              if (in_columns[begin] == col_node_lid) {
                for (auto i = 0; i < NB_DOF; ++i) {
                  for (auto j = 0; j < NB_DOF; ++j) {
                    double value = element_matrix(NB_DOF * cur_row_node_idx + i, NB_DOF * cur_col_node_idx + j);
                    auto l_block_start = g_block_start + (NB_DOF * x);
                    if (i != 0)
                      l_block_start += (NB_DOF * in_nz_per_row[row_node_lid]);
                    Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(inout_values[l_block_start + j], value);
                  }
                }
              }
              ++begin;
              ++x;
            }
          }
          ++cur_col_node_idx;
        }
        ++cur_row_node_idx;
      }
    };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Assembles the global BSR matrix for a bilinear operator, cell-wise.
   *
   * This function constructs the Block Sparse Row (BSR) matrix by iterating over mesh cells
   * and accumulating contributions from element-level matrices. It uses a user-defined
   * callback to compute the local matrix for each cell and updates the global matrix
   * with atomic operations to ensure parallel safety.
   *
   * ### Parameters:
   * - `compute_element_matrix`: A callback function that computes the local matrix
   *   for a given cell. It must return a matrix compatible with the expected block size.
   *
   * ### Key Details:
   * - Uses cell-to-node connectivity from the mesh and DoF mappings for assembly.
   * - Employs atomic operations to safely handle concurrent updates in parallel execution.
   * - Handles block structure updates for `blockSize x blockSize` dimensions.
   *
   */
  /*---------------------------------------------------------------------------*/

  template <class Function> void assembleBilinearAtomic(Function compute_element_matrix)
  {
    info() << "BSRFormat(assembleBilinearAtomic): Integrating over elements...";
    auto startTime = platform::getRealTime();

    if (m_bsr_matrix.orderValuePerBlock())
      assembleBilinearOrderedPerBlock(compute_element_matrix);
    else
      assembleBilinearOrderedPerRow(compute_element_matrix);

    info() << "[ArcaneFem-Timer] Time to assemble (atomic implementation) BSR matrix = " << (platform::getRealTime() - startTime);
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  template <class Function>
  void assembleBilinearOrderedPerBlockAtomicFree(Function compute_element_matrix)
  {
    UnstructuredMeshConnectivityView m_connectivity_view(m_mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();
    auto node_cell_cv = m_connectivity_view.nodeCell();

    ItemGenericInfoListView nodes_infos(m_mesh->nodeFamily());

    constexpr int NB_DOF_SQ = NB_DOF * NB_DOF;
    auto matrix_nb_row = m_bsr_matrix.nbRow();
    auto matrix_nb_column = m_bsr_matrix.nbCol();
    auto matrix_nb_nz = m_bsr_matrix.nbNz();

    auto command = makeCommand(m_queue);
    auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
    auto in_columns = viewIn(command, m_bsr_matrix.columns());
    auto inout_values = viewInOut(command, m_bsr_matrix.values());

    command << RUNCOMMAND_ENUMERATE(Node, row_node, m_mesh->allNodes())
    {
      auto cur_row_node_idx = 0;
      for (auto cell : node_cell_cv.cells(row_node)) {

        // Find the index of the node in the current cell
        for (Int32 i = 1; i <= 3; ++i) {
          if (row_node == cell_node_cv.nodeId(cell, i)) {
            cur_row_node_idx = i;
            break;
          }
          cur_row_node_idx = 0;
        }

        auto element_matrix = compute_element_matrix(cell);

        auto cur_col_node_idx = 0;
        for (NodeLocalId col_node_lid : cell_node_cv.nodes(cell)) {
          if (nodes_infos.isOwn(row_node)) {

            auto begin = in_row_index[row_node];
            auto end = (row_node == matrix_nb_row - 1) ? matrix_nb_column : in_row_index[row_node + 1];

            while (begin < end) {
              if (in_columns[begin] == col_node_lid) {
                auto block_start = begin * NB_DOF_SQ;

                for (auto i = 0; i < NB_DOF; ++i) {
                  for (auto j = 0; j < NB_DOF; ++j) {
                    double value = element_matrix(NB_DOF * cur_row_node_idx + i, NB_DOF * cur_col_node_idx + j);
                    inout_values[block_start + (i * NB_DOF + j)] += value;
                  }
                }
                break;
              }
              ++begin;
            }
          }
          ++cur_col_node_idx;
        }
      }
    };
  };

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  template <class Function>
  void assembleBilinearOrderedPerRowAtomicFree(Function compute_element_matrix)
  {
    UnstructuredMeshConnectivityView m_connectivity_view(m_mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();
    auto node_cell_cv = m_connectivity_view.nodeCell();

    ItemGenericInfoListView nodes_infos(m_mesh->nodeFamily());

    constexpr int NB_DOF_SQ = NB_DOF * NB_DOF;
    auto matrix_nb_row = m_bsr_matrix.nbRow();
    auto matrix_nb_column = m_bsr_matrix.nbCol();
    auto matrix_nb_nz = m_bsr_matrix.nbNz();

    auto command = makeCommand(m_queue);
    auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
    auto in_columns = viewIn(command, m_bsr_matrix.columns());
    auto inout_values = viewInOut(command, m_bsr_matrix.values());
    auto in_nz_per_row = viewIn(command, m_bsr_matrix.nbNzPerRow());

    command << RUNCOMMAND_ENUMERATE(Node, row_node, m_mesh->allNodes())
    {
      auto cur_row_node_idx = 0;
      for (auto cell : node_cell_cv.cells(row_node)) {

        // Find the index of the node in the current cell
        for (Int32 i = 1; i <= 3; ++i) {
          if (row_node == cell_node_cv.nodeId(cell, i)) {
            cur_row_node_idx = i;
            break;
          }
          cur_row_node_idx = 0;
        }

        auto element_matrix = compute_element_matrix(cell);

        auto cur_col_node_idx = 0;
        for (NodeLocalId col_node_lid : cell_node_cv.nodes(cell)) {
          if (nodes_infos.isOwn(row_node)) {

            auto g_block_start = in_row_index[row_node] * NB_DOF_SQ;
            auto begin = in_row_index[row_node];
            auto end = (row_node == matrix_nb_row - 1) ? matrix_nb_column : in_row_index[row_node + 1];

            auto x = 0;
            while (begin < end) {
              if (in_columns[begin] == col_node_lid) {
                for (auto i = 0; i < NB_DOF; ++i) {
                  for (auto j = 0; j < NB_DOF; ++j) {
                    double value = element_matrix(NB_DOF * cur_row_node_idx + i, NB_DOF * cur_col_node_idx + j);
                    auto l_block_start = g_block_start + (NB_DOF * x);
                    if (i != 0)
                      l_block_start += (NB_DOF * in_nz_per_row[row_node]);
                    inout_values[l_block_start + j] += value;
                  }
                }
              }
              ++begin;
              ++x;
            }
            ++cur_col_node_idx;
          }
        }
      };
    };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Assembles the global BSR matrix for a bilinear operator, node-wise.
   *
   * This function constructs the Block Sparse Row (BSR) matrix by iterating over mesh nodes
   * and accumulating contributions from element-level matrices. It uses a user-defined
   * callback to compute the local matrix for each cell and updates the global matrix.
   *
   * ### Parameters:
   * - `compute_element_matrix`: A callback function that computes the local matrix
   *   for a given cell. It must return a matrix compatible with the expected block size.
   *
   * ### Key Details:
   * - Uses cell-to-node, node-to-cell connectivity from the mesh and DoF mappings for assembly.
   * - Handles block structure updates for `blockSize x blockSize` dimensions.
   *
   */
  /*---------------------------------------------------------------------------*/

  template <class Function> void assembleBilinearAtomicFree(Function compute_element_matrix)
  {
    info() << "BSRFormat(assembleBilinearAtomicFree): Integrating over nodes...";
    auto startTime = platform::getRealTime();

    if (m_bsr_matrix.orderValuePerBlock())
      assembleBilinearOrderedPerBlockAtomicFree(compute_element_matrix);
    else
      assembleBilinearOrderedPerRowAtomicFree(compute_element_matrix);

    info() << "[ArcaneFem-Timer] Time to assemble (atomic-free implementation) BSR matrix = " << (platform::getRealTime() - startTime);
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  template <class Function> void assembleBilinear(Function compute_element_matrix)
  {
    if (m_use_atomic_free)
      assembleBilinearAtomicFree(compute_element_matrix);
    else
      assembleBilinearAtomic(compute_element_matrix);
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  BSRMatrix<NB_DOF>& matrix()
  {
    return m_bsr_matrix;
  };
  void resetMatrixValues()
  {
    m_bsr_matrix.values().fill(0, m_queue);
  };
  void dumpMatrix(std::string filename) const
  {
    m_bsr_matrix.dump(filename);
  };

 private:

  bool m_use_csr_in_linear_system = false;
  bool m_use_atomic_free = false;

  BSRMatrix<NB_DOF> m_bsr_matrix;
  CsrFormat m_csr_matrix;

  IMesh* m_mesh;
  RunQueue& m_queue;
  const FemDoFsOnNodes& m_dofs_on_nodes;
};
}; // namespace Arcane::FemUtils

#endif // ! BSRFORMAT_H
