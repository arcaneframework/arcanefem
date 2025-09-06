// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BSRFormat.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Matrix format using Block Sparse Row.                                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "BSRFormat.h"

// #include <arccore/trace/TraceAccessor.h>

#include <arcane/utils/ArgumentException.h>

#include <arcane/core/IIndexedIncrementalItemConnectivityMng.h>
#include <arcane/core/IIndexedIncrementalItemConnectivity.h>
#include <arcane/core/MeshUtils.h>

#include <arcane/accelerator/RunCommandLoop.h>
#include <arcane/accelerator/GenericSorter.h>
#include <arcane/accelerator/Scan.h>

#include "DoFLinearSystem.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

BSRMatrix::
BSRMatrix(ITraceMng* tm, const eMemoryRessource& mem_ressource, const RunQueue& queue)
: TraceAccessor(tm)
, m_values(mem_ressource)
, m_columns(mem_ressource)
, m_row_index(mem_ressource)
, m_queue(queue) {};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRMatrix::
initialize(Int32 nb_non_zero_value, Int32 nb_col, Int32 nb_row, Int8 nb_block, bool order_values_per_block)
{
  if (nb_block <= 0 || nb_non_zero_value <= 0 || nb_row <= 0)
    ARCANE_THROW(ArgumentException, "BSRMatrix(initialize): arguments should be positive and not null (block_size={0}, nb_non_zero_value={1} and nb_row={2})", nb_block, nb_non_zero_value, nb_row);

  if (nb_block > nb_row)
    ARCANE_THROW(ArgumentException, "BSRMatrix(initialize): block_size should be less than nb_row");

  info() << "BSRMatrix(initialize): Initialize BSRMatrix with block_size=" << nb_block
         << ", nb_non_zero_value=" << nb_non_zero_value
         << ", nb_col=" << nb_col << ", nb_row=" << nb_row
         << ", order_values_per_block=" << std::boolalpha << order_values_per_block;

  m_order_values_per_block = order_values_per_block;

  m_nb_non_zero_value = nb_non_zero_value;
  m_nb_col = nb_col;
  m_nb_row = nb_row;
  m_nb_block = nb_block;

  m_values.resize(nb_non_zero_value);
  m_values.fill(0, &m_queue);
  m_columns.resize(nb_col);
  m_row_index.resize(nb_row);
  m_nb_nz_per_row.resize(nb_row);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Int32 BSRMatrix::
findValueIndex(DoFLocalId row, DoFLocalId col) const
{
  auto block_row = row / m_nb_block;
  auto block_col = col / m_nb_block;

  auto block_start = m_row_index[block_row];
  auto block_end = (block_row == m_nb_row - 1) ? m_nb_col : m_row_index[block_row + 1];

  auto row_offset = row % m_nb_block;
  auto col_offset = col % m_nb_block;

  auto nb_block_sq = m_nb_block * m_nb_block;
  auto block_start_in_value = block_start * nb_block_sq;
  auto col_index = 0;
  while (block_start < block_end) {
    if (m_columns[block_start] == block_col) {
      if (!m_order_values_per_block)
        return block_start_in_value + (m_nb_block * col_index) + (row_offset * m_nb_block * m_nb_nz_per_row[block_row]) + col_offset;
      else
        return (block_start * nb_block_sq) + ((row_offset * m_nb_block) + col_offset);
    }
    ++block_start;
    ++col_index;
  }

  ARCANE_FATAL("BSRMatrix(findValueIndex): Value not found");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRMatrix::
toCsr(CsrFormat* csr_matrix)
{
  info() << "BSRMatrix(toCsr): Convert matrix to CSR";
  auto startTime = platform::getRealTime();

  auto nb_block_rows = m_row_index.extent0();
  auto nb_rows = nb_block_rows * m_nb_block;
  auto total_non_zero_elements = m_nb_non_zero_value;

  csr_matrix->initialize(nullptr, total_non_zero_elements, nb_rows, m_queue);

  if (m_nb_block == 1) {
    csr_matrix->m_matrix_row = m_row_index;
    csr_matrix->m_matrix_column = m_columns;
    csr_matrix->m_matrix_rows_nb_column = m_nb_nz_per_row;
  }
  else {
    // Translate `row_index`
    csr_matrix->m_matrix_row[0] = 0;
    auto offset = 1;
    for (auto i = 0; i < m_row_index.extent0() - 1; ++i) {
      for (auto j = 0; j < m_nb_block; ++j) {
        auto start = m_row_index[i];
        auto end = i == m_row_index.extent0() - 1 ? m_nb_col : m_row_index[i + 1];
        auto nombre_de_dof_dans_la_rangee = (end - start) * m_nb_block;
        csr_matrix->m_matrix_row[offset] = csr_matrix->m_matrix_row[offset - 1] + nombre_de_dof_dans_la_rangee;
        offset++;
      }
    }

    // Translate `columns`
    offset = 0;
    for (auto i = 0; i < m_row_index.extent0(); ++i) {
      for (auto j = 0; j < m_nb_block; ++j) {
        auto start = m_row_index[i];
        auto end = i == m_row_index.extent0() - 1 ? m_nb_col : m_row_index[i + 1];
        for (auto block_index = start; block_index < end; ++block_index) {
          for (auto k = 0; k < m_nb_block; ++k) {
            csr_matrix->m_matrix_column[offset] = m_columns[block_index] * m_nb_block + k;
            offset++;
          }
        }
      }
    }

    // Translate `nb_nz_per_row`
    offset = 0;
    for (auto i = 0; i < m_nb_nz_per_row.extent0(); ++i) {
      for (auto j = 0; j < m_nb_block; ++j) {
        csr_matrix->m_matrix_rows_nb_column[offset++] = m_nb_nz_per_row[i] * m_nb_block;
      }
    }

    // Fill last entry of CSR matrix row_index array
    auto row_index_size = csr_matrix->m_matrix_row.extent0();
    csr_matrix->m_matrix_row[row_index_size - 1] = csr_matrix->m_matrix_row[row_index_size - 2] + csr_matrix->m_matrix_rows_nb_column[csr_matrix->m_matrix_rows_nb_column.extent0() - 1];
  }

  // NOTE: If we don't want to keep bsr matrix coefficients we could move the data instead of copying it.
  csr_matrix->m_matrix_value = m_values;
  info() << std::left << std::setw(40) << "[BsrMatrix-Timer] convert-bsr-to-csr" << " = " << (platform::getRealTime() - startTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRMatrix::
toLinearSystem(DoFLinearSystem& linear_system)
{
  info() << "BSRMatrix(toLinearSystem): Translate matrix to linear system using `matrixAddValue`";

  for (auto row = 0; row < m_nb_row; ++row) {
    auto row_start = m_row_index[row];
    auto row_end = (row + 1 < m_nb_row) ? m_row_index[row + 1] : m_nb_col;
    for (auto block_idx = row_start; block_idx < row_end; ++block_idx) {
      auto col = m_columns[block_idx];
      for (auto i = 0; i < m_nb_block; ++i) {
        for (auto j = 0; j < m_nb_block; ++j) {
          auto global_row = (row * m_nb_block) + i;
          auto global_col = (col * m_nb_block) + j;
          auto value = getValue(DoFLocalId(global_row), DoFLocalId(global_col));
          linear_system.matrixAddValue(DoFLocalId(global_row), DoFLocalId(global_col), value);
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRMatrix::
dump(const String& filename)
{
  info() << "BSRMatrix(dump): Dump matrix in '" << filename << "'";
  ofstream file(filename.localstr());

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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

BSRFormat::
BSRFormat(ITraceMng* tm, RunQueue& queue, const FemDoFsOnNodes& dofs_on_nodes)
: TraceAccessor(tm)
, m_dofs_on_nodes(dofs_on_nodes)
, m_queue(queue)
, m_bsr_matrix(tm, queue.memoryRessource(), queue)
, m_csr_matrix(tm)
{}

Int64 BSRFormat::
computeNbEdge(IMesh* mesh)
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::
initialize(IMesh* mesh, Int8 nb_dof, bool does_linear_system_use_csr, bool use_atomic_free)
{
  ARCANE_CHECK_POINTER(mesh);

  if (mesh->dimension() != 2 && mesh->dimension() != 3)
    ARCANE_THROW(NotImplementedException, "BSRFormat(initialize): Only supports 2D and 3D");

  auto startTime = platform::getRealTime();

  m_mesh = mesh;
  m_nb_dof = nb_dof;
  Int64 nb_edge = computeNbEdge(mesh);
  Int32 nb_node = m_mesh->nbNode();
  Int32 nb_col = 2 * nb_edge + nb_node;
  Int32 nb_non_zero_value = (m_nb_dof * m_nb_dof) * (2 * nb_edge + nb_node);

  m_use_csr_in_linear_system = does_linear_system_use_csr;
  bool order_values_per_block = !does_linear_system_use_csr;
  m_bsr_matrix.initialize(nb_non_zero_value, nb_col, nb_node, m_nb_dof, order_values_per_block);
  m_use_atomic_free = use_atomic_free;

  info() << std::left << std::setw(40) << "[BsrMatrix-Timer] initialize-bsr-matrix" << " = " << (platform::getRealTime() - startTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::
toLinearSystem(DoFLinearSystem& linear_system)
{
  auto startTime = platform::getRealTime();
  if (m_use_csr_in_linear_system) {
    if (!linear_system.hasSetCSRValues())
      ARCANE_THROW(ArgumentException, "BSRFormat(toLinearSystem): Linear system was set to use CSR but is incompatible");

    m_bsr_matrix.toCsr(&m_csr_matrix);

    info() << "BSRFormat(toLinearSystem): Set CSR values into linear system";
    CSRFormatView csr_view(m_csr_matrix.view());
    linear_system.setCSRValues(csr_view);
  }
  else
    m_bsr_matrix.toLinearSystem(linear_system);

  info() << std::left << std::setw(40) << "[BsrMatrix-Timer] setup-linear-systm" << " = " << (platform::getRealTime() - startTime);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::
computeNzPerRowArray()
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

  info() << std::left << std::setw(40) << "[BsrMatrix-Timer] compute-nnz-bsr" << " = " << (platform::getRealTime() - startTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::
computeNeighborsAtomicFree(SmallSpan<Int32>& neighbors_ss)
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

void BSRFormat::
computeRowIndexAtomicFree()
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

void BSRFormat::
computeColumnsAtomicFree()
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

void BSRFormat::
computeSparsityAtomicFree()
{
  info() << "BSRFormat(computeSparsityAtomicFree): Compute sparsity of BSR matrix using Arcane connectivities (e.g node-face for 2D and node-node for 3D)";
  auto startTime = platform::getRealTime();

  computeRowIndexAtomicFree();
  computeColumnsAtomicFree();

  info() << std::left << std::setw(40) << "[BsrMatrix-Timer] build-sparsity-af_bsr" << " = " << (platform::getRealTime() - startTime);

  if (m_use_csr_in_linear_system)
    computeNzPerRowArray();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::
computeSortedEdges(Int8 edges_per_element, Int64 nb_edge_total, SmallSpan<UInt64>& sorted_edges_ss)
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

void BSRFormat::
computeNeighbors(Int8 edges_per_element, Int64 nb_edge_total, NumArray<Int32, MDDim1>& neighbors, SmallSpan<UInt64>& sorted_edges_ss)
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

void BSRFormat::
computeRowIndex(Int8 edges_per_element, Int64 nb_edge_total, SmallSpan<UInt64>& sorted_edges_ss)
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

void BSRFormat::
computeColumns(Int8 edges_per_element, Int64 nb_edge_total, SmallSpan<uint64_t>& sorted_edges_ss)
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

void BSRFormat::
computeSparsityAtomic()
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

  info() << std::left << std::setw(40) << "[BsrMatrix-Timer] build-sparsity-bsr" << " = " << (platform::getRealTime() - startTime);

  if (m_use_csr_in_linear_system)
    computeNzPerRowArray();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::
computeSparsity()
{
  if (m_use_atomic_free)
    computeSparsityAtomicFree();
  else
    computeSparsityAtomic();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

BSRMatrix& BSRFormat::
matrix()
{
  return m_bsr_matrix;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::
resetMatrixValues()
{
  m_bsr_matrix.values().fill(0, m_queue);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::
dumpMatrix(const String& filename)
{
  m_bsr_matrix.dump(filename);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
