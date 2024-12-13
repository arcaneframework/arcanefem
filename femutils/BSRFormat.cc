// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BSRFormat.cc                                                (C) 2022-2024 */
/*                                                                           */
/* BSRFormat class implementation.                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "BSRFormat.h"
#include <arccore/base/ArccoreGlobal.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

void BSRMatrix::initialize(Int32 block_size, Int32 nb_non_zero_value, Int32 nb_col, Int32 nb_row)
{
  if (block_size <= 0 || nb_non_zero_value <= 0 || nb_row <= 0)
    ARCANE_THROW(ArgumentException, "BSRMatrix(initialize): arguments should be positive and not null (block_size={0}, nb_non_zero_value={1} and nb_row={2})", block_size, nb_non_zero_value, nb_row);

  if (block_size > nb_row)
    ARCANE_THROW(ArgumentException, "BSRMatrix(initialize): block_size should be less than nb_row");

  info() << "BSRMatrix(initialize): Initialize BSRMatrix with block_size=" << block_size << ", nb_non_zero_value=" << nb_non_zero_value << ", nb_row=" << nb_row;

  m_block_size = block_size;
  m_nb_non_zero_value = nb_non_zero_value;
  m_nb_col = nb_col;
  m_nb_row = nb_row;

  m_values.resize(nb_non_zero_value);
  m_values.fill(0, &m_queue);
  m_columns.resize(nb_col);
  m_row_index.resize(nb_row);
  m_nb_nz_per_row.resize(nb_row);
}

/**
 * Retrieves the value of a specific degree of freedom (DoF) from the BSR matrix.
 *
 * @param row The row index of the DoF.
 * @param col The column index of the DoF.
 * @param isValueArrayCsrLike A flag indicating if the value array uses a CSR-like layout.
 *                            - true: The value array is structured in a CSR-like manner.
 *                            - false: The value array is structured in a standard BSR layout.
 *
 * @return The value corresponding to the specified row and column indices.
 *
 * @throws std::runtime_error If the value is not found in the BSR matrix.
 */

Real BSRMatrix::getValue(DoFLocalId row, DoFLocalId col, bool isValueArrayCsrLike)
{
  auto block_row = row / 2;
  auto block_col = col / 2;

  auto block_start = m_row_index[block_row];
  auto block_end = (block_row == m_nb_row - 1) ? m_nb_col : m_row_index[block_row + 1];

  auto row_offset = row % 2;
  auto col_offset = col % 2;

  auto block_start_in_value = block_start * (m_block_size * m_block_size);
  auto col_index = 0;
  while (block_start < block_end) {
    if (m_columns[block_start] == block_col) {
      if (isValueArrayCsrLike) {
        auto value_index = block_start_in_value + (m_block_size * col_index) + ((row_offset) * (m_block_size * m_nb_nz_per_row[block_row]));
        return m_values[value_index + (col_offset)];
      }
      else {
        block_start = block_start * (m_block_size * m_block_size);
        auto value_index = block_start + ((row_offset)*m_block_size + (col_offset));
        return m_values[value_index];
      }
    }
    ++block_start;
    ++col_index;
  }

  throw std::runtime_error("BSRMatrix(getValue): Value not found");
}

/**
 * Sets the value of a specific degree of freedom (DoF) in the BSR matrix.
 *
 * @param row The row index of the DoF.
 * @param col The column index of the DoF.
 * @param value The value to set at the specified row and column indices.
 * @param isValueArrayCsrLike A flag indicating if the value array uses a CSR-like layout.
 *                            - true: The value array is structured in a CSR-like manner.
 *                            - false: The value array is structured in a standard BSR layout.
 */

void BSRMatrix::setValue(DoFLocalId row, DoFLocalId col, Real value, bool isValueArrayCsrLike)
{
  auto block_row = row / 2;
  auto block_col = col / 2;

  auto block_start = m_row_index[block_row];
  auto block_end = (block_row == m_nb_row - 1) ? m_nb_col : m_row_index[block_row + 1];

  auto row_offset = row % 2;
  auto col_offset = col % 2;

  auto block_start_in_value = block_start * (m_block_size * m_block_size);
  auto col_index = 0;
  while (block_start < block_end) {
    if (m_columns[block_start] == block_col) {
      if (isValueArrayCsrLike) {
        auto value_index = block_start_in_value + (m_block_size * col_index) + ((row_offset) * (m_block_size * m_nb_nz_per_row[block_row]));
        m_values[value_index + (col_offset)] = value;
        return;
      }
      else {
        block_start = block_start * (m_block_size * m_block_size);
        auto value_index = block_start + ((row_offset)*m_block_size + (col_offset));
        m_values[value_index] = value;
        return;
      }
    }
    ++block_start;
    ++col_index;
  }

  throw std::runtime_error("BSRMatrix(setValue): Failed to set value in BSR matrix");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRMatrix::toCsr(CsrFormat* csr_matrix)
{
  ARCCORE_CHECK_POINTER(csr_matrix);
  info() << "BSRMatrix(toCsr): Convert matrix to CSR";

   auto block_size = m_block_size;
   auto b_squared = block_size * block_size;

   auto nb_block_rows = m_row_index.extent0() - 1;
   auto nb_rows = nb_block_rows * block_size;
   auto total_non_zero_elements = m_columns.extent0() * b_squared;

  csr_matrix->initialize(nullptr, total_non_zero_elements, nb_rows, m_queue);

  csr_matrix->m_matrix_row.resize(nb_rows + 1); // Initialize all entries to 0
  csr_matrix->m_matrix_row.fill(0);
  csr_matrix->m_matrix_column.resize(total_non_zero_elements);

  // Compute csr_matrix.m_matrix_row
  for (auto block_row = 0; block_row < nb_block_rows; ++block_row) {
     auto start_block = m_row_index[block_row];
     auto end_block = m_row_index[block_row + 1];
     auto nb_blocks = end_block - start_block;

    for (auto offset = 0; offset < block_size; ++offset) {
       auto row = block_row * block_size + offset;
      csr_matrix->m_matrix_row[row + 1] = csr_matrix->m_matrix_row[row] + nb_blocks * block_size;
    }
  }

  // Compute csr_columns
  size_t csr_index = 0;
  for (auto block_row = 0; block_row < nb_block_rows; ++block_row) {
     auto start_block = m_row_index[block_row];
     auto end_block = block_row == nb_block_rows - 1 ? m_nb_col :  m_row_index[block_row + 1];

    for (auto row_offset = 0; row_offset < block_size; ++row_offset) {
      for (auto block_index = start_block; block_index < end_block; ++block_index) {
         auto block_col = m_columns[block_index];
        for (Int32 col_offset = 0; col_offset < block_size; ++col_offset) {
          assert(csr_index < total_non_zero_elements); // Ensure no out-of-bounds
          csr_matrix->m_matrix_column[csr_index++] = block_col * block_size + col_offset;
        }
      }
    }
  }

  // Copy bsr_values into csr_matrix.m_matrix_value
  csr_matrix->m_matrix_value.resize(m_values.extent0());
  for (auto i = 0; i < m_values.extent0(); ++i)
    csr_matrix->m_matrix_value[i] = m_values[i];

  // Copy bsr_nb_nz_per_row into csr_matrix
  csr_matrix->m_matrix_rows_nb_column.resize(m_nb_nz_per_row.extent0());
  for (auto i = 0; i < m_nb_nz_per_row.extent0(); ++i)
    csr_matrix->m_matrix_rows_nb_column[i] = m_nb_nz_per_row[i];
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRMatrix::toLinearSystem(DoFLinearSystem& linear_system, CsrFormat* csr_matrix)
{
  info() << "BSRMatrix(toLinearSystem): Translate matrix to linear system";

  if (linear_system.hasSetCSRValues()) {
    ARCCORE_CHECK_POINTER(csr_matrix);
    this->toCsr(csr_matrix);
    CSRFormatView csr_view(csr_matrix->m_matrix_row.to1DSpan(), csr_matrix->m_matrix_column.to1DSpan(), csr_matrix->m_matrix_rows_nb_column.to1DSpan(), m_values.to1DSpan());
    linear_system.setCSRValues(csr_view);
  }
  else {
    for (auto row = 0; row < m_nb_row; ++row) {
      auto row_start = m_row_index[row];
      auto row_end = (row + 1 < m_nb_row) ? m_row_index[row + 1] : m_nb_col;

      for (auto block_idx = row_start; block_idx < row_end; ++block_idx) {
        auto col = m_columns[block_idx];
        for (auto i = 0; i < m_block_size; ++i) {
          for (auto j = 0; j < m_block_size; ++j) {
            auto global_row = row * m_block_size + i;
            auto global_col = col * m_block_size + j;
            auto value = getValue(DoFLocalId(global_row), DoFLocalId(global_col), true);
            linear_system.matrixAddValue(DoFLocalId(global_row), DoFLocalId(global_col), value);
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRMatrix::dump(std::string filename)
{
  info() << "BSRMatrix(dump): Dump matrix in \"" << filename << "\"";
  ofstream file(filename);

  file << "size :" << nbNz() << "\n";
  for (auto i = 0; i < nbRow(); ++i) {
    file << m_row_index(i) << " ";
    for (Int32 j = m_row_index(i) + 1; (i + 1 < m_row_index.dim1Size() && j < m_row_index(i + 1)) || (i + 1 == m_row_index.dim1Size() && j < m_columns.dim1Size()); j++)
      file << "  ";
  }
  file << "\n";

  for (auto i = 0; i < nbNz(); ++i)
    file << m_columns(i) << " ";
  file << "\n";

  for (auto i = 0; i < nbNz(); ++i)
    file << m_values(i) << " ";
  file << "\n";

  file.close();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::initialize(Int32 nb_edge)
{
  Int32 nb_node = m_mesh.nbNode(); // Is the number of row
  Int32 nb_dof = m_dofs_on_nodes.nbDofPerNode();
  Int32 nb_col = 2 * nb_edge + nb_node;
  Int32 nb_non_zero_value = (nb_dof * nb_dof) * (2 * nb_edge + nb_node);

  m_bsr_matrix.initialize(nb_dof, nb_non_zero_value, nb_col, nb_node);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::computeSparsityRowIndex2D(Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> copy_out_data)
{
  auto command = makeCommand(m_queue);
  UnstructuredMeshConnectivityView connectivity_view(&m_mesh);
  auto node_face_cv = connectivity_view.nodeFace();
  command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, m_mesh.allNodes())
  {
    copy_out_data[node_lid.asInt32()] = node_face_cv.nbFace(node_lid) + 1;
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::computeSparsityRowIndex3D(Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> copy_out_data)
{
  auto command = makeCommand(m_queue);
  auto connectivity_mng = m_mesh.indexedConnectivityMng();
  auto connectivity_ptr = connectivity_mng->findOrCreateConnectivity(m_mesh.nodeFamily(), m_mesh.nodeFamily(), "NodeNodeViaEdge");
  IndexedNodeNodeConnectivityView node_node_cv = connectivity_ptr->view();
  command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, m_mesh.allNodes())
  {
    copy_out_data[node_lid.asInt32()] = node_node_cv.nbItem(node_lid) + 1;
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::computeSparsityRowIndex()
{
  info() << "BSRFormat(computeSparsityRowIndex): Compute row index sparsity of BSRMatrix";
  auto command = makeCommand(m_queue);

  NumArray<Int32, MDDim1> out_data;
  out_data.resize(m_bsr_matrix.nbRow());
  auto copy_out_data = viewInOut(command, out_data);

  m_mesh.dimension() == 2 ? computeSparsityRowIndex2D(copy_out_data) : computeSparsityRowIndex3D(copy_out_data);
  m_queue.barrier();

  Accelerator::Scanner<Int32> scanner;
  scanner.exclusiveSum(&m_queue, out_data, m_bsr_matrix.rowIndex());
  m_queue.barrier();

  // NOTE: Compute "m_matrix_rows_nb_column" array
  // Better to do it later, like in the translate to linear
  // system loop but needed during the assembly with BSR...
  auto nb_row = m_bsr_matrix.nbRow();
  {
    auto command = makeCommand(m_queue);

    m_bsr_matrix.nbNzPerRow().copy(m_bsr_matrix.rowIndex());
    auto inout_nb_nz_per_row = viewInOut(command, m_bsr_matrix.nbNzPerRow());

    auto nb_nz_per_row = m_bsr_matrix.nbNzPerRow();
    NumArray<Int32, MDDim1> nb_nz_per_row_copy = nb_nz_per_row;
    auto in_nb_nz_per_row_copy = viewIn(command, nb_nz_per_row_copy);

    command << RUNCOMMAND_LOOP1(iter, nb_row - 1)
    {
      auto [i] = iter();
      auto x = inout_nb_nz_per_row[i];
      inout_nb_nz_per_row[i] = in_nb_nz_per_row_copy[i + 1] - x;
    };

    m_bsr_matrix.nbNzPerRow()[nb_row - 1] = m_bsr_matrix.nbCol() - nb_nz_per_row[nb_row - 1];
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::computeSparsityColumns2D(Accelerator::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> in_row_index, Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> inout_columns)
{
  auto command = makeCommand(m_queue);
  UnstructuredMeshConnectivityView connectivity_view(&m_mesh);
  auto node_face_cv = connectivity_view.nodeFace();
  auto face_node_cv = connectivity_view.faceNode();
  command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, m_mesh.allNodes())
  {
    auto offset = in_row_index[node_lid.asInt32()];
    for (auto face_lid : node_face_cv.faceIds(node_lid)) {
      auto nodes = face_node_cv.nodes(face_lid);
      inout_columns[offset] = nodes[0] == node_lid ? nodes[1] : nodes[0];
      ++offset;
    }
    inout_columns[offset] = node_lid.asInt32();
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::computeSparsityColumns3D(Accelerator::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> in_row_index, Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> inout_columns)
{
  auto command = makeCommand(m_queue);
  auto connectivity_mng = m_mesh.indexedConnectivityMng();
  auto connectivity_ptr = connectivity_mng->findOrCreateConnectivity(m_mesh.nodeFamily(), m_mesh.nodeFamily(), "NodeNodeViaEdge");
  IndexedNodeNodeConnectivityView node_node_cv = connectivity_ptr->view();
  command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, m_mesh.allNodes())
  {
    auto offset = in_row_index[node_lid.asInt32()];
    for (auto neighbor_lid : node_node_cv.items(node_lid)) {
      inout_columns[offset] = neighbor_lid.asInt32();
      ++offset;
    }
    inout_columns[offset] = node_lid.asInt32();
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::computeSparsityColumns()
{
  info() << "BSRFormat(computeSparsityColumns): Compute columns sparsity of BSRMatrix";
  auto command = makeCommand(m_queue);
  auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
  auto inout_columns = viewInOut(command, m_bsr_matrix.columns());
  m_mesh.dimension() == 2 ? computeSparsityColumns2D(in_row_index, inout_columns) : computeSparsityColumns3D(in_row_index, inout_columns);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void BSRFormat::computeSparsity()
{
  info() << "BSRFormat(computeSparsity): Compute sparsity of BSRMatrix";
  computeSparsityRowIndex();
  computeSparsityColumns();
}

}; // namespace Arcane::FemUtils
