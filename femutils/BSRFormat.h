// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BSRFormat.h                                                 (C) 2022-2024 */
/*                                                                           */
/* BSRFormat class definition.                                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef BSRFORMAT_H
#define BSRFORMAT_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arccore/trace/TraceAccessor.h>

#include <arccore/base/NotImplementedException.h>
#include <arccore/base/ArgumentException.h>
#include <arccore/base/ArccoreGlobal.h>

#include <arcane/core/IIndexedIncrementalItemConnectivityMng.h>
#include <arcane/core/IIndexedIncrementalItemConnectivity.h>
#include <arcane/core/UnstructuredMeshConnectivity.h>
#include <arcane/core/IIncrementalItemConnectivity.h>
#include <arcane/core/IndexedItemConnectivityView.h>
#include <arcane/core/ItemEnumerator.h>
#include <arcane/core/ItemTypes.h>
#include <arcane/core/MeshUtils.h>
#include <arcane/core/DataView.h>
#include <arcane/core/IMesh.h>

#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/ArrayLayout.h>
#include <arcane/utils/UtilsTypes.h>
#include <arcane/utils/NumArray.h>
#include <arcane/utils/MDDim.h>

#include <arcane/accelerator/RunCommandEnumerate.h>
#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/accelerator/NumArrayViews.h>
#include <arcane/accelerator/Atomic.h>
#include <arcane/accelerator/Scan.h>

#include "DoFLinearSystem.h"
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
 *       - Values:    [ 1, 2, 3, 4, 5, 6, 9, 10, 7, 8, 11, 12, 13, 14, 15, 16 ]
 *       - Columns:   [ 0, 1, 2, 2 ]  // Column indices of each block
 *       - Row Index: [ 0, 1, 3 ]     // Starting indices of blocks per row
 *
 * @note The class supports flexible memory management for computations on 
 * different hardware backends (e.g., CPU, GPU).
 */
/*---------------------------------------------------------------------------*/

class BSRMatrix : public TraceAccessor
{
 public:

  BSRMatrix(ITraceMng* tm, const eMemoryRessource& mem_ressource)
  : TraceAccessor(tm)
  , m_values(mem_ressource)
  , m_columns(mem_ressource)
  , m_row_index(mem_ressource) {};

  void initialize(Int32 block_size, Int32 nb_non_zero_value, Int32 nb_col, Int32 nb_row, const RunQueue& queue);

  Int32 blockSize() { return m_block_size; };
  Int32 nbNz() { return m_nb_non_zero_value; };
  Int32 nbCol() { return m_nb_col; };
  Int32 nbRow() { return m_nb_row; };

  NumArray<Real, MDDim1>& values() { return m_values; }
  NumArray<Int32, MDDim1>& columns() { return m_columns; }
  NumArray<Int32, MDDim1>& rowIndex() { return m_row_index; }

  void toLinearSystem(DoFLinearSystem& linear_system); // TODO: Make it use GPU ?
  void dump(std::string filename);

 private:

  Int32 m_block_size;
  Int32 m_nb_non_zero_value;
  Int32 m_nb_col;
  Int32 m_nb_row;

  NumArray<Real, MDDim1> m_values;
  NumArray<Int32, MDDim1> m_columns;
  NumArray<Int32, MDDim1> m_row_index;
};

/*---------------------------------------------------------------------------*/
/**index
 * @brief A class for assembling Block Sparse Row (BSR) matrices from mesh data.
 *
 * The `BSRFormat` class is designed to construct and manage BSR matrices using 
 * mesh connectivity and degree-of-freedom (DoF) information. It supports operations 
 * such as sparsity pattern computation, exporting matrices to linear systems, and 
 * assembling bilinear operators for finite element applications.
 *
 * In 2D, the sparsity is computed based on node-face connectivity, while in 3D, 
 * node-node connectivity is used. Below is an example workflow:
 * - **2D Connectivity**:
 *   - Row indices are determined by node-to-face relationships.
 *   - Column indices are computed based on face adjacency.
 * - **3D Connectivity**:
 *   - Both row and column indices are computed from node-node relationships.
 *
 * @note This class uses GPU-accelerated compute kernels when applicable via `RunQueue`.
 * It relies on a pre-initialized `BSRMatrix` for matrix storage  and operations.
 */
/*---------------------------------------------------------------------------*/

class BSRFormat : public TraceAccessor
{
 public:

  BSRFormat(ITraceMng* tm, RunQueue& queue, IMesh& mesh, const FemDoFsOnNodes& dofs_on_nodes)
  : TraceAccessor(tm)
  , m_queue(queue)
  , m_mesh(mesh)
  , m_dofs_on_nodes(dofs_on_nodes)
  , m_bsr_matrix(tm, queue.memoryRessource())
  {
    if (m_mesh.dimension() != 2 && m_mesh.dimension() != 3)
      ARCANE_THROW(NotImplementedException, "BSRFormat(Ctor): Only supports 2D and 3D");
  };

  void initialize(Int32 nb_edge); // NOTE: Could compute nb_edge inside function via m_mesh
  void toLinearSystem(DoFLinearSystem& linear_system) { m_bsr_matrix.toLinearSystem(linear_system); };
  void dumpMatrix(std::string filename) { m_bsr_matrix.dump(filename); };
  void resetMatrixValues() { m_bsr_matrix.values().fill(0, &m_queue); };

  void computeSparsity();

  void computeSparsityRowIndex2D(Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> copy_out_data);
  void computeSparsityRowIndex3D(Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> copy_out_data);
  void computeSparsityRowIndex();

  void computeSparsityColumns2D(Accelerator::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> in_row_index, Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> inout_columns);
  void computeSparsityColumns3D(Accelerator::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> in_row_index, Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> inout_columns);
  void computeSparsityColumns();

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

  template <class Function> void assembleCellWise(Function compute_element_matrix)
  {
    info() << "BSRFormat(assembleCellWise): Assemble bilinear operator cell-wise";

    UnstructuredMeshConnectivityView m_connectivity_view(&m_mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();

    ItemGenericInfoListView nodes_infos(m_mesh.nodeFamily());

    auto block_size = m_bsr_matrix.blockSize();
    auto matrix_nb_row = m_bsr_matrix.nbRow();
    auto matrix_nb_column = m_bsr_matrix.nbCol();
    auto matrix_nb_nz = m_bsr_matrix.nbNz();

    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    auto command = makeCommand(m_queue);
    auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
    auto in_columns = viewIn(command, m_bsr_matrix.columns());
    auto inout_values = viewInOut(command, m_bsr_matrix.values());

    command << RUNCOMMAND_ENUMERATE(Cell, cell, m_mesh.allCells())
    {
      auto element_matrix = compute_element_matrix(cell);

      auto cur_row_node_idx = 0;
      for (NodeLocalId row_node_lid : cell_node_cv.nodes(cell)) {
        auto cur_col_node_idx = 0;
        for (NodeLocalId col_node_lid : cell_node_cv.nodes(cell)) {
          if (nodes_infos.isOwn(row_node_lid)) {

            auto begin = in_row_index[row_node_lid];
            auto end = (row_node_lid == matrix_nb_row - 1) ? matrix_nb_column : in_row_index[row_node_lid + 1];

            // Find the column index
            while (begin < end) {
              if (in_columns[begin] == col_node_lid) {
                auto block_start = begin * (block_size * block_size);

                // Update the block values
                for (auto i = 0; i < block_size; ++i) {
                  for (auto j = 0; j < block_size; ++j) {
                    double value = element_matrix(block_size * cur_row_node_idx + i, block_size * cur_col_node_idx + j);
                    ARCANE_ASSERT((block_start + (i * block_size) + j) < matrix_nb_nz, ("BSRFormat(assembleCellWise): Index out of bounds in inout_values"));
                    Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(inout_values[block_start + (i * block_size + j)], value);
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

  template <class Function> void assembleNodeWise(Function compute_element_matrix)
  {
    info() << "BSRFormat(assembleNodeWise): Assemble bilinear operator node-wise";

    UnstructuredMeshConnectivityView m_connectivity_view(&m_mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();
    auto node_cell_cv = m_connectivity_view.nodeCell();

    ItemGenericInfoListView nodes_infos(m_mesh.nodeFamily());

    auto block_size = m_bsr_matrix.blockSize();
    auto matrix_nb_row = m_bsr_matrix.nbRow();
    auto matrix_nb_column = m_bsr_matrix.nbCol();
    auto matrix_nb_nz = m_bsr_matrix.nbNz();

    auto command = makeCommand(m_queue);
    auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
    auto in_columns = viewIn(command, m_bsr_matrix.columns());
    auto inout_values = viewInOut(command, m_bsr_matrix.values());

    command << RUNCOMMAND_ENUMERATE(Node, row_node, m_mesh.allNodes())
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

            // Find the column index
            while (begin < end) {
              if (in_columns[begin] == col_node_lid) {
                auto block_start = begin * (block_size * block_size);

                // Update the block values
                for (auto i = 0; i < block_size; ++i) {
                  for (auto j = 0; j < block_size; ++j) {
                    double value = element_matrix(block_size * cur_row_node_idx + i, block_size * cur_col_node_idx + j);
                    ARCANE_ASSERT((block_start + (i * block_size) + j) < matrix_nb_nz, ("BSRFormat(assembleNodeWise): Index out of bounds in inout_values"));
                    inout_values[block_start + (i * block_size + j)] += value;
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
  }

 private:

  BSRMatrix m_bsr_matrix;
  RunQueue& m_queue;
  IMesh& m_mesh;
  const FemDoFsOnNodes& m_dofs_on_nodes;
};
}; // namespace Arcane::FemUtils

#endif // ! BSRFORMAT_H
