#include "DoFLinearSystem.h"
#include "FemUtils.h"
#include "FemDoFsOnNodes.h"
#include "arcane/core/DataView.h"
#include "arcane/core/IndexedItemConnectivityView.h"
#include "arcane/utils/ArrayLayout.h"
#include "arcane/utils/MDDim.h"
#include "arccore/base/ArccoreGlobal.h"
#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/core/IMesh.h>
#include <arcane/core/ItemEnumerator.h>
#include <arcane/core/ItemTypes.h>
#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/UtilsTypes.h>
#include <arccore/base/ArgumentException.h>
#include <arccore/base/NotImplementedException.h>
#include <arccore/trace/TraceAccessor.h>
#include <arcane/utils/NumArray.h>
#include <arcane/accelerator/NumArrayViews.h>
#include <arcane/accelerator/RunCommandEnumerate.h>
#include <arcane/accelerator/Scan.h>
#include <arcane/core/UnstructuredMeshConnectivity.h>
#include <arcane/accelerator/Atomic.h>
#include <arcane/core/IIndexedIncrementalItemConnectivity.h>
#include <arcane/core/IIndexedIncrementalItemConnectivityMng.h>
#include <arcane/core/IIncrementalItemConnectivity.h>
#include <arcane/core/IndexedItemConnectivityView.h>
#include <arcane/core/MeshUtils.h>

namespace Arcane::FemUtils
{

class BSRMatrix : public TraceAccessor
{
 public:

  BSRMatrix(ITraceMng* tm, const eMemoryRessource& mem_ressource);

  void initialize(Int32 block_size, Int32 nb_non_zero_value, Int32 nb_col, Int32 nb_row, const RunQueue& queue);

  Int32 blockSize() { return m_block_size; };
  Int32 nbNz() { return m_nb_non_zero_value; };
  Int32 nbCol() { return m_nb_col; };
  Int32 nbRow() { return m_nb_row; };

  NumArray<Real, MDDim1>& values() { return m_values; }
  NumArray<Int32, MDDim1>& columns() { return m_columns; }
  NumArray<Int32, MDDim1>& rowIndex() { return m_row_index; }

  void toLinearSystem(DoFLinearSystem& linear_system);
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
    if (m_mesh.dimension() != 2 && m_mesh.dimension() != 3) // TODO: Why dimension can't be called on a const mesh ?
      ARCANE_THROW(NotImplementedException, "BSRFormat(Ctor): Only supports 2D and 3D");
  };

  void initialize(Int32 nb_edge); // TODO: Un peu dommage de devoir passer un argument...
  // Could remove the argument by passing via an UnstructuedMeshConnectivityMessh initialized with m_mesh
  void computeSparsity();

  void computeSparsityRowIndex2D(Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> copy_out_data);
  void computeSparsityRowIndex3D(Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> copy_out_data);
  void computeSparsityRowIndex();

  void computeSparsityColumns2D(Accelerator::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> in_row_index, Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> inout_columns);
  void computeSparsityColumns3D(Accelerator::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> in_row_index, Accelerator::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> inout_columns);
  void computeSparsityColumns();

  // TODO: Be able to call the .dump() method of bsr_matrix from here ?
  // TODO: Be able to access bsr_matrix with a getter ?

  // TODO: try to make it less than 40 loc
  template <class Function> void assembleBilinear(Function compute_element_matrix)
  {
    info() << "BSRFormat(assembleBilinear): Assemble bilinear operator";
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
            Int32 begin = in_row_index[row_node_lid];
            auto end = (row_node_lid == matrix_nb_row - 1) ? matrix_nb_column : in_row_index[row_node_lid + 1];
            while (begin < end) {
              if (in_columns[begin] == col_node_lid) {
                auto block_start = begin * (block_size * block_size);
                for (Int32 i = 0; i < block_size; ++i) {
                  for (Int32 j = 0; j < block_size; ++j) {
                    double value = element_matrix(block_size * cur_row_node_idx + i, block_size * cur_col_node_idx + j);
                    ARCANE_ASSERT((block_start + (i * block_size) + j) < matrix_nb_nz, ("Index out of bounds in inout_values"));
                    if ((block_start + (i * block_size) + j) < matrix_nb_nz)
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

  // TODO: make it private
  BSRMatrix m_bsr_matrix;

 private:

  RunQueue& m_queue;
  IMesh& m_mesh;
  const FemDoFsOnNodes& m_dofs_on_nodes;
};
}; // namespace Arcane::FemUtils
