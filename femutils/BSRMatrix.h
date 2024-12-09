#include "FemUtils.h"
#include "FemDoFsOnNodes.h"
#include "arcane/accelerator/VariableViews.h"
#include "arcane/core/DataView.h"
#include "arcane/core/IndexedItemConnectivityView.h"
#include "arcane/utils/ArrayLayout.h"
#include "arcane/utils/MDDim.h"
#include "arccore/base/ArccoreGlobal.h"
#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/core/IMesh.h>
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
#include <variant>

namespace Arcane::FemUtils
{

class BSRMatrix : public TraceAccessor
{
 public:

  BSRMatrix(ITraceMng* tm, const eMemoryRessource& mem_ressource);

  void initialize(Int32 block_size, Int32 nb_non_zero_value, Int32 nb_row, const RunQueue& queue);
  Int32 blockSize() { return m_block_size; };
  Int32 nbRow() { return m_row_index.extent0(); };
  Int32 nbNz() { return m_values.extent0(); };

  NumArray<Real, MDDim1>& values() { return m_values; }
  NumArray<Int32, MDDim1>& columns() { return m_columns; }
  NumArray<Int32, MDDim1>& rowIndex() { return m_row_index; }

  void dump(std::string filename);

 private:

  Int32 m_block_size; // TODO: What size for this int ?
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
  template <int N, class Function>
  void assembleBilinear(Function compute_element_matrix)
  {
    info() << "BSRFormat(assembleBilinear): Assemble bilinear operator";
    UnstructuredMeshConnectivityView m_connectivity_view(&m_mesh);
    auto cell_node_cv = m_connectivity_view.cellNode();

    ItemGenericInfoListView nodes_infos(m_mesh.nodeFamily());

    Int32 matrix_nb_row = m_bsr_matrix.nbRow();
    Int32 matrix_nb_column = m_bsr_matrix.nbNz();

    auto command = makeCommand(m_queue);
    auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
    auto in_columns = viewIn(command, m_bsr_matrix.columns());
    auto inout_values = viewInOut(command, m_bsr_matrix.values());

    command << RUNCOMMAND_ENUMERATE(Cell, cell, m_mesh.allCells())
    {
      auto element_matrix = compute_element_matrix(cell);

      Int32 cur_src_node_index = 0; // index on cell !
      for (NodeLocalId src_node_lid : cell_node_cv.nodes(cell)) { // TODO: cell should be of type ItemLocalId and not Cell ??
        Int32 cur_dst_node_index = 0;
        for (NodeLocalId dst_node_lid : cell_node_cv.nodes(cell)) {
          if (nodes_infos.isOwn(src_node_lid)) {
            // TODO: Existe-t-il une conversion implicite entre NodeLocalId et Int32 (pour ne pas faire .localId()) ?
            double value = element_matrix(cur_src_node_index, cur_dst_node_index);

            Int32 row_in_matrix = src_node_lid.localId();
            Int32 col_in_matrix = dst_node_lid.localId();
            Int32 begin = in_row_index[row_in_matrix];
            Int32 end = (row_in_matrix == matrix_nb_row - 1) ? matrix_nb_column : in_row_index[row_in_matrix + 1];

            while (begin < end) {
              if (in_columns[begin] == col_in_matrix) {
                Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(inout_values(begin), value);
                break;
              }
              ++begin;
            }
          }
          ++cur_dst_node_index;
        }
        ++cur_src_node_index;
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
