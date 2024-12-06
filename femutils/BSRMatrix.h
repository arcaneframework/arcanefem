#include "FemDoFsOnNodes.h"
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

namespace Arcane::FemUtils
{

class BSRMatrix : public TraceAccessor
{
 public:

  BSRMatrix(ITraceMng* tm, const eMemoryRessource& mem_ressource);

  void initialize(Int8 block_size, Int32 nb_non_zero_value, Int32 nb_row, const RunQueue& queue);
  Int32 blockSize() { return m_block_size; };
  Int32 nbRow() { return m_row_index.extent0(); };
  Int32 nbNz() { return m_values.extent0(); };

  NumArray<Int32, MDDim1>& columns() { return m_columns; }
  NumArray<Int32, MDDim1>& rowIndex() { return m_row_index; }

  void dump(std::string filename);

 private:

  Int8 m_block_size; // TODO: What size for this int ?
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
    if (m_mesh.dimension() != 3) // TODO: Why dimension can't be called on const ?
      ARCANE_THROW(NotImplementedException, "BSRFormat(Ctor): Only supports 3D");
  };

  void initialize();
  void computeSparsity(const IndexedNodeNodeConnectivityView& node_node_cv);

  // TODO: I would like to make those 2 methods private
  void computeSparsityRowIndex(const IndexedNodeNodeConnectivityView& node_node_cv);
  void computeSparsityColumns(const IndexedNodeNodeConnectivityView& node_node_cv);

 private:

  RunQueue& m_queue;
  IMesh& m_mesh;
  const FemDoFsOnNodes& m_dofs_on_nodes;
  BSRMatrix m_bsr_matrix;
};
}; // namespace Arcane::FemUtils
