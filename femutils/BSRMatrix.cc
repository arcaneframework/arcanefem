#include "BSRMatrix.h"
#include "DoFLinearSystem.h"
#include <arcane/core/IndexedItemConnectivityView.h>

namespace Arcane::FemUtils
{

BSRMatrix::BSRMatrix(ITraceMng* tm, const eMemoryRessource& mem_ressource)
: TraceAccessor(tm)
, m_values(mem_ressource)
, m_columns(mem_ressource)
, m_row_index(mem_ressource) {};

void BSRMatrix::initialize(Int32 block_size, Int32 nb_non_zero_value, Int32 nb_col, Int32 nb_row, const RunQueue& queue)
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
  m_values.fill(0, &queue);
  m_columns.resize(nb_col);
  m_row_index.resize(nb_row);
}

void BSRMatrix::toLinearSystem(DoFLinearSystem& linear_system)
{
  info() << "BSRMatrix(toLinearSystem): Translate matrix to linear system";
  for (auto row = 0; row < m_nb_row; ++row) {
    auto row_start = m_row_index[row];
    auto row_end = (row + 1 < m_nb_row) ? m_row_index[row + 1] : m_nb_col;

    if (row_start == row_end)
      continue; // Skip empty rows

    for (auto block_idx = row_start; block_idx < row_end; ++block_idx) {
      auto col = m_columns[block_idx];
      auto block_start = block_idx * (m_block_size * m_block_size);

      for (auto i = 0; i < m_block_size; ++i) {
        for (auto j = 0; j < m_block_size; ++j) {
          auto value = m_values[block_start + i * m_block_size + j];
          auto global_row = row * m_block_size + i;
          auto global_col = col * m_block_size + j;
          linear_system.matrixAddValue(DoFLocalId(global_row), DoFLocalId(global_col), value);
        }
      }
    }
  }
}

void BSRMatrix::dump(std::string filename)
{
  info() << "BSRMatrix(dump): Dump BSRMatrix in \"" << filename << "\"";
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

void BSRFormat::initialize(Int32 nb_edge)
{
  Int32 nb_node = m_mesh.nbNode();
  Int32 nb_dof = m_dofs_on_nodes.nbDofPerNode();
  Int32 nb_col = 2 * nb_edge + nb_node;
  Int32 nb_non_zero_value = (nb_dof * nb_dof) * (2 * nb_edge + nb_node);
  Int32 nb_row = nb_node;

  m_bsr_matrix.initialize(nb_dof, nb_non_zero_value, nb_col, nb_row, m_queue);
}

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
}

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

void BSRFormat::computeSparsityColumns()
{
  info() << "BSRFormat(computeSparsityColumns): Compute columns sparsity of BSRMatrix";
  auto command = makeCommand(m_queue);
  auto in_row_index = viewIn(command, m_bsr_matrix.rowIndex());
  auto inout_columns = viewInOut(command, m_bsr_matrix.columns());
  m_mesh.dimension() == 2 ? computeSparsityColumns2D(in_row_index, inout_columns) : computeSparsityColumns3D(in_row_index, inout_columns);
}

void BSRFormat::computeSparsity()
{
  info() << "BSRFormat(computeSparsity): Compute sparsity of BSRMatrix";
  computeSparsityRowIndex();
  computeSparsityColumns();
}

}; // namespace Arcane::FemUtils
