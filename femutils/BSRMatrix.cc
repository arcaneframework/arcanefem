#include "BSRMatrix.h"

namespace Arcane::FemUtils
{

BSRMatrix::BSRMatrix(ITraceMng* tm, const eMemoryRessource& mem_ressource)
: TraceAccessor(tm)
, m_values(mem_ressource)
, m_columns(mem_ressource)
, m_row_index(mem_ressource) {};

void BSRMatrix::initialize(Int8 block_size, Int32 nb_non_zero_value, Int32 nb_row, const RunQueue& queue)
{
  if (block_size <= 0 || nb_non_zero_value <= 0 || nb_row <= 0)
    ARCANE_THROW(ArgumentException, "BSRMatrix(initialize): block_size, nb_non_zero_value and nb_row should be positive and not nul");

  if (block_size > nb_row)
    ARCANE_THROW(ArgumentException, "BSRMatrix(initialize): block_size should be less than nb_row");

  info() << "BSRMatrix(initialize): Initiliaze BSRMatrix with block_size=" << block_size << ", nb_non_zero_value=" << nb_non_zero_value << ", nb_row=" << nb_row;

  m_block_size = block_size;
  m_values.resize(nb_non_zero_value);
  m_values.fill(0, &queue);
  m_columns.resize(nb_non_zero_value);
  m_row_index.resize(nb_row);
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

void BSRFormat::initialize()
{
  Int32 nb_node = m_mesh.nbNode();
  Int32 nb_dof = m_dofs_on_nodes.dofFamily()->nbItem();
  Int32 nb_non_zero_value = (nb_dof * nb_dof) * (2 * m_mesh.nbEdge() + nb_node);
  Int32 nb_row = nb_node * nb_dof;

  m_bsr_matrix.initialize(nb_dof, nb_non_zero_value, nb_row, m_queue);
}

void BSRFormat::computeSparsityRowIndex(const IndexedNodeNodeConnectivityView& node_node_cv)
{
  info() << "BSRFormat(computeSparsityRowIndex): Compute row index sparsity of BSRMatrix";
  auto command = makeCommand(m_queue);

  NumArray<Int32, MDDim1> out_data;
  out_data.resize(m_bsr_matrix.nbRow());
  auto copy_out_data = viewInOut(command, out_data);

  {
    auto command = makeCommand(m_queue);
    command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, m_mesh.allNodes())
    {
      copy_out_data[node_lid.asInt32()] = node_node_cv.nbNode(node_lid) + 1; // Add one to count the node itself
    };
  }
  m_queue.barrier();

  Accelerator::Scanner<Int32> scanner;
  scanner.exclusiveSum(&m_queue, out_data, m_bsr_matrix.rowIndex());
}

void BSRFormat::computeSparsityColumns(const IndexedNodeNodeConnectivityView& node_node_cv)
{
  info() << "BSRFormat(computeSparsityColumns): Compute columns sparsity of BSRMatrix";
  auto command = makeCommand(m_queue);

  auto out_row_index = viewIn(command, m_bsr_matrix.rowIndex());
  auto inout_columns = viewInOut(command, m_bsr_matrix.columns());

  {
    auto command = makeCommand(m_queue);
    command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, m_mesh.allNodes())
    {
      auto offset = out_row_index[node_lid.asInt32()];
      for (auto neighbor_lid : node_node_cv.nodeIds(node_lid)) {
        inout_columns[offset] = neighbor_lid.asInt32();
        ++offset;
      }
      inout_columns[offset] = node_lid.asInt32();
    };
  }
}

void BSRFormat::computeSparsity(const IndexedNodeNodeConnectivityView& node_node_cv)
{
  info() << "BSRFormat(computeSparsity): Compute sparsity of BSRMatrix";
  computeSparsityRowIndex(node_node_cv);
  computeSparsityColumns(node_node_cv);
  m_bsr_matrix.dump("bsr_matrix_dump.txt");
}

}; // namespace Arcane::FemUtils
