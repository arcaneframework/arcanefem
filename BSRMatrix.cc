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

  info() << "BSRMatrix(initialize): Create BSRMatrix with block_size=" << block_size << ", nb_non_zero_value=" << nb_non_zero_value << ", nb_row=" << nb_row;

  m_block_size = block_size;
  m_values.resize(nb_non_zero_value);
  m_values.fill(0, &queue);
  m_columns.resize(nb_non_zero_value);
  m_row_index.resize(nb_row);
}

void BSRFormat::initialize()
{
  Int32 nb_node = m_mesh.nbNode();
  Int32 nb_dof = m_dofs_on_nodes.dofFamily()->nbItem();
  Int32 nb_non_zero_value = (nb_dof * nb_dof) * (2 * m_mesh.nbEdge() + nb_node);
  Int32 nb_row = nb_node * nb_dof;

  m_bsr_matrix.initialize(nb_dof, nb_non_zero_value, nb_row, m_queue);
}

void BSRFormat::computeSparsity()
{
  auto command = makeCommand(m_queue);

  NumArray<Int32, MDDim1> out_data;
  out_data.resize(m_bsr_matrix.nbRow());

  auto copy_out_data = viewInOut(command, out_data);
}

}; // namespace Arcane::FemUtils
