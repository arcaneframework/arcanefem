// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrFormatMatrix.cc                                          (C) 2022-2026 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/NumArray.h>

#include <arcane/core/VariableTypes.h>
#include <arcane/core/IItemFamily.h>

#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/accelerator/NumArrayViews.h>
#include <arcane/accelerator/RunCommandLoop.h>

#include "CsrFormatMatrix.h"
#include "DoFLinearSystem.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrFormat::
initialize(IItemFamily* dof_family, Int32 nnz, Int32 nb_row, RunQueue& queue)
{
  info() << "Initialize CsrFormat: nb_non_zero=" << nnz << " nb_row=" << nb_row;

  eMemoryRessource mem_resource = queue.memoryResource();
  m_matrix_row = NumArray<Int32, MDDim1>(mem_resource);
  m_matrix_column = NumArray<Int32, MDDim1>(mem_resource);
  m_matrix_value = NumArray<Real, MDDim1>(mem_resource);
  m_matrix_rows_nb_column = NumArray<Int32, MDDim1>(mem_resource);

  m_matrix_row.resize(nb_row);
  m_matrix_column.resize(nnz);
  m_matrix_value.resize(nnz);
  m_matrix_row.fill(-1, &queue);
  m_matrix_column.fill(-1, &queue);
  m_matrix_value.fill(0, &queue);
  m_matrix_rows_nb_column.resize(nb_row);
  m_matrix_rows_nb_column.fill(0, &queue);
  m_dof_family = dof_family;
  m_nnz = nnz;
  info() << "Filling CSR Matrix with zeros";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrFormat::
initialize(IItemFamily* dof_family, NumArray<Int32, MDDim1>&& rows_index, NumArray<Int32, MDDim1>&& columns, RunQueue& queue)
{
  info() << "Initialize CsrFormat with rows_index : nb_row=" << (rows_index.extent0() - 1);

  if (rows_index.extent0() == 0)
    ARCANE_FATAL("rows_index is empty");

  eMemoryRessource mem_resource = queue.memoryResource();
  if (rows_index.memoryResource() != mem_resource)
    ARCANE_FATAL("Bad memory resource '{0}' for 'rown_index' (expected={1})", rows_index.memoryResource(), mem_resource);
  if (columns.memoryResource() != mem_resource)
    ARCANE_FATAL("Bad memory resource '{0}' for 'columns' (expected={1})", columns.memoryResource(), mem_resource);

  Int32 nb_row = rows_index.extent0() - 1;
  Int32 nnz = rows_index[nb_row];
  if (nnz != columns.extent0())
    ARCANE_FATAL("Incoherent sizes for columns (from_rows={0} from_columns={1})", nnz, columns.extent0());

  m_matrix_row = rows_index;
  m_matrix_column = columns;

  m_matrix_value = NumArray<Real, MDDim1>(mem_resource);
  m_matrix_rows_nb_column = NumArray<Int32, MDDim1>(mem_resource);

  // TODO: Make the filling optional
  m_matrix_value.resize(nnz);
  m_matrix_value.fill(0, &queue);

  m_matrix_rows_nb_column.resize(nb_row);
  m_matrix_rows_nb_column.fill(0, &queue);
  for (Int32 i = 0; i < nb_row; ++i)
    m_matrix_rows_nb_column[i] = rows_index[i + 1] - rows_index[i];

  m_dof_family = dof_family;
  m_nnz = nnz;

  checkValid();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrFormat::
translateToLinearSystem(DoFLinearSystem& linear_system, const RunQueue& queue)
{
  bool do_set_csr = linear_system.hasSetCSRValues();
  info() << "TranslateToLinearSystem this=" << this << " is_csr=" << do_set_csr;

  const Int32 nb_row = m_matrix_row.dim1Size();
  const Int32 matrix_column_size = m_matrix_column.dim1Size();

  // When using CSR format, we need to know the number of non zero values for
  // each row.
  // NOTE: it should be possible to compute that in setCoordinates().
  // and this value is constant if the structure of the matrix do not change
  // so we can store these values instead of recomputing them.
  if (do_set_csr) {
    m_matrix_rows_nb_column.resize(nb_row);
    auto command = makeCommand(queue);
    auto out_matrix_rows_nb_column = viewOut(command, m_matrix_rows_nb_column);
    auto in_matrix_rows = viewIn(command, m_matrix_row);
    command << RUNCOMMAND_LOOP1(iter, nb_row)
    {
      auto [i] = iter();
      Int32 nb_column = 0;
      if (((i + 1) < nb_row) && (in_matrix_rows(i) == in_matrix_rows(i + 1))) {
        out_matrix_rows_nb_column[0];
        return;
      }
      for (Int32 j = in_matrix_rows(i); ((i + 1) < nb_row && j < in_matrix_rows(i + 1)) || ((i + 1) == nb_row && j < matrix_column_size); j++) {
        ++nb_column;
      }
      out_matrix_rows_nb_column[i] = nb_column;
    };
    CSRFormatView csr_view(view());
    linear_system.setCSRValues(csr_view);
    return;
  }

  for (Int32 i = 0; i < nb_row; i++) {
    m_matrix_rows_nb_column[i] = 0;
    if (((i + 1) < nb_row) && (m_matrix_row(i) == m_matrix_row(i + 1)))
      continue;
    for (Int32 j = m_matrix_row(i); ((i + 1) < nb_row && j < m_matrix_row(i + 1)) || ((i + 1) == nb_row && j < matrix_column_size); j++) {
      if (DoFLocalId(m_matrix_column(j)).isNull())
        continue;
      //info() << "Add: (" << i << ", " << m_matrix_column(j) << " v=" << m_matrix_value(j);
      linear_system.matrixAddValue(DoFLocalId(i), DoFLocalId(m_matrix_column(j)), m_matrix_value(j));
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

CsrFormatMatrixView CsrFormat::
view()
{
  return CSRFormatView(m_matrix_row.to1DSmallSpan(), m_matrix_rows_nb_column.to1DSmallSpan(),
                       m_matrix_column.to1DSmallSpan(), m_matrix_value.to1DSmallSpan());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrFormat::
checkValid(bool force) const
{
  if (!arcaneIsCheck() && !force)
    return;
  Int32 nb_row = nbRow();
  Int32 n1 = m_matrix_row.extent0();
  if (n1 != (nb_row + 1))
    ARCCORE_FATAL("Bad size '{0}' for rowIndexes() (expected value = {1})", n1, nb_row + 1);
  Int32 nb_value = m_matrix_row[nb_row];
  Int32 n2 = m_matrix_column.extent0();
  Int32 n3 = m_matrix_value.extent0();
  if (n2 != nb_value)
    ARCCORE_FATAL("Bad size '{0}' for columns() (expected value = {1})", n2, nb_value);
  if (n3 != nb_value)
    ARCCORE_FATAL("Bad size '{0}' for values() (expected value = {1})", n3, nb_value);
  auto mem_resource = m_matrix_row.memoryResource();
  bool can_compare = (mem_resource != eMemoryResource::Device);
  if (can_compare) {
    for (Int32 i = 0; i < nb_row; ++i) {
      Int32 x0 = m_matrix_rows_nb_column[i];
      Int32 x1 = m_matrix_row[i + 1] - m_matrix_row[i];
      if (x0 != x1)
        ARCCORE_FATAL("Bad number of column for row='{0}' v={1} expected={2}", i, x1, x0);
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrFormat::
printMatrix(std::string fileName)
{
  ofstream file(fileName);
  file << "size :" << m_nnz << "\n";
  for (auto i = 0; i < m_matrix_row.dim1Size(); i++) {
    file << m_matrix_row(i) << " ";
    for (Int32 j = m_matrix_row(i) + 1; (i + 1 < m_matrix_row.dim1Size() && j < m_matrix_row(i + 1)) || (i + 1 == m_matrix_row.dim1Size() && j < m_matrix_column.dim1Size()); j++) {
      file << "  ";
    }
  }
  file << "\n";
  for (auto i = 0; i < m_nnz; i++) {
    file << m_matrix_column(i) << " ";
  }
  file << "\n";
  for (auto i = 0; i < m_nnz; i++) {
    file << m_matrix_value(i) << " ";
  }
  file << "\n";
  file.close();
}

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane
{
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Convert CSR format rows into COO format rows
 */
void FemUtils::
_translateCSRToCOO(Span<const Int32> csr_rows, SmallSpan<Int32> coo_rows,
                   const RunQueue& queue)
{
  const Int32 nb_value = coo_rows.size();
  const Int32 nb_row = csr_rows.size();

  {
    auto command = makeCommand(queue);
    command << RUNCOMMAND_LOOP1(iter, nb_row)
    {
      auto [i] = iter();
      if (i != (nb_row - 1)) {
        for (int j = csr_rows[i]; j < csr_rows[i + 1]; j++)
          coo_rows[j] = i;
      }
      else {
        // The last iteration fill the remaining values
        for (int j = csr_rows[nb_row - 1]; j < nb_value; j++)
          coo_rows[j] = nb_row - 1;
      }
    };
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
