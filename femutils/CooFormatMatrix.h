// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrFormatMatrix.cc                                     (C) 2022-2024 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/NumArray.h>

#include <arcane/VariableTypes.h>
#include <arcane/IItemFamily.h>

#include <arcane/aleph/AlephTypesSolver.h>
#include <arcane/aleph/Aleph.h>

#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "arcane_version.h"

#include <iostream>
#include <fstream>

#include "arcane/accelerator/NumArrayViews.h"

namespace Arcane::FemUtils
{
using namespace Arcane;

namespace ax = Arcane::Accelerator;
class CooFormat : TraceAccessor
{
 public:

  CooFormat(ISubDomain* sd)
  : TraceAccessor(sd->traceMng())
  {
    info() << "Creating COO Matrix";
  }

  void initialize(IItemFamily* dof_family, Int32 nnz)
  {
    m_matrix_row.resize(nnz);
    m_matrix_column.resize(nnz);
    m_matrix_value.resize(nnz);
    m_matrix_row.fill(0);
    m_matrix_column.fill(0);
    m_matrix_value.fill(0);
    m_dof_family = dof_family;
    m_last_value = 0;
    m_nnz = nnz;
    info() << "Filling COO Matrix with zeros";
  }

  /**
  * @brief
  *
  * @param row
  * @param column
  * @param value
  */
  void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value)
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (column.isNull())
      ARCANE_FATAL("Column is null");
    if (value == 0.0)
      return;
    m_matrix_value(indexValue(row, column)) += value;
  }

  /**
 * @brief  Translate to Arcane linear system
 *
 * @param linear_system
 */
  void translateToLinearSystem(DoFLinearSystem& linear_system)
  {
    for (Int32 i = 0; i < m_nnz; i++) {
      linear_system.matrixAddValue(DoFLocalId(m_matrix_row(i)), DoFLocalId(m_matrix_column(i)), m_matrix_value(i));
    }
  }

  /**
 * @brief function to print the current content of the csr matrix
 *
 * @param fileName
 * @param nonzero if set to true, print only the non zero values
 */
  void
  printMatrix(std::string fileName, bool nonzero)
  {
    ofstream file(fileName);
    file << "size :" << m_matrix_row.extent0() << "\n";
    for (auto i = 0; i < m_matrix_row.extent0(); i++) {
      if (nonzero && m_matrix_value(i) == 0)
        continue;
      file << m_matrix_row(i) << " ";
    }
    file << "\n";
    for (auto i = 0; i < m_nnz; i++) {
      if (nonzero && m_matrix_value(i) == 0)
        continue;
      file << m_matrix_column(i) << " ";
    }
    file << "\n";
    for (auto i = 0; i < m_nnz; i++) {
      if (nonzero && m_matrix_value(i) == 0)
        continue;
      file << m_matrix_value(i) << " ";
    }
    file.close();
  }

  void setCoordinates(DoFLocalId row, DoFLocalId column)
  {
    m_matrix_row(m_last_value) = row.localId();
    m_matrix_column(m_last_value) = column.localId();
    m_last_value++;
  }

  void sort()
  {
    sortMatrix(true, 0, m_matrix_row.extent0() - 1);
    Int32 begin = 0;
    for (Int32 i = 0; i < m_matrix_row.extent0(); i++) {
      if (i + 1 == m_matrix_row.extent0() || m_matrix_row(i + 1) != m_matrix_row(begin)) {
        sortMatrix(false, begin, i);
        begin = i + 1;
      }
    }
  }

 public:

  Int32 m_nnz;
  // To become parallelizable, have all the index
  // inside a queue that would gradually pop ?
  // or link the idnex to the index of the core ?
  Int32 m_last_value;
  NumArray<Int32, MDDim1> m_matrix_row;
  NumArray<Int32, MDDim1> m_matrix_column;
  NumArray<Real, MDDim1> m_matrix_value;
  IItemFamily* m_dof_family = nullptr;

  /*
  getValue return the Value at the (row, column) coordinates.
*/
  Int32 getValue(DoFLocalId row, DoFLocalId column)
  {
    return m_matrix_value(indexValue(row, column));
  }

  /**
  * @brief binSearchRow is a binary search through the row to get the
  * leftmost corresponding index.
  *
  * @param row
  * @return Int32
  */
  Int32 binSearchRow(Int32 row)
  {
    Int32 begin = 0;
    Int32 end = m_matrix_row.totalNbElement() - 1;
    while (begin <= end) {
      Int32 mid = begin + (end - begin) / 2;
      if (row == m_matrix_row(mid)) {
        while (mid - 1 >= 0 && m_matrix_row(mid - 1) == row) {
          mid--;
        }
        return mid;
      }
      if (row > m_matrix_row(mid)) {
        begin = mid + 1;
      }
      if (row < m_matrix_row(mid)) {
        end = mid - 1;
      }
    }
    return -1;
  }

  /**
  * @brief indexValue is a Binsearch through the row and then the column
  *  to get the index of the corresponding value.
  *
  * @param row
  * @param column
  * @return Int32
  */
  Int32 indexValue(Int32 row, Int32 column)
  {

    Int32 i = binSearchRow(row);
    while (i != m_matrix_row.totalNbElement() && m_matrix_row(i) == row) {
      if (m_matrix_column(i) == column)
        return i;
      i++;
    }
    //binsearch only on the row and iterate through the column
    /*
    while (begin <= end) {
      /*
      Int32 mid = begin + (end - begin) / 2;
      if (column == m_matrix_column(mid)) {
        return mid;
      }
      if (column > m_matrix_column(mid)) {
        begin = mid + 1;
      }
      if (column < m_matrix_column(mid)) {
        end = mid - 1;
      }
  }
      */
    return -1;
  }

  /**
 * @brief Quicksort algorithm for the CSR Matrix
 *
 * @param is_row
 * @param start
 * @param end
 */
  void
  sortMatrix(bool is_row, Int32 start, Int32 end)
  {
    if (start >= end) {
      return;
    }

    int pivot = partition(is_row, start, end);

    sortMatrix(is_row, start, pivot - 1);

    sortMatrix(is_row, pivot + 1, end);
  }

  /**
 * @brief Partition helper for the quickSort
 *
 * @param is_row
 * @param start
 * @param end
 * @return Int32
 */
  Int32 partition(bool is_row, Int32 start, Int32 end)
  {
    Int32 pivot;
    if (is_row)
      pivot = m_matrix_row[end];
    else
      pivot = m_matrix_column[end];

    Int32 pIndex = start;

    for (Int32 i = start; i < end; i++) {
      if ((is_row && m_matrix_row[i] <= pivot) || (!is_row && m_matrix_column[i] <= pivot)) {

        swap(is_row, i, pIndex);
        pIndex++;
      }
    }

    swap(is_row, pIndex, end);

    return pIndex;
  }

  /**
 * @brief Swap helper for the quickSort
 *
 * @param is_row
 * @param i
 * @param j
 */
  void swap(bool is_row, Int32 i, Int32 j)
  {
    if (is_row) {
      Int32 tmp = m_matrix_row(i);
      m_matrix_row(i) = m_matrix_row(j);
      m_matrix_row(j) = tmp;
    }
    Int32 tmp = m_matrix_column(i);
    m_matrix_column(i) = m_matrix_column(j);
    m_matrix_column(j) = tmp;
    Real tmp_val = m_matrix_value(i);
    m_matrix_value(i) = m_matrix_value(j);
    m_matrix_value(j) = tmp_val;
  }
};

ARCCORE_HOST_DEVICE Int32 findIndexBinarySearch(Int32 row, Int32 col,
                                                const auto& in_row_coo,
                                                const auto& in_col_coo,
                                                Int32 row_length)
{
  Int32 begin = 0;
  Int32 end = row_length - 1;
  Int32 i = -1; // Default value in case the index is not found

  // Binary search to find the first occurrence of `row`
  while (begin <= end) {
    Int32 mid = begin + (end - begin) / 2;

    if (row == in_row_coo(mid)) {
      // Move back to the first occurrence of the row
      while (mid > 0 && in_row_coo(mid - 1) == row) {
        mid--;
      }
      i = mid;
      break;
    }

    if (row > in_row_coo(mid)) {
      begin = mid + 1;
    }
    else {
      end = mid - 1;
    }
  }

  // If no occurrence of the row was found, return -1
  if (i == -1) {
    return -1;
  }

  // Search for the matching column in the found row's block
  while (i < row_length && in_row_coo(i) == row) {
    if (in_col_coo(i) == col) {
      return i;
    }
    i++;
  }

  // If no matching column is found, return -1
  return -1;
}

} // namespace Arcane::FemUtils
