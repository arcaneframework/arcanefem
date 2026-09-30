// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrFormatMatrix.h                                           (C) 2022-2026 */
/*                                                                           */
/* Container for Matrix using CSR format.                                    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_CSRFORMATMATRIX_H
#define ARCANEFEM_FEMUTILS_CSRFORMATMATRIX_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/NumArray.h>

#include <arcane/VariableTypes.h>
#include <arcane/IItemFamily.h>

#include "FemUtils.h"
#include "CsrFormatMatrixView.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
class DoFLinearSystem;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Matrix using Compressed Sparse Row (CSR) format.
 *
 * You have to call initialize() before using this matrix.
 */
class CsrFormat
: public TraceAccessor
{
 public:

  explicit CsrFormat(Arcane::ITraceMng* tm)
  : TraceAccessor(tm)
  {
  }

 public:

  void initialize(IItemFamily* dof_family, Int32 nnz, Int32 nbRow, RunQueue& queue);

  /*!
   * \brief Initialize the matrix.
   *
   * The values of the arrays \a rows_index and \a columns will be moved into its class
   * and should not be user after. It should have be allocated using queue.memoryRessource().
   *
   * The number of rows of the matrix will be equal to `rows_index.size()-1` and
   * the number of column for a given \a row is `rows_index[row+1]-rows_index[row]`.
   * So the number of non-zero is equal to `rows_index[nb_row]` and should be
   * equal to `columns.size()`.
   */
  void initialize(IItemFamily* dof_family, NumArray<Int32, MDDim1>&& rows_index, NumArray<Int32, MDDim1>&& columns, RunQueue& queue);

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

  Int32 indexValue(DoFLocalId row, DoFLocalId column)
  {
    Int32 begin = m_matrix_row(row.localId());
    Int32 end = 0;
    if (row.localId() == m_matrix_row.extent0() - 1) {

      end = m_matrix_column.extent0();
    }
    else {

      end = m_matrix_row(row + 1);
    }
    for (Int32 i = begin; i < end; i++) {
      if (m_matrix_column(i) == column.localId()) {
        return i;
      }
    }
    return -1;
  }

  //! Number of rows in the matrix
  constexpr Int32 nbRow() const { return m_matrix_rows_nb_column.extent0(); }

  /**
   * @brief
   *
   * @param linear_system
   */
  void translateToLinearSystem(DoFLinearSystem& linear_system, const RunQueue& queue);

  /**
   * @brief function to print the current content of the csr matrix
   *
   * @param fileName
   * @param nonzero if set to true, print only the non zero values
   */
  void printMatrix(std::string fileName);

  // Warning : does not support empty row (or does it ?)
  void setCoordinates(DoFLocalId row, DoFLocalId column)
  {
    Int32 row_lid = row.localId();
    if (m_matrix_row(row_lid) == -1) {
      m_matrix_row(row_lid) = m_last_value;
    }
    m_matrix_column(m_last_value) = column.localId();
    m_last_value++;
  }

  void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value)
  {
    m_matrix_value(indexValue(row, column)) = value;
  }

  //! View of the matrix
  CsrFormatMatrixView view();

  /*!
   * \brief Check that sizes are valid:
   * - rowIndexes().size() = nbRow() + 1;
   * - columns().size() = rowIndexes[nbRow()];
   * - values().size() = rowIndexes[nbRow()];
   *
   * \note: At the moment theses properties are not always verified.
   */
  void checkValid(bool force = false) const;

 public:

  Int32 m_nnz = 0;
  // To become parallelizable, have all the index
  // inside a queue that would gradually pop ?
  // or link the idnex to the index of the core ?
  Int32 m_last_value = 0;
  NumArray<Int32, MDDim1> m_matrix_row;
  NumArray<Int32, MDDim1> m_matrix_column;
  NumArray<Real, MDDim1> m_matrix_value;
  //! Nombre de colonnes de chaque lignes.
  NumArray<Int32, MDDim1> m_matrix_rows_nb_column;
  IItemFamily* m_dof_family = nullptr;

  //! Return the Value at the (row, column) coordinates.
  Int32 getValue(DoFLocalId row, DoFLocalId column)
  {
    return m_matrix_value(indexValue(row, column));
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
