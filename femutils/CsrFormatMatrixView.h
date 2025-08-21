// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrFormatMatrixView.h                                       (C) 2022-2025 */
/*                                                                           */
/* View of a Matrix with CSR format.                                         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_CSRFORMATMATRIXVIEW_H
#define ARCANEFEM_FEMUTILS_CSRFORMATMATRIXVIEW_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/ArrayView.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class CsrRowColumnIterator;
class CsrFormatMatrixView;
class CsrRow;
class CsrFormat;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Index in the RowColumn list of a CSR Matrix.
 */
class CsrRowColumnIndex
{
 public:

  using IndexType = Int32;

 public:

  CsrRowColumnIndex() = default;
  explicit constexpr ARCCORE_HOST_DEVICE CsrRowColumnIndex(IndexType index)
  : m_index(index)
  {}

 public:

  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE IndexType value() const { return m_index; }
  constexpr ARCCORE_HOST_DEVICE operator IndexType() const { return m_index; }

 private:

  IndexType m_index = -1;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Represents an iterator over the indexes of the columns of a CsrRow.
 */
class CsrRowColumnIterator
{
  friend CsrRow;

 public:

  CsrRowColumnIterator() = default;

 private:

  explicit constexpr ARCCORE_HOST_DEVICE CsrRowColumnIterator(Int32 index)
  : m_index(index)
  {}

 public:

  constexpr ARCCORE_HOST_DEVICE CsrRowColumnIndex operator*() const { return CsrRowColumnIndex(m_index); }
  constexpr ARCCORE_HOST_DEVICE void operator++() { ++m_index; }
  friend constexpr ARCCORE_HOST_DEVICE bool operator!=(const CsrRowColumnIterator& lhs, const CsrRowColumnIterator& rhs)
  {
    return lhs.m_index != rhs.m_index;
  }

 private:

  Int32 m_index = -1;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Represents a row of a CSR Matrix.
 */
class CsrRow
{
  friend CsrFormatMatrixView;

 public:

  CsrRow() = default;

 private:

  constexpr ARCCORE_HOST_DEVICE CsrRow(Int32 begin, Int32 end)
  : m_begin(begin)
  , m_end(end)
  {}

 public:

  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE CsrRowColumnIterator begin() const { return CsrRowColumnIterator(m_begin); }
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE CsrRowColumnIterator end() const { return CsrRowColumnIterator(m_end); }

 private:

  Int32 m_begin = -1;
  Int32 m_end = -1;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief View of a Matrix with CSR format.
 *
 * This view is a temporary object and is invalided when the underlying matrix
 * structure is modified.
 */
class CsrFormatMatrixView
{
  friend CsrFormat;

 public:

  CsrFormatMatrixView() = default;

 private:

  CsrFormatMatrixView(SmallSpan<const Int32> rows,
                      SmallSpan<const Int32> matrix_rows_nb_column,
                      SmallSpan<const Int32> columns,
                      SmallSpan<Real> values)
  : m_matrix_rows(rows)
  , m_matrix_rows_nb_column(matrix_rows_nb_column)
  , m_matrix_columns(columns)
  , m_values(values)
  {}

 public:

  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Span<const Int32> rows() const { return m_matrix_rows; }
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Span<const Int32> rowsNbColumn() const { return m_matrix_rows_nb_column; }
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Span<const Int32> columns() const { return m_matrix_columns; }
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Span<Real> values() const { return m_values; }

  //! Number of the rows in the matrix
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Int32 nbRow() const { return m_matrix_rows.size(); }
  //! Number of the values in the matrix
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Int32 nbColumn() const { return m_matrix_columns.size(); }
  //! Number of the values in the matrix
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Int32 nbValue() const { return m_values.size(); }

  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Int32 row(Int32 index) const { return m_matrix_rows[index]; }

  //! Local index of the column for the given RowColumnIndex \a rc_index
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Int32 column(CsrRowColumnIndex rc_index) const { return m_matrix_columns[rc_index]; }
  //! Value of the matrix for the given RowColumnIndex \a rc_index
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE Real& value(CsrRowColumnIndex rc_index) const { return m_values[rc_index]; }

  //! Range of CsrRowColumnIndex for the given row \a row
  [[nodiscard]] constexpr ARCCORE_HOST_DEVICE CsrRow rowRange(Int32 row) const
  {
    // TODO: look if we can use rowsNbColumn() to compute 'end'.
    auto begin = m_matrix_rows[row];
    auto end = (row == (nbRow() - 1)) ? nbColumn() : m_matrix_rows[row + 1];
    return { begin, end };
  }

 private:

  SmallSpan<const Int32> m_matrix_rows;
  SmallSpan<const Int32> m_matrix_rows_nb_column;
  SmallSpan<const Int32> m_matrix_columns;
  SmallSpan<Real> m_values;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

//! Old name to keep compatibility with existing code.
using CSRFormatView = CsrFormatMatrixView;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
