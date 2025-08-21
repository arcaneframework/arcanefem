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
/*!
 * \brief View of a Matrix with CSR format.
 *
 * This view is a temporary object and is invalided when the underlying matrix
 * structure is modified.
 */
class CSRFormatView
{
 public:

  CSRFormatView() = default;
  CSRFormatView(Span<const Int32> rows,
                Span<const Int32> matrix_rows_nb_column,
                Span<const Int32> columns,
                Span<Real> values)
  : m_matrix_rows(rows)
  , m_matrix_rows_nb_column(matrix_rows_nb_column)
  , m_matrix_columns(columns)
  , m_values(values)
  {}

 public:

  constexpr ARCCORE_HOST_DEVICE Span<const Int32> rows() const { return m_matrix_rows; }
  constexpr ARCCORE_HOST_DEVICE Span<const Int32> rowsNbColumn() const { return m_matrix_rows_nb_column; }
  constexpr ARCCORE_HOST_DEVICE Span<const Int32> columns() const { return m_matrix_columns; }
  constexpr ARCCORE_HOST_DEVICE Span<Real> values() { return m_values; }

  constexpr ARCCORE_HOST_DEVICE Int32 nbRow() { return m_matrix_rows.size(); }
  constexpr ARCCORE_HOST_DEVICE Int32 nbColumn() { return m_matrix_columns.size(); }
  constexpr ARCCORE_HOST_DEVICE Int32 nbValue() { return m_values.size(); }

  constexpr ARCCORE_HOST_DEVICE Int32 row(Int32 index) { return m_matrix_rows[index]; }

 private:

  Span<const Int32> m_matrix_rows;
  Span<const Int32> m_matrix_rows_nb_column;
  Span<const Int32> m_matrix_columns;
  Span<Real> m_values;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
