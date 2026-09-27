// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoKDoFLinearSystemImpl.cc                                   (C) 2000-2026 */
/*                                                                           */
/* Implementation of IDoFLinearSystemImpl using a matrix with DoK format.    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/Array.h>

#include "internal/DoKDoFLinearSystemImpl.h"
#include "CsrFormatMatrixView.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class DoKDoFLinearSystemImpl::InternalCSRMatrix
{
 public:

  UniqueArray<Real> matrix_values;
  UniqueArray<Int32> matrix_column_indexes;
  UniqueArray<Int32> row_indexes;
  UniqueArray<Int32> nb_value_per_row;

 public:

  CsrFormatMatrixView view() const
  {
    return CsrFormatMatrixView(row_indexes.view(), nb_value_per_row.view(),
                               matrix_column_indexes.view(), matrix_values.view());
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

DoKDoFLinearSystemImpl::
~DoKDoFLinearSystemImpl()
{
  delete m_internal_csr_matrix;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoKDoFLinearSystemImpl::
fillRowColumnEliminationInfos()
{
  OrderedRowColumnMap& rc_elimination_map = _rowColumnEliminationMap();
  rc_elimination_map.clear();
  DoFInfoListView item_list_view(dofFamily());

  auto& dof_elimination_info = getEliminationInfo();
  auto& dof_elimination_value = getEliminationValue();

  for (const auto& rc_value : m_dok_matrix.m_values_map) {
    RowColumn rc = rc_value.first;
    Real value = rc_value.second;
    DoF dof_row = item_list_view[rc.row_id];
    DoF dof_column = item_list_view[rc.column_id];
    Byte row_elimination_info = dof_elimination_info[dof_row];
    Byte column_elimination_info = dof_elimination_info[dof_column];
    if (row_elimination_info == ELIMINATE_ROW_COLUMN || column_elimination_info == ELIMINATE_ROW_COLUMN)
      rc_elimination_map.setValue({ rc.row_id, rc.column_id }, value);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoKDoFLinearSystemImpl::
applyRHSTransformation()
{
  const bool do_print_filling = m_do_print_filling;

  // Apply Row+Column elimination
  // Phase 1:
  // - subtract values of the RHS vector if Row+Column elimination
  _applyRowColumnEliminationToRHS(do_print_filling);

  IItemFamily* dof_family = dofFamily();

  auto& dof_elimination_info = getEliminationInfo();
  auto& dof_elimination_value = getEliminationValue();
  auto& rhs_variable = rhsVariable();

  // Apply Row or Row+Column elimination on RHS
  ENUMERATE_ (DoF, idof, dof_family->allItems()) {
    DoF dof = *idof;
    if (!dof.isOwn())
      continue;
    Byte elimination_info = dof_elimination_info[dof];
    if (elimination_info == ELIMINATE_ROW || elimination_info == ELIMINATE_ROW_COLUMN) {
      Real elimination_value = dof_elimination_value[dof];
      rhs_variable[dof] = elimination_value;
      if (do_print_filling)
        info() << "EliminateRHS info=" << static_cast<int>(elimination_info) << " row="
               << std::setw(4) << dof.localId() << " value=" << elimination_value;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoKDoFLinearSystemImpl::
convertToCSRMatrix()
{
  delete m_internal_csr_matrix;
  m_internal_csr_matrix = new InternalCSRMatrix();

  IItemFamily* dof_family = dofFamily();
  Int32 nb_row = dof_family->maxLocalId();
  DoFInfoListView item_list_view(dof_family);

  info() << "Convert To CSR Matrix nb_row=" << nb_row;

  UniqueArray<Int32> csr_matrix_nb_row(nb_row, 0);
  Int32 nb_value = 0;
  auto count_row = [&](DoF row, DoF, Real) {
    ++csr_matrix_nb_row[row.localId()];
    ++nb_value;
  };
  visitDoKMatrix(count_row);

  info() << "Convert To CSR Matrix nb_row=" << nb_row << " nb_value=" << nb_value;
  m_internal_csr_matrix->matrix_values.resize(nb_value);
  m_internal_csr_matrix->matrix_column_indexes.resize(nb_value);
  m_internal_csr_matrix->row_indexes.resize(nb_row + 1);

  // Now we know the number of columns per row and total number of non zeros.
  SmallSpan<Real> csr_matrix_values = m_internal_csr_matrix->matrix_values.view();
  SmallSpan<Int32> csr_matrix_column_indexes = m_internal_csr_matrix->matrix_column_indexes.view();
  SmallSpan<Int32> csr_row_indexes = m_internal_csr_matrix->row_indexes.view();
  Int32 current_index = 0;
  for (Int32 i = 0; i < nb_row; ++i) {
    csr_row_indexes[i] = current_index;
    //info() << "ROW i=" << i << " nb_row=" << csr_matrix_nb_row[i] << " index=" << csr_row_indexes[i];
    current_index += csr_matrix_nb_row[i];
  }
  csr_row_indexes[nb_row] = current_index;

  // Fill the column indexes and the values of the CSR Matrix
  m_internal_csr_matrix->nb_value_per_row.resize(nb_row, 0);
  SmallSpan<Int32> work_nb_value_per_row = m_internal_csr_matrix->nb_value_per_row.view();
  auto set_csr_matrix_value = [&](DoF row, DoF column, Real value) {
    Int32 row_id = row.localId();
    Int32 index = csr_row_indexes[row_id] + work_nb_value_per_row[row_id];
    //info() << "SET_VALUE (" << row_id << ", " << column.localId() << ") = " << value << " index=" << index;
    csr_matrix_column_indexes[index] = column.localId();
    csr_matrix_values[index] = value;
    ++work_nb_value_per_row[row_id];
  };
  visitDoKMatrix(set_csr_matrix_value);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

CsrFormatMatrixView DoKDoFLinearSystemImpl::
getCsrFormatMatrixView() const
{
  if (!m_internal_csr_matrix)
    return {};
  return m_internal_csr_matrix->view();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
