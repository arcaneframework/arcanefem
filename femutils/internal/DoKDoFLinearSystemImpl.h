// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoKDoFLinearSystemImpl.h                                    (C) 2000-2026 */
/*                                                                           */
/* Implementation of IDoFLinearSystemImpl using a matrix with DoK format.    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_INTERNAL_DOKDOFLINEARSYSTEMIMPL_H
#define ARCANEFEM_FEMUTILS_INTERNAL_DOKDOFLINEARSYSTEMIMPL_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/VariableTypes.h>

#include "internal/DoFLinearSystemImplBase.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class DoKMatrix
{
 public:

  /*!
   * \brief Map to store values by Row/Column.
   *
   * This map has to be sorted if we want to reuse the internal structure because
   * the matrix filled has to be in the same order when we reuse it.
   */
  using RowColumnMap = OrderedRowColumnMap;
  using RowColumn = RowColumnMap::RowColumn;

 public:

  void addValue(Int32 row, Int32 column, Real value)
  {
    m_values_map.addValue({ row, column }, value);
  };
  void setValue(Int32 row, Int32 column, Real value)
  {
    m_forced_set_values_map.setValue({ row, column }, value);
  };
  void clearValues()
  {
    m_values_map.clear();
    m_forced_set_values_map.clear();
  }

 public:

  //! List of (i,j) values added to the matrix
  RowColumnMap m_values_map;
  //! List of (i,j) whose value is fixed. This will override added values in m_values_map.
  RowColumnMap m_forced_set_values_map;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Linear system implementation using a Dictionary of Keys (DoK) matrix.
 */
class DoKDoFLinearSystemImpl
: public DoFLinearSystemImplBase
{
  using RowColumn = DoKMatrix::RowColumn;

  class InternalCSRMatrix;

 public:

  DoKDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
  : DoFLinearSystemImplBase(dof_family, solver_name)
  {}
  ~DoKDoFLinearSystemImpl() override;

 public:

  void clearValues() override
  {
    DoFLinearSystemImplBase::clearValues();
    m_dok_matrix.clearValues();
  }

  void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (column.isNull())
      ARCANE_FATAL("Column is null");
    m_dok_matrix.addValue(row.localId(), column.localId(), value);
  }

  void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (column.isNull())
      ARCANE_FATAL("Column is null");
    m_dok_matrix.setValue(row.localId(), column.localId(), value);
  }

  void eliminateRow(DoFLocalId row, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    getEliminationInfo()[row] = ELIMINATE_ROW;
    getEliminationValue()[row] = value;
    info() << "EliminateRow row=" << row.localId() << " v=" << value;
  }

  void eliminateRowColumn(DoFLocalId row, Real value) override
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    getEliminationInfo()[row] = ELIMINATE_ROW_COLUMN;
    getEliminationValue()[row] = value;
    info() << "EliminateRowColumn row=" << row.localId() << " v=" << value;
  }

  void applyRHSTransformation() override;

  void setCSRValues(const CSRFormatView& csr_view) override
  {
    ARCANE_THROW(NotSupportedException, "DoKDoFLinearSystem does not support CSR format");
  }
  CSRFormatView& getCSRValues() override
  {
    ARCANE_THROW(NotSupportedException, "DoKDoFLinearSystem does not support CSR format");
  }
  bool hasSetCSRValues() const override { return false; }

 public:

  void setPrintFilling(bool v) { m_do_print_filling = v; }

  void fillRowColumnEliminationInfos();

  void convertToCSRMatrix();

  /*!
   * \brief Return the view of this matrix converted to CSR.
   *
   * You have to call convertToCSRMatrix() before, otherwise it will return
   * a null view.
   */
  CsrFormatMatrixView getCsrFormatMatrixView() const;

  /*!
   * \brief Visit all the non zero elements of the matrix and apply \a func.
   *
   * You need to make sure fillRowColumnEliminationInfos() has been called
   * before calling this method
   */
  template <typename Lambda> void
  visitDoKMatrix(const Lambda& func)
  {
    OrderedRowColumnMap& rc_elimination_map = _rowColumnEliminationMap();

    IItemFamily* dof_family = dofFamily();
    DoFInfoListView item_list_view(dof_family);

    auto& dof_elimination_info = getEliminationInfo();
    auto& dof_elimination_value = getEliminationValue();

    // Fill the matrix from the values of \a m_values_map
    // Skip (row,column) values which are part of an elimination.
    for (const auto& rc_value : m_dok_matrix.m_values_map) {
      RowColumn rc = rc_value.first;
      Real value = rc_value.second;
      DoF dof_row = item_list_view[rc.row_id];
      DoF dof_column = item_list_view[rc.column_id];

      Byte row_elimination_info = dof_elimination_info[dof_row];

      if (row_elimination_info == ELIMINATE_ROW)
        // Will be computed after this loop
        continue;
      if (rc_elimination_map.contains(rc))
        continue;

      // Check if value is forced for current RowColumn
      auto x = m_dok_matrix.m_forced_set_values_map.find(rc);
      if (x != m_dok_matrix.m_forced_set_values_map.end()) {
        info(4) << "FORCED VALUE R=" << rc.row_id << " C=" << rc.column_id
                << " old=" << value << " new=" << x->second;
        value = x->second;
      }

      func(dof_row, dof_column, value);
    }

    const bool do_print_filling = m_do_print_filling;

    // Apply Row or Row+Column elimination on Matrix
    // Phase 2: set the diagonal value for elimination row to 1.0
    ENUMERATE_ (DoF, idof, dof_family->allItems()) {
      DoF dof = *idof;
      if (!dof.isOwn())
        continue;
      Byte elimination_info = dof_elimination_info[dof];
      if (elimination_info == ELIMINATE_ROW || elimination_info == ELIMINATE_ROW_COLUMN) {
        Real elimination_value = dof_elimination_value[dof];
        if (do_print_filling)
          info() << "EliminateMatrix info=" << static_cast<int>(elimination_info) << " row="
                 << std::setw(4) << dof.localId() << " value=" << elimination_value;
        func(dof, dof, 1.0);
      }
    }
  }

 private:

  //! Container to store matrix values
  DoKMatrix m_dok_matrix;
  //! True is we want to print values during filling
  bool m_do_print_filling = false;
  // CSR Matrix used for conversion
  InternalCSRMatrix* m_internal_csr_matrix = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
