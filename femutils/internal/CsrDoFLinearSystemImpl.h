// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrDoFLinearSystemImpl.h                                    (C) 2022-2025 */
/*                                                                           */
/* Implementation of IDoFLinearSystemImpl using a matrix with CSR format.    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_INTERNAL_CSRDOFLINEARSYSTEMIMPL_H
#define ARCANEFEM_FEMUTILS_INTERNAL_CSRDOFLINEARSYSTEMIMPL_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/VariableTypes.h>
#include <arcane/accelerator/core/Runner.h>

#include "internal/DoFLinearSystemImplBase.h"
#include "CsrFormatMatrixView.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \internal
 * \brief Implementation of IDoFLinearSystemImpl using a matrix with CSR format.
 */
class CsrDoFLinearSystemImpl
: public DoFLinearSystemImplBase
{
 public:

  CsrDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name);

 public:

  Int32 indexValue(DoFLocalId row_lid, DoFLocalId column_lid)
  {
    auto begin = m_csr_view.rows()[row_lid];
    auto end = row_lid == m_csr_view.nbRow() - 1 ? m_csr_view.nbColumn() : m_csr_view.row(row_lid + 1);
    for (auto i = begin; i < end; ++i)
      if (m_csr_view.columns()[i] == column_lid)
        return i;
    return -1;
  }

 public:

  void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    m_csr_view.values()[indexValue(row, column)] += value;
  }

  void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    m_csr_view.values()[indexValue(row, column)] = value;
  }

  void eliminateRow(DoFLocalId row, Real value) override
  {
    ARCANE_THROW(NotImplementedException, "");
  }

  void eliminateRowColumn(DoFLocalId row, Real value) override
  {
    ARCANE_THROW(NotImplementedException, "");
  }

  void applyMatrixTransformation() override;
  void applyRHSTransformation() override;

  void clearValues() override
  {
    info() << "[CsrImpl]: Clear values";
    DoFLinearSystemImplBase::clearValues();
    m_csr_view = {};
  }

  CSRFormatView& getCSRValues() override { return m_csr_view; };

  void setCSRValues(const CSRFormatView& csr_view) override
  {
    m_csr_view = csr_view;
  }

  bool hasSetCSRValues() const override { return true; }

 private:

  CSRFormatView m_csr_view;
  bool m_has_row_column_elimination = false;

 public:

  // These methods should be private but has to be public because of NVidia compiler
  void _applyForcedValuesToLhs();
  void _fillRowColumnEliminationInfos();
  void _applyRowEliminationOnMatrix();
  void _applyRowColumnEliminationOnMatrix();
  void _applyRowOrRowColumnEliminationOnRHS();

 protected:

  void _applyRowOrRowColumnEliminationOnMatrix();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
