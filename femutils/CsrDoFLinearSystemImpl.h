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
#ifndef ARCANEFEM_FEMUTILS_CSRDOFLINEARSYSTEMIMPL_H
#define ARCANEFEM_FEMUTILS_CSRDOFLINEARSYSTEMIMPL_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/VariableTypes.h>
#include <arcane/accelerator/core/Runner.h>

#include "DoFLinearSystem.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Implementation of IDoFLinearSystemImpl using a matrix with CSR format.
 */
class CsrDoFLinearSystemImpl
: public TraceAccessor
, public IDoFLinearSystemImpl
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

  VariableDoFReal& solutionVariable() override
  {
    return m_dof_variable;
  }

  VariableDoFReal& rhsVariable() override
  {
    return m_rhs_variable;
  }

  void clearValues() override
  {
    info() << "[Hypre-Info]: Clear values";
    m_csr_view = {};
    m_dof_forced_info.fill(false);
    m_dof_elimination_info.fill(ELIMINATE_NONE);
    m_dof_elimination_value.fill(0);
  }

  CSRFormatView& getCSRValues() override { return m_csr_view; };
  VariableDoFBool& getForcedInfo() override { return m_dof_forced_info; }
  VariableDoFReal& getForcedValue() override { return m_dof_forced_value; }
  VariableDoFByte& getEliminationInfo() override { return m_dof_elimination_info; }
  VariableDoFReal& getEliminationValue() override { return m_dof_elimination_value; }

  void setCSRValues(const CSRFormatView& csr_view) override
  {
    m_csr_view = csr_view;
  }

  bool hasSetCSRValues() const override { return true; }

  void setRunner(const Runner& r) override { m_runner = r; }
  Runner runner() const override { return m_runner; }

 public:

  IItemFamily* dofFamily() const { return m_dof_family; }

 private:

  // TODO: make all these fields private

  IItemFamily* m_dof_family = nullptr;
  VariableDoFReal m_rhs_variable;
  VariableDoFReal m_dof_variable;
  Runner m_runner;

  CSRFormatView m_csr_view;

  VariableDoFBool m_dof_forced_info;
  VariableDoFReal m_dof_forced_value;
  static constexpr Byte ELIMINATE_NONE = 0;
  static constexpr Byte ELIMINATE_ROW = 1;
  VariableDoFByte m_dof_elimination_info;
  VariableDoFReal m_dof_elimination_value;

 public:

  // These methods should be private but has to be public because of NVidia compiler
  void _applyRowElimination();
  void _applyForcedValuesToLhs();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
