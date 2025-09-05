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

#include "internal/CsrDoFLinearSystemImpl.h"

#include <arcane/utils/ITraceMng.h>
#include <arcane/core/IItemFamily.h>

#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/accelerator/core/Runner.h>

#include <arcane/accelerator/VariableViews.h>
#include <arcane/accelerator/RunCommandLoop.h>

#include "ArcaneFemFunctionsGpu.h"

// NOTE:
// DoF family must be compacted (i.e maxLocalId()==nbItem()) and sorted
// for this implementation to works.

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

CsrDoFLinearSystemImpl::
CsrDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
: DoFLinearSystemImplBase(dof_family, solver_name)
{}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrDoFLinearSystemImpl::
_applyForcedValuesToLhs()
{
  IItemFamily* dof_family = dofFamily();
  auto nb_dof = dof_family->nbItem();

  RunQueue queue = makeQueue(runner());
  auto command = makeCommand(queue);

  auto in_out_forced_info = Accelerator::viewInOut(command, getForcedInfo());
  auto in_out_forced_value = Accelerator::viewInOut(command, getForcedValue());
  auto csr_view = m_csr_view;

  command << RUNCOMMAND_LOOP1(iter, nb_dof)
  {
    auto [i] = iter();
    DoFLocalId dof_id(i);
    if (in_out_forced_info[dof_id]) {
      auto index = csr_view.tryFindColumnInRow(dof_id, dof_id);
      csr_view.value(index) = in_out_forced_value[dof_id];
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrDoFLinearSystemImpl::
_applyRowOrRowColumnElimination()
{
  _applyRowElimination();
  if (m_has_row_column_elimination)
    _applyRowColumnElimination();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrDoFLinearSystemImpl::
_applyRowColumnElimination()
{
  IItemFamily* dof_family = dofFamily();

  auto nb_dof = dof_family->nbItem();

  RunQueue queue = makeQueue(runner());
  auto command = makeCommand(queue);

  auto in_elimination_info = Accelerator::viewIn(command, getEliminationInfo());
  auto in_elimination_value = Accelerator::viewIn(command, getEliminationValue());

  auto in_out_rhs_variable = Accelerator::viewInOut(command, rhsVariable());
  auto csr_view = m_csr_view;
  command << RUNCOMMAND_LOOP1(iter, nb_dof)
  {
    auto [row_index] = iter();
    DoFLocalId dof_row(row_index);
    auto row_elimination_info = in_elimination_info[dof_row];
    bool is_row_elimination = (row_elimination_info == ELIMINATE_ROW) || (row_elimination_info == ELIMINATE_ROW_COLUMN);
    for (CsrRowColumnIndex csr_index : csr_view.rowRange(dof_row)) {
      Int32 column_index = csr_view.column(csr_index);
      if (column_index > 0) {
        DoFLocalId dof_column(column_index);
        auto column_elimination_info = in_elimination_info[dof_column];
        bool is_column_elimination = (column_elimination_info == ELIMINATE_ROW) || (column_elimination_info == ELIMINATE_ROW_COLUMN);
        if (is_row_elimination || is_column_elimination) {
          csr_view.value(csr_index) = (dof_column == dof_row) ? 1.0 : 0.0;
        }
      }
    }
    if (is_row_elimination) {
      auto elimination_value = in_elimination_value[dof_row];
      in_out_rhs_variable[dof_row] = elimination_value;
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrDoFLinearSystemImpl::
_applyRowElimination()
{
  IItemFamily* dof_family = dofFamily();

  auto nb_dof = dof_family->nbItem();

  RunQueue queue = makeQueue(runner());
  auto command = makeCommand(queue);

  auto in_elimination_info = Accelerator::viewIn(command, getEliminationInfo());
  auto in_elimination_value = Accelerator::viewIn(command, getEliminationValue());

  auto in_out_rhs_variable = Accelerator::viewInOut(command, rhsVariable());
  auto csr_view = m_csr_view;
  command << RUNCOMMAND_LOOP1(iter, nb_dof)
  {
    auto [thread_id] = iter();
    DoFLocalId dof_id(thread_id);
    auto elimination_info = in_elimination_info[dof_id];
    if (elimination_info == ELIMINATE_ROW) {
      auto elimination_value = in_elimination_value[dof_id];
      for (CsrRowColumnIndex csr_index : csr_view.rowRange(dof_id))
        csr_view.value(csr_index) = (csr_view.column(csr_index) == dof_id) ? 1.0 : 0.0;
      in_out_rhs_variable[dof_id] = elimination_value;
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrDoFLinearSystemImpl::
_fillRowColumnEliminationInfos()
{
  // TODO: add a way to detect if user has added row/column elimination
  // to avoid doing this loop if this is not needed (because at the moment it
  // is only on CPU because OrderedRowColumnMap need to be ordered and so is
  // not available on GPU)
  OrderedRowColumnMap& rc_elimination_map = _rowColumnEliminationMap();
  rc_elimination_map.clear();
  DoFInfoListView item_list_view(dofFamily());

  auto& dof_elimination_info = getEliminationInfo();
  auto& dof_elimination_value = getEliminationValue();
  auto csr_view = m_csr_view;
  const Int32 nb_dof = dofFamily()->nbItem();

  m_has_row_column_elimination = false;
  for (Int32 i = 0; i < nb_dof; ++i) {
    DoFLocalId dof_row(i);
    Byte row_elimination_info = dof_elimination_info[dof_row];
    if (row_elimination_info == ELIMINATE_ROW_COLUMN) {
      m_has_row_column_elimination = true;
      break;
    }
  }
  if (!m_has_row_column_elimination)
    return;

  for (Int32 i = 0; i < nb_dof; ++i) {
    DoFLocalId dof_row(i);
    for (CsrRowColumnIndex csr_index : csr_view.rowRange(dof_row)) {
      Int32 col_index = csr_view.column(csr_index);
      if (col_index >= 0) {
        DoFLocalId dof_column(col_index);
        Byte row_elimination_info = dof_elimination_info[dof_row];
        Byte column_elimination_info = dof_elimination_info[dof_column];
        if (row_elimination_info == ELIMINATE_ROW_COLUMN || column_elimination_info == ELIMINATE_ROW_COLUMN) {
          Real value = csr_view.value(csr_index);
          rc_elimination_map[{ i, col_index }] = value;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
