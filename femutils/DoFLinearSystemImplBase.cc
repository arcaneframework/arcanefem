// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoFDoFLinearSystemImplBase.cc                                (C) 2022-2025 */
/*                                                                           */
/* Base class of implementation of IDoFLinearSystemImpl.                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "internal/DoFLinearSystemImplBase.h"

#include <arcane/core/IItemFamily.h>
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

DoFLinearSystemImplBase::
DoFLinearSystemImplBase(IItemFamily* dof_family, const String& solver_name)
: TraceAccessor(dof_family->traceMng())
, m_dof_family(dof_family)
, m_rhs_variable(VariableBuildInfo(dof_family, solver_name + "RHSVariable"))
, m_solution_variable(VariableBuildInfo(dof_family, solver_name + "SolutionVariable"))
, m_dof_forced_info(VariableBuildInfo(dof_family, solver_name + "DoFForcedInfo"))
, m_dof_forced_value(VariableBuildInfo(dof_family, solver_name + "DoFForcedValue"))
, m_dof_elimination_info(VariableBuildInfo(dof_family, solver_name + "DoFEliminationInfo"))
, m_dof_elimination_value(VariableBuildInfo(dof_family, solver_name + "DoFEliminationValue"))
{}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystemImplBase::
clearValues()
{
  m_dof_forced_info.fill(false);
  m_dof_elimination_info.fill(ELIMINATE_NONE);
  m_dof_elimination_value.fill(0.0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystemImplBase::
_applyRowColumnEliminationToRHS(bool is_verbose)
{
  DoFInfoListView item_list_view(m_dof_family);

  auto& dof_elimination_info = getEliminationInfo();
  auto& dof_elimination_value = getEliminationValue();
  auto& rhs_variable = rhsVariable();

  const bool do_print_filling = true;
  for (const auto& rc_value : m_row_column_elimination_map) {
    auto rc = rc_value.first;
    Real matrix_value = rc_value.second;
    DoF dof_row = item_list_view[rc.row_id];
    DoF dof_column = item_list_view[rc.column_id];
    if (dof_row == dof_column)
      continue;
    if (!dof_column.isOwn())
      continue;
    Byte row_elimination_info = dof_elimination_info[dof_row];
    Real elimination_value = dof_elimination_value[dof_row];
    // Subtract the value of RHS vector for current column.
    if (row_elimination_info == ELIMINATE_ROW_COLUMN) {
      Real v = rhs_variable[dof_column];
      rhs_variable[dof_column] = v - matrix_value * elimination_value;
      if (is_verbose)
        info() << "EliminateRowColumn (" << std::setw(4) << rc.row_id
               << "," << std::setw(4) << rc.column_id << ")"
               << " elimination_value=" << std::setw(25) << elimination_value
               << "  old_rhs=" << std::setw(25) << v
               << "  new_rhs=" << std::setw(25) << rhs_variable[dof_column];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
