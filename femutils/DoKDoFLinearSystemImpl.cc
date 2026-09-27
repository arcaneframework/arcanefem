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
/*---------------------------------------------------------------------------*/

#include "internal/DoKDoFLinearSystemImpl.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

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

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
