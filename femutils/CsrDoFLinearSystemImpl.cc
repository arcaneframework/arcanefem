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

#include "CsrDoFLinearSystemImpl.h"

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
: TraceAccessor(dof_family->traceMng())
, m_dof_family(dof_family)
, m_rhs_variable(VariableBuildInfo(dof_family, solver_name + "RHSVariable"))
, m_dof_variable(VariableBuildInfo(dof_family, solver_name + "SolutionVariable"))
, m_dof_forced_info(VariableBuildInfo(dof_family, solver_name + "DoFForcedInfo"))
, m_dof_forced_value(VariableBuildInfo(dof_family, solver_name + "DoFForcedValue"))
, m_dof_elimination_info(VariableBuildInfo(dof_family, solver_name + "DoFEliminationInfo"))
, m_dof_elimination_value(VariableBuildInfo(dof_family, solver_name + "DoFEliminationValue"))
{}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrDoFLinearSystemImpl::
_applyForcedValuesToLhs()
{
  auto nb_dof = m_dof_family->nbItem();

  RunQueue queue = makeQueue(m_runner);
  auto command = makeCommand(queue);

  auto in_out_forced_info = Accelerator::viewInOut(command, m_dof_forced_info);
  auto in_out_forced_value = Accelerator::viewInOut(command, m_dof_forced_value);

  auto csr_row_size = m_csr_view.nbRow();
  auto csr_columns_size = m_csr_view.nbColumn();
  auto in_csr_row = m_csr_view.rows();
  auto in_csr_columns = m_csr_view.columns();
  auto in_out_csr_values = m_csr_view.values();

  command << RUNCOMMAND_LOOP1(iter, nb_dof)
  {
    auto [dof_id] = iter();
    if (in_out_forced_info[(DoFLocalId)dof_id]) {
      auto begin = in_csr_row[dof_id];
      auto end = dof_id == csr_row_size - 1 ? csr_columns_size : in_csr_row[dof_id + 1];
      auto index = FemUtils::Gpu::Csr::findIndex(begin, end, dof_id, in_csr_columns);
      in_out_csr_values[index] = in_out_forced_value[(DoFLocalId)dof_id];
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrDoFLinearSystemImpl::
_applyRowElimination()
{
  auto nb_dof = m_dof_family->nbItem();

  RunQueue queue = makeQueue(m_runner);
  auto command = makeCommand(queue);

  auto in_elimination_info = Accelerator::viewIn(command, m_dof_elimination_info);
  auto in_elimination_value = Accelerator::viewIn(command, m_dof_elimination_value);

  auto in_out_rhs_variable = Accelerator::viewInOut(command, m_rhs_variable);
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

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
