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

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
