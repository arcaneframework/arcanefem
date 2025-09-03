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

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
