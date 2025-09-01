// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoFLinearSystemImplBase.h                                   (C) 2022-2025 */
/*                                                                           */
/* Base class of implementation of IDoFLinearSystemImpl.                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_INTERNAL_DOFLINEARSYSTEMIMPLBASE_H
#define ARCANEFEM_FEMUTILS_INTERNAL_DOFLINEARSYSTEMIMPLBASE_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/VariableTypes.h>
#include <arcane/accelerator/core/Runner.h>

#include "internal/IDoFLinearSystemImpl.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \internal
 * \brief Base class of implementation of IDoFLinearSystemImpl.
 *
 * This class will handle parts common to all implementation of IDoFLinearSystemImpl.
 */
class DoFLinearSystemImplBase
: public TraceAccessor
, public IDoFLinearSystemImpl
{
 public:

  static constexpr Byte ELIMINATE_NONE = 0;
  static constexpr Byte ELIMINATE_ROW = 1;
  static constexpr Byte ELIMINATE_ROW_COLUMN = 2;

 public:

  DoFLinearSystemImplBase(IItemFamily* dof_family, const String& solver_name);

 public:

  void clearValues() override;

  VariableDoFReal& solutionVariable() final  {    return m_solution_variable;  }

  VariableDoFReal& rhsVariable() final  {    return m_rhs_variable;  }

  VariableDoFBool& getForcedInfo() final { return m_dof_forced_info; }
  VariableDoFReal& getForcedValue() final { return m_dof_forced_value; }
  VariableDoFByte& getEliminationInfo() final { return m_dof_elimination_info; }
  VariableDoFReal& getEliminationValue() final { return m_dof_elimination_value; }

  void setRunner(const Runner& r) override { m_runner = r; }
  Runner runner() const override { return m_runner; }

  IItemFamily* dofFamily() const { return m_dof_family; }

 private:


  IItemFamily* m_dof_family = nullptr;
  Runner m_runner;

  VariableDoFReal m_rhs_variable;
  VariableDoFReal m_solution_variable;
  VariableDoFBool m_dof_forced_info;
  VariableDoFReal m_dof_forced_value;
  VariableDoFByte m_dof_elimination_info;
  VariableDoFReal m_dof_elimination_value;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
