// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* IDoFLinearSystemImpl.h                                      (C) 2022-2025 */
/*                                                                           */
/* Interface of the implementation of a DoFLinearSystem.                     */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_INTERNAL_IDOFLINEARSYSTEMIMPL_H
#define ARCANEFEM_FEMUTILS_INTERNAL_IDOFLINEARSYSTEMIMPL_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/VariableTypedef.h>
#include <arcane/core/ItemTypes.h>

#include "FemUtilsGlobal.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
class IDoFLinearSystemFactory;
class DoFLinearSystem;
class CsrFormatMatrixView;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \internal.
 * \brief Interface of the implementation of a DoFLinearSystem.
 */
class IDoFLinearSystemImpl
{
 public:

  virtual ~IDoFLinearSystemImpl() = default;

 public:

  virtual void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value) = 0;
  virtual void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value) = 0;
  virtual void eliminateRow(DoFLocalId row, Real value) = 0;
  virtual void eliminateRowColumn(DoFLocalId row, Real value) = 0;
  virtual void solve() = 0;
  virtual VariableDoFReal& solutionVariable() = 0;
  virtual VariableDoFReal& rhsVariable() = 0;
  virtual void setSolverCommandLineArguments(const CommandLineArguments& args) = 0;
  virtual void clearValues() = 0;
  virtual void setCSRValues(const CSRFormatView& csr_view) = 0;
  virtual CSRFormatView& getCSRValues() = 0;
  virtual VariableDoFReal& getForcedValue() = 0;
  virtual VariableDoFBool& getForcedInfo() = 0;
  virtual VariableDoFByte& getEliminationInfo() = 0;
  virtual VariableDoFReal& getEliminationValue() = 0;
  [[nodiscard]] virtual bool hasSetCSRValues() const = 0;
  virtual void setRunner(const Runner& r) = 0;
  [[nodiscard]] virtual Runner runner() const = 0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
