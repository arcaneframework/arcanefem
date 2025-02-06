// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* IDoFLinearSystemFactory.h                                   (C) 2022-2025 */
/*                                                                           */
/* Interface to a factory to build a linear system implementation.           */
/*---------------------------------------------------------------------------*/
#ifndef FEMTEST_IDOFLINEARSYSTEMFACTORY_H
#define FEMTEST_IDOFLINEARSYSTEMFACTORY_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/ItemTypes.h>
#include <arcane/VariableTypedef.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class DoFLinearSystemImpl;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class IDoFLinearSystemFactory
{
 public:

  virtual ~IDoFLinearSystemFactory() = default;

 public:

  virtual DoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name) = 0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
