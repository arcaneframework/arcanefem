// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* IDoFLinearSystemFactory.h                                   (C) 2000-2026 */
/*                                                                           */
/* Interface to a factory to build a linear system implementation.           */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_IDOFLINEARSYSTEMFACTORY_H
#define ARCANEFEM_FEMUTILS_IDOFLINEARSYSTEMFACTORY_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/ItemTypes.h>

#include "FemUtilsGlobal.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Interface to a factory to build a linear system implementation.
 */
class IDoFLinearSystemFactory
{
 public:

  virtual ~IDoFLinearSystemFactory() = default;

 public:

  //! Whether the FEM module should provide an AMG near-null-space basis.
  virtual bool amgNearNullSpace() { return false; }

  //! Create an instance using the default of the factory implementation format for matrix
  virtual IDoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name) = 0;

  //! Create an instance using the format 'matrix_format' for the matrix
  virtual IDoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name,
                 eLinearSystemMatrixFormat matrix_format);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
