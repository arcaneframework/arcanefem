// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemUtilsGlobal.h                                            (C) 2022-2025 */
/*                                                                           */
/* Defines types for FemUtils component.                                     */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_FEMUTILSGLOBAL_H
#define ARCANEFEM_FEMUTILS_FEMUTILSGLOBAL_H

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemDoFsOnNodes.h"
#include <arcane/utils/ArcaneGlobal.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class CsrRowColumnIterator;
class CsrFormatMatrixView;
class CsrRow;
class CsrFormat;
class CsrRowColumnIndex;
class IDoFLinearSystemFactory;
class DoFLinearSystem;
class IDoFLinearSystemImpl;
class FemDoFsOnNodes;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

//! Old name to keep compatibility with existing code.
using CSRFormatView = CsrFormatMatrixView;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief List of matrix elimination type
 */
enum class MatrixEliminationType : Byte
{
  //! No elimination
  None = 0,
  //! Row elimination
  Row = 1,
  //! RowColumn elimination
  RowColumn = 2,
};
static constexpr Byte ELIMINATE_NONE = static_cast<Byte>(MatrixEliminationType::None);
static constexpr Byte ELIMINATE_ROW = static_cast<Byte>(MatrixEliminationType::Row);
static constexpr Byte ELIMINATE_ROW_COLUMN = static_cast<Byte>(MatrixEliminationType::RowColumn);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
