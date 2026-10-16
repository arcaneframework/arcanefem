// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemUtilsGlobal.h                                            (C) 2000-2026 */
/*                                                                           */
/* Defines types for FemUtils component.                                     */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_FEMUTILSGLOBAL_H
#define ARCANEFEM_FEMUTILS_FEMUTILSGLOBAL_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/UtilsTypes.h>

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

//! Vector of Real of size 4.
using Real4 = NumVector<Arcane::Real, 4>;

//! Vector of Real of size N.
template <int N> using RealVector = NumVector<Arcane::Real, N>;

//! Matrix of Real of size (Row, Column)
template <int Row, int Column> using RealMatrix = NumMatrix<Arcane::Real, Row, Column>;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief List of supported matrix format for linear systems.
 */
enum class eLinearSystemMatrixFormat
{
  //! Format Column Sparse Row
  Csr,
  //! Format Dictionary of Keys
  DoK
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

namespace ArcaneFemFunctions
{

class FemOperation2D;
class FemOperation3D;

//! To keep compatibility with existing code
using FeOperation2D = FemOperation2D;
//! To keep compatibility with existing code
using FeOperation3D = FemOperation3D;

} // namespace ArcaneFemFunctions

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
