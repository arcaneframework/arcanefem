// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemUtils.cc                                                 (C) 2022-2022 */
/*                                                                           */
/* Utilitary classes for FEM.                                                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemUtils.h"

#include <arcane/utils/FatalErrorException.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
using namespace Arcane::MatVec;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void _convertNumArrayToCSRMatrix(Matrix& out_matrix, MDSpan<const Real, MDDim2> in_matrix)
{
  const Int32 matrix_size = in_matrix.extent0();
  if (matrix_size != in_matrix.extent1())
    ARCANE_FATAL("Not a square matrix ({0},{1})", matrix_size, in_matrix.extent1());
  UniqueArray<Int32> nb_non_zero_in_row(matrix_size);

  Int32 nb_non_zero = 0;
  for (Int32 i = 0; i < matrix_size; ++i) {
    Int32 row_non_zero = 0;
    for (Int32 j = 0; j < matrix_size; ++j)
      if (in_matrix(i, j) != 0)
        ++row_non_zero;
    nb_non_zero_in_row[i] = row_non_zero;
    nb_non_zero += row_non_zero;
  }

  UniqueArray<Int32> columns(nb_non_zero);
  UniqueArray<Real> values(nb_non_zero);

  {
    Int32 fill_index = 0;
    for (Int32 i = 0; i < matrix_size; ++i) {
      for (Int32 j = 0; j < matrix_size; ++j) {
        Real v = in_matrix(i, j);
        if (v != 0) {
          columns[fill_index] = j;
          values[fill_index] = v;
          ++fill_index;
        }
      }
    }
  }

  out_matrix.setRowsSize(nb_non_zero_in_row);
  out_matrix.setValues(columns, values);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
