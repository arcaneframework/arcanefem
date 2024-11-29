// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrFormatMatrix.h                                           (C) 2022-2024 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/NumArray.h>

#include <arcane/VariableTypes.h>
#include <arcane/IItemFamily.h>

#include <arcane/aleph/AlephTypesSolver.h>
#include <arcane/aleph/Aleph.h>

#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "arcane_version.h"

#include <iostream>
#include <fstream>

namespace Arcane::FemUtils
{
using namespace Arcane;

class CsrFormat
: public TraceAccessor
{
 public:

  explicit CsrFormat(ITraceMng* tm)
  : TraceAccessor(tm)
  {
  }

 public:

  void initialize(IItemFamily* dof_family, Int32 nnz, Int32 nbRow, RunQueue& queue);

  /**
   * @brief
   *
   * @param row
   * @param column
   * @param value
   */
  void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value)
  {
    if (row.isNull())
      ARCANE_FATAL("Row is null");
    if (column.isNull())
      ARCANE_FATAL("Column is null");
    if (value == 0.0)
      return;
    m_matrix_value(indexValue(row, column)) += value;
  }

  Int32 indexValue(DoFLocalId row, DoFLocalId column)
  {
    Int32 begin = m_matrix_row(row.localId());
    Int32 end = 0;
    if (row.localId() == m_matrix_row.extent0() - 1) {

      end = m_matrix_column.extent0();
    }
    else {

      end = m_matrix_row(row + 1);
    }
    for (Int32 i = begin; i < end; i++) {
      if (m_matrix_column(i) == column.localId()) {
        return i;
      }
    }
    return -1;
  }

  /**
   * @brief
   *
   * @param linear_system
   */
  void translateToLinearSystem(DoFLinearSystem& linear_system, const RunQueue& queue);

  /**
   * @brief function to print the current content of the csr matrix
   *
   * @param fileName
   * @param nonzero if set to true, print only the non zero values
   */
  void printMatrix(std::string fileName);

  // Warning : does not support empty row (or does it ?)
  void setCoordinates(DoFLocalId row, DoFLocalId column)
  {
    Int32 row_lid = row.localId();
    if (m_matrix_row(row_lid) == -1) {
      m_matrix_row(row_lid) = m_last_value;
    }
    m_matrix_column(m_last_value) = column.localId();
    m_last_value++;
  }

  void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value)
  {
    m_matrix_value(indexValue(row, column)) = value;
  }

 public:

  Int32 m_nnz = 0;
  // To become parallelizable, have all the index
  // inside a queue that would gradually pop ?
  // or link the idnex to the index of the core ?
  Int32 m_last_value = 0;
  NumArray<Int32, MDDim1> m_matrix_row;
  NumArray<Int32, MDDim1> m_matrix_column;
  NumArray<Real, MDDim1> m_matrix_value;
  //! Nombre de colonnes de chaque lignes.
  NumArray<Int32, MDDim1> m_matrix_rows_nb_column;
  IItemFamily* m_dof_family = nullptr;

  //! Return the Value at the (row, column) coordinates.
  Int32 getValue(DoFLocalId row, DoFLocalId column)
  {
    return m_matrix_value(indexValue(row, column));
  }
};

} // namespace Arcane::FemUtils
