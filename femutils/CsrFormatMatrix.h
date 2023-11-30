// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrFormatMatrix.cc                                     (C) 2022-2023 */
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

class CsrFormat : TraceAccessor
{
 public:

  CsrFormat(ISubDomain* sd)
  : TraceAccessor(sd->traceMng())
  {
    info() << "Creating CSR Matrix";
  }

  void initialize(IItemFamily* dof_family, Int32 nnz, Int32 nbRow)
  {
    m_matrix_row.resize(nbRow);
    m_matrix_column.resize(nnz);
    m_matrix_value.resize(nnz);
    m_matrix_row.fill(-1);
    m_matrix_column.fill(0);
    m_matrix_value.fill(0);
    m_dof_family = dof_family;
    m_last_value = 0;
    m_nnz = nnz;
    info() << "Filling CSR Matrix with zeros";
  }

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
    Int32 end;
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
  void translateToLinearSystem(DoFLinearSystem& linear_system)
  {
    for (Int32 i = 0; i < m_matrix_row.dim1Size(); i++) {
      for (Int32 j = m_matrix_row(i); (i + 1 < m_matrix_row.dim1Size() && j < m_matrix_row(i + 1)) || (i + 1 == m_matrix_row.dim1Size() && j < m_matrix_column.dim1Size()); j++)
        linear_system.matrixSetValue(DoFLocalId(i), DoFLocalId(m_matrix_column(j)), m_matrix_value(j));
    }
  }

  /**
 * @brief function to print the current content of the csr matrix 
 * 
 * @param fileName 
 * @param nonzero if set to true, print only the non zero values 
 */
  void
  printMatrix(std::string fileName)
  {
    ofstream file(fileName);
    file << "size :" << m_nnz << "\n";
    for (auto i = 0; i < m_matrix_row.dim1Size(); i++) {
      file << m_matrix_row(i) << " ";
      for (Int32 j = m_matrix_row(i) + 1; (i + 1 < m_matrix_row.dim1Size() && j < m_matrix_row(i + 1)) || (i + 1 == m_matrix_row.dim1Size() && j < m_matrix_column.dim1Size()); j++) {
        file << "  ";
      }
    }
    file << "\n";
    for (auto i = 0; i < m_nnz; i++) {
      file << m_matrix_column(i) << " ";
    }
    file << "\n";
    for (auto i = 0; i < m_nnz; i++) {
      file << m_matrix_value(i) << " ";
    }
    file << "\n";
    file.close();
  }

  // Warning : does not support empty row
  void
  setCoordinates(DoFLocalId row, DoFLocalId column)
  {
    if (m_matrix_row(row.localId()) == -1) {
      m_matrix_row(row.localId()) = m_last_value;
    }
    m_matrix_column(m_last_value) = column.localId();
    m_last_value++;
  }

 public:

  Int32 m_nnz;
  // To become parallelizable, have all the index
  // inside a queue that would gradually pop ?
  // or link the idnex to the index of the core ?
  Int32 m_last_value;
  NumArray<Int32, MDDim1> m_matrix_row;
  NumArray<Int32, MDDim1> m_matrix_column;
  NumArray<Real, MDDim1> m_matrix_value;
  IItemFamily* m_dof_family = nullptr;

  /*
  getValue return the Value at the (row, column) coordinates.
*/
  Int32 getValue(DoFLocalId row, DoFLocalId column)
  {
    return m_matrix_value(indexValue(row, column));
  }
};
} // namespace Arcane::FemUtils