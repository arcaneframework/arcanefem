// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* LinkedCsrFormatMatrix.h                                     (C) 2022-2023 */
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

class LinkedCsrFormat : TraceAccessor
{
 public:

  LinkedCsrFormat(ISubDomain* sd)
  : TraceAccessor(sd->traceMng())
  {
    info() << "Creating LinkedCSR Matrix";
  }

 public:

  void initialize(IItemFamily* dof_family, Int32 nb_nodes)
  {
    m_row_to_col.resize(nb_nodes);
    m_row_to_val.resize(nb_nodes);
    m_dof_family = dof_family;
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
    if (value == 0.0)
      return;
    m_row_to_val(row.localId())(indexValue(row, column)) += value;
  }

  /**
 * @brief 
 * 
 * @param row 
 * @param column 
 * @return Int64
 */
  Int64 indexValue(DoFLocalId row, DoFLocalId column)
  {
    for (Int64 i = 0; i < m_row_to_col(row.localId()).capacity(); i++) {
      if (m_row_to_col(row.localId())(i) == column.localId())
        return i;
    }
    return -1;
  }

  /**
 * @brief  Translate to Arcane linear system
 * 
 * @param linear_system 
 *
  */
  void translateToLinearSystem(DoFLinearSystem& linear_system)
  {
    for (Int32 i = 0; i < m_row_to_col.capacity(); i++) {
      for (Int32 j = 0; m_row_to_col(i).dim1Size(); j++)
        linear_system.matrixAddValue(DoFLocalId(i), DoFLocalId(m_row_to_col(i)(j)), m_row_to_val(i)(j));
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
    file << "size :" << m_row_to_col.capacity() << "\n";
    for (auto i = 0; i < m_row_to_col.capacity(); i++) {
      file << i << " ";
      for (auto j = 1; j < m_row_to_col(i).dim1Size(); j++) {
        file << "  ";
      }
    }
    file << "\n";
    for (auto i = 0; i < m_row_to_col.capacity(); i++) {
      for (auto j = 0; j < m_row_to_col(i).dim1Size(); j++) {
        file << m_row_to_col(i)(j) << " ";
      }
    }
    file << "\n";
    for (auto i = 0; i < m_row_to_val.capacity(); i++) {
      for (auto j = 0; j < m_row_to_val(i).dim1Size(); j++) {
        file << m_row_to_val(i)(j) << " ";
      }
    }
    file.close();
  }

  SharedArray<NumArray<Int32, MDDim1>> m_row_to_col;
  SharedArray<NumArray<Real, MDDim1>> m_row_to_val;
  IItemFamily* m_dof_family = nullptr;
};
} // namespace Arcane::FemUtils