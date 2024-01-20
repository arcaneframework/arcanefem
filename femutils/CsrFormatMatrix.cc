// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrFormatMatrix.cc                                          (C) 2022-2024 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/NumArray.h>

#include <arcane/VariableTypes.h>
#include <arcane/IItemFamily.h>

#include "CsrFormatMatrix.h"

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrFormat::
initialize(IItemFamily* dof_family, Int32 nnz, Int32 nbRow)
{
  info() << "Initialize CsrFormat: nb_non_zero=" << nnz << " nb_row=" << nbRow;

  m_matrix_row.resize(nbRow);
  m_matrix_column.resize(nnz);
  m_matrix_value.resize(nnz);
  m_matrix_row.fill(-1);
  m_matrix_column.fill(-1);
  m_matrix_value.fill(0);
  m_matrix_rows_nb_column.resize(nbRow);
  m_matrix_rows_nb_column.fill(0);
  m_dof_family = dof_family;
  m_last_value = 0;
  m_nnz = nnz;
  info() << "Filling CSR Matrix with zeros";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CsrFormat::
translateToLinearSystem(DoFLinearSystem& linear_system)
{
  info() << "TranslateToLinearSystem this=" << this;
  bool do_set_csr = linear_system.hasSetCSRValues();
  // When using CSR format, we need to know the number of non zero values for
  // each row.
  // NOTE: it should be possible to compute that in setCoordinates().
  // and this value is constant if the structure of the matrix do not change
  // so we can store these values instead of recomputing them.
  if (do_set_csr){
    m_matrix_rows_nb_column.resize(m_matrix_row.extent0());
    m_matrix_rows_nb_column.fill(0);
  }
  Int32 nb_row = m_matrix_row.dim1Size();
  for (Int32 i = 0; i < nb_row; i++) {
    if (((i + 1) < nb_row) && (m_matrix_row(i) == m_matrix_row(i + 1)))
      continue;
    for (Int32 j = m_matrix_row(i); ((i + 1) < nb_row && j < m_matrix_row(i + 1)) || ((i + 1) == nb_row && j < m_matrix_column.dim1Size()); j++) {
      if (do_set_csr){
        ++m_matrix_rows_nb_column[i];
        continue;
      }
      if (DoFLocalId(m_matrix_column(j)).isNull())
        continue;
      //info() << "Add: (" << i << ", " << m_matrix_column(j) << " v=" << m_matrix_value(j);
      linear_system.matrixAddValue(DoFLocalId(i), DoFLocalId(m_matrix_column(j)), m_matrix_value(j));
    }
  }

  if (do_set_csr){
    CSRFormatView csr_view(m_matrix_row.to1DSpan(),m_matrix_rows_nb_column.to1DSpan(),
                           m_matrix_column.to1DSpan(),m_matrix_value.to1DSpan());
    linear_system.setCSRValues(csr_view);
  }
}

void CsrFormat::
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

} // namespace Arcane::FemUtils
