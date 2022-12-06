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
#include <arcane/utils/PlatformUtils.h>
#include <arcane/utils/ValueConvert.h>

#include <arcane/VariableTypes.h>
#include <arcane/IItemFamily.h>

#include <map>

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

void checkNodeResultFile(ITraceMng* tm, const String& filename,
                         const VariableNodeReal& node_values, double epsilon)
{
  std::cout << "CheckNodeResultFile filename=" << filename << "\n";
  if (filename.empty())
    ARCANE_FATAL("Invalid empty filename");
  IItemFamily* node_family = node_values.variable()->itemFamily();
  if (!node_family)
    ARCANE_FATAL("Variable '{0}' is not allocated", node_values.name());

  std::map<Int64, double> item_reference_values;
  {
    std::ifstream sbuf(filename.localstr());
    double read_value = 0.0;
    Int64 read_uid = 0;
    if (!sbuf.eof())
      sbuf >> ws;
    while (!sbuf.eof()) {
      sbuf >> read_uid >> ws >> read_value;
      if (sbuf.fail() || sbuf.bad())
        ARCANE_FATAL("Error during parsing of file '{0}'", filename);
      item_reference_values.insert(std::make_pair(read_uid, read_value));
      sbuf >> ws;
    }
  }

  std::cout << "NB_Values=" << item_reference_values.size() << "\n";

  // Get Max UID
  Int64 max_uid = 0;
  ENUMERATE_ (Node, inode, node_family->allItems()) {
    Node node = *inode;
    if (node.uniqueId() > max_uid)
      max_uid = node.uniqueId();
  }

  std::map<Int64, double> item_current_values;

  ENUMERATE_ (Node, inode, node_family->allItems()) {
    Node node = *inode;
    item_current_values[node.uniqueId()] = node_values[node];
  }

  Int32 nb_ref_value = item_reference_values.size();
  Int32 nb_current_value = item_current_values.size();
  if (nb_ref_value != nb_current_value)
    ARCANE_FATAL("Can not compare files because there is not the same number of values nb_ref={0} nb_current={1}",
                 nb_ref_value, nb_current_value);

  Int64 nb_error = 0;
  for (const auto& x : item_current_values) {
    Int64 uid = x.first;
    Real v = x.second;
    auto x_ref = item_reference_values.find(uid);
    if (x_ref != item_reference_values.end()) {
      Real ref_v = x_ref->second;
      if (!TypeEqualT<double>::isNearlyEqualWithEpsilon(ref_v, v, epsilon)) {
        ++nb_error;
        if (nb_error < 15)
          std::cout << String::format("ERROR: ref={0} v={1} diff={2}", ref_v, v, ref_v - v) << "\n";
      }
    }
  }
  if (nb_error > 0)
    ARCANE_FATAL("Error checking values nb_error={0}", nb_error);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
