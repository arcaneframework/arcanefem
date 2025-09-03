// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* OrderedRowColumndMap.h                                      (C) 2022-2025 */
/*                                                                           */
/* Ordered map to keep a set a values indexed by (row,column).               */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_INTERNAL_ORDEREDROWCOLUMNMAP_H
#define ARCANEFEM_FEMUTILS_INTERNAL_ORDEREDROWCOLUMNMAP_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/ItemTypes.h>

#include "FemUtilsGlobal.h"

#include <map>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Ordered map to keep a set a values indexed by (row,column).
 */
class OrderedRowColumnMap
{
 public:

  //! Index in a RowColumnMap
  struct RowColumn
  {
    Int32 row_id = 0;
    Int32 column_id = 0;
    friend bool operator==(RowColumn rc1, RowColumn rc2)
    {
      if (rc1.row_id != rc2.row_id)
        return false;
      return rc1.column_id == rc2.column_id;
    }
    friend bool operator<(RowColumn rc1, RowColumn rc2)
    {
      if (rc1.row_id == rc2.row_id)
        return rc1.column_id < rc2.column_id;
      return rc1.row_id < rc2.row_id;
    }
  };

 private:

  using MapType = std::map<RowColumn, Real>;

 public:

  using iterator = typename MapType::iterator;
  using const_iterator = typename MapType::const_iterator;

 public:

  void addValue(RowColumn rc, Real value)
  {
    auto x = m_values_map.find(rc);
    if (x == m_values_map.end())
      m_values_map.insert(std::make_pair(rc, value));
    else
      x->second += value;
  }

  Real& operator[](RowColumn rc)
  {
    return m_values_map[rc];
  }

  void clear() { m_values_map.clear(); }
  iterator begin() { return m_values_map.begin(); }
  iterator end() { return m_values_map.end(); }
  const_iterator begin() const { return m_values_map.begin(); }
  const_iterator end() const { return m_values_map.end(); }
  iterator find(RowColumn rc) { return m_values_map.find(rc); }
  const_iterator find(RowColumn rc) const { return m_values_map.find(rc); }

 private:

  MapType m_values_map;
};


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
