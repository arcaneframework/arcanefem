// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* OrderedRowColumnMap.h                                       (C) 2000-2026 */
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

  using iterator = MapType::iterator;
  using const_iterator = MapType::const_iterator;

 public:

  //! Add the value at index \a rc.
  void addValue(RowColumn rc, Real value)
  {
    auto x = m_values_map.find(rc);
    if (x == m_values_map.end())
      m_values_map.insert(std::make_pair(rc, value));
    else
      x->second += value;
  }
  //! Set the value at index \a rc, replacing current value if it exists
  void setValue(RowColumn rc, Real value)
  {
    auto x = m_values_map.find(rc);
    if (x == m_values_map.end())
      m_values_map.insert(std::make_pair(rc, value));
    else
      x->second = value;
  }
  /*!
   * \brief Read-only value of the matrix at index \a rc.
   *
   * Return zero if there is no value at the current index.
   */
  const Real operator[](RowColumn rc)
  {
    auto x = m_values_map.find(rc);
    if (x == m_values_map.end())
      return {};
    return x->second;
  }

  void clear() { m_values_map.clear(); }
  iterator begin() { return m_values_map.begin(); }
  iterator end() { return m_values_map.end(); }
  [[nodiscard]] const_iterator begin() const { return m_values_map.begin(); }
  [[nodiscard]] const_iterator end() const { return m_values_map.end(); }
  iterator find(RowColumn rc) { return m_values_map.find(rc); }
  [[nodiscard]] const_iterator find(RowColumn rc) const { return m_values_map.find(rc); }
  [[nodiscard]] bool contains(RowColumn rc) const { return find(rc) != end(); }
  [[nodiscard]] Int32 size() const { return static_cast<Int32>(m_values_map.size()); }

 private:

  MapType m_values_map;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
