// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemUtils.h                                                  (C) 2022-2022 */
/*                                                                           */
/* Classes utilisateurs pour FEM.                                            */
/*---------------------------------------------------------------------------*/
#ifndef FEMTEST_FEMUTILS_H
#define FEMTEST_FEMUTILS_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/ArcaneTypes.h>

#include <array>
#include <iostream>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Matrice NxM de taille fixe.
 */
template <int N, int M>
class FixedMatrix
{
  using ThatClass = FixedMatrix<N, M>;

 public:

  static constexpr Arcane::Int32 totalNbElement() { return N * M; }

 public:

  Arcane::Real& operator()(Arcane::Int32 i, Arcane::Int32 j)
  {
    ARCANE_CHECK_AT(i,N);
    ARCANE_CHECK_AT(j,M);
    return m_values[i * M + j];
  }

  Arcane::Real operator()(Arcane::Int32 i, Arcane::Int32 j) const
  {
    ARCANE_CHECK_AT(i,N);
    ARCANE_CHECK_AT(j,M);
    return m_values[i * M + j];
  }

 public:

  //! Multiply all the components by \a v
  void multInPlace(Arcane::Real v)
  {
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] *= v;
  }

  //! Dump matrix values
  void dump(std::ostream& o) const
  {
    const ThatClass& values = *this;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      o << "[ ";
      for (Arcane::Int32 j = 0; j < M; ++j) {
        if (j != 0)
          o << ' ';
        o << values(i, j);
      }
      o << "]\n";
    }
  }

 private:

  std::array<Arcane::Real, totalNbElement()> m_values = {};
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

