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

 private:

  std::array<Arcane::Real, totalNbElement()> m_values = {};
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

