// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemUtils.h                                                  (C) 2022-2022 */
/*                                                                           */
/* Utilitary classes for FEM.                                                */
/*---------------------------------------------------------------------------*/
#ifndef FEMTEST_FEMUTILS_H
#define FEMTEST_FEMUTILS_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/ArcaneTypes.h>
#include <arcane/utils/MDSpan.h>
#include <arcane/matvec/Matrix.h>
#include <arcane/VariableTypedef.h>

#include <array>
#include <iostream>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

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
    ARCANE_CHECK_AT(i, N);
    ARCANE_CHECK_AT(j, M);
    return m_values[i * M + j];
  }

  Arcane::Real operator()(Arcane::Int32 i, Arcane::Int32 j) const
  {
    ARCANE_CHECK_AT(i, N);
    ARCANE_CHECK_AT(j, M);
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

template <int N, int M> inline FixedMatrix<N, N>
matrixMultiplication(const FixedMatrix<N, M>& a, const FixedMatrix<M, N>& b)
{
  using namespace Arcane;
  FixedMatrix<N, N> new_matrix;

  for (Int32 i = 0; i < N; ++i) {
    for (Int32 j = 0; j < N; ++j) {
      Real x = 0.0;
      for (Int32 k = 0; k < M; ++k) {
        x += a(i, k) * b(k, j);
      }
      new_matrix(i, j) += x;
    }
  }
  return new_matrix;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <int N, int M> inline FixedMatrix<M, N>
matrixTranspose(const FixedMatrix<N, M>& a)
{
  using namespace Arcane;

  FixedMatrix<M, N> t_matrix;
  for (Int32 i = 0; i < N; ++i) {
    for (Int32 j = 0; j < M; ++j) {
      t_matrix(j, i) = a(i, j);
    }
  }
  return t_matrix;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Convert a dense matrix to an Arcane sequential CSR Matrix.
 */
extern "C++" void
_convertNumArrayToCSRMatrix(Arcane::MatVec::Matrix& out_matrix,
                            Arcane::MDSpan<const Arcane::Real, Arcane::MDDim2> in_matrix);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Check the values of the variable against a reference file.
 *
 * The reference file \a filename is an ASCII file containing N values where
 * \a N is the number of nodes in the mesh.
 *
 * For each item the current value of \a node_values is compared to the
 * reference value and if the relative difference is greater than \a epsilon
 * this is an error.
 */
extern "C++" void
checkNodeResultFile(ITraceMng* tm, const String& filename,
                    const VariableNodeReal& node_values, double epsilon);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
