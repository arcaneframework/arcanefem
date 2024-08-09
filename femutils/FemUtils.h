// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemUtils.h                                                  (C) 2022-2024 */
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
#include <arcane/Parallel.h>
#include <arcane/IIOMng.h>
#include <arcane/CaseTable.h>

#include <arcane/utils/Real3.h>
#include <arcane/utils/Real3x3.h>

#include <array>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

struct Real4
{
    Arcane::Real data[4];
    
    Arcane::Real& operator[](std::size_t i) { return data[i]; }
    const Arcane::Real& operator[](std::size_t i) const { return data[i]; }
};

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

  //! Define the addition operator
  FixedMatrix<N, M> operator+(const FixedMatrix<N, M>& other) const
  {
    FixedMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = (*this)(i, j) + other(i, j);
      }
    }
    return result;
  }

  //! Define the subtraction operator
  FixedMatrix<N, M> operator-(const FixedMatrix<N, M>& other) const
  {
    FixedMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = (*this)(i, j) - other(i, j);
      }
    }
    return result;
  }

  //! Define the unary negation operator
  FixedMatrix<N, M> operator-() const
  {
    FixedMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = -(*this)(i, j);
      }
    }
    return result;
  }

  //! Scalar multiplication: FixedMatrix * scalar
  FixedMatrix<N, M> operator*(Real scalar) const
  {
    FixedMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = (*this)(i, j) * scalar;
      }
    }
    return result;
  }

  //! Friend function for scalar multiplication: scalar * FixedMatrix
  friend FixedMatrix<N, M> operator*(Real scalar, const FixedMatrix<N, M>& matrix)
  {
    FixedMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = scalar * matrix(i, j);
      }
    }
    return result;
  }

 private:

  std::array<Arcane::Real, totalNbElement()> m_values = {};
};

/*---------------------------------------------------------------------------*/
//  Outer product of two Real3 vectors to produce a FixedMatrix<3, 3>
/*---------------------------------------------------------------------------*/
inline FixedMatrix<3, 3> operator^(const Arcane::Real3& lhs, const Arcane::Real3& rhs)
{
  FixedMatrix<3, 3> result;
  for (Arcane::Int32 i = 0; i < 3; ++i) {
    for (Arcane::Int32 j = 0; j < 3; ++j) {
      result(i, j) = lhs[i] * rhs[j];
    }
  }
  return result;
}

/*---------------------------------------------------------------------------*/
//  Outer product of two Real4 vectors to produce a FixedMatrix<4, 4>
/*---------------------------------------------------------------------------*/
inline FixedMatrix<4, 4> operator^(const Real4& lhs, const Real4& rhs)
{
    FixedMatrix<4, 4> result;
    for (Arcane::Int32 i = 0; i < 4; ++i) {
        for (Arcane::Int32 j = 0; j < 4; ++j) {
            result(i, j) = lhs[i] * rhs[j];
        }
    }
    return result;
}

/*---------------------------------------------------------------------------*/
// Define the conversion from Real3x3 to FixedMatrix<3, 3>
/*---------------------------------------------------------------------------*/
inline FixedMatrix<3, 3> convertReal3x3ToFixedMatrix(const Arcane::Real3x3& rhs)
{
  FixedMatrix<3, 3> result;
  for (Arcane::Int32 i = 0; i < 3; ++i) {
    for (Arcane::Int32 j = 0; j < 3; ++j) {
      result(i, j) = rhs(i, j);
    }
  }
  return result;
}

/*---------------------------------------------------------------------------*/
// Overload operator+ to handle FixedMatrix<3, 3> + Real3x3
/*---------------------------------------------------------------------------*/
inline FixedMatrix<3, 3> operator+(const FixedMatrix<3, 3>& lhs, const Arcane::Real3x3& rhs)
{
  FixedMatrix<3, 3> converted_rhs = convertReal3x3ToFixedMatrix(rhs);
  return lhs + converted_rhs;
}

/*---------------------------------------------------------------------------*/
// Overload operator+ to handle Real3x3 + FixedMatrix<3, 3>
/*---------------------------------------------------------------------------*/
inline FixedMatrix<3, 3> operator+(const Arcane::Real3x3& lhs, const FixedMatrix<3, 3>& rhs)
{
  FixedMatrix<3, 3> converted_lhs = convertReal3x3ToFixedMatrix(lhs);
  return converted_lhs + rhs;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <int N, int M> inline FixedMatrix<N, N>
matrixAddition(const FixedMatrix<N, M>& a, const FixedMatrix<M, N>& b)
{
  using namespace Arcane;
  FixedMatrix<N, N> new_matrix;

  for (Int32 i = 0; i < N; ++i) {
    for (Int32 j = 0; j < N; ++j) {
      new_matrix(i, j) = a(i, j) + b(i, j);
    }
  }
  return new_matrix;
}

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
 *
 * Only values greater whose absolute value is greater than \a min_value are
 * checked.
 */
extern "C++" void
checkNodeResultFile(ITraceMng* tm, const String& filename,
                    const VariableNodeReal& node_values, double epsilon,
                    double min_value = 0.0);

/*!
 * \brief Check the values of the variable against a reference file.
 *
 * The reference file \a filename is an ASCII file containing \a N \a Real2
 * values where \a N is the number of nodes in the mesh (so there will be
 * \a N*2 real in the file). The file should have the following format:
 * uid1 x1 y1 uid2 x2 y2 ... uidN xN yN.
 *
 * For each item the current value of \a node_values is compared to the
 * reference value and if the relative difference is greater than \a epsilon
 * this is an error.
 *
 * Only values greater whose absolute value is greater than \a min_value are
 * checked.
 */
extern "C++" void
checkNodeResultFile(ITraceMng* tm, const String& filename,
                    const VariableNodeReal2& node_values, double epsilon,
                    double min_value = 0.0);

/*!
 * \brief Check the values of the variable against a reference file.
 *
 * The reference file \a filename is an ASCII file containing \a N \a Real3
 * values where \a N is the number of nodes in the mesh (so there will be
 * \a N*3 real in the file). The file should have the following format:
 * uid1 x1 y1 z1 uid2 x2 y2 z2 ... uidN xN yN zN.
 *
 * For each item the current value of \a node_values is compared to the
 * reference value and if the relative difference is greater than \a epsilon
 * this is an error.
 *
 * Only values greater whose absolute value is greater than \a min_value are
 * checked.
 */
extern "C++" void
checkNodeResultFile(ITraceMng* tm, const String& filename,
                    const VariableNodeReal3& node_values, double epsilon,
                    double min_value = 0.0);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Sample to read value from a file and create an associated CaseTable.
 *
 * The file should contains 3 values for each time step (so the number of
 * values should be a multiple of 4).
 */
extern "C++" CaseTable*
readFileAsCaseTable(IParallelMng* pm, const String& filename, const Int32& ndim);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
