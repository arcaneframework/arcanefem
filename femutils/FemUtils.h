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

#include <arccore/base/ArccoreGlobal.h>
#include <array>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

struct Real4
{
  Arcane::Real data[4];

  ARCCORE_HOST_DEVICE Arcane::Real& operator[](std::size_t i) { return data[i]; }
  ARCCORE_HOST_DEVICE const Arcane::Real& operator[](std::size_t i) const { return data[i]; }
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

  ARCCORE_HOST_DEVICE Arcane::Real& operator()(Arcane::Int32 i, Arcane::Int32 j)
  {
    ARCANE_CHECK_AT(i, N);
    ARCANE_CHECK_AT(j, M);
    return m_values[i * M + j];
  }

  ARCCORE_HOST_DEVICE Arcane::Real operator()(Arcane::Int32 i, Arcane::Int32 j) const
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
  ARCCORE_HOST_DEVICE FixedMatrix<N, M> operator+(const FixedMatrix<N, M>& other) const
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
  ARCCORE_HOST_DEVICE friend FixedMatrix<N, M> operator*(Real scalar, const FixedMatrix<N, M>& matrix)
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
ARCCORE_HOST_DEVICE inline FixedMatrix<3, 3> operator^(const Arcane::Real3& lhs, const Arcane::Real3& rhs)
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
ARCCORE_HOST_DEVICE inline FixedMatrix<4, 4> operator^(const Real4& lhs, const Real4& rhs)
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
 * \brief Vector N de taille fixe.
 */
template <int N>
class FixedVector
{
  using ThatClass = FixedVector<N>;

 public:

  static constexpr Arcane::Int32 totalNbElement() { return N; }

 public:

  Arcane::Real& operator()(Arcane::Int32 i)
  {
    ARCANE_CHECK_AT(i, N);
    return m_values[i];
  }

  ARCCORE_HOST_DEVICE Arcane::Real operator()(Arcane::Int32 i) const
  {
    ARCANE_CHECK_AT(i, N);
    return m_values[i];
  }

 public:

  //! Multiply all the components by \a v
  void multInPlace(Arcane::Real v)
  {
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] *= v;
  }

  //! Add \a v to all the components
  void addInPlace(Arcane::Real v)
  {
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] += v;
  }

  //! Dump values
  void dump(std::ostream& o) const
  {
    const ThatClass& values = *this;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      o << "[ ";
      o << values(i);
      o << "]\n";
    }
  }

  //! Set this vector equal to b
  void setEqualTo(const FixedVector<N>& b)
  {
    ARCANE_CHECK_AT(totalNbElement(), N);
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] = b[i];
  }

  //! Add b to this vector
  void add(const FixedVector<N>& b)
  {
    ARCANE_CHECK_AT(totalNbElement(), N);
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] += b[i];
  }

  //! Substract b to this vector
  void sub(const FixedVector<N>& b)
  {
    ARCANE_CHECK_AT(totalNbElement(), N);
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] -= b[i];
  }

 private:

  std::array<Arcane::Real, totalNbElement()> m_values = {};
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Tensor: used for symmetric 2nd-order tensors (useful for stresses, strains)
// Storage in vectorial form (xx yy zz xy yz zx)
using Tensor = FixedVector<6>;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Real
trace(const Tensor& b)
{
  return (b(0) + b(1) + b(3));
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Tensor
operator+(const Tensor& t1, const Tensor& t2)
{
  Tensor new_vector;

  for (Int32 i = 0; i < 6; ++i) {
    new_vector(i) = t1(i) + t2(i);
  }
  return new_vector;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Tensor
operator-(const Tensor& t1, const Tensor& t2)
{
  Tensor new_vector;

  for (Int32 i = 0; i < 6; ++i) {
    new_vector(i) = t1(i) - t2(i);
  }
  return new_vector;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Real3
tensorDiagonal(const Tensor& m)
{ // xx yy zz
  return { m(0), m(1), m(2) };
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Real3
tensorOutDiagonal(const Tensor& m)
{ // xy yz xz
  return { m(3), m(4), m(5) };
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Real3x3
tensorToMatrix3x3(const Tensor& m)
{
  Real3x3 mat;
  for (Int32 i = 0; i < 3; ++i)
    mat[i][i] = m(i);
  mat[0][1] = mat[1][0] = m(3);
  mat[0][2] = mat[2][0] = m(4);
  mat[1][2] = mat[2][1] = m(5);
  return mat;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Tensor
matrix3x3ToTensor(const Real3x3& m)
{
  Tensor t;
  for (Int32 i = 0; i < 3; ++i)
    t(i) = m[i][i];
  t(3) = m[0][1];
  t(4) = m[0][2];
  t(5) = m[2][1];
  return t;
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
