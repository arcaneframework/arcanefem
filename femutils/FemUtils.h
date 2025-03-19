// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemUtils.h                                                  (C) 2022-2025 */
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
  // Vector addition: Real4 + Real4
  ARCCORE_HOST_DEVICE Real4 operator+(const Real4& other) const
  {
    Real4 result;
    for (std::size_t i = 0; i < 4; ++i)
      result[i] = data[i] + other[i];
    return result;
  }

  // Vector subtraction: Real4 - Real4
  ARCCORE_HOST_DEVICE Real4 operator-(const Real4& other) const
  {
    Real4 result;
    for (std::size_t i = 0; i < 4; ++i)
      result[i] = data[i] - other[i];
    return result;
  }
  // Scalar multiplication: Real4 * scalar
  ARCCORE_HOST_DEVICE Real4 operator*(Arcane::Real scalar) const
  {
    Real4 result;
    for (std::size_t i = 0; i < 4; ++i)
      result[i] = data[i] * scalar;
    return result;
  }
  friend ARCCORE_HOST_DEVICE Real4 operator*(Arcane::Real scalar, const Real4& vec)
  {
    return vec * scalar;
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*!
 * \brief Matrix of size NxM.
 */
/*---------------------------------------------------------------------------*/
template <int N, int M>
class RealMatrix
{
  using ThatClass = RealMatrix<N, M>;

 public:

  static constexpr Arcane::Int32 totalNbElement() { return N * M; }
  ARCCORE_HOST_DEVICE RealMatrix() {};
  ARCCORE_HOST_DEVICE RealMatrix(std::initializer_list<Real> init_list)
  {
    auto i = 0;
    for (auto it = init_list.begin(); it != init_list.end(); it++) {
      m_values[i] = *it;
      i++;
    }
  };

  ARCCORE_HOST_DEVICE RealMatrix(std::initializer_list<std::initializer_list<Real>> init_list)
  {
    Arcane::Int32 i = 0;
    for (const auto& row : init_list) {
      Arcane::Int32 j = 0;
      for (const auto& value : row) {
          m_values[i * M + j] = value;
        ++j;
      }
      ++i;
    }
  }

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
  ARCCORE_HOST_DEVICE void multInPlace(Arcane::Real v)
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
  ARCCORE_HOST_DEVICE RealMatrix<N, M> operator+(const RealMatrix<N, M>& other) const
  {
    RealMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = (*this)(i, j) + other(i, j);
      }
    }
    return result;
  }

  //! Define the subtraction operator
  RealMatrix<N, M> operator-(const RealMatrix<N, M>& other) const
  {
    RealMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = (*this)(i, j) - other(i, j);
      }
    }
    return result;
  }

  //! Define the unary negation operator
  RealMatrix<N, M> operator-() const
  {
    RealMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = -(*this)(i, j);
      }
    }
    return result;
  }

  //! Scalar multiplication: RealMatrix * scalar
  ARCCORE_HOST_DEVICE RealMatrix<N, M> operator*(Real scalar) const
  {
    RealMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = (*this)(i, j) * scalar;
      }
    }
    return result;
  }

  //! Scalar division: RealMatrix / scalar
  ARCCORE_HOST_DEVICE RealMatrix<N, M> operator/(Real scalar) const
  {
    RealMatrix<N, M> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      for (Arcane::Int32 j = 0; j < M; ++j) {
        result(i, j) = (*this)(i, j) / scalar;
      }
    }
    return result;
  }

  //! Friend function for scalar multiplication: scalar * RealMatrix
  ARCCORE_HOST_DEVICE friend RealMatrix<N, M> operator*(Real scalar, const RealMatrix<N, M>& matrix)
  {
    RealMatrix<N, M> result;
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
//  Outer product of two Real3 vectors to produce a RealMatrix<3, 3>
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE inline RealMatrix<3, 3> operator^(const Arcane::Real3& lhs, const Arcane::Real3& rhs)
{
  RealMatrix<3, 3> result;
  for (Arcane::Int32 i = 0; i < 3; ++i) {
    for (Arcane::Int32 j = 0; j < 3; ++j) {
      result(i, j) = lhs[i] * rhs[j];
    }
  }
  return result;
}

/*---------------------------------------------------------------------------*/
//  Outer product of two Real4 vectors to produce a RealMatrix<4, 4>
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE inline RealMatrix<4, 4> operator^(const Real4& lhs, const Real4& rhs)
{
  RealMatrix<4, 4> result;
  for (Arcane::Int32 i = 0; i < 4; ++i) {
    for (Arcane::Int32 j = 0; j < 4; ++j) {
      result(i, j) = lhs[i] * rhs[j];
    }
  }
  return result;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <int N, int M> inline RealMatrix<N, N>
ARCCORE_HOST_DEVICE matrixAddition(const RealMatrix<N, M>& a, const RealMatrix<M, N>& b)
{
  using namespace Arcane;
  RealMatrix<N, N> new_matrix;

  for (Int32 i = 0; i < N; ++i) {
    for (Int32 j = 0; j < N; ++j) {
      new_matrix(i, j) = a(i, j) + b(i, j);
    }
  }
  return new_matrix;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <int N, int M> inline RealMatrix<N, N>
ARCCORE_HOST_DEVICE matrixMultiplication(const RealMatrix<N, M>& a, const RealMatrix<M, N>& b)
{
  using namespace Arcane;
  RealMatrix<N, N> new_matrix;

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

template <int N, int M> inline RealMatrix<M, N>
ARCCORE_HOST_DEVICE matrixTranspose(const RealMatrix<N, M>& a)
{
  using namespace Arcane;

  RealMatrix<M, N> t_matrix;
  for (Int32 i = 0; i < N; ++i) {
    for (Int32 j = 0; j < M; ++j) {
      t_matrix(j, i) = a(i, j);
    }
  }
  return t_matrix;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <int N> inline RealMatrix<N, N>
ARCCORE_HOST_DEVICE massMatrix(const RealMatrix<1, N>& lhs, const RealMatrix<1, N>& rhs)
{
  using namespace Arcane;

  RealMatrix<N, N> m_matrix;
  for (Arcane::Int32 i = 0; i < N; ++i)
    for (Arcane::Int32 j = 0; j < N; ++j)
      m_matrix(i, j) = lhs(0, i) * rhs(0, j);

  for (Arcane::Int32 i = 0; i < N; ++i)
    m_matrix(i, i) *= 2.;

  return m_matrix;
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Vector of size N.
 */
/*---------------------------------------------------------------------------*/
template <int N>
class RealVector
{
  using ThatClass = RealVector<N>;

 public:

  static constexpr Arcane::Int32 totalNbElement() { return N; }

  ARCCORE_HOST_DEVICE RealVector() {};

  ARCCORE_HOST_DEVICE RealVector(std::initializer_list<Arcane::Real> init_list)
  {
    Arcane::Int32 i = 0;
    for (auto it = init_list.begin(); it != init_list.end() && i < N; ++it, ++i) {
      m_values[i] = *it;
    }
  }

 public:

  ARCCORE_HOST_DEVICE Arcane::Real& operator()(Arcane::Int32 i)
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
  ARCCORE_HOST_DEVICE void multInPlace(Arcane::Real v)
  {
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] *= v;
  }

  //! Add \a v to all the components
  ARCCORE_HOST_DEVICE void addInPlace(Arcane::Real v)
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

  //! Define the addition operator
  ARCCORE_HOST_DEVICE RealVector<N> operator+(const RealVector<N>& other) const
  {
    RealVector<N> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(i) = (*this)(i) + other(i);
    }
    return result;
  }

  //! Define the subtraction operator
  RealVector<N> operator-(const RealVector<N>& other) const
  {
    RealVector<N> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(i) = (*this)(i)-other(i);
    }
    return result;
  }

  //! Define the unary negation operator
  RealVector<N> operator-() const
  {
    RealVector<N> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(i) = -(*this)(i);
    }
    return result;
  }

  //! Set this vector equal to b
  ARCCORE_HOST_DEVICE void setEqualTo(const RealVector<N>& b)
  {
    ARCANE_CHECK_AT(totalNbElement(), N);
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] = b[i];
  }

  //! Add b to this vector
  ARCCORE_HOST_DEVICE void add(const RealVector<N>& b)
  {
    ARCANE_CHECK_AT(totalNbElement(), N);
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] += b[i];
  }

  //! Substract b to this vector
  ARCCORE_HOST_DEVICE void sub(const RealVector<N>& b)
  {
    ARCANE_CHECK_AT(totalNbElement(), N);
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] -= b[i];
  }

  //! Scalar multiplication: RealMatrix * scalar
  ARCCORE_HOST_DEVICE RealVector<N> operator*(Real scalar) const
  {
    RealVector<N> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(i) = (*this)(i)*scalar;
    }
    return result;
  }

  //! Scalar multiplication: RealMatrix * scalar
  ARCCORE_HOST_DEVICE RealVector<N> operator/(Real scalar) const
  {
    RealVector<N> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(i) = (*this)(i) / scalar;
    }
    return result;
  }

  //! Friend function for scalar multiplication: scalar * RealMatrix
  ARCCORE_HOST_DEVICE friend RealVector<N> operator*(Real scalar, const RealVector<N>& vector)
  {
    RealVector<N> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(i) = scalar * vector(i);
    }
    return result;
  }

 private:

  std::array<Arcane::Real, totalNbElement()> m_values = {};
};

/*---------------------------------------------------------------------------*/
//  Matrix vector product of RealVector<N> vectors and matrix RealMatrix<N, N>
/*---------------------------------------------------------------------------*/

template <int N> inline RealVector<N>
ARCCORE_HOST_DEVICE operator*(const RealVector<N>& lhs, const RealMatrix<N, N>& rhs)
{
  RealVector<N> result;
  for (Arcane::Int32 j = 0; j < N; ++j) {
    result(j) = 0; // Initialize result element
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(j) += lhs(i) * rhs(i, j);
    }
  }
  return result;
}

/*---------------------------------------------------------------------------*/
//  Outer product of two RealVector<N> vectors to produce a RealMatrix<N, N>
/*---------------------------------------------------------------------------*/
template <int N> inline RealMatrix<N, N>
ARCCORE_HOST_DEVICE operator^(const RealVector<N>& lhs, const RealVector<N>& rhs)
{
  RealMatrix<N, N> result;
  for (Arcane::Int32 i = 0; i < N; ++i) {
    for (Arcane::Int32 j = 0; j < N; ++j) {
      result(i, j) = lhs(i) * rhs(j);
    }
  }
  return result;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <int N> inline RealMatrix<N, N>
ARCCORE_HOST_DEVICE massMatrix(const RealVector<N>& lhs, const RealVector<N>& rhs)
{
  using namespace Arcane;

  RealMatrix<N, N> m_matrix;
  for (Arcane::Int32 i = 0; i < N; ++i)
    for (Arcane::Int32 j = 0; j < N; ++j)
      m_matrix(i, j) = lhs(i) * rhs(j);

  for (Arcane::Int32 i = 0; i < N; ++i)
    m_matrix(i, i) *= 2.;

  return m_matrix;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Tensor: used for symmetric 2nd-order tensors (useful for stresses, strains)
// Storage in vectorial form (xx yy zz xy yz zx)
using Tensor = RealVector<6>;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Real
trace(const Tensor& b)
{
  return (b(0) + b(1) + b(3));
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
