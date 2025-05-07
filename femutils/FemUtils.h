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
#include <arcane/MeshVariableArrayRef.h>
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

  ARCCORE_HOST_DEVICE RealVector() = default;

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

  //! Scalar multiplication: RealVector * scalar
  ARCCORE_HOST_DEVICE RealVector<N> operator*(Real scalar) const
  {
    RealVector<N> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(i) = (*this)(i)*scalar;
    }
    return result;
  }

  //! Scalar division: RealVector / scalar
  ARCCORE_HOST_DEVICE RealVector<N> operator/(Real scalar) const
  {
    RealVector<N> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(i) = (*this)(i) / scalar;
    }
    return result;
  }

  //! Friend function for scalar multiplication: scalar * RealVector
  ARCCORE_HOST_DEVICE friend RealVector<N> operator*(Real scalar, const RealVector<N>& vector)
  {
    RealVector<N> result;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result(i) = scalar * vector(i);
    }
    return result;
  }

  //! Friend function for dot product of RealVector
  ARCCORE_HOST_DEVICE friend Real dot(const RealVector<N>& u, const RealVector<N>& v)
  {
    Real result{0.};
    for (Arcane::Int32 i = 0; i < N; ++i) {
      result += u(i) * v(i);
    }
    return result;
  }

  // private:
 protected:

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
class Tensor2
{
  RealVector<6> m_vec;

 public:

  static constexpr Int32 totalNbElement() { return 6; }

  ARCCORE_HOST_DEVICE Tensor2() = default;

  ARCCORE_HOST_DEVICE Tensor2(std::initializer_list<Real> init_list):m_vec(init_list)
  {}

  ARCCORE_HOST_DEVICE Tensor2(const Real3& a,const Real3& b) {
    m_vec = {a.x, a.y, a.z, b.x, b.y, b.z};
  }
  ARCCORE_HOST_DEVICE explicit Tensor2(const RealVector<6>& vec) {
    m_vec = {vec(0), vec(1), vec(2), vec(3), vec(4), vec(5)};
  }

  ARCCORE_HOST_DEVICE explicit Tensor2(const Real3x3& mat) {
    ARCANE_ASSERT(real3x3IsSym(mat), true);
    for (Arcane::Int32 i = 0; i < 3; i++) m_vec(i) = mat[i][i];
    for (Arcane::Int32 i = 3; i < 5; i++) m_vec(i) = mat[0][i - 2];
    m_vec(5) = mat[1][2];
  }

  ARCCORE_HOST_DEVICE ~Tensor2() = default;

 public:

  ARCCORE_HOST_DEVICE Real& operator()(Int32 i) {
    return m_vec(i);
  }

  ARCCORE_HOST_DEVICE [[nodiscard]] RealVector<6> getVec() const {  return m_vec; }

  ARCCORE_HOST_DEVICE Real operator()(Int32 i) const {
    return m_vec(i);
  }

  ARCCORE_HOST_DEVICE Real& operator()(Int32 i, Int32 j)
  {
    int ij = get_index(i, j);
    return m_vec(ij);
  }

  ARCCORE_HOST_DEVICE Real operator()(Int32 i, Int32 j) const
  {
    int ij = get_index(i, j);
    return m_vec(ij);
  }

  //! Convert this Tensor to Real3x3 matrix
  ARCCORE_HOST_DEVICE explicit operator Real3x3() const {
    Real3x3 mat;
    for (Int32 i = 0; i < 3; i++) mat[i][i] = (*this)(i);
    for (Int32 i = 3; i < 5; i++) {
      mat[0][i - 2] = (*this)(i);
      mat[i - 2][0] = (*this)(i);
    }
    mat[1][2] = mat[2][1] = (*this)(5);

    return mat;
  }

  //! Return Identity Tensor2
  static Tensor2 identity() { return {Real3(1.,1.,1.), Real3::zero()};  }
  //! Return zero Tensor2
  static Tensor2 zero() { return {};  }

 public:
  //! Multiply all the components by \a v
  ARCCORE_HOST_DEVICE void multInPlace(Real v) { m_vec.multInPlace(v); }

  //! Add \a v to all the components
  ARCCORE_HOST_DEVICE void addInPlace(Real v) { m_vec.addInPlace(v); }

  //! Dump values
  void dump(std::ostream& o) const { m_vec.dump(o); }

  //! Define the == operator
  ARCCORE_HOST_DEVICE bool operator==(const Tensor2& vec) {
    Real eps{1.0e-15};
    for (Arcane::Int32 i = 0; i < 6; ++i) {
      if (fabs(m_vec(i) - vec(i)) > eps)
        return false;
    }
    return true;
  }

  //! Define the = operator
  ARCCORE_HOST_DEVICE Tensor2& operator=(const Tensor2& vec) {
    m_vec = {vec(0), vec(1), vec(2), vec(3), vec(4), vec(5)};
    return (*this);
  }

  //! Define the = operator
  ARCCORE_HOST_DEVICE Tensor2& operator=(const Real3x3& mat) {
    ARCANE_ASSERT(real3x3IsSym(mat), true);
    for (Arcane::Int32 i = 0; i < 3; i++) m_vec(i) = mat[i][i];
    for (Arcane::Int32 i = 3; i < 5; i++) m_vec(i) = mat[0][i - 2];
    m_vec(5) = mat[1][2];
    return (*this);
  }

  //! Define the addition operator
  ARCCORE_HOST_DEVICE Tensor2 operator+(const Tensor2& other) const {
    Tensor2 result;
    for (Arcane::Int32 i = 0; i < 6; ++i) {
      result(i) = m_vec(i) + other(i);
    }
    return result;
  }

  //! Define the += operator
  ARCCORE_HOST_DEVICE Tensor2& operator+=(const Tensor2& vec) {
    *this = this->operator+(vec);
    return (*this);
  }

  //! Define the subtraction operator
  ARCCORE_HOST_DEVICE Tensor2 operator-(const Tensor2& other) const {
    Tensor2 result;
    for (Arcane::Int32 i = 0; i < 6; ++i) {
      result(i) = m_vec(i)-other(i);
    }
    return result;
  }

  //! Define the -= operator
  ARCCORE_HOST_DEVICE Tensor2& operator-=(const Tensor2& vec) {
    *this = this->operator-(vec);
    return (*this);
  }

  //! Define the unary negation operator
  ARCCORE_HOST_DEVICE Tensor2 operator-() const {
    Tensor2 result;
    for (Arcane::Int32 i = 0; i < 6; ++i) {
      result(i) = -m_vec(i);
    }
    return result;
  }

  //! Set this Tensor2 equal to b
  ARCCORE_HOST_DEVICE void setEqualTo(const Tensor2& b) {
    for (Arcane::Int32 i = 0; i < 6; ++i)
      m_vec(i) = b(i);
  }

  ARCCORE_HOST_DEVICE void setVec(const Real3& d, const Real3& s) {
    m_vec = {};
    Arcane::Int32 i{0};
    for (; i < 3; i++) m_vec(i) = d[i]; // xx yy zz
    for (; i < 6; i++) m_vec(i) = s[i-3]; // xy xz yz
  }

  //! Add b to this Tensor2
  ARCCORE_HOST_DEVICE void add(const Tensor2& b) {
    for (Arcane::Int32 i = 0; i < 6; ++i)
      m_vec(i) += b(i);
  }
  //! Substract b to this Tensor2
  ARCCORE_HOST_DEVICE void add(const Real3x3& b) {
    Tensor2 tb(b);
    this->add(tb);
  }

  //! Substract b to this Tensor2
  ARCCORE_HOST_DEVICE void sub(const Tensor2& b) {
    for (Arcane::Int32 i = 0; i < 6; ++i)
      m_vec(i) -= b(i);
  }
  //! Substract b to this Tensor2
  ARCCORE_HOST_DEVICE void sub(const Real3x3& b) {
    Tensor2 tb(b);
    this->sub(tb);
  }

  //! Scalar multiplication: Tensor * scalar
  ARCCORE_HOST_DEVICE Tensor2 operator*(Real scalar) const {
    Tensor2 result;
    for (Arcane::Int32 i = 0; i < 6; ++i) {
      result(i) = m_vec(i)*scalar;
    }
    return result;
  }

  //! Define the *= operator
  ARCCORE_HOST_DEVICE Tensor2& operator*=(Real scalar) {
    *this = this->operator*(scalar);
    return (*this);
  }

  //! Scalar division: Tensor / scalar
  ARCCORE_HOST_DEVICE Tensor2 operator/(Real scalar) const {
    Tensor2 result;
    for (Arcane::Int32 i = 0; i < 6; ++i) {
      result(i) = m_vec(i) / scalar;
    }
    return result;
  }

  //! Define the *= operator
  ARCCORE_HOST_DEVICE Tensor2& operator/=(Real scalar) {
    *this = this->operator/(scalar);
    return (*this);
  }

  static Arcane::Int32 get_index(Arcane::Int32 i, Arcane::Int32 j) {
    ARCANE_CHECK_AT(i, 3);
    ARCANE_CHECK_AT(j, 3);
    if (j == i) return i;

    Arcane::Int32 ij{0}, nrow{3};

    if (j > i) ij = (i + 1) * (2 * nrow - i) / 2 + j - i - 1;
    else ij = (j + 1) * (2 * nrow - j) / 2 + i - j - 1;
    return ij;
  }

  ARCCORE_HOST_DEVICE [[nodiscard]] Real3 get_diagonal() const { // xx yy zz
    return {m_vec(0), m_vec(1), m_vec(2)};
  }

  ARCCORE_HOST_DEVICE void set_diagonal(const Real3& d) {
    for (Int32 i = 0; i < 3; i++) m_vec(i) = d[i];
  }

  ARCCORE_HOST_DEVICE [[nodiscard]] Real3 get_outdiagonal() const { // xy xz yz
    return {m_vec(3), m_vec(4), m_vec(5)};
  }
  ARCCORE_HOST_DEVICE void set_outdiagonal(const Real3& s) {
    for (Int32 i = 3; i < 6; i++) m_vec(i) = s[i-3];
  }

  ARCCORE_HOST_DEVICE [[nodiscard]] Real trace() const {
    return m_vec(0) + m_vec(1) + m_vec(2);
  }

  ARCCORE_HOST_DEVICE [[nodiscard]] Real norm() const {
    return sqrt(dot(m_vec,m_vec));
  }

  //! Function to convert Real3x3 matrix (symmetric) to Tensor
  ARCCORE_HOST_DEVICE void fromReal3x3ToTensor2(const Real3x3& mat) {
    for (Arcane::Int32 i = 0; i < 3; i++) (*this)(i) = mat[i][i];
    for (Arcane::Int32 i = 3; i < 5; i++) (*this)(i) = mat[0][i - 2];
    (*this)(5) = mat[1][2];
  }


  //! Friend function for scalar multiplication: scalar * Tensor
  ARCCORE_HOST_DEVICE friend Tensor2 operator*(Real scalar, const Tensor2& vector) {
    Tensor2 result;
    for (Arcane::Int32 i = 0; i < 6; ++i) {
      result(i) = scalar * vector(i);
    }
    return result;
  }

  //! Friend function to convert Tensor to Real3x3 matrix
  ARCCORE_HOST_DEVICE friend Real3x3 fromTensor2Real3x3(const Tensor2& vector) {
    Real3x3 mat;
    for (Arcane::Int32 i = 0; i < 3; i++) mat[i][i] = vector(i);
    for (Arcane::Int32 i = 3; i < 5; i++) {
      mat[0][i - 2] = vector(i);
      mat[i - 2][0] = vector(i);
    }
    mat[1][2] = mat[2][1] = vector(5);

    return mat;
  }
  //! Friend function to convert Real3x3 matrix (symmetric) to Tensor
  ARCCORE_HOST_DEVICE friend Tensor2 fromReal3x3ToTensor2(const Real3x3& mat) {
    Tensor2 vector;
    for (Arcane::Int32 i = 0; i < 3; i++) vector(i) = mat[i][i];
    for (Arcane::Int32 i = 3; i < 5; i++) vector(i) = mat[0][i - 2];
    vector(5) = mat[1][2];
    return vector;
  }
};

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
