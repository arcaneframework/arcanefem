// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* utilFEM.h                                                   (C) 2022-2024 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef PASSMO_UTILFEM_H_
#define PASSMO_UTILFEM_H_

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "FemUtils.h"
#include <arcane/Assertion.h>
#include "ArcaneFemFunctions.h"

using namespace Arcane;
using namespace Arcane::FemUtils;


/*---------------------------------------------------------------------------*/
Real trace(const Real3x3& mat) {
  return mat[0][0] + mat[1][1] + mat[2][2];
}
/*---------------------------------------------------------------------------*/
Real3 real3x3GetSupOutdiagonal(const Real3x3& mat) {
  return {mat[0][1], mat[0][2], mat[1][2]};
}
/*---------------------------------------------------------------------------*/
Real3 real3x3GetLowOutdiagonal(const Real3x3& mat){
  return {mat[1][0], mat[2][0], mat[2][1]};
}
/*---------------------------------------------------------------------------*/
Real3x3 diagonalReal3x3(const Real3x3& mat){
  Real3x3 newmat;
  newmat[0][0] = mat[0][0];
  newmat[1][1] = mat[1][1];
  newmat[2][2] = mat[2][2];
  return newmat;
}
/*---------------------------------------------------------------------------*/
Real3x3 outdiagonalReal3x3(const Real3x3& mat){
  return (mat - diagonalReal3x3(mat));
}

bool	real3x3IsSym(const Real3x3& mat)
{
  Real3 matsup = real3x3GetSupOutdiagonal(mat);
  Real3 matlow = real3x3GetLowOutdiagonal(mat);
  return (matsup == matlow);
}


/*---------------------------------------------------------------------------*/
/*!
 * \brief Symmetric Tensor used for stress and strain tensors
 */
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
  ARCCORE_HOST_DEVICE void multInPlace(Real v) { return m_vec.multInPlace(v); }

  //! Add \a v to all the components
  ARCCORE_HOST_DEVICE void addInPlace(Real v) { return m_vec.addInPlace(v); }

  //! Dump values
  void dump(std::ostream& o) const { m_vec.dump(o); }

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
  Tensor2 operator-(const Tensor2& other) const {
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
  Tensor2 operator-() const {
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
  ARCCORE_HOST_DEVICE void sub(const Tensor2& b) {
    for (Arcane::Int32 i = 0; i < 6; ++i)
      m_vec(i) -= b(i);
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

  ARCCORE_HOST_DEVICE Real trace(const Tensor2& vec)  {
    return m_vec(0) + m_vec(1) + m_vec(2);
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
typedef MeshVariableArrayRefT<Arcane::DoF,Tensor2> VariableDoFTensor;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
inline RealUniqueArray2
bothMultiply(const ElastTensor& a, const RealUniqueArray2& b){
#ifdef _DEBUG
  assert(b.dim1Size() == 6);
#endif
  Int32 n{b.dim2Size()};
  RealUniqueArray2 bt_ab(n,n);
  RealUniqueArray2 bt(n,6);
  RealUniqueArray2 ab(6,n);
  for (Int32 i = 0; i < 6; ++i) {
    for (Int32 j = 0; j < n; ++j) {
        bt[j][i] = b[i][j];
    }
  }
  for (Int32 i = 0; i < 6; ++i) {
    for (Int32 j = 0; j < n; ++j) {
      for (Int32 k = 0; k < 6; ++k) {
        ab[i][j] += a(i, k)*b[k][j];
      }
    }
  }
  for (Int32 i = 0; i < n; ++i) {
    for (Int32 j = 0; j < n; ++j) {
      for (Int32 k = 0; k < 6; ++k) {
        bt_ab[i][j] += bt[i][k]*ab[k][j];
      }
    }
  }
  return bt_ab;
}
*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
inline void
addArray2(RealUniqueArray2& a, const RealUniqueArray2& b, const Real& wt){
#ifdef _DEBUG
  assert(a.dim1Size() == b.dim1Size() && a.dim2Size() == b.dim2Size());
#endif
  for (Int32 i = 0; i < a.dim1Size(); ++i) {
    for (Int32 j = 0; j <  a.dim2Size(); ++j) {
      a[i][j] += wt*b[i][j];
    }
  }
}
*/


#endif // PASSMO_UTILFEM_H_
