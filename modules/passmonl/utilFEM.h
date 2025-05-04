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
Real real3x3Trace(const Real3x3& mat) {
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
 * \brief 2nd-order symmetric tensor used for stress and strain tensors
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
  ARCCORE_HOST_DEVICE void multInPlace(Real v) { m_vec.multInPlace(v); }

  //! Add \a v to all the components
  ARCCORE_HOST_DEVICE void addInPlace(Real v) { m_vec.addInPlace(v); }

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
/*!
 * \brief 4th-order tensor used for constitutive operators
//      _          _
//     |  D     Sup |
// C = |             | with D, S, Sup, Slow = Real3x3 blocks
//     |             |
//     |  Slow   S  |
//     |_           _|
//
// example: for elastic Tensor4
//      _                                      _
//     | lambda + 2mu   lambda           lambda |
// D = |                                        |
//     | lambda       lambda + 2mu       lambda |
//     |                                        |
//     | lambda         lambda     lambda + 2mu |
//     |_                                      _|
//      _           _
//     | mu   0    0 |
// S = |             |  Sup = Slow = 0 (Real3x3)
//     | 0    mu   0 |
//     |             |
//     | 0    0   mu |
//     |_           _|
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class Tensor4
{
 // ***** ATTRIBUTES
 protected:
  Real3x3 m_tab[4];
  bool 	m_constitutive;// introduced if some special operations are required for constitutive models
  bool 	m_sym; 		// all Tensor2 are symmetric in this case

 public:

  ARCCORE_HOST_DEVICE explicit Tensor4(bool constitutive = true, bool sym = true): m_constitutive(constitutive), m_sym(sym)
  {}

  ARCCORE_HOST_DEVICE Tensor4(const Real3x3& D,const Real3x3& S){
    m_constitutive = true;
    m_tab[0] = D;
    m_tab[1] = S;
    m_sym = (real3x3IsSym(D) && real3x3IsSym(S));
  }

  ARCCORE_HOST_DEVICE Tensor4(const Real3x3& D,const Real3x3& S,const Real3x3& Sup,const Real3x3& Slow,bool constitutive = true)  {
    m_constitutive = constitutive;
    m_tab[0] = D;
    m_tab[1] = S;
    m_tab[2] = Sup;
    m_tab[3] = Slow;
    m_sym = (real3x3IsSym(D) && real3x3IsSym(S) && Sup == Slow);
  }

  ARCCORE_HOST_DEVICE Tensor4(Real lambda, Real mu) {
    m_sym = true;
    m_constitutive = true;
    m_tab[0] = lambda;
    for (int i = 0; i < 3; i++) m_tab[0][i][i] += 2*mu;
    m_tab[1] = mu*Real3x3::identity();
  }

  ARCCORE_HOST_DEVICE Tensor4(const Tensor4& tens)  {
    for (int i = 0; i < 4; i++) m_tab[i] = tens.m_tab[i];
    m_sym = tens.m_sym;
    m_constitutive = tens.m_constitutive;
  }

  ARCCORE_HOST_DEVICE ~Tensor4() = default;

 public:

  ARCCORE_HOST_DEVICE Real3x3 operator[](Int32& i) const {
    ARCANE_CHECK_AT(i,4);
    return m_tab[i];
  }
  ARCCORE_HOST_DEVICE Real3x3& operator[](Int32& i) {
    ARCANE_CHECK_AT(i,4);
    return m_tab[i];
  }

  //! Return Identity Tensor4
  static Tensor4 identity() {
    return { Real3x3::identity(), Real3x3::identity() };
  }

  //! Return zero Tensor2
  static Tensor4 zero() {
    return { Real3x3::zero(), Real3x3::zero() };
  }

 public:
  ARCCORE_HOST_DEVICE Tensor4& operator=(const Tensor4& tens) {
    for (int i = 0; i < 4; i++) m_tab[i] = tens.m_tab[i];
    m_sym = tens.m_sym;
    m_constitutive = tens.m_constitutive;
    return (*this);
  }

  //! Define the addition operator
  ARCCORE_HOST_DEVICE Tensor4 operator+(const Tensor4& tens) const {
    Tensor4 result;
    for (Arcane::Int32 i = 0; i < 4; ++i) {
      result.m_tab[i] = m_tab[i] + tens[i];
    }
    return result;
  }

  //! Define the += operator
  ARCCORE_HOST_DEVICE Tensor4& operator+=(const Tensor4& tens) {
    for (Int32 i = 0; i < 4; i++) m_tab[i] += tens[i];
    return (*this);
  }

  //! Define the subtraction operator
  ARCCORE_HOST_DEVICE Tensor4 operator-(const Tensor4& tens) const {
    Tensor4 result;
    for (Arcane::Int32 i = 0; i < 4; ++i) {
      result.m_tab[i] = m_tab[i] + tens[i];
    }
    return result;
  }

  //! Define the -= operator
  ARCCORE_HOST_DEVICE Tensor4& operator-=(const Tensor4& tens) {
    for (Int32 i = 0; i < 4; i++) m_tab[i] -= tens[i];
    return (*this);
  }

  //! Multiply all the components by \a v
  ARCCORE_HOST_DEVICE void multInPlace(Real scalar) {
    for (auto & i : m_tab) {
      i *= scalar;
    }
  }

  //! Define the unary negation operator
  ARCCORE_HOST_DEVICE Tensor4 operator-() const {
    Tensor4 tens(*this);
    for (Int32 i = 0; i < 4; i++) tens[i] *= (-1);
    return tens;
  }

  //! Add \a v to all the components
  ARCCORE_HOST_DEVICE void addInPlace(Real scalar) {
    for (auto & i : m_tab) {
      i.addSame(Real3(scalar, scalar, scalar));
    }
  }

  //! Substract \a v to all the components
  ARCCORE_HOST_DEVICE void subInPlace(Real v) {
    for (auto & i : m_tab) {
      i.subSame(Real3(v, v, v));
    }
  }

  //! Scalar multiplication: Tensor4 * scalar
  ARCCORE_HOST_DEVICE Tensor4 operator*(Real scalar) const {
    Tensor4 result;
    for (auto & i : result.m_tab) {
      i *= scalar;
    }
    return result;
  }

  //! Define the *= operator
  ARCCORE_HOST_DEVICE Tensor4& operator*=(Real scalar) {
    *this = this->operator*(scalar);
    return (*this);
  }

  //! Scalar division: Tensor / scalar
  ARCCORE_HOST_DEVICE Tensor4 operator/(Real scalar) const {
    Tensor4 result;
    for (auto & i : result.m_tab) {
      i /= scalar;
    }
    return result;
  }

  //! Define the *= operator
  ARCCORE_HOST_DEVICE Tensor4& operator/=(Real scalar) {
    *this = this->operator/(scalar);
    return (*this);
  }

  ARCCORE_HOST_DEVICE Real trace(const Tensor4& vec)  {
    return real3x3Trace(m_tab[0]) + real3x3Trace(m_tab[1]);
  }

  //! Friend function for scalar multiplication: scalar * Tensor4
  ARCCORE_HOST_DEVICE friend Tensor4 operator*(Real scalar, const Tensor4& tens) {
    Tensor4 result;
    for (Arcane::Int32 i = 0; i < 4; ++i) {
      result.m_tab[i] = scalar * tens.m_tab[i];
    }
    return result;
  }
};
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
