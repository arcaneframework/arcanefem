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
#include <arcane/Assertion.h>
#include <arcane/VariableTypes.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
extern Real real3x3Trace(const Real3x3& mat);
/*---------------------------------------------------------------------------*/
extern Real3 real3x3GetSupOutdiagonal(const Real3x3& mat);
/*---------------------------------------------------------------------------*/
extern Real3 real3x3GetLowOutdiagonal(const Real3x3& mat);
/*---------------------------------------------------------------------------*/
extern Real3x3 diagonalReal3x3(const Real3x3& mat);
/*---------------------------------------------------------------------------*/
extern Real3x3 outdiagonalReal3x3(const Real3x3& mat);

extern bool	real3x3IsSym(const Real3x3& mat);
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
  Real3x3 m_tab[4];//0 = D, 1 = S, 2 = Sup, 3 = Slow
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

  ARCCORE_HOST_DEVICE Real3x3 operator[](Int32 i) const {
    ARCANE_CHECK_AT(i,4);
    return m_tab[i];
  }
  ARCCORE_HOST_DEVICE Real3x3& operator[](Int32 i) {
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

  //! Define the trace of the Tensor
  ARCCORE_HOST_DEVICE [[nodiscard]] Real trace() const {
    return real3x3Trace(m_tab[0]) + real3x3Trace(m_tab[1]);
  }

  //! Get the (i,j) value of the Tensor
  ARCCORE_HOST_DEVICE [[nodiscard]] Real value_at(Int32 i,Int32 j) const
  {
    ARCANE_CHECK_AT(i, 6);
    ARCANE_CHECK_AT(j, 6);

    if (i < 3) {

      if (j < 3) return m_tab[0][i][j]; // D(i,j)
      return m_tab[2][i][j-3]; // Sup(i,j-3);
    }

    if (j < 3) {

      if (m_tab[3] == Real3x3::zero())// if Slow == 0
        return m_tab[2][j][i-3]; // Sup(j,i-3)
      return m_tab[3][i-3][j]; // Slow(i-3,j)
    }

    return m_tab[1][i-3][j-3];// S(i-3,j-3)
  }

  //! Set the (i,j) value of the Tensor
  ARCCORE_HOST_DEVICE void value_at(Int32 i,Int32 j, Real x)
  {
    ARCANE_CHECK_AT(i, 6);
    ARCANE_CHECK_AT(j, 6);

    if (i < 3) {

      if (j < 3)
        m_tab[0][i][j] = x; // D(i,j)
      else
        m_tab[2][i][j-3] = x; // Sup(i,j-3);

    } else {

      if (j < 3) {

        if (m_tab[3] == Real3x3::zero()) // if Slow == 0
          m_tab[2][j][i - 3] = x; // Sup(j,i-3)
        else
          m_tab[3][i - 3][j] = x; // Slow(i-3,j)

      } else {
        m_tab[1][i - 3][j - 3] = x; // S(i-3,j-3)
      }
    }
  }

  //! Get the (i,j,k,l) value of the Tensor4
  ARCCORE_HOST_DEVICE Real value_at(Int32 i,Int32 j, Int32 k, Int32 l)
  {
    ARCANE_CHECK_AT(i, 3);
    ARCANE_CHECK_AT(j, 3);
    ARCANE_CHECK_AT(k, 3);
    ARCANE_CHECK_AT(l, 3);

    int p = math::abs(i - j),q = math::abs(k - l);
    if (!p)
    {
      if (!q) return m_tab[0][i][k]; //D(i,k)
      if (q == 1 && k) q++;
      else q--;
      return m_tab[2][i][q]; //Sup(i,q)
    }

    if (p == 1 && i) p++;
    else p--;

    if (!q) {

      if (m_tab[3] == Real3x3::zero()) // if Slow == 0
        return m_tab[2][p][k];//Sup(p,k)

      return m_tab[3][p][k];//Slow(p,k)
    }
    if (q == 1 && k) q++;
    else q--;
    return m_tab[1][p][q];//S(p,q)
  }

  //! Set the (i,j,k,l) value of the Tensor4
  ARCCORE_HOST_DEVICE void value_at(Int32 i,Int32 j, Int32 k, Int32 l, Real x)
  {
    ARCANE_CHECK_AT(i, 3);
    ARCANE_CHECK_AT(j, 3);
    ARCANE_CHECK_AT(k, 3);
    ARCANE_CHECK_AT(l, 3);

    int p = math::abs(i - j),q = math::abs(k - l);
    if (!p)
    {
      if (!q) m_tab[0][i][k] = x;// D(i,k)
      else
      {
        if (q == 1 && k) q++;
        else q--;
        m_tab[2][i][q] = x;//Sup(i,q) = x;
      }
    }
    else
    {
      if (p == 1 && i) p++;
      else p--;

      if (!q) {

        if (m_tab[3] == Real3x3::zero()) // if Slow == 0
          m_tab[2][p][k] = x;//Sup(p,k);
        else
          m_tab[3][p][k] = x;//Slow(p,k)
      }
      else
      {
        if (q == 1 && k) q++;
        else q--;
        m_tab[1][p][q] = x;//S(p,q);
      }
    }
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

#endif // PASSMO_UTILFEM_H_
