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
#include "GaussQuadrature.h"

using namespace Arcane;
using namespace Arcane::FemUtils;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Functions useful for class CellFEMDispatcher
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
extern void DirVectors(const Face& face,const VariableNodeReal3& n, const Int32& ndim, Real3& e1, Real3& e2, Real3& e3);

extern Real Line2Length(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Line2ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Line2ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Line2Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Line3Length(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Line3ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Line3ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Line3Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Tri3Surface(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Tri3ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Tri3ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Tri3Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Tri6Surface(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Tri6ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Tri6ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Tri6Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Quad4Surface(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Quad4ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Quad4ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Quad4Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Quad8Surface(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Quad8ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Quad8ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Quad8Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Hexa8Volume(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Hexa8ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Hexa8ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Hexa8Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Hexa20Volume(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Hexa20ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Hexa20ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Hexa20Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Tetra4Volume(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Tetra4ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Tetra4ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Tetra4Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Tetra10Volume(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Tetra10ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Tetra10ShapeFuncDeriv(const Integer& inod, const Real3& coord);
//extern Integer3 Tetra10Orientation(const ItemWithNodes& item, const VariableNodeReal3& n);

extern Real Penta6Volume(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Penta6ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Penta6ShapeFuncDeriv(const Integer& inod, const Real3& coord);

extern Real Pyramid5Volume(const ItemWithNodes& item, const VariableNodeReal3& n);
extern Real Pyramid5ShapeFuncVal(const Integer& inod, const Real3& coord);
extern Real3 Pyramid5ShapeFuncDeriv(const Integer& inod, const Real3& coord);

extern Int32 getGeomDimension(const ItemWithNodes& item);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class CellFEMDispatcher
{
 public:
  CellFEMDispatcher();

 public:

//  void set_node_coords(VariableNodeReal3& /*node_coords*/);

//  Real3 getBarycenter(const ItemWithNodes& /*cell*/);

  Real getShapeFuncVal(const Int16& /*item_type*/, const Int32& /*inod*/, const Real3& /*ref coord*/);
  Real3 getShapeFuncDeriv(const Int16& /*item_type*/, const Int32& /*inod*/, const Real3& /*ref coord*/);

//  Integer3 getOrientation(const ItemWithNodes& /*cell*/);
  RealUniqueArray getGaussData(const ItemWithNodes& item, const Integer& nint, Int32& ngauss);

 private:

//  std::function<Real(const ItemWithNodes& item, const VariableNodeReal3& n)> m_geomfunc[NB_BASIC_ITEM_TYPE];
  std::function<Real(const Int32& inod, const Real3& coord)> m_shapefunc[NB_BASIC_ITEM_TYPE];
  std::function<Real3(const Int32& inod, const Real3& coord)> m_shapefuncderiv[NB_BASIC_ITEM_TYPE];
//  std::function<Integer3(const ItemWithNodes& item, const VariableNodeReal3& n)> m_orientfunc[NB_BASIC_ITEM_TYPE];
//  VariableNodeReal3 *m_node_coords{ nullptr};
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// class ElastTensor: 4th-order elastic tensor
//
//      _          _
//     |  D      0  |
// C = |            |
//     |            |
//     |  0      S  |
//     |_          _|
//
//      _                                      _
//     | lambda + 2mu   lambda           lambda |  (Real3x3)
// D = |                                        |
//     | lambda       lambda + 2mu       lambda |
//     |                                        |
//     | lambda         lambda     lambda + 2mu |
//     |_                                      _|
//      _           _
//     | mu   0    0 | (Real3x3)
// S = |             |
//     | 0    mu   0 |
//     |             |
//     | 0    0   mu |
//     |_           _|


class ElastTensor{
 public:
  using ThatClass = ElastTensor;
 private:
//  FixedMatrix<6,6>	m_values;
  Real lambda{};
  Real mu{};
  Real lamba2mu{};
  Int32 dim{};

 public:
  ElastTensor() = default;
  ElastTensor(const Real& young, const Real& nu,const Int32& ndim): dim(ndim){
    mu = young/2./(1 + nu);
    lambda = 2.*mu*nu/(1 - 2.*nu);
    lamba2mu = lambda + 2.*mu;
/*    for (Arcane::Int32 i = 0; i < 3; ++i) {
      m_values(i, i) = lambda + 2 * mu;
      m_values(i+3, i+3) = mu;

      for (Arcane::Int32 j = i + 1; j < 3; ++j)
        m_values(i, j) = m_values(j, i) = lambda;
    }*/

  }
  ~ElastTensor() = default;

 public:

//  Arcane::Real& operator()(Arcane::Int32 i,Arcane::Int32 j){
//    return m_values(i,j);
//  }
  [[nodiscard]] auto getDim() const { return dim; }

  Arcane::Real operator()(Arcane::Int32 i,Arcane::Int32 j) const{
//    return m_values(i,j);
    if (i < dim){
      if (i == j)
        return lamba2mu;
      if (j < dim)
        return lambda;
    }
    if (i == j)
      return mu;
    return 0.;
  }
};
/*---------------------------------------------------------------------------*/
inline Real
trace(ElastTensor m)
{
  auto ndim = m.getDim();
  Real x{0.};
  for (Arcane::Int32 i = 0; i <= ndim; ++i)
    x += m(i, i);
  return x;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
template <int N>
class FixedVector{
  using ThatClass = FixedVector<N>;

 public:

  static constexpr Arcane::Int32 totalNbElement() { return N; }

 public:

  Arcane::Real& operator()(Arcane::Int32 i)  {
    ARCANE_CHECK_AT(i, N);
    return m_values[i];
  }

  Arcane::Real operator()(Arcane::Int32 i) const  {
    ARCANE_CHECK_AT(i, N);
    return m_values[i];
  }

 public:

  //! Multiply all the components by \a v
  void multInPlace(Arcane::Real v)  {
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] *= v;
  }

  //! Add \a v to all the components
  void addInPlace(Arcane::Real v)  {
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] += v;
  }

  //! Dump values
  void dump(std::ostream& o) const  {
    const ThatClass& values = *this;
    for (Arcane::Int32 i = 0; i < N; ++i) {
      o << "[ ";
      o << values(i);
      o << "]\n";
    }
  }

  //! Set this vector equal to b
  void setEqualTo(const FixedVector<N>& b)  {
    ARCANE_CHECK_AT(totalNbElement(), N);
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] = b[i];
  }

  //! Add b to this vector
  void add(const FixedVector<N>& b)  {
    ARCANE_CHECK_AT(totalNbElement(), N);
    for (Arcane::Int32 i = 0, n = totalNbElement(); i < n; ++i)
      m_values[i] += b[i];
  }

  //! Substract b to this vector
  void sub(const FixedVector<N>& b)  {
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
using Tensor= FixedVector<6>;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Tensor
tensorMultiply(const ElastTensor& a, const Tensor& b){
  Tensor new_vector;

  for (Int32 i = 0; i < 6; ++i) {
    Real x = 0.0;
    for (Int32 j = 0; j < 6; ++j) {
        x += a(i, j) * b(j);
    }
    new_vector(i) += x;
  }
  return new_vector;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Real
trace(const Tensor& b){
  return  (b(0) + b(1) + b(3));
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Tensor
operator+(const Tensor& t1,const Tensor& t2) {
  Tensor new_vector;

  for (Int32 i = 0; i < 6; ++i) {
    new_vector(i) = t1(i) + t2(i);
  }
  return new_vector;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Tensor
operator-(const Tensor& t1,const Tensor& t2) {
  Tensor new_vector;

  for (Int32 i = 0; i < 6; ++i) {
    new_vector(i) = t1(i) - t2(i);
  }
  return new_vector;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Real3
tensorDiagonal(const Tensor& m) { // xx yy zz
  return {m(0), m(1), m(2) };
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Real3
tensorOutDiagonal(const Tensor& m) { // xy yz xz
  return {m(3), m(4), m(5) };
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Real3x3
tensorToMatrix3x3(const Tensor& m){
  Real3x3 mat;
  for (Int32 i = 0; i < 3; ++i) mat[i][i] = m(i);
  mat[0][1] = mat[1][0] = m(3);
  mat[0][2] = mat[2][0] = m(4);
  mat[1][2] = mat[2][1] = m(5);
  return mat;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline Tensor
matrix3x3ToTensor(const Real3x3& m){
  Tensor t;
  for (Int32 i = 0; i < 3; ++i) t(i) = m[i][i];
  t(3) = m[0][1];
  t(4) = m[0][2];
  t(5) = m[2][1];
  return t;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
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


#endif // PASSMO_UTILFEM_H_
