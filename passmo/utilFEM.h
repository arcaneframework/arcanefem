/*
 * PASSMO : Performant Assessment for Seismic Site Modelling
 *
 * Definition of classes to implement finite-element cells and related shape functions
 * and their derivatives
 *
 * utilFEM.h: declarations
 *
 *  Created on: October 2019
 *      Author: E. Foerster
 */

#ifndef PASSMO_UTILFEM_H_
#define PASSMO_UTILFEM_H_

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "FemUtils.h"
#include <arcane/Assertion.h>

using namespace Arcane;
using namespace Arcane::FemUtils;

extern Real	REL_PREC;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Tables of abscissa & integration weights for Gauss-Legendre quadrature
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// For edges & quadrangles only
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// order 1
const Real xgauss1 = 0.;
const Real wgauss1 = 2.;

// order 2
const Real xgauss2[2] = { -1. / sqrt(3.), 1. / sqrt(3.) };
const Real wgauss2[2] = { 1., 1. };

// order 3
const Real xgauss3[3] = { -sqrt(0.6), 0., sqrt(0.6) };
const Real wgauss3[3] = { 5. / 9., 8. / 9., 5. / 9. };

// order 4
const Real xgauss4[4] = {
  -sqrt((3. + 2. * sqrt(1.2)) / 7.), -sqrt((3. - 2. * sqrt(1.2)) / 7.),
  sqrt((3. - 2. * sqrt(1.2)) / 7.), sqrt((3. + 2. * sqrt(1.2)) / 7.)
};

const Real wgauss4[4] = {
  0.5 - sqrt(5. / 6.) / 6., 0.5 + sqrt(5. / 6.) / 6.,
  0.5 + sqrt(5. / 6.) / 6., 0.5 - sqrt(5. / 6.) / 6.
};

// order 5
const Real xgauss5[5] = {
  -sqrt(245. + 14. * sqrt(70.)) / 21., -sqrt(245. - 14. * sqrt(70.)) / 21., 0.,
  sqrt(245. - 14. * sqrt(70.)) / 21., sqrt(245. + 14. * sqrt(70.)) / 21.
};

const Real wgauss5[5] = {
  (322. - 13 * sqrt(70.)) / 900., (322. + 13 * sqrt(70.)) / 900., 128. / 225.,
  (322. + 13 * sqrt(70.)) / 900., (322. - 13 * sqrt(70.)) / 900.
};

// order 6
const Real xgauss6[6] = {
  -0.932469514203152, -0.661209386466265,
  -0.238619186083197, 0.238619186083197,
  0.661209386466265, 0.932469514203152
};

const Real wgauss6[6] = {
  0.171324492379170, 0.360761573048139,
  0.467913934572691, 0.467913934572691,
  0.360761573048139, 0.171324492379170
};

// order 7
const Real xgauss7[7] = {
  -0.949107912342759, -0.741531185599394,
  -0.405845151377397, 0., 0.405845151377397,
  0.741531185599394, 0.949107912342759
};

const Real wgauss7[7] = {
  0.129484966168870, 0.279705391489277,
  0.381830050505119, 0.417959183673469, 0.381830050505119,
  0.279705391489277, 0.129484966168870
};

// order 8
const Real xgauss8[8] = {
  -0.960289856497536, -0.796666477413627,
  -0.525532409916329, -0.183434642495650,
  0.183434642495650, 0.525532409916329,
  0.796666477413627, 0.960289856497536
};

const Real wgauss8[8] = {
  0.101228536290376, 0.222381034453374,
  0.313706645877887, 0.362683783378362,
  0.362683783378362, 0.313706645877887,
  0.222381034453374, 0.101228536290376
};

// order 9
const Real xgauss9[9] = {
  -0.968160239507626, -0.836031107326636,
  -0.613371432700590, -0.324253423403809, 0.,
  0.324253423403809, 0.613371432700590,
  0.836031107326636, 0.968160239507626
};

const Real wgauss9[9] = {
  0.081274388361574, 0.180648160694857,
  0.260610696402935, 0.312347077040003, 0.330239355001260,
  0.312347077040003, 0.260610696402935,
  0.180648160694857, 0.081274388361574
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// For triangles only
// integration order max = 4 => number of points max = 7
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Local (reference) coordinates on the 1st edge of the triangle
const Real xg1[4][7] = {
  { 1./3., 0., 0., 0., 0., 0., 0. },
  { 0.5, 0., 0.5, 0., 0., 0., 0. },
  { 1./3., 0.6, 0.2, 0.2, 0., 0., 0. },
  { 1./3., (9. - 2. * sqrt(15.)) / 21., (6. + sqrt(15.)) / 21., (6. + sqrt(15.)) / 21.,
    (9. + 2. * sqrt(15.)) / 21., (6. - sqrt(15.)) / 21., (6. - sqrt(15.)) / 21. }
};

// Local (reference) coordinates on the 2nd edge of the triangle
const Real xg2[4][7] = {
  { 1./3., 0., 0., 0., 0., 0., 0. },
  { 0.5, 0.5, 0., 0., 0., 0., 0. },
  { 1./3., 0.2, 0.6, 0.2, 0., 0., 0. },
  { 1./3., (6. + sqrt(15.)) / 21., (9. - 2. * sqrt(15.)) / 21., (6. + sqrt(15.)) / 21., (6. - sqrt(15.)) / 21.,
    (9. + 2. * sqrt(15.)) / 21., (6. - sqrt(15.)) / 21. }
};

// Local (reference) coordinates on the 3rd edge of the triangle
/*
const Real xg3[4][7] = {
  { 1./3., 0., 0., 0., 0., 0., 0. },
  { 0., 0.5, 0.5, 0., 0., 0., 0. },
  { 1./3., 0.2, 0.2, 0.6, 0., 0., 0. },
  { 1./3., (6. + sqrt(15.)) / 21., (6. + sqrt(15.)) / 21., (9. - 2. * sqrt(15.)) / 21.,
    (6. - sqrt(15.)) / 21., (6. - sqrt(15.)) / 21., (9. + 2. * sqrt(15.)) / 21. }
};
*/

// Integration weights
const Real wg[4][7] = {
  { 0.5, 0., 0., 0., 0., 0., 0. },
  { 1. / 6., 1. / 6., 1. / 6., 0., 0., 0., 0. },
  { -27./96., 25./96., 25./96., 25./96., 0., 0., 0. },
  { 0.1125, (155. + sqrt(15.)) / 2400., (155. + sqrt(15.)) / 2400., (155. + sqrt(15.)) / 2400.,
    (155. - sqrt(15.)) / 2400., (155. - sqrt(15.)) / 2400., (155. - sqrt(15.)) / 2400. }
};

//  Correspondence between integration order & number of integration points
const Integer nptg[4] = { 1, 3, 4, 7 };

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// For tetrahedrons from 4 to 10 nodes only
// integration order = ? => number of integration points = 4
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Local (reference) coordinates along x axis
const Real xit[4] = { 0.58541020, 0.13819660, 0.13819660, 0.13819660 };

// Local (reference) coordinates along y axis
const Real yit[4] = { 0.13819660, 0.58541020, 0.13819660, 0.13819660 };

// Local (reference) coordinates along z axis
const Real zit[4] = { 0.13819660, 0.13819660, 0.58541020, 0.13819660 };

// Integration weight
const Real wgtetra = 1./24.;

const Integer maxnint = 9;

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
//  explicit CellFEMDispatcher(VariableNodeReal3& /*node_coords*/);

 public:

  void set_node_coords(VariableNodeReal3& /*node_coords*/);

//  Real getMeasure(const ItemWithNodes& /*cell*/);

//  Real3 getBarycenter(const ItemWithNodes& /*cell*/);

  Real getShapeFuncVal(const Int16& /*item_type*/, const Int32& /*inod*/, const Real3& /*ref coord*/);

  //  RealUniqueArray getShapeFuncDeriv(const Int16& /*item_type*/,const Integer& /*idir*/,const Real3& /*ref coord*/);
  Real3 getShapeFuncDeriv(const Int16& /*item_type*/, const Int32& /*inod*/, const Real3& /*ref coord*/);

//  Integer3 getOrientation(const ItemWithNodes& /*cell*/);
  RealUniqueArray getGaussData(const ItemWithNodes& item, const Integer3& nint, Int32& ngauss);

 private:

  std::function<Real(const ItemWithNodes& item, const VariableNodeReal3& n)> m_geomfunc[NB_BASIC_ITEM_TYPE];
  std::function<Real(const Int32& inod, const Real3& coord)> m_shapefunc[NB_BASIC_ITEM_TYPE];
  //  std::function<RealUniqueArray(const Integer& idir,const Real3& coord)> m_shapefuncderiv[NB_BASIC_ITEM_TYPE];
  std::function<Real3(const Int32& inod, const Real3& coord)> m_shapefuncderiv[NB_BASIC_ITEM_TYPE];
  std::function<Integer3(const ItemWithNodes& item, const VariableNodeReal3& n)> m_orientfunc[NB_BASIC_ITEM_TYPE];
  VariableNodeReal3 *m_node_coords{ nullptr};
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

extern Real getRefPosition(const Integer& indx, const Integer& ordre);
extern Real getWeight(const Integer& indx, const Integer& ordre);

extern Real3 LineRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real LineWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 TriRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real TriWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 QuadRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real QuadWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 HexaRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real HexaWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 TetraRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real TetraWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 PentaRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real PentaWeight(const Integer3& indices, const Integer3& ordre);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class GaussPointDispatcher
{
 public:

  GaussPointDispatcher();
//  GaussPointDispatcher(const Integer3& indices, const Integer3& int_order);

 public:

  void init_order(const Integer3& /*int_order*/);

  Real3 getRefPosition(const ItemWithNodes&, const Integer3& /*indices*/);

  Real getWeight(const ItemWithNodes&, const Integer3& /*indices*/);

 private:

  std::function<Real(const Integer3& indices, const Integer3& ordre)> m_weightfunc[NB_BASIC_ITEM_TYPE];
  std::function<Real3(const Integer3& indices, const Integer3& ordre)> m_refpositionfunc[NB_BASIC_ITEM_TYPE];
  Integer3 m_integ_order;
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
