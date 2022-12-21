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


using namespace Arcane;

extern Real	REL_PREC;
///////////////////////////////////////////////////////////////////////////////
// Tables of abscissa & integration weights for Gauss-Legendre quadrature

////////////////////////////////////////
// For edges & quadrangles only

// order 1
const Real xgauss1 = 0.;
const Real wgauss1 = 2.;

// order 2
const Real xgauss2[2] = {-1./sqrt(3.),1./sqrt(3.)};
const Real wgauss2[2] = {1.,1.};


// order 3
const Real xgauss3[3] = {-sqrt(0.6),0.,sqrt(0.6)};
const Real wgauss3[3] = {5./9.,8./9.,5./9.};

// order 4
const Real xgauss4[4] = {
-sqrt((3.+2.*sqrt(1.2))/7.),-sqrt((3.-2.*sqrt(1.2))/7.),
sqrt((3.-2.*sqrt(1.2))/7.),sqrt((3.+2.*sqrt(1.2))/7.)};

const Real wgauss4[4] = {
0.5-sqrt(5./6.)/6.,0.5+sqrt(5./6.)/6.,
0.5+sqrt(5./6.)/6.,0.5-sqrt(5./6.)/6.};

// order 5
const Real xgauss5[5] = {
-sqrt(245.+14.*sqrt(70.))/21.,-sqrt(245.-14.*sqrt(70.))/21.,0.,
sqrt(245.-14.*sqrt(70.))/21.,sqrt(245.+14.*sqrt(70.))/21.};

const Real wgauss5[5] = {
(322.-13*sqrt(70.))/900.,(322.+13*sqrt(70.))/900.,128./225.,
(322.+13*sqrt(70.))/900.,(322.-13*sqrt(70.))/900.};

// order 6
const Real xgauss6[6] = {
-0.932469514203152,-0.661209386466265,
-0.238619186083197,0.238619186083197,
0.661209386466265,0.932469514203152};

const Real wgauss6[6] = {
0.171324492379170,0.360761573048139,
0.467913934572691,0.467913934572691,
0.360761573048139,0.171324492379170};

// order 7
const Real xgauss7[7] = {
-0.949107912342759,-0.741531185599394,
-0.405845151377397,0.,0.405845151377397,
0.741531185599394,0.949107912342759};

const Real wgauss7[7] = {
0.129484966168870,0.279705391489277,
0.381830050505119,0.417959183673469,0.381830050505119,
0.279705391489277,0.129484966168870};

// order 8
const Real xgauss8[8] = {
-0.960289856497536,-0.796666477413627,
-0.525532409916329,-0.183434642495650,
0.183434642495650,0.525532409916329,
0.796666477413627,0.960289856497536};

const Real wgauss8[8] = {
0.101228536290376,0.222381034453374,
0.313706645877887,0.362683783378362,
0.362683783378362,0.313706645877887,
0.222381034453374,0.101228536290376};

// order 9
const Real xgauss9[9] = {
-0.968160239507626,-0.836031107326636,
-0.613371432700590,-0.324253423403809,0.,
0.324253423403809,0.613371432700590,
0.836031107326636,0.968160239507626};

const Real wgauss9[9] = {
0.081274388361574,0.180648160694857,
0.260610696402935,0.312347077040003,0.330239355001260,
0.312347077040003,0.260610696402935,
0.180648160694857,0.081274388361574};

////////////////////////////////////////
// For triangles only
// integration order max = 4 => number of points max = 7

// Local (reference) coordinates on the 1st edge of the triangle
const Real xg1[4][7] = {
	{1./3.,0.,0.,0.,0.,0.,0.},
	{0.5,0.,0.5,0.,0.,0.,0.},
	{1./3.,0.6,0.2,0.2,0.,0.,0.},
	{1./3.,(9.-2.*sqrt(15.))/21.,(6.+sqrt(15.))/21.,(6.+sqrt(15.))/21.,
			(9.+2.*sqrt(15.))/21.,(6.-sqrt(15.))/21.,(6.-sqrt(15.))/21.}
};

// Local (reference) coordinates on the 2nd edge of the triangle
const Real xg2[4][7] = {
	{1./3., 0., 0., 0., 0., 0., 0.},
	{0.5, 0.5, 0., 0., 0., 0., 0.},
	{1./3., 0.2, 0.6, 0.2, 0., 0., 0.},
	{1./3.,(6.+sqrt(15.))/21., (9.-2.*sqrt(15.))/21., (6.+sqrt(15.))/21., (6.-sqrt(15.))/21.,
		(9.+2.*sqrt(15.))/21., (6.-sqrt(15.))/21.}
};

// Local (reference) coordinates on the 3rd edge of the triangle
const Real xg3[4][7] = {
	{1./3., 0., 0., 0., 0., 0., 0.},
	{0., 0.5, 0.5, 0., 0., 0., 0.},
	{1./3., 0.2, 0.2, 0.6, 0., 0., 0.},
	{1./3.,(6.+sqrt(15.))/21.,(6.+sqrt(15.))/21.,(9.-2.*sqrt(15.))/21.,
	(6.-sqrt(15.))/21.,(6.-sqrt(15.))/21.,(9.+2.*sqrt(15.))/21.}
};

// Integration weights
const Real wg[4][7] = {
	{0.5,0.,0.,0.,0.,0.,0.},
	{1./6.,1./6.,1./6.,0.,0.,0.,0.},
	{-0.28125,0.2604166667,0.2604166667,0.2604166667,0.,0.,0.},
	{0.1125,(155.+sqrt(15.))/2400.,(155.+sqrt(15.))/2400.,(155.+sqrt(15.))/2400.,
	(155.-sqrt(15.))/2400.,(155.-sqrt(15.))/2400.,(155.-sqrt(15.))/2400.}
};

//  Correspondence between integration order & number of integration points
const Integer nptg[4] = {1,3,4,7};

///////////////////////////////////////////
// For tetrahedrons from 4 to 10 nodes only
// integration order = ? => number of integration points = 4

// Local (reference) coordinates along x axis
const Real xit[4] = {0.58541020,0.13819660,0.13819660,0.13819660};

// Local (reference) coordinates along y axis
const Real yit[4] = {0.13819660,0.58541020,0.13819660,0.13819660};

// Local (reference) coordinates along z axis
const Real zit[4] = {0.13819660,0.13819660,0.58541020,0.13819660};

// Integration weight
const Real wgtetra = 0.04166666666667;

const Integer maxnint = 9;

/////////////////////////////////////////////////////////////////////////////
// Functions useful for class CellFEMDispatcher

extern Real Line2Length(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Line2ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Line2ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Line2Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Real Line3Length(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Line3ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Line3ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Line3Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Real Tri3Surface(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Tri3ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Tri3ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Tri3Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Real Tri6Surface(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Tri6ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Tri6ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Tri6Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Real Quad4Surface(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Quad4ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Quad4ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Quad4Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Real Quad8Surface(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Quad8ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Quad8ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Quad8Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Real Hexa8Volume(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Hexa8ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Hexa8ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Hexa8Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Real Hexa20Volume(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Hexa20ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Hexa20ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Hexa20Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Real Tetra4Volume(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Tetra4ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Tetra4ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Tetra4Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Real Tetra10Volume(const ItemWithNodes& item,const VariableNodeReal3& n);
extern Real Tetra10ShapeFuncVal(const Integer& inod,const Real3& coord);
extern RealUniqueArray Tetra10ShapeFuncDeriv(const Integer& idir,const Real3& coord);
extern Integer3 Tetra10Orientation(const ItemWithNodes& item,const VariableNodeReal3& n);

extern Integer getGeomDimension(const ItemWithNodes& item);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class CellFEMDispatcher
{
 public:
	explicit CellFEMDispatcher(VariableNodeReal3& /*node_coords*/);

 public:
  Real getGeomSize(const ItemWithNodes&);

  Real3 getBarycenter(const ItemWithNodes&);

  Real getShapeFuncVal(const ItemWithNodes&,const Integer& /*inod*/,const Real3& /*coord*/);

  RealUniqueArray getShapeFuncDeriv(const ItemWithNodes&,const Integer& /*idir*/,const Real3& /*coord*/);

  Integer3 getOrientation(const ItemWithNodes&);

// Computing the Jacobian of this finite-element at an integration point (e.g. Gauss) given by its local coordinates in the element
  Real computeJacobian(const ItemWithNodes&,const Real3& /*coord*/);

 private:
  std::function<Real(const ItemWithNodes& item,const VariableNodeReal3& n)> m_geomfunc[NB_BASIC_ITEM_TYPE];
  std::function<Real(const Integer& inod,const Real3& coord)> m_shapefunc[NB_BASIC_ITEM_TYPE];
  std::function<RealUniqueArray(const Integer& idir,const Real3& coord)> m_shapefuncderiv[NB_BASIC_ITEM_TYPE];
  std::function<Integer3(const ItemWithNodes& item,const VariableNodeReal3& n)> m_orientfunc[NB_BASIC_ITEM_TYPE];
  VariableNodeReal3 m_node_coords;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/////////////////////////////////////////////////////////////////////////////
// Functions useful for class GaussPointDispatcher

extern Real getRefPosition(const Integer& indx,const Integer& ordre);
extern Real getWeight(const Integer& indx,const Integer& ordre);

extern Real3 LineRefPosition(const Integer3& indices,const Integer3& ordre);
extern Real LineWeight(const Integer3& indices,const Integer3& ordre);

extern Real3 TriRefPosition(const Integer3& indices,const Integer3& ordre);
extern Real TriWeight(const Integer3& indices,const Integer3& ordre);

extern Real3 QuadRefPosition(const Integer3& indices,const Integer3& ordre);
extern Real QuadWeight(const Integer3& indices,const Integer3& ordre);

extern Real3 HexaRefPosition(const Integer3& indices,const Integer3& ordre);
extern Real HexaWeight(const Integer3& indices,const Integer3& ordre);

extern Real3 TetraRefPosition(const Integer3& indices,const Integer3& ordre);
extern Real TetraWeight(const Integer3& indices,const Integer3& ordre);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class GaussPointDispatcher
{
 public:
	GaussPointDispatcher(const Integer3& indices,const Integer3& int_order);

 public:
  Real3 getRefPosition(const ItemWithNodes&);

  Real getWeight(const ItemWithNodes&);

 private:
  std::function<Real(const Integer3& indices,const Integer3& ordre)> m_weightfunc[NB_BASIC_ITEM_TYPE];
  std::function<Real3(const Integer3& indices,const Integer3& ordre)> m_refpositionfunc[NB_BASIC_ITEM_TYPE];
  Integer3	m_indices;
  Integer3	m_integ_order;
};

#endif // PASSMO_UTILFEM_H_
