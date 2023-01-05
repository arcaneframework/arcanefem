/*
 * PASSMO : Performant Assessment for Seismic Site Modelling
 *
 * Definition of classes to implement finite-element cells and related shape functions
 * and their derivatives
 *
 * utilFEM.cc: declarations
 *
 *  Created on: December 2019
 *      Author: E. Foerster
 */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "arcane/MathUtils.h"
#include <arcane/utils/NumArray.h>
#include "arcane/utils/ArgumentException.h"
#include <arcane/IParallelMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/geometry/IGeometry.h>
#include <arcane/VariableTypes.h>

#include "Integer3std.h"
#include "utilFEM.h"

using namespace Arcane;
using namespace Arcane::FemUtils;
/////////////////////////////////////////////////////////////////////////////
// class CellFEMDispatcher: constructor

CellFEMDispatcher::CellFEMDispatcher(VariableNodeReal3& node_coords): m_node_coords(node_coords)
{
        // Setting to null default value
        for(int i = 0; i < NB_BASIC_ITEM_TYPE; ++i )
        {
            m_geomfunc[i] = nullptr;
            m_shapefunc[i] = nullptr;
            m_shapefuncderiv[i] = nullptr;
            m_orientfunc[i] = nullptr;
        }

        // Gives functions to compute geometric properties (length, surface or volume)
        // Linear elements
        m_geomfunc[IT_Line2] = Line2Length;
        m_geomfunc[IT_Triangle3] = Tri3Surface;
        m_geomfunc[IT_Quad4] = Quad4Surface;
        m_geomfunc[IT_Tetraedron4] = Tetra4Volume;
        m_geomfunc[IT_Hexaedron8] = Hexa8Volume;

        // Quadratic elements
        m_geomfunc[IT_Line3] = Line3Length;
        m_geomfunc[IT_Triangle6] = Tri6Surface;
        m_geomfunc[IT_Quad8] = Quad8Surface;
        m_geomfunc[IT_Tetraedron10] = Tetra10Volume;
        m_geomfunc[IT_Hexaedron20] = Hexa20Volume;

        // Gives functions to compute shape function value in finite-element reference coordinate system
        // Linear elements
        m_shapefunc[IT_Line2] = Line2ShapeFuncVal;
        m_shapefunc[IT_Triangle3] = Tri3ShapeFuncVal;
        m_shapefunc[IT_Quad4] = Quad4ShapeFuncVal;
        m_shapefunc[IT_Tetraedron4] = Tetra4ShapeFuncVal;
        m_shapefunc[IT_Hexaedron8] = Hexa8ShapeFuncVal;

        // Quadratic elements
        m_shapefunc[IT_Line3] = Line3ShapeFuncVal;
        m_shapefunc[IT_Triangle6] = Tri6ShapeFuncVal;
        m_shapefunc[IT_Quad8] = Quad8ShapeFuncVal;
        m_shapefunc[IT_Tetraedron10] = Tetra10ShapeFuncVal;
        m_shapefunc[IT_Hexaedron20] = Hexa20ShapeFuncVal;

        // Gives functions to compute shape function derivate vector at all nodes of a finite-element
        // along a local direction (in reference coordinate system)
        // Linear elements
        m_shapefuncderiv[IT_Line2] = Line2ShapeFuncDeriv;
        m_shapefuncderiv[IT_Triangle3] = Tri3ShapeFuncDeriv;
        m_shapefuncderiv[IT_Quad4] = Quad4ShapeFuncDeriv;
        m_shapefuncderiv[IT_Tetraedron4] = Tetra4ShapeFuncDeriv;
        m_shapefuncderiv[IT_Hexaedron8] = Hexa8ShapeFuncDeriv;

        // Quadratic elements
        m_shapefuncderiv[IT_Line3] = Line3ShapeFuncDeriv;
        m_shapefuncderiv[IT_Triangle6] = Tri6ShapeFuncDeriv;
        m_shapefuncderiv[IT_Quad8] = Quad8ShapeFuncDeriv;
        m_shapefuncderiv[IT_Tetraedron10] = Tetra10ShapeFuncDeriv;
        m_shapefuncderiv[IT_Hexaedron20] = Hexa20ShapeFuncDeriv;

        // Gives functions to compute orientation vector of the finite-element in global axes (x, y or z)
        // Linear elements
        m_orientfunc[IT_Line2] = Line2Orientation;
        m_orientfunc[IT_Triangle3] = Tri3Orientation;
        m_orientfunc[IT_Quad4] = Quad4Orientation;
        m_orientfunc[IT_Tetraedron4] = Tetra4Orientation;
        m_orientfunc[IT_Hexaedron8] = Hexa8Orientation;

        // Quadratic elements
        m_orientfunc[IT_Line3] = Line3Orientation;
        m_orientfunc[IT_Triangle6] = Tri6Orientation;
        m_orientfunc[IT_Quad8] = Quad8Orientation;
        m_orientfunc[IT_Tetraedron10] = Tetra10Orientation;
        m_orientfunc[IT_Hexaedron20] = Hexa20Orientation;
}

/////////////////////////////////////////////////////////////////////////////
// class CellFEMDispatcher: implementation methods

// ! Computes the Jacobian = geometric size (length, surface or volume) of any finite-element
// aligned with the coordinate axes (otherwise, jacobian to be computed by the
// general Jacobian matrix)
Real CellFEMDispatcher::getJacobian(const ItemWithNodes& item)
{
    Int32 item_type = item.type();
    auto f = m_geomfunc[item_type];
    if (f!=nullptr)
      return f(item,m_node_coords);
    return (-1.0);
}

// ! Computes the barycenter of a Cell
Real3 CellFEMDispatcher::getBarycenter(const ItemWithNodes& item)
{

 	Real3 sum = Real3::zero();
  	Integer nnod = item.nbNode();

  	for (Integer i = 0; i < nnod; i++) {
  		sum += m_node_coords[item.node(i)];
  	}

  	if (nnod) sum /= nnod;
  	return sum;
}

// ! Computes the value of the nodal shape functions of a given FE element at a given Gauss point (allowing to the element type)
Real CellFEMDispatcher::getShapeFuncVal(const Int16& item_type,const Integer& inod,const Real3& coord)
{
    auto f = m_shapefunc[item_type];
    if (f!=nullptr)
      return f(inod,coord);
    return 0.;
}

// ! Gives the values of the derivatives for the function of node inod for a given FE element at a given Gauss point
// ! Derivation is performed along 3 directions x, y,z
Real3 CellFEMDispatcher::getShapeFuncDeriv(const Int16& item_type,const Integer& inod,const Real3& ref_coord){
    auto f = m_shapefuncderiv[item_type];
    if (f!=nullptr)
      return f(inod,ref_coord);
    return {};
}

// ! Computes the cartesian directions of a finite-element (used for derivation of shape functions)
// ! For instance, assuming the returned Integer3 vector=idir, idir = (0,-1,-1) means there will derivation along x only
Integer3 CellFEMDispatcher::getOrientation(const ItemWithNodes& cell)
{
     Int32 item_type = cell.type();
     auto f = m_orientfunc[item_type];
     if (f!=nullptr)
       return f(cell,m_node_coords);
     return Integer3::zero();
}

/////////////////////////////////////////////////////////////////////////////
// Functions useful for class CellFEMDispatcher

/*---------------------------------------------------------------------------*/
// Line2: linear edge finite-element
//
// 1           0
// o-----------o---> x
//-1           1
// direct : 0->1 (local numbering)
/*---------------------------------------------------------------------------*/

Real Line2Length(const ItemWithNodes& item,const VariableNodeReal3& n){
  const Real3& n0 = n[item.node(0)];
  const Real3& n1 = n[item.node(1)];
  return (n0-n1).normL2();
}

Real Line2ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
	assert(inod >= 0 && inod < 2);
#endif

	Real r = ref_coord[0];
	if (!inod) return (0.5*(1 + r));
	return (0.5*(1 - r));
}

Real3 Line2ShapeFuncDeriv(const Integer& inod,const Real3&){
  if (!inod) return { 0.5,0.,0. };
  return { -0.5,0.,0. };
}

Integer3 Line2Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
	const Real3& n0 = n[item.node(0)];
	const Real3& n1 = n[item.node(1)];
	Integer3 dir;
	Real3 pt = n0-n1;

	Integer j = -1;
	for (Integer i = 0; i < 3; i++)	{
		if (fabs(pt[i]) != 0.) dir[i] = ++j;
	}
	return dir;
}

/*---------------------------------------------------------------------------*/
// Line3: quadratic edge finite-element
//
// 1     2     0
// o-----o-----o---> x
//-1           1
// direct : 0->1 (local numbering)
/*---------------------------------------------------------------------------*/

Real Line3Length(const ItemWithNodes& item,const VariableNodeReal3& n){
    return Line2Length(item,n);
}

Real Line3ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
        assert(inod >= 0 && inod < 3);
#endif

    Real ri = ref_coord[0];
    if (inod == 1) ri *= -1;

    if (inod < 2) return 0.5*ri*(1 + ri); // nodes 0 or 1
    return (1 - ri*ri); // middle node
}

Real3 Line3ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
    if (!inod) return {0.5 + ref_coord[0], 0.,0.};
    if (inod == 1) return {-0.5 + ref_coord[0], 0.,0.};
    return {-2.*ref_coord[0], 0.,0.};
}

Integer3 Line3Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
    return Line2Orientation(item,n);
}

/*---------------------------------------------------------------------------*/
// Tri3: linear triangle finite-element
//
//        0 o
//         . .
//        .   .
//       .     .
//      .       .
//     .         .
//    .           .
//   o-------------o
//  1               2
// direct : 0,1,2 (local numbering)
/*---------------------------------------------------------------------------*/

Real Tri3Surface(const ItemWithNodes& item,const VariableNodeReal3& n){
  Real3 n0 = n[item.node(0)];
  Real3 n1 = n[item.node(1)];
  Real3 n2 = n[item.node(2)];

  Real x1 = n1.x - n0.x;
  Real y1 = n1.y - n0.y;
  Real x2 = n2.x - n1.x;
  Real y2 = n2.y - n1.y;

  return x1 * y2 - y1 * x2;
}

Real Tri3ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
	assert(inod >= 0 && inod < 3);
#endif

	return ref_coord[inod];
}

Real3 Tri3ShapeFuncDeriv(const Integer& inod,const Real3&){
  if (!inod) return {1., 0.,0.};
  if (inod==1) return {0., 1., 0.};
  return {-1., -1.,0.};
}

Integer3 Tri3Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
	Real3 p0 = n[item.node(0)];
	Real3 u1 = p0 - n[item.node(1)];
	Real3 u2 = p0 - n[item.node(2)];
	Integer3 dir;

	for (Integer i = 0, j = -1; i < 3; i++)
	{
		Real3 e;
		e[i] = 1.;
		if (math::dot(u1,e) != 0. || math::dot(u2,e) != 0.) dir[i] = ++j;
	}
	return dir;
}

/*---------------------------------------------------------------------------*/
// Tri6: quadratic triangle finite-element
//
//        0 o
//         . .
//        .   .
//      3o     o5
//      .       .
//     .         .
//    .           .
//   o------o4-----o
//  1               2
// direct : 0,1,...,5 (local numbering)
/*---------------------------------------------------------------------------*/

Real Tri6Surface(const ItemWithNodes& item,const VariableNodeReal3& n){
    return Tri3Surface(item,n);
}

Real Tri6ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
    assert(inod >= 0 && inod < 6);
#endif

    if (inod < 3) return ref_coord[inod];

    auto	wi = 0.,ri = ref_coord[0],si = ref_coord[1],ti = ref_coord[2];
    switch(inod) {
        default: break;

        case 3:	wi = 4.*ri*si; break;

        case 4:	wi = 4.*si*ti; break;

        case 5:	wi = 4.*ri*ti; break;
    }
    return wi;
}

Real3 Tri6ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
  if (!inod) return {1., 0.,0.};
  if (inod==1) return {0., 1., 0.};
  if (inod == 2) return {-1., -1.,0.};

  auto	ri = ref_coord[0],si = ref_coord[1],ti = ref_coord[2];
  if (inod == 3) return {4.*si, 4.*ri,0.};
  if (inod == 4) return {-4.*si, 4.*(ti - si),0.};
  return {4.*(ti - ri), -4.*ri,0.};
}

Integer3 Tri6Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
    return Tri3Orientation(item,n);
}

/*---------------------------------------------------------------------------*/
// Quad4: linear quadrangle finite-element
//
//        ^y
//        1|
//  1o-----|-----o0
//   |     |     |
//   |     |     |quadratic
//   |     |     |
//   |-----|---------> x
// -1|     |     |1
//   |     |     |
//   |     |     |
//  2o-----|-----o3
//         -1
// direct : 0,1,2,3 (local numbering)
/*---------------------------------------------------------------------------*/

Real Quad4Surface(const ItemWithNodes& item,const VariableNodeReal3& n){
  Real3 n0 = n[item.node(0)];
  Real3 n1 = n[item.node(1)];
  Real3 n2 = n[item.node(2)];
  Real3 n3 = n[item.node(3)];

  Real x1 = n1.x - n0.x;
  Real y1 = n1.y - n0.y;
  Real x2 = n2.x - n1.x;
  Real y2 = n2.y - n1.y;
  Real surface = x1 * y2 - y1 * x2;

  x1 = n2.x - n0.x;
  y1 = n2.y - n0.y;
  x2 = n3.x - n2.x;
  y2 = n3.y - n2.y;

  surface += x1 * y2 - y1 * x2;
  return surface;
}

Real Quad4ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
	assert(inod >= 0 && inod < 4);
#endif

	Real	ri = ref_coord[0],si = ref_coord[1]; // default is first node (index 0)

	switch(inod){
		default: break;
		case 1:	ri *= -1; break;
		case 2:	ri *= -1; si *= -1; break;
		case 3:	si *= -1; break;
	}
	return ( 0.25*(1 + ri)*(1 + si) );
}

Real3 Quad4ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
	assert (idir >= 0 && idir < 2);
#endif

  auto ri{ref_coord[0]}, si{ref_coord[1]};
  if (!inod) return {0.25*(1 + si),0.25*(1 + ri),0.};
  if (inod == 1) return {-0.25*(1 + si),0.25*(1 - ri),0.};
  if (inod == 2) return {-0.25*(1 - si),-0.25*(1 - ri),0.};
  return {0.25*(1 - si),-0.25*(1 + ri),0.};
}

Integer3 Quad4Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
	const Real3& p0 = n[item.node(0)];
	Real3 u1 = p0 - n[item.node(1)];
	Real3 u2 = p0 - n[item.node(2)];
	Integer3 dir;

	for (Integer i = 0,j = -1; i < 3; i++)	{
		Real3 e;
		e[i] = 1.;
		if (math::dot(u1,e) != 0. || math::dot(u2,e) != 0.) dir[i] = ++j;
	}
	return dir;
}

/*---------------------------------------------------------------------------*/
// Quad8: quadratic quadrangle finite-element
//
//        ^y
//        1|
//  1o-----o-----o0
//   |     |4    |
//   |     |     |
//   |     |     |
//   o5----|----7o---> x
// -1|     |     |1
//   |     |     |
//   |     6     |
//  2o-----o-----o3
//         -1
// direct : 0,1,2,3,...,7 (local numbering)
/*---------------------------------------------------------------------------*/

Real Quad8Surface(const ItemWithNodes& item,const VariableNodeReal3& n){
    return Quad4Surface(item,n);
}


Real Quad8ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
    assert(inod >= 0 && inod < 8);
#endif

    auto	wi{0.},
          r = ref_coord[0],s = ref_coord[1],
          ri{ r },si{ s }; // default is first node (index 0)

    switch(inod){
        default: break;

        case 1:	ri *= -1; break;

        case 2:	ri *= -1; si *= -1; break;

        case 3:	si *= -1; break;

        case 4:	ri = 0; break;

        case 5:	ri *= -1; si = 0; break;

        case 6:	ri = 0; si *= -1; break;

        case 7:	si = 0; break;
    }

    if (inod < 4) wi = 0.25*(1 + ri)*(1 + si)*(ri + si - 1);
    else if (ri == 0.) wi = 0.5*(1 - r*r)*(1 + si);
    else if (si == 0.) wi = 0.5*(1 - s*s)*(1 + ri);

    return wi;
}

Real3 Quad8ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

    RealUniqueArray r(8),s(8);
    r[0] = 1; s[0] = 1;
    r[1] = -1; s[1] = 1;
    r[2] = -1; s[2] = -1;
    r[3] = 1; s[3] = -1;
    r[4] = 0; s[4] = 1;
    r[5] = -1; s[5] = 0;
    r[6] = 0; s[6] = -1;
    r[7] = 1; s[7] = 0;
    auto wi = Quad4ShapeFuncVal(inod,ref_coord),
         ri = r[inod]*ref_coord[0],
         si = s[inod]*ref_coord[1],
         r2 = ref_coord[0]*ref_coord[0],
         s2 = ref_coord[1]*ref_coord[1];

    // derivatives on corner nodes
    if (inod < 4) return {r[inod]*(0.25*(1 + si)*(ri + si - 1) + wi),s[inod]*(0.25*(1 + ri)*(ri + si - 1) + wi), 0.};

    // derivatives on middle nodes
    if (r[inod] != 0.) return {-ref_coord[0]*(1 + si), 0.5*si*(1 - r2),0.};
    return {0.5*ri*(1 - s2), -ref_coord[1]*(1 + ri),0.};
}

Integer3 Quad8Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
    return Quad4Orientation(item,n);
}

/*---------------------------------------------------------------------------*/
// Hexa8: linear hexaedron finite-element
//
//         1-----------0
//        /|          /|y=1
//      /  |  z=1   /  |   z
//    2----|-- ---3    |   |    y
//    |    |      |    |   |  /
//    |    |      |    |   |---->x
//    |    |      |    |
//x=-1|    |      |x=1 |
//    |    5------|----4
//    |  /        |  /
//    |/          |/
//    6-----------7
//      y=-1,z=-1
// direct : 0,1,2,3,...,7 (local numbering)
/*---------------------------------------------------------------------------*/

Real Hexa8Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
  const Real3& n0 = n[item.node(0)];
  const Real3& n1 = n[item.node(1)];
  const Real3& n2 = n[item.node(2)];
  const Real3& n3 = n[item.node(3)];
  const Real3& n4 = n[item.node(4)];
  const Real3& n5 = n[item.node(5)];
  const Real3& n6 = n[item.node(6)];
  const Real3& n7 = n[item.node(7)];

  Real v1 = math::matDet((n6 - n1) + (n7 - n0), n6 - n3, n2 - n0);
  Real v2 = math::matDet(n7 - n0, (n6 - n3) + (n5 - n0), n6 - n4);
  Real v3 = math::matDet(n6 - n1, n5 - n0, (n6 - n4) + (n2 - n0));
  return (v1 + v2 + v3) / 12.0;
}

Real Hexa8ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
	assert(inod >= 0 && inod < 8);
#endif

	auto	ri = ref_coord[0],si = ref_coord[1],ti = ref_coord[2]; // default is first node (index 0)

	switch(inod){
		default: break;
		case 1:	ri *= -1; break;
		case 2:	ri *= -1; si *= -1; break;
		case 3:	si *= -1; break;
		case 4:	ti *= -1; break;
		case 5:	ri *= -1; ti *= -1; break;
		case 6:	ri *= -1; si *= -1; ti *= -1; break;
		case 7:	si *= -1; ti *= -1; break;
	}
	return (0.125*(1 + ri)*(1 + si)*(1 + ti));
}

Real3 Hexa8ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
	Real	rp = 1.+ ref_coord[0],sp = 1. + ref_coord[1],tp = 1. + ref_coord[2],
			rm = 1.- ref_coord[0],sm = 1. - ref_coord[1],tm = 1. - ref_coord[2];

  if (!inod) return {0.125*sp*tp,0.125*rp*tp, 0.125*rp*sp};
  if (inod == 1) return {-0.125*sp*tp,0.125*rm*tp, 0.125*rm*sp};
  if (inod == 2) return {-0.125*sm*tp,-0.125*rm*tp, 0.125*rm*sm};
  if (inod == 3) return {0.125*sm*tp,-0.125*rp*tp, 0.125*rp*sm};
  if (inod == 4) return {0.125*sp*tm,-0.125*rp*tm, -0.125*rp*sp};
  if (inod == 5) return {-0.125*sp*tm,0.125*rm*tm, -0.125*rm*sp};
  if (inod == 6) return {-0.125*sm*tm,-0.125*rm*tm, -0.125*rm*sm};
  return {0.125*sm*tm,-0.125*rp*tm, -0.125*rp*sm};
}

Integer3 Hexa8Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
	return {1,1,1};
}

/*---------------------------------------------------------------------------*/
// Hexa20: quadratic hexaedron finite-element
//
//         1-----------0
//        /|          /|y=1
//      /  |  z=1   /  |   z
//    2----|-- ---3    |   |    y
//    |    |      |    |   |  /
//    |    |      |    |   |---->x
//    |    |      |    |
//x=-1|    |      |x=1 |
//    |    5------|----4
//    |  /        |  /
//    |/          |/
//    6-----------7
//      y=-1,z=-1
// direct : 0,1,2,3,...,7 (local numbering)
/*---------------------------------------------------------------------------*/

Real Hexa20Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
    return Hexa8Volume(item, n);
}

Real Hexa20ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
    assert(inod >= 0 && inod < 20);
#endif

    auto	wi = 0.,ri = ref_coord[0],si = ref_coord[1],ti = ref_coord[2],
            rm = 1. - ri,rp = 1. + ri,
            sm = 1. - si,sp = 1. + si,
            tm = 1. - ti,tp = 1. + ti,
            r2 = 1. - ri*ri,s2 = 1. - si*si,t2 = 1. - ti*ti;

    if (inod < 8) return Hexa8ShapeFuncVal(inod,ref_coord);

    switch(inod){
        default: break;
        case 8: wi = 0.25*r2*sp*tp; break;
        case 9: wi = 0.25*rm*s2*tp; break;
        case 10: wi = 0.25*r2*sm*tp; break;
        case 11: wi = 0.25*rp*s2*tp; break;
        case 12: wi = 0.25*r2*sp*tm; break;
        case 13: wi = 0.25*rm*s2*tm; break;
        case 14: wi = 0.25*r2*sm*tm; break;
        case 15: wi = 0.25*rp*s2*tm; break;
        case 16: wi = 0.25*rp*sp*t2; break;
        case 17: wi = 0.25*rm*sp*t2; break;
        case 18: wi = 0.25*rm*sm*t2; break;
        case 19: wi = 0.25*rp*sm*t2; break;
    }
    return wi;
}

Real3 Hexa20ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

  if (inod < 8) return Hexa8ShapeFuncDeriv(inod,ref_coord);

  auto ri = ref_coord[0],si = ref_coord[1],ti = ref_coord[2],
       rm = 1. - ri,rp = 1. + ri,
       sm = 1. - si,sp = 1. + si,
       tm = 1. - ti,tp = 1. + ti,
       r2 = 1. - ri*ri,
       s2 = 1. - si*si,
       t2 = 1. - ti*ti;

  if (inod == 8) return { -0.5*ri*sp*tp,0.25*r2*tp,0.25*r2*sp };
  if (inod == 9) return { -0.25*s2*tp,-0.5*si*rm*tp,0.25*s2*rm };
  if (inod == 10) return { -0.5*ri*sm*tp,-0.25*r2*tp,0.25*r2*sm };
  if (inod == 11) return { 0.25*s2*tp,-0.5*si*rp*tp,0.25*s2*rp };
  if (inod == 12) return { -0.5*ri*sp*tm,0.25*r2*tm,-0.25*r2*sp };
  if (inod == 13) return { -0.25*s2*tm,-0.5*si*rm*tm,-0.25*s2*rm };
  if (inod == 14) return { -0.5*ri*sm*tm,-0.25*r2*tm,-0.25*r2*sm };
  if (inod == 15) return { 0.25*s2*tm,-0.5*si*rp*tm,-0.25*s2*rp };
  if (inod == 16) return { 0.25*t2*sp,0.25*t2*rp,-0.5*ti*rp*sp };
  if (inod == 17) return { -0.25*t2*sp,0.25*t2*rm,-0.5*ti*rm*sp };
  if (inod == 18) return { -0.25*t2*sm,-0.25*t2*rm,-0.5*ti*rm*sm };
  if (inod == 19) return { 0.25*t2*sm,-0.25*t2*rp,-0.5*ti*rp*sm };
}

Integer3 Hexa20Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
    return {1,1,1};
}

/*---------------------------------------------------------------------------*/
// Tetra4: linear tetrahedral finite-element
//
//    (0,0,1)                     0
//       .                        *.*
//       .                        * . *
//       .                        *  .  *
//       .                        *   .   *
//       Z   (0,1,0)              *    .    *
//       .    .                   *     2     *
//       .   .                    *   .    .    *
//       .  Y                     *  .        .   *
//       . .                      * .            .  *
//       ..           (1,0,0)     *.                . *
//       --------X------>         3********************1
//
// direct : 0,1,2,3 (local numbering)
/*---------------------------------------------------------------------------*/

Real Tetra4Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
  const Real3& n0 = n[item.node(0)];
  const Real3& n1 = n[item.node(1)];
  const Real3& n2 = n[item.node(2)];
  const Real3& n3 = n[item.node(3)];
  return math::matDet(n1 - n0, n2 - n0, n3 - n0)/6.;
}

Real Tetra4ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
	assert(inod >= 0 && inod < 4);
#endif

	auto	ri = ref_coord[0],si = ref_coord[1],ti = ref_coord[2]; // default is first node (index 0)

	switch(inod){
		default: break;
		case 1:	return ri;
		case 2:	return si;
		case 3:	return (1. - ri - si - ti);
	}
	return ti;
}

Real3 Tetra4ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

  if (!inod) return {0.,0.,1.};
  if (inod == 1) return {1.,0.,0.};
  if (inod == 2) return {0.,1.,0.};
  return {-1.,-1.,-1.};
}

Integer3 Tetra4Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
	return {1,1,1};
}

/*---------------------------------------------------------------------------*/
// Tetra10: quadratic tetrahedral finite-element
//
//
//    (0,0,1)                     x0
//       .                        *.*
//       .                        * . *
//       .                        *  .  *
//       .                        *   o8  *
//       Z   (0,1,0)              *    .    o7
//       .    .                  9o     x2    *
//       .   .                    *   .    .    *
//       .  Y                     *  o5      4o   *
//       . .                      * .            .  *
//       ..           (1,0,0)     *.                . *
//       --------X------>        3x*******o6***********x1
//
// direct : 0,1,2,...,9 (local numbering)
/*---------------------------------------------------------------------------*/

Real Tetra10Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
    return Tetra4Volume(item,n);
}

Real Tetra10ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
    assert(inod >= 0 && inod < 10);
#endif

    auto x = ref_coord[0],y = ref_coord[1],z = ref_coord[2],
         t = 1. - x - y - z,
         wi{0.};

    switch(inod){
        default: break;

        // Corner nodes
        case 0:	wi = (1. - 2*t - 2*x - 2*y)*z;break;
        case 1:	wi = (1. - 2*t - 2*y - 2*z)*x;break;
        case 2:	wi = (1. - 2*x - 2*t - 2*z)*y;break;
        case 3:	wi = (1. - 2*x - 2*y - 2*z)*t;break;

        // Middle nodes
        case 4:	wi = 4*x*y;break;
        case 5:	wi = 4*y*t;break;
        case 6:	wi = 4*x*t;break;
        case 7:	wi = 4*z*x;break;
        case 8:	wi = 4*z*y;break;
        case 9:	wi = 4*z*t;break;
    }
    return wi;
}

Real3 Tetra10ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
    auto x = ref_coord[0],y = ref_coord[1],z = ref_coord[2],
         t = 1. - x - y - z;

    // Corner nodes
    if (!inod) return {0.,0.,1. - 2*t - 2*x - 2*y + 2*z};
    if (inod == 1) return {1. - 2*t - 2*y - 2*z + 2*x,0.,0.};
    if (inod == 2) return {0.,1. - 2*x - 2*t - 2*z + 2*y,0.};
    if (inod == 3) return {-1. - 2*t + 2*x + 2*y + 2*z,-1. - 2*t+ 2*x + 2*y + 2*z,-1. - 2*t + 2*x + 2*y + 2*z};

    // Middle nodes
    if (inod == 4) return {4*y,4*x,0.};
    if (inod == 5) return {-4*y,4*(t - y),-4*y};
    if (inod == 6) return {4*(t - x),-4*x,-4*x};
    if (inod == 7) return {4*z,0.,4*x};
    if (inod == 8) return {0.,4*z,4*y};
    return {-4*z,-4*z,4*(t - z)};
}

Integer3 Tetra10Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
    return {1,1,1};
}

/*---------------------------------------------------------------------------*/

Integer getGeomDimension(const ItemWithNodes& item){
    Int32 item_type = item.type();
    Integer dim = 1; // default geometric dimension is 1D (Line2 and Line3 finite-elements)

    switch(item_type) {

        // 2D elements
        case IT_Triangle3:
        case IT_Quad4: dim = 2; break;

/*        case IT_Triangle6:
        case IT_Quad8: dim = 2; break;*/

            // 3D elements
        case IT_Tetraedron4:
        case IT_Hexaedron8: dim = 3; break;

/*        case IT_Tetraedron10:
        case IT_Hexaedron20: dim = 3; break;*/

        default: break;

    }
    return dim;
}

/*---------------------------------------------------------------------------*/
/////////////////////////////////////////////////////////////////////////////
// class GaussPointDispatcher: construction methods

GaussPointDispatcher::GaussPointDispatcher(const Integer3& indices,const Integer3& int_order):
        m_indices(indices),m_integ_order(int_order){
    // Setting to null default value
    for(Integer i = 0; i < NB_BASIC_ITEM_TYPE; ++i ){
        m_weightfunc[i] = nullptr;
        m_refpositionfunc[i] = nullptr;
    }

    // Gives functions to compute weight for a Gauss Point in finite-element reference coordinate system
    // Linear elements
    m_weightfunc[IT_Line2] = LineWeight;
    m_weightfunc[IT_Triangle3] = TriWeight;
    m_weightfunc[IT_Quad4] = QuadWeight;
    m_weightfunc[IT_Tetraedron4] = TetraWeight;
    m_weightfunc[IT_Hexaedron8] = HexaWeight;

    // Quadratic elements
    m_weightfunc[IT_Line3] = LineWeight;
    m_weightfunc[IT_Triangle6] = TriWeight;
    m_weightfunc[IT_Quad8] = QuadWeight;
    m_weightfunc[IT_Tetraedron10] = TetraWeight;
    m_weightfunc[IT_Hexaedron20] = HexaWeight;

    // Gives functions to compute position of a Gauss Point in finite-element reference coordinate system
    // Linear elements
    m_refpositionfunc[IT_Line2] = LineRefPosition;
    m_refpositionfunc[IT_Triangle3] = TriRefPosition;
    m_refpositionfunc[IT_Quad4] = QuadRefPosition;
    m_refpositionfunc[IT_Tetraedron4] = TetraRefPosition;
    m_refpositionfunc[IT_Hexaedron8] = HexaRefPosition;

    // Quadratic elements
    m_refpositionfunc[IT_Line3] = LineRefPosition;
    m_refpositionfunc[IT_Triangle6] = TriRefPosition;
    m_refpositionfunc[IT_Quad8] = QuadRefPosition;
    m_refpositionfunc[IT_Tetraedron10] = TetraRefPosition;
    m_refpositionfunc[IT_Hexaedron20] = HexaRefPosition;

}

/////////////////////////////////////////////////////////////////////////////
// class GaussPointDispatcher: implementation methods

Real3 GaussPointDispatcher::getRefPosition(const ItemWithNodes& item){
	  Int32 item_type = item.type();
	  auto f = m_refpositionfunc[item_type];
	  if (f!=nullptr)
		  return f(m_indices,m_integ_order);
	  return Real3::zero();
}

Real GaussPointDispatcher::getWeight(const ItemWithNodes& item){
	  Int32 item_type = item.type();
	  auto f = m_weightfunc[item_type];
	  if (f!=nullptr)
		  return f(m_indices,m_integ_order);
	  return 0.;
}

/*---------------------------------------------------------------------------*/
/////////////////////////////////////////////////////////////////////////////
// Functions useful for class GaussPointDispatcher

Real getRefPosition(const Integer& indx,const Integer& ordre){
		Real x = xgauss1; // default is order 1

		switch(ordre){
			case 2: x = xgauss2[indx]; break;
			case 3: x = xgauss3[indx]; break;
			case 4: x = xgauss4[indx]; break;
			case 5: x = xgauss5[indx]; break;
			case 6: x = xgauss6[indx]; break;
			case 7: x = xgauss7[indx]; break;
			case 8: x = xgauss8[indx]; break;
			case 9: x = xgauss9[indx]; break;
			default: break;
		}
		return x;
}

Real getWeight(const Integer& indx,const Integer& ordre){
		Real w = wgauss1; // default is order 1

		switch(ordre){
			case 2: w = wgauss2[indx]; break;
			case 3: w = wgauss3[indx]; break;
			case 4: w = wgauss4[indx]; break;
			case 5: w = wgauss5[indx]; break;
			case 6: w = wgauss6[indx]; break;
			case 7: w = wgauss7[indx]; break;
			case 8: w = wgauss8[indx]; break;
			case 9: w = wgauss9[indx]; break;
			default: break;
		}
		return w;
}

/*---------------------------------------------------------------------------*/

Real3 LineRefPosition(const Integer3& indices,const Integer3& ordre){
	return {getRefPosition(indices[0],ordre[0]),0.,0.};
}

Real LineWeight(const Integer3& indices,const Integer3& ordre){
	return getWeight(indices[0],ordre[0]);
}

/*---------------------------------------------------------------------------*/

Real3 TriRefPosition(const Integer3& indices,const Integer3& ordre){
	Integer o = ordre[0];
	Integer i = indices[0];
	return {xg1[o][i],xg2[o][i],xg3[o][i]};
}

Real TriWeight(const Integer3& indices,const Integer3& ordre){
	return wg[ordre[0]][indices[0]];
}

/*---------------------------------------------------------------------------*/

Real3 QuadRefPosition(const Integer3& indices,const Integer3& ordre){
	Real3 pos;
	for (Integer i = 0; i < 2; i++) pos[i] = getRefPosition(indices[i],ordre[i]);
	return pos;
}

Real QuadWeight(const Integer3& indices,const Integer3& ordre){
	Real w = 1.;
	for (Integer i = 0; i < 2; i++) w *= getWeight(indices[i],ordre[i]);
	return w;
}

/*---------------------------------------------------------------------------*/

Real3 HexaRefPosition(const Integer3& indices,const Integer3& ordre){
	Real3 pos;
	for (Integer i = 0; i < 3; i++) pos[i] = getRefPosition(indices[i],ordre[i]);
	return pos;
}

Real HexaWeight(const Integer3& indices,const Integer3& ordre){
	Real w = 1.;
	for (Integer i = 0; i < 3; i++) w *= getWeight(indices[i],ordre[i]);
	return w;
}

/*---------------------------------------------------------------------------*/

Real3 TetraRefPosition(const Integer3& indices,const Integer3& ordre){
	Integer i = indices[0];
	return {xit[i],yit[i],zit[i]};
}

Real TetraWeight(const Integer3& indices,const Integer3& ordre){
	return wgtetra;
}

/*---------------------------------------------------------------------------*/
