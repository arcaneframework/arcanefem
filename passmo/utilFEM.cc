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

// ! Computes the geometric size (length, surface or volume) of any finite-element
Real CellFEMDispatcher::getGeomSize(const ItemWithNodes& item)
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
Real CellFEMDispatcher::getShapeFuncVal(const ItemWithNodes& item,const Integer& inod,const Real3& coord)
{
    Int32 item_type = item.type();
    auto f = m_shapefunc[item_type];
    if (f!=nullptr)
      return f(inod,coord);
    return 0.;
}

// ! Gives the values of the derivatives for all nodal shape functions for a given FE element at a given Gauss point
// ! Size of the returned vector = number of nodes allowing to the element type
// ! Derivation is performed along idir direction:
// ! idir = 0 -> along x
//        = 1 -> along y
//        = 2 -> along z
RealUniqueArray CellFEMDispatcher::getShapeFuncDeriv(const ItemWithNodes& item,const Integer& idir,const Real3& coord)
{
    Int32 item_type = item.type();
    auto f = m_shapefuncderiv[item_type];
    if (f!=nullptr)
      return f(idir,coord);
    return RealUniqueArray();
}

// ! Computes the cartesian directions of a finite-element (used for derivation of shape functions)
// ! For instance, assuming the returned Integer3 vector=idir, idir = (0,-1,-1) means there will derivation along x only
Integer3 CellFEMDispatcher::getOrientation(const ItemWithNodes& item)
{
     Int32 item_type = item.type();
     auto f = m_orientfunc[item_type];
     if (f!=nullptr)
       return f(item,m_node_coords);
     return Integer3::zero();
}

// ! Computes the Jacobian of a finite-element at a Gauss point given by its local coordinates in the element
Real CellFEMDispatcher::computeJacobian(const ItemWithNodes& item,const Real3& coord)
{
	  Integer	n = item.nbNode();
//	  Integer	ndim = getGeomDimension(item);
	  Integer3	dir = getOrientation(item);

	  // Jacobian matrix computed at the integration point
	  Real3x3	jac;

	  for (int i = 0, ii = 0; i < 3; i++) {

	      int idir = dir[i];
		  if (idir != -1)
		  {
			  // vector of local derivatives at this integration point, for all nodes of the cell
			  RealUniqueArray	dwi = getShapeFuncDeriv(item,idir,coord);

			  for (int j = 0, jj = 0; j < 3; j++) {

				  Real temp = 0.;

                  if (dir[j] != -1) {
                      for (int k = 0; k < n; k++) {

                          Real3 coord_nodk = m_node_coords[item.node(k)];
                          temp += dwi[k] * coord_nodk[j];
                      }
                      jac[ii][jj++] = temp;
                  }
			  }
			  ii++;
		  }
	  }

	  // Jacobian = determinant of matrix jac
      Real jacobian = math::matrixDeterminant(jac);

	  if (fabs(jacobian) < REL_PREC)
	  {
		  // Error
	  }
	  return jacobian;
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

Real Line2Length(const ItemWithNodes& item,const VariableNodeReal3& n)
{
  Real3 n0 = n[item.node(0)];
  Real3 n1 = n[item.node(1)];

  return (n0-n1).normL2();
}

Real Line2ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
	assert(inod >= 0 && inod < 2);
#endif

	Real r = coord[0];
	if (!inod) return (0.5*(1 + r));
	return (0.5*(1 - r));
}

RealUniqueArray Line2ShapeFuncDeriv(const Integer& idir,const Real3& coord)
{
	RealUniqueArray vec(2);

	// idir = local direction for derivation
	// derivation vector at all nodes of the finite-element
	if (idir >= 0)
	{
		vec[0] = 0.5;
		vec[1] = -0.5;
	}
	return vec;
}

Integer3 Line2Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
	Real3 n0 = n[item.node(0)];
	Real3 n1 = n[item.node(1)];
	Integer3 dir;

	Real3 pt = n0-n1;

	Integer j = -1;
	for (Integer i = 0; i < 3; i++)
	{
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

Real Line3Length(const ItemWithNodes& item,const VariableNodeReal3& n)
{
    return Line2Length(item,n);
}

Real Line3ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
        assert(inod >= 0 && inod < 3);
#endif

    Real ri = coord[0];
    if (inod == 1) ri *= -1;

    if (inod < 2) return 0.5*ri*(1 + ri); // nodes 0 or 1
    return (1 - ri*ri); // middle node
}

RealUniqueArray Line3ShapeFuncDeriv(const Integer& idir,const Real3& coord)
{
    RealUniqueArray vec(3);

    // idir = local direction for derivation
    // derivation vector at all nodes of the finite-element
    if (idir >= 0)
    {
        vec[0] = 0.5 + coord[0];
        vec[1] = -0.5 + coord[0];
        vec[2] = -2.*coord[0];
    }
    return vec;
}

Integer3 Line3Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
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

Real Tri3Surface(const ItemWithNodes& item,const VariableNodeReal3& n)
{
  Real3 n0 = n[item.node(0)];
  Real3 n1 = n[item.node(1)];
  Real3 n2 = n[item.node(2)];

  Real x1 = n1.x - n0.x;
  Real y1 = n1.y - n0.y;
  Real x2 = n2.x - n1.x;
  Real y2 = n2.y - n1.y;

  return x1 * y2 - y1 * x2;
}

Real Tri3ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
	assert(inod >= 0 && inod < 3);
#endif

	return coord[inod];
}

RealUniqueArray Tri3ShapeFuncDeriv(const Integer& idir,const Real3&)
{
#ifdef _DEBUG
	assert (idir >= 0 && idir < 2);
#endif

	// idir = local direction for derivation
	// derivation vector at all nodes of the finite-element
	RealUniqueArray vec(3);

	vec[2] = -1.;
	vec[idir] = 1.;
	return vec;
}

Integer3 Tri3Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
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

Real Tri6Surface(const ItemWithNodes& item,const VariableNodeReal3& n)
{
    return Tri3Surface(item,n);
}

Real Tri6ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
    assert(inod >= 0 && inod < 6);
#endif

    if (inod < 3) return coord[inod];

    double	wi = 0.,ri = coord[0],si = coord[1],ti = coord[2];
    switch(inod)
    {
        default: break;

        case 3:	wi = 4.*ri*si; break;

        case 4:	wi = 4.*si*ti; break;

        case 5:	wi = 4.*ri*ti; break;
    }
    return wi;
}

RealUniqueArray Tri6ShapeFuncDeriv(const Integer& idir,const Real3& coord)
{
#ifdef _DEBUG
    assert (idir >= 0 && idir < 2);
#endif

    // idir = local direction for derivation
    // derivation vector at all nodes of the finite-element
    RealUniqueArray vec(6);

    vec[idir] = 1.;
    vec[2] = -1.;

    double	ri = coord[0],si = coord[1],ti = coord[2];
    if (!idir)
    {
        vec[3] = 4.*si;
        vec[4] = -4.*si;
        vec[5] = 4.*(ti - ri);
    }
    else
    {
        vec[3] = 4.*ri;
        vec[4] = 4.*(ti - si);
        vec[5] = -4.*ri;
    }
    return vec;
}

Integer3 Tri6Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
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

Real Quad4Surface(const ItemWithNodes& item,const VariableNodeReal3& n)
{
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


Real Quad4ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
	assert(inod >= 0 && inod < 4);
#endif

	Real	ri = coord[0],si = coord[1]; // default is first node (index 0)

	switch(inod)
	{
		default: break;
		case 1:	ri *= -1; break;
		case 2:	ri *= -1; si *= -1; break;
		case 3:	si *= -1; break;
	}

	return ( 0.25*(1 + ri)*(1 + si) );
}

RealUniqueArray Quad4ShapeFuncDeriv(const Integer& idir,const Real3& coord)
{
#ifdef _DEBUG
	assert (idir >= 0 && idir < 2);
#endif

	// idir = local direction for derivation
	// derivation vector at all nodes of the finite-element
	RealUniqueArray vec(4);

	if (!idir)
	{
		Real val = 0.25*(1 + coord[1]);
		vec[0] = val;
		vec[1] = -val;
		val = -0.25*(1 - coord[1]);
		vec[2] = val;
		vec[3] = -val;
	}
	else
	{
		Real val = 0.25*(1 + coord[0]);
		vec[0] = val;
		vec[3] = -val;
		val = 0.25*(1 - coord[0]);
		vec[1] = val;
		vec[2] = -val;
	}
	return vec;
}

Integer3 Quad4Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
	Real3 p0 = n[item.node(0)];
	Real3 u1 = p0 - n[item.node(1)];
	Real3 u2 = p0 - n[item.node(2)];
	Integer3 dir;

	for (Integer i = 0,j = -1; i < 3; i++)
	{
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

Real Quad8Surface(const ItemWithNodes& item,const VariableNodeReal3& n)
{
    return Quad4Surface(item,n);
}


Real Quad8ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
    assert(inod >= 0 && inod < 8);
#endif

    Real	wi = 0.,
            r = coord[0],s = coord[1],
            ri = r,si = s; // default is first node (index 0)

    switch(inod)
    {
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

RealUniqueArray Quad8ShapeFuncDeriv(const Integer& idir,const Real3& coord)
{
#ifdef _DEBUG
    assert (idir >= 0 && idir < 2);
#endif

    // idir = local direction for derivation
    // derivation vector at all nodes of the finite-element
    RealUniqueArray vec(8);
    RealUniqueArray r(8),s(8);

    r[0] = 1; s[0] = 1;
    r[1] = -1; s[1] = 1;
    r[2] = -1; s[2] = -1;
    r[3] = 1; s[3] = -1;
    r[4] = 0; s[4] = 1;
    r[5] = -1; s[5] = 0;
    r[6] = 0; s[6] = -1;
    r[7] = 1; s[7] = 0;

    if (!idir) // r direction
    {
        // derivatives on 4 corner nodes
        for (int j = 0; j < 4; j++)
        {
            Real    wj = Quad4ShapeFuncVal(j,coord),
                    rj = r[j]*coord[0],
                    sj = s[j]*coord[1];

            vec[j] = r[j]*(0.25*(1 + sj)*(rj + sj - 1) + wj);
        }

        Real s2 = coord[1]*coord[1];
        // derivatives on 4 middle nodes
        for (int j = 4; j < 8; j++)
        {
            Real    rj = r[j]*coord[0],
                    sj = s[j]*coord[1];

            if (!rj) vec[j] = -coord[0]*(1 + sj);
            else if (!sj) vec[j] = 0.5*rj*(1 - s2);
        }
    }
    else // s direction
    {
        // derivatives on 4 corner nodes
        for (int j = 0; j < 4; j++)
        {
            Real    wj = Quad4ShapeFuncVal(j,coord),
                    rj = r[j]*coord[0],
                    sj = s[j]*coord[1];

            vec[j] = s[j]*(0.25*(1 + rj)*(rj + sj - 1) + wj);
        }

        Real r2 = coord[0]*coord[0];
        // derivatives on 4 middle nodes
        for (int j = 4; j < 8; j++)
        {
            Real rj = r[j]*coord[0], sj = s[j]*coord[1];
            if (!rj) vec[j] = 0.5*sj*(1 - r2);
            else if (!sj) vec[j] = -coord[1]*(1 + rj);
        }
    }
    return vec;
}

Integer3 Quad8Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
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

Real Hexa8Volume(const ItemWithNodes& item,const VariableNodeReal3& n)
{
  Real3 n0 = n[item.node(0)];
  Real3 n1 = n[item.node(1)];
  Real3 n2 = n[item.node(2)];
  Real3 n3 = n[item.node(3)];
  Real3 n4 = n[item.node(4)];
  Real3 n5 = n[item.node(5)];
  Real3 n6 = n[item.node(6)];
  Real3 n7 = n[item.node(7)];

  Real v1 = math::matDet((n6 - n1) + (n7 - n0), n6 - n3, n2 - n0);
  Real v2 = math::matDet(n7 - n0, (n6 - n3) + (n5 - n0), n6 - n4);
  Real v3 = math::matDet(n6 - n1, n5 - n0, (n6 - n4) + (n2 - n0));

  Real res = (v1 + v2 + v3) / 12.0;

  return res;
}

Real Hexa8ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
	assert(inod >= 0 && inod < 8);
#endif

	Real	ri = coord[0],si = coord[1],ti = coord[2]; // default is first node (index 0)

	switch(inod)
	{
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

RealUniqueArray Hexa8ShapeFuncDeriv(const Integer& idir,const Real3& coord)
{
#ifdef _DEBUG
	assert(idir >= 0 && idir < 3);
#endif

	// idir = local direction for derivation
	// derivation vector at all nodes of the finite-element
	RealUniqueArray	vec(8);
	Real	rp = 1.+ coord[0],sp = 1. + coord[1],tp = 1. + coord[2],
			rm = 1.- coord[0],sm = 1. - coord[1],tm = 1. - coord[2];

	if (!idir)
	{
		auto val = 0.125*sp;
		vec[0] = val*tp;
		vec[1] = -val*tp;
		vec[4] = val*tm;
		vec[5] = -val*tm;

		val = -0.125*sm;
		vec[2] = val*tp;
		vec[3] = -val*tp;
		vec[6] = val*tm;
		vec[7] = -val*tm;
	}
	else if (idir == 1)
	{
		auto val = 0.125*rp;
		vec[0] = val*tp;
		vec[3] = -val*tp;
		vec[4] = val*tm;
		vec[7] = -val*tm;

		val = 0.125*rm;
		vec[1] = val*tp;
		vec[2] = -val*tp;
		vec[5] = val*tm;
		vec[6] = -val*tm;
	}
	else if (idir == 2)
	{
		auto val = 0.125*rp;
		vec[0] = val*sp;
		vec[4] = -val*sp;
		vec[3] = val*sm;
		vec[7] = -val*sm;

		val = 0.125*rm;
		vec[1] = val*sp;
		vec[5] = -val*sp;
		vec[2] = val*sm;
		vec[6] = -val*sm;
	}
	return vec;
}

Integer3 Hexa8Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
	return Integer3(1,1,1);
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

Real Hexa20Volume(const ItemWithNodes& item,const VariableNodeReal3& n)
{
    return Hexa8Volume(item, n);
}

Real Hexa20ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
    assert(inod >= 0 && inod < 20);
#endif

    Real	wi = 0.,ri = coord[0],si = coord[1],ti = coord[2],
            rm = 1. - ri,rp = 1. + ri,
            sm = 1. - si,sp = 1. + si,
            tm = 1. - ti,tp = 1. + ti,
            r2 = 1. - ri*ri,s2 = 1. - si*si,t2 = 1. - ti*ti;

    if (inod < 8) return Hexa8ShapeFuncVal(inod,coord);

    switch(inod)
    {
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

RealUniqueArray Hexa20ShapeFuncDeriv(const Integer& idir,const Real3& coord)
{
#ifdef _DEBUG
    assert(idir >= 0 && idir < 3);
#endif

    // idir = local direction for derivation
    RealUniqueArray	vec(20),    // derivation vector at all nodes of the finite-element
            vec_h8(Hexa8ShapeFuncDeriv(idir,coord));  // derivation vector at the 8 corner nodes

    Real	ri = coord[0],si = coord[1],ti = coord[2],
            rm = 1. - ri,rp = 1. + ri,
            sm = 1. - si,sp = 1. + si,
            tm = 1. - ti,tp = 1. + ti,
            r2 = 1. - ri*ri,s2 = 1. - si*si,t2 = 1. - ti*ti;

    // corner nodes
    for (int j = 0; j < 8; j++)
    {
        vec[j] = vec_h8[j];
    }

    if (!idir) // r direction
    {
        auto val = 0.25*s2;
        vec[9] = -val*tp;
        vec[11] = val*tp;
        vec[13] = -val*tm;
        vec[15] = val*tm;

        val = 0.25*t2;
        vec[16] = val*sp;
        vec[17] = -val*sp;
        vec[18] = -val*sm;
        vec[19] = val*sm;

        val = -0.5*ri;
        vec[8] = val*sp*tp;
        vec[10] = val*sm*tp;
        vec[12] = val*sp*tm;
        vec[14] = val*sm*tm;
    }
    else if (idir == 1) // s direction
    {
        auto val = -0.5*si;
        vec[9] = val*rm*tp;
        vec[11] = val*rp*tp;
        vec[13] = val*rm*tm;
        vec[15] = val*rp*tm;

        val = 0.25*t2;
        vec[16] = val*rp;
        vec[17] = val*rm;
        vec[18] = -val*rm;
        vec[19] = -val*rp;

        val = 0.25*r2;
        vec[8] = val*tp;
        vec[10] = -val*tp;
        vec[12] = val*tm;
        vec[14] = -val*tm;
    }
    else if (idir == 2) // t direction
    {
        auto val = 0.25*s2;
        vec[9] = val*rm;
        vec[11] = val*rp;
        vec[13] = -val*rm;
        vec[15] = -val*rp;

        val = -0.5*ti;
        vec[16] = val*rp*sp;
        vec[17] = val*rm*sp;
        vec[18] = val*rm*sm;
        vec[19] = val*rp*sm;

        val = 0.25*r2;
        vec[8] = val*sp;
        vec[10] = val*sm;
        vec[12] = -val*sp;
        vec[14] = -val*sm;
    }
    return vec;
}

Integer3 Hexa20Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
    return Integer3(1,1,1);
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

Real Tetra4Volume(const ItemWithNodes& item,const VariableNodeReal3& n)
{
  Real3 n0 = n[item.node(0)];
  Real3 n1 = n[item.node(1)];
  Real3 n2 = n[item.node(2)];
  Real3 n3 = n[item.node(3)];

  return math::matDet(n1 - n0, n2 - n0, n3 - n0)/6.;
}

Real Tetra4ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
	assert(inod >= 0 && inod < 4);
#endif

	Real	ri = coord[0],si = coord[1],ti = coord[2]; // default is first node (index 0)

	switch(inod)
	{
		default: break;

		case 1:	return ri;
		case 2:	return si;
		case 3:	return (1. - ri - si - ti);
	}

	return ti;
}

RealUniqueArray Tetra4ShapeFuncDeriv(const Integer& idir,const Real3& coord)
{
#ifdef _DEBUG
	assert(idir >= 0 && idir < 3);
#endif

	// idir = local direction for derivation
	// derivation vector at all nodes of the finite-element
	RealUniqueArray	vec(4);

	if (!idir) // x
	{
		vec[0] = 0.;
		vec[1] = 1.;
		vec[2] = 0.;
		vec[3] = -1.;
	}
	else if (idir == 1) // y
	{
		vec[0] = 0.;
		vec[1] = 0.;
		vec[2] = 1.;
		vec[3] = -1.;
	}
	else if (idir == 2) // z
	{
		vec[0] = 1.;
		vec[1] = 0.;
		vec[2] = 0.;
		vec[3] = -1.;
	}
	return vec;
}

Integer3 Tetra4Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
	return Integer3(1,1,1);
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

Real Tetra10Volume(const ItemWithNodes& item,const VariableNodeReal3& n)
{
    return Tetra4Volume(item,n);
}

Real Tetra10ShapeFuncVal(const Integer& inod,const Real3& coord)
{
#ifdef _DEBUG
    assert(inod >= 0 && inod < 10);
#endif

    Real	x = coord[0],y = coord[1],z = coord[2],
            t = 1. - x - y - z,
            wi = 0.;

    switch(inod)
    {
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

RealUniqueArray Tetra10ShapeFuncDeriv(const Integer& idir,const Real3& coord)
{
#ifdef _DEBUG
    assert(idir >= 0 && idir < 3);
#endif

    // idir = local direction for derivation
    // derivation vector at all nodes of the finite-element
    RealUniqueArray	vec(10);
    Real	x = coord[0],y = coord[1],z = coord[2],
            t = 1. - x - y - z;

    if (!idir) // x direction
    {
        // Corner nodes
        vec[0] = 0.;
        vec[1] = 1. - 2*t - 2*y - 2*z + 2*x;
        vec[2] = 0.;
        vec[3] = -1. - 2*t + 2*x + 2*y + 2*z;

        // Middle nodes
        vec[4] = 4*y;
        vec[5] = -4*y;
        vec[6] = 4*(t - x);
        vec[7] = 4*z;
        vec[8] = 0.;
        vec[9] = -4*z;
    }
    else if (idir == 1) // y direction
    {
        // Corner nodes
        vec[0] = 0.;
        vec[1] = 0.;
        vec[2] = 1. - 2*x - 2*t - 2*z + 2*y;
        vec[3] = -1. - 2*t+ 2*x + 2*y + 2*z;

        // Middle nodes
        vec[4] = 4*x;
        vec[5] = 4*(t - y);
        vec[6] = -4*x;
        vec[7] = 0.;
        vec[8] = 4*z;
        vec[9] = -4*z;
    }
    else if (idir == 2) // z direction
    {
        // Corner nodes
        vec[0] = 1. - 2*t - 2*x - 2*y + 2*z;
        vec[1] = 0.;
        vec[2] = 0.;
        vec[3] = -1. - 2*t + 2*x + 2*y + 2*z;

        // Middle nodes
        vec[4] = 0.;
        vec[5] = -4*y;
        vec[6] = -4*x;
        vec[7] = 4*x;
        vec[8] = 4*y;
        vec[9] = 4*(t - z);
    }
    return vec;
}

Integer3 Tetra10Orientation(const ItemWithNodes& item,const VariableNodeReal3& n)
{
    return Integer3(1,1,1);
}

/*---------------------------------------------------------------------------*/

Integer getGeomDimension(const ItemWithNodes& item)
{
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
        m_indices(indices),m_integ_order(int_order)
{
    // Setting to null default value
    for(Integer i = 0; i < NB_BASIC_ITEM_TYPE; ++i )
    {
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

Real3 GaussPointDispatcher::getRefPosition(const ItemWithNodes& item)
{
	  Int32 item_type = item.type();
	  auto f = m_refpositionfunc[item_type];
	  if (f!=nullptr)
		  return f(m_indices,m_integ_order);
	  return Real3::zero();
}

Real GaussPointDispatcher::getWeight(const ItemWithNodes& item)
{
	  Int32 item_type = item.type();
	  auto f = m_weightfunc[item_type];
	  if (f!=nullptr)
		  return f(m_indices,m_integ_order);
	  return 0.;
}

/*---------------------------------------------------------------------------*/
/////////////////////////////////////////////////////////////////////////////
// Functions useful for class GaussPointDispatcher

Real getRefPosition(const Integer& indx,const Integer& ordre)
{
		Real x = xgauss1; // default is order 1

		switch(ordre)
		{
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

Real getWeight(const Integer& indx,const Integer& ordre)
{
		Real w = wgauss1; // default is order 1

		switch(ordre)
		{
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

Real3 LineRefPosition(const Integer3& indices,const Integer3& ordre)
{
	return Real3(getRefPosition(indices[0],ordre[0]),0.,0.);
}

Real LineWeight(const Integer3& indices,const Integer3& ordre)
{
	return getWeight(indices[0],ordre[0]);
}

/*---------------------------------------------------------------------------*/

Real3 TriRefPosition(const Integer3& indices,const Integer3& ordre)
{
	Integer o = ordre[0];
	Integer i = indices[0];
	return Real3(xg1[o][i],xg2[o][i],xg3[o][i]);
}

Real TriWeight(const Integer3& indices,const Integer3& ordre)
{
	return wg[ordre[0]][indices[0]];
}

/*---------------------------------------------------------------------------*/

Real3 QuadRefPosition(const Integer3& indices,const Integer3& ordre)
{
	Real3 pos;
	for (Integer i = 0; i < 2; i++) pos[i] = getRefPosition(indices[i],ordre[i]);
	return pos;
}

Real QuadWeight(const Integer3& indices,const Integer3& ordre)
{
	Real w = 1.;
	for (Integer i = 0; i < 2; i++) w *= getWeight(indices[i],ordre[i]);
	return w;
}

/*---------------------------------------------------------------------------*/

Real3 HexaRefPosition(const Integer3& indices,const Integer3& ordre)
{
	Real3 pos;
	for (Integer i = 0; i < 3; i++) pos[i] = getRefPosition(indices[i],ordre[i]);
	return pos;
}

Real HexaWeight(const Integer3& indices,const Integer3& ordre)
{
	Real w = 1.;
	for (Integer i = 0; i < 3; i++) w *= getWeight(indices[i],ordre[i]);
	return w;
}

/*---------------------------------------------------------------------------*/

Real3 TetraRefPosition(const Integer3& indices,const Integer3& ordre)
{
	Integer i = indices[0];
	return Real3(xit[i],yit[i],zit[i]);
}

Real TetraWeight(const Integer3& indices,const Integer3& ordre)
{
	return wgtetra;
}

/*---------------------------------------------------------------------------*/
