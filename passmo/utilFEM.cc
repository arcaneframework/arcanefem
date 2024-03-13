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
//CellFEMDispatcher::CellFEMDispatcher(VariableNodeReal3& node_coords): m_node_coords(node_coords)
CellFEMDispatcher::CellFEMDispatcher(){
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
 m_geomfunc[IT_Pentaedron6] = Penta6Volume;
 m_geomfunc[IT_Pyramid5] = Pyramid5Volume;

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
 m_shapefunc[IT_Pentaedron6] = Penta6ShapeFuncVal;
 m_shapefunc[IT_Pyramid5] = Pyramid5ShapeFuncVal;

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
 m_shapefuncderiv[IT_Pentaedron6] = Penta6ShapeFuncDeriv;
 m_shapefuncderiv[IT_Pyramid5] = Pyramid5ShapeFuncDeriv;

 // Quadratic elements
 m_shapefuncderiv[IT_Line3] = Line3ShapeFuncDeriv;
 m_shapefuncderiv[IT_Triangle6] = Tri6ShapeFuncDeriv;
 m_shapefuncderiv[IT_Quad8] = Quad8ShapeFuncDeriv;
 m_shapefuncderiv[IT_Tetraedron10] = Tetra10ShapeFuncDeriv;
 m_shapefuncderiv[IT_Hexaedron20] = Hexa20ShapeFuncDeriv;

 // Gives functions to compute orientation vector of the finite-element in global axes (x, y or z)
 // Linear elements
/*
 m_orientfunc[IT_Line2] = Line2Orientation;
 m_orientfunc[IT_Triangle3] = Tri3Orientation;
 m_orientfunc[IT_Quad4] = Quad4Orientation;
 m_orientfunc[IT_Tetraedron4] = Tetra4Orientation;
 m_orientfunc[IT_Hexaedron8] = Hexa8Orientation;
*/

 // Quadratic elements
/*
 m_orientfunc[IT_Line3] = Line3Orientation;
 m_orientfunc[IT_Triangle6] = Tri6Orientation;
 m_orientfunc[IT_Quad8] = Quad8Orientation;
 m_orientfunc[IT_Tetraedron10] = Tetra10Orientation;
 m_orientfunc[IT_Hexaedron20] = Hexa20Orientation;
*/
}

/////////////////////////////////////////////////////////////////////////////
// class CellFEMDispatcher: implementation methods

// ! Computes the Jacobian = geometric size (length, surface or volume) of any finite-element
// aligned with the coordinate axes (otherwise, jacobian to be computed by the
// general Jacobian matrix)
void CellFEMDispatcher::
set_node_coords(VariableNodeReal3& node_coords){  m_node_coords = &node_coords; }

RealUniqueArray CellFEMDispatcher::
getGaussData(const ItemWithNodes& item, const Integer3& nint, Int32& ngauss){

 Int32 ndim = getGeomDimension(item);
 const Int32& nnod = item.nbNode();
 ngauss = nint.m_i;
 auto cell_type = item.type();
 auto nint1 {nint.m_i};
 auto nint2 {nint.m_j};
 auto nint3 {nint.m_k};

 if (ndim >= 2) {
   if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
     nint2 = 1;
     nint1 = nptg[nint.m_i-1];
   }
   else if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10){
     nint2 = 1;
     nint3 = 1;
     nint1 = 4;
   }

   ngauss = nint1*nint2;

   if (ndim == 3)
     ngauss *= nint3;
 }

 // Vector of double containing:
 // ngauss points * [weight, gauss ref coord [Real3], nnod * (shapefunc values, 3*shapefunc deriv
 // in ref. coord system)]
 Int32 nsize = ngauss * 4 * (1 + nnod);
 RealUniqueArray vec(nsize);

 GaussPointDispatcher gausspt;
 gausspt.init_order(nint);

 Int32 index{ 0 };

 if (ndim == 2) {
   for (Int32 i1 = 0; i1 < nint1; ++i1) {
     for (Int32 i2 = 0; i2 < nint2; ++i2) {

       Integer3 indices{ i1, i2, -1 };
       auto wt = gausspt.getWeight(item, indices);
       auto pos = gausspt.getRefPosition(item, indices);
       vec[index++] = wt;
       vec[index++] = pos.x;
       vec[index++] = pos.y;
       vec[index++] = 0.;

       for (Int32 inod = 0; inod < item.nbNode(); ++inod) {
         auto Phi_i = getShapeFuncVal(item.type(), inod, pos);
         vec[index++] = Phi_i;
         auto dPhi = getShapeFuncDeriv(item.type(), inod, pos);
         vec[index++] = dPhi.x;
         vec[index++] = dPhi.y;
         vec[index++] = 0.;
       }
     }
   }
 } else if (ndim == 3) {

   for (Int32 i1 = 0; i1 < nint.m_i; ++i1) {
     for (Int32 i2 = 0; i2 < nint.m_j; ++i2) {
       for (Int32 i3 = 0; i3 < nint.m_k; ++i3) {
         Integer3 indices{ i1, i2, i3 };
         auto wt = gausspt.getWeight(item, indices);
         auto pos = gausspt.getRefPosition(item, indices);
         vec[index++] = wt;
         vec[index++] = pos.x;
         vec[index++] = pos.y;
         vec[index++] = pos.z;

         for (Int32 inod = 0; inod < item.nbNode(); ++inod) {
           auto num = item.node(inod).uniqueId();
           auto Phi_i = getShapeFuncVal(item.type(), inod, pos);
           vec[index++] = Phi_i;
           auto dPhi = getShapeFuncDeriv(item.type(), inod, pos);
           vec[index++] = dPhi.x;
           vec[index++] = dPhi.y;
           vec[index++] = dPhi.z;
         }
       }
     }
   }
 } else {

   for (Int32 i1 = 0; i1 < nint.m_i; ++i1) {

     Integer3 indices{ i1, -1, -1 };
     auto wt = getWeight(i1, nint.m_i);
     Real pos{ getRefPosition(i1, nint.m_i) };
     vec[index++] = wt;
     vec[index++] = pos;
     vec[index++] = 0.;
     vec[index++] = 0.;

     for (Int32 inod = 0; inod < item.nbNode(); ++inod) {
       auto Phi_i = getShapeFuncVal(item.type(), inod, { pos, 0., 0. });
       vec[index++] = Phi_i;
       auto dPhi = getShapeFuncDeriv(item.type(), inod, { pos, 0., 0. });
       vec[index++] = dPhi.x;
       vec[index++] = 0.;
       vec[index++] = 0.;
     }
   }
 }
 return vec;
}

// ! Compute Length, Area or Volume of an element
/*
Real CellFEMDispatcher::
getMeasure(const ItemWithNodes& item)
{
 Int32 item_type = item.type();
 auto f = m_geomfunc[item_type];
 if (f!=nullptr)
   return f(item,*m_node_coords);
 return (0.0);
}
*/

// ! Computes the barycenter of a Cell
/*
Real3 CellFEMDispatcher::
getBarycenter(const ItemWithNodes& item)
{

 Real3 sum = Real3::zero();
 Integer nnod = item.nbNode();

 for (Integer i = 0; i < nnod; i++) {
   sum += (*m_node_coords)[item.node(i)];
 }

 if (nnod) sum /= nnod;
 return sum;
}
*/

// ! Computes the value of the nodal shape functions of a given FE element at a given Gauss point (allowing to the element type)
Real CellFEMDispatcher::
getShapeFuncVal(const Int16& item_type,const Int32& inod,const Real3& coord)
{
 auto f = m_shapefunc[item_type];
 if (f!=nullptr)
   return f(inod,coord);
 return 0.;
}

// ! Gives the values of the derivatives for the function of node inod for a given FE element at a given Gauss point
// ! Derivation is performed along 3 directions x, y,z
Real3 CellFEMDispatcher::
getShapeFuncDeriv(const Int16& item_type,const Int32& inod,const Real3& ref_coord){
 auto f = m_shapefuncderiv[item_type];
 if (f!=nullptr)
   return f(inod,ref_coord);
 return {};
}

// ! Computes the cartesian directions of a finite-element (used for derivation of shape functions)
// ! For instance, assuming the returned Integer3 vector=idir,
// idir = (0,-1,-1) means there will be derivation along x only
/*
Integer3 CellFEMDispatcher::
getOrientation(const ItemWithNodes& cell){
 Int32 item_type = cell.type();
 auto f = m_orientfunc[item_type];
 if (f!=nullptr)
   return f(cell,*m_node_coords);
 return Integer3::zero();
}
*/

/*---------------------------------------------------------------------------*/
/* Functions used for geometric transformations (rotations, projections, ...)*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Edge & Face normal & tangent vectors (normalized, direct oriented)        */
/* 2D: assuming the edge lies in x-y plane (z coord = 0.)                    */
/*---------------------------------------------------------------------------*/
void DirVectors(const Face& face,const VariableNodeReal3& n, const Int32& ndim, Real3& e1, Real3& e2, Real3& e3){

 Real3 n0 = n[face.node(0)];
 Real3 n1 = n[face.node(1)];

 if (!face.isSubDomainBoundaryOutside())
   std::swap(n0, n1);

 // 1st in-plane vector/along edge
 e1 = n1 - n0;

 if (ndim == 3) {

   const Real3& n2 = n[face.node(2)];

   // out Normal to the face plane
   e3 = math::cross(e1, n2 - n0);

   // 2nd in-plane vector
   e2 = math::cross(e3,e1);

   e3.normalize();
 }
 else {

   Cell cell {face.boundaryCell()};
   Node nod;
   for (Node node : cell.nodes()) {
     if (node != face.node(0) && node != face.node(1)) {
       nod = node;
       break;
     }
   }

   // Out Normal to the edge
   e2 = { -e1.y, e1.x, 0. };

   const Real3& n2 = n[nod];
   auto sgn = math::dot(e2,n2 - n0);
   if (sgn > 0.) e2 *= -1.;
 }
 e1.normalize();
 e2.normalize();
}

/*---------------------------------------------------------------------------*/
/* Functions useful for class CellFEMDispatcher                              */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
// Line2: linear edge finite-element
//
// 0           1
// o-----------o---> x
//-1           1
// direct : 0->1 (local numbering)
/*---------------------------------------------------------------------------*/

Real Line2Length(const ItemWithNodes& item,const VariableNodeReal3& n){
 const Real3& n0 = n[item.node(0)];
 const Real3& n1 = n[item.node(1)];
 return (n1-n0).normL2();
}

Real Line2ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
 assert(inod >= 0 && inod < 2);
#endif

 Real r = ref_coord[0];
 if (inod == 1) return (0.5*(1 + r));
 return (0.5*(1 - r));
}

Real3 Line2ShapeFuncDeriv(const Integer& inod,const Real3&){
 if (inod == 1) return { 0.5,0.,0. };
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
 const Real3& n0 = n[item.node(0)];
 const Real3& n1 = n[item.node(1)];
 const Real3& n2 = n[item.node(2)];

 auto v1 = n1 - n0;
 auto v2 = n2 - n0;
 auto v = math::cross(v1,v2);

 auto norm = v.normL2();
 return 0.5*norm;
}

Real Tri3ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
 assert(inod >= 0 && inod < 3);
#endif
 Real r = ref_coord[0];
 Real s = ref_coord[1];
 if (!inod) return (1 - r - s);
 if (inod == 1) return r;
 return s;
}

Real3 Tri3ShapeFuncDeriv(const Integer& inod,const Real3&){
 if (!inod) return {-1., -1.,0.};
 if (inod==1) return {1., 0.,0.};
 return {0., 1., 0.};
}

/*
Integer3 Tri3Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
 const Real3& p0 = n[item.node(0)];
 Real3 u1 = p0 - n[item.node(1)];
 Real3 u2 = p0 - n[item.node(2)];
 Integer3 dir;

 for (Int32 i = 0, j = -1; i < 3; i++)	{
   Real3 e;
   e[i] = 1.;
   if (math::dot(u1,e) != 0. || math::dot(u2,e) != 0.) dir[i] = ++j;
 }
 return dir;
}
*/

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

/*
Integer3 Tri6Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
 return Tri3Orientation(item,n);
}
*/

/*---------------------------------------------------------------------------*/
// Quad4: linear quadrangle finite-element
//
//        ^y
//        1|
//  1o-----|-----o0
//   |     |     |
//   |     |     |
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
 const Real3& n0 = n[item.node(0)];
 const Real3& n1 = n[item.node(1)];
 const Real3& n2 = n[item.node(2)];

 auto v = math::cross(n1 - n0,n2 - n0);
 return v.normL2();
}

Real Quad4ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
 assert(inod >= 0 && inod < 4);
#endif

 auto	r{ ref_coord[0] },s{ ref_coord[1] };
 auto	ri{1.},si{1.};

 switch(inod){
 default: break;// default is first node (index 0)
 case 2:	si = -1;
 case 1:	ri = -1; break;

 case 3:	si = -1; break;
 }
 return ( (1 + ri*r)*(1 + si*s) / 4. );
}

Real3 Quad4ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
 assert(inod >= 0 && inod < 4);
#endif

 auto	r{ ref_coord[0] },s{ ref_coord[1] };
 auto	ri{1.},si{1.}; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

 switch(inod){
 default: break;// default is first node (index 0)
 case 2:	si = -1;
 case 1:	ri = -1; break;

 case 3:	si = -1; break;
 }
 return {0.25 * ri * (1 + si*s), 0.25 * si * (1 + ri*r), 0.};
}

/*
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
*/

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

 auto	r{ ref_coord[0] },s{ ref_coord[1] };
 auto	ri{1.},si{1.};

 switch(inod){
 default: break;// default is first node (index 0)
 case 2:	si = -1;
 case 1:	ri = -1; break;

 case 3:	si = -1; break;

 case 6:	si = -1;
 case 4:	ri = 0; break;

 case 5:	ri = -1;
 case 7:	si = 0; break;
 }

 auto r0 {r*ri}, s0 {s*si};
 Real Phi{0.};
 auto t0{ r0 + s0 - 1. };

 if (inod < 4) // Corner nodes
   Phi = (1 + r0) * (1 + s0) * t0 / 4.;

 else { // Middle nodes
   if (fabs(ri) < REL_PREC)
     Phi = (1 - r * r) * (1 + s0) / 2.;
   else if (fabs(si) < REL_PREC)
     Phi = (1 - s * s) * (1 + r0) / 2.;
 }
 return Phi;
}

Real3 Quad8ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

 auto	r{ ref_coord[0] },s{ ref_coord[1] };
 auto	ri{1.},si{1.};

 switch(inod){
 default: break;// default is first node (index 0)
 case 2:	si = -1;
 case 1:	ri = -1; break;

 case 3:	si = -1; break;

 case 6:	si = -1;
 case 4:	ri = 0; break;

 case 5:	ri = -1;
 case 7:	si = 0; break;
 }

 auto r0 {r*ri}, s0 {s*si};
 Real Phi{0.};
 Real3 dPhi;
 auto t0{ r0 + s0 - 1. };

 if (inod < 4) { // Corner nodes
   dPhi.x = ri * (1 + s0) * (t0 + 1. + r0) / 4.;
   dPhi.y = si * (1 + r0) * (t0 + 1. + s0) / 4.;
 }
 else { // Middle nodes
   if (fabs(ri) < REL_PREC) {
     dPhi.x = -r * (1 + s0);
     dPhi.y = si * (1 - r * r) / 2.;
   }
   else if (fabs(si) < REL_PREC) {
     dPhi.x = -s * (1 + r0);
     dPhi.y = ri * (1 - s * s) / 2.;
   }
 }
 dPhi.z = 0.;
 return dPhi;
}

/*
Integer3 Quad8Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
 return Quad4Orientation(item,n);
}
*/

/*---------------------------------------------------------------------------*/
// Hexa8: linear hexaedron finite-element
//  Normalized coordinates (triplets): x, y, z varying between -1/+1
//     (-1, 1,1)
//         1-------------0 (1,1,1)
//        /|            /|
//       / |           / |
//     /   |          /  |
//    2----|---------3   |   z   y
//  (-1,-1,1)        |   |   | /
//    |    |         |   |   |/--->x
//    |    |         |   |
//    |    |         |   |
//    |    5---------|---4 (1,1,-1)
//    |  /           |  /
//    | /            | /
//    |/             |/
//    6--------------7 (1,-1,-1)
// (-1,-1,-1)
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
 return (v1 + v2 + v3) / 12.;
}

Real Hexa8ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
 assert(inod >= 0 && inod < 8);
#endif
 auto	x{ ref_coord[0] },y{ ref_coord[1] },z{ ref_coord[2] };
 auto	ri{1.},si{1.}, ti{1.}; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

 switch(inod){
 default: break;
 case 3:
 case 2:	ri = -1; break;
 case 0:
 case 1: ri = -1; si = -1;
   break;
 case 4:
 case 5: si = -1;
   break;
 }
 if (inod == 1 || inod == 2 || inod == 5 || inod == 6) ti = -1;

 auto r0 {x*ri}, s0 {y*si}, t0 {z*ti};
 auto Phi = (1 + r0) * (1 + s0) * (1 + t0) / 8.;

 return Phi;

}

Real3 Hexa8ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

 auto	x{ ref_coord[0] },y{ ref_coord[1] },z{ ref_coord[2] };
 auto	ri{1.},si{1.}, ti{1.}; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

 switch(inod){
 default: break;
 case 3:
 case 2:	ri = -1; break;
 case 0:
 case 1: ri = -1; si = -1;
   break;
 case 4:
 case 5: si = -1;
   break;
 }
 if (inod == 1 || inod == 2 || inod == 5 || inod == 6) ti = -1;

 auto r0 {x*ri}, s0 {y*si}, t0 {z*ti};
 Real3 dPhi;
 dPhi.x = ri * (1 + s0) * (1 + t0) / 8.;
 dPhi.y = si * (1 + r0) * (1 + t0) / 8.;
 dPhi.z = ti * (1 + r0) * (1 + s0) / 8.;
 return dPhi;
}

/*
Integer3 Hexa8Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
 return {1,1,1};
}
*/

/*---------------------------------------------------------------------------*/
// Hexa20: quadratic hexaedron finite-element
//  Normalized coordinates (triplets): x, y, z varying between -1/+1
//    (-1, 1, 1)
//         1------8------0 (1, 1,1)
//        /|            /|
//       9 |           11|
//     /   |          /  |
//    2----|-10------3   |   z   y
//  (-1,-1,1)        |   |   | /
//    |   16         |  19   |/--->x
//    |    |         |   |
//   17    |        18   |
//    |    5-----12--|---4 (1,1,-1)
//    |   /          |  /
//    | 13           | 15
//    |/             |/
//    6------14------7 (1,-1,-1)
// (-1,-1,-1)
// direct : 0,1,2,3,...,7 (local numbering)
/*---------------------------------------------------------------------------*/

Real Hexa20Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
 return Hexa8Volume(item, n);
}

Real Hexa20ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
 assert(inod >= 0 && inod < 20);
#endif
 auto	x{ ref_coord[0] },y{ ref_coord[1] },z{ ref_coord[2] };
 auto	ri{1.},si{1.}, ti{1.}; // Normalized coordinates (=+-1) =>node index 0 = (1,1,1)

 switch(inod){
 default: break;

 case 5:	ti = -1.;
 case 1: ri = -1; break;

 case 6: ti = -1.;
 case 2: ri = -1; si = -1; break;

 case 7: ti = -1.;
 case 3: si = -1; break;

 case 4:	ti = -1.; break;

 case 9: ri = -1.;
 case 11: si = 0.; break;

 case 10: si = -1.;
 case 8: ri = 0.; break;

 case 14: si = -1.;
 case 12: ri = 0.; ti = -1.; break;

 case 17: si = -1.;
 case 16: ri = -1.; ti = 0.; break;

 case 18: si = -1.;
 case 19: ti = 0.; break;
 }

 auto r0 {x*ri}, s0 {y*si}, t0 {z*ti};
 Real Phi{0.};
 auto t{ r0 + s0 + t0 - 2. };

 if (inod < 8)  // Corner nodes
   Phi = (1 + r0) * (1 + s0) * (1 + t0) * t / 8.;

 else{  // Middle nodes
   if (fabs(ri) < REL_PREC)
     Phi = (1 - x*x) * (1 + s0) * (1 + t0) / 4.;
   else if (fabs(si) < REL_PREC)
     Phi = (1 - y*y) * (1 + r0) * (1 + t0) / 4.;
   else if (fabs(ti) < REL_PREC)
     Phi = (1 - z*z) * (1 + r0) * (1 + s0) / 4.;
 }
 return Phi;
}

Real3 Hexa20ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

 auto	x{ ref_coord[0] },y{ ref_coord[1] },z{ ref_coord[2] };
 auto	ri{1.},si{1.}, ti{1.}; // Normalized coordinates (=+-1) =>node index 0 = (1,1,1)

 switch(inod){
 default: break;

 case 5:	ti = -1.;
 case 1: ri = -1; break;

 case 6: ti = -1.;
 case 2: ri = -1; si = -1; break;

 case 7: ti = -1.;
 case 3: si = -1; break;

 case 4:	ti = -1.; break;

 case 9: ri = -1.;
 case 11: si = 0.; break;

 case 10: si = -1.;
 case 8: ri = 0.; break;

 case 14: si = -1.;
 case 12: ri = 0.; ti = -1.; break;

 case 17: si = -1.;
 case 16: ri = -1.; ti = 0.; break;

 case 18: si = -1.;
 case 19: ti = 0.; break;
 }

 auto r0 {x*ri}, s0 {y*si}, t0 {z*ti};
 auto t{ r0 + s0 + t0 - 2. };
 Real3 dPhi;

 if (inod < 8) { // Corner nodes
   dPhi = Hexa8ShapeFuncDeriv(inod, ref_coord);
   dPhi.x *= (t + 1. + r0);
   dPhi.y *= (t + 1. + s0);
   dPhi.z *= (t + 1. + t0);
 }
 else { // Middle nodes
   auto x2{ x * x }, y2{ y * y }, z2{ z * z };
   if (fabs(ri) < REL_PREC) {
     dPhi.x = -x * (1 + s0) * (1 + t0) / 2.;
     dPhi.y = si * (1 - x2) * (1 + t0) / 4.;
     dPhi.z = ti * (1 - x2) * (1 + s0) / 4.;
   }
   else if (fabs(si) < REL_PREC) {
     dPhi.x = ri * (1 - y2) * (1 + t0) / 4.;
     dPhi.y = -y * (1 + r0) * (1 + t0) / 2.;
     dPhi.z = ti * (1 - y2) * (1 + r0) / 4.;
   }
   else if (fabs(ti) < REL_PREC) {
     dPhi.x = ri * (1 - z2) * (1 + s0) / 4.;
     dPhi.y = si * (1 - z2) * (1 + r0) / 4.;
     dPhi.z = -z * (1 + r0) * (1 + s0) / 2.;
   }
 }
 return dPhi;
}

/*
Integer3 Hexa20Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
 return {1,1,1};
}
*/

/*---------------------------------------------------------------------------*/
// Tetra4: linear tetrahedral finite-element
//
//    (0,0,1)                     3
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
//       --------X------>         0********************1
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

 auto	ri = ref_coord[0],si = ref_coord[1],ti = ref_coord[2]; // default is first node (index 3)

 switch(inod){
 default: break;
 case 1:	return ri;
 case 2:	return si;
 case 0:	return (1. - ri - si - ti);
 }
 return ti;
}

Real3 Tetra4ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

 if (inod == 3) return {0.,0.,1.};
 if (inod == 1) return {1.,0.,0.};
 if (inod == 2) return {0.,1.,0.};
 return {-1.,-1.,-1.};
}

/*
Integer3 Tetra4Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
 return {1,1,1};
}
*/

/*---------------------------------------------------------------------------*/
// Tetra10: quadratic tetrahedral finite-element
//
//
//    (0,0,1)                     x3
//       .                        *.*
//       .                        * . *
//       .                        *  .  *
//       .                        *   9   *
//       Z   (0,1,0)              *    .     8
//       .    .                   7      x2    *
//       .   .                    *   .    .    *
//       .  Y                     *  6        5   *
//       . .                      * .            .  *
//       ..           (1,0,0)     *.                . *
//       --------X------>        0x****** 4 ***********x1
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
 case 0:	wi = t * (2*t - 1.); break;//=(1. - 2*x - 2*y - 2*z) * t
 case 1:	wi = x * (2*x - 1.); break;//=(1. - 2*t - 2*y - 2*z)*x
 case 2:	wi = y * (2*y - 1.); break;//=(1. - 2*x - 2*t - 2*z)*y
 case 3:	wi = z * (2*z - 1.);break;//=(1. - 2*t - 2*x - 2*y)*z

 // Middle nodes
 case 4:	wi = 4*x*t;break;
 case 5:	wi = 4*x*y;break;
 case 6:	wi = 4*y*t;break;
 case 7:	wi = 4*z*t;break;
 case 8:	wi = 4*z*x;break;
 case 9:	wi = 4*z*y;break;
 }
 return wi;
}

Real3 Tetra10ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
 auto  x{ ref_coord[0] },y{ ref_coord[1] },z{ ref_coord[2] },
 t{ 1. - x - y - z },
 x4{ 4 * x },
 y4{ 4 * y },
 z4{ 4 * z },
 t4{ 4 * t };

 // Corner nodes
 /*
   if (inod == 3) return {0.,0.,1. + 2*t - 2*x - 2*y + 2*z};
   if (inod == 1) return {1. - 2*t - 2*y - 2*z + 2*x,0.,0.};
   if (inod == 2) return {0.,1. - 2*x - 2*t - 2*z + 2*y,0.};
   if (!inod) return {-1. - 2*t + 2*x + 2*y + 2*z,-1. - 2*t+ 2*x + 2*y + 2*z,-1. - 2*t + 2*x + 2*y + 2*z};
*/
 if (!inod) return {1. - t4, 1. - t4, 1. - t4};
 if (inod == 1) return {x4 - 1., 0.,0.};
 if (inod == 2) return {0., y4 - 1., 0.};
 if (inod == 3) return {0., 0., z4 - 1.};

 // Middle nodes
 if (inod == 4) return {t4 - x4, -x4, -x4};
 if (inod == 5) return {y4, x4, 0.};
 if (inod == 6) return {-y4, t4 - y4, -y4};
 if (inod == 8) return {z4, 0., x4};
 if (inod == 9) return {0., z4, y4};
 return {-z4, -z4, t4 - z4};//inod == 7
}

/*
Integer3 Tetra10Orientation(const ItemWithNodes& item,const VariableNodeReal3& n){
 return {1,1,1};
}
*/

/*---------------------------------------------------------------------------*/
// Penta6: linear pentaedron (wedge or triangular prism) finite-element
//  Normalized coordinates (triplets): x, y, z varying between -1/+1
//
//                     5 (0,1,1)
//                   . |  .
//                  .  |     .
//                 .   Z        .
//                .    |           .
//               .     |             .
//       (0,0,1) 3 ------------------ 4 (1,0,1)
//               |     |              |
//               |     |              |
//               |     |              |
//               |     |              |
//               |     2 (0,1,-1)     |
//               |   .    .           |
//               |  Y        .        |
//               | .            .     |
//               |.                .  |
//      (0,0,-1) 0 -------- X ------- 1 (1,0,-1)
//
// direct : 0,1,2,3,...,6 (local numbering)
/*---------------------------------------------------------------------------*/

Real Penta6Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
 const Real3& n0 = n[item.node(0)];
 const Real3& n1 = n[item.node(1)];
 const Real3& n2 = n[item.node(2)];
 const Real3& n3 = n[item.node(3)];
 const Real3& n4 = n[item.node(4)];
 const Real3& n5 = n[item.node(5)];

 auto v = math::cross(n1 - n0,n2 - n0);
 auto base = 0.5*v.normL2();
 auto h1 = (n3 - n0).normL2();
 auto h2 = (n4 - n1).normL2();
 auto h3 = (n5 - n2).normL2();

 return base * (h1 + h2 + h3)/3.0;
}

Real Penta6ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
 assert(inod >= 0 && inod < 6);
#endif
 auto	r{ ref_coord[0] },s{ ref_coord[1] },t{ ref_coord[2] };
 auto	r0{1.},s0{1.}, ti{-1.};
 auto rs {1.- r - s};

 if (inod >= 3) ti = 1.;
 auto t0{1 + ti*t};

 switch(inod){
   default: break;// Node 0
   case 4:
   case 1:	r0 = r; rs = 1.; break;
   case 5:
   case 2:	s0 = s; rs = 1.; break;
 }

 return 0.5 * r0 * s0 * rs * t0;
}

Real3 Penta6ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

#ifdef _DEBUG
 assert(inod >= 0 && inod < 6);
#endif
 auto	r{ ref_coord[0] },s{ ref_coord[1] },t{ ref_coord[2] };
 auto	ri{1.},si{1.};
 auto	r0{1.},s0{1.}, ti{-1.};
 auto rs {1.- r - s};

 if (inod >= 3) ti = 1.;
 auto t0{1 + ti*t};

 switch(inod){
 default: break;
 case 3:
 case 0:	ri = -1.; si = -1.; break;
 case 4:
 case 1:	r0 = r; si = 0.; rs = 1.; break;
 case 5:
 case 2:	s0 = s; rs = 1.; break;
 }

 Real3 dPhi;
 dPhi.x =  0.5 * ri * t0;
 dPhi.y =  0.5 * si * t0;
 dPhi.z =  0.5 * ti * rs * r0 * s0;
 return dPhi;
}

/*---------------------------------------------------------------------------*/
// Pyramid5: linear pyramid finite-element
//  Normalized coordinates (triplets): x, y, z varying between -1/+1
//
//                          ^
//                          |
//                          Z
//                          |
//                          4 (0,0,1)
//                          *
//                         ** *
//                        *  *  *
//                       *    *   *
//                      *      *     *
//                     *        *       *
//          (-1,1,-1) 3 ---------*-------- 2 (1,1,-1)
//                   .            *     .
//                  Y              *    .
//                 .                *  .
//                .                  *.
//   (-1,-1,-1)  0 -------- X ------- 1 (1,-1,-1)
//
// direct : 0,1,2,3,...,6 (local numbering)
/*---------------------------------------------------------------------------*/

Real Pyramid5Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
 const Real3& n0 = n[item.node(0)];
 const Real3& n1 = n[item.node(1)];
 const Real3& n2 = n[item.node(2)];
 const Real3& n4 = n[item.node(4)];

 auto v = math::cross(n1 - n0,n2 - n0);
 auto base = v.normL2();
 v.normalize();
 auto h = fabs(math::dot(v, n4 - n0));
 return base * h / 3.0;
}

Real Pyramid5ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
 assert(inod >= 0 && inod < 5);
#endif
 auto	r{ ref_coord[0] },s{ ref_coord[1] },t{ ref_coord[2] } ;
 auto	ri{-1.},si{-1.};

 switch(inod){
 default: break;// default is first node (index 0)
 case 2:	si = 1;
 case 1:	ri = 1; break;
 case 3:	si = 1; break;
 }

 if (inod == 4)
   return (1. + t) / 2.;

 return (1. + ri*r) * (1. + si*s) * (1. - t) / 8.;
}

Real3 Pyramid5ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

#ifdef _DEBUG
 assert(inod >= 0 && inod < 5);
#endif
 auto	r{ ref_coord[0] },s{ ref_coord[1] },t{ ref_coord[2] } ;
 auto	ri{-1.},si{-1.};

 switch(inod){
 default: break;// default is first node (index 0)
 case 2:	si = 1;
 case 1:	ri = 1; break;
 case 3:	si = 1; break;
 }

 Real3 dPhi;
 if (inod == 4) {
   dPhi.x = dPhi.y = 0.;
   dPhi.z = 0.5;
 }
 else{
   dPhi.x = ri * (1. + si*s) * (1. - t) / 8.;
   dPhi.y = si * (1. + ri*r) * (1. - t) / 8.;
   dPhi.z = - (1. + ri*r) * (1. + si*s) / 8.;
 }

 return dPhi;
}

/*---------------------------------------------------------------------------*/

Int32 getGeomDimension(const ItemWithNodes& item){
 Int32 item_type = item.type();
 Int32 dim = 1; // default geometric dimension is 1D (Line2 and Line3 finite-elements)

 switch(item_type) {

 // 2D elements
 case IT_Triangle3:
 case IT_Quad4:
 case IT_Triangle6:
 case IT_Quad8: dim = 2; break;

   // 3D elements
 case IT_Tetraedron4:
 case IT_Hexaedron8:
 case IT_Tetraedron10:
 case IT_Hexaedron20: dim = 3; break;

 default: break;

 }
 return dim;
}

/*---------------------------------------------------------------------------*/
/////////////////////////////////////////////////////////////////////////////
// class GaussPointDispatcher: construction methods

//GaussPointDispatcher::GaussPointDispatcher(const Integer3& indices,const Integer3& int_order):
// m_indices(indices),m_integ_order(int_order){
GaussPointDispatcher::GaussPointDispatcher(){
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
 m_weightfunc[IT_Pentaedron6] = PentaWeight;

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
 m_refpositionfunc[IT_Pentaedron6] = PentaRefPosition;

 // Quadratic elements
 m_refpositionfunc[IT_Line3] = LineRefPosition;
 m_refpositionfunc[IT_Triangle6] = TriRefPosition;
 m_refpositionfunc[IT_Quad8] = QuadRefPosition;
 m_refpositionfunc[IT_Tetraedron10] = TetraRefPosition;
 m_refpositionfunc[IT_Hexaedron20] = HexaRefPosition;

}

/////////////////////////////////////////////////////////////////////////////
// class GaussPointDispatcher: implementation methods
void GaussPointDispatcher::init_order(const Integer3& int_order) { m_integ_order = int_order; }

Real3 GaussPointDispatcher::getRefPosition(const ItemWithNodes& item, const Integer3& indices){
 Int32 item_type = item.type();
 auto f = m_refpositionfunc[item_type];
 if (f!=nullptr)
   return f(indices,m_integ_order);
 return Real3::zero();
}

Real GaussPointDispatcher::getWeight(const ItemWithNodes& item, const Integer3& indices){
 Int32 item_type = item.type();
 auto f = m_weightfunc[item_type];
 if (f!=nullptr)
   return f(indices,m_integ_order);
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
 Integer o = ordre[0]-1;
 Integer i = indices[0];
 return {xg1[o][i],xg2[o][i],0.};
}

Real TriWeight(const Integer3& indices,const Integer3& ordre){
 return wg[ordre[0]-1][indices[0]];
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

Real3 TetraRefPosition(const Integer3& indices,const Integer3& /*ordre*/){
 Integer i = indices[0];
 return {xit[i],yit[i],zit[i]};
}

Real TetraWeight(const Integer3& indices,const Integer3& /*ordre*/){
 return wgtetra;
}

/*---------------------------------------------------------------------------*/

Real3 PentaRefPosition(const Integer3& indices,const Integer3& ordre){

 // Same as TriRefPosition on reference coordinate plane (r,s)
 // and LineRefPosition along reference coordinate t (vertical)
 auto pos = TriRefPosition(indices,ordre);
 pos.z = getRefPosition(indices[2],ordre[2]);

 return pos;
}

Real PentaWeight(const Integer3& indices,const Integer3& ordre){

 // Same as TriWeight on reference coordinate plane (r,s)
 // and LineWeight with ordre[2] to account for reference coordinate t (vertical)
 Real wgpenta = TriWeight(indices,ordre)*getWeight(indices[2],ordre[2]);
 return wgpenta;
}

/*---------------------------------------------------------------------------*/