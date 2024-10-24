#ifndef ARCANE_FEM_FUNCTIONS_H
#define ARCANE_FEM_FUNCTIONS_H

#include <arcane/core/IStandardFunction.h>
#include <arcane/IndexedItemConnectivityView.h>
#include <arcane/VariableTypes.h>
#include <arcane/IMesh.h>
#include "IArcaneFemBC.h"
#include "GaussQuadrature.h"
#include "DoFLinearSystem.h"

using namespace Arcane;
using namespace Arcane::FemUtils;
Real REL_PREC {1.0e-15};

/*---------------------------------------------------------------------------*/
/**
 * @brief Contains various functions & operations related to FEM calculations.
 *
 * The class provides methods organized into different nested classes for:
 * - MeshOperation: Mesh related operations.
 * - FeOperation2D: Finite element operations at element level.
 * - BoundaryConditions2D: Boundary condition related operations.
 */
/*---------------------------------------------------------------------------*/

class ArcaneFemFunctions
{
 public:

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for various mesh-related operations.
   *
   * This class includes static methods for computing geometric properties
   * of mesh elements, such as the area of triangles, the length of edges,
   * and the normal vectors of edges.
   */
  /*---------------------------------------------------------------------------*/
  class MeshOperation
  {

   public:

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the volume of a tetrahedra defined by four nodes.
     *
     * This method calculates the volume using the scalar triple product formula.
     * We do the following:
     *   1. get the four nodes
     *   2. ge the vector representing the edges of tetrahedron
     *   3. compute volume using scalar triple product
     */
    /*---------------------------------------------------------------------------*/

    /* ef: needing to compute with different entity inputs (face or cell) */
    /* static inline Real computeVolumeTetra4(Cell cell, const VariableNodeReal3& node_coord)*/
    static inline Real computeVolumeTetra4(const ItemWithNodes& item, const VariableNodeReal3& node_coord)
    {
/*
      Real3 vertex0 = node_coord[cell.node(0)];
      Real3 vertex1 = node_coord[cell.node(1)];
      Real3 vertex2 = node_coord[cell.node(2)];
      Real3 vertex3 = node_coord[cell.node(3)];
*/
      const Real3& vertex0 = node_coord[item.node(0)];
      const Real3& vertex1 = node_coord[item.node(1)];
      const Real3& vertex2 = node_coord[item.node(2)];
      const Real3& vertex3 = node_coord[item.node(3)];

      Real3 v0 = vertex1 - vertex0;
      Real3 v1 = vertex2 - vertex0;
      Real3 v2 = vertex3 - vertex0;

      return std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2))) / 6.0;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the area of a triangle defined by three nodes.
     *
     * This method calculates the area using the determinant formula for a triangle.
     * The area is computed as half the value of the determinant of the matrix
     * formed by the coordinates of the triangle's vertices.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real computeAreaTria3(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.node(0)];
      Real3 vertex1 = node_coord[cell.node(1)];
      Real3 vertex2 = node_coord[cell.node(2)];

      return 0.5 * ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));
    }

    // ef: needing to compute with different entity inputs (face or cell)
    static inline Real Tri3Surface(const ItemWithNodes& item,const VariableNodeReal3& node_coord){
      const Real3& n0 = node_coord[item.node(0)];
      const Real3& n1 = node_coord[item.node(1)];
      const Real3& n2 = node_coord[item.node(2)];

      auto v = math::cross(n1 - n0,n2 - n0);
      return  0.5*v.normL2();
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the area of a quadrilateral defined by four nodes.
     *
     * This method calculates the area of a quadrilateral by breaking it down 
     * into two triangles and using the determinant formula. The area is computed 
     * as half the value of the determinant of the matrix formed by the coordinates 
     * of the quadrilateral's vertices.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real computeAreaQuad4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];
      Real3 vertex3 = node_coord[cell.nodeId(3)];

      return 0.5 * ((vertex1.x * vertex2.y + vertex2.x * vertex3.y + vertex3.x * vertex0.y + vertex0.x * vertex1.y) - (vertex2.x * vertex1.y + vertex3.x * vertex2.y + vertex0.x * vertex3.y + vertex1.x * vertex0.y));
    }

    // ef: needing to compute with different entity inputs (face or cell)
    static inline Real Quad4Surface(const ItemWithNodes& item,const VariableNodeReal3& n){
      const Real3& n0 = n[item.node(0)];
      const Real3& n1 = n[item.node(1)];
      const Real3& n2 = n[item.node(2)];

      auto v = math::cross(n1 - n0,n2 - n0);
      return v.normL2();
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the volume of a hexaedron defined by eight nodes.
     *
     */
    /*---------------------------------------------------------------------------*/
    static inline Real Hexa8Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
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

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the volume of a pentaedron (wedge or triangular prism)
     * defined by six nodes.
     *
     */
    /*---------------------------------------------------------------------------*/

    static inline Real Penta6Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
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

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the volume of a pyramid defined by five nodes.
     *
     */
    /*---------------------------------------------------------------------------*/

    static inline Real Pyramid5Volume(const ItemWithNodes& item,const VariableNodeReal3& n){
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

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the barycenter (centroid) of a triangle.
     *
     * This method calculates the barycenter of a triangle defined by three nodes.
     * The barycenter is computed as the average of the vertices' coordinates.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeBaryCenterTria3(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];

      Real Center_x = (vertex0.x + vertex1.x + vertex2.x) / 3.;
      Real Center_y = (vertex0.y + vertex1.y + vertex2.y) / 3.;

      return { Center_x, Center_y, 0 };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the length of the edge defined by a given face.
     *
     * This method calculates Euclidean distance between the two nodes of the face.
     */
    /*---------------------------------------------------------------------------*/

    // ef:
    // - needing to compute edge length with different entity inputs (edge, face or cell)
    // - simplifying using available arcane functions
    /* ef    static inline Real computeLengthEdge2(Face face, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[face.nodeId(0)];
      Real3 vertex1 = node_coord[face.nodeId(1)];

      Real dx = vertex1.x - vertex0.x;
      Real dy = vertex1.y - vertex0.y;

      return math::sqrt(dx * dx + dy * dy);
      }*/
    static inline Real computeLengthEdge2(const ItemWithNodes& item, const VariableNodeReal3& node_coord)
    {
      const Real3& vertex0 = node_coord[item.node(0)];
      const Real3& vertex1 = node_coord[item.node(1)];
      return (vertex1-vertex0).normL2();
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the normalized edge normal for a given face.
     *
     * This method calculates normal vector to the edge defined by nodes of the face,
     * normalizes it, and ensures the correct orientation.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real2 computeNormalEdge2(const Face& face, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[face.node(0)];
      Real3 vertex1 = node_coord[face.node(1)];

      if (!face.isSubDomainBoundaryOutside())
        std::swap(vertex0, vertex1);

      Real dx = vertex1.x - vertex0.x;
      Real dy = vertex1.y - vertex0.y;
      Real norm_N = math::sqrt(dx * dx + dy * dy);

      return { dy / norm_N, -dx / norm_N };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the factor used for integration of 1D, 2D and 3D finite elements
     *
     * This method calculates length or surface for a given finite-element (P1, P2, ...) and returns
     * the associated factor for elementary integrals.
     */
    /*---------------------------------------------------------------------------*/
    static inline Real computeFacLengthOrArea(const Face& face, const VariableNodeReal3& node_coord)
    {
      Int32 item_type = face.type();
      Real fac_el{0.};

      switch (item_type) {

      // Lines
      case IT_Line2:
      case IT_Line3:
        fac_el = computeLengthEdge2(face, node_coord) / 2.;
        break;

      // Faces
      case IT_Triangle3:
      case IT_Triangle6:
        fac_el = Tri3Surface(face, node_coord) / 3.;
        break;

      case IT_Quad4:
      case IT_Quad8:
        fac_el = Quad4Surface(face, node_coord) / 4.;
        break;

      default:
        break;
      }
      return fac_el;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes finite-element entity (Edge, Face or Cell) geometric dimension
     * This method is used for the FEM 2D & 3D needs (coming from PASSMO)
     * for jacobians & elementary matrices computations
     */
    /*---------------------------------------------------------------------------*/
    static inline Int32 getGeomDimension(const ItemWithNodes& item){
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
      case IT_Hexaedron20:
      case IT_Pentaedron6:
      case IT_Pyramid5:
        dim = 3; break;

      default: break;

      }
      return dim;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes edge & Face normal & tangent vectors (normalized, direct oriented)
     * This method is used for the FEM 2D & 3D needs (coming from PASSMO)
     * for geometric transformations (rotations, projections, ...)
     * In 2D, it assumes the edge lies in x-y plane (z coord = 0.)
     */
    /*---------------------------------------------------------------------------*/
    static inline void DirVectors(const Face& face,const VariableNodeReal3& n, const Int32& ndim, Real3& e1, Real3& e2, Real3& e3){

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
  };

  /*---------------------------------------------------------------------------*/
  /**
  * @brief Provides methods for finite element operations in 2D.
  *
  * This class includes static methods for calculating gradients of basis
  * functions and integrals for P1 triangles in 2D finite element analysis.
  *
  * Reference Tri3 and Quad4 element in Arcane
  *
  *               0 o                    1 o . . . . o 0
  *                . .                     .         .
  *               .   .                    .         .
  *              .     .                   .         .
  *           1 o  . .  o 2              2 o . . . . o 3
  */
  /*---------------------------------------------------------------------------*/
  class FeOperation2D
  {

   public:

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the gradients of given function U for P1 triangles.
     *
     * This method calculates gradient operator ∇ Ui for a given P1 cell
     * with i = 1,..,3 for the three values of Ui hence  at  cell nodes.
     * The output is ∇ Ui is P0 (piece-wise constant) hence Real3  value
     * per cell
     *
     *         ∇ Ui = [ ∂U/∂x   ∂U/∂y   ∂U/∂z ]
     *
     *         ∂U/∂x = ( u1*(y2 − y3) + u2*(y3 − y1) + u3*(y1 − y2) ) / (2*A)
     *         ∂U/∂y = ( u1*(x3 − x2) + u2*(x1 − x3) + u3*(x2 − x1) ) / (2*A)
     *         ∂U/∂z = 0
     *
     * @note we can adapt the same for 3D by filling the third component
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeGradientTria3(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal& u)
    {
      Real3 n0 = node_coord[cell.nodeId(0)];
      Real3 n1 = node_coord[cell.nodeId(1)];
      Real3 n2 = node_coord[cell.nodeId(2)];

      Real u0 = u[cell.nodeId(0)];
      Real u1 = u[cell.nodeId(1)];
      Real u2 = u[cell.nodeId(2)];

      Real A2 = ((n1.x - n0.x) * (n2.y - n0.y) - (n2.x - n0.x) * (n1.y - n0.y));

      return { (u0*(n1.y - n2.y) + u1*(n2.y - n0.y) + u2*(n0.y - n1.y)) / A2 , (u0*(n2.x - n1.x) + u1*(n0.x - n2.x) + u2*(n1.x - n0.x)) / A2 , 0};
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the gradients of basis functions N for P1 triangles.
     *
     * This method calculates gradient operator ∇ Ni for a given P1 cell
     * with i = 1,..,3 for the three shape function Ni hence output is a
     * matrix of size 3x3
     *
     *              [ ∂N1/∂x     ∂N1/∂y   ∂N1/∂z  ]
     *         ∇N = [ ∂N2/∂x     ∂N2/∂y   ∂N2/∂z  ]
     *              [ ∂N3/∂x     ∂N3/∂y   ∂N3/∂z  ]
     *
     *                     [ y2​−y3     x3​−x2    0  ]
     *         ∇N = 1/(2A) [ y3−y1     x1−x3    0  ]
     *                     [ y1−y2     x2−x1    0  ]
     *
     * @note we can adapt the same for 3D by filling the third component
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3x3 computeGradientTria3(const Cell& cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];

      Real A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

      return {Real3((vertex1.y - vertex2.y) / A2, (vertex2.x - vertex1.x) / A2, 0),
                     Real3((vertex2.y - vertex0.y) / A2, (vertex0.x - vertex2.x) / A2, 0),
                     Real3((vertex0.y - vertex1.y) / A2, (vertex1.x - vertex0.x) / A2, 0)};
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the X gradients of basis functions N for P1 triangles.
     *
     * This method calculates gradient operator ∂/∂x of Ni for a given P1
     * cell with i = 1,..,3 for the three shape function  Ni  hence output
     * is a vector of size 3
     *
     *         ∂N/∂x = [ ∂N1/∂x  ∂N1/∂x  ∂N3/∂x ]
     *
     *         ∂N/∂x = 1/(2A) [ y2​−y3  y3−y1  y1−y2 ]
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeGradientXTria3(const Cell& cell, const VariableNodeReal3& node_coord)
    {
      const Real3& vertex0 = node_coord[cell.node(0)];
      const Real3& vertex1 = node_coord[cell.node(1)];
      const Real3& vertex2 = node_coord[cell.node(2)];

      auto A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

      return {(vertex1.y - vertex2.y) / A2, (vertex2.y - vertex0.y) / A2, (vertex0.y - vertex1.y) / A2};
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the Y gradients of basis functions N for P1 triangles.
     *
     * This method calculates gradient operator ∂/∂y of Ni for a given P1
     * cell with i = 1,..,3 for the three shape function  Ni  hence output
     * is a vector of size 3
     *
     *         ∂N/∂x = [ ∂N1/∂y  ∂N1/∂y  ∂N3/∂y ]
     *
     *         ∂N/∂x = 1/(2A) [ x3​−x2  x1−x3  x2−x1 ]
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeGradientYTria3(const Cell& cell, const VariableNodeReal3& node_coord)
    {
      const Real3& vertex0 = node_coord[cell.node(0)];
      const Real3& vertex1 = node_coord[cell.node(1)];
      const Real3& vertex2 = node_coord[cell.node(2)];

      auto A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

      return {(vertex2.x - vertex1.x) / A2, (vertex0.x - vertex2.x) / A2, (vertex1.x - vertex0.x) / A2};
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the X gradients of basis functions N for P1 Quad.
     *
     * This method calculates gradient operator ∂/∂x of Ni for a given P1
     * cell with i = 1,..,4 for the four shape function  Ni  hence output
     * is a vector of size 4
     *
     *         ∂N/∂x = [ ∂N1/∂x  ∂N1/∂x  ∂N3/∂x  ∂N4/∂x ]
     *
     *         ∂N/∂x = 1/(2A) [ y2​−y3  y3−y0  y0−y1  y1−y2 ]
     */
    /*---------------------------------------------------------------------------*/

    static inline  Real4 computeGradientXQuad4(const Cell& cell, const VariableNodeReal3& node_coord)
    {
      const Real3& vertex0 = node_coord[cell.node(0)];
      const Real3& vertex1 = node_coord[cell.node(1)];
      const Real3& vertex2 = node_coord[cell.node(2)];
      const Real3& vertex3 = node_coord[cell.node(3)];

      auto A2 = ((vertex1.x * vertex2.y + vertex2.x * vertex3.y + vertex3.x * vertex0.y + vertex0.x * vertex1.y) - (vertex2.x * vertex1.y + vertex3.x * vertex2.y + vertex0.x * vertex3.y + vertex1.x * vertex0.y));

      Real4 dx{};

      dx[0] = (vertex2.y - vertex3.y) / A2;
      dx[1] = (vertex3.y - vertex0.y) / A2;
      dx[2] = (vertex0.y - vertex1.y) / A2;
      dx[3] = (vertex1.y - vertex2.y) / A2;  

      return dx;
    }

    /*---------------------------------------------------------------------------*
    /**
     * @brief Computes the Y gradients of basis functions N for P1 Quad.
     *
     * This method calculates gradient operator ∂/∂y of Ni for a given P1
     * cell with i = 1,..,4 for the three shape function  Ni  hence output
     * is a vector of size 4
     *
     *         ∂N/∂x = [ ∂N1/∂y  ∂N1/∂y  ∂N3/∂y  ∂N4/∂y ]
     *
     *         ∂N/∂x = 1/(2A) [ x3​−x2  x0−x3  x1−x0  x2−x1 ]
     */
    /*---------------------------------------------------------------------------*/

    static inline  Real4 computeGradientYQuad4(const Cell& cell, const VariableNodeReal3& node_coord)
    {
      const Real3& vertex0 = node_coord[cell.nodeId(0)];
      const Real3& vertex1 = node_coord[cell.nodeId(1)];
      const Real3& vertex2 = node_coord[cell.nodeId(2)];
      const Real3& vertex3 = node_coord[cell.nodeId(3)];

      Real A2 = ((vertex1.x * vertex2.y + vertex2.x * vertex3.y + vertex3.x * vertex0.y + vertex0.x * vertex1.y) - (vertex2.x * vertex1.y + vertex3.x * vertex2.y + vertex0.x * vertex3.y + vertex1.x * vertex0.y));

      Real4 dy{};
      
      dy[0] = (vertex3.x - vertex2.x) / A2;
      dy[1] = (vertex0.x - vertex3.x) / A2;
      dy[2] = (vertex1.x - vertex0.x) / A2;
      dy[3] = (vertex2.x - vertex1.x) / A2;

      return dy;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the integral (u*v) for P1 triangles.
     *
     * here the element matrix will read
     *
     *              [ 1/6     1/12   1/12 ]
     *         a  = [ 1/12    1/6    1/12 ]
     *              [ 1/12    1/12   1/6  ]
     *
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3x3 computeUVTria3(const Cell& cell, const VariableNodeReal3& node_coord)
    {
      Real aii = 1. / 6.;
      Real aij = 1. / 12.;
      return {Real3(aii, aij, aij), Real3(aij, aii, aij), Real3(aij, aij, aii)};
    }
  };

  /*---------------------------------------------------------------------------*/
  /**
  * @brief Provides methods for finite element operations in 3D.
  *
  * This class includes static methods for calculating gradients of basis
  * functions and integrals for P1 triangles in 3D finite element analysis.
  *
  * Reference Tetra4 element in Arcane
  *
  *               3 o
  *                /|\
  *               / | \
  *              /  |  \
  *             /   o 2 \
  *            / .    .  \
  *         0 o-----------o 1
  */
  /*---------------------------------------------------------------------------*/
  class FeOperation3D
  {
   public:
    /*-------------------------------------------------------------------------*/
    /**
     * @brief Computes the X gradients of basis functions N for P1 Tetrahedron.
     *
     * This method calculates gradient operator ∂/∂x of Ni for a given P1
     * cell with i = 1,..,4 for the four shape function  Ni  hence output
     * is a vector of size 4
     *
     *         ∂N/∂x = [ ∂N1/∂x  ∂N2/∂x  ∂N3/∂x  ∂N4/∂x ]
     *
     *         ∂N/∂x = 1/(6V) [ b0  b1  b2  b3 ]
     *
     * where:
     *    b0 = (m1.y * (m3.z - m2.z) + m2.y * (m1.z - m3.z) + m3.y * (m2.z - m1.z)),
     *    b1 = (m0.y * (m2.z - m3.z) + m2.y * (m3.z - m0.z) + m3.y * (m0.z - m2.z)),
     *    b2 = (m0.y * (m3.z - m1.z) + m1.y * (m0.z - m3.z) + m3.y * (m1.z - m0.z)),
     *    b3 = (m0.y * (m1.z - m2.z) + m1.y * (m2.z - m0.z) + m2.y * (m0.z - m1.z)).
     *
     */
    /*-------------------------------------------------------------------------*/

    static inline  Real4 computeGradientXTetra4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];
      Real3 vertex3 = node_coord[cell.nodeId(3)];

      Real3 v0 = vertex1 - vertex0;
      Real3 v1 = vertex2 - vertex0;
      Real3 v2 = vertex3 - vertex0;

      // 6 x Volume of tetrahedron
      Real V6 =  std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

      Real4 dx{};

      dx[0] = (vertex1.y * (vertex3.z - vertex2.z) + vertex2.y * (vertex1.z - vertex3.z) + vertex3.y * (vertex2.z - vertex1.z))/V6;
      dx[1] = (vertex0.y * (vertex2.z - vertex3.z) + vertex2.y * (vertex3.z - vertex0.z) + vertex3.y * (vertex0.z - vertex2.z))/V6;
      dx[2] = (vertex0.y * (vertex3.z - vertex1.z) + vertex1.y * (vertex0.z - vertex3.z) + vertex3.y * (vertex1.z - vertex0.z))/V6;
      dx[3] = (vertex0.y * (vertex1.z - vertex2.z) + vertex1.y * (vertex2.z - vertex0.z) + vertex2.y * (vertex0.z - vertex1.z))/V6; 

      return dx;
    }

    /*-------------------------------------------------------------------------*/
/**
 * @brief Computes the Y gradients of basis functions N for P1 Tetrahedron.
 *
 * This method calculates gradient operator ∂/∂y of Ni for a given P1
 * cell with i = 1,..,4 for the four shape functions Ni, hence the output
 * is a vector of size 4.
 *
 *         ∂N/∂y = [ ∂N1/∂y  ∂N2/∂y  ∂N3/∂y  ∂N4/∂y ]
 *
 *         ∂N/∂y = 1/(6V) [ c0  c1  c2  c3 ]
 *
 * where:
 *    c0 = (m1.z * (m3.x - m2.x) + m2.z * (m1.x - m3.x) + m3.z * (m2.x - m1.x)),
 *    c1 = (m0.z * (m2.x - m3.x) + m2.z * (m3.x - m0.x) + m3.z * (m0.x - m2.x)),
 *    c2 = (m0.z * (m3.x - m1.x) + m1.z * (m0.x - m3.x) + m3.z * (m1.x - m0.x)),
 *    c3 = (m0.z * (m1.x - m2.x) + m1.z * (m2.x - m0.x) + m2.z * (m0.x - m1.x)).
 *
 */
/*-------------------------------------------------------------------------*/

static inline Real4 computeGradientYTetra4(Cell cell, const VariableNodeReal3& node_coord)
{
  Real3 vertex0 = node_coord[cell.nodeId(0)];
  Real3 vertex1 = node_coord[cell.nodeId(1)];
  Real3 vertex2 = node_coord[cell.nodeId(2)];
  Real3 vertex3 = node_coord[cell.nodeId(3)];

  Real3 v0 = vertex1 - vertex0;
  Real3 v1 = vertex2 - vertex0;
  Real3 v2 = vertex3 - vertex0;

  // 6 x Volume of tetrahedron
  Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

  Real4 dy{};

  dy[0] = (vertex1.z * (vertex3.x - vertex2.x) + vertex2.z * (vertex1.x - vertex3.x) + vertex3.z * (vertex2.x - vertex1.x))/V6;
  dy[1] = (vertex0.z * (vertex2.x - vertex3.x) + vertex2.z * (vertex3.x - vertex0.x) + vertex3.z * (vertex0.x - vertex2.x))/V6;
  dy[2] = (vertex0.z * (vertex3.x - vertex1.x) + vertex1.z * (vertex0.x - vertex3.x) + vertex3.z * (vertex1.x - vertex0.x))/V6;
  dy[3] = (vertex0.z * (vertex1.x - vertex2.x) + vertex1.z * (vertex2.x - vertex0.x) + vertex2.z * (vertex0.x - vertex1.x))/V6;

  return dy;
}

/*-------------------------------------------------------------------------*/
/**
 * @brief Computes the Z gradients of basis functions N for P1 Tetrahedron.
 *
 * This method calculates gradient operator ∂/∂z of Ni for a given P1
 * cell with i = 1,..,4 for the four shape functions Ni, hence the output
 * is a vector of size 4.
 *
 *         ∂N/∂z = [ ∂N1/∂z  ∂N2/∂z  ∂N3/∂z  ∂N4/∂z ]
 *
 *         ∂N/∂z = 1/(6V) [ d0  d1  d2  d3 ]
 *
 * where:
 *    d0 = (m1.x * (m3.y - m2.y) + m2.x * (m1.y - m3.y) + m3.x * (m2.y - m1.y)),
 *    d1 = (m0.x * (m2.y - m3.y) + m2.x * (m3.y - m0.y) + m3.x * (m0.y - m2.y)),
 *    d2 = (m0.x * (m3.y - m1.y) + m1.x * (m0.y - m3.y) + m3.x * (m1.y - m0.y)),
 *    d3 = (m0.x * (m1.y - m2.y) + m1.x * (m2.y - m0.y) + m2.x * (m0.y - m1.y)).
 *
 */
/*-------------------------------------------------------------------------*/

static inline Real4 computeGradientZTetra4(const Cell& cell, const VariableNodeReal3& node_coord)
{
  const Real3& vertex0 = node_coord[cell.nodeId(0)];
  const Real3& vertex1 = node_coord[cell.nodeId(1)];
  const Real3& vertex2 = node_coord[cell.nodeId(2)];
  const Real3& vertex3 = node_coord[cell.nodeId(3)];

  auto v0 = vertex1 - vertex0;
  auto v1 = vertex2 - vertex0;
  auto v2 = vertex3 - vertex0;

  // 6 x Volume of tetrahedron
  Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

  Real4 dz{};

  dz[0] = (vertex1.x * (vertex3.y - vertex2.y) + vertex2.x * (vertex1.y - vertex3.y) + vertex3.x * (vertex2.y - vertex1.y))/V6;
  dz[1] = (vertex0.x * (vertex2.y - vertex3.y) + vertex2.x * (vertex3.y - vertex0.y) + vertex3.x * (vertex0.y - vertex2.y))/V6;
  dz[2] = (vertex0.x * (vertex3.y - vertex1.y) + vertex1.x * (vertex0.y - vertex3.y) + vertex3.x * (vertex1.y - vertex0.y))/V6;
  dz[3] = (vertex0.x * (vertex1.y - vertex2.y) + vertex1.x * (vertex2.y - vertex0.y) + vertex2.x * (vertex0.y - vertex1.y))/V6;

  return dz;
}

  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for applying boundary conditions in 2D FEM problems.
   *
   * This class includes static methods for applying Neumann boundary conditions
   * to the right-hand side (RHS) of finite element method equations in 2D.
   */
  /*---------------------------------------------------------------------------*/
  class BoundaryConditions2D
  {
   public:

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies a constant source term to the RHS vector.
     *
     * This method adds a constant source term `qdot` to the RHS vector for each 
     * node in the mesh. The contribution to each node is weighted by the area of 
     * the cell and evenly distributed among the number of nodes of the cell.
     *
     * @param [IN]  qdot       : The constant source term.
     * @param [IN]  mesh       : The mesh containing all cells.
     * @param [IN]  node_dof   : DOF connectivity view.
     * @param [IN]  node_coord : The coordinates of the nodes.
     * @param [OUT] rhs_values : The RHS values to update.
     */
    /*---------------------------------------------------------------------------*/

    static inline void applyConstantSourceToRhs(const Real& qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
    {
      ENUMERATE_ (Cell, icell, mesh->allCells()) {
        Cell cell = *icell;
        Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, node_coord);
        for (Node node : cell.nodes()) {
          if (node.isOwn())
            rhs_values[node_dof.dofId(node, 0)] += qdot * area / cell.nbNode();
        }
      }
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies a manufactured source term to the RHS vector.
     *
     * This method adds a manufactured source term to the RHS vector for each 
     * node in the mesh. The contribution to each node is weighted by the area of 
     * the cell and evenly distributed among the nodes of the cell.
     *
     * @param [IN]  qdot       : The constant source term.
     * @param [IN]  mesh       : The mesh containing all cells.
     * @param [IN]  node_dof   : DOF connectivity view.
     * @param [IN]  node_coord : The coordinates of the nodes.
     * @param [OUT] rhs_values : The RHS values to update.
     */
    /*---------------------------------------------------------------------------*/

    static inline void applyManufacturedSourceToRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_source, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
    {
      ENUMERATE_ (Cell, icell, mesh->allCells()) {
        Cell cell = *icell;
        Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, node_coord);
        Real3 bcenter = ArcaneFemFunctions::MeshOperation::computeBaryCenterTria3(cell, node_coord);

        for (Node node : cell.nodes()) {
          if (node.isOwn())
            rhs_values[node_dof.dofId(node, 0)] += manufactured_source->apply(area / cell.nbNode(), bcenter);
        }
      }
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies Neumann conditions to the right-hand side (RHS) values.
     *
     * This method updates the RHS values of the finite element method equations
     * based on the provided Neumann boundary condition. The boundary condition
     * can specify a value or its components along the x and y directions.
     *
     * @param [IN]  bs         : The Neumann boundary condition values.
     * @param [IN]  node_dof   : Connectivity view for degrees of freedom at nodes.
     * @param [IN]  node_coord : Coordinates of the nodes in the mesh.
     * @param [OUT] rhs_values : The right-hand side values to be updated.
     */
    /*---------------------------------------------------------------------------*/

    static inline void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
    {
      FaceGroup group = bs->getSurface();

      Real value = 0.0;
      Real valueX = 0.0;
      Real valueY = 0.0;
      bool hasValue = bs->hasValue();
      bool hasValueX = bs->getValueX();
      bool hasValueY = bs->getValueY();

      if (hasValue) {
        value = bs->getValue();
      }
      else {
        if (hasValueX)
          valueX = bs->getValueX();
        if (hasValueY)
          valueY = bs->getValueY();
      }

      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;

        Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, node_coord);
        Real2 normal = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, node_coord);

        for (Node node : iface->nodes()) {
          if (!node.isOwn())
            continue;
          Real rhs_value = 0.0;

          if (hasValue) {
            rhs_value = value * length / 2.0;
          }
          else {
            rhs_value = (normal.x * valueX + normal.y * valueY) * length / 2.0;
          }

          rhs_values[node_dof.dofId(node, 0)] += rhs_value;
        }
      }
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies Dirichlet boundary conditions to RHS and LHS.
     *
     * Updates the LHS matrix and RHS vector to enforce Dirichlet conditions.
     *
     * - For LHS matrix `A`, the diagonal term for the Dirichlet DOF is set to `P`.
     * - For RHS vector `b`, the Dirichlet DOF term is scaled by `P`.
     *
     * @param [IN]  bs              : Boundary condition values.
     * @param [IN]  node_dof        : DOF connectivity view.
     * @param [IN]  node_coord      : Node coordinates.
     * @param [OUT] m_linear_system : Linear system for LHS.
     * @param [OUT] rhs_values RHS  : RHS values to update.
     */
    /*---------------------------------------------------------------------------*/
    static inline void applyDirichletToLhsAndRhs(BC::IDirichletBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, DoFLinearSystem& m_linear_system, VariableDoFReal& rhs_values)
    {
      FaceGroup group = bs->getSurface();
      Real value = bs->getValue();
      Real Penalty = bs->getPenalty();

      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          if (node.isOwn()) {
            m_linear_system.matrixSetValue(node_dof.dofId(node, 0), node_dof.dofId(node, 0), Penalty);
            Real u_g = Penalty * value;
            rhs_values[node_dof.dofId(node, 0)] = u_g;
          }
        }
      }
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies Point Dirichlet boundary conditions to RHS and LHS.
     *
     * Updates the LHS matrix and RHS vector to enforce the Dirichlet.
     *
     * - For LHS matrix `A`, the diagonal term for the Dirichlet DOF is set to `P`.
     * - For RHS vector `b`, the Dirichlet DOF term is scaled by `P`.
     *
     * @param [IN]  bs              : Boundary condition values.
     * @param [IN]  node_dof        : DOF connectivity view.
     * @param [IN]  node_coord      : Node coordinates.
     * @param [OUT] m_linear_system : Linear system for LHS.
     * @param [OUT] rhs_values RHS  : RHS values to update.
     */
    /*---------------------------------------------------------------------------*/
    static inline void applyPointDirichletToLhsAndRhs(BC::IDirichletPointCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, DoFLinearSystem& m_linear_system, VariableDoFReal& rhs_values)
    {
      NodeGroup group = bs->getNode();
      Real value = bs->getValue();
      Real Penalty = bs->getPenalty();

      ENUMERATE_ (Node, inode, group) {
        Node node = *inode;
        if (node.isOwn()) {
          m_linear_system.matrixSetValue(node_dof.dofId(node, 0), node_dof.dofId(node, 0), Penalty);
          Real u_g = Penalty * value;
          rhs_values[node_dof.dofId(node, 0)] = u_g;
        }
      }
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies Manufactured Dirichlet boundary conditions to RHS and LHS.
     *
     * Updates the LHS matrix and RHS vector to enforce the Dirichlet.
     *
     * - For LHS matrix `A`, the diagonal term for the Dirichlet DOF is set to `P`.
     * - For RHS vector `b`, the Dirichlet DOF term is scaled by `P`.
     *
     * @param [IN]  manufactured_dirichlet   : External function for Dirichlet.
     * @param [IN]  group           : Group of all external faces.
     * @param [IN]  bs              : Boundary condition values.
     * @param [IN]  node_dof        : DOF connectivity view.
     * @param [IN]  node_coord      : Node coordinates.
     * @param [OUT] m_linear_system : Linear system for LHS.
     * @param [OUT] rhs_values RHS  : RHS values to update.
     */
    /*---------------------------------------------------------------------------*/
    static inline void applyManufacturedDirichletToLhsAndRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_dirichlet, const Real& lambda, const FaceGroup& group, BC::IManufacturedSolution* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, DoFLinearSystem& m_linear_system, VariableDoFReal& rhs_values)
    {
      Real Penalty = bs->getPenalty();

      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          if (node.isOwn()) {
            m_linear_system.matrixSetValue(node_dof.dofId(node, 0), node_dof.dofId(node, 0), Penalty);
            double tt = 1.;
            Real u_g = Penalty * manufactured_dirichlet->apply(tt, node_coord[node]);
            rhs_values[node_dof.dofId(node, 0)] = u_g;
          }
        }
      }
    }
  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods based on the Dispatcher mechanism available in Arcane,
   * allowing to compute FEM methods without declaring the FE entity type
   * (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  class CellFEMDispatcher
  {

   public:

    Real getShapeFuncVal(const Int16& /*item_type*/, const Int32& /*inod*/, const Real3& /*ref coord*/);
    Real3 getShapeFuncDeriv(const Int16& /*item_type*/, const Int32& /*inod*/, const Real3& /*ref coord*/);

    RealUniqueArray getGaussData(const ItemWithNodes& item, const Integer& nint, Int32& ngauss);

    CellFEMDispatcher();

   private:

    std::function<Real(const Int32& inod, const Real3& coord)> m_shapefunc[NB_BASIC_ITEM_TYPE];
    std::function<Real3(const Int32& inod, const Real3& coord)> m_shapefuncderiv[NB_BASIC_ITEM_TYPE];
  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for various FEM-related operations.
   *
   * This class includes static methods for computing shape functions and their
   * derivatives depending on finite element types. These methods are used within
   * the Dispatcher mechanism available in Arcane through the class
   * CellFEMDispatcher (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  class FemShapeMethods {
   public:
    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) edge finite-element
     * The "Line2" reference element is assumed as follows:
     *  0           1
     *  o-----------o---> x
     * -1           1
     * direct local numbering : 0->1
     */
    /*---------------------------------------------------------------------------*/
    static inline Real Line2ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
      assert(inod >= 0 && inod < 2);
#endif

      Real r = ref_coord[0];
      if (inod == 1) return (0.5*(1 + r));
      return (0.5*(1 - r));
    }

    static inline Real3 Line2ShapeFuncDeriv(const Integer& inod,const Real3&){
      if (inod == 1) return { 0.5,0.,0. };
      return { -0.5,0.,0. };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) edge finite-element
     * The "Line3" reference element is assumed as follows:
     *  0     2      1
     *  o-----o------o---> x
     * -1     0      1
     * direct local numbering : 0->1->2
     */
    /*---------------------------------------------------------------------------*/
    static inline Real Line3ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
      assert(inod >= 0 && inod < 3);
#endif

      Real ri = ref_coord[0];
      if (inod == 0) ri *= -1;

      if (inod < 2) return 0.5*ri*(1 + ri); // nodes 0 or 1
      return (1 - ri*ri); // middle node
    }

    static inline Real3 Line3ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
      Real ri = ref_coord[0];
      if (!inod) return {-0.5 + ri, 0.,0.};
      if (inod == 1) return {0.5 + ri, 0.,0.};
      return {-2.*ri, 0.,0.};
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) triangle finite-element
     * The "Tri3" reference element is assumed as follows:
     *        0 o
     *         . .
     *        .   .
     *       .     .
     *      .       .
     *     .         .
     *    .           .
     *   o-------------o
     *  1               2
     * direct local numbering : 0->1->2
     */
    /*---------------------------------------------------------------------------*/
    static inline Real Tri3ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
      assert(inod >= 0 && inod < 3);
#endif
      Real r = ref_coord[0];
      Real s = ref_coord[1];
      if (!inod) return (1 - r - s);
      if (inod == 1) return r;
      return s;
    }

    static inline Real3 Tri3ShapeFuncDeriv(const Integer& inod,const Real3&){
      if (!inod) return {-1., -1.,0.};
      if (inod==1) return {1., 0.,0.};
      return {0., 1., 0.};
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) triangle finite-element
     * The "Tri6" reference element is assumed as follows:
     *        0 o
     *         . .
     *        .   .
     *     3 o     o 5
     *      .       .
     *     .         .
     *    .           .
     *   o------o------o
     *  1       4       2
     * direct local numbering : 0->1->2->3->4->5
     */
    /*---------------------------------------------------------------------------*/
    static inline Real Tri6ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
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

    static inline Real3 Tri6ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
      if (!inod) return {1., 0.,0.};
      if (inod==1) return {0., 1., 0.};
      if (inod == 2) return {-1., -1.,0.};

      auto	ri = ref_coord[0],si = ref_coord[1],ti = ref_coord[2];
      if (inod == 3) return {4.*si, 4.*ri,0.};
      if (inod == 4) return {-4.*si, 4.*(ti - si),0.};
      return {4.*(ti - ri), -4.*ri,0.};
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) quadrangle finite-element
     * The "Quad4" reference element is assumed as follows:
     *         ^y
     *          |
     *  1 o-----1-----o 0
     *    |     |     |
     *    |     |     |
     *    |     |     |
     *   -1 ----|---- 1 ---> x
     *    |     |     |
     *    |     |     |
     *    |     |     |
     *  2 o--- -1 ----o 3
     *
     * direct local numbering : 0->1->2->3
     */
    /*---------------------------------------------------------------------------*/
    static inline Real Quad4ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
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

    static inline Real3 Quad4ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
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

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) quadrangle finite-element
     * The "Quad8" reference element is assumed as follows:
     *         ^y
     *          |
     *          4
     *  1 o-----o-----o 0
     *    |     |     |
     *    |     |     |
     *    |     |     |
     *  5 o ----|---- o 7 ---> x
     *    |     |     |
     *    |     |     |
     *    |     |     |
     *  2 o-----o-----o 3
     *          6
     * Normalized coordinates (x, y) vary between -1/+1
     * Nodes 4, 6 are on line (x = 0)
     * direct local numbering :  0->1->2->...->5->6->7
     */
    /*---------------------------------------------------------------------------*/
    static inline Real Quad8ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
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

    static inline Real3 Quad8ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

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

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) hexaedron finite-element
     * The "Hexa8" reference element is assumed as follows:
     *     (-1, 1,1)
     *         1-------------0 (1,1,1)
     *        /|            /|
     *       / |           / |
     *     /   |          /  |
     *    2----|---------3   |   z   y
     *  (-1,-1,1)        |   |   | /
     *    |    |         |   |   |/--->x
     *    |    |         |   |
     *    |    |         |   |
     *    |    5---------|---4 (1,1,-1)
     *    |  /           |  /
     *    | /            | /
     *    |/             |/
     *    6--------------7 (1,-1,-1)
     * (-1,-1,-1)
     * Normalized coordinates (x, y, z) vary between -1/+1
     * direct local numbering : 0->1->2->3->4->5->6->7
     */
    /*---------------------------------------------------------------------------*/
    static inline Real Hexa8ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
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

    static inline Real3 Hexa8ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

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

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) hexaedron finite-element
     * The "Hexa20" reference element is assumed as follows:
     *     (-1, 1,1)
     *         1------8------0 (1,1,1)
     *        /|            /|
     *      9  |          11 |
     *     /   |          /  |
     *    2----|--10-----3   |   z   y
     *  (-1,-1,1)        |   |   | /
     *    |   16         |  19   |/--->x
     *    |    |         |   |
     *   17    |        18   |
     *    |    5----12---|---4 (1,1,-1)
     *    |  /           |  /
     *    | 13           | 15
     *    |/             |/
     *    6-----14-------7 (1,-1,-1)
     * (-1,-1,-1)
     * Normalized coordinates (x, y, z) vary between -1/+1
     * direct local numbering : 0->1->2->3->...->18->19
     */
    /*---------------------------------------------------------------------------*/

    static inline Real Hexa20ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
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

    static inline Real3 Hexa20ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

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

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) tetrahedral finite-element
     * The "Tetra4" reference element is assumed as follows:
     *
     *    (0,0,1)                     3
     *       .                        *.*
     *       .                        * . *
     *       .                        *  .  *
     *       .                        *   .   *
     *       Z   (0,1,0)              *    .    *
     *       .    .                   *     2     *
     *       .   .                    *   .    .    *
     *       .  Y                     *  .        .   *
     *       . .                      * .            .  *
     *       ..           (1,0,0)     *.                . *
     *       --------X------>         0********************1
     *
     * direct local numbering : 0->1->2->3
     */
    /*---------------------------------------------------------------------------*/

    static inline Real Tetra4ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
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

    static inline Real3 Tetra4ShapeFuncDeriv(const Integer& inod,const Real3& /*ref_coord*/){

      if (inod == 3) return {0.,0.,1.};
      if (inod == 1) return {1.,0.,0.};
      if (inod == 2) return {0.,1.,0.};
      return {-1.,-1.,-1.};
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) tetrahedral finite-element
     * The "Tetra10" reference element is assumed as follows:
     *
     *    (0,0,1)                     x 3
     *       .                        *.*
     *       .                        * . *
     *       .                        *  .  *
     *       .                        *  9   *
     *       Z   (0,1,0)              *    .    *
     *       .    .                   *     x     8
     *       .   .                    7   . 2  .    *
     *       .  Y                     *  6       5    *
     *       . .                      * .            .  *
     *       ..           (1,0,0)     *.                . *
     *       --------X------>       0 x ******* 4 ******** x 1
     *
     * direct local numbering : 0->1->2->3->...->8->9
     */
    /*---------------------------------------------------------------------------*/

    static inline Real Tetra10ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
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

    static inline Real3 Tetra10ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){
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

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) pentaedron finite-element
     * The "Penta6" reference element is assumed as follows:
     *
     *                     5 (0,1,1)
     *                   . |  .
     *                  .  |     .
     *                 .   Z        .
     *                .    |           .
     *               .     |             .
     *       (0,0,1) 3 ------------------ 4 (1,0,1)
     *               |     |              |
     *               |     |              |
     *               |     |              |
     *               |     |              |
     *               |     2 (0,1,-1)     |
     *               |   .    .           |
     *               |  Y        .        |
     *               | .            .     |
     *               |.                .  |
     *      (0,0,-1) 0 -------- X ------- 1 (1,0,-1)
     *
     * direct local numbering : 0->1->2->3->4->5->6
     */
    /*---------------------------------------------------------------------------*/

    static inline Real Penta6ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
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

    static inline Real3 Penta6ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

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
    /**
     * @brief Provides methods for reference linear (P1) pyramid finite-element
     * The "Pyramid5" reference element is assumed as follows:
     *
     *                               ^
     *                               |
     *                               Z
     *                               |
     *                               4 (0,0,1)
     *                              *
     *                             * **
     *                            ** *  *
     *                           * * |*   *
     *                          * *  | *    *
     *                         * *   |  *     *
     *                        * *    |   *      *      .Y
     *                       * *     |    *       *  .
     *            (-1,0,0)  * 2 -----|-----*------ 1 (0,1,0)
     *                     * .       |      *  .  .
     *                    * .        |     . *   .
     *                   *.             .     * .
     *                  *                  X . *
     *        (0,-1,0) 3 --------------------- 0 (1,0,0)
     *
     * direct local numbering : 0->1->2->3->4
     */
    /*---------------------------------------------------------------------------*/

    static inline Real Pyramid5ShapeFuncVal(const Integer& inod,const Real3& ref_coord){
#ifdef _DEBUG
      assert(inod >= 0 && inod < 5);
#endif
      auto	r{ ref_coord[0] },s{ ref_coord[1] },t{ ref_coord[2] } ;
      auto	r1{-1.},s1{1.}, r2{-1.},s2{-1.};

      if (inod == 4) return t;
      auto ti{t - 1.};
      auto t0{0.};

      if (fabs(ti) < REL_PREC)
        ti = 0.;
      else
        t0 = -1./ti / 4.;

      switch(inod){
      case 1:	s1 = -1.; r2 = 1.; break;
      case 2:	r1 = 1.; r2 = 1.; break;
      case 3:	r1 = 1.; s2 = 1.; break;
      default: break;// default is for node 0
      }

      return (r1*r + s1*s + ti) * (r2*r + s2*s + ti) * t0;
    }

    static inline Real3 Pyramid5ShapeFuncDeriv(const Integer& inod,const Real3& ref_coord){

#ifdef _DEBUG
      assert(inod >= 0 && inod < 5);
#endif
      auto	r{ ref_coord[0] },s{ ref_coord[1] },t{ ref_coord[2] } ;
      auto	r1{-1.},s1{1.}, r2{-1.},s2{-1.};

      auto ti{t - 1.};
      auto t0{0.};

      if (fabs(ti) < REL_PREC)
        ti = 0.;
      else
        t0 = -1./ti / 4.;

      switch(inod){
      case 1:	s1 = -1.; r2 = 1.; break;
      case 2:	r1 = 1.; r2 = 1.; break;
      case 3:	r1 = 1.; s2 = 1.; break;
      default: break;// default is for node 0
      }

      if (inod == 4) return {0.,0.,1.};

      Real3 dPhi;
      auto r12{r1+r2}, rr{2.*r1*r2}, s12{s1+s2}, ss{2.*s1*s2}, rs{r1*s2 + r2*s1}, t02{4.*t0*t0};

      dPhi.x = t0 * (rr*r + rs*s + r12*ti);
      dPhi.y = t0 * (rs*r + ss*s + s12*ti);

      if (fabs(ti) < REL_PREC) dPhi.z = 0.;
      else
        dPhi.z = t0 * (r12*r + s12*s + 2.*ti) + t02 * (r1*r + s1*s + ti) * (r2*r + s2*s + ti);

      return dPhi;
    }

    /*---------------------------------------------------------------------------*/

  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for Gauss quadrature.
   *
   * This class includes static methods for computing Gauss-Legendre integration
   * depending on finite element types (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  class FemGaussQuadrature {

   public:
  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the number of Gauss Points for a given finite element type,
   * depending on the integration order chosen by user (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  static inline Integer getNbGaussPointsfromOrder(const Int16& cell_type, const Integer& ninteg){
      Integer nbgauss{0};
      auto ninteg2{ ninteg * ninteg };
      auto ninteg3{ ninteg2 * ninteg };

      if (ninteg <= 1)
        nbgauss = 1;
      else if (cell_type == IT_Line2 || cell_type == IT_Line3)
        nbgauss = ninteg;
      else if (cell_type == IT_Quad4 || cell_type == IT_Quad8)
        nbgauss = ninteg2;
      else if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20)
        nbgauss = ninteg3;

      else if (ninteg == 2) {
        switch (cell_type) {
        default: break;

        case IT_Triangle3:
        case IT_Triangle6: nbgauss = 3;
          break;

        case IT_Tetraedron4:
        case IT_Tetraedron10: nbgauss = 4;
          break;

        case IT_Pentaedron6: nbgauss = 6;
          break;

        case IT_Pyramid5: nbgauss = 5;
          break;
        }
      }
      else if (ninteg == 3) {
        switch (cell_type) {
        default: break;

        case IT_Triangle3:
        case IT_Triangle6: nbgauss = 4;
          break;

        case IT_Tetraedron4:
        case IT_Tetraedron10: nbgauss = 5;
          break;

        case IT_Pentaedron6: nbgauss = 8;
          break;

        case IT_Pyramid5: nbgauss = 6;
          break;
        }
      }
      else if (ninteg >= 4) {
        switch (cell_type) {
        default: break;

        case IT_Triangle3:
        case IT_Triangle6: nbgauss = 7;
          break;

        case IT_Tetraedron4:
        case IT_Tetraedron10: nbgauss = 15;
          break;

        case IT_Pentaedron6: nbgauss = 21;
          break;

        case IT_Pyramid5: nbgauss = 27;
          break;
        }
      }
      return nbgauss;
    }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the position of a Gauss Point in the reference (local) element
   * for a given finite element, depending on the rank of the point in the
   * element loop (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 getGaussRefPosition(const ItemWithNodes& cell, const Integer& ninteg, const Int32& rank)
    {
      auto cell_type = cell.type();
      Integer nint{ninteg};

      if (nint < 1)
        nint = 1;

      if (cell_type == IT_Line2 || cell_type == IT_Line3)
        return LineRefPosition({rank,-1,-1},{nint,0,0});

      if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
        auto in{0};
        for (Int32 i1 = 0; i1 < nint; ++i1) {
          for (Int32 i2 = 0; i2 < nint; ++i2) {

            if (rank == in)
              return QuadRefPosition({i1,i2,-1},{nint,nint,0});

            ++in;
          }
        }
      }

      if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
        auto in{0};
        for (Int32 i1 = 0; i1 < nint; ++i1) {
          for (Int32 i2 = 0; i2 < nint; ++i2) {
            for (Int32 i3 = 0; i3 < nint; ++i3) {

              if (rank == in)
                return HexaRefPosition({ i1, i2, i3 }, { nint, nint, nint });

              ++in;
            }
          }
        }
      }

      if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
        auto o {3};
        if (nint <= 3)
          o = nint-1;

        return {xg1[o][rank],xg2[o][rank],0.};
      }

      if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10) {
        auto o {3};
        if (nint <= 3)
          o = nint-1;

        return {xtet[o][rank],ytet[o][rank],ztet[o][rank]};
      }

      if (cell_type == IT_Pyramid5) {
        auto o {1};
        if (nint <= 2)
          o = 0;

        return {xpyr[o][rank],ypyr[o][rank],zpyr[o][rank]};
      }

      if (cell_type == IT_Pentaedron6) {
        auto o {1};
        if (nint <= 2)
          o = 0;

        return {xpent[o][rank],ypent[o][rank],zpent[o][rank]};
      }
      return {};
    }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for a given
   * finite element, depending on the rank of the point in the element loop
   * (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  static inline Real getGaussWeight(const ItemWithNodes& cell, const Integer& ninteg, const Int32& rank)
    {
      auto cell_type = cell.type();
      Integer nint{ninteg};

      if (nint < 1)
        nint = 1;

      if (cell_type == IT_Line2 || cell_type == IT_Line3)
        return LineWeight({rank,-1,-1},{nint,0,0});

      if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
        auto in{0};
        for (Int32 i1 = 0; i1 < nint; ++i1) {
          for (Int32 i2 = 0; i2 < nint; ++i2) {

            if (rank == in)
              return QuadWeight({i1,i2,-1},{nint,nint,0});

            ++in;
          }
        }
      }

      if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
        auto in{0};
        for (Int32 i1 = 0; i1 < nint; ++i1) {
          for (Int32 i2 = 0; i2 < nint; ++i2) {
            for (Int32 i3 = 0; i3 < nint; ++i3) {

              if (rank == in)
                return HexaWeight({ i1, i2, i3 }, { nint, nint, nint });

              ++in;
            }
          }
        }
      }

      if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
        auto o {3};
        if (nint <= 3)
          o = nint-1;

        return wg[o][rank];
      }

      if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10) {
        if (nint == 1)
          return wgtet1;
        if (nint == 2)
          return wgtet2;
        if (nint == 3) {
          auto i = npwgtet3[rank];
          return wgtet3[i];
        }
        // nint >= 4
        auto i = npwgtet4[rank];
        return wgtet4[i];
      }

      if (cell_type == IT_Pyramid5) {
        if (nint <= 2)
          return wgpyr2;

        // nint >= 3
        auto i = npwgpyr3[rank];
        return wgpyr3[i];
      }

      if (cell_type == IT_Pentaedron6) {
        if (nint <= 2)
          return wgpent2;

        // nint >= 3
        auto i = npwgpent3[rank];
        return wgpent3[i];
      }

      return 1.;
    }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the position of a Gauss Point in the reference (local) element
   * coordinates, depending on the index of the Gauss Point and integration order
   * chosen by user (coming fro PASSMO).
   * This method is generic (not depending on the finite element type) and is called
   * by specialized methods (depending on FE type)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real getRefPosition(const Integer& indx,const Integer& ordre){
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

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight depending on the index of the Gauss Point
   * and integration order chosen by user (coming fro PASSMO).
   * This method is generic (not depending on the finite element type) and is called
   * by specialized methods (depending on FE type)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real getWeight(const Integer& indx,const Integer& ordre){
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
  /**
   * @brief Provides the reference coordinates of a Gauss Point for edge elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/

  static inline Real3 LineRefPosition(const Integer3& indices,const Integer3& ordre){
    return {getRefPosition(indices[0],ordre[0]),0.,0.};
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for triangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 TriRefPosition(const Integer3& indices,const Integer3& ordre){
    Integer o = ordre[0]-1;
    Integer i = indices[0];
    return {xg1[o][i],xg2[o][i],0.};
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for quadrangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 QuadRefPosition(const Integer3& indices,const Integer3& ordre){
    Real3 pos;
    for (Integer i = 0; i < 2; i++) pos[i] = getRefPosition(indices[i],ordre[i]);
    return pos;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for hexaedron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 HexaRefPosition(const Integer3& indices,const Integer3& ordre){
    Real3 pos;
    for (Integer i = 0; i < 3; i++) pos[i] = getRefPosition(indices[i],ordre[i]);
    return pos;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for tetraedron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 TetraRefPosition(const Integer3& indices,const Integer3& /*ordre*/){
    Integer i = indices[0];
    return {xit[i],yit[i],zit[i]};
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for wedge (pentaedron)
   * elements (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 PentaRefPosition(const Integer3& indices,const Integer3& ordre){

    // Same as TriRefPosition on reference coordinate plane (r,s)
    // and LineRefPosition along reference coordinate t (vertical)
    auto pos = TriRefPosition(indices,ordre);
    pos.z = getRefPosition(indices[2],ordre[2]);

    return pos;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for edge elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real LineWeight(const Integer3& indices,const Integer3& ordre){
    return getWeight(indices[0],ordre[0]);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for triangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real TriWeight(const Integer3& indices,const Integer3& ordre){
    return wg[ordre[0]-1][indices[0]];
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for quadrangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real QuadWeight(const Integer3& indices,const Integer3& ordre){
    Real w = 1.;
    for (Integer i = 0; i < 2; i++) w *= getWeight(indices[i],ordre[i]);
    return w;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for hexadron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real HexaWeight(const Integer3& indices,const Integer3& ordre){
    Real w = 1.;
    for (Integer i = 0; i < 3; i++) w *= getWeight(indices[i],ordre[i]);
    return w;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for tetraedron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real TetraWeight(const Integer3& indices,const Integer3& /*ordre*/){
    return wgtetra;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for wedge (pentaedron)
   * elements (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real PentaWeight(const Integer3& indices,const Integer3& ordre){

    // Same as TriWeight on reference coordinate plane (r,s)
    // and LineWeight with ordre[2] to account for reference coordinate t (vertical)
    Real wgpenta = TriWeight(indices,ordre)*getWeight(indices[2],ordre[2]);
    return wgpenta;
  }

  /*---------------------------------------------------------------------------*/

  };


};

#endif // ARCANE_FEM_FUNCTIONS_H
