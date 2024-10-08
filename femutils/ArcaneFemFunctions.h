#ifndef ARCANE_FEM_FUNCTIONS_H
#define ARCANE_FEM_FUNCTIONS_H

#include <arcane/core/IStandardFunction.h>

using namespace Arcane;

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

    static inline Real computeVolumeTetra4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];
      Real3 vertex3 = node_coord[cell.nodeId(3)];

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
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];

      return 0.5 * ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));
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

    static inline Real computeLengthEdge2(Face face, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[face.nodeId(0)];
      Real3 vertex1 = node_coord[face.nodeId(1)];

      Real dx = vertex1.x - vertex0.x;
      Real dy = vertex1.y - vertex0.y;

      return math::sqrt(dx * dx + dy * dy);
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the normalized edge normal for a given face.
     *
     * This method calculates normal vector to the edge defined by nodes of the face,
     * normalizes it, and ensures the correct orientation.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real2 computeNormalEdge2(Face face, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[face.nodeId(0)];
      Real3 vertex1 = node_coord[face.nodeId(1)];

      if (!face.isSubDomainBoundaryOutside())
        std::swap(vertex0, vertex1);

      Real dx = vertex1.x - vertex0.x;
      Real dy = vertex1.y - vertex0.y;
      Real norm_N = math::sqrt(dx * dx + dy * dy);

      return { dy / norm_N, -dx / norm_N };
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

      return Real3 ( (u0*(n1.y - n2.y) + u1*(n2.y - n0.y) + u2*(n0.y - n1.y)) / A2 , (u0*(n2.x - n1.x) + u1*(n0.x - n2.x) + u2*(n1.x - n0.x)) / A2 , 0);
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

    static inline Real3x3 computeGradientTria3(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];

      Real A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

      return Real3x3(Real3((vertex1.y - vertex2.y) / A2, (vertex2.x - vertex1.x) / A2, 0),
                     Real3((vertex2.y - vertex0.y) / A2, (vertex0.x - vertex2.x) / A2, 0),
                     Real3((vertex0.y - vertex1.y) / A2, (vertex1.x - vertex0.x) / A2, 0));
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

    static inline Real3 computeGradientXTria3(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];

      Real A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

      return Real3((vertex1.y - vertex2.y) / A2, (vertex2.y - vertex0.y) / A2, (vertex0.y - vertex1.y) / A2);
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

    static inline Real3 computeGradientYTria3(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];

      Real A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

      return Real3((vertex2.x - vertex1.x) / A2, (vertex0.x - vertex2.x) / A2, (vertex1.x - vertex0.x) / A2);
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

    static inline  Real4 computeGradientXQuad4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];
      Real3 vertex3 = node_coord[cell.nodeId(3)];

      Real A2 = ((vertex1.x * vertex2.y + vertex2.x * vertex3.y + vertex3.x * vertex0.y + vertex0.x * vertex1.y) - (vertex2.x * vertex1.y + vertex3.x * vertex2.y + vertex0.x * vertex3.y + vertex1.x * vertex0.y));

      Real4 dx;

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

    static inline  Real4 computeGradientYQuad4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];
      Real3 vertex3 = node_coord[cell.nodeId(3)];

      Real A2 = ((vertex1.x * vertex2.y + vertex2.x * vertex3.y + vertex3.x * vertex0.y + vertex0.x * vertex1.y) - (vertex2.x * vertex1.y + vertex3.x * vertex2.y + vertex0.x * vertex3.y + vertex1.x * vertex0.y));

      Real4 dy;
      
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

    static inline Real3x3 computeUVTria3(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real aii = 1. / 6.;
      Real aij = 1. / 12.;
      return Real3x3(Real3(aii, aij, aij), Real3(aij, aii, aij), Real3(aij, aij, aii));
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

      Real4 dx;

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

  Real4 dy;

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

static inline Real4 computeGradientZTetra4(Cell cell, const VariableNodeReal3& node_coord)
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

  Real4 dz;

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

    static inline void applyConstantSourceToRhs(const Real& qdot, IMesh* mesh, const Arcane::IndexedNodeDoFConnectivityView& node_dof, const Arcane::VariableNodeReal3& node_coord, Arcane::VariableDoFReal& rhs_values)
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

    static inline void applyManufacturedSourceToRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_source, IMesh* mesh, const Arcane::IndexedNodeDoFConnectivityView& node_dof, const Arcane::VariableNodeReal3& node_coord, Arcane::VariableDoFReal& rhs_values)
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

    static inline void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const Arcane::IndexedNodeDoFConnectivityView& node_dof, const Arcane::VariableNodeReal3& node_coord, Arcane::VariableDoFReal& rhs_values)
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
    static inline void applyDirichletToLhsAndRhs(BC::IDirichletBoundaryCondition* bs, const Arcane::IndexedNodeDoFConnectivityView& node_dof, const Arcane::VariableNodeReal3& node_coord, FemUtils::DoFLinearSystem& m_linear_system, Arcane::VariableDoFReal& rhs_values)
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
    static inline void applyPointDirichletToLhsAndRhs(BC::IDirichletPointCondition* bs, const Arcane::IndexedNodeDoFConnectivityView& node_dof, const Arcane::VariableNodeReal3& node_coord, FemUtils::DoFLinearSystem& m_linear_system, Arcane::VariableDoFReal& rhs_values)
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
    static inline void applyManufacturedDirichletToLhsAndRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_dirichlet, const Arcane::Real& lambda, const Arcane::FaceGroup& group, BC::IManufacturedSolution* bs, const Arcane::IndexedNodeDoFConnectivityView& node_dof, const Arcane::VariableNodeReal3& node_coord, FemUtils::DoFLinearSystem& m_linear_system, Arcane::VariableDoFReal& rhs_values)
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
};

#endif // ARCANE_FEM_FUNCTIONS_H
