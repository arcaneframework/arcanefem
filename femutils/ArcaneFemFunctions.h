#ifndef ARCANE_FEM_FUNCTIONS_H
#define ARCANE_FEM_FUNCTIONS_H

#include <arcane/VariableTypes.h>

using namespace Arcane;

class ArcaneFemFunctions
{
 public:

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the area of a triangle defined by three nodes.
   *
   * This method calculates the area using the determinant formula for a triangle.
   * The area is computed as half the value of the determinant of the matrix
   * formed by the coordinates of the triangle's vertices.
   */
  /*---------------------------------------------------------------------------*/

  static inline Real computeAreaTriangle3(Cell cell, VariableNodeReal3& node_coord)
  {
    Real3 vertex0 = node_coord[cell.nodeId(0)];
    Real3 vertex1 = node_coord[cell.nodeId(1)];
    Real3 vertex2 = node_coord[cell.nodeId(2)];

    return 0.5 * ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the length of the edge defined by a given face.
   *
   * This method calculates Euclidean distance between the two nodes of the face.
   */
  /*---------------------------------------------------------------------------*/

  static inline Real computeEdgeLength2(Face face, VariableNodeReal3& node_coord)
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

  static inline Real2 computeEdgeNormal2(Face face, VariableNodeReal3& node_coord)
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

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the gradients of basis functions N for P1 triangles.
   *
   * This method calculates gradient operator for a given P1 cell
   *
   * @note we are using 3x3 matrix for this even in 2D hence the z gradient is 0.
   */
  /*---------------------------------------------------------------------------*/

  static inline Real3x3 computeGradientTria3(Cell cell, VariableNodeReal3& node_coord)
  {
    Real3 vertex0 = node_coord[cell.nodeId(0)];
    Real3 vertex1 = node_coord[cell.nodeId(1)];
    Real3 vertex2 = node_coord[cell.nodeId(2)];

    return Real3x3(Real3(vertex1.y - vertex2.y, vertex2.x - vertex1.x, 0),
                   Real3(vertex2.y - vertex0.y, vertex0.x - vertex2.x, 0),
                   Real3(vertex0.y - vertex1.y, vertex1.x - vertex0.x, 0));
  }
};

#endif // ARCANE_FEM_FUNCTIONS_H
