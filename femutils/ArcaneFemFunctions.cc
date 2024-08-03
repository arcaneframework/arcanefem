#include "ArcaneFemFunctions.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the area of a triangle defined by three nodes.
 *
 * This method calculates the area using the determinant formula for a triangle.
 * The area is computed as half the value of the determinant of the matrix
 * formed by the coordinates of the triangle's vertices.
 */
/*---------------------------------------------------------------------------*/

Real ArcaneFemFunctions::computeAreaTriangle3(Cell cell, VariableNodeReal3& node_coord)
{
  Real3 m0 = node_coord[cell.nodeId(0)];
  Real3 m1 = node_coord[cell.nodeId(1)];
  Real3 m2 = node_coord[cell.nodeId(2)];

  return 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the length of the edge defined by a given face.
 *
 * This method calculates Euclidean distance between the two nodes of the face.
 */
/*---------------------------------------------------------------------------*/

Real ArcaneFemFunctions::computeEdgeLength2(Face face, VariableNodeReal3& node_coord)
{
  Real3 m0 = node_coord[face.nodeId(0)];
  Real3 m1 = node_coord[face.nodeId(1)];

  Real dx = m1.x - m0.x;
  Real dy = m1.y - m0.y;

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

Real2 ArcaneFemFunctions::computeEdgeNormal2(Face face, VariableNodeReal3& node_coord)
{
  Real3 m0 = node_coord[face.nodeId(0)];
  Real3 m1 = node_coord[face.nodeId(1)];

  if (!face.isSubDomainBoundaryOutside())
    std::swap(m0, m1);

  Real dx = m1.x - m0.x;
  Real dy = m1.y - m0.y;
  Real norm_N = math::sqrt(dx * dx + dy * dy);

  return { dy / norm_N, -dx / norm_N };
}
