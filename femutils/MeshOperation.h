// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* MeshOperation.h                                             (C) 2000-2026 */
/*                                                                           */
/* Various geometric operations on mesh items.                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMTUTILS_MESHOPERATION_H
#define ARCANEFEM_FEMTUTILS_MESHOPERATION_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemUtilsGlobal.h"

#include <arcane/core/VariableTypes.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace ArcaneFemFunctions
{
using namespace Arcane;
using namespace Arcane::FemUtils;

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
   * @brief Computes the area of a polygon (strictly 2D)
   *
   * This method calculates the area of a polygon using the shoelace formula.
   */
  /*---------------------------------------------------------------------------*/

  static inline Real computeAreaPolygon2D(ItemWithNodes item, const VariableNodeReal3& node_coord)
  {
    Int8 n = item.nbNode();
    Real area = 0.0;

    for (Int8 i = 0; i < n; ++i) {
      Int8 j = (i + 1) % n;

      const Real3& pi = node_coord[item.nodeId(i)];
      const Real3& pj = node_coord[item.nodeId(j)];

      area += pi.x * pj.y - pj.x * pi.y;
    }

    return 0.5 * math::abs(area);
  }

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

  static inline Real computeVolumeTetra4(ItemWithNodes item, const VariableNodeReal3& node_coord)
  {
    Real3 vertex0 = node_coord[item.nodeId(0)];
    Real3 vertex1 = node_coord[item.nodeId(1)];
    Real3 vertex2 = node_coord[item.nodeId(2)];
    Real3 vertex3 = node_coord[item.nodeId(3)];

    Real3 v0 = vertex1 - vertex0;
    Real3 v1 = vertex2 - vertex0;
    Real3 v2 = vertex3 - vertex0;

    return math::abs(math::dot(v0, math::cross(v1, v2))) / 6.0;
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

  static inline Real computeAreaTria3(ItemWithNodes item, const VariableNodeReal3& node_coord)
  {
    Real3 n0 = node_coord[item.nodeId(0)];
    Real3 n1 = node_coord[item.nodeId(1)];
    Real3 n2 = node_coord[item.nodeId(2)];

    auto v = math::cross(n1 - n0, n2 - n0);

    return v.normL2() / 2.0;
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

  static inline Real computeAreaQuad4(ItemWithNodes item, const VariableNodeReal3& node_coord)
  {
    Real3 n0 = node_coord[item.nodeId(0)];
    Real3 n1 = node_coord[item.nodeId(1)];
    Real3 n2 = node_coord[item.nodeId(2)];
    Real3 n3 = node_coord[item.nodeId(3)];

    auto tri1x2 = math::cross(n2 - n1, n0 - n1);
    auto tri2x2 = math::cross(n0 - n3, n2 - n3);

    return 0.5 * (tri1x2.normL2() + tri2x2.normL2());
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the volume of a hexahedron defined by eight nodes.
   */
  static inline Real computeVolumeHexa8(ItemWithNodes item, const VariableNodeReal3& node_coord)
  {
    Real3 n0 = node_coord[item.nodeId(0)];
    Real3 n1 = node_coord[item.nodeId(1)];
    Real3 n2 = node_coord[item.nodeId(2)];
    Real3 n3 = node_coord[item.nodeId(3)];
    Real3 n4 = node_coord[item.nodeId(4)];
    Real3 n5 = node_coord[item.nodeId(5)];
    Real3 n6 = node_coord[item.nodeId(6)];
    Real3 n7 = node_coord[item.nodeId(7)];

    Real v1 = math::matDet((n6 - n1) + (n7 - n0), n6 - n3, n2 - n0);
    Real v2 = math::matDet(n7 - n0, (n6 - n3) + (n5 - n0), n6 - n4);
    Real v3 = math::matDet(n6 - n1, n5 - n0, (n6 - n4) + (n2 - n0));
    return (v1 + v2 + v3) / 12.;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the volume of a pentahedron (wedge or triangular prism)
   * defined by six nodes.
   */
  static inline Real penta6Volume(ItemWithNodes item, const VariableNodeReal3& node_coord)
  {
    Real3 n0 = node_coord[item.nodeId(0)];
    Real3 n1 = node_coord[item.nodeId(1)];
    Real3 n2 = node_coord[item.nodeId(2)];
    Real3 n3 = node_coord[item.nodeId(3)];
    Real3 n4 = node_coord[item.nodeId(4)];
    Real3 n5 = node_coord[item.nodeId(5)];

    auto v = math::cross(n1 - n0, n2 - n0);
    auto base = 0.5 * v.normL2();
    auto h1 = (n3 - n0).normL2();
    auto h2 = (n4 - n1).normL2();
    auto h3 = (n5 - n2).normL2();

    return base * (h1 + h2 + h3) / 3.0;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the volume of a pyramid defined by five nodes.
   */
  static inline Real pyramid5Volume(ItemWithNodes item, const VariableNodeReal3& node_coord)
  {
    Real3 n0 = node_coord[item.nodeId(0)];
    Real3 n1 = node_coord[item.nodeId(1)];
    Real3 n2 = node_coord[item.nodeId(2)];
    Real3 n3 = node_coord[item.nodeId(3)];
    Real3 n4 = node_coord[item.nodeId(4)];

    // Compute the area of the base triangle
    auto v = math::cross(n1 - n0, n2 - n0);

    auto tri1x2 = math::cross(n2 - n1, n0 - n1);
    auto tri2x2 = math::cross(n0 - n3, n2 - n3);
    auto base = 0.5 * (tri1x2.normL2() + tri2x2.normL2());

    // Compute the height of the pyramid
    // The height is the distance from the apex (n4) to the base plane
    // formed by the triangle (n0, n1, n2)
    auto h = math::dot(n4 - n0, v) / v.normL2();

    return (base * h) / 3.0;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the centroid of an item with nodes.
   *
   * This method calculates the centroid (geometric center) of an item by
   * averaging the coordinates of its nodes. The centroid is computed as the
   * mean of the vertex coordinates, providing a central point that represents
   * the item's position in space.
   */
  static inline Real3 computeCentroid(ItemWithNodes item, const VariableNodeReal3& node_coord)
  {
    Int8 nb_node = item.nbNode();
    Real3 centroid = { 0., 0., 0. };

    for (Int8 i = 0; i < nb_node; ++i) {
      Real3 vertex = node_coord[item.nodeId(i)];
      centroid.x += vertex.x;
      centroid.y += vertex.y;
      centroid.z += vertex.z;
    }

    centroid.x /= nb_node;
    centroid.y /= nb_node;
    centroid.z /= nb_node;

    return centroid;
  }

  /*---------------------------------------------------------------------------*/
  /** @brief Returns the DG penalty length 2|K|/|F| in two dimensions.
   *
   * This method computes the Discontinuous Galerkin (DG) penalty length for a
   * 2D cell and face. The penalty length is calculated  as  twice the area of
   * the cell divided by the length of the face.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real computeDGPenaltyLength2D(Cell cell, Face face, const VariableNodeReal3& node_coord)
  {
    return 2.0 * computeAreaPolygon2D(cell, node_coord) / computeLengthEdge2(face, node_coord);
  }

  /*---------------------------------------------------------------------------*/
  /** @brief Computes the area of a planar polygon embedded in three dimensions.
   *
   * This method calculates the area of a polygon in 3D space by decomposing it
   * into triangles. Area is  computed as the sum of the areas of the triangles
   * formed by the polygon's vertices and an anchor point.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real computeAreaPolygon3D(Face face, const VariableNodeReal3& node_coord)
  {
    Real area = 0.0;
    Real3 anchor = node_coord[face.nodeId(0)];
    for (Int32 i = 1; i + 1 < face.nbNode(); ++i) {
      Real3 b = node_coord[face.nodeId(i)];
      Real3 c = node_coord[face.nodeId(i + 1)];
      area += 0.5 * math::cross(b - anchor, c - anchor).normL2();
    }
    return area;
  }

  /*---------------------------------------------------------------------------*/
  /** @brief Computes convex polyhedron volume by a face-tetrahedra decomposition.
   *
   * This method calculates the volume of a convex polyhedron by decomposing it
   * into tetrahedra. The volume  is  computed  as the sum of the volumes of the
   * tetrahedra formed by the polyhedron's faces and a center point.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real computeVolumePolyhedron3D(Cell cell, const VariableNodeReal3& node_coord)
  {
    Real volume = 0.0;
    Real3 center = computeCentroid(cell, node_coord);
    for (Face face : cell.faces()) {
      Real3 anchor = node_coord[face.nodeId(0)];
      for (Int32 i = 1; i + 1 < face.nbNode(); ++i) {
        Real3 b = node_coord[face.nodeId(i)];
        Real3 c = node_coord[face.nodeId(i + 1)];
        volume += math::abs(math::dot(anchor - center, math::cross(b - center, c - center))) / 6.0;
      }
    }
    return volume;
  }

  /*---------------------------------------------------------------------------*/
  /** @brief Computes cell diameter used as the DG penalty length in three dimensions.
   *
   * This method computes the Discontinuous Galerkin (DG) penalty length for a
   * 3D cell. The penalty length is calculated as the maximum distance between
   * any two nodes of the cell.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real computeDGPenaltyLength3D(Cell cell, Face, const VariableNodeReal3& node_coord)
  {
    Real diameter = 0.0;
    for (Int32 i = 0; i < cell.nbNode(); ++i) {
      Real3 a = node_coord[cell.nodeId(i)];
      for (Int32 j = i + 1; j < cell.nbNode(); ++j)
        diameter = math::max(diameter, (node_coord[cell.nodeId(j)] - a).normL2());
    }
    return diameter;
  }

  /*---------------------------------------------------------------------------*/
  /** @brief Computes the unit face normal directed out of a 2D cell.
   *
   * This method calculates the unit normal vector of a face in 2D space, ensuring
   * that it is directed outward from the  specified  cell. The  normal  vector is
   * computed based on  the  edge  defined by the  face's  nodes and is normalized
   * to have a length of one. The orientation is adjusted to ensure it points away
   * from the cell's centroid.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 computeUnitNormal2D(Face face, Cell cell, const VariableNodeReal3& node_coord)
  {
    Real3 edge = node_coord[face.nodeId(1)] - node_coord[face.nodeId(0)];
    Real3 normal = { edge.y, -edge.x, 0.0 };
    Real norm = normal.normL2();
#ifdef _DEBUG
    if (norm <= 1.0e-30)
      ARCANE_FATAL("Degenerate face {0} has a zero normal", face.uniqueId());
#endif
    normal = normal / norm;
    if (math::dot(normal, computeCentroid(face, node_coord) - computeCentroid(cell, node_coord)) < 0.0)
      normal = -normal;
    return normal;
  }

  /*---------------------------------------------------------------------------*/
  /** @brief Computes the Newell unit normal directed out of a 3D cell.
   *
   * This method calculates the unit normal vector of a face in 3D space using
   * Newell's method, ensuring that it is directed outward  from the specified
   * cell. The normal vector is computed based on the vertices of the face and
   * is normalized to have a  length  of one. The  orientation  is adjusted to
   * ensure it points away from the cell's centroid.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 computeUnitNormal3D(Face face, Cell cell, const VariableNodeReal3& node_coord)
  {
    Real3 normal = { 0.0, 0.0, 0.0 };
    for (Int32 i = 0; i < face.nbNode(); ++i) {
      Real3 p = node_coord[face.nodeId(i)];
      Real3 q = node_coord[face.nodeId((i + 1) % face.nbNode())];
      normal.x += (p.y - q.y) * (p.z + q.z);
      normal.y += (p.z - q.z) * (p.x + q.x);
      normal.z += (p.x - q.x) * (p.y + q.y);
    }
    Real norm = normal.normL2();
#ifdef _DEBUG
    if (norm <= 1.0e-30)
      ARCANE_FATAL("Degenerate face {0} has a zero normal", face.uniqueId());
#endif
    normal = normal / norm;
    if (math::dot(normal, computeCentroid(face, node_coord) - computeCentroid(cell, node_coord)) < 0.0)
      normal = -normal;
    return normal;
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
   * @brief Computes the barycenter (centroid) of a tetrahedron.
   *
   * This method calculates the barycenter of a tetrahedron defined by its nodes.
   * The barycenter is computed as the average of the vertices' coordinates.
   */
  /*---------------------------------------------------------------------------*/

  static inline Real3 computeBaryCenterTetra4(Cell cell, const VariableNodeReal3& node_coord)
  {
    Real3 vertex0 = node_coord[cell.nodeId(0)];
    Real3 vertex1 = node_coord[cell.nodeId(1)];
    Real3 vertex2 = node_coord[cell.nodeId(2)];
    Real3 vertex3 = node_coord[cell.nodeId(3)];

    Real Center_x = (vertex0.x + vertex1.x + vertex2.x + vertex3.x) / 4.;
    Real Center_y = (vertex0.y + vertex1.y + vertex2.y + vertex3.x) / 4.;
    Real Center_z = (vertex0.z + vertex1.z + vertex2.z + vertex3.z) / 4.;

    return { Center_x, Center_y, Center_z };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the length of the edge defined by a given face.
   *
   * This method calculates Euclidean distance between the two nodes of the face.
   */
  /*---------------------------------------------------------------------------*/

  static inline Real computeLengthEdge2(ItemWithNodes item, const VariableNodeReal3& node_coord)
  {
    Real3 vertex0 = node_coord[item.nodeId(0)];
    Real3 vertex1 = node_coord[item.nodeId(1)];

    return (vertex1 - vertex0).normL2();
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

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the factor used for integration of 1D, 2D and 3D finite elements
   *
   * This method calculates length or surface for a given finite-element (P1, P2, ...) and returns
   * the associated factor for elementary integrals.
   */
  /*---------------------------------------------------------------------------*/
  static inline Real computeFacLengthOrArea(Face face, const VariableNodeReal3& node_coord)
  {
    Int32 item_type = face.type();
    Real fac_el{ 0. };

    switch (item_type) {

    // Lines
    case IT_Line2:
    case IT_Line3:
      fac_el = computeLengthEdge2(face, node_coord) / 2.;
      break;

    // Faces
    case IT_Triangle3:
    case IT_Triangle6:
      fac_el = computeAreaTria3(face, node_coord) / 3.;
      break;

    case IT_Quad4:
    case IT_Quad8:
      fac_el = computeAreaQuad4(face, node_coord) / 4.;
      break;

    default:
      break;
    }
    return fac_el;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the normalized triangle normal for a given face.
   *
   * This method calculates normal vector to the triangle defined by nodes,
   * of the face and normalizes it, and ensures the correct orientation.
   */
  /*---------------------------------------------------------------------------*/

  static inline Real3 computeNormalTriangle(Face face, const VariableNodeReal3& node_coord)
  {
    // Get the vertices of the triangle
    Real3 vertex0 = node_coord[face.nodeId(0)];
    Real3 vertex1 = node_coord[face.nodeId(1)];
    Real3 vertex2 = node_coord[face.nodeId(2)];

    if (!face.isSubDomainBoundaryOutside())
      std::swap(vertex0, vertex1);

    // Calculate two edges of the triangle
    Real3 edge1 = { vertex1.x - vertex0.x, vertex1.y - vertex0.y, vertex1.z - vertex0.z };
    Real3 edge2 = { vertex2.x - vertex0.x, vertex2.y - vertex0.y, vertex2.z - vertex0.z };

    // Compute the cross product of the two edges
    Real3 normal = {
      edge1.y * edge2.z - edge1.z * edge2.y,
      edge1.z * edge2.x - edge1.x * edge2.z,
      edge1.x * edge2.y - edge1.y * edge2.x
    };

    // Calculate the magnitude of the normal vector
    Real norm = math::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

    // Normalize the vector to unit length and return
    return { normal.x / norm, normal.y / norm, normal.z / norm };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes the normalized quad normal for a given face.
   *
   * This method calculates normal vector to the quad defined by nodes,
   * of the face and normalizes it, and ensures the correct orientation.
   * Newell's Method is used to compute the normal vector for the quadrilateral.
   */
  /*---------------------------------------------------------------------------*/

  static inline Real3 computeNormalQuad(Face face, const VariableNodeReal3& node_coord)
  {
    Real3 normal = { 0.0, 0.0, 0.0 };

    Int32 n = face.nbNode(); // Should be 4 for quad

    for (Int32 i = 0; i < n; ++i) {
      Real3 v_curr = node_coord[face.nodeId(i)];
      Real3 v_next = node_coord[face.nodeId((i + 1) % n)];

      normal.x += (v_curr.y - v_next.y) * (v_curr.z + v_next.z);
      normal.y += (v_curr.z - v_next.z) * (v_curr.x + v_next.x);
      normal.z += (v_curr.x - v_next.x) * (v_curr.y + v_next.y);
    }

    // Normalize
    Real norm = math::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    return { normal.x / norm, normal.y / norm, normal.z / norm };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes finite-element entity (Edge, Face or Cell) geometric dimension
   * This method is used for the FEM 2D & 3D needs (coming from PASSMO)
   * for Jacobian & elementary matrices computations
   */
  /*---------------------------------------------------------------------------*/
  static inline Int32 getGeomDimension(ItemWithNodes item)
  {
    return item.typeInfo()->dimension();
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Computes edge & Face normal & tangent vectors (normalized, direct oriented)
   * This method is used for the FEM 2D & 3D needs (coming from PASSMO)
   * for geometric transformations (rotations, projections, ...)
   * In 2D, it assumes the edge lies in x-y plane (z coord = 0.)
   */
  /*---------------------------------------------------------------------------*/
  static inline void dirVectors(Face face, const VariableNodeReal3& node_coord, Integer ndim, Real3& e1, Real3& e2, Real3& e3)
  {
    Real3 n0 = node_coord[face.nodeId(0)];
    Real3 n1 = node_coord[face.nodeId(1)];

    if (!face.isSubDomainBoundaryOutside())
      std::swap(n0, n1);

    // 1st in-plane vector/along edge
    e1 = n1 - n0;

    if (ndim == 3) {

      Real3 n2 = node_coord[face.nodeId(2)];

      // out Normal to the face plane
      e3 = math::cross(e1, n2 - n0);

      // 2nd in-plane vector
      e2 = math::cross(e3, e1);
      e3 = math::mutableNormalize(e3);
    }
    else {

      Cell cell{ face.boundaryCell() };
      Node nod;
      for (Node node : cell.nodes()) {
        if (node != face.node(0) && node != face.node(1)) {
          nod = node;
          break;
        }
      }

      // Out Normal to the edge
      e2 = { -e1.y, e1.x, 0. };

      Real3 n2 = node_coord[nod];
      auto sgn = math::dot(e2, n2 - n0);
      if (sgn > 0.)
        e2 *= -1.;
    }
    e1 = math::mutableNormalize(e1);
    e2 = math::mutableNormalize(e2);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace ArcaneFemFunctions

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
