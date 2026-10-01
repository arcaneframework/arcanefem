// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DgQuadrature.cc                                             (C) 2000-2026 */
/*                                                                           */
/* Provides quadrature rules used by discontinuous Galerkin methods.         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "DgQuadrature.h"

#include "MeshOperation.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace ArcaneFemFunctions
{

UniqueArray<DgQuadraturePoint> DgQuadrature::
computeCellQuadrature2D(Cell cell, const VariableNodeReal3& node_coord)
{
  return { { MeshOperation::computeCentroid(cell, node_coord),
             MeshOperation::computeAreaPolygon2D(cell, node_coord) } };
}

/*---------------------------------------------------------------------------*/

UniqueArray<DgQuadraturePoint> DgQuadrature::
computeFaceQuadrature2D(Face face, const VariableNodeReal3& node_coord)
{
  return { { MeshOperation::computeCentroid(face, node_coord),
             MeshOperation::computeLengthEdge2(face, node_coord) } };
}

/*---------------------------------------------------------------------------*/

UniqueArray<DgQuadraturePoint> DgQuadrature::
computeFaceQuadrature3D(Face face, const VariableNodeReal3& node_coord)
{
  UniqueArray<DgQuadraturePoint> quadrature;
  Real3 center = MeshOperation::computeCentroid(face, node_coord);
  for (Int32 i = 0; i < face.nbNode(); ++i) {
    Real3 b = node_coord[face.nodeId(i)];
    Real3 c = node_coord[face.nodeId((i + 1) % face.nbNode())];
    Real triangle_area = 0.5 * math::cross(b - center, c - center).normL2();
    // Skip degenerate triangles
    if (triangle_area <= 1.0e-30)
      continue;
    Real weight = triangle_area / 3.0;
    quadrature.add({ { (4.0 * center.x + b.x + c.x) / 6.0,
                       (4.0 * center.y + b.y + c.y) / 6.0,
                       (4.0 * center.z + b.z + c.z) / 6.0 },
                     weight });
    quadrature.add({ { (center.x + 4.0 * b.x + c.x) / 6.0,
                       (center.y + 4.0 * b.y + c.y) / 6.0,
                       (center.z + 4.0 * b.z + c.z) / 6.0 },
                     weight });
    quadrature.add({ { (center.x + b.x + 4.0 * c.x) / 6.0,
                       (center.y + b.y + 4.0 * c.y) / 6.0,
                       (center.z + b.z + 4.0 * c.z) / 6.0 },
                     weight });
  }
  return quadrature;
}

/*---------------------------------------------------------------------------*/

UniqueArray<DgQuadraturePoint> DgQuadrature::
computeCellQuadrature3D(Cell cell, const VariableNodeReal3& node_coord)
{
  UniqueArray<DgQuadraturePoint> quadrature;
  Real3 center = MeshOperation::computeCentroid(cell, node_coord);
  for (Face face : cell.faces()) {
    Real3 anchor = node_coord[face.nodeId(0)];
    for (Int32 i = 1; i + 1 < face.nbNode(); ++i) {
      Real3 b = node_coord[face.nodeId(i)];
      Real3 c = node_coord[face.nodeId(i + 1)];
      Real volume = math::abs(math::dot(anchor - center,
                                        math::cross(b - center, c - center))) /
      6.0;
      // Skip degenerate tetrahedra
      if (volume <= 1.0e-30)
        continue;

      Real3 vertices[4] = { center, anchor, b, c };
      Real3 tetra_center = (center + anchor + b + c) / 4.0;
      quadrature.add({ tetra_center, -4.0 * volume / 5.0 });
      for (Int32 vertex = 0; vertex < 4; ++vertex) {
        Real3 point = { 0.0, 0.0, 0.0 };
        for (Int32 j = 0; j < 4; ++j)
          point += vertices[j] * ((j == vertex) ? 0.5 : 1.0 / 6.0);
        quadrature.add({ point, 9.0 * volume / 20.0 });
      }
    }
  }
  return quadrature;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace ArcaneFemFunctions
