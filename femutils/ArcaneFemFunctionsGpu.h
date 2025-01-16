// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArcaneFemFunctionsGpu.h                                     (C) 2022-2025 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef ARCANE_FEM_FUNCTIONS_GPU_H
#define ARCANE_FEM_FUNCTIONS_GPU_H

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/accelerator/VariableViews.h>
#include <arccore/base/ArccoreGlobal.h>
#include <arcane/utils/ArcaneGlobal.h>

#include <arcane/core/IndexedItemConnectivityView.h>
#include <arcane/core/MathUtils.h>
#include <arcane/core/Item.h>

#include "FemUtils.h"
#include "IArcaneFemBC.h"

namespace Arcane::FemUtils::Gpu::MeshOperation
{

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the area of a triangle defined by three nodes.
 *
 * This method calculates the area using the determinant formula for a triangle.
 * The area is computed as half the value of the determinant of the matrix
 * formed by the coordinates of the triangle's vertices.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static inline Real computeAreaTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 vertex0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 vertex1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 vertex2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];

  return 0.5 * ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));
};

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

ARCCORE_HOST_DEVICE static inline Real computeVolumeTetra4(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 vertex0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 vertex1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 vertex2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real3 vertex3 = in_node_coord[cn_cv.nodeId(cell_lid, 3)];

  Real3 v0 = vertex1 - vertex0;
  Real3 v1 = vertex2 - vertex0;
  Real3 v2 = vertex3 - vertex0;

  return math::abs(math::dot(v0, math::cross(v1, v2))) / 6.0;
};

} // namespace Arcane::FemUtils::Gpu::MeshOperation

namespace Arcane::FemUtils::Gpu::FeOperation2D
{

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
 *         ∂N/∂x = 1/(2A) [ y2-y3  y3−y1  y1−y2 ]
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static inline Real3 computeGradientXTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 vertex0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 vertex1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 vertex2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];

  auto A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

  return { (vertex1.y - vertex2.y) / A2, (vertex2.y - vertex0.y) / A2, (vertex0.y - vertex1.y) / A2 };
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
 *         ∂N/∂x = 1/(2A) [ x3−x2  x1−x3  x2−x1 ]
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static inline Real3 computeGradientYTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 vertex0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 vertex1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 vertex2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];

  auto A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

  return { (vertex2.x - vertex1.x) / A2, (vertex0.x - vertex2.x) / A2, (vertex1.x - vertex0.x) / A2 };
}

} // namespace Arcane::FemUtils::Gpu::FeOperation2D

namespace Arcane::FemUtils::Gpu::FeOperation3D
{

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

ARCCORE_HOST_DEVICE static inline Real4 computeGradientXTetra4(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 vertex0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 vertex1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 vertex2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real3 vertex3 = in_node_coord[cn_cv.nodeId(cell_lid, 3)];

  Real3 v0 = vertex1 - vertex0;
  Real3 v1 = vertex2 - vertex0;
  Real3 v2 = vertex3 - vertex0;

  // 6 x Volume of tetrahedron
  Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

  Real4 dx{};

  dx[0] = (vertex1.y * (vertex3.z - vertex2.z) + vertex2.y * (vertex1.z - vertex3.z) + vertex3.y * (vertex2.z - vertex1.z)) / V6;
  dx[1] = (vertex0.y * (vertex2.z - vertex3.z) + vertex2.y * (vertex3.z - vertex0.z) + vertex3.y * (vertex0.z - vertex2.z)) / V6;
  dx[2] = (vertex0.y * (vertex3.z - vertex1.z) + vertex1.y * (vertex0.z - vertex3.z) + vertex3.y * (vertex1.z - vertex0.z)) / V6;
  dx[3] = (vertex0.y * (vertex1.z - vertex2.z) + vertex1.y * (vertex2.z - vertex0.z) + vertex2.y * (vertex0.z - vertex1.z)) / V6;

  return dx;
};

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

ARCCORE_HOST_DEVICE static inline Real4 computeGradientYTetra4(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 vertex0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 vertex1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 vertex2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real3 vertex3 = in_node_coord[cn_cv.nodeId(cell_lid, 3)];

  Real3 v0 = vertex1 - vertex0;
  Real3 v1 = vertex2 - vertex0;
  Real3 v2 = vertex3 - vertex0;

  // 6 x Volume of tetrahedron
  Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

  Real4 dy{};

  dy[0] = (vertex1.z * (vertex3.x - vertex2.x) + vertex2.z * (vertex1.x - vertex3.x) + vertex3.z * (vertex2.x - vertex1.x)) / V6;
  dy[1] = (vertex0.z * (vertex2.x - vertex3.x) + vertex2.z * (vertex3.x - vertex0.x) + vertex3.z * (vertex0.x - vertex2.x)) / V6;
  dy[2] = (vertex0.z * (vertex3.x - vertex1.x) + vertex1.z * (vertex0.x - vertex3.x) + vertex3.z * (vertex1.x - vertex0.x)) / V6;
  dy[3] = (vertex0.z * (vertex1.x - vertex2.x) + vertex1.z * (vertex2.x - vertex0.x) + vertex2.z * (vertex0.x - vertex1.x)) / V6;

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

ARCCORE_HOST_DEVICE static inline Real4 computeGradientZTetra4(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 vertex0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 vertex1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 vertex2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real3 vertex3 = in_node_coord[cn_cv.nodeId(cell_lid, 3)];

  auto v0 = vertex1 - vertex0;
  auto v1 = vertex2 - vertex0;
  auto v2 = vertex3 - vertex0;

  // 6 x Volume of tetrahedron
  Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

  Real4 dz{};

  dz[0] = (vertex1.x * (vertex3.y - vertex2.y) + vertex2.x * (vertex1.y - vertex3.y) + vertex3.x * (vertex2.y - vertex1.y)) / V6;
  dz[1] = (vertex0.x * (vertex2.y - vertex3.y) + vertex2.x * (vertex3.y - vertex0.y) + vertex3.x * (vertex0.y - vertex2.y)) / V6;
  dz[2] = (vertex0.x * (vertex3.y - vertex1.y) + vertex1.x * (vertex0.y - vertex3.y) + vertex3.x * (vertex1.y - vertex0.y)) / V6;
  dz[3] = (vertex0.x * (vertex1.y - vertex2.y) + vertex1.x * (vertex2.y - vertex0.y) + vertex2.x * (vertex0.y - vertex1.y)) / V6;

  return dz;
};

} // namespace Arcane::FemUtils::Gpu::FeOperation3D

namespace Arcane::FemUtils::Gpu::BoundaryConditions2D
{

static inline void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values) {}

} // namespace Arcane::FemUtils::Gpu::BoundaryConditions2D

#endif // ! ARCANE_FEM_FUNCTIONS_GPU_H
