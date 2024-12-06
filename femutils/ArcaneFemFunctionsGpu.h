#ifndef ARCANE_FEM_FUNCTIONS_GPU_H
#define ARCANE_FEM_FUNCTIONS_GPU_H

#include "FemUtils.h"
#include <arcane/accelerator/VariableViews.h>
#include <arcane/core/IndexedItemConnectivityView.h>
#include <arcane/core/Item.h>
#include <arcane/core/MathUtils.h>
#include <arcane/utils/ArcaneGlobal.h>
#include <arccore/base/ArccoreGlobal.h>

namespace Arcane::FemUtils::Gpu::MeshOperation
{

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

} // namespace Arcane::FemUtils::Gpu::MeshOperation

#endif // ! ARCANE_FEM_FUNCTIONS_GPU_H
