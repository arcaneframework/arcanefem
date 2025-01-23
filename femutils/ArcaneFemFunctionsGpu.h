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

#include <arcane/VariableTypes.h>
#include <arcane/IMesh.h>

#include <arcane/accelerator/RunCommandEnumerate.h>
#include <arcane/accelerator/AcceleratorGlobal.h>
#include <arcane/accelerator/NumArrayViews.h>
#include <arcane/accelerator/VariableViews.h>
#include <arcane/accelerator/RunCommand.h>
#include <arcane/accelerator/RunQueue.h>
#include <arcane/accelerator/Atomic.h>

#include <arccore/base/ArccoreGlobal.h>

#include <arcane/core/UnstructuredMeshConnectivity.h>
#include <arcane/core/IndexedItemConnectivityView.h>
#include <arcane/core/ItemInfoListView.h>
#include <arcane/core/ItemEnumerator.h>
#include <arcane/core/MathUtils.h>
#include <arcane/core/ItemTypes.h>
#include <arcane/core/DataView.h>
#include <arcane/core/Item.h>

#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/ArrayLayout.h>
#include "arcane/core/VariableTypedef.h"
#include "arcane/utils/UtilsTypes.h"
#include <arcane/utils/NumArray.h>

#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "IArcaneFemBC.h"
#include "FemUtils.h"

namespace Arcane::FemUtils::Gpu::Csr
{
ARCCORE_HOST_DEVICE static Int32 findIndex(Int32 begin, Int32 end, Int32 column_lid, Span<const Int32> in_csr_columns)
{
  for (auto i = begin; i < end; ++i)
    if (in_csr_columns[i] == column_lid)
      return i;
  return -1;
}
} // namespace Arcane::FemUtils::Gpu::Csr

namespace Arcane::FemUtils::Gpu::MeshOperation
{

/*---------------------------------------------------------------------------*/
/**
 * @brief computes the area of a triangle defined by three nodes.
 *
 * this method calculates the area using the determinant formula for a triangle.
 * the area is computed as half the value of the determinant of the matrix
 * formed by the coordinates of the triangle's vertices.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static inline Real computeAreaTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 vertex0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 vertex1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 vertex2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];

  auto v = math::cross(vertex1 - vertex0, vertex2 - vertex0);
  return v.normL2() / 2.0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static inline Real computeAreaTria(FaceLocalId face_lid, const IndexedFaceNodeConnectivityView& fn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 vertex0 = in_node_coord[fn_cv.nodeId(face_lid, 0)];
  Real3 vertex1 = in_node_coord[fn_cv.nodeId(face_lid, 1)];
  Real3 vertex2 = in_node_coord[fn_cv.nodeId(face_lid, 2)];

  auto v = math::cross(vertex1 - vertex0, vertex2 - vertex0);
  return v.normL2() / 2.0;
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

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the length of the edge defined by a given face.
 *
 * This method calculates Euclidean distance between the two nodes of the face.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static inline Real computeLengthFace(FaceLocalId face_lid, IndexedFaceNodeConnectivityView fn_cv, Accelerator::VariableNodeReal3InView in_node_coord)
{

  Real3 m0 = in_node_coord[fn_cv.nodeId(face_lid, 0)];
  Real3 m1 = in_node_coord[fn_cv.nodeId(face_lid, 1)];
  return math::sqrt((m1.x - m0.x) * (m1.x - m0.x) + (m1.y - m0.y) * (m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the normalized edge normal for a given face.
 *
 * This method calculates normal vector to the edge defined by nodes of the face,
 * normalizes it, and ensures the correct orientation.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static inline Real2 computeNormalFace(FaceLocalId face_lid, IndexedFaceNodeConnectivityView fn_cv, Accelerator::VariableNodeReal3InView in_node_coord, FaceInfoListView faces_infos)
{
  Real3 m0 = in_node_coord[fn_cv.nodeId(face_lid, 0)];
  Real3 m1 = in_node_coord[fn_cv.nodeId(face_lid, 1)];

  if (!faces_infos.isSubDomainBoundaryOutside(face_lid)) {
    Real3 tmp = m0;
    m0 = m1;
    m1 = tmp;
  }

  Real2 N;
  Real norm_N = math::sqrt((m1.y - m0.y) * (m1.y - m0.y) + (m1.x - m0.x) * (m1.x - m0.x));
  N.x = (m1.y - m0.y) / norm_N;
  N.y = (m0.x - m1.x) / norm_N;
  return N;
}

/*---------------------------------------------------------------------------*/
/**
     * @brief Computes the normalized triangle normal for a given face.
     *
     * This method calculates normal vector to the triangle defined by nodes,
     * of the face and normalizes it, and ensures the correct orientation.
     */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE static inline Real3 computeNormalTriangle(FaceLocalId face_lid, IndexedFaceNodeConnectivityView fn_cv, Accelerator::VariableNodeReal3InView in_node_coord, FaceInfoListView faces_infos)
{

  Real3 vertex0 = in_node_coord[fn_cv.nodeId(face_lid, 0)];
  Real3 vertex1 = in_node_coord[fn_cv.nodeId(face_lid, 1)];
  Real3 vertex2 = in_node_coord[fn_cv.nodeId(face_lid, 2)];

  if (!faces_infos.isSubDomainBoundaryOutside(face_lid))
    std::swap(vertex0, vertex1);

  Real3 edge1 = { vertex1.x - vertex0.x, vertex1.y - vertex0.y, vertex1.z - vertex0.z };
  Real3 edge2 = { vertex2.x - vertex0.x, vertex2.y - vertex0.y, vertex2.z - vertex0.z };

  Real3 normal = {
    edge1.y * edge2.z - edge1.z * edge2.y,
    edge1.z * edge2.x - edge1.x * edge2.z,
    edge1.x * edge2.y - edge1.y * edge2.x
  };

  Real norm = math::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
  return { normal.x / norm, normal.y / norm, normal.z / norm };
}

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

namespace Arcane::FemUtils::Gpu
{

namespace BoundaryCondtionsHelpers
{
  template <class Function>
  inline void applyConstantSourceToRhsBase(Function computeSpatialIntegral, Real qdot, Accelerator::RunQueue* queue, VariableDoFReal& rhs_variable_na, IMesh* mesh, const FemDoFsOnNodes& dofs_on_nodes, VariableNodeReal3 node_coord)
  {
    ARCANE_CHECK_PTR(queue);
    ARCANE_CHECK_PTR(mesh);
    ARCANE_CHECK_PTR(computeSpatialIntegral);

    UnstructuredMeshConnectivityView connectivity_view;
    connectivity_view.setMesh(mesh);
    NodeInfoListView nodes_infos(mesh->nodeFamily());
    auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
    auto cn_cv = connectivity_view.cellNode();

    auto command = Accelerator::makeCommand(queue);

    auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
    auto in_node_coord = Accelerator::viewIn(command, node_coord);

    command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh->allCells())
    {
      Real domain = computeSpatialIntegral(cell_lid, cn_cv, in_node_coord);
      for (NodeLocalId node_lid : cn_cv.nodes(cell_lid)) {
        if (nodes_infos.isOwn(node_lid)) {
          Real value = qdot * domain / cn_cv.nbItem(cell_lid);
          Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], value);
        }
      }
    };
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  inline void applyDirichletToNodeGroupViaPenalty(Real value, Real penalty, Accelerator::RunQueue* queue, IMesh* mesh, DoFLinearSystem& linear_system, const FemDoFsOnNodes& dofs_on_nodes, NodeGroup& node_group, Int32 dof_index = 0)
  {
    ARCANE_CHECK_PTR(queue);
    ARCANE_CHECK_PTR(mesh);

    NodeInfoListView nodes_infos(mesh->nodeFamily());
    auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());

    auto command = Accelerator::makeCommand(queue);
    auto in_out_forced_info = Accelerator::viewInOut(command, linear_system.getForcedInfo());
    auto in_out_forced_value = Accelerator::viewInOut(command, linear_system.getForcedValue());
    auto in_out_rhs_variable = Accelerator::viewInOut(command, linear_system.rhsVariable());

    command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, node_group)
    {
      if (nodes_infos.isOwn(node_lid)) {
        DoFLocalId dof_id = node_dof.dofId(node_lid, dof_index);
        in_out_forced_info[dof_id] = true;
        in_out_forced_value[dof_id] = penalty;
        in_out_rhs_variable[dof_id] = penalty * value;
      }
    };
  }

} // namespace BoundaryCondtionsHelpers

namespace BoundaryConditions
{

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Applies Dirichlet boundary conditions to RHS and LHS.
   *
   * Updates the LHS matrix and RHS vector to enforce Dirichlet conditions.
   *
   * - For LHS matrix `A`, the diagonal term for the Dirichlet DOF is set to `P`.
   * - For RHS vector `b`, the Dirichlet DOF term is scaled by `P`.
   *
   */
  /*---------------------------------------------------------------------------*/

  inline void applyDirichletViaPenalty(BC::IDirichletBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes, DoFLinearSystem& linear_system, IMesh* mesh, Accelerator::RunQueue* queue)
  {
    ARCANE_CHECK_PTR(bs);
    Real value = bs->getValue();
    Real penalty = bs->getPenalty();
    FaceGroup face_group = bs->getSurface();
    NodeGroup node_group = face_group.nodeGroup();
    BoundaryCondtionsHelpers::applyDirichletToNodeGroupViaPenalty(value, penalty, queue, mesh, linear_system, dofs_on_nodes, node_group);
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  inline void applyDirichletViaPenaltyVectorial(const FemDoFsOnNodes& dofs_on_nodes, DoFLinearSystem& linear_system, IMesh* mesh, Accelerator::RunQueue* queue, FaceGroup face_group, Real penalty, const UniqueArray<String> u_dirichlet_str)
  {
    NodeGroup node_group = face_group.nodeGroup();
    for (auto dof_index = 0; dof_index < u_dirichlet_str.size(); ++dof_index) {
      if (u_dirichlet_str[dof_index] != "NULL") {
        Real u_dirichlet = std::stod(u_dirichlet_str[dof_index].localstr());
        BoundaryCondtionsHelpers::applyDirichletToNodeGroupViaPenalty(u_dirichlet, penalty, queue, mesh, linear_system, dofs_on_nodes, node_group, dof_index);
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
   */
  /*---------------------------------------------------------------------------*/

  static inline void applyPointDirichletViaPenalty(BC::IDirichletPointCondition* bs, const FemDoFsOnNodes& dofs_on_nodes, DoFLinearSystem& linear_system, IMesh* mesh, Accelerator::RunQueue* queue)
  {
    ARCANE_CHECK_PTR(bs);
    Real value = bs->getValue();
    Real penalty = bs->getPenalty();
    NodeGroup node_group = bs->getNode();
    BoundaryCondtionsHelpers::applyDirichletToNodeGroupViaPenalty(value, penalty, queue, mesh, linear_system, dofs_on_nodes, node_group);
  }

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

  static inline void applyPointDirichletViaPenaltyVectorial(const FemDoFsOnNodes& dofs_on_nodes, DoFLinearSystem& linear_system, IMesh* mesh, Accelerator::RunQueue* queue, NodeGroup node_group, Real penalty, const UniqueArray<String> u_dirichlet_str)
  {
    for (auto dof_index = 0; dof_index < u_dirichlet_str.size(); ++dof_index) {
      if (u_dirichlet_str[dof_index] != "NULL") {
        Real u_dirichlet = std::stod(u_dirichlet_str[dof_index].localstr());
        BoundaryCondtionsHelpers::applyDirichletToNodeGroupViaPenalty(u_dirichlet, penalty, queue, mesh, linear_system, dofs_on_nodes, node_group, dof_index);
      }
    }
  }

}; // namespace BoundaryConditions

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
   */
  /*---------------------------------------------------------------------------*/

  static inline void applyConstantSourceToRhs(Real qdot, const FemDoFsOnNodes& dofs_on_nodes, VariableNodeReal3 node_coord, VariableDoFReal& rhs_variable_na, IMesh* mesh, Accelerator::RunQueue* queue)
  {
    BoundaryCondtionsHelpers::applyConstantSourceToRhsBase([] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord) { return MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord); }, qdot, queue, rhs_variable_na, mesh, dofs_on_nodes, node_coord);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Applies Neumann conditions to the right-hand side (RHS) values.
   *
   * This method updates the RHS values of the finite element method equations
   * based on the provided Neumann boundary condition. The boundary condition
   * can specify a value or its components along the x and y directions.
   *
   */
  /*---------------------------------------------------------------------------*/

  static inline void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes, VariableNodeReal3 node_coord, VariableDoFReal& rhs_variable_na, IMesh* mesh, Accelerator::RunQueue* queue)
  {
    ARCANE_CHECK_PTR(bs);
    ARCANE_CHECK_PTR(queue);
    ARCANE_CHECK_PTR(mesh);

    FaceGroup group = bs->getSurface();
    bool has_value = bs->hasValue();

    UnstructuredMeshConnectivityView connectivity_view;
    connectivity_view.setMesh(mesh);
    NodeInfoListView nodes_infos(mesh->nodeFamily());
    auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
    auto fn_cv = connectivity_view.faceNode();

    if (has_value) {
      {
        Real value = bs->getValue();

        auto command = Accelerator::makeCommand(queue);
        auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
        auto in_node_coord = Accelerator::viewIn(command, node_coord);
        command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
        {
          Real length = MeshOperation::computeLengthFace(face_lid, fn_cv, in_node_coord);
          for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
            if (nodes_infos.isOwn(node_lid)) {
              Real rhs_value = value * length / 2.0;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        };
      }
    }
    else {
      {
        Real value_x = bs->hasValueX() ? bs->getValueX() : 0.0;
        Real value_y = bs->hasValueY() ? bs->getValueY() : 0.0;

        auto command = Accelerator::makeCommand(queue);
        auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
        auto in_node_coord = Accelerator::viewIn(command, node_coord);
        FaceInfoListView faces_infos(mesh->faceFamily());
        command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
        {
          Real length = MeshOperation::computeLengthFace(face_lid, fn_cv, in_node_coord);
          Real2 normal = MeshOperation::computeNormalFace(face_lid, fn_cv, in_node_coord, faces_infos);
          for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
            if (nodes_infos.isOwn(node_lid)) {
              Real rhs_value = (normal.x * value_x + normal.y * value_y) * length / 2.0;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        };
      }
    }
  }
};

class BoundaryConditions3D
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
   */
  /*---------------------------------------------------------------------------*/

  static inline void applyConstantSourceToRhs(Real qdot, const FemDoFsOnNodes& dofs_on_nodes, VariableNodeReal3 node_coord, VariableDoFReal& rhs_variable_na, IMesh* mesh, Accelerator::RunQueue* queue)
  {
    BoundaryCondtionsHelpers::applyConstantSourceToRhsBase([] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord) { return MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord); }, qdot, queue, rhs_variable_na, mesh, dofs_on_nodes, node_coord);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Applies Neumann conditions to the right-hand side (RHS) values.
   *
   * This method updates the RHS values of the finite element method equations
   * based on the provided Neumann boundary condition. The boundary condition
   * can specify a value or its components along the x and y directions.
   *
   */
  /*---------------------------------------------------------------------------*/

  static inline void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes, VariableNodeReal3 node_coord, VariableDoFReal& rhs_variable_na, IMesh* mesh, Accelerator::RunQueue* queue)
  {
    ARCANE_CHECK_PTR(bs);
    ARCANE_CHECK_PTR(queue);
    ARCANE_CHECK_PTR(mesh);

    FaceGroup group = bs->getSurface();
    bool has_value = bs->hasValue();

    UnstructuredMeshConnectivityView connectivity_view;
    connectivity_view.setMesh(mesh);
    NodeInfoListView nodes_infos(mesh->nodeFamily());
    auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
    auto fn_cv = connectivity_view.faceNode();

    if (has_value) {
      {
        Real value = bs->getValue();

        auto command = Accelerator::makeCommand(queue);
        auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
        auto in_node_coord = Accelerator::viewIn(command, node_coord);
        command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
        {
          Real area = MeshOperation::computeAreaTria(face_lid, fn_cv, in_node_coord);
          for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
            if (nodes_infos.isOwn(node_lid)) {
              Real rhs_value = value * area / 3.0;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        };
      }
    }
    else {
      {
        Real value_x = bs->hasValueX() ? bs->getValueX() : 0.0;
        Real value_y = bs->hasValueY() ? bs->getValueY() : 0.0;
        Real value_z = bs->hasValueZ() ? bs->getValueZ() : 0.0;

        auto command = Accelerator::makeCommand(queue);
        auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
        auto in_node_coord = Accelerator::viewIn(command, node_coord);
        FaceInfoListView faces_infos(mesh->faceFamily());
        command << RUNCOMMAND_ENUMERATE(FaceLocalId, face_lid, group)
        {
          Real area = MeshOperation::computeAreaTria(face_lid, fn_cv, in_node_coord);
          Real3 normal = MeshOperation::computeNormalTriangle(face_lid, fn_cv, in_node_coord, faces_infos);
          for (NodeLocalId node_lid : fn_cv.nodes(face_lid)) {
            if (nodes_infos.isOwn(node_lid)) {
              Real rhs_value = (normal.x * value_x + normal.y * value_y + normal.z * value_z) * area / 3.0;
              Accelerator::doAtomic<Accelerator::eAtomicOperation::Add>(in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)], rhs_value);
            }
          }
        };
      }
    }
  }
};

} // namespace Arcane::FemUtils::Gpu

#endif // ! ARCANE_FEM_FUNCTIONS_GPU_H
