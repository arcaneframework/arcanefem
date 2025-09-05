// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArcaneFemFunctionsGpu.h                                     (C) 2022-2025 */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMTUILS_ARCANEFEMFUNCTIONSGPU_H
#define ARCANEFEM_FEMTUILS_ARCANEFEMFUNCTIONSGPU_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/UnstructuredMeshConnectivity.h>
#include <arcane/core/IndexedItemConnectivityView.h>
#include <arcane/core/ItemInfoListView.h>
#include <arcane/core/ItemEnumerator.h>
#include <arcane/core/MathUtils.h>
#include <arcane/core/ItemTypes.h>
#include <arcane/core/DataView.h>
#include <arcane/core/Item.h>
#include <arcane/core/VariableTypes.h>
#include <arcane/core/IMesh.h>

#include <arcane/accelerator/RunCommandEnumerate.h>
#include <arcane/accelerator/VariableViews.h>
#include <arcane/accelerator/RunCommand.h>
#include <arcane/accelerator/RunQueue.h>
#include <arcane/accelerator/Atomic.h>

#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "IArcaneFemBC.h"
#include "FemUtils.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils::Gpu::Csr
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Int32
findIndex(Int32 begin, Int32 end, Int32 column_lid, Span<const Int32> in_csr_columns)
{
  for (auto i = begin; i < end; ++i)
    if (in_csr_columns[i] == column_lid)
      return i;
  return -1;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils::Gpu::Csr

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils::Gpu::MeshOperation
{

/*---------------------------------------------------------------------------*/
/**
 * @brief computes the area of a triangle defined by three nodes.
 *
 * The area is determined as half the magnitude of the cross-product 
 * of two edge vectors originating from the first node.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real
computeAreaTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv,
                 const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 n0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 n1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 n2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];

  auto v = math::cross(n1 - n0, n2 - n0);
  return v.normL2() / 2.0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real
computeAreaTria(FaceLocalId face_lid,
                const IndexedFaceNodeConnectivityView& fn_cv,
                const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 n0 = in_node_coord[fn_cv.nodeId(face_lid, 0)];
  Real3 n1 = in_node_coord[fn_cv.nodeId(face_lid, 1)];
  Real3 n2 = in_node_coord[fn_cv.nodeId(face_lid, 2)];

  auto v = math::cross(n1 - n0, n2 - n0);
  return v.normL2() / 2.0;
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the volume of a tetrahedra defined by four nodes.
 *
 * The volume is given by one-sixth of the absolute value of the determinant
 * formed by three edge vectors originating from the first node.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real
computeVolumeTetra4(CellLocalId cell_lid,
                    const IndexedCellNodeConnectivityView& cn_cv,
                    const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 n0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 n1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 n2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real3 n3 = in_node_coord[cn_cv.nodeId(cell_lid, 3)];

  Real3 v0 = n1 - n0;
  Real3 v1 = n2 - n0;
  Real3 v2 = n3 - n0;

  return math::abs(math::dot(v0, math::cross(v1, v2))) / 6.0;
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the length of the edge defined by a given face.
 *
 * This method calculates Euclidean distance between the two nodes of the face.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real
computeLengthFace(FaceLocalId face_lid,
                  IndexedFaceNodeConnectivityView fn_cv,
                  Accelerator::VariableNodeReal3InView in_node_coord)
{
  Real3 n0 = in_node_coord[fn_cv.nodeId(face_lid, 0)];
  Real3 n1 = in_node_coord[fn_cv.nodeId(face_lid, 1)];

  return math::sqrt((n1.x - n0.x) * (n1.x - n0.x) + (n1.y - n0.y) * (n1.y - n0.y));
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the normalized edge normal for a given face.
 *
 * This method calculates normal vector to the edge defined by nodes of the face,
 * normalizes it, and ensures the correct orientation.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real2
computeNormalFace(FaceLocalId face_lid,
                  IndexedFaceNodeConnectivityView fn_cv,
                  Accelerator::VariableNodeReal3InView in_node_coord,
                  FaceInfoListView faces_infos)
{
  Real3 n0 = in_node_coord[fn_cv.nodeId(face_lid, 0)];
  Real3 n1 = in_node_coord[fn_cv.nodeId(face_lid, 1)];

  if (!faces_infos.isSubDomainBoundaryOutside(face_lid)) {
    Real3 tmp = n0;
    n0 = n1;
    n1 = tmp;
  }

  Real2 N;
  Real norm_N = math::sqrt((n1.y - n0.y) * (n1.y - n0.y) + (n1.x - n0.x) * (n1.x - n0.x));
  N.x = (n1.y - n0.y) / norm_N;
  N.y = (n0.x - n1.x) / norm_N;
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

ARCCORE_HOST_DEVICE inline Real3
computeNormalTriangle(FaceLocalId face_lid,
                      IndexedFaceNodeConnectivityView fn_cv,
                      Accelerator::VariableNodeReal3InView in_node_coord,
                      FaceInfoListView faces_infos)
{

  Real3 n0 = in_node_coord[fn_cv.nodeId(face_lid, 0)];
  Real3 n1 = in_node_coord[fn_cv.nodeId(face_lid, 1)];
  Real3 n2 = in_node_coord[fn_cv.nodeId(face_lid, 2)];

  if (!faces_infos.isSubDomainBoundaryOutside(face_lid))
    std::swap(n0, n1);

  Real3 edge1 = { n1.x - n0.x, n1.y - n0.y, n1.z - n0.z };
  Real3 edge2 = { n2.x - n0.x, n2.y - n0.y, n2.z - n0.z };

  Real3 normal = {
    edge1.y * edge2.z - edge1.z * edge2.y,
    edge1.z * edge2.x - edge1.x * edge2.z,
    edge1.x * edge2.y - edge1.y * edge2.x
  };

  Real norm = math::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
  return { normal.x / norm, normal.y / norm, normal.z / norm };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils::Gpu::MeshOperation

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils::Gpu::FeOperation2D
{

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the 𝑥 gradients of basis functions 𝐍 for ℙ1 triangles.
 *
 * This method calculates gradient operator ∂/∂𝑥 of 𝑁ᵢ for a given ℙ1
 * cell with i = 1,..,3 for the three shape function  𝑁ᵢ  hence output
 * is a vector of size 3
 *
 *         ∂𝐍/∂𝑥 = [ ∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ∂𝑁₃/∂𝑥 ]
 *
 *         ∂𝐍/∂𝑥 = 1/(2𝐴) [ y2-y3  y3−y1  y1−y2 ]
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real3
computeGradientXTria3(CellLocalId cell_lid,
                      const IndexedCellNodeConnectivityView& cn_cv,
                      const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 n0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 n1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 n2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];

  auto A2 = ((n1.x - n0.x) * (n2.y - n0.y) - (n2.x - n0.x) * (n1.y - n0.y));

  return { (n1.y - n2.y) / A2, (n2.y - n0.y) / A2, (n0.y - n1.y) / A2 };
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the 𝑦 gradients of basis functions 𝐍 for ℙ1 triangles.
 *
 * This method calculates gradient operator ∂/∂𝑦 of 𝑁ᵢ for a given ℙ1
 * cell with i = 1,..,3 for the three shape function  𝑁ᵢ  hence output
 * is a vector of size 3
 *
 *         ∂𝐍/∂𝑥 = [ ∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ∂𝑁₃/∂𝑦 ]
 *
 *         ∂𝐍/∂𝑥 = 1/(2𝐴) [ 𝑥₃−𝑥₂  𝑥₁−𝑥₃  𝑥₂−𝑥₁ ]
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real3
computeGradientYTria3(CellLocalId cell_lid,
                      const IndexedCellNodeConnectivityView& cn_cv,
                      const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 n0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 n1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 n2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];

  auto A2 = ((n1.x - n0.x) * (n2.y - n0.y) - (n2.x - n0.x) * (n1.y - n0.y));

  return { (n2.x - n1.x) / A2, (n0.x - n2.x) / A2, (n1.x - n0.x) / A2 };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils::Gpu::FeOperation2D

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils::Gpu::FeOperation3D
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
/**
 * @brief Computes the X gradients of basis functions 𝐍 for ℙ1 Tetrahedron.
 *
 * This method calculates gradient operator ∂/∂𝑥 of 𝑁ᵢ for a given ℙ1
 * cell with i = 1,..,4 for the four shape function  𝑁ᵢ  hence output
 * is a vector of size 4
 *
 *         ∂𝐍/∂𝑥 = [ ∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ∂𝑁₃/∂𝑥  ∂𝑁₄/∂𝑥 ]
 *
 *         ∂𝐍/∂𝑥 = 1/(6𝑉) [ b0  b1  b2  b3 ]
 *
 * where:
 *    b0 = (n1.y * (n3.z - n2.z) + n2.y * (n1.z - n3.z) + n3.y * (n2.z - n1.z)),
 *    b1 = (n0.y * (n2.z - n3.z) + n2.y * (n3.z - n0.z) + n3.y * (n0.z - n2.z)),
 *    b2 = (n0.y * (n3.z - n1.z) + n1.y * (n0.z - n3.z) + n3.y * (n1.z - n0.z)),
 *    b3 = (n0.y * (n1.z - n2.z) + n1.y * (n2.z - n0.z) + n2.y * (n0.z - n1.z)).
 *
 */
/*-------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real4
computeGradientXTetra4(CellLocalId cell_lid,
                       const IndexedCellNodeConnectivityView& cn_cv,
                       const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 n0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 n1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 n2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real3 n3 = in_node_coord[cn_cv.nodeId(cell_lid, 3)];

  Real3 v0 = n1 - n0;
  Real3 v1 = n2 - n0;
  Real3 v2 = n3 - n0;

  // 6 x Volume of tetrahedron
  Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

  Real4 dx{};

  dx[0] = (n1.y * (n3.z - n2.z) + n2.y * (n1.z - n3.z) + n3.y * (n2.z - n1.z)) / V6;
  dx[1] = (n0.y * (n2.z - n3.z) + n2.y * (n3.z - n0.z) + n3.y * (n0.z - n2.z)) / V6;
  dx[2] = (n0.y * (n3.z - n1.z) + n1.y * (n0.z - n3.z) + n3.y * (n1.z - n0.z)) / V6;
  dx[3] = (n0.y * (n1.z - n2.z) + n1.y * (n2.z - n0.z) + n2.y * (n0.z - n1.z)) / V6;

  return dx;
};

/*-------------------------------------------------------------------------*/
/**
 * @brief Computes the Y gradients of basis functions 𝐍 for ℙ1 Tetrahedron.
 *
 * This method calculates gradient operator ∂/∂𝑦 of 𝑁ᵢ for a given ℙ1
 * cell with i = 1,..,4 for the four shape functions 𝑁ᵢ, hence the output
 * is a vector of size 4.
 *
 *         ∂𝐍/∂𝑦 = [ ∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ∂𝑁₃/∂𝑦  ∂𝑁₄/∂𝑦 ]
 *
 *         ∂𝐍/∂𝑦 = 1/(6𝑉) [ c0  c1  c2  c3 ]
 *
 * where:
 *    c0 = (n1.z * (n3.x - n2.x) + n2.z * (n1.x - n3.x) + n3.z * (n2.x - n1.x)),
 *    c1 = (n0.z * (n2.x - n3.x) + n2.z * (n3.x - n0.x) + n3.z * (n0.x - n2.x)),
 *    c2 = (n0.z * (n3.x - n1.x) + n1.z * (n0.x - n3.x) + n3.z * (n1.x - n0.x)),
 *    c3 = (n0.z * (n1.x - n2.x) + n1.z * (n2.x - n0.x) + n2.z * (n0.x - n1.x)).
 *
 */
/*-------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real4
computeGradientYTetra4(CellLocalId cell_lid,
                       const IndexedCellNodeConnectivityView& cn_cv,
                       const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 n0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 n1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 n2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real3 n3 = in_node_coord[cn_cv.nodeId(cell_lid, 3)];

  Real3 v0 = n1 - n0;
  Real3 v1 = n2 - n0;
  Real3 v2 = n3 - n0;

  // 6 x Volume of tetrahedron
  Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

  Real4 dy{};

  dy[0] = (n1.z * (n3.x - n2.x) + n2.z * (n1.x - n3.x) + n3.z * (n2.x - n1.x)) / V6;
  dy[1] = (n0.z * (n2.x - n3.x) + n2.z * (n3.x - n0.x) + n3.z * (n0.x - n2.x)) / V6;
  dy[2] = (n0.z * (n3.x - n1.x) + n1.z * (n0.x - n3.x) + n3.z * (n1.x - n0.x)) / V6;
  dy[3] = (n0.z * (n1.x - n2.x) + n1.z * (n2.x - n0.x) + n2.z * (n0.x - n1.x)) / V6;

  return dy;
}

/*-------------------------------------------------------------------------*/
/**
 * @brief Computes the 𝑧 gradients of basis functions 𝐍 for ℙ1 Tetrahedron.
 *
 * This method calculates gradient operator ∂/∂𝑧 of 𝑁ᵢ for a given ℙ1
 * cell with i = 1,..,4 for the four shape functions 𝑁ᵢ, hence the output
 * is a vector of size 4.
 *
 *         ∂𝐍/∂𝑧 = [ ∂𝑁₁/∂𝑧  ∂𝑁₂/∂𝑧  ∂𝑁₃/∂𝑧  ∂𝑁₄/∂𝑧 ]
 *
 *         ∂𝐍/∂𝑧 = 1/(6𝑉) [ d0  d1  d2  d3 ]
 *
 * where:
 *    d0 = (n1.x * (n3.y - n2.y) + n2.x * (n1.y - n3.y) + n3.x * (n2.y - n1.y)),
 *    d1 = (n0.x * (n2.y - n3.y) + n2.x * (n3.y - n0.y) + n3.x * (n0.y - n2.y)),
 *    d2 = (n0.x * (n3.y - n1.y) + n1.x * (n0.y - n3.y) + n3.x * (n1.y - n0.y)),
 *    d3 = (n0.x * (n1.y - n2.y) + n1.x * (n2.y - n0.y) + n2.x * (n0.y - n1.y)).
 *
 */
/*-------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Real4
computeGradientZTetra4(CellLocalId cell_lid,
                       const IndexedCellNodeConnectivityView& cn_cv,
                       const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real3 n0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 n1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 n2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real3 n3 = in_node_coord[cn_cv.nodeId(cell_lid, 3)];

  auto v0 = n1 - n0;
  auto v1 = n2 - n0;
  auto v2 = n3 - n0;

  // 6 x Volume of tetrahedron
  Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

  Real4 dz{};

  dz[0] = (n1.x * (n3.y - n2.y) + n2.x * (n1.y - n3.y) + n3.x * (n2.y - n1.y)) / V6;
  dz[1] = (n0.x * (n2.y - n3.y) + n2.x * (n3.y - n0.y) + n3.x * (n0.y - n2.y)) / V6;
  dz[2] = (n0.x * (n3.y - n1.y) + n1.x * (n0.y - n3.y) + n3.x * (n1.y - n0.y)) / V6;
  dz[3] = (n0.x * (n1.y - n2.y) + n1.x * (n2.y - n0.y) + n2.x * (n0.y - n1.y)) / V6;

  return dz;
};

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils::Gpu::FeOperation3D

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

namespace Arcane::FemUtils::Gpu
{

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

class BoundaryConditionsHelpers
{
 public:

  template <class Function>
  static void applyConstantSourceToRhsBase(Function computeSpatialIntegral, Real qdot,
                                           RunQueue* queue, VariableDoFReal& rhs_variable_na, IMesh* mesh,
                                           const FemDoFsOnNodes& dofs_on_nodes, const VariableNodeReal3& node_coord)
  {
    ARCANE_CHECK_PTR(queue);
    ARCANE_CHECK_PTR(mesh);
    ARCANE_CHECK_PTR(computeSpatialIntegral);

    UnstructuredMeshConnectivityView connectivity_view;
    connectivity_view.setMesh(mesh);
    NodeInfoListView nodes_infos(mesh->nodeFamily());
    auto node_dof(dofs_on_nodes.nodeDoFConnectivityView());
    auto cn_cv = connectivity_view.cellNode();
    auto nc_cv = connectivity_view.nodeCell();

    auto command = Accelerator::makeCommand(queue);

    auto in_out_rhs_variable_na = Accelerator::viewInOut(command, rhs_variable_na);
    auto in_node_coord = Accelerator::viewIn(command, node_coord);

    // Iterate over all nodes and compute the contribution to the RHS
    command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, mesh->allNodes())
    {
      if (nodes_infos.isOwn(node_lid)) {
        Real sum = 0.0;
        for (CellLocalId cell_lid : nc_cv.cells(node_lid)) {
          Real domain = computeSpatialIntegral(cell_lid, cn_cv, in_node_coord);
          sum += qdot * domain / cn_cv.nbItem(cell_lid);
        }
        in_out_rhs_variable_na[node_dof.dofId(node_lid, 0)] = sum;
      }
    };
  }

  static void
  applyDirichletToNodeGroupRhsOnly(Int32 dof_index, Real value, RunQueue* queue,
                                   IMesh* mesh, DoFLinearSystem& linear_system, const FemDoFsOnNodes& dofs_on_nodes,
                                   const NodeGroup& node_group);

  static void
  applyDirichletToNodeGroupViaPenalty(Int32 dof_index, Real value, Real penalty, RunQueue* queue,
                                      IMesh* mesh, DoFLinearSystem& linear_system, const FemDoFsOnNodes& dofs_on_nodes,
                                      const NodeGroup& node_group);

  static void
  applyDirichletToNodeGroupViaRowOrRowColumnElimination(Byte elimination_type, Int32 dof_index, Real value, RunQueue* queue,
                                                        DoFLinearSystem& linear_system, const FemDoFsOnNodes& dofs_on_nodes,
                                                        const NodeGroup& node_group);
};

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

class BoundaryConditions
{
 public:

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Applies Dirichlet boundary conditions to RHS and LHS.
   *
   * - For LHS matrix `𝐀`, the diagonal term for the Dirichlet DOF is set to `𝑃`.
   * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
   */
  /*---------------------------------------------------------------------------*/

  static void
  applyDirichletToLhsAndRhs(BC::IDirichletBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                            DoFLinearSystem& linear_system,
                            IMesh* mesh, Accelerator::RunQueue* queue);

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Applies Point Dirichlet boundary conditions to RHS and LHS.
   *
   * - For LHS matrix `𝐀`, the diagonal term for the Dirichlet DOF is set to `𝑃`.
   * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
   */
  /*---------------------------------------------------------------------------*/

  static void
  applyPointDirichletToLhsAndRhs(BC::IDirichletPointCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                                 DoFLinearSystem& linear_system,
                                 IMesh* mesh, RunQueue* queue);

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Applies Dirichlet boundary conditions to RHS.
   *
   * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
   */
  /*---------------------------------------------------------------------------*/

  static void
  applyDirichletToRhs(BC::IDirichletBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                      DoFLinearSystem& linear_system,
                      IMesh* mesh, Accelerator::RunQueue* queue);

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Applies Point Dirichlet boundary conditions to RHS.
   *
   * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
   */
  /*---------------------------------------------------------------------------*/

  static void
  applyPointDirichletToRhs(BC::IDirichletPointCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                           DoFLinearSystem& linear_system,
                           IMesh* mesh, Accelerator::RunQueue* queue);
};

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

  static void applyConstantSourceToRhs(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
                                       const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                                       IMesh* mesh, Accelerator::RunQueue* queue);

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

  static void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                                const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                                IMesh* mesh, Accelerator::RunQueue* queue);
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

  static void applyConstantSourceToRhs(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
                                       const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                                       IMesh* mesh, Accelerator::RunQueue* queue);

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

  static void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                                const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                                IMesh* mesh, Accelerator::RunQueue* queue);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils::Gpu

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
