// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArcaneFemFunctionsGpu.h                                     (C) 2000-2026 */
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
struct Quad4GaussPointInfo
{
  RealVector<4> dN_dx; // Derivatives of shape functions in x {∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ∂𝑁₃/∂𝑥  ∂𝑁₄/∂𝑥}
  RealVector<4> dN_dy; // Derivatives of shape functions in y {∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ∂𝑁₃/∂𝑦  ∂𝑁₄/∂𝑦}
  Real det_j; // Determinant of the Jacobian matrix at the Gauss point.
};

ARCCORE_HOST_DEVICE inline Quad4GaussPointInfo
computeGradientsAndJacobianQuad4Gpu(CellLocalId cell_lid,
                                    const IndexedCellNodeConnectivityView& cn_cv,
                                    const Accelerator::VariableNodeReal3InView& in_node_coord,
                                    Real xi, Real eta)
{
  // Shape function derivatives ∂𝐍/∂ξ and ∂𝐍/∂η
  //     ∂𝐍/∂ξ = [ ∂C₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ ]
  //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η ]
  const Real dN_dxi[4]  = { -0.25 * (1.0 - eta),  0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta) };
  const Real dN_deta[4] = { -0.25 * (1.0 - xi),  -0.25 * (1.0 + xi),  0.25 * (1.0 + xi),   0.25 * (1.0 - xi) };

  // Jacobian calculation 𝑱
  //    𝑱 = [ 𝒋₀₀  𝒋₀₁ ] = [ ∂𝑥/∂ξ  ∂𝑦/∂ξ ]
  //        [ 𝒋₁₀  𝒋₁₁ ]   [ ∂𝑥/∂η  ∂𝑦/∂η ]
  //
  // The Jacobian is computed as follows:
  //   𝒋₀₀ = ∑ (∂𝑁ᵢ/∂ξ * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟒
  //   𝒋₀₁ = ∑ (∂𝑁ᵢ/∂ξ * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟒
  //   𝒋₁₀ = ∑ (∂𝑁ᵢ/∂η * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟒
  //   𝒋₁₁ = ∑ (∂𝑁ᵢ/∂η * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟒
  Real2x2 J;
  J[0][0] = 0.0; J[0][1] = 0.0;
  J[1][0] = 0.0; J[1][1] = 0.0;

  for (Int8 a = 0; a < 4; ++a) {
    const auto& coord = in_node_coord[cn_cv.nodeId(cell_lid, a)];
    J[0][0] += dN_dxi[a] * coord.x;
    J[0][1] += dN_dxi[a] * coord.y;
    J[1][0] += dN_deta[a] * coord.x;
    J[1][1] += dN_deta[a] * coord.y;
  }

  // Determinant of the Jacobian
  const Real detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

  // TODO : Handle non-positive Jacobian determinant case appropriately for GPU execution.
  //if (detJ <= 0.0) {
  //  ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
  //}

  // Inverse of the Jacobian
  //    𝑱⁻¹ = [ invJ00 invJ01 ]
  //          [ invJ10 invJ11 ]
  const Real invJ00 = J[1][1] / detJ;
  const Real invJ01 = -J[0][1] / detJ;
  const Real invJ10 = -J[1][0] / detJ;
  const Real invJ11 = J[0][0] / detJ;

  //   Gradients in physical space (∂𝐍/∂𝑥, ∂𝐍/∂𝑦)
  //    {∂𝐍/∂𝑥} = [J]⁻¹ {∂𝐍/∂ξ}
  //    {∂𝐍/∂𝑦}         {∂𝐍/∂η}
  RealVector<4> dN_dx_result, dN_dy_result;
  for (Int8 a = 0; a < 4; ++a) {
    dN_dx_result(a) = invJ00 * dN_dxi[a] + invJ01 * dN_deta[a];
    dN_dy_result(a) = invJ10 * dN_dxi[a] + invJ11 * dN_deta[a];
  }

  return { dN_dx_result, dN_dy_result, detJ };
}

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
  Real V6 = math::abs(math::dot(v0, math::cross(v1, v2)));

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
  Real V6 = math::abs(math::dot(v0, math::cross(v1, v2)));

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
  Real V6 = math::abs(math::dot(v0, math::cross(v1, v2)));

  Real4 dz{};

  dz[0] = (n1.x * (n3.y - n2.y) + n2.x * (n1.y - n3.y) + n3.x * (n2.y - n1.y)) / V6;
  dz[1] = (n0.x * (n2.y - n3.y) + n2.x * (n3.y - n0.y) + n3.x * (n0.y - n2.y)) / V6;
  dz[2] = (n0.x * (n3.y - n1.y) + n1.x * (n0.y - n3.y) + n3.x * (n1.y - n0.y)) / V6;
  dz[3] = (n0.x * (n1.y - n2.y) + n1.x * (n2.y - n0.y) + n2.x * (n0.y - n1.y)) / V6;

  return dz;
};

/*---------------------------------------------------------------------------*/
/**
  * @brief Holds information for a Hexa8 element at a single Gauss point.
  *
  * This includes the gradients of the shape functions in the physical space (𝑥,𝑦,𝑧)
  * and the determinant of the Jacobian matrix.
  */
/*---------------------------------------------------------------------------*/
struct Hexa8GaussPointInfo
{
  RealVector<8> dN_dx; // Derivatives of shape functions in x {∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ...  ∂𝑁₈/∂𝑥}
  RealVector<8> dN_dy; // Derivatives of shape functions in y {∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ...  ∂𝑁₈/∂𝑦}
  RealVector<8> dN_dz; // Derivatives of shape functions in z {∂𝑁₁/∂𝑧  ∂𝑁₂/∂𝑧  ...  ∂𝑁₈/∂𝑧}
  Real det_j; // Determinant of the Jacobian matrix at the Gauss point.
};

/*---------------------------------------------------------------------------*/
/**
  * @brief Computes shape function gradients and the Jacobian determinant for a Hexa8 element.
  *
  * @param cell_lid The local ID of the cell in the mesh.
  * @param cn_cv The connectivity view for the cell nodes.
  * @param in_node_coord The coordinates of the mesh nodes.
  * @param xi The ξ coordinate of the evaluation point (-1 to 1).
  * @param eta The η coordinate of the evaluation point (-1 to 1).
  * @param zeta The ζ coordinate of the evaluation point (-1 to 1).
  * @return A Hexa8GaussPointInfo struct containing the gradients and Jacobian determinant.
  */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline Hexa8GaussPointInfo
computeGradientsAndJacobianHexa8Gpu(CellLocalId cell_lid,
                                    const IndexedCellNodeConnectivityView& cn_cv,
                                    const Accelerator::VariableNodeReal3InView& in_node_coord,
                                    Real xi, Real eta, Real zeta)
{
  // Shape function derivatives ∂𝐍/∂ξ, ∂𝐍/∂η, ∂𝐍/∂ζ
  //     ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ  ∂𝑁₅/∂ξ  ∂𝑁₆/∂ξ  ∂𝑁₇/∂ξ  ∂𝑁₈/∂ξ ]
  //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η  ∂𝑁₅/∂η  ∂𝑁₆/∂η  ∂𝑁₇/∂η  ∂𝑁₈/∂η ]
  //     ∂𝐍/∂ζ = [ ∂𝑁₁/∂ζ  ∂𝑁₂/∂ζ  ∂𝑁₃/∂ζ  ∂𝑁₄/∂ζ  ∂𝑁₅/∂ζ  ∂𝑁₆/∂ζ  ∂𝑁₇/∂ζ  ∂𝑁₈/∂ζ ]
  Real dN_dxi[8], dN_deta[8], dN_dzeta[8];
  const Real one_minus_eta = 1.0 - eta;
  const Real one_plus_eta = 1.0 + eta;
  const Real one_minus_xi = 1.0 - xi;
  const Real one_plus_xi = 1.0 + xi;
  const Real one_minus_zeta = 1.0 - zeta;
  const Real one_plus_zeta = 1.0 + zeta;

  dN_dxi[0] = -0.125 * one_minus_eta * one_minus_zeta;
  dN_dxi[1] = 0.125 * one_minus_eta * one_minus_zeta;
  dN_dxi[2] = 0.125 * one_plus_eta * one_minus_zeta;
  dN_dxi[3] = -0.125 * one_plus_eta * one_minus_zeta;
  dN_dxi[4] = -0.125 * one_minus_eta * one_plus_zeta;
  dN_dxi[5] = 0.125 * one_minus_eta * one_plus_zeta;
  dN_dxi[6] = 0.125 * one_plus_eta * one_plus_zeta;
  dN_dxi[7] = -0.125 * one_plus_eta * one_plus_zeta;

  dN_deta[0] = -0.125 * one_minus_xi * one_minus_zeta;
  dN_deta[1] = -0.125 * one_plus_xi * one_minus_zeta;
  dN_deta[2] = 0.125 * one_plus_xi * one_minus_zeta;
  dN_deta[3] = 0.125 * one_minus_xi * one_minus_zeta;
  dN_deta[4] = -0.125 * one_minus_xi * one_plus_zeta;
  dN_deta[5] = -0.125 * one_plus_xi * one_plus_zeta;
  dN_deta[6] = 0.125 * one_plus_xi * one_plus_zeta;
  dN_deta[7] = 0.125 * one_minus_xi * one_plus_zeta;

  dN_dzeta[0] = -0.125 * one_minus_xi * one_minus_eta;
  dN_dzeta[1] = -0.125 * one_plus_xi * one_minus_eta;
  dN_dzeta[2] = -0.125 * one_plus_xi * one_plus_eta;
  dN_dzeta[3] = -0.125 * one_minus_xi * one_plus_eta;
  dN_dzeta[4] = 0.125 * one_minus_xi * one_minus_eta;
  dN_dzeta[5] = 0.125 * one_plus_xi * one_minus_eta;
  dN_dzeta[6] = 0.125 * one_plus_xi * one_plus_eta;
  dN_dzeta[7] = 0.125 * one_minus_xi * one_plus_eta;

  // Jacobian matrix (default-initialized to zero see Real3x3.h)
  //    𝑱 = [ 𝒋₀₀  𝒋₀₁  𝒋₀₂ ]
  //        [ 𝒋₁₀  𝒋₁₁  𝒋₁₂ ]
  //        [ 𝒋₂₀  𝒋₂₁  𝒋₂₂ ]
  //
  // The Jacobian is computed as follows:
  //   𝒋₀₀ = ∑ (∂𝑥/∂ξ * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟖
  //   𝒋₀₁ = ∑ (∂𝑥/∂ξ * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟖
  //   𝒋₀₂ = ∑ (∂𝑥/∂ξ * 𝑧ᵢ) ∀ 𝑖= 𝟏,……,𝟖
  //   𝒋₁₀ = ∑ (∂𝑦/∂η * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟖
  //   𝒋₁₁ = ∑ (∂𝑦/∂η * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟖
  //   𝒋₁₂ = ∑ (∂𝑦/∂η * 𝑧ᵢ) ∀ 𝑖= 𝟏,……,𝟖
  //   𝒋₂₀ = ∑ (∂𝑧/∂ζ * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟖
  //   𝒋₂₁ = ∑ (∂𝑧/∂ζ * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟖
  //   𝒋₂₂ = ∑ (∂𝑧/∂ζ * 𝑧ᵢ) ∀ 𝑖= 𝟏,……,𝟖
  Real3x3 J;
  for (Int8 a = 0; a < 8; ++a) {
    const Real3& n_coord = in_node_coord[cn_cv.nodeId(cell_lid, a)];
    J[0][0] += dN_dxi[a] * n_coord.x; // ∂𝑥/∂ξ
    J[0][1] += dN_dxi[a] * n_coord.y; // ∂𝑦/∂ξ
    J[0][2] += dN_dxi[a] * n_coord.z; // ∂𝑧/∂ξ
    J[1][0] += dN_deta[a] * n_coord.x; // ∂𝑥/∂η
    J[1][1] += dN_deta[a] * n_coord.y; // ∂𝑦/∂η
    J[1][2] += dN_deta[a] * n_coord.z; // ∂𝑧/∂η
    J[2][0] += dN_dzeta[a] * n_coord.x; // ∂𝑥/∂ζ
    J[2][1] += dN_dzeta[a] * n_coord.y; // ∂𝑦/∂ζ
    J[2][2] += dN_dzeta[a] * n_coord.z; // ∂𝑧/∂ζ
  }

  // Determinant and Inverse of the Jacobian
  const Real detJ = math::matrixDeterminant(J);
  // if (detJ <= 0.0) {
  //   ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
  // }
  const Real3x3 invJ = math::inverseMatrix(J, detJ);

  // Gradients in physical space (∂𝐍/∂𝑥, ∂𝐍/∂𝑦, ∂𝐍/∂𝑧)
  //    {∂𝐍/∂𝑥} = [J]⁻¹ {∂𝐍/∂ξ}
  //    {∂𝐍/∂𝑦}         {∂𝐍/∂η}
  //    {∂𝐍/∂𝑧}         {∂𝐍/∂ζ}
  RealVector<8> dN_dx_result, dN_dy_result, dN_dz_result;
  for (Int8 a = 0; a < 8; ++a) {
    dN_dx_result(a) = invJ[0][0] * dN_dxi[a] + invJ[0][1] * dN_deta[a] + invJ[0][2] * dN_dzeta[a];
    dN_dy_result(a) = invJ[1][0] * dN_dxi[a] + invJ[1][1] * dN_deta[a] + invJ[1][2] * dN_dzeta[a];
    dN_dz_result(a) = invJ[2][0] * dN_dxi[a] + invJ[2][1] * dN_deta[a] + invJ[2][2] * dN_dzeta[a];
  }

  return { dN_dx_result, dN_dy_result, dN_dz_result, detJ };
}

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

  static void applyConstantSourceToRhsQuad4(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
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

  static void applyNeumannToRhsQuad4(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
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

  static void applyConstantSourceToRhsHexa8(Real qdot, const FemDoFsOnNodes& dofs_on_nodes,
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

  static void applyNeumannToHexa8(BC::INeumannBoundaryCondition* bs, const FemDoFsOnNodes& dofs_on_nodes,
                       const VariableNodeReal3& node_coord, VariableDoFReal& rhs_variable_na,
                       IMesh* mesh, Accelerator::RunQueue* queue);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils::Gpu

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
