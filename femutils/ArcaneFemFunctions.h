// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArcaneFemFunctions.h                                        (C) 2000-2026 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef ARCANE_FEM_FUNCTIONS_H
#define ARCANE_FEM_FUNCTIONS_H

#define M_SQRT1_3 0.57735026918962576451 /* 1/sqrt(3) */

#include <arcane/utils/ITraceMng.h>
#include <arcane/utils/StringList.h>
#include <arcane/utils/CommandLineArguments.h>

#include <arcane/core/IStandardFunction.h>
#include <arcane/core/UnstructuredMeshConnectivity.h>
#include <arcane/core/IndexedItemConnectivityView.h>
#include <arcane/core/VariableTypes.h>
#include <arcane/core/IMesh.h>

#include "IArcaneFemBC.h"
#include "GaussQuadrature.h"
#include "DoFLinearSystem.h"
#include "ShapeFunctions.h"

#include "MeshOperation.h"
#include "FeOperation.h"

using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/**
 * @brief Contains various functions & operations related to FEM calculations.
 *
 * The class provides methods organized into different nested classes for:
 * - MeshOperation: Mesh related operations.
 * - FeOperation2D/3D: Finite element operations at element level.
 * - BoundaryConditions2D/3D: Boundary condition related operations.
 */
/*---------------------------------------------------------------------------*/
namespace ArcaneFemFunctions
{
/*---------------------------------------------------------------------------*/
/**
 * @brief Provides general purpose fuctions
 */
/*---------------------------------------------------------------------------*/
class GeneralFunctions
{
 public:

  static void printArcaneFemTime(ITraceMng* tm, const String& label, const Real& value);

  static CommandLineArguments getPetscFlagsFromCommandline(const String& petsc_flags);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods to help build boundary conditions in 2D/3D.
 */
/*---------------------------------------------------------------------------*/
class BoundaryConditionsHelpers
{
 public:

  static void applyDirichletToNodeGroupRhsOnly(const Int32 dof_index, Real value,
                                               const IndexedNodeDoFConnectivityView& node_dof,
                                               VariableDoFReal& rhs_values,
                                               NodeGroup& node_group);

  static void applyDirichletToNodeGroupViaPenalty(const Int32 dof_index, Real value, Real penalty,
                                                  const IndexedNodeDoFConnectivityView& node_dof,
                                                  DoFLinearSystem& linear_system,
                                                  VariableDoFReal& rhs_values,
                                                  NodeGroup& node_group);

  static void applyDirichletToNodeGroupViaRowElimination(const Int32 dof_index, Real value,
                                                         const IndexedNodeDoFConnectivityView& node_dof,
                                                         DoFLinearSystem& linear_system,
                                                         VariableDoFReal& rhs_values,
                                                         NodeGroup& node_group);

  static void applyDirichletToNodeGroupViaRowColumnElimination(const Int32 dof_index, Real value,
                                                               const IndexedNodeDoFConnectivityView& node_dof,
                                                               DoFLinearSystem& linear_system,
                                                               VariableDoFReal& rhs_values,
                                                               NodeGroup& node_group);
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods for applying boundary conditions in 2D/3D FEM.
 *
 * This class includes static methods for applying:
 * - Dirichlet boundary condition
 * - Point Dirichlet condition
 */
/*---------------------------------------------------------------------------*/
class BoundaryConditions
{
 public:

  static void applyDirichletToLhsAndRhs(BC::IDirichletBoundaryCondition* bs,
                                        const IndexedNodeDoFConnectivityView& node_dof,
                                        DoFLinearSystem& m_linear_system,
                                        VariableDoFReal& rhs_values);

  static void applyPointDirichletToLhsAndRhs(BC::IDirichletPointCondition* bs,
                                             const IndexedNodeDoFConnectivityView& node_dof,
                                             DoFLinearSystem& m_linear_system,
                                             VariableDoFReal& rhs_values);

  static void applyDirichletToRhs(BC::IDirichletBoundaryCondition* bs,
                                  const IndexedNodeDoFConnectivityView& node_dof,
                                  VariableDoFReal& rhs_values);

  static void applyPointDirichletToRhs(BC::IDirichletPointCondition* bs,
                                       const IndexedNodeDoFConnectivityView& node_dof,
                                       VariableDoFReal& rhs_values);

  static void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, IMesh* mesh,
                                const IndexedNodeDoFConnectivityView& node_dof,
                                const VariableNodeReal3& node_coord,
                                VariableDoFReal& rhs_values);

  static void applyConstantSourceToRhs(Real qdot, IMesh* mesh,
                                       const IndexedNodeDoFConnectivityView& node_dof,
                                       const VariableNodeReal3& node_coord,
                                       VariableDoFReal& rhs_values);
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods for applying boundary conditions in 3D FEM problems.
 *
 * This class includes static methods for applying Neumann boundary conditions
 * to the right-hand side (RHS) of finite element method equations in 3D.
 */
/*---------------------------------------------------------------------------*/
class BoundaryConditions3D
{
 public:

  static void applyConstantSourceToRhsTetra4(Real qdot, IMesh* mesh,
                                             const IndexedNodeDoFConnectivityView& node_dof,
                                             const VariableNodeReal3& node_coord,
                                             VariableDoFReal& rhs_values);

  static void applyConstantSourceToRhsHexa8(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values);
  static void applyConstantSourceToRhsHexa20(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof,
                                             const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values);
  static void applyConstantSourceToRhsHexa27(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof,
                                             const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values);

  static void integrateNodalFieldToRhsTetra4(VariableNodeReal& field, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof,
                                             const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values);

  static void integrateNodalFieldToRhsHexa8(VariableNodeReal& field, IMesh* mesh,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);

  static void applyManufacturedSourceToRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_source,
                                           IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof,
                                           const VariableNodeReal3& node_coord,
                                           VariableDoFReal& rhs_values);

  static void applyNeumannToRhsTetra4(BC::INeumannBoundaryCondition* bs,
                                      const IndexedNodeDoFConnectivityView& node_dof,
                                      const VariableNodeReal3& node_coord,
                                      VariableDoFReal& rhs_values);

  static void applyNeumannToRhsHexa8(BC::INeumannBoundaryCondition* bs,
                                     const IndexedNodeDoFConnectivityView& node_dof,
                                     const VariableNodeReal3& node_coord,
                                     VariableDoFReal& rhs_values);

  static void applyNeumannToRhsHexa20(BC::INeumannBoundaryCondition* bs,
                                      const IndexedNodeDoFConnectivityView& node_dof,
                                      const VariableNodeReal3& node_coord,
                                      VariableDoFReal& rhs_values);

  static void applyNeumannToRhsHexa27(BC::INeumannBoundaryCondition* bs,
                                      const IndexedNodeDoFConnectivityView& node_dof,
                                      const VariableNodeReal3& node_coord,
                                      VariableDoFReal& rhs_values);

  static void applyTractionToRhsTetra4(BC::ITractionBoundaryCondition* bs,
                                       const IndexedNodeDoFConnectivityView& node_dof,
                                       const VariableNodeReal3& node_coord,
                                       VariableDoFReal& rhs_values);

  static void applyTractionToRhsHexa8(BC::ITractionBoundaryCondition* bs,
                                      const IndexedNodeDoFConnectivityView& node_dof,
                                      const VariableNodeReal3& node_coord,
                                      VariableDoFReal& rhs_values);

  static void applyTractionTableToRhsTetra4(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
                                            const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);

  static void applyTractionTableToRhsHexa8(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
                                           const ConstArrayView<CaseTableInfo>& traction_case_table_list,
                                           const IndexedNodeDoFConnectivityView& node_dof,
                                           const VariableNodeReal3& node_coord,
                                           VariableDoFReal& rhs_values);

  static void applyManufacturedDirichletToLhsAndRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_dirichlet,
                                                    Real /*lambda*/, const FaceGroup& group,
                                                    BC::IManufacturedSolution* bs,
                                                    const IndexedNodeDoFConnectivityView& node_dof,
                                                    const VariableNodeReal3& node_coord,
                                                    DoFLinearSystem& m_linear_system,
                                                    VariableDoFReal& rhs_values);
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods for applying boundary conditions in 3D FEM problems.
 *
 * This class includes static methods for applying Neumann boundary conditions
 * to the right-hand side (RHS) of finite element method equations in 3D.
 */
/*---------------------------------------------------------------------------*/
class BoundaryConditions2D
{
 public:

  static void applyConstantSourceToRhsTria3(Real qdot, IMesh* mesh,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);

  static void applyConstantSourceToRhsQuad4(Real qdot, IMesh* mesh,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);
  static void applyConstantSourceToRhsQuad8(Real qdot, IMesh* mesh,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);
  static void applyConstantSourceToRhsQuad9(Real qdot, IMesh* mesh,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);

  static void integrateNodalFieldToRhsTria3(VariableNodeReal& field, IMesh* mesh,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);

  static void integrateNodalFieldToRhsQuad4(VariableNodeReal& field, IMesh* mesh,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);

  static void applyManufacturedSourceToRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_source, IMesh* mesh,
                                           const IndexedNodeDoFConnectivityView& node_dof,
                                           const VariableNodeReal3& node_coord,
                                           VariableDoFReal& rhs_values);

  static void applyNeumannToRhsTria3(BC::INeumannBoundaryCondition* bs,
                                     const IndexedNodeDoFConnectivityView& node_dof,
                                     const VariableNodeReal3& node_coord,
                                     VariableDoFReal& rhs_values);

  static void applyNeumannToRhsQuad4(BC::INeumannBoundaryCondition* bs,
                                     const IndexedNodeDoFConnectivityView& node_dof,
                                     const VariableNodeReal3& node_coord,
                                     VariableDoFReal& rhs_values);

  static void applyNeumannToRhsLine3(BC::INeumannBoundaryCondition* bs,
                                     const IndexedNodeDoFConnectivityView& node_dof,
                                     const VariableNodeReal3& node_coord,
                                     VariableDoFReal& rhs_values);

  static void applyNeumannToRhsQuad8(BC::INeumannBoundaryCondition* bs,
                                     const IndexedNodeDoFConnectivityView& node_dof,
                                     const VariableNodeReal3& node_coord,
                                     VariableDoFReal& rhs_values);

  static void applyNeumannToRhsQuad9(BC::INeumannBoundaryCondition* bs,
                                     const IndexedNodeDoFConnectivityView& node_dof,
                                     const VariableNodeReal3& node_coord,
                                     VariableDoFReal& rhs_values);

  static void applyTractionToRhsTria3(BC::ITractionBoundaryCondition* bs,
                                      const IndexedNodeDoFConnectivityView& node_dof,
                                      const VariableNodeReal3& node_coord,
                                      VariableDoFReal& rhs_values);

  static void applyTractionToRhsQuad4(BC::ITractionBoundaryCondition* bs,
                                      const IndexedNodeDoFConnectivityView& node_dof,
                                      const VariableNodeReal3& node_coord,
                                      VariableDoFReal& rhs_values);

  static void applyTractionTableToRhsTria3(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
                                           const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list,
                                           const IndexedNodeDoFConnectivityView& node_dof,
                                           const VariableNodeReal3& node_coord,
                                           VariableDoFReal& rhs_values);

  static void applyTractionTableToRhsQuad4(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
                                           const UniqueArray<CaseTableInfo>& traction_case_table_list,
                                           const IndexedNodeDoFConnectivityView& node_dof,
                                           const VariableNodeReal3& node_coord,
                                           VariableDoFReal& rhs_values);

  static void applyManufacturedDirichletToLhsAndRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_dirichlet, Real /*lambda*/,
                                                    const FaceGroup& group, BC::IManufacturedSolution* bs,
                                                    const IndexedNodeDoFConnectivityView& node_dof,
                                                    const VariableNodeReal3& node_coord,
                                                    DoFLinearSystem& m_linear_system,
                                                    VariableDoFReal& rhs_values);
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods based on the Dispatcher mechanism available in Arcane,
 * allowing to compute FEM methods without declaring the FE entity type
 * (coming from PASSMO).
 */
/*---------------------------------------------------------------------------*/
class CellFEMDispatcher
{

 public:

  Real getShapeFuncVal(Int16 /*item_type*/, Integer /*inod*/, Real3 /*ref coord*/);
  Real3 getShapeFuncDeriv(Int16 /*item_type*/, Integer /*inod*/, Real3 /*ref coord*/);

  RealUniqueArray getGaussData(ItemWithNodes item, Integer nint, Integer ngauss);

  CellFEMDispatcher();

 private:

  std::function<Real(Integer inod, Real3 coord)> m_shapefunc[NB_BASIC_ITEM_TYPE];
  std::function<Real3(Integer inod, Real3 coord)> m_shapefuncderiv[NB_BASIC_ITEM_TYPE];
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods for various FEM-related operations.
 *
 * This class includes static methods for computing shape functions and their
 * derivatives depending on finite element types. These methods are used within
 * the Dispatcher mechanism available in Arcane through the class
 * CellFEMDispatcher (coming from PASSMO).
 */
/*---------------------------------------------------------------------------*/
class FemShapeMethods
{
 public:

  /*---------------------------------------------------------------------------*/
  /**
     * @brief Provides methods for reference linear (P1) edge finite-element
     * The "Line2" reference element is assumed as follows:
     *  0           1
     *  o-----------o---> x
     * -1           1
     * direct local numbering : 0->1
     */
  /*---------------------------------------------------------------------------*/
  static inline Real line2ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 2);
#endif

    Real r = ref_coord[0];
    if (inod == 1)
      return (0.5 * (1 + r));
    return (0.5 * (1 - r));
  }

  static inline Real3 line2ShapeFuncDeriv(Integer inod, Real3)
  {
    if (inod == 1)
      return { 0.5, 0., 0. };
    return { -0.5, 0., 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
     * @brief Provides methods for reference quadratic (P2) edge finite-element
     * The "Line3" reference element is assumed as follows:
     *  0     2      1
     *  o-----o------o---> x
     * -1     0      1
     * direct local numbering : 0->1->2
     */
  /*---------------------------------------------------------------------------*/
  static inline Real line3ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 3);
#endif

    Real ri = ref_coord[0];
    if (inod == 0)
      ri *= -1;

    if (inod < 2)
      return 0.5 * ri * (1 + ri); // nodes 0 or 1
    return (1 - ri * ri); // middle node
  }

  static inline Real3 line3ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    Real ri = ref_coord[0];
    if (!inod)
      return { -0.5 + ri, 0., 0. };
    if (inod == 1)
      return { 0.5 + ri, 0., 0. };
    return { -2. * ri, 0., 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
     * @brief Provides methods for reference linear (P1) triangle finite-element
     * The "Tri3" reference element is assumed as follows:
     *
     *   ^ s
     *   |
     *  2 (1,0)
     *   o
     *   . .
     *   .   .
     *   .     .
     *   .       .
     *   .         .
     *   .           .
     *   o-------------o---------> r
     *  0 (0,0)         1 (1,0)
     * direct local numbering : 0->1->2
     */
  /*---------------------------------------------------------------------------*/
  static inline Real tri3ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 3);
#endif
    Real r = ref_coord[0];
    Real s = ref_coord[1];
    if (!inod)
      return (1 - r - s);
    if (inod == 1)
      return r;
    return s;
  }

  static inline Real3 tri3ShapeFuncDeriv(Integer inod, Real3)
  {
    if (!inod)
      return { -1., -1., 0. };
    if (inod == 1)
      return { 1., 0., 0. };
    return { 0., 1., 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
     * @brief Provides methods for reference quadratic (P2) triangle finite-element
     * The "Tri6" reference element is assumed as follows:
     *   ^ s
     *   |
     *  2 (1,0)
     *   o
     *   .  .
     *   .    .
     *   .      .
     *   o 6      o 5(0.5;0.5)
     *   .(0;0.5)   .
     *   .            .
     *   .              .
     *   o-------o-------o---------> r
     * 0(0,0)  4(0.5;0)  1(1,0)
     * direct local numbering : 0->1->2->3->4->5
     */
  /*---------------------------------------------------------------------------*/
  static inline Real tri6ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 6);
#endif
    auto wi = 0., ri = ref_coord[0], si = ref_coord[1];
    auto ri2 = 2. * ri - 1.;
    auto si2 = 2. * si - 1.;
    auto ti = 1. - ri - si, ti2 = 2. * ti - 1.;

    switch (inod) {
    default:
      break;
    case 0:
      wi = ti * ti2;
      break;
    case 1:
      wi = ri * ri2;
      break;
    case 2:
      wi = si * si2;
      break;
    case 3:
      wi = 4. * ri * ti;
      break;
    case 4:
      wi = 4. * ri * si;
      break;
    case 5:
      wi = 4. * si * ti;
      break;
    }
    return wi;
  }

  static inline Real3 tri6ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    auto ri = ref_coord[0], si = ref_coord[1];
    auto ti = 1. - ri - si;

    if (!inod) {
      auto wi = -3. + 4. * (ri + si);
      return { wi, wi, 0. };
    }
    if (inod == 1)
      return { -1. + 4. * ri, 0., 0. };
    if (inod == 2)
      return { 0., -1. + 4. * si, 0. };

    if (inod == 3)
      return { 4. * (ti - ri), -4. * ri, 0. };
    if (inod == 4)
      return { 4. * si, 4. * ri, 0. };
    return { -4. * si, 4. * (ti - si), 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) quadrangle finite-element
   * The "Quad4" reference element is assumed as follows:
   *         ^y
   *          |
   *  1 o-----1-----o 0
   *    |     |     |
   *    |     |     |
   *    |     |     |
   *   -1 ----|---- 1 ---> x
   *    |     |     |
   *    |     |     |
   *    |     |     |
   *  2 o--- -1 ----o 3
   *
   * direct local numbering : 0->1->2->3
   */
  /*---------------------------------------------------------------------------*/
  static inline Real quad4ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 4);
#endif

    auto r{ ref_coord[0] }, s{ ref_coord[1] };
    auto ri{ 1. }, si{ 1. };

    switch (inod) {
    default:
      break; // default is first node (index 0)
    case 2:
      si = -1;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 3:
      si = -1;
      break;
    }
    return ((1 + ri * r) * (1 + si * s) / 4.);
  }

  static inline Real3 quad4ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 4);
#endif

    auto r{ ref_coord[0] }, s{ ref_coord[1] };
    auto ri{ 1. }, si{ 1. }; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

    switch (inod) {
    default:
      break; // default is first node (index 0)
    case 2:
      si = -1;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 3:
      si = -1;
      break;
    }
    return { 0.25 * ri * (1 + si * s), 0.25 * si * (1 + ri * r), 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference quadratic (P2) quadrangle finite-element
   * The "Quad8" reference element is assumed as follows:
   *         ^y
   *          |
   *          4
   *  1 o-----o-----o 0
   *    |     |     |
   *    |     |     |
   *    |     |     |
   *  5 o ----|---- o 7 ---> x
   *    |     |     |
   *    |     |     |
   *    |     |     |
   *  2 o-----o-----o 3
   *          6
   * Normalized coordinates (x, y) vary between -1/+1
   * Nodes 4, 6 are on line (x = 0)
   * direct local numbering :  0->1->2->...->5->6->7
   */
  /*---------------------------------------------------------------------------*/
  static inline Real quad8ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 8);
#endif
    Real tol{ 1.0e-15 };

    auto r{ ref_coord[0] }, s{ ref_coord[1] };
    auto ri{ 1. }, si{ 1. };

    switch (inod) {
    default:
      break; // default is first node (index 0)
    case 2:
      si = -1;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 3:
      si = -1;
      break;

    case 6:
      si = -1;
      [[fallthrough]];
    case 4:
      ri = 0;
      break;

    case 5:
      ri = -1;
      [[fallthrough]];
    case 7:
      si = 0;
      break;
    }

    auto r0{ r * ri }, s0{ s * si };
    Real Phi{ 0. };
    auto t0{ r0 + s0 - 1. };

    if (inod < 4) // Corner nodes
      Phi = (1 + r0) * (1 + s0) * t0 / 4.;

    else { // Middle nodes
      if (fabs(ri) < tol)
        Phi = (1 - r * r) * (1 + s0) / 2.;
      else if (fabs(si) < tol)
        Phi = (1 - s * s) * (1 + r0) / 2.;
    }
    return Phi;
  }

  static inline Real3 quad8ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    Real tol{ 1.0e-15 };

    auto r{ ref_coord[0] }, s{ ref_coord[1] };
    auto ri{ 1. }, si{ 1. };

    switch (inod) {
    default:
      break; // default is first node (index 0)
    case 2:
      si = -1;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 3:
      si = -1;
      break;

    case 6:
      si = -1;
      [[fallthrough]];
    case 4:
      ri = 0;
      break;

    case 5:
      ri = -1;
      [[fallthrough]];
    case 7:
      si = 0;
      break;
    }

    auto r0{ r * ri }, s0{ s * si };
    Real3 dPhi;
    auto t0{ r0 + s0 - 1. };

    if (inod < 4) { // Corner nodes
      dPhi.x = ri * (1 + s0) * (t0 + 1. + r0) / 4.;
      dPhi.y = si * (1 + r0) * (t0 + 1. + s0) / 4.;
    }
    else { // Middle nodes
      if (fabs(ri) < tol) {
        dPhi.x = -r * (1 + s0);
        dPhi.y = si * (1 - r * r) / 2.;
      }
      else if (fabs(si) < tol) {
        dPhi.x = -s * (1 + r0);
        dPhi.y = ri * (1 - s * s) / 2.;
      }
    }
    dPhi.z = 0.;
    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) hexaedron finite-element
   * The "Hexa8" reference element is assumed as follows:
   *     (-1, 1,1)
   *         1-------------0 (1,1,1)
   *        /|            /|
   *       / |           / |
   *     /   |          /  |
   *    2----|---------3   |   z   y
   *  (-1,-1,1)        |   |   | /
   *    |    |         |   |   |/--->x
   *    |    |         |   |
   *    |    |         |   |
   *    |    5---------|---4 (1,1,-1)
   *    |  /           |  /
   *    | /            | /
   *    |/             |/
   *    6--------------7 (1,-1,-1)
   * (-1,-1,-1)
   * Normalized coordinates (x, y, z) vary between -1/+1
   * direct local numbering : 0->1->2->3->4->5->6->7
   */
  /*---------------------------------------------------------------------------*/
  static inline Real hexa8ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 8);
#endif
    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

    switch (inod) {
    default:
      break;
    case 3:
    case 2:
      ri = -1;
      break;
    case 0:
    case 1:
      ri = -1;
      si = -1;
      break;
    case 4:
    case 5:
      si = -1;

      break;
    }
    if (inod == 1 || inod == 2 || inod == 5 || inod == 6)
      ti = -1;

    auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
    auto Phi = (1 + r0) * (1 + s0) * (1 + t0) / 8.;

    return Phi;
  }

  static inline Real3 hexa8ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {

    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

    switch (inod) {
    default:
      break;
    case 3:
    case 2:
      ri = -1;
      break;
    case 0:
    case 1:
      ri = -1;
      si = -1;
      break;
    case 4:
    case 5:
      si = -1;
      break;
    }
    if (inod == 1 || inod == 2 || inod == 5 || inod == 6)
      ti = -1;

    auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
    Real3 dPhi;
    dPhi.x = ri * (1 + s0) * (1 + t0) / 8.;
    dPhi.y = si * (1 + r0) * (1 + t0) / 8.;
    dPhi.z = ti * (1 + r0) * (1 + s0) / 8.;
    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference quadratic (P2) hexaedron finite-element
   * The "Hexa20" reference element is assumed as follows:
   *     (-1, 1,1)
   *         1------8------0 (1,1,1)
   *        /|            /|
   *      9  |          11 |
   *     /   |          /  |
   *    2----|--10-----3   |   z   y
   *  (-1,-1,1)        |   |   | /
   *    |   17         |  16   |/--->x
   *    |    |         |   |
   *   18    |        19   |
   *    |    5----12---|---4 (1,1,-1)
   *    |  /           |  /
   *    | 13           | 15
   *    |/             |/
   *    6-----14-------7 (1,-1,-1)
   * (-1,-1,-1)
   * Normalized coordinates (x, y, z) vary between -1/+1
   * direct local numbering : 0->1->2->3->...->18->19
   */
  /*---------------------------------------------------------------------------*/

  static inline Real hexa20ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 20);
#endif
    Real tol{ 1.0e-15 };

    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 0 = (1,1,1)

    switch (inod) {
    default:
      break;

    case 5:
      ti = -1.;
      [[fallthrough]];
    case 1:
      ri = -1;
      break;

    case 6:
      ti = -1.;
      [[fallthrough]];
    case 2:
      ri = -1;
      si = -1;
      break;

    case 7:
      ti = -1.;
      [[fallthrough]];
    case 3:
      si = -1;
      break;

    case 4:
      ti = -1.;
      break;

    case 9:
      ri = -1.;
      [[fallthrough]];
    case 11:
      si = 0.;
      break;

    case 10:
      si = -1.;
      [[fallthrough]];
    case 8:
      ri = 0.;
      break;

    case 14:
      si = -1.;
      [[fallthrough]];
    case 12:
      ri = 0.;
      ti = -1.;
      break;

    case 17:
      ri = -1.;
      [[fallthrough]];
    case 16:
      ti = 0.;
      break;

    case 18:
      ri = -1.;
      [[fallthrough]];
    case 19:
      si = -1.;
      ti = 0.;
      break;
    }

    auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
    Real Phi{ 0. };
    auto t{ r0 + s0 + t0 - 2. };

    if (inod < 8) // Corner nodes
      Phi = (1 + r0) * (1 + s0) * (1 + t0) * t / 8.;

    else { // Middle nodes
      if (math::abs(ri) < tol)
        Phi = (1 - x * x) * (1 + s0) * (1 + t0) / 4.;
      else if (math::abs(si) < tol)
        Phi = (1 - y * y) * (1 + r0) * (1 + t0) / 4.;
      else if (math::abs(ti) < tol)
        Phi = (1 - z * z) * (1 + r0) * (1 + s0) / 4.;
    }
    return Phi;
  }

  static inline Real3 hexa20ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    Real tol{ 1.0e-15 };

    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 0 = (1,1,1)

    switch (inod) {
    default:
      break;

    case 5:
      ti = -1.;
    case 1:
      ri = -1;
      break;

    case 6:
      ti = -1.;
    case 2:
      ri = -1;
      si = -1;
      break;

    case 7:
      ti = -1.;
    case 3:
      si = -1;
      break;

    case 4:
      ti = -1.;
      break;

    case 9:
      ri = -1.;
    case 11:
      si = 0.;
      break;

    case 10:
      si = -1.;
    case 8:
      ri = 0.;
      break;

    case 14:
      si = -1.;
    case 12:
      ri = 0.;
      ti = -1.;
      break;

    case 17:
      ri = -1.;
      [[fallthrough]];
    case 16:
      ti = 0.;
      break;

    case 18:
      ri = -1.;
      [[fallthrough]];
    case 19:
      si = -1.;
      ti = 0.;
      break;
    }

    auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
    auto t{ r0 + s0 + t0 - 2. };
    Real3 dPhi;

    if (inod < 8) { // Corner nodes
      dPhi = hexa8ShapeFuncDeriv(inod, ref_coord);
      dPhi.x *= (t + 1. + r0);
      dPhi.y *= (t + 1. + s0);
      dPhi.z *= (t + 1. + t0);
    }
    else { // Middle nodes
      auto x2{ x * x }, y2{ y * y }, z2{ z * z };
      if (math::abs(ri) < tol) {
        dPhi.x = -x * (1 + s0) * (1 + t0) / 2.;
        dPhi.y = si * (1 - x2) * (1 + t0) / 4.;
        dPhi.z = ti * (1 - x2) * (1 + s0) / 4.;
      }
      else if (math::abs(si) < tol) {
        dPhi.x = ri * (1 - y2) * (1 + t0) / 4.;
        dPhi.y = -y * (1 + r0) * (1 + t0) / 2.;
        dPhi.z = ti * (1 - y2) * (1 + r0) / 4.;
      }
      else if (math::abs(ti) < tol) {
        dPhi.x = ri * (1 - z2) * (1 + s0) / 4.;
        dPhi.y = si * (1 - z2) * (1 + r0) / 4.;
        dPhi.z = -z * (1 + r0) * (1 + s0) / 2.;
      }
    }
    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) tetrahedral finite-element
   * The "Tetra4" reference element is assumed as follows:
   *
   *    (0,0,1)                     3
   *       .                        *.*
   *       .                        * . *
   *       .                        *  .  *
   *       .                        *   .   *
   *       Z   (0,1,0)              *    .    *
   *       .    .                   *     2     *
   *       .   .                    *   .    .    *
   *       .  Y                     *  .        .   *
   *       . .                      * .            .  *
   *       ..           (1,0,0)     *.                . *
   *       --------X------>         0********************1
   *
   * direct local numbering : 0->1->2->3
   */
  /*---------------------------------------------------------------------------*/

  static inline Real tetra4ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 4);
#endif

    auto ri = ref_coord[0], si = ref_coord[1], ti = ref_coord[2]; // default is first node (index 3)

    switch (inod) {
    default:
      break;
    case 1:
      return ri;
    case 2:
      return si;
    case 0:
      return (1. - ri - si - ti);
    }
    return ti;
  }

  static inline Real3 tetra4ShapeFuncDeriv(Integer inod, Real3 /*ref_coord*/)
  {

    if (inod == 3)
      return { 0., 0., 1. };
    if (inod == 1)
      return { 1., 0., 0. };
    if (inod == 2)
      return { 0., 1., 0. };
    return { -1., -1., -1. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference quadratic (P2) tetrahedral finite-element
   * The "Tetra10" reference element is assumed as follows:
   *
   *    (0,0,1)                     x 3
   *       .                        *.*
   *       .                        * . *
   *       .                        *  .  *
   *       .                        *  9   *
   *       Z   (0,1,0)              *    .    *
   *       .    .                   *     x     8
   *       .   .                    7   . 2  .    *
   *       .  Y                     *  6       5    *
   *       . .                      * .            .  *
   *       ..           (1,0,0)     *.                . *
   *       --------X------>       0 x ******* 4 ******** x 1
   *
   * direct local numbering : 0->1->2->3->...->8->9
   */
  /*---------------------------------------------------------------------------*/

  static inline Real tetra10ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 10);
#endif

    auto x = ref_coord[0], y = ref_coord[1], z = ref_coord[2],
         t = 1. - x - y - z,
         wi{ 0. };

    switch (inod) {
    default:
      break;

    // Corner nodes
    case 0:
      wi = t * (2 * t - 1.);
      break; //=(1. - 2*x - 2*y - 2*z) * t
    case 1:
      wi = x * (2 * x - 1.);
      break; //=(1. - 2*t - 2*y - 2*z)*x
    case 2:
      wi = y * (2 * y - 1.);
      break; //=(1. - 2*x - 2*t - 2*z)*y
    case 3:
      wi = z * (2 * z - 1.);
      break; //=(1. - 2*t - 2*x - 2*y)*z

    // Middle nodes
    case 4:
      wi = 4 * x * t;
      break;
    case 5:
      wi = 4 * x * y;
      break;
    case 6:
      wi = 4 * y * t;
      break;
    case 7:
      wi = 4 * z * t;
      break;
    case 8:
      wi = 4 * z * x;
      break;
    case 9:
      wi = 4 * z * y;
      break;
    }
    return wi;
  }

  static inline Real3 tetra10ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {
    auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] },
    t{ 1. - x - y - z },
    x4{ 4 * x },
    y4{ 4 * y },
    z4{ 4 * z },
    t4{ 4 * t };

    // Corner nodes
    /*
   if (inod == 3) return {0.,0.,1. + 2*t - 2*x - 2*y + 2*z};
   if (inod == 1) return {1. - 2*t - 2*y - 2*z + 2*x,0.,0.};
   if (inod == 2) return {0.,1. - 2*x - 2*t - 2*z + 2*y,0.};
   if (!inod) return {-1. - 2*t + 2*x + 2*y + 2*z,-1. - 2*t+ 2*x + 2*y + 2*z,-1. - 2*t + 2*x + 2*y + 2*z};
*/
    if (!inod)
      return { 1. - t4, 1. - t4, 1. - t4 };
    if (inod == 1)
      return { x4 - 1., 0., 0. };
    if (inod == 2)
      return { 0., y4 - 1., 0. };
    if (inod == 3)
      return { 0., 0., z4 - 1. };

    // Middle nodes
    if (inod == 4)
      return { t4 - x4, -x4, -x4 };
    if (inod == 5)
      return { y4, x4, 0. };
    if (inod == 6)
      return { -y4, t4 - y4, -y4 };
    if (inod == 8)
      return { z4, 0., x4 };
    if (inod == 9)
      return { 0., z4, y4 };
    return { -z4, -z4, t4 - z4 }; //inod == 7
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) pentaedron finite-element
   * The "Penta6" reference element is assumed as follows:
   *
   *                     5 (0,1,1)
   *                   . |  .
   *                  .  |     .
   *                 .   Z        .
   *                .    |           .
   *               .     |             .
   *       (0,0,1) 3 ------------------ 4 (1,0,1)
   *               |     |              |
   *               |     |              |
   *               |     |              |
   *               |     |              |
   *               |     2 (0,1,-1)     |
   *               |   .    .           |
   *               |  Y        .        |
   *               | .            .     |
   *               |.                .  |
   *      (0,0,-1) 0 -------- X ------- 1 (1,0,-1)
   *
   * direct local numbering : 0->1->2->3->4->5->6
   */
  /*---------------------------------------------------------------------------*/

  static inline Real penta6ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 6);
#endif
    auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
    auto r0{ 1. }, s0{ 1. }, ti{ -1. };
    auto rs{ 1. - r - s };

    if (inod >= 3)
      ti = 1.;
    auto t0{ 1 + ti * t };

    switch (inod) {
    default:
      break; // Node 0
    case 4:
    case 1:
      r0 = r;
      rs = 1.;
      break;
    case 5:
    case 2:
      s0 = s;
      rs = 1.;
      break;
    }

    return 0.5 * r0 * s0 * rs * t0;
  }

  static inline Real3 penta6ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {

#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 6);
#endif
    auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
    auto ri{ 1. }, si{ 1. };
    auto r0{ 1. }, s0{ 1. }, ti{ -1. };
    auto rs{ 1. - r - s };

    if (inod >= 3)
      ti = 1.;
    auto t0{ 1 + ti * t };

    switch (inod) {
    default:
      break;
    case 3:
    case 0:
      ri = -1.;
      si = -1.;
      break;
    case 4:
    case 1:
      r0 = r;
      si = 0.;
      rs = 1.;
      break;
    case 5:
    case 2:
      s0 = s;
      rs = 1.;
      break;
    }

    Real3 dPhi;
    dPhi.x = 0.5 * ri * t0;
    dPhi.y = 0.5 * si * t0;
    dPhi.z = 0.5 * ti * rs * r0 * s0;
    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for reference linear (P1) pyramid finite-element
   * The "Pyramid5" reference element is assumed as follows:
   *
   *                               ^
   *                               |
   *                               Z
   *                               |
   *                               4 (0,0,1)
   *                              *
   *                             * **
   *                            ** *  *
   *                           * * |*   *
   *                          * *  | *    *
   *                         * *   |  *     *
   *                        * *    |   *      *      .Y
   *                       * *     |    *       *  .
   *            (-1,0,0)  * 2 -----|-----*------ 1 (0,1,0)
   *                     * .       |      *  .  .
   *                    * .        |     . *   .
   *                   *.             .     * .
   *                  *                  X . *
   *        (0,-1,0) 3 --------------------- 0 (1,0,0)
   *
   * direct local numbering : 0->1->2->3->4
   */
  /*---------------------------------------------------------------------------*/

  static inline Real pyramid5ShapeFuncVal(Integer inod, Real3 ref_coord)
  {
#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 5);
#endif
    Real tol{ 1.0e-15 };

    auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
    auto r1{ -1. }, s1{ 1. }, r2{ -1. }, s2{ -1. };

    if (inod == 4)
      return t;
    auto ti{ t - 1. };
    auto t0{ 0. };

    if (math::abs(ti) < tol)
      ti = 0.;
    else
      t0 = -1. / ti / 4.;

    switch (inod) {
    case 1:
      s1 = -1.;
      r2 = 1.;
      break;
    case 2:
      r1 = 1.;
      r2 = 1.;
      break;
    case 3:
      r1 = 1.;
      s2 = 1.;
      break;
    default:
      break; // default is for node 0
    }

    return (r1 * r + s1 * s + ti) * (r2 * r + s2 * s + ti) * t0;
  }

  static inline Real3 pyramid5ShapeFuncDeriv(Integer inod, Real3 ref_coord)
  {

#ifdef _DEBUG
    ARCANE_ASSERT(inod >= 0 && inod < 5);
#endif
    Real tol{ 1.0e-15 };

    auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
    auto r1{ -1. }, s1{ 1. }, r2{ -1. }, s2{ -1. };

    auto ti{ t - 1. };
    auto t0{ 0. };

    if (math::abs(ti) < tol)
      ti = 0.;
    else
      t0 = -1. / ti / 4.;

    switch (inod) {
    case 1:
      s1 = -1.;
      r2 = 1.;
      break;
    case 2:
      r1 = 1.;
      r2 = 1.;
      break;
    case 3:
      r1 = 1.;
      s2 = 1.;
      break;
    default:
      break; // default is for node 0
    }

    if (inod == 4)
      return { 0., 0., 1. };

    Real3 dPhi;
    auto r12{ r1 + r2 }, rr{ 2. * r1 * r2 }, s12{ s1 + s2 }, ss{ 2. * s1 * s2 }, rs{ r1 * s2 + r2 * s1 }, t02{ 4. * t0 * t0 };

    dPhi.x = t0 * (rr * r + rs * s + r12 * ti);
    dPhi.y = t0 * (rs * r + ss * s + s12 * ti);

    if (math::abs(ti) < tol)
      dPhi.z = 0.;
    else
      dPhi.z = t0 * (r12 * r + s12 * s + 2. * ti) + t02 * (r1 * r + s1 * s + ti) * (r2 * r + s2 * s + ti);

    return dPhi;
  }

  /*---------------------------------------------------------------------------*/
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides methods for Gauss quadrature.
 *
 * This class includes static methods for computing Gauss-Legendre integration
 * depending on finite element types (coming from PASSMO).
 */
/*---------------------------------------------------------------------------*/
class FemGaussQuadrature
{

 public:

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the number of Gauss Points for a given finite element type,
   * depending on the integration order chosen by user (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  static inline Integer getNbGaussPointsfromOrder(Int16 cell_type, Integer ninteg)
  {
    Integer nbgauss{ 0 };
    auto ninteg2{ ninteg * ninteg };
    auto ninteg3{ ninteg2 * ninteg };

    if (ninteg <= 1)
      nbgauss = 1;
    else if (cell_type == IT_Line2 || cell_type == IT_Line3)
      nbgauss = ninteg;
    else if (cell_type == IT_Quad4 || cell_type == IT_Quad8)
      nbgauss = ninteg2;
    else if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20)
      nbgauss = ninteg3;

    else if (ninteg == 2) {
      switch (cell_type) {
      default:
        break;

      case IT_Triangle3:
      case IT_Triangle6:
        nbgauss = 3;
        break;

      case IT_Tetraedron4:
      case IT_Tetraedron10:
        nbgauss = 4;
        break;

      case IT_Pentaedron6:
        nbgauss = 6;
        break;

      case IT_Pyramid5:
        nbgauss = 5;
        break;
      }
    }
    else if (ninteg == 3) {
      switch (cell_type) {
      default:
        break;

      case IT_Triangle3:
      case IT_Triangle6:
        nbgauss = 4;
        break;

      case IT_Tetraedron4:
      case IT_Tetraedron10:
        nbgauss = 5;
        break;

      case IT_Pentaedron6:
        nbgauss = 8;
        break;

      case IT_Pyramid5:
        nbgauss = 6;
        break;
      }
    }
    else if (ninteg >= 4) {
      switch (cell_type) {
      default:
        break;

      case IT_Triangle3:
      case IT_Triangle6:
        nbgauss = 7;
        break;

      case IT_Tetraedron4:
      case IT_Tetraedron10:
        nbgauss = 15;
        break;

      case IT_Pentaedron6:
        nbgauss = 21;
        break;

      case IT_Pyramid5:
        nbgauss = 27;
        break;
      }
    }
    return nbgauss;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the position of a Gauss Point in the reference (local) element
   * for a given finite element, depending on the rank of the point in the
   * element loop (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 getGaussRefPosition(ItemWithNodes cell, Integer ninteg, Integer rank)
  {
    auto cell_type = cell.type();
    Integer nint{ ninteg };

    if (nint < 1)
      nint = 1;

    if (cell_type == IT_Line2 || cell_type == IT_Line3)
      return lineRefPosition({ rank, -1, -1 }, { nint, 0, 0 });

    if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
      auto in{ 0 };
      for (Int32 i1 = 0; i1 < nint; ++i1) {
        for (Int32 i2 = 0; i2 < nint; ++i2) {

          if (rank == in)
            return quadRefPosition({ i1, i2, -1 }, { nint, nint, 0 });

          ++in;
        }
      }
    }

    if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
      auto in{ 0 };
      for (Int32 i1 = 0; i1 < nint; ++i1) {
        for (Int32 i2 = 0; i2 < nint; ++i2) {
          for (Int32 i3 = 0; i3 < nint; ++i3) {

            if (rank == in)
              return hexaRefPosition({ i1, i2, i3 }, { nint, nint, nint });

            ++in;
          }
        }
      }
    }

    if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
      auto o{ 3 };
      if (nint <= 3)
        o = nint - 1;

      return { xg1[o][rank], xg2[o][rank], 0. };
    }

    if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10) {
      auto o{ 3 };
      if (nint <= 3)
        o = nint - 1;

      return { xtet[o][rank], ytet[o][rank], ztet[o][rank] };
    }

    if (cell_type == IT_Pyramid5) {
      auto o{ 1 };
      if (nint <= 2)
        o = 0;

      return { xpyr[o][rank], ypyr[o][rank], zpyr[o][rank] };
    }

    if (cell_type == IT_Pentaedron6) {
      auto o{ 1 };
      if (nint <= 2)
        o = 0;

      return { xpent[o][rank], ypent[o][rank], zpent[o][rank] };
    }
    return {};
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for a given
   * finite element, depending on the rank of the point in the element loop
   * (coming fro PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  static inline Real getGaussWeight(ItemWithNodes cell, Integer ninteg, Integer rank)
  {
    auto cell_type = cell.type();
    Integer nint{ ninteg };

    if (nint < 1)
      nint = 1;

    if (cell_type == IT_Line2 || cell_type == IT_Line3)
      return lineWeight({ rank, -1, -1 }, { nint, 0, 0 });

    if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
      auto in{ 0 };
      for (Int32 i1 = 0; i1 < nint; ++i1) {
        for (Int32 i2 = 0; i2 < nint; ++i2) {

          if (rank == in)
            return quadWeight({ i1, i2, -1 }, { nint, nint, 0 });

          ++in;
        }
      }
    }

    if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
      auto in{ 0 };
      for (Int32 i1 = 0; i1 < nint; ++i1) {
        for (Int32 i2 = 0; i2 < nint; ++i2) {
          for (Int32 i3 = 0; i3 < nint; ++i3) {

            if (rank == in)
              return hexaWeight({ i1, i2, i3 }, { nint, nint, nint });

            ++in;
          }
        }
      }
    }

    if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
      auto o{ 3 };
      if (nint <= 3)
        o = nint - 1;

      return wg[o][rank];
    }

    if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10) {
      if (nint == 1)
        return wgtet1;
      if (nint == 2)
        return wgtet2;
      if (nint == 3) {
        auto i = npwgtet3[rank];
        return wgtet3[i];
      }
      // nint >= 4
      auto i = npwgtet4[rank];
      return wgtet4[i];
    }

    if (cell_type == IT_Pyramid5) {
      if (nint <= 2)
        return wgpyr2;

      // nint >= 3
      auto i = npwgpyr3[rank];
      return wgpyr3[i];
    }

    if (cell_type == IT_Pentaedron6) {
      if (nint <= 2)
        return wgpent2;

      // nint >= 3
      auto i = npwgpent3[rank];
      return wgpent3[i];
    }

    return 1.;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the position of a Gauss Point in the reference (local) element
   * coordinates, depending on the index of the Gauss Point and integration order
   * chosen by user (coming fro PASSMO).
   * This method is generic (not depending on the finite element type) and is called
   * by specialized methods (depending on FE type)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real getRefPosition(Integer indx, Integer ordre)
  {
    Real x = xgauss1; // default is order 1

    switch (ordre) {
    case 2:
      x = xgauss2[indx];
      break;
    case 3:
      x = xgauss3[indx];
      break;
    case 4:
      x = xgauss4[indx];
      break;
    case 5:
      x = xgauss5[indx];
      break;
    case 6:
      x = xgauss6[indx];
      break;
    case 7:
      x = xgauss7[indx];
      break;
    case 8:
      x = xgauss8[indx];
      break;
    case 9:
      x = xgauss9[indx];
      break;
    default:
      break;
    }
    return x;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight depending on the index of the Gauss Point
   * and integration order chosen by user (coming fro PASSMO).
   * This method is generic (not depending on the finite element type) and is called
   * by specialized methods (depending on FE type)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real getWeight(Integer indx, Integer ordre)
  {
    Real w = wgauss1; // default is order 1

    switch (ordre) {
    case 2:
      w = wgauss2[indx];
      break;
    case 3:
      w = wgauss3[indx];
      break;
    case 4:
      w = wgauss4[indx];
      break;
    case 5:
      w = wgauss5[indx];
      break;
    case 6:
      w = wgauss6[indx];
      break;
    case 7:
      w = wgauss7[indx];
      break;
    case 8:
      w = wgauss8[indx];
      break;
    case 9:
      w = wgauss9[indx];
      break;
    default:
      break;
    }
    return w;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for edge elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/

  static inline Real3 lineRefPosition(Integer3 indices, Integer3 ordre)
  {
    return { getRefPosition(indices[0], ordre[0]), 0., 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for triangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 triRefPosition(Integer3 indices, Integer3 ordre)
  {
    Integer o = ordre[0] - 1;
    Integer i = indices[0];
    return { xg1[o][i], xg2[o][i], 0. };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for quadrangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 quadRefPosition(Integer3 indices, Integer3 ordre)
  {
    Real3 pos;
    for (Integer i = 0; i < 2; i++)
      pos[i] = getRefPosition(indices[i], ordre[i]);
    return pos;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for hexaedron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real3 hexaRefPosition(Integer3 indices, Integer3 ordre)
  {
    Real3 pos;
    for (Integer i = 0; i < 3; i++)
      pos[i] = getRefPosition(indices[i], ordre[i]);
    return pos;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for tetraedron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  [[maybe_unused]] static inline Real3 tetraRefPosition(Integer3 indices, Integer3 /*ordre*/)
  {
    Integer i = indices[0];
    return { xit[i], yit[i], zit[i] };
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the reference coordinates of a Gauss Point for wedge (pentaedron)
   * elements (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  [[maybe_unused]] static inline Real3 pentaRefPosition(Integer3 indices, Integer3 ordre)
  {
    // Same as TriRefPosition on reference coordinate plane (r,s)
    // and LineRefPosition along reference coordinate t (vertical)
    auto pos = triRefPosition(indices, ordre);
    pos.z = getRefPosition(indices[2], ordre[2]);

    return pos;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for edge elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real lineWeight(Integer3 indices, Integer3 ordre)
  {
    return getWeight(indices[0], ordre[0]);
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for triangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real triWeight(Integer3 indices, Integer3 ordre)
  {
    return wg[ordre[0] - 1][indices[0]];
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for quadrangle elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real quadWeight(Integer3 indices, Integer3 ordre)
  {
    Real w = 1.;
    for (Integer i = 0; i < 2; i++)
      w *= getWeight(indices[i], ordre[i]);
    return w;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for hexadron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  static inline Real hexaWeight(Integer3 indices, Integer3 ordre)
  {
    Real w = 1.;
    for (Integer i = 0; i < 3; i++)
      w *= getWeight(indices[i], ordre[i]);
    return w;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for tetraedron elements
   * (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  [[maybe_unused]] static inline Real tetraWeight(Integer3 /*indices*/, Integer3 /*ordre*/)
  {
    return wgtetra;
  }

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides the integration weight of a Gauss Point for wedge (pentaedron)
   * elements (coming fro PASSMO).
   * This method takes the indices and integration orders chosen
   * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
   */
  /*---------------------------------------------------------------------------*/
  [[maybe_unused]] static inline Real pentaWeight(Integer3 indices, Integer3 ordre)
  {

    // Same as TriWeight on reference coordinate plane (r,s)
    // and LineWeight with ordre[2] to account for reference coordinate t (vertical)
    Real wgpenta = triWeight(indices, ordre) * getWeight(indices[2], ordre[2]);
    return wgpenta;
  }

  /*---------------------------------------------------------------------------*/
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace ArcaneFemFunctions

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif // ARCANE_FEM_FUNCTIONS_H
