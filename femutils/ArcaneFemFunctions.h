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

#include "FemUtilsGlobal.h"
#include "IArcaneFemBC.h"
#include "GaussQuadrature.h"
#include "DoFLinearSystem.h"
#include "ShapeFunctions.h"

#include "MeshOperation.h"
#include "FemOperation.h"
#include "FemShapeMethods.h"
#include "FemGaussQuadrature.h"

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
/*---------------------------------------------------------------------------*/

} // namespace ArcaneFemFunctions

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif // ARCANE_FEM_FUNCTIONS_H
