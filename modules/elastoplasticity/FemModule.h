// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.h                                                (C) 2000-2026  */
/*                                                                           */
/* FemModuleElastoplasticity class definition.                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANFEM_ELASTOPLATICITY_FEMMODULE
#define ARCANFEM_ELASTOPLATICITY_FEMMODULE
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/ParameterList.h>
#include <arcane/utils/ApplicationInfo.h>

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/core/ItemTypes.h>
#include <arccore/base/ArccoreGlobal.h>
#include "arccore/base/NotImplementedException.h"

#include "IArcaneFemBC.h"
#include "IDoFLinearSystemFactory.h"
#include "Fem_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "BSRFormat.h"

#include "ArcaneFemFunctions.h"
#include "ArcaneFemFunctionsGpu.h"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
/**
 * @brief A module for finite element method.
 *
 * This class handles the initialization and computation for finite element
 * method (FEM) simulations, providing methods to  set  up and solve linear
 * systems, assemble FEM operators, and perform result checks.
 */
/*---------------------------------------------------------------------------*/

class FemModuleElastoplasticity
: public ArcaneFemObject
{
 public:

  explicit FemModuleElastoplasticity(const ModuleBuildInfo& mbi)
  : ArcaneFemObject(mbi)
  , m_dofs_on_nodes(mbi.subDomain()->traceMng())
  , m_bsr_format(mbi.subDomain()->traceMng(), *(mbi.subDomain()->acceleratorMng()->defaultQueue()), m_dofs_on_nodes)
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
  }
  ~FemModuleElastoplasticity()
  {
    for (const CaseTableInfo& t : m_traction_case_table_list)
      delete t.case_table;
    for (const CaseTableInfo& t : m_dirichlet_case_table_list)
      delete t.case_table;
  }

  void startInit() override; //! Method called at the beginning of the simulation
  void compute() override; //! Method called at each iteration
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

  void _doStationarySolve();
  void _assembleBilinearOperatorGlobal();
  void _assembleBilinearOperatorLocal(bool elastic_assembly = false);
  void _assembleDirichletsNewtonGpu();
  void _assembleZeroRHSOnConstrainedDOFsGpu();

  inline void _applyInternalBodyForceTria3Gpu(VariableDoFReal& rhs_values, const FemDoFsOnNodes& dofs_on_nodes, const VariableNodeReal3& node_coord, IMesh* mesh, RunQueue* queue);

  inline void _updateGlobalTangentMaterialTensorVonMisesTria3Gpu();
  inline void _updateGlobalTangentMaterialTensorDruckerPragerTria3Gpu();

 private:

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;
  BSRFormat m_bsr_format;

  // List of CaseTable for traction boundary conditions
  UniqueArray<CaseTableInfo> m_traction_case_table_list;
  // List of CaseTable for Dirichlet boundary conditions
  UniqueArray<CaseTableInfo> m_dirichlet_case_table_list;
  Real t = 0.;
  Real dt = 0.;
  Real tmax = 0.;
  Real E; // Youngs modulus
  Real nu; // Poisson ratio
  Real sig0; // Yield strength
  Real cohesion; // Yield strength
  Real friction_angle; // Yield strength
  Real mu;
  Real lambda;
  Real Et; // Tangent modulus
  Real H; // Hardening modulus
  Real Qlim; // Limiting pressure
  Real bulk; // Bulk modulus
  Real dpEta; //
  Real dpC; //
  Real max_settlement; // Limiting settlement
  Real footing_width; // Footing width
  Real alg_reaction; // Algebraic reaction
  Real m_newton_atol;
  Real m_newton_rtol;
  Real m_residual_norm0 = 0.0;
  Real m_increment_norm = 0.0;
  Real m_residual_norm = 0.0;

  Real3 f;

  RealMatrix<3, 3> m_C_elas_2d;
  RealMatrix<6, 6> m_C_elas_3d;

  Int8 m_dof_per_node;
  Int8 m_nGP = 1;
  Int32 m_newton_iter;
  Int32 m_newton_max_iters;

  String m_petsc_flags;
  String m_matrix_format = "DOK";
  String m_constitutive_law = "VonMises";
  String m_gp_material_tensor_strategy = "local";
  String m_newton_converged_reason = "";

  bool m_use_gpu_functions = true;
  bool m_assemble_linear_system = true;
  bool m_solve_linear_system = true;
  bool m_solve_nonlinear_system = true;
  bool m_cross_validation = false;
  bool m_hex_quad_mesh = false;

  bool m_material_initialized = false;
  bool m_newton_solver_converged = false;
  bool m_check_with_bilinear_operator = false;

  void _updateTime();
  void _getMaterialParameters();
  void _setElasticMaterialTensorAtGPs();
  void _solveNewton();
  void _checkNewtonConvergence();
  void _incrementVariables();
  void _solve();
  void _assembleLinearOperator();
  void _validateResults();
  void _readCaseTables();
  void _updateNewtonIncrements();
  void _updateTimeVariables();
  void _initBsr();
  void _initConstitutiveLaw();

  // Von Mises Law
  inline void _restoreConvergedStateVonMises();
  inline void _commitInternalVariablesVonMises();
  inline void _updateGlobalTangentMaterialTensorVonMises();
  inline void _updateGlobalTangentMaterialTensorVonMisesTria3Cpu();
  inline RealMatrix<3,3> _updateGlobalTangentMaterialTensorVonMisesTria3CpuBase(const Cell& cell, Int8& iGP);

  // Drucker Prager Law
  inline void _restoreConvergedStateDruckerPrager();
  inline void _commitInternalVariablesDruckerPrager();
  inline void _updateGlobalTangentMaterialTensorDruckerPrager();
  inline void _updateGlobalTangentMaterialTensorDruckerPragerTria3Cpu();

  // RHS assembly helper functions
  inline void _applyInternalBodyForce(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  inline void _applyInternalBodyForceTria3Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);

  inline void _applyExternalBodyForce(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);

  inline void _applyTraction(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  static inline void _applyPressureTableToRhsTria3(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index, const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values);

  inline void _applyDirichletNewton(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  inline void _applyZeroRHSOnConstrainedDOFs(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);

  inline Real _normL2(VariableNodeReal3& u);
  inline Real _normL2(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  inline Real _normL1(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);

  RealMatrix<6, 6> _computeElementMatrixTria3(Cell cell);
  RealMatrix<12, 12> _computeElementMatrixTetra4(Cell cell);
  RealMatrix<8, 8> _computeElementMatrixQuad4(Cell cell);
  RealMatrix<24, 24> _computeElementMatrixHexa8(Cell cell);

  RealMatrix<6, 6> _computeLocalVonMisesElementMatrixTria3Cpu(Cell cell, bool elastic_assembly = false);

  IBinaryMathFunctor<Real, Real3, Real>* m_prescribed_settlement = nullptr;

  template <int N>
  void _assembleBilinearOperatorCpu(const std::function<RealMatrix<N, N>(const Cell&)>& compute_element_matrix);

};
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/

#endif