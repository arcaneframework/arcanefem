// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Elastoplasticity2Module.h                                   (C) 2000-2026 */
/*                                                                           */
/* Elastoplasticity2Module class definition.                                 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANFEM_ELASTOPLATICITY2_ELASTOPLATICITY2MODULE
#define ARCANFEM_ELASTOPLATICITY2_ELASTOPLATICITY2MODULE
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

//#include <arcane/utils/CommandLineArguments.h>
//#include <arcane/utils/ParameterList.h>
//#include <arcane/utils/ApplicationInfo.h>
//#include <arcane/utils/NumArray.h>

//#include <arcane/ITimeLoopMng.h>
//#include <arcane/IMesh.h>
//#include <arcane/IItemFamily.h>
//#include <arcane/ItemGroup.h>
//#include <arcane/accelerator/core/IAcceleratorMng.h>
//#include <arcane/accelerator/core/RunQueue.h>
//#include <arcane/core/ItemTypes.h>
//#include <arccore/base/ArccoreGlobal.h>
//#include "arccore/base/NotImplementedException.h"

#include "femutils/IArcaneFemBC.h"
#include "femutils/IDoFLinearSystemFactory.h"
// #include "femutils/FemUtils.h"
#include "femutils/DoFLinearSystem.h"
// #include "femutils/FemDoFsOnNodes.h"
#include "femutils/BSRFormat.h"

#include "modules/elastoplasticity2/Elastoplasticity2_axl.h"

//#include "femutils/ArcaneFemFunctions.h"
//#include "femutils/ArcaneFemFunctionsGpu.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::ArcaneFem
{

using namespace FemUtils;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief A module for finite element method.
 *
 * This class handles the initialization and computation for finite element
 * method (FEM) simulations, providing methods to  set  up and solve linear
 * systems, assemble FEM operators, and perform result checks.
 */
class Elastoplasticity2Module
: public ArcaneElastoplasticity2Object
{
public:

  explicit Elastoplasticity2Module(const ModuleBuildInfo& mbi);
  ~Elastoplasticity2Module() override;

  void startInit() override; //! Method called at the beginning of the simulation
  void compute() override; //! Method called at each iteration
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

  void _doStationarySolve();
  void _assembleBilinearOperatorGlobal();
  void _assembleBilinearOperatorLocalVonMises(bool elastic_assembly = false);
  void _assembleBilinearOperatorLocalDruckerPrager(bool elastic_assembly = false);
  void _assembleDirichletsNewtonGpu();
  void _assembleZeroRHSOnConstrainedDOFsGpu();

  inline void _applyInternalBodyForceTria3Gpu(VariableDoFReal& rhs_values, const FemDoFsOnNodes& dofs_on_nodes, const VariableNodeReal3& node_coord, IMesh* mesh, RunQueue* queue);

  inline void _updateGlobalTangentMaterialTensorVonMisesTria3Gpu();
  inline void _updateStressAndInVarsVonMisesTria3Gpu();
  inline void _updateGlobalTangentMaterialTensorDruckerPragerTria3Gpu();
  inline void _updateStressAndInVarsDruckerPragerTria3Gpu();


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
  Real E = 0.0; // Youngs modulus
  Real nu = 0.0; // Poisson ratio
  Real sig0 = 0.0; // Yield strength
  Real cohesion = 0.0; // Yield strength
  Real friction_angle = 0.0; // Yield strength
  Real mu = 0.0;
  Real lambda = 0.0;
  Real Et = 0.0; // Tangent modulus
  Real H = 0.0; // Hardening modulus
  Real Qlim = 0.0; // Limiting pressure
  Real bulk = 0.0; // Bulk modulus
  Real dpEta = 0.0; //
  Real dpC = 0.0; //
  Real max_settlement = 0.0; // Limiting settlement
  Real footing_width = 0.0; // Footing width
  Real alg_reaction = 0.0; // Algebraic reaction
  Real m_newton_atol = 0.0;
  Real m_newton_rtol = 0.0;
  Real m_residual_norm0 = 0.0;
  Real m_increment_norm = 0.0;
  Real m_residual_norm = 0.0;

  Real3 f;

  RealMatrix<3, 3> m_C_elas_2d;
  RealMatrix<6, 6> m_C_elas_3d;

  Int8 m_dof_per_node = 0;
  Int8 m_nGP = 1;
  Int8 m_nodes_per_cell = 0;
  Int32 m_newton_iter = 0;
  Int32 m_newton_max_iters = 0;

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
  bool m_use_rigid_body_near_null_space = false;
  bool m_hex_quad_mesh = false;

  bool m_material_initialized = false;
  bool m_newton_solver_converged = false;
  bool m_check_with_bilinear_operator = false;

  NumArray<Real, MDDim2> m_near_null_space_vectors;

  void _updateTime();
  void _getMaterialParameters();
  void _setGlobalElasticMaterialTensorAtGPs();
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
  void _buildRigidBodyNearNullSpace();

  // Von Mises Law
  void _restoreConvergedStateVonMises();
  void _commitInternalVariablesVonMises();
  void _updateGlobalTangentMaterialTensorVonMises();
  void _updateGlobalTangentMaterialTensorVonMisesTria3Cpu();
  void _updateGlobalTangentMaterialTensorVonMisesQuad4Cpu();
  void _updateGlobalTangentMaterialTensorVonMisesQuad8Cpu();
  void _updateGlobalTangentMaterialTensorVonMisesQuad9Cpu();
  void _updateStressAndInVarsVonMises();
  void _updateStressAndInVarsVonMisesQuad4Cpu();
  void _updateStressAndInVarsVonMisesQuad8Cpu();
  void _updateStressAndInVarsVonMisesQuad9Cpu();

  // Drucker Prager Law
  void _restoreConvergedStateDruckerPrager();
  void _commitInternalVariablesDruckerPrager();
  void _updateGlobalTangentMaterialTensorDruckerPrager();
  void _updateGlobalTangentMaterialTensorDruckerPragerTria3Cpu();
  void _updateStressAndInVarsDruckerPrager();

  // RHS assembly helper functions
  void _applyInternalBodyForce(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  void _applyInternalBodyForceTria3Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  void _applyInternalBodyForceQuad4Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  void _applyInternalBodyForceQuad8Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  void _applyInternalBodyForceQuad9Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);

  void _applyExternalBodyForce(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  void _applyExternalBodyForceQuad4Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  void _applyExternalBodyForceQuad8Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  void _applyExternalBodyForceQuad9Cpu(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);

  void _applyTraction(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  static void _applyPressureTableToRhsTria3(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
                                            const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);
  static void _applyPressureTableToRhsLine3(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
                                            const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list,
                                            const IndexedNodeDoFConnectivityView& node_dof,
                                            const VariableNodeReal3& node_coord,
                                            VariableDoFReal& rhs_values);
  static void _applyTractionToRhsLine3(BC::ITractionBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof,
                                       const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values);

  void _applyDirichletNewton(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  void _applyZeroRHSOnConstrainedDOFs(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);

  Real _normL2(VariableNodeReal3& u);
  Real _normL2(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  Real _normL1(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);

  inline RealMatrix<6, 6> _computeElementMatrixTria3(Cell cell);
  inline RealMatrix<12, 12> _computeElementMatrixTetra4(Cell cell);
  inline RealMatrix<8, 8> _computeElementMatrixQuad4(Cell cell);
  inline RealMatrix<16, 16> _computeElementMatrixQuad8(Cell cell);
  inline RealMatrix<18, 18> _computeElementMatrixQuad9(Cell cell);
  inline RealMatrix<24, 24> _computeElementMatrixHexa8(Cell cell);

  inline RealMatrix<6, 6> _computeLocalVonMisesElementMatrixTria3Cpu(Cell cell, bool elastic_assembly = false);
  inline RealMatrix<8, 8> _computeLocalVonMisesElementMatrixQuad4Cpu(Cell cell, bool elastic_assembly = false);
  inline RealMatrix<16, 16> _computeLocalVonMisesElementMatrixQuad8Cpu(Cell cell, bool elastic_assembly = false);
  inline RealMatrix<18, 18> _computeLocalVonMisesElementMatrixQuad9Cpu(Cell cell, bool elastic_assembly = false);
  inline RealMatrix<6, 6> _computeLocalDruckerPragerElementMatrixTria3Cpu(Cell cell, bool elastic_assembly = false);

  IBinaryMathFunctor<Real, Real3, Real>* m_prescribed_settlement = nullptr;

  template <int N>
  void _assembleBilinearOperatorCpu(const std::function<RealMatrix<N, N>(const Cell&)>& compute_element_matrix);

};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
