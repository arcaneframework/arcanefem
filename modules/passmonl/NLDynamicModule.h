// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-20245CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* NLDModule.h                                                 (C) 2022-2025 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef PASSMO_NLDYNAMICMODULE_H
#define PASSMO_NLDYNAMICMODULE_H

#include "TypesNLdynamic.h"
#include "NLDynamic_axl.h"
#include "FemUtils.h"
#include "analytical_func.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "GaussDoFsOnCells.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Class to define the seismic inputmotion features
 */
/*class Inputmotion{
 public:
  Real3	m_ampli_factors{1.,1.,1.};// Amplification factors to apply on (X, Y, Z) components (default = 1. => no amplification)
  Real  m_max_frequency{10.};// Max frequency for the input signal
  CaseTable* m_acc{ nullptr};
  CaseTable* m_vel{ nullptr};
  CaseTable* m_displ{ nullptr};
  bool m_rigid_base{true};
  bool m_is_vel{ false};
  bool m_is_displ{false};
  bool m_component[3]{true,true,true};//= true when the "ith" component has an input motion curve set
  NodeGroup m_node_group{};

  Inputmotion() = default;
  ~Inputmotion() = default;
};*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Module NLDynamic.
 */

class NLDynamicModule
: public ArcaneNLDynamicObject {
 public:
  explicit NLDynamicModule(const ModuleBuildInfo &mbi);
  ~NLDynamicModule();

 public:

  //! Method called at the beginning of the simulation
  void startInit() override;

  //! Method called at each step
  void compute() override;

  VersionInfo versionInfo() const override;

 private:

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;
  GaussDoFsOnCells m_gauss_on_cells;
//  IntegPoints m_integ_points;

  // Struct to make sure we are using a CaseTable associated
  // to the right file
  struct CaseTableInfo
  {
    String file_name;
    CaseTable* case_table = nullptr;
  };
  // List of CaseTable for traction boundary conditions
  UniqueArray<CaseTableInfo> m_traction_case_table_list;

  // List of CaseTable for dirichlet boundary (surface) conditions
  UniqueArray<CaseTableInfo> m_sacc_case_table_list;
  UniqueArray<CaseTableInfo> m_sdispl_case_table_list;
  UniqueArray<CaseTableInfo> m_svel_case_table_list;
  UniqueArray<CaseTableInfo> m_sforce_case_table_list;

  // List of CaseTable for dirichlet point conditions
  UniqueArray<CaseTableInfo> m_acc_case_table_list;
  UniqueArray<CaseTableInfo> m_displ_case_table_list;
  UniqueArray<CaseTableInfo> m_vel_case_table_list;
  UniqueArray<CaseTableInfo> m_force_case_table_list;

  // List of CaseTable for input motions set on paraxial boundaries
  UniqueArray<CaseTableInfo> m_ain_case_table_list;
  UniqueArray<CaseTableInfo> m_vin_case_table_list;
  UniqueArray<CaseTableInfo> m_uin_case_table_list;

  // List of CaseTable for double couple conditions (seismic moments or loadings)
  UniqueArray<CaseTableInfo> m_dc_case_table_list;

  Integer ninteg{2};
  Int32 NDIM{2};
  ArcaneFemFunctions::CellFEMDispatcher cell_fem{};
  Real3 gravity{0.,0.,0.};
  Real penalty{1.e64};
  Real gamma{0.5};
  Real beta{0.25};
  Real alfam{0.};
  Real alfaf{0.};
  Real utol{0.001}, ftol{0.01}, etol{0.};
  Real m_norm_R0{0.};

  bool is_alfa_method{false};
  bool keep_constop{false};
  bool is_linear{true};
  bool m_converge{false};
  bool m_law_stop{ false};
  bool m_ref{false};
  bool m_deseq{false}; // indicator for initial unbalance due to stresses (set by user)
  Integer m_nb_law_param{2};//max nb law constitutive parameters at Gauss points
  Integer m_nb_law_hist_param{1};//max nb law history parameters at Gauss points
  String  m_law_param_file{};

  Real dt2{0.};
  Int32 linop_nstep{1000}, linop_nstep_counter{0}, ite_max{15};
  TypesNLDynamic::eElastType elast_type{TypesNLDynamic::NoElastPropType};
  TypesNLDynamic::eAnalysisType analysis_type{TypesNLDynamic::PlaneStrain};
  TypesNLDynamic::eIntegType integ_type{TypesNLDynamic::FemCell};
  TypesNLDynamic::eAlgoType algo_type{TypesNLDynamic::ModNewtonRaphson};
  AnalyticFunc m_inputfunc{};
 private:

  void _initDofs();
  void _startInitGauss();
  void _initGaussStep();
  void _initCells();
  void _applyInitialNodeConditions();
  void _applyInitialCellConditions();
  void _assembleLinearLHS();
  void _assembleLinearRHS();
  void _assembleNonLinRHS(bool init);
  void _doSolve();
  void _initBoundaryConditions();
  void _initDCConditions();
  void _applyDirichletBoundaryConditions();
  void _applyParaxialBoundaryConditions();
  void _applyNeumannBoundaryConditions();
  void _getParaxialContribution(VariableDoFReal& rhs_values);
  void _assembleLHSParaxialContribution();
  void _getTractionContribution(Arcane::VariableDoFReal& rhs_values);
  Real3x3 _computeJacobian(const ItemWithNodes& cell, const Int32& ig, const RealUniqueArray& vec, Real& jac);

  /*  Predict nodal dofs vector for the Newmark or Generalized-alfa time integration schemes */
  void _predictNewmark();

  /*  Update nodal dofs vector for the Newmark or Generalized-alfa time integration schemes */
  void _updateNewmark();

  void _checkResultFile();

  void _computeK(const Real& lambda, const Real& mu, const DoFLocalId& igauss, const Int32& nb_nodes, RealUniqueArray2& Ke);
  void _computeElemMass(const Real& rho, const DoFLocalId& igauss, const Int32& nb_nodes, RealUniqueArray2& Me);
  void _computeKParax(const Face& face, const Int32& ig, const RealUniqueArray& vec, const Real& jacobian,
                      RealUniqueArray2& Ke, const Real3& RhoC);

  void _iterate();
  void _check_convergence(Int32 iter);

  RealUniqueArray2 _getB(const DoFLocalId& igauss, const Int32& nb_nodes);

  void _compute_stress(bool init, bool store);
  void _stress_prediction(bool init);
  void _stress_correction();

};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif // PASSMO_NLDYNAMICMODULE_H
