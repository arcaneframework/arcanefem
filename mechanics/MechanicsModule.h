// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* MechanicsModule.h                                           (C) 2022-2024 */
/*                                                                           */
/* Quasi-static implicit mechanical FEM solver                               */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef MECHANICSMODULE_H
#define MECHANICSMODULE_H

#include "Types.h"
#include "Mechanics_axl.h"
#include "FemUtils.h"
#include "utilFEM.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "GaussDoFsOnCells.h"


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;
Real REL_PREC{1.0e-15};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Class to define the seismic inputmotion features
 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Module Mechanics.
 */

class MechanicsModule
: public ArcaneMechanicsObject {
 public:
  explicit MechanicsModule(const ModuleBuildInfo &mbi);
  ~MechanicsModule();

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
  UniqueArray<CaseTableInfo> m_sdispl_case_table_list;
  UniqueArray<CaseTableInfo> m_sforce_case_table_list;

  // List of CaseTable for dirichlet point conditions
  UniqueArray<CaseTableInfo> m_displ_case_table_list;
  UniqueArray<CaseTableInfo> m_force_case_table_list;

  Integer ninteg{2};
  Int32 NDIM{2};
  CellFEMDispatcher cell_fem{};
  Real3 gravity{0.,0.,0.};
  Real penalty{1.e64};
  Real theta{0.5};
  bool keep_constop{false};

  Real dt2{0.};
  Int32 linop_nstep{1000}, linop_nstep_counter{0};
  TypesMechanics::eElastType elast_type{TypesMechanics::NoElastPropType};
  TypesMechanics::eAnalysisType analysis_type{TypesMechanics::PlaneStrain};

 private:

  void _initDofs();
  void _startInitGauss();
  void _initGaussStep();
  void _initCells();
  void _applyInitialNodeConditions();
  void _applyInitialCellConditions();
  void _assembleLinearLHS();
  void _assembleLinearRHS();
  void _doSolve();
  void _initBoundaryConditions();
  void _applyDirichletBoundaryConditions();
  void _applyNeumannBoundaryConditions();
  void _getTractionContribution(Arcane::VariableDoFReal& rhs_values);
  Real3x3 _computeJacobian(const ItemWithNodes& cell, const Int32& ig, const RealUniqueArray& vec, Real& jac);
  Real _computeFacLengthOrArea(const Face& face);

  /*  Update nodal dofs vector and stress/strains at end of time step (convergence) */
  void _update();

  void _computeK(const Real& lambda, const Real& mu, const DoFLocalId& igauss, const Int32& nb_nodes, RealUniqueArray2& Ke);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif // MECHANICSMODULE_H
