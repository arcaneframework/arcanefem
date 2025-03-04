// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.h                                                (C) 2022-2025 */
/*                                                                           */
/* FemModule class definition.                                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef FEMMODULES_H
#define FEMMODULES_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/CaseTable.h>

#include "IDoFLinearSystemFactory.h"
#include "Fem_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "ArcaneFemFunctions.h"


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Module Fem.
 */
class FemModule
: public ArcaneFemObject
{
 public:

  explicit FemModule(const ModuleBuildInfo& mbi)
  : ArcaneFemObject(mbi)
  , m_dofs_on_nodes(mbi.subDomain()->traceMng())
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
  }
  ~FemModule()
  {
    for( const CaseTableInfo&  t : m_traction_case_table_list )
      delete t.case_table;
  }

 void startInit() override; //! Method called at the beginning of the simulation
 void compute() override; //! Method called at each iteration
 VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 private:

  Real t; // time variable 𝑡
  Real dt; // time step δ𝑡
  Real tmax; // max time 𝑡ₘₐₓ
  Real alpm; // time discretization param αᵐ
  Real alpf; // time discretization param αᶠ
  Real beta; // time discretization param β
  Real gamma; // time discretization param γ

  Real etam; // damping parameter ηₘ
  Real etak; // damping parameter ηₖ

  Real E; // Youngs modulus 𝐸
  Real nu; // Poissons ratio ν
  Real rho; // Density ρ
  Real mu; // Lame parameter μ
  Real lambda; // Lame parameter λ

  Real3 f = {0,0,0}; // body force 𝐟

  Real c0; // constant c₀
  Real c1; // constant c₁
  Real c2; // constant c₂
  Real c3; // constant c₃
  Real c4; // constant c₄
  Real c5; // constant c₅
  Real c6; // constant c₆
  Real c7; // constant c₇
  Real c8; // constant c₈
  Real c9; // constant c₉
  Real c10; // constant c₁₀

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;

  // Struct to make sure we are using a CaseTable associated
  // to the right file
  struct CaseTableInfo
  {
    String file_name;
    CaseTable* case_table = nullptr;
  };
  // List of CaseTable for traction boundary conditions
  UniqueArray<CaseTableInfo> m_traction_case_table_list;

 private:

  void _doStationarySolve();
  void _getParameters();
  void _updateVariables();
  void _updateTime();
  void _assembleBilinearOperatorTRIA3();
  void _solve();
  void _assembleLinearOperator();
  FixedMatrix<6, 6> _computeElementMatrixTRIA3(Cell cell);
  Real2 _computeDxDyOfRealTRIA3(Cell cell);
  void _checkResultFile();
  void _readCaseTables();
};

#endif