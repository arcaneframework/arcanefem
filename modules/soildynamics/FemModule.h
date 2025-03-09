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

#include <arcane/utils/PlatformUtils.h>

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
    for( const CaseTableInfo&  t : m_double_couple_case_table_list )
      delete t.case_table;
  }

 public:

  //! Method called at each iteration
  void compute() override;

  //! Method called at the beginning of the simulation
  void startInit() override;

  VersionInfo versionInfo() const override
  {
    return VersionInfo(1, 0, 0);
  }

 private:

  Real t;                     // time variable
  Real dt;                    // time step
  Real tmax;                  // max time

  Real etam;                  // time discretization param etam
  Real etak;                  // time discretization param etak
  Real alpm;                  // time discretization param alpm
  Real alpf;                  // time discretization param alpf
  Real beta;                  // time discretization param beta
  Real gamma;                 // time discretization param gamma

  Real3 f = {0,0,0}; // body force 𝐟

  Real E;                     // Youngs modulus
  Real nu;                    // Poissons ratio
  Real rho;                   // Density
  Real cp;                    // Primary wave velocity of soil
  Real cs;                    // Secondary wave velocity of soil
  Real f1;                    // Body force in x
  Real f2;                    // Body force in y
  Real mu;                    // Lame parameter mu
  Real mu2;                   // Lame parameter mu * 2
  Real lambda;                // Lame parameter lambda

  Real c0;                    // constant
  Real c1;                    // constant
  Real c2;                    // constant
  Real c3;                    // constant
  Real c4;                    // constant
  Real c5;                    // constant
  Real c6;                    // constant
  Real c7;                    // constant
  Real c8;                    // constant
  Real c9;                    // constant

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
  UniqueArray<CaseTableInfo> m_double_couple_case_table_list;

 private:

  void _doStationarySolve();
  void _getParameters();
  void _updateVariables();
  void _updateTime();
  void _assembleBilinearOperatorTria3();
  void _assembleBilinearOperatorEdge2();
  void _solve();
  void _assembleLinearOperator();
  void _applyDirichletBoundaryConditions();
  void _validateResults();
  void _readCaseTables();

  FixedMatrix<4, 4> _computeElementMatrixEdge2(Face face);
  FixedMatrix<6, 6> _computeElementMatrixTria3(Cell cell);
};

#endif