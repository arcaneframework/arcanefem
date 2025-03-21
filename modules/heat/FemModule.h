﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
/*---------------------------------------------------------------------------*/
/* FemModule.h                                                 (C) 2022-2025 */
/*                                                                           */
/* FemModule class definition.                                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef FEMMODULES_H
#define FEMMODULES_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/PlatformUtils.h>

#include <arcane/utils/NumArray.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>

#include "IDoFLinearSystemFactory.h"
#include "Fem_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"

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

  //! Time variables
  Real t   ,
       dt  ,
       tmax;
  //! Temperature
  Real Tinit ,
       Text  ;
  //! Material parameters
  Real lambda ,
       h      ,
       qdot   ;
  //! FEM parameter
  Real ElementNodes;

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnNodes m_dofs_on_nodes;

 private:

  void _initTime();
  void _updateTime();
  void _updateVariables();
  void _initTemperature();
  void _doStationarySolve();
  void _getParameters();
  void _updateBoundayConditions();
  void _assembleBilinearOperator();
  void _assembleBilinearOperatorTRIA3();
  void _assembleBilinearOperatorEDGE2();
  void _solve();
  void _initBoundaryconditions();
  void _assembleLinearOperator();
  void _applyDirichletBoundaryConditions();
  void _checkResultFile();

  void _printArcaneFemTime(const String label, const Real value);

  RealMatrix<2, 2> _computeElementMatrixEDGE2(Face face);
  RealMatrix<3, 3> _computeElementMatrixTRIA3(Cell cell);
  Real  _computeDxOfRealTRIA3(Cell cell);
  Real  _computeDyOfRealTRIA3(Cell cell);
  Real2 _computeDxDyOfRealTRIA3(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
  Real _computeEdgeLength2(Face face);

};

#endif