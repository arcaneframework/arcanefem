// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2024 */
/*                                                                           */
/* FemModule class definition.                                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef FEMMODULES_H
#define FEMMODULES_H

#include <arcane/utils/NumArray.h>
#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/StringList.h>

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>

#include "IArcaneFemBC.h"
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

  Real rho;
  Real epsilon;
  Real ElementNodes;

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnNodes m_dofs_on_nodes;

 private:

  void _doStationarySolve();
  void _getMaterialParameters();
  void _assembleBilinearOperator();
  void _solve();
  void _getE();
  void _assembleLinearOperator();
  void _validateResults();

  RealMatrix<3, 3> _computeElementMatrixTria3(Cell cell);

  template<int N>
  void _assembleBilinear( const std::function<RealMatrix<N, N>(const Cell&)>& compute_element_matrix);
};

#endif