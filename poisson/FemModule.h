// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
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
/**
 * @brief A module for finite element method.
 *
 * This class handles the initialization and computation for finite element
 * method (FEM) simulations, providing methods to  set  up and solve linear
 * systems, assemble FEM operators, and perform result checks.
 */
/*---------------------------------------------------------------------------*/

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

  void startInit() override; //! Method called at the beginning of the simulation
  void compute() override; //! Method called at each iteration
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 private:

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnNodes m_dofs_on_nodes;

  Real f;

  void _doStationarySolve();
  void _getMaterialParameters();
  void _assembleBilinearOperator();
  void _solve();
  void _assembleLinearOperator();
  void _validateResults();

  FixedMatrix<3, 3> _computeElementMatrixTria3(Cell cell);
  FixedMatrix<4, 4> _computeElementMatrixTetra4(Cell cell);

  template<int N>
  void _assembleBilinear( const std::function<FixedMatrix<N, N>(const Cell&)>& compute_element_matrix);
};

#endif
