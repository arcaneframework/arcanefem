﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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

#include <arcane/utils/PlatformUtils.h>
#include <arcane/utils/NumArray.h>
#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/StringList.h>

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>

#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/accelerator/VariableViews.h>

#include "IArcaneFemBC.h"
#include "IDoFLinearSystemFactory.h"
#include "Fem_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "ArcaneFemFunctions.h"

#include "ArcaneFemFunctionsGpu.h"
#include "BSRFormat.h"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
using namespace Arcane::FemUtils;
namespace ax = Arcane::Accelerator;

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
  , m_bsr_format(mbi.subDomain()->traceMng(), *(mbi.subDomain()->acceleratorMng()->defaultQueue()), m_dofs_on_nodes)
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
  }

 public:

  void startInit() override; //! Method called at the beginning of the simulation
  void compute() override; //! Method called at each iteration
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

  void _assembleBilinearOperator();
  void _assembleLinearOperatorGpu();

  private:

  Real rho;
  Real epsilon;

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnNodes m_dofs_on_nodes;
  BSRFormat m_bsr_format;

  String m_petsc_flags;
  String m_matrix_format = "DOK";

  bool m_assemble_linear_system = true;
  bool m_solve_linear_system = true;
  bool m_cross_validation = true;

  void _doStationarySolve();
  void _getMaterialParameters();
  void _assembleLinearOperatorCpu();
  void _solve();
  void _assembleLinearOperator();
  void _updateVariables();
  void _validateResults();

  RealMatrix<3, 3> _computeElementMatrixTria3(Cell cell);
  RealMatrix<4, 4> _computeElementMatrixTetra4(Cell cell);

  template <int N>
  void _assembleBilinear(const std::function<RealMatrix<N, N>(const Cell&)>& compute_element_matrix);
};

#endif