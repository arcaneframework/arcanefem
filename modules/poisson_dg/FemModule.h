// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.h                                                 (C) 2022-2026 */
/*                                                                           */
/* FemModulePoisson class definition.                                        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef FEMMODULES_H
#define FEMMODULES_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/PlatformUtils.h>
#include <arcane/utils/ParameterList.h>
#include <arcane/utils/ApplicationInfo.h>
#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/NumArray.h>

#include <arcane/ITimeLoopMng.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/IMesh.h>

#include <arccore/base/NotImplementedException.h>
#include <arcane/mesh/IncrementalItemConnectivity.h>

#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/accelerator/VariableViews.h>

#include "IDoFLinearSystemFactory.h"
#include "ArcaneFemFunctionsGpu.h"
#include "ArcaneFemFunctions.h"
#include "CsrFormatMatrix.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnCells.h"
#include "IArcaneFemBC.h"
#include "FemUtils.h"

#include "Fem_axl.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
using namespace Arcane::FemUtils;
namespace ax = Arcane::Accelerator;

/*---------------------------------------------------------------------------*/
/**
 * @brief A module for finite element method.
 *
 * This class handles the initialization and computation for finite element
 * method (FEM) simulations, providing methods to  set  up and solve linear
 * systems, assemble FEM operators, and perform result checks.
 */
/*---------------------------------------------------------------------------*/

class FemModulePoisson
: public ArcaneFemObject
{
 public:

  explicit FemModulePoisson(const ModuleBuildInfo& mbi)
  : ArcaneFemObject(mbi)
  , m_dofs_on_cells(mbi.subDomain()->traceMng())
  , m_csr_matrix(mbi.subDomain()->traceMng())
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
  }

  void startInit() override; //! Method called at the beginning of the simulation
  void compute() override; //! Method called at each iteration
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

  void _assembleBilinearOperator();

 private:

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnCells m_dofs_on_cells;
  CsrFormat m_csr_matrix;

  Real f;

  String m_petsc_flags;
  String m_matrix_format = "DOK";

  IndexedCellCellConnectivityView m_cell_cell_connectivity_view;

  bool m_assemble_linear_system = true;
  bool m_solve_linear_system = true;
  bool m_cross_validation = false;

  void _doStationarySolve();
  void _getMaterialParameters();
  void _solve();
  void _assembleLinearSystem();
  void _buildCsrSparsity();
  void _updateVariables();
  void _validateResults();

};

#endif