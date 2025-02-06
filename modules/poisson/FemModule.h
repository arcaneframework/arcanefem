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

#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/PlatformUtils.h>
#include <arcane/utils/ParameterList.h>
#include <arcane/utils/ApplicationInfo.h>
#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/StringList.h>
#include <arcane/utils/NumArray.h>

#include <arcane/ITimeLoopMng.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/IMesh.h>

#include <arccore/base/NotImplementedException.h>

#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/accelerator/VariableViews.h>

#include "IDoFLinearSystemFactory.h"
#include "ArcaneFemFunctionsGpu.h"
#include "ArcaneFemFunctions.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "IArcaneFemBC.h"
#include "BSRFormat.h"
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

  void startInit() override; //! Method called at the beginning of the simulation
  void compute() override; //! Method called at each iteration
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

  void _assembleBilinearOperator();
  void _assembleLinearOperatorGpu();

 private:

  BSRFormat<1> m_bsr_format;
  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnNodes m_dofs_on_nodes;

  Real f;

  String m_petsc_flags;
  String m_matrix_format = "DOK";

  bool m_assemble_linear_system = true;
  bool m_solve_linear_system = true;
  bool m_cross_validation = true;
  bool m_use_bsr_matrix_format = false;
  bool m_use_atomic_free_bsr_matrix_format = false;


  void _doStationarySolve();
  void _getMaterialParameters();
  void _solve();
  void _assembleLinearOperatorCpu();
  void _assembleLinearOperator();

  void _updateVariables();
  void _validateResults();
  void _handleCommandLineFlags();
  void _setPetscFlagsFromCommandline();

  void _printArcaneFemTime(const String label, const Real value);

  FixedMatrix<3, 3> _computeElementMatrixTria3(Cell cell);
  FixedMatrix<4, 4> _computeElementMatrixTetra4(Cell cell);

  template <int N>
  void _assembleBilinear(const std::function<FixedMatrix<N, N>(const Cell&)>& compute_element_matrix);
};

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (P1 FE).
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
 *
 * This function calculates the integral of the expression:
 * integral2D (u.dx * v.dx + u.dy * v.dy)
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the gradients of the shape functions.
 * 3. Return (u.dx * v.dx + u.dy * v.dy);
 */
/*---------------------------------------------------------------------------*/

FixedMatrix<3, 3> FemModule::_computeElementMatrixTria3(Cell cell)
{
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  Real3 dxU = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyU = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

  return area * (dxU ^ dxU) + area * (dyU ^ dyU);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<3, 3> _computeElementMatrixTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord)
{
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  Real3 dxU = FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyU = FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);

  return area * (dxU ^ dxU) + area * (dyU ^ dyU);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a tetrahedral element (P1 FE).
 *
 * This function calculates the integral of the expression:
 * integral3D (u.dx * v.dx + u.dy * v.dy)
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the gradients of the shape functions.
 * 3. Return (u.dx * v.dx + u.dy * v.dy);
 */
/*---------------------------------------------------------------------------*/

FixedMatrix<4, 4> FemModule::_computeElementMatrixTetra4(Cell cell)
{
  Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);

  Real4 dxU = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
  Real4 dyU = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
  Real4 dzU = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);

  return volume * (dxU ^ dxU) + volume * (dyU ^ dyU) + volume * (dzU ^ dzU);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<4, 4> _computeElementMatrixTetra4Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord)
{
  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  Real4 dxU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);

  return volume * (dxU ^ dxU) + volume * (dyU ^ dyU) + volume * (dzU ^ dzU);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
