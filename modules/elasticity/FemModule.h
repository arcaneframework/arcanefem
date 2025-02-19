﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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

#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/ParameterList.h>
#include <arcane/utils/ApplicationInfo.h>
#include <arcane/utils/StringList.h>

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/core/ItemTypes.h>
#include <arccore/base/ArccoreGlobal.h>
#include "arccore/base/NotImplementedException.h"

#include "IDoFLinearSystemFactory.h"
#include "Fem_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "BSRFormat.h"

#include "ArcaneFemFunctions.h"
#include "ArcaneFemFunctionsGpu.h"

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
  , m_bsr_format(mbi.subDomain()->traceMng(), *(mbi.subDomain()->acceleratorMng()->defaultQueue()), m_dofs_on_nodes)
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
  }

  void startInit() override; //! Method called at the beginning of the simulation
  void compute() override; //! Method called at each iteration
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

  void _doStationarySolve();
  void _assembleBilinearOperator();
  void _assembleDirichletsGpu();

 private:

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;
  BSRFormat<2> m_bsr_format;

  Real E;
  Real nu;
  Real mu2;
  Real lambda;

  Real3 f;
  Real3 t;

  String m_petsc_flags;
  String m_matrix_format = "DOK";

  bool m_assemble_linear_system = true;
  bool m_solve_linear_system = true;
  bool m_cross_validation = true;
  bool m_use_bsr_matrix_format = false;
  bool m_use_atomic_free_bsr_matrix_format = false;

  bool m_use_bsr = false;
  bool m_use_bsr_atomic_free = false;

  void _getMaterialParameters();
  void _assembleBilinearOperatorTRIA3();
  void _solve();
  void _assembleLinearOperator();
  void _validateResults();
  void _updateVariables();
  void _initBsr();

  void _handleCommandLineFlags();
  void _setPetscFlagsFromCommandline();
  void _printArcaneFemTime(const String label, const Real value);

  Real _computeEdgeLength2(Face face);
  Real _computeAreaTriangle3(Cell cell);

  FixedMatrix<6, 6> _computeElementMatrixTRIA3(Cell cell);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<6, 6> computeElementMatrixTRIA3Base(Real3 m0, Real3 m1, Real3 m2, Real area, Real lambda, Real mu2)
{
  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  Real6 dxU1 = { dPhi0.x, 0., dPhi1.x, 0., dPhi2.x, 0. };
  Real6 dyU1 = { dPhi0.y, 0., dPhi1.y, 0., dPhi2.y, 0. };
  Real6 dxU2 = { 0., dPhi0.x, 0., dPhi1.x, 0., dPhi2.x };
  Real6 dyU2 = { 0., dPhi0.y, 0., dPhi1.y, 0., dPhi2.y };

  // -----------------------------------------------------------------------------
  //  step1 = (dx(u1)dx(v1) + dy(u2)dx(v1) + dx(u1)dy(v2) + dy(u2)dy(v2)) * lambda
  //------------------------------------------------------------------------------
  FixedMatrix<6, 6> step1_res = (((dxU1 ^ dxU1) + (dxU1 ^ dyU2) + (dyU2 ^ dxU1) + (dyU2 ^ dyU2)) * lambda) / (4 * area);

  // -----------------------------------------------------------------------------
  //  step2 = 2*mu * (dx(u1)dx(v1) + dy(u2)dy(v2) + 0.5 *
  //                  (dy(u1)dy(v1) + dx(u2)dy(v1) + dy(u1)dx(v2) + dx(u2)dx(v2)))
  //------------------------------------------------------------------------------
  FixedMatrix<6, 6> step2_res = ((((dxU1 ^ dxU1) + (dyU2 ^ dyU2)) + ((dyU1 ^ dyU1) + (dyU1 ^ dxU2) + (dxU2 ^ dyU1) + (dxU2 ^ dxU2)) * 0.5) * mu2) / (4 * area);

  return (step1_res + step2_res);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<6, 6> computeElementMatrixTRIA3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu2)
{
  Real3 m0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 m1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 m2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);
  return computeElementMatrixTRIA3Base(m0, m1, m2, area, lambda, mu2);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<2, 6> computeElementVectorTRIA3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu2, Int32 node_lid)
{
  Real3 m0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 m1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 m2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  FixedMatrix<1, 6> dxU1 = { dPhi0.x, 0., dPhi1.x, 0., dPhi2.x, 0. };
  FixedMatrix<1, 6> dyU1 = { dPhi0.y, 0., dPhi1.y, 0., dPhi2.y, 0. };
  FixedMatrix<1, 6> dxU2 = { 0., dPhi0.x, 0., dPhi1.x, 0., dPhi2.x };
  FixedMatrix<1, 6> dyU2 = { 0., dPhi0.y, 0., dPhi1.y, 0., dPhi2.y };

  // -----------------------------------------------------------------------------
  //  step1 = (dx(u1)dx(v1) + dy(u2)dx(v1) + dx(u1)dy(v2) + dy(u2)dy(v2)) * lambda
  //------------------------------------------------------------------------------
  FixedMatrix<1, 6> step1_res1 = (((dxU1(0,node_lid*2) * dxU1) + (dxU1(0,node_lid*2) * dyU2) + (dyU2(0,node_lid*2) * dxU1) + (dyU2(0,node_lid*2) * dyU2)) * lambda) / (4 * area);
  FixedMatrix<1, 6> step1_res2 = (((dxU1(0,node_lid*2+1) * dxU1) + (dxU1(0,node_lid*2+1) * dyU2) + (dyU2(0,node_lid*2+1) * dxU1) + (dyU2(0,node_lid*2+1) * dyU2)) * lambda) / (4 * area);

  // -----------------------------------------------------------------------------
  //  step2 = 2*mu * (dx(u1)dx(v1) + dy(u2)dy(v2) + 0.5 *
  //                  (dy(u1)dy(v1) + dx(u2)dy(v1) + dy(u1)dx(v2) + dx(u2)dx(v2)))
  //------------------------------------------------------------------------------
  FixedMatrix<1, 6> step2_res1 = ((((dxU1(0,node_lid*2) * dxU1) + (dyU2(0,node_lid*2) * dyU2)) + ((dyU1(0,node_lid*2) * dyU1) + (dyU1(0,node_lid*2) * dxU2) + (dxU2(0,node_lid*2) * dyU1) + (dxU2(0,node_lid*2) * dxU2)) * 0.5) * mu2) / (4 * area);
  FixedMatrix<1, 6> step2_res2 = ((((dxU1(0,node_lid*2+1) * dxU1) + (dyU2(0,node_lid*2+1) * dyU2)) + ((dyU1(0,node_lid*2+1) * dyU1) + (dyU1(0,node_lid*2+1) * dxU2) + (dxU2(0,node_lid*2+1) * dyU1) + (dxU2(0,node_lid*2+1) * dxU2)) * 0.5) * mu2) / (4 * area);

  FixedMatrix <1, 6> result1 = step1_res1 + step2_res1;
  FixedMatrix <1, 6> result2 = step1_res2 + step2_res2;

  FixedMatrix <2, 6> result;
  result(0,0)=result1(0,0); result(0,1)=result1(0,1); result(0,2)=result1(0,2); result(0,3)=result1(0,3); result(0,4)=result1(0,4); result(0,5)=result1(0,5);
  result(1,0)=result2(0,0); result(1,1)=result2(0,1); result(1,2)=result2(0,2); result(1,3)=result2(0,3); result(1,4)=result2(0,4); result(1,5)=result2(0,5);

  return result;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<6, 6> FemModule::_computeElementMatrixTRIA3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
  return computeElementMatrixTRIA3Base(m0, m1, m2, area, lambda, mu2);
}

#endif