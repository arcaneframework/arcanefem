// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.h                                                (C) 2000-2026  */
/*                                                                           */
/* FemModuleElastodynamics class definition.                                 */
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
#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/accelerator/core/RunQueue.h>

#include "IArcaneFemBC.h"
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
namespace ax = Arcane::Accelerator;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Module Fem.
 */
class FemModuleElastodynamics
: public ArcaneFemObject
{
 public:

  explicit FemModuleElastodynamics(const ModuleBuildInfo& mbi)
  : ArcaneFemObject(mbi)
  , m_dofs_on_nodes(mbi.subDomain()->traceMng())
  , m_bsr_format(mbi.subDomain()->traceMng(), *(mbi.subDomain()->acceleratorMng()->defaultQueue()), m_dofs_on_nodes)
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
  }
  ~FemModuleElastodynamics()
  {
    for( const CaseTableInfo&  t : m_traction_case_table_list )
      delete t.case_table;
  }

 void startInit() override; //! Method called at the beginning of the simulation
 void compute() override; //! Method called at each iteration
 VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 void _assembleBilinearOperatorTria3Gpu();
 void _assembleBilinearOperatorTetra4Gpu();
 void _assembleDirichletsGpu();

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

  String m_petsc_flags;
  String m_matrix_format = "DOK";

  bool m_assemble_linear_system = true;
  bool m_solve_linear_system = true;
  bool m_cross_validation = false;
  bool m_hex_quad_mesh = false;

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;
  BSRFormat m_bsr_format;

  // List of CaseTable for traction boundary conditions
  UniqueArray<CaseTableInfo> m_traction_case_table_list;

 private:

  void _doStationarySolve();
  void _getParameters();
  void _updateVariables();
  void _updateTime();
  void _solve();
  void _assembleLinearOperator();
  void _validateResults();
  void _readCaseTables();
  void _assembleBilinearOperator();

  inline void _applyTraction(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  inline void _applyDirichlet(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  inline void _applySourceTerm(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  inline void _applySourceTermTria3(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  inline void _applySourceTermQuad4(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  inline void _applySourceTermTetra4(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);
  inline void _applySourceTermHexa8(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof);

  RealMatrix<6, 6> _computeElementMatrixTria3(Cell cell);
  RealMatrix<12, 12> _computeElementMatrixTetra4(Cell cell);
  RealMatrix<8, 8> _computeElementMatrixQuad4(Cell cell);
  RealMatrix<24, 24> _computeElementMatrixHexa8(Cell cell);

  template <int N>
  void _assembleBilinearOperatorCpu(const std::function<RealMatrix<N, N>(const Cell&)>& compute_element_matrix);
};

#endif