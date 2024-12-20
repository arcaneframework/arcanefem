// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.h                                                 (C) 2022-2024 */
/*                                                                           */
/* Simple module to test simple FEM mechanism.                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_POISSON_FEMMODULE_H
#define ARCANEFEM_POISSON_FEMMODULE_H
/*---------------------------------------------------------------------------*/

#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/StringList.h>

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/Connectivity.h>

#include "CooFormatMatrix.h"

#include "CsrFormatMatrix.h"

#include "IDoFLinearSystemFactory.h"
#include "Fem_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"

#include <fstream>
#include <iostream>
#include <chrono>
#include <filesystem>

#include "arcane/core/IIndexedIncrementalItemConnectivityMng.h"
#include "arcane/core/IIndexedIncrementalItemConnectivity.h"

//include for GPU use
#include "arcane/accelerator/core/IAcceleratorMng.h"
#include "arcane/accelerator/core/ProfileRegion.h"
#include "arcane/accelerator/Accelerator.h"
#include "arcane/accelerator/core/RunQueue.h"
#include "arcane/accelerator/Atomic.h"

#include <arcane/utils/NumArray.h>

//include for connectivity view
#include "arcane/UnstructuredMeshConnectivity.h"
#include "arcane/ItemGenericInfoListView.h"

// Pour avoir les vues sur les variables
#include "arcane/accelerator/VariableViews.h"

// Pour avoir les vues sur les NumArray
#include "arcane/accelerator/NumArrayViews.h"

// Fichier à inclure pour avoir RUNCOMMAND_ENUMERATE
#include "arcane/accelerator/RunCommandEnumerate.h"

// Fichier à inclure pour avoir RUNCOMMAND_LOOP
#include "arcane/accelerator/RunCommandLoop.h"

#include "arcane/accelerator/Reduce.h"
#include "arcane/accelerator/Accelerator.h"

// Fichier à inclure afin d'avoir le scan
#include "arcane/accelerator/Scan.h"

#include "arcane/AcceleratorRuntimeInitialisationInfo.h"

// Fichier à inclure afin de pouvoir effectuer des connectivités personnalisées
#include "arcane/core/IndexedItemConnectivityView.h"
#include "arcane/core/IIndexedIncrementalItemConnectivityMng.h"
#include "arcane/core/IIndexedIncrementalItemConnectivity.h"
#include "arcane/core/IIncrementalItemConnectivity.h"

//Fichier à inclure afin d'avoir le timer
#include "arcane/Timer.h"
#include "arcane/utils/JSONWriter.h"
#include "arcane/utils/JSONReader.h"
#include "arcane/core/ITimeStats.h"

#include "arcane/utils/ApplicationInfo.h"

#include "arcane/utils/CommandLineArguments.h"
#include "arcane/utils/ParameterList.h"
#include "arcane/core/IApplication.h"

#include "arcane/utils/ValueConvert.h"

#include <arcane/core/MeshUtils.h>
#include "BSRFormat.h"
#include "ArcaneFemFunctionsGpu.h"
#include "arcane/core/VariableTypedef.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
using namespace Arcane::FemUtils;
namespace ax = Arcane::Accelerator;
namespace fs = std::filesystem;

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
  , m_coo_matrix(mbi.subDomain()->traceMng())
  , m_csr_matrix(mbi.subDomain()->traceMng())
  , m_time_stats(mbi.subDomain()->timeStats())
  , m_bsr_format(mbi.subDomain()->traceMng(), *(mbi.subDomain()->acceleratorMng()->defaultQueue()), m_dofs_on_nodes)
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);

    addEntryPoint(this, "Build",
                  &FemModule::_build,
                  IEntryPoint::WBuild,
                  IEntryPoint::PAutoLoadBegin);
  }

 public:

  //! Method called at each iteration
  void compute() override;

  //! Method called at the beginning of the simulation
  void startInit() override;

  //! Method called at the end of the simulation
  void endModule() override;

  VersionInfo versionInfo() const override
  {
    return VersionInfo(1, 0, 0);
  }

 private:

  Real f;
  Real ElementNodes;

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnNodes m_dofs_on_nodes;
  bool m_register_time = false;
  bool m_arcane_timer = false;
  Integer m_cache_warming = 1;
  bool m_use_coo = false;
  bool m_use_coo_sort = false;
  bool m_use_coo_gpu = false;
  bool m_use_coo_sort_gpu = false;
  bool m_use_csr = false;
  bool m_use_csr_gpu = false;
  bool m_use_nodewise_csr = false;
  bool m_use_buildless_csr = false;
  bool m_use_legacy = true;
  bool m_use_bsr = false;
  bool m_running_on_gpu = false;
  bool m_solve_linear_system = true;
  bool m_cross_validation = true;

  ITimeStats* m_time_stats;

  CooFormat m_coo_matrix;

  CsrFormat m_csr_matrix;

  BSRFormat<1> m_bsr_format;

  NumArray<Real, MDDim1> m_rhs_vect;

  Integer cache_index;

  /*!
   * \brief List of nodes connected to an another node via a edge of the a cell.
   */
  Ref<IIndexedIncrementalItemConnectivity> m_node_node_via_edge_connectivity;

  //! Number of edges (only for 3D meshes)
  Int64 m_nb_edge = -1;

  //! Default queue used for computation.
  RunQueue m_queue;

 public:

  void _doStationarySolve();

 private:

  void _handleFlags();
  void _getMaterialParameters();
  void _updateBoundayConditions();
  void _checkCellType();
  void _assembleBilinearOperatorTRIA3();
  void _assembleBilinearOperatorTETRA4();
  void _solve();
  void _initBoundaryconditions();
  void _assembleLinearOperator(BSRMatrix<1>* bsr_matrix = nullptr);
  void _applyDirichletBoundaryConditions();
  void _checkResultFile();
  void _dumpTimeStats();
  FixedMatrix<3, 3> _computeElementMatrixTRIA3(Cell cell);
  FixedMatrix<4, 4> _computeElementMatrixTETRA4(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
  Real _computeAreaTetra4(Cell cell);
  Real _computeEdgeLength2(Face face);
  Real2 _computeEdgeNormal2(Face face);

 public:

  void _applyDirichletBoundaryConditionsGpu();
  void _assembleCsrGpuLinearOperator();

  static ARCCORE_HOST_DEVICE Int32
  _getValIndexCsrGpu(Int32 begin, Int32 end, DoFLocalId col, ax::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> csr_col);

  static ARCCORE_HOST_DEVICE Real
  _computeAreaTetra4Gpu(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                        ax::VariableNodeReal3InView in_node_coord);
  static ARCCORE_HOST_DEVICE Real
  _computeAreaTriangle3Gpu(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                           ax::VariableNodeReal3InView in_node_coord);
  static ARCCORE_HOST_DEVICE Real
  _computeEdgeLength2Gpu(FaceLocalId iface,
                         IndexedFaceNodeConnectivityView fnc,
                         ax::VariableNodeReal3InView in_node_coord);
  static ARCCORE_HOST_DEVICE Real2
  _computeEdgeNormal2Gpu(FaceLocalId iface, IndexedFaceNodeConnectivityView fnc,
                         ax::VariableNodeReal3InView in_node_coord,
                         FaceInfoListView faces_infos);

 private:

  void _buildMatrixCoo();
  void _assembleCooBilinearOperatorTRIA3();
  void _assembleCooBilinearOperatorTETRA4();

  void _buildMatrixCooSort();
  void _assembleCooSortBilinearOperatorTRIA3();
  void _assembleCooSortBilinearOperatorTETRA4();

 public:

  static ARCCORE_HOST_DEVICE void
  _computeElementMatrixTRIA3GPU(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                                ax::VariableNodeReal3InView in_node_coord, Real K_e[9]);

  static ARCCORE_HOST_DEVICE inline void
  _computeElementMatrixTETRA4GPU(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                                 ax::VariableNodeReal3InView in_node_coord, Real K_e[16]);

 private:
 public:

  void _buildMatrixCsrGPU();
  void _buildOffsets(const SmallSpan<uint>& offsets_smallspan);
  void _fillDiagonal(Int64 nb_edge, NodeGroup nodes);
  void _assembleCsrGPUBilinearOperatorTRIA3();
  void _assembleCsrGPUBilinearOperatorTETRA4();

 private:
 public:

  void _buildMatrixCooGPU();
  void _assembleCooGPUBilinearOperatorTRIA3();
  void _assembleCooGPUBilinearOperatorTETRA4();

  void _buildMatrixCooSortGPU();
  void _assembleCooSortGPUBilinearOperatorTRIA3();
  void _assembleCooSortGPUBilinearOperatorTETRA4();

 private:

  void _assembleCsrBilinearOperatorTRIA3();
  void _assembleCsrBilinearOperatorTETRA4();
  void _buildMatrixCsr();

 public:

  void _buildMatrixNodeWiseCsr();
  void _buildOffsetsNodeWiseCsr(const SmallSpan<uint>& offsets_smallspan);
  void _assembleNodeWiseCsrBilinearOperatorTria3();
  void _assembleNodeWiseCsrBilinearOperatorTetra4();

  void _buildMatrixBuildLessCsr();
  void _buildMatrixGpuBuildLessCsr();
  static ARCCORE_HOST_DEVICE Real
  _computeCellMatrixGpuTETRA4(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                              ax::VariableNodeReal3InView in_node_coord, Real b_matrix[12]);
  static ARCCORE_HOST_DEVICE inline Real
  _computeCellMatrixGpuTRIA3(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                             ax::VariableNodeReal3InView in_node_coord, Real b_matrix[6]);
  static ARCCORE_HOST_DEVICE void
  _addValueToGlobalMatrixGpu(Int32 begin, Int32 end, Int32 col,
                             ax::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> in_out_col_csr,
                             ax::NumArrayView<DataViewGetterSetter<Real>, MDDim1, DefaultLayout> in_out_val_csr, Real x);
  void _assembleBuildLessCsrBilinearOperatorTria3();
  void _assembleBuildLessCsrBilinearOperatorTetra4();

 private:

  void _build();
  void _assembleCsrLinearOperator();
  void _translateRhs();
  bool _isMasterRank() const;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline void FemModule::
_computeElementMatrixTRIA3GPU(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                              ax::VariableNodeReal3InView in_node_coord, Real K_e[9])
{
  // Get coordiantes of the triangle element  TRI3
  //------------------------------------------------
  //                  0 o
  //                   . .
  //                  .   .
  //                 .     .
  //              1 o . . . o 2
  //------------------------------------------------
  // We might want to replace the next 4 lines of codes with _computeAreaTriangle3Gpu()
  Real3 m0 = in_node_coord[cnc.nodeId(icell, 0)];
  Real3 m1 = in_node_coord[cnc.nodeId(icell, 1)];
  Real3 m2 = in_node_coord[cnc.nodeId(icell, 2)];

  Real area = 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y)); //_computeAreaTriangle3Gpu(icell, cnc, in_node_coord);

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  //We will want to replace fixed matrix by some numarray ? Will not work because NumArray function are host functions
  //NumArray<Real, ExtentsV<2, 3>> b_matrix(eMemoryRessource::Device);

  Real A2 = 2.0 * area;
  Real b_matrix[2][3] = { { dPhi0.x / A2, dPhi1.x / A2, dPhi2.x / A2 },
                          { dPhi0.y / A2, dPhi1.y / A2, dPhi2.y / A2 } };

  //NumArray<Real, ExtentsV<3, 3>> int_cdPi_dPj;

  //Multiplying b_matrix by its transpose, and doing the mult in place in the same loop
  // Compute the upper triangular part of the matrix
  for (Int32 i = 0; i < 3; ++i) {
    for (Int32 j = i; j < 3; ++j) {
      for (Int32 k = 0; k < 2; ++k) {
        K_e[i * 3 + j] += b_matrix[k][i] * b_matrix[k][j];
      }
      // Multiply by A2 to complete the matrix
      K_e[i * 3 + j] *= area;

      // Mirror to the lower triangular part
      K_e[j * 3 + i] = K_e[i * 3 + j];
    }
  }

  //info() << "Cell=" << cell.localId();
  //std::cout << " int_cdPi_dPj=";
  //int_cdPi_dPj.dump(std::cout);
  //std::cout << "\n";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline void FemModule::
_computeElementMatrixTETRA4GPU(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                               ax::VariableNodeReal3InView in_node_coord, Real K_e[16])
{
  // Get coordinates of the triangle element  TETRA4
  //------------------------------------------------
  //                3 o
  //                 /|\
  //                / | \
  //               /  |  \
  //              /   o 2 \
  //             / .    .  \
  //            o-----------o
  //            0           1
  //------------------------------------------------
  Real3 m0 = in_node_coord[cnc.nodeId(icell, 0)];
  Real3 m1 = in_node_coord[cnc.nodeId(icell, 1)];
  Real3 m2 = in_node_coord[cnc.nodeId(icell, 2)];
  Real3 m3 = in_node_coord[cnc.nodeId(icell, 3)];

  Real3 v0 = m1 - m0;
  Real3 v1 = m2 - m0;
  Real3 v2 = m3 - m0;

  Real volume = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2))) / 6.0;

  // Compute gradients of shape functions
  Real3 dPhi0 = Arcane::math::cross(m2 - m1, m1 - m3);
  Real3 dPhi1 = Arcane::math::cross(m3 - m0, m0 - m2);
  Real3 dPhi2 = Arcane::math::cross(m1 - m0, m0 - m3);
  Real3 dPhi3 = Arcane::math::cross(m0 - m1, m1 - m2);

  Real mul = 1.0 / (6.0 * volume);

  // Construct the B-matrix
  Real b_matrix[3][4];
  b_matrix[0][0] = dPhi0.x;
  b_matrix[1][0] = dPhi0.y;
  b_matrix[2][0] = dPhi0.z;

  b_matrix[0][1] = dPhi1.x;
  b_matrix[1][1] = dPhi1.y;
  b_matrix[2][1] = dPhi1.z;

  b_matrix[0][2] = dPhi2.x;
  b_matrix[1][2] = dPhi2.y;
  b_matrix[2][2] = dPhi2.z;

  b_matrix[0][3] = dPhi3.x;
  b_matrix[1][3] = dPhi3.y;
  b_matrix[2][3] = dPhi3.z;

  for (Int32 i = 0; i < 3; ++i)
    for (Int32 j = 0; j < 4; ++j)
      b_matrix[i][j] *= mul;

  // Compute the element matrix
  for (Int32 i = 0; i < 4; ++i) {
    for (Int32 j = i; j < 4; ++j) {
      for (Int32 k = 0; k < 3; ++k)
        K_e[i * 4 + j] += b_matrix[k][i] * b_matrix[k][j];
      K_e[i * 4 + j] *= volume;
      K_e[j * 4 + i] = K_e[i * 4 + j];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE
inline Real FemModule::_computeCellMatrixGpuTRIA3(CellLocalId icell, IndexedCellNodeConnectivityView cnc, ax::VariableNodeReal3InView in_node_coord, Real b_matrix[6])
{
  Real3 m0 = in_node_coord[cnc.nodeId(icell, 0)];
  Real3 m1 = in_node_coord[cnc.nodeId(icell, 1)];
  Real3 m2 = in_node_coord[cnc.nodeId(icell, 2)];

  Real area = 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y)); //_computeAreaTriangle3Gpu(icell, cnc, in_node_coord);

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  Real mul = 0.5 / area;
  b_matrix[0] = dPhi0.x * mul;
  b_matrix[1] = dPhi0.y * mul;

  b_matrix[2] = dPhi1.x * mul;
  b_matrix[3] = dPhi1.y * mul;

  b_matrix[4] = dPhi2.x * mul;
  b_matrix[5] = dPhi2.y * mul;

  return area;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE
inline Real FemModule::_computeCellMatrixGpuTETRA4(CellLocalId icell, IndexedCellNodeConnectivityView cnc, ax::VariableNodeReal3InView in_node_coord, Real b_matrix[12])
{
  Real3 m0 = in_node_coord[cnc.nodeId(icell, 0)];
  Real3 m1 = in_node_coord[cnc.nodeId(icell, 1)];
  Real3 m2 = in_node_coord[cnc.nodeId(icell, 2)];
  Real3 m3 = in_node_coord[cnc.nodeId(icell, 3)];

  // Calculate vectors representing edges of the tetrahedron
  Real3 v0 = m1 - m0;
  Real3 v1 = m2 - m0;
  Real3 v2 = m3 - m0;

  // Compute volume using scalar triple product
  Real volume = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2))) / 6.0;

  // Compute gradients of shape functions
  Real3 dPhi0 = Arcane::math::cross(m2 - m1, m1 - m3);
  Real3 dPhi1 = Arcane::math::cross(m3 - m0, m0 - m2);
  Real3 dPhi2 = Arcane::math::cross(m1 - m0, m0 - m3);
  Real3 dPhi3 = Arcane::math::cross(m0 - m1, m1 - m2);

  // Construct the B-matrix as a vector
  Real mul = 1.0 / (6.0 * volume);
  b_matrix[0] = dPhi0.x * mul;
  b_matrix[1] = dPhi0.y * mul;
  b_matrix[2] = dPhi0.z * mul;

  b_matrix[3] = dPhi1.x * mul;
  b_matrix[4] = dPhi1.y * mul;
  b_matrix[5] = dPhi1.z * mul;

  b_matrix[6] = dPhi2.x * mul;
  b_matrix[7] = dPhi2.y * mul;
  b_matrix[8] = dPhi2.z * mul;

  b_matrix[9] = dPhi3.x * mul;
  b_matrix[10] = dPhi3.y * mul;
  b_matrix[11] = dPhi3.z * mul;

  return volume;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
