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

#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/StringList.h>

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>

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

#if defined(USE_CUSPARSE_ADD)
//include for cusparse
#include <cusparse_v2.h>
#endif

//include for GPU use
#include "arcane/accelerator/core/IAcceleratorMng.h"
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
using namespace Arcane::FemUtils;
namespace ax = Arcane::Accelerator;
namespace fs = std::filesystem;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifdef USE_CUSPARSE_ADD
/**
 * @brief Macro for use of cusparse
 * 
 */
#define CHECK_CUSPARSE(func) \
  { \
    cusparseStatus_t status = (func); \
    if (status != CUSPARSE_STATUS_SUCCESS) { \
      printf("CUSPARSE API failed at line %d with error: %s (%d)\n", \
             __LINE__, cusparseGetErrorString(status), status); \
      printf("CUSPARSE STATUS INTERNAL ERROR : %d\n", status == CUSPARSE_STATUS_INTERNAL_ERROR); \
      printf("CUSPARSE STATUS EXECUTION FAILED : %d\n", status == CUSPARSE_STATUS_EXECUTION_FAILED); \
      printf("CUSPARSE STATUS INSUFFICENT RESOURCES: %d\n", status == CUSPARSE_STATUS_INSUFFICIENT_RESOURCES); \
      return; \
    } \
  }

#define CHECK_CUDA(func) \
  { \
    cudaError_t status = (func); \
    if (status != cudaSuccess) { \
      printf("CUDA API failed at line %d with error: %s (%d)\n", \
             __LINE__, cudaGetErrorString(status), status); \
      return; \
    } \
  }

/**
 * @brief struct for the csr of cusparse 
 * 
 */
struct cusparseCsr
{
  cusparseMatDescr_t desc;
  Int32 nnz = 0;
  Int32* csrRow;
  Int32* csrCol;
  float* csrVal;
};

#endif

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
  , m_coo_matrix(mbi.subDomain())
  , m_csr_matrix(mbi.subDomain())
  , m_time_stats(mbi.subDomain()->timeStats())
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
  bool m_use_csr = false;
  bool m_use_csr_gpu = false;
  bool m_use_nodewise_csr = false;
  bool m_use_buildless_csr = false;
  bool m_use_cusparse_add = false;
  bool m_use_legacy = true;
  bool m_running_on_gpu = false;
  ITimeStats* m_time_stats;

  CooFormat m_coo_matrix;

  CsrFormat m_csr_matrix;

  NumArray<Real, MDDim1> m_rhs_vect;

  std::ofstream logger;
  std::ofstream wbuild;
  std::ofstream timer;

  Integer cache_index;

 private:

  void fileNumArray(bool ref, NumArray<Real, MDDim1> numarray);

  void _handleFlags();
  void _doStationarySolve();
  void _getMaterialParameters();
  void _updateBoundayConditions();
  void _checkCellType();
  void _assembleBilinearOperatorTRIA3();
  void _assembleBilinearOperatorQUAD4();
  void _assembleBilinearOperatorTETRA4();
  void _solve();
  void _initBoundaryconditions();
  void _assembleLinearOperator();
  void _applyDirichletBoundaryConditions();
  void _checkResultFile();
  void _writeInJson();
  void _saveTimeInCSV();
  void _saveNoBuildTimeInCSV();
  void _benchBuildRow();
  Real _readTimeFromJson(String main_time, String sub_time);
  FixedMatrix<3, 3> _computeElementMatrixTRIA3(Cell cell);
  FixedMatrix<4, 4> _computeElementMatrixTETRA4(Cell cell);
  FixedMatrix<4, 4> _computeElementMatrixQUAD4(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
  Real _computeAreaTetra4(Cell cell);
  Real _computeAreaQuad4(Cell cell);
  Real _computeEdgeLength2(Face face);
  Real2 _computeEdgeNormal2(Face face);

 public:

  void _applyDirichletBoundaryConditionsGpu();
  void _assembleCsrGpuLinearOperator();
  static ARCCORE_HOST_DEVICE Int32
  _getValIndexCsrGpu(Int32 begin, Int32 end, DoFLocalId col, ax::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> csr_col);

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


#ifdef USE_CUSPARSE_ADD
  void printCsrMatrix(std::string fileName, cusparseCsr csr, bool is_coo);
  void _computeCusparseElementMatrix(cusparseCsr& result, cusparseCsr& global, Cell icell, cusparseHandle_t handle, IndexedNodeDoFConnectivityView node_dof);
  void _assembleCusparseBilinearOperatorTRIA3();
#endif
  void _buildMatrix();
  void _assembleCooBilinearOperatorTRIA3();

  void _buildMatrixSort();
  void _assembleCooSortBilinearOperatorTRIA3();

 public:

  static ARCCORE_HOST_DEVICE void
  _computeElementMatrixTRIA3GPU(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                                ax::VariableNodeReal3InView in_node_coord, Real K_e[9]);

 private:


 public:

  void _buildMatrixCsrGPU();
  void _assembleCsrGPUBilinearOperatorTRIA3();

 private:

#ifdef USE_COO_GPU

 public:

  void _buildMatrixGPU();
  void _assembleCooGPUBilinearOperatorTRIA3();

 private:

#endif
  void _assembleCsrBilinearOperatorTRIA3();
  void _assembleCsrBilinearOperatorTETRA4();
  void _buildMatrixCsr();
 public:

  void _buildMatrixNodeWiseCsr();
  void _assembleNodeWiseCsrBilinearOperatorTria3();

  void _buildMatrixBuildLessCsr();
  void _buildMatrixGpuBuildLessCsr();
  static ARCCORE_HOST_DEVICE Real
  _computeCellMatrixGpuTRIA3(CellLocalId icell, IndexedCellNodeConnectivityView cnc,
                             ax::VariableNodeReal3InView in_node_coord, Real b_matrix[6]);
  static ARCCORE_HOST_DEVICE void
  _addValueToGlobalMatrixTria3Gpu(Int32 begin, Int32 end, Int32 col,
                                  ax::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> in_out_col_csr,
                                  ax::NumArrayView<DataViewGetterSetter<Real>, MDDim1, DefaultLayout> in_out_val_csr, Real x);
  void _assembleBuildLessCsrBilinearOperatorTria3();

 private:

  void _assembleCsrLinearOperator();
  void _translateRhs();
  bool _isMasterRank() const;
};
