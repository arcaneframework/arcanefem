// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.h                                                 (C) 2022-2023 */
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

#include "arcane/AcceleratorRuntimeInitialisationInfo.h"

// Fichier à inclure afin de pouvoir effectuer des connectivités personnalisées
#include "arcane/core/IndexedItemConnectivityView.h"
#include "arcane/core/IIndexedIncrementalItemConnectivityMng.h"
#include "arcane/core/IIndexedIncrementalItemConnectivity.h"
#include "arcane/core/IIncrementalItemConnectivity.h"

//Fichier à inclure afin d'avoir le timer
#include "arcane/Timer.h"

#include "arcane/utils/ApplicationInfo.h"

#include "arcane/utils/CommandLineArguments.h"
#include "arcane/utils/ParameterList.h"
#include "arcane/core/IApplication.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
using namespace Arcane::FemUtils;
namespace ax = Arcane::Accelerator;

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

struct computeTimer
{
  double add_glob = 0;
  double compute_el = 0;
  double sort_coo = 0;
  double convert_coo = 0;
  double convert_coo_tot = 0;
  double convert_csr_tot = 0;
  double convert_tot = 0;
  double iter_time = 0;
  double compute_tot = 0;
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
  bool m_use_blcsr = false;
  bool m_use_cusparse_add = false;
  bool m_use_legacy = true;

  CooFormat m_coo_matrix;

  CsrFormat m_csr_matrix;

  std::ofstream logger;
  std::ofstream wbuild;
  std::ofstream timer;

  double lhs_time = 0;
  double rhs_time = 0;
  double solver_time = 0;

  Integer cache_index;

 private:

  void _handleFlags();
  void _doStationarySolve();
  void _getMaterialParameters();
  void _updateBoundayConditions();
  void _checkCellType();
  void _assembleBilinearOperatorTRIA3();
  void _assembleBilinearOperatorQUAD4();
  void _solve();
  void _initBoundaryconditions();
  void _assembleLinearOperator();
  void _applyDirichletBoundaryConditions();
  void _checkResultFile();
  FixedMatrix<3, 3> _computeElementMatrixTRIA3(Cell cell);
  FixedMatrix<4, 4> _computeElementMatrixQUAD4(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
  Real _computeAreaQuad4(Cell cell);
  Real _computeEdgeLength2(Face face);
  Real2 _computeEdgeNormal2(Face face);
#ifdef USE_CUSPARSE_ADD
  void printCsrMatrix(std::string fileName, cusparseCsr csr, bool is_coo);
  void _computeCusparseElementMatrix(cusparseCsr& result, cusparseCsr& global, Cell icell, cusparseHandle_t handle, IndexedNodeDoFConnectivityView node_dof,
                                     computeTimer& timer);
  void _assembleCusparseBilinearOperatorTRIA3();
#endif
  void _buildMatrix();
  void _assembleCooBilinearOperatorTRIA3();

  void _buildMatrixSort();
  void _assembleCooSortBilinearOperatorTRIA3();

 public:

  void _computeElementMatrixTRIA3GPU(CellLocalId icell, IndexedCellNodeConnectivityView cnc, ax::VariableNodeReal3InView in_node_coord, Real K_e[9]);

 private:

#ifdef ARCANE_HAS_CUDA

 public:

  void _buildMatrixCsrGPU();
  void _assembleCsrGPUBilinearOperatorTRIA3();

 private:

#endif
#ifdef USE_COO_GPU

 public:

  void _buildMatrixGPU();
  void _assembleCooGPUBilinearOperatorTRIA3();

 private:

#endif
  void _assembleCsrBilinearOperatorTRIA3();
  void _buildMatrixCsr();
#ifdef ARCANE_HAS_CUDA
 public:

  void _buildMatrixBuildLessCsr();
  void _assembleBLCsrBilinearOperatorTria3();

 private:

#endif
};