// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2023 */
/*                                                                           */
/* Simple module to test simple FEM mechanism.                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

//Having H files for each steps of the FEM

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

#include <iostream>
#include <fstream>
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

  ofstream logger;
  ofstream wbuild;
  ofstream timer;

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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "Module Fem COMPUTE";

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());

  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");
  // Test for adding parameters for PETSc.
  // This is only used for the first call.
  {
    StringList string_list;
    string_list.add("-trmalloc");
    string_list.add("-log_trace");

    string_list.add("-ksp_monitor");
    string_list.add("-ksp_view");
    string_list.add("-math_view");
    string_list.add("draw");
    string_list.add("-draw_pause");
    string_list.add("-10");

    CommandLineArguments args(string_list);
    m_linear_system.setSolverCommandLineArguments(args);
  }
  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();

  _doStationarySolve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "Module Fem INIT";

  if (m_register_time) {
    logger = ofstream("timer.txt");
    wbuild = ofstream("with_build.csv", std::ios_base::app);
    wbuild << nbNode() << ",";
    timer = ofstream("timer.csv", std::ios_base::app);
    timer << nbNode() << ",";
  }

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();

  //_buildDoFOnNodes();
  //Int32 nb_node = allNodes().size();
  //m_k_matrix.resize(nb_node, nb_node);
  //m_k_matrix.fill(0.0);

  //m_rhs_vector.resize(nb_node);
  //m_rhs_vector.fill(0.0);

  // # init mesh
  // # init behavior
  // # init behavior on mesh entities
  // # init BCs
  _handleFlags();
  _initBoundaryconditions();

  _checkCellType();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule ::
_handleFlags()
{
  ParameterList parameter_list = this->subDomain()->application()->applicationInfo().commandLineArguments().parameters();
  info() << "-----------------------------------------------------------------------------------------";
  info() << "The time will be registered by arcane in the output/listing/logs.0 file";
  if (parameter_list.getParameterOrNull("REGISTER_TIME") == "TRUE") {
    m_register_time = true;
    info() << "The time will also be registered in the timer.txt, with_build.csv and timer.csv file";
  }
  if (parameter_list.getParameterOrNull("CACHE_WARMING") != NULL) {
    std::stringstream test("");
    test << parameter_list.getParameterOrNull("CACHE_WARMING");
    m_cache_warming = std::stoi(test.str());
    info() << "A cache warming of " << m_cache_warming << " iterations will happen";
  }
  if (parameter_list.getParameterOrNull("COO") == "TRUE") {
    m_use_coo = true;
    info() << "The COO datastructure and its associated methods will be used";
  }
  if (parameter_list.getParameterOrNull("COO_SORT") == "TRUE") {
    m_use_coo_sort = true;
    info() << "The COO with sorting datastructure and its associated methods will be used";
  }
  if (parameter_list.getParameterOrNull("CSR") == "TRUE") {
    m_use_csr = true;
    info() << "The CSR datastructure and its associated methods will be used";
  }
#ifdef ARCANE_HAS_CUDA
  if (parameter_list.getParameterOrNull("CSR_GPU") == "TRUE") {
    m_use_csr_gpu = true;
    info() << "The CSR datastructure GPU compatible and its associated methods will be used";
  }
  if (parameter_list.getParameterOrNull("BLCSR") == "TRUE") {
    m_use_blcsr = true;
    info() << "The BLCSR datastructure (GPU compatible) and its associated methods will be used";
  }
  if (parameter_list.getParameterOrNull("CUSPARSE_ADD") == "TRUE") {
    m_use_cusparse_add = true;
    info() << "CUSPARSE and its associated methods will be used";
  }
#endif
  if (parameter_list.getParameterOrNull("LEGACY") == "FALSE") {
    m_use_cusparse_add = false;
  }
  else {
    info() << "The Legacy datastructure and its associated methods will be used";
  }
  info() << "-----------------------------------------------------------------------------------------";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  Timer::Action timer_action(this->subDomain(), "StationarySolve");
  std::chrono::_V2::system_clock::time_point fem_start;
  if (m_register_time && m_cache_warming)
    fem_start = std::chrono::high_resolution_clock::now();

  // # get material parameters
  _getMaterialParameters();

  // # update BCs
  _updateBoundayConditions();

  // Assemble the FEM bilinear operator (LHS - matrix A)
  if (options()->meshType == "QUAD4")
    _assembleBilinearOperatorQUAD4();
  else {

#ifdef ARCANE_HAS_CUDA
    if (m_use_cusparse_add) {
      for (cache_index = 0; cache_index < m_cache_warming; cache_index++)
        _assembleCusparseBilinearOperatorTRIA3();
    }
  }
#endif
  if (m_use_coo) {
    for (cache_index = 0; cache_index < m_cache_warming; cache_index++)
      _assembleCooBilinearOperatorTRIA3();
    //m_coo_matrix.translateToLinearSystem(m_linear_system);
  }
  if (m_use_coo_sort) {
    for (cache_index = 0; cache_index < m_cache_warming; cache_index++)
      _assembleCooSortBilinearOperatorTRIA3();
    //m_coo_matrix.translateToLinearSystem(m_linear_system);
  }
#ifdef USE_COO_GPU
#ifdef m_cache_warming
  for (i = 0; i < 3; i++)
#endif
    _assembleCooGPUBilinearOperatorTRIA3();
    //m_coo_matrix.translateToLinearSystem(m_linear_system);
#endif
  if (m_use_csr) {
    for (cache_index = 0; cache_index < m_cache_warming; cache_index++)
      _assembleCsrBilinearOperatorTRIA3();
    //m_csr_matrix.translateToLinearSystem(m_linear_system);
  }
  if (m_use_legacy) {
    for (cache_index = 0; cache_index < m_cache_warming; cache_index++)
      _assembleBilinearOperatorTRIA3();
  }

#ifdef ARCANE_HAS_CUDA
  if (m_use_csr_gpu) {
    for (cache_index = 0; cache_index < m_cache_warming; cache_index++)
      _assembleCsrGPUBilinearOperatorTRIA3();
    //m_csr_matrix.translateToLinearSystem(m_linear_system);
  }
#endif
  if (m_use_blcsr) {
    for (cache_index = 0; cache_index < m_cache_warming; cache_index++)
      _assembleBLCsrBilinearOperatorTria3();
    //m_csr_matrix.translateToLinearSystem(m_linear_system);
  }

  // Assemble the FEM linear operator (RHS - vector b)
  //_assembleLinearOperator();

  // # T=linalg.solve(K,RHS)
  //_solve();

  // Check results
  //_checkResultFile();

  if (m_register_time) {
    auto fem_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> fem_duration = fem_stop - fem_start;
    double total_duration = fem_duration.count();
    logger << "FEM total duration : " << fem_duration.count() << "\n"
           << "LHS time in total duration : " << lhs_time / total_duration * 100 << "%\n"
           << "RHS time in total duration : " << rhs_time / total_duration * 100 << "%\n"
           << "Solver time in total duration : " << solver_time / total_duration * 100 << "%\n";
    logger.close();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "Get material parameters...";
  f = options()->f();
  ElementNodes = 3.;

  if (options()->meshType == "QUAD4")
    ElementNodes = 4.;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_initBoundaryconditions()
{
  info() << "Init boundary conditions...";

  info() << "Apply boundary conditions";
  _applyDirichletBoundaryConditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_applyDirichletBoundaryConditions()
{
  // Handle all the Dirichlet boundary conditions.
  // In the 'arc' file, there are in the following format:
  //   <dirichlet-boundary-condition>
  //   <surface>Haut</surface>
  //   <value>21.0</value>
  // </dirichlet-boundary-condition>

  for (const auto& bs : options()->dirichletBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real value = bs->value();
    info() << "Apply Dirichlet boundary condition surface=" << group.name() << " v=" << value;
    ENUMERATE_ (Face, iface, group) {
      for (Node node : iface->nodes()) {
        m_u[node] = value;
        m_u_dirichlet[node] = true;
      }
    }
  }

  for (const auto& bs : options()->dirichletPointCondition()) {
    NodeGroup group = bs->node();
    Real value = bs->value();
    info() << "Apply Dirichlet point condition node=" << group.name() << " v=" << value;
    ENUMERATE_ (Node, inode, group) {
      Node node = *inode;
      m_u[node] = value;
      m_u_dirichlet[node] = true;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_checkCellType()
{
  Int16 type = 0;
  if (options()->meshType == "QUAD4") {
    type = IT_Quad4;
  }
  else {
    type = IT_Triangle3;
  }
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    if (cell.type() != type)
      ARCANE_FATAL("Only Triangle3 cell type is supported");
  }
}

void FemModule::
_updateBoundayConditions()
{
  info() << "TODO " << A_FUNCINFO;
}

/*---------------------------------------------------------------------------*/
// Assemble the FEM linear operator
//  - This function enforces a Dirichlet boundary condition in a weak sense
//    via the penalty method
//  - The method also adds source term
//  - TODO: external fluxes
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "Assembly of FEM linear operator ";
  info() << "Applying Dirichlet boundary condition via  penalty method ";

  // time registration
  std::chrono::_V2::system_clock::time_point rhs_start;
  double penalty_time = 0;
  double wpenalty_time = 0;
  double sassembly_time = 0;
  double fassembly_time = 0;
  if (m_register_time) {
    rhs_start = std::chrono::high_resolution_clock::now();
  }

  Timer::Action timer_action(this->subDomain(), "AssembleLinearOperator");

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (options()->enforceDirichletMethod() == "Penalty") {

    Timer::Action timer_action(this->subDomain(), "Penalty");

    //----------------------------------------------
    // penalty method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'P' be the penalty term and let 'i' be the set of DOF for which
    //  Dirichlet condition needs to be applied
    //
    //  - For LHS matrix A the diag term corresponding to the Dirichlet DOF
    //           a_{i,i} = 1. * P
    //
    //  - For RHS vector b the term that corresponds to the Dirichlet DOF
    //           b_{i} = b_{i} * P
    //----------------------------------------------

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    Real Penalty = options()->penalty(); // 1.0e30 is the default
    std::chrono::_V2::system_clock::time_point penalty_start;
    if (m_register_time) {
      penalty_start = std::chrono::high_resolution_clock::now();
    }

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u_dirichlet[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_linear_system.matrixSetValue(dof_id, dof_id, Penalty);
        Real u_g = Penalty * m_u[node_id];
        rhs_values[dof_id] = u_g;
      }
    }
    std::chrono::_V2::system_clock::time_point penalty_stop;
    if (m_register_time) {
      penalty_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> penalty_duration = penalty_stop - penalty_start;
      penalty_time = penalty_duration.count();
      logger << "Penalty duration : " << penalty_time << "\n";
    }
  }
  else if (options()->enforceDirichletMethod() == "WeakPenalty") {
    Timer::Action timer_action(this->subDomain(), "WeakPenalty");

    //----------------------------------------------
    // weak penalty method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'P' be the penalty term and let 'i' be the set of DOF for which
    //  Dirichlet condition needs to be applied
    //
    //  - For LHS matrix A the diag term corresponding to the Dirichlet DOF
    //           a_{i,i} = a_{i,i} + P
    //
    //  - For RHS vector b the term that corresponds to the Dirichlet DOF
    //           b_{i} = b_{i} * P
    //----------------------------------------------

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    Real Penalty = options()->penalty(); // 1.0e30 is the default
    std::chrono::_V2::system_clock::time_point wpenalty_start;
    if (m_register_time) {
      wpenalty_start = std::chrono::high_resolution_clock::now();
    }

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u_dirichlet[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_linear_system.matrixAddValue(dof_id, dof_id, Penalty);
        Real u_g = Penalty * m_u[node_id];
        rhs_values[dof_id] = u_g;
      }
    }
    std::chrono::_V2::system_clock::time_point wpenalty_stop;
    if (m_register_time) {
      wpenalty_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> wpenalty_duration = wpenalty_stop - wpenalty_start;
      wpenalty_time = wpenalty_duration.count();
      logger << "Weak Penalty duration : " << wpenalty_time << "\n";
    }
  }
  else if (options()->enforceDirichletMethod() == "RowElimination") {

    //----------------------------------------------
    // Row elimination method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'I' be the set of DOF for which  Dirichlet condition needs to be applied
    //
    //  to apply the Dirichlet on 'i'th DOF
    //  - For LHS matrix A the row terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j
    //           a_{i,j} = 1.  : i==j
    //----------------------------------------------

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    // TODO
  }
  else if (options()->enforceDirichletMethod() == "RowColumnElimination") {

    //----------------------------------------------
    // Row elimination method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'I' be the set of DOF for which  Dirichlet condition needs to be applied
    //
    //  to apply the Dirichlet on 'i'th DOF
    //  - For LHS matrix A the row terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j  for all j
    //           a_{i,j} = 1.  : i==j
    //    also the column terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j  for all i
    //----------------------------------------------

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    // TODO
  }
  else {

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " is not supported \n"
           << "enforce-Dirichlet-method only supports:\n"
           << "  - Penalty\n"
           << "  - WeakPenalty\n"
           << "  - RowElimination\n"
           << "  - RowColumnElimination\n";
  }
  std::chrono::_V2::system_clock::time_point sassemby_start;
  if (m_register_time)
    sassemby_start = std::chrono::high_resolution_clock::now();

  std::chrono::_V2::system_clock::time_point fassemby_start;
  {
    Timer::Action timer_action(this->subDomain(), "ConstantSourceTermAssembly");
    //----------------------------------------------
    // Constant source term assembly
    //----------------------------------------------
    //
    //  $int_{Omega}(f*v^h)$
    //  only for noded that are non-Dirichlet
    //----------------------------------------------
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real area = _computeAreaTriangle3(cell);
      for (Node node : cell.nodes()) {
        if (!(m_u_dirichlet[node]) && node.isOwn())
          rhs_values[node_dof.dofId(node, 0)] += f * area / ElementNodes;
      }
    }
    std::chrono::_V2::system_clock::time_point sassemby_stop;
    if (m_register_time) {
      sassemby_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> sassembly_duration = sassemby_stop - sassemby_start;
      sassembly_time = sassembly_duration.count();
      logger << "Constant source term assembly duration : " << sassembly_time << "\n";
      fassemby_start = std::chrono::high_resolution_clock::now();
    }
  }
  {
    Timer::Action timer_action(this->subDomain(), "ConstantSourceTermAssembly");

    //----------------------------------------------
    // Constant flux term assembly
    //----------------------------------------------
    //
    //  only for noded that are non-Dirichlet
    //  $int_{dOmega_N}((q.n)*v^h)$
    // or
    //  $int_{dOmega_N}((n_x*q_x + n_y*q_y)*v^h)$
    //----------------------------------------------
    for (const auto& bs : options()->neumannBoundaryCondition()) {
      FaceGroup group = bs->surface();

      if (bs->value.isPresent()) {
        Real value = bs->value();
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real length = _computeEdgeLength2(face);
          for (Node node : iface->nodes()) {
            if (!(m_u_dirichlet[node]) && node.isOwn())
              rhs_values[node_dof.dofId(node, 0)] += value * length / 2.;
          }
        }
        continue;
      }

      if (bs->valueX.isPresent() && bs->valueY.isPresent()) {
        Real valueX = bs->valueX();
        Real valueY = bs->valueY();
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real length = _computeEdgeLength2(face);
          Real2 Normal = _computeEdgeNormal2(face);
          for (Node node : iface->nodes()) {
            if (!(m_u_dirichlet[node]) && node.isOwn())
              rhs_values[node_dof.dofId(node, 0)] += (Normal.x * valueX + Normal.y * valueY) * length / 2.;
          }
        }
        continue;
      }

      if (bs->valueX.isPresent()) {
        Real valueX = bs->valueX();
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real length = _computeEdgeLength2(face);
          Real2 Normal = _computeEdgeNormal2(face);
          for (Node node : iface->nodes()) {
            if (!(m_u_dirichlet[node]) && node.isOwn())
              rhs_values[node_dof.dofId(node, 0)] += (Normal.x * valueX) * length / 2.;
          }
        }
        continue;
      }

      if (bs->valueY.isPresent()) {
        Real valueY = bs->valueY();
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real length = _computeEdgeLength2(face);
          Real2 Normal = _computeEdgeNormal2(face);
          for (Node node : iface->nodes()) {
            if (!(m_u_dirichlet[node]) && node.isOwn())
              rhs_values[node_dof.dofId(node, 0)] += (Normal.y * valueY) * length / 2.;
          }
        }
        continue;
      }
    }
    std::chrono::_V2::system_clock::time_point fassemby_stop;
    if (m_register_time) {

      fassemby_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> fassembly_duration = fassemby_stop - fassemby_start;
      fassembly_time = fassembly_duration.count();
      logger << "Constant flux term assembly duration : " << fassembly_time << "\n";
      auto rhs_end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration = rhs_end - rhs_start;
      rhs_time = duration.count();
      logger << "RHS total duration : " << duration.count() << "\n";
      if (penalty_time != 0)
        logger << "Penalty time in rhs : " << penalty_time / rhs_time * 100 << "%\n";
      else
        logger << "Weak Penalty time in rhs : " << wpenalty_time / rhs_time * 100 << "%\n";
      logger << "Constant source term assembly time in rhs : " << sassembly_time / rhs_time * 100 << "%\n"
             << "Constant flux term assembly time in rhs : " << fassembly_time / rhs_time * 100 << "%\n\n"
             << "-------------------------------------------------------------------------------------\n\n";
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeAreaQuad4(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  Real3 m3 = m_node_coord[cell.nodeId(3)];
  return 0.5 * ((m1.x * m2.y + m2.x * m3.y + m3.x * m0.y + m0.x * m1.y) - (m2.x * m1.y + m3.x * m2.y + m0.x * m3.y + m1.x * m0.y));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeAreaTriangle3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  return 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/*----------------------------#endif-----------------------------------------------*/

Real FemModule::
_computeEdgeLength2(Face face)
{
  Real3 m0 = m_node_coord[face.nodeId(0)];
  Real3 m1 = m_node_coord[face.nodeId(1)];
  return math::sqrt((m1.x - m0.x) * (m1.x - m0.x) + (m1.y - m0.y) * (m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real2 FemModule::
_computeEdgeNormal2(Face face)
{
  Real3 m0 = m_node_coord[face.nodeId(0)];
  Real3 m1 = m_node_coord[face.nodeId(1)];
  if (!face.isSubDomainBoundaryOutside())
    std::swap(m0, m1);
  Real2 N;
  Real norm_N = math::sqrt((m1.y - m0.y) * (m1.y - m0.y) + (m1.x - m0.x) * (m1.x - m0.x)); // for normalizing
  N.x = (m1.y - m0.y) / norm_N;
  N.y = (m0.x - m1.x) / norm_N;
  return N;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<3, 3> FemModule::
_computeElementMatrixTRIA3(Cell cell)
{
  // Get coordiantes of the triangle element  TRI3
  //------------------------------------------------
  //                  0 o
  //                   . .
  //                  .   .
  //                 .     .
  //              1 o . . . o 2
  //------------------------------------------------
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real area = _computeAreaTriangle3(cell); // calculate area

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  FixedMatrix<2, 3> b_matrix;
  b_matrix(0, 0) = dPhi0.x;
  b_matrix(0, 1) = dPhi1.x;
  b_matrix(0, 2) = dPhi2.x;

  b_matrix(1, 0) = dPhi0.y;
  b_matrix(1, 1) = dPhi1.y;
  b_matrix(1, 2) = dPhi2.y;

  b_matrix.multInPlace(1.0 / (2.0 * area));

  FixedMatrix<3, 3> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(area);

  //info() << "Cell=" << cell.localId();
  //std::cout << " int_cdPi_dPj=";
  //int_cdPi_dPj.dump(std::cout);
  //std::cout << "\n";

  return int_cdPi_dPj;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<4, 4> FemModule::
_computeElementMatrixQUAD4(Cell cell)
{
  // Get coordiantes of the quadrangular element  QUAD4
  //------------------------------------------------
  //             1 o . . . . o 0
  //               .         .
  //               .         .
  //               .         .
  //             2 o . . . . o 3
  //------------------------------------------------
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  Real3 m3 = m_node_coord[cell.nodeId(3)];

  Real area = _computeAreaQuad4(cell); // calculate area

  Real2 dPhi0(m2.y - m3.y, m3.x - m2.x);
  Real2 dPhi1(m3.y - m0.y, m0.x - m3.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);
  Real2 dPhi3(m1.y - m2.y, m2.x - m1.x);

  FixedMatrix<2, 4> b_matrix;
  b_matrix(0, 0) = dPhi0.x;
  b_matrix(0, 1) = dPhi1.x;
  b_matrix(0, 2) = dPhi2.x;
  b_matrix(0, 3) = dPhi3.x;

  b_matrix(1, 0) = dPhi0.y;
  b_matrix(1, 1) = dPhi1.y;
  b_matrix(1, 2) = dPhi2.y;
  b_matrix(1, 3) = dPhi3.y;

  b_matrix.multInPlace(1.0 / (2.0 * area));

  FixedMatrix<4, 4> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(area);

  //info() << "Cell=" << cell.localId();
  //std::cout << " int_cdPi_dPj=";
  //int_cdPi_dPj.dump(std::cout);
  //std::cout << "\n";

  return int_cdPi_dPj;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorQUAD4()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    if (cell.type() != IT_Quad4)
      ARCANE_FATAL("Only Quad4 cell type is supported");

    auto K_e = _computeElementMatrixQUAD4(cell); // element stifness matrix
    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e(n1_index, n2_index);
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (node1.isOwn()) {
          m_linear_system.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifdef ARCANE_HAS_CUDA

ARCCORE_HOST_DEVICE void FemModule::
_computeElementMatrixTRIA3GPU(CellLocalId icell, IndexedCellNodeConnectivityView cnc, ax::VariableNodeReal3InView in_node_coord, Real K_e[9])
{
  // Get coordiantes of the triangle element  TRI3
  //------------------------------------------------
  //                  0 o
  //                   . .
  //                  .   .
  //                 .     .
  //              1 o . . . o 2
  //------------------------------------------------
  Real3 m0 = in_node_coord[cnc.nodeId(icell, 0)];
  Real3 m1 = in_node_coord[cnc.nodeId(icell, 1)];
  Real3 m2 = in_node_coord[cnc.nodeId(icell, 2)];

  Real area = 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y)); // calculate area

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  //We will want to replace fixed matrix by some numarray ? Will not work because NumArray function are host functions
  //NumArray<Real, ExtentsV<2, 3>> b_matrix(eMemoryRessource::Device);
  Real b_matrix[2][3] = { 0 };
  b_matrix[0][0] = dPhi0.x * (1.0 / (2.0 * area));
  b_matrix[0][1] = dPhi1.x * (1.0 / (2.0 * area));
  b_matrix[0][2] = dPhi2.x * (1.0 / (2.0 * area));

  b_matrix[1][0] = dPhi0.y * (1.0 / (2.0 * area));
  b_matrix[1][1] = dPhi1.y * (1.0 / (2.0 * area));
  b_matrix[1][2] = dPhi2.y * (1.0 / (2.0 * area));

  //NumArray<Real, ExtentsV<3, 3>> int_cdPi_dPj;

  //Multiplying b_matrix by its transpose, and doing the mult in place in the same loop
  for (Int32 i = 0; i < 3; i++) {
    for (Int32 j = 0; j < 3; j++) {
      Real x = 0.0;
      for (Int32 k = 0; k < 2; k++) {
        x += b_matrix[k][i] * b_matrix[k][j];
      }
      K_e[i * 3 + j] = x * area;
    }
  }

  //info() << "Cell=" << cell.localId();
  //std::cout << " int_cdPi_dPj=";
  //int_cdPi_dPj.dump(std::cout);
  //std::cout << "\n";

  //No need to return anymore
  //return int_cdPi_dPj;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* USE_CSR_GPU methods                                                       */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Initialization of the coo matrix. It only works for p=1 since there is
 * one node per Edge. Currently, there is no difference between buildMatrixCsr and this method.
 * 
 * 
 */
void FemModule::
_buildMatrixCsrGPU()
{
  //Initialization of the csr matrix;
  //This formula only works in p=1

  /*
  //Create a connection between nodes through the faces
  //Useless here because we only need this information once
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  auto idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  auto* cn = idx_cn->connectivity();
  ENUMERATE_NODE (node, allNodes()) {
  }
  */

  Int32 nnz = nbFace() * 2 + nbNode();
  m_csr_matrix.initialize(m_dof_family, nnz, nbNode());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  //We iterate through the node, and we do not sort anymore : we assume the nodes ID are sorted, and we will iterate throught the column to avoid making < and > comparison
  ENUMERATE_NODE (inode, allNodes()) {
    Node node = *inode;

    m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

    for (Face face : node.faces()) {
      if (face.nodeId(0) == node.localId())
        m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
      else
        m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrGPUBilinearOperatorTRIA3()
{

  Timer::Action timer_gpu_bili(this->subDomain(), "CsrGpuBilinearOperatorTria3");

  std::chrono::_V2::system_clock::time_point lhs_start;
  double global_build_average = 0;
  double build_time = 0;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using GPU csr with NumArray format\n";
    lhs_start = std::chrono::high_resolution_clock::now();
  }
  {
    Timer::Action timer_gpu_build(this->subDomain(), "CsrGpuBuildMatrix");
    // Build the csr matrix
    _buildMatrixCsrGPU();
  }

  std::chrono::_V2::system_clock::time_point build_stop;
  std::chrono::_V2::system_clock::time_point var_init_start;
  if (m_register_time) {
    build_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> build_duration = build_stop - lhs_start;
    build_time = build_duration.count();
    var_init_start = std::chrono::high_resolution_clock::now();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();
  // Boucle sur les mailles déportée sur accélérateur
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
  Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
  auto in_out_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);
  UnstructuredMeshConnectivityView m_connectivity_view;
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_connectivity_view.setMesh(this->mesh());
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
  Arcane::ItemGenericInfoListView cells_infos(this->mesh()->cellFamily());

  std::chrono::_V2::system_clock::time_point var_init_stop;
  std::chrono::_V2::system_clock::time_point loop_start;
  double var_init_time = 0;
  if (m_register_time) {
    var_init_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> var_init_duration = var_init_stop - var_init_start;
    var_init_time = var_init_duration.count();
    loop_start = std::chrono::high_resolution_clock::now();
  }

  Timer::Action timer_add_compute(this->subDomain(), "CsrGpuAddComputeLoop");

  command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
  {

    Real K_e[9] = { 0 };

    _computeElementMatrixTRIA3GPU(icell, cnc, in_node_coord, K_e); // element stifness matrix

    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
    Int32 n1_index = 0;
    for (NodeLocalId node1 : cnc.nodes(icell)) {
      Int32 n2_index = 0;
      for (NodeLocalId node2 : cnc.nodes(icell)) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        double v = K_e[n1_index * 3 + n2_index];
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (nodes_infos.isOwn(node1)) {
          Int32 row = node_dof.dofId(node1, 0).localId();
          Int32 col = node_dof.dofId(node2, 0).localId();
          Int32 begin = in_row_csr[row];
          Int32 end;
          if (row == row_csr_size - 1) {
            end = col_csr_size;
          }
          else {
            end = in_row_csr[row + 1];
          }
          while (begin < end) {
            if (in_col_csr[begin] == col) {
              // t is necessary to get the right type for the atomicAdd (but that means that we have more operations ?)
              // The Macro is there to avoid compilation error if not in c++ 20
              ax::doAtomic<ax::eAtomicOperation::Add>(in_out_val_csr(begin), v);
              break;
            }
            begin++;
          }
        }
        ++n2_index;
      }
      ++n1_index;
    }
  };
  if (m_register_time) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - lhs_start;
    std::chrono::duration<double> loop_duration = lhs_end - loop_start;

    double loop_time = loop_duration.count();
    double lhs_loc_time = duration.count();
    logger << "Building time of the coo matrix :" << build_time << "\n"
           << "Variable initialisation time : " << var_init_time << "\n"
           << "Computation and Addition time : " << loop_time << "\n"
           << "LHS Total time : " << lhs_loc_time << "\n"
           << "Build matrix time in lhs :" << build_time / lhs_loc_time * 100 << "%\n"
           << "Variable initialisation time in lhs : " << var_init_time / lhs_loc_time * 100 << "%\n"
           << "Computation and Addition time in lhs : " << loop_time / lhs_loc_time * 100 << "%\n\n"
           << "-------------------------------------------------------------------------------------\n\n";
    lhs_time += lhs_loc_time;
    wbuild << lhs_loc_time << ",";
    timer << loop_time << ",";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
#ifdef USE_COO_GPU
//Currently, this code does not work
/**
 * @brief Initialization of the coo matrix. It only works for p=1 since there is
 * one node per Edge. There is no difference with buildMatrix() and this method currently.
 * 
 * 
 */
void FemModule::
_buildMatrixGPU()
{
  //Initialization of the coo matrix;
  //This formula only works in p=1

  /*
  //Create a connection between nodes through the faces
  //Useless here because we only need this information once
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  auto idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  auto* cn = idx_cn->connectivity();
  ENUMERATE_NODE (node, allNodes()) {
  }
  */

  Int32 nnz = nbFace() * 2 + nbNode();
  m_coo_matrix.initialize(m_dof_family, nnz);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  //We iterate through the node, and we do not sort anymore : we assume the nodes ID are sorted, and we will iterate throught the column to avoid making < and > comparison
  ENUMERATE_NODE (inode, allNodes()) {
    Node node = *inode;

    m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

    for (Face face : node.faces()) {
      if (face.nodeId(0) == node.localId())
        m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
      else
        m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
    }
  }

  //In this commented code, we begin by filling the diagonal before filling what's left by iterating through the nodes. It corresponds to the COO-sort method in the diagrams
  /*
  //Fill the diagonal
  ENUMERATE_NODE (inode, allNodes()) {
    Node node = *inode;
    m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));
  }

  //Fill what is left
  ENUMERATE_FACE (iface, allFaces()) {
    Face face = *iface;
    auto nodes = face.nodes();
    for (Int32 i = 0; i < nodes.size() - i - 1; i++) {
      m_coo_matrix.setCoordinates(node_dof.dofId(nodes[i], 0), node_dof.dofId(nodes[nodes.size() - i - 1], 0));
      m_coo_matrix.setCoordinates(node_dof.dofId(nodes[nodes.size() - i - 1], 0), node_dof.dofId(nodes[i], 0));
    }
  }

  //Sort the matrix
  m_coo_matrix.sort();
  */
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCooGPUBilinearOperatorTRIA3()
{
#if defined(m_register_time) && defined(m_cache_warming)
  logger << "-------------------------------------------------------------------------------------\n"
         << "Using CPU coo with NumArray format\n";
  auto lhs_start = std::chrono::high_resolution_clock::now();
  double compute_average = 0;
  double global_build_average = 0;
  double build_time = 0;
#endif
  // Build the coo matrix
  _buildMatrixGPU();
#if defined(m_register_time) && defined(m_cache_warming)
  auto build_stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> build_duration = build_stop - lhs_start;
  build_time = build_duration.count();
#endif

  RunQueue* queue = acceleratorMng()->defaultQueue();
  // Boucle sur les mailles déportée sur accélérateur
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_coo = ax::viewIn(command, m_coo_matrix.m_matrix_row);
  auto in_col_coo = ax::viewIn(command, m_coo_matrix.m_matrix_column);
  auto in_out_val_coo = ax::viewInOut(command, m_coo_matrix.m_matrix_value);
  UnstructuredMeshConnectivityView m_connectivity_view;
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_connectivity_view.setMesh(this->mesh());
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
  Arcane::ItemGenericInfoListView cells_infos(this->mesh()->cellFamily());

  command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
  {

    Real K_e[9] = { 0 };

    _computeElementMatrixTRIA3GPU(icell, cnc, in_node_coord, K_e); // element stifness matrix

    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
    Int32 n1_index = 0;
    for (NodeLocalId node1 : cnc.nodes(icell)) {
      Int32 n2_index = 0;
      for (NodeLocalId node2 : cnc.nodes(icell)) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e[n1_index * 3 + n2_index];
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        //replacing the isOwn (probably with a nice view)
        if (nodes_infos.isOwn(node1)) {
          //m_coo_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }
  };

#if defined(m_register_time) && defined(m_cache_warming)
  if (i == 3) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - lhs_start;
    double lhs_loc_time = duration.count();
    logger << "Building time of the coo matrix :" << build_time << "\n"
           << "Compute Elements average time : " << compute_average / nbCell() << "\n"
           << "Compute Elements total time : " << compute_average << "\n"
           << "Add in global matrix average time : " << global_build_average / nbCell() << "\n"
           << "Add in global matrix total time : " << global_build_average << "\n"
           << "LHS Total time : " << lhs_loc_time << "\n"
           << "Build matrix time in lhs :" << build_time / lhs_loc_time * 100 << "%\n"
           << "Compute element time in lhs : " << compute_average / lhs_loc_time * 100 << "%\n"
           << "Add in global matrix time in lhs : " << global_build_average / lhs_loc_time * 100 << "%\n\n"
           << "-------------------------------------------------------------------------------------\n\n";
    lhs_time += lhs_loc_time;
  }
#endif
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* USE_CSR methods                                                           */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Initialization of the coo matrix. It only works for p=1 since there is
 * one node per Edge.
 * 
 * 
 */
void FemModule::
_buildMatrixCsr()
{
  //Initialization of the coo matrix;
  //This formula only works in p=1

  /*
  //Create a connection between nodes through the faces
  //Useless here because we only need this information once
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  auto idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  auto* cn = idx_cn->connectivity();
  ENUMERATE_NODE (node, allNodes()) {
  }
  */

  Int32 nnz = nbFace() * 2 + nbNode();
  m_csr_matrix.initialize(m_dof_family, nnz, nbNode());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  //We iterate through the node, and we do not sort anymore : we assume the nodes ID are sorted, and we will iterate throught the column to avoid making < and > comparison
  ENUMERATE_NODE (inode, allNodes()) {
    Node node = *inode;

    m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

    for (Face face : node.faces()) {
      if (face.nodeId(0) == node.localId())
        m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
      else
        m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrBilinearOperatorTRIA3()
{

  Timer::Action timer_csr_bili(this->subDomain(), "AssembleCsrBilinearOperatorTria3");

  double compute_average = 0;
  double global_build_average = 0;
  double build_time = 0;
  std::chrono::_V2::system_clock::time_point lhs_start;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using CPU CSR with NumArray format\n";
    lhs_start = std::chrono::high_resolution_clock::now();
  }
  {
    Timer::Action timer_csr_build(this->subDomain(), "CsrBuildMatrix");
    // Build the csr matrix
    _buildMatrixCsr();
  }

  std::chrono::_V2::system_clock::time_point build_stop;
  if (m_register_time) {
    build_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> build_duration = build_stop - lhs_start;
    build_time = build_duration.count();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    std::chrono::_V2::system_clock::time_point compute_El_start;
    if (m_register_time)
      compute_El_start = std::chrono::high_resolution_clock::now();

    FixedMatrix<3, 3> K_e;
    {
      Timer::Action timer_csr_compute_add(this->subDomain(), "CsrComputeElementMatrixTria3");
      K_e = _computeElementMatrixTRIA3(cell); // element stifness matrix
    }

    if (m_register_time) {
      auto compute_El_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> compute_duration = compute_El_stop - compute_El_start;
      compute_average += compute_duration.count();
    }

    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
    std::chrono::_V2::system_clock::time_point global_build_start;
    if (m_register_time) {
      global_build_start = std::chrono::high_resolution_clock::now();
    }
    Timer::Action timer_action(this->subDomain(), "CsrAddToGlobalMatrix");
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e(n1_index, n2_index);
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (node1.isOwn()) {
          m_csr_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }

    if (m_register_time) {
      auto global_build_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> global_build_duration = global_build_stop - global_build_start;
      global_build_average += global_build_duration.count();
    }
  }

  if (m_register_time) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - lhs_start;
    double lhs_loc_time = duration.count();
    logger << "Building time of the csr matrix :" << build_time << "\n"
           << "Compute Elements average time : " << compute_average / nbCell() << "\n"
           << "Compute Elements total time : " << compute_average << "\n"
           << "Add in global matrix average time : " << global_build_average / nbCell() << "\n"
           << "Add in global matrix total time : " << global_build_average << "\n"
           << "LHS Total time : " << lhs_loc_time << "\n"
           << "Build matrix time in lhs :" << build_time / lhs_loc_time * 100 << "%\n"
           << "Compute element time in lhs : " << compute_average / lhs_loc_time * 100 << "%\n"
           << "Add in global matrix time in lhs : " << global_build_average / lhs_loc_time * 100 << "%\n\n"
           << "-------------------------------------------------------------------------------------\n\n";
    lhs_time += lhs_loc_time;
    wbuild << lhs_loc_time << ",";
    timer << compute_average + global_build_average << ",";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* USE_COO_SORT methods                                                      */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Initialization of the coo matrix. It only works for p=1 since there is
 * one node per Edge.
 * 
 * 
 */
void FemModule::
_buildMatrixSort()
{
  //Initialization of the coo matrix;
  //This formula only works in p=1

  /*
  //Create a connection between nodes through the faces
  //Useless here because we only need this information once
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  auto idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  auto* cn = idx_cn->connectivity();
  ENUMERATE_NODE (node, allNodes()) {
  }
  */

  Int32 nnz = nbFace() * 2 + nbNode();
  m_coo_matrix.initialize(m_dof_family, nnz);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  //In this commented code, we begin by filling the diagonal before filling what's left by iterating through the nodes. It corresponds to the COO-sort method in the diagrams

  //In this one, we begin by filling the diagonal before filling what's left by iterating through the nodes

  //Fill the diagonal
  ENUMERATE_NODE (inode, allNodes()) {
    Node node = *inode;
    m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));
  }

  //Fill what is left
  ENUMERATE_FACE (iface, allFaces()) {
    Face face = *iface;
    auto nodes = face.nodes();
    for (Int32 i = 0; i < nodes.size() - i - 1; i++) {
      m_coo_matrix.setCoordinates(node_dof.dofId(nodes[i], 0), node_dof.dofId(nodes[nodes.size() - i - 1], 0));
      m_coo_matrix.setCoordinates(node_dof.dofId(nodes[nodes.size() - i - 1], 0), node_dof.dofId(nodes[i], 0));
    }
  }

  //Sort the matrix
  m_coo_matrix.sort();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCooSortBilinearOperatorTRIA3()
{

  Timer::Action timer_coosort_bili(this->subDomain(), "AssembleCooSortBilinearOperatorTria3");

  double compute_average = 0;
  double global_build_average = 0;
  double build_time = 0;
  std::chrono::_V2::system_clock::time_point lhs_start;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using CPU coo with NumArray format\n";
    lhs_start = std::chrono::high_resolution_clock::now();
  }
  {
    Timer::Action timer_build_coosort(this->subDomain(), "BuildMatrixCooSort");
    // Build the coo matrix
    _buildMatrixSort();
  }
  if (m_register_time) {
    auto build_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> build_duration = build_stop - lhs_start;
    build_time = build_duration.count();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    std::chrono::_V2::system_clock::time_point compute_El_start;
    if (m_register_time) {
      compute_El_start = std::chrono::high_resolution_clock::now();
    }

    FixedMatrix<3, 3> K_e;
    {
      Timer::Action timer_element_coosort(this->subDomain(), "CooSortComputeElementMatrixTria3");
      K_e = _computeElementMatrixTRIA3(cell); // element stifness matrix
    }

    if (m_register_time) {
      auto compute_El_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> compute_duration = compute_El_stop - compute_El_start;
      compute_average += compute_duration.count();
    }

    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
    std::chrono::_V2::system_clock::time_point global_build_start;
    if (m_register_time) {
      global_build_start = std::chrono::high_resolution_clock::now();
    }

    Timer::Action timer_coosort_add_compute(this->subDomain(), "CooSortAddToGlobalMatrix");
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e(n1_index, n2_index);
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (node1.isOwn()) {
          m_coo_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }
    if (m_register_time) {
      auto global_build_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> global_build_duration = global_build_stop - global_build_start;
      global_build_average += global_build_duration.count();
    }
  }

  if (m_register_time) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - lhs_start;
    double lhs_loc_time = duration.count();
    logger << "Building time of the coo matrix :" << build_time << "\n"
           << "Compute Elements average time : " << compute_average / nbCell() << "\n"
           << "Compute Elements total time : " << compute_average << "\n"
           << "Add in global matrix average time : " << global_build_average / nbCell() << "\n"
           << "Add in global matrix total time : " << global_build_average << "\n"
           << "LHS Total time : " << lhs_loc_time << "\n"
           << "Build matrix time in lhs :" << build_time / lhs_loc_time * 100 << "%\n"
           << "Compute element time in lhs : " << compute_average / lhs_loc_time * 100 << "%\n"
           << "Add in global matrix time in lhs : " << global_build_average / lhs_loc_time * 100 << "%\n\n"
           << "-------------------------------------------------------------------------------------\n\n";
    lhs_time += lhs_loc_time;
    wbuild << lhs_loc_time << ",";
    timer << compute_average + global_build_average << ",";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* USE_COO methods                                                           */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Initialization of the coo matrix. It only works for p=1 since there is
 * one node per Edge.
 * 
 * 
 */
void FemModule::
_buildMatrix()
{
  //Initialization of the coo matrix;
  //This formula only works in p=1

  /*
  //Create a connection between nodes through the faces
  //Useless here because we only need this information once
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  auto idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  auto* cn = idx_cn->connectivity();
  ENUMERATE_NODE (node, allNodes()) {
  }
  */

  Int32 nnz = nbFace() * 2 + nbNode();
  m_coo_matrix.initialize(m_dof_family, nnz);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  //We iterate through the node, and we do not sort anymore : we assume the nodes ID are sorted, and we will iterate throught the column to avoid making < and > comparison
  ENUMERATE_NODE (inode, allNodes()) {
    Node node = *inode;

    m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

    for (Face face : node.faces()) {
      if (face.nodeId(0) == node.localId())
        m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
      else
        m_coo_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCooBilinearOperatorTRIA3()
{

  Timer::Action timer_coo_bili(this->subDomain(), "AssembleCooBilinearOperatorTria3");

  double compute_average = 0;
  double global_build_average = 0;
  double build_time = 0;
  std::chrono::_V2::system_clock::time_point lhs_start;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using CPU coo with NumArray format\n";
    lhs_start = std::chrono::high_resolution_clock::now();
  }
  {
    Timer::Action timer_coo_build(this->subDomain(), "CooBuildMatrix");
    // Build the coo matrix
    _buildMatrix();
  }
  if (m_register_time) {
    auto build_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> build_duration = build_stop - lhs_start;
    build_time = build_duration.count();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    std::chrono::_V2::system_clock::time_point compute_El_start;
    if (m_register_time) {
      compute_El_start = std::chrono::high_resolution_clock::now();
    }

    FixedMatrix<3, 3> K_e;
    {
      Timer::Action timer_coo_compute(this->subDomain(), "CooComputeElementMatrix");
      K_e = _computeElementMatrixTRIA3(cell); // element stifness matrix
    }

    if (m_register_time) {
      auto compute_El_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> compute_duration = compute_El_stop - compute_El_start;
      compute_average += compute_duration.count();
    }

    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]

    std::chrono::_V2::system_clock::time_point global_build_start;
    if (m_register_time) {
      auto global_build_start = std::chrono::high_resolution_clock::now();
    }
    Timer::Action timer_coo_add(this->subDomain(), "CooAddToGlobalMatrix");
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e(n1_index, n2_index);
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (node1.isOwn()) {
          m_coo_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }
    if (m_register_time) {
      auto global_build_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> global_build_duration = global_build_stop - global_build_start;
      global_build_average += global_build_duration.count();
    }
  }

  if (m_register_time) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - lhs_start;
    double lhs_loc_time = duration.count();
    logger << "Building time of the coo matrix :" << build_time << "\n"
           << "Compute Elements average time : " << compute_average / nbCell() << "\n"
           << "Compute Elements total time : " << compute_average << "\n"
           << "Add in global matrix average time : " << global_build_average / nbCell() << "\n"
           << "Add in global matrix total time : " << global_build_average << "\n"
           << "LHS Total time : " << lhs_loc_time << "\n"
           << "Build matrix time in lhs :" << build_time / lhs_loc_time * 100 << "%\n"
           << "Compute element time in lhs : " << compute_average / lhs_loc_time * 100 << "%\n"
           << "Add in global matrix time in lhs : " << global_build_average / lhs_loc_time * 100 << "%\n\n"
           << "-------------------------------------------------------------------------------------\n\n";
    lhs_time += lhs_loc_time;
    wbuild << lhs_loc_time << ",";
    timer << compute_average + global_build_average << ",";
  }
}

#ifdef USE_CUSPARSE_ADD

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* USE_CUSPARSE_ADD methods                                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
printCsrMatrix(std::string fileName, cusparseCsr csr, bool is_coo)
{
  ofstream file(fileName);
  file << "size :" << csr.nnz << "\n";
  for (auto i = 0; i < (is_coo ? csr.nnz : nbNode()); i++) {
    file << csr.csrRow[i] << " ";
  }
  file << "\n";
  for (auto i = 0; i < csr.nnz; i++) {
    file << csr.csrCol[i] << " ";
  }
  file << "\n";
  for (auto i = 0; i < csr.nnz; i++) {
    file << csr.csrVal[i] << " ";
  }
  file.close();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_computeCusparseElementMatrix(cusparseCsr& result, cusparseCsr& global, Cell cell, cusparseHandle_t handle, IndexedNodeDoFConnectivityView node_dof,
                              computeTimer& timer)
{

  Timer::Action timer_action(this->subDomain(), "ComputeCusparseElementMatrix");

  std::chrono::_V2::system_clock::time_point compute_start;
  if (m_register_time) {
    compute_start = std::chrono::high_resolution_clock::now();
  }
  /*-------------------------------------------------------------------------------------------------------------------------------*/
  //First part : compute element matrix
  // Get coordiantes of the triangle element  TRI3
  //------------------------------------------------
  //                  0 o
  //                   . .
  //                  .   .
  //                 .     .
  //              1 o . . . o 2
  //------------------------------------------------
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real area = _computeAreaTriangle3(cell); // calculate area

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  FixedMatrix<2, 3> b_matrix;
  b_matrix(0, 0) = dPhi0.x;
  b_matrix(0, 1) = dPhi1.x;
  b_matrix(0, 2) = dPhi2.x;

  b_matrix(1, 0) = dPhi0.y;
  b_matrix(1, 1) = dPhi1.y;
  b_matrix(1, 2) = dPhi2.y;

  b_matrix.multInPlace(1.0 / (2.0 * area));

  FixedMatrix<3, 3> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(area);

  std::chrono::_V2::system_clock::time_point convert_coo;
  if (m_register_time) {
    auto compute_el_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> compute_el_duration = compute_el_stop - compute_start;
    timer.compute_el += compute_el_duration.count();
    convert_coo = std::chrono::high_resolution_clock::now();
  }
  /*-------------------------------------------------------------------------------------------------------------------------------*/
  //Second part : putting the matrix in COO format (might want to optimsie that part by doing it earlier) before converting it to csr

  //Must change int_cdPi_dPj in a COO matrix (before converting it to csr);
  void* row_void;
  CHECK_CUDA(cudaMallocManaged(&row_void, 9 * sizeof(Int32)));
  Int32* row_indexes = (Int32*)row_void;
  void* col_void;
  CHECK_CUDA(cudaMallocManaged(&col_void, 9 * sizeof(Int32)));
  Int32* col_indexes = (Int32*)col_void;
  void* vals_void;
  CHECK_CUDA(cudaMallocManaged(&vals_void, 9 * sizeof(float)));
  float* vals = (float*)vals_void;

  cusparseMatDescr_t local_mat;
  CHECK_CUSPARSE(cusparseCreateMatDescr(&local_mat));
  cusparseCsr local;
  local.desc = local_mat;
  local.csrRow = row_indexes;
  local.csrCol = col_indexes;
  local.csrVal = (float*)vals;
  local.nnz = 9;

  int i = 0;
  int j = 0;
  for (NodeLocalId node1 : cell.nodes()) {
    j = 0;
    for (NodeLocalId node2 : cell.nodes()) {
      vals[i * 3 + j] = int_cdPi_dPj(i, j);
      row_indexes[i * 3 + j] = node_dof.dofId(node1, 0);
      col_indexes[i * 3 + j] = node_dof.dofId(node2, 0);
      j++;
    }
    i++;
  }

  std::chrono::_V2::system_clock::time_point sort_coo;
  std::chrono::duration<double> convert_coo_time;
  if (m_register_time) {
    auto convert_coo_stop = std::chrono::high_resolution_clock::now();
    convert_coo_time = convert_coo_stop - convert_coo;
    timer.convert_coo += convert_coo_time.count();
    sort_coo = std::chrono::high_resolution_clock::now();
  }
  //Sorting of the COO values with an insertion sort
  Int32 rj = 0;
  Int32 cj = 0;
  float vj = 0;
  for (i = 1; i < 9; i++) {
    rj = row_indexes[i];
    cj = col_indexes[i];
    vj = vals[i];
    j = i - 1;
    while (j >= 0 && row_indexes[j] > rj) {
      row_indexes[j + 1] = row_indexes[j];
      col_indexes[j + 1] = col_indexes[j];
      vals[j + 1] = vals[j];
      j--;
    }
    row_indexes[j + 1] = rj;
    col_indexes[j + 1] = cj;
    vals[j + 1] = vj;
    Int32 k = j - 1;
    Int32 rk, ck;
    float vk;
    if (j > 0) {
      rk = row_indexes[j];
      ck = col_indexes[j];
      vk = vals[j];
      while (k >= 0 && row_indexes[k] == rk && col_indexes[k] > ck) {
        col_indexes[k + 1] = col_indexes[k];
        vals[k + 1] = vals[k];
        k--;
      }
      col_indexes[k + 1] = ck;
      vals[k + 1] = vk;
    }
  }
  std::chrono::_V2::system_clock::time_point convert_csr;
  if (m_register_time) {
    auto sort_coo_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sort_coo_time = sort_coo_stop - sort_coo;
    std::chrono::duration<double> convert_coo_tot = sort_coo_stop - convert_coo;
    timer.sort_coo += sort_coo_time.count();
    timer.convert_coo_tot += convert_coo_tot.count();
    convert_csr = std::chrono::high_resolution_clock::now();
  }

  //conversion from COO to CSR
  void* csrRowPtr_void;
  CHECK_CUDA(cudaMallocManaged(&csrRowPtr_void, nbNode() * sizeof(Int32)));
  Int32* csrRowPtr = (Int32*)csrRowPtr_void;
  CHECK_CUSPARSE(cusparseXcoo2csr(handle, row_indexes, 9, nbNode(), csrRowPtr, CUSPARSE_INDEX_BASE_ZERO));
  local.csrRow = csrRowPtr;
  CHECK_CUDA(cudaFree(row_indexes));

  std::chrono::_V2::system_clock::time_point adding_global;
  if (m_register_time) {
    auto convert_csr_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> convert_csr_time = convert_csr_stop - convert_csr;
    std::chrono::duration<double> convert_tot = convert_csr_stop - convert_coo;
    timer.convert_csr_tot += convert_coo_time.count();
    timer.convert_tot += convert_tot.count();
    adding_global = std::chrono::high_resolution_clock::now();
  }
  /*-------------------------------------------------------------------------------------------------------------------------------*/
  // Third part : adding the local and global, storing result in the res

  //Adding the CSR local matrix to the global one using cusparsecsrgeam2
  //see https://docs.nvidia.com/cuda/cusparse/index.html?highlight=cusparseScsrgeam#cusparse-t-csrgeam2 for the example code
  Int32 baseC,
  nnzC;
  size_t bufferSizeInBytes;
  char* buffer = NULL;
  void* buffer_void = NULL;
  Int32* nnzTotalDevHostPtr = &nnzC;
  CHECK_CUSPARSE(cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST));
  //Int32* csrRowPtrC;
  //CHECK_CUDA(cudaMallocManaged(&csrRowPtrC, nbNode() + 1 * sizeof(Int32)));
  float alpha = 1.0;
  float beta = 1.0;
  Int32 m = nbNode();
  CHECK_CUSPARSE(cusparseScsrgeam2_bufferSizeExt(handle, m, m,
                                                 &alpha,
                                                 global.desc, global.nnz,
                                                 global.csrVal, global.csrRow, global.csrCol,
                                                 &beta,
                                                 local.desc, local.nnz,
                                                 local.csrVal, local.csrRow, local.csrCol,
                                                 result.desc,
                                                 result.csrVal, result.csrRow, result.csrCol, &bufferSizeInBytes));
  CHECK_CUDA(cudaMallocManaged(&buffer_void, bufferSizeInBytes * sizeof(char)));
  buffer = (char*)buffer_void;

  CHECK_CUSPARSE(cusparseXcsrgeam2Nnz(handle, m, m,
                                      local.desc, local.nnz, local.csrRow, local.csrCol,
                                      global.desc, global.nnz, global.csrRow, global.csrCol,
                                      result.desc, result.csrRow, nnzTotalDevHostPtr,
                                      buffer));
  if (NULL != nnzTotalDevHostPtr)
    nnzC = *nnzTotalDevHostPtr;
  else {
    CHECK_CUDA(cudaMemcpy(&nnzC, result.csrRow + m, sizeof(int), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(&baseC, result.csrRow, sizeof(int), cudaMemcpyDeviceToHost));
    nnzC -= baseC;
  }
  result.nnz = nnzC;
  void* res_col_void;
  void* res_val_void;

  CHECK_CUDA(cudaMallocManaged(&res_col_void, sizeof(Int32) * nnzC));
  CHECK_CUDA(cudaMallocManaged(&res_val_void, sizeof(float) * nnzC));
  result.csrCol = (Int32*)res_col_void;
  result.csrVal = (float*)res_val_void;
  CHECK_CUSPARSE(cusparseScsrgeam2(handle, m, m,
                                   &alpha,
                                   local.desc, local.nnz,
                                   local.csrVal, local.csrRow, local.csrCol,
                                   &beta,
                                   global.desc, global.nnz,
                                   global.csrVal, global.csrRow, global.csrCol,
                                   result.desc,
                                   result.csrVal, result.csrRow, result.csrCol,
                                   buffer));

  CHECK_CUDA(cudaFree(local.csrVal));
  CHECK_CUDA(cudaFree(local.csrCol));
  CHECK_CUDA(cudaFree(local.csrRow));
  CHECK_CUSPARSE(cusparseDestroyMatDescr(local.desc));

  if (m_register_time) {
    auto adding_global_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> adding_tot = adding_global_stop - adding_global;
    std::chrono::duration<double> compute_tot = adding_global_stop - compute_start;
    timer.add_glob += adding_tot.count();
    timer.compute_tot += compute_tot.count();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Assemble Bilinear TRIA3 with cusparse help. It only works for p=1 since there is
 * one node per Edge.
 * 
 * 
 */
void FemModule::
_assembleCusparseBilinearOperatorTRIA3()
{

  Timer::Action timer_action(this->subDomain(), "AssembleCusparseBilinearOperator");

  computeTimer t = {};
  std::chrono::_V2::system_clock::time_point lhs_s;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using Cusparse for Bilinear assembly\n";
    lhs_s = std::chrono::high_resolution_clock::now();
  }
  //Initialization of the CSR matrix;
  //This formula only works in p=1
  std::chrono::_V2::system_clock::time_point cuda_init_start;
  if (m_register_time) {
    cuda_init_start = std::chrono::high_resolution_clock::now();
  }
  CHECK_CUDA(cudaFree(0));

  double cuda_init = 0;
  std::chrono::_V2::system_clock::time_point cuda_init_stop;
  if (m_register_time) {
    cuda_init_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> cuda_init_time = cuda_init_stop - cuda_init_start;
    cuda_init = cuda_init_time.count();
  }

  Int32 nnz = nbFace() * 2 + nbNode();
  //Initialize the global matrix. Everything is in the unified memory
  void* res1_row_void;
  void* res2_row_void;
  CHECK_CUDA(cudaMallocManaged(&res1_row_void, sizeof(Int32) * nbNode()));
  Int32* res1_row = (Int32*)res1_row_void;
  CHECK_CUDA(cudaMallocManaged(&res2_row_void, sizeof(Int32) * nbNode()));
  Int32* res2_row = (Int32*)res2_row_void;

  cusparseHandle_t handle;
  CHECK_CUSPARSE(cusparseCreate(&handle));
  //The number of Node must be changed when p != 1

  //init result matrix
  cusparseCsr res1;
  cusparseCsr res2;

  cusparseMatDescr_t res1_desc;
  cusparseMatDescr_t res2_desc;
  CHECK_CUSPARSE(cusparseCreateMatDescr(&res1_desc));
  CHECK_CUSPARSE(cusparseCreateMatDescr(&res2_desc));

  res1.desc = res1_desc;
  res2.desc = res2_desc;
  res1.csrRow = res1_row;
  res2.csrRow = res2_row;
  res1.csrCol = NULL;
  res2.csrCol = NULL;

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Int32 i = 0;

  ENUMERATE_CELL (icell, allCells()) {
    Cell cell = *icell;

    if (i % 2 == 0)
      //computation of the local matrix and adding it in the global one
      _computeCusparseElementMatrix(res1, res2, cell, handle, node_dof,
                                    t);
    else
      _computeCusparseElementMatrix(res2, res1, cell, handle, node_dof,
                                    t);
    i++;
  }
  /*
  if (nbNode() % 2 == 0)
    printCsrMatrix("csrTest.txt", res1, false);
  else
    printCsrMatrix("csrTest.txt", res2, false);

*/
  CHECK_CUSPARSE(cusparseDestroyMatDescr(res1.desc));
  CHECK_CUSPARSE(cusparseDestroyMatDescr(res2.desc));

  //Must free the resulting vectors, not done there yet
  if (m_register_time) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - cuda_init_stop;
    double lhs_loc_time = duration.count();
    logger << "Average Time to compute element matrix : " << t.compute_el / nbCell() << "\n"
           << "Total Time to compute element matrix : " << t.compute_el << "\n"
           << "Percentage time to compute element matrix : " << t.compute_el / lhs_loc_time * 100 << "%\n"
           << "Average Time to convert to coo : " << t.convert_coo / nbCell() << "\n"
           << "Total Time to convert to coo : " << t.convert_coo << "\n"
           << "Percentage Time to convert to coo : " << t.convert_coo / lhs_loc_time * 100 << "%\n"
           << "Average Time to sort the coo : " << t.sort_coo / nbCell() << "\n"
           << "Total Time to sort the coo : " << t.sort_coo << "\n"
           << "Percentage Time to sort the coo : " << t.sort_coo / lhs_loc_time * 100 << "%\n"
           << "Average Time to convert and sort to coo : " << t.convert_coo_tot / nbCell() << "\n"
           << "Total Time to convert and sort to coo : " << t.convert_coo_tot << "\n"
           << "Percentage Time to convert and sort to coo : " << t.convert_coo_tot / lhs_loc_time * 100 << "%\n"
           << "Average Time to convert to csr : " << t.convert_csr_tot / nbCell() << "\n"
           << "Total Time to convert to csr : " << t.convert_csr_tot << "\n"
           << "Percentage Time to convert to csr : " << t.convert_csr_tot / lhs_loc_time * 100 << "%\n"
           << "Average Time to convert the computed matrix : " << t.convert_tot / nbCell() << "\n"
           << "Total Time to convert the computed matrix : " << t.convert_tot << "\n"
           << "Percentage Time to convert the computed matrix : " << t.convert_tot / lhs_loc_time * 100 << "%\n"
           << "Average Time to add to the global matrix : " << t.add_glob / nbCell() << "\n"
           << "Total Time to add to the global matrix : " << t.add_glob << "\n"
           << "Percentage Time to add to the global matrix : " << t.add_glob / lhs_loc_time * 100 << "%\n"
           << "Average Time to make the computation operation : " << t.compute_tot / nbCell() << "\n"
           << "Total Time to make the computation operation : " << t.compute_tot << "\n"
           << "Percentage Time to make the computation operation : " << t.compute_tot / lhs_loc_time * 100 << "%\n"
           << "Total time for the lhs computation : " << lhs_loc_time << "\n"
           << "Total time of the cuda init : " << cuda_init << "\n"
           << "Total time of lhs with the init : " << cuda_init + lhs_loc_time << "\n"
           << "-------------------------------------------------------------------------------------\n\n";
    wbuild << lhs_loc_time << ",";
    lhs_time += lhs_loc_time;
  }
}

#endif

#ifdef ARCANE_HAS_CUDA

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* USE_BLCSR methods                                                         */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_buildMatrixBuildLessCsr()
{

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // Compute the number of nnz and initialize the memory space
  Int32 nnz = nbFace() * 2 + nbNode();
  m_csr_matrix.initialize(m_dof_family, nnz, nbNode());

  /*removing the neoighbouring currently as it is useless
  // Creating a connectivity from node to their neighbouring nodes
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  cn = idx_cn->connectivity();
  */
  ENUMERATE_NODE (inode, allNodes()) {

    //Since we compute the neighbouring connectivity here, we also fill the csr matrix

    Node node = *inode;

    m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

    for (Face face : node.faces()) {
      if (face.nodeId(0) == node.localId()) {
        //    cn->addConnectedItem(node, face.node(0));
        m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
      }
      else {
        //  cn->addConnectedItem(node, face.node(1));
        m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_assembleBLCsrBilinearOperatorTria3()
{
  Timer::Action timer_blcsr_bili(this->subDomain(), "AssembleBlcsrBilinearOperatorTria3");

  std::chrono::_V2::system_clock::time_point lhs_start;
  double global_build_average = 0;
  double build_time = 0;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using GPU BLCSR with NumArray format\n";
    lhs_start = std::chrono::high_resolution_clock::now();
  }

  {
    Timer::Action timer_blcsr_build(this->subDomain(), "BlcsrBuildMatrix");
    // Build the csr matrix
    _buildMatrixBuildLessCsr();
  }

  std::chrono::_V2::system_clock::time_point var_init_start;
  if (m_register_time) {
    auto build_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> build_duration = build_stop - lhs_start;
    build_time = build_duration.count();
    var_init_start = std::chrono::high_resolution_clock::now();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();

  // Boucle sur les noeuds déportée sur accélérateur
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
  Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
  auto in_out_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);

  auto in_node_coord = ax::viewIn(command, m_node_coord);

  UnstructuredMeshConnectivityView m_connectivity_view;
  m_connectivity_view.setMesh(this->mesh());
  auto ncc = m_connectivity_view.nodeCell();
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
  Arcane::ItemGenericInfoListView cells_infos(this->mesh()->cellFamily());

  std::chrono::_V2::system_clock::time_point loop_start;
  double var_init_time = 0;
  if (m_register_time) {
    auto var_init_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> var_init_duration = var_init_stop - var_init_start;
    var_init_time = var_init_duration.count();
    loop_start = std::chrono::high_resolution_clock::now();
  }

  Timer::Action timer_blcsr_add_compute(this->subDomain(), "BlcsrAddAndCompute");
  command << RUNCOMMAND_ENUMERATE(Node, inode, allNodes())
  {
    Int32 inode_index = 0;
    for (auto cell : ncc.cells(inode)) {

      // How can I know the right index ?
      // By checking in the global id ?
      // Working currently, but maybe only because p = 1 ?
      if (inode == cnc.nodeId(cell, 1)) {
        inode_index = 1;
      }
      else if (inode == cnc.nodeId(cell, 2)) {
        inode_index = 2;
      }
      else {
        inode_index = 0;
      }
      Real3 m0 = in_node_coord[cnc.nodeId(cell, 0)];
      Real3 m1 = in_node_coord[cnc.nodeId(cell, 1)];
      Real3 m2 = in_node_coord[cnc.nodeId(cell, 2)];

      Real area = 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y)); // calculate area

      Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
      Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
      Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

      Real b_matrix[3][2] = { 0 };
      Real mul = (1.0 / (2.0 * area));
      b_matrix[0][0] = dPhi0.x * mul;
      b_matrix[0][1] = dPhi0.y * mul;

      b_matrix[1][0] = dPhi1.x * mul;
      b_matrix[1][1] = dPhi1.y * mul;

      b_matrix[2][0] = dPhi2.x * mul;
      b_matrix[2][1] = dPhi2.y * mul;

      Int32 i = 0;
      for (NodeLocalId node2 : cnc.nodes(cell)) {
        Real x = 0.0;
        for (Int32 k = 0; k < 2; k++) {
          x += b_matrix[inode_index][k] * b_matrix[i][k];
        }
        if (nodes_infos.isOwn(inode)) {

          Int32 row = node_dof.dofId(inode, 0).localId();
          Int32 col = node_dof.dofId(node2, 0).localId();
          Int32 begin = in_row_csr[row];
          Int32 end;
          if (row == row_csr_size - 1) {
            end = col_csr_size;
          }
          else {
            end = in_row_csr[row + 1];
          }
          while (begin < end) {
            if (in_col_csr[begin] == col) {
              in_out_val_csr[begin] += x * area;
              break;
            }
            begin++;
          }
        }
        i++;
      }
    }
  };

  if (m_register_time) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - lhs_start;
    std::chrono::duration<double> loop_duration = lhs_end - loop_start;

    double loop_time = loop_duration.count();
    double lhs_loc_time = duration.count();
    logger << "Building time of the coo matrix :" << build_time << "\n"
           << "Variable initialisation time : " << var_init_time << "\n"
           << "Computation and Addition time : " << loop_time << "\n"
           << "LHS Total time : " << lhs_loc_time << "\n"
           << "Build matrix time in lhs :" << build_time / lhs_loc_time * 100 << "%\n"
           << "Variable initialisation time in lhs : " << var_init_time / lhs_loc_time * 100 << "%\n"
           << "Computation and Addition time in lhs : " << loop_time / lhs_loc_time * 100 << "%\n\n"
           << "-------------------------------------------------------------------------------------\n\n";
    lhs_time += lhs_loc_time;
    wbuild << lhs_loc_time << "\n";
    timer << loop_time << "\n";
  }
}

#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* USE_LEGACY methods                                                        */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTRIA3()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Timer::Action timer_action(this->subDomain(), "AssembleLegacyBilinearOperatorTria3");

  std::chrono::_V2::system_clock::time_point lhs_start;
  double compute_average = 0;
  double global_build_average = 0;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using hashmap legacy format\n";
    lhs_start = std::chrono::high_resolution_clock::now();
  }

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    std::chrono::_V2::system_clock::time_point compute_El_start;
    if (m_register_time) {
      compute_El_start = std::chrono::high_resolution_clock::now();
    }

    FixedMatrix<3, 3> K_e;
    {
      Timer::Action timer_action(this->subDomain(), "LegacyComputeElementMatrix");
      K_e = _computeElementMatrixTRIA3(cell); // element stifness matrix
    }

    if (m_register_time) {
      auto compute_El_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> compute_duration = compute_El_stop - compute_El_start;
      compute_average += compute_duration.count();
    }

    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]

    std::chrono::_V2::system_clock::time_point global_build_start;
    if (m_register_time) {
      global_build_start = std::chrono::high_resolution_clock::now();
    }

    {
      Timer::Action timer_action(this->subDomain(), "LegacyAddToGlobalMatrix");
      Int32 n1_index = 0;
      for (Node node1 : cell.nodes()) {
        Int32 n2_index = 0;
        for (Node node2 : cell.nodes()) {
          // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
          Real v = K_e(n1_index, n2_index);
          // m_k_matrix(node1.localId(), node2.localId()) += v;
          if (node1.isOwn()) {
            m_linear_system.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
          }
          ++n2_index;
        }
        ++n1_index;
      }
    }
    if (m_register_time) {
      auto global_build_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> global_build_duration = global_build_stop - global_build_start;
      global_build_average += global_build_duration.count();
    }
  }

  if (m_register_time) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - lhs_start;
    double lhs_loc_time = duration.count();
    logger << "compute elements average time : " << compute_average / nbCell() << "\n"
           << "compute elements total time : " << compute_average << "\n"
           << "add in global matrix average time : " << global_build_average / nbCell() << "\n"
           << "add in global matrix total time : " << global_build_average << "\n"
           << "lhs total time : " << lhs_loc_time << "\n"
           << "compute element time in lhs : " << compute_average / lhs_loc_time * 100 << "%\n"
           << "add in global matrix time in lhs : " << global_build_average / lhs_loc_time * 100 << "%\n\n"
           << "-------------------------------------------------------------------------------------\n\n";
    lhs_time += lhs_loc_time;
    wbuild << lhs_loc_time << ",";
    timer << lhs_loc_time << ",";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{

  Timer::Action timer_action(this->subDomain(), "Solving");

  std::chrono::_V2::system_clock::time_point solve_start;
  if (m_register_time) {
    solve_start = std::chrono::high_resolution_clock::now();
  }

  m_linear_system.solve();

  // Re-Apply boundary conditions because the solver has modified the value
  // of u on all nodes
  _applyDirichletBoundaryConditions();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    // Copy RHS DoF to Node u
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_u[node_dof.dofId(node, 0)];
      m_u[node] = v;
    }
  }

  //test
  m_u.synchronize();
  // def update_T(self,T):
  //     """Update u value on nodes after the FE resolution"""
  //     for i in range(0,len(self.mesh.nodes)):
  //         node=self.mesh.nodes[i]
  //         # don't update T imposed by Dirichlet BC
  //         if not node.is_T_fixed:
  //             self.mesh.nodes[i].T=T[i]

  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "T[" << node.localId() << "][" << node.uniqueId() << "] = "
             << m_u[node];
      //info() << "T[]" << node.uniqueId() << " "
      //       << m_u[node];
    }
  }

  if (m_register_time) {
    auto solve_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> solve_duration = solve_end - solve_start;
    solver_time = solve_duration.count();
    logger << "Solver duration : " << solver_time << "\n";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_checkResultFile()
{
  String filename = options()->resultFile();
  info() << "CheckResultFile filename=" << filename;
  if (filename.empty())
    return;
  const double epsilon = 1.0e-4;
  checkNodeResultFile(traceMng(), filename, m_u, epsilon);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
