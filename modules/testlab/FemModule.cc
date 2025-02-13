// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Testlab is testing ground for algorithms in ArcaneFEM.                    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_dumpTimeStats()
{
  Int64 nb_node = mesh()->ownNodes().size();
  Int64 total_nb_node = mesh()->parallelMng()->reduce(Parallel::ReduceSum, nb_node);

  Int64 nb_face = mesh()->outerFaces().size();
  Int64 total_nb_boundary_elt = mesh()->parallelMng()->reduce(Parallel::ReduceSum, nb_face);

  Int64 nb_cell = mesh()->ownCells().size();
  Int64 total_nb_elt = mesh()->parallelMng()->reduce(Parallel::ReduceSum, nb_cell);

  // Only master sub domain values are representative for time statistics
  if (!_isMasterRank())
    return;

  ofstream dump_file("./output/listing/time_stats.json");
  JSONWriter json_writer(JSONWriter::FormatFlags::None);

  json_writer.beginObject();

  json_writer.write("cacheWarming", m_cache_warming);
  json_writer.write("nbParallelInstance", parallelMng()->commSize());

  ParameterList parameter_list = this->subDomain()->application()->applicationInfo().commandLineArguments().parameters();
  if (m_running_on_gpu)
    json_writer.write("acceleratorRuntime", parameter_list.getParameterOrNull("AcceleratorRuntime"));

  json_writer.write("meshDim", defaultMesh()->dimension());
  json_writer.write("nbNode", total_nb_node);
  json_writer.write("nbBoundaryElement", total_nb_boundary_elt);
  json_writer.write("nbElement", total_nb_elt);

  m_time_stats->dumpStatsJSON(json_writer);

  json_writer.endObject();

  dump_file << json_writer.getBuffer();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
endModule()
{
  _dumpTimeStats();
}

void FemModule::
compute()
{
  info() << "[ArcaneFem-Info] Started module compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");
  if (m_petsc_flags != NULL)
  _setPetscFlagsFromCommandline();
  _printArcaneFemTime("[ArcaneFem-Timer] init-linear-system", (platform::getRealTime() - elapsedTime));

  Int64 nb_node = mesh()->ownNodes().size();
  Int64 total_nb_node = mesh()->parallelMng()->reduce(Parallel::ReduceSum, nb_node);

  Int64 nb_face = mesh()->outerFaces().size();
  Int64 total_nb_boundary_elt = mesh()->parallelMng()->reduce(Parallel::ReduceSum, nb_face);

  Int64 nb_cell = mesh()->ownCells().size();
  Int64 total_nb_elt = mesh()->parallelMng()->reduce(Parallel::ReduceSum, nb_cell);

  info() << "[ArcaneFem-Info] mesh dimension " << defaultMesh()->dimension();
  info() << "[ArcaneFem-Info] mesh boundary elements " << total_nb_boundary_elt;
  info() << "[ArcaneFem-Info] mesh cells " << total_nb_elt;
  info() << "[ArcaneFem-Info] mesh nodes " << total_nb_node;

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "[ArcaneFem-Info] Started module startInit()";
  Real elapsedTime = platform::getRealTime();
  Real TimeStart;

  m_queue = *(acceleratorMng()->defaultQueue());
  // When everything will be available on the GPU we will be able to
  // use device memory.
  //if (m_queue.isAcceleratorPolicy())
  //m_queue.setMemoryRessource(eMemoryRessource::Device);

  {
    IMesh* mesh = defaultMesh();
    if (mesh->dimension() == 3) {
      TimeStart = platform::getRealTime();
      m_node_node_via_edge_connectivity = MeshUtils::computeNodeNodeViaEdgeConnectivity(defaultMesh(), "NodeNodeViaEdge");
      m_node_node_via_edge_connectivity->connectivity()->dumpStats(std::cout);
      IndexedNodeNodeConnectivityView nn_cv = m_node_node_via_edge_connectivity->view();
      Int64 nb_edge = 0;
      ENUMERATE_NODE (inode, allNodes()) {
        Node node = *inode;
        nb_edge += nn_cv.nbNode(node);
      }
      m_nb_edge = nb_edge / 2;
      _printArcaneFemTime("[ArcaneFem-Timer] init-nde-nde-contvty", (platform::getRealTime() - TimeStart));
    }
    else {
      m_nb_edge = mesh->nbEdge();
    }

    if (options()->bsr || options()->bsrAtomicFree()) {
      bool use_csr_in_linear_system = options()->linearSystem.serviceName() == "HypreLinearSystem";
      m_bsr_format.initialize(mesh, use_csr_in_linear_system, options()->bsrAtomicFree);
    }
  }

  TimeStart = platform::getRealTime();
  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();
  _printArcaneFemTime("[ArcaneFem-Timer] initialize-DOFs", (platform::getRealTime() - TimeStart));

  _handleFlags();

  TimeStart = platform::getRealTime();
  _initBoundaryconditions();
  _printArcaneFemTime("[ArcaneFem-Timer] initialize-bc", (platform::getRealTime() - TimeStart));

  _checkCellType();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] initialize-total", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule ::
_handleFlags()
{
  ParameterList parameter_list = this->subDomain()->application()->applicationInfo().commandLineArguments().parameters();
  String cache_warm = parameter_list.getParameterOrNull("cache_warming");
  if (cache_warm != NULL) {
    auto tmp = Convert::Type<Integer>::tryParse(cache_warm);
    m_cache_warming = *tmp;
    info() << "[ArcaneFem-Info] cache warming activated via CLI";
    info() << "[ArcaneFem-Info] A cache warming of " << m_cache_warming << " iterations will happen";
  }
  if (cache_warm == NULL) {
    m_cache_warming = options()->cacheWarming();
    if (m_cache_warming != 1)
      info() << "[ArcaneFem-Info] A cache warming of " << m_cache_warming << " iterations will happen";
  }
  if (parameter_list.getParameterOrNull("COO") == "TRUE" || options()->coo()) {
    m_use_coo = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format COO activated via CLI";
  }
  if (parameter_list.getParameterOrNull("S-COO") == "TRUE" || options()->cooSorting()) {
    m_use_coo_sort = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format S_COO activated via CLI";
  }
  if (parameter_list.getParameterOrNull("COO_GPU") == "TRUE" || options()->cooGpu()) {
    m_use_coo_gpu = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format COO_GPU activated via CLI";
  }
  if (parameter_list.getParameterOrNull("S-COO_GPU") == "TRUE" || options()->cooSortingGpu()) {
    m_use_coo_sort_gpu = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format S_COO_GPU activated via CLI";
  }
  if (parameter_list.getParameterOrNull("CSR") == "TRUE" || options()->csr()) {
    m_use_csr = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format CSR activated via CLI";
  }
  if (parameter_list.getParameterOrNull("CSR_GPU") == "TRUE" || options()->csrGpu()) {
    m_use_csr_gpu = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format CSR_GPU activated via CLI";
  }
  if (parameter_list.getParameterOrNull("AF-CSR_GPU") == "TRUE" || options()->nwcsr()) {
    m_use_nodewise_csr = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format AF_CSR_GPU activated via CLI";
  }
  if (parameter_list.getParameterOrNull("BL-CSR_GPU") == "TRUE" || options()->blcsr()) {
    m_use_buildless_csr = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format BL_CSR_GPU activated via CLI";
  }
  if (parameter_list.getParameterOrNull("BSR_GPU") == "TRUE" || options()->bsr) {
    m_use_bsr = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format BSR_GPU activated via CLI";
  }
  if (parameter_list.getParameterOrNull("AF-BSR_GPU") == "TRUE" || options()->bsrAtomicFree()) {
    m_use_bsr_atomic_free = true;
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format atomic free BSR_GPU activated via CLI";
  }
  if (parameter_list.getParameterOrNull("DOK") == "TRUE" || m_use_legacy || options()->legacy()) {
    m_use_legacy = true;
    info() << "[ArcaneFem-Info] Matrix format DOK LEGACYactivated via CLI";
  }
  else if (parameter_list.getParameterOrNull("DOK") == "FALSE" || options()->legacy()) {
    m_use_legacy = false;
    info() << "[ArcaneFem-Info] Matrix format DOK disabled via CLI";
  }
  if (parameter_list.getParameterOrNull("AcceleratorRuntime") == "cuda") {
    m_running_on_gpu = true;
    info() << "[ArcaneFem-Info] CUDA Accelerator Runtime for GPU";
  }
  if (parameter_list.getParameterOrNull("AcceleratorRuntime") == "hip") {
    m_running_on_gpu = true;
    info() << "[ArcaneFem-Info] HIP Accelerator Runtime for GPU";
  }
  if (parameter_list.getParameterOrNull("solve_linear_system") == "FALSE") {
    m_solve_linear_system = false;
    info() << "[ArcaneFem-Info] Linear system assembled not solved (solve_linear_system = FALSE)";
  }
  if (parameter_list.getParameterOrNull("cross_validation") == "FALSE") {
    m_cross_validation = false;
    info() << "[ArcaneFem-Info] Cross validation disabled (cross_validation = FALSE)";
  }
  m_petsc_flags = parameter_list.getParameterOrNull("petsc_flags");
  if (m_petsc_flags != NULL) {
    info() << "[ArcaneFem-Info] PETSc flags the user provided will be used (petsc_flags != NULL)";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<3, 3> computeElementMatrixTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dxU = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyU = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  return area * (dxU ^ dxU) + area * (dyU ^ dyU);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<4, 4> computeElementMatrixTetra4(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dxU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzU = Arcane::FemUtils::Gpu::FeOperation3D::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);
  return volume * (dxU ^ dxU) + volume * (dyU ^ dyU) + volume * (dzU ^ dzU);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  Real assemblyTimeStart; // Timer variable

  Timer::Action timer_action(m_time_stats, "StationarySolve");
  Accelerator::ProfileRegion ps1(m_queue, "StationarySolve", 0xFF00FF);

  _getMaterialParameters();

  auto dim = mesh()->dimension();

  // Assemble the FEM bilinear operator (LHS - matrix A)
  if (m_use_legacy) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleBilinearOperatorTRIA3 : &FemModule::_assembleBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-DOK-mat", (platform::getRealTime() - assemblyTimeStart));
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Legacy");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-DOK-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_bsr) {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(m_queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    assemblyTimeStart = platform::getRealTime();

    m_bsr_format.computeSparsityAtomic();
    if (dim == 2)
      m_bsr_format.assembleBilinear([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3(cell_lid, cn_cv, in_node_coord); });
    else
      m_bsr_format.assembleBilinear([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4(cell_lid, cn_cv, in_node_coord); });
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-BSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));

    for (auto i = 1; i < m_cache_warming; ++i) {
      assemblyTimeStart = platform::getRealTime();
      m_bsr_format.resetMatrixValues();
      m_bsr_format.computeSparsityAtomic();
      if (dim == 2)
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3(cell_lid, cn_cv, in_node_coord); });
      else
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4(cell_lid, cn_cv, in_node_coord); });
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-BSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_bsr_atomic_free) {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(m_queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    assemblyTimeStart = platform::getRealTime();

    m_bsr_format.computeSparsityAtomicFree();
    if (dim == 2)
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3(cell_lid, cn_cv, in_node_coord); });
    else
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4(cell_lid, cn_cv, in_node_coord); });
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-AF_BSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));

    for (auto i = 1; i < m_cache_warming; ++i) {
      assemblyTimeStart = platform::getRealTime();
      m_bsr_format.resetMatrixValues();
      m_bsr_format.computeSparsityAtomicFree();
      if (dim == 2)
        m_bsr_format.assembleBilinear([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3(cell_lid, cn_cv, in_node_coord); });
      else
        m_bsr_format.assembleBilinear([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4(cell_lid, cn_cv, in_node_coord); });
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-AF_BSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_csr) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCsrBilinearOperatorTRIA3 : &FemModule::_assembleCsrBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-CSR-mat", (platform::getRealTime() - assemblyTimeStart));
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Csr");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-CSR-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_coo) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooBilinearOperatorTRIA3 : &FemModule::_assembleCooBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-COO-mat", (platform::getRealTime() - assemblyTimeStart));
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Coo");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-COO-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_coo_sort) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooSortBilinearOperatorTRIA3 : &FemModule::_assembleCooSortBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-COO-mat", (platform::getRealTime() - assemblyTimeStart));
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_CooSort");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-S_COO-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_coo_gpu) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooGPUBilinearOperatorTRIA3 : &FemModule::_assembleCooGPUBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-COO_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Coo_Gpu");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-COO_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_coo_sort_gpu) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooSortGPUBilinearOperatorTRIA3 : &FemModule::_assembleCooSortGPUBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-S_COO_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_CooSort_Gpu");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-S_COO_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_csr_gpu) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCsrGPUBilinearOperatorTRIA3 : &FemModule::_assembleCsrGPUBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-CSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Csr_Gpu");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-CSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_nodewise_csr) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleNodeWiseCsrBilinearOperatorTria3 : &FemModule::_assembleNodeWiseCsrBilinearOperatorTetra4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-AF_CSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_CsrNodeWise");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-AF_CSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  if (m_use_buildless_csr) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleBuildLessCsrBilinearOperatorTria3 : &FemModule::_assembleBuildLessCsrBilinearOperatorTetra4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    _printArcaneFemTime("[ArcaneFem-Timer] assemble-BL_CSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_CsrBuildLess");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      _printArcaneFemTime("[ArcaneFem-Timer] assemble-BL_CSR_GPU-mat", (platform::getRealTime() - assemblyTimeStart));
    }
  }

  // Assemble the FEM linear operator (RHS - vector b)
  if (m_use_buildless_csr || m_use_csr_gpu || m_use_nodewise_csr || m_use_csr) {
    //_assembleCsrLinearOperator();
    _assembleCsrGpuLinearOperator();
    {
      Timer::Action timer_action(m_time_stats, "TranslateToLinearSystem");
      m_csr_matrix.translateToLinearSystem(m_linear_system, m_queue);
    }
    _translateRhs();
  }
  else if (m_use_bsr || m_use_bsr_atomic_free) {
    _assembleLinearOperator(&(m_bsr_format.matrix()));
    Timer::Action timer_action(m_time_stats, "TranslateToLinearSystem");
    m_bsr_format.toLinearSystem(m_linear_system);
  }
  else {
    if (m_use_coo || m_use_coo_sort || m_use_coo_gpu || m_use_coo_sort_gpu) {
      Timer::Action timer_action(m_time_stats, "TranslateToLinearSystem");
      m_coo_matrix.translateToLinearSystem(m_linear_system);
    }
    _assembleLinearOperator();
  }

  // solve linear system
  if (m_solve_linear_system)
    _solve();

  // Check results
  if (m_solve_linear_system && m_cross_validation)
    _checkResultFile();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  f = options()->f();
  ElementNodes = 3.;

  if (options()->meshType == "TETRA4")
    ElementNodes = 4.;

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_initBoundaryconditions()
{
  info() << "[ArcaneFem-Info] Started module _initBoundaryconditions()";

  _applyDirichletBoundaryConditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_applyDirichletBoundaryConditionsGpu()
{
  info() << "[ArcaneFem-Info] Started module _applyDirichletBoundaryConditionsGpu()";

  for (const auto& bs : options()->dirichletBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real value = bs->value();
    info() << "Apply Dirichlet boundary condition surface=" << group.name() << " v=" << value
           << " via " << options()->enforceDirichletMethod() << " method ";

    RunQueue* queue = acceleratorMng()->defaultQueue();
    auto command = makeCommand(queue);

    UnstructuredMeshConnectivityView m_connectivity_view;
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    m_connectivity_view.setMesh(this->mesh());
    auto fnc = m_connectivity_view.faceNode();
    auto out_u_dirichlet = ax::viewOut(command, m_u_dirichlet);
    auto out_u = ax::viewOut(command, m_u);

    command << RUNCOMMAND_ENUMERATE(Face, iface, group)
    {
      for (NodeLocalId node : fnc.nodes(iface)) {
        out_u[node] = value;
        out_u_dirichlet[node] = true;
      }
    };
  }

  for (const auto& bs : options()->dirichletPointCondition()) {

    RunQueue* queue = acceleratorMng()->defaultQueue();
    auto command = makeCommand(queue);

    auto out_u = ax::viewOut(command, m_u);
    auto out_u_dirichlet = ax::viewOut(command, m_u_dirichlet);

    NodeGroup group = bs->node();
    Real value = bs->value();
    info() << "Apply Dirichlet point condition node=" << group.name() << " v=" << value
           << " via " << options()->enforceDirichletMethod() << " method ";
    command << RUNCOMMAND_ENUMERATE(Node, inode, group)
    {
      out_u[inode] = value;
      out_u_dirichlet[inode] = true;
    };
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_applyDirichletBoundaryConditions()
{
  info() << "[ArcaneFem-Info] Started module _applyDirichletBoundaryConditions()";

  for (const auto& bs : options()->dirichletBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real value = bs->value();
    info() << "Apply Dirichlet boundary condition surface=" << group.name() << " v=" << value
           << " via " << options()->enforceDirichletMethod() << " method ";
    ENUMERATE_ (Face, iface, group) {
      for (Node node : iface->nodes()) {
        //Original Code
        m_u[node] = value;
        m_u_dirichlet[node] = true;
      }
    }
  }

  for (const auto& bs : options()->dirichletPointCondition()) {
    NodeGroup group = bs->node();
    Real value = bs->value();
    info() << "Apply Dirichlet point condition node=" << group.name() << " v=" << value
           << " via " << options()->enforceDirichletMethod() << " method ";
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
  info() << "[ArcaneFem-Info] Started module _checkCellType()";
  Real elapsedTime = platform::getRealTime();

  Int16 type = 0;
  if (options()->meshType == "TETRA4") {
    type = IT_Tetraedron4;
  }
  else {
    type = IT_Triangle3;
  }
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    if (cell.type() != type)
      ARCANE_FATAL("Only Triangle3 cell type is supported");
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] mesh-check", elapsedTime);
}

/*---------------------------------------------------------------------------*/
// Assemble the FEM linear operator
//  - This function enforces a Dirichlet boundary condition in a weak sense
//    via the penalty method
//  - The method also adds source term
//  - TODO: external fluxes
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator(BSRMatrix<1>* bsr_matrix)
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  // time registration
  Timer::Action timer_action(m_time_stats, "AssembleLinearOperator");

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (options()->enforceDirichletMethod() == "Penalty") {

    Timer::Action timer_action(m_time_stats, "Penalty");

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

    Real Penalty = options()->penalty(); // 1.0e30 is the default

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u_dirichlet[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        // This SetValue should be updated in the acoording format we have (such as COO or CSR)
        if (bsr_matrix)
          bsr_matrix->setValue(dof_id, dof_id, Penalty);
        else
          m_linear_system.matrixSetValue(dof_id, dof_id, Penalty);
        Real u_g = Penalty * m_u[node_id];
        // This should be changed for a numArray
        rhs_values[dof_id] = u_g;
      }
    }
  }
  else if (options()->enforceDirichletMethod() == "WeakPenalty") {
    Timer::Action timer_action(m_time_stats, "WeakPenalty");

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

    Real Penalty = options()->penalty(); // 1.0e30 is the default

    // The same as before
    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u_dirichlet[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_linear_system.matrixAddValue(dof_id, dof_id, Penalty);
        Real u_g = Penalty * m_u[node_id];
        rhs_values[dof_id] = u_g;
      }
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
    ARCANE_FATAL("RowElimination Not Implemented for CSR GPU.");
    // The same as before
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
    ARCANE_FATAL("RowColumnElimination Not Implemented for CSR GPU.");
    // The same as before
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

  {
    Timer::Action timer_action(m_time_stats, "ConstantSourceTermAssembly");
    //----------------------------------------------
    // Constant source term assembly
    //----------------------------------------------
    //
    //  $int_{Omega}(f*v^h)$
    //  only for noded that are non-Dirichlet
    //----------------------------------------------
    if (options()->meshType == "TRIA3") {
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real area = _computeAreaTriangle3(cell);
        for (Node node : cell.nodes()) {
          if (!(m_u_dirichlet[node]) && node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += f * area / ElementNodes;
          }
        }
      }
    }

    if (options()->meshType == "TETRA4") {
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real area = _computeAreaTetra4(cell);
        for (Node node : cell.nodes()) {
          if (!(m_u_dirichlet[node]) && node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += f * area / ElementNodes;
          }
        }
      }
    }
  }
  {
    Timer::Action timer_action(m_time_stats, "ConstantFluxTermAssembly");

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
              // must replace rhs_values with numArray
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
              // must replace rhs_values with numArray
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
              // must replace rhs_values with numArray
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
              // must replace rhs_values with numArray
              rhs_values[node_dof.dofId(node, 0)] += (Normal.y * valueY) * length / 2.;
          }
        }
        continue;
      }
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] rhs-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleCsrLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  // time registration
  Timer::Action timer_action(m_time_stats, "CsrAssembleLinearOperator");

  m_rhs_vect.resize(nbNode());
  m_rhs_vect.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (options()->enforceDirichletMethod() == "Penalty") {

    Timer::Action timer_action(m_time_stats, "CsrPenalty");

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

    Real Penalty = options()->penalty(); // 1.0e30 is the default

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u_dirichlet[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_csr_matrix.matrixSetValue(dof_id, dof_id, Penalty);
        Real u_g = Penalty * m_u[node_id];
        m_rhs_vect[dof_id] = u_g;
      }
    }
  }
  else if (options()->enforceDirichletMethod() == "WeakPenalty") {
    Timer::Action timer_action(m_time_stats, "CsrWeakPenalty");

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

    Real Penalty = options()->penalty(); // 1.0e30 is the default

    // The same as before
    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u_dirichlet[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_csr_matrix.matrixAddValue(dof_id, dof_id, Penalty);
        Real u_g = Penalty * m_u[node_id];
        m_rhs_vect[dof_id] = u_g;
      }
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
    ARCANE_FATAL("RowElimination Not Implemented for CSR GPU.");
    // The same as before
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
    ARCANE_FATAL("RowColumnElimination Not Implemented for CSR GPU.");
    // The same as before
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

  {
    Timer::Action timer_action(m_time_stats, "CsrConstantSourceTermAssembly");
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
        if (!(m_u_dirichlet[node]) && node.isOwn()) {
          // Original code
          m_rhs_vect[node_dof.dofId(node, 0)] += f * area / ElementNodes;
        }
      }
    }
  }
  {
    Timer::Action timer_action(m_time_stats, "CsrConstantFluxTermAssembly");

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
              m_rhs_vect[node_dof.dofId(node, 0)] += value * length / 2.;
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
              m_rhs_vect[node_dof.dofId(node, 0)] += (Normal.x * valueX + Normal.y * valueY) * length / 2.;
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
              m_rhs_vect[node_dof.dofId(node, 0)] += (Normal.x * valueX) * length / 2.;
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
              m_rhs_vect[node_dof.dofId(node, 0)] += (Normal.y * valueY) * length / 2.;
          }
        }
        continue;
      }
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] rhs-assembly-csr", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE
Int32 FemModule::
_getValIndexCsrGpu(Int32 begin, Int32 end, DoFLocalId col, ax::NumArrayView<DataViewGetter<Int32>, MDDim1, DefaultLayout> csr_col)
{
  Int32 i = begin;
  while (i < end && col != csr_col(i)) {
    i++;
  }
  // The value has not been found
  if (i == end) {
    return -1;
  }
  // The value as been found
  return i;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrGpuLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleCsrGpuLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  Timer::Action timer_action(m_time_stats, "CsrGpuAssembleLinearOperator");

  m_rhs_vect.resize(nbNode());
  m_rhs_vect.fill(0.0);

  if (options()->enforceDirichletMethod() == "Penalty") {

    Timer::Action timer_action(m_time_stats, "CsrGpuPenalty");

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

    Real Penalty = options()->penalty(); // 1.0e30 is the default

    RunQueue* queue = acceleratorMng()->defaultQueue();
    auto command = makeCommand(queue);

    auto in_out_rhs_vect = ax::viewInOut(command, m_rhs_vect);
    auto in_csr_row = ax::viewIn(command, m_csr_matrix.m_matrix_row);
    auto in_csr_col = ax::viewIn(command, m_csr_matrix.m_matrix_column);
    auto in_out_csr_val = ax::viewInOut(command, m_csr_matrix.m_matrix_value);
    Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
    Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    auto in_m_u_dirichlet = ax::viewIn(command, m_u_dirichlet);
    auto in_m_u = ax::viewIn(command, m_u);
    // In this loop :
    // m_u_dirichlet must be adapted -> module variable, just need a view
    // m_u must be adapted
    // Set value must be replaced
    // m_rhs_vect must also be replaced
    command << RUNCOMMAND_ENUMERATE(Node, inode, ownNodes())
    {
      if (in_m_u_dirichlet(inode)) {
        DoFLocalId dof_id = node_dof.dofId(inode, 0);
        Int32 begin = in_csr_row(dof_id);
        Int32 end = 0;
        if (begin == row_csr_size - 1) {
          end = col_csr_size;
        }
        else {
          end = in_csr_row(dof_id + 1);
        }
        Int32 index = _getValIndexCsrGpu(begin, end, dof_id, in_csr_col);
        in_out_csr_val(index) = Penalty;
        Real u_g = Penalty * in_m_u(inode);
        in_out_rhs_vect(dof_id) = u_g;
      }
    };
  }
  else if (options()->enforceDirichletMethod() == "WeakPenalty") {
    Timer::Action timer_action(m_time_stats, "CsrGpuWeakPenalty");

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

    Real Penalty = options()->penalty(); // 1.0e30 is the default

    RunQueue* queue = acceleratorMng()->defaultQueue();
    auto command = makeCommand(queue);

    auto in_out_rhs_vect = ax::viewInOut(command, m_rhs_vect);
    auto in_csr_row = ax::viewIn(command, m_csr_matrix.m_matrix_row);
    auto in_csr_col = ax::viewIn(command, m_csr_matrix.m_matrix_column);
    auto in_out_csr_val = ax::viewInOut(command, m_csr_matrix.m_matrix_value);
    Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
    Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    auto in_m_u_dirichlet = ax::viewIn(command, m_u_dirichlet);
    auto in_m_u = ax::viewIn(command, m_u);
    // In this loop :
    // m_u_dirichlet must be adapted
    // m_u must have a view
    // Set value must be replaced
    // m_rhs_vect must also be replaced
    command << RUNCOMMAND_ENUMERATE(Node, inode, ownNodes())
    {
      if (in_m_u_dirichlet(inode)) {
        DoFLocalId dof_id = node_dof.dofId(inode, 0);
        Int32 begin = in_csr_row(dof_id);
        Int32 end = 0;
        if (begin == row_csr_size - 1) {
          end = col_csr_size;
        }
        else {
          end = in_csr_row(dof_id + 1);
        }
        Int32 index = _getValIndexCsrGpu(begin, end, dof_id, in_csr_col);
        ax::doAtomic<ax::eAtomicOperation::Add>(in_out_csr_val(index), Penalty);
        //in_out_csr_val(index) += Penalty;

        Real u_g = Penalty * in_m_u(inode);
        in_out_rhs_vect(dof_id) = u_g;
      }
    };
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
    ARCANE_FATAL("RowElimination Not Implemented for CSR GPU.");
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
    ARCANE_FATAL("RowColumnElimination Not Implemented for CSR GPU.");
  }
  else {

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " is not supported \n"
           << "enforce-Dirichlet-method only supports:\n"
           << "  - Penalty\n"
           << "  - WeakPenalty\n";
  }

  if (options()->meshType == "TRIA3") {
    Timer::Action timer_action(m_time_stats, "CsrGpuConstantSourceTermAssembly");
    //----------------------------------------------
    // Constant source term assembly
    //----------------------------------------------
    //
    //  $int_{Omega}(f*v^h)$
    //  only for noded that are non-Dirichlet
    //----------------------------------------------

    RunQueue* queue = acceleratorMng()->defaultQueue();
    auto command = makeCommand(queue);

    auto in_out_rhs_vect = ax::viewInOut(command, m_rhs_vect);

    auto in_m_u_dirichlet = ax::viewIn(command, m_u_dirichlet);

    Real tmp_f = f;
    Real tmp_ElementNodes = ElementNodes;

    UnstructuredMeshConnectivityView m_connectivity_view;
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    m_connectivity_view.setMesh(this->mesh());
    auto cnc = m_connectivity_view.cellNode();
    Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    // In this loop :
    // m_u_dirichlet must be adapted
    // node.isOwn must be adapted
    // m_rhs_vect must also be replaced
    // f and Element nodes must be put in local variable
    // computeArea must be replaced

    command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
    {
      Real area = _computeAreaTriangle3Gpu(icell, cnc, in_node_coord);
      for (NodeLocalId node : cnc.nodes(icell)) {
        if (!(in_m_u_dirichlet(node)) && nodes_infos.isOwn(node)) {
          // Original code
          Real val = tmp_f * area / tmp_ElementNodes;
          ax::doAtomic<ax::eAtomicOperation::Add>(in_out_rhs_vect(node_dof.dofId(node, 0)), val);
        }
      }
    };
  }

  if (options()->meshType == "TETRA4") {
    Timer::Action timer_action(m_time_stats, "CsrGpuConstantSourceTermAssembly");
    //----------------------------------------------
    // Constant source term assembly
    //----------------------------------------------
    //
    //  $int_{Omega}(f*v^h)$
    //  only for noded that are non-Dirichlet
    //----------------------------------------------

    RunQueue* queue = acceleratorMng()->defaultQueue();
    auto command = makeCommand(queue);

    auto in_out_rhs_vect = ax::viewInOut(command, m_rhs_vect);

    auto in_m_u_dirichlet = ax::viewIn(command, m_u_dirichlet);

    Real tmp_f = f;
    Real tmp_ElementNodes = ElementNodes;

    UnstructuredMeshConnectivityView m_connectivity_view;
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    m_connectivity_view.setMesh(this->mesh());
    auto cnc = m_connectivity_view.cellNode();
    Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    // In this loop :
    // m_u_dirichlet must be adapted
    // node.isOwn must be adapted
    // m_rhs_vect must also be replaced
    // f and Element nodes must be put in local variable
    // computeArea must be replaced

    command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
    {
      Real area = _computeAreaTetra4Gpu(icell, cnc, in_node_coord);
      for (NodeLocalId node : cnc.nodes(icell)) {
        if (!(in_m_u_dirichlet(node)) && nodes_infos.isOwn(node)) {
          // Original code
          Real val = tmp_f * area / tmp_ElementNodes;
          ax::doAtomic<ax::eAtomicOperation::Add>(in_out_rhs_vect(node_dof.dofId(node, 0)), val);
        }
      }
    };
  }

  {
    Timer::Action timer_action(m_time_stats, "CsrGpuConstantFluxTermAssembly");

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

        RunQueue* queue = acceleratorMng()->defaultQueue();
        auto command = makeCommand(queue);

        auto in_out_rhs_vect = ax::viewInOut(command, m_rhs_vect);

        auto in_m_u_dirichlet = ax::viewIn(command, m_u_dirichlet);

        UnstructuredMeshConnectivityView m_connectivity_view;
        auto in_node_coord = ax::viewIn(command, m_node_coord);
        m_connectivity_view.setMesh(this->mesh());
        auto fnc = m_connectivity_view.faceNode();
        Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
        auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

        // In this loop :
        // m_u_dirichlet must be adapted
        // node.isOwn must be adapted
        // m_rhs_vect must also be replaced
        // computeEdgeLength2 must be reimplemented
        command << RUNCOMMAND_ENUMERATE(Face, iface, group)
        {
          Real length = _computeEdgeLength2Gpu(iface, fnc, in_node_coord);
          for (NodeLocalId node : fnc.nodes(iface)) {
            if (!(in_m_u_dirichlet[node]) && nodes_infos.isOwn(node))
              ax::doAtomic<ax::eAtomicOperation::Add>(in_out_rhs_vect[node_dof.dofId(node, 0)], value * length / 2.);
            //in_out_rhs_vect[node_dof.dofId(node, 0)] += value * length / 2.;
          }
        };
        continue;
      }

      if (bs->valueX.isPresent() && bs->valueY.isPresent()) {
        Real valueX = bs->valueX();
        Real valueY = bs->valueY();

        RunQueue* queue = acceleratorMng()->defaultQueue();
        auto command = makeCommand(queue);

        auto in_out_rhs_vect = ax::viewInOut(command, m_rhs_vect);

        auto in_m_u_dirichlet = ax::viewIn(command, m_u_dirichlet);

        UnstructuredMeshConnectivityView m_connectivity_view;
        auto in_node_coord = ax::viewIn(command, m_node_coord);
        m_connectivity_view.setMesh(this->mesh());
        auto fnc = m_connectivity_view.faceNode();
        Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
        Arcane::FaceInfoListView faces_infos(this->mesh()->nodeFamily());
        auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

        // In this loop :
        // m_u_dirichlet must be adapted
        // node.isOwn must be adapted
        // m_rhs_vect must also be replaced
        // computeEdgeLength2 must be reimplemented
        // computeEdgeNormal2 must be reimplemented
        command << RUNCOMMAND_ENUMERATE(Face, iface, group)
        {
          Real length = _computeEdgeLength2Gpu(iface, fnc, in_node_coord);
          Real2 Normal = _computeEdgeNormal2Gpu(iface, fnc, in_node_coord, faces_infos);
          for (NodeLocalId node : fnc.nodes(iface)) {
            if (!(in_m_u_dirichlet[node]) && nodes_infos.isOwn(node)) {
              Real value = (Normal.x * valueX + Normal.y * valueY) * length / 2.;
              ax::doAtomic<ax::eAtomicOperation::Add>(in_out_rhs_vect[node_dof.dofId(node, 0)], value);
              //in_out_rhs_vect[node_dof.dofId(node, 0)] += value;
            }
          }
        };
        continue;
      }

      if (bs->valueX.isPresent()) {
        Real valueX = bs->valueX();

        RunQueue* queue = acceleratorMng()->defaultQueue();
        auto command = makeCommand(queue);

        auto in_out_rhs_vect = ax::viewInOut(command, m_rhs_vect);

        auto in_m_u_dirichlet = ax::viewIn(command, m_u_dirichlet);

        UnstructuredMeshConnectivityView m_connectivity_view;
        auto in_node_coord = ax::viewIn(command, m_node_coord);
        m_connectivity_view.setMesh(this->mesh());
        auto fnc = m_connectivity_view.faceNode();
        Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
        Arcane::FaceInfoListView faces_infos(this->mesh()->nodeFamily());
        auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

        // In this loop :
        // m_u_dirichlet must be adapted
        // node.isOwn must be adapted
        // m_rhs_vect must also be replaced
        // computeEdgeLength2 must be reimplemented
        // computeEdgeNormal2 must be reimplemented
        command << RUNCOMMAND_ENUMERATE(Face, iface, group)
        {
          Real length = _computeEdgeLength2Gpu(iface, fnc, in_node_coord);
          Real2 Normal = _computeEdgeNormal2Gpu(iface, fnc, in_node_coord, faces_infos);
          for (NodeLocalId node : fnc.nodes(iface)) {
            if (!(in_m_u_dirichlet[node]) && nodes_infos.isOwn(node)) {
              Real value = (Normal.x * valueX) * length / 2.;
              ax::doAtomic<ax::eAtomicOperation::Add>(in_out_rhs_vect[node_dof.dofId(node, 0)], value);
              //in_out_rhs_vect[node_dof.dofId(node, 0)] += value;
            }
          }
        };
        continue;
      }

      if (bs->valueY.isPresent()) {
        Real valueY = bs->valueY();

        RunQueue* queue = acceleratorMng()->defaultQueue();
        auto command = makeCommand(queue);

        auto in_out_rhs_vect = ax::viewInOut(command, m_rhs_vect);

        auto in_m_u_dirichlet = ax::viewIn(command, m_u_dirichlet);

        UnstructuredMeshConnectivityView m_connectivity_view;
        auto in_node_coord = ax::viewIn(command, m_node_coord);
        m_connectivity_view.setMesh(this->mesh());
        auto fnc = m_connectivity_view.faceNode();
        Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
        Arcane::FaceInfoListView faces_infos(this->mesh()->nodeFamily());
        auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

        // In this loop :
        // m_u_dirichlet must be adapted
        // node.isOwn must be adapted
        // m_rhs_vect must also be replaced
        // computeEdgeLength2 must be reimplemented
        // computeEdgeNormal2 must be reimplemented
        command << RUNCOMMAND_ENUMERATE(Face, iface, group)
        {
          Real length = _computeEdgeLength2Gpu(iface, fnc, in_node_coord);
          Real2 Normal = _computeEdgeNormal2Gpu(iface, fnc, in_node_coord, faces_infos);
          for (NodeLocalId node : fnc.nodes(iface)) {
            if (!(in_m_u_dirichlet[node]) && nodes_infos.isOwn(node)) {
              Real value = (Normal.y * valueY) * length / 2.;
              ax::doAtomic<ax::eAtomicOperation::Add>(in_out_rhs_vect[node_dof.dofId(node, 0)], value);
              //in_out_rhs_vect[node_dof.dofId(node, 0)] += (Normal.y * valueY) * length / 2.;
            }
          }
        };
        continue;
      }
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] rhs-assembly-csr-gpu", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void FemModule::
_translateRhs()
{
  info() << "[ArcaneFem-Info] Started module _translateRhs()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);
  for (Int32 i = 0; i < m_rhs_vect.dim1Size(); i++) {
    rhs_values[DoFLocalId(i)] = m_rhs_vect[DoFLocalId(i)];
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] rhs-translation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE
Real FemModule::
_computeAreaTetra4Gpu(CellLocalId icell, IndexedCellNodeConnectivityView cnc, ax::VariableNodeReal3InView in_node_coord)
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
  return std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2))) / 6.0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE
Real FemModule::
_computeAreaTriangle3Gpu(CellLocalId icell, IndexedCellNodeConnectivityView cnc, ax::VariableNodeReal3InView in_node_coord)
{
  Real3 m0 = in_node_coord[cnc.nodeId(icell, 0)];
  Real3 m1 = in_node_coord[cnc.nodeId(icell, 1)];
  Real3 m2 = in_node_coord[cnc.nodeId(icell, 2)];

  return 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y));
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
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE
Real FemModule::
_computeEdgeLength2Gpu(FaceLocalId iface, IndexedFaceNodeConnectivityView fnc, ax::VariableNodeReal3InView in_node_coord)
{
  Real3 m0 = in_node_coord[fnc.nodeId(iface, 0)];
  Real3 m1 = in_node_coord[fnc.nodeId(iface, 1)];
  return math::sqrt((m1.x - m0.x) * (m1.x - m0.x) + (m1.y - m0.y) * (m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeEdgeLength2(Face face)
{
  Real3 m0 = m_node_coord[face.nodeId(0)];
  Real3 m1 = m_node_coord[face.nodeId(1)];
  return math::sqrt((m1.x - m0.x) * (m1.x - m0.x) + (m1.y - m0.y) * (m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeAreaTetra4(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  Real3 m3 = m_node_coord[cell.nodeId(3)];

  // Calculate vectors representing edges of the tetrahedron
  Real3 v0 = m1 - m0;
  Real3 v1 = m2 - m0;
  Real3 v2 = m3 - m0;

  // Compute volume using scalar triple product
  return std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2))) / 6.0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE
Real2 FemModule::
_computeEdgeNormal2Gpu(FaceLocalId iface, IndexedFaceNodeConnectivityView fnc,
                       ax::VariableNodeReal3InView in_node_coord, Arcane::FaceInfoListView faces_infos)
{
  Real3 m0 = in_node_coord[fnc.nodeId(iface, 0)];
  Real3 m1 = in_node_coord[fnc.nodeId(iface, 1)];
  // We need to access this information on GPU
  if (!faces_infos.isSubDomainBoundaryOutside(iface)) {
    Real3 tmp = m0;
    m0 = m1;
    m1 = tmp;
  }
  Real2 N;
  Real norm_N = math::sqrt((m1.y - m0.y) * (m1.y - m0.y) + (m1.x - m0.x) * (m1.x - m0.x)); // for normalizing
  N.x = (m1.y - m0.y) / norm_N;
  N.y = (m0.x - m1.x) / norm_N;
  return N;
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
_computeElementMatrixTETRA4(Cell cell)
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
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  Real3 m3 = m_node_coord[cell.nodeId(3)];

  Real volume = _computeAreaTetra4(cell);

  // Compute gradients of shape functions
  Real3 dPhi0 = Arcane::math::cross(m2 - m1, m1 - m3);
  Real3 dPhi1 = Arcane::math::cross(m3 - m0, m0 - m2);
  Real3 dPhi2 = Arcane::math::cross(m1 - m0, m0 - m3);
  Real3 dPhi3 = Arcane::math::cross(m0 - m1, m1 - m2);

  // Construct the B-matrix
  FixedMatrix<3, 4> b_matrix;
  b_matrix(0, 0) = dPhi0.x;
  b_matrix(1, 0) = dPhi0.y;
  b_matrix(2, 0) = dPhi0.z;

  b_matrix(0, 1) = dPhi1.x;
  b_matrix(1, 1) = dPhi1.y;
  b_matrix(2, 1) = dPhi1.z;

  b_matrix(0, 2) = dPhi2.x;
  b_matrix(1, 2) = dPhi2.y;
  b_matrix(2, 2) = dPhi2.z;

  b_matrix(0, 3) = dPhi3.x;
  b_matrix(1, 3) = dPhi3.y;
  b_matrix(2, 3) = dPhi3.z;

  b_matrix.multInPlace(1.0 / (6.0 * volume));

  // Compute the element matrix
  FixedMatrix<4, 4> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(volume);

  /*
  cout << " Ae \n"
       << "\t" << int_cdPi_dPj(0,0)<<"\t"<< int_cdPi_dPj(0,1)<<"\t"<< int_cdPi_dPj(0,2)<<"\t"<< int_cdPi_dPj(0,3)<<"\n"
       << "\t" << int_cdPi_dPj(1,0)<<"\t"<< int_cdPi_dPj(1,1)<<"\t"<< int_cdPi_dPj(1,2)<<"\t"<< int_cdPi_dPj(1,3)<<"\n"
       << "\t" << int_cdPi_dPj(2,0)<<"\t"<< int_cdPi_dPj(2,1)<<"\t"<< int_cdPi_dPj(2,2)<<"\t"<< int_cdPi_dPj(2,3)<<"\n"
       << "\t" << int_cdPi_dPj(3,0)<<"\t"<< int_cdPi_dPj(3,1)<<"\t"<< int_cdPi_dPj(3,2)<<"\t"<< int_cdPi_dPj(3,3)<<"\n"
       << endl;
*/

  return int_cdPi_dPj;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  info() << "[ArcaneFem-Info] Started module _solve()";
  Real elapsedTime = platform::getRealTime();

  Real TimeStart;
  ITimeStats* tstat = m_time_stats;
  Timer::Action timer_action(tstat, "Solving");

  {
    TimeStart = platform::getRealTime();
    Timer::Action ta1(tstat, "LinearSystemSolve");
    m_linear_system.solve();
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] solving-linear-system", elapsedTime);

  // Re-Apply boundary conditions because the solver has modified the value
  {
    TimeStart = platform::getRealTime();
    Timer::Action ta1(tstat, "ApplyBoundaryConditions");
    _applyDirichletBoundaryConditions();
    _printArcaneFemTime("[ArcaneFem-Timer] reapply-bc", (platform::getRealTime() - TimeStart));
  }

  {
    TimeStart = platform::getRealTime();
    Timer::Action ta1(tstat, "CopySolution");
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    // Copy RHS DoF to Node u
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_u[node_dof.dofId(node, 0)];
      m_u[node] = v;
    }
    _printArcaneFemTime("[ArcaneFem-Timer] update-solution", (platform::getRealTime() - TimeStart));
  }

  TimeStart = platform::getRealTime();
  m_u.synchronize();
  _printArcaneFemTime("[ArcaneFem-Timer] synchronize", (platform::getRealTime() - TimeStart));

  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "u[" << node.localId() << "][" << node.uniqueId() << "] = "
             << m_u[node];
      //info() << "u[]" << node.uniqueId() << " "
      //       << m_u[node];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_build()
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_checkResultFile()
{
  info() << "[ArcaneFem-Info] Started module _checkResultFile()";
  Real elapsedTime = platform::getRealTime();

  String filename = options()->resultFile();
  if (filename.empty())
    return;
  const double epsilon = 1.0e-4;
  const double skipValuesMinLim = 1.0e-16;
  checkNodeResultFile(traceMng(), filename, m_u, epsilon, skipValuesMinLim);

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] cross-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

bool FemModule::
_isMasterRank() const
{
  return parallelMng()->isMasterIO();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Function to prints the execution time `value` of phase `label`
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_printArcaneFemTime(const String label, const Real value)
{
  info() << std::left << std::setw(45) << label << " = " << value;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Function to set PETSc flags from commandline
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_setPetscFlagsFromCommandline()
{
  StringList string_list;
  std::string petsc_flags_std = m_petsc_flags.localstr();
  // Use a string stream to split the string by spaces
  std::istringstream iss(petsc_flags_std);
  String token;
  while (iss >> token) {
    string_list.add(token);
  }

  CommandLineArguments args(string_list);
  m_linear_system.setSolverCommandLineArguments(args);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/