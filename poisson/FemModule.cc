// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2024 */
/*                                                                           */
/* Simple module to test simple FEM mechanism.                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_dumpTimeStats()
{
  Int64 nb_node = mesh()->ownNodes().size();
  Int64 total_nb_node = mesh()->parallelMng()->reduce(Parallel::ReduceSum, nb_node);

  Int64 nb_face = mesh()->ownFaces().size(); // Face in 3D, edge in 2D
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
  info() << "Module Fem COMPUTE";

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  Real TimeStart = platform::getRealTime();
  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());

  m_linear_system.initialize(subDomain(), acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");
  // Test for adding parameters for PETSc.
  // This is only used for the first call.
  {
    StringList string_list;
    /*
    string_list.add("-trmalloc");
    string_list.add("-log_trace");
    string_list.add("-ksp_monitor");
    string_list.add("-ksp_view");
    string_list.add("-math_view");
    string_list.add("draw");
    string_list.add("-draw_pause");
    string_list.add("-10");
    */
    CommandLineArguments args(string_list);
    m_linear_system.setSolverCommandLineArguments(args);
  }
  info() << "[ArcaneFem-Timer] Time to initialize linear-system = " << (platform::getRealTime() - TimeStart);

  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();

  _doStationarySolve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  Real TimeStart;
  info() << "Module Fem INIT";

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
      std::cout << "\n";
      IndexedNodeNodeConnectivityView nn_cv = m_node_node_via_edge_connectivity->view();
      Int64 nb_edge = 0;
      ENUMERATE_NODE (inode, allNodes()) {
        Node node = *inode;
        nb_edge += nn_cv.nbNode(node);
      }
      m_nb_edge = nb_edge / 2;
      info() << "Using custom node-node via edge connectivity: nb_edge=" << m_nb_edge;
      info() << "[ArcaneFem-Timer] Time to initialize node-node connectivity view = " << (platform::getRealTime() - TimeStart);
    }
    else {
      m_nb_edge = mesh->nbEdge();
      info() << "Number of edge: nb_edge=" << m_nb_edge;
    }

    if (options()->bsr) {
      bool use_csr_in_linear_system = options()->linearSystem.serviceName() == "HypreLinearSystem";
      m_bsr_format.initialize(mesh, mesh->dimension() == 2 ? nbFace() : m_nb_edge, use_csr_in_linear_system);
      m_bsr_format.computeSparsity(); // Need to be done just once.
    }
  }

  TimeStart = platform::getRealTime();
  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();
  info() << "[ArcaneFem-Timer] Time to initialize DOFs = " << (platform::getRealTime() - TimeStart);

  _handleFlags();

  TimeStart = platform::getRealTime();
  _initBoundaryconditions();
  info() << "[ArcaneFem-Timer] Time to initialize boundary conditions = " << (platform::getRealTime() - TimeStart);

  _checkCellType();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule ::
_handleFlags()
{
  ParameterList parameter_list = this->subDomain()->application()->applicationInfo().commandLineArguments().parameters();
  info() << "-----------------------------------------------------------------------------------------";
  String cache_warm = parameter_list.getParameterOrNull("CACHE_WARMING");
  if (cache_warm != NULL) {
    auto tmp = Convert::Type<Integer>::tryParse(cache_warm);
    m_cache_warming = *tmp;
    info() << "CACHE_WARMING: A cache warming of " << m_cache_warming << " iterations will happen";
  }
  if (cache_warm == NULL) {
    m_cache_warming = options()->cacheWarming();
    if (m_cache_warming != 1)
      info() << "CACHE_WARMING: A cache warming of " << m_cache_warming << " iterations will happen";
  }
  if (parameter_list.getParameterOrNull("COO") == "TRUE" || options()->coo()) {
    m_use_coo = true;
    m_use_legacy = false;
    info() << "COO: The COOrdinate data structure is used for sparse matrices";
  }
  if (parameter_list.getParameterOrNull("COO_SORT") == "TRUE" || options()->cooSorting()) {
    m_use_coo_sort = true;
    m_use_legacy = false;
    info() << "COO_SORT: The COOrdinate data structure with SORTing is used for sparse matrices";
  }
  if (parameter_list.getParameterOrNull("COO_GPU") == "TRUE" || options()->cooGpu()) {
    m_use_coo_gpu = true;
    m_use_legacy = false;
    info() << "COO_GPU: The GPU-compatible COOrdinate data structure is used for sparse matrices";
  }
  if (parameter_list.getParameterOrNull("COO_SORT_GPU") == "TRUE" || options()->cooSortingGpu()) {
    m_use_coo_sort_gpu = true;
    m_use_legacy = false;
    info() << "COO_SORT_GPU: The GPU-compatible COOrdinate data structure with SORTing is used for sparse matrices";
  }
  if (parameter_list.getParameterOrNull("CSR") == "TRUE" || options()->csr()) {
    m_use_csr = true;
    m_use_legacy = false;
    info() << "CSR: The Compressed Sparse Row data structure is used for sparse matrices";
  }
  if (parameter_list.getParameterOrNull("CSR_GPU") == "TRUE" || options()->csrGpu()) {
    m_use_csr_gpu = true;
    m_use_legacy = false;
    info() << "CSR_GPU: The GPU-compatible Compressed Sparse Row data structure is used for sparse matrices";
  }
  if (parameter_list.getParameterOrNull("NWCSR") == "TRUE" || options()->nwcsr()) {
    m_use_nodewise_csr = true;
    m_use_legacy = false;
    info() << "NWCSR: The GPU-compatible Compressed Sparse Row data structure is used for sparse matrices with Node-Wise computation";
  }
  if (parameter_list.getParameterOrNull("BLCSR") == "TRUE" || options()->blcsr()) {
    m_use_buildless_csr = true;
    m_use_legacy = false;
    info() << "BLCSR: The GPU-compatible Compressed Sparse Row (CSR) data structure is used for sparse matrices with Node-Wise computation in a Build Less manner";
  }
  if (parameter_list.getParameterOrNull("LEGACY") == "TRUE" || m_use_legacy || options()->legacy()) {
    m_use_legacy = true;
    info() << "DOK: The Dictionary Of Key ata structure is used for sparse matrices";
  }
  else if (parameter_list.getParameterOrNull("LEGACY") == "FALSE" || options()->legacy()) {
    m_use_legacy = false;
  }
  if (parameter_list.getParameterOrNull("AcceleratorRuntime") == "cuda") {
    m_running_on_gpu = true;
    info() << "CUDA: The methods able to use GPU will use it";
  }
  if (parameter_list.getParameterOrNull("AcceleratorRuntime") == "hip") {
    m_running_on_gpu = true;
    info() << "HIP: The methods able to use GPU will use it";
  }
  if (parameter_list.getParameterOrNull("SOLVE_LINEAR_SYSTEM") == "FALSE") {
    m_solve_linear_system = false;
    info() << "Linear system assembled but not solved (SOLVE_LINEAR_SYSTEM = FALSE)";
  }
  if (parameter_list.getParameterOrNull("CROSS_VALIDATION") == "FALSE") {
    m_cross_validation = false;
    info() << "Cross validation disabled (CROSS_VALIDATION = FALSE)";
  }
  if (options()->bsr) {
    m_use_bsr = true;
    m_use_legacy = false;
  }
  info() << "-----------------------------------------------------------------------------------------";
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

  assemblyTimeStart = platform::getRealTime();
  _getMaterialParameters();
  info() << "[ArcaneFem-Timer] Time to get material parameters = " << (platform::getRealTime() - assemblyTimeStart);

  auto dim = mesh()->dimension();

  if (m_use_bsr) {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(m_queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);

    if (dim == 2)
      m_bsr_format.assembleCellWise([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3(cell_lid, cn_cv, in_node_coord); });
    else
      m_bsr_format.assembleCellWise([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4(cell_lid, cn_cv, in_node_coord); });

    _assembleLinearOperator(&(m_bsr_format.matrix()));
    m_bsr_format.toLinearSystem(m_linear_system);
    _solve();
    _checkResultFile();
    return;
  }

  // Assemble the FEM bilinear operator (LHS - matrix A)
  if (m_use_legacy) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleBilinearOperatorTRIA3 : &FemModule::_assembleBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    info() << "[ArcaneFem-Timer] Time to assemble DOK matrix = " << (platform::getRealTime() - assemblyTimeStart);
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Legacy");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      info() << "[ArcaneFem-Timer] Time to assemble DOK matrix = " << (platform::getRealTime() - assemblyTimeStart);
    }
  }

  if (m_use_csr) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCsrBilinearOperatorTRIA3 : &FemModule::_assembleCsrBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    info() << "[ArcaneFem-Timer] Time to assemble CSR matrix = " << (platform::getRealTime() - assemblyTimeStart);
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Csr");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      info() << "[ArcaneFem-Timer] Time to assemble CSR matrix = " << (platform::getRealTime() - assemblyTimeStart);
    }
  }

  if (m_use_coo) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooBilinearOperatorTRIA3 : &FemModule::_assembleCooBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    info() << "[ArcaneFem-Timer] Time to assemble COO matrix = " << (platform::getRealTime() - assemblyTimeStart);
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Coo");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      info() << "[ArcaneFem-Timer] Time to assemble COO matrix = " << (platform::getRealTime() - assemblyTimeStart);
    }
  }

  if (m_use_coo_sort) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooSortBilinearOperatorTRIA3 : &FemModule::_assembleCooSortBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    info() << "[ArcaneFem-Timer] Time to assemble S-COO matrix = " << (platform::getRealTime() - assemblyTimeStart);
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_CooSort");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      info() << "[ArcaneFem-Timer] Time to assemble S-COO matrix = " << (platform::getRealTime() - assemblyTimeStart);
    }
  }

  if (m_use_coo_gpu) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooGPUBilinearOperatorTRIA3 : &FemModule::_assembleCooGPUBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    info() << "[ArcaneFem-Timer] Time to assemble COO_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Coo_Gpu");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      info() << "[ArcaneFem-Timer] Time to assemble COO_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    }
  }

  if (m_use_coo_sort_gpu) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooSortGPUBilinearOperatorTRIA3 : &FemModule::_assembleCooSortGPUBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    info() << "[ArcaneFem-Timer] Time to assemble S-COO_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_CooSort_Gpu");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      info() << "[ArcaneFem-Timer] Time to assemble S-COO_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    }
  }

  if (m_use_csr_gpu) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCsrGPUBilinearOperatorTRIA3 : &FemModule::_assembleCsrGPUBilinearOperatorTETRA4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    info() << "[ArcaneFem-Timer] Time to assemble CSR_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_Csr_Gpu");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      info() << "[ArcaneFem-Timer] Time to assemble CSR_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    }
  }

  if (m_use_nodewise_csr) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleNodeWiseCsrBilinearOperatorTria3 : &FemModule::_assembleNodeWiseCsrBilinearOperatorTetra4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    info() << "[ArcaneFem-Timer] Time to assemble NW-CSR_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_CsrNodeWise");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      info() << "[ArcaneFem-Timer] Time to assemble NW-CSR_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    }
  }

  if (m_use_buildless_csr) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleBuildLessCsrBilinearOperatorTria3 : &FemModule::_assembleBuildLessCsrBilinearOperatorTetra4;
    m_linear_system.clearValues();
    assemblyTimeStart = platform::getRealTime();
    (this->*assembly_fun)();
    info() << "[ArcaneFem-Timer] Time to assemble BL-CSR_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    if (m_cache_warming != 1)
      m_time_stats->resetStats("AssembleBilinearOperator_CsrBuildLess");
    for (auto i = 1; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      assemblyTimeStart = platform::getRealTime();
      (this->*assembly_fun)();
      info() << "[ArcaneFem-Timer] Time to assemble BL-CSR_GPU matrix = " << (platform::getRealTime() - assemblyTimeStart);
    }
  }

  // Assemble the FEM linear operator (RHS - vector b)
  assemblyTimeStart = platform::getRealTime();
  if (m_use_buildless_csr || m_use_csr_gpu || m_use_nodewise_csr || m_use_csr) {
    //_assembleCsrLinearOperator();
    _assembleCsrGpuLinearOperator();
    {
      Timer::Action timer_action(m_time_stats, "TranslateToLinearSystem");
      m_csr_matrix.translateToLinearSystem(m_linear_system, m_queue);
    }
    _translateRhs();
  }
  else {
    if (m_use_coo || m_use_coo_sort || m_use_coo_gpu || m_use_coo_sort_gpu) {
      Timer::Action timer_action(m_time_stats, "TranslateToLinearSystem");
      m_coo_matrix.translateToLinearSystem(m_linear_system);
    }
    _assembleLinearOperator();
  }
  info() << "[ArcaneFem-Timer] Time to assemble RHS vector = " << (platform::getRealTime() - assemblyTimeStart);

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
  info() << "Get material parameters...";
  f = options()->f();
  ElementNodes = 3.;

  if (options()->meshType == "TETRA4")
    ElementNodes = 4.;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_initBoundaryconditions()
{

  info() << "Apply boundary conditions";
  _applyDirichletBoundaryConditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_applyDirichletBoundaryConditionsGpu()
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
    info() << "Apply Dirichlet point condition node=" << group.name() << " v=" << value;
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
        //Original Code
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
  info() << "Assembly of FEM linear operator  ";
  info() << "Applying Dirichlet boundary condition via  penalty method ";

  // time registration
  Timer::Action timer_action(m_time_stats, "AssembleLinearOperator");

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (m_use_csr) {
  m_rhs_vect.resize(nbNode());
  m_rhs_vect.fill(0.0);

  std::array<VariableNodeByte, 1> u_dirichlet_arr = { m_u_dirichlet };

  auto method = options()->enforceDirichletMethod();
  m_bsr_format.applyDirichlet(options()->enforceDirichletMethod(), options()->penalty(), m_rhs_vect, u_dirichlet_arr, m_u, &m_linear_system);

  _translateRhs();
  }
  else {
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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";
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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrLinearOperator()
{
  info() << "Assembly of FEM linear operator ";
  info() << "Applying Dirichlet boundary condition via  penalty method for Csr";

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";
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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

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
  info() << "Assembly of FEM linear operator ";
  info() << "Applying Dirichlet boundary condition via penalty method for Csr, designed for GPU";

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void FemModule::
_translateRhs()
{
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);
  for (Int32 i = 0; i < m_rhs_vect.dim1Size(); i++) {

    rhs_values[DoFLocalId(i)] = m_rhs_vect[DoFLocalId(i)];
  }
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
/*----------------------------#endif-----------------------------------------------*/

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
  Real TimeStart;
  ITimeStats* tstat = m_time_stats;
  Timer::Action timer_action(tstat, "Solving");

  {
    TimeStart = platform::getRealTime();
    Timer::Action ta1(tstat, "LinearSystemSolve");
    m_linear_system.solve();
    info() << "[ArcaneFem-Timer] Time to solve linear system = " << (platform::getRealTime() - TimeStart);
  }

  // Re-Apply boundary conditions because the solver has modified the value
  {
    TimeStart = platform::getRealTime();
    Timer::Action ta1(tstat, "ApplyBoundaryConditions");
    _applyDirichletBoundaryConditions();
    info() << "[ArcaneFem-Timer] Time to Re-Apply boundary conditions = " << (platform::getRealTime() - TimeStart);
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
    info() << "[ArcaneFem-Timer] Time to prepare solution for post-process = " << (platform::getRealTime() - TimeStart);
  }

  TimeStart = platform::getRealTime();
  m_u.synchronize();
  info() << "[ArcaneFem-Timer] Time to synchronize solution across subdomains = " << (platform::getRealTime() - TimeStart);

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
  String filename = options()->resultFile();
  info() << "CheckResultFile filename=" << filename;
  if (filename.empty())
    return;
  const double epsilon = 1.0e-4;
  const double skipValuesMinLim = 1.0e-16;
  checkNodeResultFile(traceMng(), filename, m_u, epsilon, skipValuesMinLim);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

bool FemModule::
_isMasterRank() const
{
  return parallelMng()->isMasterIO();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
