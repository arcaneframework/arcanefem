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
#include <arcane/core/MeshUtils.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_writeInJson()
{
  if (!_isMasterRank())
    return;
  ofstream jsonFile("time.json");
  JSONWriter json_writer(JSONWriter::FormatFlags::None);
  json_writer.beginObject();
  {
    JSONWriter::Object jo(json_writer, "Timer");
    m_time_stats->dumpStatsJSON(json_writer);
  }
  json_writer.endObject();
  jsonFile << json_writer.getBuffer();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_readTimeFromJson(String main_time, String sub_time)
{
  UniqueArray<Byte> bytes;
  IParallelMng* pm = this->mesh()->parallelMng();
  pm->ioMng()->collectiveRead("time.json", bytes, false);
  JSONDocument json_doc;
  json_doc.parse(bytes, "time.json");
  //Parsing through the JSON
  JSONValue root = json_doc.root();
  //From root to the list of subactions in Main
  JSONValueList main = root.child("Timer").child("Current").child("Main").child("SubActions").valueAsArray();
  //From the list of subactions in Main to the list of subactions in Loop
  JSONValueList loop = (main.begin() + 3)->child("SubActions").valueAsArray();
  //From the list of subactions in Loop to the list of subactions in LoopEntryPoints
  JSONValueList loopEntryPoint = (loop.begin() + 1)->child("SubActions").valueAsArray();
  //From the list of subactions in LoopEntryPoints to the list of subactions in Fem
  JSONValueList fem = (loopEntryPoint.begin() + 7)->child("SubActions").valueAsArray();
  //From the list of subactions in Fem to the list of subactions in Compute
  JSONValueList compute = (fem.begin() + 1)->child("SubActions").valueAsArray();
  //From the list of subactions in Compute to the list of subactions in StationarySolve
  String prev = "";
  JSONValueList stationarySolve;
  for (JSONValue el : compute) {
    if (prev == "StationarySolve") {
      stationarySolve = el.child("SubActions").valueAsArray();
      break;
    }
    prev = el.valueAsStringView();
  }
  //Selecting the right 'main' action
  JSONValue function;
  prev = "";
  for (JSONValue el : stationarySolve) {
    if (prev == main_time) {
      function = el;
      break;
    }
    prev = el.valueAsStringView();
  }
  //Selecting the sub action if we want it
  if (sub_time != "") {
    prev = "";
    for (JSONValue el : function.child("SubActions").valueAsArray()) {
      if (prev == sub_time) {
        function = el;
        break;
      }
      prev = el.valueAsStringView();
    }
  }
  // The timer has not been found
  if (prev == "") {
    return 0;
  }
  String val;
  std::stringstream ss;
  //Get only the Cumulative value
  ss << function.child("Cumulative").value();
  ss >> val;
  return *Convert::Type<Real>::tryParse(val);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_saveTimeInCSV()
{
  std::ofstream csv_save;
  String csv_file_name = String::format("time.{0}.csv",parallelMng()->commRank());
  if (!fs::exists(csv_file_name.localstr())) {
    csv_save.open(csv_file_name.localstr());
    csv_save << "Number of Nodes,Legacy,COO with sorting,COO,CSR,CSR made for GPU,Node Wise CSR made for GPU,BLCSR made for GPU,CSR GPU,Node Wise CSR GPU,BLCSR GPU\n";
  }
  else {
    csv_save.open(csv_file_name.localstr(), std::ios_base::app);
  }
  Integer denume = m_cache_warming;
  if (denume > 1)
    denume--;
  csv_save << nbNode() << ",";
  csv_save << _readTimeFromJson("AssembleLegacyBilinearOperatorTria3", "") / denume << ",";
  csv_save << _readTimeFromJson("AssembleCooSortBilinearOperatorTria3", "") / denume << ",";
  csv_save << _readTimeFromJson("AssembleCooBilinearOperatorTria3", "") / denume << ",";
  csv_save << _readTimeFromJson("AssembleCsrBilinearOperatorTria3", "") / denume << ",";
  if (m_running_on_gpu) {
    csv_save << "0,0,0,";
    csv_save << _readTimeFromJson("AssembleCsrGpuBilinearOperatorTria3", "") / denume << ",";
    csv_save << _readTimeFromJson("AssembleNodeWiseCsrBilinearOperatorTria3", "") / denume << ",";
    csv_save << _readTimeFromJson("AssembleBuildLessCsrBilinearOperatorTria3", "") / denume;
  }
  else {
    csv_save << _readTimeFromJson("AssembleCsrGpuBilinearOperatorTria3", "") / (m_cache_warming == 1 ? 1 : m_cache_warming - 1) << ",";
    csv_save << _readTimeFromJson("AssembleNodeWiseCsrBilinearOperatorTria3", "") / denume << ",";
    csv_save << _readTimeFromJson("AssembleBuildLessCsrBilinearOperatorTria3", "") / denume << ",";
    csv_save << "0,0,0";
  }
  csv_save << "\n";
  csv_save.close();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_saveNoBuildTimeInCSV()
{
  std::ofstream csv_save;
  String csv_file_name = String::format("timeNoBuild.{0}.csv",parallelMng()->commRank());
  if (!fs::exists(csv_file_name.localstr())) {
    csv_save.open(csv_file_name.localstr());
    csv_save << "Number of Nodes,Legacy,COO with sorting,COO,CSR,CSR made for GPU,Node Wise CSR made for GPU,BLCSR made for GPU,CSR GPU,Node Wise CSR GPU,BLCSR GPU\n";
  }
  else {
    csv_save.open(csv_file_name.localstr(), std::ios_base::app);
  }
  Integer denume = m_cache_warming;
  if (denume > 1)
    denume--;
  csv_save << nbNode() << ",";
  csv_save << _readTimeFromJson("AssembleLegacyBilinearOperatorTria3", "") / denume << ",";
  csv_save << (_readTimeFromJson("AssembleCooSortBilinearOperatorTria3", "CooSortComputeElementMatrixTria3") + _readTimeFromJson("AssembleCooSortBilinearOperatorTria3", "CooSortAddToGlobalMatrix")) / denume << ",";
  csv_save << (_readTimeFromJson("AssembleCooBilinearOperatorTria3", "CooComputeElementMatrixTria3") + _readTimeFromJson("AssembleCooBilinearOperatorTria3", "CooAddToGlobalMatrix")) / denume << ",";
  csv_save << (_readTimeFromJson("AssembleCsrBilinearOperatorTria3", "CsrComputeElementMatrixTria3") + _readTimeFromJson("AssembleCsrBilinearOperatorTria3", "CsrAddToGlobalMatrix")) / denume << ",";
  if (m_running_on_gpu) {
    csv_save << "0,0,0,";
    csv_save << _readTimeFromJson("AssembleCsrGpuBilinearOperatorTria3", "CsrGpuAddComputeLoop") / denume << ",";
    csv_save << _readTimeFromJson("AssembleNodeWiseCsrBilinearOperatorTria3", "NodeWiseCsrAddAndCompute") / denume << ",";
    csv_save << _readTimeFromJson("AssembleBuildLessCsrBilinearOperatorTria3", "BuildLessCsrAddAndCompute") / denume;
  }
  else {
    csv_save << _readTimeFromJson("AssembleCsrGpuBilinearOperatorTria3", "CsrGpuAddComputeLoop") / denume << ",";
    csv_save << _readTimeFromJson("AssembleNodeWiseCsrBilinearOperatorTria3", "NodeWiseCsrAddAndCompute") / denume << ",";
    csv_save << _readTimeFromJson("AssembleBuildLessCsrBilinearOperatorTria3", "BuildLessCsrAddAndCompute") / denume << ",";
    csv_save << "0,0,0";
  }
  csv_save << "\n";
  csv_save.close();
}

void FemModule::
_benchBuildRow()
{
  std::ofstream csv_save;
  String csv_file_name = String::format("buildRow.{0}.csv",parallelMng()->commRank());
  if (!fs::exists(csv_file_name.localstr())) {
    csv_save.open(csv_file_name.localstr());
    csv_save << "Number of Nodes,Build on CPU,Build on GPU\n";
  }
  else {
    csv_save.open(csv_file_name.localstr(), std::ios_base::app);
  }
  csv_save << nbNode() << ",";
  csv_save << _readTimeFromJson("AssembleBuildLessCsrBilinearOperatorTria3", "BuildLessCsrBuildMatrix") / m_cache_warming << ",";
  csv_save << _readTimeFromJson("AssembleBuildLessCsrBilinearOperatorTria3", "BuildLessCsrBuildMatrixGPU") / m_cache_warming << "\n";
  csv_save.close();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
endModule()
{
  _writeInJson();
  _saveTimeInCSV();
  _saveNoBuildTimeInCSV();
  //_benchBuildRow();
}

void FemModule::
compute()
{
  info() << "Module Fem COMPUTE";

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

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
  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();

  _doStationarySolve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "Module Fem INIT";

  if (m_register_time && _isMasterRank()) {
    logger = ofstream("timer.txt");
    wbuild = ofstream("with_build.csv", std::ios_base::app);
    wbuild << nbNode() << ",";
    timer = ofstream("timer.csv", std::ios_base::app);
    timer << nbNode() << ",";
  }

  {
    IMesh* mesh = defaultMesh();
    // If we do not create edges, we need to create custom connectivity
    // to store the neighbours node of a node
    if (mesh->dimension()==3 && !options()->createEdges()){
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
    }
    else{
      m_nb_edge = mesh->nbEdge();
      info() << "Number of edge: nb_edge=" << m_nb_edge;
    }
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
  info() << "The time will be registered by arcane in the output/listing/logs.0 file, and will be added to (or will create) the time.csv (with time for the various bilinear assembly phases) and timeNoBuild.csv (with time without the building part of COO and CSR for the various bilinear assembly phases) fil";
  if (parameter_list.getParameterOrNull("REGISTER_TIME") == "TRUE") {
    m_register_time = true;
    info() << "REGISTER_TIME: The time will also be registered in the timer.txt, with_build.csv and timer.csv file";
  }
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
    info() << "COO: The COO datastructure and its associated methods will be used";
  }
  if (parameter_list.getParameterOrNull("COO_SORT") == "TRUE" || options()->cooSorting()) {
    m_use_coo_sort = true;
    m_use_legacy = false;
    info() << "COO_SORT: The COO with sorting datastructure and its associated methods will be used";
  }
  if (parameter_list.getParameterOrNull("COO_GPU") == "TRUE" || options()->cooGpu()) {
    m_use_coo_gpu = true;
    m_use_legacy = false;
    info() << "COO_GPU: The COO datastructure GPU comptaible and its associated methods will be used";
  }
  if (parameter_list.getParameterOrNull("COO_SORT_GPU") == "TRUE" || options()->cooSortingGpu()) {
    m_use_coo_sort_gpu = true;
    m_use_legacy = false;
    info() << "COO_SORT_GPU: The COO with sorting datastructure GPU comptaible and its associated methods will be used";
  }
  if (parameter_list.getParameterOrNull("CSR") == "TRUE" || options()->csr()) {
    m_use_csr = true;
    m_use_legacy = false;
    info() << "CSR: The CSR datastructure and its associated methods will be used";
  }
#ifdef ARCANE_HAS_ACCELERATOR
  if (parameter_list.getParameterOrNull("CSR_GPU") == "TRUE" || options()->csrGpu()) {
    m_use_csr_gpu = true;
    m_use_legacy = false;
    info() << "CSR_GPU: The CSR datastructure GPU compatible and its associated methods will be used";
  }
#endif
  if (parameter_list.getParameterOrNull("NWCSR") == "TRUE" || options()->nwcsr()) {
    m_use_nodewise_csr = true;
    m_use_legacy = false;
    info() << "NWCSR: The Csr datastructure (GPU compatible) and its associated methods will be used with computation in a nodewise manner";
  }
  if (parameter_list.getParameterOrNull("BLCSR") == "TRUE" || options()->blcsr()) {
    m_use_buildless_csr = true;
    m_use_legacy = false;
    info() << "BLCSR: The Csr datastructure (GPU compatible) and its associated methods will be used with computation in a nodewise manner with the building phases incorporated in the computation";
  }
  if (parameter_list.getParameterOrNull("LEGACY") == "TRUE" || m_use_legacy || options()->legacy()) {
    m_use_legacy = true;
    info() << "LEGACY: The Legacy datastructure and its associated methods will be used";
  }
  else if (parameter_list.getParameterOrNull("LEGACY") == "FALSE" || options()->legacy()) {
    m_use_legacy = false;
  }
  if (parameter_list.getParameterOrNull("AcceleratorRuntime") == "cuda") {
    m_running_on_gpu = true;
    info() << "CUDA: The methods able to use GPU will use it";
  }
  info() << "-----------------------------------------------------------------------------------------";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  Timer::Action timer_action(m_time_stats, "StationarySolve");

  _getMaterialParameters();

  auto dim = mesh()->dimension();

  // Assemble the FEM bilinear operator (LHS - matrix A)
  if (m_use_legacy) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleBilinearOperatorTRIA3 : &FemModule::_assembleBilinearOperatorTETRA4;
    for (auto i = 0; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      m_time_stats->resetStats("AssembleBilinearOperator_Legacy");
      (this->*assembly_fun)();
    }
  }

  if (m_use_csr) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCsrBilinearOperatorTRIA3 : &FemModule::_assembleCsrBilinearOperatorTETRA4;
    for (auto i = 0; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      m_time_stats->resetStats("AssembleBilinearOperator_Csr");
      (this->*assembly_fun)();
    }
  }

  if (m_use_coo) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooBilinearOperatorTRIA3 : &FemModule::_assembleCooBilinearOperatorTETRA4;
    for (auto i = 0; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      m_time_stats->resetStats("AssembleBilinearOperator_Coo");
      (this->*assembly_fun)();
    }
    m_coo_matrix.translateToLinearSystem(m_linear_system);
  }

  if (m_use_coo_sort) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooSortBilinearOperatorTRIA3 : &FemModule::_assembleCooSortBilinearOperatorTETRA4;
    for (auto i = 0; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      m_time_stats->resetStats("AssembleBilinearOperator_CooSort");
      (this->*assembly_fun)();
    }
    m_coo_matrix.translateToLinearSystem(m_linear_system);
  }

  if (m_use_coo_gpu || m_use_coo_sort_gpu) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCooGPUBilinearOperatorTRIA3 : &FemModule::_assembleCooGPUBilinearOperatorTETRA4;
    for (auto i = 0; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      m_time_stats->resetStats("AssembleBilinearOperator_Coo_GPU");
      m_time_stats->resetStats("AssembleBilinearOperator_CooSort_GPU");
      (this->*assembly_fun)();
    }
    m_coo_matrix.translateToLinearSystem(m_linear_system);
  }

  if (m_use_csr_gpu) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleCsrGPUBilinearOperatorTRIA3 : &FemModule::_assembleCsrGPUBilinearOperatorTETRA4;
    for (auto i = 0; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      m_time_stats->resetStats("AssembleBilinearOperator_Csr_GPU");
      (this->*assembly_fun)();
    }
  }

  if (m_use_nodewise_csr) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleNodeWiseCsrBilinearOperatorTria3 : &FemModule::_assembleNodeWiseCsrBilinearOperatorTetra4;
    for (auto i = 0; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      m_time_stats->resetStats("AssembleBilinearOperator_CsrNodeWise");
      (this->*assembly_fun)();
    }
  }

  if (m_use_buildless_csr) {
    void (FemModule::*assembly_fun)() = dim == 2 ? &FemModule::_assembleBuildLessCsrBilinearOperatorTria3 : &FemModule::_assembleBuildLessCsrBilinearOperatorTetra4;
    for (auto i = 0; i < m_cache_warming; ++i) {
      m_linear_system.clearValues();
      m_time_stats->resetStats("AssembleBilinearOperator_CsrBuildLess");
      (this->*assembly_fun)();
    }
  }

  // Assemble the FEM linear operator (RHS - vector b)
  if (m_use_buildless_csr || m_use_csr_gpu || m_use_nodewise_csr || m_use_csr) {
    _assembleCsrGpuLinearOperator();
    //_assembleCsrLinearOperator();
    m_csr_matrix.translateToLinearSystem(m_linear_system);
    _translateRhs();
  }
  else {
    _assembleLinearOperator();
  }

  _solve();

  // Check results
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
  info() << "Init boundary conditions...";

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
_assembleLinearOperator()
{
  info() << "Assembly of FEM linear operator  ";
  info() << "Applying Dirichlet boundary condition via  penalty method ";

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    Real Penalty = options()->penalty(); // 1.0e30 is the default

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u_dirichlet[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        // This SetValue should be updated in the acoording format we have (such as COO or CSR)
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

  {
    Timer::Action timer_action(m_time_stats, "ConstantSourceTermAssembly");
    //----------------------------------------------
    // Constant source term assembly
    //----------------------------------------------
    //
    //  $int_{Omega}(f*v^h)$
    //  only for noded that are non-Dirichlet
    //----------------------------------------------
    if (options()->meshType == "TRIA3"){
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

    if (options()->meshType == "TETRA4"){
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
    Timer::Action timer_action(m_time_stats, "ConstantSourceTermAssembly");

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

  if (options()->meshType == "TRIA3")
  {
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

  if (options()->meshType == "TETRA4")
  {
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
  Real3 dPhi0 = Arcane::math::cross(m2 - m1, m1 - m3) ;
  Real3 dPhi1 = Arcane::math::cross(m3 - m0, m0 - m2) ;
  Real3 dPhi2 = Arcane::math::cross(m1 - m0, m0 - m3) ;
  Real3 dPhi3 = Arcane::math::cross(m0 - m1, m1 - m2) ;

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
_assembleBilinearOperatorTETRA4()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    auto K_e = _computeElementMatrixTETRA4(cell);  // element stiffness matrix

    //             # assemble elementary matrix into the global one
    //             # elementary terms are positioned into K according
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

void FemModule::
_solve()
{
  ITimeStats* tstat = m_time_stats;
  Timer::Action timer_action(tstat, "Solving");

  {
    Timer::Action ta1(tstat, "LinearSystemSolve");
    // # T=linalg.solve(K,RHS)
    m_linear_system.solve();
  }

  // Re-Apply boundary conditions because the solver has modified the value
  // of u on all nodes
  {
    Timer::Action ta1(tstat, "ApplyBoundaryConditions");
    _applyDirichletBoundaryConditions();
  }

  {
    Timer::Action ta1(tstat, "CopySolution");
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
  Connectivity c(mesh()->connectivity());
  if (options()->meshType == "TETRA4" && options()->createEdges()){
    info() << "Adding edge connectivity";
    c.enableConnectivity(Connectivity::CT_HasEdge);
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

void FemModule::fileNumArray(bool ref, NumArray<Real, MDDim1> numarray)
{
  ofstream file;
  if (ref)
    file.open("ref.txt");
  else
    file.open("test.txt");
  for (auto i = 0; i < numarray.dim1Size(); i++) {
    file << numarray(i) << " ";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
