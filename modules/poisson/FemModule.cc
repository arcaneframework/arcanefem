// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Poisson solver module of ArcaneFEM.                                       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModule at the start of the simulation.
 *
 * This method initializes degrees of freedom (DoFs) on nodes.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "[ArcaneFem-Info] Started module startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();

  m_matrix_format = options()->matrixFormat();
  _handleCommandLineFlags();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs the main computation for the FemModule.
 *
 * This method:
 *   1. Stops the time loop after 1 iteration since the equation is steady state.
 *   2. Resets, configures, and initializes the linear system.
 *   3. Executes the stationary solve.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "[ArcaneFem-Info] Started module compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  info() << "[ArcaneFem-Info] Matrix format used " << m_matrix_format;
  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");

  if (m_petsc_flags != NULL)
    _setPetscFlagsFromCommandline();

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") {
    auto use_csr_in_linear_system = options()->linearSystem.serviceName() == "HypreLinearSystem";
    if (m_matrix_format == "BSR")
      m_bsr_format.initialize(mesh(), use_csr_in_linear_system, 0);
    else
      m_bsr_format.initialize(mesh(), use_csr_in_linear_system, 1);
    m_bsr_format.computeSparsity();
  }

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows a sequence of steps to solve FEM system:
 *
 *   1. _getMaterialParameters()     Retrieves material parameters via
 *   2. _assembleBilinearOperator()  Assembles the FEM  matrix A
 *   3. _assembleLinearOperator()    Assembles the FEM RHS vector b
 *   4. _solve()                     Solves for solution vector u = A^-1*b
 *   5. _updateVariables()           Updates FEM variables u = x
 *   6. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  _getMaterialParameters();

  if(m_assemble_linear_system){
    _assembleBilinearOperator();
    _assembleLinearOperator();
  }

  if(m_solve_linear_system){
    _solve();
    _updateVariables();
  }

  if(m_cross_validation){
    _validateResults();
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Retrieves and sets the material parameters for the simulation.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  f = options()->f();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the FEM linear operator for the current simulation step.
 *
 * This method constructs the right-hand side (RHS) vector by calling theq
 * appropriate assembly function based on the execution context:
 *  - CPU-exclusive execution: Calls _assembleLinearOperatorCpu().
 *  - CPU or GPU execution: Calls _assembleLinearOperatorGpu().
 *
 */
 /*---------------------------------------------------------------------------*/

void FemModule::_assembleLinearOperator()
{
  if (options()->linearSystem.serviceName() == "HypreLinearSystem")
    _assembleLinearOperatorGpu();
  else
    _assembleLinearOperatorCpu();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief FEM linear operator for the current simulation step.
 * GPU compatible. Currently working with HypreDoFLinearSystem.
 *
 * This method constructs the linear  system by  assembling the LHS matrix
 * and  RHS vector, applying various boundary conditions and source terms.
 *
 * Steps involved:
 *  1. The RHS vector is initialized to zero before applying any conditions.
 *  2. If Neumann BC are specified applied to the RHS.
 *  3. If Dirichlet BC/Point are specified apply to the LHS & RHS.
 */
/*---------------------------------------------------------------------------*/

void FemModule::_assembleLinearOperatorGpu()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperatorGpu()";
  Real elapsedTime = platform::getRealTime();

  auto& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  // Because of Dirichlet (Penalty) implementation in Hypre.
  m_linear_system.getForcedInfo().fill(false);

  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  auto mesh_ptr = mesh();

  auto applyBoundaryConditions = [&](auto BCFunctions) {
    if (options()->f.isPresent())
      BCFunctions.applyConstantSourceToRhs(f, m_dofs_on_nodes, m_node_coord, rhs_values, mesh_ptr, queue);

    BC::IArcaneFemBC* bc = options()->boundaryConditions();

    if (bc) {
      for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions())
        BCFunctions.applyNeumannToRhs(bs, m_dofs_on_nodes, m_node_coord, rhs_values, mesh_ptr, queue);

      for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
        FemUtils::Gpu::BoundaryConditions::applyDirichletViaPenalty(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);

      for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions())
        FemUtils::Gpu::BoundaryConditions::applyPointDirichletViaPenalty(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);
    }
  };

  if (mesh()->dimension() == 3)
    applyBoundaryConditions(FemUtils::Gpu::BoundaryConditions3D());
  else
    applyBoundaryConditions(FemUtils::Gpu::BoundaryConditions2D());

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] rhs-vector-assembly-gpu", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief FEM linear operator for the current simulation step.
 *
 * This method constructs the linear  system by  assembling the LHS matrix
 * and  RHS vector, applying various boundary conditions and source terms.
 *
 * Steps involved:
 *  1. The RHS vector is initialized to zero before applying any conditions.
 *  2. If Neumann BC are specified applied to the RHS.
 *  3. If Dirichlet BC/Point are specified apply to the LHS & RHS.
 */
/*---------------------------------------------------------------------------*/

void FemModule::_assembleLinearOperatorCpu()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperatorCpu()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // Helper lambda to apply boundary conditions
  auto applyBoundaryConditions = [&](auto BCFunctions) {
    if (options()->f.isPresent())
      BCFunctions.applyConstantSourceToRhs(f, mesh(), node_dof, m_node_coord, rhs_values);

    BC::IArcaneFemBC* bc = options()->boundaryConditions();
    if (bc) {
      for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions())
        BCFunctions.applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);

      for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
        BCFunctions.applyDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);

      for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions())
        BCFunctions.applyPointDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);
    }
  };

  // Apply the correct boundary conditions based on mesh dimension
  if (mesh()->dimension() == 3) {
    using BCFunctions = ArcaneFemFunctions::BoundaryConditions3D;
    applyBoundaryConditions(BCFunctions());
  }
  else {
    using BCFunctions = ArcaneFemFunctions::BoundaryConditions2D;
    applyBoundaryConditions(BCFunctions());
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly given as mesh type.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto queue = subDomain()->acceleratorMng()->defaultQueue();
    auto command = makeCommand(queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);

    if (mesh()->dimension() == 2)
      m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return _computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord); });
    else
      m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return _computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord); });

    m_bsr_format.toLinearSystem(m_linear_system);
  }

  if (m_matrix_format == "AF-BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto queue = subDomain()->acceleratorMng()->defaultQueue();
    auto command = makeCommand(queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);

    if (mesh()->dimension() == 2)
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return _computeElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, node_lid); });
    else
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return _computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, node_lid); });

    m_bsr_format.toLinearSystem(m_linear_system);
  }

  if (m_matrix_format == "DOK") {
    if (mesh()->dimension() == 3)
      _assembleBilinear<4>([this](const Cell& cell) {
        return _computeElementMatrixTetra4(cell);
      });
    else
      _assembleBilinear<3>([this](const Cell& cell) {
        return _computeElementMatrixTria3(cell);
      });
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] lhs-matrix-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system.
 *
 * The method performs the following steps:
 *   1. For each cell, retrieves the cell-specific constant `lambda`.
 *   2. Computes element matrix using provided `compute_element_matrix` function.
 *   3. Assembles global matrix by adding contributions from each cell's element
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

template <int N>
void FemModule::
_assembleBilinear(const std::function<FixedMatrix<N, N>(const Cell&)>& compute_element_matrix)
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    auto K_e = compute_element_matrix(cell); // element matrix based on the provided function
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        Real v = K_e(n1_index, n2_index);
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
/**
 * @brief Solves the linear system.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  info() << "[ArcaneFem-Info] Started module _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.solve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] solve-linear-system", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Update the FEM variables.
 *
 * This method performs the following actions:
 *   1. Fetches values of solution from solved linear system to FEM variables,
 *      i.e., it copies RHS DOF to u.
 *   2. Performs synchronize of FEM variables across subdomains.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_updateVariables()
{
  info() << "[ArcaneFem-Info] Started module _updateVariables()";
  Real elapsedTime = platform::getRealTime();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      m_u[node] = dof_u[node_dof.dofId(node, 0)];
    }
  }

  m_u.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] update-variables", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Validates and prints the results of the FEM computation.
 *
 * This method performs the following actions:
 *   1. If number of nodes < 200, prints the computed values for each node.
 *   2. Retrieves the filename for the result file from options.
 *   3. If a filename is provided, checks the computed results against result file.
 *
 * @note The result comparison uses a tolerance of 1.0e-4.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_validateResults()
{
  info() << "[ArcaneFem-Info] Started module _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200)
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "u[" << node.localId() << "][" << node.uniqueId() << "] = " << m_u[node];
    }

  String filename = options()->resultFile();
  info() << "ValidateResultFile filename=" << filename;

  if (!filename.empty())
    checkNodeResultFile(traceMng(), filename, m_u, 1.0e-4);

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Function to prints the execution time `value` of phase `label`
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_printArcaneFemTime(const String label, const Real value)
{
  info() << std::left << std::setw(40) << label << " = " << value;
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
/**
 * @brief Function to hande commandline flags
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_handleCommandLineFlags()
{
  info() << "[ArcaneFem-Info] Started module _handleCommandLineFlags()";
  ParameterList parameter_list = this->subDomain()->application()->applicationInfo().commandLineArguments().parameters();

  if (parameter_list.getParameterOrNull("assemble_linear_system") == "FALSE") {
    m_assemble_linear_system = false;
    info() << "[ArcaneFem-Info] Linear system not assembled (assemble_linear_system = FALSE)";
  }

  if (parameter_list.getParameterOrNull("solve_linear_system") == "FALSE") {
    m_solve_linear_system = false;
    m_cross_validation = false;
    info() << "[ArcaneFem-Info] Linear system assembled but not solved (solve_linear_system = FALSE)";
  }

  if (parameter_list.getParameterOrNull("cross_validation") == "FALSE") {
    m_cross_validation = false;
    info() << "[ArcaneFem-Info] Cross validation disabled (cross_validation = FALSE)";
  }

  m_petsc_flags = parameter_list.getParameterOrNull("petsc_flags");
  if (m_petsc_flags != NULL) {
    info() << "[ArcaneFem-Info] PETSc flags the user provided will be used (petsc_flags != NULL)";
  }

  String matrix_format_from_commandline = parameter_list.getParameterOrNull("matrix_format");
  if (matrix_format_from_commandline != NULL){
    m_matrix_format = matrix_format_from_commandline;
    info() << "[ArcaneFem-Info] Using commandline format for matrix format (matrix_format != NULL)";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
