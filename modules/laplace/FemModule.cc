// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* FEM code to solve Laplace problem.                                        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include "ElementMatrix.h"
#include "ElementMatrixHexQuad.h"

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
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->hasSolutionComparisonFile();
  m_petsc_flags = options()->petscFlags();
  m_hex_quad_mesh = options()->hexQuadMesh();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"initialize", elapsedTime);
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

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");

  if (m_petsc_flags != NULL){
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") {
    bool use_csr_in_linear_system =
    options()->linearSystem.serviceName() == "HypreLinearSystem" ||
    options()->linearSystem.serviceName() == "PETScLinearSystem";
    if (m_matrix_format == "BSR")
      m_bsr_format.initialize(mesh(), 1, use_csr_in_linear_system, 0);
    else
      m_bsr_format.initialize(mesh(), 1, use_csr_in_linear_system, 1);
    m_bsr_format.computeSparsity();
  }

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows via the following steps:
 *   1. _assembleBilinearOperator()  Assembles the FEM  matrix 𝐀
 *   2. _assembleLinearOperator()    Assembles the FEM RHS vector 𝐛
 *   3. _solve()                     Solves for solution vector 𝐮 = 𝐀⁻¹𝐛
 *   4. _updateVariables()           Updates FEM variables 𝐮 = 𝐱
 *   5. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
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
 * @brief Assembles the FEM linear operator for the current simulation step.
 *
 * This method constructs the right-hand side (RHS) vector by calling then
 * appropriate assembly function based on the execution context:
 *  - CPU-exclusive execution: Calls _assembleLinearOperatorCpu().
 *  - CPU or GPU execution: Calls _assembleLinearOperatorGpu().
 *
 */
/*---------------------------------------------------------------------------*/

void FemModule::_assembleLinearOperator()
{
  if (options()->linearSystem.serviceName() == "HypreLinearSystem" ||
      options()->linearSystem.serviceName() == "PETScLinearSystem")
    _assembleLinearOperatorGpu();
  else
    _assembleLinearOperatorCpu();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief FEM linear operator for the current simulation step.
 *
 * This method constructs the FEM linear  systems RHS vector by applying
 * various boundary conditions and source terms.
 *
 * Steps involved:
 *  1. The RHS vector is initialized to zero before applying any conditions.
 *  2. If Neumann BC are specified applied to the RHS.
 *  3. If Dirichlet BC/Point are specified apply to the LHS & RHS.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperatorCpu()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR")
    m_bsr_format.toLinearSystem(m_linear_system);

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  BC::IArcaneFemBC* bc = options()->boundaryConditions();
  if (bc) {
    for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions()) {
      if (mesh()->dimension() == 2) {
        if (m_hex_quad_mesh)
          ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhsQuad4(bs, node_dof, m_node_coord, rhs_values);
        else
          ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);
      }
      if (mesh()->dimension() == 3) {
        if (m_hex_quad_mesh)
          ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhsHexa8(bs, node_dof, m_node_coord, rhs_values);
        else
          ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);
      }
    }

    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
      ArcaneFemFunctions::BoundaryConditions::applyDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);

    for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions())
      ArcaneFemFunctions::BoundaryConditions::applyPointDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"rhs-vector-assembly", elapsedTime);
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

  m_bsr_format.toLinearSystem(m_linear_system);

  auto& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  // Because of Dirichlet (Penalty) implementation in Hypre.
  m_linear_system.getForcedInfo().fill(false);

  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  auto mesh_ptr = mesh();

  auto applyBoundaryConditions = [&](auto BCFunctions) {
    BC::IArcaneFemBC* bc = options()->boundaryConditions();

    if (bc) {
      for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions())
        BCFunctions.applyNeumannToRhs(bs, m_dofs_on_nodes, m_node_coord, rhs_values, mesh_ptr, queue);

      for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
        FemUtils::Gpu::BoundaryConditions::applyDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);

      for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions())
        FemUtils::Gpu::BoundaryConditions::applyPointDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);
    }
  };

  if (mesh()->dimension() == 3)
    applyBoundaryConditions(FemUtils::Gpu::BoundaryConditions3D());
  else
    applyBoundaryConditions(FemUtils::Gpu::BoundaryConditions2D());

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"rhs-vector-assembly-gpu", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto m_queue = subDomain()->acceleratorMng()->defaultQueue();
    auto command = makeCommand(m_queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    if (m_matrix_format == "BSR") {
      if (mesh()->dimension() == 2)
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord); });
      if (mesh()->dimension() == 3)
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord); });
    }
    if (m_matrix_format == "AF-BSR") {
      if (mesh()->dimension() == 2)
        m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, node_lid); });
      if (mesh()->dimension() == 3)
        m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, node_lid); });
    }
  }

  if (m_matrix_format == "DOK") {
    if (mesh()->dimension() == 3)
      if(m_hex_quad_mesh)
        _assembleBilinear<8>([this](const Cell& cell) { return _computeElementMatrixHexa8(cell); });
      else
        _assembleBilinear<4>([this](const Cell& cell) { return _computeElementMatrixTetra4(cell); });
    else
      if(m_hex_quad_mesh)
        _assembleBilinear<4>([this](const Cell& cell) { return _computeElementMatrixQuad4(cell); });
      else
        _assembleBilinear<3>([this](const Cell& cell) { return _computeElementMatrixTria3(cell); });
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"lhs-matrix-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system.
 *
 * The method performs the following steps:
 *   1. Computes element matrix using provided `compute_element_matrix` function.
 *   2. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

template <int N>
void FemModule::
_assembleBilinear(const std::function<RealMatrix<N, N>(const Cell&)>& compute_element_matrix)
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
 * @brief Solves the linear system and updates the solution vector.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  info() << "[ArcaneFem-Info] Started module _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.applyLinearSystemTransformationAndSolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"solve-linear-system", elapsedTime);
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
  info() << "[ArcaneFem-Module] _updateVariables()";
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
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"update-variables", elapsedTime);
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
      info() << "u[" << node.uniqueId() << "] = " << m_u[node];
    }

  String filename = options()->solutionComparisonFile();
  const double epsilon = 1.0e-4;
  const double min_value_to_test = 1.0e-16;

  checkNodeResultFile(traceMng(), filename, m_u, epsilon, min_value_to_test);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/