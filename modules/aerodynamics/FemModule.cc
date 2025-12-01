// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Simple module to solve aerodynamics using FEM.                            */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include "ElementMatrix.h"
#include "ElementMatrixHexQuad.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModule at the start of the simulation.
 *
 * - This method initializes degrees of freedom (DoFs) on nodes.
 * - It also gets values of some solver parameters.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "[ArcaneFem-Info] Started module startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dofs_on_nodes.initialize(mesh(), 1);

  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->hasSolutionComparisonFile();
  m_petsc_flags = options()->petscFlags();
  m_hex_quad_mesh = options()->hexQuadMesh();

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
 *   3. Applies PETSc commandline flags to the solver (if PETSc is used).
 *   4. Executes the stationary solve and extracts psi.
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

  if (m_petsc_flags != NULL) {
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  _doStationarySolve();
  _getPsi();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getPsi()
{
  info() << "[ArcaneFem-Info] Started module _getPsi()";
  Real elapsedTime = platform::getRealTime();

  if (mesh()->dimension() == 2)
    if (m_hex_quad_mesh)
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real3 grad = ArcaneFemFunctions::FeOperation2D::computeGradientQuad4(cell, m_node_coord, m_u);
        m_psi[cell] = -grad.x * grad.x - grad.y * grad.y;
      }
    else
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real3 grad = ArcaneFemFunctions::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_u);
        m_psi[cell] = -grad.x * grad.x - grad.y * grad.y;
      }

  if (mesh()->dimension() == 3)
    if (m_hex_quad_mesh)
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real3 grad = ArcaneFemFunctions::FeOperation3D::computeGradientHexa8(cell, m_node_coord, m_u);
        m_psi[cell] = -grad.x * grad.x - grad.y * grad.y - grad.z * grad.z;
      }
    else
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real3 grad = ArcaneFemFunctions::FeOperation3D::computeGradientTetra4(cell, m_node_coord, m_u);
        m_psi[cell] = -grad.x * grad.x - grad.y * grad.y - grad.z * grad.z;
      }

  m_psi.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "update-psi", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows a sequence of steps to solve FEM system:
 *   1. _assembleBilinearOperator()       Assembles the FEM  matrix 𝑨
 *   2. _assembleLinearOperator()         Assembles the FEM RHS vector 𝒃
 *   3. _solve()                          Solves for solution vector 𝒖 = 𝑨⁻¹𝒃
 *   4. _updateVariables()                Updates FEM variables 𝒖 = 𝒙
 *   5. _validateResults()                Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  if (m_assemble_linear_system) {
    _assembleBilinearOperator();
    _assembleLinearOperator();
  }
  if (m_solve_linear_system) {
    _solve();
    _updateVariables();
  }
  if (m_cross_validation) {
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
  if (options()->linearSystem.serviceName() == "HypreLinearSystem")
    _assembleLinearOperatorGpu();
  else
    _assembleLinearOperatorCpu();
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
 *  2. Dirichlet BC/Point are specified apply to the LHS & RHS.
 *  3. Far-field BC are applied to the LHS & RHS.
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
  if (bc){
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
      ArcaneFemFunctions::BoundaryConditions::applyDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);
  }

  for (const auto& bs : options()->farfieldBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real angle = bs->angle();
    Real penalty = options()->penalty();

    ENUMERATE_ (Face, iface, group) {
      for (Node node : iface->nodes()) {
        if (node.isOwn()) {
          m_linear_system.matrixSetValue(node_dof.dofId(node, 0), node_dof.dofId(node, 0), penalty);
          Real u_g = 0;
          if (mesh()->dimension() == 2)
            u_g = (m_node_coord[node].y - angle * m_node_coord[node].x) * penalty;
          if (mesh()->dimension() == 3)
            u_g = (m_node_coord[node].z - angle * m_node_coord[node].x) * penalty;
          rhs_values[node_dof.dofId(node, 0)] = u_g;
        }
      }
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "assemble-rhs", elapsedTime);
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

void FemModule::
_assembleLinearOperatorGpu()
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
  auto dim = mesh()->dimension();
  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  if (bc) {
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
      FemUtils::Gpu::BoundaryConditions::applyDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);
  }

  for (const auto& bs : options()->farfieldBoundaryCondition()) {
    Real angle = bs->angle();
    Real penalty = options()->penalty();
    FaceGroup face_group = bs->getSurface();
    NodeGroup node_group = face_group.nodeGroup();

    NodeInfoListView nodes_infos(mesh_ptr->nodeFamily());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    auto command = Accelerator::makeCommand(queue);
    auto in_m_node_coord = ax::viewIn(command, m_node_coord);
    auto in_out_forced_info = Accelerator::viewInOut(command, m_linear_system.getForcedInfo());
    auto in_out_forced_value = Accelerator::viewInOut(command, m_linear_system.getForcedValue());
    auto in_out_rhs_variable = Accelerator::viewInOut(command, m_linear_system.rhsVariable());

    command << RUNCOMMAND_ENUMERATE(NodeLocalId, node_lid, node_group)
    {
      if (nodes_infos.isOwn(node_lid)) {
        DoFLocalId dof_id = node_dof.dofId(node_lid, 0);
        in_out_forced_info[dof_id] = true;
        in_out_forced_value[dof_id] = penalty;
        if (dim == 2)
          in_out_rhs_variable[dof_id] = (in_m_node_coord[node_lid].y - angle * in_m_node_coord[node_lid].x) * penalty;
        else
          in_out_rhs_variable[dof_id] = (in_m_node_coord[node_lid].z - angle * in_m_node_coord[node_lid].x) * penalty;
      }
    };
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "rhs-vector-assembly-gpu", elapsedTime);
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
  else {
    if (mesh()->dimension() == 3)
      if (m_hex_quad_mesh)
        _assembleBilinear<8>([this](const Cell& cell) { return _computeElementMatrixHexa8(cell); });
      else
        _assembleBilinear<4>([this](const Cell& cell) { return _computeElementMatrixTetra4(cell); });
    if (mesh()->dimension() == 2)
      if (m_hex_quad_mesh)
        _assembleBilinear<4>([this](const Cell& cell) { return _computeElementMatrixQuad4(cell); });
      else
        _assembleBilinear<3>([this](const Cell& cell) { return _computeElementMatrixTria3(cell); });
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "lhs-matrix-assembly", elapsedTime);
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
        if (node1.isOwn())
          m_linear_system.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
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
  info() << "[ArcaneFem-Info] Started module _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.applyLinearSystemTransformationAndSolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "solve-linear-system", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Update the FEM variables.
 *
 * This method performs the following actions:
 *   1. Fetches values of solution from solved linear system to FEM variables,
 *      i.e., it copies RHS DOF to 𝒖.
 *   2. Performs synchronize of FEM variables across subdomains.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_updateVariables()
{
  info() << "[ArcaneFem-Info] Started module _updateVariables()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& dof_u(m_linear_system.solutionVariable());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Node, inode, ownNodes()) {
    Node node = *inode;
    m_u[node] = dof_u[node_dof.dofId(node, 0)];
  }

  m_u.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "update-variables", elapsedTime);
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

  checkNodeResultFile(traceMng(), filename, m_u, 1.0e-4);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
