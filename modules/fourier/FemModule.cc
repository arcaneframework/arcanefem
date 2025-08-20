// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Simple module to solve Fourier's equation using FEM.                      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include "ElementMatrix.h"
#include "ElementMatrixHexQuad.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModule at the start of the simulation.
 *
 *  - initializes degrees of freedom (DoFs) on nodes.
 *  - builds support for manufactured test case (optional).
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
  m_cross_validation = options()->crossValidation();
  m_petsc_flags = options()->petscFlags();
  m_hex_quad_mesh = options()->hexQuadMesh();

  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  for (BC::IManufacturedSolution* bs : bc->manufacturedSolutions()) {
    if (bs->getManufacturedSource()) {
      ICaseFunction* opt_function_source = bs->getManufacturedSourceFunction();
      IStandardFunction* scf_source = bs->getManufacturedSourceStandardFunction();
      if (!scf_source)
        ARCANE_FATAL("No standard case function for option 'manufactured-source-condition'");
      auto* functorS = scf_source->getFunctorRealReal3ToReal();
      if (!functorS)
        ARCANE_FATAL("Standard function '{0}' is not convertible to f(Real,Real3) -> Real", opt_function_source->name());
      m_manufactured_source = functorS;
    }

    if (bs->getManufacturedDirichlet()) {
      ICaseFunction* opt_function = bs->getManufacturedDirichletFunction();
      IStandardFunction* scf = bs->getManufacturedDirichletStandardFunction();
      if (!scf)
        ARCANE_FATAL("No standard case function for option 'manufactured-dirichlet-condition'");
      auto* functor = scf->getFunctorRealReal3ToReal();
      if (!functor)
        ARCANE_FATAL("Standard function '{0}' is not convertible to f(Real,Real3) -> Real", opt_function->name());
      m_manufactured_dirichlet = functor;
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs the main computation for the FemModule.
 *
 * - Stops the time loop after 1 iteration since the equation is steady state.
 * - Resets, configures, and initializes the linear system.
 * - Executes the stationary solve.
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
  m_linear_system.initialize(subDomain(),  acceleratorMng()->defaultRunner(),  m_dofs_on_nodes.dofFamily(), "Solver");

  if (m_petsc_flags != NULL){
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") {
    auto use_csr_in_linear_system = options()->linearSystem.serviceName() == "HypreLinearSystem";
    if (m_matrix_format == "BSR")
      m_bsr_format.initialize(mesh(), 1, use_csr_in_linear_system, 0);
    else
      m_bsr_format.initialize(mesh(), 1, use_csr_in_linear_system, 1);
    m_bsr_format.computeSparsity();
  }

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows via the following steps:
 *   1. _getMaterialParameters()     Retrieves material parameters
 *   2. _assembleBilinearOperator()  Assembles the FEM  matrix 𝐀
 *   3. _assembleLinearOperator()    Assembles the FEM RHS vector 𝐛
 *   4. _solve()                     Solves for solution vector 𝐮 = 𝐀⁻¹𝐛
 *   5. _updateVariables()           Updates FEM variables 𝐮 = 𝐱
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
 *
 * This method initializes:
 *  - material properties:
 *       # thermal conductivity coefficient (`lambda`)
 *       # heat source term (`qdot`)
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  lambda = options()->lambda();
  qdot = options()->qdot();

  m_cell_lambda.fill(lambda);

  for (const auto& bs : options()->materialProperty()) {
    CellGroup group = bs->volume();
    Real value = bs->lambda();
    info() << "Lambda for group= " << group.name() << " v=" << value;

    ENUMERATE_ (Cell, icell, group) {
      Cell cell = *icell;
      m_cell_lambda[cell] = value;
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "get-material-params", elapsedTime);
}

void FemModule::
_assembleLinearOperator()
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
 *  2. If a constant source term is specified (`qdot`), apply it to the RHS.
 *  3. If Neumann BC are specified applied to the RHS.
 *  4. If Dirichlet BC are specified apply to the LHS & RHS. 
 *  5. If manufactured source conditions is specified, apply to RHS.
 *  6. If  manufactured Dirichlet BC are specified apply to the LHS & RHS. 
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

  if (options()->qdot.isPresent()) {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh)
        ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhsQuad4(qdot, mesh(), node_dof, m_node_coord, rhs_values);
      else
        ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhs(qdot, mesh(), node_dof, m_node_coord, rhs_values);
    }
    else {
      if (m_hex_quad_mesh)
        ArcaneFemFunctions::BoundaryConditions3D::applyConstantSourceToRhsHexa8(qdot, mesh(), node_dof, m_node_coord, rhs_values);
      else
        ArcaneFemFunctions::BoundaryConditions3D::applyConstantSourceToRhs(qdot, mesh(), node_dof, m_node_coord, rhs_values);
    }
  }

  // Helper lambda to apply boundary conditions
  auto applyBoundaryConditions = [&](auto BCFunctions) {
    BC::IArcaneFemBC* bc = options()->boundaryConditions();
    if (bc) {

      // Neumann
      for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions())
        if (mesh()->dimension() == 2)
          if (m_hex_quad_mesh)
            ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhsQuad4(bs, node_dof, m_node_coord, rhs_values);
          else
            ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);
        else if (m_hex_quad_mesh)
          ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhsHexa8(bs, node_dof, m_node_coord, rhs_values);
        else
          ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);

      // Dirichlet
      for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
        BCFunctions.applyDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);

      // Manufactured boundary conditions
      for (BC::IManufacturedSolution* bs : bc->manufacturedSolutions()) {
        if (bs->getManufacturedSource()) {
          ARCANE_CHECK_POINTER(m_manufactured_source);
          info() << "Apply manufactured Source condition to all cells";
          BCFunctions.applyManufacturedSourceToRhs(m_manufactured_source, mesh(), node_dof, m_node_coord, rhs_values);
        }
        if (bs->getManufacturedDirichlet()) {
          ARCANE_CHECK_POINTER(m_manufactured_dirichlet);
          info() << "Apply manufactured dirichlet condition to all borders";
          FaceGroup group = mesh()->outerFaces();
          BCFunctions.applyManufacturedDirichletToLhsAndRhs(m_manufactured_dirichlet, lambda, group, bs, node_dof, m_node_coord, m_linear_system, rhs_values);
        }
      }
    }
  };

  if (mesh()->dimension() == 3)
    applyBoundaryConditions(ArcaneFemFunctions::BoundaryConditions3D());
  else
    applyBoundaryConditions(ArcaneFemFunctions::BoundaryConditions2D());

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "rhs-vector-assembly", elapsedTime);
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
    if (options()->qdot.isPresent())
      BCFunctions.applyConstantSourceToRhs(qdot, m_dofs_on_nodes, m_node_coord, rhs_values, mesh_ptr, queue);

    BC::IArcaneFemBC* bc = options()->boundaryConditions();
    if (bc) {
      for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions())
        BCFunctions.applyNeumannToRhs(bs, m_dofs_on_nodes, m_node_coord, rhs_values, mesh_ptr, queue);

      for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
        FemUtils::Gpu::BoundaryConditions::applyDirichletViaPenalty(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);

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
 * @brief Calls the right function for LHS assembly given as mesh type.
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
    auto in_cell_lambda = ax::viewIn(command, m_cell_lambda);

    if (mesh()->dimension() == 2)
      if (m_matrix_format == "BSR")
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord, in_cell_lambda); });
      else
        m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, in_cell_lambda, node_lid); });
    else
      if (m_matrix_format == "BSR")
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_cell_lambda); });
      else
        m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_cell_lambda, node_lid); });
  }

  if (m_matrix_format == "DOK") {
    if (mesh()->dimension() == 3)
      if(m_hex_quad_mesh)
        _assembleBilinear<8>([this](const Cell& cell) { return _computeElementMatrixHexa8(cell); });
      else
        _assembleBilinear<4>([this](const Cell& cell) { return _computeElementMatrixTetra4(cell); });

    if (mesh()->dimension() == 2)
      if(m_hex_quad_mesh)
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
 *   1. For each cell, retrieves the cell-specific constant `lambda`.
 *   2. Computes element matrix using provided `compute_element_matrix` function.
 *   3. Assembles global matrix by adding contributions from each cell's element 
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

    lambda = m_cell_lambda[cell]; // lambda is always considered cell constant
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
 *
 * This method performs the following actions:
 *   1. Solves the linear system to compute the solution.
 *   2. Copies the computed solution from the DoF to the node values.
 *   3. Synchronizes the updated node values.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  info() << "[ArcaneFem-Info] Started module _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.solve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "solve-linear-system", elapsedTime);
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

  { // copies solution (and optionally exact solution) to FEM output
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_u[node_dof.dofId(node, 0)];
      m_u[node] = v;
    }

    for (BC::IManufacturedSolution* bs : options()->boundaryConditions()->manufacturedSolutions())
      ENUMERATE_ (Node, inode, ownNodes()) {
        Node node = *inode;
        m_u_exact[node] = m_manufactured_dirichlet->apply(lambda, m_node_coord[node]);
      }
  }

  m_u.synchronize();

  for (BC::IManufacturedSolution* bs : options()->boundaryConditions()->manufacturedSolutions())
    m_u_exact.synchronize();

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

  String filename = options()->resultFile();

  if (!filename.empty())
    checkNodeResultFile(traceMng(), filename, m_u, 1.0e-4);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
