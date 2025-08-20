// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Simple FEM solver to solve acoustics.                                     */
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
  m_cross_validation = options()->crossValidation();
  m_petsc_flags = options()->petscFlags();
  m_hex_quad_mesh = options()->hexQuadMesh();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "initialize", elapsedTime);
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

  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  if (m_petsc_flags != NULL) {
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows a sequence of steps to solve FEM system:
 *   1. _getMaterialParameters()          Retrieves material parameters via
 *   2. _assembleBilinearOperator()       Assembles the FEM  matrix 𝑨
 *   3. _assembleLinearOperator()         Assembles the FEM RHS vector 𝒃
 *   4.  _solve()                         Solves for solution vector 𝒖 = 𝑨⁻¹𝒃
 *   5. _updateVariables()                Updates FEM variables 𝒖 = 𝒙
 *   6. _validateResults()                Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  _getMaterialParameters();

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
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  m_kc2 = options()->kc2();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "get-material-parms", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the FEM linear operator.
 *
 * This method assembles the RHS 'b' of the linear system by:
 *   1. Initializing the RHS values to zero.
 *   2. Iterating over Neumann boundary conditions to apply it to RHS.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable for RHS values
  rhs_values.fill(0.0);

  const auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  BC::IArcaneFemBC* bc = options()->boundaryConditions();
  if (bc) {
    for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions())
      if (mesh()->dimension() == 2)
        if (m_hex_quad_mesh)
          ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhsQuad4(bs, node_dof, m_node_coord, rhs_values);
        else
          ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);
      else
        if (m_hex_quad_mesh)
          ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhsHexa8(bs, node_dof, m_node_coord, rhs_values);
        else
          ArcaneFemFunctions::BoundaryConditions3D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "rhs-vector-assembly", elapsedTime);
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
  const double epsilon = 1.0e-4;
  const double min_value_to_test = 1.0e-16;

  info() << "[ArcaneFem-Info] Validating results filename=" << filename << " epsilon =" << epsilon;

  if (!filename.empty())
    checkNodeResultFile(traceMng(), filename, m_u, epsilon, min_value_to_test);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
