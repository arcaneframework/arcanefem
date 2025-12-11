// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* FEM code to test vectorial FE for Elasticity problem.                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/IParallelMng.h>

#include "FemModule.h"
#include "ElementMatrix.h"
#include "ElementMatrixHexQuad.h"
#include "BodyForce.h"
#include "Traction.h"
#include "Dirichlet.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModuleElasticity at the start of the simulation.
 *
 * This method initializes degrees of freedom (DoFs) on nodes.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElasticity::
startInit()
{
  info() << "[ArcaneFem-Info] Started module  startInit()";
  Real elapsedTime = platform::getRealTime();

  _getMaterialParameters();

  m_dofs_on_nodes.initialize(defaultMesh(), m_dof_per_node);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs the main computation for the FemModuleElasticity.
 *
 * This method:
 *   1. Stops the time loop after 1 iteration since the equation is steady state.
 *   2. Resets, configures, and initializes the linear system.
 *   3. Sets PETSc flags if user has provided them.
 *   4. Executes the stationary solve.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElasticity::
compute()
{
  info() << "[ArcaneFem-Info] Started module  compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");
  m_linear_system.clearValues();

  if (m_petsc_flags != NULL){
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR")
    _initBsr();

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
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes BSR matrix.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElasticity::_initBsr()
{
  info() << "[ArcaneFem-Info] Started module  _initBsr()";
  Real elapsedTime = platform::getRealTime();

  bool use_csr_in_linearsystem = options()->linearSystem.serviceName() == "HypreLinearSystem" ||
    options()->linearSystem.serviceName() == "PETScLinearSystem";

  if (m_matrix_format == "BSR")
    m_bsr_format.initialize(defaultMesh(), m_dof_per_node, use_csr_in_linearsystem, 0);
  else
    m_bsr_format.initialize(defaultMesh(), m_dof_per_node, use_csr_in_linearsystem, 1);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"initialize-bsr-matrix", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows a sequence of steps to solve FEM system:
 *
 *   1. _getMaterialParameters()     Retrieves material parameters via
 *   2. _assembleBilinearOperator()  Assembles the FEM  matrix 𝐀
 *   3. _assembleLinearOperator()    Assembles the FEM RHS vector 𝐛
 *   4. _solve()                     Solves for solution vector 𝐮 = 𝐀⁻¹𝐛
 *   5. _updateVariables()           Updates FEM variables 𝐮 = 𝐱
 *   6. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModuleElasticity::
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
 * @brief Retrieves and sets the material parameters for the simulation.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElasticity::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module  _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  E = options()->E(); // Youngs modulus
  nu = options()->nu(); // Poission ratio ν

  mu = (E / (2 * (1 + nu))); // lame parameter μ
  lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter λ

  m_dof_per_node = defaultMesh()->dimension();
  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->hasSolutionComparisonFile();
  m_petsc_flags = options()->petscFlags();
  m_hex_quad_mesh = options()->hexQuadMesh();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/**
 * @brief Assemble the FEM linear operator.
 *
 * This method follows a sequence of steps to assemble RHS of FEM linear system:
 *
 *   1. assembles the bodyforce contribution (source term) ∫∫∫ (𝐟.𝐯) on Ω
 *   2. assembles the traction contribution (Neumann term) ∫∫ (𝐭.𝐯)  on ∂Ω
 *   3. apply Dirichlet contributions to LHS and RHS
 */
/*---------------------------------------------------------------------------*/

void FemModuleElasticity::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  _applyBodyForce(rhs_values, node_dof);
  _applyTraction(rhs_values, node_dof);
  _applyDirichlet(rhs_values, node_dof);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly given as mesh type.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElasticity::
_assembleBilinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());
    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto lambda_copy = lambda;
    auto mu_copy = mu;

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2)
      m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord, lambda_copy, mu_copy); });
    if (mesh()->dimension() == 3)
      m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord, lambda_copy, mu_copy); });
    m_bsr_format.toLinearSystem(m_linear_system);
  }
  else if (m_matrix_format == "AF-BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());
    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto lambda_copy = lambda;
    auto mu_copy = mu;

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2)
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, lambda_copy, mu_copy, node_lid); });
    if (mesh()->dimension() == 3)
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, lambda_copy, mu_copy, node_lid); });
    m_bsr_format.toLinearSystem(m_linear_system);
  }
  else if (m_matrix_format == "DOK") {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {
        _assembleBilinearOperatorCpu<8>([this](const Cell& cell) { return _computeElementMatrixQuad4(cell); });
      }
      else {
        _assembleBilinearOperatorCpu<6>([this](const Cell& cell) { return _computeElementMatrixTria3(cell); });
      }
    }
    if (mesh()->dimension() == 3) {
      if (m_hex_quad_mesh) {
        _assembleBilinearOperatorCpu<24>([this](const Cell& cell) { return _computeElementMatrixHexa8(cell); });
      }
      else {
        _assembleBilinearOperatorCpu<12>([this](const Cell& cell) { return _computeElementMatrixTetra4(cell); });
      }
    }
  }
  else {
    ARCANE_FATAL("Unsupported matrix type, only DOK| BSR|AF-BSR is supported.");
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"lhs-matrix-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the FEM bilinear operator on CPU.
 *
 * This method assembles the FEM stiffness matrix by iterating over each cell,
 * computing the element stiffness matrix using the provided function, and
 * populating the global stiffness matrix accordingly.
 *
 * @tparam N Total DOF size (nodes_per_element × dimensions).
 * @param compute_element_matrix function computing cell's element stiffness matrix.
 */
/*---------------------------------------------------------------------------*/

template <int N>
void FemModuleElasticity::
_assembleBilinearOperatorCpu(const std::function<RealMatrix<N, N>(const Cell&)>& compute_element_matrix)
{
  const Int32 dim = mesh()->dimension();
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    auto K_e = compute_element_matrix(cell);

    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      if (node1.isOwn()) {
        Int32 n2_index = 0;
        for (Node node2 : cell.nodes()) {
          for (Int32 i = 0; i < dim; ++i) {
            DoFLocalId dof1 = node_dof.dofId(node1, i);
            for (Int32 j = 0; j < dim; ++j) {
              DoFLocalId dof2 = node_dof.dofId(node2, j);
              Real value = K_e(dim * n1_index + i, dim * n2_index + j);
              m_linear_system.matrixAddValue(dof1, dof2, value);
            }
          }
          ++n2_index;
        }
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

void FemModuleElasticity::
_solve()
{
  info() << "[ArcaneFem-Info] Started module  _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.applyLinearSystemTransformationAndSolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"solve-linear-system", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModuleElasticity::
_validateResults()
{
  info() << "[ArcaneFem-Info] Started module  _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200) {
    int p = std::cout.precision();
    std::cout.precision(17);
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      std::cout << "( N_id, u1, u2, u3 ) = ( "
                << node.uniqueId() << ", " << m_U[node].x << ", " << m_U[node].y << ", " << m_U[node].z
                << ")\n";
    }
    std::cout.precision(p);
  }

  String filename = options()->solutionComparisonFile();
  const double epsilon = options()->resultEpsilon();
  const double min_value_to_test = 1.0e-16;

  Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_U, epsilon, min_value_to_test);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"result-validation", elapsedTime);
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

void FemModuleElasticity::
_updateVariables()
{
  info() << "[ArcaneFem-Info] Started module  _updateVariables()";
  Real elapsedTime = platform::getRealTime();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    if (mesh()->dimension() == 3)
      ENUMERATE_ (Node, inode, ownNodes()) {
        Node node = *inode;
        Real u1_val = dof_u[node_dof.dofId(node, 0)];
        Real u2_val = dof_u[node_dof.dofId(node, 1)];
        Real u3_val = dof_u[node_dof.dofId(node, 2)];
        m_U[node] = Real3(u1_val, u2_val, u3_val);
      }
    else
      ENUMERATE_ (Node, inode, ownNodes()) {
        Node node = *inode;
        Real u1_val = dof_u[node_dof.dofId(node, 0)];
        Real u2_val = dof_u[node_dof.dofId(node, 1)];
        m_U[node] = Real3(u1_val, u2_val, 0.);
      }
  }

  m_U.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"update-variables", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModuleElasticity);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
