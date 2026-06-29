// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2000-2026 */
/*                                                                           */
/* FEM code to test vectorial FE for Elastoplasticity problem.                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/IParallelMng.h>
// #include <arcane/VariableTypes.h>

#include "FemModule.h"
#include "ElementMatrix.h"
#include "ElementMatrixHexQuad.h"
#include "BodyForce.h"
#include "Traction.h"
#include "Dirichlet.h"
#include "NLResidualRHS.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModuleElastoplasticity at the start of the simulation.
 *
 * This method initializes degrees of freedom (DoFs) on nodes.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
startInit()
{
  info() << "[ArcaneFem-Info] Started module  startInit()";
  Real elapsedTime = platform::getRealTime();

  E = options()->E(); // Youngs modulus
  nu = options()->nu(); // Poission ratio ν

  m_dof_per_node = defaultMesh()->dimension();
  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->hasSolutionComparisonFile();
  m_petsc_flags = options()->petscFlags();
  m_hex_quad_mesh = options()->hexQuadMesh();

  m_dofs_on_nodes.initialize(defaultMesh(), m_dof_per_node);

  m_nonlinear_law = options()->nonlinearLaw();
  m_newton_max_iters = options()->newtonMaxIters();
  m_newton_atol = options()->newtonAtol();
  m_newton_rtol = options()->newtonRtol();

  m_gp_material_tensor_strategy = options()->gpMaterialTensorStrategy();

  if (m_gp_material_tensor_strategy == "global") {
    if (mesh()->dimension() == 2) {
      m_C_2d_cell.reshape({ 3, 3 });
      m_C_tang_2d_cell.reshape({ 3, 3 });
    } else {
      m_C_3d_cell.reshape({ 6, 6 });
      m_C_tang_3d_cell.reshape({ 6, 6 });
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs the main computation for the FemModuleElastoplasticity.
 *
 * This method:
 *   1. Stops the time loop after 1 iteration since the equation is steady state.
 *   2. Resets, configures, and initializes the linear system.
 *   3. Sets Petsc flags if user has provided them.
 *   4. Executes the stationary solve.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
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

void FemModuleElastoplasticity::_initBsr()
{
  info() << "[ArcaneFem-Info] Started module  _initBsr()";
  Real elapsedTime = platform::getRealTime();

  bool use_csr_in_linearsystem =
  options()->linearSystem.serviceName() == "HypreLinearSystem" ||
  options()->linearSystem.serviceName() == "AlienLinearSystem" ||
  options()->linearSystem.serviceName() == "PetscLinearSystem";

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
 * This method does the following
 *   1. Solves the FEM system either with:
 *      _solveNewton()
 *      or
 *      _solveLinear()
 *      based on the nature of the constitutive law
 *
 *   2. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_doStationarySolve()
{
  if (m_nonlinear_law) {
    _solveNewton();
  }
  else {
    _solveLinear();
  }

  if(m_cross_validation){
    _validateResults();
  }

}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a linear solve for the FEM system using Newton method.
 *
 * This method follows a sequence of steps to solve FEM system:
 *
 *   1. _getMaterialParameters()     Updates nonlinear material parameters
 *   2. _assembleBilinearOperator()  Assembles the FEM  matrix 𝐀
 *   3. _assembleLinearOperator()    Assembles the FEM RHS vector 𝐛
 *   4. _solve()                     Solves for solution vector 𝐮 = 𝐀⁻¹𝐛
 *   5. _updateVariables()           Updates FEM variables 𝐮 = 𝐱
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_solveLinear()
{
  info() << "[ArcaneFem-Info] Started module  _solveLinear()";
  _getMaterialParameters();

  if(m_assemble_linear_system){
    _assembleBilinearOperator();
    _assembleLinearOperator();
  }
  if(m_solve_linear_system){
    _solve();
    _updateVariables();
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a nonlinear solve for the FEM system using Newton method.
 *
 * This method follows a sequence of steps to solve FEM system:
 *
 *   1. _getMaterialParameters()     Updates nonlinear material parameters
 *   2. _assembleBilinearOperator()  Assembles the FEM  matrix 𝐀ʹ
 *   3. _assembleLinearOperator()    Assembles the FEM RHS vector 𝐛
 *   4. _solve()                     Solves for solution vector 𝐝𝐮 = 𝐀ʹ⁻¹𝐛
 *   5. _updateNewtonIncrements()    Updates FEM variables 𝐝𝐮 = 𝐱
 *   5. _checkNewtonConvergence()    Check convergence on norm of 𝐝𝐮 /𝐮
 *   6. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_solveNewton()
{
  info() << "[ArcaneFem-Info] Started module  _solveNewton()";

  while (m_newton_iter < m_newton_max_iters) {
    _getMaterialParameters();

    if(m_assemble_nonlinear_system) {
      if (m_linear_system.isInitialized() && m_newton_iter != 0) {
        m_linear_system.clearValues();

        if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR")
          m_bsr_format.resetMatrixValues();

        _assembleBilinearOperator(); // assembles Jacobian
        _assembleLinearOperator(); // assembles Residuals(m_U) + BCs
      } else {
        _assembleBilinearOperator(); // iter 0: initialisation contd.
        _assembleLinearOperator();   // iter 0: initialisation contd.
      }
    }

    if(m_solve_nonlinear_system){
      _solve();
      _updateNewtonIncrements();
    }

    ++m_newton_iter;
    _checkNewtonConvergence();

    // m_U.add(m_dU);
    _incrementVariables();

    if (m_newton_solver_converged) {
      info() << "[ArcaneFem-Info] Newton solver converged after " << m_newton_iter << " iterations.";
      m_newton_iter = 0;
      m_newton_solver_converged = false;
      break;
    } else {
      _updateGuessFromIncrement(); //TODO check if m_linear_solve keeps the previous solution
    }
  }

  if (m_newton_iter == m_newton_max_iters && !m_newton_solver_converged) {
    info() << "[ArcaneFem-Info] Newton iterations did not converge after maximum (" << m_newton_max_iters << ") iterations";
    ARCANE_FATAL("Newton iterations diverged after max iters");
  }

}

/*---------------------------------------------------------------------------*/
/**
 * @brief Retrieves and sets the material parameters for the simulation.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module  _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  mu = (E / (2 * (1 + nu))); // lame parameter μ
  lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter λ

  /*
   {{lambda + 2. * mu, lambda,           0.},
    {lambda,           lambda + 2. * mu, 0.},
    {0.,               0.,               mu}}
  */ // 2D elastic material tensor

  /*
   {{lambda + 2.*mu, lambda,         lambda,         0.,  0.,  0.},
    {lambda,         lambda + 2.*mu, lambda,         0.,  0.,  0.},
    {lambda,         lambda,         lambda + 2.*mu, 0.,  0.,  0.},
    {0.,             0.,             0.,             mu,  0.,  0.},
    {0.,             0.,             0.,             0.,  mu,  0.},
    {0.,             0.,             0.,             0.,  0.,  mu}}
  */ // 3D elastic material tensor

  if (m_gp_material_tensor_strategy == "local") {
     if (mesh()->dimension() == 2) {
       m_C_tang_2d.fill(0.);
       m_C_tang_2d(0, 0) = lambda + 2. * mu;
       m_C_tang_2d(1, 1) = lambda + 2. * mu;
       m_C_tang_2d(2, 2) = mu;
       m_C_tang_2d(0, 1) = lambda;
       m_C_tang_2d(1, 0) = lambda;

       m_C_2d = m_C_tang_2d; // elasticity
     } else {
       m_C_tang_3d.fill(0.);
       m_C_tang_3d(0, 0) = lambda + 2. * mu;
       m_C_tang_3d(1, 1) = lambda + 2. * mu;
       m_C_tang_3d(2, 2) = lambda + 2. * mu;
       m_C_tang_3d(3, 3) = mu;
       m_C_tang_3d(4, 4) = mu;
       m_C_tang_3d(5, 5) = mu;
       m_C_tang_3d(0, 1) = lambda;
       m_C_tang_3d(1, 0) = lambda;
       m_C_tang_3d(0, 2) = lambda;
       m_C_tang_3d(2, 0) = lambda;
       m_C_tang_3d(1, 2) = lambda;
       m_C_tang_3d(2, 1) = lambda;

       m_C_3d = m_C_tang_3d; // elasticity
     }
  } else {
    if (mesh()->dimension() == 2) {
      ENUMERATE_ (Cell, icell, allCells()) {
        for (Int32 ix = 0; ix < 3; ++ix) {
          for (Int32 iy = 0; iy < 3; ++iy) {
            m_C_tang_2d_cell(icell, ix, iy) = 0.0;
          }
        }
        m_C_tang_2d_cell(icell, 0, 0) = lambda + 2. * mu;
        m_C_tang_2d_cell(icell, 1, 1) = lambda + 2. * mu;
        m_C_tang_2d_cell(icell, 2, 2) = mu;
        m_C_tang_2d_cell(icell, 0, 1) = lambda;
        m_C_tang_2d_cell(icell, 1, 0) = lambda;

        for (Int32 ix = 0; ix < 3; ++ix) {
          for (Int32 iy = 0; iy < 3; ++iy) {
            m_C_2d_cell(icell, ix, iy) = m_C_tang_2d_cell(icell, ix, iy);
          }
        } // elasticity
      }
    } else {
      ENUMERATE_ (Cell, icell, allCells()) {
        for (Int32 ix = 0; ix < 6; ++ix) {
          for (Int32 iy = 0; iy < 6; ++iy) {
            m_C_tang_3d_cell(icell, ix, iy) = 0.0;
          }
        }
        m_C_tang_3d_cell(icell,0, 0) = lambda + 2. * mu;
        m_C_tang_3d_cell(icell,1, 1) = lambda + 2. * mu;
        m_C_tang_3d_cell(icell,2, 2) = lambda + 2. * mu;
        m_C_tang_3d_cell(icell,3, 3) = mu;
        m_C_tang_3d_cell(icell,4, 4) = mu;
        m_C_tang_3d_cell(icell,5, 5) = mu;
        m_C_tang_3d_cell(icell,0, 1) = lambda;
        m_C_tang_3d_cell(icell,1, 0) = lambda;
        m_C_tang_3d_cell(icell,0, 2) = lambda;
        m_C_tang_3d_cell(icell,2, 0) = lambda;
        m_C_tang_3d_cell(icell,1, 2) = lambda;
        m_C_tang_3d_cell(icell,2, 1) = lambda;

        for (Int32 ix = 0; ix < 6; ++ix) {
          for (Int32 iy = 0; iy < 6; ++iy) {
            m_C_3d_cell(icell, ix, iy) = m_C_tang_3d_cell(icell, ix, iy);
          }
        }
      }
    }
  }
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

void FemModuleElastoplasticity::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (m_nonlinear_law) {
    _applyResidualRHS(rhs_values, node_dof);
  }

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

void FemModuleElastoplasticity::
_assembleBilinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());
    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    RealMatrix<3, 3> C_tang_2d = m_C_tang_2d;
    RealMatrix<6, 6> C_tang_3d = m_C_tang_3d;

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2)
      m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord, C_tang_2d); });
    if (mesh()->dimension() == 3)
      m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord, C_tang_3d); });
    m_bsr_format.toLinearSystem(m_linear_system);
  }
  else if (m_matrix_format == "AF-BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());
    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    RealMatrix<3, 3> C_tang_2d = m_C_tang_2d;
    RealMatrix<6, 6> C_tang_3d = m_C_tang_3d;

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2)
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, C_tang_2d, node_lid); });
    if (mesh()->dimension() == 3)
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, C_tang_3d, node_lid); });
    m_bsr_format.toLinearSystem(m_linear_system);
  }
  else if (m_matrix_format == "DOK") {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {
        _assembleBilinearOperatorCpu<8>([this](const Cell& cell) { return _computeElementMatrixQuad4(cell); });
      }
      else {
        _assembleBilinearOperatorCpu<6>([&](const Cell& cell) { return _computeElementMatrixTria3(cell); });
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
void FemModuleElastoplasticity::
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

void FemModuleElastoplasticity::
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

void FemModuleElastoplasticity::
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
  const double min_value_to_test = 1.0e-10;

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

void FemModuleElastoplasticity::
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
/**
 * @brief Update the FEM Newton increment.
 *
 * This method performs the following actions:
 *   1. Fetches values of solution from solved linear system to FEM variables,
 *      i.e., it copies RHS DOF to du.
 *   2. Performs synchronize of FEM variables across subdomains.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_updateNewtonIncrements()
{
  info() << "[ArcaneFem-Info] Started module  _updateNewtonIncrements()";
  Real elapsedTime = platform::getRealTime();

  {
    VariableDoFReal& dof_du(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    if (mesh()->dimension() == 3)
      ENUMERATE_ (Node, inode, ownNodes()) {
        Node node = *inode;
        Real du1_val = dof_du[node_dof.dofId(node, 0)];
        Real du2_val = dof_du[node_dof.dofId(node, 1)];
        Real du3_val = dof_du[node_dof.dofId(node, 2)];
        m_dU[node] = Real3(du1_val, du2_val, du3_val);
      }
    else
      ENUMERATE_ (Node, inode, ownNodes()) {
        Node node = *inode;
        Real du1_val = dof_du[node_dof.dofId(node, 0)];
        Real du2_val = dof_du[node_dof.dofId(node, 1)];
        m_dU[node] = Real3(du1_val, du2_val, 0.);
      }
  }
  m_dU.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"update-Newton-increments", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Reinitialize the solution vector of the linear solve with the FEM variables.
 *
 * This method performs the following actions:
 *   1. Performs synchronization of FEM increment variables across subdomains.
 *   2. Fetches the FEM increment variables to the solution vector of the
 *      linear solver for the next nonlinear solver iteration.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_updateGuessFromIncrement()
{
  info() << "[ArcaneFem-Info] Started module _updateGuessFromIncrement()";
  Real elapsedTime = platform::getRealTime();

  m_dU.synchronize();

  {
    VariableDoFReal& dof_du(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    if (mesh()->dimension() == 3)
      ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      dof_du[node_dof.dofId(node, 0)] = m_dU[node][0];
      dof_du[node_dof.dofId(node, 1)] = m_dU[node][1];
      dof_du[node_dof.dofId(node, 2)] = m_dU[node][2];
    }
    else
      ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      dof_du[node_dof.dofId(node, 0)] = m_dU[node][0];
      dof_du[node_dof.dofId(node, 1)] = m_dU[node][1];
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "_update-guess-from-increment", elapsedTime);
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/**
 * @brief Increments the FEM variables.
 *
 * This method updates the FEM solutions with the increment of
 * the current Newton iteration
 *
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_incrementVariables()
{
  info() << "[ArcaneFem-Info] Started module _incrementVariables()";
  Real elapsedTime = platform::getRealTime();

  m_dU.synchronize();
  m_U.synchronize();
  {
      ENUMERATE_ (Node, inode, ownNodes()) {
      m_U[inode] += m_dU[inode];
    }
  }
  m_U.synchronize();

  IParallelMng* pm = defaultMesh()->parallelMng();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "increment-fem-variables", elapsedTime);
}
/*---------------------------------------------------------------------------*/
/**
 * @brief Check for the convergence of nonlinear solver.
 *
 * This method performs the following actions:
 *   1. Evaluates the convergence norm with Newton increment FEM variables
 *   2. Checks for convergence
 *
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_checkNewtonConvergence()
{
  info() << "[ArcaneFem-Info] Started module _checkNewtonConvergence()";
  Real elapsedTime = platform::getRealTime();

  m_dU.synchronize();
  m_U.synchronize();

  Real l2_norm_du = 0.0;
  Real l2_norm_u = 0.0;
  {
    ENUMERATE_ (Node, inode, ownNodes()) {
      const Real norm_du = math::pow(m_dU[inode][0], 2.0) + math::pow(m_dU[inode][1], 2.0) + math::pow(m_dU[inode][2], 2.0);
      const Real norm_u = math::pow(m_U[inode][0], 2.0) + math::pow(m_U[inode][1], 2.0) + math::pow(m_U[inode][2], 2.0);
      l2_norm_du += norm_du;
      l2_norm_u += norm_u;
    }
  }
  IParallelMng* pm = defaultMesh()->parallelMng();
  l2_norm_du = pm->reduce(Parallel::ReduceSum, l2_norm_du);
  l2_norm_u = pm->reduce(Parallel::ReduceSum, l2_norm_u);

  l2_norm_du = math::sqrt(l2_norm_du);
  l2_norm_u = math::sqrt(l2_norm_u);

  convergence_norm = l2_norm_du / (m_newton_rtol * l2_norm_u  + m_newton_atol);

  if ( convergence_norm <= 1.0){
    VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
    rhs_values.fill(0.0);
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    _applyResidualRHS(rhs_values, node_dof, true);
    info() << "[ArcaneFem-Info] At newton iteration "<< m_newton_iter <<": rel. conv. norm = " << convergence_norm << " and residual norm = " << rhs_norm;
    m_newton_solver_converged = true;
  } else {
    m_newton_solver_converged = false;
    info() << "[ArcaneFem-Info] At newton iteration "<< m_newton_iter <<": rel. conv. norm = " << convergence_norm;

  }


  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "check-newton-convergence", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModuleElastoplasticity);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
