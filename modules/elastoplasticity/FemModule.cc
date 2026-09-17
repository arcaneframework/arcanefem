// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2000-2026 */
/*                                                                           */
/* FEM code for Elastoplasticity problem.                                    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/IParallelMng.h>
#include <arcane/accelerator/MDVariableViews.h>
#include <arcane/utils/ValueConvert.h>

#include "FemModule.h"
#include "ElementMatrix.h"
#include "ElementMatrixHexQuad.h"
#include "ExternalBodyForce.h"
#include "Traction.h"
#include "Dirichlet.h"
#include "InternalBodyForce.h"
#include "VonMisesLaw.h"
#include "DruckerPragerLaw.h"

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

  tmax = options()->tmax(); // max time 𝑡ₘₐₓ
  dt = options()->dt(); // time step δ𝑡

  m_dof_per_node = defaultMesh()->dimension();
  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_solve_nonlinear_system = options()->solveNonlinearSystem();
  m_cross_validation = options()->hasSolutionComparisonFile();
  m_petsc_flags = options()->petscFlags();
  m_hex_quad_mesh = options()->hexQuadMesh();

  m_dofs_on_nodes.initialize(defaultMesh(), m_dof_per_node);

  m_newton_max_iters = options()->newtonMaxIters();
  m_newton_atol = options()->newtonAtol();
  m_newton_rtol = options()->newtonRtol();

  m_use_gpu_functions =   (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") &&
                          (options()->linearSystem.serviceName() == "HypreLinearSystem" ||
                           options()->linearSystem.serviceName() == "PetscLinearSystem" ||
                           options()->linearSystem.serviceName() == "AlephLinearSystem");

  m_gp_material_tensor_strategy = options()->gpMaterialTensorStrategy();
  m_check_with_bilinear_operator = options()->checkBilinearOperatorForResidual();

  if (mesh()->dimension() == 2) {
    if (m_hex_quad_mesh)
      m_nGP = 4;
    else
      m_nGP = 1;

    if (m_gp_material_tensor_strategy == "global")
      m_C_tang_gp.reshape({m_nGP, 3, 3 });

    m_sigma_gp.reshape({m_nGP, 3});
    m_sigma_old_gp.reshape({m_nGP, 3}); // TODO move to von mises only
    m_sigma_zz_gp.reshape({m_nGP});
    m_sigma_zz_old_gp.reshape({m_nGP}); // TODO move to von mises only
  } else {
    if (m_hex_quad_mesh)
      m_nGP = 8;
    else
      m_nGP = 1;

    if (m_gp_material_tensor_strategy == "global")
      m_C_tang_gp.reshape({ m_nGP, 6, 6 });

    m_sigma_gp.reshape({m_nGP, 6});
    m_sigma_old_gp.reshape({m_nGP, 6});
  }

  t = dt;
  tmax = tmax - dt;
  m_global_deltat.assign(dt);

  _initConstitutiveLaw();

  // The Drucker-Prager return mapping is currently implemented on the CPU.
  // BSR may still assemble the matrix on an accelerator, but constitutive,
  // internal-force, and Dirichlet updates must use their CPU implementations.
  if (m_constitutive_law == "DruckerPrager")
    m_use_gpu_functions = false;

  _readCaseTables();

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

  // Stop the computation loop if the maximum time is reached
  if (t >= tmax)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  info() << "[ArcaneFem-Info] Time iteration at t : " << t << " (s) ";
  bool keep_struct = true;
  if (m_linear_system.isInitialized() && keep_struct) {
    m_linear_system.clearValues();
  } else {
    m_linear_system.reset();
    m_linear_system.setLinearSystemFactory(options()->linearSystem());
    m_linear_system.initialize(subDomain(), acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");
  }

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
  _updateTime();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"compute", elapsedTime);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Updates time.
 */
/*---------------------------------------------------------------------------*/
void FemModuleElastoplasticity::
_updateTime()
{
  t += dt;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes constitutive law parameters.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_initConstitutiveLaw()
{
  info() << "[ArcaneFem-Info] Started module  _initConstitutiveLaw()";
  Real elapsedTime = platform::getRealTime();

  for (const auto& constitutive_law : options()->constitutiveLaw()) {
    String law_name = constitutive_law->law();
    m_constitutive_law = law_name;
    if (law_name == "VonMises") {
      for (const auto von_mises : constitutive_law->vonMises()) {
        E = von_mises->E(); // Youngs modulus
        nu = von_mises->nu(); // Poission ratio ν
        sig0 = von_mises->sig0(); // Yield Strength
      }
    } else if (law_name == "DruckerPrager") {
      for (const auto drucker_prager : constitutive_law->druckerPrager()) {
        E = drucker_prager->E(); // Youngs modulus
        nu = drucker_prager->nu(); // Poission ratio ν
        cohesion = drucker_prager->cohesion(); // Cohesion
        friction_angle = drucker_prager->frictionAngle(); // Friction angle
      }
    } else {
      ARCANE_FATAL("Undefined constitutive law");
    }
  }

  if (m_constitutive_law == "VonMises" || m_constitutive_law == "DruckerPrager") {
    if (mesh()->dimension() != 2 || m_hex_quad_mesh)
      ARCANE_FATAL("Native von Mises plasticity currently supports only 2D Tria3 elements");

    if (m_constitutive_law == "VonMises") {
      m_p_old_gp.reshape({m_nGP});
      m_dp_gp.reshape({m_nGP});
    }

    if (m_constitutive_law == "DruckerPrager") {
      m_eps_p_gp.reshape({m_nGP, 3});
      m_eps_p_old_gp.reshape({m_nGP, 3});
      m_eps_p_zz_gp.reshape({m_nGP});
      m_eps_p_zz_old_gp.reshape({m_nGP});
    }

  }
  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"initialize-constitutive-law", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes BSR matrix.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_initBsr()
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
 *   1. Solves the nonlinear FEM system with _solveNewton()
 *   2. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_doStationarySolve()
{
    if (m_solve_nonlinear_system)
      _solveNewton();

   if(m_cross_validation)
     if (t > 0. && t==tmax)
       _validateResults();

}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a nonlinear solve for the FEM system using Newton method.
 *
 * This method follows a sequence of steps to solve FEM system:
 *
 *   1. _getMaterialParameters()     Updates nonlinear material parameters
 *   2. _assembleBilinearOperatorGlobal()  Assembles the FEM  matrix 𝐀ʹ
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

  _getMaterialParameters();

  // --- initialize_increment ---- //
  m_DUn.fill({0., 0., 0.});
  m_DUk.fill({0., 0., 0.});
  m_newton_iter = 0;

  if (m_constitutive_law == "VonMises") {
    _restoreConvergedStateVonMises();
  } else if (m_constitutive_law == "DruckerPrager") {
    _restoreConvergedStateDruckerPrager();
  }

  if (m_gp_material_tensor_strategy == "global")
    _setGlobalElasticMaterialTensorAtGPs();

  // --- assemble_linear_system ---- //
  if (m_assemble_linear_system) {
    if (m_gp_material_tensor_strategy == "global")
      _assembleBilinearOperatorGlobal();
    else
      if (m_constitutive_law == "VonMises")
        _assembleBilinearOperatorLocalVonMises(true); //checks law an assembles corresponding matrix
      else if (m_constitutive_law == "DruckerPrager")
        _assembleBilinearOperatorLocalDruckerPrager(true);
      else
        ARCANE_FATAL("Constitutive law not supported");

    _assembleLinearOperator();
  }

  // --- calculate_residual ---- //
  VariableDoFReal& residual_values(m_linear_system.rhsVariable());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  m_residual_norm0 = _normL2(residual_values, node_dof);
  info() << "[ArcaneFem-Info] Initial residual norm = " << m_residual_norm0;

  // --- start_newton_loop ---- //
  while (m_newton_iter < m_newton_max_iters && !m_newton_solver_converged) {
    m_newton_iter++;

    // --- solve_linear_system ---- //
    if(m_solve_linear_system){
      _solve();
      _updateNewtonIncrements();
    }

    // --- update_increment ---- //
    _incrementVariables();

    if (m_gp_material_tensor_strategy == "global") {
      if (m_constitutive_law == "VonMises") {
        _updateGlobalTangentMaterialTensorVonMises();
      } else if (m_constitutive_law == "DruckerPrager") {
        _updateGlobalTangentMaterialTensorDruckerPrager();
      }
    }

    // --- assemble_linear_system ---- //
    if (m_linear_system.isInitialized()) {
      m_linear_system.clearValues();

      if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR")
        m_bsr_format.resetMatrixValues();
      if (m_gp_material_tensor_strategy == "global")
        _assembleBilinearOperatorGlobal(); // assembles Jacobian
      else
        if (m_constitutive_law == "VonMises")
          _assembleBilinearOperatorLocalVonMises(); // assembles Jacobian for Von Mises
        else if (m_constitutive_law == "DruckerPrager")
          _assembleBilinearOperatorLocalDruckerPrager(); // assembles Jacobian for Drucker Prager
        else
          ARCANE_FATAL("Constitutive law not supported");

      _assembleLinearOperator(); // assembles Residuals(m_DUn) + BCs
    }

    if (m_newton_iter == 1) {
      if (m_constitutive_law == "DruckerPrager") {
        // The first assembled rhs contains the large algebraic enforcement of the
        // non-zero footing displacement. Reset the Newton reference norm after that
        // correction so convergence is measured using the physical equilibrium
        // residual, not the Dirichlet penalty.
        _applyZeroRHSOnConstrainedDOFs(residual_values, node_dof);
        m_residual_norm0 = _normL2(residual_values, node_dof);
      }
    }

    // --- calculate_and_check_residual ---- //
    _checkNewtonConvergence();

  }

  if (m_newton_solver_converged) {
    info() << "[ArcaneFem-Info] Newton solver converged after " << m_newton_iter << " iterations.";
    //-- update global displacement after newton convergence -- //
    _updateTimeVariables();

    if (m_constitutive_law == "VonMises") {
      //-- commit increment for von mises -- //
      _commitInternalVariablesVonMises();

      if (t == dt) {
        Real Ri = 1.0;
        Real Re = 1.3;
        Qlim = 2./math::sqrt(3.) * math::log( Re/Ri) * sig0;
      }
      Real tl = math::sqrt(1.1 / tmax * (t));
      info() << "[ArcaneFem-Info] At Time Step "
             << t - 1 << ":\tPressure applied: " << Qlim * tl
             << "\tNewton iters: " << m_newton_iter
             << "\tResidual norm: " << m_residual_norm;
    }

    if (m_constitutive_law == "DruckerPrager") {
      // -- commit increment for Drucker Prager --//
      _commitInternalVariablesDruckerPrager();

      if (t == dt) {
        max_settlement = 0.03;
        footing_width = 1.0;
      }

      VariableDoFReal& algebraic_reaction(m_linear_system.rhsVariable());
      algebraic_reaction.fill(0.);

      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;

        Int8 iGP = 0; // for tria P1 elements nGP=1
        Real sigma_xx = m_sigma_gp(cell , iGP, 0);
        Real sigma_yy = m_sigma_gp(cell , iGP, 1);
        Real sigma_xy = m_sigma_gp(cell , iGP, 2);

        RealVector<6> u = {0., 0., 0., 0., 0., 0.};
        for (Int8 i = 0; i < 3; ++i) {
          Real3 vertex = m_node_coord[cell.nodeId(i)];
          u[2*i + 1] = (math::abs(vertex.y - 10.) < 1.e-8 && vertex.x <= footing_width + 1.e-8);
        }

        Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
        Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
        Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

        RealVector<6> epsxx = { dxu[0] * u[0], 0., dxu[1] * u[2], 0., dxu[2] * u[4], 0. };
        RealVector<6> epsyy = { 0., dyu[0] * u[1], 0., dyu[1] * u[3], 0., dyu[2] * u[5] };
        RealVector<6> epsxy = { dyu[0] * u[0], dxu[0] * u[1], dyu[1] * u[2], dxu[1] * u[3], dyu[2] * u[4], dxu[2] * u[5] };
        epsxy = 0.70710678118654746172 * epsxy;

        RealVector<6> rhs = area * (sigma_xx * epsxx + sigma_yy * epsyy + sigma_xy * epsxy);

        algebraic_reaction[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
        algebraic_reaction[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
        algebraic_reaction[node_dof.dofId(cell.nodeId(1), 0)] += rhs(2);
        algebraic_reaction[node_dof.dofId(cell.nodeId(1), 1)] += rhs(3);
        algebraic_reaction[node_dof.dofId(cell.nodeId(2), 0)] += rhs(4);
        algebraic_reaction[node_dof.dofId(cell.nodeId(2), 1)] += rhs(5);
      }

      alg_reaction = _normL1(algebraic_reaction, node_dof);

      Real normalized_pressure = - alg_reaction / (footing_width * cohesion);
      Real settlement = t / tmax * max_settlement;

      info() << "[ArcaneFem-Info] At Time Step "
             << t - 1 << ":\tSettlement: " << settlement
             << "\tNormalised pressure: " << normalized_pressure
             << "\tNewton iters: " << m_newton_iter
             << "\tResidual norm: " << m_residual_norm;
    }

    m_newton_solver_converged = false;
    m_newton_iter = 0;
  }

  if (m_newton_iter == m_newton_max_iters && !m_newton_solver_converged) {
    info() << "[ArcaneFem-Info] Newton iterations did not converge after maximum (" << m_newton_max_iters << ") iterations";
    ARCANE_FATAL("Newton iterations diverged after max iters");
  }
}


/*---------------------------------------------------------------------------*/
/**
 * @brief Sets the material tensor at Gauss points to Elastic Law.
 */
/*---------------------------------------------------------------------------*/
void FemModuleElastoplasticity::
_setGlobalElasticMaterialTensorAtGPs()
{
  Int8 nDim = (mesh()->dimension() == 2) ? 3 : 6;
  ENUMERATE_ (Cell, icell, allCells())
  {
    Cell cell = *icell;
    for (Int8 iGP = 0; iGP < m_nGP; ++iGP)
      for (Int8 ix = 0; ix < nDim; ++ix)
        for (Int8 iy = 0; iy < nDim; ++iy)
          m_C_tang_gp(cell, iGP, ix, iy) = m_C_elas_2d(ix, iy); // set tangent C equal to elastic C
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

  if (m_material_initialized)
    return;

  if (m_constitutive_law == "VonMises") {
    mu = (E / (2 * (1 + nu))); // lame parameter μ
    lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter λ

    Et = E / 100.;
    H = E * Et / (E - Et);

    ENUMERATE_ (Cell, icell, allCells())
    {
      for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
        m_p_old_gp(icell, iGP) = 0.;
        m_dp_gp(icell, iGP) = 0.;
      }
    }

  } else if (m_constitutive_law == "DruckerPrager") {

    mu = (E / (2 * (1 + nu))); // lame parameter μ
    lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter λ

    bulk   = E/(3.*(1.-2.*nu));
    dpEta = 3. * tan(friction_angle) / (math::sqrt(9. + 12. * tan(friction_angle) * tan(friction_angle)));
    dpC = 3. * cohesion / (math::sqrt(9. + 12. * tan(friction_angle) * tan(friction_angle)));

    ENUMERATE_ (Cell, icell, allCells())
    {
      for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
        m_eps_p_gp(icell, iGP, 0) = 0.;
        m_eps_p_gp(icell, iGP, 1) = 0.;
        m_eps_p_gp(icell, iGP, 2) = 0.;
        m_eps_p_zz_gp(icell, iGP) = 0.;
        m_eps_p_old_gp(icell, iGP, 0) = 0.;
        m_eps_p_old_gp(icell, iGP, 1) = 0.;
        m_eps_p_old_gp(icell, iGP, 2) = 0.;
        m_eps_p_zz_old_gp(icell, iGP) = 0.;
      }
    }
  }

  if (m_constitutive_law == "VonMises" || m_constitutive_law == "DruckerPrager") {
    if (mesh()->dimension() == 2) {

      // Initialize elastic part of the material tensor
      /*
        lambda + 2mu  lambda         0
        lambda        lambda + 2mu   0
          0             0           2mu
      */
      m_C_elas_2d.fill(0.);
      m_C_elas_2d(0, 0) = lambda + 2. * mu;
      m_C_elas_2d(1, 1) = lambda + 2. * mu;
      m_C_elas_2d(2, 2) = 2. * mu;
      m_C_elas_2d(0, 1) = lambda;
      m_C_elas_2d(1, 0) = lambda;

      // Initialize constitutive history
      ENUMERATE_ (Cell, icell, allCells()) // TODO check if MDMeshVars provide initialisation method
      {
        for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
          m_sigma_gp(icell, iGP, 0) = 0.;
          m_sigma_gp(icell, iGP, 1) = 0.;
          m_sigma_gp(icell, iGP, 2) = 0.;
          m_sigma_zz_gp(icell, iGP) = 0.;

          m_sigma_old_gp(icell, iGP, 0) = 0.;
          m_sigma_old_gp(icell, iGP, 1) = 0.;
          m_sigma_old_gp(icell, iGP, 2) = 0.;
          m_sigma_zz_old_gp(icell, iGP) = 0.;
        }
      }

      // Initialize the tangent material tensor
      if (m_gp_material_tensor_strategy == "global") {
        _setGlobalElasticMaterialTensorAtGPs();
      }
    } else {
      ARCANE_FATAL("Not implemented for 3D yet");

      // Initialize elastic part of the material tensor
        /*
          lambda+2mu  lambda     lambda     0    0    0
          lambda      lambda+2mu lambda     0    0    0
          lambda      lambda     lambda+2mu 0    0    0
          0           0          0          2mu  0    0
          0           0          0          0    2mu  0
          0           0          0          0    0    2mu
      */
      m_C_elas_3d.fill(0.);
      m_C_elas_3d(0, 0) = lambda + 2. * mu;
      m_C_elas_3d(1, 1) = lambda + 2. * mu;
      m_C_elas_3d(2, 2) = lambda + 2. * mu;
      m_C_elas_3d(3, 3) = 2 * mu;
      m_C_elas_3d(4, 4) = 2 * mu;
      m_C_elas_3d(5, 5) = 2 * mu;
      m_C_elas_3d(0, 1) = lambda;
      m_C_elas_3d(1, 0) = lambda;
      m_C_elas_3d(0, 2) = lambda;
      m_C_elas_3d(2, 0) = lambda;
      m_C_elas_3d(1, 2) = lambda;
      m_C_elas_3d(2, 1) = lambda;

      // Initialize the tangent material tensor
      if (m_gp_material_tensor_strategy == "global") {
        _setGlobalElasticMaterialTensorAtGPs();
      }
    }
  }
  m_material_initialized = true;
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

  _applyExternalBodyForce(rhs_values, node_dof);
  _applyTraction(rhs_values, node_dof);

  _applyInternalBodyForce(rhs_values, node_dof);

  _applyDirichletNewton(rhs_values, node_dof);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly given as mesh type
 * following the global material tensor strategy.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_assembleBilinearOperatorGlobal()
{
  info() << "[ArcaneFem-Info] Started module  _assembleBilinearOperatorGlobal()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());
    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto in_C_tang = Accelerator::viewIn(command, m_C_tang_gp);

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2) {
      m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord, in_C_tang); });
    }
    else {
      m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_C_tang); });
    }
    m_bsr_format.toLinearSystem(m_linear_system);

  } else if (m_matrix_format == "AF-BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());
    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto in_C_tang = Accelerator::viewIn(command, m_C_tang_gp);

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2) {
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, in_C_tang, node_lid); });
    } else {
      m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_C_tang, node_lid); });
    }
    m_bsr_format.toLinearSystem(m_linear_system);

  } else if (m_matrix_format == "DOK") {
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
  } else {
    ARCANE_FATAL("Unsupported matrix type, only DOK| BSR|AF-BSR is supported.");
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"lhs-matrix-assembly", elapsedTime);
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly for a given mesh type,
 * it follows the local material tensor strategy for Von Mises Law.
 */
/*---------------------------------------------------------------------------*/
void FemModuleElastoplasticity::
_assembleBilinearOperatorLocalVonMises(bool elastic_assembly)
{
  info() << "[ArcaneFem-Info] Started module  _assembleBilinearOperatorLocalVonMises()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());

    auto in_out_dp_gp = Accelerator::viewInOut(command, m_dp_gp);
    auto in_out_sigma_gp = Accelerator::viewInOut(command, m_sigma_gp);
    auto in_out_sigma_zz_gp = Accelerator::viewInOut(command, m_sigma_zz_gp);

    auto in_sigma_old_gp = Accelerator::viewIn(command, m_sigma_old_gp);
    auto in_sigma_zz_old_gp = Accelerator::viewIn(command, m_sigma_zz_old_gp);
    auto in_p_old_gp = Accelerator::viewIn(command, m_p_old_gp);

    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto in_DUn = Accelerator::viewIn(command, m_DUn);

    auto C_elas_2d = m_C_elas_2d;
    auto in_sig0 = sig0;
    auto in_H = H;
    auto in_mu = mu;
    auto assemble_elastic = elastic_assembly;

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2) {
      m_bsr_format.assembleBilinearAtomic(
      [=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) {
        return computeLocalVonMisesElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord, in_DUn,
                                                         in_out_sigma_gp, in_out_sigma_zz_gp,
                                                         in_out_dp_gp,
                                                         in_sigma_old_gp, in_sigma_zz_old_gp,
                                                         in_p_old_gp,
                                                         C_elas_2d, in_sig0,
                                                         in_H, in_mu,
                                                         assemble_elastic); });
    }
    else {
      // m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_C_tang); });
      ARCANE_FATAL("3D not supported");
    }
    m_bsr_format.toLinearSystem(m_linear_system);
  }
  else if (m_matrix_format == "AF-BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());

    auto in_out_dp_gp = Accelerator::viewInOut(command, m_dp_gp);
    auto in_out_sigma_gp = Accelerator::viewInOut(command, m_sigma_gp);
    auto in_out_sigma_zz_gp = Accelerator::viewInOut(command, m_sigma_zz_gp);

    auto in_sigma_old_gp = Accelerator::viewIn(command, m_sigma_old_gp);
    auto in_sigma_zz_old_gp = Accelerator::viewIn(command, m_sigma_zz_old_gp);
    auto in_p_old_gp = Accelerator::viewIn(command, m_p_old_gp);

    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto in_DUn = Accelerator::viewIn(command, m_DUn);

    auto C_elas_2d = m_C_elas_2d;
    auto in_sig0 = sig0;
    auto in_H = H;
    auto in_mu = mu;
    auto assemble_elastic = elastic_assembly;

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2) {
      m_bsr_format.assembleBilinearAtomicFree(
    [=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) {
      return computeLocalVonMisesElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, in_DUn,
                                                       in_out_sigma_gp, in_out_sigma_zz_gp,
                                                       in_out_dp_gp,
                                                       in_sigma_old_gp, in_sigma_zz_old_gp,
                                                       in_p_old_gp,
                                                       C_elas_2d, in_sig0,
                                                       in_H, in_mu,
                                                       assemble_elastic,
                                                       node_lid); });
  } else {
      ARCANE_FATAL("3D not supported");
      // m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_C_tang, node_lid); });
    }
    m_bsr_format.toLinearSystem(m_linear_system);

  } else if (m_matrix_format == "DOK") {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Unsupported 2D quad DOK for local assembly");
        // _assembleBilinearOperatorCpu<8>([this](const Cell& cell) { return _computeElementMatrixQuad4(cell); });
      }
      else {
        _assembleBilinearOperatorCpu<6>([&](const Cell& cell) { return _computeLocalVonMisesElementMatrixTria3Cpu(cell, elastic_assembly); });
      }
    }
    if (mesh()->dimension() == 3) {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Unsupported 3D hexa DOK type for local assembly");
        // _assembleBilinearOperatorCpu<24>([this](const Cell& cell) { return _computeElementMatrixHexa8(cell); });
      }
      else {
        ARCANE_FATAL("Unsupported 3D tetra DOK for local assembly");
        // _assembleBilinearOperatorCpu<12>([this](const Cell& cell) { return _computeElementMatrixTetra4(cell); });
      }
    }
  } else {
    ARCANE_FATAL("Unsupported matrix type, only DOK| BSR|AF-BSR is supported.");
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "lhs-matrix-assembly", elapsedTime);
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly for a given mesh type,
 * it follows the local material tensor strategy for Von Mises Law.
 */
/*---------------------------------------------------------------------------*/
void FemModuleElastoplasticity::
_assembleBilinearOperatorLocalDruckerPrager(bool elastic_assembly)
{
  info() << "[ArcaneFem-Info] Started module  _assembleBilinearOperatorLocalDruckerPrager()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());

    auto in_out_sigma_gp = Accelerator::viewInOut(command, m_sigma_gp);
    auto in_out_sigma_zz_gp = Accelerator::viewInOut(command, m_sigma_zz_gp);
    auto in_out_eps_p_gp = Accelerator::viewInOut(command, m_eps_p_gp);
    auto in_out_eps_p_zz_gp = Accelerator::viewInOut(command, m_eps_p_zz_gp);

    auto in_eps_p_old_gp = Accelerator::viewIn(command, m_eps_p_old_gp);
    auto in_eps_p_zz_old_gp = Accelerator::viewIn(command, m_eps_p_zz_old_gp);

    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto in_DUn = Accelerator::viewIn(command, m_DUn);
    auto in_U = Accelerator::viewIn(command, m_U);

    auto C_elas_2d = m_C_elas_2d;
    auto in_bulk = bulk;
    auto in_dpEta = dpEta;
    auto in_dpC = dpC;
    auto in_mu = mu;
    auto assemble_elastic = elastic_assembly;

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2) {
    m_bsr_format.assembleBilinearAtomic(
    [=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) {
      return computeLocalDruckerPragerElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord,
                                                             in_DUn, in_U,
                                                             in_out_sigma_gp, in_out_sigma_zz_gp,
                                                             in_out_eps_p_gp, in_out_eps_p_zz_gp,
                                                             in_eps_p_old_gp, in_eps_p_zz_old_gp,
                                                             C_elas_2d,
                                                             in_bulk, in_dpEta, in_dpC, in_mu,
                                                             assemble_elastic); });
    } else {
      // m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_C_tang); });
      ARCANE_FATAL("3D not supported");
    }
    m_bsr_format.toLinearSystem(m_linear_system);
  } else if (m_matrix_format == "AF-BSR") {
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());

    auto in_out_sigma_gp = Accelerator::viewInOut(command, m_sigma_gp);
    auto in_out_sigma_zz_gp = Accelerator::viewInOut(command, m_sigma_zz_gp);
    auto in_out_eps_p_gp = Accelerator::viewInOut(command, m_eps_p_gp);
    auto in_out_eps_p_zz_gp = Accelerator::viewInOut(command, m_eps_p_zz_gp);

    auto in_eps_p_old_gp = Accelerator::viewIn(command, m_eps_p_old_gp);
    auto in_eps_p_zz_old_gp = Accelerator::viewIn(command, m_eps_p_zz_old_gp);

    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto in_DUn = Accelerator::viewIn(command, m_DUn);
    auto in_U = Accelerator::viewIn(command, m_U);

    auto C_elas_2d = m_C_elas_2d;
    auto in_bulk = bulk;
    auto in_dpEta = dpEta;
    auto in_dpC = dpC;
    auto in_mu = mu;
    auto assemble_elastic = elastic_assembly;

    m_bsr_format.computeSparsity();
    if (mesh()->dimension() == 2) {
    m_bsr_format.assembleBilinearAtomicFree(
    [=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) {
      return computeLocalDruckerPragerElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord,
                                                             in_DUn, in_U,
                                                             in_out_sigma_gp, in_out_sigma_zz_gp,
                                                             in_out_eps_p_gp, in_out_eps_p_zz_gp,
                                                             in_eps_p_old_gp, in_eps_p_zz_old_gp,
                                                             C_elas_2d,
                                                             in_bulk, in_dpEta, in_dpC, in_mu,
                                                             assemble_elastic, node_lid); });
    } else {
      ARCANE_FATAL("3D not supported");
      // m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_C_tang, node_lid); });
    }
    m_bsr_format.toLinearSystem(m_linear_system);
  }
  else if (m_matrix_format == "DOK") {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Unsupported 2D quad DOK for local assembly");
        // _assembleBilinearOperatorCpu<8>([this](const Cell& cell) { return _computeElementMatrixQuad4(cell); });
      }
      else {
        _assembleBilinearOperatorCpu<6>([&](const Cell& cell) { return _computeLocalDruckerPragerElementMatrixTria3Cpu(cell, elastic_assembly); });
      }
    }
    if (mesh()->dimension() == 3) {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Unsupported 3D hexa DOK type for local assembly");
        // _assembleBilinearOperatorCpu<24>([this](const Cell& cell) { return _computeElementMatrixHexa8(cell); });
      }
      else {
        ARCANE_FATAL("Unsupported 3D tetra DOK for local assembly");
        // _assembleBilinearOperatorCpu<12>([this](const Cell& cell) { return _computeElementMatrixTetra4(cell); });
      }
    }
  }
  else {
    ARCANE_FATAL("Unsupported matrix type, only DOK| BSR|AF-BSR is supported.");
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "lhs-matrix-assembly", elapsedTime);
}
/*---------------------------------------------------------------------------*/

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

  ENUMERATE_ (Node, inode, allNodes()) {
    Node node = *inode;
    info() << "U["<< node.uniqueId() << "] = " << m_DUn[node].x << " " << m_DUn[node].y << " " << m_DUn[node].z;
  }

  String filename = options()->solutionComparisonFile();
  const double epsilon = options()->resultEpsilon();
  const double min_value_to_test = 1.0e-10;

  Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_DUn, epsilon, min_value_to_test);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*
  * @brief Reads case tables for traction and Dirichlet boundary conditions.
  *
  * This method reads the case tables specified in the options and stores
  * them in a list for later use.
  */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_readCaseTables()
{
  IParallelMng* pm = subDomain()->parallelMng();
  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  // Keep one entry per boundary condition so that the table list and the
  // boundary-condition list always have identical indices.
  for (BC::ITractionBoundaryCondition* bs : bc->tractionBoundaryConditions()) {
    CaseTable* case_table = nullptr;
    auto traction_table_file_name = bs->getTractionInputFile();
    bool getTractionFromTable = !traction_table_file_name.empty();
    if (getTractionFromTable)
      case_table = readFileAsCaseTable(pm, traction_table_file_name, 3);
    m_traction_case_table_list.add(CaseTableInfo{ traction_table_file_name, case_table });
  }

  for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions()) {
    auto dirichlet_table_file_name = bs->getDirichletInputFile();
    if (!dirichlet_table_file_name.empty())
      m_dirichlet_case_table_list.add(readDirichletFileAsCaseTable(pm, dirichlet_table_file_name));
    else
      m_dirichlet_case_table_list.add(CaseTableInfo{ dirichlet_table_file_name, nullptr });
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Update the FEM variables in time.
 *
 * This method performs the following actions:
 *   1. Fetches values of FEM variable DU for solved time step and,
 *      adds it to global time FEM variable U.
 *   2. Performs synchronize of FEM variables across subdomains.
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::
_updateTimeVariables()
{
  info() << "[ArcaneFem-Info] Started module  _updateTimeVariables()";
  Real elapsedTime = platform::getRealTime();

  m_U.synchronize();
  m_DUn.synchronize();
  ENUMERATE_ (Node, inode, ownNodes()) {
    m_U[inode] += m_DUn[inode];
  }
  m_U.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"update-time-variables", elapsedTime);
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
        m_DUk[node] = Real3(du1_val, du2_val, du3_val);
      }
    else
      ENUMERATE_ (Node, inode, ownNodes()) {
        Node node = *inode;
        Real du1_val = dof_du[node_dof.dofId(node, 0)];
        Real du2_val = dof_du[node_dof.dofId(node, 1)];
        m_DUk[node] = Real3(du1_val, du2_val, 0.);
      }
  }
  m_DUk.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"update-Newton-increments", elapsedTime);
}

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

  m_DUk.synchronize();
  m_DUn.synchronize();
  {
      ENUMERATE_ (Node, inode, ownNodes()) {
      m_DUn[inode] += m_DUk[inode];
    }
  }
  m_DUn.synchronize();

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

  m_DUk.synchronize();
  m_DUn.synchronize();

  Real l2_norm_du = _normL2(m_DUk);
  Real l2_norm_u = _normL2(m_DUn);

  m_increment_norm = l2_norm_u != 0.0 ? l2_norm_du / l2_norm_u : 1.0;
  Real convergence_error_increment = l2_norm_du / (m_newton_rtol * l2_norm_u  + m_newton_atol);

  VariableDoFReal& residual_values(m_linear_system.rhsVariable());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  _applyZeroRHSOnConstrainedDOFs(residual_values, node_dof);
  Real l2_norm_rhs = _normL2(residual_values, node_dof);

  m_residual_norm = m_residual_norm0 !=0. ? l2_norm_rhs / (m_residual_norm0 + 1e-30) : l2_norm_rhs / (1.0 + 1e-30);
  Real convergence_error_residual = m_residual_norm;

  // The OR criterion follows petsc SNES
  if (convergence_error_residual <= m_newton_rtol) {
    m_newton_solver_converged = true;
    m_newton_converged_reason = "RESIDUAL_CONVERGED";
    info() << "[ArcaneFem-Info] At Newton iteration " << m_newton_iter
           << ": ||X_k+1 - X_k||/||X_k+1|| = " << m_increment_norm
           << " ||(F_ext - F_int(X_k))||/||(F_ext - F_int(X_0))|| = " << m_residual_norm
           << " => " << "CONVERGED with " << m_newton_converged_reason;
  } else if (convergence_error_increment <= 1.0) {
    m_newton_solver_converged = true;
    m_newton_converged_reason = "INCREMENT_CONVERGED";
    info() << "[ArcaneFem-Info] At Newton iteration " << m_newton_iter
           << ": ||X_k+1 - X_k||/||X_k+1|| = " << m_increment_norm
           << " ||(F_ext - F_int(X_k))||/||(F_ext - F_int(X_0))|| = " << m_residual_norm
           << " => " << "CONVERGED with " << m_newton_converged_reason;
  } else {
    m_newton_solver_converged = false;
    info() << "[ArcaneFem-Info] At Newton iteration " << m_newton_iter
           << ": ||X_k+1 - X_k||/||X_k+1|| = " << m_increment_norm
           << " ||(F_ext - F_int(X_k))||/||(F_ext - F_int(X_0))|| = " << m_residual_norm
           << " => " << "NOT CONVERGED";
  }

    elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "check-newton-convergence", elapsedTime);
}

inline Real FemModuleElastoplasticity::
_normL2(VariableNodeReal3& u) {
  Real l2_norm_u = 0.0;
  {
    ENUMERATE_ (Node, inode, ownNodes()) {
      const Real norm_u = math::pow(u[inode][0], 2.0) + math::pow(u[inode][1], 2.0) + math::pow(u[inode][2], 2.0);
      l2_norm_u += norm_u;
    }
  }
  IParallelMng* pm = defaultMesh()->parallelMng();
  l2_norm_u = pm->reduce(Parallel::ReduceSum, l2_norm_u);

  return math::sqrt(l2_norm_u);
}

inline Real FemModuleElastoplasticity::
_normL2(VariableDoFReal& u, const IndexedNodeDoFConnectivityView& node_dof) {
  Real l2_norm_u = 0.0;
  Int32 mesh_dimension = mesh()->dimension();
  {
    ENUMERATE_ (Node, inode, ownNodes()) {
      Real norm_residual = 0.0;
      if (mesh_dimension == 2) {
        norm_residual =  math::pow(u[node_dof.dofId(inode, 0)], 2.0)
                  + math::pow(u[node_dof.dofId(inode, 1)], 2.0);
      } else {
        norm_residual =  math::pow(u[node_dof.dofId(inode, 0)], 2.0)
                  + math::pow(u[node_dof.dofId(inode, 1)], 2.0)
                  + math::pow(u[node_dof.dofId(inode, 2)], 2.0);
      }
      l2_norm_u += norm_residual;
    }
  }
  IParallelMng* pm = defaultMesh()->parallelMng();
  l2_norm_u = pm->reduce(Parallel::ReduceSum, l2_norm_u);
  return math::sqrt(l2_norm_u);
}

inline Real FemModuleElastoplasticity::
_normL1(VariableDoFReal& u, const IndexedNodeDoFConnectivityView& node_dof) {
  Real l1_norm_u = 0.0;
  Int32 mesh_dimension = mesh()->dimension();
  {
    ENUMERATE_ (Node, inode, ownNodes()) {
      Real norm_residual = 0.0;
      if (mesh_dimension == 2) {
        norm_residual = u[node_dof.dofId(inode, 0)]
                      + u[node_dof.dofId(inode, 1)];
      } else {
        norm_residual = u[node_dof.dofId(inode, 0)]
                      + u[node_dof.dofId(inode, 1)]
                      + u[node_dof.dofId(inode, 2)];
      }
      l1_norm_u += norm_residual;
    }
  }
  IParallelMng* pm = defaultMesh()->parallelMng();
  l1_norm_u = pm->reduce(Parallel::ReduceSum, l1_norm_u);
  return l1_norm_u;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModuleElastoplasticity);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
