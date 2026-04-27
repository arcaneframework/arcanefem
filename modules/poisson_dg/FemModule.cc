// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2026 */
/*                                                                           */
/* Poisson solver module of ArcaneFEM.                                       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModulePoisson at the start of the simulation.
 *
 * This method initializes degrees of freedom (DoFs) on nodes.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
startInit()
{
  info() << "[ArcaneFem-Info] Started module startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dofs_on_cells.initialize(mesh(), 3); // P1 DG elements have 3 DoFs per cell
  m_dof_family = m_dofs_on_cells.dofFamily();

  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->hasSolutionComparisonFile();
  m_petsc_flags = options()->petscFlags();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs the main computation for the FemModulePoisson.
 *
 * This method:
 *   1. Stops the time loop after 1 iteration since the equation is steady state.
 *   2. Resets, configures, and initializes the linear system.
 *   3. Executes the stationary solve.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
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
  m_linear_system.initialize(subDomain(), acceleratorMng()->defaultRunner(), m_dof_family, "Solver");

  if (m_petsc_flags != NULL) {
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows a sequence of steps to solve FEM system:
 *
 *   1. _getMaterialParameters()     Retrieves material parameters via
 *   2. _assembleLinearSystem()      Assembles the FEM  matrix A RHS vector b
 *   3. _solve()                     Solves for solution vector u = A^-1*b
 *   4. _updateVariables()           Updates FEM variables u = x
 *   5. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
_doStationarySolve()
{
  _getMaterialParameters();

  if (m_assemble_linear_system) {
    _assembleLinearSystem();
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
 * @brief Retrieves and sets the material parameters for the simulation.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  f = options()->f();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the FEM linear system (matrix and RHS) for the Poisson 
 * problem using DG formulation.
 */
 /*---------------------------------------------------------------------------*/

void FemModulePoisson::
_assembleLinearSystem()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearSystem()";
  Real elapsedTime = platform::getRealTime();

  // For DG, we use DOK format
  auto cell_dof = m_dofs_on_cells.cellDoFConnectivityView();
  VariableDoFReal& rhs_values = m_linear_system.rhsVariable();
  rhs_values.fill(0.0);

  Real penalty = 10.0; // Can be made configurable later

  // Volume terms: (∇v, ∇u)_K
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
    // Compute centroid
    Real3 centroid = {0.0, 0.0, 0.0};
    for (Node node : cell.nodes()) {
      centroid += m_node_coord[node];
    }
    centroid /= cell.nbNode();

    // Evaluate basis at centroid for volume integral approximation
    // Basis: phi_0 = 1, phi_1 = x - x_c, phi_2 = y - y_c
    Real x_rel = centroid.x - centroid.x;
    Real y_rel = centroid.y - centroid.y;
    Real phi[3] = {1.0, x_rel, y_rel};
    Real grad_x[3] = {0.0, 1.0, 0.0};
    Real grad_y[3] = {0.0, 0.0, 1.0};

    // RHS: (f, phi)_K
    for (Int32 i = 0; i < 3; ++i) {
      DoFLocalId dof_i = cell_dof.dofId(cell, i);
      rhs_values[dof_i] += f * phi[i] * area;
    }

    // Stiffness: ∫ ∇phi_i · ∇phi_j dK
    for (Int32 i = 0; i < 3; ++i) {
      for (Int32 j = 0; j < 3; ++j) {
        DoFLocalId dof_i = cell_dof.dofId(cell, i);
        DoFLocalId dof_j = cell_dof.dofId(cell, j);
        Real stiff = (grad_x[i] * grad_x[j] + grad_y[i] * grad_y[j]) * area;
        m_linear_system.matrixAddValue(dof_i, dof_j, stiff);
      }
    }
  }

  // Face terms (SIPG)
  ENUMERATE_ (Face, iface, allFaces()) {
    Face face = *iface;
    // Compute length
    Real3 p0 = m_node_coord[face.nodeId(0)];
    Real3 p1 = m_node_coord[face.nodeId(1)];
    Real length = (p1 - p0).normL2();
    // Compute centroid
    Real3 face_center = (p0 + p1) * 0.5;

    if (face.nbCell() == 2) { // Interior face
      Cell cell_i = face.cell(0);
      Cell cell_j = face.cell(1);

      // Compute normal (assuming 2D, outward from cell_i)
      Real3 edge_vec = p1 - p0;
      Real3 normal = {edge_vec.y, -edge_vec.x, 0.0};
      normal = normal / normal.normL2();
      // Ensure pointing from i to j
      Real3 cell_i_cent = {0.0, 0.0, 0.0};
      for (Node node : cell_i.nodes()) cell_i_cent += m_node_coord[node];
      cell_i_cent /= cell_i.nbNode();
      Real3 cell_j_cent = {0.0, 0.0, 0.0};
      for (Node node : cell_j.nodes()) cell_j_cent += m_node_coord[node];
      cell_j_cent /= cell_j.nbNode();
      if (math::dot(normal, cell_j_cent - cell_i_cent) < 0) normal = -normal;

      // Evaluate basis at face center
      Real3 cent_i = {0.0, 0.0, 0.0};
      for (Node node : cell_i.nodes()) cent_i += m_node_coord[node];
      cent_i /= cell_i.nbNode();
      Real3 cent_j = {0.0, 0.0, 0.0};
      for (Node node : cell_j.nodes()) cent_j += m_node_coord[node];
      cent_j /= cell_j.nbNode();
      Real x_rel_i = face_center.x - cent_i.x;
      Real y_rel_i = face_center.y - cent_i.y;
      Real phi_i[3] = {1.0, x_rel_i, y_rel_i};
      Real grad_x_i[3] = {0.0, 1.0, 0.0};
      Real grad_y_i[3] = {0.0, 0.0, 1.0};
      Real grad_n_i[3] = {grad_x_i[0] * normal.x + grad_y_i[0] * normal.y,
                          grad_x_i[1] * normal.x + grad_y_i[1] * normal.y,
                          grad_x_i[2] * normal.x + grad_y_i[2] * normal.y};

      Real x_rel_j = face_center.x - cent_j.x;
      Real y_rel_j = face_center.y - cent_j.y;
      Real phi_j[3] = {1.0, x_rel_j, y_rel_j};
      Real grad_x_j[3] = {0.0, 1.0, 0.0};
      Real grad_y_j[3] = {0.0, 0.0, 1.0};
      Real grad_n_j[3] = {grad_x_j[0] * normal.x + grad_y_j[0] * normal.y,
                          grad_x_j[1] * normal.x + grad_y_j[1] * normal.y,
                          grad_x_j[2] * normal.x + grad_y_j[2] * normal.y};

      // Penalty parameter σ = γ/h
      // Compute diameter as max distance between nodes
      Real diam_i = 0.0;
      for (Int32 n1 = 0; n1 < cell_i.nbNode(); ++n1) {
        for (Int32 n2 = n1 + 1; n2 < cell_i.nbNode(); ++n2) {
          Real dist = (m_node_coord[cell_i.nodeId(n1)] - m_node_coord[cell_i.nodeId(n2)]).normL2();
          diam_i = math::max(diam_i, dist);
        }
      }
      Real diam_j = 0.0;
      for (Int32 n1 = 0; n1 < cell_j.nbNode(); ++n1) {
        for (Int32 n2 = n1 + 1; n2 < cell_j.nbNode(); ++n2) {
          Real dist = (m_node_coord[cell_j.nodeId(n1)] - m_node_coord[cell_j.nodeId(n2)]).normL2();
          diam_j = math::max(diam_j, dist);
        }
      }
      Real h = math::min(diam_i, diam_j);
      Real sigma = penalty / h;

      // SIPG terms
      for (Int32 i = 0; i < 3; ++i) {
        for (Int32 j = 0; j < 3; ++j) {
          DoFLocalId dof_i_i = cell_dof.dofId(cell_i, i);
          DoFLocalId dof_j_i = cell_dof.dofId(cell_i, j);
          DoFLocalId dof_i_j = cell_dof.dofId(cell_j, i);
          DoFLocalId dof_j_j = cell_dof.dofId(cell_j, j);

          // - <[[v]], {∇u·n}> = - <v_i - v_j, 0.5(∇u_i·n + ∇u_j·n)>
          Real term1_ii = -0.5 * phi_i[i] * grad_n_i[j] * length;
          Real term1_ij = -0.5 * phi_i[i] * grad_n_j[j] * length;
          Real term1_ji = 0.5 * phi_j[i] * grad_n_i[j] * length;
          Real term1_jj = 0.5 * phi_j[i] * grad_n_j[j] * length;

          // - <{∇v·n}, [[u]]> = - <0.5(∇v_i·n + ∇v_j·n), u_i - u_j>
          Real term2_ii = -0.5 * grad_n_i[i] * phi_i[j] * length;
          Real term2_ij = 0.5 * grad_n_i[i] * phi_j[j] * length;
          Real term2_ji = -0.5 * grad_n_j[i] * phi_i[j] * length;
          Real term2_jj = 0.5 * grad_n_j[i] * phi_j[j] * length;

          // + <σ[[v]], [[u]]> = <σ(v_i - v_j), u_i - u_j>
          Real term3_ii = sigma * phi_i[i] * phi_i[j] * length;
          Real term3_ij = -sigma * phi_i[i] * phi_j[j] * length;
          Real term3_ji = -sigma * phi_j[i] * phi_i[j] * length;
          Real term3_jj = sigma * phi_j[i] * phi_j[j] * length;

          m_linear_system.matrixAddValue(dof_i_i, dof_j_i, term1_ii + term2_ii + term3_ii);
          m_linear_system.matrixAddValue(dof_i_i, dof_j_j, term1_ij + term2_ij + term3_ij);
          m_linear_system.matrixAddValue(dof_i_j, dof_j_i, term1_ji + term2_ji + term3_ji);
          m_linear_system.matrixAddValue(dof_i_j, dof_j_j, term1_jj + term2_jj + term3_jj);
        }
      }
    } else { // Boundary face
      Cell cell_i = face.cell(0);
      // Compute normal (outward)
      Real3 edge_vec = p1 - p0;
      Real3 normal = {edge_vec.y, -edge_vec.x, 0.0};
      normal = normal / normal.normL2();
      // Ensure outward
      Real3 cell_cent = {0.0, 0.0, 0.0};
      for (Node node : cell_i.nodes()) cell_cent += m_node_coord[node];
      cell_cent /= cell_i.nbNode();
      Real3 face_cent = (p0 + p1) * 0.5;
      if (math::dot(normal, face_cent - cell_cent) < 0) normal = -normal;

      Real3 cent_i = {0.0, 0.0, 0.0};
      for (Node node : cell_i.nodes()) cent_i += m_node_coord[node];
      cent_i /= cell_i.nbNode();
      Real x_rel_i = face_center.x - cent_i.x;
      Real y_rel_i = face_center.y - cent_i.y;
      Real phi_i[3] = {1.0, x_rel_i, y_rel_i};
      Real grad_x_i[3] = {0.0, 1.0, 0.0};
      Real grad_y_i[3] = {0.0, 0.0, 1.0};
      Real grad_n_i[3] = {grad_x_i[0] * normal.x + grad_y_i[0] * normal.y,
                          grad_x_i[1] * normal.x + grad_y_i[1] * normal.y,
                          grad_x_i[2] * normal.x + grad_y_i[2] * normal.y};

      Real diam_i = 0.0;
      for (Int32 n1 = 0; n1 < cell_i.nbNode(); ++n1) {
        for (Int32 n2 = n1 + 1; n2 < cell_i.nbNode(); ++n2) {
          Real dist = (m_node_coord[cell_i.nodeId(n1)] - m_node_coord[cell_i.nodeId(n2)]).normL2();
          diam_i = math::max(diam_i, dist);
        }
      }
      Real h = diam_i;
      Real sigma = penalty / h;

      // Check if this face has Dirichlet BC
      // TODO: Add Dirchelt and Neumann support later
      Real g = 0.0; // default homogeneous

      // No Dirichlet BC - use homogeneous Neumann (natural BC)
      for (Int32 i = 0; i < 3; ++i) {
        for (Int32 j = 0; j < 3; ++j) {
          DoFLocalId dof_i_i = cell_dof.dofId(cell_i, i);
          DoFLocalId dof_j_i = cell_dof.dofId(cell_i, j);

          // - <v, ∇u·n>
          m_linear_system.matrixAddValue(dof_i_i, dof_j_i, -phi_i[i] * grad_n_i[j] * length);

          // - <∇v·n, u>
          m_linear_system.matrixAddValue(dof_i_i, dof_j_i, -grad_n_i[i] * phi_i[j] * length);

          // + <σv, u>
          m_linear_system.matrixAddValue(dof_i_i, dof_j_i, sigma * phi_i[i] * phi_i[j] * length);
        }

        // RHS: - <∇v·n, g> + <σv, g> (g=0 for homogeneous Neumann)
        DoFLocalId dof_i = cell_dof.dofId(cell_i, i);
        rhs_values[dof_i] -= grad_n_i[i] * g * length;
        rhs_values[dof_i] += sigma * phi_i[i] * g * length;
      }
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"lhs-matrix-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Solves the linear system.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
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
 *      i.e., it copies RHS DOF to u.
 *   2. Performs synchronize of FEM variables across subdomains.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
_updateVariables()
{
  info() << "[ArcaneFem-Info] Started module _updateVariables()";
  Real elapsedTime = platform::getRealTime();

  {
    // For DG, interpolate cell-based solution to node-based P1 field for visualization
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto cell_dof = m_dofs_on_cells.cellDoFConnectivityView();

    // Use vectors for temporary node values
    Int32 max_node_id = mesh()->nodeFamily()->maxLocalId();
    std::vector<Real> node_values(max_node_id, 0.0);
    std::vector<Int32> node_count(max_node_id, 0);

    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;

      // Get cell centroid
      Real3 centroid = { 0.0, 0.0, 0.0 };
      for (Node node : cell.nodes()) {
        centroid += m_node_coord[node];
      }
      centroid /= cell.nbNode();

      // Get DG coefficients for this cell
      Real a0 = dof_u[cell_dof.dofId(cell, 0)]; // constant term
      Real a1 = dof_u[cell_dof.dofId(cell, 1)]; // x-gradient
      Real a2 = dof_u[cell_dof.dofId(cell, 2)]; // y-gradient

      // Evaluate solution at each node of this cell
      for (Node node : cell.nodes()) {
        Real3 node_pos = m_node_coord[node];
        Real x_rel = node_pos.x - centroid.x;
        Real y_rel = node_pos.y - centroid.y;
        Real u_value = a0 + a1 * x_rel + a2 * y_rel;

        node_values[node.localId()] += u_value;
        node_count[node.localId()] += 1;
      }
    }

    // Average the values at nodes
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      if (node_count[node.localId()] > 0) {
        m_u[node] = node_values[node.localId()] / node_count[node.localId()];
      }
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

void FemModulePoisson::
_validateResults()
{
  info() << "[ArcaneFem-Info] Started module _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200)
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "u["  << node.uniqueId() << "] = " << m_u[node];
    }

  String filename = options()->solutionComparisonFile();

  checkNodeResultFile(traceMng(), filename, m_u, 1.0e-4);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModulePoisson);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/