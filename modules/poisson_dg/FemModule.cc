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
 * This method is run at the beginning of simulation and performs the
 * following tasks:
 *  1. Initializes the DoF connectivity for P1 DG elements (3 per cell).
 *  2. Reads configuration options for matrix format, assembly, solving, 
 *     and validation.
 *  3. Builds cell-cell connectivity if using CSR format.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
startInit()
{
  info() << "[ArcaneFem-Info] Started module startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dimension = mesh()->dimension();
  if (m_dimension == 2) {
    m_compute_dg_penalty_length = &ArcaneFemFunctions::MeshOperation::computeDGPenaltyLength2D;
    m_compute_face_normal = &ArcaneFemFunctions::MeshOperation::computeOutwardUnitNormalEdge2D;
    m_compute_cell_quadrature = &ArcaneFemFunctions::DgQuadrature::computeCellQuadrature2D;
    m_compute_face_quadrature = &ArcaneFemFunctions::DgQuadrature::computeFaceQuadrature2D;
  }
  else if (m_dimension == 3) {
    m_compute_dg_penalty_length = &ArcaneFemFunctions::MeshOperation::computeDGPenaltyLength3D;
    m_compute_face_normal = &ArcaneFemFunctions::MeshOperation::computeOutwardUnitNormalPolygon3D;
    m_compute_cell_quadrature = &ArcaneFemFunctions::DgQuadrature::computeCellQuadrature3D;
    m_compute_face_quadrature = &ArcaneFemFunctions::DgQuadrature::computeFaceQuadrature3D;
  }

  m_nb_dof_per_cell = m_dimension + 1;
  m_dofs_on_cells.initialize(mesh(), m_nb_dof_per_cell);
  m_dof_family = m_dofs_on_cells.dofFamily();

  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->hasSolutionComparisonFile();
  m_petsc_flags = options()->petscFlags();

  // Build cell-cell connectivity if using CSR format for efficient assembly
  if (m_matrix_format == "CSR") {
    IItemFamily* cell_family = mesh()->cellFamily();
    auto* cn = new mesh::IncrementalItemConnectivity(cell_family, cell_family, "NeighbourCellCell");
    ENUMERATE_CELL (icell, allCells()) {
      Cell cell = *icell;
      cn->notifySourceItemAdded(cell);
      for (Face face : cell.faces()) {
        if (face.nbCell() != 2)
          continue;
        Cell opposite_cell = face.oppositeCell(cell);
        if (!opposite_cell.null())
          cn->addConnectedItem(cell, opposite_cell);
      }
    }
    m_cell_cell_connectivity_view = cn->connectivityView();
  }

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
  m_penalty = options()->penalty();

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

  auto cell_dof = m_dofs_on_cells.cellDoFConnectivityView();
  VariableDoFReal& rhs_values = m_linear_system.rhsVariable();
  rhs_values.fill(0.0);

  if (m_matrix_format == "CSR")
    _buildCsrSparsity();

  auto add_matrix_value = [&](DoFLocalId row, DoFLocalId column, Real value) {
    if (m_matrix_format == "CSR")
      m_csr_matrix.matrixAddValue(row, column, value);
    else
      m_linear_system.matrixAddValue(row, column, value);
  };

  const Real3 gradients[4] = { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 },
                               { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

  // Volume terms: (∇v, ∇u)_K and (v,f)_K.
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real3 centroid = ArcaneFemFunctions::MeshOperation::computeCentroid(cell, m_node_coord);
    UniqueArray<QuadraturePoint> quadrature = m_compute_cell_quadrature(cell, m_node_coord);
    Real measure = 0.0;

    for (const QuadraturePoint& qp : quadrature) {
      measure += qp.weight;
      Real3 relative = qp.point - centroid;
      Real phi[4] = { 1.0, relative.x, relative.y, relative.z };
      for (Int32 i = 0; i < m_nb_dof_per_cell; ++i)
        rhs_values[cell_dof.dofId(cell, i)] += f * phi[i] * qp.weight;
    }

    for (Int32 i = 0; i < m_nb_dof_per_cell; ++i) {
      for (Int32 j = 0; j < m_nb_dof_per_cell; ++j) {
        Real stiff = math::dot(gradients[i], gradients[j]) * measure;
        add_matrix_value(cell_dof.dofId(cell, i), cell_dof.dofId(cell, j), stiff);
      }
    }
  }

  // Interior-face SIPG terms.
  ENUMERATE_ (Face, iface, allFaces()) {
    Face face = *iface;
    if (face.nbCell() != 2)
      continue;

    Cell cell_i = face.cell(0);
    Cell cell_j = face.cell(1);
    Real3 center_i = ArcaneFemFunctions::MeshOperation::computeCentroid(cell_i, m_node_coord);
    Real3 center_j = ArcaneFemFunctions::MeshOperation::computeCentroid(cell_j, m_node_coord);
    Real3 normal = m_compute_face_normal(face, cell_i, m_node_coord);
    Real grad_n[4];
    for (Int32 i = 0; i < m_nb_dof_per_cell; ++i)
      grad_n[i] = math::dot(gradients[i], normal);

    Real h_i = m_compute_dg_penalty_length(cell_i, face, m_node_coord);
    Real h_j = m_compute_dg_penalty_length(cell_j, face, m_node_coord);
    Real sigma = m_penalty / math::min(h_i, h_j);

    for (const QuadraturePoint& qp : m_compute_face_quadrature(face, m_node_coord)) {
      Real3 relative_i = qp.point - center_i;
      Real3 relative_j = qp.point - center_j;
      Real phi_i[4] = { 1.0, relative_i.x, relative_i.y, relative_i.z };
      Real phi_j[4] = { 1.0, relative_j.x, relative_j.y, relative_j.z };

      for (Int32 i = 0; i < m_nb_dof_per_cell; ++i) {
        for (Int32 j = 0; j < m_nb_dof_per_cell; ++j) {
          DoFLocalId row_i = cell_dof.dofId(cell_i, i);
          DoFLocalId row_j = cell_dof.dofId(cell_j, i);
          DoFLocalId column_i = cell_dof.dofId(cell_i, j);
          DoFLocalId column_j = cell_dof.dofId(cell_j, j);
          Real weight = qp.weight;

          // - <[[v]], {∇u·n}> = - <v_i - v_j, 0.5(∇u_i·n + ∇u_j·n)>
          // - <{∇v·n}, [[u]]> = - <0.5(∇v_i·n + ∇v_j·n), u_i - u_j>
          // + <σ[[v]], [[u]]> = <σ(v_i - v_j), u_i - u_j>
          add_matrix_value(row_i, column_i, weight * (
            -0.5 * phi_i[i] * grad_n[j] - 0.5 * grad_n[i] * phi_i[j]
            + sigma * phi_i[i] * phi_i[j]));
          add_matrix_value(row_i, column_j, weight * (
            -0.5 * phi_i[i] * grad_n[j] + 0.5 * grad_n[i] * phi_j[j]
            - sigma * phi_i[i] * phi_j[j]));
          add_matrix_value(row_j, column_i, weight * (
            0.5 * phi_j[i] * grad_n[j] - 0.5 * grad_n[i] * phi_i[j]
            - sigma * phi_j[i] * phi_i[j]));
          add_matrix_value(row_j, column_j, weight * (
            0.5 * phi_j[i] * grad_n[j] + 0.5 * grad_n[i] * phi_j[j]
            + sigma * phi_j[i] * phi_j[j]));
        }
      }
    }
  }

  BC::IArcaneFemBC* bc = options()->boundaryConditions();
  if (bc) { // Only process if BCs are defined
    // Loop over Dirichlet BCs and apply them as penalty terms in SIPG
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions()) {
      FaceGroup face_group = bs->getSurface();

      // Retrieve the Dirichlet value and convert to Real
      const StringConstArrayView value = bs->getValue();
      Real g = 0.0;
      if (builtInGetValue(g,value[0]))
          ARCANE_FATAL("Can not convert '{0}' to real",value[0]);

      ENUMERATE_ (Face, iface, face_group) {
        Face face = *iface;
        Cell cell = face.cell(0);
        Real3 center = ArcaneFemFunctions::MeshOperation::computeCentroid(cell, m_node_coord);
        Real3 normal = m_compute_face_normal(face, cell, m_node_coord);
        Real grad_n[4];
        for (Int32 i = 0; i < m_nb_dof_per_cell; ++i)
          grad_n[i] = math::dot(gradients[i], normal);

        Real h = m_compute_dg_penalty_length(cell, face, m_node_coord);
        Real sigma = m_penalty / h;

        for (const QuadraturePoint& qp : m_compute_face_quadrature(face, m_node_coord)) {
          Real3 relative = qp.point - center;
          Real phi[4] = { 1.0, relative.x, relative.y, relative.z };
          for (Int32 i = 0; i < m_nb_dof_per_cell; ++i) {
            DoFLocalId row = cell_dof.dofId(cell, i);
            for (Int32 j = 0; j < m_nb_dof_per_cell; ++j) {
              DoFLocalId column = cell_dof.dofId(cell, j);
              add_matrix_value(row, column, qp.weight * (
                -phi[i] * grad_n[j] - grad_n[i] * phi[j] + sigma * phi[i] * phi[j]));
            }
            rhs_values[row] += qp.weight * (-grad_n[i] * g + sigma * phi[i] * g);
          }
        }
      }
    }

    for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions()) {
      FaceGroup face_group = bs->getSurface();

      // Retrieve the Dirichlet value and convert to Real
      const StringConstArrayView value = bs->getValue();
      Real g = 0.0;
      if (builtInGetValue(g,value[0]))
          ARCANE_FATAL("Can not convert '{0}' to real",value[0]);

      ENUMERATE_ (Face, iface, face_group) {
        Face face = *iface;
        Cell cell = face.cell(0);
        Real3 center = ArcaneFemFunctions::MeshOperation::computeCentroid(cell, m_node_coord);

        for (const QuadraturePoint& qp : m_compute_face_quadrature(face, m_node_coord)) {
          Real3 relative = qp.point - center;
          Real phi[4] = { 1.0, relative.x, relative.y, relative.z };
          for (Int32 i = 0; i < m_nb_dof_per_cell; ++i)
            rhs_values[cell_dof.dofId(cell, i)] += phi[i] * g * qp.weight;
        }
      }
    }
  }

  if (m_matrix_format == "CSR") {
    RunQueue* queue = subDomain()->acceleratorMng()->defaultQueue();
    m_csr_matrix.translateToLinearSystem(m_linear_system, *queue);
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "lhs-matrix-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Builds the CSR sparsity pattern for the DG matrix.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
_buildCsrSparsity()
{
  auto cell_dof = m_dofs_on_cells.cellDoFConnectivityView();
  CellInfoListView cells(mesh()->cellFamily());
  const Int32 nb_dof = m_dof_family->nbItem();
  Int32 nnz = 0;

  ENUMERATE_CELL (icell, allCells()) {
    Int32 nb_connected_cell = m_cell_cell_connectivity_view.nbCell(icell);
    nnz += m_nb_dof_per_cell * m_nb_dof_per_cell * (1 + nb_connected_cell);
  }

  RunQueue queue = subDomain()->acceleratorMng()->queue();

  NumArray<Int32, MDDim1> rows_index(queue.memoryResource());
  NumArray<Int32, MDDim1> columns(queue.memoryResource());
  rows_index.resize(nb_dof + 1);
  columns.resize(nnz);

  Int32 index = 0;
  Int32 column_index = 0;
  rows_index[0] = 0;
  ENUMERATE_CELL (icell, allCells()) {
    Cell cell = *icell;
    for (Int32 i = 0; i < m_nb_dof_per_cell; ++i) {
      DoFLocalId row_dof = cell_dof.dofId(cell, i);
      for (Int32 j = 0; j < m_nb_dof_per_cell; ++j) {
        DoFLocalId column_dof = cell_dof.dofId(cell, j);
        columns[column_index++] = column_dof;
      }
      for (CellLocalId neighbor_cell_id : m_cell_cell_connectivity_view.cells(icell)) {
        Cell neighbor_cell = cells[neighbor_cell_id];
        for (Int32 j = 0; j < m_nb_dof_per_cell; ++j) {
          DoFLocalId column_dof = cell_dof.dofId(neighbor_cell, j);
          columns[column_index++] = column_dof;
        }
      }
      rows_index[++index] = column_index;
    }
  }
  m_csr_matrix.initialize(m_dof_family, std::move(rows_index), std::move(columns), queue);
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
 *   1. Fetches values of solution from solved linear system to FEM variables.
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

    dof_u.synchronize(); // Ensure solution is up to date across subdomains before interpolation
    m_u.fill(0.0);

    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real3 centroid =  ArcaneFemFunctions::MeshOperation::computeCentroid(cell, m_node_coord);

      Real coefficients[4] = {};
      for (Int32 i = 0; i < m_nb_dof_per_cell; ++i)
        coefficients[i] = dof_u[cell_dof.dofId(cell, i)];

      // Evaluate solution at each node of this cell
      for (Node node : cell.nodes()) {
        Real3 relative = m_node_coord[node] - centroid;
        Real u_value = coefficients[0];
        for (Int32 d = 0; d < m_dimension; ++d)
          u_value += coefficients[d + 1] * relative[d];

        m_u[node] += u_value / node.nbCell(); // Average contributions from all cells sharing this node
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
 *   1. Prints the computed values for each node.
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

  ENUMERATE_ (Node, inode, allNodes()) {
    Node node = *inode;
    info() << "u[" << node.uniqueId() << "] = " << m_u[node];
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
