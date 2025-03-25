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
  m_cross_validation = options()->crossValidation();
  m_petsc_flags = options()->petscFlags();

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") {
    auto use_csr_in_linear_system = options()->linearSystem.serviceName() == "HypreLinearSystem";
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
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

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
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real2 grad = _computeGradientOfRealTria3(cell);
      m_psi[cell] = -grad.x * grad.x - grad.y * grad.y;
    }

  if (mesh()->dimension() == 3)
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real3 grad = _computeGradientOfRealTetra4(cell);
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
 *   3.  _solve()                         Solves for solution vector 𝒖 = 𝑨⁻¹𝒃
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
 * @brief FEM linear operator for the current simulation step.
 *
 * This method constructs the FEM linear  systems RHS vector by applying
 * various boundary conditions and source terms.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR")
    m_bsr_format.toLinearSystem(m_linear_system);

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  BC::IArcaneFemBC* bc = options()->boundaryConditions();
  if (bc)
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
      ArcaneFemFunctions::BoundaryConditions2D::applyDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);

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
/*---------------------------------------------------------------------------*/

Real3 FemModule::
_computeGradientOfRealTetra4(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  Real3 m3 = m_node_coord[cell.nodeId(3)];

  Real f0 = m_u[cell.nodeId(0)];
  Real f1 = m_u[cell.nodeId(1)];
  Real f2 = m_u[cell.nodeId(2)];
  Real f3 = m_u[cell.nodeId(3)];

  Real3 v0 = m1 - m0;
  Real3 v1 = m2 - m0;
  Real3 v2 = m3 - m0;

  // 6 x Volume of tetrahedron
  Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

  // Compute gradient components
  Real3 grad;
  grad.x = (f0 * (m1.y * m2.z + m2.y * m3.z + m3.y * m1.z - m3.y * m2.z - m2.y * m1.z - m1.y * m3.z)
           - f1 * (m0.y * m2.z + m2.y * m3.z + m3.y * m0.z - m3.y * m2.z - m2.y * m0.z - m0.y * m3.z)
           + f2 * (m0.y * m1.z + m1.y * m3.z + m3.y * m0.z - m3.y * m1.z - m1.y * m0.z - m0.y * m3.z)
           - f3 * (m0.y * m1.z + m1.y * m2.z + m2.y * m0.z - m2.y * m1.z - m1.y * m0.z - m0.y * m2.z)) / V6;
  grad.y = (f0 * (m1.z * m2.x + m2.z * m3.x + m3.z * m1.x - m3.z * m2.x - m2.z * m1.x - m1.z * m3.x)
           - f1 * (m0.z * m2.x + m2.z * m3.x + m3.z * m0.x - m3.z * m2.x - m2.z * m0.x - m0.z * m3.x)
           + f2 * (m0.z * m1.x + m1.z * m3.x + m3.z * m0.x - m3.z * m1.x - m1.z * m0.x - m0.z * m3.x)
           - f3 * (m0.z * m1.x + m1.z * m2.x + m2.z * m0.x - m2.z * m1.x - m1.z * m0.x - m0.z * m2.x)) / V6;
  grad.z = (f0 * (m1.x * m2.y + m2.x * m3.y + m3.x * m1.y - m3.x * m2.y - m2.x * m1.y - m1.x * m3.y)
           - f1 * (m0.x * m2.y + m2.x * m3.y + m3.x * m0.y - m3.x * m2.y - m2.x * m0.y - m0.x * m3.y)
           + f2 * (m0.x * m1.y + m1.x * m3.y + m3.x * m0.y - m3.x * m1.y - m1.x * m0.y - m0.x * m3.y)
           - f3 * (m0.x * m1.y + m1.x * m2.y + m2.x * m0.y - m2.x * m1.y - m1.x * m0.y - m0.x * m2.y)) / V6;
  return grad;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real2 FemModule::
_computeGradientOfRealTria3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real f0 = m_u[cell.nodeId(0)];
  Real f1 = m_u[cell.nodeId(1)];
  Real f2 = m_u[cell.nodeId(2)];

  Real detA = (m0.x * (m1.y - m2.y) - m0.y * (m1.x - m2.x) + (m1.x * m2.y - m2.x * m1.y));

  Real2 grad;
  grad.x = (m0.x * (f1 - f2) - f0 * (m1.x - m2.x) + (f2 * m1.x - f1 * m2.x)) / detA;
  grad.y = (f0 * (m1.y - m2.y) - m0.y * (f1 - f2) + (f1 * m2.y - f2 * m1.y)) / detA;

  return grad;
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
      _assembleBilinear<4>([this](const Cell& cell) {
        return _computeElementMatrixTetra4(cell);
      });
    if (mesh()->dimension() == 2)
      _assembleBilinear<3>([this](const Cell& cell) {
        return _computeElementMatrixTria3(cell);
      });
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
      info() << "u[" << node.localId() << "][" << node.uniqueId() << "] = " << m_u[node];
    }

  String filename = options()->resultFile();
  info() << "ValidateResultFile filename=" << filename;

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
