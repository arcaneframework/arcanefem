// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Heat equation solver module of ArcaneFEM.                                 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include "ElementMatrix.h"

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

  bool use_csr_in_linearsystem = options()->linearSystem.serviceName() == "HypreLinearSystem";
  if (m_matrix_format == "BSR"){
    m_bsr_format.initialize(mesh(), 1, use_csr_in_linearsystem, 0);
    m_bsr_format.computeSparsity();
  }

  _initTime(); // initialize time
  _getParameters(); // get material parameters
  _initTemperature(); // initialize temperature
  m_global_deltat.assign(dt);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "[ArcaneFem-Info] Started module compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop code after computations
  if (t >= tmax)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  // Set if we want to keep the matrix structure between calls.
  // The matrix has to have the same structure (same structure for non-zero)
  bool keep_struct = true;
  if (m_linear_system.isInitialized() && keep_struct) {
    m_linear_system.clearValues();
  }
  else {
    m_linear_system.reset();
    m_linear_system.setLinearSystemFactory(options()->linearSystem());
    m_linear_system.initialize(subDomain(), acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");
  }

  if (m_petsc_flags != NULL){
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  _doStationarySolve();
  _updateTime();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_initTime()
{
  info() << "[ArcaneFem-Info] Started module _initTime()";

  tmax = options()->tmax();
  dt = options()->dt();
  tmax = tmax;
  t = 0.0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateTime()
{
  info() << "[ArcaneFem-Info] Started module _updateTime()";

  t += dt;
  info() << "Time t is :" << t << " (s)";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateVariables()
{
  info() << "[ArcaneFem-Info] Started module _updateVariables()";
  Real elapsedTime = platform::getRealTime();

  {
    VariableDoFReal& dof_temperature(m_linear_system.solutionVariable());
    // Copy RHS DoF to Node temperature
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_temperature[node_dof.dofId(node, 0)];
      m_node_temperature[node] = v;
      m_node_temperature_old[node] = m_node_temperature[node];
    }

    m_node_temperature.synchronize();
    m_node_temperature_old.synchronize();
    if (mesh()->dimension() == 2)
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;

        Real3 grad = ArcaneFemFunctions::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_node_temperature);
        m_flux[cell].x = -m_cell_lambda[cell] * grad.x;
        m_flux[cell].y = -m_cell_lambda[cell] * grad.y;
        m_flux[cell].z = 0.;
      }

    if (mesh()->dimension() == 3)
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;

        Real3 grad = ArcaneFemFunctions::FeOperation3D::computeGradientTetra4(cell, m_node_coord, m_node_temperature);
        m_flux[cell].x = -m_cell_lambda[cell] * grad.x;
        m_flux[cell].y = -m_cell_lambda[cell] * grad.y;
        m_flux[cell].z = 0.;
      }

    m_flux.synchronize();
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "update-variables", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_initTemperature()
{
  info() << "[ArcaneFem-Info] Started module _initTemperature()";

  Tinit = options()->Tinit();
  m_node_temperature_old.fill(Tinit);
}

/*---------------------------------------------------------------------------*/
/*
 * @brief Solve a single time step of the FEM problem.
 *
 * This function performs the following steps: 
 * 1. Assembles the bilinear operator (if required).
 * 2. Assembles the linear operator (if required).
 * 3. Solves the linear system (if required).
 * 4. Updates the variables (temperature and flux).
 * 5. Validates the results (if cross-validation is enabled).
 *
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
  if (m_cross_validation && (t >= tmax)) {
    _validateResults();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getParameters()
{
  info() << "[ArcaneFem-Info] Started module _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  lambda = options()->lambda();
  qdot = options()->qdot();

  // Initialize cell lambda with default value
  m_cell_lambda.fill(lambda);

  for (const auto& bs : options()->materialProperty()) {
    CellGroup group = bs->volume();
    Real value = bs->lambda();
    info() << "Lambda for group=" << group.name() << " v=" << value;

    ENUMERATE_ (Cell, icell, group) {
      Cell cell = *icell;
      m_cell_lambda[cell] = value;
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*
  * @brief Assembles the linear operator for the FEM linear system.
  *
  * This method performs the following steps:
  *   1. Initializes the right-hand side (RHS) vector.
  *   2. Assembles the convection term for each boundary condition.
  *   3. Applies the variable source term to the RHS vector.
  *   4. Applies Neumann and Dirichlet boundary conditions.
  *   5. Applies point Dirichlet boundary conditions.
  */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR")
    m_bsr_format.toLinearSystem(m_linear_system);

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // Convective term ∫(𝒉 ⋅ 𝑇 ⋅ 𝑣ʰ ) for domain ∂Ωₐ
  for (const auto& bs : options()->convectionBoundaryCondition()) {
    FaceGroup group = bs->surface();
    h = bs->h();
    Text = bs->Text();

    auto processConvectionBoundaryCondition = [&](auto U, Real factor, auto computeMeasure, Integer num_nodes) {
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real measure = computeMeasure(face, m_node_coord);

        if (t <= dt - 1e-8) {
          auto K_e = factor * massMatrix(U, U) * measure;

          Int32 n1_index = 0;
          for (Node node1 : face.nodes()) {
            Int32 n2_index = 0;
            for (Node node2 : face.nodes()) {
              Real v = K_e(n1_index, n2_index);
              if (node1.isOwn())
                m_linear_system.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
              ++n2_index;
            }
            ++n1_index;
          }
        }

        for (Node node : iface->nodes()) {
          if (node.isOwn())
            rhs_values[node_dof.dofId(node, 0)] += h * Text * measure / num_nodes;
        }
      }
    };

    if (mesh()->dimension() == 2)
      processConvectionBoundaryCondition(RealVector<2>{1., 1.}, 1 / 6., ArcaneFemFunctions::MeshOperation::computeLengthEdge2, 2);
    else if (mesh()->dimension() == 3)
      processConvectionBoundaryCondition(RealVector<3>{1., 1., 1.}, 1 / 12., ArcaneFemFunctions::MeshOperation::computeAreaTria3, 3);
  }

  // Helper lambda to apply boundary conditions
  auto applyBoundaryConditions = [&](auto BCFunctions) {
    m_node_temperature_old.mult(1.0 / dt);
    BCFunctions.applyVariableSourceToRhs(m_node_temperature_old, mesh(), node_dof, m_node_coord, rhs_values);

    BC::IArcaneFemBC* bc = options()->boundaryConditions();
    if (bc) {
      for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions())
        BCFunctions.applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);

      for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions()) {
        //  BCFunctions.applyDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);
        FaceGroup face_group = bs->getSurface();
        NodeGroup node_group = face_group.nodeGroup();
        Real value = bs->getValue();
        if (bs->getEnforceDirichletMethod() == "Penalty") {
          Real penalty = bs->getPenalty();
          ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(value, penalty, node_dof, m_linear_system, rhs_values, node_group);
        }
        else if (bs->getEnforceDirichletMethod() == "RowElimination") {
          if (t <= dt - 1e-8)
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowElimination(value, node_dof, m_linear_system, rhs_values, node_group);
          else
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(value, node_dof, m_linear_system, rhs_values, node_group);
        }
        // else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
        //     ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowColumnElimination(value, node_dof, m_linear_system, rhs_values, node_group);
        // }
        else {
          ARCANE_FATAL("Unknown Dirichlet method");
        }
      }

      for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions()){
      //  BCFunctions.applyPointDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);
      NodeGroup node_group = bs->getNode();
      Real value = bs->getValue();
      if (bs->getEnforceDirichletMethod() == "Penalty") {
        Real penalty = bs->getPenalty();
        ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(value, penalty, node_dof, m_linear_system, rhs_values, node_group);
      }
      else if (bs->getEnforceDirichletMethod() == "RowElimination") {
        if (t <= dt - 1e-8)
          ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowElimination(value, node_dof, m_linear_system, rhs_values, node_group);
        else
          ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(value, node_dof, m_linear_system, rhs_values, node_group);
      }
      // else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
      //     ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowColumnElimination(value, node_dof, m_linear_system, rhs_values, node_group);
      // }
      else {
        ARCANE_FATAL("Unknown Dirichlet method");
      }
      }
    }
  };

  // Apply the correct boundary conditions based on mesh dimension
  if (mesh()->dimension() == 3) {
    using BCFunctions = ArcaneFemFunctions::BoundaryConditions3D;
    applyBoundaryConditions(BCFunctions());
  }
  else {
    using BCFunctions = ArcaneFemFunctions::BoundaryConditions2D;
    applyBoundaryConditions(BCFunctions());
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "rhs-vector-assembly", elapsedTime);
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

  if (t <= dt - 1e-8) {
    if (m_matrix_format == "BSR") {
      UnstructuredMeshConnectivityView m_connectivity_view(mesh());
      auto cn_cv = m_connectivity_view.cellNode();
      auto queue = subDomain()->acceleratorMng()->defaultQueue();
      auto command = makeCommand(queue);
      auto in_node_coord = ax::viewIn(command, m_node_coord);
      auto in_cell_lambda = ax::viewIn(command, m_cell_lambda);
      auto in_dt = dt;

      if (mesh()->dimension() == 2)
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return _computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord, in_cell_lambda, in_dt); });
      else
        m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return _computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord, in_cell_lambda, in_dt); });
    }
    else {
      if (mesh()->dimension() == 3)
        _assembleBilinear<4>([this](const Cell& cell) {
          return _computeElementMatrixTetra4(cell);
        });
      else
        _assembleBilinear<3>([this](const Cell& cell) {
          return _computeElementMatrixTria3(cell);
        });
    }
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
/*---------------------------------------------------------------------------*/

void FemModule::
_validateResults()
{
  info() << "[ArcaneFem-Info] Started module _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200)
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "T[" << node.uniqueId() << "] = " << m_node_temperature[node];
    }

  String filename = options()->resultFile();
  info() << "ValidateResultFile filename=" << filename;

  if (!filename.empty())
    checkNodeResultFile(traceMng(), filename, m_node_temperature, 1.0e-4);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
