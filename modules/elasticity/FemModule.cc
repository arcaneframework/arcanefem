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
  info() << "[ArcaneFem-Info] Started module  startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dof_per_node = defaultMesh()->dimension();
  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->crossValidation();
  m_petsc_flags = options()->petscFlags();

  m_dofs_on_nodes.initialize(defaultMesh(), m_dof_per_node);

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
 *   3. Sets PETSc flags if user has provided them.
 *   4. Executes the stationary solve.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
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

void FemModule::_initBsr()
{
  info() << "[ArcaneFem-Info] Started module  _initBsr()";
  Real elapsedTime = platform::getRealTime();

  bool use_csr_in_linearsystem = options()->linearSystem.serviceName() == "HypreLinearSystem";

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
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "[ArcaneFem-Info] Started module  _getMaterialParameters()";
  Real elapsedTime = platform::getRealTime();

  E = options()->E(); // Youngs modulus
  nu = options()->nu(); // Poission ratio ν

  mu = (E / (2 * (1 + nu))); // lame parameter μ
  lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter λ

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_assembleDirichletsGpu()
{
  info() << "[ArcaneFem-Info] Started module  _assembleLinearOperatorGpu()";

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  auto mesh_ptr = mesh();

  // Dirichlet conditions to LHS and RHS
  for (const auto& bs : options()->dirichletBoundaryCondition()) {
    FaceGroup group = bs->surface();
    const UniqueArray<String> u_dirichlet_string = bs->u();

    auto method = options()->enforceDirichletMethod();

    info() << "[ArcaneFem-Info] Applying Dirichlet " << u_dirichlet_string << " on Gpu";
    info() << "[ArcaneFem-Info] Dirichlet surface '" << bs->surface().name() << "'";
    info() << "[ArcaneFem-Info] Dirichlet method '" << method << "'";

    if (method == "Penalty")
      FemUtils::Gpu::BoundaryConditions::applyDirichletViaPenaltyVectorial(m_dofs_on_nodes, m_linear_system, mesh_ptr, queue, group, options()->penalty(), u_dirichlet_string);
    else if (method == "RowElimination")
      FemUtils::Gpu::BoundaryConditions::applyDirichletViaRowEliminationVectorial(m_dofs_on_nodes, m_linear_system, mesh_ptr, queue, group, u_dirichlet_string);
    else
      ARCANE_THROW(Arccore::NotImplementedException, "Method is not supported.");
  }

  for (const auto& bs : options()->dirichletPointCondition()) {
    NodeGroup group = bs->node();
    const UniqueArray<String> u_dirichlet_string = bs->u();

    auto method = options()->enforceDirichletMethod();

    info() << "[ArcaneFem-Info] Applying point Dirichlet " << u_dirichlet_string << " on Gpu";
    info() << "[ArcaneFem-Info] Dirichlet points '" << group.name() << "'";
    info() << "[ArcaneFem-Info] Dirichlet method '" << method << "'";

    if (method == "Penalty")
      FemUtils::Gpu::BoundaryConditions::applyPointDirichletViaPenaltyVectorial(m_dofs_on_nodes, m_linear_system, mesh(), queue, group, options()->penalty(), u_dirichlet_string);
    else if (method == "RowElimination")
      FemUtils::Gpu::BoundaryConditions::applyPointDirichletViaRowEliminationVectorial(m_dofs_on_nodes, m_linear_system, mesh_ptr, queue, group, u_dirichlet_string);
    else
      ARCANE_THROW(Arccore::NotImplementedException, "Method is not supported.");
  }
}

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

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  //----------------------------------------------------------------------
  // body force ∫∫∫ (𝐟.𝐯)  with 𝐟 = (𝑓𝑥, 𝑓𝑦, 𝑓𝑧) = (f[0], f[1], f[2])
  //----------------------------------------------------------------------
  const UniqueArray<String> f_string = options()->f();
  info() << "[ArcaneFem-Info] Applying Bodyforce " << f_string;
  for (Int32 i = 0; i < f_string.size(); ++i) {
    f[i] = 0.0;
    if (f_string[i] != "NULL") {
      f[i] = std::stod(f_string[i].localstr());
    }
  }

  if (mesh()->dimension() == 2)
    if (f_string[0] != "NULL" || f_string[1] != "NULL")
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
        for (Node node : cell.nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += f[0] * area / 3;
            rhs_values[node_dof.dofId(node, 1)] += f[1] * area / 3;
          }
        }
      }
  if (mesh()->dimension() == 3)
    if (f_string[0] != "NULL" || f_string[1] != "NULL" || f_string[1] != "NULL")
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;
        Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);
        for (Node node : cell.nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += f[0] * volume / 4;
            rhs_values[node_dof.dofId(node, 1)] += f[1] * volume / 4;
            rhs_values[node_dof.dofId(node, 2)] += f[2] * volume / 4;
          }
        }
      }

  //----------------------------------------------------------------------
  // traction term ∫∫ (𝐭.𝐯)  with 𝐭 = (𝑡𝑥, 𝑡𝑦, 𝑡𝑧) = (t[0], t[1], t[2])
  //----------------------------------------------------------------------
  for (const auto& bs : options()->tractionBoundaryCondition()) {
    FaceGroup group = bs->surface();
    const UniqueArray<String> t_string = bs->t();

    info() << "[ArcaneFem-Info] Applying Traction " << t_string;
    info() << "[ArcaneFem-Info] Traction surface '" << bs->surface().name() << "'";

    for (Int32 i = 0; i < t_string.size(); ++i) {
      t[i] = 0.0;
      if (t_string[i] != "NULL") {
        t[i] = std::stod(t_string[i].localstr());
      }
    }

    if (mesh()->dimension() == 2)
      if (t_string[0] != "NULL" || t_string[1] != "NULL")
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);
          for (Node node : iface->nodes()) {
            if (node.isOwn()) {
              rhs_values[node_dof.dofId(node, 0)] += t[0] * length / 2.;
              rhs_values[node_dof.dofId(node, 1)] += t[1] * length / 2.;
            }
          }
        }

    if (mesh()->dimension() == 3)
      if (t_string[0] != "NULL" || t_string[1] != "NULL" || t_string[2] != "NULL")
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, m_node_coord);
          for (Node node : iface->nodes()) {
            if (node.isOwn()) {
              rhs_values[node_dof.dofId(node, 0)] += t[0] * area / 3.;
              rhs_values[node_dof.dofId(node, 1)] += t[1] * area / 3.;
              rhs_values[node_dof.dofId(node, 2)] += t[2] * area / 3.;
            }
          }
        }
  }

  auto use_hypre = options()->linearSystem.serviceName() == "HypreLinearSystem";
  if (use_hypre) {
    // The rest of the assembly can be handled on Gpu because of Hypre solver.
    _assembleDirichletsGpu();
    return;
  }

  //----------------------------------------------
  // Dirichlet conditions to LHS and RHS
  //----------------------------------------------

  for (const auto& bs : options()->dirichletBoundaryCondition()) {
    FaceGroup group = bs->surface();
    const UniqueArray<String> u_dirichlet_string = bs->u();

    info() << "[ArcaneFem-Info] Applying Dirichlet " << u_dirichlet_string;
    info() << "[ArcaneFem-Info] Dirichlet surface '" << bs->surface().name() << "'";
    info() << "[ArcaneFem-Info] Dirichlet method '" << options()->enforceDirichletMethod() << "'";

    if (options()->enforceDirichletMethod() == "Penalty") {

      Real Penalty = options()->penalty();

      for (Int32 i = 0; i < u_dirichlet_string.size(); ++i) {
        if (u_dirichlet_string[i] != "NULL") {
          Real u_dirichlet = std::stod(u_dirichlet_string[i].localstr());
          ENUMERATE_ (Face, iface, group) {
            for (Node node : iface->nodes()) {
              DoFLocalId dof_id = node_dof.dofId(node, i);
              if (node.isOwn()) {
                m_linear_system.matrixSetValue(dof_id, dof_id, Penalty);
                rhs_values[dof_id] = Penalty * u_dirichlet;
              }
            }
          }
        }
      }
    }
    else if (options()->enforceDirichletMethod() == "RowElimination") {

      for (Int32 i = 0; i < u_dirichlet_string.size(); ++i) {
        if (u_dirichlet_string[i] != "NULL") {
          Real u_dirichlet = std::stod(u_dirichlet_string[i].localstr());
          ENUMERATE_ (Face, iface, group) {
            for (Node node : iface->nodes()) {
              DoFLocalId dof_id = node_dof.dofId(node, i);
              if (node.isOwn()) {
                m_linear_system.eliminateRow(dof_id, u_dirichlet);
              }
            }
          }
        }
      }
    }
    else if (options()->enforceDirichletMethod() == "RowColumnElimination") {

      for (Int32 i = 0; i < u_dirichlet_string.size(); ++i) {
        if (u_dirichlet_string[i] != "NULL") {
          Real u_dirichlet = std::stod(u_dirichlet_string[i].localstr());
          ENUMERATE_ (Face, iface, group) {
            for (Node node : iface->nodes()) {
              DoFLocalId dof_id = node_dof.dofId(node, i);
              if (node.isOwn()) {
                m_linear_system.eliminateRowColumn(dof_id, u_dirichlet);
              }
            }
          }
        }
      }
    }
  }

  for (const auto& bs : options()->dirichletPointCondition()) {
    NodeGroup group = bs->node();
    const UniqueArray<String> u_dirichlet_string = bs->u();

    info() << "[ArcaneFem-Info] Applying point Dirichlet " << u_dirichlet_string;
    info() << "[ArcaneFem-Info] Dirichlet points '" << group.name() << "'";
    info() << "[ArcaneFem-Info] Dirichlet method '" << options()->enforceDirichletMethod() << "'";

    if (options()->enforceDirichletMethod() == "Penalty") {
      Real Penalty = options()->penalty();

      for (Int32 i = 0; i < u_dirichlet_string.size(); ++i) {
        if (u_dirichlet_string[i] != "NULL") {
          Real u_dirichlet = std::stod(u_dirichlet_string[i].localstr());
          ENUMERATE_ (Node, inode, group) {
            Node node = *inode;
            DoFLocalId dof_id = node_dof.dofId(node, i);
            if (node.isOwn()) {
              m_linear_system.matrixSetValue(dof_id, dof_id, Penalty);
              rhs_values[dof_id] = Penalty * u_dirichlet;
            }
          }
        }
      }
    }
    else if (options()->enforceDirichletMethod() == "RowElimination") {
      for (Int32 i = 0; i < u_dirichlet_string.size(); ++i) {
        if (u_dirichlet_string[i] != "NULL") {
          Real u_dirichlet = std::stod(u_dirichlet_string[i].localstr());
          ENUMERATE_ (Node, inode, group) {
            Node node = *inode;
            DoFLocalId dof_id = node_dof.dofId(node, i);
            if (node.isOwn()) {
              m_linear_system.eliminateRow(dof_id, u_dirichlet);
            }
          }
        }
      }
    }
    else if (options()->enforceDirichletMethod() == "RowColumnElimination") {
      for (Int32 i = 0; i < u_dirichlet_string.size(); ++i) {
        if (u_dirichlet_string[i] != "NULL") {
          Real u_dirichlet = std::stod(u_dirichlet_string[i].localstr());
          ENUMERATE_ (Node, inode, group) {
            Node node = *inode;
            DoFLocalId dof_id = node_dof.dofId(node, i);
            if (node.isOwn()) {
              m_linear_system.eliminateRowColumn(dof_id, u_dirichlet);
            }
          }
        }
      }
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly given as mesh type.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
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
    if (mesh()->dimension() == 2)
      _assembleBilinearOperatorTria3();
    if (mesh()->dimension() == 3)
    _assembleBilinearOperatorTetra4();
  }
  else {
    ARCANE_FATAL("Unsupported matrix type, only DOK| BSR|AF-BSR is supported.");
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"lhs-matrix-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTetra4()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    auto K_e = _computeElementMatrixTetra4(cell);
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        Real v1 = K_e(3 * n1_index, 3 * n2_index);
        Real v2 = K_e(3 * n1_index, 3 * n2_index + 1);
        Real v3 = K_e(3 * n1_index, 3 * n2_index + 2);

        Real v4 = K_e(3 * n1_index + 1, 3 * n2_index);
        Real v5 = K_e(3 * n1_index + 1, 3 * n2_index + 1);
        Real v6 = K_e(3 * n1_index + 1, 3 * n2_index + 2);

        Real v7 = K_e(3 * n1_index + 2, 3 * n2_index);
        Real v8 = K_e(3 * n1_index + 2, 3 * n2_index + 1);
        Real v9 = K_e(3 * n1_index + 2, 3 * n2_index + 2);
        if (node1.isOwn()) {
          DoFLocalId node1_dof1 = node_dof.dofId(node1, 0);
          DoFLocalId node1_dof2 = node_dof.dofId(node1, 1);
          DoFLocalId node1_dof3 = node_dof.dofId(node1, 2);
          DoFLocalId node2_dof1 = node_dof.dofId(node2, 0);
          DoFLocalId node2_dof2 = node_dof.dofId(node2, 1);
          DoFLocalId node2_dof3 = node_dof.dofId(node2, 2);

          m_linear_system.matrixAddValue(node1_dof1, node2_dof1, v1);
          m_linear_system.matrixAddValue(node1_dof1, node2_dof2, v2);
          m_linear_system.matrixAddValue(node1_dof1, node2_dof3, v3);

          m_linear_system.matrixAddValue(node1_dof2, node2_dof1, v4);
          m_linear_system.matrixAddValue(node1_dof2, node2_dof2, v5);
          m_linear_system.matrixAddValue(node1_dof2, node2_dof3, v6);

          m_linear_system.matrixAddValue(node1_dof3, node2_dof1, v7);
          m_linear_system.matrixAddValue(node1_dof3, node2_dof2, v8);
          m_linear_system.matrixAddValue(node1_dof3, node2_dof3, v9);
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
_assembleBilinearOperatorTria3()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    if (cell.type() != IT_Triangle3)
      ARCANE_FATAL("Only Triangle3 cell type is supported");

    auto K_e = _computeElementMatrixTria3(cell); // element stiffness matrix
    // assemble elementary matrix into the global one elementary terms are
    // positioned into K according to the rank  of  associated  node in the
    // mesh.nodes list and according the dof  number. Here  for  each  node
    // two dofs exists [u1,u2].  For each TRIA3 there are 3 nodes hence the
    // elementary stiffness matrix size is (3*2 x 3*2)=(6x6). We will  fill
    // this below in 4 at a time.
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        Real v1 = K_e(2 * n1_index, 2 * n2_index);
        Real v2 = K_e(2 * n1_index, 2 * n2_index + 1);
        Real v3 = K_e(2 * n1_index + 1, 2 * n2_index);
        Real v4 = K_e(2 * n1_index + 1, 2 * n2_index + 1);

        if (node1.isOwn()) {
          DoFLocalId node1_dof1 = node_dof.dofId(node1, 0);
          DoFLocalId node1_dof2 = node_dof.dofId(node1, 1);
          DoFLocalId node2_dof1 = node_dof.dofId(node2, 0);
          DoFLocalId node2_dof2 = node_dof.dofId(node2, 1);

          m_linear_system.matrixAddValue(node1_dof1, node2_dof1, v1);
          m_linear_system.matrixAddValue(node1_dof1, node2_dof2, v2);
          m_linear_system.matrixAddValue(node1_dof2, node2_dof1, v3);
          m_linear_system.matrixAddValue(node1_dof2, node2_dof2, v4);
        }
        ++n2_index;
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

void FemModule::
_solve()
{
  info() << "[ArcaneFem-Info] Started module  _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.solve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"solve-linear-system", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
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

  String filename = options()->resultFile();
  const double epsilon = 1.0e-3;
  const double min_value_to_test = 1.0e-16;

  info() << "[ArcaneFem-Info] Validating results filename=" << filename << " epsilon =" << epsilon;

  if (!filename.empty())
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

void FemModule::
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

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
