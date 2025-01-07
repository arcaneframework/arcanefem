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

  m_dofs_on_nodes.initialize(defaultMesh(), 2);

  _initBoundaryconditions();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] initialize", elapsedTime);
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
  info() << "[ArcaneFem-Info] Started module  compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  _doStationarySolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Initilizes BSR matrix.
 */
/*---------------------------------------------------------------------------*/

void FemModule::_initBsr()
{
  info() << "[ArcaneFem-Info] Started module  compute()";
  Real elapsedTime = platform::getRealTime();

  bool use_csr_in_linearsystem = options()->linearSystem.serviceName() == "HypreLinearSystem";
  m_bsr_format.initialize(defaultMesh(), nbFace(), use_csr_in_linearsystem);
  m_bsr_format.computeSparsity();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows a sequence of steps to solve FEM system:
 *
 *   1. _getMaterialParameters()     Retrieves material parameters via
 *   2. _assembleBilinearOperator()  Assembles the FEM  matrix A
 *   3. _assembleLinearOperator()    Assembles the FEM RHS vector b
 *   4. _solve()                     Solves for solution vector u = A^-1*b
 *   5. _updateVariables()           Updates FEM variables u = x
 *   6. _validateResults()           Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  _getMaterialParameters();
  _assembleBilinearOperator();
  _assembleLinearOperator();
  _solve();
  _updateVariables();
  _validateResults();
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

  f1 = options()->f1(); // body force in Y
  f2 = options()->f2(); // body force in Y
  E = options()->E(); // Youngs modulus
  nu = options()->nu(); // Poission ratio

  mu2 = (E / (2 * (1 + nu))) * 2; // lame parameter mu * 2
  lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter lambda

  m_use_bsr = options()->bsr;

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// TODO : To be removed

void FemModule::
_initBoundaryconditions()
{
  info() << "[ArcaneFem-Info] Started module  _initBoundaryconditions()";
  Real elapsedTime = platform::getRealTime();

  _applyDirichletBoundaryConditions();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] applying-boundary-cond", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_applyDirichletBoundaryConditions()
{
  for (const auto& bs : options()->dirichletBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real u1_val = bs->u1();
    Real u2_val = bs->u2();

    if (bs->u1.isPresent() && bs->u2.isPresent()) {
      info() << "[ArcaneFem-Info] Applying Dirichlet on group = " << group.name() << " u1= " << u1_val << " u2= " << u2_val;
      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          m_U[node].x = u1_val;
          m_U[node].y = u2_val;
          m_u1_fixed[node] = true;
          m_u2_fixed[node] = true;
        }
      }
      continue;
    }

    if (bs->u1.isPresent()) {
      info() << "[ArcaneFem-Info] Applying Dirichlet on group = " << group.name() << " u1=" << u1_val;
      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          m_U[node].x = u1_val;
          m_u1_fixed[node] = true;
        }
      }
      continue;
    }

    if (bs->u2.isPresent()) {
      info() << "[ArcaneFem-Info] Applying Dirichlet on group = " << group.name() << " u2=" << u2_val;
      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          m_U[node].y = u2_val;
          m_u2_fixed[node] = true;
        }
      }
      continue;
    }
  }

  for (const auto& bs : options()->dirichletPointCondition()) {
    NodeGroup group = bs->node();
    Real u1_val = bs->u1();
    Real u2_val = bs->u2();

    if (bs->u1.isPresent() && bs->u2.isPresent()) {
      info() << "[ArcaneFem-Info] Applying point Dirichlet on group = " << group.name() << " u1= " << u1_val << " u2= " << u2_val;
      ENUMERATE_ (Node, inode, group) {
        Node node = *inode;
        m_U[node].x = u1_val;
        m_U[node].y = u2_val;
        m_u1_fixed[node] = true;
        m_u2_fixed[node] = true;
      }
      continue;
    }

    if (bs->u1.isPresent()) {
      info() << "[ArcaneFem-Info] Applying point Dirichlet on group = " << group.name() << " u1=" << u1_val;
      ENUMERATE_ (Node, inode, group) {
        Node node = *inode;
        m_U[node].x = u1_val;
        m_u1_fixed[node] = true;
      }
      continue;
    }

    if (bs->u2.isPresent()) {
      info() << "[ArcaneFem-Info] Applying point Dirichlet on group = " << group.name() << " u2=" << u2_val;
      ENUMERATE_ (Node, inode, group) {
        Node node = *inode;
        m_U[node].y = u2_val;
        m_u2_fixed[node] = true;
      }
      continue;
    }
  }
}

/*---------------------------------------------------------------------------*/
// Assemble the FEM linear operator
//  - This function enforces a Dirichlet boundary condition in a weak sense
//    via the penalty method
//  - The method also adds source term
//  - The method also adds external fluxes
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (options()->enforceDirichletMethod() == "Penalty") {

    //----------------------------------------------
    // penalty method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'P' be the penalty term and let 'i' be the set of DOF for which
    //  Dirichlet condition needs to be applied
    //
    //  - For LHS matrix A the diag term corresponding to the Dirichlet DOF
    //           a_{i,i} = 1. * P
    //
    //  - For RHS vector b the term that corresponds to the Dirichlet DOF
    //           b_{i} = b_{i} * P
    //----------------------------------------------

    info() << "[ArcaneFem-Info] Applying Dirichlet via "
           << options()->enforceDirichletMethod() << " method ";

    Real Penalty = options()->penalty(); // 1.0e30 is the default

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u1_fixed[node_id]) {
        DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);
        if (m_use_bsr)
          m_bsr_format.matrix().setValue(dof_id1, dof_id1, Penalty);
        else
          m_linear_system.matrixSetValue(dof_id1, dof_id1, Penalty);
        {
          Real u1_dirichlet = Penalty * m_U[node_id].x;
          rhs_values[dof_id1] = u1_dirichlet;
        }
      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);
        if (m_use_bsr)
          m_bsr_format.matrix().setValue(dof_id2, dof_id2, Penalty);
        else
          m_linear_system.matrixSetValue(dof_id2, dof_id2, Penalty);
        {
          Real u2_dirichlet = Penalty * m_U[node_id].y;
          rhs_values[dof_id2] = u2_dirichlet;
        }
      }
    }
  }
  else if (options()->enforceDirichletMethod() == "WeakPenalty") {

    //----------------------------------------------
    // weak penalty method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'P' be the penalty term and let 'i' be the set of DOF for which
    //  Dirichlet condition needs to be applied
    //
    //  - For LHS matrix A the diag term corresponding to the Dirichlet DOF
    //           a_{i,i} = a_{i,i} + P
    //
    //  - For RHS vector b the term that corresponds to the Dirichlet DOF
    //           b_{i} = b_{i} * P
    //----------------------------------------------

    info() << "[ArcaneFem-Info] Applying Dirichlet via "
           << options()->enforceDirichletMethod() << " method ";

    Real Penalty = options()->penalty(); // 1.0e30 is the default

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u1_fixed[node_id]) {
        DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);
        m_linear_system.matrixAddValue(dof_id1, dof_id1, Penalty);
        {
          Real u1_dirichlet = Penalty * m_U[node_id].x;
          rhs_values[dof_id1] = u1_dirichlet;
        }
      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);
        m_linear_system.matrixAddValue(dof_id2, dof_id2, Penalty);
        {
          Real u2_dirichlet = Penalty * m_U[node_id].y;
          rhs_values[dof_id2] = u2_dirichlet;
        }
      }
    }
  }
  else if (options()->enforceDirichletMethod() == "RowElimination") {

    //----------------------------------------------
    // Row elimination method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'i' be the DOF for which  Dirichlet condition 'g_i' needs to be applied
    //
    //  to apply the Dirichlet on 'i'th DOF
    //  - For LHS matrix A the row terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j
    //           a_{i,j} = 1.  : i==j
    //  - For RHS vector b the terms corresponding to the Dirichlet DOF
    //           b_i = g_i
    //----------------------------------------------

    info() << "[ArcaneFem-Info] Applying Dirichlet via "
           << options()->enforceDirichletMethod() << " method ";

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u1_fixed[node_id]) {
        DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);

        Real u1_dirichlet = m_U[node_id].x;
        m_linear_system.eliminateRow(dof_id1, u1_dirichlet);
      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);

        Real u2_dirichlet = m_U[node_id].y;
        m_linear_system.eliminateRow(dof_id2, u2_dirichlet);
      }
    }
  }
  else if (options()->enforceDirichletMethod() == "RowColumnElimination") {

    //----------------------------------------------
    // Row elimination method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'I' be the set of DOF for which  Dirichlet condition needs to be applied
    //
    //  to apply the Dirichlet on 'i'th DOF
    //  - For LHS matrix A the row terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j  for all j
    //           a_{i,j} = 1.  : i==j
    //    also the column terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j  for all i
    //----------------------------------------------

    info() << "[ArcaneFem-Info] Applying Dirichlet via "
           << options()->enforceDirichletMethod() << " method ";

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u1_fixed[node_id]) {
        DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);

        Real u1_dirichlet = m_U[node_id].x;
        m_linear_system.eliminateRowColumn(dof_id1, u1_dirichlet);
      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);

        Real u2_dirichlet = m_U[node_id].y;
        m_linear_system.eliminateRowColumn(dof_id2, u2_dirichlet);
      }
    }
  }
  else {

    info() << "[ArcaneFem-Info] Applying Dirichlet via "
           << options()->enforceDirichletMethod() << " is not supported \n"
           << "enforce-Dirichlet-method only supports:\n"
           << "  - Penalty\n"
           << "  - WeakPenalty\n"
           << "  - RowElimination\n"
           << "  - RowColumnElimination\n";
  }

  //----------------------------------------------
  // Body force term assembly
  //----------------------------------------------
  //
  //  $int_{Omega}(f1*v1^h)$
  //  $int_{Omega}(f2*v2^h)$
  //  only for noded that are non-Dirichlet
  //----------------------------------------------

  if (options()->f1.isPresent()) {
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real area = _computeAreaTriangle3(cell);
      for (Node node : cell.nodes()) {
        if (!(m_u1_fixed[node]) && node.isOwn()) {
          DoFLocalId dof_id1 = node_dof.dofId(node, 0);
          rhs_values[dof_id1] += f1 * area / 3;
        }
      }
    }
  }

  if (options()->f2.isPresent()) {
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real area = _computeAreaTriangle3(cell);
      for (Node node : cell.nodes()) {
        if (!(m_u2_fixed[node]) && node.isOwn()) {
          DoFLocalId dof_id2 = node_dof.dofId(node, 1);
          rhs_values[dof_id2] += f2 * area / 3;
        }
      }
    }
  }

  //----------------------------------------------
  // Traction term assembly
  //----------------------------------------------
  //
  //  $int_{dOmega_N}((tx.nx)*v1^h)$
  //  $int_{dOmega_N}((ty.ny)*v1^h)$
  //  only for noded that are non-Dirichlet
  //----------------------------------------------
  for (const auto& bs : options()->tractionBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real t1_val = bs->t1();
    Real t2_val = bs->t2();

    if (bs->t1.isPresent() && bs->t2.isPresent()) {
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real length = _computeEdgeLength2(face);
        for (Node node : iface->nodes()) {
          if (!(m_u1_fixed[node]) && node.isOwn()) {
            DoFLocalId dof_id1 = node_dof.dofId(node, 0);
            rhs_values[dof_id1] += t1_val * length / 2.;
          }
          if (!(m_u2_fixed[node]) && node.isOwn()) {
            DoFLocalId dof_id2 = node_dof.dofId(node, 1);
            rhs_values[dof_id2] += t2_val * length / 2.;
          }
        }
      }
      continue;
    }

    if (bs->t1.isPresent()) {
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real length = _computeEdgeLength2(face);
        for (Node node : iface->nodes()) {
          if (!(m_u1_fixed[node]) && node.isOwn()) {
            DoFLocalId dof_id1 = node_dof.dofId(node, 0);
            rhs_values[dof_id1] += t1_val * length / 2.;
          }
        }
      }
      continue;
    }

    if (bs->t2.isPresent()) {
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real length = _computeEdgeLength2(face);
        for (Node node : iface->nodes()) {
          if (!(m_u2_fixed[node]) && node.isOwn()) {
            DoFLocalId dof_id2 = node_dof.dofId(node, 1);
            rhs_values[dof_id2] += t2_val * length / 2.;
          }
        }
      }
      continue;
    }
  }

  if (m_use_bsr)
    m_bsr_format.toLinearSystem(m_linear_system);

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeAreaTriangle3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  return 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeEdgeLength2(Face face)
{
  Real3 m0 = m_node_coord[face.nodeId(0)];
  Real3 m1 = m_node_coord[face.nodeId(1)];
  return math::sqrt((m1.x - m0.x) * (m1.x - m0.x) + (m1.y - m0.y) * (m1.y - m0.y));
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

  if (m_use_bsr) {
    _initBsr();

    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(acceleratorMng()->defaultQueue());
    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto lambda_copy = lambda;
    auto mu2_copy = mu2;

    m_bsr_format.assembleCellWise([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTRIA3Gpu(cell_lid, cn_cv, in_node_coord, lambda_copy, mu2_copy); });
  }
  else {
    _assembleBilinearOperatorTRIA3();
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] lhs-matrix-assembly", elapsedTime);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTRIA3()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    if (cell.type() != IT_Triangle3)
      ARCANE_FATAL("Only Triangle3 cell type is supported");

    auto K_e = _computeElementMatrixTRIA3(cell); // element stiffness matrix
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
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (node1.isOwn()) {
          DoFLocalId node1_dof1 = node_dof.dofId(node1, 0);
          DoFLocalId node1_dof2 = node_dof.dofId(node1, 1);
          DoFLocalId node2_dof1 = node_dof.dofId(node2, 0);
          DoFLocalId node2_dof2 = node_dof.dofId(node2, 1);
          //          m_linear_system.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
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
  _printArcaneFemTime("[ArcaneFem-Timer] solve-linear-system", elapsedTime);
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
  const double epsilon = 1.0e-4;
  const double min_value_to_test = 1.0e-16;

  info() << "[ArcaneFem-Info] Validating results filename=" << filename << " epsilon =" << epsilon;

  if (!filename.empty())
    Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_U, epsilon, min_value_to_test);

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] result-validation", elapsedTime);
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

  // Re-Apply boundary conditions because the solver has modified the value
  _applyDirichletBoundaryConditions();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real u1_val = dof_u[node_dof.dofId(node, 0)];
      Real u2_val = dof_u[node_dof.dofId(node, 1)];
      m_U[node] = Real3(u1_val, u2_val, 0.0);
    }
  }

  m_U.synchronize();

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] update-variables", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Function to prints the execution time `value` of phase `label`
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_printArcaneFemTime(const String label, const Real value)
{
  info() << std::left << std::setw(40) << label << " = " << value;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
