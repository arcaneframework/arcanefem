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

  E = options()->E(); // Youngs modulus
  nu = options()->nu(); // Poission ratio

  mu2 = (E / (2 * (1 + nu))) * 2; // lame parameter mu * 2
  lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter lambda

  m_use_bsr = options()->bsr;

  elapsedTime = platform::getRealTime() - elapsedTime;
  _printArcaneFemTime("[ArcaneFem-Timer] get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
// Assemble the FEM linear operator
//  - The method also adds source term
//  - The method also adds external fluxes
//  - This function enforces a Dirichlet boundary condition in a weak sense
//    via the penalty method
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  //----------------------------------------------
  // Body force term assembly $int_{Omega}(f.v)$
  //----------------------------------------------

  // body force (fx, fy, fz) = (f[0], f[1], f[2])
  const UniqueArray<String> f_string = options()->f();
  info() << "[ArcaneFem-Info] Applying Bodyforce " << f_string;
  for (Int32 i = 0; i < f_string.size(); ++i) {
    f[i] = 0.0;
    if (f_string[i] != "NULL") {
      f[i] = std::stod(f_string[i].localstr());
    }
  }

  if (f_string[0] != "NULL" || f_string[1] != "NULL") {
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real area = _computeAreaTriangle3(cell);
      for (Node node : cell.nodes()) {
        if (node.isOwn()) {
          rhs_values[node_dof.dofId(node, 0)] += f[0] * area / 3;
          rhs_values[node_dof.dofId(node, 1)] += f[1] * area / 3;
        }
      }
    }
  }

  //----------------------------------------------
  // Traction term assembly  $int_{dOmega_N}((t.v)$
  //----------------------------------------------

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

    if (t_string[0] != "NULL" || t_string[1] != "NULL") {
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real length = _computeEdgeLength2(face);
        for (Node node : iface->nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += t[0] * length / 2.;
            rhs_values[node_dof.dofId(node, 1)] += t[1] * length / 2.;
          }
        }
      }
    }
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
                if (m_use_bsr)
                  m_bsr_format.matrix().setValue(dof_id, dof_id, Penalty);
                else
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

    if (options()->enforceDirichletMethod() == "Penalty") {

      info() << "[ArcaneFem-Info] Applying Dirichlet via "
             << options()->enforceDirichletMethod() << " method "
             << "with Dirichlet " << u_dirichlet_string;

      Real Penalty = options()->penalty();

      for (Int32 i = 0; i < u_dirichlet_string.size(); ++i) {
        if (u_dirichlet_string[i] != "NULL") {
          Real u_dirichlet = std::stod(u_dirichlet_string[i].localstr());
          ENUMERATE_ (Node, inode, group) {
            Node node = *inode;
            DoFLocalId dof_id = node_dof.dofId(node, i);
            if (node.isOwn()) {
              if (m_use_bsr)
                m_bsr_format.matrix().setValue(dof_id, dof_id, Penalty);
              else
                m_linear_system.matrixSetValue(dof_id, dof_id, Penalty);
              rhs_values[dof_id] = Penalty * u_dirichlet;
            }
          }
        }
      }
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
