// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* FEM code to test vectorial FE for Soildynamics problem.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include "ElementMatrix.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "[ArcaneFem-Info] Started module  compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop code after computations
  if (t >= tmax)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  info() << "Time iteration at t : " << t << " (s) ";

  // Set if we want to keep the matrix structure between calls.
  // The matrix has to have the same structure (same structure for non-zero)
  bool keep_struct = true;
  if (m_linear_system.isInitialized() && keep_struct) {
    m_linear_system.clearValues();
  }
  else {
    m_linear_system.reset();
    m_linear_system.setLinearSystemFactory(options()->linearSystem());
    m_linear_system.initialize(subDomain(),  acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");
  }

  if (m_petsc_flags != NULL){
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  _doStationarySolve();
  if (t >= tmax  && m_cross_validation)
    _validateResults();
  _updateTime();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "[ArcaneFem-Info] Started module  startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dofs_on_nodes.initialize(mesh(), mesh()->dimension());

  _getParameters();

  bool use_csr_in_linearsystem = options()->linearSystem.serviceName() == "HypreLinearSystem";
  if (m_matrix_format == "BSR")
    m_bsr_format.initialize(defaultMesh(), mesh()->dimension(), use_csr_in_linearsystem, 0);
  else if (m_matrix_format == "AF-BSR")
    m_bsr_format.initialize(defaultMesh(), mesh()->dimension(), use_csr_in_linearsystem, 1);

  t = dt;
  tmax = tmax;
  m_global_deltat.assign(dt);

  _readCaseTables();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateTime()
{
  t += dt;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getParameters()
{
  info() << "[ArcaneFem-Info] Started module  _getParameters()";
  Real elapsedTime = platform::getRealTime();

  //--------- time parameters -----------//
  tmax = options()->tmax(); // max time
  dt = options()->dt(); // time step

  //--------- material parameter ---------//
  E = options()->E(); // Youngs modulus
  nu = options()->nu(); // Poission ratio
  rho = options()->rho(); // Density
  cp = options()->cp(); // Wave velocity primary
  cs = options()->cs(); // Wave velocity secondary

  if (options()->E.isPresent() && options()->nu.isPresent()) {
    mu = E / (2 * (1 + nu));
    lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    cs = math::sqrt(mu / rho);
    cp = math::sqrt((lambda + (2. * mu)) / rho);
  }

  if (options()->mu.isPresent() && options()->lambda.isPresent()) {
    mu = options()->mu;
    lambda = options()->lambda;
    cs = math::sqrt(mu / rho);
    cp = math::sqrt((lambda + (2. * mu)) / rho);
  }

  if ((options()->cp.isPresent()) && (options()->cs.isPresent())) {
    mu = cs * cs * rho;
    lambda = cp * cp * rho - 2 * mu;
  }

  if (options()->f.isPresent()) {
    const UniqueArray<String> f_string = options()->f();
    info() << "[ArcaneFem-Info] Applying Bodyforce " << f_string;
    for (Int32 i = 0; i < f_string.size(); ++i) {
      if (f_string[i] != "NULL") {
        f[i] = std::stod(f_string[i].localstr());
      }
    }
  }

  gamma = 0.5;
  beta = (1. / 4.) * (gamma + 0.5) * (gamma + 0.5);

  c0 = rho / (beta * dt * dt);
  c1 = lambda;
  c2 = mu;
  c3 = rho / (beta * dt);
  c4 = rho * (1. / 2. / beta - 1.);
  c5 = 0.;
  c6 = 0.;
  c7 = rho * gamma / beta / dt;
  c8 = rho * (1. - gamma / beta);
  c9 = rho * dt * (1. - gamma / (2. * beta));

  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->crossValidation();
  m_petsc_flags = options()->petscFlags();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_readCaseTables()
{
  IParallelMng* pm = subDomain()->parallelMng();

  for (const auto& bs : options()->tractionBoundaryCondition()) {
    CaseTable* case_table = nullptr;

    if (bs->tractionInputFile.isPresent())
      case_table = readFileAsCaseTable(pm, bs->tractionInputFile(), 3);

    m_traction_case_table_list.add(CaseTableInfo{ bs->tractionInputFile(), case_table });
  }

  for (const auto& bs : options()->doubleCouple()) {
    CaseTable* case_table = nullptr;

    if (bs->doubleCoupleInputFile.isPresent())
      case_table = readFileAsCaseTable(pm, bs->doubleCoupleInputFile(), 1);

    m_double_couple_case_table_list.add(CaseTableInfo{ bs->doubleCoupleInputFile(), case_table });
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateVariables()
{

  VariableDoFReal& dof_u(m_linear_system.solutionVariable());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  if (mesh()->dimension() == 2)
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real u1_val = dof_u[node_dof.dofId(node, 0)];
      Real u2_val = dof_u[node_dof.dofId(node, 1)];
      Real3 u_disp = {u1_val, u2_val, 0.};
      m_dU[node] = u_disp;
    }
  if (mesh()->dimension() == 3)
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real u1_val = dof_u[node_dof.dofId(node, 0)];
      Real u2_val = dof_u[node_dof.dofId(node, 1)];
      Real u3_val = dof_u[node_dof.dofId(node, 2)];
      Real3 u_disp = {u1_val, u2_val, u3_val};
      m_dU[node] = u_disp;
    }

  m_dU.synchronize();
  m_U.synchronize();
  m_V.synchronize();
  m_A.synchronize();

  // Note at this stage we already have calculated dU
  Real alocX;
  Real alocY;
  Real alocZ;

  if (mesh()->dimension() == 2)
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;

      alocX = (m_dU[node].x - m_U[node].x - dt * m_V[node].x) / beta / (dt * dt) - (1. - 2. * beta) / 2. / beta * m_A[node].x;
      alocY = (m_dU[node].y - m_U[node].y - dt * m_V[node].y) / beta / (dt * dt) - (1. - 2. * beta) / 2. / beta * m_A[node].y;

      m_V[node].x = m_V[node].x + dt * ((1. - gamma) * m_A[node].x + gamma * alocX);
      m_V[node].y = m_V[node].y + dt * ((1. - gamma) * m_A[node].y + gamma * alocY);

      m_A[node].x = alocX;
      m_A[node].y = alocY;

      m_U[node].x = m_dU[node].x;
      m_U[node].y = m_dU[node].y;
    }
  if (mesh()->dimension() == 3)
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;

      alocX = (m_dU[node].x - m_U[node].x - dt * m_V[node].x) / beta / (dt * dt) - (1. - 2. * beta) / 2. / beta * m_A[node].x;
      alocY = (m_dU[node].y - m_U[node].y - dt * m_V[node].y) / beta / (dt * dt) - (1. - 2. * beta) / 2. / beta * m_A[node].y;
      alocZ = (m_dU[node].z - m_U[node].z - dt * m_V[node].z) / beta / (dt * dt) - (1. - 2. * beta) / 2. / beta * m_A[node].z;

      m_V[node].x = m_V[node].x + dt * ((1. - gamma) * m_A[node].x + gamma * alocX);
      m_V[node].y = m_V[node].y + dt * ((1. - gamma) * m_A[node].y + gamma * alocY);
      m_V[node].z = m_V[node].z + dt * ((1. - gamma) * m_A[node].z + gamma * alocZ);

      m_A[node].x = alocX;
      m_A[node].y = alocY;
      m_A[node].z = alocZ;

      m_U[node].x = m_dU[node].x;
      m_U[node].y = m_dU[node].y;
      m_U[node].z = m_dU[node].z;
    }
}

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (mesh()->dimension() == 2)
    if (m_matrix_format == "BSR")
      _assembleLinearOperator2d(&(m_bsr_format.matrix()));
    else
      _assembleLinearOperator2d();
  if (mesh()->dimension() == 3)
    _assembleLinearOperator3d();

  if(m_matrix_format == "BSR")
    m_bsr_format.toLinearSystem(m_linear_system);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
// Assemble the FEM linear operator
//  - This function enforces a Dirichlet boundary condition in a weak sense
//    via the penalty method
//  - The method also adds source term
//  - The method also adds external fluxes
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator2d(BSRMatrix* bsr_matrix)
{
  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
    Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
    Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

    RealVector<6> Uy = { 0., 1., 0., 1., 0., 1. };
    RealVector<6> Ux = { 1., 0., 1., 0., 1., 0. };

    RealVector<6> F = { f[0], f[1], f[0], f[1], f[0], f[1] };

    RealVector<6> Un = { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y,
                             m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y,
                             m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y };
    RealVector<6> Vn = { m_V[cell.nodeId(0)].x, m_V[cell.nodeId(0)].y,
                             m_V[cell.nodeId(1)].x, m_V[cell.nodeId(1)].y,
                             m_V[cell.nodeId(2)].x, m_V[cell.nodeId(2)].y };
    RealVector<6> An = { m_A[cell.nodeId(0)].x, m_A[cell.nodeId(0)].y,
                             m_A[cell.nodeId(1)].x, m_A[cell.nodeId(1)].y,
                             m_A[cell.nodeId(2)].x, m_A[cell.nodeId(2)].y };
    //----------------------------------------------------------------------
    //  ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯)
    //----------------------------------------------------------------------
    RealVector<6> rhs = (   F * (1 / 3.)
                              + Un * (massMatrix(Ux, Ux) + massMatrix(Uy, Uy)) * (c0 * 1 / 12.)
                              + Vn * (massMatrix(Ux, Ux) + massMatrix(Uy, Uy)) * (c3 * 1 / 12.)
                              + An * (massMatrix(Ux, Ux) + massMatrix(Uy, Uy)) * (c4 * 1 / 12.)
                            ) * area;

    rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
    rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
    rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(2);
    rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(3);
    rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(4);
    rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(5);
  }

  //----------------------------------------------------------------------
  // traction term ∫∫ (𝐭.𝐯)  with 𝐭 = (𝑡𝑥, 𝑡𝑦, 𝑡𝑧) = (t[0], t[1], t[2])
  //----------------------------------------------------------------------
  // Index of the boundary condition. Needed to associate a CaseTable
  Int32 boundary_condition_index = 0;
  for (const auto& bs : options()->tractionBoundaryCondition()) {
    FaceGroup group = bs->surface();
    const CaseTableInfo& case_table_info = m_traction_case_table_list[boundary_condition_index];
    ++boundary_condition_index;

    Real3 trac; // traction in x, y and z

    if (bs->tractionInputFile.isPresent()) {

      String file_name = bs->tractionInputFile();
      info() << "Applying traction boundary conditions for surface " << group.name()
             << " via CaseTable" << file_name;
      CaseTable* inn = case_table_info.case_table;
      if (!inn)
        ARCANE_FATAL("CaseTable is null. Maybe there is a missing call to _readCaseTables()");
      if (file_name != case_table_info.file_name)
        ARCANE_FATAL("Incoherent CaseTable. The current CaseTable is associated to file '{0}'", case_table_info.file_name);
      inn->value(t, trac);

      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);
        for (Node node : iface->nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += trac.x * length / 2.;
            rhs_values[node_dof.dofId(node, 1)] += trac.y * length / 2.;
          }
        }
      }
      continue;
    }
    else {
      const UniqueArray<String> t_string = bs->t();
      Real3 trac;

      info() << "[ArcaneFem-Info] Applying Traction " << t_string;
      info() << "[ArcaneFem-Info] Traction surface '" << bs->surface().name() << "'";

      for (Int32 i = 0; i < t_string.size(); ++i) {
        trac[i] = 0.0;
        if (t_string[i] != "NULL") {
          trac[i] = std::stod(t_string[i].localstr());
        }
      }

      if (t_string[0] != "NULL" || t_string[1] != "NULL")
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);
          for (Node node : iface->nodes()) {
            if (node.isOwn()) {
              rhs_values[node_dof.dofId(node, 0)] += trac[0] * length / 2.;
              rhs_values[node_dof.dofId(node, 1)] += trac[1] * length / 2.;
            }
          }
        }
    }
  }

  //----------------------------------------------
  // Paraxial term assembly for LHS and RHS
  //----------------------------------------------
  for (const auto& bs : options()->paraxialBoundaryCondition()) {
    FaceGroup group = bs->surface();

    info() << "[ArcaneFem-Info] Applying constant Paraxial boundary conditions for surface " << group.name();

    ENUMERATE_ (Face, iface, group) {
      Face face = *iface;
      Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);

      RealVector<4> Uy = { 0., 1., 0., 1. };
      RealVector<4> Ux = { 1., 0., 1., 0. };

      RealVector<4> Un = { m_U[face.nodeId(0)].x, m_U[face.nodeId(0)].y, m_U[face.nodeId(1)].x, m_U[face.nodeId(1)].y };
      RealVector<4> Vn = { m_V[face.nodeId(0)].x, m_V[face.nodeId(0)].y, m_V[face.nodeId(1)].x, m_V[face.nodeId(1)].y };
      RealVector<4> An = { m_A[face.nodeId(0)].x, m_A[face.nodeId(0)].y, m_A[face.nodeId(1)].x, m_A[face.nodeId(1)].y };

      RealMatrix<4, 4>  ParaxialElementMatrixEdge2 = _computeParaxialElementMatrixEdge2(face);
      RealVector<4> rhs = length * (c7 * Un * ParaxialElementMatrixEdge2 - c8 * Vn * ParaxialElementMatrixEdge2  + c9 * An * ParaxialElementMatrixEdge2) ;

      rhs_values[node_dof.dofId(face.nodeId(0), 0)] += rhs(0);
      rhs_values[node_dof.dofId(face.nodeId(0), 1)] += rhs(1);
      rhs_values[node_dof.dofId(face.nodeId(1), 0)] += rhs(2);
      rhs_values[node_dof.dofId(face.nodeId(1), 1)] += rhs(3);

      if (t <= dt) {
        ParaxialElementMatrixEdge2 = c7 * length * ParaxialElementMatrixEdge2;
        Int32 n1_index = 0;
        for (Node node1 : face.nodes()) {
          Int32 n2_index = 0;
          for (Node node2 : face.nodes()) {
            Real v1 = ParaxialElementMatrixEdge2(2 * n1_index, 2 * n2_index);
            Real v2 = ParaxialElementMatrixEdge2(2 * n1_index, 2 * n2_index + 1);
            Real v3 = ParaxialElementMatrixEdge2(2 * n1_index + 1, 2 * n2_index);
            Real v4 = ParaxialElementMatrixEdge2(2 * n1_index + 1, 2 * n2_index + 1);
            if (node1.isOwn()) {
              DoFLocalId node1_dof1 = node_dof.dofId(node1, 0);
              DoFLocalId node1_dof2 = node_dof.dofId(node1, 1);
              DoFLocalId node2_dof1 = node_dof.dofId(node2, 0);
              DoFLocalId node2_dof2 = node_dof.dofId(node2, 1);

              if (bsr_matrix) {
                bsr_matrix->addValue(node1_dof1, node2_dof1, v1);
                bsr_matrix->addValue(node1_dof1, node2_dof2, v2);
                bsr_matrix->addValue(node1_dof2, node2_dof1, v3);
                bsr_matrix->addValue(node1_dof2, node2_dof2, v4);
              }
              else {
                m_linear_system.matrixAddValue(node1_dof1, node2_dof1, v1);
                m_linear_system.matrixAddValue(node1_dof1, node2_dof2, v2);
                m_linear_system.matrixAddValue(node1_dof2, node2_dof1, v3);
                m_linear_system.matrixAddValue(node1_dof2, node2_dof2, v4);
              }
            }
            ++n2_index;
          }
          ++n1_index;
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
                m_linear_system.matrixSetValue(dof_id, dof_id, Penalty);
                rhs_values[dof_id] = Penalty * u_dirichlet;
              }
            }
          }
        }
      }
    }
    else {
      ARCANE_FATAL("Dirichlet enforcing method not implemented use Penalty");
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
    else {
      ARCANE_FATAL("Dirichlet enforcing Method Not implemented use Penalty");
    }
  }

  //----------------------------------------------
  // Double-couple term assembly
  //----------------------------------------------

  // Index of the boundary condition. Needed to associate a CaseTable
  Int32 boundary_condition_index_dc = 0;

  for (const auto& bs : options()->doubleCouple()) {

    const CaseTableInfo& case_table_dc_info = m_double_couple_case_table_list[boundary_condition_index_dc];

    ++boundary_condition_index_dc;

    Real dc_force; // double-couple force

    String file_name = bs->doubleCoupleInputFile();
    info() << "Applying boundary conditions for surface via CaseTable" << file_name;

    CaseTable* dc_case_table_inn = case_table_dc_info.case_table;

    dc_case_table_inn->value(t, dc_force);

    NodeGroup north = bs->northNodeName();
    NodeGroup south = bs->southNodeName();
    NodeGroup east = bs->eastNodeName();
    NodeGroup west = bs->westNodeName();

    ENUMERATE_ (Node, inode, north) {
      Node node = *inode;
      DoFLocalId dof_id1 = node_dof.dofId(node, 0);
      rhs_values[dof_id1] = dc_force;
    }
    ENUMERATE_ (Node, inode, south) {
      Node node = *inode;
      DoFLocalId dof_id1 = node_dof.dofId(node, 0);
      rhs_values[dof_id1] = -dc_force;
    }
    ENUMERATE_ (Node, inode, east) {
      Node node = *inode;
      DoFLocalId dof_id2 = node_dof.dofId(node, 1);
      rhs_values[dof_id2] = -dc_force;
    }
    ENUMERATE_ (Node, inode, west) {
      Node node = *inode;
      DoFLocalId dof_id2 = node_dof.dofId(node, 1);
      rhs_values[dof_id2] = dc_force;
    }
  }
}

void FemModule::
_assembleLinearOperator3d(BSRMatrix* bsr_matrix)
{
  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);

    RealVector<12> Ux = { 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0. };
    RealVector<12> Uy = { 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0. };
    RealVector<12> Uz = { 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1. };

    RealVector<12> F = { f[0], f[1], f[2], f[0], f[1], f[2], f[0], f[1], f[2] };

    RealVector<12> Un = { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y, m_U[cell.nodeId(0)].z,
                              m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y, m_U[cell.nodeId(1)].z,
                              m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y, m_U[cell.nodeId(2)].z,
                              m_U[cell.nodeId(3)].x, m_U[cell.nodeId(3)].y, m_U[cell.nodeId(3)].z };

    RealVector<12> Vn = { m_V[cell.nodeId(0)].x, m_V[cell.nodeId(0)].y, m_V[cell.nodeId(0)].z,
                              m_V[cell.nodeId(1)].x, m_V[cell.nodeId(1)].y, m_V[cell.nodeId(1)].z,
                              m_V[cell.nodeId(2)].x, m_V[cell.nodeId(2)].y, m_V[cell.nodeId(2)].z,
                              m_V[cell.nodeId(3)].x, m_V[cell.nodeId(3)].y, m_V[cell.nodeId(3)].z };

    RealVector<12> An = { m_A[cell.nodeId(0)].x, m_A[cell.nodeId(0)].y, m_A[cell.nodeId(0)].z,
                              m_A[cell.nodeId(1)].x, m_A[cell.nodeId(1)].y, m_A[cell.nodeId(1)].z,
                              m_A[cell.nodeId(2)].x, m_A[cell.nodeId(2)].y, m_A[cell.nodeId(2)].z,
                              m_A[cell.nodeId(3)].x, m_A[cell.nodeId(3)].y, m_A[cell.nodeId(3)].z };

    //----------------------------------------------------------------------
    //  ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯)
    //----------------------------------------------------------------------
    RealVector<12> rhs = (  F * (1 / 4.)
                              + Un * (massMatrix(Ux, Ux) + massMatrix(Uy, Uy) + massMatrix(Uz, Uz)) * (c0 * 1 / 20.)
                              + Vn * (massMatrix(Ux, Ux) + massMatrix(Uy, Uy) + massMatrix(Uz, Uz)) * (c3 * 1 / 20.)
                              + An * (massMatrix(Ux, Ux) + massMatrix(Uy, Uy) + massMatrix(Uz, Uz)) * (c4 * 1 / 20.)
                             ) * volume;

    rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
    rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
    rhs_values[node_dof.dofId(cell.nodeId(0), 2)] += rhs(2);
    rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(3);
    rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(4);
    rhs_values[node_dof.dofId(cell.nodeId(1), 2)] += rhs(5);
    rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(6);
    rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(7);
    rhs_values[node_dof.dofId(cell.nodeId(2), 2)] += rhs(8);
    rhs_values[node_dof.dofId(cell.nodeId(3), 0)] += rhs(9);
    rhs_values[node_dof.dofId(cell.nodeId(3), 1)] += rhs(10);
    rhs_values[node_dof.dofId(cell.nodeId(3), 2)] += rhs(11);
  }

  //----------------------------------------------------------------------
  // traction term ∫∫ (𝐭.𝐯)  with 𝐭 = (𝑡𝑥, 𝑡𝑦, 𝑡𝑧) = (t[0], t[1], t[2])
  //----------------------------------------------------------------------
  // Index of the boundary condition. Needed to associate a CaseTable
  Int32 boundary_condition_index = 0;
  for (const auto& bs : options()->tractionBoundaryCondition()) {
    FaceGroup group = bs->surface();
    const CaseTableInfo& case_table_info = m_traction_case_table_list[boundary_condition_index];
    ++boundary_condition_index;

    Real3 trac; // traction in x, y and z

    if (bs->tractionInputFile.isPresent()) {

      String file_name = bs->tractionInputFile();
      info() << "Applying traction boundary conditions for surface " << group.name()
             << " via CaseTable" << file_name;
      CaseTable* inn = case_table_info.case_table;
      if (!inn)
        ARCANE_FATAL("CaseTable is null. Maybe there is a missing call to _readCaseTables()");
      if (file_name != case_table_info.file_name)
        ARCANE_FATAL("Incoherent CaseTable. The current CaseTable is associated to file '{0}'", case_table_info.file_name);
      inn->value(t, trac);

      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, m_node_coord);
        for (Node node : iface->nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += trac.x * area / 3.;
            rhs_values[node_dof.dofId(node, 1)] += trac.y * area / 3.;
            rhs_values[node_dof.dofId(node, 2)] += trac.z * area / 3.;
          }
        }
      }
      continue;
    }
    else {
      const UniqueArray<String> t_string = bs->t();
      Real3 trac;

      info() << "[ArcaneFem-Info] Applying Traction " << t_string;
      info() << "[ArcaneFem-Info] Traction surface '" << bs->surface().name() << "'";

      for (Int32 i = 0; i < t_string.size(); ++i) {
        trac[i] = 0.0;
        if (t_string[i] != "NULL") {
          trac[i] = std::stod(t_string[i].localstr());
        }
      }

      if (t_string[0] != "NULL" || t_string[1] != "NULL" || t_string[2] != "NULL")
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, m_node_coord);
          for (Node node : iface->nodes()) {
            if (node.isOwn()) {
              rhs_values[node_dof.dofId(node, 0)] += trac[0] * area / 3.;
              rhs_values[node_dof.dofId(node, 1)] += trac[1] * area / 3.;
              rhs_values[node_dof.dofId(node, 2)] += trac[2] * area / 3.;
            }
          }
        }
    }
  }

  //----------------------------------------------
  // Paraxial term assembly
  //----------------------------------------------
  for (const auto& bs : options()->paraxialBoundaryCondition()) {
    FaceGroup group = bs->surface();

    info() << "Applying constant Paraxial boundary conditions for surface " << group.name();

    ENUMERATE_ (Face, iface, group) {
      Face face = *iface;
      Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, m_node_coord);

      RealVector<9> Un = { m_U[face.nodeId(0)].x, m_U[face.nodeId(0)].y, m_U[face.nodeId(0)].z,
                               m_U[face.nodeId(1)].x, m_U[face.nodeId(1)].y, m_U[face.nodeId(1)].z,
                               m_U[face.nodeId(2)].x, m_U[face.nodeId(2)].y, m_U[face.nodeId(2)].z, };

      RealVector<9> Vn = { m_V[face.nodeId(0)].x, m_V[face.nodeId(0)].y, m_V[face.nodeId(0)].z,
                               m_V[face.nodeId(1)].x, m_V[face.nodeId(1)].y, m_V[face.nodeId(1)].z,
                               m_V[face.nodeId(2)].x, m_V[face.nodeId(2)].y, m_V[face.nodeId(2)].z, };

      RealVector<9> An = { m_A[face.nodeId(0)].x, m_A[face.nodeId(0)].y, m_A[face.nodeId(0)].z,
                               m_A[face.nodeId(1)].x, m_A[face.nodeId(1)].y, m_A[face.nodeId(1)].z,
                               m_A[face.nodeId(2)].x, m_A[face.nodeId(2)].y, m_A[face.nodeId(2)].z, };

      RealMatrix<9, 9>  ParaxialElementMatrixTria3 = _computeParaxialElementMatrixTria3(face);
      RealVector<9> rhs = area * (c7 * Un * ParaxialElementMatrixTria3 - c8 * Vn * ParaxialElementMatrixTria3  + c9 * An * ParaxialElementMatrixTria3) ;

      rhs_values[node_dof.dofId(face.nodeId(0), 0)] += rhs(0);
      rhs_values[node_dof.dofId(face.nodeId(0), 1)] += rhs(1);
      rhs_values[node_dof.dofId(face.nodeId(0), 2)] += rhs(2);
      rhs_values[node_dof.dofId(face.nodeId(1), 0)] += rhs(3);
      rhs_values[node_dof.dofId(face.nodeId(1), 1)] += rhs(4);
      rhs_values[node_dof.dofId(face.nodeId(1), 2)] += rhs(5);
      rhs_values[node_dof.dofId(face.nodeId(2), 0)] += rhs(6);
      rhs_values[node_dof.dofId(face.nodeId(2), 1)] += rhs(7);
      rhs_values[node_dof.dofId(face.nodeId(2), 2)] += rhs(8);

      if (t <= dt) {
        ParaxialElementMatrixTria3 = area * c7 * ParaxialElementMatrixTria3;
        Int32 n1_index = 0;
        for (Node node1 : face.nodes()) {
          Int32 n2_index = 0;
          for (Node node2 : face.nodes()) {
            Real v1 = ParaxialElementMatrixTria3(3 * n1_index, 3 * n2_index);
            Real v2 = ParaxialElementMatrixTria3(3 * n1_index, 3 * n2_index + 1);
            Real v3 = ParaxialElementMatrixTria3(3 * n1_index, 3 * n2_index + 2);

            Real v4 = ParaxialElementMatrixTria3(3 * n1_index + 1, 3 * n2_index);
            Real v5 = ParaxialElementMatrixTria3(3 * n1_index + 1, 3 * n2_index + 1);
            Real v6 = ParaxialElementMatrixTria3(3 * n1_index + 1, 3 * n2_index + 2);

            Real v7 = ParaxialElementMatrixTria3(3 * n1_index + 2, 3 * n2_index);
            Real v8 = ParaxialElementMatrixTria3(3 * n1_index + 2, 3 * n2_index + 1);
            Real v9 = ParaxialElementMatrixTria3(3 * n1_index + 2, 3 * n2_index + 2);
            if (node1.isOwn()) {
              DoFLocalId node1_dof1 = node_dof.dofId(node1, 0);
              DoFLocalId node1_dof2 = node_dof.dofId(node1, 1);
              DoFLocalId node1_dof3 = node_dof.dofId(node1, 2);

              DoFLocalId node2_dof1 = node_dof.dofId(node2, 0);
              DoFLocalId node2_dof2 = node_dof.dofId(node2, 1);
              DoFLocalId node2_dof3 = node_dof.dofId(node2, 2);

              if (bsr_matrix) {
                bsr_matrix->addValue(node1_dof1, node2_dof1, v1);
                bsr_matrix->addValue(node1_dof1, node2_dof2, v2);
                bsr_matrix->addValue(node1_dof1, node2_dof3, v3);

                bsr_matrix->addValue(node1_dof2, node2_dof1, v4);
                bsr_matrix->addValue(node1_dof2, node2_dof2, v5);
                bsr_matrix->addValue(node1_dof2, node2_dof3, v6);

                bsr_matrix->addValue(node1_dof3, node2_dof1, v7);
                bsr_matrix->addValue(node1_dof3, node2_dof2, v8);
                bsr_matrix->addValue(node1_dof3, node2_dof3, v9);
              }
              else {
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
            }
            ++n2_index;
          }
          ++n1_index;
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
                m_linear_system.matrixSetValue(dof_id, dof_id, Penalty);
                rhs_values[dof_id] = Penalty * u_dirichlet;
              }
            }
          }
        }
      }
    }
    else {
      ARCANE_FATAL("Dirichlet enforcing Method Not implemented use Penalty");
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
    else {
      ARCANE_FATAL("Dirichlet enforcing Method Not implemented use Penalty");
    }
  }

  //----------------------------------------------
  // Double-couple term assembly
  //----------------------------------------------

  // Index of the boundary condition. Needed to associate a CaseTable
  Int32 boundary_condition_index_dc = 0;

  for (const auto& bs : options()->doubleCouple()) {

    const CaseTableInfo& case_table_dc_info = m_double_couple_case_table_list[boundary_condition_index_dc];

    ++boundary_condition_index_dc;

    Real dc_force; // double-couple force

    String file_name = bs->doubleCoupleInputFile();
    info() << "Applying boundary conditions for surface via CaseTable" << file_name;

    CaseTable* dc_case_table_inn = case_table_dc_info.case_table;

    dc_case_table_inn->value(t, dc_force);

    NodeGroup north = bs->northNodeName();
    NodeGroup south = bs->southNodeName();
    NodeGroup east = bs->eastNodeName();
    NodeGroup west = bs->westNodeName();

    ENUMERATE_ (Node, inode, north) {
      Node node = *inode;
      DoFLocalId dof_id1 = node_dof.dofId(node, 0);
      rhs_values[dof_id1] = dc_force;
    }
    ENUMERATE_ (Node, inode, south) {
      Node node = *inode;
      DoFLocalId dof_id1 = node_dof.dofId(node, 0);
      rhs_values[dof_id1] = -dc_force;
    }
    ENUMERATE_ (Node, inode, east) {
      Node node = *inode;
      DoFLocalId dof_id2 = node_dof.dofId(node, 2);
      rhs_values[dof_id2] = -dc_force;
    }
    ENUMERATE_ (Node, inode, west) {
      Node node = *inode;
      DoFLocalId dof_id2 = node_dof.dofId(node, 2);
      rhs_values[dof_id2] = dc_force;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (t <= dt) {
    if (mesh()->dimension() == 2)
      if(m_matrix_format == "DOK")
        _assembleBilinearOperatorTria3();
      else if(m_matrix_format == "BSR")
        _assembleBilinearOperatorTria3Gpu();
    if (mesh()->dimension() == 3)
      _assembleBilinearOperatorTetra4();
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "lhs-matrix-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTria3Gpu()
{
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();
  auto fn_cv = m_connectivity_view.faceNode();
  auto command = makeCommand(acceleratorMng()->defaultQueue());
  FaceInfoListView faces_infos(mesh()->faceFamily());
  auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
  auto c0_copy = c0;
  auto c1_copy = c1;
  auto c2_copy = c2;

  m_bsr_format.computeSparsity();
  m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord, c0_copy, c1_copy, c2_copy); });
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTria3()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    auto K_e = _computeElementMatrixTria3(cell); // element stiffness matrix
    // assemble elementary matrix into  the global one elementary terms are
    // positioned into  K according  to the rank of associated  node in the
    // mesh.nodes list  and according the dof number. Here  for  each  node
    // two dofs  exists [u1,u2]. For each TRIA3 there are 3 nodes hence the
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
_solve()
{
  info() << "[ArcaneFem-Info] Started module  _solve()";
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
  info() << "[ArcaneFem-Info] Started module  _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200) {
    int p = std::cout.precision();
    std::cout.precision(17);
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      std::cout << "( N_id, u1, u2, u3 ) = ( "
                << node.uniqueId() << ", " << m_dU[node].x << ", " << m_dU[node].y << ", " << m_dU[node].z
                << ")\n";
    }
    std::cout.precision(p);
  }

  String filename = options()->resultFile();
  const double epsilon = 1.0e-4;
  const double min_value_to_test = 1.0e-16;

  info() << "[ArcaneFem-Info] Validating results filename=" << filename << " epsilon =" << epsilon;

  if (!filename.empty())
    Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_dU, epsilon, min_value_to_test);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "cross-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
