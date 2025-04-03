// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* FEM code to test vectorial FE for Elastodynamics problem.                 */
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

  m_dofs_on_nodes.initialize(mesh(), mesh()->dimension());

  _getParameters();

  t = dt;
  tmax = tmax - dt;
  m_global_deltat.assign(dt);

  _readCaseTables();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "initialize", elapsedTime);
}

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
    m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");
  }

  if (m_petsc_flags != NULL){
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  _doStationarySolve();
  _updateTime();

  if ((t > tmax + dt - 1e-8) && m_cross_validation)
    _validateResults();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "compute", elapsedTime);
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
/*
 * @brief Reads the parameters from the options.
 *
 * This method retrieves various parameters such as time step, material properties,
 * and time discretization methods from the options provided in the simulation.
 * It also initializes the material properties based on the retrieved parameters.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_getParameters()
{
  info() << "[ArcaneFem-Info] Started module  _getParameters()";
  Real elapsedTime = platform::getRealTime();

  //--------- time parameters -----------//
  tmax = options()->tmax(); // max time
  dt = options()->dt(); // time step

  //--- damping term parameter ---//
  etam = options()->etam(); // damping param etam
  etak = options()->etak(); // damping param etak

  //--- time discretization parameter ---//
  alpm = options()->alpm(); // time discretization param alpm
  alpf = options()->alpf(); // time discretization param alpf

  //--------- material parameter ---------//
  E = options()->E(); // Youngs modulus
  nu = options()->nu(); // Poission ratio
  rho = options()->rho(); // Density

  mu = E / (2 * (1 + nu)); // lame parameter mu
  lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter lambda

  if( options()->mu.isPresent())
    mu = options()->mu;

  if( options()->lambda.isPresent())
    lambda = options()->lambda;

  //--------- body force ---------//
  if (options()->f.isPresent()) {
    const UniqueArray<String> f_string = options()->f();
    info() << "[ArcaneFem-Info] Applying Bodyforce " << f_string;
    for (Int32 i = 0; i < f_string.size(); ++i) {
      if (f_string[i] != "NULL") {
        f[i] = std::stod(f_string[i].localstr());
      }
    }
  }

  //----- time discretization Newmark-Beta or Generalized-alpha  -----//
  if (options()->timeDiscretization == "Newmark-beta") {

    info() << "[ArcaneFem-Info] Apply time discretization via Newmark-beta ";

    gamma = 0.5;
    beta  = (1./4.)*(gamma+0.5)*(gamma+0.5)  ;

    c0 = rho / (beta * dt * dt) + etam * rho * gamma / beta / dt;
    c1 = lambda + lambda * etak * gamma / beta / dt;
    c2 = mu + mu * etak * gamma / beta / dt;
    c3 = rho / beta / dt - etam * rho * (1 - gamma / beta);
    c4 = rho * ((1. - 2. * beta) / 2. / beta - etam * dt * (1. - gamma / 2 / beta));
    c5 = -lambda * etak * gamma / beta / dt;
    c6 = -mu * etak * gamma / beta / dt;
    c7 = etak * lambda * (gamma / beta - 1);
    c8 = etak * lambda * dt * ((1. - 2 * beta) / 2. / beta - (1. - gamma));
    c9 = etak * mu * (gamma / beta - 1);
    c10 = etak * mu * dt * ((1. - 2 * beta) / 2. / beta - (1. - gamma));
  }
  else if (options()->timeDiscretization == "Generalized-alpha") {

    info() << "[ArcaneFem-Info] Apply time discretization via Generalized-alpha ";

    gamma = 0.5 + alpf - alpm                ;
    beta  = (1./4.)*(gamma+0.5)*(gamma+0.5)  ;

    c0 = rho * (1. - alpm) / (beta * dt * dt) + etam * rho * gamma * (1 - alpf) / beta / dt;
    c1 = lambda * (1. - alpf) + lambda * etak * gamma * (1. - alpf) / beta / dt;
    c2 = mu * (1. - alpf) + mu * etak * gamma * (1. - alpf) / beta / dt;
    c3 = rho * (1. - alpm) / beta / dt - etam * rho * (1 - gamma * (1 - alpf) / beta);
    c4 = rho * ((1. - alpm) * (1. - 2. * beta) / 2. / beta - alpm - etam * dt * (1. - alpf) * (1. - gamma / 2 / beta));
    c5 = lambda * alpf - lambda * etak * gamma * (1. - alpf) / beta / dt;
    c6 = mu * alpf - mu * etak * gamma * (1. - alpf) / beta / dt;
    c7 = etak * lambda * (gamma * (1. - alpf) / beta - 1);
    c8 = etak * lambda * dt * (1. - alpf) * ((1. - 2 * beta) / 2. / beta - (1. - gamma));
    c9 = etak * mu * (gamma * (1. - alpf) / beta - 1);
    c10 = etak * mu * dt * (1. - alpf) * ((1. - 2 * beta) / 2. / beta - (1. - gamma));
  }
  else {
    ARCANE_FATAL("Only Newmark-beta | Generalized-alpha are supported for time-discretization ");
  }

  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->crossValidation();
  m_petsc_flags = options()->petscFlags();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*
  * @brief Reads case tables for traction boundary conditions.
  *
  * This method reads the case tables specified in the options and stores
  * them in a list for later use.
  */
/*---------------------------------------------------------------------------*/

void FemModule::
_readCaseTables()
{
  IParallelMng* pm = subDomain()->parallelMng();
  for (const auto& bs : options()->tractionBoundaryCondition()) {
    CaseTable* case_table = nullptr;
    String file_name;
    if (bs->tractionInputFile.isPresent()) {
      file_name = bs->tractionInputFile();
      case_table = readFileAsCaseTable(pm, file_name, 3);
    }
    m_traction_case_table_list.add(CaseTableInfo{ file_name, case_table });
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateVariables()
{
  info() << "[ArcaneFem-Info] Started module  _updateVariables()";
  Real elapsedTime = platform::getRealTime();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    if (mesh()->dimension() == 2)
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real u1_val = dof_u[node_dof.dofId(node, 0)];
      Real u2_val = dof_u[node_dof.dofId(node, 1)];
      Real3 u_disp;
      u_disp.x = u1_val;
      u_disp.y = u2_val;
      u_disp.z = 0.0;
      m_dU[node] = u_disp;
    }
    if (mesh()->dimension() == 3)
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real u1_val = dof_u[node_dof.dofId(node, 0)];
      Real u2_val = dof_u[node_dof.dofId(node, 1)];
      Real u3_val = dof_u[node_dof.dofId(node, 2)];
      Real3 u_disp;
      u_disp.x = u1_val;
      u_disp.y = u2_val;
      u_disp.z = u3_val;
      m_dU[node] = u_disp;
    }
  }

  m_dU.synchronize();
  m_U.synchronize();
  m_V.synchronize();
  m_A.synchronize();

  // Note at this stage we already have calculated dU
  Real alocX;
  Real alocY;
  Real alocZ;

  VariableDoFReal& dof_u(m_linear_system.solutionVariable());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

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

    elapsedTime = platform::getRealTime() - elapsedTime;
    ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"update-variables", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*
 * @brief Assembles the linear operator for the FEM simulation.
 *
 * This method computes the linear operator for the FEM simulation
 * by assembling the linear parts of the system. It uses the
 * material properties and the mesh to compute the contributions
 * from each element. The results are stored in the right-hand side
 * of the linear system.
 * 
 * @note This method assumes that the mesh and material properties
 *  are already set up.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (mesh()->dimension() == 2)
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
      Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
      Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

      RealVector<6> Uy = { 0., 1., 0., 1., 0., 1. };
      RealVector<6> Ux = { 1., 0., 1., 0., 1., 0. };
      RealVector<6> F = { f[0], f[1], f[0], f[1], f[0], f[1] };
      RealVector<6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
      RealVector<6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
      RealVector<6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
      RealVector<6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

      RealVector<6> Un = { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y,
                               m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y,
                               m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y };

                               RealVector<6>Vn = { m_V[cell.nodeId(0)].x, m_V[cell.nodeId(0)].y,
                               m_V[cell.nodeId(1)].x, m_V[cell.nodeId(1)].y,
                               m_V[cell.nodeId(2)].x, m_V[cell.nodeId(2)].y };

                               RealVector<6> An = { m_A[cell.nodeId(0)].x, m_A[cell.nodeId(0)].y,
                               m_A[cell.nodeId(1)].x, m_A[cell.nodeId(1)].y,
                               m_A[cell.nodeId(2)].x, m_A[cell.nodeId(2)].y };

      //----------------------------------------------------------------------
      //  ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯) +
      //  ∫∫∫ (c₅)(∇𝐮ₙ.∇𝐯) + ∫∫∫ (c₆)(ε(𝐮ₙ):ε(𝐯)) +
      //  ∫∫∫ (c₇)(∇𝐮ᵗₙ.∇𝐯) + ∫∫∫ (c₉)(ε(𝐮ᵗₙ):ε(𝐯)) +
      //  ∫∫∫ (c₈)(∇𝐮ᵗᵗₙ.∇𝐯) + ∫∫∫ (c₁₀)(ε(𝐮ᵗᵗₙ):ε(𝐯))
      //----------------------------------------------------------------------
      RealVector<6> rhs = ( F * (1/3.)
                              + Un * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy))*(c0*1/12.)
                              + Vn * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy))*(c3*1/12.)
                              + An * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy))*(c4*1/12.)
                              - Un * ((dyUy ^ dxUx) + (dxUx ^ dyUy)  +  (dxUx ^ dxUx) + (dyUy ^ dyUy)) * c5
                              - Un * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy)) +  ((dxUy + dyUx) ^ (dyUx + dxUy)))*c6
                              + Vn * ((dyUy ^ dxUx) + (dxUx ^ dyUy)  +  (dxUx ^ dxUx) + (dyUy ^ dyUy)) * c7
                              + Vn * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy)) +  ((dxUy + dyUx) ^ (dyUx + dxUy)))*c9
                              + An * ((dyUy ^ dxUx) + (dxUx ^ dyUy)  +  (dxUx ^ dxUx) + (dyUy ^ dyUy)) * c8
                              + An * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy)) +  ((dxUy + dyUx) ^ (dyUx + dxUy)))*c10
                              ) * area;

      rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
      rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
      rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(2);
      rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(3);
      rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(4);
      rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(5);
    }

  if (mesh()->dimension() == 3)
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);
      Real4 dxu = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
      Real4 dyu = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
      Real4 dzu = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);

      RealVector<12> Uy = { 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0. };
      RealVector<12> Ux = { 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0. };
      RealVector<12> Uz = { 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1. };

      RealVector<12> F = { f[0], f[1], f[2], f[0], f[1], f[2], f[0], f[1], f[2] };
      RealVector<12> dxUx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
      RealVector<12> dyUx = { dyu[0], 0., 0.,    dyu[1], 0., 0.,    dyu[2], 0., 0.,    dyu[3], 0., 0. };
      RealVector<12> dzUx = { dzu[0], 0., 0.,    dzu[1], 0., 0.,    dzu[2], 0., 0.,    dzu[3], 0., 0. };

      RealVector<12> dxUy = { 0., dxu[0], 0.,    0., dxu[1], 0.,    0., dxu[2], 0.,    0., dxu[3], 0. };
      RealVector<12> dyUy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
      RealVector<12> dzUy = { 0., dzu[0], 0.,    0., dzu[1], 0.,    0., dzu[2], 0.,    0., dzu[3], 0. };

      RealVector<12> dxUz = { 0., 0., dxu[0],    0., 0., dxu[1],    0., 0., dxu[2],    0., 0., dxu[3] };
      RealVector<12> dyUz = { 0., 0., dyu[0],    0., 0., dyu[1],    0., 0., dyu[2],    0., 0., dyu[3] };
      RealVector<12> dzUz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

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
      //  ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯) +
      //  ∫∫∫ (c₅)(∇𝐮ₙ.∇𝐯) + ∫∫∫ (c₆)(ε(𝐮ₙ):ε(𝐯)) +
      //  ∫∫∫ (c₇)(∇𝐮ᵗₙ.∇𝐯) + ∫∫∫ (c₉)(ε(𝐮ᵗₙ):ε(𝐯)) +
      //  ∫∫∫ (c₈)(∇𝐮ᵗᵗₙ.∇𝐯) + ∫∫∫ (c₁₀)(ε(𝐮ᵗᵗₙ):ε(𝐯))
      //----------------------------------------------------------------------
      RealVector<12> rhs = ( F * (1/4.)
                              + Un * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy) + massMatrix(Uz,Uz))*(c0*1/20.)
                              + Vn * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy) + massMatrix(Uz,Uz))*(c3*1/20.)
                              + An * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy) + massMatrix(Uz,Uz))*(c4*1/20.)
                              - Un * ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) +
                                      (dyUy ^ dxUx) + (dxUx ^ dyUy) +
                                      (dzUz ^ dxUx) + (dxUx ^ dzUz) +
                                      (dyUy ^ dzUz) + (dzUz ^ dyUy)) * c5
                              - Un * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) ) +
                                         ( ((dxUy + dyUx) ^ (dyUx + dxUy)) +
                                           ((dzUy + dyUz) ^ (dyUz + dzUy)) +
                                           ((dxUz + dzUx) ^ (dzUx + dxUz)) ))*c6
                              + Vn * ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) +
                                      (dyUy ^ dxUx) + (dxUx ^ dyUy) +
                                      (dzUz ^ dxUx) + (dxUx ^ dzUz) +
                                      (dyUy ^ dzUz) + (dzUz ^ dyUy)) * c7
                              + Vn * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) ) +
                                      ( ((dxUy + dyUx) ^ (dyUx + dxUy)) +
                                        ((dzUy + dyUz) ^ (dyUz + dzUy)) +
                                        ((dxUz + dzUx) ^ (dzUx + dxUz)) ))*c9
                                        + An * ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) +
                                        (dyUy ^ dxUx) + (dxUx ^ dyUy) +
                                        (dzUz ^ dxUx) + (dxUx ^ dzUz) +
                                        (dyUy ^ dzUz) + (dzUz ^ dyUy)) * c8
                                + An * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) ) +
                                        ( ((dxUy + dyUx) ^ (dyUx + dxUy)) +
                                          ((dzUy + dyUz) ^ (dyUz + dzUy)) +
                                          ((dxUz + dzUx) ^ (dzUx + dxUz)) ))*c10
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

      if (mesh()->dimension() == 2)
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
      if (mesh()->dimension() == 3)
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

      if (mesh()->dimension() == 2)
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

      if (mesh()->dimension() == 3)
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
                if (t <= dt)
                  m_linear_system.eliminateRow(dof_id, u_dirichlet);
                else
                  rhs_values[dof_id] =  u_dirichlet;
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
                if (t <= dt)
                  m_linear_system.eliminateRowColumn(dof_id, u_dirichlet);
                else
                  rhs_values[dof_id] =  u_dirichlet;
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
              if (t <= dt)
                m_linear_system.eliminateRow(dof_id, u_dirichlet);
              else
                rhs_values[dof_id] =  u_dirichlet;
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
              if (t <= dt)
                m_linear_system.eliminateRowColumn(dof_id, u_dirichlet);
              else
                rhs_values[dof_id] =  u_dirichlet;
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
_assembleBilinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (t <= dt) {
    if (mesh()->dimension() == 2)
      _assembleBilinearOperatorTria3();
    else
      _assembleBilinearOperatorTetra4();
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "lhs-matrix-assembly", elapsedTime);
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
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"cross-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
