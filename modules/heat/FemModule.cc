﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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

  //_initBoundaryconditions();    // initialize boundary conditions
  _initTime();                  // initialize time
  _getParameters();             // get material parameters
  _initTemperature();           // initialize temperature
  m_global_deltat.assign(dt);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"[ArcaneFem-Timer] initialize", elapsedTime);
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

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  _doStationarySolve();
  _updateVariables();
  _updateTime();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"[ArcaneFem-Timer] compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_initTime()
{
  info() << "[ArcaneFem-Info] Started module _initTime()";

  tmax   = options()->tmax();
  dt     = options()->dt();

  tmax = tmax ;
  t    = 0.0;
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

  {
    // Copy Node temperature to Node temperature old
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      m_node_temperature_old[node] = m_node_temperature[node];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_initTemperature()
{
  info() << "[ArcaneFem-Info] Started module _initTemperature()";

  Tinit   = options()->Tinit();

  {
    // Copy Node temperature to Node temperature old
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      m_node_temperature_old[node] = Tinit;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{

  // Assemble the FEM bilinear operator (LHS - matrix A)
  _assembleBilinearOperator();

  // Assemble the FEM linear operator (RHS - vector b)
  _assembleLinearOperator();

  // # T=linalg.solve(K,RHS)
  _solve();

  // # update time
  if (t >= tmax)
    _checkResultFile();
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
  ElementNodes = 3.;

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    m_cell_lambda[cell] = lambda;
  }

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
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"[ArcaneFem-Timer] get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
// Assemble the FEM linear operator
//  - This function enforces a Dirichlet boundary condition in a weak sense
//    via the penalty method
//  - The method also adds source term
//  - Also adds external fluxes
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  //----------------------------------------------
  // Convection term assembly
  //----------------------------------------------
  //
  //  $int_{dOmega_C}( h*Text*v^h)$
  //  only for nodes that are non-Dirichlet
  //----------------------------------------------
  for (const auto& bs : options()->convectionBoundaryCondition()) {
    FaceGroup group = bs->surface();
    h = bs->h();
    Text = bs->Text();
    ENUMERATE_ (Face, iface, group) {
      Face face = *iface;
      Real length =  ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);
      for (Node node : iface->nodes()) {
        if (node.isOwn())
          rhs_values[node_dof.dofId(node, 0)] += h * Text * length / 2.;
      }
    }
  }

  // Helper lambda to apply boundary conditions
  auto applyBoundaryConditions = [&](auto BCFunctions) {
    m_node_temperature_old.mult(1.0 / dt);
    BCFunctions.applyVariableSourceToRhs(m_node_temperature_old, mesh(), node_dof, m_node_coord, rhs_values);

    BC::IArcaneFemBC* bc = options()->boundaryConditions();
    if (bc) {
      for (BC::INeumannBoundaryCondition* bs : bc->neumannBoundaryConditions())
        BCFunctions.applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);

      for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
        BCFunctions.applyDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);

      for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions())
        BCFunctions.applyPointDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);
    }
  };

  using BCFunctions = ArcaneFemFunctions::BoundaryConditions2D;
  applyBoundaryConditions(BCFunctions());

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"[ArcaneFem-Timer] rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeDyOfRealTRIA3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real f0 = m_node_temperature[cell.nodeId(0)];
  Real f1 = m_node_temperature[cell.nodeId(1)];
  Real f2 = m_node_temperature[cell.nodeId(2)];

  // Using Cramer's rule  det (adj (A)) / det (A)
  return ( m0.x*(f1 - f2) - f0*(m1.x - m2.x) + (f2*m1.x - f1*m2.x) ) /
         ( m0.x*(m1.y - m2.y) - m0.y*(m1.x - m2.x) + (m1.x*m2.y - m2.x*m1.y) );
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeDxOfRealTRIA3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real f0 = m_node_temperature[cell.nodeId(0)];
  Real f1 = m_node_temperature[cell.nodeId(1)];
  Real f2 = m_node_temperature[cell.nodeId(2)];

  // Using Cramer's rule  det (adj (A)) / det (A)
  return ( f0*(m1.y - m2.y) - m0.y*(f1 - f2) + (f1*m2.y - f2*m1.y) ) /
         ( m0.x*(m1.y - m2.y) - m0.y*(m1.x - m2.x) + (m1.x*m2.y - m2.x*m1.y) );
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real2 FemModule::
_computeDxDyOfRealTRIA3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real f0 = m_node_temperature[cell.nodeId(0)];
  Real f1 = m_node_temperature[cell.nodeId(1)];
  Real f2 = m_node_temperature[cell.nodeId(2)];

  Real detA = ( m0.x*(m1.y - m2.y) - m0.y*(m1.x - m2.x) + (m1.x*m2.y - m2.x*m1.y) );

  Real2 DX;
        DX.x = ( m0.x*(f1 - f2) - f0*(m1.x - m2.x) + (f2*m1.x - f1*m2.x) ) / detA;
        DX.y = ( f0*(m1.y - m2.y) - m0.y*(f1 - f2) + (f1*m2.y - f2*m1.y) ) / detA;

  return DX ;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

RealMatrix<2, 2> FemModule::
_computeElementMatrixEDGE2(Face face)
{
  // Get coordinates of the triangle element  EDGE2
  //------------------------------------------------
  //                   1         0
  //                   o . . . . o
  //
  //------------------------------------------------
  Real3 m0 = m_node_coord[face.nodeId(0)];
  Real3 m1 = m_node_coord[face.nodeId(1)];

  Real area =  ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);

  Real2 dPhi0(m0.y - m1.y, m1.x - m0.x);
  Real2 dPhi1(m1.y - m0.y, m0.x - m1.x);

  RealMatrix<1, 2> b_matrix;
  RealMatrix<2, 1> bT_matrix;
  RealMatrix<2, 2> int_DOmega_i;

  for (Int32 i = 0; i<2; i++)
    for (Int32 j = 0; j<2; j++)
      int_DOmega_i(i,j) = 0.;

  // uv //
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 1.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = area/3;
  bT_matrix(1, 0) = area/3;

  bT_matrix.multInPlace(1.);

  RealMatrix<2, 2> int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<2; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(h);
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  return int_DOmega_i;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

RealMatrix<3, 3> FemModule::
_computeElementMatrixTRIA3(Cell cell)
{
  // Get coordinates of the triangle element  TRI3
  //------------------------------------------------
  //                  0 o
  //                   . .
  //                  .   .
  //                 .     .
  //              1 o . . . o 2
  //------------------------------------------------
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real area =  ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  RealMatrix<1, 3> b_matrix;
  RealMatrix<3, 1> bT_matrix;
  RealMatrix<3, 3> int_Omega_i;

  for (Int32 i = 0; i<3; i++)
    for (Int32 j = 0; j<3; j++)
      int_Omega_i(i,j) = 0.;

// -----------------------------------------------------------------------------
//  lambda*(dx(u)dx(v) + dy(u)dy(v)) + uv/dt
//------------------------------------------------------------------------------


  // dx(u)dx(v) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = dPhi1.x/area;
  b_matrix(0, 2) = dPhi2.x/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = dPhi1.x;
  bT_matrix(2, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.5f);

  RealMatrix<3, 3> int_dxUdxV = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxUdxV);


  // dy(u)dy(v) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = dPhi1.y/area;
  b_matrix(0, 2) = dPhi2.y/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = dPhi1.y;
  bT_matrix(2, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5f);

  RealMatrix<3, 3> int_dyUdyV = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyUdyV);

  int_Omega_i.multInPlace(lambda);

  // uv //
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 1.;
  b_matrix(0, 2) = 1.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = area/3.;
  bT_matrix(1, 0) = area/3.;
  bT_matrix(2, 0) = area/3.;

  bT_matrix.multInPlace(0.5f);

  RealMatrix<3, 3> int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<3; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(1./dt);

  int_Omega_i = matrixAddition( int_Omega_i, int_UV);

  return int_Omega_i;
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

  _assembleBilinearOperatorTRIA3();
  _assembleBilinearOperatorEDGE2();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"[ArcaneFem-Timer] lhs-matrix-assembly", elapsedTime);
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

    lambda = m_cell_lambda[cell];                 // lambda is always considered cell constant
    auto K_e = _computeElementMatrixTRIA3(cell);  // element stiffness matrix
    // assemble elementary matrix into the global one elementary terms are
    // positioned into K according to  the rank of associated  node in the
    // mesh.nodes list and according the dof number. For each TRIA3  there
    // are 3 nodes hence the elementary stiffness matrix  size  is (3x3)=6
    // will be  filled.
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e(n1_index, n2_index);
        // m_k_matrix(node1.localId(), node2.localId()) += v;
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
_assembleBilinearOperatorEDGE2()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  for (const auto& bs : options()->convectionBoundaryCondition()) {
    FaceGroup group = bs->surface();
    h = bs->h();
    ENUMERATE_ (Face, iface, group) {
      Face face = *iface;

      auto K_e = _computeElementMatrixEDGE2(face);  // element stiffness matrix

      Int32 n1_index = 0;
      for (Node node1 : face.nodes() ) {
        Int32 n2_index = 0;
        for (Node node2 : face.nodes()) {
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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  info() << "[ArcaneFem-Info] Started module _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.solve();

  {
    VariableDoFReal& dof_temperature(m_linear_system.solutionVariable());
    // Copy RHS DoF to Node temperature
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_temperature[node_dof.dofId(node, 0)];
      m_node_temperature[node] = v;
    }

    m_node_temperature.synchronize();
    m_node_temperature_old.synchronize();

    if(m_flux.tagValue("PostProcessing")=="1") {
      ENUMERATE_ (Cell, icell, allCells()) {
        Cell cell = *icell;

        Real2 DX = _computeDxDyOfRealTRIA3(cell);
        m_flux[cell].x = -m_cell_lambda[cell] * DX.x;
        m_flux[cell].y = -m_cell_lambda[cell] * DX.y;
        m_flux[cell].z = 0.;
      }

      m_flux.synchronize();
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"[ArcaneFem-Timer] solve-linear-system", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_checkResultFile()
{
  info() << "[ArcaneFem-Info] Started module _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200) {
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "T[" << node.uniqueId() << "] = "
             << m_node_temperature[node];
    }
  }

  String filename = options()->resultFile();
  info() << "ValidateResultFile filename=" << filename;

  if (!filename.empty())
    checkNodeResultFile(traceMng(), filename, m_node_temperature, 1.0e-4);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(),"[ArcaneFem-Timer] result-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
