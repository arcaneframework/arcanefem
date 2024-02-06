﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2024 */
/*                                                                           */
/* FEM code to test vectorial FE for Soildynamics problem.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/CaseTable.h>

#include "IDoFLinearSystemFactory.h"
#include "Fem_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Module Fem.
 */
class FemModule
: public ArcaneFemObject
{
 public:

  explicit FemModule(const ModuleBuildInfo& mbi)
  : ArcaneFemObject(mbi)
  , m_dofs_on_nodes(mbi.subDomain()->traceMng())
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
  }
  ~FemModule()
  {
    for( const CaseTableInfo&  t : m_traction_case_table_list )
      delete t.case_table;
    for( const CaseTableInfo&  t : m_double_couple_case_table_list_north )
      delete t.case_table;
    for( const CaseTableInfo&  t : m_double_couple_case_table_list_south )
      delete t.case_table;
    for( const CaseTableInfo&  t : m_double_couple_case_table_list_east )
      delete t.case_table;
    for( const CaseTableInfo&  t : m_double_couple_case_table_list_west )
      delete t.case_table;
  }

 public:

  //! Method called at each iteration
  void compute() override;

  //! Method called at the beginning of the simulation
  void startInit() override;

  VersionInfo versionInfo() const override
  {
    return VersionInfo(1, 0, 0);
  }

 private:

  Real t;                     // time variable
  Real dt;                    // time step
  Real tmax;                  // max time

  Real etam;                  // time discretization param etam
  Real etak;                  // time discretization param etak
  Real alpm;                  // time discretization param alpm
  Real alpf;                  // time discretization param alpf
  Real beta;                  // time discretization param beta
  Real gamma;                 // time discretization param gamma

  Real E;                     // Youngs modulus
  Real nu;                    // Poissons ratio
  Real rho;                   // Density
  Real cp;                    // Primary wave velocity of soil
  Real cs;                    // Secondary wave velocity of soil
  Real f1;                    // Body force in x
  Real f2;                    // Body force in y
  Real mu;                    // Lame parameter mu
  Real mu2;                   // Lame parameter mu * 2
  Real lambda;                // Lame parameter lambda

  Real c0;                    // constant
  Real c1;                    // constant
  Real c2;                    // constant
  Real c3;                    // constant
  Real c4;                    // constant
  Real c5;                    // constant
  Real c6;                    // constant
  Real c7;                    // constant
  Real c8;                    // constant
  Real c9;                    // constant

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;

  // Struct to make sure we are using a CaseTable associated
  // to the right file
  struct CaseTableInfo
  {
    String file_name;
    CaseTable* case_table = nullptr;
  };

  // List of CaseTable for traction boundary conditions
  UniqueArray<CaseTableInfo> m_traction_case_table_list;
  UniqueArray<CaseTableInfo> m_double_couple_case_table_list_north;
  UniqueArray<CaseTableInfo> m_double_couple_case_table_list_south;
  UniqueArray<CaseTableInfo> m_double_couple_case_table_list_east;
  UniqueArray<CaseTableInfo> m_double_couple_case_table_list_west;

 private:

  void _doStationarySolve();
  void _getParameters();
  void _updateVariables();
  void _updateTime();
  void _assembleBilinearOperatorTRIA3();
  void _assembleBilinearOperatorEDGE2();
  void _solve();
  void _assembleLinearOperator();
  void _applyDirichletBoundaryConditions();
  void _applyDoubleCoupleLinear();
  void _applyDoubleCoupleBilinear();
  void _checkResultFile();
  void _readCaseTables();
  FixedMatrix<4, 4> _computeElementMatrixEDGE2(Face face);
  FixedMatrix<6, 6> _computeElementMatrixTRIA3(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
  Real _computeEdgeLength2(Face face);
  Real2 _computeDxDyOfRealTRIA3(Cell cell);
  Real2 _computeEdgeNormal2(Face face);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "Module Fem COMPUTE";
  info() << " \n\n***[WIP] this is module is not working yet please dont trust the results***[\n\n";

  // Stop code after computations
  if (t >= tmax)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  info() << "Time iteration at t : " << t << " (s) ";

  // Set if we want to keep the matrix structure between calls.
  // The matrix has to have the same structure (same structure for non-zero)
  bool keep_struct = true;
  if (m_linear_system.isInitialized() && keep_struct){
    m_linear_system.clearValues();
  }
  else{
    m_linear_system.reset();
    m_linear_system.setLinearSystemFactory(options()->linearSystem());
    m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");
  }

  info() << " \n\n***[WIP] this is module is not working yet please dont trust the results***[\n\n";

  _doStationarySolve();

  info() << " \n\n***[WIP] this is module is not working yet please dont trust the results***[\n\n";

  _updateVariables();

  info() << " \n\n***[WIP] this is module is not working yet please dont trust the results***[\n\n";

  _updateTime();

  info() << " \n\n***[WIP] this is module is not working yet please dont trust the results***[\n\n";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "Module Fem INIT";

  m_dofs_on_nodes.initialize(mesh(), 2);

  _applyDirichletBoundaryConditions();

  // # get parameters
  _getParameters();

  t    = dt;
  tmax = tmax - dt;
  m_global_deltat.assign(dt);

  _readCaseTables();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateTime()
{
  info() << "Update time";
  t += dt;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{

  // Assemble the FEM bilinear operator (LHS - matrix A)
  _assembleBilinearOperatorTRIA3();

  _assembleBilinearOperatorEDGE2();

  // Assemble the FEM linear operator (RHS - vector b)
  _assembleLinearOperator();

  // Solve for [u1,u2]
  _solve();

  // Check results TODO
  // _checkResultFile();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getParameters()
{
  info() << "Get material parameters...";

  //--------- time parameters -----------//
  tmax = options()->tmax();                // max time
  dt   = options()->dt();                  // time step

  //--- time discretization parameter ---//
  alpm = options()->alpm();                // time discretization param alpm
  alpf = options()->alpf();                // time discretization param alpf

  //--------- material parameter ---------//
  f1   = options()->f1();                  // body force in X
  f2   = options()->f2();                  // body force in Y
  E    = options()->E();                   // Youngs modulus
  nu   = options()->nu();                  // Poission ratio
  rho  = options()->rho();                 // Density
  cp   = options()->cp();                  // Wave velocity primary
  cs   = options()->cs();                  // Wave velocity secondary

  if( options()->E.isPresent() && options()->nu.isPresent()) {
    mu     = E/(2*(1+nu));                   // lame parameter mu
    lambda = E*nu/((1+nu)*(1-2*nu));         // lame parameter lambda
    cs     = math::sqrt(mu/rho);
    cp     = math::sqrt((lambda+(2.*mu))/rho) ;
  }

  if( options()->mu.isPresent() && options()->lambda.isPresent()) {
    mu     = options()->mu;
    lambda = options()->lambda;
    cs     = math::sqrt(mu/rho);
    cp     = math::sqrt((lambda+(2.*mu))/rho) ;
  }

  if( (options()->cp.isPresent()) && (options()->cs.isPresent()) ) {
    mu     =  cs*cs*rho;
    lambda =  cp*cp*rho - 2*mu;
  }

  mu2 =  mu*2;                             // lame parameter mu * 2

  //----- time discretization Newmark-Beta or Generalized-alpha  -----//
  if (options()->timeDiscretization == "Newmark-beta") {

    info() << "Apply time discretization via Newmark-beta ";

    gamma = 0.5;
    beta  = (1./4.)*(gamma+0.5)*(gamma+0.5)  ;

    c0 =   rho/(beta*dt*dt)                                ;
    c1 =   lambda                                          ;
    c2 =   2.*mu                                           ;
    c3 =   rho/(beta*dt)                                   ;
    c4 =   rho*(1./2./beta -1.)                            ;
    c5 =   0.                                              ;
    c6 =   0.                                              ;
    c7 =   rho*gamma/beta/dt                               ;
    c8 =   rho*(1.-gamma/beta)                             ;
    c9 =   rho*dt*(1.-gamma/(2.*beta))                     ;

    }

  else if (options()->timeDiscretization == "Generalized-alpha") {

    info() << "Apply time discretization via Generalized-alpha ";

    gamma = 0.5 + alpf - alpm                ;
    beta  = (1./4.)*(gamma+0.5)*(gamma+0.5)  ;

    c0 =   rho*(1.-alpm)/(beta*dt*dt)                      ;
    c1 =   lambda*(1.-alpf)                                ;
    c2 =   2.*mu*(1.-alpf)                                 ;
    c3 =   rho*(1.-alpm)/(beta*dt)                         ;
    c4 =   rho*((1.-alpm)/2./beta -1.)                     ;
    c5 =   lambda*alpf                                     ;
    c6 =   2*mu*alpf                                       ;
    c7 =   rho*(1.-alpf)*gamma/beta/dt                     ;
    c8 =   rho*(1.-gamma*(1-alpf)/beta)                    ;
    c9 =   rho*(1.-alpf)*dt*(1.-gamma/(2.*beta))           ;

    ARCANE_FATAL("Only Newmark-beta works for time-discretization Generalized-alpha WIP ");

    }

  else {

    ARCANE_FATAL("Only Newmark-beta | Generalized-alpha are supported for time-discretization ");

    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_readCaseTables()
{
  IParallelMng* pm = subDomain()->parallelMng();
  for (const auto& bs : options()->tractionBoundaryCondition()) {
    CaseTable* case_table = nullptr;
    String file_name;
    if(bs->tractionInputFile.isPresent()){
      file_name = bs->tractionInputFile();
      case_table = readFileAsCaseTable(pm, file_name, 3);
    }
    m_traction_case_table_list.add(CaseTableInfo{file_name,case_table});
  }

  for (const auto& bs : options()->doubleCouple()) {
    CaseTable* case_table_north = nullptr;
    CaseTable* case_table_south = nullptr;
    CaseTable* case_table_east  = nullptr;
    CaseTable* case_table_west = nullptr;

    String file_name_north;
    String file_name_south;
    String file_name_east;
    String file_name_west;

    if(bs->doubleCoupleInputFile.isPresent()){

      file_name_north = bs->doubleCoupleInputFile()+"_north.txt";
      file_name_south = bs->doubleCoupleInputFile()+"_south.txt";
      file_name_east  = bs->doubleCoupleInputFile()+"_east.txt";
      file_name_west  = bs->doubleCoupleInputFile()+"_west.txt";

      case_table_north = readFileAsCaseTable(pm, file_name_north, 3);
      case_table_south = readFileAsCaseTable(pm, file_name_south, 3);
      case_table_east  = readFileAsCaseTable(pm, file_name_east, 3);
      case_table_west  = readFileAsCaseTable(pm, file_name_west, 3);
    }

    m_double_couple_case_table_list_north.add(CaseTableInfo{file_name_north,case_table_north});
    m_double_couple_case_table_list_south.add(CaseTableInfo{file_name_south,case_table_south});
    m_double_couple_case_table_list_east.add(CaseTableInfo{file_name_east,case_table_east});
    m_double_couple_case_table_list_west.add(CaseTableInfo{file_name_west,case_table_west});
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_applyDirichletBoundaryConditions()
{

  info() << "Apply boundary conditions";

  // Handle all the Dirichlet boundary conditions.
  // In the 'arc' file, there are in the following format:
  //   <dirichlet-boundary-condition>
  //   <surface>Left</surface>
  //   <u1>1.0</u1>
  //   <u1>0.0</u2>
  // </dirichlet-boundary-condition>

  for (const auto& bs : options()->dirichletBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real u1_val = bs->u1();
    Real u2_val = bs->u2();

    if( bs->u1.isPresent() && bs->u2.isPresent()) {
      info() << "Apply Dirichlet boundary condition surface=" << group.name() << " u1= " << u1_val << " u2= " << u2_val;
      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          m_dU[node].x = u1_val;
          m_dU[node].y = u2_val;
          m_u1_fixed[node] = true;
          m_u2_fixed[node] = true;
        }
      }
      continue;
    }

    if(bs->u1.isPresent()) {
      info() << "Apply Dirichlet boundary condition surface=" << group.name() << " u1=" << u1_val;
      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          m_dU[node].x = u1_val;
          m_u1_fixed[node] = true;
        }
      }
      continue;
    }

    if(bs->u2.isPresent()) {
      info() << "Apply Dirichlet boundary condition surface=" << group.name() << " u2=" << u2_val;
      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          m_dU[node].y = u2_val;
          m_u2_fixed[node] = true;
        }
      }
      continue;
    }
  }

  // Handle all the Dirichlet point conditions.
  // In the 'arc' file, there are in the following format:
  //   <dirichlet-point-condition>
  //   <surface>Left</surface>
  //   <u1>1.0</u1>
  //   <u1>0.0</u2>
  // </dirichlet-point-condition>

  for (const auto& bs : options()->dirichletPointCondition()) {
    NodeGroup group = bs->node();
    Real u1_val = bs->u1();
    Real u2_val = bs->u2();

    if( bs->u1.isPresent() && bs->u2.isPresent()) {
      info() << "Apply Dirichlet point condition on node=" << group.name() << " u1= " << u1_val << " u2= " << u2_val;
      ENUMERATE_ (Node, inode, group) {
        Node node = *inode;
        m_dU[node].x = u1_val;
        m_dU[node].y = u2_val;
        m_u1_fixed[node] = true;
        m_u2_fixed[node] = true;
      }
      continue;
    }

    if(bs->u1.isPresent()) {
      info() << "Apply Dirichlet point condition on node=" << group.name() << " u1=" << u1_val;
      ENUMERATE_ (Node, inode, group) {
        Node node = *inode;
        m_dU[node].x = u1_val;
        m_u1_fixed[node] = true;
      }
      continue;
    }

    if(bs->u2.isPresent()) {
      info() << "Apply Dirichlet point condition on node=" << group.name() << " u2=" << u2_val;
      ENUMERATE_ (Node, inode, group) {
        Node node = *inode;
        m_dU[node].y = u2_val;
        m_u2_fixed[node] = true;
      }
      continue;
    }
  }

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateVariables()
{
  // Note at this stage we already have calculated dU
  Real alocX;
  Real alocY;

  VariableDoFReal& dof_u(m_linear_system.solutionVariable());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Node, inode, allNodes()) {
    Node node = *inode;

    alocX = (m_dU[node].x - m_U[node].x - dt*m_V[node].x)/beta/(dt*dt)
                  - (1.-2.*beta)/2./beta*m_A[node].x;
    alocY = (m_dU[node].y - m_U[node].y - dt*m_V[node].y)/beta/(dt*dt)
                  - (1.-2.*beta)/2./beta*m_A[node].y;

    m_V[node].x = m_V[node].x + dt*((1.-gamma)*m_A[node].x + gamma*alocX);
    m_V[node].y = m_V[node].y + dt*((1.-gamma)*m_A[node].y + gamma*alocY);

    m_A[node].x = alocX;
    m_A[node].y = alocY;

    m_U[node].x = m_dU[node].x;
    m_U[node].y = m_dU[node].y;
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
  info() << "Assembly of FEM linear operator ";

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    Real Penalty = options()->penalty();        // 1.0e30 is the default

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u1_fixed[node_id]) {
        DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);
        m_linear_system.matrixSetValue(dof_id1, dof_id1, Penalty);
        {
          Real u1_dirichlet = Penalty * m_dU[node_id].x;
          rhs_values[dof_id1] = u1_dirichlet;
        }
      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);
        m_linear_system.matrixSetValue(dof_id2, dof_id2, Penalty);
        {
          Real u2_dirichlet = Penalty * m_dU[node_id].y;
          rhs_values[dof_id2] = u2_dirichlet;
        }
      }
    }
  }else if (options()->enforceDirichletMethod() == "WeakPenalty") {

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    Real Penalty = options()->penalty();        // 1.0e30 is the default

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u1_fixed[node_id]) {
        DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);
        m_linear_system.matrixAddValue(dof_id1, dof_id1, Penalty);
        {
          Real u1_dirichlet = Penalty * m_dU[node_id].x;
          rhs_values[dof_id1] = u1_dirichlet;
        }
      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);
        m_linear_system.matrixAddValue(dof_id2, dof_id2, Penalty);
        {
          Real u2_dirichlet = Penalty * m_dU[node_id].y;
          rhs_values[dof_id2] = u2_dirichlet;
        }
      }
    }
  }else if (options()->enforceDirichletMethod() == "RowElimination") {

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u1_fixed[node_id]) {
        DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);

        Real u1_dirichlet = m_dU[node_id].x;
        m_linear_system.eliminateRow(dof_id1, u1_dirichlet);

      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);

        Real u2_dirichlet = m_dU[node_id].y;
        m_linear_system.eliminateRow(dof_id2, u2_dirichlet);

      }
    }
  }else if (options()->enforceDirichletMethod() == "RowColumnElimination") {

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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_u1_fixed[node_id]) {
        DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);

        Real u1_dirichlet = m_dU[node_id].x;
        m_linear_system.eliminateRowColumn(dof_id1, u1_dirichlet);

      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);

        Real u2_dirichlet = m_dU[node_id].y;
        m_linear_system.eliminateRowColumn(dof_id2, u2_dirichlet);

      }
    }
  }else {

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " is not supported \n"
           << "enforce-Dirichlet-method only supports:\n"
           << "  - Penalty\n"
           << "  - WeakPenalty\n"
           << "  - RowElimination\n"
           << "  - RowColumnElimination\n";

    ARCANE_FATAL( "Dirichlet boundary conditions were not applied " );
  }
  //----------------------------------------------
  // Body force term assembly
  //----------------------------------------------
  //
  //  $int_{Omega}(f1*v1^h)$
  //  $int_{Omega}(f2*v2^h)$
  //  only for nodes that are non-Dirichlet
  //----------------------------------------------

  if ( options()->f1.isPresent()) {
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

  if ( options()->f2.isPresent()) {
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

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real area = _computeAreaTriangle3(cell);

    Real3 m0 = m_node_coord[cell.nodeId(0)];
    Real3 m1 = m_node_coord[cell.nodeId(1)];
    Real3 m2 = m_node_coord[cell.nodeId(2)];

    Real Uold1 = m_U[cell.nodeId(0)].x + m_U[cell.nodeId(1)].x + m_U[cell.nodeId(2)].x;
    Real Uold2 = m_U[cell.nodeId(0)].y + m_U[cell.nodeId(1)].y + m_U[cell.nodeId(2)].y;

    Real Vold1 = m_V[cell.nodeId(0)].x + m_V[cell.nodeId(1)].x + m_V[cell.nodeId(2)].x;
    Real Vold2 = m_V[cell.nodeId(0)].y + m_V[cell.nodeId(1)].y + m_V[cell.nodeId(2)].y;

    Real Aold1 = m_A[cell.nodeId(0)].x + m_A[cell.nodeId(1)].x + m_A[cell.nodeId(2)].x;
    Real Aold2 = m_A[cell.nodeId(0)].y + m_A[cell.nodeId(1)].y + m_A[cell.nodeId(2)].y;

/*
$$
\int_{\Omega}(
                    (U \cdot v) c_0
                  + (V \cdot v) c_3
                  + (A \cdot v) c_4
               )
$$
*/

    // Info: We could also use the following logic
    //    cell.node(i).isOwn();
    //    cell.node(i);
    //    Node node = cell.node(i);
    // Then we can loop 1:3
    int i = 0;
    for (Node node : cell.nodes()) {
      if (node.isOwn()) {
        DoFLocalId dof_id1 = node_dof.dofId(node, 0);
        DoFLocalId dof_id2 = node_dof.dofId(node, 1);
        
        if (!(m_u1_fixed[node]))
        rhs_values[dof_id1] +=   (Uold1 + m_U[node].x) * (area / 12.) * c0
                               + (Vold1 + m_V[node].x) * (area / 12.) * c3
                               + (Aold1 + m_A[node].x) * (area / 12.) * c4   // TODO add c5 and c6 contribution for Galpha
                               ;

        if (!(m_u2_fixed[node]))
        rhs_values[dof_id2] +=   (Uold2 + m_U[node].y) * (area / 12.) * c0
                               + (Vold2 + m_V[node].y) * (area / 12.) * c3
                               + (Aold2 + m_A[node].y) * (area / 12.) * c4  // TODO add c5 and c6 contribution for Galpha
                               ;
      }
      i++;
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

  // Index of the boundary condition. Needed to associate a CaseTable
  Int32 boundary_condition_index = 0;

  for (const auto& bs : options()->tractionBoundaryCondition()) {

    FaceGroup group = bs->surface();

    const CaseTableInfo& case_table_info = m_traction_case_table_list[boundary_condition_index];
    ++boundary_condition_index;

    Real3 trac; // traction in x, y and z


    if (bs->tractionInputFile.isPresent()){

      String file_name = bs->tractionInputFile();
      info() << "Applying traction boundary conditions for surface "<< group.name()
             << " via CaseTable" <<  file_name;
      CaseTable* inn = case_table_info.case_table;

      if (!inn)
        ARCANE_FATAL("CaseTable is null. Maybe there is a missing call to _readCaseTables()");
      if (file_name!=case_table_info.file_name)
        ARCANE_FATAL("Incoherent CaseTable. The current CaseTable is associated to file '{0}'",case_table_info.file_name);

      inn->value(t, trac);


      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real length = _computeEdgeLength2(face);
        for (Node node : iface->nodes()) {
          if (!(m_u1_fixed[node]) && node.isOwn()) {
            DoFLocalId dof_id1 = node_dof.dofId(node, 0);

            rhs_values[dof_id1] += trac.x * length / 2.;
          }
          if (!(m_u2_fixed[node]) && node.isOwn()) {
            DoFLocalId dof_id2 = node_dof.dofId(node, 1);
            rhs_values[dof_id2] += trac.y * length / 2.;

          }
        }
      }

      continue;
    }
    else {

      info() << "Applying constant traction boundary conditions for surface "<< group.name();

      trac.x = bs->t1();
      trac.y = bs->t2();

      if( bs->t1.isPresent() && bs->t2.isPresent()) {
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real length = _computeEdgeLength2(face);
          for (Node node : iface->nodes()) {
            if (!(m_u1_fixed[node]) && node.isOwn()) {
              DoFLocalId dof_id1 = node_dof.dofId(node, 0);
              rhs_values[dof_id1] += trac.x * length / 2.;
            }
            if (!(m_u2_fixed[node]) && node.isOwn()) {
              DoFLocalId dof_id2 = node_dof.dofId(node, 1);
              rhs_values[dof_id2] += trac.y * length / 2.;
            }
          }
        }
        continue;
      }

      if( bs->t1.isPresent()) {
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
        Real length = _computeEdgeLength2(face);
          for (Node node : iface->nodes()) {
            if (!(m_u1_fixed[node]) && node.isOwn()) {
              DoFLocalId dof_id1 = node_dof.dofId(node, 0);
              rhs_values[dof_id1] += trac.x * length / 2.;
            }
          }
        }
        continue;
      }

      if( bs->t2.isPresent()) {
        ENUMERATE_ (Face, iface, group) {
          Face face = *iface;
          Real length = _computeEdgeLength2(face);
          for (Node node : iface->nodes()) {
            if (!(m_u2_fixed[node]) && node.isOwn()) {
              DoFLocalId dof_id2 = node_dof.dofId(node, 1);
              rhs_values[dof_id2] += trac.y * length / 2. ;
            }
          }
        }
        continue;
      }

    }
  }

  //----------------------------------------------
  // Paraxial term assembly
  //----------------------------------------------
  //
  //  $int_{dOmega_P}  c7 (.....)
  //  only for noded that are non-Dirichlet
  //----------------------------------------------

  for (const auto& bs : options()->paraxialBoundaryCondition()) {
    FaceGroup group = bs->surface();

    info() << "Applying constant paraxial boundary conditions for surface "<< group.name();

    ENUMERATE_ (Face, iface, group) {
      Face face = *iface;

      Real  length = _computeEdgeLength2(face);
      Real2 Normal = _computeEdgeNormal2(face);

      for (Node node : iface->nodes()) {
        if (node.isOwn()) {
          if (!(m_u1_fixed[node])) {
            DoFLocalId dof_id1 = node_dof.dofId(node, 0);
            rhs_values[dof_id1] += (  c7*( cp*( Normal.x*Normal.x*m_U[node].x + Normal.x*Normal.y*m_U[node].y ) +
                                           cs*( Normal.y*Normal.y*m_U[node].x - Normal.x*Normal.y*m_U[node].y )
                                         )
                                    - c8*( cp*( Normal.x*Normal.x*m_V[node].x + Normal.x*Normal.y*m_V[node].y ) +
                                           cs*( Normal.y*Normal.y*m_V[node].x - Normal.x*Normal.y*m_V[node].y )
                                         )
                                    - c9*( cp*( Normal.x*Normal.x*m_A[node].x + Normal.x*Normal.y*m_A[node].y ) +
                                           cs*( Normal.y*Normal.y*m_A[node].x - Normal.x*Normal.y*m_A[node].y )
                                         )
                                    ) * length / 2.;
          }
          if (!(m_u2_fixed[node])) {
            DoFLocalId dof_id2 = node_dof.dofId(node, 1);
            rhs_values[dof_id2] += (  c7*( cp*( Normal.x*Normal.y*m_U[node].x + Normal.y*Normal.y*m_U[node].y ) +
                                           cs*(-Normal.x*Normal.y*m_U[node].x + Normal.x*Normal.x*m_U[node].y )
                                         )
                                    - c8*( cp*( Normal.x*Normal.y*m_V[node].x + Normal.y*Normal.y*m_V[node].y ) +
                                           cs*(-Normal.x*Normal.y*m_V[node].x + Normal.x*Normal.x*m_V[node].y )
                                         )
                                    - c9*( cp*( Normal.x*Normal.y*m_A[node].x + Normal.y*Normal.y*m_A[node].y ) +
                                           cs*(-Normal.x*Normal.y*m_A[node].x + Normal.x*Normal.x*m_A[node].y )
                                         )
                                   ) * length / 2.;
          }
        }
      }
    }
  }

  //----------------------------------------------
  // Double-couple term assembly
  //----------------------------------------------
  //----------------------------------------------

  // Index of the boundary condition. Needed to associate a CaseTable
  Int32 boundary_condition_index_dc = 0;

  for (const auto& bs : options()->doubleCouple()) {


    const CaseTableInfo& case_table_dc_info_north = m_double_couple_case_table_list_north[boundary_condition_index_dc];
    const CaseTableInfo& case_table_dc_info_south = m_double_couple_case_table_list_south[boundary_condition_index_dc];
    const CaseTableInfo& case_table_dc_info_east  = m_double_couple_case_table_list_east[boundary_condition_index_dc];
    const CaseTableInfo& case_table_dc_info_west  = m_double_couple_case_table_list_west[boundary_condition_index_dc];
    ++boundary_condition_index_dc;

      Real3 trac_north; // values in x, y and z
      Real3 trac_south; // values in x, y and z
      Real3 trac_east;  // values in x, y and z
      Real3 trac_west;  // values in x, y and z

      String file_name = bs->doubleCoupleInputFile();
      info() << "Applying boundary conditions for surface via CaseTable" <<  file_name;

      CaseTable* inn_north = case_table_dc_info_north.case_table;
      CaseTable* inn_south = case_table_dc_info_south.case_table;
      CaseTable* inn_east  = case_table_dc_info_east.case_table;
      CaseTable* inn_west  = case_table_dc_info_west.case_table;

      //if (!inn_north)
      //  ARCANE_FATAL("CaseTable is null. Maybe there is a missing call to _readCaseTables()");

      //if (file_name!=case_table_dc_info_north.file_name)
      //  ARCANE_FATAL("Incoherent CaseTable. The current CaseTable is associated to file '{0}'",case_table_dc_info_north.file_name);

      inn_north->value(t, trac_north);
      inn_south->value(t, trac_south);
      inn_east->value(t, trac_east);
      inn_west->value(t, trac_west);

      //cout << "ArcFemDebug north "<< trac_north.x <<  "  " << trac_north.y <<  "  " << trac_north.z << endl;
      //cout << "ArcFemDebug south "<< trac_south.x <<  "  " << trac_south.y <<  "  " << trac_south.z << endl;
      //cout << "ArcFemDebug east "<< trac_east.x  <<  "  " << trac_east.y <<  "  " << trac_east.z << endl;
      //cout << "ArcFemDebug west "<< trac_west.x  <<  "  " << trac_west.y <<  "  " << trac_west.z << endl;

    NodeGroup north = bs->northNodeName();
    NodeGroup south = bs->southNodeName();
    NodeGroup east  = bs->eastNodeName();
    NodeGroup west  = bs->westNodeName();



       ENUMERATE_ (Node, inode, north) {
        Node node = *inode;
        DoFLocalId dof_id1 = node_dof.dofId(node, 0);
        rhs_values[dof_id1] = trac_north.x;
      }
       ENUMERATE_ (Node, inode, south) {
        Node node = *inode;
        DoFLocalId dof_id1 = node_dof.dofId(node, 0);
        rhs_values[dof_id1] = trac_south.x;
      }
       ENUMERATE_ (Node, inode, east) {
        Node node = *inode;
        DoFLocalId dof_id2 = node_dof.dofId(node, 1);
        rhs_values[dof_id2] = trac_east.y;
      }
       ENUMERATE_ (Node, inode, west) {
        Node node = *inode;
        DoFLocalId dof_id2 = node_dof.dofId(node, 1);
        rhs_values[dof_id2] = trac_west.y;
      }

  }
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
  return  math::sqrt((m1.x-m0.x)*(m1.x-m0.x) + (m1.y-m0.y)*(m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real2 FemModule::
_computeEdgeNormal2(Face face)
{
  Real3 m0 = m_node_coord[face.nodeId(0)];
  Real3 m1 = m_node_coord[face.nodeId(1)];
  if (!face.isSubDomainBoundaryOutside())
    std::swap(m0,m1);
  Real2 N;
  Real norm_N = math::sqrt( (m1.y - m0.y)*(m1.y - m0.y) + (m1.x - m0.x)*(m1.x - m0.x) );   // for normalizing
  N.x = (m1.y - m0.y)/ norm_N;
  N.y = (m0.x - m1.x)/ norm_N;
  return  N;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<6, 6> FemModule::
_computeElementMatrixTRIA3(Cell cell)
{
  // Get coordiantes of the triangle element  TRI3
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

  Real area = _computeAreaTriangle3(cell);    // calculate area

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  FixedMatrix<1, 6> b_matrix;
  FixedMatrix<6, 1> bT_matrix;
  FixedMatrix<6, 6> int_Omega_i;

  for (Int32 i = 0; i<6; i++)
    for (Int32 j = 0; j<6; j++)
      int_Omega_i(i,j) = 0.;

// -----------------------------------------------------------------------------
//  c1( dx(du1)dx(v1) + dy(du2)dx(v1) + dx(du1)dy(v2) + dy(du2)dy(v2) )
//------------------------------------------------------------------------------


  // dx(u1)dx(v1) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_dxU1dxV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxU1dxV1);

  // dy(u2)dx(v1) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_dyU1dyV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyU1dyV1);

  // dx(u1)dy(v2) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_dxU2dxV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxU2dxV1);

  // dy(u2)dy(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_dyU2dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyU2dyV1);

  // lambda * (.....)
  int_Omega_i.multInPlace(c1);


// -----------------------------------------------------------------------------
//  c2( dx(u1)dx(v1) + dy(u2)dy(v2) + 0.5*(   dy(u1)dy(v1) + dx(u2)dy(v1)
//                                            + dy(u1)dx(v2) + dx(u2)dx(v2) )
//      )
//------------------------------------------------------------------------------

  // mu*dx(u1)dx(v1) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudxU1dxV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_mudxU1dxV1.multInPlace(c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU1dxV1);


  // mu*dy(u2)dy(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudyU2dyV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudyU2dyV2.multInPlace(c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU2dyV2);

  // 0.5*0.5*mu*dy(u1)dy(v1) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.y/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.y/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.y;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.y;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudyU1dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudyU1dyV1.multInPlace(0.5*c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU1dyV1);

  // 0.5*mu*dx(u2)dy(v1) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.x/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.x/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.x/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.y;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.y;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudxU2dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudxU2dyV1.multInPlace(0.5*c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU2dyV1);

  // 0.5*mu*dy(u1)dx(v2) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.y/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.y/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.x;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.x;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudyU1dxV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudyU1dxV2.multInPlace(0.5*c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU1dxV2);

  // 0.5*mu*dx(u2)dx(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.x/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.x/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.x/area;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.x;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.x;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_mudxU2dxV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_mudxU2dxV2.multInPlace(0.5*c2);

  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU2dxV2);

// -----------------------------------------------------------------------------
//   c0(du1v1 + du2v2)
// -----------------------------------------------------------------------------
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = 1.;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = 1.;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = area/3.;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = area/3.;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = area/3.;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5);

  FixedMatrix<6, 6> int_U2V2   = matrixMultiplication(bT_matrix, b_matrix);
  int_U2V2.multInPlace(c0);

  for (Int32 i = 0; i<6; i++)
    int_U2V2(i,i) *= 2.;

  int_Omega_i = matrixAddition( int_Omega_i, int_U2V2);



  // du2v2 //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = 1.;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = 1.;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = 1.;

  b_matrix.multInPlace(0.5);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = area/3.;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = area/3.;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = area/3.;

  bT_matrix.multInPlace(0.5);

  int_U2V2   = matrixMultiplication(bT_matrix, b_matrix);
  int_U2V2.multInPlace(c0);

  for (Int32 i = 0; i<6; i++)
    int_U2V2(i,i) *= 2.;


  int_Omega_i = matrixAddition( int_Omega_i, int_U2V2);

  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<4, 4> FemModule::
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

  Real area = _computeEdgeLength2(face);    // calculate length
  Real2 N   = _computeEdgeNormal2(face);

  FixedMatrix<1, 4> b_matrix;
  FixedMatrix<4, 1> bT_matrix;
  FixedMatrix<4, 4> int_DOmega_i;

  for (Int32 i = 0; i<4; i++)
    for (Int32 j = 0; j<4; j++)
      int_DOmega_i(i,j) = 0.;

  // c7*u1v1 //
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = 1.;
  b_matrix(0, 3) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = area/3;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = area/3;
  bT_matrix(3, 0) = 0.;

  bT_matrix.multInPlace(1.);

  FixedMatrix<4, 4> int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<4; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(c7*(N.x*N.x*cp + N.y*N.y*cs));
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  // c7*u2v2 //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = 1.;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = 1.;

  b_matrix.multInPlace(0.5f);


  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = area/3;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = area/3;

  bT_matrix.multInPlace(1.);

  int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<4; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(c7*(N.y*N.y*cp + N.x*N.x*cs));
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  // c7*u1v2 //
  b_matrix(0, 0) = 1.;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = 1.;
  b_matrix(0, 3) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = area/3;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = area/3;

  bT_matrix.multInPlace(1.);

  int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  int_UV.multInPlace(c7*(N.x*N.y*(cp - cs)) );
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  // c7*u2v1 //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = 1.;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = 1.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = area/3;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = area/3;
  bT_matrix(3, 0) = 0.;

  bT_matrix.multInPlace(1.);

  int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  int_UV.multInPlace( c7*(N.x*N.y*(cp - cs)) );
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  return int_DOmega_i;
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

    auto K_e = _computeElementMatrixTRIA3(cell);  // element stiffness matrix
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
        Real v1 = K_e(2 * n1_index    , 2 * n2_index    );
        Real v2 = K_e(2 * n1_index    , 2 * n2_index + 1);
        Real v3 = K_e(2 * n1_index + 1, 2 * n2_index    );
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
_assembleBilinearOperatorEDGE2()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  for (const auto& bs : options()->paraxialBoundaryCondition()) {
    FaceGroup group = bs->surface();

    ENUMERATE_ (Face, iface, group) {
      Face face = *iface;

      auto K_e = _computeElementMatrixEDGE2(face);  // element stiffness matrix

      Int32 n1_index = 0;
      for (Node node1 : face.nodes() ) {
        Int32 n2_index = 0;
        for (Node node2 : face.nodes()) {
          Real v1 = K_e(2 * n1_index    , 2 * n2_index    );
          Real v2 = K_e(2 * n1_index    , 2 * n2_index + 1);
          Real v3 = K_e(2 * n1_index + 1, 2 * n2_index    );
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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  info() << "Solving Linear system";
  m_linear_system.solve();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real u1_val = dof_u[node_dof.dofId(node, 0)];
      Real u2_val = dof_u[node_dof.dofId(node, 1)];
      Real3 u_disp;
      u_disp.x = u1_val;
      u_disp.y = u2_val;
      u_disp.z = 0.0;
      m_dU[node] = u_disp;
      info() << "Node: " << node.localId() << " U1=" << u1_val << " U2=" << u2_val;
    }
  }

  // Re-Apply boundary conditions because the solver has modified the value
  _applyDirichletBoundaryConditions();  // ************ CHECK

  m_dU.synchronize();
  m_U.synchronize();
  m_V.synchronize();
  m_A.synchronize();

  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    int p = std::cout.precision();
    std::cout.precision(17);
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      std::cout << "U1[" << node.localId() << "][" << node.uniqueId() << "] = "
                << m_dU[node].x << " U2[" << node.localId() << "][" << node.uniqueId() << "] = "
                << m_dU[node].y << "\n";
      //std::cout << "U1[]" << node.uniqueId() << " " << m_u1[node] << "\n";
    }
    std::cout.precision(p);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_checkResultFile()
{
  String filename = options()->resultFile();
  info() << "CheckResultFile filename=" << filename;
  if (filename.empty())
    return;
  const double epsilon = 1.0e-4;
  // TODO
  //Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_U.x, epsilon);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
