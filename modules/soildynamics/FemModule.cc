﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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
  info() << "Module Fem COMPUTE";

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

  _doStationarySolve();
  _updateVariables();
  _updateTime();

  // At the last time stepp check error
  if (t > tmax + dt - 1e-8){
    info() << "Perfroming check";
    _checkResultFile();
  }

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "Module Fem INIT";

  m_dofs_on_nodes.initialize(mesh(), 2);

  // # get parameters
  _getParameters();

  t    = dt;
  tmax = tmax;
  m_global_deltat.assign(dt);

  _readCaseTables();

  _applyDirichletBoundaryConditions();
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
  if(t<=dt){
    _assembleBilinearOperatorTRIA3();
    _assembleBilinearOperatorEDGE2();
  }

  // Assemble the FEM linear operator (RHS - vector b)
  _assembleLinearOperator();

  // Solve for [u1,u2]
  _solve();

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
    CaseTable* case_table = nullptr;

    if(bs->doubleCoupleInputFile.isPresent()){
      case_table = readFileAsCaseTable(pm, bs->doubleCoupleInputFile(), 1);
    }

    m_double_couple_case_table_list.add(CaseTableInfo{bs->doubleCoupleInputFile(),case_table});
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_applyDirichletBoundaryConditions()
{

  info() << "Apply boundary conditions";

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

      Real f0 = m_U[face.nodeId(0)].x;
      Real f1 = m_U[face.nodeId(1)].x;
      Real Uold1 = f0 + f1;

      f0 = m_U[face.nodeId(0)].y;
      f1 = m_U[face.nodeId(1)].y;
      Real Uold2 = f0 + f1;

      f0 = m_V[face.nodeId(0)].x;
      f1 = m_V[face.nodeId(1)].x;
      Real Vold1 = f0 + f1;

      f0 = m_V[face.nodeId(0)].y;
      f1 = m_V[face.nodeId(1)].y;
      Real Vold2 = f0 + f1;

      f0 = m_A[face.nodeId(0)].x;
      f1 = m_A[face.nodeId(1)].x;
      Real Aold1 = f0 + f1;

      f0 = m_A[face.nodeId(0)].y;
      f1 = m_A[face.nodeId(1)].y;
      Real Aold2 = f0 + f1;

      for (Node node : iface->nodes()) {
        if (node.isOwn()) {
          if (!(m_u1_fixed[node])) {
            DoFLocalId dof_id1 = node_dof.dofId(node, 0);
             rhs_values[dof_id1] += (  c7*( cp*( Normal.x*Normal.x*(Uold1+m_U[node].x) + Normal.x*Normal.y*(Uold2+m_U[node].y) ) +
                                           cs*( Normal.y*Normal.y*(Uold1+m_U[node].x) - Normal.x*Normal.y*(Uold2+m_U[node].y) )
                                         )
                                    - c8*( cp*( Normal.x*Normal.x*(Vold1+m_V[node].x) + Normal.x*Normal.y*(Vold2+m_V[node].y) ) +
                                           cs*( Normal.y*Normal.y*(Vold1+m_V[node].x) - Normal.x*Normal.y*(Vold2+m_V[node].y) )
                                         )
                                    - c9*( cp*( Normal.x*Normal.x*(Aold1+m_A[node].x) + Normal.x*Normal.y*(Aold2+m_A[node].y) ) +
                                           cs*( Normal.y*Normal.y*(Aold1+m_A[node].x) - Normal.x*Normal.y*(Aold2+m_A[node].y) )
                                         )
                                    ) * length / 6.;
          }
          if (!(m_u2_fixed[node])) {
            DoFLocalId dof_id2 = node_dof.dofId(node, 1);
             rhs_values[dof_id2] += (  c7*( cp*( Normal.x*Normal.y*(Uold1+m_U[node].x) + Normal.y*Normal.y*(Uold2+m_U[node].y)) +
                                           cs*(-Normal.x*Normal.y*(Uold1+m_U[node].x) + Normal.x*Normal.x*(Uold2+m_U[node].y) )
                                         )
                                    - c8*( cp*( Normal.x*Normal.y*(Vold1+m_V[node].x) + Normal.y*Normal.y*(Vold2+m_V[node].y) ) +
                                           cs*(-Normal.x*Normal.y*(Vold1+m_V[node].x) + Normal.x*Normal.x*(Vold2+m_V[node].y) )
                                         )
                                    - c9*( cp*( Normal.x*Normal.y*(Aold1+m_A[node].x) + Normal.y*Normal.y*(Aold2+m_A[node].y) ) +
                                           cs*(-Normal.x*Normal.y*(Aold1+m_A[node].x) + Normal.x*Normal.x*(Aold2+m_A[node].y) )
                                         )
                                   ) * length / 6.;
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

    const CaseTableInfo& case_table_dc_info = m_double_couple_case_table_list[boundary_condition_index_dc];

    ++boundary_condition_index_dc;

    Real dc_force; // double-couple force

    String file_name = bs->doubleCoupleInputFile();
    info() << "Applying boundary conditions for surface via CaseTable" <<  file_name;

    CaseTable* dc_case_table_inn = case_table_dc_info.case_table;

    dc_case_table_inn->value(t, dc_force);

    NodeGroup north = bs->northNodeName();
    NodeGroup south = bs->southNodeName();
    NodeGroup east  = bs->eastNodeName();
    NodeGroup west  = bs->westNodeName();

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
      Real  u1_val = dof_u[node_dof.dofId(node, 0)];
      Real  u2_val = dof_u[node_dof.dofId(node, 1)];
      Real3 u_disp;
      u_disp.x = u1_val;
      u_disp.y = u2_val;
      u_disp.z = 0.;
      m_dU[node] = u_disp;
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
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      //std::cout.precision(17);
      //std::cout << "U1[" << node.localId() << "][" << node.uniqueId() << "] = "
      //          << m_dU[node].x << " U2[" << node.localId() << "][" << node.uniqueId() << "] = "
      //          << m_dU[node].y << "\n";
      info() << "Node: " << node.uniqueId() << " U1=" << m_dU[node].x << " U2=" << m_dU[node].y << " U3=0.0";
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
  Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_dU, epsilon);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
