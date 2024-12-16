// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2024 */
/*                                                                           */
/* FEM code to test vectorial FE for Elasticity problem.                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/core/ItemTypes.h>
#include <arccore/base/ArccoreGlobal.h>

#include "IDoFLinearSystemFactory.h"
#include "Fem_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"
#include "BSRFormat.h"
#include "ArcaneFemFunctionsGpu.h"

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
  , m_bsr_format(mbi.subDomain()->traceMng(), *(mbi.subDomain()->acceleratorMng()->defaultQueue()), m_dofs_on_nodes)
  {
    ICaseMng* cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
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

  Real E;
  Real nu;
  Real f1;
  Real f2;
  Real mu2;
  Real lambda;

  bool m_use_bsr = false;

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;

  BSRFormat m_bsr_format;

 private:

  void _doStationarySolve();
  void _getMaterialParameters();
  void _solve();
  void _initBoundaryconditions();
  void _assembleLinearOperator();
  FixedMatrix<6, 6> _computeElementMatrixTRIA3(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
  Real _computeEdgeLength2(Face face);
  void _applyDirichletBoundaryConditions();
  void _checkResultFile();

 public:

  void _assembleBilinearOperatorTRIA3();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "Module Fem COMPUTE";

  // Stop code after computations
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();
  _doStationarySolve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "Module Fem INIT";

  IMesh* mesh = defaultMesh();
  m_dofs_on_nodes.initialize(mesh, 2);

  bool use_csr_in_linearsystem = options()->linearSystem.serviceName() == "HypreLinearSystem";
  m_bsr_format.initialize(mesh, nbFace(), use_csr_in_linearsystem);
  m_bsr_format.computeSparsity();

  _initBoundaryconditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  // # get material parameters
  _getMaterialParameters();

  // Assemble the FEM bilinear operator (LHS - matrix A)
  _assembleBilinearOperatorTRIA3();

  // Assemble the FEM linear operator (RHS - vector b)
  _assembleLinearOperator();

  // Solve for [u1,u2]
  _solve();

  // Check results
  _checkResultFile();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "Get material parameters...";
  f1   = options()->f1();                  // body force in Y
  f2   = options()->f2();                  // body force in Y
  E    = options()->E();                   // Youngs modulus
  nu   = options()->nu();                  // Poission ratio

  mu2 = ( E/(2*(1+nu)) )*2;                // lame parameter mu * 2
  lambda = E*nu/((1+nu)*(1-2*nu));         // lame parameter lambda
  
  if (options()->bsr) {
    m_use_bsr = true;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_initBoundaryconditions()
{
  info() << "Init boundary conditions...";

  info() << "Apply boundary conditions";
  _applyDirichletBoundaryConditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_applyDirichletBoundaryConditions()
{
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
          m_U[node].x = u1_val;
          m_U[node].y = u2_val;
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
          m_U[node].x = u1_val;
          m_u1_fixed[node] = true;
        }
      }
      continue;
    }

    if(bs->u2.isPresent()) {
      info() << "Apply Dirichlet boundary condition surface=" << group.name() << " u2=" << u2_val;
      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          m_U[node].y = u2_val;
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
        m_U[node].x = u1_val;
        m_U[node].y = u2_val;
        m_u1_fixed[node] = true;
        m_u2_fixed[node] = true;
      }
      continue;
    }

    if(bs->u1.isPresent()) {
      info() << "Apply Dirichlet point condition on node=" << group.name() << " u1=" << u1_val;
      ENUMERATE_ (Node, inode, group) {
        Node node = *inode;
        m_U[node].x = u1_val;
        m_u1_fixed[node] = true;
      }
      continue;
    }

    if(bs->u2.isPresent()) {
      info() << "Apply Dirichlet point condition on node=" << group.name() << " u2=" << u2_val;
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

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    Real Penalty = options()->penalty();        // 1.0e30 is the default

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

        Real u1_dirichlet = m_U[node_id].x;
        m_linear_system.eliminateRow(dof_id1, u1_dirichlet);

      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);

        Real u2_dirichlet = m_U[node_id].y;
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

        Real u1_dirichlet = m_U[node_id].x;
        m_linear_system.eliminateRowColumn(dof_id1, u1_dirichlet);

      }
      if (m_u2_fixed[node_id]) {
        DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);

        Real u2_dirichlet = m_U[node_id].y;
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
  }

  //----------------------------------------------
  // Body force term assembly
  //----------------------------------------------
  //
  //  $int_{Omega}(f1*v1^h)$
  //  $int_{Omega}(f2*v2^h)$
  //  only for noded that are non-Dirichlet
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

    if( bs->t1.isPresent() && bs->t2.isPresent()) {
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


    if( bs->t1.isPresent()) {
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



    if( bs->t2.isPresent()) {
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
//  lambda( dx(u1)dx(v1) + dy(u2)dx(v1) + dx(u1)dy(v2) + dy(u2)dy(v2) ) + u2v2
//------------------------------------------------------------------------------


  // dx(u1)dx(v1) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5f);

  FixedMatrix<6, 6> int_dxU1dxV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxU1dxV1);

  // dy(u2)dx(v1) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5f);

  FixedMatrix<6, 6> int_dyU1dyV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyU1dyV1);



  // dx(u1)dy(v2) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5f);

  FixedMatrix<6, 6> int_dxU2dxV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxU2dxV1);

  // dy(u2)dy(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5f);

  FixedMatrix<6, 6> int_dyU2dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyU2dyV1);

  // lambda * (.....)
  int_Omega_i.multInPlace(lambda);


// -----------------------------------------------------------------------------
//  2*mu( dx(u1)dx(v1) + dy(u2)dy(v2) + 0.5*(   dy(u1)dy(v1) + dx(u2)dy(v1)
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

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5f*mu2);

  FixedMatrix<6, 6> int_mudxU1dxV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU1dxV1);


  // mu*dy(u2)dy(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5f*mu2);

  FixedMatrix<6, 6> int_mudyU2dyV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU2dyV2);

  // 0.5*0.5*mu*dy(u1)dy(v1) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.y/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.y/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.y;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.y;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.25f*mu2);

  FixedMatrix<6, 6> int_mudyU1dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU1dyV1);


  // 0.5*mu*dx(u2)dy(v1) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.x/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.x/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.x/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.y;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.y;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.25f*mu2);

  FixedMatrix<6, 6> int_mudxU2dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU2dyV1);

  // 0.5*mu*dy(u1)dx(v2) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.y/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.y/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.x;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.x;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.25f*mu2);

  FixedMatrix<6, 6> int_mudyU1dxV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU1dxV2);

  // 0.5*mu*dx(u2)dx(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.x/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.x/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.x/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.x;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.x;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.25f*mu2);

  FixedMatrix<6, 6> int_mudxU2dxV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU2dxV2);

  //info() << "Cell=" << cell.localId();
  //std::cout << " int_cdPi_dPj=";
  //int_cdPi_dPj.dump(std::cout);
  //std::cout << "\n";

  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<6, 6> computeElementMatrixTRIA3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord, Real lambda, Real mu2)
{
  // Get coordiantes of the triangle element  TRI3
  //------------------------------------------------
  //                  0 o
  //                   . .
  //                  .   .
  //                 .     .
  //              1 o . . . o 2
  //------------------------------------------------
  Real3 m0 = in_node_coord[cn_cv.nodeId(cell_lid, 0)];
  Real3 m1 = in_node_coord[cn_cv.nodeId(cell_lid, 1)];
  Real3 m2 = in_node_coord[cn_cv.nodeId(cell_lid, 2)];

  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

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
//  lambda( dx(u1)dx(v1) + dy(u2)dx(v1) + dx(u1)dy(v2) + dy(u2)dy(v2) ) + u2v2
//------------------------------------------------------------------------------


  // dx(u1)dx(v1) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5f);

  FixedMatrix<6, 6> int_dxU1dxV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxU1dxV1);

  // dy(u2)dx(v1) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5f);

  FixedMatrix<6, 6> int_dyU1dyV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyU1dyV1);



  // dx(u1)dy(v2) //
  b_matrix(0, 0) = dPhi0.x/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.x/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.x/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5f);

  FixedMatrix<6, 6> int_dxU2dxV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dxU2dxV1);

  // dy(u2)dy(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5f);

  FixedMatrix<6, 6> int_dyU2dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_dyU2dyV1);

  // lambda * (.....)
  int_Omega_i.multInPlace(lambda);


// -----------------------------------------------------------------------------
//  2*mu( dx(u1)dx(v1) + dy(u2)dy(v2) + 0.5*(   dy(u1)dy(v1) + dx(u2)dy(v1)
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

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.x;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.x;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.x;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.5f*mu2);

  FixedMatrix<6, 6> int_mudxU1dxV1 = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU1dxV1);


  // mu*dy(u2)dy(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.y/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.y/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.y/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.y;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.y;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.y;

  bT_matrix.multInPlace(0.5f*mu2);

  FixedMatrix<6, 6> int_mudyU2dyV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU2dyV2);

  // 0.5*0.5*mu*dy(u1)dy(v1) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.y/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.y/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.y;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.y;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.25f*mu2);

  FixedMatrix<6, 6> int_mudyU1dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU1dyV1);


  // 0.5*mu*dx(u2)dy(v1) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.x/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.x/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.x/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = dPhi0.y;
  bT_matrix(1, 0) = 0.;
  bT_matrix(2, 0) = dPhi1.y;
  bT_matrix(3, 0) = 0.;
  bT_matrix(4, 0) = dPhi2.y;
  bT_matrix(5, 0) = 0.;

  bT_matrix.multInPlace(0.25f*mu2);

  FixedMatrix<6, 6> int_mudxU2dyV1  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU2dyV1);

  // 0.5*mu*dy(u1)dx(v2) //
  b_matrix(0, 0) = dPhi0.y/area;
  b_matrix(0, 1) = 0.;
  b_matrix(0, 2) = dPhi1.y/area;
  b_matrix(0, 3) = 0.;
  b_matrix(0, 4) = dPhi2.y/area;
  b_matrix(0, 5) = 0.;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.x;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.x;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.25f*mu2);

  FixedMatrix<6, 6> int_mudyU1dxV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudyU1dxV2);

  // 0.5*mu*dx(u2)dx(v2) //
  b_matrix(0, 0) = 0.;
  b_matrix(0, 1) = dPhi0.x/area;
  b_matrix(0, 2) = 0.;
  b_matrix(0, 3) = dPhi1.x/area;
  b_matrix(0, 4) = 0.;
  b_matrix(0, 5) = dPhi2.x/area;

  b_matrix.multInPlace(0.5f);

  bT_matrix(0, 0) = 0.;
  bT_matrix(1, 0) = dPhi0.x;
  bT_matrix(2, 0) = 0.;
  bT_matrix(3, 0) = dPhi1.x;
  bT_matrix(4, 0) = 0.;
  bT_matrix(5, 0) = dPhi2.x;

  bT_matrix.multInPlace(0.25f*mu2);

  FixedMatrix<6, 6> int_mudxU2dxV2  = matrixMultiplication(bT_matrix, b_matrix);
  int_Omega_i = matrixAddition( int_Omega_i, int_mudxU2dxV2);

  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTRIA3()
{
  if (m_use_bsr) {
    auto m_queue = acceleratorMng()->defaultQueue();

    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(m_queue);
    auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
    auto lambda_copy = lambda;
    auto mu2_copy = mu2;

    m_bsr_format.assembleCellWise([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTRIA3Gpu(cell_lid, cn_cv, in_node_coord, lambda_copy, mu2_copy); });
  }
  else {
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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  info() << "Solving Linear system";
  m_bsr_format.toLinearSystem(m_linear_system);
  m_linear_system.solve();

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
      info() << "Node: " << node.uniqueId() << " " << u1_val << " " << u2_val;
    }
  }

  m_U.synchronize();

  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    int p = std::cout.precision();
    std::cout.precision(17);
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      std::cout << "U1[" << node.localId() << "][" << node.uniqueId() << "] = "
                << m_U[node].x << " U2[" << node.localId() << "][" << node.uniqueId() << "] = "
                << m_U[node].y << "\n";
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
  const double min_value_to_test = 1.0e-16;
  Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_U, epsilon, min_value_to_test);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
