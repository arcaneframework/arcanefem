// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2023 */
/*                                                                           */
/* FEM code to test vectorial FE for Elasticity problem.                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>

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
  Real fy;
  Real mu2;
  Real lambda;

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;

 private:

  void _doStationarySolve();
  void _getMaterialParameters();
  void _updateBoundayConditions();
  void _assembleBilinearOperatorTRIA3();
  void _assembleBilinearOperatorQUAD4();
  void _solve();
  void _initBoundaryconditions();
  void _assembleLinearOperator();
  FixedMatrix<6, 6> _computeElementMatrixTRIA3(Cell cell);
  FixedMatrix<4, 4> _computeElementMatrixQUAD4(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
  Real _computeAreaQuad4(Cell cell);
  Real _computeEdgeLength2(Face face);
  void _applyDirichletBoundaryConditions();
  void _checkResultFile();
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

  m_dofs_on_nodes.initialize(mesh(), 2);

  _initBoundaryconditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  // # get material parameters
  _getMaterialParameters();

  // # update BCs
  _updateBoundayConditions();

  // Assemble the FEM bilinear operator (LHS - matrix A)
  if (options()->meshType == "QUAD4")
    _assembleBilinearOperatorQUAD4();
  else
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
  fy   = options()->fy();                  // body force in Y
  E    = options()->E();                   // Youngs modulus
  nu   = options()->nu();                  // Poission ratio
  
  mu2 = ( E/(2*(1+nu)) )*2;                // lame parameter mu * 2
  lambda = E*nu/((1+nu)*(1-2*nu));         // lame parameter lambda
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
          m_u1[node] = u1_val; 
          m_u2[node] = u2_val;
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
          m_u1[node] = u1_val; 
          m_u1_fixed[node] = true;
        }
      }
      continue;
    }

    if(bs->u2.isPresent()) {
      info() << "Apply Dirichlet boundary condition surface=" << group.name() << " u2=" << u2_val;
      ENUMERATE_ (Face, iface, group) {
        for (Node node : iface->nodes()) {
          m_u2[node] = u2_val;
          m_u2_fixed[node] = true;
        }
      }
      continue;
    }
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateBoundayConditions()
{
  info() << "TODO " << A_FUNCINFO;
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
  info() << "Applying Dirichlet boundary condition via  penalty method ";

  // Temporary variable to keep values for the RHS part of the linear system

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  //----------------------------------------------
  // penelty method for assembly of Dirichlet BC
  //----------------------------------------------
  // TODO: 1.0e6 is a user value, moreover we should use something like 1e31
  //----------------------------------------------

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Node, inode, ownNodes()) {
    NodeLocalId node_id = *inode;
    if (m_u1_fixed[node_id]) {
      DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);
      m_linear_system.matrixAddValue(dof_id1, dof_id1, 1.0e30);
      {
        Real u1_dirichlet = 1.0e30 * m_u1[node_id];
        rhs_values[dof_id1] = u1_dirichlet;
      }
    }
    if (m_u2_fixed[node_id]) {
      DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);
      m_linear_system.matrixAddValue(dof_id2, dof_id2, 1.0e30);
      {
        Real u2_dirichlet = 1.0e30 * m_u2[node_id];
        rhs_values[dof_id2] = u2_dirichlet;
      }
    }
  }

  //----------------------------------------------
  // Constant source term assembly
  //----------------------------------------------
  //
  //  $int_{Omega}(fy*v2^h)$
  //  only for noded that are non-Dirichlet
  //----------------------------------------------
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real area = _computeAreaTriangle3(cell);
    for (Node node : cell.nodes()) {
      if (!(m_u1_fixed[node]) && node.isOwn()) {
        DoFLocalId dof_id2 = node_dof.dofId(node, 1);
        rhs_values[dof_id2] += fy * area / 3;
      }
    }
  }

  //----------------------------------------------
  // Constant flux term assembly
  //----------------------------------------------
  //
  //  $int_{dOmega_N}((q.n)*v^h)$
  //  only for noded that are non-Dirichlet
  //  TODO : take flux vector and use normals at boundaries
  //----------------------------------------------
  for (const auto& bs : options()->neumannBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real value = bs->value();
    ENUMERATE_ (Face, iface, group) {
      Face face = *iface;
      Real length = _computeEdgeLength2(face);
      for (Node node : iface->nodes()) {
        if (!(m_u1_fixed[node]) && node.isOwn()) {
          DoFLocalId dof_id1 = node_dof.dofId(node, 0);
          rhs_values[dof_id1] += value * length / 2.;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real FemModule::
_computeAreaQuad4(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  Real3 m3 = m_node_coord[cell.nodeId(3)];
  return 0.5 * (  (m1.x*m2.y + m2.x*m3.y + m3.x*m0.y + m0.x*m1.y)
                 -(m2.x*m1.y + m3.x*m2.y + m0.x*m3.y + m1.x*m0.y) );
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

FixedMatrix<4, 4> FemModule::
_computeElementMatrixQUAD4(Cell cell)
{
  // Get coordiantes of the quadrangular element  QUAD4
  //------------------------------------------------
  //             1 o . . . . o 0
  //               .         .
  //               .         .
  //               .         .
  //             2 o . . . . o 3
  //------------------------------------------------
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  Real3 m3 = m_node_coord[cell.nodeId(3)];

  Real area = _computeAreaQuad4(cell);    // calculate area

  Real2 dPhi0(m2.y - m3.y, m3.x - m2.x);
  Real2 dPhi1(m3.y - m0.y, m0.x - m3.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);
  Real2 dPhi3(m1.y - m2.y, m2.x - m1.x);

  FixedMatrix<2, 4> b_matrix;
  b_matrix(0, 0) = dPhi0.x;
  b_matrix(0, 1) = dPhi1.x;
  b_matrix(0, 2) = dPhi2.x;
  b_matrix(0, 3) = dPhi3.x;

  b_matrix(1, 0) = dPhi0.y;
  b_matrix(1, 1) = dPhi1.y;
  b_matrix(1, 2) = dPhi2.y;
  b_matrix(1, 3) = dPhi3.y;

  b_matrix.multInPlace(1.0 / (2.0 * area));

  FixedMatrix<4, 4> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(area);

  //info() << "Cell=" << cell.localId();
  //std::cout << " int_cdPi_dPj=";
  //int_cdPi_dPj.dump(std::cout);
  //std::cout << "\n";

  return int_cdPi_dPj;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorQUAD4()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    if (cell.type() != IT_Quad4)
      ARCANE_FATAL("Only Quad4 cell type is supported");

    auto K_e = _computeElementMatrixQUAD4(cell);  // element stifness matrix
    // assemble elementary matrix into the global one elementary terms are
    // positionned into K according to the rank of associated  node in the
    // mesh.nodes list and acoording the dof number. Here  for  each  node
    // two dofs exists [u1,u2]. For each TRIA3 there are 3 nodes hence the
    // elementary stifness matrix size is (3*2 x 3*2)=(6x6). We will  fill
    // this below in 4 at a time.
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        //Real v = K_e(n1_index, n2_index);
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
_assembleBilinearOperatorTRIA3()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    if (cell.type() != IT_Triangle3)
      ARCANE_FATAL("Only Triangle3 cell type is supported");

    auto K_e = _computeElementMatrixTRIA3(cell);  // element stifness matrix
    // assemble elementary matrix into the global one elementary terms are
    // positionned into K according to the rank of associated  node in the
    // mesh.nodes list and acoording the dof number. Here  for  each  node
    // two dofs exists [u1,u2]. For each TRIA3 there are 3 nodes hence the
    // elementary stifness matrix size is (3*2 x 3*2)=(6x6). We will  fill
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
_solve()
{
  info() << "Solving Linear system";
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
      m_u1[node] = u1_val;
      m_u2[node] = u2_val;
      Real3 u_disp;  
      u_disp.x = u1_val;
      u_disp.y = u2_val; 
      u_disp.z = 0.0;
      m_U[node] = u_disp;
      info() << "Node: " << node.localId() << " U1=" << u1_val << " U2=" << u2_val;
    }
  }

  m_u1.synchronize();
  m_u2.synchronize();
  m_U.synchronize();

  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    int p = std::cout.precision();
    std::cout.precision(17);
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      std::cout << "U1[" << node.localId() << "][" << node.uniqueId() << "] = "
                << m_u1[node] << " U2[" << node.localId() << "][" << node.uniqueId() << "] = "
                << m_u2[node] << "\n";
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
  Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_u1, epsilon);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
