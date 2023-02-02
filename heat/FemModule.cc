// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2023 */
/*                                                                           */
/* Simple module to test simple FEM mechanism.                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/NumArray.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>

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

  //! Time variables
  Real t   ,
       dt  ,
       tmax;
  //! Temperature
  Real Tinit ,
       Text  ;
  //! Material parameters
  Real lambda ,
       h      ,
       qdot   ;
  //! FEM parameter
  Real ElementNodes;

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnNodes m_dofs_on_nodes;

 private:

  void _initTime();
  void _updateTime();
  void _updateVariables();
  void _initTemperature();
  void _doStationarySolve();
  void _getParameters();
  void _updateBoundayConditions();
  void _assembleBilinearOperatorTRIA3();
  void _assembleBilinearOperatorQUAD4();
  void _assembleBilinearOperatorEDGE2();
  void _solve();
  void _initBoundaryconditions();
  void _assembleLinearOperator();
  FixedMatrix<2, 2> _computeElementMatrixEDGE2(Face face);
  FixedMatrix<3, 3> _computeElementMatrixTRIA3(Cell cell);
  FixedMatrix<4, 4> _computeElementMatrixQUAD4(Cell cell);
  Real  _computeDxOfRealTRIA3(Cell cell);
  Real  _computeDyOfRealTRIA3(Cell cell);
  Real2 _computeDxDyOfRealTRIA3(Cell cell);
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
  if (t >= tmax)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();
  _doStationarySolve();
  _updateVariables();
  _updateTime();

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "Module Fem INIT";

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();

  _initBoundaryconditions();    // initialize boundary conditions
  _initTime();                  // initialize time
  _getParameters();             // get material parameters
  _initTemperature();           // initialize temperature
  m_global_deltat.assign(dt);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_initTime()
{
  info() << "Initiate time";

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
  info() << "Update time";

  t += dt;
  info() << "Time t is :" << t << " (s)";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateVariables()
{
  info() << "Update FEM variables";

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
  info() << "Init Temperature";

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

  // # update BCs
  _updateBoundayConditions();

  // Assemble the FEM bilinear operator (LHS - matrix A)
  if (options()->meshType == "QUAD4")
    _assembleBilinearOperatorQUAD4();
  else
    _assembleBilinearOperatorTRIA3();

  _assembleBilinearOperatorEDGE2();

  // Assemble the FEM linear operator (RHS - vector b)
  _assembleLinearOperator();

  // # T=linalg.solve(K,RHS)
  _solve();

  // Check results
  _checkResultFile();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getParameters()
{
  info() << "Get material parameters...";
  lambda = options()->lambda();
  qdot   = options()->qdot();
  ElementNodes = 3.;

  if (options()->meshType == "QUAD4")
    ElementNodes = 4.;

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
  //   <surface>Haut</surface>
  //   <value>21.0</value>
  // </dirichlet-boundary-condition>

  for (const auto& bs : options()->dirichletBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real value = bs->value();
    info() << "Apply Dirichlet boundary condition surface=" << group.name() << " v=" << value;
    ENUMERATE_ (Face, iface, group) {
      for (Node node : iface->nodes()) {
        m_node_temperature[node] = value;
        m_node_is_temperature_fixed[node] = true;
      }
    }
  }

  // Handle all Dirichlet point boundary condition
  // In the 'arc' file, these are recognized by
  //     <dirichlet-point-condition>
  //    <node>botLeftCorner</node>
  //    <value>20.0</value>
  //  </dirichlet-point-condition>
  for (const auto& bs : options()->dirichletPointCondition()) {
    NodeGroup group = bs->node();
    Real value = bs->value();
    info() << "Apply Dirichlet point condition node=" << group.name() << " v=" << value;
    ENUMERATE_ (Node, inode, group) {
      Node node = *inode;
      m_node_temperature[node] = value;
      m_node_is_temperature_fixed[node] = true;
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
//  - TODO: external fluxes
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
      if (m_node_is_temperature_fixed[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_linear_system.matrixSetValue(dof_id, dof_id, Penalty);
        Real temperature = Penalty * m_node_temperature[node_id];
        rhs_values[dof_id] = temperature;
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
      if (m_node_is_temperature_fixed[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_linear_system.matrixAddValue(dof_id, dof_id, Penalty);
        Real temperature = Penalty * m_node_temperature[node_id];
        rhs_values[dof_id] = temperature;
      }
    }
  }else if (options()->enforceDirichletMethod() == "RowElimination") {

    //----------------------------------------------
    // Row elimination method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'I' be the set of DOF for which  Dirichlet condition needs to be applied
    //
    //  to apply the Dirichlet on 'i'th DOF
    //  - For LHS matrix A the row terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j
    //           a_{i,j} = 1.  : i==j
    //----------------------------------------------

    info() << "Applying Dirichlet boundary condition via "
           << options()->enforceDirichletMethod() << " method ";

    ENUMERATE_ (Node, inode, ownNodes()) {
      NodeLocalId node_id = *inode;
      if (m_node_is_temperature_fixed[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_linear_system.matrixEliminateRow(dof_id);
        m_linear_system.matrixSetValue(dof_id, dof_id, 1.);
        Real temperature = m_node_temperature[node_id];
        rhs_values[dof_id] = temperature;
      }
    }
    // TODO
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

    // TODO
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
  // Constant source term assembly
  //----------------------------------------------
  //
  //  $int_{Omega}(qdot*v^h)$
  //  only for nodes that are non-Dirichlet
  //----------------------------------------------
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real area = _computeAreaTriangle3(cell);
    for (Node node : cell.nodes()) {
      if (!(m_node_is_temperature_fixed[node]) && node.isOwn())
        rhs_values[node_dof.dofId(node, 0)] += (m_node_temperature_old[node]/dt) * area / ElementNodes;
    }
  }

  //----------------------------------------------
  // Constant flux term assembly
  //----------------------------------------------
  //
  //  $int_{dOmega_N}((q.n)*v^h)$
  //  only for nodes that are non-Dirichlet
  //  TODO : take flux vector and use normal at boundaries
  //----------------------------------------------
  for (const auto& bs : options()->neumannBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real value = bs->value();
    ENUMERATE_ (Face, iface, group) {
      Face face = *iface;
      Real length = _computeEdgeLength2(face);
      for (Node node : iface->nodes()) {
        if (!(m_node_is_temperature_fixed[node]) && node.isOwn())
          rhs_values[node_dof.dofId(node, 0)] += value * length / 2.;
      }
    }
  }

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
      Real length = _computeEdgeLength2(face);
      for (Node node : iface->nodes()) {
        if (!(m_node_is_temperature_fixed[node]) && node.isOwn())
          rhs_values[node_dof.dofId(node, 0)] += h * Text * length / 2.;
      }
    }
  }

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

FixedMatrix<2, 2> FemModule::
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

  Real2 dPhi0(m0.y - m1.y, m1.x - m0.x);
  Real2 dPhi1(m1.y - m0.y, m0.x - m1.x);

  FixedMatrix<1, 2> b_matrix;
  FixedMatrix<2, 1> bT_matrix;
  FixedMatrix<2, 2> int_DOmega_i;

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

  FixedMatrix<2, 2> int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<2; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(h);
  int_DOmega_i = matrixAddition( int_DOmega_i, int_UV);

  return int_DOmega_i;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<3, 3> FemModule::
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

  Real area = _computeAreaTriangle3(cell);    // calculate area

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  FixedMatrix<1, 3> b_matrix;
  FixedMatrix<3, 1> bT_matrix;
  FixedMatrix<3, 3> int_Omega_i;

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

  FixedMatrix<3, 3> int_dxUdxV = matrixMultiplication(bT_matrix, b_matrix);
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

  FixedMatrix<3, 3> int_dyUdyV = matrixMultiplication(bT_matrix, b_matrix);
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

  FixedMatrix<3, 3> int_UV   = matrixMultiplication(bT_matrix, b_matrix);

  for (Int32 i = 0; i<3; i++)
    int_UV(i,i) *= 2.;

  int_UV.multInPlace(1./dt);

  int_Omega_i = matrixAddition( int_Omega_i, int_UV);

  return int_Omega_i;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<4, 4> FemModule::
_computeElementMatrixQUAD4(Cell cell)
{
  // Get coordinates of the quadrangular element  QUAD4
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
  int_cdPi_dPj.multInPlace(area * lambda);

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

    lambda = m_cell_lambda[cell];                 // lambda is always considered cell constant
    auto K_e = _computeElementMatrixQUAD4(cell);  // element stiffness matrix
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
  m_linear_system.solve();

  // Re-Apply boundary conditions because the solver has modified the value
  // of node_temperature on all nodes
  _applyDirichletBoundaryConditions();

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

        //m_dx_node_temperature[cell] = _computeDyOfRealTRIA3(cell);

        Real2 DX = _computeDxDyOfRealTRIA3(cell);
        m_flux[cell].x = -m_cell_lambda[cell] * DX.x;
        m_flux[cell].y = -m_cell_lambda[cell] * DX.y;
        m_flux[cell].z = 0.;
        /*
        Real3 x0 = m_node_coord[cell.nodeId(0)];
        Real3 x1 = m_node_coord[cell.nodeId(1)];
        Real3 x2 = m_node_coord[cell.nodeId(2)];

        Real f0 = m_node_temperature[cell.nodeId(0)];
        Real f1 = m_node_temperature[cell.nodeId(1)];
        Real f2 = m_node_temperature[cell.nodeId(2)];

        // Using Cramer's rule  det (adj (A)) / det (A)
        m_dx_node_temperature[cell]  =    f0*(x1.y - x2.y) - x0.y*(f1 - f2) + (f1*x2.y - f2*x1.y);
        m_dx_node_temperature[cell] /=  ( x0.x*(x1.y - x2.y) - x0.y*(x1.x - x2.x) + (x1.x*x2.y - x2.x*x1.y) );
        */
      }

      m_flux.synchronize();
    }

  }
  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "T[" << node.localId() << "][" << node.uniqueId() << "] = "
             << m_node_temperature[node];
      //info() << "T[]" << node.uniqueId() << " "
      //       << m_node_temperature[node];
    }

/*
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;

      info() << "DX[" << cell.localId() << "] = "
             << m_dx_node_temperature[cell];
    }
*/
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
  checkNodeResultFile(traceMng(), filename, m_node_temperature, epsilon);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
