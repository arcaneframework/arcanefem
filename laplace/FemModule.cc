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
#include <arcane/utils/CommandLineArguments.h>
#include <arcane/utils/StringList.h>

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

  Real ElementNodes;

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
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
  void _applyDirichletBoundaryConditions();
  void _checkResultFile();
  FixedMatrix<3, 3> _computeElementMatrixTRIA3(Cell cell);
  FixedMatrix<4, 4> _computeElementMatrixQUAD4(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
  Real _computeAreaQuad4(Cell cell);
  Real _computeEdgeLength2(Face face);
  Real2 _computeEdgeNormal2(Face face);

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
  // Test for adding parameters for PETSc.
  // This is only used for the first call.
  {
    StringList string_list;
    string_list.add("-trmalloc");
    string_list.add("-log_trace");
    CommandLineArguments args(string_list);
    m_linear_system.setSolverCommandLineArguments(args);
  }
  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();
  _doStationarySolve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "Module Fem INIT";

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();

  //_buildDoFOnNodes();
  //Int32 nb_node = allNodes().size();
  //m_k_matrix.resize(nb_node, nb_node);
  //m_k_matrix.fill(0.0);

  //m_rhs_vector.resize(nb_node);
  //m_rhs_vector.fill(0.0);

  // # init mesh
  // # init behavior
  // # init behavior on mesh entities
  // # init BCs
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

  // # T=linalg.solve(K,RHS)
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
  ElementNodes = 3.;

  if (options()->meshType == "QUAD4")
    ElementNodes = 4.;
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
        m_u[node] = value;
        m_u_dirichlet[node] = true;
      }
    }
  }

  for (const auto& bs : options()->dirichletPointCondition()) {
    NodeGroup group = bs->node();
    Real value = bs->value();
    info() << "Apply Dirichlet point condition node=" << group.name() << " v=" << value;
    ENUMERATE_ (Node, inode, group) {
      Node node = *inode;
      m_u[node] = value;
      m_u_dirichlet[node] = true;
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
  info() << "Applying Dirichlet boundary condition via  penalty method ";

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
      if (m_u_dirichlet[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_linear_system.matrixSetValue(dof_id, dof_id, Penalty);
        Real u_g = Penalty * m_u[node_id];
        rhs_values[dof_id] = u_g;
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
      if (m_u_dirichlet[node_id]) {
        DoFLocalId dof_id = node_dof.dofId(*inode, 0);
        m_linear_system.matrixAddValue(dof_id, dof_id, Penalty);
        Real u_g = Penalty * m_u[node_id];
        rhs_values[dof_id] = u_g;
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
  // Constant flux term assembly
  //----------------------------------------------
  //
  //  only for noded that are non-Dirichlet
  //  $int_{dOmega_N}((q.n)*v^h)$
  // or
  //  $int_{dOmega_N}((n_x*q_x + n_y*q_y)*v^h)$
  //----------------------------------------------
  for (const auto& bs : options()->neumannBoundaryCondition()) {
    FaceGroup group = bs->surface();

    if(bs->value.isPresent()) {
      Real value = bs->value();
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real length = _computeEdgeLength2(face);
        for (Node node : iface->nodes()) {
          if (!(m_u_dirichlet[node]) && node.isOwn())
            rhs_values[node_dof.dofId(node, 0)] += value * length / 2.;
        }
      }
      continue;
    }


    if(bs->valueX.isPresent()  && bs->valueY.isPresent()) {
      Real valueX = bs->valueX();
      Real valueY = bs->valueY();
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real  length = _computeEdgeLength2(face);
        Real2 Normal = _computeEdgeNormal2(face);
        for (Node node : iface->nodes()) {
          if (!(m_u_dirichlet[node]) && node.isOwn())
            rhs_values[node_dof.dofId(node, 0)] += (Normal.x*valueX + Normal.y*valueY) * length / 2.;
        }
      }
      continue;
    }

    if(bs->valueX.isPresent()) {
      Real valueX = bs->valueX();
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real  length = _computeEdgeLength2(face);
        Real2 Normal = _computeEdgeNormal2(face);
        for (Node node : iface->nodes()) {
          if (!(m_u_dirichlet[node]) && node.isOwn())
            rhs_values[node_dof.dofId(node, 0)] += (Normal.x*valueX) * length / 2.;
        }
      }
      continue;
    }

    if(bs->valueY.isPresent()) {
      Real valueY = bs->valueY();
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real  length = _computeEdgeLength2(face);
        Real2 Normal = _computeEdgeNormal2(face);
        for (Node node : iface->nodes()) {
          if (!(m_u_dirichlet[node]) && node.isOwn())
            rhs_values[node_dof.dofId(node, 0)] += (Normal.y*valueY) * length / 2.;
        }
      }
      continue;
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

FixedMatrix<3, 3> FemModule::
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

  FixedMatrix<2, 3> b_matrix;
  b_matrix(0, 0) = dPhi0.x;
  b_matrix(0, 1) = dPhi1.x;
  b_matrix(0, 2) = dPhi2.x;

  b_matrix(1, 0) = dPhi0.y;
  b_matrix(1, 1) = dPhi1.y;
  b_matrix(1, 2) = dPhi2.y;

  b_matrix.multInPlace(1.0 / (2.0 * area));

  FixedMatrix<3, 3> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(area);

  //info() << "Cell=" << cell.localId();
  //std::cout << " int_cdPi_dPj=";
  //int_cdPi_dPj.dump(std::cout);
  //std::cout << "\n";

  return int_cdPi_dPj;
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
    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
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

    auto K_e = _computeElementMatrixTRIA3(cell);  // element stifness matrix
    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
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
_solve()
{
  m_linear_system.solve();

  // Re-Apply boundary conditions because the solver has modified the value
  // of u on all nodes
  _applyDirichletBoundaryConditions();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    // Copy RHS DoF to Node u
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_u[node_dof.dofId(node, 0)];
      m_u[node] = v;
    }
  }

  m_u.synchronize();

  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "u[" << node.localId() << "][" << node.uniqueId() << "] = "
             << m_u[node];
      //info() << "u[]" << node.uniqueId() << " "
      //       << m_u[node];
    }
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
  checkNodeResultFile(traceMng(), filename, m_u, epsilon);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
