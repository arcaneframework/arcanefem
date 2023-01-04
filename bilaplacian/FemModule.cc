// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2022 */
/*                                                                           */
/* FEM code to test vectorial FE for bilaplacian problem.                    */
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

  Real f;

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;

 private:

  void _doStationarySolve();
  void _getMaterialParameters();
  void _updateBoundayConditions();
  void _assembleBilinearOperator();
  void _solve();
  void _initBoundaryconditions();
  void _assembleLinearOperator();
  FixedMatrix<6, 6> _computeIntCDPhiiDPhij(Cell cell);
  FixedMatrix<6, 6> _computeIntPhiiDPhij(Cell cell);
  FixedMatrix<2, 6> _computeBMatrix(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
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
  _assembleBilinearOperator();

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
  f   = options()->f();
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
        m_u1[node] = value;
        m_u1_fixed[node] = true;
      }
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
  //VariableNodeReal rhs1_values(VariableBuildInfo(defaultMesh(), "NodeRHS1Values"));
  //rhs1_values.fill(0.0);

  //VariableNodeReal rhs2_values(VariableBuildInfo(defaultMesh(), "NodeRHS2Values"));
  //rhs2_values.fill(0.0);

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  //----------------------------------------------
  // penelty method for assembly of Dirichlet BC
  //----------------------------------------------
  //
  // # adapt K and RHS to take into account Dirichlet BCs
  //         for node in self.mesh.nodes:
  //             if node.is_T_fixed:
  //                 K[node.rank,node.rank]=K[node.rank,node.rank]+10**6
  //                 RHS[node.rank]=RHS[node.rank]+(10**6)*node.T
  // TODO: 1.0e6 is a user value, moreover we should use something like 1e31
  //----------------------------------------------

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Node, inode, ownNodes()) {
    NodeLocalId node_id = *inode;
    if (m_u1_fixed[node_id]) {
      DoFLocalId dof_id1 = node_dof.dofId(node_id, 0);
      DoFLocalId dof_id2 = node_dof.dofId(node_id, 1);
      //m_k_matrix(node_id, node_id) += 1.0e6;
      //                                             u1   ,NA,NA, u2
      // m_linear_system.matrixAddValue(*inode, *inode, 1.0e30, 0, 0, 0);
      m_linear_system.matrixAddValue(dof_id1, dof_id1, 1.0e30);
      m_linear_system.matrixAddValue(dof_id2, dof_id2, 0.0);
      //m_rhs_vector[node_id] += 1.0e6 * m_node_temperature[node_id];
      {
        Real temperature = 1.0e30 * m_u1[node_id];
        rhs_values[dof_id1] = temperature;
      }
    }
  }

  info() << "Applying Dirichlet boundary condition via  penalty method 1 ";
  //----------------------------------------------
  // Constant source term assembly
  //----------------------------------------------
  //
  //  $int_{Omega}(f*v^h)$
  //  only for noded that are non-Dirichlet
  //----------------------------------------------
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real area = _computeAreaTriangle3(cell);
    for (Node node : cell.nodes()) {
      if (!(m_u1_fixed[node]) && node.isOwn()) {
        DoFLocalId dof_id1 = node_dof.dofId(node, 0);
        rhs_values[dof_id1] += f * area / 3;
      }
    }
  }
  info() << "Applying Dirichlet boundary condition via  penalty method 2 ";
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
  info() << "Applying Dirichlet boundary condition via  penalty method 3 ";
#if 0
  {
    // For the LinearSystem class we need an array
    // with only the values for the ownNodes().
    // The values of 'rhs_values' should not be updated after
    // this call.
    UniqueArray<Real> rhs_values_for_linear_system;
    ENUMERATE_ (Node, inode, ownNodes()) {
      rhs_values_for_linear_system.add(rhs1_values[inode]);
      rhs_values_for_linear_system.add(rhs2_values[inode]);
    }
    for (Int32 i = 0; i < rhs_values_for_linear_system.size(); ++i) {
      cout << "VECT[" << i << "] = " << rhs_values_for_linear_system[i] << endl;
    }
    m_linear_system.setRHSValues(rhs_values_for_linear_system);
  }
#endif
  info() << "Applying Dirichlet boundary condition via  penalty method 4 ";
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

//     """Compute matrix of gradient of FE shape functions for current element
//     B=[grad(Phi_0) grad(Phi1) grad(Phi2)] and return a numpy array
//     """
FixedMatrix<2, 6> FemModule::
_computeBMatrix(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  //     (M0,M1,M2)=(self.nodes[0],self.nodes[1],self.nodes[2])

  //     area=self.compute_area()
  Real area = _computeAreaTriangle3(cell); //m0, m1, m2);
  //     dPhi0=[M1.y-M2.y,M2.x-M1.x]
  //     dPhi1=[M2.y-M0.y,M0.x-M2.x]
  //     dPhi2=[M0.y-M1.y,M1.x-M0.x]
  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  FixedMatrix<2, 6> b_matrix;
  b_matrix(0, 0) = dPhi0.x +0;
  b_matrix(0, 1) = dPhi0.x +0;  
  b_matrix(0, 2) = dPhi1.x +0;
  b_matrix(0, 3) = dPhi1.x +0;
  b_matrix(0, 4) = dPhi2.x +0;
  b_matrix(0, 5) = dPhi2.x +0;

  b_matrix(1, 0) = dPhi0.y +0;
  b_matrix(1, 1) = dPhi0.y +0;
  b_matrix(1, 2) = dPhi1.y +0;
  b_matrix(1, 3) = dPhi1.y +0;
  b_matrix(1, 4) = dPhi2.y +0;
  b_matrix(1, 5) = dPhi2.y +0;

  //     B=1/(2*area)*array([[dPhi0[0],dPhi1[0],dPhi2[0]],
  //                             [dPhi0[1],dPhi1[1],dPhi2[1]]])
  b_matrix.multInPlace(1.0 / (2.0 * area));
  //     return(B)

  //std::cout << "B=";
  //b_matrix.dump(std::cout);
  //std::cout << "\n";

  return b_matrix;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<6, 6> FemModule::
_computeIntCDPhiiDPhij(Cell cell)
{
  //const Real c = 1.75;
  FixedMatrix<2, 6> b_matrix = _computeBMatrix(cell);
  //         B=self.compute_B_matrix()
  //         print("B=",B, "\nT=",B.T)
  //         area=self.compute_area()
  Real area = _computeAreaTriangle3(cell); //m0, m1, m2);
  //         #z = dot(B.T,B)
  //         #print("Z=",z)
  //         #z2 = matmul(B.T,B)
  //         #print("Z2=",z2)
  //         int_cdPi_dPj=area*c*dot(B.T,B)
  FixedMatrix<6, 6> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(area);
  //         #print(int_cdPi_dPj)
  //        return int_cdPi_dPj

  //info() << "Cell=" << cell.localId();
  //std::cout << " int_cdPi_dPj=";
  //int_cdPi_dPj.dump(std::cout);
  //std::cout << "\n";

  return int_cdPi_dPj;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<6, 6> FemModule::
_computeIntPhiiDPhij(Cell cell)
{

  Real area = _computeAreaTriangle3(cell); //m0, m1, m2);
  
  //const Real c = 1.75;
  FixedMatrix<1, 6> b_matrix ;// = _computeBMatrix(cell);
  
  b_matrix(0, 0) = 1/12.;//dPhi0.x;
  b_matrix(0, 1) = 1/12.;//dPhi1.x;
  b_matrix(0, 2) = 1/12.;//dPhi2.x;
  b_matrix(0, 3) = 1/12.;//dPhi0.x;
  b_matrix(0, 4) = 1/12.;//dPhi1.x;
  b_matrix(0, 5) = 1/12.;//dPhi2.x;
    
  FixedMatrix<6, 1> bT_matrix ;// = _computeBMatrix(cell);
  
  bT_matrix(0, 0) = 1.;//dPhi0.x;
  bT_matrix(1, 0) = 1.;//dPhi1.x;
  bT_matrix(2, 0) = 1.;//dPhi2.x;  
  bT_matrix(3, 0) = 1.;//dPhi0.x;
  bT_matrix(4, 0) = 1.;//dPhi1.x;
  bT_matrix(5, 0) = 1.;//dPhi2.x; 
  
  //bT_matrix.multInPlace(1.0 / (2.0 * area));
        
  //         B=self.compute_B_matrix()
  //         print("B=",B, "\nT=",B.T)
  //         area=self.compute_area()

  //         #z = dot(B.T,B)
  //         #print("Z=",z)
  //         #z2 = matmul(B.T,B)
  //         #print("Z2=",z2)
  //         int_cdPi_dPj=area*c*dot(B.T,B)
  FixedMatrix<6, 6> int_cdPi_dPj = matrixMultiplication(bT_matrix, b_matrix);
  
  for (Int32 i = 0; i<6; i++)
    int_cdPi_dPj(i,i) *= 2.; 
  int_cdPi_dPj.multInPlace(area);
  //         #print(int_cdPi_dPj)
  //        return int_cdPi_dPj

  //info() << "Cell=" << cell.localId();
  //std::cout << " int_cdPi_dPj=";
  //int_cdPi_dPj.dump(std::cout);
  //std::cout << "\n";

  return int_cdPi_dPj;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperator()
{
  // for elem in self.mesh.elements:

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    if (cell.type() != IT_Triangle3)
      ARCANE_FATAL("Only Triangle3 cell type is supported");
    //             # first compute elementary thermal conductivity matrix
    //             K_e=elem.int_c_dPhii_dPhij(elem.k)
    auto K_e = _computeIntCDPhiiDPhij(cell);
    auto K_e2 = _computeIntPhiiDPhij(cell);
    Real area = _computeAreaTriangle3(cell); //m0, m1, m2);
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
        //                 K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v1 = 0; //K_e(2*n1_index  ,   2*n2_index);
        Real v2 = K_e(2 * n1_index, 2 * n2_index + 1);
        Real v3 = K_e(2 * n1_index + 1, 2 * n2_index);
        Real v4 = K_e2(2 * n1_index + 1, 2 * n2_index + 1); //area/12.;//K_e2(2*n1_index+1,  2*n2_index+1 );
        //m_k_matrix(node1.localId(), node2.localId()) += v;
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
_solve()
{
  info() << "Solving Linear system";
  m_linear_system.solve();

  // Re-Apply boundary conditions because the solver has modified the value
  _applyDirichletBoundaryConditions();

  {
    VariableDoFReal& dof_temperature(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real u1_val = dof_temperature[node_dof.dofId(node, 0)];
      Real u2_val = dof_temperature[node_dof.dofId(node, 1)];
      m_u1[node] = u1_val;
      m_u2[node] = u2_val;
      info() << "Node: " << node.localId() << " U1=" << u1_val << " U2=" << u2_val;
    }
  }

  m_u1.synchronize();
  m_u2.synchronize();

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
