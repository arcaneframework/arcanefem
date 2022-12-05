// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Fem1Module.cc                                               (C) 2022-2022 */
/*                                                                           */
/* Simple module to test simple FEM mechanism.                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/NumArray.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>

#include "Fem1_axl.h"
#include "./FemUtils.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*!
 * \brief Module Fem1.
 */
class Fem1Module
: public ArcaneFem1Object
{
 public:

  explicit Fem1Module(const ModuleBuildInfo& mbi)
  : ArcaneFem1Object(mbi)
  {}

 public:

  /*!
   * \brief Méthode appelée à chaque itération.
   */
  void compute() override;
  /*!
   * \brief Méthode appelée lors de l'initialisation.
   */
  void startInit() override;

  /** Retourne le numéro de version du module */
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 private:

  //! K matrix
  NumArray<Real, MDDim2> m_k_matrix;
  //! RHS (Right Hand Side) vector
  NumArray<Real, MDDim1> m_rhs_vector;

  Real lambda;

 private:

  void _doStationarySolve();
  void _getMaterialParameters();
  void _updateBoundayConditions();
  void _computeConductivity();
  void _computeGeneralizedFluxes();
  void _computeSourceTerm();
  void _solve();
  void _initBoundaryconditions();
  void _applyPenaltyDirichletBC();
  FixedMatrix<3, 3> _computeIntCDPhiiDPhij(Cell cell);
  FixedMatrix<2, 3> _computeBMatrix(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
  void _applyDirichletBoundaryConditions();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
compute()
{
  info() << "Module Fem1 COMPUTE";

  // Stop code after computations
  subDomain()->timeLoopMng()->stopComputeLoop(true);

  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();
  _doStationarySolve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
startInit()
{
  info() << "Module Fem1 INIT";

  Int32 nb_node = allNodes().size();
  m_k_matrix.resize(nb_node, nb_node);
  m_k_matrix.fill(0.0);

  m_rhs_vector.resize(nb_node);
  m_rhs_vector.fill(0.0);

  // # init mesh
  // # init behavior
  // # init behavior on mesh entities
  // # init BCs
  _initBoundaryconditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_doStationarySolve()
{

  // # get material parameters
  _getMaterialParameters();

  // # update BCs
  _updateBoundayConditions();

  // K=self.compute_conductivity()
  _computeConductivity();

  //      # compute flux component RHS
  _computeGeneralizedFluxes();

  //      # compute source component RHS
  _computeSourceTerm();

  // Plenalty method to enforce Dirichlet conditions
  _applyPenaltyDirichletBC();

  // # T=linalg.solve(K,RHS)
  _solve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_getMaterialParameters()
{
  info() << "Get material parameters...";
  lambda = options()->lambda();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_initBoundaryconditions()
{
  info() << "Init boundary conditions...";

  info() << "Apply boundary conditions";
  _applyDirichletBoundaryConditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_updateBoundayConditions()
{
  info() << "TODO " << A_FUNCINFO;
}

/*---------------------------------------------------------------------------*/
// This function enforces a Dirichlet boundary condition in a weak sense by
// penalty method
/*---------------------------------------------------------------------------*/

void Fem1Module::
_applyPenaltyDirichletBC()
{
  info() << "Applying Dirichlet boundary condition via  penalty method ";

  // # adapt K and RHS to take into account Dirichlet BCs
  //         for node in self.mesh.nodes:
  //             if node.is_T_fixed:
  //                 K[node.rank,node.rank]=K[node.rank,node.rank]+10**6
  //                 RHS[node.rank]=RHS[node.rank]+(10**6)*node.T

  // TODO: 1.0e6 is a user value, moreover we should use seomthing like 1e31
  ENUMERATE_ (Node, inode, allNodes()) {
    NodeLocalId node_id = *inode;
    if (m_node_is_temperature_fixed[node_id]) {
      m_k_matrix(node_id, node_id) += 1.0e6;
      m_rhs_vector[node_id] += 1.0e6 * m_node_temperature[node_id];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real Fem1Module::
_computeAreaTriangle3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];
  return 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

//     """Compute matrix of gradient of FE shape functions for current element
//     B=[grad(Phi_0) grad(Phi1) grad(Phi2)] and return a numpy array
//     """
FixedMatrix<2, 3> Fem1Module::
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

  FixedMatrix<2, 3> b_matrix;
  b_matrix(0, 0) = dPhi0.x;
  b_matrix(0, 1) = dPhi1.x;
  b_matrix(0, 2) = dPhi2.x;

  b_matrix(1, 0) = dPhi0.y;
  b_matrix(1, 1) = dPhi1.y;
  b_matrix(1, 2) = dPhi2.y;

  //     B=1/(2*area)*array([[dPhi0[0],dPhi1[0],dPhi2[0]],
  //                             [dPhi0[1],dPhi1[1],dPhi2[1]]])
  b_matrix.multInPlace(1.0 / (2.0 * area));
  //     return(B)

  std::cout << "B=";
  b_matrix.dump(std::cout);
  std::cout << "\n";

  return b_matrix;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<3, 3> Fem1Module::
_computeIntCDPhiiDPhij(Cell cell)
{
  //const Real c = 1.75;
  FixedMatrix<2, 3> b_matrix = _computeBMatrix(cell);
  //         B=self.compute_B_matrix()
  //         print("B=",B, "\nT=",B.T)
  //         area=self.compute_area()
  Real area = _computeAreaTriangle3(cell); //m0, m1, m2);
  //         #z = dot(B.T,B)
  //         #print("Z=",z)
  //         #z2 = matmul(B.T,B)
  //         #print("Z2=",z2)
  //         int_cdPi_dPj=area*c*dot(B.T,B)
  FixedMatrix<3, 3> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(area * lambda);
  //         #print(int_cdPi_dPj)
  //        return int_cdPi_dPj

  info() << "Cell=" << cell.localId();
  std::cout << " int_cdPi_dPj=";
  int_cdPi_dPj.dump(std::cout);
  std::cout << "\n";

  return int_cdPi_dPj;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_computeConductivity()
{
  // for elem in self.mesh.elements:

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    if (cell.type() != IT_Triangle3)
      ARCANE_FATAL("Only Triangle3 cell type is supported");
    //             # first compute elementary thermal conductivity matrix
    //             K_e=elem.int_c_dPhii_dPhij(elem.k)
    auto K_e = _computeIntCDPhiiDPhij(cell);
    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
    Int32 n1_index = 0;
    for (NodeLocalId node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (NodeLocalId node2 : cell.nodes()) {
        //                 K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        m_k_matrix(node1.localId(), node2.localId()) += K_e(n1_index, n2_index);
        ++n2_index;
      }
      ++n1_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_computeGeneralizedFluxes()
{
  // TODO: Loop over all faces on the border instead
  m_rhs_vector.fill(0.0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_computeSourceTerm()
{
  // TODO: Loop over all cells and fill the source term
  m_rhs_vector.fill(0.0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_solve()
{
  Int32 matrix_size = m_k_matrix.extent0();
  Arcane::MatVec::Matrix matrix(matrix_size, matrix_size);
  _convertNumArrayToCSRMatrix(matrix, m_k_matrix.span());

  info() << "TODO " << A_FUNCINFO;
  int p = std::cout.precision();
  std::cout.precision(12);
  matrix.dump(std::cout);
  std::cout.precision(p);

  Arcane::MatVec::Vector vector_b(matrix_size);
  Arcane::MatVec::Vector vector_x(matrix_size);
  {
    auto vector_b_view = vector_b.values();
    auto vector_x_view = vector_x.values();
    for (Int32 i = 0; i < matrix_size; ++i) {
      vector_b_view(i) = m_rhs_vector[i];
      vector_x_view(i) = 0.0;
    }
  }

  {
    Real epsilon = 1.0e-15;
    Arcane::MatVec::DiagonalPreconditioner p(matrix);
    Arcane::MatVec::ConjugateGradientSolver solver;
    solver.solve(matrix, vector_b, vector_x, epsilon, &p);
  }

  // def update_T(self,T):
  //     """Update temperature value on nodes after the FE resolution"""
  //     for i in range(0,len(self.mesh.nodes)):
  //         node=self.mesh.nodes[i]
  //         # don't update T imposed by Dirichlet BC
  //         if not node.is_T_fixed:
  //             self.mesh.nodes[i].T=T[i]

  {
    auto vector_x_view = vector_x.values();
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      if (!m_node_is_temperature_fixed[node]) {
        m_node_temperature[node] = vector_x_view[node.localId()];
      }
      info() << "T[" << node.localId() << "] = " << m_node_temperature[node];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM1(Fem1Module);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
