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

#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>

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

  void _doStationarySolve();
  void _updateBoundayConditions();
  void _computeConductivity();
  void _computeGeneralizedFluxes();
  void _solve();
  void _initBoundaryconditions();
  void _applyOneBoundaryCondition(const String& group_name, Real value);
  FixedMatrix<3, 3> _computeIntCDPhiiDPhij(Cell cell);
  FixedMatrix<2, 3> _computeBMatrix(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
compute()
{
  info() << "Module Fem1 COMPUTE";

  // Stop code after computatio
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
  // # update BCs
  _updateBoundayConditions();

  // K=self.compute_conductivity()
  _computeConductivity();
  //      # compute RHS
  _computeGeneralizedFluxes();

  // # adapt K and RHS to take into account Dirichlet BCs
  //         for node in self.mesh.nodes:
  //             if node.is_T_fixed:
  //                 K[node.rank,node.rank]=K[node.rank,node.rank]+10**6
  //                 RHS[node.rank]=RHS[node.rank]+(10**6)*node.T

  info() << "TODO: adapt K and RHS";

  // # T=linalg.solve(K,RHS)
  _solve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_initBoundaryconditions()
{
  info() << "Init boundary conditions...";

  info() << "TODO: Init boundary conditions";

  // For this test apply fixed boundary conditions.
  // TODO: Adapt 'axl' file to allow the user to specifiy Boundary conditions

  // BC T=50.0 GROUP=Cercle
  // BC T=5.0 GROUP=Bas
  // BC T=21.0 GROUP=Haut
  _applyOneBoundaryCondition("Cercle", 50.0);
  _applyOneBoundaryCondition("Bas", 5.0);
  _applyOneBoundaryCondition("Haut", 21.0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_applyOneBoundaryCondition(const String& group_name, Real value)
{
  FaceGroup group = mesh()->faceFamily()->findGroup(group_name);
  if (group.null())
    ARCANE_FATAL("Can not find FaceGroup '{0}'", group_name);
  ENUMERATE_ (Face, iface, group) {
    for (Node node : iface->nodes()) {
      m_node_temperature[node] = value;
      m_node_is_temperature_fixed[node] = true;
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

  info() << "TODO: _computeBMatrix() m0=" << m0 << " m1=" << m1 << " m2=" << m2;

  // Copy 3D coordinates in 2D
  //Real2 m0(n0.x, n0.y);
  //Real2 m1(n1.x, n1.y);
  //Real2 m2(n2.x, n2.y);

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
  info() << "TODO: _compute int_c_dPhii_dPhij()";
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
  //         #print(int_cdPi_dPj)
  //        return int_cdPi_dPj
  FixedMatrix<3, 3> int_cdPi_dPj;
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
    _computeIntCDPhiiDPhij(cell);
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
        //C++                    K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        ++n2_index;
      }
      ++n1_index;
    }
  }

  info() << "TODO " << A_FUNCINFO;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_computeGeneralizedFluxes()
{
  info() << "TODO " << A_FUNCINFO;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_solve()
{
  info() << "TODO " << A_FUNCINFO;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM1(Fem1Module);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
