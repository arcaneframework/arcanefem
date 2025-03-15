// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* Simple module to solve aerodynamics using FEM.                            */
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
  if (m_global_iteration() > 0)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();
  _doStationarySolve();
  _getPsi();
 
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "Module Fem INIT";

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_getPsi()
{
  info() << "Postprocessing PSI";

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Real2 DX = _computeDxDyOfRealTria3(cell);
    m_psi[cell] = -DX.x * DX.x - DX.y * DX.y;
  }

  m_psi.synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  _assembleBilinearOperatorTria3();
  _assembleLinearOperator();
  _solve();
  _checkResultFile();
}

/*---------------------------------------------------------------------------*/
// Assemble the FEM linear operator
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "Assembly of FEM linear operator ";

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  BC::IArcaneFemBC* bc = options()->boundaryConditions();
  if (bc)
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
      ArcaneFemFunctions::BoundaryConditions2D::applyDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);

  for (const auto& bs : options()->farfieldBoundaryCondition()) {
    FaceGroup group = bs->surface();
    Real value = bs->angle();

    Real Penalty = options()->penalty();

    ENUMERATE_ (Face, iface, group) {
      for (Node node : iface->nodes()) {
        if (node.isOwn()) {
          m_linear_system.matrixSetValue(node_dof.dofId(node, 0), node_dof.dofId(node, 0), Penalty);
          Real u_g = (m_node_coord[node].y - value * m_node_coord[node].x) * Penalty;
          rhs_values[node_dof.dofId(node, 0)] = u_g;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real2 FemModule::
_computeDxDyOfRealTria3(Cell cell)
{
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real f0 = m_u[cell.nodeId(0)];
  Real f1 = m_u[cell.nodeId(1)];
  Real f2 = m_u[cell.nodeId(2)];

  Real detA = ( m0.x*(m1.y - m2.y) - m0.y*(m1.x - m2.x) + (m1.x*m2.y - m2.x*m1.y) );

  Real2 DX;
        DX.x = ( m0.x*(f1 - f2) - f0*(m1.x - m2.x) + (f2*m1.x - f1*m2.x) ) / detA;
        DX.y = ( f0*(m1.y - m2.y) - m0.y*(f1 - f2) + (f1*m2.y - f2*m1.y) ) / detA;

  return DX ;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTria3()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    auto K_e = _computeElementMatrixTria3(cell);  // element stifness matrix
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  m_linear_system.solve();

  {
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_u[node_dof.dofId(node, 0)];
      m_u[node] = v;
    }
  }

  m_u.synchronize();

  if (allNodes().size() < 200)
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "u[" << node.localId() << "][" << node.uniqueId() << "] = " << m_u[node];
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
