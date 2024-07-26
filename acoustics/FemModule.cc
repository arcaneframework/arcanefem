// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2024 */
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

  Real kc2;
  Real ElementNodes;

  DoFLinearSystem m_linear_system;
  IItemFamily* m_dof_family = nullptr;
  FemDoFsOnNodes m_dofs_on_nodes;

 private:

  void _doStationarySolve();
  void _getMaterialParameters();
  void _assembleBilinearOperatorTRIA3();
  void _solve();
  void _assembleLinearOperator();
  void _checkResultFile();
  FixedMatrix<3, 3> _computeElementMatrixTRIA3(Cell cell);
  Real _computeAreaTriangle3(Cell cell);
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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  // get material parameters
  _getMaterialParameters();

  // Assemble the FEM bilinear operator (LHS - matrix A)
  _assembleBilinearOperatorTRIA3();

  // Assemble the FEM linear operator (RHS - vector b)
  _assembleLinearOperator();

  // Solve the FEM linear system u= A^1*b
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
  kc2 = options()->kc2();

  ElementNodes = 3.;
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

  //----------------------------------------------
  // Constant flux term assembly
  //----------------------------------------------
  //
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
          if (node.isOwn())
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
          if (node.isOwn())
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
          if (node.isOwn())
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
          if (node.isOwn())
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
// Compute the element matrix for a triangular element (P1 finite elements)
// - This function calculates the integral of ( -(u.dx*v.dx + u.dy*v.dy) + sigma*u*v)
//   over the triangular element using linear shape functions (P1)
/*---------------------------------------------------------------------------*/

FixedMatrix<3, 3> FemModule::
_computeElementMatrixTRIA3(Cell cell)
{
  // Get coordinates of the triangle vertices
  Real3 vertex0 = m_node_coord[cell.nodeId(0)];
  Real3 vertex1 = m_node_coord[cell.nodeId(1)];
  Real3 vertex2 = m_node_coord[cell.nodeId(2)];

  // Calculate the area of the triangle
  Real area = _computeAreaTriangle3(cell);

  // Compute the gradients of the shape functions
  Real2 gradN0(vertex1.y - vertex2.y, vertex2.x - vertex1.x);
  Real2 gradN1(vertex2.y - vertex0.y, vertex0.x - vertex2.x);
  Real2 gradN2(vertex0.y - vertex1.y, vertex1.x - vertex0.x);

  // Normalize the gradients by the area
  Real A2 = 2.0 * area;
  gradN0 /= A2;
  gradN1 /= A2;
  gradN2 /= A2;


  // Compute the element matrix
  FixedMatrix<3, 3> elementMatrix;

  // -(u.dx*v.dx + u.dy*v.dy) terms
  // For P1 elements, dx and dy are symmetric in the integral evaluation
  // So these terms are the same as dx(u) * dx(v) terms
  elementMatrix(0, 0) = -area * Arcane::math::dot(gradN0, gradN0);
  elementMatrix(0, 1) = -area * Arcane::math::dot(gradN0, gradN1);
  elementMatrix(0, 2) = -area * Arcane::math::dot(gradN0, gradN2);

  elementMatrix(1, 0) = -area * Arcane::math::dot(gradN1, gradN0);
  elementMatrix(1, 1) = -area * Arcane::math::dot(gradN1, gradN1);
  elementMatrix(1, 2) = -area * Arcane::math::dot(gradN1, gradN2);

  elementMatrix(2, 0) = -area * Arcane::math::dot(gradN2, gradN0);
  elementMatrix(2, 1) = -area * Arcane::math::dot(gradN2, gradN1);
  elementMatrix(2, 2) = -area * Arcane::math::dot(gradN2, gradN2);

  // Subtract uv term
  // For P1 elements, the integral of uv over the element is:
  // Diagonal terms: area / 6
  // Off-diagonal terms: area / 12
  Real uvTerm = kc2*area / 12.0;
  elementMatrix(0, 0) += uvTerm*2.;
  elementMatrix(0, 1) += uvTerm;
  elementMatrix(0, 2) += uvTerm;

  elementMatrix(1, 0) += uvTerm;
  elementMatrix(1, 1) += uvTerm*2.;
  elementMatrix(1, 2) += uvTerm;

  elementMatrix(2, 0) += uvTerm;
  elementMatrix(2, 1) += uvTerm;
  elementMatrix(2, 2) += uvTerm*2.;

  return elementMatrix;
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

    // element stiffness matrix
    auto K_e = _computeElementMatrixTRIA3(cell);
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
