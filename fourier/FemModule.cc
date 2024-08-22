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

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModule at the start of the simulation.
 *
 *  - initializes degrees of freedom (DoFs) on nodes.
 *  - builds support for manufactured test case (optional).
 */
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "Module Fem INIT";

  m_dofs_on_nodes.initialize(mesh(), 1);
  m_dof_family = m_dofs_on_nodes.dofFamily();

  if (options()->manufacturedSolution.isPresent()) {
    const auto& bs = options()->manufacturedSolution()[0];

    if (bs->manufacturedSource.isPresent()) {
      ICaseFunction* opt_function_source = bs->manufacturedSource.function();
      IStandardFunction* scf_source = bs->manufacturedSource.standardFunction();
      if (!scf_source)
        ARCANE_FATAL("No standard case function for option 'manufactured-source-condition'");
      auto* functorS = scf_source->getFunctorRealReal3ToReal();
      if (!functorS)
        ARCANE_FATAL("Standard function '{0}' is not convertible to f(Real,Real3) -> Real", opt_function_source->name());
      m_manufactured_source = functorS;
    }

    if (bs->manufacturedDirichlet.isPresent()) {
      ICaseFunction* opt_function = bs->manufacturedDirichlet.function();
      IStandardFunction* scf = bs->manufacturedDirichlet.standardFunction();
      if (!scf)
        ARCANE_FATAL("No standard case function for option 'manufactured-dirichlet-condition'");
      auto* functor = scf->getFunctorRealReal3ToReal();
      if (!functor)
        ARCANE_FATAL("Standard function '{0}' is not convertible to f(Real,Real3) -> Real", opt_function->name());
      m_manufactured_dirichlet = functor;
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs the main computation for the FemModule.
 *
 * - Stops the time loop after 1 iteration since the equation is steady state.
 * - Resets, configures, and initializes the linear system.
 * - Executes the stationary solve.
 */
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
    string_list.add("-ksp_monitor");
    CommandLineArguments args(string_list);
    m_linear_system.setSolverCommandLineArguments(args);
  }
  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();
  _doStationarySolve();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows a sequence of steps to solve FEM system:
 *   1. _getMaterialParameters()          Retrieves material parameters via
 *   2. _assembleBilinearOperator()       Assembles the FEM  matrix A
 *   3. _assembleLinearOperator()         Assembles the FEM RHS vector b
 *   4.  _solve()                         Solves for solution vector u = A^-1*b
 *   5. _validateResults()                Regression test
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  _getMaterialParameters();
  _assembleBilinearOperator();
  _assembleLinearOperator();
  _solve();
  _validateResults();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Retrieves and sets the material parameters for the simulation.
 *
 * This method initializes:
 *  - material properties:
 *       # thermal conductivity coefficient (`lambda`)
 *       # heat source term (`qdot`) 
 *  - mesh properties: 
 *       # number of cell nodes per element based on the mesh (quad or tria)
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_getMaterialParameters()
{
  info() << "Get material parameters ";

  lambda = options()->lambda();
  qdot = options()->qdot();

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    m_cell_lambda[cell] = lambda;
  }

  for (const auto& bs : options()->materialProperty()) {
    CellGroup group = bs->volume();
    Real value = bs->lambda();
    info() << "Lambda for group= " << group.name() << " v=" << value;

    ENUMERATE_ (Cell, icell, group) {
      Cell cell = *icell;
      m_cell_lambda[cell] = value;
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief FEM linear operator for the current simulation step.      
 *
 * This method constructs the linear  system by  assembling the LHS matrix 
 * and  RHS vector, applying various boundary conditions and source terms.
 *
 * Steps involved:
 *  1. The RHS vector is initialized to zero before applying any conditions.
 *  2. If a constant source term is specified (`qdot`), apply it to the RHS.
 *  3. If Neumann BC are specified applied to the RHS.
 *  4. If Dirichlet BC are specified apply to the LHS & RHS. 
 *  5. If manufactured source conditions is specified, apply to RHS.
 *  6. If  manufactured Dirichlet BC are specified apply to the LHS & RHS. 
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "Assembly of FEM linear operator";

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  if (options()->qdot.isPresent())
    ArcaneFemFunctions::BoundaryConditions2D::applyConstantSourceToRhs(qdot, mesh(), node_dof, m_node_coord, rhs_values);

  for (const auto& bs : options()->neumannBoundaryCondition())
    ArcaneFemFunctions::BoundaryConditions2D::applyNeumannToRhs(bs, node_dof, m_node_coord, rhs_values);

  for (const auto& bs : options()->dirichletBoundaryCondition())
    ArcaneFemFunctions::BoundaryConditions2D::applyDirichletToLhsAndRhs(bs, node_dof, m_node_coord, m_linear_system, rhs_values);

  if (options()->manufacturedSolution.isPresent()) {
    const auto& bs = options()->manufacturedSolution()[0];

    if (bs->manufacturedSource.isPresent()) {
      ARCANE_CHECK_POINTER(m_manufactured_source);
      info() << "Apply manufactured Source condition to all cells";
      ArcaneFemFunctions::BoundaryConditions2D::applyManufacturedSourceToRhs(m_manufactured_source, mesh(), node_dof, m_node_coord, rhs_values);
    }

    if (bs->manufacturedDirichlet.isPresent()) {
      ARCANE_CHECK_POINTER(m_manufactured_dirichlet);
      info() << "Apply manufactured dirichlet condition to all borders";
      FaceGroup group = mesh()->outerFaces();
      ArcaneFemFunctions::BoundaryConditions2D::applyManufacturedDirichletToLhsAndRhs(m_manufactured_dirichlet, lambda, group, bs, node_dof, m_node_coord, m_linear_system, rhs_values);
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the element matrix for a triangular element (P1 FE).
 *
 * This function calculates the integral of the expression:
 * lambda*(u.dx * v.dx + u.dy * v.dy)
 *
 * Steps involved:
 * 1. Calculate the area of the triangle.
 * 2. Compute the gradients of the shape functions.
 * 3. Return lambda*(u.dx * v.dx + u.dy * v.dy);
 */
/*---------------------------------------------------------------------------*/

FixedMatrix<3, 3> FemModule::
_computeElementMatrixTria3(Cell cell)
{
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  Real3 dxU = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyU = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

  return area * lambda * (dxU ^ dxU) + area * lambda * (dyU ^ dyU);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FixedMatrix<4, 4> FemModule::
_computeElementMatrixQuad4(Cell cell)
{
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaQuad4(cell, m_node_coord);

  Real4 dxU = ArcaneFemFunctions::FeOperation2D::computeGradientXQuad4(cell, m_node_coord);
  Real4 dyU = ArcaneFemFunctions::FeOperation2D::computeGradientYQuad4(cell, m_node_coord);

  return area * lambda * (dxU ^ dxU) + area * lambda * (dyU ^ dyU);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Calls the right function for LHS assembly given as mesh type.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperator()
{
  if (options()->meshType == "QUAD4")
    _assembleBilinear<4>([this](const Cell& cell) {
      return _computeElementMatrixQuad4(cell);
    });
  else if (options()->meshType == "TRIA3")
    _assembleBilinear<3>([this](const Cell& cell) {
      return _computeElementMatrixTria3(cell);
    });
  else
    ARCANE_FATAL("Non supported meshType");
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system.
 *
 * The method performs the following steps:
 *   1. For each cell, retrieves the cell-specific constant `lambda`.
 *   2. Computes element matrix using provided `compute_element_matrix` function.
 *   3. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

template <int N>
void FemModule::
_assembleBilinear(const std::function<FixedMatrix<N, N>(const Cell&)>& compute_element_matrix)
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    lambda = m_cell_lambda[cell]; // lambda is always considered cell constant
    auto K_e = compute_element_matrix(cell); // element matrix based on the provided function
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
/**
 * @brief Solves the linear system and updates the solution vector.
 *
 * This method performs the following actions:
 *   1. Solves the linear system to compute the solution.
 *   2. Copies the computed solution from the DoF to the node values.
 *   3. Synchronizes the updated node values.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  m_linear_system.solve();

  { // copies solution (and optionally exact solution) to FEM output
    VariableDoFReal& dof_u(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      Real v = dof_u[node_dof.dofId(node, 0)];
      m_u[node] = v;
    }

    if (options()->manufacturedSolution.isPresent()) {
      ENUMERATE_ (Node, inode, ownNodes()) {
        Node node = *inode;
        m_u_exact[node] = m_manufactured_dirichlet->apply(lambda, m_node_coord[node]);
      }
    }
  }

  m_u.synchronize();
  if (options()->manufacturedSolution.isPresent())
    m_u_exact.synchronize();
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Validates and prints the results of the FEM computation.
 *
 * This method performs the following actions:
 *   1. If number of nodes < 200, prints the computed values for each node.
 *   2. Retrieves the filename for the result file from options.
 *   3. If a filename is provided, checks the computed results against result file.
 *
 * @note The result comparison uses a tolerance of 1.0e-4.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_validateResults()
{
  if (allNodes().size() < 200)
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "u[" << node.localId() << "][" << node.uniqueId() << "] = " << m_u[node];
    }

  String filename = options()->resultFile();
  info() << "ValidateResultFile filename=" << filename;

  if (!filename.empty())
    checkNodeResultFile(traceMng(), filename, m_u, 1.0e-4);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
