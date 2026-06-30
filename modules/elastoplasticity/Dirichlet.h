// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Dirichlet.h                                                 (C) 2000-2026 */
/*                                                                           */
/* Contains functions to compute and assemble dirichlet contribution to RHS  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies dirichlet to LHS matrix and RHS vector of the linear system.
 * 
 * This function applies Dirichlet boundary conditions to both the LHS matrix
 * and RHS vector of the linear system.
 * 
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 */
/*---------------------------------------------------------------------------*/

inline void FemModuleElastoplasticity::
_applyDirichlet(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  // check if Hypre|Petsc solver is used and delegate to GPU for dirichlet assembly
  auto use_gpu = options()->linearSystem.serviceName() == "HypreLinearSystem" ||
    options()->linearSystem.serviceName() == "PetscLinearSystem";
  if (use_gpu) {
    _assembleDirichletsGpu();
    return;
  }

  info() << "[ArcaneFem-Info] Started module _assembleDirichletsCpu()";

  BC::IArcaneFemBC* bc = options()->boundaryConditions();
  if (bc) {
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
      ArcaneFemFunctions::BoundaryConditions::applyDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);

    for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions())
      ArcaneFemFunctions::BoundaryConditions::applyPointDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies dirichlet to LHS matrix and RHS vector of the linear system on Gpu.
 * 
 * This function applies Dirichlet boundary conditions to both the LHS matrix
 * and RHS vector of the linear system using GPU acceleration.
 * 
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::_assembleDirichletsGpu()
{
  info() << "[ArcaneFem-Info] Started module  _assembleDirichletsGpu()";

  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  auto mesh_ptr = mesh();

  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  if (bc) {
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions())
      FemUtils::Gpu::BoundaryConditions::applyDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);

    for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions())
      FemUtils::Gpu::BoundaryConditions::applyPointDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);
  }
}


/*---------------------------------------------------------------------------*/
/**
 * @brief Applies dirichlet to LHS matrix and RHS vector of the nonlinear system.
 *
 * This function applies Dirichlet boundary conditions to both the LHS matrix
 * and RHS vector of the nonlinear system.
 *
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::
_applyDirichletNewton(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  // check if Hypre|Petsc solver is used and delegate to GPU for dirichlet assembly
  auto use_gpu = options()->linearSystem.serviceName() == "HypreLinearSystem" ||
    options()->linearSystem.serviceName() == "PetscLinearSystem";
  if (use_gpu) {
    _assembleDirichletsNewtonGpu();
    return;
  }

  info() << "[ArcaneFem-Info] Started module _assembleDirichletsNewtonCpu()";

  BC::IArcaneFemBC* bc = options()->boundaryConditions();
  if (bc) {
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions()) {
      FaceGroup face_group = bs->getSurface();
      NodeGroup node_group = face_group.nodeGroup();
      const StringConstArrayView u_dirichlet_string = bs->getValue();
      for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
        if (u_dirichlet_string[dof_index] != "NULL") {

          Real value = 0.0;
          if (m_newton_iter == 0) {
            value = std::stod(u_dirichlet_string[dof_index].localstr());
          }

          if (bs->getEnforceDirichletMethod() == "Penalty") {
            Real penalty = bs->getPenalty();
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, node_dof, m_linear_system, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowColumnElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
          }
          else {
            ARCANE_FATAL("Unknown Dirichlet method");
          }
        }
      }
    }
      // ArcaneFemFunctions::BoundaryConditions::applyDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);

    for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions()) {
      NodeGroup node_group = bs->getNode();
      const StringConstArrayView u_dirichlet_string = bs->getValue();
      for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
        if (u_dirichlet_string[dof_index] != "NULL") {
          Real value = 0.0;
          if (m_newton_iter == 0) {
            value = std::stod(u_dirichlet_string[dof_index].localstr());
          }
          if (bs->getEnforceDirichletMethod() == "Penalty") {
            Real penalty = bs->getPenalty();
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, node_dof, m_linear_system, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowColumnElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
          }
          else {
            ARCANE_FATAL("Unknown Dirichlet method");
          }
        }
      }
    }
      // ArcaneFemFunctions::BoundaryConditions::applyPointDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies dirichlet to LHS matrix and RHS vector of the nonlinear system on Gpu.
 *
 * This function applies Dirichlet boundary conditions to both the LHS matrix
 * and RHS vector of the nonlinear system using GPU acceleration.
 *
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 */
/*---------------------------------------------------------------------------*/

void FemModuleElastoplasticity::_assembleDirichletsNewtonGpu()
{
  info() << "[ArcaneFem-Info] Started module  _assembleDirichletsNewtonGpu()";

  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  auto mesh_ptr = mesh();

  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  if (bc) {
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions()) {
      ARCANE_CHECK_PTR(bs);

      FaceGroup face_group = bs->getSurface();
      NodeGroup node_group = face_group.nodeGroup();

      const StringConstArrayView u_dirichlet_string = bs->getValue();

      for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
        if (u_dirichlet_string[dof_index] != "NULL") {
          Real value = 0.0;
          if (m_newton_iter == 0) {
            value = std::stod(u_dirichlet_string[dof_index].localstr());
          }
          if (bs->getEnforceDirichletMethod() == "Penalty") {
            Real penalty = bs->getPenalty();
            Gpu::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, queue, mesh_ptr, m_linear_system, m_dofs_on_nodes, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowElimination") {
            Gpu::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowOrRowColumnElimination(ELIMINATE_ROW, dof_index, value, queue, m_linear_system, m_dofs_on_nodes, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
            Gpu::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowOrRowColumnElimination(ELIMINATE_ROW_COLUMN, dof_index, value, queue, m_linear_system, m_dofs_on_nodes, node_group);
          }
          else {
            ARCANE_FATAL("Unknown method to enforce Dirichlet BC: '{0}'", bs->getEnforceDirichletMethod());
          }
        }
      }
    }
      // FemUtils::Gpu::BoundaryConditions::applyDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);

    for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions()) {
      ARCANE_CHECK_PTR(bs);
      NodeGroup node_group = bs->getNode();

      const StringConstArrayView u_dirichlet_str = bs->getValue();

      for (Int32 dof_index = 0; dof_index < u_dirichlet_str.size(); ++dof_index) {
        if (u_dirichlet_str[dof_index] != "NULL") {
          Real value = 0.0;
          if (m_newton_iter == 0) {
            value = std::stod(u_dirichlet_str[dof_index].localstr());
          }

          if (bs->getEnforceDirichletMethod() == "Penalty") {
            Real penalty = bs->getPenalty();
            Gpu::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, queue, mesh_ptr, m_linear_system, m_dofs_on_nodes, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowElimination") {
            Gpu::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowOrRowColumnElimination(ELIMINATE_ROW, dof_index, value, queue, m_linear_system, m_dofs_on_nodes, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
            Gpu::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowOrRowColumnElimination(ELIMINATE_ROW_COLUMN, dof_index, value, queue, m_linear_system, m_dofs_on_nodes, node_group);
          }
          else {
            ARCANE_FATAL("Unknown method to enforce Dirichlet BC: '{0}'", bs->getEnforceDirichletMethod());
          }
        }
      }
    }
      // FemUtils::Gpu::BoundaryConditions::applyPointDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);
  }
}