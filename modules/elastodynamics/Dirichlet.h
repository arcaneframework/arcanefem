﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Dirichlet.h                                                 (C) 2022-2025 */
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

inline void FemModule::
_applyDirichlet(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{

  // check if Hypre solver is used and delegate to GPU for dirichlet assembly
  auto use_hypre = options()->linearSystem.serviceName() == "HypreLinearSystem";
  if (use_hypre) {
    _assembleDirichletsGpu();
    return;
  }

  info() << "[ArcaneFem-Info] Started module _assembleDirichletsCpu()";

  BC::IArcaneFemBC* bc = options()->boundaryConditions();
  if (bc) {
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions()) {
      if (bs->getEnforceDirichletMethod() == "Penalty") {
        ArcaneFemFunctions::BoundaryConditions::applyDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);
      }
      else {
        if (t <= dt)
          ArcaneFemFunctions::BoundaryConditions::applyDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);
        else
          ArcaneFemFunctions::BoundaryConditions::applyDirichletToRhs(bs, node_dof, rhs_values);
      }
    }

    for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions()) {
      if (bs->getEnforceDirichletMethod() == "Penalty") {
        ArcaneFemFunctions::BoundaryConditions::applyPointDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);
      }
      else {
        if (t <= dt)
          ArcaneFemFunctions::BoundaryConditions::applyPointDirichletToLhsAndRhs(bs, node_dof, m_linear_system, rhs_values);
        else
          ArcaneFemFunctions::BoundaryConditions::applyPointDirichletToRhs(bs, node_dof, rhs_values);
      }
    }
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

void FemModule::_assembleDirichletsGpu()
{
  info() << "[ArcaneFem-Info] Started module  _assembleDirichletsGpu()";

  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  auto mesh_ptr = mesh();

  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  if (bc) {
    for (BC::IDirichletBoundaryCondition* bs : bc->dirichletBoundaryConditions()) {
      if (t <= dt)
        FemUtils::Gpu::BoundaryConditions::applyDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);
      else
        FemUtils::Gpu::BoundaryConditions::applyDirichletToRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);
    }

    for (BC::IDirichletPointCondition* bs : bc->dirichletPointConditions()) {
      if (t <= dt)
        FemUtils::Gpu::BoundaryConditions::applyPointDirichletToLhsAndRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);
      else
        FemUtils::Gpu::BoundaryConditions::applyPointDirichletToRhs(bs, m_dofs_on_nodes, m_linear_system, mesh_ptr, queue);
    }
  }
}