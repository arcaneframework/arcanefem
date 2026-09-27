// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BoundaryConditionsHelpers.cc                                (C) 2000-2026 */
/*                                                                           */
/* Helper functions for handling boundary conditions.                        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "ArcaneFemFunctions.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditionsHelpers::
applyDirichletToNodeGroupRhsOnly(const Int32 dof_index, Real value,
                                 const IndexedNodeDoFConnectivityView& node_dof,
                                 VariableDoFReal& rhs_values,
                                 NodeGroup& node_group)
{
  ENUMERATE_ (Node, inode, node_group) {
    Node node = *inode;
    if (node.isOwn()) {
      rhs_values[node_dof.dofId(node, dof_index)] = value;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditionsHelpers::
applyDirichletToNodeGroupViaPenalty(const Int32 dof_index, Real value, Real penalty,
                                    const IndexedNodeDoFConnectivityView& node_dof,
                                    DoFLinearSystem& linear_system,
                                    VariableDoFReal& rhs_values,
                                    NodeGroup& node_group)
{
  ENUMERATE_ (Node, inode, node_group) {
    Node node = *inode;
    if (node.isOwn()) {
      linear_system.matrixSetValue(node_dof.dofId(node, dof_index), node_dof.dofId(node, dof_index), penalty);
      Real u_g = penalty * value;
      rhs_values[node_dof.dofId(node, dof_index)] = u_g;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditionsHelpers::
applyDirichletToNodeGroupViaRowElimination(const Int32 dof_index, Real value,
                                           const IndexedNodeDoFConnectivityView& node_dof,
                                           DoFLinearSystem& linear_system,
                                           VariableDoFReal& rhs_values,
                                           NodeGroup& node_group)
{
  DoFLinearSystemRowEliminationHelper elimination_helper(linear_system.rowEliminationHelper());
  ENUMERATE_ (Node, inode, node_group) {
    Node node = *inode;
    if (node.isOwn()) {
      elimination_helper.addElimination(node_dof.dofId(*inode, dof_index), value);
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditionsHelpers::
applyDirichletToNodeGroupViaRowColumnElimination(const Int32 dof_index, Real value,
                                                 const IndexedNodeDoFConnectivityView& node_dof,
                                                 DoFLinearSystem& linear_system,
                                                 VariableDoFReal& rhs_values,
                                                 NodeGroup& node_group)
{
  DoFLinearSystemRowColumnEliminationHelper elimination_helper(linear_system.rowColumnEliminationHelper());

  ENUMERATE_ (Node, inode, node_group) {
    Node node = *inode;
    if (node.isOwn()) {
      elimination_helper.addElimination(node_dof.dofId(*inode, dof_index), value);
    }
  }
}

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/
