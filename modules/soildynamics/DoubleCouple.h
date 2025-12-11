// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoubleCouple.h                                               (C) 2022-2025*/
/*                                                                           */
/* Functions to compute & assemble double-couple contribution to LHS/RHS     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies double-couple boundary conditions to the RHS values.
 * 
 * This function iterates over the defined double-couple boundary conditions,
 * retrieves the corresponding force values from the associated CaseTable,
 * and updates the RHS values for the specified nodes accordingly.
 * 
 * @param rhs_values The RHS values to be updated.
 * @param node_dof The node DoF connectivity view.
 * 
/*---------------------------------------------------------------------------*/

void FemModuleSoildynamics::
_applyDoubleCouple(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  // Index of the boundary condition. Needed to associate a CaseTable
  Int32 boundary_condition_index = 0;

  const Int8 NorthSouthId = 0; // X-Direction DoF imposed by double-couple
  const Int8 EastWestId = (mesh()->dimension() == 3) ? 2 : 1; // Y or Z-Direction DoF imposed by double-couple

  for (const auto& bs : options()->doubleCouple()) {

    const CaseTableInfo& case_table_dc_info = m_double_couple_case_table_list[boundary_condition_index];

    ++boundary_condition_index;

    Real dc_force; // double-couple force

    CaseTable* dc_case_table = case_table_dc_info.case_table;
    dc_case_table->value(t, dc_force);

    NodeGroup north = bs->northNodeName();
    NodeGroup south = bs->southNodeName();
    NodeGroup east = bs->eastNodeName();
    NodeGroup west = bs->westNodeName();

    ENUMERATE_ (Node, inode, north.own()) {
      Node node = *inode;
      DoFLocalId dof_id1 = node_dof.dofId(node, NorthSouthId);
      rhs_values[dof_id1] = dc_force;
    }
    ENUMERATE_ (Node, inode, south.own()) {
      Node node = *inode;
      DoFLocalId dof_id1 = node_dof.dofId(node, NorthSouthId);
      rhs_values[dof_id1] = -dc_force;
    }
    ENUMERATE_ (Node, inode, east.own()) {
      Node node = *inode;
      DoFLocalId dof_id2 = node_dof.dofId(node, EastWestId);
      rhs_values[dof_id2] = -dc_force;
    }
    ENUMERATE_ (Node, inode, west.own()) {
      Node node = *inode;
      DoFLocalId dof_id2 = node_dof.dofId(node, EastWestId);
      rhs_values[dof_id2] = dc_force;
    }
  }
}