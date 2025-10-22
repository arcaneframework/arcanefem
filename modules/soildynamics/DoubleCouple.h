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
 * @brief Applies paraxial boundary conditions to the RHS and LHS
 * of the linear system.
 * 
 * @param rhs_values The RHS values to be updated.
 * @param node_dof The node DoF connectivity view.
 * @param bsr_matrix The BSR matrix to be updated (optional).
 * 
/*---------------------------------------------------------------------------*/

inline void FemModule::
_applyDoubleCouple(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  // Index of the boundary condition. Needed to associate a CaseTable
  Int32 boundary_condition_index_dc = 0;

  Int8 NorthSouthId = 0;
  Int8 EastWestId = 1;

  if (mesh()->dimension() == 3) {
    EastWestId = 2;
  }

  for (const auto& bs : options()->doubleCouple()) {

    const CaseTableInfo& case_table_dc_info = m_double_couple_case_table_list[boundary_condition_index_dc];

    ++boundary_condition_index_dc;

    Real dc_force; // double-couple force

    String file_name = bs->doubleCoupleInputFile();
    CaseTable* dc_case_table_inn = case_table_dc_info.case_table;
    dc_case_table_inn->value(t, dc_force);

    NodeGroup north = bs->northNodeName();
    NodeGroup south = bs->southNodeName();
    NodeGroup east = bs->eastNodeName();
    NodeGroup west = bs->westNodeName();

    ENUMERATE_ (Node, inode, north) {
      Node node = *inode;
      DoFLocalId dof_id1 = node_dof.dofId(node, NorthSouthId);
      rhs_values[dof_id1] = dc_force;
    }
    ENUMERATE_ (Node, inode, south) {
      Node node = *inode;
      DoFLocalId dof_id1 = node_dof.dofId(node, NorthSouthId);
      rhs_values[dof_id1] = -dc_force;
    }
    ENUMERATE_ (Node, inode, east) {
      Node node = *inode;
      DoFLocalId dof_id2 = node_dof.dofId(node, EastWestId);
      rhs_values[dof_id2] = -dc_force;
    }
    ENUMERATE_ (Node, inode, west) {
      Node node = *inode;
      DoFLocalId dof_id2 = node_dof.dofId(node, EastWestId);
      rhs_values[dof_id2] = dc_force;
    }
  }
}