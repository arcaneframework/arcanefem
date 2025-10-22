// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Traction.h                                                  (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute and assemble traction contribution to RHS   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies traction to the RHS vector of the linear system.
 * 
 * This function  computes  the  contribution of traction to the RHS vector 
 * of the linear system. It iterates over all cells in the mesh, calculates 
 * the appropriate force  contributions  based on the element type and mesh 
 * dimension, and updates the RHS vector accordingly.
 * 
 * traction term ∫∫ (𝐭.𝐯)  with 𝐭 = (𝑡𝑥, 𝑡𝑦, 𝑡𝑧) = (t[0], t[1], t[2])
 * 
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding 
 *                 degrees of freedom (DoFs).
 * 
/*---------------------------------------------------------------------------*/

void FemModule::
_applyTraction(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  Int32 boundary_condition_index = 0;
  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  if (bc) {
    for (BC::ITractionBoundaryCondition* bs : bc->tractionBoundaryConditions()) {
      const auto traction_table_file_name = bs->getTractionInputFile();
      const bool is_transient_traction = !traction_table_file_name.empty();

      auto transientTraction = [&](auto fn) { fn(bs, t, boundary_condition_index, m_traction_case_table_list, node_dof, m_node_coord, rhs_values); };
      auto constantTraction = [&](auto fn) { fn(bs, node_dof, m_node_coord, rhs_values); };

      if (mesh()->dimension() == 2) {
        is_transient_traction ? transientTraction(ArcaneFemFunctions::BoundaryConditions2D::applyTractionTableToRhsTria3)
                              : constantTraction(ArcaneFemFunctions::BoundaryConditions2D::applyTractionToRhsTria3);
      }
      else if (mesh()->dimension() == 3) {
        is_transient_traction ? transientTraction(ArcaneFemFunctions::BoundaryConditions3D::applyTractionTableToRhsTetra4)
                              : constantTraction(ArcaneFemFunctions::BoundaryConditions3D::applyTractionToRhsTetra4);
      }
    }
  }
}