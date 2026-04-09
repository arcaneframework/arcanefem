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

inline void FemModuleElasticity::
_applyTraction(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  for (BC::ITractionBoundaryCondition* bs : bc->tractionBoundaryConditions()) {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh)
        ArcaneFemFunctions::BoundaryConditions2D::applyTractionToRhsQuad4(bs, node_dof, m_node_coord, rhs_values);
      else
        ArcaneFemFunctions::BoundaryConditions2D::applyTractionToRhsTria3(bs, node_dof, m_node_coord, rhs_values);
    }
    else if (mesh()->dimension() == 3) {
      if (m_hex_quad_mesh)
        ArcaneFemFunctions::BoundaryConditions3D::applyTractionToRhsHexa8(bs, node_dof, m_node_coord, rhs_values);
      else
        ArcaneFemFunctions::BoundaryConditions3D::applyTractionToRhsTetra4(bs, node_dof, m_node_coord, rhs_values);
    }
  }
}