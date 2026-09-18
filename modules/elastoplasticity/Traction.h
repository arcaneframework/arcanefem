// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Traction.h                                                  (C) 2000-2026 */
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

inline void FemModuleElastoplasticity::
_applyTraction(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  Int32 boundary_condition_index = 0;
  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  for (BC::ITractionBoundaryCondition* bs : bc->tractionBoundaryConditions()) {
    const auto traction_table_file_name = bs->getTractionInputFile();
    const bool is_transient_traction = !traction_table_file_name.empty();

    auto transientTraction = [&](auto fn) { fn(bs, t, boundary_condition_index, m_traction_case_table_list, node_dof, m_node_coord, rhs_values); };
    auto constantTraction = [&](auto fn) { fn(bs, node_dof, m_node_coord, rhs_values); };

    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh)
        is_transient_traction ? transientTraction(ArcaneFemFunctions::BoundaryConditions2D::applyTractionTableToRhsQuad4)
                              : constantTraction(ArcaneFemFunctions::BoundaryConditions2D::applyTractionToRhsQuad4);
      else
        is_transient_traction ? transientTraction(_applyPressureTableToRhsTria3)
                              : constantTraction(ArcaneFemFunctions::BoundaryConditions2D::applyTractionToRhsTria3);
    }
    else if (mesh()->dimension() == 3) {
      if (m_hex_quad_mesh)
        is_transient_traction ? transientTraction(ArcaneFemFunctions::BoundaryConditions3D::applyTractionTableToRhsHexa8)
                              : constantTraction(ArcaneFemFunctions::BoundaryConditions3D::applyTractionToRhsHexa8);
      else
        is_transient_traction ? transientTraction(ArcaneFemFunctions::BoundaryConditions3D::applyTractionTableToRhsTetra4)
                              : constantTraction(ArcaneFemFunctions::BoundaryConditions3D::applyTractionToRhsTetra4);
    }
  }
  ++boundary_condition_index;
}


inline void FemModuleElastoplasticity::
_applyPressureTableToRhsTria3(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index, const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
{
  // mesh boundary group on which traction is applied
  FaceGroup group = bs->getSurface();

  bool applyTraction = false;
  Real3 trac;
  auto traction_table_file_name = bs->getTractionInputFile();
  bool getTractionFromTable = !traction_table_file_name.empty();

  if (getTractionFromTable) {

    const Arcane::FemUtils::CaseTableInfo& case_table_info = traction_case_table_list[boundary_condition_index];
    applyTraction = true;

    CaseTable* ct = case_table_info.case_table;
    if (!ct)
      ARCANE_FATAL("CaseTable is null. Maybe there is a missing call to _readCaseTables()");
    if (traction_table_file_name != case_table_info.file_name)
      ARCANE_FATAL("Incoherent CaseTable. The current CaseTable is associated to file '{0}'", case_table_info.file_name);

    ct->value(t, trac);
  }

  // no traction to apply hence return
  if (!applyTraction)
    return;

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;
    Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, node_coord);
    Real2 normal = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, node_coord);
    for (Node node : iface->nodes()) {
      if (node.isOwn()) {
        // The table stores the pressure magnitude in its x/y columns. Internal
        // pressure acts opposite to the outward normal of the hollow domain.
        rhs_values[node_dof.dofId(node, 0)] -= trac[0] * normal.x * length / 2.;
        rhs_values[node_dof.dofId(node, 1)] -= trac[1] * normal.y * length / 2.;
      }
    }
  }
}