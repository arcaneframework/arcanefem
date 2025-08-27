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

inline void FemModule::
_applyTraction(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
  // Index of the boundary condition. Needed to associate a CaseTable
  Int32 boundary_condition_index = 0;

  for (const auto& bs : options()->tractionBoundaryCondition()) {

    FaceGroup group = bs->surface();
    

    // get traction components
    Real3 trac = { 0, 0, 0 };
    {
      bool getTractionFromTable = bs->tractionInputFile.isPresent();

      if (getTractionFromTable) {

        const CaseTableInfo& case_table_info = m_traction_case_table_list[boundary_condition_index++];
        String file_name = bs->tractionInputFile();
        info() << "Applying traction boundary conditions for surface " << group.name() << " via CaseTable" << file_name;

        CaseTable* ct = case_table_info.case_table;
        if (!ct)
          ARCANE_FATAL("CaseTable is null. Maybe there is a missing call to _readCaseTables()");
        if (file_name != case_table_info.file_name)
          ARCANE_FATAL("Incoherent CaseTable. The current CaseTable is associated to file '{0}'", case_table_info.file_name);

        ct->value(t, trac);
      }
      else {

        const UniqueArray<String> t_string = bs->t();

        for (Int32 i = 0; i < t_string.size(); ++i) {
          trac[i] = 0.0;
          if (t_string[i] != "NULL") {
            trac[i] = std::stod(t_string[i].localstr());
          }
        }
      }
    }

    // Only apply traction if any component is active
    if (trac.x == 0.0 && trac.y == 0.0 && trac.z == 0.0)
      return;

    info() << "[ArcaneFem-Info] Applying Traction " << trac;
    info() << "[ArcaneFem-Info] Traction surface '" << bs->surface().name() << "'";

    if (mesh()->dimension() == 2){
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, m_node_coord);
        for (Node node : iface->nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += trac[0] * length / 2.;
            rhs_values[node_dof.dofId(node, 1)] += trac[1] * length / 2.;
          }
        }
      }
    }

    if (mesh()->dimension() == 3){
      ENUMERATE_ (Face, iface, group) {
        Face face = *iface;
        Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, m_node_coord);
        for (Node node : iface->nodes()) {
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += trac[0] * area / 3.;
            rhs_values[node_dof.dofId(node, 1)] += trac[1] * area / 3.;
            rhs_values[node_dof.dofId(node, 2)] += trac[2] * area / 3.;
          }
        }
      }
    }
  }
}