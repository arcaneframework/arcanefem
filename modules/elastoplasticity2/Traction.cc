// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Traction.cc                                                 (C) 2000-2026 */
/*                                                                           */
/* Contains functions to compute and assemble traction contribution to RHS   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "modules/elastoplasticity2/Elastoplasticity2Module.h"

#include "femutils/ArcaneFemFunctions.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::ArcaneFem
{

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
 */
void Elastoplasticity2Module::
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
      if (m_hex_quad_mesh) {
        if (is_transient_traction) {
          if (m_nodes_per_cell == 4)
            transientTraction(_applyPressureTableToRhsTria3);
          else
            transientTraction(_applyPressureTableToRhsLine3);
        }
        else {
          if (m_nodes_per_cell == 4)
            constantTraction(ArcaneFemFunctions::BoundaryConditions2D::applyTractionToRhsQuad4);
          else
            constantTraction(_applyTractionToRhsLine3);
        }
      }
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
    ++boundary_condition_index;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies traction to the RHS vector for a Line3 element.
 *
 * This function computes the contribution of traction to the RHS vector
 * specifically for Line3 elements. It uses Gaussian quadrature to evaluate
 * the integral of the traction over the element and updates the RHS vector
 * accordingly.
 *
 * @param bs The traction boundary condition object containing traction values.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 *                 degrees of freedom (DoFs).
 * @param node_coord The coordinates of the mesh nodes.
 * @param rhs_values The variable representing the RHS vector to be updated.
 */
void Elastoplasticity2Module::
_applyTractionToRhsLine3(BC::ITractionBoundaryCondition* bs,
                         const IndexedNodeDoFConnectivityView& node_dof,
                         const VariableNodeReal3& node_coord,
                         VariableDoFReal& rhs_values)
{
  const StringConstArrayView values = bs->getValue();
  Real3 traction{};
  for (Int32 i = 0; i < values.size() && i < 2; ++i)
    if (values[i] != "NULL")
      traction[i] = std::stod(values[i].localstr());

  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
  ENUMERATE_ (Face, iface, bs->getSurface()) {
    Face face = *iface;
    Node nodes[3] = { face.node(0), face.node(1), face.node(2) };
    const Real3 coords[3] = { node_coord[nodes[0]], node_coord[nodes[1]], node_coord[nodes[2]] };
    for (Int8 igauss = 0; igauss < 3; ++igauss) {
      const Real xi = gp[igauss];
      const RealVector<3> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsLine3(xi);
      const Real dN[3] = { xi - 0.5, xi + 0.5, -2.0 * xi };
      Real dx_dxi = 0.0;
      Real dy_dxi = 0.0;
      for (Int8 i = 0; i < 3; ++i) {
        dx_dxi += dN[i] * coords[i].x;
        dy_dxi += dN[i] * coords[i].y;
      }
      const Real jacobian = math::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);
      if (jacobian <= 0.0)
        ARCANE_FATAL("Invalid (non-positive) Line3 Jacobian: {0}", jacobian);
      const Real integration_weight = weights[igauss] * jacobian;
      for (Int8 i = 0; i < 3; ++i) {
        if (!nodes[i].isOwn())
          continue;
        rhs_values[node_dof.dofId(nodes[i], 0)] += traction[0] * N[i] * integration_weight;
        rhs_values[node_dof.dofId(nodes[i], 1)] += traction[1] * N[i] * integration_weight;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/** 
 * @brief Applies pressure table to the RHS vector for a Line3 element.
 *
 * This function computes the contribution of a pressure table to the RHS vector
 * specifically for Line3 elements. It uses Gaussian quadrature to evaluate the
 * integral of the pressure over the element and updates the RHS vector accordingly.
 *
 * @param bs The traction boundary condition object containing pressure values.
 * @param t The current time for evaluating the pressure table.
 * @param boundary_condition_index The index of the boundary condition in the list.
 * @param traction_case_table_list The list of case tables for traction boundary conditions.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 *                degrees of freedom (DoFs).
 * @param node_coord The coordinates of the mesh nodes.
 * @param rhs_values The variable representing the RHS vector to be updated.
 *
/*---------------------------------------------------------------------------*/

void Elastoplasticity2Module::
_applyPressureTableToRhsLine3(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
                              const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list,
                              const IndexedNodeDoFConnectivityView& node_dof,
                              const VariableNodeReal3& node_coord,
                              VariableDoFReal& rhs_values)
{
  const auto& case_table_info = traction_case_table_list[boundary_condition_index];
  CaseTable* ct = case_table_info.case_table;
  if (!ct)
    ARCANE_FATAL("CaseTable is null. Maybe there is a missing call to _readCaseTables()");
  if (bs->getTractionInputFile() != case_table_info.file_name)
    ARCANE_FATAL("Incoherent CaseTable. The current CaseTable is associated to file '{0}'", case_table_info.file_name);

  Real3 traction;
  ct->value(t, traction);
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  ENUMERATE_ (Face, iface, bs->getSurface()) {
    Face face = *iface;
    Node nodes[3] = { face.node(0), face.node(1), face.node(2) };
    const Real3 coords[3] = { node_coord[nodes[0]], node_coord[nodes[1]], node_coord[nodes[2]] };
    const Real orientation = face.isSubDomainBoundaryOutside() ? 1.0 : -1.0;

    for (Int8 igauss = 0; igauss < 3; ++igauss) {
      const Real xi = gp[igauss];
      const RealVector<3> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsLine3(xi);
      const Real dN[3] = { xi - 0.5, xi + 0.5, -2.0 * xi };
      Real dx_dxi = 0.0;
      Real dy_dxi = 0.0;
      for (Int8 i = 0; i < 3; ++i) {
        dx_dxi += dN[i] * coords[i].x;
        dy_dxi += dN[i] * coords[i].y;
      }
      const Real jacobian = math::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);
      if (jacobian <= 0.0)
        ARCANE_FATAL("Invalid (non-positive) Line3 Jacobian: {0}", jacobian);
      const Real normal_x = orientation * dy_dxi / jacobian;
      const Real normal_y = orientation * -dx_dxi / jacobian;
      const Real integration_weight = weights[igauss] * jacobian;
      for (Int8 i = 0; i < 3; ++i) {
        if (!nodes[i].isOwn())
          continue;
        rhs_values[node_dof.dofId(nodes[i], 0)] -= traction[0] * normal_x * N[i] * integration_weight;
        rhs_values[node_dof.dofId(nodes[i], 1)] -= traction[1] * normal_y * N[i] * integration_weight;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies pressure table to the RHS vector for a Tri3 element.
 *
 * This function computes the contribution of a pressure table to the RHS vector
 * specifically for Tri3 elements. It uses Gaussian quadrature to evaluate the
 * integral of the pressure over the element and updates the RHS vector accordingly.
 *
 * @param bs The traction boundary condition object containing pressure values.
 * @param t The current time for evaluating the pressure table.
 * @param boundary_condition_index The index of the boundary condition in the list.
 * @param traction_case_table_list The list of case tables for traction boundary conditions.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 *                degrees of freedom (DoFs).
 * @param node_coord The coordinates of the mesh nodes.
 * @param rhs_values The variable representing the RHS vector to be updated.
 *
/*---------------------------------------------------------------------------*/

void Elastoplasticity2Module::
_applyPressureTableToRhsTria3(BC::ITractionBoundaryCondition* bs, const Real t,
                              Int32 boundary_condition_index,
                              const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list,
                              const IndexedNodeDoFConnectivityView& node_dof,
                              const VariableNodeReal3& node_coord,
                              VariableDoFReal& rhs_values)
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::ArcaneFem

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
