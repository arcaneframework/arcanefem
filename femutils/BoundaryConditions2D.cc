// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BoundaryConditions2D.cc                                     (C) 2000-2026 */
/*                                                                           */
/* Helper functions for handling 2D boundary conditions.                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "ArcaneFemFunctions.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector.
 *
 * This method adds a constant source term `qdot` to the RHS vector for each
 * node in the mesh. The contribution to each node is weighted by the area of
 * the cell and evenly distributed among the number of nodes of the cell.
 *
 * @param qdot The constant source term.
 * @param mesh The mesh containing all cells.
 * @param node_dof DOF connectivity view.
 * @param node_coord The coordinates of the nodes.
 * @param rhs_values The RHS values to update.
 */
void ArcaneFemFunctions::BoundaryConditions2D::
applyConstantSourceToRhsTria3(Real qdot, IMesh* mesh,
                              const IndexedNodeDoFConnectivityView& node_dof,
                              const VariableNodeReal3& node_coord,
                              VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, node_coord);
    for (Node node : cell.nodes()) {
      if (node.isOwn())
        rhs_values[node_dof.dofId(node, 0)] += qdot * area / cell.nbNode();
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector for Quad4 elements.
 *
 * Uses a 2x2 Gauss rule to integrate the biquadratic (serendipity) shape
 * functions exactly. For each Gauss point the shape functions, their
 * derivatives, the Jacobian and its determinant are computed, then the
 * weighted contribution N[i]*qdot*detJ is scattered onto the owned nodes.
 */
void ArcaneFemFunctions::BoundaryConditions2D::
applyConstantSourceToRhsQuad4(Real qdot, IMesh* mesh,
                              const IndexedNodeDoFConnectivityView& node_dof,
                              const VariableNodeReal3& node_coord,
                              VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;
    // Real area = ArcaneFemFunctions::MeshOperation::computeAreaQuad4(cell, node_coord);
    // for (Node node : cell.nodes()) {
    //   if (node.isOwn())
    //     rhs_values[node_dof.dofId(node, 0)] += qdot * area / cell.nbNode();
    // }

    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 };
    constexpr Real weights[2] = { 1.0, 1.0 };

    for (Int32 ixi = 0; ixi < 2; ++ixi) {
      for (Int32 ieta = 0; ieta < 2; ++ieta) {

        // Get the coordinates of the Gauss point
        Real xi = gp[ixi]; // Get the ξ coordinate of the Gauss point
        Real eta = gp[ieta]; // Get the η coordinate of the Gauss point
        Real weight = weights[ixi] * weights[ieta];

        // Shape functions 𝐍 for Quad4
        RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);

        // Shape function derivatives ∂𝐍/∂ξ and ∂𝐍/∂η
        //     ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ ]
        //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η ]
        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad4(xi, eta);

        // Jacobian calculation 𝑱
        //    𝑱 = [ 𝒋₀₀  𝒋₀₁ ] = [ ∂𝑥/∂ξ  ∂𝑦/∂ξ ]
        //        [ 𝒋₁₀  𝒋₁₁ ]   [ ∂𝑥/∂η  ∂𝑦/∂η ]
        //
        // The Jacobian is computed as follows:
        //   𝒋₀₀ = ∑ (∂𝑁ᵢ/∂ξ * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟒
        //   𝒋₀₁ = ∑ (∂𝑁ᵢ/∂ξ * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟒
        //   𝒋₁₀ = ∑ (∂𝑁ᵢ/∂η * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟒
        //   𝒋₁₁ = ∑ (∂𝑁ᵢ/∂η * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟒

        Real J00 = 0, J01 = 0, J10 = 0, J11 = 0;
        for (Int8 a = 0; a < 4; ++a) {
          J00 += reference_gradients.dN_dxi[a] * node_coord[cell.nodeId(a)].x;
          J01 += reference_gradients.dN_dxi[a] * node_coord[cell.nodeId(a)].y;
          J10 += reference_gradients.dN_deta[a] * node_coord[cell.nodeId(a)].x;
          J11 += reference_gradients.dN_deta[a] * node_coord[cell.nodeId(a)].y;
        }

        // Determinant of the Jacobian
        Real detJ = J00 * J11 - J01 * J10;

        // Compute integration weight
        Real integration_weight = weight * detJ;

        // Assemble RHS
        for (Int32 i = 0; i < 4; ++i) {
          Node node = cell.node(i);
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += N[i] * qdot * integration_weight;
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector for Quad8 elements.
 *
 * Uses a 3x3 Gauss rule to integrate the biquadratic (serendipity) shape
 * functions exactly. For each Gauss point the shape functions, their
 * derivatives, the Jacobian and its determinant are computed, then the
 * weighted contribution N[i]*qdot*detJ is scattered onto the owned nodes.
 */
void ArcaneFemFunctions::BoundaryConditions2D::
applyConstantSourceToRhsQuad8(Real qdot, IMesh* mesh,
                              const IndexedNodeDoFConnectivityView& node_dof,
                              const VariableNodeReal3& node_coord,
                              VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;

    // 3-point Gauss rule per direction (needed for exact integration of quadratic Quad8 shape functions)
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(5/9) , 0 , sqrt(5/9)]
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {

        // Get the coordinates of the Gauss point
        Real xi = gp[ixi]; // Get the ξ coordinate of the Gauss point
        Real eta = gp[ieta]; // Get the η coordinate of the Gauss point
        Real weight = weights[ixi] * weights[ieta];

        // Shape functions 𝐍 for Quad8 (serendipity)
        RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad8(xi, eta);

        // Shape function derivatives ∂𝐍/∂ξ and ∂𝐍/∂η
        //     ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ  ∂𝑁₅/∂ξ  ∂𝑁₆/∂ξ  ∂𝑁₇/∂ξ  ∂𝑁₈/∂ξ ]
        //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η  ∂𝑁₅/∂η  ∂𝑁₆/∂η  ∂𝑁₇/∂η  ∂𝑁₈/∂η ]
        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad8(xi, eta);
        // Jacobian calculation 𝑱
        //    𝑱 = [ 𝒋₀₀  𝒋₀₁ ] = [ ∂𝑥/∂ξ  ∂𝑦/∂ξ ]
        //        [ 𝒋₁₀  𝒋₁₁ ]   [ ∂𝑥/∂η  ∂𝑦/∂η ]
        //
        // The Jacobian is computed as follows:
        //   𝒋₀₀ = ∑ (∂𝑁ᵢ/∂ξ * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₀₁ = ∑ (∂𝑁ᵢ/∂ξ * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₁₀ = ∑ (∂𝑁ᵢ/∂η * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟖
        //   𝒋₁₁ = ∑ (∂𝑁ᵢ/∂η * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟖

        Real J00 = 0, J01 = 0, J10 = 0, J11 = 0;
        for (Int8 a = 0; a < 8; ++a) {
          J00 += reference_gradients.dN_dxi[a] * node_coord[cell.nodeId(a)].x;
          J01 += reference_gradients.dN_dxi[a] * node_coord[cell.nodeId(a)].y;
          J10 += reference_gradients.dN_deta[a] * node_coord[cell.nodeId(a)].x;
          J11 += reference_gradients.dN_deta[a] * node_coord[cell.nodeId(a)].y;
        }

        // Determinant of the Jacobian
        Real detJ = J00 * J11 - J01 * J10;

        // Compute integration weight
        Real integration_weight = weight * detJ;

        // Assemble RHS
        for (Int32 i = 0; i < 8; ++i) {
          Node node = cell.node(i);
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += N[i] * qdot * integration_weight;
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector for Quad9 elements.
 *
 * Uses a 3x3 Gauss rule to integrate the biquadratic (Lagrange) shape
 * functions exactly. For each Gauss point the shape functions, their
 * derivatives, the Jacobian and its determinant are computed, then the
 * weighted contribution N[i]*qdot*detJ is scattered onto the owned nodes.
 */
/*---------------------------------------------------------------------------*/
void ArcaneFemFunctions::BoundaryConditions2D::
applyConstantSourceToRhsQuad9(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;

    // 3-point Gauss rule per direction (needed for exact integration of quadratic Quad9 shape functions)
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(5/9) , 0 , sqrt(5/9)]
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {

        // Get the coordinates of the Gauss point
        Real xi = gp[ixi]; // Get the ξ coordinate of the Gauss point
        Real eta = gp[ieta]; // Get the η coordinate of the Gauss point
        Real weight = weights[ixi] * weights[ieta];

        // Shape functions 𝐍 for Quad9 (Lagrange)
        RealVector<9> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad9(xi, eta);

        // Shape function derivatives ∂𝐍/∂ξ and ∂𝐍/∂η
        //     ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ  ∂𝑁₅/∂ξ  ∂𝑁₆/∂ξ  ∂𝑁₇/∂ξ  ∂𝑁₈/∂ξ  ∂𝑁₉/∂ξ ]
        //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η  ∂𝑁₅/∂η  ∂𝑁₆/∂η  ∂𝑁₇/∂η  ∂𝑁₈/∂η  ∂𝑁₉/∂η ]
        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad9(xi, eta);
        // Jacobian calculation 𝑱
        //    𝑱 = [ 𝒋₀₀  𝒋₀₁ ] = [ ∂𝑥/∂ξ  ∂𝑦/∂ξ ]
        //        [ 𝒋₁₀  𝒋₁₁ ]   [ ∂𝑥/∂η  ∂𝑦/∂η ]
        //
        // The Jacobian is computed as follows:
        //   𝒋₀₀ = ∑ (∂𝑁ᵢ/∂ξ * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟗
        //   𝒋₀₁ = ∑ (∂𝑁ᵢ/∂ξ * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟗
        //   𝒋₁₀ = ∑ (∂𝑁ᵢ/∂η * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟗
        //   𝒋₁₁ = ∑ (∂𝑁ᵢ/∂η * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟗

        Real J00 = 0, J01 = 0, J10 = 0, J11 = 0;
        for (Int8 a = 0; a < 9; ++a) {
          J00 += reference_gradients.dN_dxi[a] * node_coord[cell.nodeId(a)].x;
          J01 += reference_gradients.dN_dxi[a] * node_coord[cell.nodeId(a)].y;
          J10 += reference_gradients.dN_deta[a] * node_coord[cell.nodeId(a)].x;
          J11 += reference_gradients.dN_deta[a] * node_coord[cell.nodeId(a)].y;
        }

        // Determinant of the Jacobian
        Real detJ = J00 * J11 - J01 * J10;

        // Compute integration weight
        Real integration_weight = weight * detJ;

        // Assemble RHS
        for (Int32 i = 0; i < 9; ++i) {
          Node node = cell.node(i);
          if (node.isOwn()) {
            rhs_values[node_dof.dofId(node, 0)] += N[i] * qdot * integration_weight;
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a nodal field to the RHS vector.
 *
 * @param field The field term defined on nodes.
 * @param mesh  The mesh containing all cells.
 * @param node_dof DOF connectivity view.
 * @param node_coord The coordinates of the nodes.
 * @param rhs_values The RHS values to update.
 */

void ArcaneFemFunctions::BoundaryConditions2D::
integrateNodalFieldToRhsTria3(VariableNodeReal& field, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, node_coord);

    // Get nodal values for this triangular cell
    const Real field_at_nodes[3] = {
      field[cell.nodeId(0)],
      field[cell.nodeId(1)],
      field[cell.nodeId(2)]
    };

    // Apply mass matrix integration
    Real node_contributions[3] = { 0.0, 0.0, 0.0 };

    for (Int8 i = 0; i < 3; ++i) {
      for (Int8 j = 0; j < 3; ++j) {
        Real mass_coeff;
        if (i == j) {
          mass_coeff = area / 6.0; // diagonal: area * (2/12) = area/6
        }
        else {
          mass_coeff = area / 12.0; // off-diagonal: area * (1/12)
        }
        node_contributions[i] += mass_coeff * field_at_nodes[j];
      }
    }

    // Add contributions to global RHS
    for (Int8 i = 0; i < 3; ++i) {
      Node node = cell.node(i);
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += node_contributions[i];
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a nodal field to the RHS vector.
 *
 * @param field The field term defined on nodes.
 * @param mesh  The mesh containing all cells.
 * @param node_dof DOF connectivity view.
 * @param node_coord The coordinates of the nodes.
 * @param rhs_values The RHS values to update.
 */
void ArcaneFemFunctions::BoundaryConditions2D::
integrateNodalFieldToRhsQuad4(VariableNodeReal& field, IMesh* mesh,
                              const IndexedNodeDoFConnectivityView& node_dof,
                              const VariableNodeReal3& node_coord,
                              VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;

    // Get nodal values of field for this cell
    const Real field_at_nodes[4] = {
      field[cell.nodeId(0)],
      field[cell.nodeId(1)],
      field[cell.nodeId(2)],
      field[cell.nodeId(3)]
    };

    // Initialize contributions for each node in this cell
    Real node_contributions[4] = { 0.0, 0.0, 0.0, 0.0 };

    // 2x2 Gauss integration for quadrilateral element
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
    constexpr Real w = 1.0;

    for (Int8 ixi = 0; ixi < 2; ++ixi) {
      for (Int8 ieta = 0; ieta < 2; ++ieta) {

        // Get the coordinates of the Gauss point
        Real xi = gp[ixi]; // Get the ξ coordinate of the Gauss point
        Real eta = gp[ieta]; // Get the η coordinate of the Gauss point
        Real weight = w * w; // Weight for 2D Gauss integration

        // Shape functions 𝐍 for Quad4
        RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);

        // compute the det(Jacobian)
        const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, node_coord, xi, eta);
        const Real detJ = gp_info.det_j;

        // compute integration weight
        const Real integration_weight = weight * detJ;

        // Interpolate qdot at the quadrature point: qdot_gp = ∑ 𝑁ᵢ * q̇
        Real qdot_gp = 0.0;
        for (Int8 a = 0; a < 4; ++a) {
          qdot_gp += N[a] * field_at_nodes[a];
        }

        // Add contribution to each test function: ∫ q̇ * 𝑁ᵢ dΩ
        for (Int8 i = 0; i < 4; ++i) {
          node_contributions[i] += qdot_gp * N[i] * integration_weight;
        }
      }
    }

    // Add contributions to global RHS vector
    for (Int8 i = 0; i < 4; ++i) {
      Node node = cell.node(i);
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += node_contributions[i];
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a manufactured source term to the RHS vector.
 *
 * This method adds a manufactured source term to the RHS vector for each
 * node in the mesh. The contribution to each node is weighted by the area of
 * the cell and evenly distributed among the nodes of the cell.
 *
 * @param qdot The constant source term.
 * @param mesh The mesh containing all cells.
 * @param node_dof DOF connectivity view.
 * @param node_coord The coordinates of the nodes.
 * @param rhs_values The RHS values to update.
 */

void ArcaneFemFunctions::BoundaryConditions2D::
applyManufacturedSourceToRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_source,
                             IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof,
                             const VariableNodeReal3& node_coord,
                             VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, node_coord);
    Real3 bcenter = ArcaneFemFunctions::MeshOperation::computeBaryCenterTria3(cell, node_coord);

    for (Node node : cell.nodes()) {
      if (node.isOwn())
        rhs_values[node_dof.dofId(node, 0)] += manufactured_source->apply(area / cell.nbNode(), bcenter);
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Neumann conditions to the right-hand side (RHS) values.
 *
 * This method updates the RHS values of the finite element method equations
 * based on the provided Neumann boundary condition. The boundary condition
 * can specify a value or its components along the x and y directions.
 *
 * @param bs The Neumann boundary condition values.
 * @param node_dof Connectivity view for degrees of freedom at nodes.
 * @param node_coord Coordinates of the nodes in the mesh.
 * @param rhs_values The right-hand side values to be updated.
 */
void ArcaneFemFunctions::BoundaryConditions2D::
applyNeumannToRhsTria3(BC::INeumannBoundaryCondition* bs,
                       const IndexedNodeDoFConnectivityView& node_dof,
                       const VariableNodeReal3& node_coord,
                       VariableDoFReal& rhs_values)
{
  FaceGroup group = bs->getSurface();

  Real value = 0.0;
  Real valueX = 0.0;
  Real valueY = 0.0;

  bool scalarNeumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalarNeumann = true;
    value = std::stod(neumann_str[0].localstr());
  }
  else {
    if (neumann_str.size() > 1) {
      if (neumann_str[0] != "NULL")
        valueX = std::stod(neumann_str[0].localstr());
      if (neumann_str[1] != "NULL")
        valueY = std::stod(neumann_str[1].localstr());
    }
  }

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;

    Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, node_coord);
    Real2 normal = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, node_coord);

    for (Node node : iface->nodes()) {
      if (!node.isOwn())
        continue;
      Real rhs_value;

      if (scalarNeumann) {
        rhs_value = value * length / 2.0;
      }
      else {
        rhs_value = (normal.x * valueX + normal.y * valueY) * length / 2.0;
      }

      rhs_values[node_dof.dofId(node, 0)] += rhs_value;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions2D::
applyNeumannToRhsQuad4(BC::INeumannBoundaryCondition* bs,
                       const IndexedNodeDoFConnectivityView& node_dof,
                       const VariableNodeReal3& node_coord,
                       VariableDoFReal& rhs_values)
{
  FaceGroup group = bs->getSurface();

  Real value = 0.0;
  Real valueX = 0.0;
  Real valueY = 0.0;

  bool scalarNeumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalarNeumann = true;
    value = std::stod(neumann_str[0].localstr());
  }
  else {
    if (neumann_str.size() > 1) {
      if (neumann_str[0] != "NULL")
        valueX = std::stod(neumann_str[0].localstr());
      if (neumann_str[1] != "NULL")
        valueY = std::stod(neumann_str[1].localstr());
    }
  }

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;

    // 2-point Gauss integration for line element
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
    constexpr Real weights[2] = { 1.0, 1.0 };

    Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, node_coord);
    Real2 normal = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, node_coord);

    Node node0 = face.node(0);
    Node node1 = face.node(1);

    for (Int32 i = 0; i < 2; ++i) {
      Real xi = gp[i];
      Real weight = weights[i];

      // Linear shape functions for Line2
      RealVector<2> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsLine2(xi);

      // Integration weight: weight * jacobian (length/2 for reference element [-1,1])
      Real integration_weight = weight * length * 0.5;

      // Apply to both nodes
      Node nodes[2] = { node0, node1 };
      for (Int32 j = 0; j < 2; ++j) {
        Node node = nodes[j];
        if (!node.isOwn())
          continue;

        Real rhs_value;
        if (scalarNeumann) {
          rhs_value = value * N[j] * integration_weight;
        }
        else {
          rhs_value = (normal.x * valueX + normal.y * valueY) * N[j] * integration_weight;
        }

        rhs_values[node_dof.dofId(node, 0)] += rhs_value;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a Neumann condition on a quadratic Line3 face.
 *
 * Uses three-point Gauss integration and an isoparametric Line3 mapping.
 * This supports both scalar fluxes and vector fluxes projected onto the
 * outward normal, including curved quadratic edges.
 */
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions2D::
applyNeumannToRhsLine3(BC::INeumannBoundaryCondition* bs,
                       const IndexedNodeDoFConnectivityView& node_dof,
                       const VariableNodeReal3& node_coord,
                       VariableDoFReal& rhs_values)
{
  FaceGroup group = bs->getSurface();

  Real value = 0.0;
  Real valueX = 0.0;
  Real valueY = 0.0;
  bool scalar_neumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalar_neumann = true;
    value = std::stod(neumann_str[0].localstr());
  }
  else if (neumann_str.size() > 1) {
    if (neumann_str[0] != "NULL")
      valueX = std::stod(neumann_str[0].localstr());
    if (neumann_str[1] != "NULL")
      valueY = std::stod(neumann_str[1].localstr());
  }

  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };
  constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;
    if (face.nbNode() != 3)
      ARCANE_FATAL("Expected a Line3 face for quadratic quadrilateral Neumann assembly, got '{0}' nodes", face.nbNode());

    Node nodes[3] = { face.node(0), face.node(1), face.node(2) };
    Real3 coords[3] = { node_coord[nodes[0]], node_coord[nodes[1]], node_coord[nodes[2]] };
    const Real orientation = face.isSubDomainBoundaryOutside() ? 1.0 : -1.0;

    for (Int32 igauss = 0; igauss < 3; ++igauss) {
      const Real xi = gp[igauss];
      RealVector<3> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsLine3(xi);

      const Real dN[3] = { xi - 0.5, xi + 0.5, -2.0 * xi };

      Real dx_dxi = 0.0;
      Real dy_dxi = 0.0;
      for (Int32 i = 0; i < 3; ++i) {
        dx_dxi += dN[i] * coords[i].x;
        dy_dxi += dN[i] * coords[i].y;
      }

      const Real jacobian = math::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);
      if (jacobian <= 0.0)
        ARCANE_FATAL("Invalid (non-positive) Line3 Jacobian: {0}", jacobian);

      const Real normal_x = orientation * dy_dxi / jacobian;
      const Real normal_y = orientation * -dx_dxi / jacobian;
      const Real flux = scalar_neumann ? value : normal_x * valueX + normal_y * valueY;
      const Real integration_weight = weights[igauss] * jacobian;

      for (Int32 i = 0; i < 3; ++i) {
        Node node = nodes[i];
        if (node.isOwn())
          rhs_values[node_dof.dofId(node, 0)] += flux * N[i] * integration_weight;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions2D::
applyNeumannToRhsQuad8(BC::INeumannBoundaryCondition* bs,
                       const IndexedNodeDoFConnectivityView& node_dof,
                       const VariableNodeReal3& node_coord,
                       VariableDoFReal& rhs_values)
{
  applyNeumannToRhsLine3(bs, node_dof, node_coord, rhs_values);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

//! Quad9 has the same quadratic Line3 boundary interpolation as Quad8.
void ArcaneFemFunctions::BoundaryConditions2D::
applyNeumannToRhsQuad9(BC::INeumannBoundaryCondition* bs,
                       const IndexedNodeDoFConnectivityView& node_dof,
                       const VariableNodeReal3& node_coord,
                       VariableDoFReal& rhs_values)
{
  applyNeumannToRhsLine3(bs, node_dof, node_coord, rhs_values);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies traction conditions to the right-hand side (RHS) values.
 *
 * @param bs The traction boundary condition values.
 * @param node_dof Connectivity view for degrees of freedom at nodes.
 * @param node_coord Coordinates of the nodes in the mesh.
 * @param rhs_values The right-hand side values to be updated.
 */
void ArcaneFemFunctions::BoundaryConditions2D::
applyTractionToRhsTria3(BC::ITractionBoundaryCondition* bs,
                        const IndexedNodeDoFConnectivityView& node_dof,
                        const VariableNodeReal3& node_coord,
                        VariableDoFReal& rhs_values)
{
  // mesh boundary group on which traction is applied
  FaceGroup group = bs->getSurface();

  Real3 t;

  // get traction force vector
  bool applyTraction = false;
  const UniqueArray<String> t_string = bs->getValue();
  for (Int32 i = 0; i < t_string.size(); ++i) {
    t[i] = 0.0;
    if (t_string[i] != "NULL") {
      applyTraction = true;
      t[i] = std::stod(t_string[i].localstr());
    }
  }

  // no traction to apply hence return
  if (!applyTraction)
    return;

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;
    Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, node_coord);

    for (Node node : iface->nodes()) {
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += t[0] * length / 2.;
        rhs_values[node_dof.dofId(node, 1)] += t[1] * length / 2.;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions2D::
applyTractionToRhsQuad4(BC::ITractionBoundaryCondition* bs,
                        const IndexedNodeDoFConnectivityView& node_dof,
                        const VariableNodeReal3& node_coord,
                        VariableDoFReal& rhs_values)
{
  // mesh boundary group on which traction is applied
  FaceGroup group = bs->getSurface();

  Real3 t;
  // get traction force vector
  bool applyTraction = false;
  const UniqueArray<String> t_string = bs->getValue();
  for (Int32 i = 0; i < t_string.size(); ++i) {
    t[i] = 0.0;
    if (t_string[i] != "NULL") {
      applyTraction = true;
      t[i] = std::stod(t_string[i].localstr());
    }
  }

  // no traction to apply hence return
  if (!applyTraction)
    return;

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;

    // 2-point Gauss integration for line element
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
    constexpr Real weights[2] = { 1.0, 1.0 };

    Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, node_coord);

    Node node0 = face.node(0);
    Node node1 = face.node(1);

    for (Int32 i = 0; i < 2; ++i) {
      Real xi = gp[i];
      Real weight = weights[i];

      // Linear shape functions for Line2
      RealVector<2> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsLine2(xi);

      // Integration weight: weight * jacobian (length/2 for reference element [-1,1])
      Real integration_weight = weight * length * 0.5;

      // Apply to both nodes
      Node nodes[2] = { node0, node1 };
      for (Int32 j = 0; j < 2; ++j) {
        Node node = nodes[j];
        if (!node.isOwn())
          continue;

        rhs_values[node_dof.dofId(node, 0)] += t[0] * N[j] * integration_weight;
        rhs_values[node_dof.dofId(node, 1)] += t[1] * N[j] * integration_weight;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies trasient traction conditions to the right-hand side (RHS) values.
 *
 * @param bs The traction boundary condition values.
 * @param t Time at which traction is applied.
 * @param boundary_condition_index Boundary condition index iterator.
 * @param traction_case_table_list Traction case table list.
 * @param node_dof Connectivity view for dof at nodes.
 * @param node_coord Coordinates of the nodes in the mesh.
 * @param rhs_values The right-hand side values to be updated.
 */
void ArcaneFemFunctions::BoundaryConditions2D::
applyTractionTableToRhsTria3(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
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

    const Arcane::FemUtils::CaseTableInfo& case_table_info = traction_case_table_list[boundary_condition_index++];
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
    for (Node node : iface->nodes()) {
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += (trac[0]) * length / 2.;
        rhs_values[node_dof.dofId(node, 1)] += (trac[1]) * length / 2.;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions2D::
applyTractionTableToRhsQuad4(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
                             const UniqueArray<CaseTableInfo>& traction_case_table_list,
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

    const CaseTableInfo& case_table_info = traction_case_table_list[boundary_condition_index++];
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

    // 2-point Gauss integration for line element
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
    constexpr Real weights[2] = { 1.0, 1.0 };

    Real length = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(face, node_coord);

    Node node0 = face.node(0);
    Node node1 = face.node(1);

    for (Int32 i = 0; i < 2; ++i) {
      Real xi = gp[i];
      Real weight = weights[i];

      // Linear shape functions for Line2
      RealVector<2> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsLine2(xi);

      // Integration weight: weight * jacobian (length/2 for reference element [-1,1])
      Real integration_weight = weight * length * 0.5;

      // Apply to both nodes
      Node nodes[2] = { node0, node1 };
      for (Int32 j = 0; j < 2; ++j) {
        Node node = nodes[j];
        if (!node.isOwn())
          continue;

        rhs_values[node_dof.dofId(node, 0)] += trac[0] * N[j] * integration_weight;
        rhs_values[node_dof.dofId(node, 1)] += trac[1] * N[j] * integration_weight;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies Manufactured Dirichlet boundary conditions to RHS and LHS.
 *
 * Updates the LHS matrix and RHS vector to enforce the Dirichlet.
 *
 * - For LHS matrix `𝐀`, the diagonal term for the Dirichlet DOF is set to `𝑃`.
 * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
 *
 * @param manufactured_dirichlet External function for Dirichlet.
 * @param group Group of all external faces.
 * @param bs Boundary condition values.
 * @param node_dof DOF connectivity view.
 * @param node_coord Node coordinates.
 * @param m_linear_system Linear system for LHS.
 * @param rhs_values RHS values to update.
 */
void ArcaneFemFunctions::BoundaryConditions2D::
applyManufacturedDirichletToLhsAndRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_dirichlet, Real /*lambda*/,
                                      const FaceGroup& group, BC::IManufacturedSolution* bs,
                                      const IndexedNodeDoFConnectivityView& node_dof,
                                      const VariableNodeReal3& node_coord,
                                      DoFLinearSystem& m_linear_system,
                                      VariableDoFReal& rhs_values)
{
  Real penalty = bs->getPenalty();
  NodeGroup node_group = group.nodeGroup();
  ENUMERATE_ (Node, inode, node_group) {
    Node node = *inode;
    if (node.isOwn()) {
      m_linear_system.matrixSetValue(node_dof.dofId(node, 0), node_dof.dofId(node, 0), penalty);
      Real u_g = penalty * manufactured_dirichlet->apply(1., node_coord[node]);
      rhs_values[node_dof.dofId(node, 0)] = u_g;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
