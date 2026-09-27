// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BoundaryConditions3D.cc                                     (C) 2000-2026 */
/*                                                                           */
/* Helper functions for handling 3D boundary conditions.                     */
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
void ArcaneFemFunctions::BoundaryConditions3D::
applyConstantSourceToRhsTetra4(Real qdot, IMesh* mesh,
                               const IndexedNodeDoFConnectivityView& node_dof,
                               const VariableNodeReal3& node_coord,
                               VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;
    Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, node_coord);
    for (Node node : cell.nodes()) {
      if (node.isOwn())
        rhs_values[node_dof.dofId(node, 0)] += qdot * volume / cell.nbNode();
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector for Hexa8 elements.
 */
void ArcaneFemFunctions::BoundaryConditions3D::
applyConstantSourceToRhsHexa8(Real qdot, IMesh* mesh,
                              const IndexedNodeDoFConnectivityView& node_dof,
                              const VariableNodeReal3& node_coord,
                              VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;

    // Gauss quadrature for Hexa8
    // Using 2x2x2 Gauss points for integration
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // {-1/sqrt(3) 1/sqrt(3)}
    constexpr Real weights[2] = { 1.0, 1.0 };

    for (Int32 ixi = 0; ixi < 2; ++ixi) {
      for (Int32 ieta = 0; ieta < 2; ++ieta) {
        for (Int32 izeta = 0; izeta < 2; ++izeta) {

          // Gauss point coordinates in reference space
          Real xi = gp[ixi]; // ξ coordinate
          Real eta = gp[ieta]; // η coordinate
          Real zeta = gp[izeta]; // ζ coordinate
          Real weight = weights[ixi] * weights[ieta] * weights[izeta];

          // Shape functions 𝐍 for Hexa8
          RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa8(xi, eta, zeta);

          // Shape function derivatives in reference space
          //  ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ  ∂𝑁₅/∂ξ  ∂𝑁₆/∂ξ  ∂𝑁₇/∂ξ  ∂𝑁₈/∂ξ ]
          //  ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η  ∂𝑁₅/∂η  ∂𝑁₆/∂η  ∂𝑁₇/∂η  ∂𝑁₈/∂η ]
          //  ∂𝐍/∂ζ = [ ∂𝑁₁/∂ζ  ∂𝑁₂/∂ζ  ∂𝑁₃/∂ζ  ∂𝑁₄/∂ζ  ∂𝑁₅/∂ζ  ∂𝑁₆/∂ζ  ∂𝑁₇/∂ζ  ∂𝑁₈/∂ζ ]
          const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsHexa8(xi, eta, zeta);
          // Jacobian for 3D (using your working stiffness matrix approach)
          Real3x3 J;
          for (Int8 a = 0; a < 8; ++a) {
            const Real3& n = node_coord[cell.nodeId(a)];
            J[0][0] += reference_gradients.dN_dxi[a] * n.x; // ∂𝑥/∂ξ
            J[0][1] += reference_gradients.dN_dxi[a] * n.y; // ∂𝑦/∂ξ
            J[0][2] += reference_gradients.dN_dxi[a] * n.z; // ∂𝑧/∂ξ
            J[1][0] += reference_gradients.dN_deta[a] * n.x; // ∂𝑥/∂η
            J[1][1] += reference_gradients.dN_deta[a] * n.y; // ∂𝑦/∂η
            J[1][2] += reference_gradients.dN_deta[a] * n.z; // ∂𝑧/∂η
            J[2][0] += reference_gradients.dN_dzeta[a] * n.x; // ∂𝑥/∂ζ
            J[2][1] += reference_gradients.dN_dzeta[a] * n.y; // ∂𝑦/∂ζ
            J[2][2] += reference_gradients.dN_dzeta[a] * n.z; // ∂𝑧/∂ζ
          }

          // Compute determinant of Jacobian
          Real detJ = math::matrixDeterminant(J);

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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector for Hexa20 elements.
 *
 * Uses a 3x3x3 Gauss rule to integrate the quadratic (serendipity) shape
 * functions exactly. Node ordering follows ItemTypeMng.cc (VTK convention):
 *   0-7 corners, 8-11 bottom edges, 12-15 top edges, 16-19 vertical edges.
 */
void ArcaneFemFunctions::BoundaryConditions3D::
applyConstantSourceToRhsHexa20(Real qdot, IMesh* mesh,
                               const IndexedNodeDoFConnectivityView& node_dof,
                               const VariableNodeReal3& node_coord,
                               VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;

    // 3-point Gauss rule per direction (needed for exact integration of quadratic Hexa20 shape functions)
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(3/5) , 0 , sqrt(3/5)]
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {
        for (Int32 izeta = 0; izeta < 3; ++izeta) {

          // Gauss point coordinates in reference space
          Real xi = gp[ixi]; // ξ coordinate
          Real eta = gp[ieta]; // η coordinate
          Real zeta = gp[izeta]; // ζ coordinate
          Real weight = weights[ixi] * weights[ieta] * weights[izeta];

          // Shape functions 𝐍 for Hexa20 (serendipity)
          RealVector<20> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa20(xi, eta, zeta);

          // Shape function derivatives ∂𝐍/∂ξ, ∂𝐍/∂η, ∂𝐍/∂ζ
          const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsHexa20(xi, eta, zeta);
          // Jacobian matrix (default-initialized to zero see Real3x3.h)
          Real3x3 J;
          for (Int8 a = 0; a < 20; ++a) {
            const Real3& n_coord = node_coord[cell.nodeId(a)];
            J[0][0] += reference_gradients.dN_dxi[a] * n_coord.x;
            J[0][1] += reference_gradients.dN_dxi[a] * n_coord.y;
            J[0][2] += reference_gradients.dN_dxi[a] * n_coord.z;
            J[1][0] += reference_gradients.dN_deta[a] * n_coord.x;
            J[1][1] += reference_gradients.dN_deta[a] * n_coord.y;
            J[1][2] += reference_gradients.dN_deta[a] * n_coord.z;
            J[2][0] += reference_gradients.dN_dzeta[a] * n_coord.x;
            J[2][1] += reference_gradients.dN_dzeta[a] * n_coord.y;
            J[2][2] += reference_gradients.dN_dzeta[a] * n_coord.z;
          }

          // Determinant of the Jacobian
          Real detJ = math::matrixDeterminant(J);
          if (detJ <= 0.0) {
            ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
          }

          // Compute integration weight
          Real integration_weight = weight * detJ;

          // Assemble RHS
          for (Int32 i = 0; i < 20; ++i) {
            Node node = cell.node(i);
            if (node.isOwn()) {
              rhs_values[node_dof.dofId(node, 0)] += N[i] * qdot * integration_weight;
            }
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a constant source term to the RHS vector for Hexa27 elements.
 *
 * Uses a 3x3x3 Gauss rule to integrate the triquadratic (Lagrange) shape
 * functions exactly. Node ordering follows ItemTypeMng.cc (VTK convention):
 *   0-7 corners, 8-19 edges, 20-25 face centers, 26 body center.
 */
void ArcaneFemFunctions::BoundaryConditions3D::
applyConstantSourceToRhsHexa27(Real qdot, IMesh* mesh,
                               const IndexedNodeDoFConnectivityView& node_dof,
                               const VariableNodeReal3& node_coord,
                               VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;

    // 3-point Gauss rule per direction (needed for exact integration of quadratic Hexa27 shape functions)
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(3/5) , 0 , sqrt(3/5)]
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {
        for (Int32 izeta = 0; izeta < 3; ++izeta) {

          // Gauss point coordinates in reference space
          Real xi = gp[ixi]; // ξ coordinate
          Real eta = gp[ieta]; // η coordinate
          Real zeta = gp[izeta]; // ζ coordinate
          Real weight = weights[ixi] * weights[ieta] * weights[izeta];

          // Shape functions 𝐍 for Hexa27 (triquadratic Lagrange)
          RealVector<27> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa27(xi, eta, zeta);

          // Shape function derivatives ∂𝐍/∂ξ, ∂𝐍/∂η, ∂𝐍/∂ζ
          const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsHexa27(xi, eta, zeta);
          // Jacobian matrix (default-initialized to zero see Real3x3.h)
          Real3x3 J;
          for (Int8 a = 0; a < 27; ++a) {
            const Real3& n_coord = node_coord[cell.nodeId(a)];
            J[0][0] += reference_gradients.dN_dxi[a] * n_coord.x;
            J[0][1] += reference_gradients.dN_dxi[a] * n_coord.y;
            J[0][2] += reference_gradients.dN_dxi[a] * n_coord.z;
            J[1][0] += reference_gradients.dN_deta[a] * n_coord.x;
            J[1][1] += reference_gradients.dN_deta[a] * n_coord.y;
            J[1][2] += reference_gradients.dN_deta[a] * n_coord.z;
            J[2][0] += reference_gradients.dN_dzeta[a] * n_coord.x;
            J[2][1] += reference_gradients.dN_dzeta[a] * n_coord.y;
            J[2][2] += reference_gradients.dN_dzeta[a] * n_coord.z;
          }

          // Determinant of the Jacobian
          Real detJ = math::matrixDeterminant(J);
          if (detJ <= 0.0) {
            ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
          }

          // Compute integration weight
          Real integration_weight = weight * detJ;

          // Assemble RHS
          for (Int32 i = 0; i < 27; ++i) {
            Node node = cell.node(i);
            if (node.isOwn()) {
              rhs_values[node_dof.dofId(node, 0)] += N[i] * qdot * integration_weight;
            }
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies a nodal field as a source term to the RHS vector.
 *
 * @param field The field values at cell nodes.
 * @param mesh The mesh containing all cells.
 * @param node_dof DOF connectivity view.
 * @param node_coord The coordinates of the nodes.
 * @param rhs_values The RHS values to update.
 */
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions3D::
integrateNodalFieldToRhsTetra4(VariableNodeReal& field, IMesh* mesh,
                               const IndexedNodeDoFConnectivityView& node_dof,
                               const VariableNodeReal3& node_coord,
                               VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;
    Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, node_coord);

    // Get nodal values for this tetrahedral cell
    const Real field_at_nodes[4] = {
      field[cell.nodeId(0)],
      field[cell.nodeId(1)],
      field[cell.nodeId(2)],
      field[cell.nodeId(3)]
    };

    // Apply mass matrix integration: ∫(T^n * N_i)dV
    Real node_contributions[4] = { 0.0, 0.0, 0.0, 0.0 };

    for (Int8 i = 0; i < 4; ++i) {
      for (Int8 j = 0; j < 4; ++j) {
        Real mass_coeff;
        if (i == j) {
          mass_coeff = volume / 10.0; // diagonal: volume * (2/20) = volume/10
        }
        else {
          mass_coeff = volume / 20.0; // off-diagonal: volume * (1/20)
        }
        node_contributions[i] += mass_coeff * field_at_nodes[j];
      }
    }

    // Add contributions to global RHS
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
 * @brief Applies a nodal field to the RHS vector.
 *
 * @param field The field term defined on nodes.
 * @param mesh The mesh containing all cells.
 * @param node_dof DOF connectivity view.
 * @param node_coord The coordinates of the nodes.
 * @param rhs_values The RHS values to update.
 */
void ArcaneFemFunctions::BoundaryConditions3D::
integrateNodalFieldToRhsHexa8(VariableNodeReal& field, IMesh* mesh,
                              const IndexedNodeDoFConnectivityView& node_dof,
                              const VariableNodeReal3& node_coord,
                              VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;

    // Get nodal values of field for this cell
    const Real field_at_nodes[8] = {
      field[cell.nodeId(0)],
      field[cell.nodeId(1)],
      field[cell.nodeId(2)],
      field[cell.nodeId(3)],
      field[cell.nodeId(4)],
      field[cell.nodeId(5)],
      field[cell.nodeId(6)],
      field[cell.nodeId(7)]
    };

    // Initialize contributions for each node in this cell
    Real node_contributions[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    // 2x2 Gauss integration for quadrilateral element
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
    constexpr Real w = 1.0;

    for (Int32 ixi = 0; ixi < 2; ++ixi) {
      for (Int32 ieta = 0; ieta < 2; ++ieta) {
        for (Int32 izeta = 0; izeta < 2; ++izeta) {

          // Get the coordinates of the Gauss point
          Real xi = gp[ixi]; // Get the ξ coordinate of the Gauss point
          Real eta = gp[ieta]; // Get the η coordinate of the Gauss point
          Real zeta = gp[izeta]; // ζ coordinate
          Real weight = w * w * w; // Weight for 3D Gauss integration

          // Shape functions 𝐍 for Hexa8
          RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa8(xi, eta, zeta);

          // compute the det(Jacobian)
          const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa8(cell, node_coord, xi, eta, zeta);
          const Real detJ = gp_info.det_j;

          // compute integration weight
          const Real integration_weight = weight * detJ;

          // Interpolate qdot at the quadrature point: qdot_gp = ∑ 𝑁ᵢ * q̇
          Real qdot_gp = 0.0;
          for (Int8 a = 0; a < 8; ++a) {
            qdot_gp += N[a] * field_at_nodes[a];
          }

          // Add contribution to each test function: ∫ q̇ * 𝑁ᵢ dΩ
          for (Int8 i = 0; i < 8; ++i) {
            node_contributions[i] += qdot_gp * N[i] * integration_weight;
          }
        }
      }
    }

    // Add contributions to global RHS vector
    for (Int8 i = 0; i < 8; ++i) {
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
void ArcaneFemFunctions::BoundaryConditions3D::
applyManufacturedSourceToRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_source,
                             IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof,
                             const VariableNodeReal3& node_coord,
                             VariableDoFReal& rhs_values)
{
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    Cell cell = *icell;
    Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, node_coord);
    Real3 bcenter = ArcaneFemFunctions::MeshOperation::computeBaryCenterTetra4(cell, node_coord);

    for (Node node : cell.nodes()) {
      if (node.isOwn())
        rhs_values[node_dof.dofId(node, 0)] += manufactured_source->apply(volume / cell.nbNode(), bcenter);
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
void ArcaneFemFunctions::BoundaryConditions3D::
applyNeumannToRhsTetra4(BC::INeumannBoundaryCondition* bs,
                        const IndexedNodeDoFConnectivityView& node_dof,
                        const VariableNodeReal3& node_coord,
                        VariableDoFReal& rhs_values)
{
  FaceGroup group = bs->getSurface();

  Real value = 0.0;
  Real valueX = 0.0;
  Real valueY = 0.0;
  Real valueZ = 0.0;

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
      if (neumann_str[2] != "NULL")
        valueZ = std::stod(neumann_str[2].localstr());
    }
  }

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;

    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, node_coord);
    Real3 normal = ArcaneFemFunctions::MeshOperation::computeNormalTriangle(face, node_coord);

    for (Node node : iface->nodes()) {
      if (!node.isOwn())
        continue;
      Real rhs_value;

      if (scalarNeumann) {
        rhs_value = value * area / 3.0;
      }
      else {
        rhs_value = (normal.x * valueX + normal.y * valueY + normal.z * valueZ) * area / 3.0;
      }

      rhs_values[node_dof.dofId(node, 0)] += rhs_value;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions3D::
applyNeumannToRhsHexa8(BC::INeumannBoundaryCondition* bs,
                       const IndexedNodeDoFConnectivityView& node_dof,
                       const VariableNodeReal3& node_coord,
                       VariableDoFReal& rhs_values)
{
  FaceGroup group = bs->getSurface();

  Real value = 0.0;
  Real valueX = 0.0;
  Real valueY = 0.0;
  Real valueZ = 0.0;

  bool scalarNeumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalarNeumann = true;
    value = std::stod(neumann_str[0].localstr());
  }
  else {
    if (neumann_str.size() > 2) {
      if (neumann_str[0] != "NULL")
        valueX = std::stod(neumann_str[0].localstr());
      if (neumann_str[1] != "NULL")
        valueY = std::stod(neumann_str[1].localstr());
      if (neumann_str[2] != "NULL")
        valueZ = std::stod(neumann_str[2].localstr());
    }
  }

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;

    // 2x2 Gauss integration for quadrilateral face
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
    constexpr Real w = 1.0;

    // Get face nodes (assuming quad4 face)
    Node node0 = face.node(0);
    Node node1 = face.node(1);
    Node node2 = face.node(2);
    Node node3 = face.node(3);
    Node nodes[4] = { node0, node1, node2, node3 };

    // Get node coordinates
    Real3 coords[4];
    for (Int32 i = 0; i < 4; ++i) {
      coords[i] = node_coord[nodes[i]];
    }

    // Loop through 2x2 Gauss points
    for (Int32 ixi = 0; ixi < 2; ++ixi) {
      for (Int32 ieta = 0; ieta < 2; ++ieta) {
        Real xi = gp[ixi];
        Real eta = gp[ieta];

        // Quad4 shape functions
        RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);

        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad4(xi, eta);

        // Compute tangent vectors
        Real3 t1(0.0, 0.0, 0.0); // ∂r/∂ξ
        Real3 t2(0.0, 0.0, 0.0); // ∂r/∂η

        for (Int32 i = 0; i < 4; ++i) {
          t1.x += reference_gradients.dN_dxi[i] * coords[i].x;
          t1.y += reference_gradients.dN_dxi[i] * coords[i].y;
          t1.z += reference_gradients.dN_dxi[i] * coords[i].z;

          t2.x += reference_gradients.dN_deta[i] * coords[i].x;
          t2.y += reference_gradients.dN_deta[i] * coords[i].y;
          t2.z += reference_gradients.dN_deta[i] * coords[i].z;
        }

        // Normal vector (cross product of tangent vectors)
        Real3 normal;
        normal.x = t1.y * t2.z - t1.z * t2.y;
        normal.y = t1.z * t2.x - t1.x * t2.z;
        normal.z = t1.x * t2.y - t1.y * t2.x;

        // Jacobian (magnitude of normal vector for surface integration)
        Real detJ = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        if (detJ <= 0.0) {
          ARCANE_FATAL("Invalid (non-positive) surface Jacobian: {0}", detJ);
        }

        // Unit normal
        normal.x /= detJ;
        normal.y /= detJ;
        normal.z /= detJ;

        // Integration weight
        Real integration_weight = w * w * detJ;

        // Apply to all four nodes of the face
        for (Int32 j = 0; j < 4; ++j) {
          Node node = nodes[j];
          if (!node.isOwn())
            continue;

          Real rhs_value;
          if (scalarNeumann) {
            rhs_value = value * N[j] * integration_weight;
          }
          else {
            rhs_value = (normal.x * valueX + normal.y * valueY + normal.z * valueZ) * N[j] * integration_weight;
          }

          rhs_values[node_dof.dofId(node, 0)] += rhs_value;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions3D::
applyNeumannToRhsHexa20(BC::INeumannBoundaryCondition* bs,
                        const IndexedNodeDoFConnectivityView& node_dof,
                        const VariableNodeReal3& node_coord,
                        VariableDoFReal& rhs_values)
{
  FaceGroup group = bs->getSurface();

  Real value = 0.0;
  Real valueX = 0.0;
  Real valueY = 0.0;
  Real valueZ = 0.0;

  bool scalarNeumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalarNeumann = true;
    value = std::stod(neumann_str[0].localstr());
  }
  else {
    if (neumann_str.size() > 2) {
      if (neumann_str[0] != "NULL")
        valueX = std::stod(neumann_str[0].localstr());
      if (neumann_str[1] != "NULL")
        valueY = std::stod(neumann_str[1].localstr());
      if (neumann_str[2] != "NULL")
        valueZ = std::stod(neumann_str[2].localstr());
    }
  }

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;

    // 3-point Gauss rule per direction (needed for quadratic Quad8 face)
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(3/5) , 0 , sqrt(3/5)]
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    // Get face nodes (Quad8 face)
    Node nodes[8];
    for (Int32 i = 0; i < 8; ++i)
      nodes[i] = face.node(i);

    // Get node coordinates
    Real3 coords[8];
    for (Int32 i = 0; i < 8; ++i)
      coords[i] = node_coord[nodes[i]];

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {
        Real xi = gp[ixi];
        Real eta = gp[ieta];
        Real weight = weights[ixi] * weights[ieta];

        // Quad8 (serendipity) shape functions
        RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad8(xi, eta);

        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad8(xi, eta);

        // Tangent vectors ∂r/∂ξ and ∂r/∂η
        Real3 t1(0.0, 0.0, 0.0);
        Real3 t2(0.0, 0.0, 0.0);
        for (Int32 i = 0; i < 8; ++i) {
          t1.x += reference_gradients.dN_dxi[i] * coords[i].x;
          t1.y += reference_gradients.dN_dxi[i] * coords[i].y;
          t1.z += reference_gradients.dN_dxi[i] * coords[i].z;
          t2.x += reference_gradients.dN_deta[i] * coords[i].x;
          t2.y += reference_gradients.dN_deta[i] * coords[i].y;
          t2.z += reference_gradients.dN_deta[i] * coords[i].z;
        }

        // Normal vector (cross product of tangent vectors)
        Real3 normal;
        normal.x = t1.y * t2.z - t1.z * t2.y;
        normal.y = t1.z * t2.x - t1.x * t2.z;
        normal.z = t1.x * t2.y - t1.y * t2.x;

        // Surface Jacobian
        Real detJ = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        if (detJ <= 0.0) {
          ARCANE_FATAL("Invalid (non-positive) surface Jacobian: {0}", detJ);
        }

        // Unit normal
        normal.x /= detJ;
        normal.y /= detJ;
        normal.z /= detJ;

        // Integration weight
        Real integration_weight = weight * detJ;

        // Apply to the eight nodes of the face
        for (Int32 j = 0; j < 8; ++j) {
          Node node = nodes[j];
          if (!node.isOwn())
            continue;

          Real rhs_value;
          if (scalarNeumann) {
            rhs_value = value * N[j] * integration_weight;
          }
          else {
            rhs_value = (normal.x * valueX + normal.y * valueY + normal.z * valueZ) * N[j] * integration_weight;
          }

          rhs_values[node_dof.dofId(node, 0)] += rhs_value;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions3D::
applyNeumannToRhsHexa27(BC::INeumannBoundaryCondition* bs,
                        const IndexedNodeDoFConnectivityView& node_dof,
                        const VariableNodeReal3& node_coord,
                        VariableDoFReal& rhs_values)
{
  FaceGroup group = bs->getSurface();

  Real value = 0.0;
  Real valueX = 0.0;
  Real valueY = 0.0;
  Real valueZ = 0.0;

  bool scalarNeumann = false;
  const StringConstArrayView neumann_str = bs->getValue();

  if (neumann_str.size() == 1 && neumann_str[0] != "NULL") {
    scalarNeumann = true;
    value = std::stod(neumann_str[0].localstr());
  }
  else {
    if (neumann_str.size() > 2) {
      if (neumann_str[0] != "NULL")
        valueX = std::stod(neumann_str[0].localstr());
      if (neumann_str[1] != "NULL")
        valueY = std::stod(neumann_str[1].localstr());
      if (neumann_str[2] != "NULL")
        valueZ = std::stod(neumann_str[2].localstr());
    }
  }

  ENUMERATE_ (Face, iface, group) {
    Face face = *iface;

    // 3-point Gauss rule per direction (needed for quadratic Quad9 face)
    constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 }; // [-sqrt(3/5) , 0 , sqrt(3/5)]
    constexpr Real weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    // Get face nodes (Quad9 face)
    Node nodes[9];
    for (Int32 i = 0; i < 9; ++i)
      nodes[i] = face.node(i);

    // Get node coordinates
    Real3 coords[9];
    for (Int32 i = 0; i < 9; ++i)
      coords[i] = node_coord[nodes[i]];

    for (Int32 ixi = 0; ixi < 3; ++ixi) {
      for (Int32 ieta = 0; ieta < 3; ++ieta) {
        Real xi = gp[ixi];
        Real eta = gp[ieta];
        Real weight = weights[ixi] * weights[ieta];

        // Quad9 (Lagrange) shape functions
        RealVector<9> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad9(xi, eta);

        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad9(xi, eta);

        // Tangent vectors ∂r/∂ξ and ∂r/∂η
        Real3 t1(0.0, 0.0, 0.0);
        Real3 t2(0.0, 0.0, 0.0);
        for (Int32 i = 0; i < 9; ++i) {
          t1.x += reference_gradients.dN_dxi[i] * coords[i].x;
          t1.y += reference_gradients.dN_dxi[i] * coords[i].y;
          t1.z += reference_gradients.dN_dxi[i] * coords[i].z;
          t2.x += reference_gradients.dN_deta[i] * coords[i].x;
          t2.y += reference_gradients.dN_deta[i] * coords[i].y;
          t2.z += reference_gradients.dN_deta[i] * coords[i].z;
        }

        // Normal vector (cross product of tangent vectors)
        Real3 normal;
        normal.x = t1.y * t2.z - t1.z * t2.y;
        normal.y = t1.z * t2.x - t1.x * t2.z;
        normal.z = t1.x * t2.y - t1.y * t2.x;

        // Surface Jacobian
        Real detJ = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        if (detJ <= 0.0) {
          ARCANE_FATAL("Invalid (non-positive) surface Jacobian: {0}", detJ);
        }

        // Unit normal
        normal.x /= detJ;
        normal.y /= detJ;
        normal.z /= detJ;

        // Integration weight
        Real integration_weight = weight * detJ;

        // Apply to the nine nodes of the face
        for (Int32 j = 0; j < 9; ++j) {
          Node node = nodes[j];
          if (!node.isOwn())
            continue;

          Real rhs_value;
          if (scalarNeumann) {
            rhs_value = value * N[j] * integration_weight;
          }
          else {
            rhs_value = (normal.x * valueX + normal.y * valueY + normal.z * valueZ) * N[j] * integration_weight;
          }

          rhs_values[node_dof.dofId(node, 0)] += rhs_value;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Applies traction conditions to the right-hand side (RHS) values.
 *
 * This method updates the RHS values of the finite element method equations
 * based on the provided traction boundary condition. The boundary condition
 * can specify a value or its components along the x and y directions.
 *
 * @param bs The traction boundary condition values.
 * @param node_dof Connectivity view for degrees of freedom at nodes.
 * @param node_coord Coordinates of the nodes in the mesh.
 * @param rhs_values The right-hand side values to be updated.
 */
void ArcaneFemFunctions::BoundaryConditions3D::
applyTractionToRhsTetra4(BC::ITractionBoundaryCondition* bs,
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
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, node_coord);
    for (Node node : iface->nodes()) {
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += t[0] * area / 3.;
        rhs_values[node_dof.dofId(node, 1)] += t[1] * area / 3.;
        rhs_values[node_dof.dofId(node, 2)] += t[2] * area / 3.;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions3D::
applyTractionToRhsHexa8(BC::ITractionBoundaryCondition* bs,
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

    // 2x2 Gauss integration for quadrilateral face
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
    constexpr Real w = 1.0;

    // Get face nodes (assuming quad4 face)
    Node node0 = face.node(0);
    Node node1 = face.node(1);
    Node node2 = face.node(2);
    Node node3 = face.node(3);
    Node nodes[4] = { node0, node1, node2, node3 };

    // Get node coordinates
    Real3 coords[4];
    for (Int32 i = 0; i < 4; ++i) {
      coords[i] = node_coord[nodes[i]];
    }

    // Loop through 2x2 Gauss points
    for (Int32 ixi = 0; ixi < 2; ++ixi) {
      for (Int32 ieta = 0; ieta < 2; ++ieta) {
        Real xi = gp[ixi];
        Real eta = gp[ieta];

        // Quad4 shape functions
        RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);

        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad4(xi, eta);

        // Compute tangent vectors
        Real3 t1(0.0, 0.0, 0.0); // ∂r/∂ξ
        Real3 t2(0.0, 0.0, 0.0); // ∂r/∂η

        for (Int32 i = 0; i < 4; ++i) {
          t1.x += reference_gradients.dN_dxi[i] * coords[i].x;
          t1.y += reference_gradients.dN_dxi[i] * coords[i].y;
          t1.z += reference_gradients.dN_dxi[i] * coords[i].z;

          t2.x += reference_gradients.dN_deta[i] * coords[i].x;
          t2.y += reference_gradients.dN_deta[i] * coords[i].y;
          t2.z += reference_gradients.dN_deta[i] * coords[i].z;
        }

        // Normal vector (cross product of tangent vectors)
        Real3 normal;
        normal.x = t1.y * t2.z - t1.z * t2.y;
        normal.y = t1.z * t2.x - t1.x * t2.z;
        normal.z = t1.x * t2.y - t1.y * t2.x;

        // Jacobian (magnitude of normal vector for surface integration)
        Real detJ = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        // Integration weight
        Real integration_weight = w * w * detJ;

        // Apply to all four nodes of the face
        for (Int32 j = 0; j < 4; ++j) {
          Node node = nodes[j];
          if (!node.isOwn())
            continue;

          rhs_values[node_dof.dofId(node, 0)] += t[0] * N[j] * integration_weight;
          rhs_values[node_dof.dofId(node, 1)] += t[1] * N[j] * integration_weight;
          rhs_values[node_dof.dofId(node, 2)] += t[2] * N[j] * integration_weight;
        }
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
void ArcaneFemFunctions::BoundaryConditions3D::
applyTractionTableToRhsTetra4(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
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
    Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(face, node_coord);
    for (Node node : iface->nodes()) {
      if (node.isOwn()) {
        rhs_values[node_dof.dofId(node, 0)] += trac[0] * area / 3.;
        rhs_values[node_dof.dofId(node, 1)] += trac[1] * area / 3.;
        rhs_values[node_dof.dofId(node, 2)] += trac[2] * area / 3.;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions3D::
applyTractionTableToRhsHexa8(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index,
                             const ConstArrayView<CaseTableInfo>& traction_case_table_list,
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

    // 2x2 Gauss integration for quadrilateral face
    constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 }; // -1/sqrt(3), 1/sqrt(3)
    constexpr Real w = 1.0;

    // Get face nodes (assuming quad4 face)
    Node node0 = face.node(0);
    Node node1 = face.node(1);
    Node node2 = face.node(2);
    Node node3 = face.node(3);
    Node nodes[4] = { node0, node1, node2, node3 };

    // Get node coordinates
    Real3 coords[4];
    for (Int32 i = 0; i < 4; ++i) {
      coords[i] = node_coord[nodes[i]];
    }

    // Loop through 2x2 Gauss points
    for (Int32 ixi = 0; ixi < 2; ++ixi) {
      for (Int32 ieta = 0; ieta < 2; ++ieta) {
        Real xi = gp[ixi];
        Real eta = gp[ieta];

        // Quad4 shape functions
        RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);

        const auto reference_gradients = Arcane::FemUtils::ShapeFunctions::computeReferenceGradientsQuad4(xi, eta);

        // Compute tangent vectors
        Real3 t1(0.0, 0.0, 0.0); // ∂r/∂ξ
        Real3 t2(0.0, 0.0, 0.0); // ∂r/∂η

        for (Int32 i = 0; i < 4; ++i) {
          t1.x += reference_gradients.dN_dxi[i] * coords[i].x;
          t1.y += reference_gradients.dN_dxi[i] * coords[i].y;
          t1.z += reference_gradients.dN_dxi[i] * coords[i].z;

          t2.x += reference_gradients.dN_deta[i] * coords[i].x;
          t2.y += reference_gradients.dN_deta[i] * coords[i].y;
          t2.z += reference_gradients.dN_deta[i] * coords[i].z;
        }

        // Normal vector (cross product of tangent vectors)
        Real3 normal;
        normal.x = t1.y * t2.z - t1.z * t2.y;
        normal.y = t1.z * t2.x - t1.x * t2.z;
        normal.z = t1.x * t2.y - t1.y * t2.x;

        // Jacobian (magnitude of normal vector for surface integration)
        Real detJ = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        // Integration weight
        Real integration_weight = w * w * detJ;

        // Apply to all four nodes of the face
        for (Int32 j = 0; j < 4; ++j) {
          Node node = nodes[j];
          if (!node.isOwn())
            continue;

          rhs_values[node_dof.dofId(node, 0)] += trac[0] * N[j] * integration_weight;
          rhs_values[node_dof.dofId(node, 1)] += trac[1] * N[j] * integration_weight;
          rhs_values[node_dof.dofId(node, 2)] += trac[2] * N[j] * integration_weight;
        }
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
void ArcaneFemFunctions::BoundaryConditions3D::
applyManufacturedDirichletToLhsAndRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_dirichlet,
                                      Real /*lambda*/,
                                      const FaceGroup& group,
                                      BC::IManufacturedSolution* bs,
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
