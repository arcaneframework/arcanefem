#include "arcane/MathUtils.h"
#include <arcane/utils/NumArray.h>
#include <arcane/IParallelMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/geometry/IGeometry.h>
#include <arcane/VariableTypes.h>
#include "IArcaneFemBC.h"
#include "ArcaneFemFunctions.h"

using namespace Arcane;
/*---------------------------------------------------------------------------*/
/**
 * @brief Initialize the CellFEMDispatcher class (coming from PASSMO)
 */
/*---------------------------------------------------------------------------*/
ArcaneFemFunctions::CellFEMDispatcher::CellFEMDispatcher(){
 // Setting to null default value
 for(int i = 0; i < NB_BASIC_ITEM_TYPE; ++i )
 {
   m_shapefunc[i] = nullptr;
   m_shapefuncderiv[i] = nullptr;
 }

 // Gives functions to compute shape function value in finite-element reference coordinate system
 // Linear elements
 m_shapefunc[IT_Line2] = ArcaneFemFunctions::FemShapeMethods::line2ShapeFuncVal;
 m_shapefunc[IT_Triangle3] = ArcaneFemFunctions::FemShapeMethods::tri3ShapeFuncVal;
 m_shapefunc[IT_Quad4] = ArcaneFemFunctions::FemShapeMethods::quad4ShapeFuncVal;
 m_shapefunc[IT_Tetraedron4] = ArcaneFemFunctions::FemShapeMethods::tetra4ShapeFuncVal;
 m_shapefunc[IT_Hexaedron8] = ArcaneFemFunctions::FemShapeMethods::hexa8ShapeFuncVal;
 m_shapefunc[IT_Pentaedron6] = ArcaneFemFunctions::FemShapeMethods::penta6ShapeFuncVal;
 m_shapefunc[IT_Pyramid5] = ArcaneFemFunctions::FemShapeMethods::pyramid5ShapeFuncVal;

 // Quadratic elements
 m_shapefunc[IT_Line3] = ArcaneFemFunctions::FemShapeMethods::line3ShapeFuncVal;
 m_shapefunc[IT_Triangle6] = ArcaneFemFunctions::FemShapeMethods::tri6ShapeFuncVal;
 m_shapefunc[IT_Quad8] = ArcaneFemFunctions::FemShapeMethods::quad8ShapeFuncVal;
 m_shapefunc[IT_Tetraedron10] = ArcaneFemFunctions::FemShapeMethods::tetra10ShapeFuncVal;
 m_shapefunc[IT_Hexaedron20] = ArcaneFemFunctions::FemShapeMethods::hexa20ShapeFuncVal;

 // Gives functions to compute shape function derivate vector at all nodes of a finite-element
 // along a local direction (in reference coordinate system)
 // Linear elements
 m_shapefuncderiv[IT_Line2] = ArcaneFemFunctions::FemShapeMethods::line2ShapeFuncDeriv;
 m_shapefuncderiv[IT_Triangle3] = ArcaneFemFunctions::FemShapeMethods::tri3ShapeFuncDeriv;
 m_shapefuncderiv[IT_Quad4] = ArcaneFemFunctions::FemShapeMethods::quad4ShapeFuncDeriv;
 m_shapefuncderiv[IT_Tetraedron4] = ArcaneFemFunctions::FemShapeMethods::tetra4ShapeFuncDeriv;
 m_shapefuncderiv[IT_Hexaedron8] = ArcaneFemFunctions::FemShapeMethods::hexa8ShapeFuncDeriv;
 m_shapefuncderiv[IT_Pentaedron6] = ArcaneFemFunctions::FemShapeMethods::penta6ShapeFuncDeriv;
 m_shapefuncderiv[IT_Pyramid5] = ArcaneFemFunctions::FemShapeMethods::pyramid5ShapeFuncDeriv;

 // Quadratic elements
 m_shapefuncderiv[IT_Line3] = ArcaneFemFunctions::FemShapeMethods::line3ShapeFuncDeriv;
 m_shapefuncderiv[IT_Triangle6] = ArcaneFemFunctions::FemShapeMethods::tri6ShapeFuncDeriv;
 m_shapefuncderiv[IT_Quad8] = ArcaneFemFunctions::FemShapeMethods::quad8ShapeFuncDeriv;
 m_shapefuncderiv[IT_Tetraedron10] = ArcaneFemFunctions::FemShapeMethods::tetra10ShapeFuncDeriv;
 m_shapefuncderiv[IT_Hexaedron20] = ArcaneFemFunctions::FemShapeMethods::hexa20ShapeFuncDeriv;

}

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides at once, all Gauss data of a given input finite element
 * (vector containing weights, ref. coordinates, nodal shape values & derivatives)
 * This method is generic (coming from PASSMO)
 */
/*---------------------------------------------------------------------------*/
RealUniqueArray ArcaneFemFunctions::CellFEMDispatcher::
getGaussData(ItemWithNodes item, Integer nint, Integer ngauss){

 auto nnod = item.nbNode();
 auto cell_type = item.type();
 ngauss = ArcaneFemFunctions::FemGaussQuadrature::getNbGaussPointsfromOrder(cell_type,nint);

 // Vector of double containing:
 // ngauss points * [weight, gauss ref coord [Real3], nnod * (shapefunc values, 3*shapefunc deriv
 // in ref. coord system)]
 Integer nsize = ngauss * 4 * (1 + nnod);
 RealUniqueArray vec(nsize);

 Integer index{ 0 };
 for (Integer ig = 0; ig < ngauss; ++ig) {
   auto wt = ArcaneFemFunctions::FemGaussQuadrature::getGaussWeight(item, nint, ig);
   auto pos = ArcaneFemFunctions::FemGaussQuadrature::getGaussRefPosition(item, nint, ig);
   vec[index++] = wt;
   vec[index++] = pos.x;
   vec[index++] = pos.y;
   vec[index++] = pos.z;

   for (Int32 inod = 0; inod < nnod; ++inod) {
     auto Phi_i = getShapeFuncVal(cell_type, inod, pos);
     vec[index++] = Phi_i;
     auto dPhi = getShapeFuncDeriv(cell_type, inod, pos);
     vec[index++] = dPhi.x;
     vec[index++] = dPhi.y;
     vec[index++] = dPhi.z;
   }
 }
 return vec;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the value of the nodal shape function at a given Gauss point
 * (for a given node within in a given element)
 * This method is generic (coming from PASSMO)
 */
/*---------------------------------------------------------------------------*/
Real ArcaneFemFunctions::CellFEMDispatcher::
getShapeFuncVal(Int16 item_type, Integer inod, Real3 coord)
{
 auto f = m_shapefunc[item_type];
 if (f!=nullptr)
   return f(inod,coord);
 return 0.;
}

/*---------------------------------------------------------------------------*/
/**
* @brief Computes the value of the nodal shape function derivatives
 * at a given Gauss point, along 1, 2 and/or 3 directions depending on space dimension
* (for a given node within in a given element)
* This method is generic (coming from PASSMO)
 */
/*---------------------------------------------------------------------------*/
Real3 ArcaneFemFunctions::CellFEMDispatcher::
getShapeFuncDeriv(Int16 item_type, Integer inod, Real3 ref_coord){
 auto f = m_shapefuncderiv[item_type];
 if (f!=nullptr)
   return f(inod,ref_coord);
 return {};
}

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
applyNeumannToRhsLine3(BC::INeumannBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
/**
 * @brief Applies a constant source term to the RHS vector for Quad4 elements.
 *
 * Uses a 2x2 Gauss rule to integrate the biquadratic (serendipity) shape
 * functions exactly. For each Gauss point the shape functions, their
 * derivatives, the Jacobian and its determinant are computed, then the
 * weighted contribution N[i]*qdot*detJ is scattered onto the owned nodes.
 */
/*---------------------------------------------------------------------------*/
void ArcaneFemFunctions::BoundaryConditions2D::
applyConstantSourceToRhsQuad4(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
/**
 * @brief Applies a constant source term to the RHS vector for Quad8 elements.
 *
 * Uses a 3x3 Gauss rule to integrate the biquadratic (serendipity) shape
 * functions exactly. For each Gauss point the shape functions, their
 * derivatives, the Jacobian and its determinant are computed, then the
 * weighted contribution N[i]*qdot*detJ is scattered onto the owned nodes.
 */
/*---------------------------------------------------------------------------*/
void ArcaneFemFunctions::BoundaryConditions2D::
applyConstantSourceToRhsQuad8(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
/**
 * @brief Applies a constant source term to the RHS vector for Hexa8 elements.
 *
 */
/*---------------------------------------------------------------------------*/

void ArcaneFemFunctions::BoundaryConditions3D::
applyConstantSourceToRhsHexa8(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
/**
 * @brief Applies a constant source term to the RHS vector for Hexa20 elements.
 *
 * Uses a 3x3x3 Gauss rule to integrate the quadratic (serendipity) shape
 * functions exactly. Node ordering follows ItemTypeMng.cc (VTK convention):
 *   0-7 corners, 8-11 bottom edges, 12-15 top edges, 16-19 vertical edges.
 */
/*---------------------------------------------------------------------------*/
void ArcaneFemFunctions::BoundaryConditions3D::
applyConstantSourceToRhsHexa20(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
/**
 * @brief Applies a constant source term to the RHS vector for Hexa27 elements.
 *
 * Uses a 3x3x3 Gauss rule to integrate the triquadratic (Lagrange) shape
 * functions exactly. Node ordering follows ItemTypeMng.cc (VTK convention):
 *   0-7 corners, 8-19 edges, 20-25 face centers, 26 body center.
 */
/*---------------------------------------------------------------------------*/
void ArcaneFemFunctions::BoundaryConditions3D::
applyConstantSourceToRhsHexa27(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
