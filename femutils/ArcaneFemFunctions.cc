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
        //   𝐍 = [𝑁₁  𝑁₂  𝑁₃  𝑁₄  𝑁₅  𝑁₆  𝑁₇  𝑁₈]
        //   𝑁₁ = 1/4 * (1-ξ)(1-η)(-ξ-η-1)     𝑁₅ = 1/2 * (1-ξ²)(1-η)
        //   𝑁₂ = 1/4 * (1+ξ)(1-η)( ξ-η-1)     𝑁₆ = 1/2 * (1+ξ)(1-η²)
        //   𝑁₃ = 1/4 * (1+ξ)(1+η)( ξ+η-1)     𝑁₇ = 1/2 * (1-ξ²)(1+η)
        //   𝑁₄ = 1/4 * (1-ξ)(1+η)(-ξ+η-1)     𝑁₈ = 1/2 * (1-ξ)(1-η²)
        Real N[8];
        N[0] = 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1);
        N[1] = 0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1);
        N[2] = 0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1);
        N[3] = 0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1);
        N[4] = 0.5 * (1 - xi * xi) * (1 - eta);
        N[5] = 0.5 * (1 + xi) * (1 - eta * eta);
        N[6] = 0.5 * (1 - xi * xi) * (1 + eta);
        N[7] = 0.5 * (1 - xi) * (1 - eta * eta);

        // Shape function derivatives ∂𝐍/∂ξ and ∂𝐍/∂η
        //     ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ  ∂𝑁₅/∂ξ  ∂𝑁₆/∂ξ  ∂𝑁₇/∂ξ  ∂𝑁₈/∂ξ ]
        //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η  ∂𝑁₅/∂η  ∂𝑁₆/∂η  ∂𝑁₇/∂η  ∂𝑁₈/∂η ]
        Real dN_dxi[8] = {
          0.25 * (1 - eta) * (2 * xi + eta),
          0.25 * (1 - eta) * (2 * xi - eta),
          0.25 * (1 + eta) * (2 * xi + eta),
          0.25 * (1 + eta) * (2 * xi - eta),
          -xi * (1 - eta),
          0.5 * (1 - eta * eta),
          -xi * (1 + eta),
          -0.5 * (1 - eta * eta)
        };

        Real dN_deta[8] = {
          0.25 * (1 - xi) * (2 * eta + xi),
          0.25 * (1 + xi) * (2 * eta - xi),
          0.25 * (1 + xi) * (2 * eta + xi),
          0.25 * (1 - xi) * (2 * eta - xi),
          -0.5 * (1 - xi * xi),
          -eta * (1 + xi),
          0.5 * (1 - xi * xi),
          -eta * (1 - xi)
        };

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
          J00 += dN_dxi[a] * node_coord[cell.nodeId(a)].x;
          J01 += dN_dxi[a] * node_coord[cell.nodeId(a)].y;
          J10 += dN_deta[a] * node_coord[cell.nodeId(a)].x;
          J11 += dN_deta[a] * node_coord[cell.nodeId(a)].y;
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
        //   𝐍 = [𝑁₁  𝑁₂  𝑁₃  𝑁₄  𝑁₅  𝑁₆  𝑁₇  𝑁₈  𝑁₉]
        //   𝑁₁ = 1/4 * ξ(ξ-1)η(η-1)     𝑁₅ = 1/2 * (1-ξ²)η(η-1)
        //   𝑁₂ = 1/4 * ξ(ξ+1)η(η-1)     𝑁₆ = 1/2 * ξ(ξ+1)(1-η²)
        //   𝑁₃ = 1/4 * ξ(ξ+1)η(η+1)     𝑁₇ = 1/2 * (1-ξ²)η(η+1)
        //   𝑁₄ = 1/4 * ξ(ξ-1)η(η+1)     𝑁₈ = 1/2 * ξ(ξ-1)(1-η²)
        //                               𝑁₉ = (1-ξ²)(1-η²)
        Real N[9];
        N[0] = 0.25 * xi * (xi - 1) * eta * (eta - 1);
        N[1] = 0.25 * xi * (xi + 1) * eta * (eta - 1);
        N[2] = 0.25 * xi * (xi + 1) * eta * (eta + 1);
        N[3] = 0.25 * xi * (xi - 1) * eta * (eta + 1);
        N[4] = 0.5 * (1 - xi * xi) * eta * (eta - 1);
        N[5] = 0.5 * xi * (xi + 1) * (1 - eta * eta);
        N[6] = 0.5 * (1 - xi * xi) * eta * (eta + 1);
        N[7] = 0.5 * xi * (xi - 1) * (1 - eta * eta);
        N[8] = (1 - xi * xi) * (1 - eta * eta);

        // Shape function derivatives ∂𝐍/∂ξ and ∂𝐍/∂η
        //     ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ  ∂𝑁₅/∂ξ  ∂𝑁₆/∂ξ  ∂𝑁₇/∂ξ  ∂𝑁₈/∂ξ  ∂𝑁₉/∂ξ ]
        //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η  ∂𝑁₅/∂η  ∂𝑁₆/∂η  ∂𝑁₇/∂η  ∂𝑁₈/∂η  ∂𝑁₉/∂η ]
        Real dN_dxi[9] = {
          0.25 * (2 * xi - 1) * eta * (eta - 1),
          0.25 * (2 * xi + 1) * eta * (eta - 1),
          0.25 * (2 * xi + 1) * eta * (eta + 1),
          0.25 * (2 * xi - 1) * eta * (eta + 1),
          -xi * eta * (eta - 1),
          0.5 * (2 * xi + 1) * (1 - eta * eta),
          -xi * eta * (eta + 1),
          0.5 * (2 * xi - 1) * (1 - eta * eta),
          -2 * xi * (1 - eta * eta)
        };

        Real dN_deta[9] = {
          0.25 * xi * (xi - 1) * (2 * eta - 1),
          0.25 * xi * (xi + 1) * (2 * eta - 1),
          0.25 * xi * (xi + 1) * (2 * eta + 1),
          0.25 * xi * (xi - 1) * (2 * eta + 1),
          0.5 * (1 - xi * xi) * (2 * eta - 1),
          -eta * xi * (xi + 1),
          0.5 * (1 - xi * xi) * (2 * eta + 1),
          -eta * xi * (xi - 1),
          -2 * eta * (1 - xi * xi)
        };

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
          J00 += dN_dxi[a] * node_coord[cell.nodeId(a)].x;
          J01 += dN_dxi[a] * node_coord[cell.nodeId(a)].y;
          J10 += dN_deta[a] * node_coord[cell.nodeId(a)].x;
          J11 += dN_deta[a] * node_coord[cell.nodeId(a)].y;
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