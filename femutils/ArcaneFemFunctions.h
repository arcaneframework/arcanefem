// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArcaneFemFunctions.h                                        (C) 2022-2025 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef ARCANE_FEM_FUNCTIONS_H
#define ARCANE_FEM_FUNCTIONS_H

# define M_SQRT1_3	0.57735026918962576451	/* 1/sqrt(3) */

#include <arcane/core/IStandardFunction.h>

#include <arcane/utils/ITraceMng.h>
#include <arcane/utils/StringList.h>
#include <arcane/utils/CommandLineArguments.h>

#include <arcane/IndexedItemConnectivityView.h>
#include <arcane/VariableTypes.h>
#include <arcane/IMesh.h>

#include "IArcaneFemBC.h"
#include "GaussQuadrature.h"
#include "DoFLinearSystem.h"

using namespace Arcane;
using namespace Arcane::FemUtils;
/*---------------------------------------------------------------------------*/
/**
 * @brief Contains various functions & operations related to FEM calculations.
 *
 * The class provides methods organized into different nested classes for:
 * - MeshOperation: Mesh related operations.
 * - FeOperation2D/3D: Finite element operations at element level.
 * - BoundaryConditions2D/3D: Boundary condition related operations.
 */
/*---------------------------------------------------------------------------*/

class ArcaneFemFunctions
{
 public:
  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides general purpose fuctions
   */
  /*---------------------------------------------------------------------------*/
  class GeneralFunctions
  {
    public:

     /*---------------------------------------------------------------------------*/
     /**
      * @brief Logs the execution time with ArcaneFem-Timer.
      *
      * @param tm The Arcane trace manager for logging.
      * @param label A short description of the event being timed.
      * @param value The elapsed time associated with the event.
      */
     /*---------------------------------------------------------------------------*/
     static inline void printArcaneFemTime(ITraceMng* tm, const String& label, const Real& value)
     {
       ARCANE_CHECK_POINTER(tm);
       tm->info() << std::left << std::setw(40) << "[ArcaneFem-Timer] " + label << " = " << value;
     }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Parses PETSc command-line flags into a CommandLineArguments object.
     *
     * This method processes a space-separated string of PETSc flags and converts
     * it into a list of individual arguments. We do the following:
     *   1. Convert the input String to a standard C++ string.
     *   2. Use a string stream to tokenize the input based on spaces.
     *   3. Store the tokens in a StringList.
     *   4. Return a CommandLineArguments object constructed from the StringList.
     *
     * @param `petsc_flags` space-separated string containing PETSc CLI flags.
     * @return A CommandLineArguments object containing the parsed flags.
     */
    /*---------------------------------------------------------------------------*/

     static inline CommandLineArguments getPetscFlagsFromCommandline(const String& petsc_flags)
     {
       StringList string_list;
       std::string petsc_flags_std = petsc_flags.localstr();

       // Use a string stream to split the string by spaces
       std::istringstream iss(petsc_flags_std);
       String token;
       while (iss >> token) {
         string_list.add(token);
       }
       return CommandLineArguments(string_list);
     }
  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for various mesh-related operations.
   *
   * This class includes static methods for computing geometric properties
   * of mesh elements, such as the area of triangles, the length of edges,
   * and the normal vectors of edges.
   */
  /*---------------------------------------------------------------------------*/
  class MeshOperation
  {

   public:

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the volume of a tetrahedra defined by four nodes.
     *
     * This method calculates the volume using the scalar triple product formula.
     * We do the following:
     *   1. get the four nodes
     *   2. ge the vector representing the edges of tetrahedron
     *   3. compute volume using scalar triple product
     */
    /*---------------------------------------------------------------------------*/

    static inline Real computeVolumeTetra4(ItemWithNodes item, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[item.nodeId(0)];
      Real3 vertex1 = node_coord[item.nodeId(1)];
      Real3 vertex2 = node_coord[item.nodeId(2)];
      Real3 vertex3 = node_coord[item.nodeId(3)];

      Real3 v0 = vertex1 - vertex0;
      Real3 v1 = vertex2 - vertex0;
      Real3 v2 = vertex3 - vertex0;

      return math::abs(math::dot(v0, math::cross(v1, v2))) / 6.0;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the area of a triangle defined by three nodes.
     *
     * This method calculates the area using the determinant formula for a triangle.
     * The area is computed as half the value of the determinant of the matrix
     * formed by the coordinates of the triangle's vertices.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real computeAreaTria3(ItemWithNodes item, const VariableNodeReal3& node_coord)
    {
      Real3 n0 = node_coord[item.nodeId(0)];
      Real3 n1 = node_coord[item.nodeId(1)];
      Real3 n2 = node_coord[item.nodeId(2)];

      auto v = math::cross(n1 - n0, n2 - n0);

      return v.normL2() / 2.0;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the area of a quadrilateral defined by four nodes.
     *
     * This method calculates the area of a quadrilateral by breaking it down
     * into two triangles and using the determinant formula. The area is computed
     * as half the value of the determinant of the matrix formed by the coordinates
     * of the quadrilateral's vertices.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real computeAreaQuad4(ItemWithNodes item, const VariableNodeReal3& node_coord)
    {
      Real3 n0 = node_coord[item.nodeId(0)];
      Real3 n1 = node_coord[item.nodeId(1)];
      Real3 n2 = node_coord[item.nodeId(2)];
      Real3 n3 = node_coord[item.nodeId(3)];

      auto tri1x2 = math::cross(n2 - n1, n0 - n1);
      auto tri2x2 = math::cross(n0 - n3, n2 - n3);
  
      return 0.5 * (tri1x2.normL2() + tri2x2.normL2());
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the volume of a hexahedron defined by eight nodes.
     */
    /*---------------------------------------------------------------------------
    static inline Real computeVolumeHexa8(ItemWithNodes item, const VariableNodeReal3& node_coord)
    {
      Real3 n0 = node_coord[item.nodeId(0)];
      Real3 n1 = node_coord[item.nodeId(1)];
      Real3 n2 = node_coord[item.nodeId(2)];
      Real3 n3 = node_coord[item.nodeId(3)];
      Real3 n4 = node_coord[item.nodeId(4)];
      Real3 n5 = node_coord[item.nodeId(5)];
      Real3 n6 = node_coord[item.nodeId(6)];
      Real3 n7 = node_coord[item.nodeId(7)];

      Real v1 = math::matDet((n6 - n1) + (n7 - n0), n6 - n3, n2 - n0);
      Real v2 = math::matDet(n7 - n0, (n6 - n3) + (n5 - n0), n6 - n4);
      Real v3 = math::matDet(n6 - n1, n5 - n0, (n6 - n4) + (n2 - n0));
      return (v1 + v2 + v3) / 12.;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the volume of a pentahedron (wedge or triangular prism)
     * defined by six nodes.
     */
    /*---------------------------------------------------------------------------

    static inline Real penta6Volume(ItemWithNodes item, const VariableNodeReal3& node_coord)
    {
      Real3 n0 = node_coord[item.nodeId(0)];
      Real3 n1 = node_coord[item.nodeId(1)];
      Real3 n2 = node_coord[item.nodeId(2)];
      Real3 n3 = node_coord[item.nodeId(3)];
      Real3 n4 = node_coord[item.nodeId(4)];
      Real3 n5 = node_coord[item.nodeId(5)];

      auto v = math::cross(n1 - n0, n2 - n0);
      auto base = 0.5 * v.normL2();
      auto h1 = (n3 - n0).normL2();
      auto h2 = (n4 - n1).normL2();
      auto h3 = (n5 - n2).normL2();

      return base * (h1 + h2 + h3) / 3.0;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the volume of a pyramid defined by five nodes.
     */
    /*---------------------------------------------------------------------------

    static inline Real pyramid5Volume(ItemWithNodes item, const VariableNodeReal3& node_coord)
    {
      Real3 n0 = node_coord[item.nodeId(0)];
      Real3 n1 = node_coord[item.nodeId(1)];
      Real3 n2 = node_coord[item.nodeId(2)];
      Real3 n3 = node_coord[item.nodeId(3)];
      Real3 n4 = node_coord[item.nodeId(4)];

      // Compute the area of the base triangle
      auto v = math::cross(n1 - n0, n2 - n0);

      auto tri1x2 = math::cross(n2 - n1, n0 - n1);
      auto tri2x2 = math::cross(n0 - n3, n2 - n3);
      auto base = 0.5 * (tri1x2.normL2() + tri2x2.normL2());


      // Compute the height of the pyramid
      // The height is the distance from the apex (n4) to the base plane
      // formed by the triangle (n0, n1, n2)
      auto h = math::dot(n4 - n0, v) / v.normL2();

      return (base * h) / 3.0;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the barycenter (centroid) of a triangle.
     *
     * This method calculates the barycenter of a triangle defined by three nodes.
     * The barycenter is computed as the average of the vertices' coordinates.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeBaryCenterTria3(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];

      Real Center_x = (vertex0.x + vertex1.x + vertex2.x) / 3.;
      Real Center_y = (vertex0.y + vertex1.y + vertex2.y) / 3.;

      return { Center_x, Center_y, 0 };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the barycenter (centroid) of a tetrahedron.
     *
     * This method calculates the barycenter of a tetrahedron defined by its nodes.
     * The barycenter is computed as the average of the vertices' coordinates.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeBaryCenterTetra4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];
      Real3 vertex3 = node_coord[cell.nodeId(3)];

      Real Center_x = (vertex0.x + vertex1.x + vertex2.x + vertex3.x) / 4.;
      Real Center_y = (vertex0.y + vertex1.y + vertex2.y + vertex3.x) / 4.;
      Real Center_z = (vertex0.z + vertex1.z + vertex2.z + vertex3.z) / 4.;

      return { Center_x, Center_y, Center_z };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the length of the edge defined by a given face.
     *
     * This method calculates Euclidean distance between the two nodes of the face.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real computeLengthEdge2(ItemWithNodes item, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[item.nodeId(0)];
      Real3 vertex1 = node_coord[item.nodeId(1)];

      return (vertex1 - vertex0).normL2();
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the normalized edge normal for a given face.
     *
     * This method calculates normal vector to the edge defined by nodes of the face,
     * normalizes it, and ensures the correct orientation.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real2 computeNormalEdge2(Face face, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[face.nodeId(0)];
      Real3 vertex1 = node_coord[face.nodeId(1)];

      if (!face.isSubDomainBoundaryOutside())
        std::swap(vertex0, vertex1);

      Real dx = vertex1.x - vertex0.x;
      Real dy = vertex1.y - vertex0.y;
      Real norm_N = math::sqrt(dx * dx + dy * dy);

      return { dy / norm_N, -dx / norm_N };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the factor used for integration of 1D, 2D and 3D finite elements
     *
     * This method calculates length or surface for a given finite-element (P1, P2, ...) and returns
     * the associated factor for elementary integrals.
     */
    /*---------------------------------------------------------------------------*/
    static inline Real computeFacLengthOrArea(Face face, const VariableNodeReal3& node_coord)
    {
      Int32 item_type = face.type();
      Real fac_el{ 0. };

      switch (item_type) {

      // Lines
      case IT_Line2:
      case IT_Line3:
        fac_el = computeLengthEdge2(face, node_coord) / 2.;
        break;

      // Faces
      case IT_Triangle3:
      case IT_Triangle6:
        fac_el = computeAreaTria3(face, node_coord) / 3.;
        break;

      case IT_Quad4:
      case IT_Quad8:
        fac_el = computeAreaQuad4(face, node_coord) / 4.;
        break;

      default:
        break;
      }
      return fac_el;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the normalized triangle normal for a given face.
     *
     * This method calculates normal vector to the triangle defined by nodes,
     * of the face and normalizes it, and ensures the correct orientation.
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeNormalTriangle(Face face, const VariableNodeReal3& node_coord)
    {

      // Get the vertices of the triangle
      Real3 vertex0 = node_coord[face.nodeId(0)];
      Real3 vertex1 = node_coord[face.nodeId(1)];
      Real3 vertex2 = node_coord[face.nodeId(2)];

      if (!face.isSubDomainBoundaryOutside())
        std::swap(vertex0, vertex1);

      // Calculate two edges of the triangle
      Real3 edge1 = { vertex1.x - vertex0.x, vertex1.y - vertex0.y, vertex1.z - vertex0.z };
      Real3 edge2 = { vertex2.x - vertex0.x, vertex2.y - vertex0.y, vertex2.z - vertex0.z };

      // Compute the cross product of the two edges
      Real3 normal = {
        edge1.y * edge2.z - edge1.z * edge2.y,
        edge1.z * edge2.x - edge1.x * edge2.z,
        edge1.x * edge2.y - edge1.y * edge2.x
      };

      // Calculate the magnitude of the normal vector
      Real norm = math::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

      // Normalize the vector to unit length and return
      return { normal.x / norm, normal.y / norm, normal.z / norm };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes finite-element entity (Edge, Face or Cell) geometric dimension
     * This method is used for the FEM 2D & 3D needs (coming from PASSMO)
     * for Jacobian & elementary matrices computations
     */
    /*---------------------------------------------------------------------------*/
    static inline Int32 getGeomDimension(ItemWithNodes item)
    {
      return item.typeInfo()->dimension();
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes edge & Face normal & tangent vectors (normalized, direct oriented)
     * This method is used for the FEM 2D & 3D needs (coming from PASSMO)
     * for geometric transformations (rotations, projections, ...)
     * In 2D, it assumes the edge lies in x-y plane (z coord = 0.)
     */
    /*---------------------------------------------------------------------------*/
    static inline void dirVectors(Face face, const VariableNodeReal3& node_coord, Integer ndim, Real3& e1, Real3& e2, Real3& e3)
    {
      Real3 n0 = node_coord[face.nodeId(0)];
      Real3 n1 = node_coord[face.nodeId(1)];

      if (!face.isSubDomainBoundaryOutside())
        std::swap(n0, n1);

      // 1st in-plane vector/along edge
      e1 = n1 - n0;

      if (ndim == 3) {

        Real3 n2 = node_coord[face.nodeId(2)];

        // out Normal to the face plane
        e3 = math::cross(e1, n2 - n0);

        // 2nd in-plane vector
        e2 = math::cross(e3, e1);
        e3 = math::mutableNormalize(e3);
      }
      else {

        Cell cell{ face.boundaryCell() };
        Node nod;
        for (Node node : cell.nodes()) {
          if (node != face.node(0) && node != face.node(1)) {
            nod = node;
            break;
          }
        }

        // Out Normal to the edge
        e2 = { -e1.y, e1.x, 0. };

        Real3 n2 = node_coord[nod];
        auto sgn = math::dot(e2, n2 - n0);
        if (sgn > 0.)
          e2 *= -1.;
      }
      e1 = math::mutableNormalize(e1);
      e2 = math::mutableNormalize(e2);
    }
  };

  /*---------------------------------------------------------------------------*/
  /**
  * @brief Provides methods for finite element operations in 2D.
  *
  * This class includes static methods for calculating gradients of basis
  * functions and integrals for P1 triangles in 2D finite element analysis.
  *
  * Reference Tri3 and Quad4 element in Arcane
  *
  *               0 o                    1 o . . . . o 0
  *                . .                     .         .
  *               .   .                    .         .
  *              .     .                   .         .
  *           1 o  . .  o 2              2 o . . . . o 3
  */
  /*---------------------------------------------------------------------------*/
  class FeOperation2D
  {

   public:

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the gradients of given function U for P1 triangles.
     *
     * This method calculates gradient operator ∇ Ui for a given P1 cell
     * with i = 1,..,3 for the three values of Ui hence  at  cell nodes.
     * The output is ∇ Ui is P0 (piece-wise constant) hence Real3  value
     * per cell
     *
     *         ∇ Ui = [ ∂U/∂𝑥   ∂U/∂𝑦   ∂U/∂𝑧 ]
     *
     *         ∂U/∂𝑥 = ( u1*(y2 − y3) + u2*(y3 − y1) + u3*(y1 − y2) ) / (2*A)
     *         ∂U/∂𝑦 = ( u1*(x3 − x2) + u2*(x1 − x3) + u3*(x2 − x1) ) / (2*A)
     *         ∂U/∂𝑧 = 0
     *
     * @note we can adapt the same for 3D by filling the third component
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeGradientTria3(Cell cell, const VariableNodeReal3& node_coord, VariableNodeReal u)
    {
      Real3 n0 = node_coord[cell.nodeId(0)];
      Real3 n1 = node_coord[cell.nodeId(1)];
      Real3 n2 = node_coord[cell.nodeId(2)];

      Real u0 = u[cell.nodeId(0)];
      Real u1 = u[cell.nodeId(1)];
      Real u2 = u[cell.nodeId(2)];

      Real A2 = ((n1.x - n0.x) * (n2.y - n0.y) - (n2.x - n0.x) * (n1.y - n0.y));

      return { (u0 * (n1.y - n2.y) + u1 * (n2.y - n0.y) + u2 * (n0.y - n1.y)) / A2, (u0 * (n2.x - n1.x) + u1 * (n0.x - n2.x) + u2 * (n1.x - n0.x)) / A2, 0 };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the 𝑥 gradients of basis functions 𝐍 for ℙ1 triangles.
     *
     * This method calculates gradient operator ∂/∂𝑥 of 𝑁ᵢ for a given ℙ1
     * cell with i = 1,..,3 for the three shape function  𝑁ᵢ  hence output
     * is a vector of size 3
     *
     *         ∂𝐍/∂𝑥 = [ ∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ∂𝑁₃/∂𝑥 ]
     *
     *         ∂𝐍/∂𝑥 = 1/(2𝐴) [ 𝑦₂-𝑦₃  𝑦₃-𝑦₁  𝑦₁-𝑦₂ ]
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeGradientXTria3(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];

      auto A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

      return { (vertex1.y - vertex2.y) / A2, (vertex2.y - vertex0.y) / A2, (vertex0.y - vertex1.y) / A2 };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the 𝑦 gradients of basis functions 𝐍 for ℙ1 triangles.
     *
     * This method calculates gradient operator ∂/∂𝑦 of 𝑁ᵢ for a given ℙ1
     * cell with i = 1,..,3 for the three shape function  𝑁ᵢ  hence output
     * is a vector of size 3
     *
     *         ∂𝐍/∂𝑦 = [ ∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ∂𝑁₃/∂𝑦 ]
     *
     *         ∂𝐍/∂𝑦 = 1/(2𝐴) [ 𝑥₃−𝑥₂  𝑥₁−𝑥₃  𝑥₂−𝑥₁ ]
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 computeGradientYTria3(Cell cell, VariableNodeReal3 node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];

      auto A2 = ((vertex1.x - vertex0.x) * (vertex2.y - vertex0.y) - (vertex2.x - vertex0.x) * (vertex1.y - vertex0.y));

      return { (vertex2.x - vertex1.x) / A2, (vertex0.x - vertex2.x) / A2, (vertex1.x - vertex0.x) / A2 };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Holds information for a Quad4 element at a single Gauss point.
     *
     * This includes the gradients of the shape functions in the physical space (x, y)
     * and the determinant of the Jacobian matrix.
     */
    /*---------------------------------------------------------------------------*/
    struct Quad4GaussPointInfo
    {
      RealVector<4> dN_dx; // Derivatives of shape functions in x {∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ∂𝑁₃/∂𝑥  ∂𝑁₄/∂𝑥}
      RealVector<4> dN_dy; // Derivatives of shape functions in y {∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ∂𝑁₃/∂𝑦  ∂𝑁₄/∂𝑦}
      Real det_j; // Determinant of the Jacobian matrix at the Gauss point.
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes shape function gradients and the Jacobian determinant for a Quad4 element.
     *
     * @param cell The Quad4 cell entity.
     * @param node_coord The coordinates of the mesh nodes.
     * @param xi The ξ coordinate of the evaluation point (-1 to 1).
     * @param eta The η coordinate of the evaluation point (-1 to 1).
     * @return A Quad4GaussPointInfo struct containing {∂𝐍/∂𝑥, ∂𝐍/∂𝑦, det(𝑱)}.
     */
    /*---------------------------------------------------------------------------*/
    static inline Quad4GaussPointInfo
    computeGradientsAndJacobianQuad4(Cell cell, const VariableNodeReal3& node_coord, Real xi, Real eta)
    {
      // Shape function derivatives ∂𝐍/∂ξ and ∂𝐍/∂η
      //     ∂𝐍/∂ξ = [ ∂C₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ ]
      //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η ]
      const Real dN_dxi[4] = { -0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta) };
      const Real dN_deta[4] = { -0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi) };

      // Jacobian calculation 𝑱
      //    𝑱 = [ 𝒋₀₀  𝒋₀₁ ] = [ ∂𝑥/∂ξ  ∂𝑦/∂ξ ]
      //        [ 𝒋₁₀  𝒋₁₁ ]   [ ∂𝑥/∂η  ∂𝑦/∂η ]
      //
      // The Jacobian is computed as follows:
      //   𝒋₀₀ = ∑ (∂𝑁ᵢ/∂ξ * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟒
      //   𝒋₀₁ = ∑ (∂𝑁ᵢ/∂ξ * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟒
      //   𝒋₁₀ = ∑ (∂𝑁ᵢ/∂η * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟒
      //   𝒋₁₁ = ∑ (∂𝑁ᵢ/∂η * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟒
      Real2x2 J;
      for (Int8 a = 0; a < 4; ++a) {
        const auto& coord = node_coord[cell.nodeId(a)];
        J[0][0] += dN_dxi[a] * coord.x;
        J[0][1] += dN_dxi[a] * coord.y;
        J[1][0] += dN_deta[a] * coord.x;
        J[1][1] += dN_deta[a] * coord.y;
      }

      // Determinant of the Jacobian
      const Real detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

      if (detJ <= 0.0) {
        ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
      }

      // Inverse of the Jacobian
      //    𝑱⁻¹ = [ invJ00 invJ01 ]
      //          [ invJ10 invJ11 ]
      const Real invJ00 = J[1][1] / detJ;
      const Real invJ01 = -J[0][1] / detJ;
      const Real invJ10 = -J[1][0] / detJ;
      const Real invJ11 = J[0][0] / detJ;

      //   Gradients in physical space (∂𝐍/∂𝑥, ∂𝐍/∂𝑦)
      //    {∂𝐍/∂𝑥} = [J]⁻¹ {∂𝐍/∂ξ}
      //    {∂𝐍/∂𝑦}         {∂𝐍/∂η}
      RealVector<4> dN_dx_result, dN_dy_result;
      for (Int8 a = 0; a < 4; ++a) {
        dN_dx_result(a) = invJ00 * dN_dxi[a] + invJ01 * dN_deta[a];
        dN_dy_result(a) = invJ10 * dN_dxi[a] + invJ11 * dN_deta[a];
      }

      return { dN_dx_result, dN_dy_result, detJ };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the gradient of a scalar field 'u' for a Quad4 element.
     *
     * @param cell The Quad4 cell entity.
     * @param node_coord The coordinates of the mesh nodes.
     * @param u The scalar variable 𝑢 defined at the nodes.
     * @return A Real3 vector of the gradient ∇𝑢 = {∂𝑢/∂𝑥, ∂𝑢/∂𝑦, 0} at the cell center.
     */
    /*---------------------------------------------------------------------------*/
    static inline Real3
    computeGradientQuad4(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal& u/*, Real xi, Real eta*/)
    {
      // get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian at (ξ,η) = (0,0)
      const auto gp_util = computeGradientsAndJacobianQuad4(cell, node_coord, 0.0, 0.0);
      const RealVector<4>& dN_dx = gp_util.dN_dx;
      const RealVector<4>& dN_dy = gp_util.dN_dy;

      // get the nodal values of the variable 𝑢ᵢ ∀ 𝑖= 1,……,4 for the cell.
      const Real u_nodes[4] = {
        u[cell.nodeId(0)],
        u[cell.nodeId(1)],
        u[cell.nodeId(2)],
        u[cell.nodeId(3)]
      };

      // Compute the gradient components using shape function
      //    ∂𝑢/∂𝑥 = Σ (∂𝑁ᵢ/∂𝑥 * uᵢ) ∀ 𝑖= 1,……,4
      //    ∂𝑢/∂𝑦 = Σ (∂𝑁ᵢ/∂𝑦 * uᵢ) ∀ 𝑖= 1,……,4
      Real grad_x = 0.0;
      Real grad_y = 0.0;
      for (Int8 a = 0; a < 4; ++a) {
        grad_x += dN_dx(a) * u_nodes[a];
        grad_y += dN_dy(a) * u_nodes[a];
      }

      return { grad_x, grad_y, 0.0 };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the shape functions for a Quad4 element at a given point (ξ, η).
     * 
     * @param xi The ξ coordinate of the evaluation point (-1 to 1).
     * @param eta The η coordinate of the evaluation point (-1 to 1).
     * @return A RealVector<4> containing the shape functions {𝑁₁, 𝑁₂, 𝑁₃, 𝑁₄}.
     */
    /*---------------------------------------------------------------------------*/
    static inline RealVector<4> computeShapeFunctionsQuad4(Real xi, Real eta)
    {
      RealVector<4> N;
      N(0) = 0.25 * (1. - xi) * (1. - eta); // 𝑁₁
      N(1) = 0.25 * (1. + xi) * (1. - eta); // 𝑁₂
      N(2) = 0.25 * (1. + xi) * (1. + eta); // 𝑁₃
      N(3) = 0.25 * (1. - xi) * (1. + eta); // 𝑁₄
      return N;
    }
  };

  /*---------------------------------------------------------------------------*/
  /**
  * @brief Provides methods for finite element operations in 3D.
  *
  * This class includes static methods for calculating gradients of basis
  * functions and integrals for P1 triangles in 3D finite element analysis.
  *
  * Reference Tetra4 element in Arcane
  *
  *               3 o
  *                /|\
  *               / | \
  *              /  |  \
  *             /   o 2 \
  *            / .    .  \
  *         0 o-----------o 1
  */
  /*---------------------------------------------------------------------------*/
  class FeOperation3D
  {
   public:

    /*-------------------------------------------------------------------------*/
    /**
     * @brief Computes the X gradients of basis functions 𝐍 for ℙ1 Tetrahedron.
     *
     * This method calculates gradient operator ∂/∂𝑥 of 𝑁ᵢ for a given ℙ1
     * cell with i = 1,..,4 for the four shape function  𝑁ᵢ  hence output
     * is a vector of size 4
     *
     *         ∂𝐍/∂𝑥 = [ ∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ∂𝑁₃/∂𝑥  ∂𝑁₄/∂𝑥 ]
     *
     *         ∂𝐍/∂𝑥 = 1/(6𝑉) [ dx0  dx1  dx2  dx3 ]
     *
     * where:
     *    dx0 = (n1.y * (n3.z - n2.z) + n2.y * (n1.z - n3.z) + n3.y * (n2.z - n1.z)),
     *    dx1 = (n0.y * (n2.z - n3.z) + n2.y * (n3.z - n0.z) + n3.y * (n0.z - n2.z)),
     *    dx2 = (n0.y * (n3.z - n1.z) + n1.y * (n0.z - n3.z) + n3.y * (n1.z - n0.z)),
     *    dx3 = (n0.y * (n1.z - n2.z) + n1.y * (n2.z - n0.z) + n2.y * (n0.z - n1.z)).
     *
     */
    /*-------------------------------------------------------------------------*/

    static inline Real4 computeGradientXTetra4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 n0 = node_coord[cell.nodeId(0)];
      Real3 n1 = node_coord[cell.nodeId(1)];
      Real3 n2 = node_coord[cell.nodeId(2)];
      Real3 n3 = node_coord[cell.nodeId(3)];

      Real3 v0 = n1 - n0;
      Real3 v1 = n2 - n0;
      Real3 v2 = n3 - n0;

      // 6 x Volume of tetrahedron
      Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

      Real4 dx{};

      dx[0] = (n1.y * (n3.z - n2.z) + n2.y * (n1.z - n3.z) + n3.y * (n2.z - n1.z)) / V6;
      dx[1] = (n0.y * (n2.z - n3.z) + n2.y * (n3.z - n0.z) + n3.y * (n0.z - n2.z)) / V6;
      dx[2] = (n0.y * (n3.z - n1.z) + n1.y * (n0.z - n3.z) + n3.y * (n1.z - n0.z)) / V6;
      dx[3] = (n0.y * (n1.z - n2.z) + n1.y * (n2.z - n0.z) + n2.y * (n0.z - n1.z)) / V6;

      return dx;
    }

    /*-------------------------------------------------------------------------*/
    /**
     * @brief Computes the Y gradients of basis functions 𝐍 for ℙ1 Tetrahedron.
     *
     * This method calculates gradient operator ∂/∂𝑦 of 𝑁ᵢ for a given ℙ1
     * cell with i = 1,..,4 for the four shape functions 𝑁ᵢ, hence the output
     * is a vector of size 4.
     *
     *         ∂𝐍/∂𝑦 = [ ∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ∂𝑁₃/∂𝑦  ∂𝑁₄/∂𝑦 ]
     *
     *         ∂𝐍/∂𝑦 = 1/(6𝑉) [ dy0  dy1  dy2  dy3 ]
     *
     * where:
     *    dy0 = (n1.z * (n3.x - n2.x) + n2.z * (n1.x - n3.x) + n3.z * (n2.x - n1.x)),
     *    dy1 = (n0.z * (n2.x - n3.x) + n2.z * (n3.x - n0.x) + n3.z * (n0.x - n2.x)),
     *    dy2 = (n0.z * (n3.x - n1.x) + n1.z * (n0.x - n3.x) + n3.z * (n1.x - n0.x)),
     *    dy3 = (n0.z * (n1.x - n2.x) + n1.z * (n2.x - n0.x) + n2.z * (n0.x - n1.x)).
     *
     */
    /*-------------------------------------------------------------------------*/

    static inline Real4 computeGradientYTetra4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 n0 = node_coord[cell.nodeId(0)];
      Real3 n1 = node_coord[cell.nodeId(1)];
      Real3 n2 = node_coord[cell.nodeId(2)];
      Real3 n3 = node_coord[cell.nodeId(3)];

      Real3 v0 = n1 - n0;
      Real3 v1 = n2 - n0;
      Real3 v2 = n3 - n0;

      // 6 x Volume of tetrahedron
      Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

      Real4 dy{};

      dy[0] = (n1.z * (n3.x - n2.x) + n2.z * (n1.x - n3.x) + n3.z * (n2.x - n1.x)) / V6;
      dy[1] = (n0.z * (n2.x - n3.x) + n2.z * (n3.x - n0.x) + n3.z * (n0.x - n2.x)) / V6;
      dy[2] = (n0.z * (n3.x - n1.x) + n1.z * (n0.x - n3.x) + n3.z * (n1.x - n0.x)) / V6;
      dy[3] = (n0.z * (n1.x - n2.x) + n1.z * (n2.x - n0.x) + n2.z * (n0.x - n1.x)) / V6;

      return dy;
    }

    /*-------------------------------------------------------------------------*/
    /**
     * @brief Computes the 𝑧 gradients of basis functions 𝐍 for ℙ1 Tetrahedron.
     *
     * This method calculates gradient operator ∂/∂𝑧 of 𝑁ᵢ for a given ℙ1
     * cell with i = 1,..,4 for the four shape functions 𝑁ᵢ, hence the output
     * is a vector of size 4.
     *
     *         ∂𝐍/∂𝑧 = [ ∂𝑁₁/∂𝑧  ∂𝑁₂/∂𝑧  ∂𝑁₃/∂𝑧  ∂𝑁₄/∂𝑧 ]
     *
     *         ∂𝐍/∂𝑧 = 1/(6𝑉) [ dz0  dz1  dz2  dz3 ]
     *
     * where:
     *    dz0 = (n1.x * (n3.y - n2.y) + n2.x * (n1.y - n3.y) + n3.x * (n2.y - n1.y)),
     *    dz1 = (n0.x * (n2.y - n3.y) + n2.x * (n3.y - n0.y) + n3.x * (n0.y - n2.y)),
     *    dz2 = (n0.x * (n3.y - n1.y) + n1.x * (n0.y - n3.y) + n3.x * (n1.y - n0.y)),
     *    dz3 = (n0.x * (n1.y - n2.y) + n1.x * (n2.y - n0.y) + n2.x * (n0.y - n1.y)).
     *
     */
    /*-------------------------------------------------------------------------*/

    static inline Real4 computeGradientZTetra4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 n0 = node_coord[cell.nodeId(0)];
      Real3 n1 = node_coord[cell.nodeId(1)];
      Real3 n2 = node_coord[cell.nodeId(2)];
      Real3 n3 = node_coord[cell.nodeId(3)];

      auto v0 = n1 - n0;
      auto v1 = n2 - n0;
      auto v2 = n3 - n0;

      // 6 x Volume of tetrahedron
      Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

      Real4 dz{};

      dz[0] = (n1.x * (n3.y - n2.y) + n2.x * (n1.y - n3.y) + n3.x * (n2.y - n1.y)) / V6;
      dz[1] = (n0.x * (n2.y - n3.y) + n2.x * (n3.y - n0.y) + n3.x * (n0.y - n2.y)) / V6;
      dz[2] = (n0.x * (n3.y - n1.y) + n1.x * (n0.y - n3.y) + n3.x * (n1.y - n0.y)) / V6;
      dz[3] = (n0.x * (n1.y - n2.y) + n1.x * (n2.y - n0.y) + n2.x * (n0.y - n1.y)) / V6;

      return dz;
    }


    static inline Real3 computeGradientTetra4(Cell cell, const VariableNodeReal3& node_coord, VariableNodeReal u)
    {
      Real3 m0 = node_coord[cell.nodeId(0)];
      Real3 m1 = node_coord[cell.nodeId(1)];
      Real3 m2 = node_coord[cell.nodeId(2)];
      Real3 m3 = node_coord[cell.nodeId(3)];

      Real f0 = u[cell.nodeId(0)];
      Real f1 = u[cell.nodeId(1)];
      Real f2 = u[cell.nodeId(2)];
      Real f3 = u[cell.nodeId(3)];

      Real3 v0 = m1 - m0;
      Real3 v1 = m2 - m0;
      Real3 v2 = m3 - m0;

      // 6 x Volume of tetrahedron
      Real V6 = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

      // Compute gradient components
      Real3 grad;
      grad.x = (f0 * (m1.y * m2.z + m2.y * m3.z + m3.y * m1.z - m3.y * m2.z - m2.y * m1.z - m1.y * m3.z)
               - f1 * (m0.y * m2.z + m2.y * m3.z + m3.y * m0.z - m3.y * m2.z - m2.y * m0.z - m0.y * m3.z)
               + f2 * (m0.y * m1.z + m1.y * m3.z + m3.y * m0.z - m3.y * m1.z - m1.y * m0.z - m0.y * m3.z)
               - f3 * (m0.y * m1.z + m1.y * m2.z + m2.y * m0.z - m2.y * m1.z - m1.y * m0.z - m0.y * m2.z)) / V6;
      grad.y = (f0 * (m1.z * m2.x + m2.z * m3.x + m3.z * m1.x - m3.z * m2.x - m2.z * m1.x - m1.z * m3.x)
               - f1 * (m0.z * m2.x + m2.z * m3.x + m3.z * m0.x - m3.z * m2.x - m2.z * m0.x - m0.z * m3.x)
               + f2 * (m0.z * m1.x + m1.z * m3.x + m3.z * m0.x - m3.z * m1.x - m1.z * m0.x - m0.z * m3.x)
               - f3 * (m0.z * m1.x + m1.z * m2.x + m2.z * m0.x - m2.z * m1.x - m1.z * m0.x - m0.z * m2.x)) / V6;
      grad.z = (f0 * (m1.x * m2.y + m2.x * m3.y + m3.x * m1.y - m3.x * m2.y - m2.x * m1.y - m1.x * m3.y)
               - f1 * (m0.x * m2.y + m2.x * m3.y + m3.x * m0.y - m3.x * m2.y - m2.x * m0.y - m0.x * m3.y)
               + f2 * (m0.x * m1.y + m1.x * m3.y + m3.x * m0.y - m3.x * m1.y - m1.x * m0.y - m0.x * m3.y)
               - f3 * (m0.x * m1.y + m1.x * m2.y + m2.x * m0.y - m2.x * m1.y - m1.x * m0.y - m0.x * m2.y)) / V6;
      return grad;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Holds information for a Hexa8 element at a single Gauss point.
     *
     * This includes the gradients of the shape functions in the physical space (𝑥,𝑦,𝑧)
     * and the determinant of the Jacobian matrix.
     */
    /*---------------------------------------------------------------------------*/
    struct Hexa8GaussPointInfo
    {
      RealVector<8> dN_dx; // Derivatives of shape functions in x {∂𝑁₁/∂𝑥  ∂𝑁₂/∂𝑥  ...  ∂𝑁₈/∂𝑥}
      RealVector<8> dN_dy; // Derivatives of shape functions in y {∂𝑁₁/∂𝑦  ∂𝑁₂/∂𝑦  ...  ∂𝑁₈/∂𝑦}
      RealVector<8> dN_dz; // Derivatives of shape functions in z {∂𝑁₁/∂𝑧  ∂𝑁₂/∂𝑧  ...  ∂𝑁₈/∂𝑧}
      Real det_j; // Determinant of the Jacobian matrix at the Gauss point.
    };

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes shape function gradients and the Jacobian determinant for a Hexa8 element.
     *
     * @param cell The Hexa8 cell entity.
     * @param node_coord The coordinates of the mesh nodes.
     * @param xi The ξ coordinate of the evaluation point (-1 to 1).
     * @param eta The η coordinate of the evaluation point (-1 to 1).
     * @param zeta The ζ coordinate of the evaluation point (-1 to 1).
     * @return A Hexa8GaussPointInfo struct containing the gradients and Jacobian determinant.
     */
    /*---------------------------------------------------------------------------*/

    static inline Hexa8GaussPointInfo
    computeGradientsAndJacobianHexa8(Cell cell, const VariableNodeReal3& node_coord, Real xi, Real eta, Real zeta)
    {
      // Shape function derivatives ∂𝐍/∂ξ, ∂𝐍/∂η, ∂𝐍/∂ζ
      //     ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ  ∂𝑁₅/∂ξ  ∂𝑁₆/∂ξ  ∂𝑁₇/∂ξ  ∂𝑁₈/∂ξ ]
      //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η  ∂𝑁₅/∂η  ∂𝑁₆/∂η  ∂𝑁₇/∂η  ∂𝑁₈/∂η ]
      //     ∂𝐍/∂ζ = [ ∂𝑁₁/∂ζ  ∂𝑁₂/∂ζ  ∂𝑁₃/∂ζ  ∂𝑁₄/∂ζ  ∂𝑁₅/∂ζ  ∂𝑁₆/∂ζ  ∂𝑁₇/∂ζ  ∂𝑁₈/∂ζ ]
      Real dN_dxi[8], dN_deta[8], dN_dzeta[8];
      const Real one_minus_eta = 1.0 - eta;
      const Real one_plus_eta = 1.0 + eta;
      const Real one_minus_xi = 1.0 - xi;
      const Real one_plus_xi = 1.0 + xi;
      const Real one_minus_zeta = 1.0 - zeta;
      const Real one_plus_zeta = 1.0 + zeta;

      dN_dxi[0] = -0.125 * one_minus_eta * one_minus_zeta;
      dN_dxi[1] = 0.125 * one_minus_eta * one_minus_zeta;
      dN_dxi[2] = 0.125 * one_plus_eta * one_minus_zeta;
      dN_dxi[3] = -0.125 * one_plus_eta * one_minus_zeta;
      dN_dxi[4] = -0.125 * one_minus_eta * one_plus_zeta;
      dN_dxi[5] = 0.125 * one_minus_eta * one_plus_zeta;
      dN_dxi[6] = 0.125 * one_plus_eta * one_plus_zeta;
      dN_dxi[7] = -0.125 * one_plus_eta * one_plus_zeta;

      dN_deta[0] = -0.125 * one_minus_xi * one_minus_zeta;
      dN_deta[1] = -0.125 * one_plus_xi * one_minus_zeta;
      dN_deta[2] = 0.125 * one_plus_xi * one_minus_zeta;
      dN_deta[3] = 0.125 * one_minus_xi * one_minus_zeta;
      dN_deta[4] = -0.125 * one_minus_xi * one_plus_zeta;
      dN_deta[5] = -0.125 * one_plus_xi * one_plus_zeta;
      dN_deta[6] = 0.125 * one_plus_xi * one_plus_zeta;
      dN_deta[7] = 0.125 * one_minus_xi * one_plus_zeta;

      dN_dzeta[0] = -0.125 * one_minus_xi * one_minus_eta;
      dN_dzeta[1] = -0.125 * one_plus_xi * one_minus_eta;
      dN_dzeta[2] = -0.125 * one_plus_xi * one_plus_eta;
      dN_dzeta[3] = -0.125 * one_minus_xi * one_plus_eta;
      dN_dzeta[4] = 0.125 * one_minus_xi * one_minus_eta;
      dN_dzeta[5] = 0.125 * one_plus_xi * one_minus_eta;
      dN_dzeta[6] = 0.125 * one_plus_xi * one_plus_eta;
      dN_dzeta[7] = 0.125 * one_minus_xi * one_plus_eta;

      // Jacobian matrix (default-initialized to zero see Real3x3.h)
      //    𝑱 = [ 𝒋₀₀  𝒋₀₁  𝒋₀₂ ]
      //        [ 𝒋₁₀  𝒋₁₁  𝒋₁₂ ]
      //        [ 𝒋₂₀  𝒋₂₁  𝒋₂₂ ]
      //
      // The Jacobian is computed as follows:
      //   𝒋₀₀ = ∑ (∂𝑥/∂ξ * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟖
      //   𝒋₀₁ = ∑ (∂𝑥/∂ξ * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟖
      //   𝒋₀₂ = ∑ (∂𝑥/∂ξ * 𝑧ᵢ) ∀ 𝑖= 𝟏,……,𝟖
      //   𝒋₁₀ = ∑ (∂𝑦/∂η * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟖
      //   𝒋₁₁ = ∑ (∂𝑦/∂η * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟖
      //   𝒋₁₂ = ∑ (∂𝑦/∂η * 𝑧ᵢ) ∀ 𝑖= 𝟏,……,𝟖
      //   𝒋₂₀ = ∑ (∂𝑧/∂ζ * 𝑥ᵢ) ∀ 𝑖= 𝟏,……,𝟖
      //   𝒋₂₁ = ∑ (∂𝑧/∂ζ * 𝑦ᵢ) ∀ 𝑖= 𝟏,……,𝟖
      //   𝒋₂₂ = ∑ (∂𝑧/∂ζ * 𝑧ᵢ) ∀ 𝑖= 𝟏,……,𝟖
      Real3x3 J;
      for (Int8 a = 0; a < 8; ++a) {
        const Real3& n_coord = node_coord[cell.nodeId(a)];
        J[0][0] += dN_dxi[a] * n_coord.x; // ∂𝑥/∂ξ
        J[0][1] += dN_dxi[a] * n_coord.y; // ∂𝑦/∂ξ
        J[0][2] += dN_dxi[a] * n_coord.z; // ∂𝑧/∂ξ
        J[1][0] += dN_deta[a] * n_coord.x; // ∂𝑥/∂η
        J[1][1] += dN_deta[a] * n_coord.y; // ∂𝑦/∂η
        J[1][2] += dN_deta[a] * n_coord.z; // ∂𝑧/∂η
        J[2][0] += dN_dzeta[a] * n_coord.x; // ∂𝑥/∂ζ
        J[2][1] += dN_dzeta[a] * n_coord.y; // ∂𝑦/∂ζ
        J[2][2] += dN_dzeta[a] * n_coord.z; // ∂𝑧/∂ζ
      }

      // Determinant and Inverse of the Jacobian
      const Real detJ = math::matrixDeterminant(J);
      if (detJ <= 0.0) {
        ARCANE_FATAL("Invalid (non-positive) Jacobian determinant: {0}", detJ);
      }
      const Real3x3 invJ = math::inverseMatrix(J, detJ);

      // Gradients in physical space (∂𝐍/∂𝑥, ∂𝐍/∂𝑦, ∂𝐍/∂𝑧)
      //    {∂𝐍/∂𝑥} = [J]⁻¹ {∂𝐍/∂ξ}
      //    {∂𝐍/∂𝑦}         {∂𝐍/∂η}
      //    {∂𝐍/∂𝑧}         {∂𝐍/∂ζ}
      RealVector<8> dN_dx_result, dN_dy_result, dN_dz_result;
      for (Int8 a = 0; a < 8; ++a) {
        dN_dx_result(a) = invJ[0][0] * dN_dxi[a] + invJ[0][1] * dN_deta[a] + invJ[0][2] * dN_dzeta[a];
        dN_dy_result(a) = invJ[1][0] * dN_dxi[a] + invJ[1][1] * dN_deta[a] + invJ[1][2] * dN_dzeta[a];
        dN_dz_result(a) = invJ[2][0] * dN_dxi[a] + invJ[2][1] * dN_deta[a] + invJ[2][2] * dN_dzeta[a];
      }

      return { dN_dx_result, dN_dy_result, dN_dz_result, detJ };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the gradient of a scalar field 'u' for a Hexa8 element.
     *
     * @param cell The Hexa8 cell entity.
     * @param node_coord The coordinates of the mesh nodes.
     * @param u The scalar variable 𝑢 defined at the nodes.
     * @return A Real3 vector of the gradient ∇𝑢 = {∂𝑢/∂𝑥, ∂𝑢/∂𝑦, ∂𝑢/∂𝑧} at the cell center.
     */
    /*---------------------------------------------------------------------------*/
    static inline Real3
    computeGradientHexa8(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal& u /*, Real xi, Real eta*/)
    {
      // get shape function gradients w.r.t (𝑥,𝑦) and determinant of Jacobian at (ξ,η,ζ) = (0,0,0)
      const auto gp_util = computeGradientsAndJacobianHexa8(cell, node_coord, 0.0, 0.0, 0.0);
      const RealVector<8>& dN_dx = gp_util.dN_dx;
      const RealVector<8>& dN_dy = gp_util.dN_dy;
      const RealVector<8>& dN_dz = gp_util.dN_dz;

      // get the nodal values of the variable 𝑢ᵢ ∀ 𝑖= 1,……,8 for the cell.
      const Real u_nodes[8] = {
        u[cell.nodeId(0)],
        u[cell.nodeId(1)],
        u[cell.nodeId(2)],
        u[cell.nodeId(3)],
        u[cell.nodeId(4)],
        u[cell.nodeId(5)],
        u[cell.nodeId(6)],
        u[cell.nodeId(7)]
      };

      // Compute the gradient components using shape function
      //    ∂𝑢/∂𝑥 = Σ (∂𝑁ᵢ/∂𝑥 * uᵢ) ∀ 𝑖= 1,……,8
      //    ∂𝑢/∂𝑦 = Σ (∂𝑁ᵢ/∂𝑦 * uᵢ) ∀ 𝑖= 1,……,8
      //    ∂𝑢/∂𝑧 = Σ (∂𝑁ᵢ/∂𝑧 * uᵢ) ∀ 𝑖= 1,……,8
      Real grad_x = 0.0;
      Real grad_y = 0.0;
      Real grad_z = 0.0;
      for (Int8 a = 0; a < 8; ++a) {
        grad_x += dN_dx(a) * u_nodes[a];
        grad_y += dN_dy(a) * u_nodes[a];
        grad_z += dN_dz(a) * u_nodes[a];
      }

      return { grad_x, grad_y, grad_z };
    }


    /*---------------------------------------------------------------------------*/
    /**
     * @brief Computes the shape functions for a Hexa8 element at a given point (ξ,η,ζ).
     * 
     * @param xi The ξ coordinate of the evaluation point (-1 to 1).
     * @param eta The η coordinate of the evaluation point (-1 to 1).
     * @param zeta The ζ coordinate of the evaluation point (-1 to 1).
     * @return A RealVector<8> containing the shape functions {𝑁₁, 𝑁₂, 𝑁₃, 𝑁₄, 𝑁₅, 𝑁₆, 𝑁₇, 𝑁₈}.
     */
    /*---------------------------------------------------------------------------*/
    static inline RealVector<8> computeShapeFunctionsHexa8(Real xi, Real eta, Real zeta)
    {
      RealVector<8> N;
      const Real one_minus_eta = 1.0 - eta;
      const Real one_plus_eta = 1.0 + eta;
      const Real one_minus_xi = 1.0 - xi;
      const Real one_plus_xi = 1.0 + xi;
      const Real one_minus_zeta = 1.0 - zeta;
      const Real one_plus_zeta = 1.0 + zeta;

      N(0) = 0.125 * one_minus_xi * one_minus_eta * one_minus_zeta; // 𝑁₁
      N(1) = 0.125 * one_plus_xi * one_minus_eta * one_minus_zeta; // 𝑁₂
      N(2) = 0.125 * one_plus_xi * one_plus_eta * one_minus_zeta; // 𝑁₃
      N(3) = 0.125 * one_minus_xi * one_plus_eta * one_minus_zeta; // 𝑁₄
      N(4) = 0.125 * one_minus_xi * one_minus_eta * one_plus_zeta; // 𝑁₅
      N(5) = 0.125 * one_plus_xi * one_minus_eta * one_plus_zeta; // 𝑁₆
      N(6) = 0.125 * one_plus_xi * one_plus_eta * one_plus_zeta; // 𝑁₇
      N(7) = 0.125 * one_minus_xi * one_plus_eta * one_plus_zeta; // 𝑁₈

      return N;
    }
  };

  /*---------------------------------------------------------------------------*/
  /**
     * @brief Provides methods to help build boundary conditions in 2D/3D.
     */
  /*---------------------------------------------------------------------------*/
  class BoundaryConditionsHelpers
  {
   public:

    static inline void applyDirichletToNodeGroupRhsOnly(const Int32 dof_index, Real value, const IndexedNodeDoFConnectivityView& node_dof, VariableDoFReal& rhs_values, NodeGroup& node_group)
    {
      ENUMERATE_ (Node, inode, node_group) {
        Node node = *inode;
        if (node.isOwn()) {
          rhs_values[node_dof.dofId(node, dof_index)] = value;
        }
      }
    }

    static inline void applyDirichletToNodeGroupViaPenalty(const Int32 dof_index, Real value, Real penalty, const IndexedNodeDoFConnectivityView& node_dof, DoFLinearSystem& linear_system, VariableDoFReal& rhs_values, NodeGroup& node_group)
    {
      ENUMERATE_ (Node, inode, node_group) {
        Node node = *inode;
        if (node.isOwn()) {
          linear_system.matrixSetValue(node_dof.dofId(node, dof_index), node_dof.dofId(node, dof_index), penalty);
          Real u_g = penalty * value;
          rhs_values[node_dof.dofId(node, dof_index)] = u_g;
        }
      }
    }

    static inline void applyDirichletToNodeGroupViaRowElimination(const Int32 dof_index, Real value, const IndexedNodeDoFConnectivityView& node_dof, DoFLinearSystem& linear_system, VariableDoFReal& rhs_values, NodeGroup& node_group)
    {
      ENUMERATE_ (Node, inode, node_group) {
        Node node = *inode;
        if (node.isOwn()) {
          linear_system.eliminateRow(node_dof.dofId(*inode, dof_index), value);
        }
      }
    }

    static inline void applyDirichletToNodeGroupViaRowColumnElimination(const Int32 dof_index, Real value, const IndexedNodeDoFConnectivityView& node_dof, DoFLinearSystem& linear_system, VariableDoFReal& rhs_values, NodeGroup& node_group)
    {
      ENUMERATE_ (Node, inode, node_group) {
        Node node = *inode;
        if (node.isOwn()) {
          linear_system.eliminateRowColumn(node_dof.dofId(*inode, dof_index), value);
        }
      }
    }
  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for applying boundary conditions in 2D/3D FEM.
   *
   * This class includes static methods for applying:
   * - Dirichlet boundary condition
   * - Point Dirichlet condition
   */
  /*---------------------------------------------------------------------------*/
  class BoundaryConditions
  {
   public:

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies Dirichlet boundary conditions to RHS and LHS.
     *
     * Updates the LHS matrix and RHS vector to enforce Dirichlet conditions.
     *
     * - For LHS matrix `𝐀`, the diagonal term for the Dirichlet DOF is set to `𝑃`.
     * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
     *
     * @param bs Boundary condition values.
     * @param node_dof DOF connectivity view.
     * @param node_coord Node coordinates.
     * @param m_linear_system Linear system for LHS.
     * @param rhs_values RHS values to update.
     */
    /*---------------------------------------------------------------------------*/
    static inline void applyDirichletToLhsAndRhs(BC::IDirichletBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, DoFLinearSystem& m_linear_system, VariableDoFReal& rhs_values)
    {
      FaceGroup face_group = bs->getSurface();
      NodeGroup node_group = face_group.nodeGroup();
      const StringConstArrayView u_dirichlet_string = bs->getValue();
      for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
        if (u_dirichlet_string[dof_index] != "NULL") {
          Real value = std::stod(u_dirichlet_string[dof_index].localstr());
          if (bs->getEnforceDirichletMethod() == "Penalty") {
            Real penalty = bs->getPenalty();
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, node_dof, m_linear_system, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowColumnElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
          }
          else {
            ARCANE_FATAL("Unknown Dirichlet method");
          }
        }
      }
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies Point Dirichlet boundary conditions to RHS and LHS.
     *
     * Updates the LHS matrix and RHS vector to enforce the Dirichlet.
     *
     * - For LHS matrix `𝐀`, the diagonal term for the Dirichlet DOF is set to `𝑃`.
     * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
     *
     * @param bs Boundary condition values.
     * @param node_dof DOF connectivity view.
     * @param node_coord Node coordinates.
     * @param m_linear_system Linear system for LHS.
     * @param rhs_values RHS values to update.
     */
    /*---------------------------------------------------------------------------*/
    static inline void applyPointDirichletToLhsAndRhs(BC::IDirichletPointCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, DoFLinearSystem& m_linear_system, VariableDoFReal& rhs_values)
    {
      NodeGroup node_group = bs->getNode();
      const StringConstArrayView u_dirichlet_string = bs->getValue();
      for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
        if (u_dirichlet_string[dof_index] != "NULL") {
          Real value = std::stod(u_dirichlet_string[dof_index].localstr());
          if (bs->getEnforceDirichletMethod() == "Penalty") {
            Real penalty = bs->getPenalty();
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaPenalty(dof_index, value, penalty, node_dof, m_linear_system, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupViaRowColumnElimination(dof_index, value, node_dof, m_linear_system, rhs_values, node_group);
          }
          else {
            ARCANE_FATAL("Unknown Dirichlet method");
          }
        }
      }
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies Dirichlet boundary conditions to RHS.
     *
     * Updates the RHS vector to enforce Dirichlet conditions.
     *
     * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
     *
     * @param bs Boundary condition values.
     * @param node_dof DOF connectivity view.
     * @param node_coord Node coordinates.
     * @param rhs_values RHS values to update.
     */
    /*---------------------------------------------------------------------------*/
    static inline void applyDirichletToRhs(BC::IDirichletBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, VariableDoFReal& rhs_values)
    {
      FaceGroup face_group = bs->getSurface();
      NodeGroup node_group = face_group.nodeGroup();
      const StringConstArrayView u_dirichlet_string = bs->getValue();
      for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
        if (u_dirichlet_string[dof_index] != "NULL") {
          Real value = std::stod(u_dirichlet_string[dof_index].localstr());
          if (bs->getEnforceDirichletMethod() == "Penalty") {
            Real penalty = bs->getPenalty();
            value = value * penalty;
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
          }
          else {
            ARCANE_FATAL("Unknown Dirichlet method");
          }
        }
      }
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Applies Point Dirichlet boundary conditions to RHS.
     *
     * Updates the RHS vector to enforce the Dirichlet.
     *
     * - For RHS vector `𝐛`, the Dirichlet DOF term is scaled by `𝑃`.
     *
     * @param bs Boundary condition values.
     * @param node_dof DOF connectivity view.
     * @param node_coord Node coordinates.
     * @param rhs_values RHS values to update.
     */
    /*---------------------------------------------------------------------------*/
    static inline void applyPointDirichletToRhs(BC::IDirichletPointCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, VariableDoFReal& rhs_values)
    {
      NodeGroup node_group = bs->getNode();
      const StringConstArrayView u_dirichlet_string = bs->getValue();
      for (Int32 dof_index = 0; dof_index < u_dirichlet_string.size(); ++dof_index) {
        if (u_dirichlet_string[dof_index] != "NULL") {
          Real value = std::stod(u_dirichlet_string[dof_index].localstr());
          if (bs->getEnforceDirichletMethod() == "Penalty") {
            Real penalty = bs->getPenalty();
            value = value * penalty;
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
          }
          else if (bs->getEnforceDirichletMethod() == "RowColumnElimination") {
            ArcaneFemFunctions::BoundaryConditionsHelpers::applyDirichletToNodeGroupRhsOnly(dof_index, value, node_dof, rhs_values, node_group);
          }
          else {
            ARCANE_FATAL("Unknown Dirichlet method");
          }
        }
      }
    }
  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for applying boundary conditions in 3D FEM problems.
   *
   * This class includes static methods for applying Neumann boundary conditions
   * to the right-hand side (RHS) of finite element method equations in 3D.
   */
  /*---------------------------------------------------------------------------*/
  class BoundaryConditions3D
  {
   public:

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
    /*---------------------------------------------------------------------------*/

    static inline void applyConstantSourceToRhs(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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

    static inline void applyConstantSourceToRhsHexa8(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
              //   𝐍 = [𝑁₁  𝑁₂  𝑁₃  𝑁₄  𝑁₅  𝑁₆  𝑁₇  𝑁₈]
              //   𝑁₁ = 1/8 * (1 - ξ) * (1 - η) * (1 - ζ)
              //   𝑁₂ = 1/8 * (1 + ξ) * (1 - η) * (1 - ζ)
              //   𝑁₃ = 1/8 * (1 + ξ) * (1 + η) * (1 - ζ)
              //   𝑁₄ = 1/8 * (1 - ξ) * (1 + η) * (1 - ζ)
              //   𝑁₅ = 1/8 * (1 - ξ) * (1 - η) * (1 + ζ)
              //   𝑁₆ = 1/8 * (1 + ξ) * (1 - η) * (1 + ζ)
              //   𝑁₇ = 1/8 * (1 + ξ) * (1 + η) * (1 + ζ)
              //   𝑁₈ = 1/8 * (1 - ξ) * (1 + η) * (1 + ζ)
              Real N[8];
              N[0] = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
              N[1] = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
              N[2] = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
              N[3] = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
              N[4] = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
              N[5] = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
              N[6] = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
              N[7] = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);

              // Shape function derivatives in reference space
              //  ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ  ∂𝑁₅/∂ξ  ∂𝑁₆/∂ξ  ∂𝑁₇/∂ξ  ∂𝑁₈/∂ξ ]
              //  ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η  ∂𝑁₅/∂η  ∂𝑁₆/∂η  ∂𝑁₇/∂η  ∂𝑁₈/∂η ]
              //  ∂𝐍/∂ζ = [ ∂𝑁₁/∂ζ  ∂𝑁₂/∂ζ  ∂𝑁₃/∂ζ  ∂𝑁₄/∂ζ  ∂𝑁₅/∂ζ  ∂𝑁₆/∂ζ  ∂𝑁₇/∂ζ  ∂𝑁₈/∂ζ ]
              Real dN_dxi[8], dN_deta[8], dN_dzeta[8];
              dN_dxi[0] = -0.125 * (1 - eta) * (1 - zeta);
              dN_dxi[1] = 0.125 * (1 - eta) * (1 - zeta);
              dN_dxi[2] = 0.125 * (1 + eta) * (1 - zeta);
              dN_dxi[3] = -0.125 * (1 + eta) * (1 - zeta);
              dN_dxi[4] = -0.125 * (1 - eta) * (1 + zeta);
              dN_dxi[5] = 0.125 * (1 - eta) * (1 + zeta);
              dN_dxi[6] = 0.125 * (1 + eta) * (1 + zeta);
              dN_dxi[7] = -0.125 * (1 + eta) * (1 + zeta);

              dN_deta[0] = -0.125 * (1 - xi) * (1 - zeta);
              dN_deta[1] = -0.125 * (1 + xi) * (1 - zeta);
              dN_deta[2] = 0.125 * (1 + xi) * (1 - zeta);
              dN_deta[3] = 0.125 * (1 - xi) * (1 - zeta);
              dN_deta[4] = -0.125 * (1 - xi) * (1 + zeta);
              dN_deta[5] = -0.125 * (1 + xi) * (1 + zeta);
              dN_deta[6] = 0.125 * (1 + xi) * (1 + zeta);
              dN_deta[7] = 0.125 * (1 - xi) * (1 + zeta);

              dN_dzeta[0] = -0.125 * (1 - xi) * (1 - eta);
              dN_dzeta[1] = -0.125 * (1 + xi) * (1 - eta);
              dN_dzeta[2] = -0.125 * (1 + xi) * (1 + eta);
              dN_dzeta[3] = -0.125 * (1 - xi) * (1 + eta);
              dN_dzeta[4] = 0.125 * (1 - xi) * (1 - eta);
              dN_dzeta[5] = 0.125 * (1 + xi) * (1 - eta);
              dN_dzeta[6] = 0.125 * (1 + xi) * (1 + eta);
              dN_dzeta[7] = 0.125 * (1 - xi) * (1 + eta);

              // Jacobian for 3D (using your working stiffness matrix approach)
              Real3x3 J;
              for (Int8 a = 0; a < 8; ++a) {
                const Real3& n = node_coord[cell.nodeId(a)];
                J[0][0] += dN_dxi[a] * n.x; // ∂𝑥/∂ξ
                J[0][1] += dN_dxi[a] * n.y; // ∂𝑦/∂ξ
                J[0][2] += dN_dxi[a] * n.z; // ∂𝑧/∂ξ
                J[1][0] += dN_deta[a] * n.x; // ∂𝑥/∂η
                J[1][1] += dN_deta[a] * n.y; // ∂𝑦/∂η
                J[1][2] += dN_deta[a] * n.z; // ∂𝑧/∂η
                J[2][0] += dN_dzeta[a] * n.x; // ∂𝑥/∂ζ
                J[2][1] += dN_dzeta[a] * n.y; // ∂𝑦/∂ζ
                J[2][2] += dN_dzeta[a] * n.z; // ∂𝑧/∂ζ
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
     * @brief Applies a nodal field as a source term to the RHS vector.
     *
     * @param field The field values at cell nodes.
     * @param mesh The mesh containing all cells.
     * @param node_dof DOF connectivity view.
     * @param node_coord The coordinates of the nodes.
     * @param rhs_values The RHS values to update.
     */
    /*---------------------------------------------------------------------------*/

    static inline void integrateNodalFieldToRhsTetra4(VariableNodeReal& field, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
    /**
     * @brief Applies a nodal field to the RHS vector.
     *
     * @param field The field term defined on nodes.
     * @param mesh The mesh containing all cells.
     * @param node_dof DOF connectivity view.
     * @param node_coord The coordinates of the nodes.
     * @param rhs_values The RHS values to update.
     */
    /*---------------------------------------------------------------------------*/

    static inline void integrateNodalFieldToRhsHexa8(VariableNodeReal& field, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
              //   𝐍 = [𝑁₁  𝑁₂  𝑁₃  𝑁₄  𝑁₅  𝑁₆  𝑁₇  𝑁₈]
              //   𝑁₁ = 1/8 * (1 - ξ) * (1 - η) * (1 - ζ)
              //   𝑁₂ = 1/8 * (1 + ξ) * (1 - η) * (1 - ζ)
              //   𝑁₃ = 1/8 * (1 + ξ) * (1 + η) * (1 - ζ)
              //   𝑁₄ = 1/8 * (1 - ξ) * (1 + η) * (1 - ζ)
              //   𝑁₅ = 1/8 * (1 - ξ) * (1 - η) * (1 + ζ)
              //   𝑁₆ = 1/8 * (1 + ξ) * (1 - η) * (1 + ζ)
              //   𝑁₇ = 1/8 * (1 + ξ) * (1 + η) * (1 + ζ)
              //   𝑁₈ = 1/8 * (1 - ξ) * (1 + η) * (1 + ζ)
              Real N[8];
              N[0] = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
              N[1] = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
              N[2] = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
              N[3] = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
              N[4] = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
              N[5] = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
              N[6] = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
              N[7] = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);

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
    /*---------------------------------------------------------------------------*/

    static inline void applyManufacturedSourceToRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_source, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
    /*---------------------------------------------------------------------------*/

    static inline void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
    {
      FaceGroup group = bs->getSurface();

      Real value = 0.0;
      Real valueX = 0.0;
      Real valueY = 0.0;
      Real valueZ = 0.0;

      bool scalarNeumann = false;
      const StringConstArrayView neumann_str = bs->getValue();

      if(neumann_str.size() == 1 && neumann_str[0] != "NULL") {
        scalarNeumann = true;
        value =  std::stod(neumann_str[0].localstr());
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

    static inline void applyNeumannToRhsHexa8(BC::INeumannBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
            Real N[4];
            N[0] = 0.25 * (1 - xi) * (1 - eta);
            N[1] = 0.25 * (1 + xi) * (1 - eta);
            N[2] = 0.25 * (1 + xi) * (1 + eta);
            N[3] = 0.25 * (1 - xi) * (1 + eta);

            // Shape function derivatives w.r.t. natural coordinates
            Real dN_dxi[4], dN_deta[4];
            dN_dxi[0] = -0.25 * (1 - eta);
            dN_dxi[1] = 0.25 * (1 - eta);
            dN_dxi[2] = 0.25 * (1 + eta);
            dN_dxi[3] = -0.25 * (1 + eta);

            dN_deta[0] = -0.25 * (1 - xi);
            dN_deta[1] = -0.25 * (1 + xi);
            dN_deta[2] = 0.25 * (1 + xi);
            dN_deta[3] = 0.25 * (1 - xi);

            // Compute tangent vectors
            Real3 t1(0.0, 0.0, 0.0); // ∂r/∂ξ
            Real3 t2(0.0, 0.0, 0.0); // ∂r/∂η

            for (Int32 i = 0; i < 4; ++i) {
              t1.x += dN_dxi[i] * coords[i].x;
              t1.y += dN_dxi[i] * coords[i].y;
              t1.z += dN_dxi[i] * coords[i].z;

              t2.x += dN_deta[i] * coords[i].x;
              t2.y += dN_deta[i] * coords[i].y;
              t2.z += dN_deta[i] * coords[i].z;
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
    /*---------------------------------------------------------------------------*/

    static inline void applyTractionToRhsTetra4(BC::ITractionBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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

    static inline void applyTractionToRhsHexa8(BC::ITractionBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
            Real N[4];
            N[0] = 0.25 * (1 - xi) * (1 - eta);
            N[1] = 0.25 * (1 + xi) * (1 - eta);
            N[2] = 0.25 * (1 + xi) * (1 + eta);
            N[3] = 0.25 * (1 - xi) * (1 + eta);

            // Shape function derivatives w.r.t. natural coordinates
            Real dN_dxi[4], dN_deta[4];
            dN_dxi[0] = -0.25 * (1 - eta);
            dN_dxi[1] = 0.25 * (1 - eta);
            dN_dxi[2] = 0.25 * (1 + eta);
            dN_dxi[3] = -0.25 * (1 + eta);

            dN_deta[0] = -0.25 * (1 - xi);
            dN_deta[1] = -0.25 * (1 + xi);
            dN_deta[2] = 0.25 * (1 + xi);
            dN_deta[3] = 0.25 * (1 - xi);

            // Compute tangent vectors
            Real3 t1(0.0, 0.0, 0.0); // ∂r/∂ξ
            Real3 t2(0.0, 0.0, 0.0); // ∂r/∂η

            for (Int32 i = 0; i < 4; ++i) {
              t1.x += dN_dxi[i] * coords[i].x;
              t1.y += dN_dxi[i] * coords[i].y;
              t1.z += dN_dxi[i] * coords[i].z;

              t2.x += dN_deta[i] * coords[i].x;
              t2.y += dN_deta[i] * coords[i].y;
              t2.z += dN_deta[i] * coords[i].z;
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
    /*---------------------------------------------------------------------------*/

    static inline void applyTractionTableToRhsTetra4(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index, const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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

    static inline void applyTractionTableToRhsHexa8(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index, const ConstArrayView<CaseTableInfo>& traction_case_table_list, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
            Real N[4];
            N[0] = 0.25 * (1 - xi) * (1 - eta);
            N[1] = 0.25 * (1 + xi) * (1 - eta);
            N[2] = 0.25 * (1 + xi) * (1 + eta);
            N[3] = 0.25 * (1 - xi) * (1 + eta);

            // Shape function derivatives w.r.t. natural coordinates
            Real dN_dxi[4], dN_deta[4];
            dN_dxi[0] = -0.25 * (1 - eta);
            dN_dxi[1] = 0.25 * (1 - eta);
            dN_dxi[2] = 0.25 * (1 + eta);
            dN_dxi[3] = -0.25 * (1 + eta);

            dN_deta[0] = -0.25 * (1 - xi);
            dN_deta[1] = -0.25 * (1 + xi);
            dN_deta[2] = 0.25 * (1 + xi);
            dN_deta[3] = 0.25 * (1 - xi);

            // Compute tangent vectors
            Real3 t1(0.0, 0.0, 0.0); // ∂r/∂ξ
            Real3 t2(0.0, 0.0, 0.0); // ∂r/∂η

            for (Int32 i = 0; i < 4; ++i) {
              t1.x += dN_dxi[i] * coords[i].x;
              t1.y += dN_dxi[i] * coords[i].y;
              t1.z += dN_dxi[i] * coords[i].z;

              t2.x += dN_deta[i] * coords[i].x;
              t2.y += dN_deta[i] * coords[i].y;
              t2.z += dN_deta[i] * coords[i].z;
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
    /*---------------------------------------------------------------------------*/
    static inline void applyManufacturedDirichletToLhsAndRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_dirichlet, Real /*lambda*/, const FaceGroup& group, BC::IManufacturedSolution* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, DoFLinearSystem& m_linear_system, VariableDoFReal& rhs_values)
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
  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for applying boundary conditions in 3D FEM problems.
   *
   * This class includes static methods for applying Neumann boundary conditions
   * to the right-hand side (RHS) of finite element method equations in 3D.
   */
  /*---------------------------------------------------------------------------*/
  class BoundaryConditions2D
  {
   public:

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
    /*---------------------------------------------------------------------------*/

    static inline void applyConstantSourceToRhs(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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

    static inline void applyConstantSourceToRhsQuad4(Real qdot, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
            //   𝐍 = [𝑁₁  𝑁₂  𝑁₃  𝑁₄]
            //   𝑁₁ = 1/4 * (1 - ξ) * (1 - η)
            //   𝑁₂ = 1/4 * (1 + ξ) * (1 - η)
            //   𝑁₃ = 1/4 * (1 + ξ) * (1 + η)
            //   𝑁₄ = 1/4 * (1 - ξ) * (1 + η)
            Real N[4];
            N[0] = 0.25 * (1 - xi) * (1 - eta);
            N[1] = 0.25 * (1 + xi) * (1 - eta);
            N[2] = 0.25 * (1 + xi) * (1 + eta);
            N[3] = 0.25 * (1 - xi) * (1 + eta);

            // Shape function derivatives ∂𝐍/∂ξ and ∂𝐍/∂η
            //     ∂𝐍/∂ξ = [ ∂𝑁₁/∂ξ  ∂𝑁₂/∂ξ  ∂𝑁₃/∂ξ  ∂𝑁₄/∂ξ ]
            //     ∂𝐍/∂η = [ ∂𝑁₁/∂η  ∂𝑁₂/∂η  ∂𝑁₃/∂η  ∂𝑁₄/∂η ]
            Real dN_dxi[4] = { -0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta) };
            Real dN_deta[4] = { -0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi) };

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
     * @brief Applies a nodal field to the RHS vector.
     *
     * @param field The field term defined on nodes.
     * @param mesh  The mesh containing all cells.
     * @param node_dof DOF connectivity view.
     * @param node_coord The coordinates of the nodes.
     * @param rhs_values The RHS values to update.
     */
    /*---------------------------------------------------------------------------*/

    static inline void integrateNodalFieldToRhsTria3(VariableNodeReal& field, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
    /**
     * @brief Applies a nodal field to the RHS vector.
     *
     * @param field The field term defined on nodes.
     * @param mesh  The mesh containing all cells.
     * @param node_dof DOF connectivity view.
     * @param node_coord The coordinates of the nodes.
     * @param rhs_values The RHS values to update.
     */
    /*---------------------------------------------------------------------------*/

    static inline void integrateNodalFieldToRhsQuad4(VariableNodeReal& field, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
        Real node_contributions[4] = {0.0, 0.0, 0.0, 0.0};

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
            Real N[4];
            N[0] = 0.25 * (1 - xi) * (1 - eta);
            N[1] = 0.25 * (1 + xi) * (1 - eta);
            N[2] = 0.25 * (1 + xi) * (1 + eta);
            N[3] = 0.25 * (1 - xi) * (1 + eta);

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
    /*---------------------------------------------------------------------------*/

    static inline void applyManufacturedSourceToRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_source, IMesh* mesh, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
    /*---------------------------------------------------------------------------*/

    static inline void applyNeumannToRhs(BC::INeumannBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
    {
      FaceGroup group = bs->getSurface();

      Real value = 0.0;
      Real valueX = 0.0;
      Real valueY = 0.0;

      bool scalarNeumann = false;
      const StringConstArrayView neumann_str = bs->getValue();

      if(neumann_str.size() == 1 && neumann_str[0] != "NULL") {
        scalarNeumann = true;
        value =  std::stod(neumann_str[0].localstr());
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

    static inline void applyNeumannToRhsQuad4(BC::INeumannBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
          Real N[2];
          N[0] = 0.5 * (1 - xi);
          N[1] = 0.5 * (1 + xi);

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
    /**
     * @brief Applies traction conditions to the right-hand side (RHS) values.
     *
     * @param bs The traction boundary condition values.
     * @param node_dof Connectivity view for degrees of freedom at nodes.
     * @param node_coord Coordinates of the nodes in the mesh.
     * @param rhs_values The right-hand side values to be updated.
     */
    /*---------------------------------------------------------------------------*/

    static inline void applyTractionToRhsTria3(BC::ITractionBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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

    static inline void applyTractionToRhsQuad4(BC::ITractionBoundaryCondition* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
        Real2 normal = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, node_coord);

        Node node0 = face.node(0);
        Node node1 = face.node(1);

        for (Int32 i = 0; i < 2; ++i) {
          Real xi = gp[i];
          Real weight = weights[i];

          // Linear shape functions for Line2
          Real N[2];
          N[0] = 0.5 * (1 - xi);
          N[1] = 0.5 * (1 + xi);

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
    /*---------------------------------------------------------------------------*/

    static inline void applyTractionTableToRhsTria3(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index, const UniqueArray<Arcane::FemUtils::CaseTableInfo>& traction_case_table_list, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
            rhs_values[node_dof.dofId(node, 0)] += trac[0] * length / 2.;
            rhs_values[node_dof.dofId(node, 1)] += trac[1] * length / 2.;
          }
        }
      }
    }

    static inline void applyTractionTableToRhsQuad4(BC::ITractionBoundaryCondition* bs, const Real t, Int32 boundary_condition_index, const UniqueArray<CaseTableInfo>& traction_case_table_list, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, VariableDoFReal& rhs_values)
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
        Real2 normal = ArcaneFemFunctions::MeshOperation::computeNormalEdge2(face, node_coord);

        Node node0 = face.node(0);
        Node node1 = face.node(1);

        for (Int32 i = 0; i < 2; ++i) {
          Real xi = gp[i];
          Real weight = weights[i];

          // Linear shape functions for Line2
          Real N[2];
          N[0] = 0.5 * (1 - xi);
          N[1] = 0.5 * (1 + xi);

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
    /*---------------------------------------------------------------------------*/
    static inline void applyManufacturedDirichletToLhsAndRhs(IBinaryMathFunctor<Real, Real3, Real>* manufactured_dirichlet, Real /*lambda*/, const FaceGroup& group, BC::IManufacturedSolution* bs, const IndexedNodeDoFConnectivityView& node_dof, const VariableNodeReal3& node_coord, DoFLinearSystem& m_linear_system, VariableDoFReal& rhs_values)
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
  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods based on the Dispatcher mechanism available in Arcane,
   * allowing to compute FEM methods without declaring the FE entity type
   * (coming from PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  class CellFEMDispatcher
  {

   public:

    Real getShapeFuncVal(Int16 /*item_type*/, Integer /*inod*/, Real3 /*ref coord*/);
    Real3 getShapeFuncDeriv(Int16 /*item_type*/, Integer /*inod*/, Real3 /*ref coord*/);

    RealUniqueArray getGaussData(ItemWithNodes item, Integer nint, Integer ngauss);

    CellFEMDispatcher();

   private:

    std::function<Real(Integer inod, Real3 coord)> m_shapefunc[NB_BASIC_ITEM_TYPE];
    std::function<Real3(Integer inod, Real3 coord)> m_shapefuncderiv[NB_BASIC_ITEM_TYPE];
  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for various FEM-related operations.
   *
   * This class includes static methods for computing shape functions and their
   * derivatives depending on finite element types. These methods are used within
   * the Dispatcher mechanism available in Arcane through the class
   * CellFEMDispatcher (coming from PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  class FemShapeMethods
  {
   public:

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) edge finite-element
     * The "Line2" reference element is assumed as follows:
     *  0           1
     *  o-----------o---> x
     * -1           1
     * direct local numbering : 0->1
     */
    /*---------------------------------------------------------------------------*/
    static inline Real line2ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 2);
#endif

      Real r = ref_coord[0];
      if (inod == 1)
        return (0.5 * (1 + r));
      return (0.5 * (1 - r));
    }

    static inline Real3 line2ShapeFuncDeriv(Integer inod, Real3)
    {
      if (inod == 1)
        return { 0.5, 0., 0. };
      return { -0.5, 0., 0. };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) edge finite-element
     * The "Line3" reference element is assumed as follows:
     *  0     2      1
     *  o-----o------o---> x
     * -1     0      1
     * direct local numbering : 0->1->2
     */
    /*---------------------------------------------------------------------------*/
    static inline Real line3ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 3);
#endif

      Real ri = ref_coord[0];
      if (inod == 0)
        ri *= -1;

      if (inod < 2)
        return 0.5 * ri * (1 + ri); // nodes 0 or 1
      return (1 - ri * ri); // middle node
    }

    static inline Real3 line3ShapeFuncDeriv(Integer inod, Real3 ref_coord)
    {
      Real ri = ref_coord[0];
      if (!inod)
        return { -0.5 + ri, 0., 0. };
      if (inod == 1)
        return { 0.5 + ri, 0., 0. };
      return { -2. * ri, 0., 0. };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) triangle finite-element
     * The "Tri3" reference element is assumed as follows:
     *
     *   ^ s
     *   |
     *  2 (1,0)
     *   o
     *   . .
     *   .   .
     *   .     .
     *   .       .
     *   .         .
     *   .           .
     *   o-------------o---------> r
     *  0 (0,0)         1 (1,0)
     * direct local numbering : 0->1->2
     */
    /*---------------------------------------------------------------------------*/
    static inline Real tri3ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 3);
#endif
      Real r = ref_coord[0];
      Real s = ref_coord[1];
      if (!inod)
        return (1 - r - s);
      if (inod == 1)
        return r;
      return s;
    }

    static inline Real3 tri3ShapeFuncDeriv(Integer inod, Real3)
    {
      if (!inod)
        return { -1., -1., 0. };
      if (inod == 1)
        return { 1., 0., 0. };
      return { 0., 1., 0. };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) triangle finite-element
     * The "Tri6" reference element is assumed as follows:
     *   ^ s
     *   |
     *  2 (1,0)
     *   o
     *   .  .
     *   .    .
     *   .      .
     *   o 6      o 5(0.5;0.5)
     *   .(0;0.5)   .
     *   .            .
     *   .              .
     *   o-------o-------o---------> r
     * 0(0,0)  4(0.5;0)  1(1,0)
     * direct local numbering : 0->1->2->3->4->5
     */
    /*---------------------------------------------------------------------------*/
    static inline Real tri6ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 6);
#endif
      auto wi = 0., ri = ref_coord[0], si = ref_coord[1];
      auto ri2 = 2. * ri - 1.;
      auto si2 = 2. * si - 1.;
      auto ti = 1. - ri - si, ti2 = 2. * ti - 1.;

      switch (inod) {
      default:
        break;
      case 0:
        wi = ti * ti2;
        break;
      case 1:
        wi = ri * ri2;
        break;
      case 2:
        wi = si * si2;
        break;
      case 3:
        wi = 4. * ri * ti;
        break;
      case 4:
        wi = 4. * ri * si;
        break;
      case 5:
        wi = 4. * si * ti;
        break;
      }
      return wi;
    }

    static inline Real3 tri6ShapeFuncDeriv(Integer inod, Real3 ref_coord)
    {
      auto ri = ref_coord[0], si = ref_coord[1];
      auto ti = 1. - ri - si;

      if (!inod) {
        auto wi = -3. + 4. * (ri + si);
        return { wi, wi, 0. };
      }
      if (inod == 1)
        return { -1. + 4. * ri, 0., 0. };
      if (inod == 2)
        return { 0., -1. + 4. * si, 0. };

      if (inod == 3)
        return { 4. * (ti - ri), -4. * ri, 0. };
      if (inod == 4)
        return { 4. * si, 4. * ri, 0. };
      return { -4. * si, 4. * (ti - si), 0. };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) quadrangle finite-element
     * The "Quad4" reference element is assumed as follows:
     *         ^y
     *          |
     *  1 o-----1-----o 0
     *    |     |     |
     *    |     |     |
     *    |     |     |
     *   -1 ----|---- 1 ---> x
     *    |     |     |
     *    |     |     |
     *    |     |     |
     *  2 o--- -1 ----o 3
     *
     * direct local numbering : 0->1->2->3
     */
    /*---------------------------------------------------------------------------*/
    static inline Real quad4ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 4);
#endif

      auto r{ ref_coord[0] }, s{ ref_coord[1] };
      auto ri{ 1. }, si{ 1. };

      switch (inod) {
      default:
        break; // default is first node (index 0)
      case 2:
        si = -1;
        [[fallthrough]];
      case 1:
        ri = -1;
        break;

      case 3:
        si = -1;
        break;
      }
      return ((1 + ri * r) * (1 + si * s) / 4.);
    }

    static inline Real3 quad4ShapeFuncDeriv(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 4);
#endif

      auto r{ ref_coord[0] }, s{ ref_coord[1] };
      auto ri{ 1. }, si{ 1. }; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

      switch (inod) {
      default:
        break; // default is first node (index 0)
      case 2:
        si = -1;
        [[fallthrough]];
      case 1:
        ri = -1;
        break;

      case 3:
        si = -1;
        break;
      }
      return { 0.25 * ri * (1 + si * s), 0.25 * si * (1 + ri * r), 0. };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) quadrangle finite-element
     * The "Quad8" reference element is assumed as follows:
     *         ^y
     *          |
     *          4
     *  1 o-----o-----o 0
     *    |     |     |
     *    |     |     |
     *    |     |     |
     *  5 o ----|---- o 7 ---> x
     *    |     |     |
     *    |     |     |
     *    |     |     |
     *  2 o-----o-----o 3
     *          6
     * Normalized coordinates (x, y) vary between -1/+1
     * Nodes 4, 6 are on line (x = 0)
     * direct local numbering :  0->1->2->...->5->6->7
     */
    /*---------------------------------------------------------------------------*/
    static inline Real quad8ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 8);
#endif
      Real tol{1.0e-15};

      auto r{ ref_coord[0] }, s{ ref_coord[1] };
      auto ri{ 1. }, si{ 1. };

      switch (inod) {
      default:
        break; // default is first node (index 0)
      case 2:
        si = -1;
        [[fallthrough]];
      case 1:
        ri = -1;
        break;

      case 3:
        si = -1;
        break;

      case 6:
        si = -1;
        [[fallthrough]];
      case 4:
        ri = 0;
        break;

      case 5:
        ri = -1;
        [[fallthrough]];
      case 7:
        si = 0;
        break;
      }

      auto r0{ r * ri }, s0{ s * si };
      Real Phi{ 0. };
      auto t0{ r0 + s0 - 1. };

      if (inod < 4) // Corner nodes
        Phi = (1 + r0) * (1 + s0) * t0 / 4.;

      else { // Middle nodes
        if (fabs(ri) < tol)
          Phi = (1 - r * r) * (1 + s0) / 2.;
        else if (fabs(si) < tol)
          Phi = (1 - s * s) * (1 + r0) / 2.;
      }
      return Phi;
    }

    static inline Real3 quad8ShapeFuncDeriv(Integer inod, Real3 ref_coord)
    {
      Real tol{1.0e-15};

      auto r{ ref_coord[0] }, s{ ref_coord[1] };
      auto ri{ 1. }, si{ 1. };

      switch (inod) {
      default:
        break; // default is first node (index 0)
      case 2:
        si = -1;
        [[fallthrough]];
      case 1:
        ri = -1;
        break;

      case 3:
        si = -1;
        break;

      case 6:
        si = -1;
        [[fallthrough]];
      case 4:
        ri = 0;
        break;

      case 5:
        ri = -1;
        [[fallthrough]];
      case 7:
        si = 0;
        break;
      }

      auto r0{ r * ri }, s0{ s * si };
      Real3 dPhi;
      auto t0{ r0 + s0 - 1. };

      if (inod < 4) { // Corner nodes
        dPhi.x = ri * (1 + s0) * (t0 + 1. + r0) / 4.;
        dPhi.y = si * (1 + r0) * (t0 + 1. + s0) / 4.;
      }
      else { // Middle nodes
        if (fabs(ri) < tol) {
          dPhi.x = -r * (1 + s0);
          dPhi.y = si * (1 - r * r) / 2.;
        }
        else if (fabs(si) < tol) {
          dPhi.x = -s * (1 + r0);
          dPhi.y = ri * (1 - s * s) / 2.;
        }
      }
      dPhi.z = 0.;
      return dPhi;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) hexaedron finite-element
     * The "Hexa8" reference element is assumed as follows:
     *     (-1, 1,1)
     *         1-------------0 (1,1,1)
     *        /|            /|
     *       / |           / |
     *     /   |          /  |
     *    2----|---------3   |   z   y
     *  (-1,-1,1)        |   |   | /
     *    |    |         |   |   |/--->x
     *    |    |         |   |
     *    |    |         |   |
     *    |    5---------|---4 (1,1,-1)
     *    |  /           |  /
     *    | /            | /
     *    |/             |/
     *    6--------------7 (1,-1,-1)
     * (-1,-1,-1)
     * Normalized coordinates (x, y, z) vary between -1/+1
     * direct local numbering : 0->1->2->3->4->5->6->7
     */
    /*---------------------------------------------------------------------------*/
    static inline Real hexa8ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 8);
#endif
      auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
      auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

      switch (inod) {
      default:
        break;
      case 3:
      case 2:
        ri = -1;
        break;
      case 0:
      case 1:
        ri = -1;
        si = -1;
        break;
      case 4:
      case 5:
        si = -1;
        break;
      }
      if (inod == 1 || inod == 2 || inod == 5 || inod == 6)
        ti = -1;

      auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
      auto Phi = (1 + r0) * (1 + s0) * (1 + t0) / 8.;

      return Phi;
    }

    static inline Real3 hexa8ShapeFuncDeriv(Integer inod, Real3 ref_coord)
    {

      auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
      auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 7 = (1,1,1)

      switch (inod) {
      default:
        break;
      case 3:
      case 2:
        ri = -1;
        break;
      case 0:
      case 1:
        ri = -1;
        si = -1;
        break;
      case 4:
      case 5:
        si = -1;
        break;
      }
      if (inod == 1 || inod == 2 || inod == 5 || inod == 6)
        ti = -1;

      auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
      Real3 dPhi;
      dPhi.x = ri * (1 + s0) * (1 + t0) / 8.;
      dPhi.y = si * (1 + r0) * (1 + t0) / 8.;
      dPhi.z = ti * (1 + r0) * (1 + s0) / 8.;
      return dPhi;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) hexaedron finite-element
     * The "Hexa20" reference element is assumed as follows:
     *     (-1, 1,1)
     *         1------8------0 (1,1,1)
     *        /|            /|
     *      9  |          11 |
     *     /   |          /  |
     *    2----|--10-----3   |   z   y
     *  (-1,-1,1)        |   |   | /
     *    |   16         |  19   |/--->x
     *    |    |         |   |
     *   17    |        18   |
     *    |    5----12---|---4 (1,1,-1)
     *    |  /           |  /
     *    | 13           | 15
     *    |/             |/
     *    6-----14-------7 (1,-1,-1)
     * (-1,-1,-1)
     * Normalized coordinates (x, y, z) vary between -1/+1
     * direct local numbering : 0->1->2->3->...->18->19
     */
    /*---------------------------------------------------------------------------*/

    static inline Real hexa20ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 20);
#endif
      Real tol{1.0e-15};

      auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
      auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 0 = (1,1,1)

      switch (inod) {
      default:
        break;

      case 5:
        ti = -1.;
        [[fallthrough]];
      case 1:
        ri = -1;
        break;

      case 6:
        ti = -1.;
        [[fallthrough]];
      case 2:
        ri = -1;
        si = -1;
        break;

      case 7:
        ti = -1.;
        [[fallthrough]];
      case 3:
        si = -1;
        break;

      case 4:
        ti = -1.;
        break;

      case 9:
        ri = -1.;
        [[fallthrough]];
      case 11:
        si = 0.;
        break;

      case 10:
        si = -1.;
        [[fallthrough]];
      case 8:
        ri = 0.;
        break;

      case 14:
        si = -1.;
        [[fallthrough]];
      case 12:
        ri = 0.;
        ti = -1.;
        break;

      case 17:
        si = -1.;
        [[fallthrough]];
      case 16:
        ri = -1.;
        ti = 0.;
        break;

      case 18:
        si = -1.;
        [[fallthrough]];
      case 19:
        ti = 0.;
        break;
      }

      auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
      Real Phi{ 0. };
      auto t{ r0 + s0 + t0 - 2. };

      if (inod < 8) // Corner nodes
        Phi = (1 + r0) * (1 + s0) * (1 + t0) * t / 8.;

      else { // Middle nodes
        if (math::abs(ri) < tol)
          Phi = (1 - x * x) * (1 + s0) * (1 + t0) / 4.;
        else if (math::abs(si) < tol)
          Phi = (1 - y * y) * (1 + r0) * (1 + t0) / 4.;
        else if (math::abs(ti) < tol)
          Phi = (1 - z * z) * (1 + r0) * (1 + s0) / 4.;
      }
      return Phi;
    }

    static inline Real3 hexa20ShapeFuncDeriv(Integer inod, Real3 ref_coord)
    {
      Real tol{1.0e-15};

      auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] };
      auto ri{ 1. }, si{ 1. }, ti{ 1. }; // Normalized coordinates (=+-1) =>node index 0 = (1,1,1)

      switch (inod) {
      default:
        break;

      case 5:
        ti = -1.;
      case 1:
        ri = -1;
        break;

      case 6:
        ti = -1.;
      case 2:
        ri = -1;
        si = -1;
        break;

      case 7:
        ti = -1.;
      case 3:
        si = -1;
        break;

      case 4:
        ti = -1.;
        break;

      case 9:
        ri = -1.;
      case 11:
        si = 0.;
        break;

      case 10:
        si = -1.;
      case 8:
        ri = 0.;
        break;

      case 14:
        si = -1.;
      case 12:
        ri = 0.;
        ti = -1.;
        break;

      case 17:
        si = -1.;
      case 16:
        ri = -1.;
        ti = 0.;
        break;

      case 18:
        si = -1.;
      case 19:
        ti = 0.;
        break;
      }

      auto r0{ x * ri }, s0{ y * si }, t0{ z * ti };
      auto t{ r0 + s0 + t0 - 2. };
      Real3 dPhi;

      if (inod < 8) { // Corner nodes
        dPhi = hexa8ShapeFuncDeriv(inod, ref_coord);
        dPhi.x *= (t + 1. + r0);
        dPhi.y *= (t + 1. + s0);
        dPhi.z *= (t + 1. + t0);
      }
      else { // Middle nodes
        auto x2{ x * x }, y2{ y * y }, z2{ z * z };
        if (math::abs(ri) < tol) {
          dPhi.x = -x * (1 + s0) * (1 + t0) / 2.;
          dPhi.y = si * (1 - x2) * (1 + t0) / 4.;
          dPhi.z = ti * (1 - x2) * (1 + s0) / 4.;
        }
        else if (math::abs(si) < tol) {
          dPhi.x = ri * (1 - y2) * (1 + t0) / 4.;
          dPhi.y = -y * (1 + r0) * (1 + t0) / 2.;
          dPhi.z = ti * (1 - y2) * (1 + r0) / 4.;
        }
        else if (math::abs(ti) < tol) {
          dPhi.x = ri * (1 - z2) * (1 + s0) / 4.;
          dPhi.y = si * (1 - z2) * (1 + r0) / 4.;
          dPhi.z = -z * (1 + r0) * (1 + s0) / 2.;
        }
      }
      return dPhi;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) tetrahedral finite-element
     * The "Tetra4" reference element is assumed as follows:
     *
     *    (0,0,1)                     3
     *       .                        *.*
     *       .                        * . *
     *       .                        *  .  *
     *       .                        *   .   *
     *       Z   (0,1,0)              *    .    *
     *       .    .                   *     2     *
     *       .   .                    *   .    .    *
     *       .  Y                     *  .        .   *
     *       . .                      * .            .  *
     *       ..           (1,0,0)     *.                . *
     *       --------X------>         0********************1
     *
     * direct local numbering : 0->1->2->3
     */
    /*---------------------------------------------------------------------------*/

    static inline Real tetra4ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 4);
#endif

      auto ri = ref_coord[0], si = ref_coord[1], ti = ref_coord[2]; // default is first node (index 3)

      switch (inod) {
      default:
        break;
      case 1:
        return ri;
      case 2:
        return si;
      case 0:
        return (1. - ri - si - ti);
      }
      return ti;
    }

    static inline Real3 tetra4ShapeFuncDeriv(Integer inod, Real3 /*ref_coord*/)
    {

      if (inod == 3)
        return { 0., 0., 1. };
      if (inod == 1)
        return { 1., 0., 0. };
      if (inod == 2)
        return { 0., 1., 0. };
      return { -1., -1., -1. };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference quadratic (P2) tetrahedral finite-element
     * The "Tetra10" reference element is assumed as follows:
     *
     *    (0,0,1)                     x 3
     *       .                        *.*
     *       .                        * . *
     *       .                        *  .  *
     *       .                        *  9   *
     *       Z   (0,1,0)              *    .    *
     *       .    .                   *     x     8
     *       .   .                    7   . 2  .    *
     *       .  Y                     *  6       5    *
     *       . .                      * .            .  *
     *       ..           (1,0,0)     *.                . *
     *       --------X------>       0 x ******* 4 ******** x 1
     *
     * direct local numbering : 0->1->2->3->...->8->9
     */
    /*---------------------------------------------------------------------------*/

    static inline Real tetra10ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 10);
#endif

      auto x = ref_coord[0], y = ref_coord[1], z = ref_coord[2],
           t = 1. - x - y - z,
           wi{ 0. };

      switch (inod) {
      default:
        break;

      // Corner nodes
      case 0:
        wi = t * (2 * t - 1.);
        break; //=(1. - 2*x - 2*y - 2*z) * t
      case 1:
        wi = x * (2 * x - 1.);
        break; //=(1. - 2*t - 2*y - 2*z)*x
      case 2:
        wi = y * (2 * y - 1.);
        break; //=(1. - 2*x - 2*t - 2*z)*y
      case 3:
        wi = z * (2 * z - 1.);
        break; //=(1. - 2*t - 2*x - 2*y)*z

      // Middle nodes
      case 4:
        wi = 4 * x * t;
        break;
      case 5:
        wi = 4 * x * y;
        break;
      case 6:
        wi = 4 * y * t;
        break;
      case 7:
        wi = 4 * z * t;
        break;
      case 8:
        wi = 4 * z * x;
        break;
      case 9:
        wi = 4 * z * y;
        break;
      }
      return wi;
    }

    static inline Real3 tetra10ShapeFuncDeriv(Integer inod, Real3 ref_coord)
    {
      auto x{ ref_coord[0] }, y{ ref_coord[1] }, z{ ref_coord[2] },
      t{ 1. - x - y - z },
      x4{ 4 * x },
      y4{ 4 * y },
      z4{ 4 * z },
      t4{ 4 * t };

      // Corner nodes
      /*
   if (inod == 3) return {0.,0.,1. + 2*t - 2*x - 2*y + 2*z};
   if (inod == 1) return {1. - 2*t - 2*y - 2*z + 2*x,0.,0.};
   if (inod == 2) return {0.,1. - 2*x - 2*t - 2*z + 2*y,0.};
   if (!inod) return {-1. - 2*t + 2*x + 2*y + 2*z,-1. - 2*t+ 2*x + 2*y + 2*z,-1. - 2*t + 2*x + 2*y + 2*z};
*/
      if (!inod)
        return { 1. - t4, 1. - t4, 1. - t4 };
      if (inod == 1)
        return { x4 - 1., 0., 0. };
      if (inod == 2)
        return { 0., y4 - 1., 0. };
      if (inod == 3)
        return { 0., 0., z4 - 1. };

      // Middle nodes
      if (inod == 4)
        return { t4 - x4, -x4, -x4 };
      if (inod == 5)
        return { y4, x4, 0. };
      if (inod == 6)
        return { -y4, t4 - y4, -y4 };
      if (inod == 8)
        return { z4, 0., x4 };
      if (inod == 9)
        return { 0., z4, y4 };
      return { -z4, -z4, t4 - z4 }; //inod == 7
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) pentaedron finite-element
     * The "Penta6" reference element is assumed as follows:
     *
     *                     5 (0,1,1)
     *                   . |  .
     *                  .  |     .
     *                 .   Z        .
     *                .    |           .
     *               .     |             .
     *       (0,0,1) 3 ------------------ 4 (1,0,1)
     *               |     |              |
     *               |     |              |
     *               |     |              |
     *               |     |              |
     *               |     2 (0,1,-1)     |
     *               |   .    .           |
     *               |  Y        .        |
     *               | .            .     |
     *               |.                .  |
     *      (0,0,-1) 0 -------- X ------- 1 (1,0,-1)
     *
     * direct local numbering : 0->1->2->3->4->5->6
     */
    /*---------------------------------------------------------------------------*/

    static inline Real penta6ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 6);
#endif
      auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
      auto r0{ 1. }, s0{ 1. }, ti{ -1. };
      auto rs{ 1. - r - s };

      if (inod >= 3)
        ti = 1.;
      auto t0{ 1 + ti * t };

      switch (inod) {
      default:
        break; // Node 0
      case 4:
      case 1:
        r0 = r;
        rs = 1.;
        break;
      case 5:
      case 2:
        s0 = s;
        rs = 1.;
        break;
      }

      return 0.5 * r0 * s0 * rs * t0;
    }

    static inline Real3 penta6ShapeFuncDeriv(Integer inod, Real3 ref_coord)
    {

#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 6);
#endif
      auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
      auto ri{ 1. }, si{ 1. };
      auto r0{ 1. }, s0{ 1. }, ti{ -1. };
      auto rs{ 1. - r - s };

      if (inod >= 3)
        ti = 1.;
      auto t0{ 1 + ti * t };

      switch (inod) {
      default:
        break;
      case 3:
      case 0:
        ri = -1.;
        si = -1.;
        break;
      case 4:
      case 1:
        r0 = r;
        si = 0.;
        rs = 1.;
        break;
      case 5:
      case 2:
        s0 = s;
        rs = 1.;
        break;
      }

      Real3 dPhi;
      dPhi.x = 0.5 * ri * t0;
      dPhi.y = 0.5 * si * t0;
      dPhi.z = 0.5 * ti * rs * r0 * s0;
      return dPhi;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides methods for reference linear (P1) pyramid finite-element
     * The "Pyramid5" reference element is assumed as follows:
     *
     *                               ^
     *                               |
     *                               Z
     *                               |
     *                               4 (0,0,1)
     *                              *
     *                             * **
     *                            ** *  *
     *                           * * |*   *
     *                          * *  | *    *
     *                         * *   |  *     *
     *                        * *    |   *      *      .Y
     *                       * *     |    *       *  .
     *            (-1,0,0)  * 2 -----|-----*------ 1 (0,1,0)
     *                     * .       |      *  .  .
     *                    * .        |     . *   .
     *                   *.             .     * .
     *                  *                  X . *
     *        (0,-1,0) 3 --------------------- 0 (1,0,0)
     *
     * direct local numbering : 0->1->2->3->4
     */
    /*---------------------------------------------------------------------------*/

    static inline Real pyramid5ShapeFuncVal(Integer inod, Real3 ref_coord)
    {
#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 5);
#endif
      Real tol{1.0e-15};

      auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
      auto r1{ -1. }, s1{ 1. }, r2{ -1. }, s2{ -1. };

      if (inod == 4)
        return t;
      auto ti{ t - 1. };
      auto t0{ 0. };

      if (math::abs(ti) < tol)
        ti = 0.;
      else
        t0 = -1. / ti / 4.;

      switch (inod) {
      case 1:
        s1 = -1.;
        r2 = 1.;
        break;
      case 2:
        r1 = 1.;
        r2 = 1.;
        break;
      case 3:
        r1 = 1.;
        s2 = 1.;
        break;
      default:
        break; // default is for node 0
      }

      return (r1 * r + s1 * s + ti) * (r2 * r + s2 * s + ti) * t0;
    }

    static inline Real3 pyramid5ShapeFuncDeriv(Integer inod, Real3 ref_coord)
    {

#ifdef _DEBUG
      ARCANE_ASSERT(inod >= 0 && inod < 5);
#endif
      Real tol{1.0e-15};

      auto r{ ref_coord[0] }, s{ ref_coord[1] }, t{ ref_coord[2] };
      auto r1{ -1. }, s1{ 1. }, r2{ -1. }, s2{ -1. };

      auto ti{ t - 1. };
      auto t0{ 0. };

      if (math::abs(ti) < tol)
        ti = 0.;
      else
        t0 = -1. / ti / 4.;

      switch (inod) {
      case 1:
        s1 = -1.;
        r2 = 1.;
        break;
      case 2:
        r1 = 1.;
        r2 = 1.;
        break;
      case 3:
        r1 = 1.;
        s2 = 1.;
        break;
      default:
        break; // default is for node 0
      }

      if (inod == 4)
        return { 0., 0., 1. };

      Real3 dPhi;
      auto r12{ r1 + r2 }, rr{ 2. * r1 * r2 }, s12{ s1 + s2 }, ss{ 2. * s1 * s2 }, rs{ r1 * s2 + r2 * s1 }, t02{ 4. * t0 * t0 };

      dPhi.x = t0 * (rr * r + rs * s + r12 * ti);
      dPhi.y = t0 * (rs * r + ss * s + s12 * ti);

      if (math::abs(ti) < tol)
        dPhi.z = 0.;
      else
        dPhi.z = t0 * (r12 * r + s12 * s + 2. * ti) + t02 * (r1 * r + s1 * s + ti) * (r2 * r + s2 * s + ti);

      return dPhi;
    }

    /*---------------------------------------------------------------------------*/
  };

  /*---------------------------------------------------------------------------*/
  /**
   * @brief Provides methods for Gauss quadrature.
   *
   * This class includes static methods for computing Gauss-Legendre integration
   * depending on finite element types (coming from PASSMO).
   */
  /*---------------------------------------------------------------------------*/
  class FemGaussQuadrature
  {

   public:

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the number of Gauss Points for a given finite element type,
     * depending on the integration order chosen by user (coming fro PASSMO).
     */
    /*---------------------------------------------------------------------------*/
    static inline Integer getNbGaussPointsfromOrder(Int16 cell_type, Integer ninteg)
    {
      Integer nbgauss{ 0 };
      auto ninteg2{ ninteg * ninteg };
      auto ninteg3{ ninteg2 * ninteg };

      if (ninteg <= 1)
        nbgauss = 1;
      else if (cell_type == IT_Line2 || cell_type == IT_Line3)
        nbgauss = ninteg;
      else if (cell_type == IT_Quad4 || cell_type == IT_Quad8)
        nbgauss = ninteg2;
      else if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20)
        nbgauss = ninteg3;

      else if (ninteg == 2) {
        switch (cell_type) {
        default:
          break;

        case IT_Triangle3:
        case IT_Triangle6:
          nbgauss = 3;
          break;

        case IT_Tetraedron4:
        case IT_Tetraedron10:
          nbgauss = 4;
          break;

        case IT_Pentaedron6:
          nbgauss = 6;
          break;

        case IT_Pyramid5:
          nbgauss = 5;
          break;
        }
      }
      else if (ninteg == 3) {
        switch (cell_type) {
        default:
          break;

        case IT_Triangle3:
        case IT_Triangle6:
          nbgauss = 4;
          break;

        case IT_Tetraedron4:
        case IT_Tetraedron10:
          nbgauss = 5;
          break;

        case IT_Pentaedron6:
          nbgauss = 8;
          break;

        case IT_Pyramid5:
          nbgauss = 6;
          break;
        }
      }
      else if (ninteg >= 4) {
        switch (cell_type) {
        default:
          break;

        case IT_Triangle3:
        case IT_Triangle6:
          nbgauss = 7;
          break;

        case IT_Tetraedron4:
        case IT_Tetraedron10:
          nbgauss = 15;
          break;

        case IT_Pentaedron6:
          nbgauss = 21;
          break;

        case IT_Pyramid5:
          nbgauss = 27;
          break;
        }
      }
      return nbgauss;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the position of a Gauss Point in the reference (local) element
     * for a given finite element, depending on the rank of the point in the
     * element loop (coming fro PASSMO).
     */
    /*---------------------------------------------------------------------------*/
    static inline Real3 getGaussRefPosition(ItemWithNodes cell, Integer ninteg, Integer rank)
    {
      auto cell_type = cell.type();
      Integer nint{ ninteg };

      if (nint < 1)
        nint = 1;

      if (cell_type == IT_Line2 || cell_type == IT_Line3)
        return lineRefPosition({ rank, -1, -1 }, { nint, 0, 0 });

      if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
        auto in{ 0 };
        for (Int32 i1 = 0; i1 < nint; ++i1) {
          for (Int32 i2 = 0; i2 < nint; ++i2) {

            if (rank == in)
              return quadRefPosition({ i1, i2, -1 }, { nint, nint, 0 });

            ++in;
          }
        }
      }

      if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
        auto in{ 0 };
        for (Int32 i1 = 0; i1 < nint; ++i1) {
          for (Int32 i2 = 0; i2 < nint; ++i2) {
            for (Int32 i3 = 0; i3 < nint; ++i3) {

              if (rank == in)
                return hexaRefPosition({ i1, i2, i3 }, { nint, nint, nint });

              ++in;
            }
          }
        }
      }

      if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
        auto o{ 3 };
        if (nint <= 3)
          o = nint - 1;

        return { xg1[o][rank], xg2[o][rank], 0. };
      }

      if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10) {
        auto o{ 3 };
        if (nint <= 3)
          o = nint - 1;

        return { xtet[o][rank], ytet[o][rank], ztet[o][rank] };
      }

      if (cell_type == IT_Pyramid5) {
        auto o{ 1 };
        if (nint <= 2)
          o = 0;

        return { xpyr[o][rank], ypyr[o][rank], zpyr[o][rank] };
      }

      if (cell_type == IT_Pentaedron6) {
        auto o{ 1 };
        if (nint <= 2)
          o = 0;

        return { xpent[o][rank], ypent[o][rank], zpent[o][rank] };
      }
      return {};
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the integration weight of a Gauss Point for a given
     * finite element, depending on the rank of the point in the element loop
     * (coming fro PASSMO).
     */
    /*---------------------------------------------------------------------------*/
    static inline Real getGaussWeight(ItemWithNodes cell, Integer ninteg, Integer rank)
    {
      auto cell_type = cell.type();
      Integer nint{ ninteg };

      if (nint < 1)
        nint = 1;

      if (cell_type == IT_Line2 || cell_type == IT_Line3)
        return lineWeight({ rank, -1, -1 }, { nint, 0, 0 });

      if (cell_type == IT_Quad4 || cell_type == IT_Quad8) {
        auto in{ 0 };
        for (Int32 i1 = 0; i1 < nint; ++i1) {
          for (Int32 i2 = 0; i2 < nint; ++i2) {

            if (rank == in)
              return quadWeight({ i1, i2, -1 }, { nint, nint, 0 });

            ++in;
          }
        }
      }

      if (cell_type == IT_Hexaedron8 || cell_type == IT_Hexaedron20) {
        auto in{ 0 };
        for (Int32 i1 = 0; i1 < nint; ++i1) {
          for (Int32 i2 = 0; i2 < nint; ++i2) {
            for (Int32 i3 = 0; i3 < nint; ++i3) {

              if (rank == in)
                return hexaWeight({ i1, i2, i3 }, { nint, nint, nint });

              ++in;
            }
          }
        }
      }

      if (cell_type == IT_Triangle3 || cell_type == IT_Triangle6) {
        auto o{ 3 };
        if (nint <= 3)
          o = nint - 1;

        return wg[o][rank];
      }

      if (cell_type == IT_Tetraedron4 || cell_type == IT_Tetraedron10) {
        if (nint == 1)
          return wgtet1;
        if (nint == 2)
          return wgtet2;
        if (nint == 3) {
          auto i = npwgtet3[rank];
          return wgtet3[i];
        }
        // nint >= 4
        auto i = npwgtet4[rank];
        return wgtet4[i];
      }

      if (cell_type == IT_Pyramid5) {
        if (nint <= 2)
          return wgpyr2;

        // nint >= 3
        auto i = npwgpyr3[rank];
        return wgpyr3[i];
      }

      if (cell_type == IT_Pentaedron6) {
        if (nint <= 2)
          return wgpent2;

        // nint >= 3
        auto i = npwgpent3[rank];
        return wgpent3[i];
      }

      return 1.;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the position of a Gauss Point in the reference (local) element
     * coordinates, depending on the index of the Gauss Point and integration order
     * chosen by user (coming fro PASSMO).
     * This method is generic (not depending on the finite element type) and is called
     * by specialized methods (depending on FE type)
     */
    /*---------------------------------------------------------------------------*/
    static inline Real getRefPosition(Integer indx, Integer ordre)
    {
      Real x = xgauss1; // default is order 1

      switch (ordre) {
      case 2:
        x = xgauss2[indx];
        break;
      case 3:
        x = xgauss3[indx];
        break;
      case 4:
        x = xgauss4[indx];
        break;
      case 5:
        x = xgauss5[indx];
        break;
      case 6:
        x = xgauss6[indx];
        break;
      case 7:
        x = xgauss7[indx];
        break;
      case 8:
        x = xgauss8[indx];
        break;
      case 9:
        x = xgauss9[indx];
        break;
      default:
        break;
      }
      return x;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the integration weight depending on the index of the Gauss Point
     * and integration order chosen by user (coming fro PASSMO).
     * This method is generic (not depending on the finite element type) and is called
     * by specialized methods (depending on FE type)
     */
    /*---------------------------------------------------------------------------*/
    static inline Real getWeight(Integer indx, Integer ordre)
    {
      Real w = wgauss1; // default is order 1

      switch (ordre) {
      case 2:
        w = wgauss2[indx];
        break;
      case 3:
        w = wgauss3[indx];
        break;
      case 4:
        w = wgauss4[indx];
        break;
      case 5:
        w = wgauss5[indx];
        break;
      case 6:
        w = wgauss6[indx];
        break;
      case 7:
        w = wgauss7[indx];
        break;
      case 8:
        w = wgauss8[indx];
        break;
      case 9:
        w = wgauss9[indx];
        break;
      default:
        break;
      }
      return w;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the reference coordinates of a Gauss Point for edge elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/

    static inline Real3 lineRefPosition(Integer3 indices, Integer3 ordre)
    {
      return { getRefPosition(indices[0], ordre[0]), 0., 0. };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the reference coordinates of a Gauss Point for triangle elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    static inline Real3 triRefPosition(Integer3 indices, Integer3 ordre)
    {
      Integer o = ordre[0] - 1;
      Integer i = indices[0];
      return { xg1[o][i], xg2[o][i], 0. };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the reference coordinates of a Gauss Point for quadrangle elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    static inline Real3 quadRefPosition(Integer3 indices, Integer3 ordre)
    {
      Real3 pos;
      for (Integer i = 0; i < 2; i++)
        pos[i] = getRefPosition(indices[i], ordre[i]);
      return pos;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the reference coordinates of a Gauss Point for hexaedron elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    static inline Real3 hexaRefPosition(Integer3 indices, Integer3 ordre)
    {
      Real3 pos;
      for (Integer i = 0; i < 3; i++)
        pos[i] = getRefPosition(indices[i], ordre[i]);
      return pos;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the reference coordinates of a Gauss Point for tetraedron elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    [[maybe_unused]] static inline Real3 tetraRefPosition(Integer3 indices, Integer3 /*ordre*/)
    {
      Integer i = indices[0];
      return { xit[i], yit[i], zit[i] };
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the reference coordinates of a Gauss Point for wedge (pentaedron)
     * elements (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    [[maybe_unused]] static inline Real3 pentaRefPosition(Integer3 indices, Integer3 ordre)
    {

      // Same as TriRefPosition on reference coordinate plane (r,s)
      // and LineRefPosition along reference coordinate t (vertical)
      auto pos = triRefPosition(indices, ordre);
      pos.z = getRefPosition(indices[2], ordre[2]);

      return pos;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the integration weight of a Gauss Point for edge elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    static inline Real lineWeight(Integer3 indices, Integer3 ordre)
    {
      return getWeight(indices[0], ordre[0]);
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the integration weight of a Gauss Point for triangle elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    static inline Real triWeight(Integer3 indices, Integer3 ordre)
    {
      return wg[ordre[0] - 1][indices[0]];
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the integration weight of a Gauss Point for quadrangle elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    static inline Real quadWeight(Integer3 indices, Integer3 ordre)
    {
      Real w = 1.;
      for (Integer i = 0; i < 2; i++)
        w *= getWeight(indices[i], ordre[i]);
      return w;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the integration weight of a Gauss Point for hexadron elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    static inline Real hexaWeight(Integer3 indices, Integer3 ordre)
    {
      Real w = 1.;
      for (Integer i = 0; i < 3; i++)
        w *= getWeight(indices[i], ordre[i]);
      return w;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the integration weight of a Gauss Point for tetraedron elements
     * (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    [[maybe_unused]] static inline Real tetraWeight(Integer3 /*indices*/, Integer3 /*ordre*/)
    {
      return wgtetra;
    }

    /*---------------------------------------------------------------------------*/
    /**
     * @brief Provides the integration weight of a Gauss Point for wedge (pentaedron)
     * elements (coming fro PASSMO).
     * This method takes the indices and integration orders chosen
     * by user as inputs (in (1, 2, 3 directions depending on the space dimension)
     */
    /*---------------------------------------------------------------------------*/
    [[maybe_unused]] static inline Real pentaWeight(Integer3 indices, Integer3 ordre)
    {

      // Same as TriWeight on reference coordinate plane (r,s)
      // and LineWeight with ordre[2] to account for reference coordinate t (vertical)
      Real wgpenta = triWeight(indices, ordre) * getWeight(indices[2], ordre[2]);
      return wgpenta;
    }

    /*---------------------------------------------------------------------------*/
  };
};

#endif // ARCANE_FEM_FUNCTIONS_H
