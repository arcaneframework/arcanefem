// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ManufacturedSolutions.cc                                    (C) 2026      */
/*                                                                           */
/* Manufactured-solution assembly and checks for the Poisson module.         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include <arcane/IParallelMng.h>

/*---------------------------------------------------------------------------*/

namespace
{

constexpr Real Quad2x2GaussPoint[2] = {
  -0.57735026918962576451,
  0.57735026918962576451
};

constexpr Real Quad2x2Weight = 1.0;

//Works for Q2 quad or 3D Q2 hexa
constexpr Real Quad3x3GaussPoint[3] = {
  -0.77459666924148337704,
  0.0,
  0.77459666924148337704
};

constexpr Real Quad3x3Weight[3] = {
  5.0 / 9.0,
  8.0 / 9.0,
  5.0 / 9.0
};

} // namespace

/*---------------------------------------------------------------------------*/
/**
 * @brief Updates the nodal exact solution for the selected manufactured case.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
_updateManufacturedExactSolution()
{
  ENUMERATE_ (Node, inode, ownNodes()) {
    Node node = *inode;
    m_u_exact[node] = ManufacturedSolutions::dirichlet(m_manufactured_solution_name, m_node_coord[node]);
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the selected built-in manufactured source to the RHS.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
_applyManufacturedSource()
{
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable()); // Temporary variable to keep values for the RHS
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    const Integer nb_nodes = cell.nbNode();

    if (mesh()->dimension() == 2 && nb_nodes == 3) {
      const Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
      const Real3 x = ArcaneFemFunctions::MeshOperation::computeBaryCenterTria3(cell, m_node_coord);
      const Real source = ManufacturedSolutions::source(m_manufactured_solution_name, x);
      for (Node node : cell.nodes())
        if (node.isOwn())
          rhs_values[node_dof.dofId(node, 0)] += source * area / cell.nbNode();
    }
    else if (mesh()->dimension() == 2 && nb_nodes == 4) {
      for (Int32 ixi = 0; ixi < 2; ++ixi) {
        for (Int32 ieta = 0; ieta < 2; ++ieta) {
          const Real xi = Quad2x2GaussPoint[ixi];
          const Real eta = Quad2x2GaussPoint[ieta];
          const Real weight = Quad2x2Weight * Quad2x2Weight;
          const RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);
          const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, xi, eta);

          Real3 x(0.0, 0.0, 0.0);
          for (Int32 i = 0; i < 4; ++i)
            x += N[i] * m_node_coord[cell.node(i)];

          const Real source = ManufacturedSolutions::source(m_manufactured_solution_name, x);
          const Real integration_weight = gp_info.det_j * weight;
          for (Int32 i = 0; i < 4; ++i) {
            Node node = cell.node(i);
            if (node.isOwn())
              rhs_values[node_dof.dofId(node, 0)] += N[i] * source * integration_weight;
          }
        }
      }
    }
    else if (mesh()->dimension() == 2 && nb_nodes == 8) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          const Real xi = Quad3x3GaussPoint[ixi];
          const Real eta = Quad3x3GaussPoint[ieta];
          const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta];
          const RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad8(xi, eta);
          const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad8(cell, m_node_coord, xi, eta);

          Real3 x(0.0, 0.0, 0.0);
          for (Int32 i = 0; i < 8; ++i)
            x += N[i] * m_node_coord[cell.node(i)];

          const Real source = ManufacturedSolutions::source(m_manufactured_solution_name, x);
          const Real integration_weight = gp_info.det_j * weight;
          for (Int32 i = 0; i < 8; ++i) {
            Node node = cell.node(i);
            if (node.isOwn())
              rhs_values[node_dof.dofId(node, 0)] += N[i] * source * integration_weight;
          }
        }
      }
    }
    else if (mesh()->dimension() == 2 && nb_nodes == 9) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          const Real xi = Quad3x3GaussPoint[ixi];
          const Real eta = Quad3x3GaussPoint[ieta];
          const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta];
          const RealVector<9> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad9(xi, eta);
          const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad9(cell, m_node_coord, xi, eta);
        
          Real3 x(0.0, 0.0, 0.0);
          for (Int32 i = 0; i < 9; ++i)
            x += N[i] * m_node_coord[cell.node(i)];
        
          const Real source = ManufacturedSolutions::source(m_manufactured_solution_name, x);
          const Real integration_weight = gp_info.det_j * weight;
          for (Int32 i = 0; i < 9; ++i) {
            Node node = cell.node(i);
            if (node.isOwn())
              rhs_values[node_dof.dofId(node, 0)] += N[i] * source * integration_weight;
          }
        }
      }
    }
    else if (mesh()->dimension() == 3 && nb_nodes == 8) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          for (Int32 izeta = 0; izeta < 3; ++izeta){
            // Gauss point coordinates in reference space
            const Real xi = Quad3x3GaussPoint[ixi]; 
            const Real eta = Quad3x3GaussPoint[ieta]; 
            const Real zeta = Quad3x3GaussPoint[izeta]; 
            const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta] * Quad3x3Weight[izeta];
            const RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa8(xi,eta,zeta);
            const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa8(cell, m_node_coord, xi, eta, zeta);

            Real3 x(0.0, 0.0, 0.0);
            for (Int32 i = 0; i < 8; ++i)
              x += N[i] * m_node_coord[cell.node(i)];

            const Real source = ManufacturedSolutions::source(m_manufactured_solution_name, x);
            const Real integration_weight = gp_info.det_j * weight;
            for (Int32 i = 0; i < 8; ++i) {
              Node node = cell.node(i);
              if (node.isOwn())
                rhs_values[node_dof.dofId(node, 0)] += N[i] * source * integration_weight;
            }
          }
        }
      }
    }else if (mesh()->dimension() == 3 && nb_nodes == 20) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          for (Int32 izeta = 0; izeta < 3; ++izeta){
            // Gauss point coordinates in reference space
            const Real xi = Quad3x3GaussPoint[ixi]; 
            const Real eta = Quad3x3GaussPoint[ieta]; 
            const Real zeta = Quad3x3GaussPoint[izeta]; 
            const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta] * Quad3x3Weight[izeta];
            const RealVector<20> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa20(xi,eta,zeta);
            const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa20(cell, m_node_coord, xi, eta, zeta);

            Real3 x(0.0, 0.0, 0.0);
            for (Int32 i = 0; i < 20; ++i)
              x += N[i] * m_node_coord[cell.node(i)];

            const Real source = ManufacturedSolutions::source(m_manufactured_solution_name, x);
            const Real integration_weight = gp_info.det_j * weight;
            for (Int32 i = 0; i < 20; ++i) {
              Node node = cell.node(i);
              if (node.isOwn())
                rhs_values[node_dof.dofId(node, 0)] += N[i] * source * integration_weight;
            }
          }
        }
      }
    }else if (mesh()->dimension() == 3 && nb_nodes == 27) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          for (Int32 izeta = 0; izeta < 3; ++izeta){
            // Gauss point coordinates in reference space
            const Real xi = Quad3x3GaussPoint[ixi]; 
            const Real eta = Quad3x3GaussPoint[ieta]; 
            const Real zeta = Quad3x3GaussPoint[izeta]; 
            const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta] * Quad3x3Weight[izeta];
            const RealVector<27> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa27(xi,eta,zeta);
            const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa27(cell, m_node_coord, xi, eta, zeta);

            Real3 x(0.0, 0.0, 0.0);
            for (Int32 i = 0; i < 27; ++i)
              x += N[i] * m_node_coord[cell.node(i)];

            const Real source = ManufacturedSolutions::source(m_manufactured_solution_name, x);
            const Real integration_weight = gp_info.det_j * weight;
            for (Int32 i = 0; i < 27; ++i) {
              Node node = cell.node(i);
              if (node.isOwn())
                rhs_values[node_dof.dofId(node, 0)] += N[i] * source * integration_weight;
            }
          }
        }
      }
    
    }else {
      ARCANE_FATAL("Manufactured source is not implemented for dimension '{0}' cells with '{1}' nodes", mesh()->dimension(), nb_nodes);
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the selected built-in manufactured Dirichlet condition.
 */
/*---------------------------------------------------------------------------*/

void FemModulePoisson::
_applyManufacturedDirichlet()
{
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  NodeGroup node_group = mesh()->outerFaces().nodeGroup();

  ENUMERATE_ (Node, inode, node_group) {
    Node node = *inode;
    if (node.isOwn()) {
      const Real value = ManufacturedSolutions::dirichlet(m_manufactured_solution_name, m_node_coord[node]);
      m_linear_system.matrixSetValue(node_dof.dofId(node, 0), node_dof.dofId(node, 0), 1.0e30);
      rhs_values[node_dof.dofId(node, 0)] = 1.0e30 * value;
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the L2 error against the manufactured solution.
 */
/*---------------------------------------------------------------------------*/

Real FemModulePoisson::
_computeManufacturedL2Error()
{
  Real l2_error_square = 0.0;

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    const Integer nb_nodes = cell.nbNode();

    if (mesh()->dimension() == 2 && nb_nodes == 3) {
      const Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
      const Real3 x = ArcaneFemFunctions::MeshOperation::computeBaryCenterTria3(cell, m_node_coord);
      const Real uh = (m_u[cell.node(0)] + m_u[cell.node(1)] + m_u[cell.node(2)]) / 3.0;
      const Real error = uh - ManufacturedSolutions::dirichlet(m_manufactured_solution_name, x);
      l2_error_square += error * error * area;
    }
    else if (mesh()->dimension() == 2 && nb_nodes == 4) {
      for (Int32 ixi = 0; ixi < 2; ++ixi) {
        for (Int32 ieta = 0; ieta < 2; ++ieta) {
          const Real xi = Quad2x2GaussPoint[ixi];
          const Real eta = Quad2x2GaussPoint[ieta];
          const Real weight = Quad2x2Weight * Quad2x2Weight;
          const RealVector<4> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad4(xi, eta);
          const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, m_node_coord, xi, eta);

          Real uh = 0.0;
          Real3 x(0.0, 0.0, 0.0);
          for (Int32 i = 0; i < 4; ++i) {
            uh += N[i] * m_u[cell.node(i)];
            x += N[i] * m_node_coord[cell.node(i)];
          }
          const Real error = uh - ManufacturedSolutions::dirichlet(m_manufactured_solution_name, x);
          l2_error_square += error * error * gp_info.det_j * weight;
        }
      }
    }
    else if (mesh()->dimension() == 2 && nb_nodes == 8) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          const Real xi = Quad3x3GaussPoint[ixi];
          const Real eta = Quad3x3GaussPoint[ieta];
          const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta];
          const RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad8(xi, eta);
          const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad8(cell, m_node_coord, xi, eta);

          Real uh = 0.0;
          Real3 x(0.0, 0.0, 0.0);
          for (Int32 i = 0; i < 8; ++i) {
            uh += N[i] * m_u[cell.node(i)];
            x += N[i] * m_node_coord[cell.node(i)];
          }
          const Real error = uh - ManufacturedSolutions::dirichlet(m_manufactured_solution_name, x);
          l2_error_square += error * error * gp_info.det_j * weight;
        }
      }
    }
    else if (mesh()->dimension() == 2 && nb_nodes == 9) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          const Real xi = Quad3x3GaussPoint[ixi];
          const Real eta = Quad3x3GaussPoint[ieta];
          const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta];
          const RealVector<9> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsQuad9(xi, eta);
          const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad9(cell, m_node_coord, xi, eta);
        
          Real uh = 0.0;
          Real3 x(0.0, 0.0, 0.0);
          for (Int32 i = 0; i < 9; ++i) {
            uh += N[i] * m_u[cell.node(i)];
            x += N[i] * m_node_coord[cell.node(i)];
          }
          const Real error = uh - ManufacturedSolutions::dirichlet(m_manufactured_solution_name, x);
          l2_error_square += error * error * gp_info.det_j * weight;
        }
      }
    }
    else if (mesh()->dimension() == 3 && nb_nodes == 8) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          for (Int32 izeta = 0; izeta < 3; ++izeta){
            // Gauss point coordinates in reference space
            const Real xi = Quad3x3GaussPoint[ixi]; 
            const Real eta = Quad3x3GaussPoint[ieta]; 
            const Real zeta = Quad3x3GaussPoint[izeta]; 
            const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta] * Quad3x3Weight[izeta];
            const RealVector<8> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa8(xi,eta,zeta);
            const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa8(cell, m_node_coord, xi, eta, zeta);

            Real uh = 0.0;
            Real3 x(0.0, 0.0, 0.0);
            for (Int32 i = 0; i < 8; ++i) {
              uh += N[i] * m_u[cell.node(i)];
              x += N[i] * m_node_coord[cell.node(i)];
            }
            const Real error = uh - ManufacturedSolutions::dirichlet(m_manufactured_solution_name, x);
            l2_error_square += error * error * gp_info.det_j * weight;            }
          }
        }
      }
    else if (mesh()->dimension() == 3 && nb_nodes == 20) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          for (Int32 izeta = 0; izeta < 3; ++izeta){
            // Gauss point coordinates in reference space
            const Real xi = Quad3x3GaussPoint[ixi]; 
            const Real eta = Quad3x3GaussPoint[ieta]; 
            const Real zeta = Quad3x3GaussPoint[izeta]; 
            const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta] * Quad3x3Weight[izeta];
            const RealVector<20> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa20(xi,eta,zeta);
            const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa20(cell, m_node_coord, xi, eta, zeta);

            Real uh = 0.0;
            Real3 x(0.0, 0.0, 0.0);
            for (Int32 i = 0; i < 20; ++i) {
              uh += N[i] * m_u[cell.node(i)];
              x += N[i] * m_node_coord[cell.node(i)];
            }
            const Real error = uh - ManufacturedSolutions::dirichlet(m_manufactured_solution_name, x);
            l2_error_square += error * error * gp_info.det_j * weight;          }
          }
        }
      }
    else if (mesh()->dimension() == 3 && nb_nodes == 27) {
      for (Int32 ixi = 0; ixi < 3; ++ixi) {
        for (Int32 ieta = 0; ieta < 3; ++ieta) {
          for (Int32 izeta = 0; izeta < 3; ++izeta){
            // Gauss point coordinates in reference space
            const Real xi = Quad3x3GaussPoint[ixi]; 
            const Real eta = Quad3x3GaussPoint[ieta]; 
            const Real zeta = Quad3x3GaussPoint[izeta]; 
            const Real weight = Quad3x3Weight[ixi] * Quad3x3Weight[ieta] * Quad3x3Weight[izeta];
            const RealVector<27> N = Arcane::FemUtils::ShapeFunctions::computeShapeFunctionsHexa27(xi,eta,zeta);
            const auto gp_info = ArcaneFemFunctions::FeOperation3D::computeGradientsAndJacobianHexa27(cell, m_node_coord, xi, eta, zeta);

            Real uh = 0.0;
            Real3 x(0.0, 0.0, 0.0);
            for (Int32 i = 0; i < 27; ++i) {
              uh += N[i] * m_u[cell.node(i)];
              x += N[i] * m_node_coord[cell.node(i)];
            }
            const Real error = uh - ManufacturedSolutions::dirichlet(m_manufactured_solution_name, x);
            l2_error_square += error * error * gp_info.det_j * weight;            }
          }
        }
      }
    else {
      ARCANE_FATAL("Manufactured L2 error is not implemented for dimension '{0}' cells with '{1}' nodes", mesh()->dimension(), nb_nodes);
    }
  }

  l2_error_square = mesh()->parallelMng()->reduce(Parallel::ReduceSum, l2_error_square);

  return Arcane::math::sqrt(l2_error_square);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
