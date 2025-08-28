// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Dirichlet.h                                                 (C) 2022-2025 */
/*                                                                           */
/* Contains functions to compute and assemble source term contribution to RHS*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies source term to RHS vector of the linear system.
 * 
 * @param rhs_values The variable representing the RHS vector to be updated.
 * @param node_dof The connectivity view mapping nodes to their corresponding
 */
/*---------------------------------------------------------------------------*/

inline void FemModule::
_applySourceTerm(VariableDoFReal& rhs_values, const IndexedNodeDoFConnectivityView& node_dof)
{
    if (mesh()->dimension() == 2)
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);
      Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
      Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);

      RealVector<6> Uy = { 0., 1., 0., 1., 0., 1. };
      RealVector<6> Ux = { 1., 0., 1., 0., 1., 0. };
      RealVector<6> F = { f[0], f[1], f[0], f[1], f[0], f[1] };
      RealVector<6> dxUx = { dxu[0], 0., dxu[1], 0., dxu[2], 0. };
      RealVector<6> dyUx = { dyu[0], 0., dyu[1], 0., dyu[2], 0. };
      RealVector<6> dxUy = { 0., dxu[0], 0., dxu[1], 0., dxu[2] };
      RealVector<6> dyUy = { 0., dyu[0], 0., dyu[1], 0., dyu[2] };

      RealVector<6> Un = { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y,
                               m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y,
                               m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y };

                               RealVector<6>Vn = { m_V[cell.nodeId(0)].x, m_V[cell.nodeId(0)].y,
                               m_V[cell.nodeId(1)].x, m_V[cell.nodeId(1)].y,
                               m_V[cell.nodeId(2)].x, m_V[cell.nodeId(2)].y };

                               RealVector<6> An = { m_A[cell.nodeId(0)].x, m_A[cell.nodeId(0)].y,
                               m_A[cell.nodeId(1)].x, m_A[cell.nodeId(1)].y,
                               m_A[cell.nodeId(2)].x, m_A[cell.nodeId(2)].y };

      //----------------------------------------------------------------------
      //  ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯) +
      //  ∫∫∫ (c₅)(∇𝐮ₙ.∇𝐯) + ∫∫∫ (c₆)(ε(𝐮ₙ):ε(𝐯)) +
      //  ∫∫∫ (c₇)(∇𝐮ᵗₙ.∇𝐯) + ∫∫∫ (c₉)(ε(𝐮ᵗₙ):ε(𝐯)) +
      //  ∫∫∫ (c₈)(∇𝐮ᵗᵗₙ.∇𝐯) + ∫∫∫ (c₁₀)(ε(𝐮ᵗᵗₙ):ε(𝐯))
      //----------------------------------------------------------------------
      RealVector<6> rhs = ( F * (1/3.)
                              + Un * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy))*(c0*1/12.)
                              + Vn * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy))*(c3*1/12.)
                              + An * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy))*(c4*1/12.)
                              - Un * ((dyUy ^ dxUx) + (dxUx ^ dyUy)  +  (dxUx ^ dxUx) + (dyUy ^ dyUy)) * c5
                              - Un * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy)) +  ((dxUy + dyUx) ^ (dyUx + dxUy)))*c6
                              + Vn * ((dyUy ^ dxUx) + (dxUx ^ dyUy)  +  (dxUx ^ dxUx) + (dyUy ^ dyUy)) * c7
                              + Vn * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy)) +  ((dxUy + dyUx) ^ (dyUx + dxUy)))*c9
                              + An * ((dyUy ^ dxUx) + (dxUx ^ dyUy)  +  (dxUx ^ dxUx) + (dyUy ^ dyUy)) * c8
                              + An * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy)) +  ((dxUy + dyUx) ^ (dyUx + dxUy)))*c10
                              ) * area;

      rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
      rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
      rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(2);
      rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(3);
      rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(4);
      rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(5);
    }

  if (mesh()->dimension() == 3)
    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      Real volume = ArcaneFemFunctions::MeshOperation::computeVolumeTetra4(cell, m_node_coord);
      Real4 dxu = ArcaneFemFunctions::FeOperation3D::computeGradientXTetra4(cell, m_node_coord);
      Real4 dyu = ArcaneFemFunctions::FeOperation3D::computeGradientYTetra4(cell, m_node_coord);
      Real4 dzu = ArcaneFemFunctions::FeOperation3D::computeGradientZTetra4(cell, m_node_coord);

      RealVector<12> Uy = { 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0. };
      RealVector<12> Ux = { 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0. };
      RealVector<12> Uz = { 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1. };

      RealVector<12> F = { f[0], f[1], f[2], f[0], f[1], f[2], f[0], f[1], f[2], f[0], f[1], f[2] };
      RealVector<12> dxUx = { dxu[0], 0., 0.,    dxu[1], 0., 0.,    dxu[2], 0., 0.,    dxu[3], 0., 0. };
      RealVector<12> dyUx = { dyu[0], 0., 0.,    dyu[1], 0., 0.,    dyu[2], 0., 0.,    dyu[3], 0., 0. };
      RealVector<12> dzUx = { dzu[0], 0., 0.,    dzu[1], 0., 0.,    dzu[2], 0., 0.,    dzu[3], 0., 0. };

      RealVector<12> dxUy = { 0., dxu[0], 0.,    0., dxu[1], 0.,    0., dxu[2], 0.,    0., dxu[3], 0. };
      RealVector<12> dyUy = { 0., dyu[0], 0.,    0., dyu[1], 0.,    0., dyu[2], 0.,    0., dyu[3], 0. };
      RealVector<12> dzUy = { 0., dzu[0], 0.,    0., dzu[1], 0.,    0., dzu[2], 0.,    0., dzu[3], 0. };

      RealVector<12> dxUz = { 0., 0., dxu[0],    0., 0., dxu[1],    0., 0., dxu[2],    0., 0., dxu[3] };
      RealVector<12> dyUz = { 0., 0., dyu[0],    0., 0., dyu[1],    0., 0., dyu[2],    0., 0., dyu[3] };
      RealVector<12> dzUz = { 0., 0., dzu[0],    0., 0., dzu[1],    0., 0., dzu[2],    0., 0., dzu[3] };

      RealVector<12> Un = { m_U[cell.nodeId(0)].x, m_U[cell.nodeId(0)].y, m_U[cell.nodeId(0)].z,
                                m_U[cell.nodeId(1)].x, m_U[cell.nodeId(1)].y, m_U[cell.nodeId(1)].z,
                                m_U[cell.nodeId(2)].x, m_U[cell.nodeId(2)].y, m_U[cell.nodeId(2)].z,
                                m_U[cell.nodeId(3)].x, m_U[cell.nodeId(3)].y, m_U[cell.nodeId(3)].z };

      RealVector<12> Vn = { m_V[cell.nodeId(0)].x, m_V[cell.nodeId(0)].y, m_V[cell.nodeId(0)].z,
                                m_V[cell.nodeId(1)].x, m_V[cell.nodeId(1)].y, m_V[cell.nodeId(1)].z,
                                m_V[cell.nodeId(2)].x, m_V[cell.nodeId(2)].y, m_V[cell.nodeId(2)].z,
                                m_V[cell.nodeId(3)].x, m_V[cell.nodeId(3)].y, m_V[cell.nodeId(3)].z };

      RealVector<12> An = { m_A[cell.nodeId(0)].x, m_A[cell.nodeId(0)].y, m_A[cell.nodeId(0)].z,
                                m_A[cell.nodeId(1)].x, m_A[cell.nodeId(1)].y, m_A[cell.nodeId(1)].z,
                                m_A[cell.nodeId(2)].x, m_A[cell.nodeId(2)].y, m_A[cell.nodeId(2)].z,
                                m_A[cell.nodeId(3)].x, m_A[cell.nodeId(3)].y, m_A[cell.nodeId(3)].z };

      //----------------------------------------------------------------------
      //  ∫∫∫ (𝐟.𝐯) + ∫∫∫ (c₀)(𝐮ₙ.𝐯) + ∫∫∫ (c₃)(𝐮ᵗₙ.𝐯) + ∫∫∫ (c₄)(𝐮ᵗᵗₙ.𝐯) +
      //  ∫∫∫ (c₅)(∇𝐮ₙ.∇𝐯) + ∫∫∫ (c₆)(ε(𝐮ₙ):ε(𝐯)) +
      //  ∫∫∫ (c₇)(∇𝐮ᵗₙ.∇𝐯) + ∫∫∫ (c₉)(ε(𝐮ᵗₙ):ε(𝐯)) +
      //  ∫∫∫ (c₈)(∇𝐮ᵗᵗₙ.∇𝐯) + ∫∫∫ (c₁₀)(ε(𝐮ᵗᵗₙ):ε(𝐯))
      //----------------------------------------------------------------------
      RealVector<12> rhs = ( F * (1/4.)
                              + Un * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy) + massMatrix(Uz,Uz))*(c0*1/20.)
                              + Vn * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy) + massMatrix(Uz,Uz))*(c3*1/20.)
                              + An * (massMatrix(Ux,Ux) + massMatrix(Uy,Uy) + massMatrix(Uz,Uz))*(c4*1/20.)
                              - Un * ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) +
                                      (dyUy ^ dxUx) + (dxUx ^ dyUy) +
                                      (dzUz ^ dxUx) + (dxUx ^ dzUz) +
                                      (dyUy ^ dzUz) + (dzUz ^ dyUy)) * c5
                              - Un * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) ) +
                                         ( ((dxUy + dyUx) ^ (dyUx + dxUy)) +
                                           ((dzUy + dyUz) ^ (dyUz + dzUy)) +
                                           ((dxUz + dzUx) ^ (dzUx + dxUz)) ))*c6
                              + Vn * ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) +
                                      (dyUy ^ dxUx) + (dxUx ^ dyUy) +
                                      (dzUz ^ dxUx) + (dxUx ^ dzUz) +
                                      (dyUy ^ dzUz) + (dzUz ^ dyUy)) * c7
                              + Vn * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) ) +
                                      ( ((dxUy + dyUx) ^ (dyUx + dxUy)) +
                                        ((dzUy + dyUz) ^ (dyUz + dzUy)) +
                                        ((dxUz + dzUx) ^ (dzUx + dxUz)) ))*c9
                                        + An * ((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) +
                                        (dyUy ^ dxUx) + (dxUx ^ dyUy) +
                                        (dzUz ^ dxUx) + (dxUx ^ dzUz) +
                                        (dyUy ^ dzUz) + (dzUz ^ dyUy)) * c8
                                + An * (2.*((dxUx ^ dxUx) + (dyUy ^ dyUy) + (dzUz ^ dzUz) ) +
                                        ( ((dxUy + dyUx) ^ (dyUx + dxUy)) +
                                          ((dzUy + dyUz) ^ (dyUz + dzUy)) +
                                          ((dxUz + dzUx) ^ (dzUx + dxUz)) ))*c10
                              ) * volume;

      rhs_values[node_dof.dofId(cell.nodeId(0), 0)] += rhs(0);
      rhs_values[node_dof.dofId(cell.nodeId(0), 1)] += rhs(1);
      rhs_values[node_dof.dofId(cell.nodeId(0), 2)] += rhs(2);
      rhs_values[node_dof.dofId(cell.nodeId(1), 0)] += rhs(3);
      rhs_values[node_dof.dofId(cell.nodeId(1), 1)] += rhs(4);
      rhs_values[node_dof.dofId(cell.nodeId(1), 2)] += rhs(5);
      rhs_values[node_dof.dofId(cell.nodeId(2), 0)] += rhs(6);
      rhs_values[node_dof.dofId(cell.nodeId(2), 1)] += rhs(7);
      rhs_values[node_dof.dofId(cell.nodeId(2), 2)] += rhs(8);
      rhs_values[node_dof.dofId(cell.nodeId(3), 0)] += rhs(9);
      rhs_values[node_dof.dofId(cell.nodeId(3), 1)] += rhs(10);
      rhs_values[node_dof.dofId(cell.nodeId(3), 2)] += rhs(11);
    }  
}