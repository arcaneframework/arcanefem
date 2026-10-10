// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DisplacementAndGradientFunctions.h                          (C) 2000-2026 */
/*                                                                           */
/* Functions to compute the displacement and the gradient.                   */
/*---------------------------------------------------------------------------*/
#ifndef ARCANFEM_ELASTOPLATICITY2_DISPLACEMENTANDGRADIENTFUNCTIONS_H
#define ARCANFEM_ELASTOPLATICITY2_DISPLACEMENTANDGRADIENTFUNCTIONS_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "femutils/FemUtils.h"
#include "femutils/ArcaneFemFunctions.h"
#include "femutils/ArcaneFemFunctionsGpu.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::ArcaneFem
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline Real3x3
computeDisplacementGradientQuad4(Cell cell, const VariableNodeReal3& node_coord,
                                 const VariableNodeReal3& displacement, Real xi, Real eta)
{
  const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad4(cell, node_coord, xi, eta);
  Real3x3 gradient{};
  for (Int32 i = 0; i < 4; ++i) {
    const Real3 u = displacement[cell.nodeId(i)];
    gradient(0, 0) += u.x * gp_info.dN_dx(i);
    gradient(0, 1) += u.x * gp_info.dN_dy(i);
    gradient(1, 0) += u.y * gp_info.dN_dx(i);
    gradient(1, 1) += u.y * gp_info.dN_dy(i);
  }
  return gradient;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline Real3x3
computeDisplacementGradientQuad8(Cell cell, const VariableNodeReal3& node_coord,
                                 const VariableNodeReal3& displacement, Real xi, Real eta)
{
  const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad8(cell, node_coord, xi, eta);
  Real3x3 gradient{};
  for (Int32 i = 0; i < 8; ++i) {
    const Real3 u = displacement[cell.nodeId(i)];
    gradient(0, 0) += u.x * gp_info.dN_dx(i);
    gradient(0, 1) += u.x * gp_info.dN_dy(i);
    gradient(1, 0) += u.y * gp_info.dN_dx(i);
    gradient(1, 1) += u.y * gp_info.dN_dy(i);
  }
  return gradient;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline Real3x3
computeDisplacementGradientQuad9(Cell cell, const VariableNodeReal3& node_coord,
                                 const VariableNodeReal3& displacement, Real xi, Real eta)
{
  const auto gp_info = ArcaneFemFunctions::FeOperation2D::computeGradientsAndJacobianQuad9(cell, node_coord, xi, eta);
  Real3x3 gradient{};
  for (Int32 i = 0; i < 9; ++i) {
    const Real3 u = displacement[cell.nodeId(i)];
    gradient(0, 0) += u.x * gp_info.dN_dx(i);
    gradient(0, 1) += u.x * gp_info.dN_dy(i);
    gradient(1, 0) += u.y * gp_info.dN_dx(i);
    gradient(1, 1) += u.y * gp_info.dN_dy(i);
  }
  return gradient;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::ArcaneFem

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
