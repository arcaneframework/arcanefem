// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DgQuadrature.h                                              (C) 2000-2026 */
/*                                                                           */
/* Provides quadrature rules used by discontinuous Galerkin methods.         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_DGQUADRATURE_H
#define ARCANEFEM_FEMUTILS_DGQUADRATURE_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemUtilsGlobal.h"

#include <arcane/core/Item.h>
#include <arcane/core/VariableTypes.h>
#include <arcane/utils/Array.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace ArcaneFemFunctions
{

using namespace Arcane;

struct DgQuadraturePoint
{
  Real3 point;
  Real weight;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Provides quadrature rules used by discontinuous Galerkin methods.
 */
class DgQuadrature
{
 public:

  /** @brief One-point centroid quadrature for a two-dimensional cell. */
  static void
  computeCellQuadrature2D(Cell cell, const VariableNodeReal3& node_coord,
                          Array<DgQuadraturePoint>& quadrature);

  /** @brief One-point midpoint quadrature for a two-dimensional edge. */
  static void
  computeFaceQuadrature2D(Face face, const VariableNodeReal3& node_coord,
                          Array<DgQuadraturePoint>& quadrature);

  /** @brief Degree-two quadrature for a planar polygonal face in 3D. */
  static void
  computeFaceQuadrature3D(Face face, const VariableNodeReal3& node_coord,
                          Array<DgQuadraturePoint>& quadrature);

  /** @brief Degree-three quadrature for a convex polyhedron. */
  static void
  computeCellQuadrature3D(Cell cell, const VariableNodeReal3& node_coord,
                          Array<DgQuadraturePoint>& quadrature);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace ArcaneFemFunctions

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
