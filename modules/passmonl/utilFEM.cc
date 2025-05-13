// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* utilFEM.h                                                   (C) 2022-2025 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "FemUtils.h"
#include "utilFEM.h"

using namespace Arcane::FemUtils;

//! Friend function for multiplication: Tensor4 * Tensor2
Arcane::FemUtils::Tensor2 operator*(const Tensor4& tens, const Arcane::FemUtils::Tensor2& vector) {
  Arcane::FemUtils::Tensor2 result;
  for (Arcane::Int32 i = 0; i < 6; ++i) {
    Real x{0};
    for (Arcane::Int32 j = 0; j < 6; ++j) {
      x += tens.value_at(i,j) * vector(j);
    }
    result(i) = x;
  }
  return result;
}
