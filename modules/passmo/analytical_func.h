// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* analytical_func.h                                           (C) 2022-2024 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_ANALYTICAL_FUNC_H
#define ARCANEFEM_ANALYTICAL_FUNC_H

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "arcane/MathUtils.h"
#include <arcane/VariableTypes.h>
using namespace Arcane;

/*!
 * \brief Class to define analytical functions useful to define incident wave motions
 * on paraxial boundaries
 */
class AnalyticFunc{

 public:
  [[nodiscard]] Real getRicker(const Real& time) const;
  [[nodiscard]] Real getHarmonic(const Real& time) const;
  [[nodiscard]] Real getTsang(const Real& time) const;
  [[nodiscard]] Real getDecay(const Real& time) const;
  [[nodiscard]] Real getDirac(const Real& time) const;

  Real m_tp{1.}, m_ts{1.}, m_phase{0.}, m_coef{1.}, m_amplit{1.};
  Int32 m_order{2};
};

#endif //ARCANEFEM_ANALYTICAL_FUNC_H
