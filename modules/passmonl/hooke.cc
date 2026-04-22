// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-20245CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Hooke.cc                                                    (C) 2022-2025 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Implementation of Drücker-Prager constitutive law
 */
#include "TypesNLdynamic.h"
//#include "FemUtils.h"
#include "utilFEM.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;

extern Arcane::FemUtils::Tensor2 operator*(const Tensor4& tens, const Arcane::FemUtils::Tensor2& vector);
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
   * @brief Implementation of Hooke's constitutive model (linear elasticity)
   */
/*---------------------------------------------------------------------------*/

//! Initialize intern useful constants
//RealUniqueArray HookeInitConsts(RealConstArrayView& law_params)
RealUniqueArray HookeInitConsts(RealUniqueArray* law_params)
{
  auto Lambda = (*law_params)[0];
  auto Mu = (*law_params)[1];

  RealUniqueArray consts(2);
  consts[0] = Lambda;
  consts[1] = Mu;

  return consts;
}

//! Computes elastic constitutive tensor
Tensor4 HookeComputeElastTensor(RealUniqueArray* law_params, const Tensor2& /*sig*/)
{
    return {(*law_params)[0]/*Lambda*/,(*law_params)[1]/*Mu*/};
}

//! Computes tangent constitutive (stiffness) tensor
Tensor4 HookeComputeTangentTensor(RealUniqueArray* law_params, RealUniqueArray* /*history_params*/,
                                  const Tensor2& /*sig*/, const Tensor2& /*deps*/)
{
    // Tangent stiffness tensor is equal to the elastic one
    return {(*law_params)[0]/*Lambda*/,(*law_params)[1]/*Mu*/};
}

//! Initializes the vector of intern (history) variables (nothing to do for this law)
void HookeInitHistoryVars(RealUniqueArray* /*history_vars*/)
{}

//! Initializes the state (nothing to do for this law)
bool HookeInitState(const Tensor2& /*sig*/, RealUniqueArray* /*history_vars*/)
{
    return true;
}

//! Read constitutive parameters from a file and initialize intern constants allowing to the material constitutive model type
void HookeReadLawParams(RealUniqueArray* lawparams, Real lambda, Real mu, bool /*default_param*/, const String& /*name*/, Integer /*ilaw*/)
{
  // =================================================================================
  // Hooke's model parameters stored in the lawparams vector:
  // 0-Lambda=1st Lame coef.  1-Mu=2nd Lame coef.
  // =================================================================================
  (*lawparams)[0] = lambda;
  (*lawparams)[1] = mu;
}

bool HookeComputeStress(RealUniqueArray* law_params, RealUniqueArray* /*history_vars*/, Tensor2& sig, Tensor2& eps, Tensor2& /*epsp*/, Tensor2& dsig,
                        const Tensor2& deps, Tensor4& /*tangent_tensor*/, bool /*isRef*/)
{
  Tensor4 elast_tensor = HookeComputeElastTensor(law_params,sig);
  elast_tensor.isSymmetric(true);
  elast_tensor.isConstitutive(true);

	sig += elast_tensor*deps;
	eps += deps;
  return true;
}

