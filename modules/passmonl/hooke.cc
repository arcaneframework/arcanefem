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
RealUniqueArray HookeInitConsts(RealConstArrayView& law_params)
{
  // auto E = law_params[0];
  // auto Nu = law_params[1];
  auto Lambda = law_params[0];
  auto Mu = law_params[1];

  RealUniqueArray consts(2);
  //  consts[0] = E * Nu / (1. - 2. * Nu) / (1. + Nu); // Lame coefficient Lambda
  //  consts[1] = E / 2. / (1. + Nu); // Lame coefficient Mu
  consts[0] = Lambda;
  consts[1] = Mu;

  return consts;
}

//! Computes elastic constitutive tensor
Tensor4 HookeComputeElastTensor(RealConstArrayView& law_params, const Tensor2& /*sig*/)
{
    RealUniqueArray consts = HookeInitConsts(law_params);
    return {consts[0]/*Lambda*/,consts[1]/*Mu*/};
}

//! Computes tangent constitutive tensor
Tensor4 HookeComputeTangentTensor(RealConstArrayView& law_params, RealArrayView& /*history_params*/,
                                  const Tensor2& /*sig*/, const Tensor2& /*deps*/)
{
    RealUniqueArray consts = HookeInitConsts(law_params);
    // Tangent stiffness tensor is equal to the elastic one
    return {consts[0]/*Lambda*/,consts[1]/*Mu*/};
}

//! Initializes the vector of intern (history) variables (nothing to do for this law)
RealUniqueArray HookeInitHistoryVars(RealConstArrayView& /*history_vars*/)
{
    return {};
}

//! Initializes the state (nothing to do for this law)
bool HookeInitState(const Tensor2& /*sig*/, RealArrayView& /*history_vars*/)
{
    return true;
}

//! Read constitutive parameters from a file and initialize intern constants allowing to the material constitutive model type
RealUniqueArray HookeReadLawParams(Real lambda, Real mu, bool /*default_param*/, const String& /*name*/, Integer /*ilaw*/)
{
  // =================================================================================
  // Hooke's model parameters stored in the lawparams vector:
  // 0-Lambda=1st Lame coef.  1-Mu=2nd Lame coef.
  // =================================================================================
  RealUniqueArray lawparams(2);
  lawparams[0] = lambda;
  lawparams[1] = mu;
  return lawparams;
}

Tensor4 HookeComputeStress(RealConstArrayView& law_params, RealArrayView& /*history_vars*/, Tensor2& sig, Tensor2& eps, Tensor2& /*epsp*/, Tensor2& dsig,
                        const Tensor2& deps, bool /*isRef*/)
{
  Tensor4 elast_tensor = HookeComputeElastTensor(law_params,sig);
	sig += elast_tensor*deps;
	eps += deps;
  return elast_tensor;
}

