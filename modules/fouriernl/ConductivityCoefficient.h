// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ElementMatrix.h                                             (C) 2022-2026 */
/*                                                                           */
/* Contains inline functions to compute the nonconstant conduction           */
/* coefficient for Nonlinear Fourier                                         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifndef LAMBDA_CONDUCTIVITY_COEFFICIENT
#define LAMBDA_CONDUCTIVITY_COEFFICIENT(u, m) ( math::pow((1.0 + u), m) )
#endif

/**
 * @brief Computes inline conduction coefficient for CPU
 *
 * This function calculates the expression:
 *       λ(𝑢)= LAMBDA_CONDUCTION_COEFFICIENT
 *
 * Steps involved:
 * 1. Calculate the inline conduction coefficient.
 * 2. Return λ(𝑢);
 */
/*---------------------------------------------------------------------------*/
inline Real FemModuleFourierNL::_lambdaCpu(const Real& u) {
  return LAMBDA_CONDUCTIVITY_COEFFICIENT(u, m_nlin_exp);
}

/**
 * @brief Computes inline conduction coefficient for GPU
 *
 * This function calculates the expression:
 *       λ(𝑢)= LAMBDA_CONDUCTION_COEFFICIENT
 *
 * Steps involved:
 * 1. Calculate the inline conduction coefficient.
 * 2. Return λ(𝑢);
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE inline Real _lambdaGpu_m2(const Real& u) {
  // return LAMBDA_CONDUCTIVITY_COEFFICIENT(u, 2.);
  return (1.0 + 2.0 * u + u * u);
}

ARCCORE_HOST_DEVICE inline Real _lambdaGpu_m0(const Real& u) {
  return 1.0;
}

ARCCORE_HOST_DEVICE inline Real _lambdaGpu_m5(const Real& u) {
  return LAMBDA_CONDUCTIVITY_COEFFICIENT(u, 5.);
}

ARCCORE_HOST_DEVICE inline Real _lambdaGpu_maybe_slow(const Real& u, const Real expNlin) {
  return LAMBDA_CONDUCTIVITY_COEFFICIENT(u, expNlin);
}

ARCCORE_HOST_DEVICE inline Real _lambdaGpu(const Real& u, const Real expNlin)
{
  if (expNlin == 0.) {
    return _lambdaGpu_m0(u);
  } else if (expNlin == 2.) {
    return _lambdaGpu_m2(u);
  } else if (expNlin == 5.) {
    return _lambdaGpu_m5(u);
  } else {
    // ARCANE_FATAL("The chosen value of m is not implemented for GPUs");
    return _lambdaGpu_maybe_slow(u, expNlin);
  }
}