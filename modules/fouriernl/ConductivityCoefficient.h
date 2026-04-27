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
#include <cmath>

#ifndef LAMBDA_CONDUCTION_COEFFICIENT
#define LAMBDA_CONDUCTION_COEFFICIENT(u, m) ( pow((1.0 + u), m) )
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
inline Real FemModuleFourierNL::_lambdaCpu(Real u) {
  return LAMBDA_CONDUCTION_COEFFICIENT(u, options()->expNlin);
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
ARCCORE_HOST_DEVICE __forceinline__ Real _lambdaGpu_m2(Real u) {
  return LAMBDA_CONDUCTION_COEFFICIENT(u, 2.);
}