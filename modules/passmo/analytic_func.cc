//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* analytical_func.cc                                          (C) 2022-2024 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "analytical_func.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Class to define analytical functions useful to define incident wave motions
 * on paraxial boundaries
 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Ricker analytical function.
 */
Real AnalyticFunc::getRicker(const Real& time) const{
  Real tol{1.0e-15};
  Real  PI{acos(-1.)},
  t{ PI * (time - m_ts) / m_tp },
  t2{ t * t },
  expt{ m_amplit * exp(-t2) },
  k{ 1.5 - sqrt(1.5) },
  val{0.};

  switch (m_order){
  case 0 : val = expt; break;

  case 1 : val = - sqrt(2.)*exp(0.5)*expt*t; break;

  case 3 : val = exp(k)*expt*(1.5 - t2)*t/sqrt(1.5)/sqrt(k); break;

  default : val = 0.5*exp(1.5)*expt*(t2 - 0.5); break; // Default = order 2
  }

  if (fabs(val) <= tol) val = 0.;
  return val;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Harmonic function: A*sin(2*PI*t/tp + phase)
 */
Real AnalyticFunc::getHarmonic(const Real& time) const{
  Real tol{1.0e-15};
  Real  PI{acos(-1.)};
  Real  val = m_amplit*sin(2*PI*time/m_tp + m_phase);
  if (fabs(val) <= tol) val = 0.;
  return val;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Decay analytical function: f(t) = A*exp(-coef*t)
 */
Real AnalyticFunc::getDecay(const Real& time) const{
  Real tol{1.0e-15};
  Real  x = -m_coef*time;
  if (x >= DBL_MAX_10_EXP) x = DBL_MAX_10_EXP;
  else if (x < DBL_MIN_10_EXP) x = DBL_MIN_10_EXP;

  Real val = m_amplit*exp(x);
  if (fabs(val) <= tol) val = 0.;
  return val;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Tsang function: : f(t) = A*sin(2*PI*t/tp + phase)*exp(-coef*t)
 */
Real AnalyticFunc::getTsang(const Real& time) const{
  Real tol{1.0e-15};
  Real  PI{acos(-1.)};
  Real  x = -m_coef*time;
  if (x >= DBL_MAX_10_EXP) x = DBL_MAX_10_EXP;
  else if (x < DBL_MIN_10_EXP) x = DBL_MIN_10_EXP;

  Real val = m_amplit*exp(x)*sin(2*PI*time/m_tp + m_phase);
  if (fabs(val) <= tol) val = 0.;
  return val;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Dirac function: : f(t) = A if t == tp, else = 0
 */
Real AnalyticFunc::getDirac(const Real& time) const{
  if (time != m_tp) return 0.;
  return m_amplit;
}
