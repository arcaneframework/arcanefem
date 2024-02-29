/*
* PASSMO : Performant Assessment for Seismic Site Modelling
*
* Definition of analytical functions for plane wave incident fields
*
* analytic_func.cpp
*
*  Created on: Feb. 2024
*      Author: E. Foerster
*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef ARCANEFEM_ANALYTICAL_FUNC_H
#define ARCANEFEM_ANALYTICAL_FUNC_H

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "arcane/MathUtils.h"
#include <arcane/VariableTypes.h>
using namespace Arcane;
extern Real	REL_PREC;

/*!
 * \brief Class to define analytical functions useful to define incident wave motions
 * on paraxial boundaries
 */
class AnalyticFunc{

 public:
  Real getRicker(const Real& time);
  Real getHarmonic(const Real& time);
  Real getTsang(const Real& time);
  Real getDecay(const Real& time);
  Real getDirac(const Real& time);

  Real m_tp{1.}, m_ts{1.}, m_phase{0.}, m_coef{1.}, m_amplit{1.};
  Int32 m_order{2};
};

#endif //ARCANEFEM_ANALYTICAL_FUNC_H
