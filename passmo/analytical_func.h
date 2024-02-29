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
  void setTp(const Real& tp);
  void setTs(const Real& ts);
  void setAmplit(const Real& amplit);
  void setCoef(const Real& coef);
  void setPhase(const Real& phase);
  void setOrder(const Int32& order);

  Real getRicker(const Real& time);
  Real getHarmonic(const Real& time);
  Real getTsang(const Real& time);
  Real getDecay(const Real& time);
  Real getDirac(const Real& time);

 private:
  Real m_tp{0.}, m_ts{0.}, m_phase{0.}, m_coef{0.}, m_amplit{0.};
  Int32 m_order{0};
};

#endif //ARCANEFEM_ANALYTICAL_FUNC_H
