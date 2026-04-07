// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-20245CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* LawDispatcher.cc                                            (C) 2022-2025 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "TypesNLdynamic.h"
#include "FemUtils.h"
#include "utilFEM.h"
#include "LawDispatcher.h"

using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/////////////////////////////////////////////////////////////////////////////
// Linear elasticity model: functions implemented in hooke.cc

extern RealUniqueArray HookeInitConsts(RealUniqueArray* /*law_params*/);
extern Tensor4 HookeComputeElastTensor(RealUniqueArray* /*law_params*/, const Tensor2& /*sig*/);
extern Tensor4 HookeComputeTangentTensor(RealUniqueArray* /*law_params*/, RealUniqueArray* /*history_vars*/, const Tensor2& /*sig*/, const Tensor2& /*deps*/);
extern void HookeInitHistoryVars(RealUniqueArray* /*history_vars*/);
extern bool HookeInitState(const Tensor2& /*sig*/, RealUniqueArray* /*history_vars*/);
extern void HookeReadLawParams(RealUniqueArray* /*law_params*/, Real /*lambda*/, Real /*mu*/, bool /*default_param*/, const String& /*name*/, Integer /*ilaw*/);
extern Tensor4 HookeComputeStress(RealUniqueArray* /*law_params*/, RealUniqueArray* /*history_vars*/, Tensor2& /*sig*/, Tensor2& /*eps*/, Tensor2& /*epsp*/, Tensor2& /*dsig*/,
                               const Tensor2& /*deps*/, bool /*isRef*/);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/////////////////////////////////////////////////////////////////////////////
// Drücker-Prager model: functions implemented in druckp.cc

extern RealUniqueArray DruckPInitConsts(RealUniqueArray* /*law_params*/);
extern Tensor4 DruckPComputeElastTensor(RealUniqueArray* /*law_params*/, const Tensor2& /*sig*/);
extern Tensor4 DruckPComputeTangentTensor(RealUniqueArray* /*law_params*/, RealUniqueArray* /*history_vars*/, const Tensor2& /*sig*/, const Tensor2& /*deps*/);
extern void DruckPInitHistoryVars(RealUniqueArray* /*history_vars*/);
extern bool DruckPInitState(const Tensor2& /*sig*/, RealUniqueArray* /*history_vars*/);
extern void DruckPReadLawParams(RealUniqueArray* /*law_params*/, Real /*lambda*/, Real /*mu*/, bool /*default_param*/, const String& /*name*/, Integer /*ilaw*/);
extern Tensor4 DruckPComputeStress(RealUniqueArray* /*law_params*/, RealUniqueArray* /*history_vars*/, Tensor2& /*sig*/, Tensor2& /*eps*/, Tensor2& /*epsp*/, Tensor2& /*dsig*/,
                               const Tensor2& /*deps*/, bool /*isRef*/);
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

LawDispatcher::LawDispatcher(TypesNLDynamic::eLawType law_type, bool default_param) :
        m_law_type(law_type), m_default(default_param) {

    // Setting to null default value
    for (Integer i = 0; i < NB_LAW_TYPE; ++i) {
        m_compute_stress[i] = nullptr;
        m_compute_elast_tensor[i] = nullptr;
        m_compute_tangent_tensor[i] = nullptr;
        m_init_state[i] = nullptr;
        m_init_history_vars[i] = nullptr;
        m_read_law_params[i] = nullptr;
        m_init_consts[i] = nullptr;
    }

    //! Computes stresses allowing to the material constitutive model type
    m_compute_stress[TypesNLDynamic::HOOKE] = HookeComputeStress;
    m_compute_stress[TypesNLDynamic::DRUCKP] = DruckPComputeStress;

    //! Computes elastic constitutive tensor allowing to the material constitutive model type
    m_compute_elast_tensor[TypesNLDynamic::HOOKE] = HookeComputeElastTensor;
    m_compute_elast_tensor[TypesNLDynamic::DRUCKP] = DruckPComputeElastTensor;

    //! Computes the tangent constitutive tensor allowing to the material constitutive model type (used for nonlinear models only)
    m_compute_tangent_tensor[TypesNLDynamic::HOOKE] = HookeComputeTangentTensor;
    m_compute_tangent_tensor[TypesNLDynamic::DRUCKP] = DruckPComputeTangentTensor;

    //! Initializes the state allowing to the material constitutive model type (used for nonlinear models only)
    m_init_state[TypesNLDynamic::HOOKE] = HookeInitState;
    m_init_state[TypesNLDynamic::DRUCKP] = DruckPInitState;

    //! Initializes the vector of intern (history) variables allowing to the material constitutive model type (used for nonlinear models only)
    m_init_history_vars[TypesNLDynamic::HOOKE] = HookeInitHistoryVars;
    m_init_history_vars[TypesNLDynamic::DRUCKP] = DruckPInitHistoryVars;

    //! Read constitutive parameters from a file and initialize intern constants allowing to the material constitutive model type
    m_read_law_params[TypesNLDynamic::HOOKE] = HookeReadLawParams;
    m_read_law_params[TypesNLDynamic::DRUCKP] = DruckPReadLawParams;

    //! Initialize intern useful constants allowing to the material constitutive model type
    m_init_consts[TypesNLDynamic::HOOKE] = HookeInitConsts;
    m_init_consts[TypesNLDynamic::DRUCKP] = DruckPInitConsts;

    switch(law_type){

    case TypesNLDynamic::HOOKE: m_nb_law_param = 2; m_nb_law_history_param = 0;
      break;

    case TypesNLDynamic::DRUCKP:
    case TypesNLDynamic::MOHRC: m_nb_law_param = 7; m_nb_law_history_param = 1;
      break;

    case TypesNLDynamic::UNKNOWN:
    default: m_nb_law_param = 2; m_nb_law_history_param = 0;
      break;
    }
}

/////////////////////////////////////////////////////////////////////////////
// class MaterialDispatcher: implementation methods

Tensor2	LawDispatcher::getStress() const { return m_sig; }
void	LawDispatcher::setStress(const Tensor2& tensor) { m_sig = tensor; }

Tensor2	LawDispatcher::getStrain() const { return m_eps; }
void	LawDispatcher::setStrain(const Tensor2& tensor) { m_eps = tensor; }

Tensor2	LawDispatcher::getPlasticStrain() const { return m_epsp; }
void	LawDispatcher::setPlasticStrain(const Tensor2& tensor) { m_epsp = tensor; }

Tensor2	LawDispatcher::getStressIncrement() const { return m_dsig; }
void	LawDispatcher::setStressIncrement(const Tensor2& tensor) { m_dsig = tensor; }

Tensor2	LawDispatcher::getStrainIncrement() const { return m_deps; }
void	LawDispatcher::setStrainIncrement(const Tensor2& tensor) { m_deps = tensor; }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LawDispatcher::computeStress(bool init, bool isRef) {

    auto f = m_compute_stress[m_law_type];

    if (f != nullptr)
    {
    	Tensor4 tangent_tensor = f(m_law_params,m_history_vars,m_sig,m_eps,m_epsp,m_dsig,m_deps,isRef);
    	if (init || isRef) m_tangent_tensor = tangent_tensor;
    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Tensor4 LawDispatcher::computeElastTensor(const Tensor2& sig) {

    auto f = m_compute_elast_tensor[m_law_type];

    if (f != nullptr)
        return f(m_law_params,sig);

    return Tensor4();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Tensor4 LawDispatcher::computeTangentTensor(const Tensor2& sig) {

    auto f = m_compute_tangent_tensor[m_law_type];

    if (f != nullptr)
      return f(m_law_params,m_history_vars,sig,m_deps);

    return Tensor4();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
bool LawDispatcher::initState(const Tensor2& sig)
{
    if (sig == Tensor2::zero()) return true;

    if (m_history_vars == nullptr)
      return true;

    initHistoryVars(m_history_vars);

    auto f = m_init_state[m_law_type];

    if (f != nullptr)
        return f(sig,m_history_vars);

    return false;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LawDispatcher::initHistoryVars(RealUniqueArray* history_vars)
{
    auto f = m_init_history_vars[m_law_type];

    if (f != nullptr)
        f(history_vars);

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LawDispatcher::readLawParams(RealUniqueArray* lawparams, Real lambda, Real mu, bool default_param, const String& name, Integer ilaw)
{
    if (lawparams == nullptr)
      return;

    if (lawparams->size() == 2){
      (*lawparams)[0] = lambda;
      (*lawparams)[1] = mu;
      return;
    }

    auto f = m_read_law_params[m_law_type];
    m_default = default_param;

    if (f != nullptr) {
      f(lawparams,lambda, mu, default_param, name, ilaw);

    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
RealUniqueArray LawDispatcher::initConsts(RealUniqueArray* law_params)
{
    auto f = m_init_consts[m_law_type];

    if (f != nullptr)
      return f(law_params);

    return {};
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ReadLawBlock(istream& is, Integer nblock) {

  char c[500];
  String str;
  bool stop{ false };
  String str_nblock = String::fromNumber(nblock);
  String sblock = String("<") + str_nblock + String(">");

  do {

    is.getline(c, 500); // read whole line
    std::istringstream iss(c);
    String token;
    while (iss >> token) {
      if (token.contains(sblock))
        stop = true;
    }
  } while(!stop);
}

