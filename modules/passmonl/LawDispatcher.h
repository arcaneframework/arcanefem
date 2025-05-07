// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-20245CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* LawDispatcher.h                                             (C) 2022-2025 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef PASSMO_LAWDISPATCHER_H
#define PASSMO_LAWDISPATCHER_H

#include "TypesNLdynamic.h"
#include "FemUtils.h"
#include "utilFEM.h"

//! Number of available constitutive models: for the moment, defining max 10
#define NB_LAW_TYPE 10

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
   * @brief Provides methods based on the Dispatcher mechanism available in Arcane,
   * allowing to implement constitutive laws implicitely
   */
/*---------------------------------------------------------------------------*/

class LawDispatcher {

public:
    LawDispatcher(TypesNLDynamic::eLawType law_type, bool default_param);
    LawDispatcher() = default;

public:
    Tensor4 computeElastTensor(const Tensor2& sig);
    Tensor4 computeTangentTensor(const Tensor2& sig);
    bool initState(const Tensor2& sig, RealConstArrayView& history_vars);
    RealUniqueArray initHistoryVars(RealConstArrayView& history_vars);
    void computeStress(bool is_converge);
    RealUniqueArray initConsts(RealConstArrayView& law_params);
    RealUniqueArray readLawParams(Real lambda, Real mu, bool default_param, const String& name, Integer ilaw);
    RealUniqueArray updateHistoryVars();
//    Tensor4 updateTangentTensor();

    [[nodiscard]] Tensor2	getStress() const;
    void	setStress(const Tensor2&);

    [[nodiscard]] Tensor2	getStrain() const;
    void	setStrain(const Tensor2&);

    [[nodiscard]] Tensor2	getPlasticStrain() const;
    void	setPlasticStrain(const Tensor2&);

    [[nodiscard]] Tensor2	getStressIncrement() const;
    void	setStressIncrement(const Tensor2&);

    [[nodiscard]] Tensor2	getStrainIncrement() const;
    void	setStrainIncrement(const Tensor2&);

    void  setLambda(Real lambda) { m_Lambda = lambda; }
    void  setMu(Real mu) { m_Mu = mu; }
    void  setName(const String& name) { m_name = name; }
    void  setDefault(bool is_default) { m_default = is_default; }

private:
    std::function<Tensor4(RealConstArrayView& law_params, RealArrayView& history_vars, Tensor2& sig, Tensor2& eps, Tensor2& epsp, Tensor2& dsig,
            const Tensor2& deps, bool is_converge)> m_compute_stress[NB_LAW_TYPE];
    std::function<Tensor4(RealConstArrayView& law_params, const Tensor2& sig)> m_compute_elast_tensor[NB_LAW_TYPE];
    std::function<Tensor4(RealConstArrayView& law_params, RealArrayView& history_vars, const Tensor2& sig, const Tensor2& deps)> m_compute_tangent_tensor[NB_LAW_TYPE];
    std::function<bool(const Tensor2& sig, RealArrayView& history_vars)> m_init_state[NB_LAW_TYPE];
    std::function<RealUniqueArray(RealConstArrayView& history_vars)> m_init_history_vars[NB_LAW_TYPE];
    std::function<RealUniqueArray(Real lambda, Real mu, bool default_param, const String& name, Integer ilaw)> m_read_law_params[NB_LAW_TYPE];
    std::function<RealUniqueArray(RealConstArrayView& law_params)> m_init_consts[NB_LAW_TYPE];

    Tensor2 m_sig{};
    Tensor2 m_eps{};
    Tensor2 m_dsig{};
    Tensor2 m_epsp{};
    Tensor2 m_deps{};
    Tensor4 m_elast_tensor{};
    Tensor4 m_tangent_tensor{};
    Real m_Lambda{0.}, m_Mu{0.};
    RealUniqueArray m_law_params{};
    RealUniqueArray m_history_vars{};
    RealUniqueArray m_law_consts{};
    TypesNLDynamic::eLawType m_law_type{TypesNLDynamic::HOOKE};
    String m_name{};
    bool m_default{true};
};

#endif //PASSMO_LAWDISPATCHER_H
