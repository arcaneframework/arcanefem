// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-20245CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* NLDModule.h                                                 (C) 2022-2025 */
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
#include "FemUtils.h"
#include "utilFEM.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;
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
    auto E = law_params[0];
    auto Nu = law_params[1];

    RealUniqueArray consts(2);
    consts[0] = E * Nu / (1. - 2. * Nu) / (1. + Nu); // Lame coefficient Lambda
    consts[1] = E / 2. / (1. + Nu); // Lame coefficient Mu

    return consts;
}

//! Computes elastic constitutive tensor
Tensor4 HookeComputeElastTensor(RealConstArrayView& law_params, const Tensor2& /*sig*/)
{
    RealUniqueArray consts = HookeInitConsts(law_params);
    return {consts[0]/*Lambda*/,consts[1]/*Mu*/};
}

//! Computes tangent constitutive tensor
Tensor4 HookeComputeTangentTensor(RealConstArrayView& law_params, RealArrayView& /*history_params*/, const Tensor2& /*sig*/, const Tensor2& /*deps*/)
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
RealUniqueArray HookeReadLawParams(const String& name)
{
#ifdef DEBUG
    if (name.empty())
	{
		info() << "\n Error filename for materials not found "
		        "--> Process stopped in HookeReadLawParams\n";
		return RealUniqueArray(0);
	}
#endif

        std::filebuf MatFile;

        if (MatFile.open(name.localstr(), ios::in) == nullptr)
            return {};

        istream isRead(&MatFile);
        char    c[500];

        // =================================================================================
        // Hooke's model parameters :
        // 0-E=Young modulus  1-Nu=Poisson ratio
        // =================================================================================
        RealUniqueArray lawparams(2);
        for (int i = 0; i < 2; i++)
            isRead >> lawparams[i];
        isRead.getline(c, 500); // "\n"

        MatFile.close();
        return lawparams;
}

Tensor4 HookeComputeStress(RealConstArrayView& law_params, RealArrayView& /*history_vars*/, Tensor2& sig, Tensor2& eps, Tensor2& /*epsp*/, Tensor2& dsig,
                        const Tensor2& deps, bool /*is_converge*/)
{
  Tensor4 elast_tensor = HookeComputeElastTensor(law_params,sig);
	sig += elast_tensor*deps;
	eps += deps;
  return elast_tensor;
}

