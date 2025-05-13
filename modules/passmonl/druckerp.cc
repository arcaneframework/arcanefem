// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-20245CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* druckerp.cc                                                 (C) 2022-2025 */
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

extern void ReadLawBlock(istream&, Integer);
extern Arcane::FemUtils::Tensor2 operator*(const Tensor4& tens, const Arcane::FemUtils::Tensor2& vector);
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
   * @brief Implementation of Drucker-Prager constitutive model
   */
/*---------------------------------------------------------------------------*/

//! Initialize intern useful constants
RealUniqueArray DruckPInitConsts(RealConstArrayView& law_params) {

  // auto E = law_params[0];
  // auto Nu = law_params[1];
  auto Lambda = law_params[0];
  auto Mu = law_params[1];
  auto phi = law_params[2]; // already in radians (converted when reading data)
  auto psi = law_params[3]; // already in radians
  auto cohes = law_params[4];
  auto sinphi = sin(phi);
  auto cosphi = cos(phi);
  auto sinpsi = sin(psi);
  auto sef = sqrt(3.) * (3. - sinphi);//compressive meridian

  RealUniqueArray consts(5);
//  consts[0] = E * Nu / (1. - 2. * Nu) / (1. + Nu); // Lame coefficient Lambda
//  consts[1] = E / 2. / (1. + Nu); // Lame coefficient Mu
  consts[0] = Lambda;
  consts[1] = Mu;
  consts[2] = 2. * sinphi / sef;//alfa
  consts[3] = 2. * sinpsi / sqrt(3.) / (3. - sinpsi);//alfpsi
  consts[4] = 6. * cohes * cosphi / sef;//xk

  return consts;
}

//! Computes elastic constitutive tensor
Tensor4 DruckPComputeElastTensor(RealConstArrayView& law_params, const Tensor2& /*sig*/) {
  RealUniqueArray consts = DruckPInitConsts(law_params);
  return {consts[0]/*Lambda*/,consts[1]/*Mu*/};
}

//! Computes tangent constitutive tensor
Tensor4 DruckPComputeTangentTensor(RealConstArrayView& law_params, RealArrayView& history_vars,
                                   const Tensor2& sig, const Tensor2& deps)
{
  RealUniqueArray consts = DruckPInitConsts(law_params);

  auto phi = law_params[2];
  auto psi = law_params[3];
  auto indaux = (int)law_params[6];
  auto lambda = consts[0];
  auto mu = consts[1];
  auto alfa = consts[2];
  auto alfpsi = consts[3];
  auto xk = consts[4];
  Tensor4 elast_tensor(lambda,mu);
  bool is_normal = (phi == psi);

  // Initializing tangent stiffness tensor with the elastic one
  Tensor4 tangent_tensor(elast_tensor);

  auto K = lambda + 2. * mu / 3.;
  auto depsv = deps.trace();
  Real lambdap{0.}; // plastic multiplier

  // 1st stress invariant
  Real I1 = sig.trace();
  Real3 I13(I1/3.,I1/3.,I1/3.);

  // Stress deviator
  Real3	sd = sig.get_diagonal(),// sii
  ss = sig.get_outdiagonal();//sij, j!=i

  sd -= I13;

  // 2nd deviatoric stress invariant
  Real	J2 = math::dot(ss,ss) + 0.5*math::dot(sd,sd),
       RJ2 = sqrt(J2);

  // Computing plastic multiplier:
  // lamdap = phiK*depsv + phiG*sij*depsij
  // with phiK = [3*K*alfa]/RH and phiG = G/RJ2/RH
  // In case of non associated law:
  // plastic potential: g = RJ2 + alfpsi*I1

  Real	H = mu + 9. * K * alfa * alfpsi, RH = sqrt(H),
       phiG = mu / RJ2 / RH,
       phiK = 3. * K * alfa / RH,
       psiK = 3. * K * alfpsi / RH;

  lambdap = phiK * depsv;
  lambdap += phiG * (math::dot(sd, deps.get_diagonal()) + math::dot(ss, deps.get_outdiagonal()));

  // if lambdap <= 0, nothing to do (remaining elastic): tangent_tensor = elast_tensor
  if (lambdap > 0.)
  {
    // Plastic yielding
    Real3	un(1,1,1),
    fh = phiG * sd + phiK * un,
    hs = phiG * ss,
    gh = phiG * sd + psiK * un;

    // tangent_tensor = symmetric if is_normal = true ou indaux = 0 (Sup = Slow)
    //     _          _
    //     |  D     Sup |
    // C = |            |
    //     |            |
    //     |  Slow   S  |
    //     |_          _|
    //
    // 3D tensorial products => matrix 3x3
    Real3x3	mat1 = math::prodTens(gh, fh), // D
    mat2 = math::prodTens(hs, hs), // S
    mat3 = math::prodTens(gh, hs), // Sup
    mat4 = math::prodTens(fh, hs); // Slow

    tangent_tensor[0] -= mat1; //D
    tangent_tensor[1] -= mat2; //S
    tangent_tensor[2] -= mat3; //Sup

    // if non associated law => indaux (user choice) to keep the tangent stiffness tensor symmetric or not
    // if yes, we assume Slow = Sup (so keeping upper part values only)
    if (is_normal || indaux == 0)
      mat4 = mat3;

    tangent_tensor[3] -= mat4; //Slow
    history_vars[0] = 1; // plastic state
  }
  else
    history_vars[0] = 0;// elastic state

  return tangent_tensor;
}

//! Initializes the vector of intern (history) variables
//! For Drücker-Prager law, the only variable is the plasticity indicator (Real value):
//! = 0 if material point is in elastic state, = 1 if in plastic state
RealUniqueArray DruckPInitHistoryVars(RealConstArrayView& history_vars)
{
  // Initializing plasticity indicator to 0
  return {history_vars};
}

//! Initializes the state (nothing to do for this law)
//! For Drücker-Prager law, it means initializing the plasticity indicator to 0 (elastic state)
bool DruckPInitState(const Tensor2& /*sig*/, RealArrayView& history_vars)
{
  history_vars[0] = 0;
  return true;
}

//! Read constitutive parameters from a file and initialize intern constants allowing to the material constitutive model type
RealUniqueArray DruckPReadLawParams(Real lambda, Real mu, bool default_param, const String& name, Integer ilaw)
{
  RealUniqueArray lawparams(7);
  std::filebuf MatFile;
  // Elastic parameters are taken from the general user data
  // Only the plastic parameters are read from input file
  lawparams[0] = lambda;
  lawparams[1] = mu;

  // =========================================================================================================
  // Drucker-Prager model parameters stored in the lawparams vector:
  // 0-Lambda=1st Lame coef.  1-Mu=2nd Lame coef.  2-phi(°)=friction angle 3-psi(°)=dilatancy angle
  // 4-cohesion (Pa) 5-incmax=max number of sub-increments for law integration
  // 6-indaux= indicator for tangent tensor type: 0 = elastic (default, symmetric),1 = plastic (unsymmetric)
  // =========================================================================================================

  bool is_file = (!name.empty() && MatFile.open(name.localstr(), ios::in) != nullptr);
  bool is_default = (default_param || !is_file);

  if (is_default) {
    // Taking default input law parameters
    lawparams[2] = 30.;
    lawparams[3] = 30.;
    lawparams[4] = 1.e3;
    lawparams[5] = 15;
    lawparams[6] = 0;

  } else {

    istream isRead(&MatFile);
    char c[500];

    // Find this law "block" in the file containing all models
    ReadLawBlock(isRead, ilaw);

    for (int i = 2; i < 7; i++)
      isRead >> lawparams[i];
    isRead.getline(c, 500); // "\n"

    MatFile.close();
  }

  const Real RAD = acos(-1.)/180.; //PI/180°
  lawparams[2] *= RAD;
  lawparams[3] *= RAD;

  return lawparams;
}

Tensor4 DruckPComputeStress(RealConstArrayView& law_params, RealArrayView& history_vars, Tensor2& sig, Tensor2& eps, Tensor2& epsp, Tensor2& dsig,
                                     const Tensor2& deps, bool /*isRef*/)
{
  RealUniqueArray consts = DruckPInitConsts(law_params);

  auto incmax = (int)law_params[5];
  auto indaux = (int)law_params[6];
  auto lambda = consts[0];
  auto mu = consts[1];
  auto alfa = consts[2];
  auto alfpsi = consts[3];
  auto xk = consts[4];
  bool is_plastic = (history_vars[0] == 1);
  auto elast_tensor = DruckPComputeElastTensor(law_params,Tensor2::zero());
  Tensor4 tangent_tensor;
  Real3	un(1,1,1);
  Real tol{1.0e-15};

  Tensor2 sign(sig);
  dsig = elast_tensor * deps;

  // stress estimation
  sig = sign + dsig;

  // 1st stress invariant
  Real	I1 = sig.trace();

  // stress deviator
  Real3	sd = sig.get_diagonal(),// sii
  ss = sig.get_outdiagonal();//sij, j!=i

  sd += (-I1/3.)*un;

  // 2nd deviatoric stress invariant
  auto	J2 = 0.5 * math::dot(sd,sd) + math::dot(ss,ss);
  auto    RJ2 = sqrt(J2);// square root of J2

  // Computing yield surface
  auto    fseuil = RJ2 + alfa*I1 - xk;

  if (fabs(fseuil) < tol)
    is_plastic = false;

  else  {

    // Plastic yielding
    is_plastic = true;

    // Evaluation of elastic & plastic parts for deformations:
    // dsig = dsigel + dsigp
    // dsigel = ratio*deps => sigel = sign + dsigel
    // dsigp = (1-ratio)*deps

    auto	I1n = sign.trace(),dI1 = dsig.trace();
    Real3	dsd = dsig.get_diagonal(),
          dss = dsig.get_outdiagonal(),
          sdn = sign.get_diagonal(),
          ssn = sign.get_outdiagonal();
    auto	dJ2 = math::dot(dss,dss),
         sds = math::dot(ssn,dss),
         J2n = math::dot(ssn,ssn);

    dsd += (-dI1/3.)*un;
    sdn += (-I1n/3.)*un;

    dJ2 += 0.5*math::dot(dsd,dsd);
    sds += 0.5*math::dot(sdn,dsd);
    J2n += 0.5*math::dot(sdn,sdn);

    // Computing ratio by solving: A*ratio*ratio + 2*B*ratio + E = 0
    Real	ratio = 0.,
         xkai = xk - alfa*I1n,
         A = dJ2 - pow(alfa*dI1,2),
         B = sds + alfa*dI1*xkai,
         E = J2n - xkai*xkai,
         Discrim = B*B - A*E;

    if (Discrim >= 0.)
    {
      if (fabs(A) <= tol && fabs(B) > tol)
      {
        ratio = -E/2./B;
      }
      else if (fabs(A) > tol)
        ratio = (-B + sqrt(Discrim))/A;
    }

    // Elastic part
    Tensor2 rdsig(ratio*dsig);

    sig = sign + rdsig;

    // Plastic part: sub-increments are used to project back on the yield surface
    auto ninc = incmax;// default nb max of sub-increments

    if (fabs(xk) < tol)
    {
      ninc = int(20.*sqrt(2.*fseuil/xk)) + 1;
      ninc = math::min(incmax,ninc);
    }
    Real dratio = (1. - ratio)/ninc;

    Tensor2 deps1(dratio*deps);

    if (deps1.norm() > tol)
    {
      // Loop on sub-increments
      for (int i = 0; i < ninc; i++)
      {
        tangent_tensor = DruckPComputeTangentTensor(law_params,history_vars,sig,deps1);

        sig += tangent_tensor * deps1;

//        if (is_plastic)
//        {
          I1 = sig.trace();
          sd = sig.get_diagonal() - (I1/3.)*un;
          ss = sig.get_outdiagonal();
          J2 = 0.5*math::dot(sd,sd) + math::dot(ss,ss);
          RJ2 = sqrt(J2);

          fseuil = RJ2 + alfa*I1 - xk;

          Real ft = fabs(fseuil);

          if (fabs(xk) > tol) {
            ft /= xk;
            tol = 0.005;
          }

          // Checking that we are on yield surface
          if (ft > tol)// Correcting as we are not yet
          {
            // Stress increment to project back on surface (f = fseuil):
            // dsigc = -f/(df/dsig)^2*(df/dsig)
            // Corrected stress: sig = sig - dsigc
            // if a = 0.5/RJ2:
            // (df/dsig) = alfa*I + a*sij (I = identite)
            // (df/dsig)^2 = 3*alfa^2 + 2*J2*a^2
            Real	alfa2 = alfa*alfa,
                 dfdsig2 = 3.*alfa2 + 0.5,
                 a = 0.5/RJ2;

            E = fseuil/dfdsig2;

            Tensor2 dsigc(E * (a * sd + alfa * un), E * a * ss);
            sig -= dsigc;
          }
          epsp += deps1;
        //}
      }
    }
  }
  history_vars[0] = (Real)is_plastic;
  return tangent_tensor;
}

