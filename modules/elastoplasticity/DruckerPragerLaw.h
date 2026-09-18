// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DruckerPragerLaw.h                                               (C) 2000-2026 */
/*                                                                           */
/* Contains functions to compute and assemble the Drucker-Prager             */
/* plasticity law                                                            */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Restores the initial or the converged solution state
 * from the previous time step solve for stress and material
 * tangent tensors at quadrature points
 *
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::_restoreConvergedStateDruckerPrager()
{
  ENUMERATE_ (Cell, icell, allCells())
  {
    Cell cell = *icell;

    for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
      m_eps_p_gp(cell, iGP, 0) = m_eps_p_old_gp(cell, iGP, 0);
      m_eps_p_gp(cell, iGP, 1) = m_eps_p_old_gp(cell, iGP, 1);
      m_eps_p_gp(cell, iGP, 2) = m_eps_p_old_gp(cell, iGP, 2);
      m_eps_p_zz_gp(cell, iGP) = m_eps_p_zz_old_gp(cell, iGP);
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Commits to the internal state variables after convergence of the
 * nonlinear solver for a given time step
 *
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::_commitInternalVariablesDruckerPrager()
{
  ENUMERATE_ (Cell, icell, allCells())
  {
    Cell cell = *icell;

    for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
      m_eps_p_old_gp(cell, iGP, 0) = m_eps_p_gp(cell, iGP, 0);
      m_eps_p_old_gp(cell, iGP, 1) = m_eps_p_gp(cell, iGP, 1);
      m_eps_p_old_gp(cell, iGP, 2) = m_eps_p_gp(cell, iGP, 2);
      m_eps_p_zz_old_gp(cell, iGP) = m_eps_p_zz_gp(cell, iGP);
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the DruckerPrager plasticity criteria to update the
 * tangent material tensor matrix at each quadrature point for each
 * element
 *
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::_updateGlobalTangentMaterialTensorDruckerPrager()
{
  auto use_gpu = options()->linearSystem.serviceName() == "HypreLinearSystem" ||
    options()->linearSystem.serviceName() == "PetscLinearSystem";

  if (use_gpu && m_use_gpu_functions) {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Not IMPLEMENTED");
      } else {
        _updateGlobalTangentMaterialTensorDruckerPragerTria3Gpu();
      }
    } else {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Not IMPLEMENTED");
      } else {
        ARCANE_FATAL("Not IMPLEMENTED");
      }
    }
  } else {
    if (mesh()->dimension() == 2) {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Not IMPLEMENTED");
      } else {
        _updateGlobalTangentMaterialTensorDruckerPragerTria3Cpu();
      }
    } else {
      if (m_hex_quad_mesh) {
        ARCANE_FATAL("Not IMPLEMENTED");
      } else {
        ARCANE_FATAL("Not IMPLEMENTED");
      }
    }
  }
}

ARCCORE_HOST_DEVICE void computeDruckerPragerLawAtGpBase(RealMatrix<3, 3>& C_tang_gp,
                                                          RealVector<3>& sigma_gp,
                                                          Real& sigma_zz_gp,
                                                          RealVector<3>& eps_p_gp,
                                                          Real& eps_p_zz_gp,
                                                          const Real3x3& grad_DU,
                                                          const Real3x3& grad_U,
                                                          const RealVector<3>& eps_p_old_gp,
                                                          const Real& eps_p_zz_old_gp,
                                                          const RealMatrix<3,3>& C_elas_2d,
                                                          const Real& in_bulk,
                                                          const Real& in_dpEta,
                                                          const Real& in_dpC,
                                                          const Real& in_mu)
{
  const Real SQRT1_2 = 0.70710678118654752440;
  const Real SQRT2	 = 1.41421356237309504880;

  Real eps_xx = grad_DU(0, 0);
  Real eps_yy = grad_DU(1, 1);
  Real eps_xy = SQRT1_2 * (grad_DU(0, 1) + grad_DU(1, 0));

  Real eps_trial_xx = grad_U(0, 0);
  Real eps_trial_yy = grad_U(1, 1);
  Real eps_trial_xy = SQRT1_2 * (grad_U(0, 1) + grad_U(1, 0));

  eps_xx += eps_trial_xx;
  eps_yy += eps_trial_yy;
  eps_xy += eps_trial_xy;
  Real eps_zz = 0.;

  // --- compute_elastic_trial_state ----
  // Elastic trial strain, deviator, trial stress, rho=||s||, and mean stress.

  // trial strain
  Real elastic_trial_strain_xx = eps_xx - eps_p_old_gp(0);
  Real elastic_trial_strain_yy = eps_yy - eps_p_old_gp(1);
  Real elastic_trial_strain_xy = eps_xy - eps_p_old_gp(2);
  Real elastic_trial_strain_zz = eps_zz - eps_p_zz_old_gp;

  // deviator
  Real mean_elastic_trial_strain = (elastic_trial_strain_xx + elastic_trial_strain_yy + elastic_trial_strain_zz) / 3.0;

  Real elastic_deviator_xx = elastic_trial_strain_xx - mean_elastic_trial_strain;
  Real elastic_deviator_yy = elastic_trial_strain_yy - mean_elastic_trial_strain;
  Real elastic_deviator_xy = elastic_trial_strain_xy;
  Real elastic_deviator_zz = elastic_trial_strain_zz - mean_elastic_trial_strain;

  // trial stress
  Real elastic_norm = math::sqrt(math::max(0.,elastic_trial_strain_xx * elastic_deviator_xx
                                            + elastic_trial_strain_yy * elastic_deviator_yy
                                            + elastic_trial_strain_zz * elastic_deviator_zz
                                            + elastic_trial_strain_xy * elastic_deviator_xy));
  Real rho_trial = 2. * in_mu * elastic_norm;
  Real pressure_trial = in_bulk * (elastic_trial_strain_xx + elastic_trial_strain_yy + elastic_trial_strain_zz);

  Real sigma_trial_xx = 2. * in_mu * elastic_deviator_xx + pressure_trial;
  Real sigma_trial_yy = 2. * in_mu * elastic_deviator_yy + pressure_trial;
  Real sigma_trial_xy = 2. * in_mu * elastic_deviator_xy;
  Real sigma_trial_zz = 2. * in_mu * elastic_deviator_zz + pressure_trial;


  // --- classify_elastic_smooth_or_apex_return ---- //
  Real denominator_apex = in_bulk * in_dpEta * in_dpEta;
  Real denominator_smooth = in_mu + denominator_apex;
  Real criterion1 = rho_trial / SQRT2 + in_dpEta * pressure_trial - in_dpC;
  Real criterion2 = in_dpEta * pressure_trial - denominator_apex * rho_trial / (in_mu * SQRT2) - in_dpC;
  Real plastic_switch = (criterion1>0. ? 1. : 0.);
  Real apex_switch = plastic_switch*(criterion2>0. ? 1. : 0.);
  Real smooth_switch = plastic_switch - apex_switch;
  Real lambda_smooth = smooth_switch * criterion1 / denominator_smooth;
  // Real lambda_apex = apex_switch * ( in_dpEta * pressure_trial - in_dpC) / denominator_apex;

  // --- update_stress_and_plastic_strain ---- //
  Real normal_xx = smooth_switch * elastic_deviator_xx / (elastic_norm + 1.e-30);
  Real normal_yy = smooth_switch * elastic_deviator_yy / (elastic_norm + 1.e-30);
  Real normal_xy = smooth_switch * elastic_deviator_xy / (elastic_norm + 1.e-30);
  Real normal_zz = smooth_switch * elastic_deviator_zz / (elastic_norm + 1.e-30);

  Real correction_xx = SQRT2 * in_mu * normal_xx + smooth_switch * in_bulk * in_dpEta;
  Real correction_yy = SQRT2 * in_mu * normal_yy + smooth_switch * in_bulk * in_dpEta;
  Real correction_xy = SQRT2 * in_mu * normal_xy;
  Real correction_zz = SQRT2 * in_mu * normal_zz + smooth_switch * in_bulk * in_dpEta;

  // --- update_consistent_tangent ---- //
  sigma_gp(0) = (1. - apex_switch) * sigma_trial_xx - lambda_smooth * correction_xx + apex_switch * in_dpC / in_dpEta;
  sigma_gp(1) = (1. - apex_switch) * sigma_trial_yy - lambda_smooth * correction_yy + apex_switch * in_dpC / in_dpEta;
  sigma_gp(2) = (1. - apex_switch) * sigma_trial_xy - lambda_smooth * correction_xy;
  sigma_zz_gp   = (1. - apex_switch) * sigma_trial_zz - lambda_smooth * correction_zz + apex_switch * in_dpC / in_dpEta;

  eps_p_gp(0) = eps_p_old_gp( 0) + lambda_smooth * (normal_xx / SQRT2 + in_dpEta/3.) + apex_switch * (eps_xx -in_dpC / (3. * in_bulk * in_dpEta) - eps_p_old_gp( 0));
  eps_p_gp(1) = eps_p_old_gp( 1) + lambda_smooth * (normal_yy / SQRT2 + in_dpEta/3.) + apex_switch * (eps_yy -in_dpC / (3. * in_bulk * in_dpEta) - eps_p_old_gp( 1));
  eps_p_gp(2) = eps_p_old_gp( 2) + lambda_smooth * normal_xy / SQRT2 + apex_switch * (eps_xy - eps_p_old_gp( 2));
  eps_p_zz_gp= eps_p_zz_old_gp + lambda_smooth * (normal_zz / SQRT2 + in_dpEta/3.) + apex_switch * (eps_zz -in_dpC / (3. * in_bulk * in_dpEta) - eps_p_zz_old_gp);

  Real curvature_factor = smooth_switch * 2. * SQRT2 * in_mu * in_mu * lambda_smooth / (rho_trial + 1.e-30);

  C_tang_gp(0, 0) = (1. - apex_switch) * (C_elas_2d(0, 0) - curvature_factor * ( 2./3. - normal_xx*normal_xx) - correction_xx*correction_xx / denominator_smooth);
  C_tang_gp(0, 1) = (1. - apex_switch) * (C_elas_2d(0, 1) - curvature_factor * ( -1./3. - normal_xx*normal_yy) - correction_xx*correction_yy / denominator_smooth);
  C_tang_gp(0, 2) = (1. - apex_switch) * (C_elas_2d(0, 2) - curvature_factor * ( 0. - normal_xx*normal_xy) - correction_xx*correction_xy / denominator_smooth);
  C_tang_gp(1, 0) = C_tang_gp(0, 1);
  C_tang_gp(1, 1) = (1. - apex_switch) * (C_elas_2d(1, 1) - curvature_factor * ( 2./3. - normal_yy*normal_yy) - correction_yy*correction_yy / denominator_smooth);
  C_tang_gp(1, 2) = (1. - apex_switch) * (C_elas_2d(1, 2) - curvature_factor * ( 0. - normal_yy*normal_xy) - correction_yy*correction_xy / denominator_smooth);
  C_tang_gp(2, 0) = C_tang_gp(0, 2);
  C_tang_gp(2, 1) = C_tang_gp(1, 2);
  C_tang_gp(2, 2) = (1. - apex_switch) * (C_elas_2d(2, 2) - curvature_factor * ( 1. - normal_xy*normal_xy) - correction_xy*correction_xy / denominator_smooth);
}

ARCCORE_HOST_DEVICE void computeTangentMaterialTensorDruckerPragerAtGp(RealMatrix<3, 3>& C_tang_gp,
                                                                        const Real3x3& grad_DU,
                                                                        const Real3x3& grad_U,
                                                                        const RealVector<3>& eps_p_old_gp,
                                                                        const Real& eps_p_zz_old_gp,
                                                                        const RealMatrix<3,3>& C_elas_2d,
                                                                        const Real& in_bulk,
                                                                        const Real& in_dpEta,
                                                                        const Real& in_dpC,
                                                                        const Real& in_mu)
{
  const Real SQRT1_2 = 0.70710678118654752440;
  const Real SQRT2	 = 1.41421356237309504880;

  Real eps_xx = grad_DU(0, 0);
  Real eps_yy = grad_DU(1, 1);
  Real eps_xy = SQRT1_2 * (grad_DU(0, 1) + grad_DU(1, 0));

  Real eps_trial_xx = grad_U(0, 0);
  Real eps_trial_yy = grad_U(1, 1);
  Real eps_trial_xy = SQRT1_2 * (grad_U(0, 1) + grad_U(1, 0));

  eps_xx += eps_trial_xx;
  eps_yy += eps_trial_yy;
  eps_xy += eps_trial_xy;
  Real eps_zz = 0.;

  // --- compute_elastic_trial_state ----
  // Elastic trial strain, deviator, trial stress, rho=||s||, and mean stress.

  // trial strain
  Real elastic_trial_strain_xx = eps_xx - eps_p_old_gp(0);
  Real elastic_trial_strain_yy = eps_yy - eps_p_old_gp(1);
  Real elastic_trial_strain_xy = eps_xy - eps_p_old_gp(2);
  Real elastic_trial_strain_zz = eps_zz - eps_p_zz_old_gp;

  // deviator
  Real mean_elastic_trial_strain = (elastic_trial_strain_xx + elastic_trial_strain_yy + elastic_trial_strain_zz) / 3.0;

  Real elastic_deviator_xx = elastic_trial_strain_xx - mean_elastic_trial_strain;
  Real elastic_deviator_yy = elastic_trial_strain_yy - mean_elastic_trial_strain;
  Real elastic_deviator_xy = elastic_trial_strain_xy;
  Real elastic_deviator_zz = elastic_trial_strain_zz - mean_elastic_trial_strain;

  // trial stress
  Real elastic_norm = math::sqrt(math::max(0.,elastic_trial_strain_xx * elastic_deviator_xx
                                            + elastic_trial_strain_yy * elastic_deviator_yy
                                            + elastic_trial_strain_zz * elastic_deviator_zz
                                            + elastic_trial_strain_xy * elastic_deviator_xy));
  Real rho_trial = 2. * in_mu * elastic_norm;
  Real pressure_trial = in_bulk * (elastic_trial_strain_xx + elastic_trial_strain_yy + elastic_trial_strain_zz);

  // --- classify_elastic_smooth_or_apex_return ---- //
  Real denominator_apex = in_bulk * in_dpEta * in_dpEta;
  Real denominator_smooth = in_mu + denominator_apex;
  Real criterion1 = rho_trial / SQRT2 + in_dpEta * pressure_trial - in_dpC;
  Real criterion2 = in_dpEta * pressure_trial - denominator_apex * rho_trial / (in_mu * SQRT2) - in_dpC;
  Real plastic_switch = (criterion1>0. ? 1. : 0.);
  Real apex_switch = plastic_switch*(criterion2>0. ? 1. : 0.);
  Real smooth_switch = plastic_switch - apex_switch;
  Real lambda_smooth = smooth_switch * criterion1 / denominator_smooth;
  // Real lambda_apex = apex_switch * ( in_dpEta * pressure_trial - in_dpC) / denominator_apex;

  // --- update_stress_and_plastic_strain ---- //
  Real normal_xx = smooth_switch * elastic_deviator_xx / (elastic_norm + 1.e-30);
  Real normal_yy = smooth_switch * elastic_deviator_yy / (elastic_norm + 1.e-30);
  Real normal_xy = smooth_switch * elastic_deviator_xy / (elastic_norm + 1.e-30);
  // Real normal_zz = smooth_switch * elastic_deviator_zz / (elastic_norm + 1.e-30);

  Real correction_xx = SQRT2 * in_mu * normal_xx + smooth_switch * in_bulk * in_dpEta;
  Real correction_yy = SQRT2 * in_mu * normal_yy + smooth_switch * in_bulk * in_dpEta;
  Real correction_xy = SQRT2 * in_mu * normal_xy;
  // Real correction_zz = SQRT2 * in_mu * normal_zz + smooth_switch * in_bulk * in_dpEta;

  // --- update_consistent_tangent ---- //
  Real curvature_factor = smooth_switch * 2. * SQRT2 * in_mu * in_mu * lambda_smooth / (rho_trial + 1.e-30);

  C_tang_gp(0, 0) = (1. - apex_switch) * (C_elas_2d(0, 0) - curvature_factor * ( 2./3. - normal_xx*normal_xx) - correction_xx*correction_xx / denominator_smooth);
  C_tang_gp(0, 1) = (1. - apex_switch) * (C_elas_2d(0, 1) - curvature_factor * ( -1./3. - normal_xx*normal_yy) - correction_xx*correction_yy / denominator_smooth);
  C_tang_gp(0, 2) = (1. - apex_switch) * (C_elas_2d(0, 2) - curvature_factor * ( 0. - normal_xx*normal_xy) - correction_xx*correction_xy / denominator_smooth);
  C_tang_gp(1, 0) = C_tang_gp(0, 1);
  C_tang_gp(1, 1) = (1. - apex_switch) * (C_elas_2d(1, 1) - curvature_factor * ( 2./3. - normal_yy*normal_yy) - correction_yy*correction_yy / denominator_smooth);
  C_tang_gp(1, 2) = (1. - apex_switch) * (C_elas_2d(1, 2) - curvature_factor * ( 0. - normal_yy*normal_xy) - correction_yy*correction_xy / denominator_smooth);
  C_tang_gp(2, 0) = C_tang_gp(0, 2);
  C_tang_gp(2, 1) = C_tang_gp(1, 2);
  C_tang_gp(2, 2) = (1. - apex_switch) * (C_elas_2d(2, 2) - curvature_factor * ( 1. - normal_xy*normal_xy) - correction_xy*correction_xy / denominator_smooth);
}

ARCCORE_HOST_DEVICE void computeStressAndInVarsDruckerPragerAtGp(RealVector<3>& sigma_gp,
                                                              Real& sigma_zz_gp,
                                                              RealVector<3>& eps_p_gp,
                                                              Real& eps_p_zz_gp,
                                                              const Real3x3& grad_DU,
                                                              const Real3x3& grad_U,
                                                              const RealVector<3>& eps_p_old_gp,
                                                              const Real& eps_p_zz_old_gp,
                                                              const Real& in_bulk,
                                                              const Real& in_dpEta,
                                                              const Real& in_dpC,
                                                              const Real& in_mu)
{
  const Real SQRT1_2 = 0.70710678118654752440;
  const Real SQRT2	 = 1.41421356237309504880;

  Real eps_xx = grad_DU(0, 0);
  Real eps_yy = grad_DU(1, 1);
  Real eps_xy = SQRT1_2 * (grad_DU(0, 1) + grad_DU(1, 0));

  Real eps_trial_xx = grad_U(0, 0);
  Real eps_trial_yy = grad_U(1, 1);
  Real eps_trial_xy = SQRT1_2 * (grad_U(0, 1) + grad_U(1, 0));

  eps_xx += eps_trial_xx;
  eps_yy += eps_trial_yy;
  eps_xy += eps_trial_xy;
  Real eps_zz = 0.;

  // --- compute_elastic_trial_state ----
  // Elastic trial strain, deviator, trial stress, rho=||s||, and mean stress.

  // trial strain
  Real elastic_trial_strain_xx = eps_xx - eps_p_old_gp(0);
  Real elastic_trial_strain_yy = eps_yy - eps_p_old_gp(1);
  Real elastic_trial_strain_xy = eps_xy - eps_p_old_gp(2);
  Real elastic_trial_strain_zz = eps_zz - eps_p_zz_old_gp;

  // deviator
  Real mean_elastic_trial_strain = (elastic_trial_strain_xx + elastic_trial_strain_yy + elastic_trial_strain_zz) / 3.0;

  Real elastic_deviator_xx = elastic_trial_strain_xx - mean_elastic_trial_strain;
  Real elastic_deviator_yy = elastic_trial_strain_yy - mean_elastic_trial_strain;
  Real elastic_deviator_xy = elastic_trial_strain_xy;
  Real elastic_deviator_zz = elastic_trial_strain_zz - mean_elastic_trial_strain;

  // trial stress
  Real elastic_norm = math::sqrt(math::max(0.,elastic_trial_strain_xx * elastic_deviator_xx
                                            + elastic_trial_strain_yy * elastic_deviator_yy
                                            + elastic_trial_strain_zz * elastic_deviator_zz
                                            + elastic_trial_strain_xy * elastic_deviator_xy));
  Real rho_trial = 2. * in_mu * elastic_norm;
  Real pressure_trial = in_bulk * (elastic_trial_strain_xx + elastic_trial_strain_yy + elastic_trial_strain_zz);

  Real sigma_trial_xx = 2. * in_mu * elastic_deviator_xx + pressure_trial;
  Real sigma_trial_yy = 2. * in_mu * elastic_deviator_yy + pressure_trial;
  Real sigma_trial_xy = 2. * in_mu * elastic_deviator_xy;
  Real sigma_trial_zz = 2. * in_mu * elastic_deviator_zz + pressure_trial;


  // --- classify_elastic_smooth_or_apex_return ---- //
  Real denominator_apex = in_bulk * in_dpEta * in_dpEta;
  Real denominator_smooth = in_mu + denominator_apex;
  Real criterion1 = rho_trial / SQRT2 + in_dpEta * pressure_trial - in_dpC;
  Real criterion2 = in_dpEta * pressure_trial - denominator_apex * rho_trial / (in_mu * SQRT2) - in_dpC;
  Real plastic_switch = (criterion1>0. ? 1. : 0.);
  Real apex_switch = plastic_switch*(criterion2>0. ? 1. : 0.);
  Real smooth_switch = plastic_switch - apex_switch;
  Real lambda_smooth = smooth_switch * criterion1 / denominator_smooth;
  // Real lambda_apex = apex_switch * ( in_dpEta * pressure_trial - in_dpC) / denominator_apex;

  // --- update_stress_and_plastic_strain ---- //
  Real normal_xx = smooth_switch * elastic_deviator_xx / (elastic_norm + 1.e-30);
  Real normal_yy = smooth_switch * elastic_deviator_yy / (elastic_norm + 1.e-30);
  Real normal_xy = smooth_switch * elastic_deviator_xy / (elastic_norm + 1.e-30);
  Real normal_zz = smooth_switch * elastic_deviator_zz / (elastic_norm + 1.e-30);

  Real correction_xx = SQRT2 * in_mu * normal_xx + smooth_switch * in_bulk * in_dpEta;
  Real correction_yy = SQRT2 * in_mu * normal_yy + smooth_switch * in_bulk * in_dpEta;
  Real correction_xy = SQRT2 * in_mu * normal_xy;
  Real correction_zz = SQRT2 * in_mu * normal_zz + smooth_switch * in_bulk * in_dpEta;

  // --- update_consistent_tangent ---- //
  sigma_gp(0) = (1. - apex_switch) * sigma_trial_xx - lambda_smooth * correction_xx + apex_switch * in_dpC / in_dpEta;
  sigma_gp(1) = (1. - apex_switch) * sigma_trial_yy - lambda_smooth * correction_yy + apex_switch * in_dpC / in_dpEta;
  sigma_gp(2) = (1. - apex_switch) * sigma_trial_xy - lambda_smooth * correction_xy;
  sigma_zz_gp   = (1. - apex_switch) * sigma_trial_zz - lambda_smooth * correction_zz + apex_switch * in_dpC / in_dpEta;

  eps_p_gp(0) = eps_p_old_gp( 0) + lambda_smooth * (normal_xx / SQRT2 + in_dpEta/3.) + apex_switch * (eps_xx -in_dpC / (3. * in_bulk * in_dpEta) - eps_p_old_gp( 0));
  eps_p_gp(1) = eps_p_old_gp( 1) + lambda_smooth * (normal_yy / SQRT2 + in_dpEta/3.) + apex_switch * (eps_yy -in_dpC / (3. * in_bulk * in_dpEta) - eps_p_old_gp( 1));
  eps_p_gp(2) = eps_p_old_gp( 2) + lambda_smooth * normal_xy / SQRT2 + apex_switch * (eps_xy - eps_p_old_gp( 2));
  eps_p_zz_gp= eps_p_zz_old_gp + lambda_smooth * (normal_zz / SQRT2 + in_dpEta/3.) + apex_switch * (eps_zz -in_dpC / (3. * in_bulk * in_dpEta) - eps_p_zz_old_gp);
  }

inline void FemModuleElastoplasticity::_updateGlobalTangentMaterialTensorDruckerPragerTria3Cpu()
{
  ENUMERATE_ (Cell, icell, allCells())
  {
    Cell cell = *icell;
    for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
      RealMatrix<3, 3> C_tang_gp;
      C_tang_gp(0, 0) = m_C_tang_gp(cell, iGP, 0, 0);
      C_tang_gp(0, 1) = m_C_tang_gp(cell, iGP, 0, 1);
      C_tang_gp(0, 2) = m_C_tang_gp(cell, iGP, 0, 2);
      C_tang_gp(1, 0) = m_C_tang_gp(cell, iGP, 1, 0);
      C_tang_gp(1, 1) = m_C_tang_gp(cell, iGP, 1, 1);
      C_tang_gp(1, 2) = m_C_tang_gp(cell, iGP, 1, 2);
      C_tang_gp(2, 0) = m_C_tang_gp(cell, iGP, 2, 0);
      C_tang_gp(2, 1) = m_C_tang_gp(cell, iGP, 2, 1);
      C_tang_gp(2, 2) = m_C_tang_gp(cell, iGP, 2, 2);
      RealVector<3> sigma_gp;
      sigma_gp(0) = m_sigma_gp(cell, iGP, 0);
      sigma_gp(1) = m_sigma_gp(cell, iGP, 1);
      sigma_gp(2) = m_sigma_gp(cell, iGP, 2);
      Real sigma_zz_gp = m_sigma_zz_gp(cell, iGP);
      RealVector<3> eps_p_gp;
      eps_p_gp(0) = m_eps_p_gp(cell, iGP, 0);
      eps_p_gp(1) = m_eps_p_gp(cell, iGP, 1);
      eps_p_gp(2) = m_eps_p_gp(cell, iGP, 2);
      Real eps_p_zz_gp = m_eps_p_zz_gp(cell, iGP);

      RealVector<3> eps_p_old_gp;
      eps_p_old_gp(0) = m_eps_p_old_gp(cell, iGP, 0);
      eps_p_old_gp(1) = m_eps_p_old_gp(cell, iGP, 1);
      eps_p_old_gp(2) = m_eps_p_old_gp(cell, iGP, 2);
      Real eps_p_zz_old_gp = m_eps_p_zz_old_gp(cell, iGP);

      // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
      Real3x3 grad_DU = ArcaneFemFunctions::FeOperation2D::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_DUn);
      Real3x3 grad_U = ArcaneFemFunctions::FeOperation2D::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_U);

      computeDruckerPragerLawAtGpBase(C_tang_gp,
                                      sigma_gp,sigma_zz_gp,
                                      eps_p_gp,eps_p_zz_gp,
                                      grad_DU, grad_U,
                                      eps_p_old_gp, eps_p_zz_old_gp,
                                      m_C_elas_2d, bulk, dpEta, dpC, mu);

      // update gp variables //
      m_C_tang_gp(cell, iGP, 0, 0) = C_tang_gp(0, 0);
      m_C_tang_gp(cell, iGP, 0, 1) = C_tang_gp(0, 1);
      m_C_tang_gp(cell, iGP, 0, 2) = C_tang_gp(0, 2);
      m_C_tang_gp(cell, iGP, 1, 0) = C_tang_gp(1, 0);
      m_C_tang_gp(cell, iGP, 1, 1) = C_tang_gp(1, 1);
      m_C_tang_gp(cell, iGP, 1, 2) = C_tang_gp(1, 2);
      m_C_tang_gp(cell, iGP, 2, 0) = C_tang_gp(2, 0);
      m_C_tang_gp(cell, iGP, 2, 1) = C_tang_gp(2, 1);
      m_C_tang_gp(cell, iGP, 2, 2) = C_tang_gp(2, 2);

      m_sigma_gp(cell, iGP, 0) = sigma_gp(0);
      m_sigma_gp(cell, iGP, 1) = sigma_gp(1);
      m_sigma_gp(cell, iGP, 2) = sigma_gp(2);
      m_sigma_zz_gp(cell, iGP) = sigma_zz_gp;

      m_eps_p_gp(cell, iGP, 0) = eps_p_gp(0);
      m_eps_p_gp(cell, iGP, 1) = eps_p_gp(1);
      m_eps_p_gp(cell, iGP, 2) = eps_p_gp(2);
      m_eps_p_zz_gp(cell, iGP) = eps_p_zz_gp;
    }
  }
}

inline void FemModuleElastoplasticity::_updateGlobalTangentMaterialTensorDruckerPragerTria3Gpu()
{
  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();
  auto command = makeCommand(acceleratorMng()->defaultQueue());

  auto in_out_C_tang_gp = Accelerator::viewInOut(command, m_C_tang_gp);
  auto in_out_sigma_gp = Accelerator::viewInOut(command, m_sigma_gp);
  auto in_out_sigma_zz_gp = Accelerator::viewInOut(command, m_sigma_zz_gp);
  auto in_out_eps_p_gp = Accelerator::viewInOut(command, m_eps_p_gp);
  auto in_out_eps_p_zz_gp = Accelerator::viewInOut(command, m_eps_p_zz_gp);

  auto in_eps_p_old_gp = Accelerator::viewIn(command, m_eps_p_old_gp);
  auto in_eps_p_zz_old_gp = Accelerator::viewIn(command, m_eps_p_zz_old_gp);

  auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
  auto in_DUn = Accelerator::viewIn(command, m_DUn);
  auto in_U = Accelerator::viewIn(command, m_U);

  auto in_nGP = m_nGP;
  auto C_elas_2d = m_C_elas_2d;
  auto in_bulk = bulk;
  auto in_dpEta = dpEta;
  auto in_dpC = dpC;
  auto in_mu = mu;

  command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh()->allCells())
  {
    for (Int8 iGP = 0; iGP < in_nGP; ++iGP ) {
      // read gp variables //
      RealMatrix<3, 3> C_tang_gp;
      C_tang_gp(0, 0) = in_out_C_tang_gp(cell_lid, iGP, 0, 0);
      C_tang_gp(0, 1) = in_out_C_tang_gp(cell_lid, iGP, 0, 1);
      C_tang_gp(0, 2) = in_out_C_tang_gp(cell_lid, iGP, 0, 2);
      C_tang_gp(1, 0) = in_out_C_tang_gp(cell_lid, iGP, 1, 0);
      C_tang_gp(1, 1) = in_out_C_tang_gp(cell_lid, iGP, 1, 1);
      C_tang_gp(1, 2) = in_out_C_tang_gp(cell_lid, iGP, 1, 2);
      C_tang_gp(2, 0) = in_out_C_tang_gp(cell_lid, iGP, 2, 0);
      C_tang_gp(2, 1) = in_out_C_tang_gp(cell_lid, iGP, 2, 1);
      C_tang_gp(2, 2) = in_out_C_tang_gp(cell_lid, iGP, 2, 2);
      RealVector<3> sigma_gp;
      sigma_gp(0) = in_out_sigma_gp(cell_lid, iGP, 0);
      sigma_gp(1) = in_out_sigma_gp(cell_lid, iGP, 1);
      sigma_gp(2) = in_out_sigma_gp(cell_lid, iGP, 2);
      Real sigma_zz_gp = in_out_sigma_zz_gp(cell_lid, iGP);
      RealVector<3> eps_p_gp;
      eps_p_gp(0) = in_out_eps_p_gp(cell_lid, iGP, 0);
      eps_p_gp(1) = in_out_eps_p_gp(cell_lid, iGP, 1);
      eps_p_gp(2) = in_out_eps_p_gp(cell_lid, iGP, 2);
      Real eps_p_zz_gp = in_out_eps_p_zz_gp(cell_lid, iGP);

      RealVector<3> eps_p_old_gp;
      eps_p_old_gp(0) = in_eps_p_old_gp(cell_lid, iGP, 0);
      eps_p_old_gp(1) = in_eps_p_old_gp(cell_lid, iGP, 1);
      eps_p_old_gp(2) = in_eps_p_old_gp(cell_lid, iGP, 2);
      Real eps_p_zz_old_gp = in_eps_p_zz_old_gp(cell_lid, iGP);

      // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
      Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);
      Real3x3 grad_U = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_U);

      computeDruckerPragerLawAtGpBase(C_tang_gp,
                                      sigma_gp,sigma_zz_gp,
                                      eps_p_gp,eps_p_zz_gp,
                                      grad_DU, grad_U,
                                      eps_p_old_gp, eps_p_zz_old_gp,
                                      C_elas_2d, in_bulk, in_dpEta, in_dpC, in_mu);

      // update gp variables //
      in_out_C_tang_gp(cell_lid, iGP, 0, 0) = C_tang_gp(0, 0);
      in_out_C_tang_gp(cell_lid, iGP, 0, 1) = C_tang_gp(0, 1);
      in_out_C_tang_gp(cell_lid, iGP, 0, 2) = C_tang_gp(0, 2);
      in_out_C_tang_gp(cell_lid, iGP, 1, 0) = C_tang_gp(1, 0);
      in_out_C_tang_gp(cell_lid, iGP, 1, 1) = C_tang_gp(1, 1);
      in_out_C_tang_gp(cell_lid, iGP, 1, 2) = C_tang_gp(1, 2);
      in_out_C_tang_gp(cell_lid, iGP, 2, 0) = C_tang_gp(2, 0);
      in_out_C_tang_gp(cell_lid, iGP, 2, 1) = C_tang_gp(2, 1);
      in_out_C_tang_gp(cell_lid, iGP, 2, 2) = C_tang_gp(2, 2);

      in_out_sigma_gp(cell_lid, iGP, 0) = sigma_gp(0);
      in_out_sigma_gp(cell_lid, iGP, 1) = sigma_gp(1);
      in_out_sigma_gp(cell_lid, iGP, 2) = sigma_gp(2);
      in_out_sigma_zz_gp(cell_lid, iGP) = sigma_zz_gp;

      in_out_eps_p_gp(cell_lid, iGP, 0) = eps_p_gp(0);
      in_out_eps_p_gp(cell_lid, iGP, 1) = eps_p_gp(1);
      in_out_eps_p_gp(cell_lid, iGP, 2) = eps_p_gp(2);
      in_out_eps_p_zz_gp(cell_lid, iGP) = eps_p_zz_gp;
    }
  };
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the DruckerPrager plasticity criteria to assemble the LHS
 * matrix for each element
 *
 */
/*---------------------------------------------------------------------------*/
RealMatrix<6, 6> FemModuleElastoplasticity::_computeLocalDruckerPragerElementMatrixTria3Cpu(Cell cell, bool assemble_elastic)
{
  Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  if (assemble_elastic) {
    return computeElementMatrixTria3Base(dxu, dyu, area, m_C_elas_2d);
  }

  Int8 iGP = 0;
  RealMatrix<3, 3> C_tang_gp;

  RealVector<3> eps_p_old_gp;
  eps_p_old_gp(0) = m_eps_p_old_gp(cell, iGP, 0);
  eps_p_old_gp(1) = m_eps_p_old_gp(cell, iGP, 1);
  eps_p_old_gp(2) = m_eps_p_old_gp(cell, iGP, 2);
  Real eps_p_zz_old_gp = m_eps_p_zz_old_gp(cell, iGP);

  // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
  Real3x3 grad_DU = ArcaneFemFunctions::FeOperation2D::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_DUn);
  Real3x3 grad_U = ArcaneFemFunctions::FeOperation2D::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_U);

  computeTangentMaterialTensorDruckerPragerAtGp(C_tang_gp,
                                                grad_DU, grad_U,
                                                eps_p_old_gp, eps_p_zz_old_gp,
                                                m_C_elas_2d, bulk, dpEta, dpC, mu);

  return computeElementMatrixTria3Base(dxu, dyu, area, C_tang_gp);
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE void computeLocalDruckerPragerStressAndInVarsTria3Gpu(CellLocalId cell_lid,
                    const IndexedCellNodeConnectivityView& cn_cv,
                    const Accelerator::VariableNodeReal3InView& in_node_coord,
                    const Accelerator::VariableNodeReal3InView& in_DUn,
                    const Accelerator::VariableNodeReal3InView& in_U,
                    const Accelerator::MeshMDVariableInOutView<Cell, double, ExtentsV<int, -1, -1>>& in_out_sigma_gp,
                    const Accelerator::MeshMDVariableInOutView<Cell, double, ExtentsV<int, -1>>& in_out_sigma_zz_gp,
                    const Accelerator::MeshMDVariableInOutView<Cell, double, ExtentsV<int, -1, -1>>& in_out_eps_p_gp,
                    const Accelerator::MeshMDVariableInOutView<Cell, double, ExtentsV<int, -1>>& in_out_eps_p_zz_gp,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1, -1>>& in_eps_p_old_gp,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1>>& in_eps_p_zz_old_gp,
                    const Real& in_bulk,
                    const Real& in_dpEta,
                    const Real& in_dpC,
                    const Real& in_mu)
{
  Int8 iGP = 0;
  RealVector<3> sigma_gp;
  sigma_gp(0) = in_out_sigma_gp(cell_lid, iGP, 0);
  sigma_gp(1) = in_out_sigma_gp(cell_lid, iGP, 1);
  sigma_gp(2) = in_out_sigma_gp(cell_lid, iGP, 2);
  Real sigma_zz_gp = in_out_sigma_zz_gp(cell_lid, iGP);
  RealVector<3> eps_p_gp;
  eps_p_gp(0) = in_out_eps_p_gp(cell_lid, iGP, 0);
  eps_p_gp(1) = in_out_eps_p_gp(cell_lid, iGP, 1);
  eps_p_gp(2) = in_out_eps_p_gp(cell_lid, iGP, 2);
  Real eps_p_zz_gp = in_out_eps_p_zz_gp(cell_lid, iGP);

  RealVector<3> eps_p_old_gp;
  eps_p_old_gp(0) = in_eps_p_old_gp(cell_lid, iGP, 0);
  eps_p_old_gp(1) = in_eps_p_old_gp(cell_lid, iGP, 1);
  eps_p_old_gp(2) = in_eps_p_old_gp(cell_lid, iGP, 2);
  Real eps_p_zz_old_gp = in_eps_p_zz_old_gp(cell_lid, iGP);

  // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
  Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);
  Real3x3 grad_U = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_U);

  computeStressAndInVarsDruckerPragerAtGp(sigma_gp,sigma_zz_gp,
                                       eps_p_gp,eps_p_zz_gp,
                                       grad_DU, grad_U,
                                       eps_p_old_gp, eps_p_zz_old_gp,
                                       in_bulk, in_dpEta, in_dpC, in_mu);

  // update gp variables //
  in_out_sigma_gp(cell_lid, iGP, 0) = sigma_gp(0);
  in_out_sigma_gp(cell_lid, iGP, 1) = sigma_gp(1);
  in_out_sigma_gp(cell_lid, iGP, 2) = sigma_gp(2);
  in_out_sigma_zz_gp(cell_lid, iGP) = sigma_zz_gp;

  in_out_eps_p_gp(cell_lid, iGP, 0) = eps_p_gp(0);
  in_out_eps_p_gp(cell_lid, iGP, 1) = eps_p_gp(1);
  in_out_eps_p_gp(cell_lid, iGP, 2) = eps_p_gp(2);
  in_out_eps_p_zz_gp(cell_lid, iGP) = eps_p_zz_gp;
}

ARCCORE_HOST_DEVICE RealMatrix<6, 6> computeLocalDruckerPragerElementMatrixTria3Gpu(CellLocalId cell_lid,
                    const IndexedCellNodeConnectivityView& cn_cv,
                    const Accelerator::VariableNodeReal3InView& in_node_coord,
                    const Accelerator::VariableNodeReal3InView& in_DUn,
                    const Accelerator::VariableNodeReal3InView& in_U,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1, -1>>& in_eps_p_old_gp,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1>>& in_eps_p_zz_old_gp,
                    const RealMatrix<3, 3>& C_elas_2d,
                    const Real& in_bulk,
                    const Real& in_dpEta,
                    const Real& in_dpC,
                    const Real& in_mu,
                    bool assemble_elastic)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);
  if (assemble_elastic) {
    return computeElementMatrixTria3Base(dxu, dyu, area, C_elas_2d);
  }

  Int8 iGP = 0;
  RealMatrix<3, 3> C_tang_gp;

  RealVector<3> eps_p_old_gp;
  eps_p_old_gp(0) = in_eps_p_old_gp(cell_lid, iGP, 0);
  eps_p_old_gp(1) = in_eps_p_old_gp(cell_lid, iGP, 1);
  eps_p_old_gp(2) = in_eps_p_old_gp(cell_lid, iGP, 2);
  Real eps_p_zz_old_gp = in_eps_p_zz_old_gp(cell_lid, iGP);

  // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
  Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);
  Real3x3 grad_U = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_U);

  computeTangentMaterialTensorDruckerPragerAtGp(C_tang_gp,
                                                grad_DU, grad_U,
                                                eps_p_old_gp, eps_p_zz_old_gp,
                                                C_elas_2d, in_bulk, in_dpEta, in_dpC, in_mu);

  return computeElementMatrixTria3Base(dxu, dyu, area, C_tang_gp);
}
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE RealMatrix<2, 6> computeLocalDruckerPragerElementVectorTria3Gpu(CellLocalId cell_lid,
                    const IndexedCellNodeConnectivityView& cn_cv,
                    const Accelerator::VariableNodeReal3InView& in_node_coord,
                    const Accelerator::VariableNodeReal3InView& in_DUn,
                    const Accelerator::VariableNodeReal3InView& in_U,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1, -1>>& in_eps_p_old_gp,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1>>& in_eps_p_zz_old_gp,
                    const RealMatrix<3, 3>& C_elas_2d,
                    const Real& in_bulk,
                    const Real& in_dpEta,
                    const Real& in_dpC,
                    const Real& in_mu,
                    bool assemble_elastic,
                    Int32 node_lid)
{
  Real3 dxu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyu = Arcane::FemUtils::Gpu::FeOperation2D::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);
  if (assemble_elastic) {
    return computeElementVectorTria3GpuBase(dxu, dyu, area, C_elas_2d, node_lid);
  }

  Int8 iGP = 0;
  RealMatrix<3, 3> C_tang_gp;

  RealVector<3> eps_p_old_gp;
  eps_p_old_gp(0) = in_eps_p_old_gp(cell_lid, iGP, 0);
  eps_p_old_gp(1) = in_eps_p_old_gp(cell_lid, iGP, 1);
  eps_p_old_gp(2) = in_eps_p_old_gp(cell_lid, iGP, 2);
  Real eps_p_zz_old_gp = in_eps_p_zz_old_gp(cell_lid, iGP);

  // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
  Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);
  Real3x3 grad_U = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_U);

  computeTangentMaterialTensorDruckerPragerAtGp(C_tang_gp,
                                                grad_DU, grad_U,
                                                eps_p_old_gp, eps_p_zz_old_gp,
                                                C_elas_2d, in_bulk, in_dpEta, in_dpC, in_mu);

  return computeElementVectorTria3GpuBase(dxu, dyu, area, C_tang_gp, node_lid);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the DruckerPrager plasticity criteria to update the
 * stress and internal state variables at each quadrature point for each
 * element
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::_updateStressAndInVarsDruckerPrager()
{
  if (mesh()->dimension() == 2) {
    if (m_hex_quad_mesh) {
      ARCANE_FATAL("Not IMPLEMENTED");
    } else {
      _updateStressAndInVarsDruckerPragerTria3Gpu();
    }
  } else {
    if (m_hex_quad_mesh) {
      ARCANE_FATAL("Not IMPLEMENTED");
    } else {
      ARCANE_FATAL("Not IMPLEMENTED");
    }
  }
}

inline void FemModuleElastoplasticity::_updateStressAndInVarsDruckerPragerTria3Gpu()
{
  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();
  auto command = makeCommand(acceleratorMng()->defaultQueue());

  auto in_out_sigma_gp = Accelerator::viewInOut(command, m_sigma_gp);
  auto in_out_sigma_zz_gp = Accelerator::viewInOut(command, m_sigma_zz_gp);
  auto in_out_eps_p_gp = Accelerator::viewInOut(command, m_eps_p_gp);
  auto in_out_eps_p_zz_gp = Accelerator::viewInOut(command, m_eps_p_zz_gp);

  auto in_eps_p_old_gp = Accelerator::viewIn(command, m_eps_p_old_gp);
  auto in_eps_p_zz_old_gp = Accelerator::viewIn(command, m_eps_p_zz_old_gp);

  auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
  auto in_DUn = Accelerator::viewIn(command, m_DUn);
  auto in_U = Accelerator::viewIn(command, m_U);

  auto in_nGP = m_nGP;
  auto in_bulk = bulk;
  auto in_dpEta = dpEta;
  auto in_dpC = dpC;
  auto in_mu = mu;

  command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh()->allCells())
  {
    for (Int8 iGP = 0; iGP < in_nGP; ++iGP ) {
      // read gp variables //
      RealVector<3> sigma_gp;
      sigma_gp(0) = in_out_sigma_gp(cell_lid, iGP, 0);
      sigma_gp(1) = in_out_sigma_gp(cell_lid, iGP, 1);
      sigma_gp(2) = in_out_sigma_gp(cell_lid, iGP, 2);
      Real sigma_zz_gp = in_out_sigma_zz_gp(cell_lid, iGP);
      RealVector<3> eps_p_gp;
      eps_p_gp(0) = in_out_eps_p_gp(cell_lid, iGP, 0);
      eps_p_gp(1) = in_out_eps_p_gp(cell_lid, iGP, 1);
      eps_p_gp(2) = in_out_eps_p_gp(cell_lid, iGP, 2);
      Real eps_p_zz_gp = in_out_eps_p_zz_gp(cell_lid, iGP);

      RealVector<3> eps_p_old_gp;
      eps_p_old_gp(0) = in_eps_p_old_gp(cell_lid, iGP, 0);
      eps_p_old_gp(1) = in_eps_p_old_gp(cell_lid, iGP, 1);
      eps_p_old_gp(2) = in_eps_p_old_gp(cell_lid, iGP, 2);
      Real eps_p_zz_old_gp = in_eps_p_zz_old_gp(cell_lid, iGP);

      // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
      Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);
      Real3x3 grad_U = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_U);

      computeStressAndInVarsDruckerPragerAtGp(sigma_gp,sigma_zz_gp,
                                           eps_p_gp,eps_p_zz_gp,
                                           grad_DU, grad_U,
                                           eps_p_old_gp, eps_p_zz_old_gp,
                                           in_bulk, in_dpEta, in_dpC, in_mu);

      // update gp variables //
      in_out_sigma_gp(cell_lid, iGP, 0) = sigma_gp(0);
      in_out_sigma_gp(cell_lid, iGP, 1) = sigma_gp(1);
      in_out_sigma_gp(cell_lid, iGP, 2) = sigma_gp(2);
      in_out_sigma_zz_gp(cell_lid, iGP) = sigma_zz_gp;

      in_out_eps_p_gp(cell_lid, iGP, 0) = eps_p_gp(0);
      in_out_eps_p_gp(cell_lid, iGP, 1) = eps_p_gp(1);
      in_out_eps_p_gp(cell_lid, iGP, 2) = eps_p_gp(2);
      in_out_eps_p_zz_gp(cell_lid, iGP) = eps_p_zz_gp;
    }
  };
}