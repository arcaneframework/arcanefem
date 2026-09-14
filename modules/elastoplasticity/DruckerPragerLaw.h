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
      m_sigma_gp(cell, iGP, 0) = m_sigma_old_gp(cell, iGP, 0);
      m_sigma_gp(cell, iGP, 1) = m_sigma_old_gp(cell, iGP, 1);
      m_sigma_gp(cell, iGP, 2) = m_sigma_old_gp(cell, iGP, 2);
      m_sigma_zz_gp(cell, iGP) = m_sigma_zz_old_gp(cell, iGP);
      
      m_eps_p_gp(cell, iGP, 0) = m_eps_p_old_gp(cell, iGP, 0);
      m_eps_p_gp(cell, iGP, 1) = m_eps_p_old_gp(cell, iGP, 1);
      m_eps_p_gp(cell, iGP, 2) = m_eps_p_old_gp(cell, iGP, 2);
      m_eps_p_zz_gp(cell, iGP) = m_eps_p_zz_old_gp(cell, iGP);
    }

    for (Int8 iGP = 0; iGP < m_nGP; ++iGP)
      for (Int8 ix = 0; ix < 3; ++ix)
        for (Int8 iy = 0; iy < 3; ++iy)
          m_C_tang_gp(cell, iGP, ix, iy) = m_C_elas_2d(ix, iy); // set tangent C equal to elastic C
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
      m_sigma_old_gp(cell, iGP, 0) = m_sigma_gp(cell, iGP, 0);
      m_sigma_old_gp(cell, iGP, 1) = m_sigma_gp(cell, iGP, 1);
      m_sigma_old_gp(cell, iGP, 2) = m_sigma_gp(cell, iGP, 2);
      m_sigma_zz_old_gp(cell, iGP) = m_sigma_zz_gp(cell, iGP);

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
    auto queue = subDomain()->acceleratorMng()->defaultQueue();
    auto mesh_ptr = mesh();
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

inline void FemModuleElastoplasticity::_updateGlobalTangentMaterialTensorDruckerPragerTria3Cpu()
{
  ENUMERATE_ (Cell, icell, allCells())
  {
    Cell cell = *icell;
    for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
      // --- compute_total_strain ---- //
      // Total plane-strain tensor at the current Newton iterate.

      // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
      Real3x3 grad_DU = ArcaneFemFunctions::FeOperation2D::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_DUn);
      Real eps_xx = grad_DU(0, 0);
      Real eps_yy = grad_DU(1, 1);
      Real eps_xy = M_SQRT1_2 * (grad_DU(0, 1) + grad_DU(1, 0));

      Real3x3 grad_U = ArcaneFemFunctions::FeOperation2D::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_U);
      Real eps_trial_xx = grad_U(0, 0);
      Real eps_trial_yy = grad_U(1, 1);
      Real eps_trial_xy = M_SQRT1_2 * (grad_U(0, 1) + grad_U(1, 0));

      eps_xx += eps_trial_xx;
      eps_yy += eps_trial_yy;
      eps_xy += eps_trial_xy;
      Real eps_zz = 0.;

      // --- compute_elastic_trial_state ----
      // Elastic trial strain, deviator, trial stress, rho=||s||, and mean stress.

      // trial strain
      Real elastic_trial_strain_xx = eps_xx - m_eps_p_old_gp(cell, iGP, 0);
      Real elastic_trial_strain_yy = eps_yy - m_eps_p_old_gp(cell, iGP, 1);
      Real elastic_trial_strain_xy = eps_xy - m_eps_p_old_gp(cell, iGP, 2);
      Real elastic_trial_strain_zz = eps_zz - m_eps_p_zz_old_gp(cell, iGP);

      // deviator
      Real mean_elastic_trial_strain = (elastic_trial_strain_xx + elastic_trial_strain_yy + elastic_trial_strain_zz) / 3.0;

      Real elastic_deviator_xx = elastic_trial_strain_xx - mean_elastic_trial_strain;
      Real elastic_deviator_yy = elastic_trial_strain_yy - mean_elastic_trial_strain;
      Real elastic_deviator_xy = elastic_trial_strain_xy;
      Real elastic_deviator_zz = elastic_trial_strain_zz - mean_elastic_trial_strain;

      // trial stress
      Real elastic_norm = math::sqrt(max(0.,elastic_trial_strain_xx * elastic_deviator_xx
                                                + elastic_trial_strain_yy * elastic_deviator_yy
                                                + elastic_trial_strain_zz * elastic_deviator_zz
                                                + elastic_trial_strain_xy * elastic_deviator_xy));
      Real rho_trial = 2. * mu * elastic_norm;
      Real pressure_trial = bulk * (elastic_trial_strain_xx + elastic_trial_strain_yy + elastic_trial_strain_zz);

      Real sigma_trial_xx = 2. * mu * elastic_deviator_xx + pressure_trial;
      Real sigma_trial_yy = 2. * mu * elastic_deviator_yy + pressure_trial;
      Real sigma_trial_xy = 2. * mu * elastic_deviator_xy;
      Real sigma_trial_zz = 2. * mu * elastic_deviator_zz + pressure_trial;


      // --- classify_elastic_smooth_or_apex_return ---- //
      Real denominator_apex = bulk * dpEta * dpEta;
      Real denominator_smooth = mu + denominator_apex;
      Real criterion1 = rho_trial / M_SQRT2 + dpEta * pressure_trial - dpC;
      Real criterion2 = dpEta * pressure_trial - denominator_apex * rho_trial / (mu * M_SQRT2) - dpC;
      Real plastic_switch = (criterion1>0. ? 1. : 0.);
      Real apex_switch = plastic_switch*(criterion2>0. ? 1. : 0.);
      Real smooth_switch = plastic_switch - apex_switch;
      Real lambda_smooth = smooth_switch * criterion1 / denominator_smooth;
      Real lambda_apex = apex_switch * ( dpEta * pressure_trial - dpC) / denominator_apex;

      // --- update_stress_and_plastic_strain ---- //
      Real normal_xx = smooth_switch * elastic_deviator_xx / (elastic_norm + 1.e-30);
      Real normal_yy = smooth_switch * elastic_deviator_yy / (elastic_norm + 1.e-30);
      Real normal_xy = smooth_switch * elastic_deviator_xy / (elastic_norm + 1.e-30);
      Real normal_zz = smooth_switch * elastic_deviator_zz / (elastic_norm + 1.e-30);

      Real correction_xx = M_SQRT2 * mu * normal_xx + smooth_switch * bulk * dpEta;
      Real correction_yy = M_SQRT2 * mu * normal_yy + smooth_switch * bulk * dpEta;
      Real correction_xy = M_SQRT2 * mu * normal_xy;
      Real correction_zz = M_SQRT2 * mu * normal_zz + smooth_switch * bulk * dpEta;

      Real sigma_xx = (1. - apex_switch) * sigma_trial_xx - lambda_smooth * correction_xx + apex_switch * dpC / dpEta;
      Real sigma_yy = (1. - apex_switch) * sigma_trial_yy - lambda_smooth * correction_yy + apex_switch * dpC / dpEta;
      Real sigma_xy = (1. - apex_switch) * sigma_trial_xy - lambda_smooth * correction_xy;
      Real sigma_zz = (1. - apex_switch) * sigma_trial_zz - lambda_smooth * correction_zz + apex_switch * dpC / dpEta;

      m_eps_p_gp(cell, iGP, 0) = m_eps_p_old_gp(cell, iGP, 0) + lambda_smooth * (normal_xx / M_SQRT2 + dpEta/3.) + apex_switch * (eps_xx -dpC / (3. * bulk * dpEta) - m_eps_p_old_gp(cell, iGP, 0));
      m_eps_p_gp(cell, iGP, 1) = m_eps_p_old_gp(cell, iGP, 1) + lambda_smooth * (normal_yy / M_SQRT2 + dpEta/3.) + apex_switch * (eps_yy -dpC / (3. * bulk * dpEta) - m_eps_p_old_gp(cell, iGP, 1));
      m_eps_p_gp(cell, iGP, 2) = m_eps_p_old_gp(cell, iGP, 2) + lambda_smooth * normal_xy / M_SQRT2 + apex_switch * (eps_xy - m_eps_p_old_gp(cell, iGP, 2));
      m_eps_p_zz_gp(cell, iGP) = m_eps_p_zz_old_gp(cell, iGP) + lambda_smooth * (normal_zz / M_SQRT2 + dpEta/3.) + apex_switch * (eps_zz -dpC / (3. * bulk * dpEta) - m_eps_p_zz_old_gp(cell, iGP));

      // --- update_consistent_tangent ---- //
      m_sigma_gp(cell, iGP, 0) = sigma_xx;
      m_sigma_gp(cell, iGP, 1) = sigma_yy;
      m_sigma_gp(cell, iGP, 2) = sigma_xy;

      m_sigma_zz_gp(cell, iGP) = sigma_zz;

      Real curvature_factor = smooth_switch * 2. * M_SQRT2 * mu * mu * lambda_smooth / (rho_trial + 1.e-30);

      m_C_tang_gp(cell, iGP, 0, 0) = (1. - apex_switch) * (m_C_elas_2d(0, 0) - curvature_factor * ( 2./3. - normal_xx*normal_xx) - correction_xx*correction_xx / denominator_smooth);
      m_C_tang_gp(cell, iGP, 0, 1) = (1. - apex_switch) * (m_C_elas_2d(0, 1) - curvature_factor * ( -1./3. - normal_xx*normal_yy) - correction_xx*correction_yy / denominator_smooth);
      m_C_tang_gp(cell, iGP, 0, 2) = (1. - apex_switch) * (m_C_elas_2d(0, 2) - curvature_factor * ( 0. - normal_xx*normal_xy) - correction_xx*correction_xy / denominator_smooth);

      m_C_tang_gp(cell, iGP, 1, 0) = m_C_tang_gp(cell, iGP, 0, 1);
      m_C_tang_gp(cell, iGP, 1, 1) = (1. - apex_switch) * (m_C_elas_2d(1, 1) - curvature_factor * ( 2./3. - normal_yy*normal_yy) - correction_yy*correction_yy / denominator_smooth);
      m_C_tang_gp(cell, iGP, 1, 2) = (1. - apex_switch) * (m_C_elas_2d(1, 2) - curvature_factor * ( 0. - normal_yy*normal_xy) - correction_yy*correction_xy / denominator_smooth);

      m_C_tang_gp(cell, iGP, 2, 0) = m_C_tang_gp(cell, iGP, 0, 2);
      m_C_tang_gp(cell, iGP, 2, 1) = m_C_tang_gp(cell, iGP, 1, 2);
      m_C_tang_gp(cell, iGP, 2, 2) = (1. - apex_switch) * (m_C_elas_2d(2, 2) - curvature_factor * ( 1. - normal_xy*normal_xy) - correction_xy*correction_xy / denominator_smooth);

    }
  }
}

inline void FemModuleElastoplasticity::_updateGlobalTangentMaterialTensorDruckerPragerTria3Gpu()
{

  ARCANE_FATAL("Not IMPLEMENTED");

  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();

  auto command = Accelerator::makeCommand(queue);

  auto in_out_C_tang_gp = Accelerator::viewInOut(command, m_C_tang_gp);
  // auto in_out_dp_gp = Accelerator::viewInOut(command, m_dp_gp);
  // auto in_out_sigma_gp = Accelerator::viewInOut(command, m_sigma_gp);
  // auto in_out_sigma_zz_gp = Accelerator::viewInOut(command, m_sigma_zz_gp);
  //
  // auto in_sigma_old_gp = Accelerator::viewIn(command, m_sigma_old_gp);
  // auto in_sigma_zz_old_gp = Accelerator::viewIn(command, m_sigma_zz_old_gp);
  // auto in_p_old_gp = Accelerator::viewIn(command, m_p_old_gp);

  auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
  auto in_DUn = Accelerator::viewIn(command, m_DUn);

  auto in_nGP = m_nGP;
  auto in_C_elas_2d = m_C_elas_2d;
  // auto in_sig0 = sig0;
  // auto in_H = H;
  // auto in_mu = mu;


  command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh()->allCells())
  {
    for (Int8 iGP = 0; iGP < in_nGP; ++iGP ) {
      Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);
      Real eps_xx = grad_DU(0, 0);
      Real eps_yy = grad_DU(1, 1);
      Real eps_xy = 0.70710678118654752440 * (grad_DU(0, 1) + grad_DU(1, 0));


      // _updateTangentMaterialTensorVM();
      in_out_C_tang_gp(cell_lid, iGP, 0, 0) = in_C_elas_2d(0, 0);
      in_out_C_tang_gp(cell_lid, iGP, 0, 1) = in_C_elas_2d(0, 1);
      in_out_C_tang_gp(cell_lid, iGP, 0, 2) = in_C_elas_2d(0, 2);
      in_out_C_tang_gp(cell_lid, iGP, 1, 0) = in_C_elas_2d(1, 0);
      in_out_C_tang_gp(cell_lid, iGP, 1, 1) = in_C_elas_2d(1, 1);
      in_out_C_tang_gp(cell_lid, iGP, 1, 2) = in_C_elas_2d(1, 2);
      in_out_C_tang_gp(cell_lid, iGP, 2, 0) = in_C_elas_2d(2, 0);
      in_out_C_tang_gp(cell_lid, iGP, 2, 1) = in_C_elas_2d(2, 1);
      in_out_C_tang_gp(cell_lid, iGP, 2, 2) = in_C_elas_2d(2, 2);
    }
  };
}