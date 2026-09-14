// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* VonMisesLaw.h                                               (C) 2000-2026 */
/*                                                                           */
/* Contains functions to compute and assemble the Von Mises plasticity law   */
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
inline void FemModuleElastoplasticity::_restoreConvergedStateVonMises()
{
  ENUMERATE_ (Cell, icell, allCells())
  {
    Cell cell = *icell;

    for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
      m_sigma_gp(cell, iGP, 0) = m_sigma_old_gp(cell, iGP, 0);
      m_sigma_gp(cell, iGP, 1) = m_sigma_old_gp(cell, iGP, 1);
      m_sigma_gp(cell, iGP, 2) = m_sigma_old_gp(cell, iGP, 2);
      m_sigma_zz_gp(cell, iGP) = m_sigma_zz_old_gp(cell, iGP);
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
inline void FemModuleElastoplasticity::_commitInternalVariablesVonMises()
{
  ENUMERATE_ (Cell, icell, allCells())
  {
    Cell cell = *icell;

    for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
      m_sigma_old_gp(cell, iGP, 0) = m_sigma_gp(cell, iGP, 0);
      m_sigma_old_gp(cell, iGP, 1) = m_sigma_gp(cell, iGP, 1);
      m_sigma_old_gp(cell, iGP, 2) = m_sigma_gp(cell, iGP, 2);

      m_sigma_zz_old_gp(cell, iGP) = m_sigma_zz_gp(cell, iGP);
      m_p_old_gp(cell, iGP) += m_dp_gp(cell, iGP);
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria to update the
 * tangent material tensor matrix at each quadrature point for each
 * element
 *
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::_updateGlobalTangentMaterialTensorVonMises()
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
        _updateGlobalTangentMaterialTensorVonMisesTria3Gpu();
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
        _updateGlobalTangentMaterialTensorVonMisesTria3Cpu();
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

inline void FemModuleElastoplasticity::_updateGlobalTangentMaterialTensorVonMisesTria3Cpu()
{
  ENUMERATE_ (Cell, icell, allCells())
  {
    Cell cell = *icell;
    for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
      // --- compute_trial_state ---- //
      // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
      Real3x3 grad_DU = ArcaneFemFunctions::FeOperation2D::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_DUn);
      Real eps_xx = grad_DU(0, 0);
      Real eps_yy = grad_DU(1, 1);
      Real eps_xy = M_SQRT1_2 * (grad_DU(0, 1) + grad_DU(1, 0));

      Real sigma_trial_xx = m_sigma_old_gp(cell, iGP, 0) + m_C_elas_2d(0, 0) * eps_xx + m_C_elas_2d(0, 1) * eps_yy + m_C_elas_2d(0, 2) * eps_xy;
      Real sigma_trial_yy = m_sigma_old_gp(cell, iGP, 1) + m_C_elas_2d(1, 0) * eps_xx + m_C_elas_2d(1, 1) * eps_yy + m_C_elas_2d(1, 2) * eps_xy;
      Real sigma_trial_xy = m_sigma_old_gp(cell, iGP, 2) + m_C_elas_2d(2, 0) * eps_xx + m_C_elas_2d(2, 1) * eps_yy + m_C_elas_2d(2, 2) * eps_xy;

      Real sigma_trial_zz = m_sigma_zz_old_gp(cell, iGP) + lambda * (eps_xx + eps_yy);

      // Plane strain retains sigma_zz in the three-dimensional deviator.
      Real sigma_trial_mean = (sigma_trial_xx + sigma_trial_yy + sigma_trial_zz) / 3.0;

      Real dev_xx = sigma_trial_xx - sigma_trial_mean;
      Real dev_yy = sigma_trial_yy - sigma_trial_mean;
      Real dev_xy = sigma_trial_xy;

      Real dev_zz = sigma_trial_zz - sigma_trial_mean;

      Real sigma_eq_trial = math::sqrt(1.5 * (dev_xx * dev_xx + dev_yy * dev_yy + dev_zz * dev_zz + dev_xy * dev_xy) );

      // --- evaluate_yield_function ---- //
      Real yield_function = sigma_eq_trial - sig0 - H * m_p_old_gp(cell, iGP);
      Real yield_positive = (yield_function + math::abs(yield_function)) / 2.;
      m_dp_gp(cell, iGP) = yield_positive/ (3. * mu + H);
      Real plastic_switch = yield_positive / (math::abs(yield_function) + 1e-14 * sig0);

      // --- radial_return_update ---- //
      Real flowN_xx = plastic_switch * dev_xx / (sigma_eq_trial + 1e-14 * sig0);
      Real flowN_yy = plastic_switch * dev_yy / (sigma_eq_trial + 1e-14 * sig0);
      Real flowN_xy = plastic_switch * dev_xy / (sigma_eq_trial + 1e-14 * sig0);
      // Real flowN_zz = plastic_switch * dev_zz / (sigma_eq_trial + 1e-14 * sig0);

      Real beta = 3. * mu * m_dp_gp(cell, iGP) / (sigma_eq_trial + 1e-14 * sig0);

      // --- update_consistent_tangent ---- //
      Real sigma_xx = sigma_trial_xx - dev_xx * beta;
      Real sigma_yy = sigma_trial_yy - dev_yy * beta;
      Real sigma_xy = sigma_trial_xy - dev_xy * beta;

      Real sigma_zz = sigma_trial_zz - dev_zz * beta;

      m_sigma_gp(cell, iGP, 0) = sigma_xx;
      m_sigma_gp(cell, iGP, 1) = sigma_yy;
      m_sigma_gp(cell, iGP, 2) = sigma_xy;

      m_sigma_zz_gp(cell, iGP) = sigma_zz;

      Real tangentA = 3.* mu * (3. * mu / (3. * mu + H) - beta);

      m_C_tang_gp(cell, iGP, 0, 0) = m_C_elas_2d(0, 0) - tangentA * flowN_xx * flowN_xx - 4. * mu * beta / 3.;
      m_C_tang_gp(cell, iGP, 0, 1) = m_C_elas_2d(0, 1) - tangentA * flowN_xx * flowN_yy + 2. * mu * beta / 3.;
      m_C_tang_gp(cell, iGP, 0, 2) = m_C_elas_2d(0, 2) - tangentA * flowN_xx * flowN_xy;

      m_C_tang_gp(cell, iGP, 1, 0) = m_C_tang_gp(cell, iGP, 0, 1);
      m_C_tang_gp(cell, iGP, 1, 1) = m_C_elas_2d(1, 1) - tangentA * flowN_yy * flowN_yy - 4. * mu * beta / 3.;
      m_C_tang_gp(cell, iGP, 1, 2) = m_C_elas_2d(1, 2) - tangentA * flowN_yy * flowN_xy;

      m_C_tang_gp(cell, iGP, 2, 0) = m_C_tang_gp(cell, iGP, 0, 2);
      m_C_tang_gp(cell, iGP, 2, 1) = m_C_tang_gp(cell, iGP, 1, 2);
      m_C_tang_gp(cell, iGP, 2, 2) = m_C_elas_2d(2, 2) - tangentA * flowN_xy * flowN_xy - 2. * mu * beta;

    }
  }
}

inline void FemModuleElastoplasticity::_updateGlobalTangentMaterialTensorVonMisesTria3Gpu()
{
  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();

  auto command = Accelerator::makeCommand(queue);

  auto in_out_C_tang_gp = Accelerator::viewInOut(command, m_C_tang_gp);
  auto in_out_dp_gp = Accelerator::viewInOut(command, m_dp_gp);
  auto in_out_sigma_gp = Accelerator::viewInOut(command, m_sigma_gp);
  auto in_out_sigma_zz_gp = Accelerator::viewInOut(command, m_sigma_zz_gp);

  auto in_sigma_old_gp = Accelerator::viewIn(command, m_sigma_old_gp);
  auto in_sigma_zz_old_gp = Accelerator::viewIn(command, m_sigma_zz_old_gp);
  auto in_p_old_gp = Accelerator::viewIn(command, m_p_old_gp);

  auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
  auto in_DUn = Accelerator::viewIn(command, m_DUn);

  auto in_nGP = m_nGP;
  auto in_C_elas_2d = m_C_elas_2d;
  auto in_sig0 = sig0;
  auto in_H = H;
  auto in_mu = mu;


  command << RUNCOMMAND_ENUMERATE(CellLocalId, cell_lid, mesh()->allCells())
  {
    for (Int8 iGP = 0; iGP < in_nGP; ++iGP ) {
      // --- compute_trial_state ---- //
      // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
      Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);
      Real eps_xx = grad_DU(0, 0);
      Real eps_yy = grad_DU(1, 1);
      Real eps_xy = 0.70710678118654752440 * (grad_DU(0, 1) + grad_DU(1, 0));

      Real sigma_trial_xx = in_sigma_old_gp(cell_lid, iGP, 0) + in_C_elas_2d(0, 0) * eps_xx + in_C_elas_2d(0, 1) * eps_yy + in_C_elas_2d(0, 2) * eps_xy;
      Real sigma_trial_yy = in_sigma_old_gp(cell_lid, iGP, 1) + in_C_elas_2d(1, 0) * eps_xx + in_C_elas_2d(1, 1) * eps_yy + in_C_elas_2d(1, 2) * eps_xy;
      Real sigma_trial_xy = in_sigma_old_gp(cell_lid, iGP, 2) + in_C_elas_2d(2, 0) * eps_xx + in_C_elas_2d(2, 1) * eps_yy + in_C_elas_2d(2, 2) * eps_xy;

      Real sigma_trial_zz = in_sigma_zz_old_gp(cell_lid, iGP) + in_C_elas_2d(0, 1) * eps_yy + in_C_elas_2d(1, 0) * eps_xx;

      // Plane strain retains sigma_zz in the three-dimensional deviator.
      Real sigma_trial_mean = (sigma_trial_xx + sigma_trial_yy + sigma_trial_zz) / 3.0;

      Real dev_xx = sigma_trial_xx - sigma_trial_mean;
      Real dev_yy = sigma_trial_yy - sigma_trial_mean;
      Real dev_xy = sigma_trial_xy;

      Real dev_zz = sigma_trial_zz - sigma_trial_mean;

      Real sigma_eq_trial = math::sqrt(1.5 * (dev_xx * dev_xx + dev_yy * dev_yy + dev_zz * dev_zz + dev_xy * dev_xy) );

      // --- evaluate_yield_function ---- //
      Real yield_function = sigma_eq_trial - in_sig0 - in_H * in_p_old_gp(cell_lid, iGP);
      Real yield_positive = (yield_function + math::abs(yield_function)) / 2.;
      in_out_dp_gp(cell_lid, iGP) = yield_positive/ (3. * in_mu + in_H);
      Real plastic_switch = yield_positive / (math::abs(yield_function) + 1e-14 * in_sig0);

      // --- radial_return_update ---- //
      Real flowN_xx = plastic_switch * dev_xx / (sigma_eq_trial + 1e-14 * in_sig0);
      Real flowN_yy = plastic_switch * dev_yy / (sigma_eq_trial + 1e-14 * in_sig0);
      Real flowN_xy = plastic_switch * dev_xy / (sigma_eq_trial + 1e-14 * in_sig0);
      // Real flowN_zz = plastic_switch * dev_zz / (sigma_eq_trial + 1e-14 * in_sig0);

      Real beta = 3. * in_mu * in_out_dp_gp(cell_lid, iGP) / (sigma_eq_trial + 1e-14 * in_sig0);

      // --- update_consistent_tangent ---- //
      Real sigma_xx = sigma_trial_xx - dev_xx * beta;
      Real sigma_yy = sigma_trial_yy - dev_yy * beta;
      Real sigma_xy = sigma_trial_xy - dev_xy * beta;

      Real sigma_zz = sigma_trial_zz - dev_zz * beta;

      in_out_sigma_gp(cell_lid, iGP, 0) = sigma_xx;
      in_out_sigma_gp(cell_lid, iGP, 1) = sigma_yy;
      in_out_sigma_gp(cell_lid, iGP, 2) = sigma_xy;

      in_out_sigma_zz_gp(cell_lid, iGP) = sigma_zz;

      Real tangentA = 3.* in_mu * (3. * in_mu / (3. * in_mu + in_H) - beta);

      in_out_C_tang_gp(cell_lid, iGP, 0, 0) = in_C_elas_2d(0, 0) - tangentA * flowN_xx * flowN_xx - 4. * in_mu * beta / 3.;
      in_out_C_tang_gp(cell_lid, iGP, 0, 1) = in_C_elas_2d(0, 1) - tangentA * flowN_xx * flowN_yy + 2. * in_mu * beta / 3.;
      in_out_C_tang_gp(cell_lid, iGP, 0, 2) = in_C_elas_2d(0, 2) - tangentA * flowN_xx * flowN_xy;

      in_out_C_tang_gp(cell_lid, iGP, 1, 0) = in_out_C_tang_gp(cell_lid, iGP, 0, 1);
      in_out_C_tang_gp(cell_lid, iGP, 1, 1) = in_C_elas_2d(1, 1) - tangentA * flowN_yy * flowN_yy - 4. * in_mu * beta / 3.;
      in_out_C_tang_gp(cell_lid, iGP, 1, 2) = in_C_elas_2d(1, 2) - tangentA * flowN_yy * flowN_xy;

      in_out_C_tang_gp(cell_lid, iGP, 2, 0) = in_out_C_tang_gp(cell_lid, iGP, 0, 2);
      in_out_C_tang_gp(cell_lid, iGP, 2, 1) = in_out_C_tang_gp(cell_lid, iGP, 1, 2);
      in_out_C_tang_gp(cell_lid, iGP, 2, 2) = in_C_elas_2d(2, 2) - tangentA * flowN_xy * flowN_xy - 2. * in_mu * beta;
    }
  };
}