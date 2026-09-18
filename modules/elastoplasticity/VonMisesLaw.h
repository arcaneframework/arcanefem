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
  }
}
/*---------------------------------------------------------------------------*/

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
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria at a quadrature point
 * to compute the updated stress, the plastic strain increment and the
 * consistent tangent material tensor
 *
 * The update starts from an elastic trial state (built from the gradient of
 * the displacement increment and the state of the previous converged step),
 * evaluates the yield function (linear isotropic hardening) and applies a
 * radial return when the yield limit is exceeded. Plane strain is assumed:
 * sigma_zz is kept in the three-dimensional deviator.
 *
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE void computeMaterialTensorVonMisesLawAtGpBase(RealMatrix<3, 3>& C_tang_gp,
                                                                  RealVector<3>& sigma_gp,
                                                                  Real& sigma_zz_gp,
                                                                  Real& dp_gp,
                                                                  const Real3x3& grad_DU,
                                                                  const RealVector<3>& sigma_old_gp,
                                                                  const Real& sigma_zz_old_gp,
                                                                  const Real& p_old_gp,
                                                                  const RealMatrix<3,3>& C_elas_2d,
                                                                  const Real& in_sig0,
                                                                  const Real& in_H,
                                                                  const Real& in_mu)
{
  // --- compute_trial_state ---- //
  Real eps_xx = grad_DU(0, 0);
  Real eps_yy = grad_DU(1, 1);
  Real eps_xy = 0.70710678118654752440 * (grad_DU(0, 1) + grad_DU(1, 0));

  Real sigma_trial_xx = sigma_old_gp(0) + C_elas_2d(0, 0) * eps_xx + C_elas_2d(0, 1) * eps_yy + C_elas_2d(0, 2) * eps_xy;
  Real sigma_trial_yy = sigma_old_gp(1) + C_elas_2d(1, 0) * eps_xx + C_elas_2d(1, 1) * eps_yy + C_elas_2d(1, 2) * eps_xy;
  Real sigma_trial_xy = sigma_old_gp(2) + C_elas_2d(2, 0) * eps_xx + C_elas_2d(2, 1) * eps_yy + C_elas_2d(2, 2) * eps_xy;

  Real sigma_trial_zz = sigma_zz_old_gp + C_elas_2d(0, 1) * eps_yy + C_elas_2d(1, 0) * eps_xx;

  // Plane strain retains sigma_zz in the three-dimensional deviator.
  Real sigma_trial_mean = (sigma_trial_xx + sigma_trial_yy + sigma_trial_zz) / 3.0;

  Real dev_xx = sigma_trial_xx - sigma_trial_mean;
  Real dev_yy = sigma_trial_yy - sigma_trial_mean;
  Real dev_xy = sigma_trial_xy;

  Real dev_zz = sigma_trial_zz - sigma_trial_mean;

  Real sigma_eq_trial = math::sqrt(1.5 * (dev_xx * dev_xx + dev_yy * dev_yy + dev_zz * dev_zz + dev_xy * dev_xy) );

  // --- evaluate_yield_function ---- //
  Real yield_function = sigma_eq_trial - in_sig0 - in_H * p_old_gp;
  Real yield_positive = (yield_function + math::abs(yield_function)) / 2.;
  Real dp_gp_local = yield_positive/ (3. * in_mu + in_H);
  Real plastic_switch = yield_positive / (math::abs(yield_function) + 1e-14 * in_sig0);

  // --- radial_return_update ---- //
  Real flowN_xx = plastic_switch * dev_xx / (sigma_eq_trial + 1e-14 * in_sig0);
  Real flowN_yy = plastic_switch * dev_yy / (sigma_eq_trial + 1e-14 * in_sig0);
  Real flowN_xy = plastic_switch * dev_xy / (sigma_eq_trial + 1e-14 * in_sig0);
  // Real flowN_zz = plastic_switch * dev_zz / (sigma_eq_trial + 1e-14 * in_sig0);

  Real beta = 3. * in_mu * dp_gp_local / (sigma_eq_trial + 1e-14 * in_sig0);

  // --- update_consistent_tangent ---- //
  dp_gp = dp_gp_local;

  sigma_gp( 0) = sigma_trial_xx - dev_xx * beta;
  sigma_gp( 1) = sigma_trial_yy - dev_yy * beta;
  sigma_gp( 2) = sigma_trial_xy - dev_xy * beta;
  sigma_zz_gp = sigma_trial_zz - dev_zz * beta;

  Real tangentA = 3.* in_mu * (3. * in_mu / (3. * in_mu + in_H) - beta);

  C_tang_gp( 0, 0) = C_elas_2d(0, 0) - tangentA * flowN_xx * flowN_xx - 4. * in_mu * beta / 3.;
  C_tang_gp( 0, 1) = C_elas_2d(0, 1) - tangentA * flowN_xx * flowN_yy + 2. * in_mu * beta / 3.;
  C_tang_gp( 0, 2) = C_elas_2d(0, 2) - tangentA * flowN_xx * flowN_xy;
  C_tang_gp( 1, 0) = C_tang_gp(0, 1);
  C_tang_gp( 1, 1) = C_elas_2d(1, 1) - tangentA * flowN_yy * flowN_yy - 4. * in_mu * beta / 3.;
  C_tang_gp( 1, 2) = C_elas_2d(1, 2) - tangentA * flowN_yy * flowN_xy;
  C_tang_gp( 2, 0) = C_tang_gp(0, 2);
  C_tang_gp( 2, 1) = C_tang_gp(1, 2);
  C_tang_gp( 2, 2) = C_elas_2d(2, 2) - tangentA * flowN_xy * flowN_xy - 2. * in_mu * beta;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria at a quadrature point
 * to compute the consistent tangent material tensor only
 *
 * Same trial state and radial return as in
 * computeMaterialTensorVonMisesLawAtGpBase(), but the updated stress and
 * the plastic strain increment are not returned.
 *
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE void computeTangentMaterialTensorVonMisesAtGp(RealMatrix<3, 3>& C_tang_gp,
                                                                  const Real3x3& grad_DU,
                                                                  const RealVector<3>& sigma_old_gp,
                                                                  const Real& sigma_zz_old_gp,
                                                                  const Real& p_old_gp,
                                                                  const RealMatrix<3,3>& C_elas_2d,
                                                                  const Real& in_sig0,
                                                                  const Real& in_H,
                                                                  const Real& in_mu)
{
  // --- compute_trial_state ---- //
  Real eps_xx = grad_DU(0, 0);
  Real eps_yy = grad_DU(1, 1);
  Real eps_xy = 0.70710678118654752440 * (grad_DU(0, 1) + grad_DU(1, 0));

  Real sigma_trial_xx = sigma_old_gp(0) + C_elas_2d(0, 0) * eps_xx + C_elas_2d(0, 1) * eps_yy + C_elas_2d(0, 2) * eps_xy;
  Real sigma_trial_yy = sigma_old_gp(1) + C_elas_2d(1, 0) * eps_xx + C_elas_2d(1, 1) * eps_yy + C_elas_2d(1, 2) * eps_xy;
  Real sigma_trial_xy = sigma_old_gp(2) + C_elas_2d(2, 0) * eps_xx + C_elas_2d(2, 1) * eps_yy + C_elas_2d(2, 2) * eps_xy;

  Real sigma_trial_zz = sigma_zz_old_gp + C_elas_2d(0, 1) * eps_yy + C_elas_2d(1, 0) * eps_xx;

  // Plane strain retains sigma_zz in the three-dimensional deviator.
  Real sigma_trial_mean = (sigma_trial_xx + sigma_trial_yy + sigma_trial_zz) / 3.0;

  Real dev_xx = sigma_trial_xx - sigma_trial_mean;
  Real dev_yy = sigma_trial_yy - sigma_trial_mean;
  Real dev_xy = sigma_trial_xy;

  Real dev_zz = sigma_trial_zz - sigma_trial_mean;

  Real sigma_eq_trial = math::sqrt(1.5 * (dev_xx * dev_xx + dev_yy * dev_yy + dev_zz * dev_zz + dev_xy * dev_xy) );

  // --- evaluate_yield_function ---- //
  Real yield_function = sigma_eq_trial - in_sig0 - in_H * p_old_gp;
  Real yield_positive = (yield_function + math::abs(yield_function)) / 2.;
  Real dp_gp_local = yield_positive/ (3. * in_mu + in_H);
  Real plastic_switch = yield_positive / (math::abs(yield_function) + 1e-14 * in_sig0);

  // --- radial_return_update ---- //
  Real flowN_xx = plastic_switch * dev_xx / (sigma_eq_trial + 1e-14 * in_sig0);
  Real flowN_yy = plastic_switch * dev_yy / (sigma_eq_trial + 1e-14 * in_sig0);
  Real flowN_xy = plastic_switch * dev_xy / (sigma_eq_trial + 1e-14 * in_sig0);
  // Real flowN_zz = plastic_switch * dev_zz / (sigma_eq_trial + 1e-14 * in_sig0);

  Real beta = 3. * in_mu * dp_gp_local / (sigma_eq_trial + 1e-14 * in_sig0);

  // --- update_consistent_tangent ---- //
  Real tangentA = 3.* in_mu * (3. * in_mu / (3. * in_mu + in_H) - beta);

  C_tang_gp( 0, 0) = C_elas_2d(0, 0) - tangentA * flowN_xx * flowN_xx - 4. * in_mu * beta / 3.;
  C_tang_gp( 0, 1) = C_elas_2d(0, 1) - tangentA * flowN_xx * flowN_yy + 2. * in_mu * beta / 3.;
  C_tang_gp( 0, 2) = C_elas_2d(0, 2) - tangentA * flowN_xx * flowN_xy;
  C_tang_gp( 1, 0) = C_tang_gp(0, 1);
  C_tang_gp( 1, 1) = C_elas_2d(1, 1) - tangentA * flowN_yy * flowN_yy - 4. * in_mu * beta / 3.;
  C_tang_gp( 1, 2) = C_elas_2d(1, 2) - tangentA * flowN_yy * flowN_xy;
  C_tang_gp( 2, 0) = C_tang_gp(0, 2);
  C_tang_gp( 2, 1) = C_tang_gp(1, 2);
  C_tang_gp( 2, 2) = C_elas_2d(2, 2) - tangentA * flowN_xy * flowN_xy - 2. * in_mu * beta;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria at a quadrature point
 * to compute the updated stress and the plastic strain increment only
 *
 * Same trial state and radial return as in
 * computeMaterialTensorVonMisesLawAtGpBase(), but the consistent tangent
 * material tensor is not computed.
 *
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE void computeStressAndInVarsVonMisesAtGp(RealVector<3>& sigma_gp,
                                                            Real& sigma_zz_gp,
                                                            Real& dp_gp,
                                                            const Real3x3& grad_DU,
                                                            const RealVector<3>& sigma_old_gp,
                                                            const Real& sigma_zz_old_gp,
                                                            const Real& p_old_gp,
                                                            const RealMatrix<3,3>& C_elas_2d,
                                                            const Real& in_sig0,
                                                            const Real& in_H,
                                                            const Real& in_mu)
{
  // --- compute_trial_state ---- //
  Real eps_xx = grad_DU(0, 0);
  Real eps_yy = grad_DU(1, 1);
  Real eps_xy = 0.70710678118654752440 * (grad_DU(0, 1) + grad_DU(1, 0));

  Real sigma_trial_xx = sigma_old_gp(0) + C_elas_2d(0, 0) * eps_xx + C_elas_2d(0, 1) * eps_yy + C_elas_2d(0, 2) * eps_xy;
  Real sigma_trial_yy = sigma_old_gp(1) + C_elas_2d(1, 0) * eps_xx + C_elas_2d(1, 1) * eps_yy + C_elas_2d(1, 2) * eps_xy;
  Real sigma_trial_xy = sigma_old_gp(2) + C_elas_2d(2, 0) * eps_xx + C_elas_2d(2, 1) * eps_yy + C_elas_2d(2, 2) * eps_xy;

  Real sigma_trial_zz = sigma_zz_old_gp + C_elas_2d(0, 1) * eps_yy + C_elas_2d(1, 0) * eps_xx;

  // Plane strain retains sigma_zz in the three-dimensional deviator.
  Real sigma_trial_mean = (sigma_trial_xx + sigma_trial_yy + sigma_trial_zz) / 3.0;

  Real dev_xx = sigma_trial_xx - sigma_trial_mean;
  Real dev_yy = sigma_trial_yy - sigma_trial_mean;
  Real dev_xy = sigma_trial_xy;

  Real dev_zz = sigma_trial_zz - sigma_trial_mean;

  Real sigma_eq_trial = math::sqrt(1.5 * (dev_xx * dev_xx + dev_yy * dev_yy + dev_zz * dev_zz + dev_xy * dev_xy) );

  // --- evaluate_yield_function ---- //
  Real yield_function = sigma_eq_trial - in_sig0 - in_H * p_old_gp;
  Real yield_positive = (yield_function + math::abs(yield_function)) / 2.;
  Real dp_gp_local = yield_positive/ (3. * in_mu + in_H);

  Real beta = 3. * in_mu * dp_gp_local / (sigma_eq_trial + 1e-14 * in_sig0);

  // --- update_consistent_tangent ---- //
  dp_gp = dp_gp_local;

  sigma_gp( 0) = sigma_trial_xx - dev_xx * beta;
  sigma_gp( 1) = sigma_trial_yy - dev_yy * beta;
  sigma_gp( 2) = sigma_trial_xy - dev_xy * beta;
  sigma_zz_gp = sigma_trial_zz - dev_zz * beta;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria on the CPU to update the
 * tangent material tensor, the stress and the plastic strain increment at
 * each quadrature point for each TRIA3 element
 *
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::_updateGlobalTangentMaterialTensorVonMisesTria3Cpu()
{
  ENUMERATE_ (Cell, icell, allCells())
  {
    Cell cell = *icell;
    for (Int8 iGP = 0; iGP < m_nGP; ++iGP ) {
      // read gp variables //
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

      RealVector<3> sigma_old_gp;
      sigma_old_gp(0) = m_sigma_old_gp(cell, iGP, 0);
      sigma_old_gp(1) = m_sigma_old_gp(cell, iGP, 1);
      sigma_old_gp(2) = m_sigma_old_gp(cell, iGP, 2);
      Real sigma_zz_old_gp = m_sigma_zz_old_gp(cell, iGP);

      Real dp_gp = m_dp_gp(cell, iGP);
      Real p_old_gp = m_p_old_gp(cell, iGP);

      // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
      Real3x3 grad_DU = ArcaneFemFunctions::FeOperation2D::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_DUn);

      computeMaterialTensorVonMisesLawAtGpBase(C_tang_gp, sigma_gp, sigma_zz_gp, dp_gp, grad_DU,
                                                sigma_old_gp, sigma_zz_old_gp, p_old_gp,
                                                m_C_elas_2d, sig0,H, mu);

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

      m_dp_gp(cell, iGP) = dp_gp;
    }
  }
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria on the accelerator (GPU)
 * to update the tangent material tensor, the stress and the plastic strain
 * increment at each quadrature point for each TRIA3 element
 *
 */
/*---------------------------------------------------------------------------*/
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

      Real dp_gp = in_out_dp_gp(cell_lid, iGP);

      RealVector<3> sigma_old_gp;
      sigma_old_gp(0) = in_sigma_old_gp(cell_lid, iGP, 0);
      sigma_old_gp(1) = in_sigma_old_gp(cell_lid, iGP, 1);
      sigma_old_gp(2) = in_sigma_old_gp(cell_lid, iGP, 2);
      Real sigma_zz_old_gp = in_sigma_zz_old_gp(cell_lid, iGP);

      Real p_old_gp = in_p_old_gp(cell_lid, iGP);

      // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
      Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);

      computeMaterialTensorVonMisesLawAtGpBase(C_tang_gp, sigma_gp, sigma_zz_gp, dp_gp, grad_DU,
                                                sigma_old_gp, sigma_zz_old_gp, p_old_gp,
                                                in_C_elas_2d, in_sig0, in_H, in_mu);

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
      in_out_dp_gp(cell_lid, iGP) = dp_gp;
    }
  };
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria to assemble the LHS
 * matrix for each element
 *
 */
/*---------------------------------------------------------------------------*/
RealMatrix<6, 6> FemModuleElastoplasticity::_computeLocalVonMisesElementMatrixTria3Cpu(Cell cell, bool assemble_elastic)
{
  Real3 dxu = ArcaneFemFunctions::FeOperation2D::computeGradientXTria3(cell, m_node_coord);
  Real3 dyu = ArcaneFemFunctions::FeOperation2D::computeGradientYTria3(cell, m_node_coord);
  Real area = ArcaneFemFunctions::MeshOperation::computeAreaTria3(cell, m_node_coord);

  if (assemble_elastic) {
    return computeElementMatrixTria3Base(dxu, dyu, area, m_C_elas_2d);
  }

  Int8 iGP = 0;
  RealMatrix<3, 3> C_tang_gp;

  RealVector<3> sigma_old_gp;
  sigma_old_gp(0) = m_sigma_old_gp(cell, iGP, 0);
  sigma_old_gp(1) = m_sigma_old_gp(cell, iGP, 1);
  sigma_old_gp(2) = m_sigma_old_gp(cell, iGP, 2);
  Real sigma_zz_old_gp = m_sigma_zz_old_gp(cell, iGP);

  Real p_old_gp = m_p_old_gp(cell, iGP);

  // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
  Real3x3 grad_DU = ArcaneFemFunctions::FeOperation2D::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_DUn);

  computeTangentMaterialTensorVonMisesAtGp(C_tang_gp, grad_DU,
                                          sigma_old_gp, sigma_zz_old_gp, p_old_gp,
                                          m_C_elas_2d, sig0,H, mu);

  return computeElementMatrixTria3Base(dxu, dyu, area, C_tang_gp);
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria to compute the local LHS
 * matrix of a TRIA3 element from accelerator views (device counterpart of
 * _computeLocalVonMisesElementMatrixTria3Cpu)
 *
 * If assemble_elastic is true, the elastic material tensor is used instead
 * of the VonMises tangent material tensor.
 *
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE RealMatrix<6, 6> computeLocalVonMisesElementMatrixTria3Gpu(CellLocalId cell_lid,
                    const IndexedCellNodeConnectivityView& cn_cv,
                    const Accelerator::VariableNodeReal3InView& in_node_coord,
                    const Accelerator::VariableNodeReal3InView& in_DUn,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1, -1>>& in_sigma_old_gp,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1>>& in_sigma_zz_old_gp,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1>>& in_p_old_gp,
                    const RealMatrix<3, 3>& C_elas_2d,
                    const Real& in_sig0,
                    const Real& in_H,
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

  RealVector<3> sigma_old_gp;
  sigma_old_gp(0) = in_sigma_old_gp(cell_lid, iGP, 0);
  sigma_old_gp(1) = in_sigma_old_gp(cell_lid, iGP, 1);
  sigma_old_gp(2) = in_sigma_old_gp(cell_lid, iGP, 2);
  Real sigma_zz_old_gp = in_sigma_zz_old_gp(cell_lid, iGP);

  Real p_old_gp = in_p_old_gp(cell_lid, iGP);

  // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
  Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);

  computeTangentMaterialTensorVonMisesAtGp(C_tang_gp, grad_DU,
                                            sigma_old_gp, sigma_zz_old_gp, p_old_gp,
                                            C_elas_2d, in_sig0, in_H, in_mu);

  return computeElementMatrixTria3Base(dxu, dyu, area, C_tang_gp);
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria to compute the part of
 * the local LHS matrix of a TRIA3 element related to the node node_lid
 * (2x6 block), from accelerator views
 *
 * Node-wise variant of computeLocalVonMisesElementMatrixTria3Gpu(). If
 * assemble_elastic is true, the elastic material tensor is used instead of
 * the VonMises tangent material tensor.
 *
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE RealMatrix<2, 6> computeLocalVonMisesElementVectorTria3Gpu(CellLocalId cell_lid,
                    const IndexedCellNodeConnectivityView& cn_cv,
                    const Accelerator::VariableNodeReal3InView& in_node_coord,
                    const Accelerator::VariableNodeReal3InView& in_DUn,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1, -1>>& in_sigma_old_gp,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1>>& in_sigma_zz_old_gp,
                    const Accelerator::MeshMDVariableInView<Cell, double, ExtentsV<int, -1>>& in_p_old_gp,
                    const RealMatrix<3, 3>& C_elas_2d,
                    const Real& in_sig0,
                    const Real& in_H,
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

  RealVector<3> sigma_old_gp;
  sigma_old_gp(0) = in_sigma_old_gp(cell_lid, iGP, 0);
  sigma_old_gp(1) = in_sigma_old_gp(cell_lid, iGP, 1);
  sigma_old_gp(2) = in_sigma_old_gp(cell_lid, iGP, 2);
  Real sigma_zz_old_gp = in_sigma_zz_old_gp(cell_lid, iGP);

  Real p_old_gp = in_p_old_gp(cell_lid, iGP);

  // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
  Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);

  computeTangentMaterialTensorVonMisesAtGp(C_tang_gp, grad_DU,
                                            sigma_old_gp, sigma_zz_old_gp, p_old_gp,
                                            C_elas_2d, in_sig0, in_H, in_mu);

  return computeElementVectorTria3GpuBase(dxu, dyu, area, C_tang_gp, node_lid);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria to update the
 * stress and internal state variables at each quadrature point for each
 * element
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::_updateStressAndInVarsVonMises()
{
  if (mesh()->dimension() == 2) {
    if (m_hex_quad_mesh) {
      ARCANE_FATAL("Not IMPLEMENTED");
    } else {
      _updateStressAndInVarsVonMisesTria3Gpu();
    }
  } else {
    if (m_hex_quad_mesh) {
      ARCANE_FATAL("Not IMPLEMENTED");
    } else {
      ARCANE_FATAL("Not IMPLEMENTED");
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Applies the VonMises plasticity criteria on the accelerator (GPU)
 * to update the stress and the plastic strain increment at each quadrature
 * point for each TRIA3 element
 *
 */
/*---------------------------------------------------------------------------*/
inline void FemModuleElastoplasticity::_updateStressAndInVarsVonMisesTria3Gpu()
{
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();
  auto command = makeCommand(acceleratorMng()->defaultQueue());

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
      // read gp variables //
      RealVector<3> sigma_gp;
      sigma_gp(0) = in_out_sigma_gp(cell_lid, iGP, 0);
      sigma_gp(1) = in_out_sigma_gp(cell_lid, iGP, 1);
      sigma_gp(2) = in_out_sigma_gp(cell_lid, iGP, 2);
      Real sigma_zz_gp = in_out_sigma_zz_gp(cell_lid, iGP);

      Real dp_gp = in_out_dp_gp(cell_lid, iGP);

      RealVector<3> sigma_old_gp;
      sigma_old_gp(0) = in_sigma_old_gp(cell_lid, iGP, 0);
      sigma_old_gp(1) = in_sigma_old_gp(cell_lid, iGP, 1);
      sigma_old_gp(2) = in_sigma_old_gp(cell_lid, iGP, 2);
      Real sigma_zz_old_gp = in_sigma_zz_old_gp(cell_lid, iGP);

      Real p_old_gp = in_p_old_gp(cell_lid, iGP);

      // epsilon(DU) // NOTE: for nGP>1 it has to evaluated and interpolated at Gauss points
      Real3x3 grad_DU = Gpu::FeOperation2D::computeGradientTria3(cell_lid, cn_cv, in_node_coord, in_DUn);

      computeStressAndInVarsVonMisesAtGp(sigma_gp, sigma_zz_gp, dp_gp, grad_DU,
                                          sigma_old_gp, sigma_zz_old_gp, p_old_gp,
                                          in_C_elas_2d, in_sig0, in_H, in_mu);

      // update gp variables //
      in_out_sigma_gp(cell_lid, iGP, 0) = sigma_gp(0);
      in_out_sigma_gp(cell_lid, iGP, 1) = sigma_gp(1);
      in_out_sigma_gp(cell_lid, iGP, 2) = sigma_gp(2);
      in_out_sigma_zz_gp(cell_lid, iGP) = sigma_zz_gp;
      in_out_dp_gp(cell_lid, iGP) = dp_gp;
    }
  };
}
/*---------------------------------------------------------------------------*/