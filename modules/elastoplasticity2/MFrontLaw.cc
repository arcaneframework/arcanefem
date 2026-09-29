// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* MFrontLaw.cc                                                     (C) 2000-2026 */
/*                                                                           */
/* Generic MFront (MGIS) constitutive law. Integrates the compiled behaviour */
/* at each Gauss point to update the stress and the consistent tangent, using */
/* the MFrontGenericInterfaceSupport (MGIS) library. CPU-only, Tria3.         */
/*---------------------------------------------------------------------------*/

#include "modules/elastoplasticity2/Elastoplasticity2Module.h"

#include "femutils/ArcaneFemFunctions.h"

#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace arcane;
using namespace ArcaneFemFunctions;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::ArcaneFem
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace
{
  //! Map a module 2D stress/tangent slot to the MFront plane-strain stensor slot.
  //! Module 2D slots: 0 = xx, 1 = yy, 2 = sqrt2xy.
  //! MFront plane strain: 0 = xx, 1 = yy, 2 = zz, 3 = sqrt2xy.
  constexpr int _mgisSlot(int module_slot)
  {
    switch (module_slot) {
    case 0:
      return 0; // xx
    case 1:
      return 1; // yy
    case 2:
      return 3; // sqrt2xy
    }
    return module_slot;
  }

  //! Row-major 4x4 index of the MFront consistent tangent for a (row, col) module pair.
  constexpr int _mgisK(int module_row, int module_col)
  {
    return _mgisSlot(module_row) * 4 + _mgisSlot(module_col);
  }
} // namespace

/*---------------------------------------------------------------------------*/
/**
 * @brief Allocates one MGIS BehaviourData per Gauss point and seeds the
 *        material properties. Called once, when the material is initialised.
 */
/*---------------------------------------------------------------------------*/
void Elastoplasticity2Module::
_initMFrontGpData()
{
  const Int32 n_cells = allCells().size();
  const Int32 n_gp = n_cells * m_nGP;

  m_mfront_gp_data.clear();
  m_mfront_gp_data.reserve(n_gp);
  for (Int32 idx = 0; idx < n_gp; ++idx)
    m_mfront_gp_data.emplace_back(*m_mfront_behaviour);

  const Int32 n_axl = m_mfront_material_properties.size();
  for (auto& d : m_mfront_gp_data) {
    const Int32 n_mp = static_cast<Int32>(d.s0.material_properties.size());
    const Int32 n = (n_mp < n_axl) ? n_mp : n_axl;
    for (Int32 i = 0; i < n; ++i) {
      d.s0.material_properties[i] = m_mfront_material_properties[i];
      d.s1.material_properties[i] = m_mfront_material_properties[i];
    }
  }

  info() << "[ArcaneFem-Info] MFront: allocated " << n_gp
         << " BehaviourData (cells=" << n_cells << ", GPs/cell=" << m_nGP << ")";
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Restores the previous converged stress into the behaviour's old state.
 *
 * The per-Gauss-point BehaviourData persists across time steps: the old
 * stress is re-synced from the module's `sigma_old` history while the
 * internal/external state variables (plastic strain, etc.) are kept.
 */
/*---------------------------------------------------------------------------*/
void Elastoplasticity2Module::
_restoreConvergedStateMFront()
{
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    const Int32 numcell = cell.localId();
    const Int32 base = numcell * m_nGP;
    for (Int8 iGP = 0; iGP < m_nGP; ++iGP) {
      auto& d = m_mfront_gp_data[base + iGP];
      // MFront plane-strain forces: [0]=xx, [1]=yy, [2]=zz, [3]=sqrt2xy.
      d.s0.thermodynamic_forces[0] = m_sigma_old_gp(cell, iGP, 0);
      d.s0.thermodynamic_forces[1] = m_sigma_old_gp(cell, iGP, 1);
      d.s0.thermodynamic_forces[2] = m_sigma_zz_old_gp(cell, iGP);
      d.s0.thermodynamic_forces[3] = m_sigma_old_gp(cell, iGP, 2);
      // Incremental formulation: the total-strain reference stays zero.
      d.s0.gradients[0] = 0.;
      d.s0.gradients[1] = 0.;
      d.s0.gradients[2] = 0.;
      d.s0.gradients[3] = 0.;
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Integrates the MFront behaviour at every Gauss point and writes back
 *        the stress and the consistent tangent.
 *
 * For each Gauss point the strain increment is built from the displacement
 * increment (small strain, Tria3), the behaviour is integrated, and the
 * resulting stress / 3x3 tangent sub-block are stored in the MD variables.
 */
/*---------------------------------------------------------------------------*/
void Elastoplasticity2Module::
_integrateAndSaveConstitutiveLawMFront()
{
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    const Int32 numcell = cell.localId();
    const Int32 base = numcell * m_nGP;

    for (Int8 iGP = 0; iGP < m_nGP; ++iGP) {
      Real3x3 grad_DU = ArcaneFemFunctions::FemOperation2D::computeGradientTria3(
          cell, m_node_coord, m_DUn);
      const Real eps_xx = grad_DU(0, 0);
      const Real eps_yy = grad_DU(1, 1);
      const Real eps_xy = 0.70710678118654752440 * (grad_DU(0, 1) + grad_DU(1, 0));

      auto& d = m_mfront_gp_data[base + iGP];
      mgis::behaviour::update(d); // s1 = s0, K = 0
      d.K[0] = 4.0;               // INTEGRATION_CONSISTENT_TANGENT_OPERATOR
      d.s1.gradients[0] = eps_xx;
      d.s1.gradients[1] = eps_yy;
      d.s1.gradients[2] = 0.0;    // zz (plane strain)
      d.s1.gradients[3] = eps_xy; // sqrt2xy

      mgis::behaviour::BehaviourDataView v = mgis::behaviour::make_view(d);
      const int rc = mgis::behaviour::integrate(v, *m_mfront_behaviour);
      if (rc < 0) {
        ARCANE_FATAL("MFront behaviour integration failed at cell '{0}', GP '{1}' "
                     "(code '{2}'): '{3}'",
                     numcell, iGP, rc, d.error_message.c_str());
      }

      const auto& sf = d.s1.thermodynamic_forces;
      m_sigma_gp(cell, iGP, 0) = sf[0]; // xx
      m_sigma_gp(cell, iGP, 1) = sf[1]; // yy
      m_sigma_gp(cell, iGP, 2) = sf[3]; // sqrt2xy
      m_sigma_zz_gp(cell, iGP) = sf[2]; // zz

      for (Int8 i = 0; i < 3; ++i)
        for (Int8 j = 0; j < 3; ++j)
          m_C_tang_gp(cell, iGP, i, j) = d.K[_mgisK(i, j)];
    }
  }
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Commits the converged state for the MFront law.
 *
 * Advances the persistent BehaviourData (old stress / isv / esv) to the
 * converged values and updates the module's `sigma_old` history.
 */
/*---------------------------------------------------------------------------*/
void Elastoplasticity2Module::
_commitInternalVariablesMFront()
{
  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    const Int32 numcell = cell.localId();
    const Int32 base = numcell * m_nGP;

    for (Int8 iGP = 0; iGP < m_nGP; ++iGP) {
      auto& d = m_mfront_gp_data[base + iGP];
      // Advance the persistent state (keep the zero total-strain reference).
      for (int i = 0; i < 4; ++i)
        d.s0.thermodynamic_forces[i] = d.s1.thermodynamic_forces[i];
      d.s0.internal_state_variables = d.s1.internal_state_variables;
      d.s0.external_state_variables = d.s1.external_state_variables;
      d.s0.gradients[0] = d.s0.gradients[1] = d.s0.gradients[2] = d.s0.gradients[3] = 0.0;

      // Update the module stress history.
      m_sigma_old_gp(cell, iGP, 0) = m_sigma_gp(cell, iGP, 0);
      m_sigma_old_gp(cell, iGP, 1) = m_sigma_gp(cell, iGP, 1);
      m_sigma_old_gp(cell, iGP, 2) = m_sigma_gp(cell, iGP, 2);
      m_sigma_zz_old_gp(cell, iGP) = m_sigma_zz_gp(cell, iGP);
    }
  }
}

} // namespace Arcane::ArcaneFem