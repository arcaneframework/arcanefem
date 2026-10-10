// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* MGISConstitutiveLaw.cc                                      (C) 2000-2026 */
/*                                                                           */
/* Implementation of 'IConstitutiveLaw' using MGIS library.                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/PlatformUtils.h>
#include <arcane/accelerator/VariableViews.h>

#include "femutils/ArcaneFemFunctions.h"

#include "modules/elastoplasticity2/Elastoplasticity2Module.h"
#include "modules/elastoplasticity2/ElementMatrixHexQuad.h"
#include "modules/elastoplasticity2/ConstitutiveLawBase.h"
#include "modules/elastoplasticity2/MGISConstitutiveLaw_axl.h"

#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

#include <memory>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Indices of the components of the plane strain MGIS tensors with respect to
// the two-dimensional storage of the module ([xx, yy, sqrt(2)xy]).
constexpr Int32 MGIS_PS_XX = 0;
constexpr Int32 MGIS_PS_YY = 1;
constexpr Int32 MGIS_PS_ZZ = 2;
constexpr Int32 MGIS_PS_XY = 3;
constexpr Int32 MGIS_PS_SIZE = 4;

// Indices of the components of the two-dimensional tensors of the module.
constexpr Int32 C2D_XX = 0;
constexpr Int32 C2D_YY = 1;
constexpr Int32 C2D_XY = 2;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Adds the strain increment associated with a displacement gradient
 * to the gradients of an integration point of the MGIS material data
 * manager.
 *
 * Plane strain is assumed: the out-of-plane strain increment vanishes.
 * The shear component follows the MGIS (TFEL) convention in which the
 * symmetric tensors are stored with a sqrt(2) scaled shear term.
 */
void
_addStrainIncrementToGradients(mgis::real* gradients, const Real3x3& grad_DU)
{
  gradients[MGIS_PS_XX] += grad_DU(0, 0);
  gradients[MGIS_PS_YY] += grad_DU(1, 1);
  gradients[MGIS_PS_ZZ] += 0.;
  gradients[MGIS_PS_XY] += M_SQRT1_2 * (grad_DU(0, 1) + grad_DU(1, 0));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::ArcaneFem
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class MGISConstitutiveLaw
: public ArcaneMGISConstitutiveLawObject
{
 public:

  explicit MGISConstitutiveLaw(const ServiceBuildInfo& sbi)
  : ArcaneMGISConstitutiveLawObject(sbi)
  {
    // This is needed because Elastoplasticity module use it
    m_law_name = "VonMisesMGIS";
  };

 public:

  void initialize(const ConstitutiveLawInitInfo& law_info) override;
  void getMaterialProperties() override;
  void integrateAndSave() override;
  void restoreConvergedState() override;
  void commitInternalVariables() override;

 private:

  Real E = 0.0; // Youngs modulus
  Real nu = 0.0; // Poisson ratio
  Real sig0 = 0.0; // Yield strength
  Real Et = 0.0; // Tangent modulus

  std::unique_ptr<mgis::behaviour::Behaviour> m_mgis_behaviour;
  std::unique_ptr<mgis::behaviour::MaterialDataManager> m_mgis_data_manager;

 private:

  void _initializeMGIS();
  void _integrateAndSaveTria3Cpu();
  void _integrateAndSaveQuad4Cpu();
  void _integrateAndSaveQuad8Cpu();
  void _integrateAndSaveQuad9Cpu();
  void _integrateAndSaveResults();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_MGISCONSTITUTIVELAW(MGIS, MGISConstitutiveLaw);


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MGISConstitutiveLaw::
initialize(const ConstitutiveLawInitInfo& law_info)
{
  _initialize(law_info);

  E = options()->E(); // Youngs modulus
  nu = options()->nu(); // Poission ratio ν
  sig0 = options()->sig0(); // Yield Strength
  H = options()->H(); // Linear isotropic hardening modulus

  // By default, the hardening modulus is derived from E as in the
  // native VonMises law.
  Et = E / 100.;
  if (H < 0.)
    H = E * Et / (E - Et);

  m_sigma_zz_gp.reshape({ m_nGP });

  _initializeMGIS();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MGISConstitutiveLaw::
getMaterialProperties()
{
  // The hardening modulus has already been computed in
  // _initConstitutiveLaw() and the internal state variables are handled
  // by the MGIS material data manager.
  mu = (E / (2 * (1 + nu))); // lame parameter μ
  lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter λ

  if (mesh()->dimension() == 2) {

    m_C_elas_2d.fill(0.);
    m_C_elas_2d(0, 0) = lambda + 2. * mu;
    m_C_elas_2d(1, 1) = lambda + 2. * mu;
    m_C_elas_2d(2, 2) = 2. * mu;
    m_C_elas_2d(0, 1) = lambda;
    m_C_elas_2d(1, 0) = lambda;

    // Initialize constitutive history
    ENUMERATE_ (Cell, icell, allCells()) // TODO check if MDMeshVars provide initialisation method
    {
      for (Int16 iGP = 0; iGP < m_nGP; ++iGP) {
        m_sigma_gp(icell, iGP, 0) = 0.;
        m_sigma_gp(icell, iGP, 1) = 0.;
        m_sigma_gp(icell, iGP, 2) = 0.;
        m_sigma_zz_gp(icell, iGP) = 0.;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * @brief Initializes the MGIS binding for the Von Mises law.
 *
 * This method loads the MFront behaviour from the library given in the
 * options, checks that it is a standard small strain behaviour supporting
 * the plane strain hypothesis, allocates the material data manager over
 * all the Gauss points of the mesh and sets the material properties
 * (Young modulus, Poisson ratio, yield strength and hardening slope).
 */
void MGISConstitutiveLaw::
_initializeMGIS()
{
  info() << "[ArcaneFem-Info] Started module  _initMgisVonMises()";
  Real elapsedTime = platform::getRealTime();

  if (m_mgis_behaviour.get() || m_mgis_data_manager.get())
    ARCANE_FATAL("The MGIS binding is already initialized");

  if (mesh()->dimension() != 2)
    ARCANE_FATAL("The VonMisesMGIS law currently supports only 2D elements");

  const auto hypothesis = mgis::behaviour::Hypothesis::PLANESTRAIN;

  String library_name = options()->library();
  String behaviour_name = options()->behaviour(); // Behaviour name
  info() << "Initializing MGIS library='" << library_name << "' and behaviour '" << behaviour_name << "'";
  if (library_name.empty())
    ARCANE_FATAL("Invalid empty value for MGIS library name");

  try {
    {
      std::string l(library_name.toStdStringView());
      std::string b(behaviour_name.toStdStringView());
      m_mgis_behaviour = std::make_unique<mgis::behaviour::Behaviour>(mgis::behaviour::load(l, b, hypothesis));
    }

    auto& behaviour = *m_mgis_behaviour;

    if (behaviour.btype != mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR)
      ARCANE_FATAL("The behaviour '{0}' of library '{1}' is not a standard small strain behaviour",
                   behaviour_name, library_name);

    // Number of integration points, one set of Gauss points per cell.
    const Int32 nb_cell = mesh()->allCells().size();
    const mgis::size_type nb_integration_points = nb_cell * m_nGP;

    m_mgis_data_manager = std::make_unique<mgis::behaviour::MaterialDataManager>(behaviour, nb_integration_points);
    auto& data_manager = *m_mgis_data_manager;
    data_manager.allocateArrayOfTangentOperatorBlocks();

    // The behaviour is expected to exchange one 4x4 tangent operator block
    // (Stress,Strain) with the plane strain hypothesis.
    if (data_manager.K_stride != (MGIS_PS_SIZE * MGIS_PS_SIZE))
      ARCANE_FATAL("Unexpected tangent operator size for the behaviour '{0}' of library '{1}'",
                   behaviour_name, library_name);

    // The material properties are uniform and are set on the state at the
    // end of the time step before being propagated to the state at the
    // beginning of the time step by the update() call.
    mgis::behaviour::setMaterialProperty(data_manager.s1, "YoungModulus", E);
    mgis::behaviour::setMaterialProperty(data_manager.s1, "PoissonRatio", nu);
    mgis::behaviour::setMaterialProperty(data_manager.s1, "YieldStrength", sig0);
    mgis::behaviour::setMaterialProperty(data_manager.s1, "HardeningSlope", H);
    // The mechanical behaviours generated by MFront declare the temperature
    // as an external state variable. It is kept constant.
    if (mgis::behaviour::contains(behaviour.esvs, "Temperature"))
      mgis::behaviour::setExternalStateVariable(data_manager.s1, "Temperature", 293.15);
    mgis::behaviour::update(data_manager);
  }
  catch (const std::exception& e) {
    ARCANE_FATAL("Failed to initialize the MGIS behaviour '{0}' of library '{1}': {2}",
                 behaviour_name, library_name, String(e.what()));
  }

  info() << "[ArcaneFem-Info] MGIS behaviour '" << behaviour_name
         << "' of library '" << library_name
         << "' loaded for hypothesis " << mgis::behaviour::toString(hypothesis)
         << " on " << (mesh()->allCells().size() * m_nGP) << " integration points";

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "initialize-mgis-von-mises", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * @brief Restores the converged state of the previous time step in the
 * MGIS material data manager (the state at the end of the time step is
 * reset to the state at the beginning of the time step).
 */
void MGISConstitutiveLaw::
restoreConvergedState()
{
  ARCANE_CHECK_PTR(m_mgis_data_manager);
  mgis::behaviour::revert(*m_mgis_data_manager);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * @brief Commits the internal state variables after convergence of the
 * nonlinear solver (the state at the beginning of the time step is updated
 * with the converged state at the end of the time step).
 */
void MGISConstitutiveLaw::
commitInternalVariables()
{
  ARCANE_CHECK_PTR(m_mgis_data_manager);
  mgis::behaviour::update(*m_mgis_data_manager);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * @brief Integrates the Von Mises behaviour through MGIS and saves the
 * stresses and the consistent tangent operators on the Gauss point
 * variables of the module.
 *
 * The MGIS integration is performed on the CPU. The Gauss point variables
 * are then used by the assembly of the linear system, either on the CPU or
 * on an accelerator.
 */
void MGISConstitutiveLaw::
integrateAndSave()
{
  info() << "[ArcaneFem-Info] Started module  _integrateAndSaveConstitutiveLawVonMisesMgis()";
  Real elapsedTime = platform::getRealTime();

  ARCANE_CHECK_PTR(m_mgis_data_manager);

  // The behaviour is integrated from the state at the beginning of the time
  // step with the strain increment associated with the current estimate of
  // the displacement increment m_DUn.
  mgis::behaviour::revert(*m_mgis_data_manager);

  if (m_hex_quad_mesh) {
    if (m_nodes_per_cell == 4)
      _integrateAndSaveQuad4Cpu();
    else if (m_nodes_per_cell == 8)
      _integrateAndSaveQuad8Cpu();
    else
      _integrateAndSaveQuad9Cpu();
  }
  else {
    _integrateAndSaveTria3Cpu();
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "integrate-and-save-constitutive-law", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * @brief Integrates the Von Mises behaviour through MGIS on the TRIA3
 * elements of the mesh.
 */
void MGISConstitutiveLaw::
_integrateAndSaveTria3Cpu()
{
  auto& s1 = m_mgis_data_manager->s1;
  const auto gradients_stride = s1.gradients_stride;

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    // The strain increment is constant on a P1 triangle.
    const Real3x3 grad_DU = ArcaneFemFunctions::FeOperation2D::computeGradientTria3(cell, m_node_coord, m_DUn);

    for (Int16 iGP = 0; iGP < m_nGP; ++iGP) {
      const mgis::size_type gp = cell.localId() * m_nGP + iGP;
      _addStrainIncrementToGradients(s1.gradients.data() + gp * gradients_stride, grad_DU);
    }
  }

  _integrateAndSaveResults();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * @brief Integrates the Von Mises behaviour through MGIS on the QUAD4
 * elements of the mesh.
 */
void MGISConstitutiveLaw::
_integrateAndSaveQuad4Cpu()
{
  constexpr Real gp[2] = { -M_SQRT1_3, M_SQRT1_3 };

  auto& s1 = m_mgis_data_manager->s1;
  const auto gradients_stride = s1.gradients_stride;

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Int16 iGP = 0;
    for (Int16 ixi = 0; ixi < 2; ++ixi) {
      for (Int16 ieta = 0; ieta < 2; ++ieta) {
        const Real3x3 grad_DU = computeDisplacementGradientQuad4(cell, m_node_coord, m_DUn, gp[ixi], gp[ieta]);
        const mgis::size_type gpid = cell.localId() * m_nGP + iGP;
        _addStrainIncrementToGradients(s1.gradients.data() + gpid * gradients_stride, grad_DU);
        ++iGP;
      }
    }
  }

  _integrateAndSaveResults();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * @brief Integrates the Von Mises behaviour through MGIS on the QUAD8
 * elements of the mesh.
 */
void MGISConstitutiveLaw::
_integrateAndSaveQuad8Cpu()
{
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };

  auto& s1 = m_mgis_data_manager->s1;
  const auto gradients_stride = s1.gradients_stride;

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Int16 iGP = 0;
    for (Int16 ixi = 0; ixi < 3; ++ixi) {
      for (Int16 ieta = 0; ieta < 3; ++ieta) {
        const Real3x3 grad_DU = computeDisplacementGradientQuad8(cell, m_node_coord, m_DUn, gp[ixi], gp[ieta]);
        const mgis::size_type gpid = cell.localId() * m_nGP + iGP;
        _addStrainIncrementToGradients(s1.gradients.data() + gpid * gradients_stride, grad_DU);
        ++iGP;
      }
    }
  }

  _integrateAndSaveResults();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * @brief Integrates the Von Mises behaviour through MGIS on the QUAD9
 * elements of the mesh.
 */
void MGISConstitutiveLaw::
_integrateAndSaveQuad9Cpu()
{
  constexpr Real gp[3] = { -0.77459666924148337704, 0.0, 0.77459666924148337704 };

  auto& s1 = m_mgis_data_manager->s1;
  const auto gradients_stride = s1.gradients_stride;

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Int16 iGP = 0;
    for (Int16 ixi = 0; ixi < 3; ++ixi) {
      for (Int16 ieta = 0; ieta < 3; ++ieta) {
        const Real3x3 grad_DU = computeDisplacementGradientQuad9(cell, m_node_coord, m_DUn, gp[ixi], gp[ieta]);
        const mgis::size_type gpid = cell.localId() * m_nGP + iGP;
        _addStrainIncrementToGradients(s1.gradients.data() + gpid * gradients_stride, grad_DU);
        ++iGP;
      }
    }
  }

  _integrateAndSaveResults();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * @brief Integrates the behaviour over all the Gauss points and transfers
 * the thermodynamic forces and the tangent operator blocks towards the
 * Gauss point variables of the module.
 */
void MGISConstitutiveLaw::
_integrateAndSaveResults()
{
  mgis::behaviour::BehaviourIntegrationOptions options;
  options.integration_type =
  mgis::behaviour::IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR;

  const auto result = mgis::behaviour::integrate(*m_mgis_data_manager, options, dt);

  if (result.exit_status < 0)
    ARCANE_FATAL("The MGIS integration failed at integration point {1} : {2}",
                 (Integer)result.n,
                 String(result.error_message.empty() ? "unknown error" : result.error_message));

  if (result.exit_status == 0)
    info() << "WARNING: The MGIS integration succeeded but the results are unreliable";

  const auto& s1 = m_mgis_data_manager->s1;
  const auto forces_stride = s1.thermodynamic_forces_stride;
  const auto K_stride = m_mgis_data_manager->K_stride;
  const auto* const forces = s1.thermodynamic_forces.data();
  const auto* const K = m_mgis_data_manager->K.data();

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    for (Int16 iGP = 0; iGP < m_nGP; ++iGP) {
      const auto gp = mgis::size_type(cell.localId() * m_nGP + iGP);
      const auto* const sigma = forces + gp * forces_stride;
      const auto* const K_block = K + gp * K_stride;

      m_sigma_gp(cell, iGP, C2D_XX) = sigma[MGIS_PS_XX];
      m_sigma_gp(cell, iGP, C2D_YY) = sigma[MGIS_PS_YY];
      m_sigma_gp(cell, iGP, C2D_XY) = sigma[MGIS_PS_XY];
      m_sigma_zz_gp(cell, iGP) = sigma[MGIS_PS_ZZ];

      for (Int16 i = 0; i < 3; ++i) {
        const auto row = (i == C2D_XY) ? MGIS_PS_XY : i;
        for (Int16 j = 0; j < 3; ++j) {
          const auto column = (j == C2D_XY) ? MGIS_PS_XY : j;
          m_C_tang_gp(cell, iGP, i, j) = K_block[row * MGIS_PS_SIZE + column];
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::ArcaneFem

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
