// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* NLDynamicModule.cc                                      (C) 2022-2025     */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "arcane/MathUtils.h"
#include <arcane/utils/MultiArray2.h>
#include "arcane/utils/ArgumentException.h"
#include <arcane/IParallelMng.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/IIOMng.h>
#include <arcane/CaseTable.h>
#include "IDoFLinearSystemFactory.h"
#include "ArcaneFemFunctions.h"
#include "NLDynamicModule.h"
#include "LawDispatcher.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
NLDynamicModule::NLDynamicModule(const ModuleBuildInfo& mbi)
: ArcaneNLDynamicObject(mbi)
, m_dofs_on_nodes(mbi.subDomain()->traceMng())
, m_gauss_on_cells(mbi.subDomain()->traceMng())
{
  ICaseMng *cm = mbi.subDomain()->caseMng();
  cm->setTreatWarningAsError(true);
  cm->setAllowUnkownRootElelement(false);
}

VersionInfo NLDynamicModule::versionInfo() const {
  return {1, 0, 0};
}

NLDynamicModule::~NLDynamicModule()
{
  for( const CaseTableInfo&  t : m_traction_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_dc_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_sacc_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_sdispl_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_sforce_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_svel_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_acc_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_displ_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_force_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_vel_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_ain_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_vin_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_uin_case_table_list )
    delete t.case_table;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
startInit(){

  info() << "Module Elastodynamic INIT";

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());

  ninteg = options()->getGaussNint();
  gravity.x = options()->getGx();
  gravity.y = options()->getGy();
  gravity.z = options()->getGz();

  auto opt_penalty = options()->enforceDirichletMethod().lower();
  if (opt_penalty.contains("penalty") || opt_penalty.contains("weak")) {
    penalty = options()->getPenalty();
  }
  gamma = options()->getGamma();
  beta = options()->getBeta();
  alfam = options()->getAlfam();
  alfaf = options()->getAlfaf();
  auto dt = options()->getDeltat();
  m_global_deltat = dt;
  dt2 = dt * dt;
  auto tf = options()->getFinalTime();
  m_global_final_time = tf;
  auto t = options()->getStart();
  m_global_time = t;
  linop_nstep = options()->getLinopNstep();
  auto szType = options()->initElastType().lower();
  if (szType.contains("young")) elast_type = TypesNLDynamic::YoungNu;
  else if (szType.contains("lame")) elast_type = TypesNLDynamic::Lame;
  else if (szType.contains("vel")) elast_type = TypesNLDynamic::Veloc;
  else {

    info() << "init-elast-type keywords must include (not case dependent):\n"
           << "  - Young\n"
           << "  - Lame\n"
           << "  - Velocity or Vel\n";
    ARCANE_FATAL("Type for elastic properties is undefined!");
  }
  algo_type = options()->getNonlinAlgoType();
  is_linear = (algo_type == TypesNLDynamic::Linear);

  if (options()->hasItemax())
    ite_max = options()->getItemax();

  if (options()->hasUtol())
    utol = options()->getUtol();

  if (options()->hasFtol())
    ftol = options()->getFtol();

  if (options()->hasEtol())
    etol = options()->getEtol();
  else
    etol = utol * ftol;

  auto nsteps = (int)((tf - t)/dt);
  if (!is_linear && linop_nstep > nsteps) keep_constop = true;

  is_alfa_method = options()->alfa_method();
  if (is_alfa_method) {
    gamma = 0.5 + alfaf - alfam;
    beta = 0.5*pow(0.5 + gamma,2);
  }
  else{
    alfam = 0.;
    alfaf = 0.;
  }

  analysis_type = options()->getAnalysisType();
  if (analysis_type == TypesNLDynamic::ThreeD)
    NDIM = 3;
  else
    NDIM = 2;

  auto dirichletMethod = options()->enforceDirichletMethod();
  auto dirichletMethodl = dirichletMethod.lower();

  if (!dirichletMethodl.contains("penalty") && !dirichletMethodl.contains("weak")
      && !dirichletMethodl.contains("rowelim")
      && !dirichletMethodl.contains("rowcolumnelim")) {
    info() << "Dirichlet boundary condition via "
           << dirichletMethod << " is not supported \n"
           << "enforce-Dirichlet-method only supports (not case dependent):\n"
           << "  - Penalty\n"
           << "  - WeakPenalty or Weak\n"
           << "  - RowElimination or RowElim\n"
           << "  - RowColumnElimination or RowColumnElim\n";

    ARCANE_FATAL("Dirichlet boundary conditions will not be applied ");
  }

  _initDofs();
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

  /* Initializing all nodal variables to zero*/
  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    m_prev_acc[node] = Real3::zero() ;
    m_prev_vel[node] = Real3::zero() ;
    m_prev_displ[node] = Real3::zero();
    m_prev_acc_iter[node] = Real3::zero() ;
    m_prev_vel_iter[node] = Real3::zero() ;
    m_prev_displ_iter[node] = Real3::zero() ;
    m_acc[node] = Real3::zero() ;
    m_vel[node] = Real3::zero() ;
    m_displ[node] = Real3::zero();
  }

//  _startInitGauss();
  _applyInitialNodeConditions();
  _initCells();
  _initBoundaryConditions();
  _initDCConditions();
  _startInitGauss();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_initDofs(){
  m_dofs_on_nodes.initialize(mesh(),NDIM);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_startInitGauss()
{
  Integer max_gauss_per_cell{ 0 };
  Integer max_nbnodes_per_cell{ 0 };
  Real lambda, mu;
  Integer max_nb_law_param;

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto cell_type = cell.type();
    auto nbgauss = ArcaneFemFunctions::FemGaussQuadrature::getNbGaussPointsfromOrder(cell_type, ninteg);
    m_nb_gauss[cell] = nbgauss;
    max_gauss_per_cell = math::max(nbgauss, max_gauss_per_cell);
    auto nbnodes = cell.nbNode();
    max_nbnodes_per_cell = math::max(nbnodes, max_nbnodes_per_cell);
  }
  // Make sure all sub-domains have the same number of maximum values
  IParallelMng* pm = defaultMesh()->parallelMng();
  max_gauss_per_cell = pm->reduce(Parallel::ReduceMax, max_gauss_per_cell);
  max_nbnodes_per_cell = pm->reduce(Parallel::ReduceMax, max_nbnodes_per_cell);

  m_gauss_on_cells.initialize(mesh(), max_gauss_per_cell);

  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());
  VariableDoFReal& gauss_weight(m_gauss_on_cells.gaussWeight());
  VariableDoFReal3& gauss_refpos(m_gauss_on_cells.gaussRefPosition());
  VariableDoFArrayReal& gauss_shape(m_gauss_on_cells.gaussShape());
  VariableDoFArrayReal3& gauss_shapederiv(m_gauss_on_cells.gaussShapeDeriv());
  VariableDoFArrayReal& gauss_lawparam(m_gauss_on_cells.gaussLawParam());
  VariableDoFArrayReal& gauss_histparam(m_gauss_on_cells.gaussLawHistoryParam());
  VariableDoFArrayReal3x3& gauss_stress(m_gauss_on_cells.gaussStress());
  VariableDoFArrayReal3x3& gauss_strain(m_gauss_on_cells.gaussStrain());
  VariableDoFArrayReal3x3& gauss_strain_plastic(m_gauss_on_cells.gaussStrainPlastic());
  VariableDoFArrayReal3x3& gauss_tangent_operator(m_gauss_on_cells.gaussTangentOperator());

  gauss_shape.resize(max_nbnodes_per_cell);
  gauss_shapederiv.resize(max_nbnodes_per_cell);
  gauss_lawparam.resize(m_nb_law_param);
  gauss_histparam.resize(m_nb_law_hist_param);

  /* gauss tensors (stress, strains) during the global computing time loop:
   * 0: values at start time (sig0, eps0, epsp0)
   * 1: values at previous step (sign, epsn, epspn)
   * 2: values at current step (sig, eps, epsp)
   * */
  gauss_stress.resize(3);
  gauss_strain.resize(3);
  gauss_strain_plastic.resize(3);

  /* Tangent operator at Gauss points is useful for nonlinear simulations
   * 0: 1st diagonal Real3x3 block (D)
   * 1: 2nd diagonal Real3x3 block (S)
   * 2: Upper out-of-diagonal Real3x3 block (Sup)
   * 3: Lower out-of-diagonal Real3x3 block (Slow) => in case of symmetry, Slow = 0
   * */
  gauss_tangent_operator.resize(4);

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto cell_type = cell.type();
    auto cell_nbnod = cell.nbNode();
    Int32 numcell = cell.localId();
    auto cell_nbgauss = ArcaneFemFunctions::FemGaussQuadrature::getNbGaussPointsfromOrder(cell_type, ninteg);
    Int32 ndim = ArcaneFemFunctions::MeshOperation::getGeomDimension(cell);

    lambda = m_lambda[cell];
    mu = m_mu[cell];

    auto is_default = m_default_law[cell];
    auto lawtyp = static_cast<TypesNLDynamic::eLawType>(m_law[cell]);
    auto ilaw = m_iparam_law[cell];
    LawDispatcher cell_law(lawtyp,is_default);
    RealUniqueArray lawparams = cell_law.readLawParams(lambda, mu, is_default, m_law_param_file,ilaw);

    auto nblaw = lawparams.size();// For debug only
    auto nbhist = cell_law.getNbLawHistoryParam();

    for (Int32 ig = 0; ig < cell_nbgauss; ++ig) {
      DoFLocalId gauss_pti = gauss_point.dofId(cell, ig);
      Int32 gaussnum = gauss_pti.localId();
      gauss_weight[gauss_pti] = ArcaneFemFunctions::FemGaussQuadrature::getGaussWeight(cell, ninteg, ig);
      Real3 refpos = ArcaneFemFunctions::FemGaussQuadrature::getGaussRefPosition(cell, ninteg, ig);
      if (ndim <= 2) {
        refpos.z = 0.;
        if (ndim == 1)
          refpos.y = 0.;
      }
      gauss_refpos[gauss_pti] = refpos;

      for (Int32 inod = 0; inod < cell_nbnod; ++inod) {
        auto coord_nod = m_node_coord[cell.node(inod)];
        auto Phi_i = cell_fem.getShapeFuncVal(cell_type, inod, refpos);
        gauss_shape[gauss_pti][inod] = Phi_i;

        auto dPhi = cell_fem.getShapeFuncDeriv(cell_type, inod, refpos);
        if (ndim <= 2) {
          dPhi.z = 0.;
          if (ndim == 1)
            dPhi.y = 0.;
        }
        gauss_shapederiv[gauss_pti][inod] = dPhi;
      }

      // Setting the law parameters on Gauss points for use during iterations
      gauss_lawparam[gauss_pti] = lawparams;

      // History parameters are set to 0. at this stage
      // (we will see later for more complex laws)
      gauss_histparam[gauss_pti] = RealUniqueArray(nbhist,0.);

      // Initializing all strains/stresses to zero
      gauss_strain[gauss_pti] = Real3x3UniqueArray(3);
      gauss_strain_plastic[gauss_pti] = Real3x3UniqueArray(3);
      gauss_stress[gauss_pti] = Real3x3UniqueArray(3);
    }
  }

  auto hasStressFile{ options()->hasStress0PerCellFile() };

/*  VariableDoFArrayTensor2& gauss_stress(m_gauss_on_cells.gaussStress());
  VariableDoFArrayTensor2& gauss_strain(m_gauss_on_cells.gaussStrain());
  VariableDoFArrayTensor2& gauss_strain_plastic(m_gauss_on_cells.gaussStrainPlastic());
*/

  Real3x3 sig0{};
  CaseTable* stress_table{ nullptr};

  if (hasStressFile) {
    String file_stress = options()->stress0PerCellFile();
    if (!file_stress.empty()) {
      stress_table = readFileAsCaseTable(pm, file_stress, 3);
      if (stress_table == nullptr)
        hasStressFile = false;
      else
        info() << "Initializing stresses per cell from user file (CaseTable) " << file_stress;
    }
  }

    if (hasStressFile) {
      Real3 sig3{};
      ENUMERATE_CELL (icell, allCells()) {
        const Cell& cell = *icell;
        auto numcell = cell.uniqueId().asInt32();
        auto cell_type = cell.type();
        auto cell_nbgauss = ArcaneFemFunctions::FemGaussQuadrature::getNbGaussPointsfromOrder(cell_type, ninteg);
        stress_table->value(numcell, sig3);
        Tensor2 t0 = { sig3.x, sig3.y, sig3.z, 0., 0., 0. };
        sig0 = fromTensor2Real3x3(t0);

        for (Int32 ig = 0; ig < cell_nbgauss; ++ig) {
          DoFLocalId gauss_pti = gauss_point.dofId(cell, ig);
          Int32 gaussnum = gauss_pti.localId();

          for (Int32 is = 0; is < 2; is++) { //sig0, then sign
            gauss_stress[gauss_pti][is] = sig0;
          }
          gauss_stress[gauss_pti][2] = Real3x3::zero(); //sig
        }
      }
    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_initCells(){
  Real vp, vs, E, nu, lambda, mu;

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto rho = m_rho[cell];
    if (elast_type == TypesNLDynamic::YoungNu) {
      E = m_young[cell];
      nu = m_nu[cell];
      lambda = nu*E/(1. + nu)/(1. - 2.*nu);
      mu = E/2./(1. + nu);
      vp = math::sqrt( (lambda + 2. * mu)/rho );
      vs = math::sqrt( mu/rho );

    } else if (elast_type == TypesNLDynamic::Lame) {
      lambda = m_lambda[cell];
      mu = m_mu[cell];
      vp = math::sqrt( (lambda + 2.*mu)/rho );
      vs = math::sqrt( mu/rho );
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);

    } else if (elast_type == TypesNLDynamic::Veloc) {
      vp = m_vp[cell];
      vs = m_vs[cell];
      mu = rho*vs*vs;
      lambda = rho*vp*vp - 2.*mu;
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);

    }
    m_vp[cell] = vp;
    m_vs[cell] = vs;
    m_lambda[cell] = lambda;
    m_mu[cell] = mu;
    m_young[cell] = E;
    m_nu[cell] = nu;
  }

  _applyInitialCellConditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_applyInitialNodeConditions(){

  for (Int32 i = 0, nb = options()->initialNodeCondition().size(); i < nb; ++i) {

    NodeGroup node_group = options()->initialNodeCondition[i]->nodeGroup();
    info() << "Applying init-node-condition for node group " << node_group.name();

    // Loop on nodes with this initial condition
    ENUMERATE_NODE(inode, node_group) {
      const Node & node = *inode;

      if (options()->initialNodeCondition[i]->hasA())
        m_prev_acc[node] = options()->initialNodeCondition[i]->A();

      if (options()->initialNodeCondition[i]->hasV())
        m_prev_vel[node] = options()->initialNodeCondition[i]->V();

      if (options()->initialNodeCondition[i]->hasU())
        m_prev_displ[node] = options()->initialNodeCondition[i]->U();

      if (options()->initialNodeCondition[i]->hasF())
        m_force[node] = options()->initialNodeCondition[i]->F();
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_applyInitialCellConditions(){

  m_nb_law_param = 2;
  m_nb_law_hist_param = 0.;
  if (options()->hasLawInputParam())
    m_law_param_file = options()->lawInputParam();

  for (Integer i = 0, nb = options()->lawModel().size(); i < nb; ++i){

    CellGroup cell_group = options()->lawModel[i]->cellGroup();
    auto ilaw = options()->lawModel[i]->iLawParam();
    auto lawtyp = options()->lawModel[i]->lawType();
    auto nblaw = options()->lawModel[i]->nbLawParam();
    auto nbhist = options()->lawModel[i]->nbLawHistParam();
    m_nb_law_param = math::max(m_nb_law_param, nblaw);
    m_nb_law_hist_param = math::max(m_nb_law_hist_param, nbhist);

    bool is_default = (m_law_param_file.empty() || nblaw == 2);

    ENUMERATE_CELL (icell, cell_group) {
      const Cell& cell = *icell;
      m_law[cell] = (int)lawtyp;
      m_default_law[cell] = is_default;
      m_iparam_law[cell] = ilaw;
    }
  }

  for (Integer i = 0, nb = options()->initElastProperties().size(); i < nb; ++i) {

    CellGroup cell_group = options()->initElastProperties[i]->cellGroup();

    // In the future, we will have to find a way to impose different initial
    // properties (stress/strain tensors, densities...) per element from a file
    // (e.g., coming from a previous computation)
    auto rho = options()->initElastProperties[i]->rho();
    Real vp, vs, E, nu, lambda, mu;

    if (elast_type == TypesNLDynamic::YoungNu) {
      E = options()->initElastProperties[i]->young();
      nu = options()->initElastProperties[i]->nu();
      lambda = nu*E/(1. + nu)/(1. - 2.*nu);
      mu = E/2./(1. + nu);
      vp = math::sqrt( (lambda + 2. * mu)/rho );
      vs = math::sqrt( mu/rho );

    } else if (elast_type == TypesNLDynamic::Lame) {
      lambda = options()->initElastProperties[i]->young();
      mu = options()->initElastProperties[i]->nu();
      vp = math::sqrt( (lambda + 2.*mu)/rho );
      vs = math::sqrt( mu/rho );
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);

    } else if (elast_type == TypesNLDynamic::Veloc) {
      vp = options()->initElastProperties[i]->vp();
      vs = options()->initElastProperties[i]->vs();
      mu = rho*vs*vs;
      lambda = rho*vp*vp - 2.*mu;
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);

    }

    ENUMERATE_CELL (icell, cell_group) {
      const Cell& cell = *icell;
      m_rho[cell] = rho;
      m_vp[cell] = vp;
      m_vs[cell] = vs;
      m_lambda[cell] = lambda;
      m_mu[cell] = mu;
      m_young[cell] = E;
      m_nu[cell] = nu;
    }
  }

  /*  for (Integer i = 0, nb = options()->initCellCondition().size(); i < nb; ++i) {

    CellGroup cell_group = options()->initCellCondition[i]->cellGroup();
    info() << "Applying init-cell-condition for cell group " << cell_group.name();

    // In the future, we will have to find a way to impose different initial
    // properties (stress/strain tensors, densities...) per element from a file
    // (e.g., coming from a previous computation)

    // Loop on cells with this initial condition
    ENUMERATE_CELL (icell, cell_group) {
      const Cell& cell = *icell;

      // Initialize the stress tensor for the concerned cell
      if (options()->initCellCondition[i]->hasDevStrain())
        m_strain_dev[cell] = options()->initCellCondition[i]->devStrain();
      if (options()->initCellCondition[i]->hasVolStrain())
        m_strain_vol[cell] = options()->initCellCondition[i]->volStrain();
      if (options()->initCellCondition[i]->hasDevStress())
        m_stress_dev[cell] = options()->initCellCondition[i]->devStress();
      if (options()->initCellCondition[i]->hasVolStrain())
        m_stress_vol[cell] = options()->initCellCondition[i]->volStress();

    }
  }*/
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_initGaussStep()
{
  Real eps{1.0e-15};
  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());

  VariableDoFReal& gauss_jacobian(m_gauss_on_cells.gaussJacobian());
  VariableDoFReal3& gauss_refpos(m_gauss_on_cells.gaussRefPosition());
  VariableDoFReal3x3& gauss_jacobmat(m_gauss_on_cells.gaussJacobMat());
  VariableDoFArrayReal& gauss_shape(m_gauss_on_cells.gaussShape());
  VariableDoFArrayReal3& gauss_shapederiv(m_gauss_on_cells.gaussShapeDeriv());
  VariableDoFArrayReal3x3& gauss_stress(m_gauss_on_cells.gaussStress());
  VariableDoFArrayReal3x3& gauss_strain(m_gauss_on_cells.gaussStrain());
  VariableDoFArrayReal3x3& gauss_strain_plastic(m_gauss_on_cells.gaussStrainPlastic());

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto cell_type = cell.type();
    auto cell_nbnod = cell.nbNode();
    Int32 numcell = cell.localId();
    Int32 ndim = ArcaneFemFunctions::MeshOperation::getGeomDimension(cell);
    auto cell_nbgauss = ArcaneFemFunctions::FemGaussQuadrature::getNbGaussPointsfromOrder(cell_type, ninteg);

    for (Int32 ig = 0; ig < cell_nbgauss; ++ig) {
      DoFLocalId gauss_pti = gauss_point.dofId(cell, ig);
      Int32 gaussnum = gauss_pti.localId();
      Real3 refpos = gauss_refpos[gauss_pti];

      Real3x3	jac;
      Real jacobian;
      for (Int32 inod = 0; inod < cell_nbnod; ++inod) {

        auto dPhi = gauss_shapederiv[gauss_pti][inod];
        auto coord_nod = m_node_coord[cell.node(inod)];
        for (int i = 0; i < NDIM; ++i){
          for (int j = 0; j < NDIM; ++j){
            jac[i][j] += dPhi[i] * coord_nod[j];
          }
        }
      }

      if (ndim == 3)
        jacobian = math::matrixDeterminant(jac);

      else if (ndim == 2)
        jacobian = jac.x.x * jac.y.y - jac.x.y * jac.y.x;
      else
        jacobian = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(cell, m_node_coord) / 2.;

      if (fabs(jacobian) < eps) {
        ARCANE_FATAL("Cell jacobian is null");
      }
      gauss_jacobian[gauss_pti] = jacobian;
      gauss_jacobmat[gauss_pti] = jac;

      // Tensors are initialized to the value of previous step
      gauss_stress[gauss_pti][1] = gauss_stress[gauss_pti][2];
      gauss_strain[gauss_pti][1] = gauss_strain[gauss_pti][2];
      gauss_strain_plastic[gauss_pti][1] = gauss_strain_plastic[gauss_pti][2];
      gauss_stress[gauss_pti][2] = Real3x3::zero();
      gauss_strain[gauss_pti][2] = Real3x3::zero();
      gauss_strain_plastic[gauss_pti][2] = Real3x3::zero();
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_iterate(){

  info() << "Module PASSMO-NL ITERATION LOOP";

  auto dt = m_global_deltat();
  dt2 = dt * dt;

  m_converge = false;
  Int32 iter{0}; // counter for imbalance iterations

  // Starting the iteration loop until convergence is reached
  do {

    // Compute stresses from the nonlinear constitutive models
    // with  displacements predicted according to the integration scheme (e.g. Newmark)
    _stress_prediction(!iter);

    // Add the nonlinear contributions to the RHS (B)
    _assembleNonLinRHS();

    // to do: add contributions of the set dofs to the RHS (B)

    // Solve the linear system AX = B
    _doSolve();

    _check_convergence(iter);

    // Update stresses from the nonlinear constitutive models and computed solution (displacements)
    _stress_correction();

    // Update nodal quantities for this iteration (no convergence) or time step (if convergence reached)
    // (according to the integration scheme, e.g. Newmark)
    _updateNewmark();

    iter++;

  } while (iter < ite_max && !m_converge);

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
compute(){

  info() << "Module PASSMO-NL COMPUTE";
  ++linop_nstep_counter;

  // Stop code at exact final time set by user
  auto tf = m_global_final_time();
  auto t = m_global_time();
  auto dt = m_global_deltat();
  auto t0 = options()->getStart();
  dt2 = dt * dt;

  info() << "Time (s) = " << t;

  // Set if we want to keep the matrix structure between calls
  // the rate is a user input (linop_nstep)
  // The matrix has to have the same structure (same structure for non-zero)

  //  if (m_linear_system.isInitialized() && (linop_nstep_counter < linop_nstep || keep_constop || is_linear)){
  if ( m_linear_system.isInitialized() &&
      (is_linear || keep_constop || (algo_type == TypesNLDynamic::ModNewtonRaphson && linop_nstep_counter < linop_nstep)) ){
    m_linear_system.clearValues();
    m_ref = false;
  }
  else {

    m_ref = true;
    m_linear_system.reset();
    m_linear_system.setLinearSystemFactory(options()->linearSystem());
    m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

    // Reset the counter when the linear operator is reset
    linop_nstep_counter = 0;
  }

  // This part is in common for linear and nonlinear simulations

  { // Update Gauss shape functions and their derivatives for this step
    _initGaussStep();

    // Apply Dirichlet/Neumann conditions if any
    _applyDirichletBoundaryConditions();
    _applyNeumannBoundaryConditions();

    // Apply Paraxial conditions if any
    _applyParaxialBoundaryConditions();

    // Assemble the FEM global LHS operator (A matrix)
    _assembleLinearLHS();

    // Compute the predicted displacements and velocities
    // at beginning of step (from previous step values)
    _predictNewmark();
  }

  if (is_linear)
  {
    // Assemble the FEM global RHS operator (B matrix)
    _assembleLinearRHS();

    // Solve the linear system AX = B
    _doSolve();

    // Update the nodal variable according to the integration scheme (e.g. Newmark)
    _updateNewmark();
    m_converge = true;
  }
  else {
    // Adding tangent_operators contributions at Gauss points
    // to global LHS operator (A matrix)
    if (m_ref){
      // to do
    }

    // Starting the iteration loop in case of nonlinear problems
    _iterate();
  }
  if (t < tf && m_converge) {
    if (t + dt > tf) {
      dt = tf - t;
      m_global_deltat = dt;
    }
  }
  else {
    // Save/Check results
    _checkResultFile();
    subDomain()->timeLoopMng()->stopComputeLoop(true);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NLDynamicModule::
_checkResultFile()
{
  String filename = options()->resultFile();
  info() << "CheckResultFile filename=" << filename;
  if (filename.empty())
    return;
  const double epsilon = 1.0e-4;
  const double min_value_to_test = 1.0e-10;
  Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_displ, epsilon, min_value_to_test);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_predictNewmark(){

  // Predicting the nodal displacements and velocities (before solve)
  auto dt = m_global_deltat();
  auto betadt2 = 1/beta/dt2;
  auto gammadt = 1./gamma/dt;

  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    auto an = m_prev_acc[node];
    auto vn = m_prev_vel[node];
    auto dn = m_prev_displ[node];

    if (!is_alfa_method) {
      for (Int32 i = 0; i < NDIM; ++i) {

        auto bu = (bool)m_imposed_displ[node][i];
        auto bv = (bool)m_imposed_vel[node][i];
        auto ba = (bool)m_imposed_acc[node][i];
        auto vi = vn[i] + dt * (1. - gamma) * an[i];
        auto di = dn[i] + dt * vn[i] + (0.5 - beta) * dt2 * an[i];

        if (ba){
          auto ai = m_acc[node][i];
          m_displ[node][i] = di + beta * dt2 * ai;
          m_vel[node][i] = vi + gamma * dt * ai;
        }
        else {
          if (!bu)
            m_displ[node][i] = di;

          else {
            auto ai = betadt2 * (m_displ[node][i] - di);
            m_acc[node][i] = ai;
            m_vel[node][i] = vi + gamma * dt * ai;
          }

          if (!bv)
            m_vel[node][i] = vi;

          else {
            auto ai = gammadt * (m_vel[node][i] - vi);
            m_acc[node][i] = ai;
            m_displ[node][i] = di + beta * dt2 * ai;
          }
        }
      }
    } else {
      // TO DO
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_updateNewmark(){

  // Updating the nodal accelerations and velocities (after solve)
  auto dt = m_global_deltat();

  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    auto an = m_prev_acc[node];
    auto vn = m_prev_vel[node];
    auto dn = m_prev_displ[node];
    auto dn1 = m_displ[node]; // displacements for the current iteration

    if (!is_alfa_method) {
      for (Int32 i = 0; i < NDIM; ++i) {

        auto ba = (bool)m_imposed_acc[node][i];
        auto bv = (bool)m_imposed_vel[node][i];
        auto ui = dn[i] + dt * vn[i] + dt2 * (0.5 - beta) * an[i];
        auto vi = vn[i] + dt * (1. - gamma) * an[i];

        if (!ba)
          m_acc[node][i] = (dn1[i] - ui)/beta/dt2;
        else
          m_displ[node][i] = ui + beta * dt2 * m_acc[node][i];

        if (!bv)
          m_vel[node][i] = vi + dt*gamma*m_acc[node][i];
      }
    } else {
      // TO DO
    }

    m_prev_acc_iter[node] = m_acc[node];
    m_prev_vel_iter[node] = m_vel[node];
    m_prev_displ_iter[node] = m_displ[node];

    if (m_converge) {
      m_prev_acc[node] = m_acc[node];
      m_prev_vel[node] = m_vel[node];
      m_prev_displ[node] = m_displ[node];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_initBoundaryConditions()
{
  IParallelMng* pm = subDomain()->parallelMng();

  for (const auto& bd : options()->dirichletSurfaceCondition()) {
    FaceGroup face_group = bd->surface();
    info() << "Initializing Dirichlet Surface Boundary Conditions for face group " << face_group.name();

    if (bd->hasACurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_sacc_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasUCurve()) {
      String file_name = bd->UCurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_sdispl_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasVCurve()) {
      String file_name = bd->VCurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_svel_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasFCurve()) {
      String file_name = bd->FCurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_sforce_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    auto hasUcurve{bd->hasUCurve()};
    auto hasVcurve{bd->hasVCurve()};
    auto hasAcurve{bd->hasACurve()};
    auto hasFcurve{bd->hasFCurve()};
    auto xdir{bd->getXAxis()};
    auto ydir{bd->getYAxis()};
    auto zdir{bd->getZAxis()};

    // Loop on faces of the surface
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;
      Int32 nb_node = face.nbNode();

      // Loop on nodes of the face
      for (Int32 k = 0; k < nb_node; ++k) {
        const Node& node = face.node(k);
        auto coord = m_node_coord[node];
        auto num = node.uniqueId();

        m_imposed_displ[node].x = (bd->hasUx() || (hasUcurve && xdir) ? 1 : 0);
        m_imposed_displ[node].y = (bd->hasUy() || (hasUcurve && ydir) ? 1 : 0);
        m_imposed_displ[node].z = (bd->hasUz() || (hasUcurve && zdir) ? 1 : 0);

        m_imposed_acc[node].x = (bd->hasAx() || (hasAcurve && xdir) ? 1 : 0);
        m_imposed_acc[node].y = (bd->hasAy() || (hasAcurve && ydir) ? 1 : 0);
        m_imposed_acc[node].z = (bd->hasAz() || (hasAcurve && zdir) ? 1 : 0);

        m_imposed_vel[node].x = (bd->hasVx() || (hasVcurve && xdir) ? 1 : 0);
        m_imposed_vel[node].y = (bd->hasVy() || (hasVcurve && ydir) ? 1 : 0);
        m_imposed_vel[node].z = (bd->hasVz() || (hasVcurve && zdir) ? 1 : 0);

        m_imposed_force[node].x = (bd->hasFx() || (hasFcurve && xdir) ? 1 : 0);
        m_imposed_force[node].y = (bd->hasFy() || (hasFcurve && ydir) ? 1 : 0);
        m_imposed_force[node].z = (bd->hasFz() || (hasFcurve && zdir) ? 1 : 0);
      }
    }
  }

  for (const auto& bd : options()->dirichletPointCondition()) {
    NodeGroup nodes = bd->node();
    info() << "Initializing Dirichlet Point Conditions for node group " << nodes.name();

    if (bd->hasACurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_acc_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasUCurve()) {
      String file_name = bd->UCurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_displ_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasVCurve()) {
      String file_name = bd->VCurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_vel_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasFCurve()) {
      String file_name = bd->FCurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_force_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    auto hasUcurve{bd->hasUCurve()};
    auto hasVcurve{bd->hasVCurve()};
    auto hasAcurve{bd->hasACurve()};
    auto hasFcurve{bd->hasFCurve()};
    auto xdir{bd->getXAxis()};
    auto ydir{bd->getYAxis()};
    auto zdir{bd->getZAxis()};

    // Loop on nodes
    ENUMERATE_NODE (inode, nodes) {
      const Node& node = *inode;
      auto coord = m_node_coord[node];
      auto num = node.uniqueId();

      m_imposed_acc[node].x = (bd->hasAx() || (hasAcurve && xdir) ? 1 : 0);
      m_imposed_acc[node].y = (bd->hasAy() || (hasAcurve && ydir) ? 1 : 0);
      m_imposed_acc[node].z = (bd->hasAz() || (hasAcurve && zdir) ? 1 : 0);

      m_imposed_vel[node].x = (bd->hasVx() || (hasVcurve && xdir) ? 1 : 0);
      m_imposed_vel[node].y = (bd->hasVy() || (hasVcurve && ydir) ? 1 : 0);
      m_imposed_vel[node].z = (bd->hasVz() || (hasVcurve && zdir) ? 1 : 0);

      m_imposed_force[node].x = (bd->hasFx() || (hasFcurve && xdir) ? 1 : 0);
      m_imposed_force[node].y = (bd->hasFy() || (hasFcurve && ydir) ? 1 : 0);
      m_imposed_force[node].z = (bd->hasFz() || (hasFcurve && zdir) ? 1 : 0);

      if ((bool)m_imposed_acc[node].x || (bool)m_imposed_vel[node].x
          || bd->hasUx() || (hasUcurve && xdir))
        m_imposed_displ[node].x = 1;

      if ((bool)m_imposed_acc[node].y || (bool)m_imposed_vel[node].y
          || bd->hasUy() || (hasUcurve && ydir))
        m_imposed_displ[node].y = 1;

      if ((bool)m_imposed_acc[node].z || (bool)m_imposed_vel[node].z
          || bd->hasUz() || (hasUcurve && zdir))
        m_imposed_displ[node].z = 1;
    }
  }

  for (const auto& bs : options()->neumannCondition()) {
    FaceGroup face_group = bs->surface();
    info() << "Initializing Neumann (traction) Conditions for face group " << face_group.name();
    String file_name = bs->getCurve();
    if (!file_name.empty()) {
      auto case_table = readFileAsCaseTable(pm, file_name, 3);
      m_traction_case_table_list.add(CaseTableInfo{ file_name, case_table });
    }
  }

  for (const auto& bs : options()->paraxialBoundaryCondition()) {

    FaceGroup face_group = bs->surface();
    info() << "Initializing Paraxial Boundary Conditions for face group " << face_group.name();

    if (bs->getInputMotionType() == 6) {
      if (bs->hasUInput()) {
        String file_name = bs->getUInput();
        if (!file_name.empty()) {
          auto case_table = readFileAsCaseTable(pm, file_name, 3);
          m_uin_case_table_list.add(CaseTableInfo{ file_name, case_table });
        }
      }
      if (bs->hasVInput()) {
        String file_name = bs->getVInput();
        if (!file_name.empty()) {
          auto case_table = readFileAsCaseTable(pm, file_name, 3);
          m_vin_case_table_list.add(CaseTableInfo{ file_name, case_table });
        }
      }
      if (bs->hasAInput()) {
        String file_name = bs->getAInput();
        if (!file_name.empty()) {
          auto case_table = readFileAsCaseTable(pm, file_name, 3);
          m_ain_case_table_list.add(CaseTableInfo{ file_name, case_table });
        }
      }
    }

    auto rho = bs->getRhopar();
    Real cs, cp;
    bool is_inner{ false };

    if (bs->hasEPar() && bs->hasNuPar()) {

      auto E = bs->getEPar();
      auto nu = bs->getNuPar();
      auto lambda = nu * E / (1. + nu) / (1. - 2. * nu);
      auto mu = E / 2. / (1. + nu);
      cp = math::sqrt((lambda + 2. * mu) / rho);
      cs = math::sqrt(mu / rho);
    }
    else if (bs->hasCp() && bs->hasCs()) {

      cp = bs->getCp();
      cs = bs->getCp();
    }
    else if (bs->hasLambdaPar() && bs->hasMuPar()) {

      auto mu = bs->getMuPar();
      cp = math::sqrt((bs->getLambdaPar() + 2. * mu) / rho);
      cs = math::sqrt(mu / rho);
    }
    else {
      info() << "Elastic properties expected for "
             << "Paraxial boundary condition on FaceGroup "
             << face_group.name() << ": \n"
             << "  - (E-par, nu-par) or\n"
             << "  - (lambda-par, mu-par) or\n"
             << "  - (cp, cs)\n";
      info() << "When not specified, taking elastic properties from inner domain. ";
      is_inner = true;
    }

    // Loop on the faces (=edges in 2D) concerned with the paraxial condition
    // Initializing the local referential per face (just done once) for further use
    ENUMERATE_FACE (iface, face_group) {

      const Face& face = *iface;

      if (face.isSubDomainBoundary() && face.isOwn()) {

        Real3 e1{ 0. }, e2{ 0. }, e3{ 0. };
        ArcaneFemFunctions::MeshOperation::dirVectors(face, m_node_coord, NDIM, e1, e2, e3);
        m_e1_boundary[face] = e1;
        m_e2_boundary[face] = e2;
        m_e3_boundary[face] = e3;

        if (is_inner) {
          const Cell& cell = face.boundaryCell();
          rho = m_rho[cell];
          cs = m_vs[cell];
          cp = m_vp[cell];
        }

        m_rho_parax[face] = rho;
        m_vel_parax[face].x = cs;

        if (NDIM == 3){
          m_vel_parax[face].y = cs;
          m_vel_parax[face].z = cp;
        }
        else{
          m_vel_parax[face].y = cp;
          m_vel_parax[face].z = 0.;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_initDCConditions()
{
  IParallelMng* pm = subDomain()->parallelMng();

  for (const auto& bd : options()->doubleCouple()) {
    NodeGroup east = bd->getEastNode();
    NodeGroup west = bd->getWestNode();
    NodeGroup north = bd->getNorthNode();
    NodeGroup south = bd->getSouthNode();

    info() << "Initializing Double-Couple Conditions for nodes:\n"
           << "  - North = " << north.name() << "\n"
           << "  - South = " << south.name() << "\n"
           << "  - East = " << east.name() << "\n"
           << "  - West = " << west.name() << "\n";

    auto hasMoment = bd->hasSeismicMomentFile();
    auto hasLoading = bd->hasLoadingFile();
    String file_name;

    if (hasMoment || hasLoading) {
       if (hasMoment) file_name = bd->getSeismicMomentFile();
       else file_name = bd->getLoadingFile();

      auto case_table = readFileAsCaseTable(pm, file_name, 1);
      m_dc_case_table_list.add(CaseTableInfo{ file_name, case_table });
    }
    else {
      info() << "Double-Couple Error: seismic moment or user loading not provided! ";
      ARCANE_FATAL("Double-Couple conditions cannot be applied");
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_applyDirichletBoundaryConditions(){

  Real time = globalTime();
  Int32 sac_index{ 0 }, svc_index{ 0 }, suc_index{ 0 }, sfc_index{ 0 };
  for (const auto& bd : options()->dirichletSurfaceCondition()) {
    FaceGroup face_group = bd->surface();

    Real3 acc{};
    bool is_acc_imp{bd->hasACurve() || bd->hasAx() || bd->hasAy() || bd->hasAz()};
    if (bd->hasACurve()) {
      const CaseTableInfo& table_info = m_sacc_case_table_list[sac_index++];
      String file_name = bd->ACurve();
      info() << "Applying acceleration boundary conditions for surface " << face_group.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, acc);
    }
    else if (is_acc_imp) {
      if (bd->hasAx())
        acc.x = bd->getAx();
      if (bd->hasVy())
        acc.y = bd->getAy();
      if (bd->hasAz())
        acc.z = bd->getAz();
    }

    Real3 vel{};
    bool is_vel_imp{bd->hasVCurve() || bd->hasVx() || bd->hasVy() || bd->hasVz()};
    if (bd->hasVCurve()) {
      const CaseTableInfo& table_info = m_svel_case_table_list[svc_index++];
      String file_name = bd->VCurve();
      info() << "Applying velocity boundary conditions for surface " << face_group.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, vel);
    }
    else if (is_vel_imp){
      if (bd->hasVx())
        vel.x = bd->getVx();
      if (bd->hasVy())
        vel.y = bd->getVy();
      if (bd->hasVz())
        vel.z = bd->getVz();
    }

    Real3 displ{};
    bool is_displ_imp{bd->hasUCurve() || bd->hasUx() || bd->hasUy() || bd->hasUz()};
    if (bd->hasUCurve()) {
      const CaseTableInfo& table_info = m_sdispl_case_table_list[suc_index++];
      String file_name = bd->UCurve();
      info() << "Applying displacement boundary conditions for surface " << face_group.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, displ);
    }
    else if (is_displ_imp){
      if (bd->hasUx())
        displ.x = bd->getUx();
      if (bd->hasUy())
        displ.y = bd->getUy();
      if (bd->hasUz())
        displ.z = bd->getUz();
    }

    Real3 force{};
    bool is_force_imp{bd->hasFCurve() || bd->hasFx() || bd->hasFy() || bd->hasFz()};
    if (bd->hasFCurve()) {
      const CaseTableInfo& table_info = m_sforce_case_table_list[sfc_index++];
      String file_name = bd->FCurve();
      info() << "Applying force boundary conditions for surface " << face_group.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, force);
    }
    else if (is_force_imp){
      if (bd->hasFx())
        force.x = bd->getFx();
      if (bd->hasFy())
        force.y = bd->getFy();
      if (bd->hasFz())
        force.z = bd->getFz();
    }

    // Loop on faces of the surface
    ENUMERATE_FACE (iface, face_group) {

      // Loop on nodes of the face
      for (Node node : iface->nodes()) {

        //--- For debug only!!!
        auto coord = m_node_coord[node];
        auto num = node.uniqueId();

        if (is_acc_imp) {
          if ((bool)m_imposed_acc[node].x)
            m_acc[node].x = acc.x;

          if ((bool)m_imposed_acc[node].y)
            m_acc[node].y = acc.y;

          if ((bool)m_imposed_acc[node].z)
            m_acc[node].z = acc.z;
        }

        if (is_vel_imp) {
          if ((bool)m_imposed_vel[node].x)
            m_vel[node].x = vel.x;

          if ((bool)m_imposed_vel[node].y)
            m_vel[node].y = vel.y;

          if ((bool)m_imposed_vel[node].z)
            m_vel[node].z = vel.z;
        }

        if (is_displ_imp) {
          if ((bool)m_imposed_displ[node].x)
            m_displ[node].x = displ.x;

          if ((bool)m_imposed_displ[node].y)
            m_displ[node].y = displ.y;

          if ((bool)m_imposed_displ[node].z)
            m_displ[node].z = displ.z;
        }

        if (is_force_imp) {
          if ((bool)m_imposed_force[node].x)
            m_force[node].x = force.x;

          if ((bool)m_imposed_force[node].y)
            m_force[node].y = force.y;

          if ((bool)m_imposed_force[node].z)
            m_force[node].z = force.z;
        }
      }
    }
  }

  Int32 ac_index{ 0 }, vc_index{ 0 }, uc_index{ 0 }, fc_index{ 0 };
  for (const auto& bd : options()->dirichletPointCondition()) {
    NodeGroup nodes = bd->node();

    Real3 acc{};
    bool is_acc_imp{bd->hasACurve() || bd->hasAx() || bd->hasAy() || bd->hasAz()};
    if (bd->hasACurve()) {
      const CaseTableInfo& table_info = m_acc_case_table_list[ac_index++];
      String file_name = bd->ACurve();
      info() << "Applying acceleration boundary conditions for nodes " << nodes.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, acc);
    }
    else if (is_acc_imp) {
      if (bd->hasAx())
        acc.x = bd->getAx();
      if (bd->hasVy())
        acc.y = bd->getAy();
      if (bd->hasAz())
        acc.z = bd->getAz();
    }

    Real3 vel{};
    bool is_vel_imp{bd->hasVCurve() || bd->hasVx() || bd->hasVy() || bd->hasVz()};
    if (bd->hasVCurve()) {
      const CaseTableInfo& table_info = m_vel_case_table_list[vc_index++];
      String file_name = bd->VCurve();
      info() << "Applying velocity boundary conditions for nodes " << nodes.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, vel);
    }
    else if (is_vel_imp){
      if (bd->hasVx())
        vel.x = bd->getVx();
      if (bd->hasVy())
        vel.y = bd->getVy();
      if (bd->hasVz())
        vel.z = bd->getVz();
    }

    Real3 displ{};
    bool is_displ_imp{bd->hasUCurve() || bd->hasUx() || bd->hasUy() || bd->hasUz()};
    if (bd->hasUCurve()) {
      const CaseTableInfo& table_info = m_displ_case_table_list[uc_index++];
      String file_name = bd->UCurve();
      info() << "Applying displacement boundary conditions for nodes " << nodes.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, displ);
    }
    else if (is_displ_imp){
      if (bd->hasUx())
        displ.x = bd->getUx();
      if (bd->hasUy())
        displ.y = bd->getUy();
      if (bd->hasUz())
        displ.z = bd->getUz();
    }

    Real3 force{};
    bool is_force_imp{bd->hasFCurve() || bd->hasFx() || bd->hasFy() || bd->hasFz()};
    if (bd->hasFCurve()) {
      const CaseTableInfo& table_info = m_force_case_table_list[fc_index++];
      String file_name = bd->FCurve();
      info() << "Applying force boundary conditions for nodes " << nodes.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, force);
    }
    else if (is_force_imp){
      if (bd->hasFx())
        force.x = bd->getFx();
      if (bd->hasFy())
        force.y = bd->getFy();
      if (bd->hasFz())
        force.z = bd->getFz();
    }

    // Loop on nodes
    ENUMERATE_NODE (inode, nodes) {
      const Node& node = *inode;

      //--- For debug only!!!
      auto coord = m_node_coord[node];
      auto num = node.uniqueId();

      if (is_acc_imp) {
        if ((bool)m_imposed_acc[node].x)
          m_acc[node].x = acc.x;

        if ((bool)m_imposed_acc[node].y)
          m_acc[node].y = acc.y;

        if ((bool)m_imposed_acc[node].z)
          m_acc[node].z = acc.z;
      }

      if (is_vel_imp) {
        if ((bool)m_imposed_vel[node].x)
          m_vel[node].x = vel.x;

        if ((bool)m_imposed_vel[node].y)
          m_vel[node].y = vel.y;

        if ((bool)m_imposed_vel[node].z)
          m_vel[node].z = vel.z;
      }

      if (is_displ_imp) {
        if ((bool)m_imposed_displ[node].x)
          m_displ[node].x = displ.x;

        if ((bool)m_imposed_displ[node].y)
          m_displ[node].y = displ.y;

        if ((bool)m_imposed_displ[node].z)
          m_displ[node].z = displ.z;
      }

      if (is_force_imp) {
        if ((bool)m_imposed_force[node].x)
          m_force[node].x = force.x;

        if ((bool)m_imposed_force[node].y)
          m_force[node].y = force.y;

        if ((bool)m_imposed_force[node].z)
          m_force[node].z = force.z;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_applyNeumannBoundaryConditions(){

  Real time = globalTime();
  Int32 bc_index{ 0 };
  for (const auto& bs : options()->neumannCondition()) {
    FaceGroup face_group = bs->surface();

    Real3 trac{};

    if (bs->hasCurve()) {

      const CaseTableInfo& case_table_info = m_traction_case_table_list[bc_index++];
      String file_name = bs->getCurve();
      info() << "Applying traction boundary conditions for surface " << face_group.name()
             << " via CaseTable" << file_name;
      CaseTable* inn = case_table_info.case_table;

      if (inn != nullptr)
        inn->value(time, trac);
    }
    else {
      if (bs->hasXVal())
        trac.x = bs->getXVal();
      if (bs->hasYVal())
        trac.y = bs->getYVal();
      if (bs->hasZVal())
        trac.z = bs->getZVal();
    }

    // Loop on faces of the surface
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;
      m_imposed_traction[face] = trac;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_applyParaxialBoundaryConditions(){

  Real time = globalTime();

  Int32 ui_index{ 0 }, vi_index{ 0 }, ai_index{ 0 };
  for (const auto& bs : options()->paraxialBoundaryCondition()) {

    FaceGroup face_group = bs->surface();

    // Looking for an input motion defined on the paraxial boundary
    // default = none (typ = 0)
    Int32 typ{bs->getInputMotionType()};

    Real3 uin{}, vin{}, ain{};
    if (typ > 0){

      bool is_u = bs->hasUInput();
      bool is_v = bs->hasVInput();
      bool is_a = bs->hasAInput();

      if (typ == 6) {
        if (is_u) {
            const CaseTableInfo& table_info = m_uin_case_table_list[ui_index++];
            String file_name = bs->getUInput();
            info() << "Applying input displacement for paraxial element " << face_group.name()
                   << " via CaseTable " << file_name;
            CaseTable* inn = table_info.case_table;

            if (inn != nullptr)
            inn->value(time, uin);
        }
        if (is_v) {
            const CaseTableInfo& table_info = m_uin_case_table_list[vi_index++];
            String file_name = bs->getVInput();
            info() << "Applying input velocity for paraxial element " << face_group.name()
                   << " via CaseTable " << file_name;
            CaseTable* inn = table_info.case_table;

            if (inn != nullptr)
            inn->value(time, vin);
        }
        if (is_a) {
            const CaseTableInfo& table_info = m_uin_case_table_list[ai_index++];
            String file_name = bs->getAInput();
            info() << "Applying input acceleration for paraxial element " << face_group.name()
                   << " via CaseTable " << file_name;
            CaseTable* inn = table_info.case_table;

            if (inn != nullptr)
            inn->value(time, ain);
        }

      }
      else{ // For analytical input motions, assuming input displacements
        is_u = true;
        m_inputfunc.m_amplit = bs->getAmplit();
        m_inputfunc.m_coef = bs->getCoef();
        m_inputfunc.m_order = bs->getOrder();
        m_inputfunc.m_tp = bs->getTp();
        m_inputfunc.m_ts = bs->getTs();
        m_inputfunc.m_phase = bs->getPhase();

        auto val{0.};

        switch (typ){
        case 1: val = m_inputfunc.getHarmonic(time); break;
        case 2: val = m_inputfunc.getRicker(time); break;
        case 3: val = m_inputfunc.getDecay(time); break;
        case 4: val = m_inputfunc.getTsang(time); break;
        case 5: val = m_inputfunc.getDirac(time); break;
        default: break; // if user-defined input, do nothing here
        }

        auto norm_angle = bs->getNormalAngle();
        auto plane_angle = bs->getInPlaneAngle();
        auto PI{acos(-1.)};
        auto RAD{PI/180.};
        auto cosan {cos(norm_angle*RAD)};
        auto sinan {sin(norm_angle*RAD)};
        auto cosat {cos(plane_angle*RAD)};
        auto sinat {sin(plane_angle*RAD)};

        if (NDIM == 3) {
            uin.x = sinan*cosat*val;
            uin.y = sinan*sinat*val;
            uin.z = cosan*val;
        }
        else {
            uin.x = sinan*val;
            uin.y = cosan*val;
            uin.z = 0.;
        }
      }

      // Loop on the faces (=edges in 2D) concerned with an input motion condition
      ENUMERATE_FACE (iface, face_group) {

        const Face& face = *iface;

        if (face.isSubDomainBoundary() && face.isOwn()) {
            if (is_u) m_uin_parax[face] = uin;
            if (is_v) m_vin_parax[face] = vin;
            if (is_a) m_ain_parax[face] = ain;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Computes the Jacobian Matrix of a 3D finite-element at Gauss Point ig
Real3x3 NLDynamicModule::
_computeJacobian(const ItemWithNodes& cell,const Int32& ig, const RealUniqueArray& vec, Real& jacobian) {

  auto	n = cell.nbNode();
  Real eps{1.0e-15};

  // Jacobian matrix computed at the integration point
  Real3x3	jac;

  for (Int32 inod = 0, indx = 4; inod < n; ++inod) {

    // vector of local derivatives at this integration point, for node inod
    Real3 dPhi {vec[ig + indx + 1], vec[ig + indx + 2], vec[ig + indx + 3]};
    auto coord_nod = m_node_coord[cell.node(inod)];

    for (int i = 0; i < NDIM; ++i){
      for (int j = 0; j < NDIM; ++j){
        jac[i][j] += dPhi[i] * coord_nod[j];
      }
    }
    indx += 4;
  }

  Int32 ndim = ArcaneFemFunctions::MeshOperation::getGeomDimension(cell);

  if (ndim == 3)
    jacobian = math::matrixDeterminant(jac);

  else if (ndim == 2)
    jacobian = jac.x.x * jac.y.y - jac.x.y * jac.y.x;
  else
    jacobian = ArcaneFemFunctions::MeshOperation::computeLengthEdge2(cell, m_node_coord) / 2.;

  if (fabs(jacobian) < eps) {
    ARCANE_FATAL("Cell jacobian is null");
  }
  return jac;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute elementary mass matrix in 2D at a given Gauss point
void NLDynamicModule::
_computeElemMass(const Real& rho, const DoFLocalId& igauss, const Int32& nb_nodes, RealUniqueArray2& Me){

  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());
  VariableDoFReal& gauss_weight(m_gauss_on_cells.gaussWeight());
  VariableDoFArrayReal& gauss_shape(m_gauss_on_cells.gaussShape());
  VariableDoFReal& gauss_jacobian(m_gauss_on_cells.gaussJacobian());

  auto jacobian = gauss_jacobian[igauss];

  auto wt = gauss_weight[igauss] * jacobian;
  for (Int32 inod = 0; inod < nb_nodes; ++inod) {

    auto Phi_i = gauss_shape[igauss][inod];
    auto rhoPhi_i = wt*rho*Phi_i;

    //----------------------------------------------
    // Elementary Mass (Me) Matrix assembly
    //----------------------------------------------
    for (Int32 jnod = inod ; jnod < nb_nodes; ++jnod) {

      auto Phi_j = gauss_shape[igauss][jnod];
      auto mij = rhoPhi_i*Phi_j;

      for (Int32 l = 0; l < NDIM; ++l){
        auto ii = NDIM*inod + l;
        auto jj = NDIM*jnod + l;
        Me(ii,jj) = mij;
        Me(jj,ii) = mij;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute elementary stiffness matrix in 3D at a given Gauss point
void NLDynamicModule::
_computeK(const Real& lambda, const Real& mu, const DoFLocalId& igauss, const Int32& nb_nodes, RealUniqueArray2& Ke){

  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());
  VariableDoFReal& gauss_weight(m_gauss_on_cells.gaussWeight());
  VariableDoFArrayReal3& gauss_shapederiv(m_gauss_on_cells.gaussShapeDeriv());
  VariableDoFReal3x3& gauss_jacobmat(m_gauss_on_cells.gaussJacobMat());
  VariableDoFReal& gauss_jacobian(m_gauss_on_cells.gaussJacobian());

  auto jacobian = gauss_jacobian[igauss];
  auto a{ lambda + 2.*mu };
  auto size{NDIM * nb_nodes};

  // Setting the "B" matrix size for the max number of nodes in 3D:
  // 8 nodes for a lin element/20 nodes for a quadratic one
  RealUniqueArray2 Bmat = _getB(igauss,nb_nodes);

/*
  RealUniqueArray2 Bmat(NDIM, size);

  for (int i = 0; i <  NDIM; ++i)
    for (int j = 0; j < size; ++j) {
      Bmat(i, j) = 0.;
    }


  // ! Computes the Inverse Jacobian Matrix of a 2D or 3D finite-element
  auto jac = gauss_jacobmat[igauss];
  Real3x3 ijac;

  if (NDIM == 3) {
    ijac = math::inverseMatrix(jac);
  }
  else {
    ijac.x.x = jac.y.y / jacobian;
    ijac.x.y = -jac.x.y / jacobian;
    ijac.y.x = -jac.y.x / jacobian;
    ijac.y.y = jac.x.x / jacobian;
  }
//------------------------------------------------------
// Elementary Derivation Matrix B at current Gauss point
//------------------------------------------------------
for (Int32 inod = 0; inod < nb_nodes; ++inod) {
  auto dPhi = gauss_shapederiv[igauss][inod];
  for (int i = 0; i < NDIM; ++i){
    auto bi{0.};
    for (int j = 0; j < NDIM; ++j) {
      bi += ijac[i][j] * dPhi[j];
    }
    Bmat(i, inod) = bi;
  }
}
  */
  auto wt = gauss_weight[igauss] * jacobian;

  //----------------------------------------------
  // Elementary Stiffness (Ke) Matrix assembly
  //----------------------------------------------
  if (NDIM == 3) {
    for (Int32 inod = 0; inod < nb_nodes; ++inod) {
      for (Int32 l = 0; l < 3; ++l) {

        auto ii = 3 * inod + l;
        RealVector<6> Bii;

        if (!l) {
          Bii(0) = Bmat(0, inod);
          Bii(1) = 0.;
          Bii(2) = 0.;
          Bii(3) = Bmat(1, inod);
          Bii(4) = Bmat(2, inod);
          Bii(5) = 0.;
        }
        else if (l == 1) {
          Bii(0) = 0.;
          Bii(1) = Bmat(1, inod);
          Bii(2) = 0.;
          Bii(3) = Bmat(0, inod);
          Bii(4) = 0.;
          Bii(5) = Bmat(2, inod);
        }
        else if (l == 2) {
          Bii(0) = 0.;
          Bii(1) = 0.;
          Bii(2) = Bmat(2, inod);
          Bii(3) = 0.;
          Bii(4) = Bmat(0, inod);
          Bii(5) = Bmat(1, inod);
        }

        for (Int32 jj = ii; jj < size; ++jj) {

          auto ll = jj % 3;
          RealVector<6> Bjj;

          if (!ll) {
            auto jnod{ (Int32)(jj / 3) };
            Bjj(0) = Bmat(0, jnod);
            Bjj(1) = 0.;
            Bjj(2) = 0.;
            Bjj(3) = Bmat(1, jnod);
            Bjj(4) = Bmat(2, jnod);
            Bjj(5) = 0.;
          }
          else if (ll == 1) {
            auto jnod{ (Int32)((jj - 1) / 3) };
            Bjj(0) = 0.;
            Bjj(1) = Bmat(1, jnod);
            Bjj(2) = 0.;
            Bjj(3) = Bmat(0, jnod);
            Bjj(4) = 0.;
            Bjj(5) = Bmat(2, jnod);
          }
          else if (ll == 2) {
            auto jnod{ (Int32)((jj - 2) / 3) };
            Bjj(0) = 0.;
            Bjj(1) = 0.;
            Bjj(2) = Bmat(2, jnod);
            Bjj(3) = 0.;
            Bjj(4) = Bmat(0, jnod);
            Bjj(5) = Bmat(1, jnod);
          }

          /*------------------------------------------------------------------------------------
// ! Stiffness term (ii,jj) at Gauss point (weight wt) is expressed as:
  Ke(ii, jj) = wt * (Bii(0) * (D(0,0) * Bjj(0) + D(0,1) * Bjj(1) + D(0,2) * Bjj(2)
                            +  D(0,3) * Bjj(3) + D(0,4) * Bjj(4) + D(0,5) * Bjj(5))
                  +  Bii(1) * (D(1,0) * Bjj(0) + D(1,1) * Bjj(1) + D(1,2) * Bjj(2)
                            +  D(1,3) * Bjj(3) + D(1,4) * Bjj(4) + D(1,5) * Bjj(5))
                  +  Bii(2) * (D(2,0) * Bjj(0) + D(2,1) * Bjj(1) + D(2,2) * Bjj(2)
                            +  D(2,3) * Bjj(3) + D(2,4) * Bjj(4) + D(2,5) * Bjj(5))
                  +  Bii(3) * (D(3,0) * Bjj(0) + D(3,1) * Bjj(1) + D(3,2) * Bjj(2)
                            +  D(3,3) * Bjj(3) + D(3,4) * Bjj(4) + D(3,5) * Bjj(5))
                  +  Bii(4) * (D(4,0) * Bjj(0) + D(4,1) * Bjj(1) + D(4,2) * Bjj(2)
                            +  D(4,3) * Bjj(3) + D(4,4) * Bjj(4) + D(4,5) * Bjj(5))
                  +  Bii(5) * (D(5,0) * Bjj(0) + D(5,1) * Bjj(1) + D(5,2) * Bjj(2)
                            +  D(5,3) * Bjj(3) + D(5,4) * Bjj(4) + D(5,5) * Bjj(5)) )

     with elastic tensor D (Dij = Dji):
      D(0,0) = D(1,1) = D(2,2) = lambda + 2.*mu (= a)
      D(3,3) = D(4,4) = D(5,5) = mu
      D(0,1) = D(0,2) = D(1,2) = lambda
      All other terms = 0.
------------------------------------------------------------------------------------*/
          auto kij = wt * (Bii(0) * (a * Bjj(0) + lambda * Bjj(1) + lambda * Bjj(2))
                           + Bii(1) * (lambda * Bjj(0) + a * Bjj(1) + lambda * Bjj(2))
                           + Bii(2) * (lambda * Bjj(0) + lambda * Bjj(1) + a * Bjj(2))
                           + Bii(3) * (mu * Bjj(3))
                           + Bii(4) * (mu * Bjj(4))
                           + Bii(5) * (mu * Bjj(5)));

          Ke(ii, jj) = kij;
          Ke(jj, ii) = kij;
        }
      }
    }
  }
  else{
    for (Int32 inod = 0; inod < nb_nodes; ++inod) {
      for (Int32 l = 0; l < 2; ++l){

        auto ii = 2*inod + l;
        Real3 Bii;

        if (!l){
          Bii.x = Bmat(0, inod);
          Bii.y = 0.;
          Bii.z = Bmat(1, inod);
        } else {
          Bii.x = 0.;
          Bii.y = Bmat(1, inod);
          Bii.z = Bmat(0, inod);
        }

        for (Int32 jj = ii ; jj < size; ++jj) {

          auto ll = jj%2;
          Real3 Bjj{};

          if (!ll){
            auto jnod{ (Int32)(jj / 2) };
            Bjj.x = Bmat(0, jnod);
            Bjj.y = 0.;
            Bjj.z = Bmat(1, jnod);
          }
          else {
            auto jnod{ (Int32)((jj-1) / 2) };
            Bjj.x = 0.;
            Bjj.y = Bmat(1, jnod);
            Bjj.z = Bmat(0, jnod);
          }

          /*------------------------------------------------------------------------------------
// ! Stiffness term (ii,jj) at Gauss point (weight wt) is expressed as:
     Ke(ii, jj) = wt * (Bii.x * (D(0,0) * Bjj.x + D(0,1) * Bjj.y + D(0,2) * Bjj.z)
                     +  Bii.y * (D(1,0) * Bjj.x + D(1,1) * Bjj.y + D(1,2) * Bjj.z)
                     +  Bii.z * (D(2,0) * Bjj.x + D(2,1) * Bjj.y + D(2,2) * Bjj.z) )

     with elastic tensor D (Dij = Dji):
      D(0,0) = D(1,1) = lambda + 2.*mu (= a)
      D(2,2) = mu
      D(0,1) = D(1,0) = lambda
      All other terms = 0.
------------------------------------------------------------------------------------*/

          auto kij = wt * ( Bii.x * (a * Bjj.x + lambda * Bjj.y)
                           +   Bii.y * (lambda * Bjj.x + a * Bjj.y)
                           +   Bii.z * mu * Bjj.z );

          Ke(ii,jj) = kij;
          Ke(jj,ii) = kij;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_computeKParax(const Face& face, const Int32& ig, const RealUniqueArray& vec, const Real& jacobian,
               RealUniqueArray2& Ke, const Real3& RhoC){

  auto dt = m_global_deltat();
  auto alfa{ gamma / beta / dt };
  auto c1{(1. - alfaf) * gamma / beta / dt};

  Real3 ex{ 1., 0., 0. }, ey{ 0., 1., 0. }, ez{ 0., 0., 1. };

  Real3 e1{ m_e1_boundary[face] }, e2{ m_e2_boundary[face] }, e3{ m_e3_boundary[face] };
  // In 2D, paraxial = edge => e1 = tangential vector, e2 = outbound normal vector
  // In 3D, paraxial = face => e1, e2 = on tangential plane, e3 = outbound normal vector
  Real3 nvec{e3};
  if (NDIM < 3)
    nvec = e2;

  Int32 ndim = ArcaneFemFunctions::MeshOperation::getGeomDimension(face);
  auto rhocp{RhoC[ndim]};
  auto rhocs{RhoC[0]};
  auto rhocpcs{rhocp - rhocs};

  // Tensorial product on normal vector nvec:
  Real3x3 nxn({ nvec.x * nvec.x, nvec.x * nvec.y, nvec.x * nvec.z },
              { nvec.y * nvec.x, nvec.y * nvec.y, nvec.y * nvec.z },
              { nvec.z * nvec.x, nvec.z * nvec.y, nvec.z * nvec.z });

  // Loop on the face/edge Gauss points to compute surface integrals terms on the boundary
  Int32 ngauss{ 0 };
  auto wt = vec[ig] * jacobian;

  Int32 nb_nodes = face.nbNode();
  auto size{NDIM * nb_nodes};

  // Loop on nodes of the face or edge
  for (Int32 inod = 0, iig = 4; inod < nb_nodes; ++inod) {

    auto wtPhi_i = wt*vec[ig + iig];
    Node node1 = face.node(inod);

    for (int l = 0; l < NDIM; ++l) {

      auto ii = NDIM * inod + l;
      //----------------------------------------------
      // Elementary contribution c1 * <A0(un+1),v>
      //----------------------------------------------
      for (Int32 jnod = 0, jig = 4; jnod < nb_nodes; ++jnod) {

        auto Phi_j = vec[ig + jig];

        for (int ll = 0; ll < NDIM; ++ll) {
          auto jj = NDIM * jnod + ll;

          auto aij{ rhocpcs * nxn[l][ll] };
          if (l == ll) aij += rhocs;
          auto kij{ c1 * aij * wtPhi_i * Phi_j };

          Ke(ii, jj) = kij;
        }
        jig += 4;
      }
    }
    iig += 4;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Assemble the 2D or 3D bilinear operator (Left Hand Side A matrix)
void NLDynamicModule::
_assembleLinearLHS()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());

  if (NDIM == 3)
    info() << "Assembly of the FEM 3D bilinear operator (LHS - matrix A) ";
  else
    info() << "Assembly of the FEM 2D bilinear operator (LHS - matrix A) ";

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    auto cell_type = cell.type();
    auto nb_nodes{ cell.nbNode() };
    auto nbgauss = m_nb_gauss[cell];

    // Setting the elementary matrices sizes for the max number of nodes * 2 or 3 dofs per node
    auto size{NDIM*nb_nodes};
    RealUniqueArray2 Me(size,size);
    RealUniqueArray2 Ke(size,size);

    for (Int32 i = 0; i < size; ++i) {
      for (Int32 j = i; j < size; ++j) {
        Me(i,j) = 0.;
        Me(j,i) = 0.;
        Ke(i,j) = 0.;
        Ke(j,i) = 0.;
      }
    }

    // Loop on the cell Gauss points to compute integrals terms
    auto cm{(1 - alfam)/beta/dt2};
    auto ck{(1 - alfaf)};
    auto lambda = m_lambda(cell);
    auto mu = m_mu(cell);
    auto rho = m_rho(cell);

    for (Int32 igauss = 0; igauss < nbgauss; ++igauss) {

      DoFLocalId gauss_pti = gauss_point.dofId(cell, igauss);

      // Computing elementary stiffness matrix at Gauss point ig
      _computeK(lambda, mu, gauss_pti, nb_nodes, Ke);

      // Computing elementary mass matrix at Gauss point ig
      _computeElemMass(rho, gauss_pti, nb_nodes, Me);

      // Considering a simple Newmark scheme here (Generalized-alfa will be done later)
      // Computing Me/beta/dt^2 + Ke
      Int32 n1_index{ 0 };
      for (Node node1 : cell.nodes()) {

        auto num1 = node1.uniqueId().asInt32();

        if (node1.isOwn()) {
          for (Int32 iddl = 0; iddl < NDIM; ++iddl) {
            DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
            auto ii = NDIM * n1_index + iddl;

            // Assemble global bilinear operator (LHS)
            Int32 n2_index{ 0 };
            for (Node node2 : cell.nodes()) {

              auto num2 = node2.uniqueId().asInt32();
              for (Int32 jddl = 0; jddl < NDIM; ++jddl) {
                auto node2_dofj = node_dof.dofId(node2, jddl);
                auto jj = NDIM * n2_index + jddl;
                auto mij = Me(ii, jj);
                auto kij = Ke(ii, jj);
                auto aij = cm * mij + ck * kij;

                //              if (node1.isOwn())
                m_linear_system.matrixAddValue(node1_dofi, node2_dofj, aij);
              }
              ++n2_index;
            }
          }
        }
        ++n1_index;
      }
    }
  }
  // Assemble paraxial mass contribution if any
  _assembleLHSParaxialContribution();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Assemble the 2D or 3D linear contributions to the Right Hand Side
//   operator (B vector) - For nonlinear loop, done only at iteration 0
void NLDynamicModule::
_assembleLinearRHS(){
  if (NDIM == 3)
    info() << "Assembly of the FEM 3D linear operator (RHS - vector B) ";
  else
    info() << "Assembly of the FEM 2D linear operator (RHS - vector B) ";

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto dt = m_global_deltat();

  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());
  VariableDoFReal& gauss_weight(m_gauss_on_cells.gaussWeight());
  VariableDoFReal& gauss_jacobian(m_gauss_on_cells.gaussJacobian());
  VariableDoFArrayReal& gauss_shape(m_gauss_on_cells.gaussShape());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    auto cell_type = cell.type();
    auto rho = m_rho(cell);
    auto nb_nodes{ cell.nbNode() };
    auto nbgauss = m_nb_gauss[cell];

    auto size{NDIM*nb_nodes};
    RealUniqueArray2 Me(size,size);

    for (Int32 i = 0; i < size; ++i) {
      for (Int32 j = i; j < size; ++j) {
        Me(i,j) = 0.;
        Me(j,i) = 0.;
      }
    }

    auto nc = cell.uniqueId().asInt32();

    // Loop on the cell Gauss points to compute integrals terms
    auto cm = (1 - alfam)/beta/dt2;

    for (Int32 igauss = 0; igauss < nbgauss; ++igauss) {

      DoFLocalId gauss_pti = gauss_point.dofId(cell, igauss);
      auto jacobian = gauss_jacobian[gauss_pti];

      // Computing elementary mass matrix at Gauss point ig
      _computeElemMass(rho, gauss_pti, nb_nodes, Me);

      // Considering a simple Newmark scheme here (Generalized-alfa will be done later)
      // Computing Me/beta/dt^2 + Ke
      Int32 n1_index{ 0 };
      auto wt = gauss_weight[gauss_pti] * jacobian;

      for (Node node1 : cell.nodes()) {

        if (node1.isOwn()) {
          //---- For debug only !!!
          auto coord1 = m_node_coord[node1];
          auto num1 = node1.uniqueId().asInt32();

          for (Int32 iddl = 0; iddl < NDIM; ++iddl) {
           DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
           auto ii = NDIM * n1_index + iddl;

           bool is_node1_dofi_set = (bool)m_imposed_displ[node1][iddl];
           auto rhs_i{ 0. };

           //          if (node1.isOwn() && !is_node1_dofi_set) {
           if (!is_node1_dofi_set) {

              /*----------------------------------------------------------
            // Mass contribution to the RHS:
            // (1 - alfm)*Mij/(beta*dt2)*uj_pred - alfm*aj_n
            //----------------------------------------------------------*/

              Int32 n2_index{ 0 };
              for (Node node2 : cell.nodes()) {
                //---- For debug only !!!
                auto coord2 = m_node_coord[node2];
                auto num2 = node2.uniqueId().asInt32();

                auto an = m_prev_acc[node2][iddl];
/* cef - Now, this is done with _predictNewmark()
                auto vn = m_prev_vel[node2][iddl];
                auto dn = m_prev_displ[node2][iddl];
                auto u_iddl_pred = dn + dt * vn + dt2 * (0.5 - beta) * an;
*/
                auto u_iddl_pred = m_displ[node2][iddl];
                auto jj = NDIM * n2_index + iddl;
                auto mij = Me(ii, jj);
                rhs_i += mij * (cm * u_iddl_pred - alfam * an);
                ++n2_index;
              }

              /*-------------------------------------------------
            // Other forces (imposed nodal forces, body forces)
            //-------------------------------------------------*/

              {
                //----------------------------------------------
                // Body force terms
                //----------------------------------------------
                auto Phi_i = gauss_shape[gauss_pti][n1_index];
                auto rhoPhi_i = wt * rho * Phi_i;
                rhs_i += rhoPhi_i * gravity[iddl];
              }

              {
                //----------------------------------------------
                // Imposed nodal forces
                //----------------------------------------------
                if ((bool)m_imposed_force[node1][iddl])
                  rhs_i += m_force[node1][iddl];
              }
              rhs_values[node1_dofi] += rhs_i;
           }
          }
        }
        ++n1_index;
      }
    }
  }

  String dirichletMethod = options()->enforceDirichletMethod();
  info() << "Applying Dirichlet boundary condition via "
         << dirichletMethod << " method ";

  dirichletMethod = dirichletMethod.lower();
  // Looking for Dirichlet boundary nodes & modify linear operators accordingly
  ENUMERATE_ (Node, inode, ownNodes()) {
    auto node = *inode;

    for (Int32 iddl = 0; iddl < NDIM; ++iddl) {
      bool is_node_dof_set = (bool)m_imposed_displ[node][iddl];

      if (is_node_dof_set) {
        /*----------------------------------------------------------
            // if Dirichlet node, modify operators (LHS+RHS) allowing to
            // Dirichlet method selected by user
            //----------------------------------------------------------*/
        auto node_dofi = node_dof.dofId(node, iddl);
        auto u_iddl = m_displ[node][iddl];
        if (dirichletMethod == "penalty") {
          m_linear_system.matrixSetValue(node_dofi, node_dofi, penalty);
          rhs_values[node_dofi] = u_iddl * penalty;
        }
        else if (dirichletMethod.contains("weak")) {
          m_linear_system.matrixAddValue(node_dofi, node_dofi, penalty);
          rhs_values[node_dofi] = u_iddl * penalty;
        }
        else if (dirichletMethod.contains("rowelim")) {
          m_linear_system.eliminateRow(node_dofi, u_iddl);
        }
        else if (dirichletMethod.contains("rowcolumnelim")) {
          m_linear_system.eliminateRowColumn(node_dofi, u_iddl);
        }
      }
    }
  }

  //----------------------------------------------
  // Traction contribution to RHS if any
  //----------------------------------------------
  _getTractionContribution(rhs_values);

  //----------------------------------------------
  // Paraxial contribution to RHS if any
  //----------------------------------------------
  _getParaxialContribution(rhs_values);

  //----------------------------------------------
  // Looking for double-couple contributions if any
  //----------------------------------------------
  {
    Real time = globalTime();
    Int32 bd_index{ 0 };
    for (const auto& bd : options()->doubleCouple()) {

      NodeGroup east = bd->getEastNode();
      NodeGroup west = bd->getWestNode();
      NodeGroup north = bd->getNorthNode();
      NodeGroup south = bd->getSouthNode();

      Real Ft{0.};
      auto hasMoment = bd->hasSeismicMomentFile();
      auto hasLoading = bd->hasLoadingFile();

      if (hasMoment || hasLoading) {
        const CaseTableInfo& table_info = m_dc_case_table_list[bd_index];

        if (hasMoment) {
          String file_name = bd->getSeismicMomentFile();
          info() << "Applying the seismic moment for double-couple condition via CaseTable " << file_name;
        }
        else if (hasLoading){
          String file_name = bd->getLoadingFile();
          info() << "Applying the user loading for double-couple condition via CaseTable " << file_name;
        }
        CaseTable* inn = table_info.case_table;
        if (inn != nullptr)
          inn->value(time, Ft);
      }

      auto iplane = bd->getSourcePlane();
      Int32 i1{0}, i2{0};

      if (!iplane) i2 = 1;
      else if (iplane == 1){
        i1 = 1;
        i2 = 2;
      }
      else
        i2 = 2;

      auto is_dew = bd->hasDistEwSeismicMoment();
      auto is_dns = bd->hasDistNsSeismicMoment();
      auto dew = bd->getDistEwSeismicMoment();// East-West distance
      auto dns = bd->getDistNsSeismicMoment();// North-South distance

      ENUMERATE_NODE (inode, west) {

        if (!inode->null() && inode->isOwn()) {
          const Node& dc_node_west = *inode;
          DoFLocalId node_dof_id = node_dof.dofId(dc_node_west, i2);

          //---- For debug only !!!
          auto coords = m_node_coord[dc_node_west];
          auto num = dc_node_west.uniqueId().asInt32();

          rhs_values[node_dof_id] = Ft; // default = hasLoading
          if (hasMoment) {
           if (is_dew && dew != 0.) {
              rhs_values[node_dof_id] /= dew;
           }
           else{
              info() << "EW distance for seismic moment implementation is missing or equal to 0.0! "
                     << "Applying the seismic moment as a user loading";

           }
          }
        }
      }
      ENUMERATE_NODE (inode, east) {
        if (!inode->null() && inode->isOwn()) {
          const Node& dc_node_east = *inode;
          DoFLocalId node_dof_id = node_dof.dofId(dc_node_east, i2);

          //---- For debug only !!!
          auto coords = m_node_coord[dc_node_east];
          auto num = dc_node_east.uniqueId().asInt32();

          rhs_values[node_dof_id] = -Ft;// default = hasLoading
          if (hasMoment) {
           if (is_dew && dew != 0.) {
              rhs_values[node_dof_id] /= dew;
           }
           else{
              info() << "EW distance for seismic moment implementation is missing or equal to 0.0! "
                     << "Applying the seismic moment as a user loading";

           }
          }
        }
      }
      ENUMERATE_NODE (inode, north) {
        if (!inode->null() && inode->isOwn()) {
          const Node& dc_node_north = *inode;
          DoFLocalId node_dof_id = node_dof.dofId(dc_node_north, i1);

          //---- For debug only !!!
          auto coords = m_node_coord[dc_node_north];
          auto num = dc_node_north.uniqueId().asInt32();

          rhs_values[node_dof_id] = Ft;// default = hasLoading
          if (hasMoment) {
           if (is_dns && dns != 0.) {
              rhs_values[node_dof_id] /= dns;
           }
           else{
              info() << "NS distance for seismic moment implementation is missing or equal to 0.0! "
                     << "Applying the seismic moment as a user loading";
           }
          }
        }
      }
      ENUMERATE_NODE (inode, south) {
        if (!inode->null() && inode->isOwn()) {
          const Node& dc_node_south = *inode;
          DoFLocalId node_dof_id = node_dof.dofId(dc_node_south, i1);

          //---- For debug only !!!
          auto coords = m_node_coord[dc_node_south];
          auto num = dc_node_south.uniqueId().asInt32();

          rhs_values[node_dof_id] = -Ft;// default = hasLoading
          if (hasMoment) {
           if (is_dns && dns != 0.) {
              rhs_values[node_dof_id] /= dns;
           }
           else{
              info() << "NS distance for seismic moment implementation is missing or equal to 0.0! "
                     << "Applying the seismic moment as a user loading";
           }
          }
        }
      }

      ++bd_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Compute Elementary Derivation Matrix B at current Gauss point of
 * a given cell with nb_nodes
 */
RealUniqueArray2 NLDynamicModule::
_getB(const DoFLocalId& igauss, const Int32& nb_nodes)
{

  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());
  VariableDoFArrayReal3& gauss_shapederiv(m_gauss_on_cells.gaussShapeDeriv());
  VariableDoFReal3x3& gauss_jacobmat(m_gauss_on_cells.gaussJacobMat());
  VariableDoFReal& gauss_jacobian(m_gauss_on_cells.gaussJacobian());

  auto jacobian = gauss_jacobian[igauss];
  auto size{ NDIM * nb_nodes };

  // The "B" matrix size for the max number of nodes in 3D:
  // 8 nodes for a lin element/20 nodes for a quadratic one
  RealUniqueArray2 Bmat(NDIM, size);

  for (int i = 0; i < NDIM; ++i)
    for (int j = 0; j < size; ++j) {
      Bmat(i, j) = 0.;
    }

  // ! Computes the Inverse Jacobian Matrix of a 2D or 3D finite-element
  auto jac = gauss_jacobmat[igauss];
  Real3x3 ijac;

  if (NDIM == 3) {
    ijac = math::inverseMatrix(jac);
  }
  else {
    ijac.x.x = jac.y.y / jacobian;
    ijac.x.y = -jac.x.y / jacobian;
    ijac.y.x = -jac.y.x / jacobian;
    ijac.y.y = jac.x.x / jacobian;
  }

  /*------------------------------------------------------*/
  /* Elementary Derivation Matrix B at current Gauss point
  */
  for (Int32 inod = 0; inod < nb_nodes; ++inod) {
    auto dPhi = gauss_shapederiv[igauss][inod];
    for (int i = 0; i < NDIM; ++i) {
      auto bi{ 0. };
      for (int j = 0; j < NDIM; ++j) {
        bi += ijac[i][j] * dPhi[j];
      }
      Bmat(i, inod) = bi;
    }
  }
  return Bmat;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Stress computation at each iteration
//
void NLDynamicModule::
_compute_stress(bool init, bool store)
{
  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());

  VariableDoFArrayReal3x3& gauss_stress(m_gauss_on_cells.gaussStress());
  VariableDoFArrayReal3x3& gauss_strain(m_gauss_on_cells.gaussStrain());
  VariableDoFArrayReal3x3& gauss_strain_plastic(m_gauss_on_cells.gaussStrainPlastic());
  VariableDoFArrayReal& gauss_lawparam(m_gauss_on_cells.gaussLawParam());
  VariableDoFArrayReal& gauss_histparam(m_gauss_on_cells.gaussLawHistoryParam());
  VariableDoFArrayReal3x3& gauss_tangent_operator(m_gauss_on_cells.gaussTangentOperator());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto cell_type = cell.type();
    auto cell_nbnod = cell.nbNode();
    Int32 numcell = cell.localId();

    auto cell_nbgauss = ArcaneFemFunctions::FemGaussQuadrature::getNbGaussPointsfromOrder(cell_type, ninteg);

    auto is_default = m_default_law[cell];
    auto lawtyp = static_cast<TypesNLDynamic::eLawType>(m_law[cell]);
    auto ilaw = m_iparam_law[cell];
    LawDispatcher cell_law(lawtyp, is_default);

    for (Int32 ig = 0; ig < cell_nbgauss; ++ig) {
      DoFLocalId gauss_pti = gauss_point.dofId(cell, ig);
      Int32 gaussnum = gauss_pti.localId();

      Tensor2 epsn, sign, epspn;
      epsn.fromReal3x3ToTensor2(gauss_strain[gauss_pti][1]);
      sign.fromReal3x3ToTensor2(gauss_stress[gauss_pti][1]);
      epspn.fromReal3x3ToTensor2(gauss_strain_plastic[gauss_pti][1]);

      Int32 n_index{ 0 };
      Tensor2 deps;

      // Get Elementary Derivation Matrix B at this Gauss point
      // for a given node: Bmat(k, inod) = dkPhi with k=x,y or z direction
      RealUniqueArray2 Bmat = _getB(gauss_pti, cell_nbnod);

      // Compute the strain increment at this Gauss point
      for (Node node : cell.nodes()) {

        if (node.isOwn()) {
          //---- For debug only !!!
          auto coord = m_node_coord[node];
          auto num = node.uniqueId().asInt32();
          Real3 du;

          // init = iter0 =>m_displ = u pred by Newmark
          // !init = iter > 0 : m_displ = current iteration value
          // m_prev_displ = un
          du = m_displ[node] - m_prev_displ[node];

          auto inod{ NDIM * n_index };
          Real3 Bnod;
          for (Int32 i = 0; i < NDIM; i++)
            Bnod[i] = Bmat(i, inod);

          RealVector<6> vdeps;
          vdeps(0) = Bnod.x * du.x; // xx
          vdeps(1) = Bnod.y * du.y; // yy
          vdeps(NDIM) = Bnod.x * du.y + Bnod.y * du.x; // xy

          if (NDIM == 3) {
            vdeps(2) = Bnod.z * du.z; // zz
            vdeps(4) = Bnod.x * du.z + Bnod.z * du.x; // xz
            vdeps(5) = Bnod.y * du.z + Bnod.z * du.y; // yz
          }
          deps.add(vdeps);
        }
        ++n_index;
      }

      cell_law.setLawParams(gauss_lawparam[gauss_pti]);
      cell_law.setLawHistoryParams(gauss_histparam[gauss_pti]);
      cell_law.setStrain(epsn);
      cell_law.setStress(sign);
      cell_law.setPlasticStrain(epspn);
      cell_law.setStrainIncrement(deps);
      cell_law.initState(sign);

      // Compute the current stress from the law at this Gauss point
      // Tangent stiffness operator is not computed at prediction stage
      cell_law.computeStress(init, store);
      gauss_histparam[gauss_pti] = cell_law.updateHistoryVars();

      // Current plastic strains and stresses have been updated by the law
      gauss_strain_plastic[gauss_pti][2] = fromTensor2Real3x3(cell_law.getPlasticStrain());
      gauss_stress[gauss_pti][2] = fromTensor2Real3x3(cell_law.getStress());
      gauss_strain[gauss_pti][2] = fromTensor2Real3x3(cell_law.getStrain());

      if (store){ // when convergence is reached, updating the "prev" quantities for the next step
        gauss_strain_plastic[gauss_pti][1] = gauss_strain_plastic[gauss_pti][2];
        gauss_stress[gauss_pti][1] = gauss_stress[gauss_pti][2];
        gauss_strain[gauss_pti][1] = gauss_strain[gauss_pti][2];

        // To do: tangent operator stored at Gauss points if we need to re-assemble A
        Tensor4 tangent_op = cell_law.getTangentTensor();

        for (Int32 it = 0; it < 3; it++) {
          gauss_tangent_operator[gauss_pti][it] = tangent_op[it];
        }
        if (!tangent_op.isSymmetric())
          gauss_tangent_operator[gauss_pti][3] = tangent_op[3];
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Stress prediction
//
void NLDynamicModule::
_stress_prediction(bool init)
{
  _compute_stress(init, false);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Stress correction
//
void NLDynamicModule::
_stress_correction(){
  _compute_stress(false, m_converge);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Assemble the 2D or 3D nonlinear contributions to the Right Hand Side
//   operator (B vector)
void NLDynamicModule::
_assembleNonLinRHS(){
  if (NDIM == 3)
    info() << "Assembly of the FEM 3D nonlinear contributions to RHS (vector B) ";
  else
    info() << "Assembly of the FEM 2D nonlinear contributions to RHS (vector B) ";

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto dt = m_global_deltat();

  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());
  VariableDoFReal& gauss_weight(m_gauss_on_cells.gaussWeight());
  VariableDoFReal& gauss_jacobian(m_gauss_on_cells.gaussJacobian());
  VariableDoFArrayReal& gauss_shape(m_gauss_on_cells.gaussShape());
  VariableDoFArrayReal3x3& gauss_stress(m_gauss_on_cells.gaussStress());
  VariableDoFArrayReal3x3& gauss_strain(m_gauss_on_cells.gaussStrain());
  VariableDoFArrayReal3x3& gauss_strain_plastic(m_gauss_on_cells.gaussStrainPlastic());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    auto cell_type = cell.type();
    auto rho = m_rho(cell);
    auto nb_nodes{ cell.nbNode() };
    auto nbgauss = m_nb_gauss[cell];

    auto size{NDIM*nb_nodes};
    RealUniqueArray2 Me(size,size);

    for (Int32 i = 0; i < size; ++i) {
      for (Int32 j = i; j < size; ++j) {
        Me(i,j) = 0.;
        Me(j,i) = 0.;
      }
    }

    auto nc = cell.uniqueId().asInt32();

    // Loop on the cell Gauss points to compute integrals terms
    //    auto cm = (1 - alfam)/beta/dt2;
    auto betadt2 = 1/beta/dt2; // Newmark coefficient

    for (Int32 igauss = 0; igauss < nbgauss; ++igauss) {

      DoFLocalId gauss_pti = gauss_point.dofId(cell, igauss);
      auto jacobian = gauss_jacobian[gauss_pti];
      auto Bmat = _getB(gauss_pti,nb_nodes);
      auto sig = gauss_stress[gauss_pti][2];
      auto sig0 = gauss_stress[gauss_pti][0];

      // Computing elementary mass matrix at Gauss point ig
      _computeElemMass(rho, gauss_pti, nb_nodes, Me);

      // Considering a simple Newmark scheme here (Generalized-alfa will be done later)
      // Computing Me/beta/dt^2 + Ke
      Int32 n1_index{ 0 };
      auto wt = gauss_weight[gauss_pti] * jacobian;

      for (Node node1 : cell.nodes()) {

        if (node1.isOwn()) {
          //---- For debug only !!!
          auto coord1 = m_node_coord[node1];
          auto num1 = node1.uniqueId().asInt32();
          Real3 Bnod;
          for (Int32 i = 0; i < NDIM; i++)
            Bnod[i] = Bmat(i, n1_index);

          Real3 Btsig, Btsig0;
          Btsig.x = Bnod.x * sig.x.x + Bnod.y * sig.x.y;
          Btsig.y = Bnod.y * sig.y.y + Bnod.x * sig.x.y;

          if (NDIM == 3){
            Btsig.x += Bnod.z * sig.x.z;
            Btsig.y += Bnod.z * sig.y.z;
            Btsig.z = Bnod.z * sig.z.z + Bnod.x * sig.x.z + Bnod.y * sig.y.z;
          }

          if (!m_deseq){

            Btsig0.x = Bnod.x * sig0.x.x + Bnod.y * sig0.x.y;
            Btsig0.y = Bnod.y * sig0.y.y + Bnod.x * sig0.x.y;

            if (NDIM == 3) {
              Btsig0.x += Bnod.z * sig0.x.z;
              Btsig0.y += Bnod.z * sig0.y.z;
              Btsig0.z = Bnod.z * sig0.z.z + Bnod.x * sig0.x.z + Bnod.y * sig0.y.z;
            }
          }

          for (Int32 iddl = 0; iddl < NDIM; ++iddl) {
            DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
            auto ii = NDIM * n1_index + iddl;

            bool is_node1_dofi_set = (bool)m_imposed_displ[node1][iddl];
            auto rhs_i{ 0. };

            if (!is_node1_dofi_set) {

              /*------------------------------------------------------------------
                 * Stress contribution to the RHS (Newmark only at this stage):
                 * -Bt*sig(ui) with Bt, the transposed derivative operator and sig(ui),
                 * the stresses obtained by strain increments due to previous
                 * iteration displacements
                 * if m_deseq is true, Btsig0 = 0
               */
              rhs_i += wt * (Btsig0[iddl] - Btsig[iddl]);

              Int32 n2_index{ 0 };
              for (Node node2 : cell.nodes()) {
                //---- For debug only !!!
                auto coord2 = m_node_coord[node2];
                auto num2 = node2.uniqueId().asInt32();

                /*------------------------------------------------------------
                 * Mass contribution to the RHS (Newmark only at this stage):
                 * Mij/(beta*dt2)*(u_pred - ui) with ui, the displacement at
                 * previous iteration and u_pred, the predicted displacements
                 * for the current time step
                */
                auto ui = m_prev_displ_iter[node2][iddl];
                auto u_pred = m_displ[node2][iddl];
                auto jj = NDIM * n2_index + iddl;
                auto mij = Me(ii, jj);
                rhs_i += mij * betadt2 * (u_pred - ui);

                ++n2_index;
              }
              rhs_values[node1_dofi] += rhs_i;
            }
          }
        }
        ++n1_index;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Check convergence for the nonlinear iteration loop:
//    Displacements U (= dofs vector): checking if norm(Un+1 - Un) / norm(Un) < Utol
//    Out-of-balance forces (= rhs vector): checking if norm(R) / norm(R0) < Ftol
//    R = residual (iteration i+1 - i) and R0 = norm of rhs at iteration 0
//    Optional: Energy (default = Utol*Ftol): checking if norm([UF]n+1 - [UF]n) / norm([UF]n) < Etol
//
void NLDynamicModule::
_check_convergence(Int32 iter){

}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_getParaxialContribution(Arcane::VariableDoFReal& rhs_values){

  auto dt = m_global_deltat();
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  auto c0{1. - alfaf};
  auto cgb{gamma / beta};
  auto c1{c0 * cgb / dt};
  auto c2{(1. - gamma) * dt};
  auto cc3{(0.5 - beta) * dt2};

  for (const auto& bs : options()->paraxialBoundaryCondition()) {

    FaceGroup face_group = bs->surface();

    bool is_u = bs->hasUInput();
    bool is_v = bs->hasVInput();
    bool is_a = bs->hasAInput();
    Int32 typ{bs->getInputMotionType()};

    // Loop on the faces (=edges in 2D) concerned with the paraxial condition
    ENUMERATE_FACE (iface, face_group) {

      const Face& face = *iface;

      Real3 uin{},vin{}, ain{};
      if (face.isSubDomainBoundary() && face.isOwn()) {

        if (typ > 0){
          // An input motion has been defined
          if (is_u) uin = m_uin_parax[face];
          if (is_v) vin = m_vin_parax[face];
          if (is_a) ain = m_ain_parax[face];
        }
        auto rho = m_rho_parax[face];
        auto RhoC{ rho * m_vel_parax[face] };

        // In 2D, paraxial = edge => e1 = tangential vector, e2 = outbound normal vector
        // In 3D, paraxial = face => e1, e2 = on tangential plane, e3 = outbound normal vector
        Real3 e1{ m_e1_boundary[face] }, e2{ m_e2_boundary[face] }, e3{ m_e3_boundary[face] };

        // Rotation matrix between global and local (paraxial) axes
        Real3x3 ROT({ e1.x, e1.y, e1.z },
                    { e2.x, e2.y, e2.z },
                    { e3.x, e3.y, e3.z });

        auto nb_nodes{ face.nbNode() };
        auto size{ NDIM * nb_nodes };

        // Loop on the cell Gauss points to compute integrals terms
        Int32 ngauss{ 0 };
        auto vec = cell_fem.getGaussData(face, ninteg, ngauss);

        Int32 ng{4 * (1 + nb_nodes)};
        for (Int32 igauss = 0, ig = 0; igauss < ngauss; ++igauss, ig += ng) {

          auto jacobian{ 0. };
          _computeJacobian(face, ig, vec, jacobian);

          // Loop on nodes of the paraxial face (with no Dirichlet condition)
          Int32 n1_index{ 0 };
          auto iig{ 4 };
          auto wt = vec[ig] * jacobian;
          Real3 a0{};

          for (Node node : face.nodes()) {

              auto Phi_i = vec[ig + iig];
              auto vi_pred = m_prev_vel[node] + c2 * m_prev_acc[node];
              auto ui_pred = m_prev_displ[node] + dt * m_prev_vel[node] + cc3 * m_prev_acc[node];
              auto vni = m_prev_vel[node];

              for (Int32 i = 0; i < NDIM; ++i) {
                auto vi{ 0. };

                for (Int32 j = 0; j < NDIM; ++j) {
                  vi += ROT[i][j] * (-c0 * vi_pred[j] + c1 * ui_pred[j] - alfaf * vni[j]);
                }
                a0[i] += Phi_i * RhoC[i] * vi;
              }
            iig += 4;
          }

          iig = 4;
          for (Node node : face.nodes()) {

            if (node.isOwn()) {
              //---- For debug only !!!
              auto num = node.uniqueId().asInt32();
              auto coords = m_node_coord[node];

              auto Phi_i = vec[ig + iig];
              auto wtPhi_i = wt * Phi_i;
              for (Int32 iddl = 0; iddl < NDIM; ++iddl) {
                DoFLocalId node_dofi = node_dof.dofId(node, iddl);

                bool is_node_dofi_set = (bool)m_imposed_displ[node][iddl];
                auto rhs_i{ 0. };

                if (!is_node_dofi_set) {

                  for (Int32 j = 0; j < NDIM; ++j) {
                    rhs_i += ROT[iddl][j] * a0[j];
                  }
                }
                rhs_values[node_dofi] += wtPhi_i * rhs_i;
              }
              iig += 4;
            }
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_assembleLHSParaxialContribution(){

  auto dt = m_global_deltat();
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto c1{(1. - alfaf) * gamma / beta / dt};

  for (const auto& bs : options()->paraxialBoundaryCondition()) {

    FaceGroup face_group = bs->surface();
    //      info() << "Applying constant paraxial boundary conditions for surface " << face_group.name();

    // Loop on the faces (=edges in 2D) concerned with the paraxial condition
    ENUMERATE_FACE (iface, face_group) {

      const Face& face = *iface;

      if (face.isSubDomainBoundary() && face.isOwn()) {

        auto rho = m_rho_parax[face];
        auto RhoC{ rho * m_vel_parax[face] };

        // In 3D, a quadratic face element has max 9 nodes (27 dofs)
        auto nb_nodes{face.nbNode()};
        auto size{ NDIM * nb_nodes};
        RealUniqueArray2 Ke(size,size);

        for (Int32 i = 0; i < size; ++i) {
          for (Int32 j = i; j < size; ++j) {
            Ke(i, j) = 0.;
            Ke(j, i) = 0.;
          }
        }

        // Loop on the cell Gauss points to compute integrals terms
        Int32 ngauss{ 0 };
        auto vec = cell_fem.getGaussData(face, ninteg, ngauss);

        Int32 ng{ 4 * (1 + nb_nodes) };
        for (Int32 igauss = 0, ig = 0; igauss < ngauss; ++igauss, ig += ng) {


          auto jacobian{ 0. };
          _computeJacobian(face, ig, vec, jacobian);
          _computeKParax(face, ig, vec, jacobian, Ke, RhoC);

          // Loop on nodes of the face (with no Dirichlet condition)
          Int32 n1_index{ 0 };
          auto iig{ 4 };
          for (Node node1 : face.nodes()) {

            for (Int32 iddl = 0; iddl < NDIM; ++iddl) {

              DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
              auto ii = NDIM * n1_index + iddl;

              if (node1.isOwn()){
                //----------------------------------------------
                // Elementary contribution to LHS
                //----------------------------------------------
                Int32 n2_index{ 0 };
                for (Node node2 : face.nodes()) {
                  for (Int32 jddl = 0; jddl < NDIM; ++jddl) {
                    auto node2_dofj = node_dof.dofId(node2, jddl);
                    auto jj = NDIM * n2_index + jddl;
                    auto mij = Ke(ii, jj);
                    m_linear_system.matrixAddValue(node1_dofi, node2_dofj, mij);
                  }
                  ++n2_index;
                }
              }
            }
            ++n1_index;
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_getTractionContribution(Arcane::VariableDoFReal& rhs_values){

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  for (const auto& bs : options()->neumannCondition()) {
    FaceGroup face_group = bs->surface();

    // Loop on the faces (=edges in 2D) concerned with the traction condition
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;

      Real3 trac = m_imposed_traction[face];
      auto facint = ArcaneFemFunctions::MeshOperation::computeFacLengthOrArea(face,m_node_coord);

      // Loop on nodes of the face or edge (with no Dirichlet condition)
      ENUMERATE_NODE (k, face.nodes()){
        const Node& node = *k;
        auto coord = m_node_coord[node];
        auto num = node.uniqueId();

        for (Int32 iddl = 0; iddl < NDIM; ++iddl)
          //          if (!(bool)m_imposed_displ[node][iddl] && node.isOwn()) {
          if (node.isOwn()) {
            DoFLocalId dof_id = node_dof.dofId(node, iddl);
            rhs_values[dof_id] += trac[iddl] * facint;
          }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void NLDynamicModule::
_doSolve(){
  info() << "Solving Linear system";
  m_linear_system.solve();

  {
    VariableDoFReal& dof_d(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;

      auto num = node.uniqueId().asInt32();
      auto coord = m_node_coord[node];

      auto ux = dof_d[node_dof.dofId(node, 0)];
      auto uy = dof_d[node_dof.dofId(node, 1)];
      auto uz{0.};

      if (NDIM == 3)
        uz = dof_d[node_dof.dofId(node, 2)];

      m_displ[node] = Real3(ux,uy,uz);

//      info() << "Node: " << node.localId() << " Ux=" << ux << " Uy=" << uy << " Uz=" << uz;
    }
  }

  // Re-Apply Dirichlet boundary conditions because the solver has modified the values
  // on all nodes
  _applyDirichletBoundaryConditions();// --- Check if it is required (re-apply paraxial conditions too?)

  m_displ.synchronize();
  m_vel.synchronize();
  m_acc.synchronize();
  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    long p = std::cout.precision();
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
//      info() << "Node: " << node.uniqueId() << " Ux=" << m_displ[node].x << " Uy=" << m_displ[node].y << " Uz=" << m_displ[node].z;
      info() << node.uniqueId() << " " << m_displ[node].x << " " << m_displ[node].y << " " << m_displ[node].z;
    }
    std::cout.precision(p);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ARCANE_REGISTER_MODULE_NLDYNAMIC(NLDynamicModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
