// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* MechanicsModule.cc                                          (C) 2022-2024 */
/*                                                                           */
/* Quasi-static implicit mechanical FEM solver                               */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "arcane/MathUtils.h"
#include <arcane/utils/MultiArray2.h>
#include "arcane/utils/ArgumentException.h"
#include <arcane/IParallelMng.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/geometry/IGeometry.h>
#include <arcane/IIOMng.h>
#include <arcane/CaseTable.h>

#include "IDoFLinearSystemFactory.h"
#include "MechanicsModule.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MechanicsModule::MechanicsModule(const ModuleBuildInfo& mbi)
: ArcaneMechanicsObject(mbi)
, m_dofs_on_nodes(mbi.subDomain()->traceMng())
, m_gauss_on_cells(mbi.subDomain()->traceMng())
{
  ICaseMng *cm = mbi.subDomain()->caseMng();
  cm->setTreatWarningAsError(true);
  cm->setAllowUnkownRootElelement(false);
}

VersionInfo MechanicsModule::versionInfo() const {
  return {1, 0, 0};
}

MechanicsModule::~MechanicsModule()
{
  for( const CaseTableInfo&  t : m_traction_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_sdispl_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_sforce_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_displ_case_table_list )
    delete t.case_table;
  for( const CaseTableInfo&  t : m_force_case_table_list )
    delete t.case_table;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
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
  theta = options()->getTheta();
  auto dt = options()->getDt();
  m_global_deltat = dt;
  dt2 = dt * dt;
  auto tf = options()->getTf();
  m_global_final_time = tf;
  auto t = options()->getT0();
  m_global_time = t;
  linop_nstep = options()->getLinopNstep();
  auto szType = options()->initElastType().lower();
  if (szType.contains("young")) elast_type = TypesMechanics::YoungNu;
  else if (szType.contains("lame")) elast_type = TypesMechanics::Lame;
  else if (szType.contains("bulk")) elast_type = TypesMechanics::Bulk;
  else {

    info() << "init-elast-type keywords must include (not case dependent):\n"
           << "  - Young\n"
           << "  - Lame\n"
           << "  - Bulk\n";
    ARCANE_FATAL("Type for elastic properties is undefined!");
  }
  auto nsteps = (int)((tf - t)/dt);
  if (linop_nstep > nsteps) keep_constop = true;

  analysis_type = options()->getAnalysisType();
  if (analysis_type == TypesMechanics::ThreeD)
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
    m_prev_u[node] = Real3::zero();
    m_u[node] = Real3::zero();
  }

  _startInitGauss();
  _applyInitialNodeConditions();
  _initCells();
  _initBoundaryConditions();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_initDofs(){
  m_dofs_on_nodes.initialize(mesh(),NDIM);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_startInitGauss(){
  Integer max_gauss_per_cell{0};

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto cell_type = cell.type();
    auto nbgauss = getNbGaussPointsfromOrder(cell_type, ninteg);
    m_nb_gauss[cell] = nbgauss;
    max_gauss_per_cell = math::max(nbgauss,max_gauss_per_cell);
  }
  m_gauss_on_cells.initialize(mesh(),max_gauss_per_cell);

  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());
  VariableDoFReal& gauss_weight(m_gauss_on_cells.gaussWeight());
  VariableDoFReal3& gauss_refpos(m_gauss_on_cells.gaussRefPosition());
  VariableDoFArrayReal& gauss_shape(m_gauss_on_cells.gaussShape());
  VariableDoFArrayReal3& gauss_shapederiv(m_gauss_on_cells.gaussShapeDeriv());

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto cell_type = cell.type();
    auto cell_nbnod = cell.nbNode();
    auto cell_nbgauss = getNbGaussPointsfromOrder(cell_type, ninteg);
    Int32 ndim = getGeomDimension(cell);

    for (Int32 ig = 0; ig < cell_nbgauss; ++ig) {
      DoFLocalId gauss_pti = gauss_point.dofId(cell,ig);
      gauss_weight[gauss_pti] = getGaussWeight(cell,ninteg,ig);
      Real3 refpos = getGaussRefPosition(cell,ninteg,ig);
      if (ndim <= 2) {
        refpos.z = 0.;
        if (ndim == 1) refpos.y = 0.;
      }
      gauss_refpos[gauss_pti] = refpos;
      gauss_shape[gauss_pti] = RealUniqueArray(cell_nbnod);
      gauss_shapederiv[gauss_pti] = Real3UniqueArray(cell_nbnod);

      for (Int32 inod = 0; inod < cell_nbnod; ++inod) {
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
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_initCells(){
  Real E, nu, lambda, mu, K;

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto rho = m_rho[cell];
    if (elast_type == TypesMechanics::YoungNu) {
      E = m_young[cell];
      nu = m_nu[cell];
      lambda = nu*E/(1. + nu)/(1. - 2.*nu);
      mu = E/2./(1. + nu);
      K = lambda + 2./3. *mu;

    } else if (elast_type == TypesMechanics::Lame) {
      lambda = m_lambda[cell];
      mu = m_mu[cell];
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);
      K = lambda + 2./3. *mu;

    } else if (elast_type == TypesMechanics::Bulk) {
      K = m_k[cell];
      mu = m_mu[cell];
      lambda = K- 2./3. *mu;
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);

    }
    m_lambda[cell] = lambda;
    m_mu[cell] = mu;
    m_young[cell] = E;
    m_nu[cell] = nu;
    m_k[cell] = K;
  }

  _applyInitialCellConditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_applyInitialNodeConditions(){

  for (Int32 i = 0, nb = options()->initialNodeCondition().size(); i < nb; ++i) {

    NodeGroup node_group = options()->initialNodeCondition[i]->nodeGroup();
    info() << "Applying init-node-condition for node group " << node_group.name();

    // Loop on nodes with this initial condition
    ENUMERATE_NODE(inode, node_group) {
      const Node & node = *inode;

      if (options()->initialNodeCondition[i]->hasU())
        m_prev_u[node] = options()->initialNodeCondition[i]->U();

      if (options()->initialNodeCondition[i]->hasF())
        m_f[node] = options()->initialNodeCondition[i]->F();
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_applyInitialCellConditions(){

  for (Integer i = 0, nb = options()->initElastProperties().size(); i < nb; ++i) {

    CellGroup cell_group = options()->initElastProperties[i]->cellGroup();

    // In the future, we will have to find a way to impose different initial
    // properties (stress/strain tensors, densities...) per element from a file
    // (e.g., coming from a previous computation)
    auto rho = options()->initElastProperties[i]->rho();
    Real vp, vs, E, nu, lambda, mu, K;

    if (elast_type == TypesMechanics::YoungNu) {
      E = options()->initElastProperties[i]->young();
      nu = options()->initElastProperties[i]->nu();
      lambda = nu*E/(1. + nu)/(1. - 2.*nu);
      mu = E/2./(1. + nu);
      K = lambda + 2./3. *mu;

    } else if (elast_type == TypesMechanics::Lame) {
      lambda = options()->initElastProperties[i]->young();
      mu = options()->initElastProperties[i]->nu();
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);
      K = lambda + 2./3. *mu;

    } else if (elast_type == TypesMechanics::Bulk) {
      K = options()->initElastProperties[i]->k();
      mu = options()->initElastProperties[i]->nu();
      lambda = K - 2./3. * mu;
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);

    }

    ENUMERATE_CELL (icell, cell_group) {
      const Cell& cell = *icell;
      m_rho[cell] = rho;
      m_k[cell] = K;
      m_lambda[cell] = lambda;
      m_mu[cell] = mu;
      m_young[cell] = E;
      m_nu[cell] = nu;
    }
  }

  for (Integer i = 0, nb = options()->initCellCondition().size(); i < nb; ++i) {

    CellGroup cell_group = options()->initCellCondition[i]->cellGroup();
    info() << "Applying init-cell-condition for cell group " << cell_group.name();

    // In the future, we will have to find a way to impose different initial
    // properties (stress/strain tensors, densities...) per element from a file
    // (e.g., coming from a previous computation)

    bool  hasepsv{options()->initCellCondition[i]->hasEpsvol()},
          hasepsd{ options()->initCellCondition[i]->hasEpsdev()},
          hassigv{ options()->initCellCondition[i]->hasSigvol()},
          hassigd{ options()->initCellCondition[i]->hasSigdev()};

    // Loop on cells with this initial condition
    ENUMERATE_CELL (icell, cell_group) {
      const Cell& cell = *icell;

      // Initialize the stress tensor for the concerned cell
      if (hasepsd)
        m_epsdev[cell] = options()->initCellCondition[i]->epsdev();
      if (hasepsv)
        m_epsvol[cell] = options()->initCellCondition[i]->epsvol();
      if (hassigd)
        m_sigdev[cell] = options()->initCellCondition[i]->sigdev();
      if (hassigv)
        m_sigvol[cell] = options()->initCellCondition[i]->sigvol();

    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_initGaussStep()
{
  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());

  VariableDoFReal& gauss_jacobian(m_gauss_on_cells.gaussJacobian());
  VariableDoFReal3& gauss_refpos(m_gauss_on_cells.gaussRefPosition());
  VariableDoFReal3x3& gauss_jacobmat(m_gauss_on_cells.gaussJacobMat());
  VariableDoFArrayReal& gauss_shape(m_gauss_on_cells.gaussShape());
  VariableDoFArrayReal3& gauss_shapederiv(m_gauss_on_cells.gaussShapeDeriv());

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto cell_type = cell.type();
    auto cell_nbnod = cell.nbNode();
    Int32 ndim = getGeomDimension(cell);
    auto cell_nbgauss = getNbGaussPointsfromOrder(cell_type, ninteg);

    for (Int32 ig = 0; ig < cell_nbgauss; ++ig) {
      DoFLocalId gauss_pti = gauss_point.dofId(cell, ig);
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
        jacobian = Line2Length(cell, m_node_coord) / 2.;

      if (fabs(jacobian) < REL_PREC) {
        ARCANE_FATAL("Cell jacobian is null");
      }
      gauss_jacobian[gauss_pti] = jacobian;
      gauss_jacobmat[gauss_pti] = jac;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
compute(){

  info() << "Module MECHANICS COMPUTE";
  ++linop_nstep_counter;

  // Stop code at exact final time set by user
  auto tf = m_global_final_time();
  auto t = m_global_time();
  auto dt = m_global_deltat();
  auto t0 = options()->getT0();
  dt2 = dt * dt;

  info() << "Time (s) = " << t;

  // Set if we want to keep the matrix structure between calls
  // the rate is a user input (linop_nstep)
  // The matrix has to have the same structure (same structure for non-zero)

  if (m_linear_system.isInitialized() && (linop_nstep_counter < linop_nstep || keep_constop)){
    m_linear_system.clearValues();
  }
  else {

    m_linear_system.reset();
    m_linear_system.setLinearSystemFactory(options()->linearSystem());
    m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

    // Reset the counter when the linear operator is reset
    linop_nstep_counter = 0;
  }

  // Update Gauss shape functions and their derivatives for this step
  _initGaussStep();

  // Apply Dirichlet/Neumann conditions if any
  _applyDirichletBoundaryConditions();
  _applyNeumannBoundaryConditions();

  // Assemble the FEM global operators (LHS matrix/RHS vector b)
  _assembleLinearLHS();

  // Nonlinear Loop
  _assembleLinearRHS();

  // Solve the linear system AX = B
  _doSolve();

  // Save/Check results
  //  _checkResultFile();

  if (t < tf) {
    if (t + dt > tf) {
      dt = tf - t;
      m_global_deltat = dt;
    }
  }
  else
    subDomain()->timeLoopMng()->stopComputeLoop(true);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_update(){

  // Updating the nodal accelerations and velocities (after solve) with
  auto dt = m_global_deltat();

  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    auto dn = m_prev_u[node];
    auto du = m_du[node];

      for (Int32 i = 0; i < NDIM; ++i) {
        if (!(bool)m_u_imp[node][i]) {
          m_u[node][i] = dn[i] + du[i];
        }
      }
      m_prev_u[node] = m_u[node];
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_initBoundaryConditions()
{
  IParallelMng* pm = subDomain()->parallelMng();

  for (const auto& bd : options()->dirichletSurfaceCondition()) {
    FaceGroup face_group = bd->surface();
    info() << "Initializing Dirichlet Surface Boundary Conditions for face group " << face_group.name();

    if (bd->hasUCurve()) {
      String file_name = bd->UCurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_sdispl_case_table_list.add(CaseTableInfo{ file_name, case_table });
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

        m_u_imp[node].x = (bd->hasUx() || (hasUcurve && xdir) ? 1 : 0);
        m_u_imp[node].y = (bd->hasUy() || (hasUcurve && ydir) ? 1 : 0);
        m_u_imp[node].z = (bd->hasUz() || (hasUcurve && zdir) ? 1 : 0);

        m_f_imp[node].x = (bd->hasFx() || (hasFcurve && xdir) ? 1 : 0);
        m_f_imp[node].y = (bd->hasFy() || (hasFcurve && ydir) ? 1 : 0);
        m_f_imp[node].z = (bd->hasFz() || (hasFcurve && zdir) ? 1 : 0);
      }
    }
  }

  for (const auto& bd : options()->dirichletPointCondition()) {
    NodeGroup nodes = bd->node();
    info() << "Initializing Dirichlet Point Conditions for node group " << nodes.name();

    if (bd->hasUCurve()) {
      String file_name = bd->UCurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_displ_case_table_list.add(CaseTableInfo{ file_name, case_table });
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
    auto hasFcurve{bd->hasFCurve()};
    auto xdir{bd->getXAxis()};
    auto ydir{bd->getYAxis()};
    auto zdir{bd->getZAxis()};

    // Loop on nodes
    ENUMERATE_NODE (inode, nodes) {
      const Node& node = *inode;
      auto coord = m_node_coord[node];
      auto num = node.uniqueId();

      m_f_imp[node].x = (bd->hasFx() || (hasFcurve && xdir) ? 1 : 0);
      m_f_imp[node].y = (bd->hasFy() || (hasFcurve && ydir) ? 1 : 0);
      m_f_imp[node].z = (bd->hasFz() || (hasFcurve && zdir) ? 1 : 0);

      if (bd->hasUx() || (hasUcurve && xdir))
        m_u_imp[node].x = 1;

      if (bd->hasUy() || (hasUcurve && ydir))
        m_u_imp[node].y = 1;

      if (bd->hasUz() || (hasUcurve && zdir))
        m_u_imp[node].z = 1;
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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_applyDirichletBoundaryConditions(){

  Real time = globalTime();
  Int32 suc_index{ 0 }, sfc_index{ 0 };
  for (const auto& bd : options()->dirichletSurfaceCondition()) {
    FaceGroup face_group = bd->surface();

    Real3 u{};
    bool is_u_imp{bd->hasUCurve() || bd->hasUx() || bd->hasUy() || bd->hasUz()};
    if (bd->hasUCurve()) {
      const CaseTableInfo& table_info = m_sdispl_case_table_list[suc_index++];
      String file_name = bd->UCurve();
      info() << "Applying displacement boundary conditions for surface " << face_group.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, u);
    }
    else if (is_u_imp){
      if (bd->hasUx())
        u.x = bd->getUx();
      if (bd->hasUy())
        u.y = bd->getUy();
      if (bd->hasUz())
        u.z = bd->getUz();
    }

    Real3 force{};
    bool is_f_imp{bd->hasFCurve() || bd->hasFx() || bd->hasFy() || bd->hasFz()};
    if (bd->hasFCurve()) {
      const CaseTableInfo& table_info = m_sforce_case_table_list[sfc_index++];
      String file_name = bd->FCurve();
      info() << "Applying force boundary conditions for surface " << face_group.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, force);
    }
    else if (is_f_imp){
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

        if (is_u_imp) {
          if ((bool)m_u_imp[node].x)
            m_u[node].x = u.x;

          if ((bool)m_u_imp[node].y)
            m_u[node].y = u.y;

          if ((bool)m_u_imp[node].z)
            m_u[node].z = u.z;
        }

        if (is_f_imp) {
          if ((bool)m_f_imp[node].x)
            m_f[node].x = force.x;

          if ((bool)m_f_imp[node].y)
            m_f[node].y = force.y;

          if ((bool)m_f_imp[node].z)
            m_f[node].z = force.z;
        }
      }
    }
  }

  Int32 uc_index{ 0 }, fc_index{ 0 };
  for (const auto& bd : options()->dirichletPointCondition()) {
    NodeGroup nodes = bd->node();

    Real3 u{};
    bool is_u_imp{bd->hasUCurve() || bd->hasUx() || bd->hasUy() || bd->hasUz()};
    if (bd->hasUCurve()) {
      const CaseTableInfo& table_info = m_displ_case_table_list[uc_index++];
      String file_name = bd->UCurve();
      info() << "Applying displacement boundary conditions for nodes " << nodes.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, u);
    }
    else if (is_u_imp){
      if (bd->hasUx())
        u.x = bd->getUx();
      if (bd->hasUy())
        u.y = bd->getUy();
      if (bd->hasUz())
        u.z = bd->getUz();
    }

    Real3 force{};
    bool is_f_imp{bd->hasFCurve() || bd->hasFx() || bd->hasFy() || bd->hasFz()};
    if (bd->hasFCurve()) {
      const CaseTableInfo& table_info = m_force_case_table_list[fc_index++];
      String file_name = bd->FCurve();
      info() << "Applying force boundary conditions for nodes " << nodes.name()
             << " via CaseTable " << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(time, force);
    }
    else if (is_f_imp){
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

       if (is_u_imp) {
        if ((bool)m_u_imp[node].x)
          m_u[node].x = u.x;

        if ((bool)m_u_imp[node].y)
          m_u[node].y = u.y;

        if ((bool)m_u_imp[node].z)
          m_u[node].z = u.z;
      }

      if (is_f_imp) {
        if ((bool)m_f_imp[node].x)
          m_f[node].x = force.x;

        if ((bool)m_f_imp[node].y)
          m_f[node].y = force.y;

        if ((bool)m_f_imp[node].z)
          m_f[node].z = force.z;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
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
      if (bs->hasTx())
        trac.x = bs->getTx();
      if (bs->hasTy())
        trac.y = bs->getTy();
      if (bs->hasTz())
        trac.z = bs->getTz();
    }

    // Loop on faces of the surface
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;
      m_t_imp[face] = trac;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Computes the Jacobian Matrix of a 3D finite-element at Gauss Point ig
Real3x3 MechanicsModule::
_computeJacobian(const ItemWithNodes& cell,const Int32& ig, const RealUniqueArray& vec, Real& jacobian) {

  auto	n = cell.nbNode();

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

  Int32 ndim = getGeomDimension(cell);
  //  if (NDIM == 3)
  if (ndim == 3)
    jacobian = math::matrixDeterminant(jac);

  //  else if (NDIM == 2) {
  else if (ndim == 2)
    jacobian = jac.x.x * jac.y.y - jac.x.y * jac.y.x;
  else
    jacobian = Line2Length(cell, m_node_coord) / 2.;

  if (fabs(jacobian) < REL_PREC) {
    ARCANE_FATAL("Cell jacobian is null");
  }
  return jac;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute elementary stiffness matrix in 3D at a given Gauss point
void MechanicsModule::
_computeK(const Real& lambda, const Real& mu, const DoFLocalId& igauss, const Int32& nb_nodes, RealUniqueArray2& Ke){

  auto gauss_point(m_gauss_on_cells.gaussCellConnectivityView());
  VariableDoFReal& gauss_weight(m_gauss_on_cells.gaussWeight());
  VariableDoFArrayReal3& gauss_shapederiv(m_gauss_on_cells.gaussShapeDeriv());
  VariableDoFReal3x3& gauss_jacobmat(m_gauss_on_cells.gaussJacobMat());
  VariableDoFReal& gauss_jacobian(m_gauss_on_cells.gaussJacobian());

  auto jacobian = gauss_jacobian[igauss];

  auto size{NDIM * nb_nodes};

  // Setting the "B" matrix size for the max number of nodes in 3D:
  // 8 nodes for a lin element/20 nodes for a quadratic one
  RealUniqueArray2 Bmat(NDIM, size);

  auto a{ lambda + 2.*mu };

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

  auto wt = gauss_weight[igauss] * jacobian;

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

  //----------------------------------------------
  // Elementary Stiffness (Ke) Matrix assembly
  //----------------------------------------------
  if (NDIM == 3) {
    for (Int32 inod = 0; inod < nb_nodes; ++inod) {
      for (Int32 l = 0; l < 3; ++l) {

        auto ii = 3 * inod + l;
        FixedVector<6> Bii;

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
          FixedVector<6> Bjj;

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
          auto kij = wt * (Bii(0) * (a * Bjj(0) + lambda * Bjj(1) + lambda * Bjj(2)) + Bii(1) * (lambda * Bjj(0) + a * Bjj(1) + lambda * Bjj(2)) + Bii(2) * (lambda * Bjj(0) + lambda * Bjj(1) + a * Bjj(2)) + Bii(3) * (mu * Bjj(3)) + Bii(4) * (mu * Bjj(4)) + Bii(5) * (mu * Bjj(5)));

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
// ! Assemble the 2D or 3D bilinear operator (Left Hand Side A matrix)
void MechanicsModule::
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
    RealUniqueArray2 Ke(size,size);

    for (Int32 i = 0; i < size; ++i) {
      for (Int32 j = i; j < size; ++j) {
        Ke(i,j) = 0.;
        Ke(j,i) = 0.;
      }
    }

    // Loop on the cell Gauss points to compute integrals terms
    auto lambda = m_lambda(cell);
    auto mu = m_mu(cell);
    auto rho = m_rho(cell);

    for (Int32 igauss = 0; igauss < nbgauss; ++igauss) {

      DoFLocalId gauss_pti = gauss_point.dofId(cell, igauss);

      // Computing elementary stiffness matrix at Gauss point ig
      _computeK(lambda, mu, gauss_pti, nb_nodes, Ke);

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
                auto kij = Ke(ii, jj);

                //              if (node1.isOwn())
                m_linear_system.matrixAddValue(node1_dofi, node2_dofj, kij);
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Assemble the 2D or 3D linear operator (Right Hand Side B vector)
void MechanicsModule::
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
    for (Int32 igauss = 0; igauss < nbgauss; ++igauss) {

      DoFLocalId gauss_pti = gauss_point.dofId(cell, igauss);
      auto jacobian = gauss_jacobian[gauss_pti];

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

           bool is_node1_dofi_set = (bool)m_u_imp[node1][iddl];
           auto rhs_i{ 0. };

           //          if (node1.isOwn() && !is_node1_dofi_set) {
           if (!is_node1_dofi_set) {

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
                if ((bool)m_f_imp[node1][iddl])
                  rhs_i += m_f[node1][iddl];
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
      bool is_node_dof_set = (bool)m_u_imp[node][iddl];

      if (is_node_dof_set) {
        /*----------------------------------------------------------
            // if Dirichlet node, modify operators (LHS+RHS) allowing to
            // Dirichlet method selected by user
            //----------------------------------------------------------*/
        auto node_dofi = node_dof.dofId(node, iddl);
        auto u_iddl = m_du[node][iddl];
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

}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
_getTractionContribution(Arcane::VariableDoFReal& rhs_values){

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  for (const auto& bs : options()->neumannCondition()) {
    FaceGroup face_group = bs->surface();

    // Loop on the faces (=edges in 2D) concerned with the traction condition
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;

      Real3 trac = m_t_imp[face];
      auto facint = _computeFacLengthOrArea(face);

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
Real MechanicsModule::
_computeFacLengthOrArea(const Face& face)
{
  Int32 item_type = face.type();
  Real fac_el{0.};

  switch (item_type) {

  // Lines
  case IT_Line2:
  case IT_Line3:
    fac_el = Line2Length(face, m_node_coord) / 2.;
    break;

  // Faces
  case IT_Triangle3:
  case IT_Triangle6:
    fac_el = Tri3Surface(face, m_node_coord) / 3.;
    break;

  case IT_Quad4:
  case IT_Quad8:
    fac_el = Quad4Surface(face, m_node_coord) / 4.;
    break;

  default:
    break;
  }
  return fac_el;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MechanicsModule::
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

      auto dux = dof_d[node_dof.dofId(node, 0)];
      auto duy = dof_d[node_dof.dofId(node, 1)];
      auto duz{0.};

      if (NDIM == 3)
        duz = dof_d[node_dof.dofId(node, 2)];

      m_du[node] = Real3(dux,duy,duz);

      info() << "Node: " << node.localId() << " dUx=" << dux << " dUy=" << duy << " dUz=" << duz;
    }
  }

  // Re-Apply Dirichlet boundary conditions because the solver has modified the values
  // on all nodes
  _applyDirichletBoundaryConditions();// --- Check if it is required (re-apply paraxial conditions too?)

  m_du.synchronize();
  // Update the nodal variable according to the integration scheme (e.g. Newmark)
  _update();

  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    long p = std::cout.precision();
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      info() << "Node: " << node.uniqueId() << " Ux=" << m_u[node].x << " Uy=" << m_u[node].y << " Uz=" << m_u[node].z;
    }
    std::cout.precision(p);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ARCANE_REGISTER_MODULE_MECHANICS(MechanicsModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
