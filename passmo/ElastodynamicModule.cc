#include <arcane/utils/FatalErrorException.h>
#include "arcane/MathUtils.h"
#include <arcane/utils/NumArray.h>
#include <arcane/utils/MultiArray2.h>
#include "arcane/utils/ArgumentException.h"
#include <arcane/IParallelMng.h>
#include <arcane/ITimeLoopMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ICaseMng.h>
#include <arcane/geometry/IGeometry.h>
#include <arcane/IIOMng.h>
#include <arcane/CaseTable.h>

#include "IDoFLinearSystemFactory.h"
#include "Integer3std.h"
#include "ElastodynamicModule.h"
#include "utilFEM.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ElastodynamicModule::ElastodynamicModule(const ModuleBuildInfo& mbi)
        : ArcaneElastodynamicObject(mbi)
        , m_dofs_on_nodes(mbi.subDomain()->traceMng())
{
    ICaseMng *cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
}

VersionInfo ElastodynamicModule::versionInfo() const {
  return {1, 0, 0};
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
startInit(){

  info() << "Module Elastodynamic INIT";

  m_linear_system.reset();
  m_linear_system.setLinearSystemFactory(options()->linearSystem());

  integ_order.m_i = options()->getNint1();
  integ_order.m_j = options()->getNint2();
  integ_order.m_k = options()->getNint3();
  gravity.x = options()->getGx();
  gravity.y = options()->getGy();
  gravity.z = options()->getGz();

  if (options()->enforceDirichletMethod() == "Penalty") {
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
  if (szType.contains("young")) elast_type = TypesElastodynamic::YoungNu;
  else if (szType.contains("lame")) elast_type = TypesElastodynamic::Lame;
  else if (szType.contains("vel")) elast_type = TypesElastodynamic::Veloc;
  else
    ARCANE_FATAL("Type for elastic properties is undefined!");

  auto nsteps = (int)((tf - t)/dt);
  if (linop_nstep > nsteps) keep_constop = true;

  is_alfa_method = options()->alfa_method();
  if (is_alfa_method) {
    gamma = 0.5 + alfaf - alfam;
    beta = 0.5*pow(0.5 + gamma,2);
  }
  else{
    alfam = 0.;
    alfaf = 0.;
  }

  if (options()->getAnalysisType() == TypesElastodynamic::ThreeD)
    NDIM = 3;
  else
    NDIM = 2;
  if (NDIM == 2) integ_order.m_k = 0;

  cell_fem.set_node_coords(m_node_coord);
  gausspt.init_order(integ_order);

  String dirichletMethod = options()->enforceDirichletMethod();

  if (dirichletMethod != "Penalty" && dirichletMethod != "WeakPenalty"
      && dirichletMethod != "RowElimination"
      && !dirichletMethod.contains("RowColumnElimination")) {
    info() << "Dirichlet boundary condition via "
           << dirichletMethod << " is not supported \n"
           << "enforce-Dirichlet-method only supports:\n"
           << "  - Penalty\n"
           << "  - WeakPenalty\n"
           << "  - RowElimination\n"
           << "  - RowColumnElimination\n";

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
    m_acc[node] = Real3::zero() ;
    m_vel[node] = Real3::zero() ;
    m_displ[node] = Real3::zero();
  }

  _applyInitialNodeConditions();
  _initCells();
  _initBoundaryConditions();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_initDofs(){
  m_dofs_on_nodes.initialize(mesh(),NDIM);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_initCells(){
  Real vp, vs, E, nu, lambda, mu;

  ENUMERATE_CELL (icell, allCells()) {
    const Cell& cell = *icell;
    auto rho = m_rho[cell];
    if (elast_type == TypesElastodynamic::YoungNu) {
      E = m_young[cell];
      nu = m_nu[cell];
      lambda = nu*E/(1. + nu)/(1. - 2.*nu);
      mu = E/2./(1. + nu);
      vp = math::sqrt( (lambda + 2. * mu)/rho );
      vs = math::sqrt( mu/rho );

    } else if (elast_type == TypesElastodynamic::Lame) {
      lambda = m_lambda[cell];
      mu = m_mu[cell];
      vp = math::sqrt( (lambda + 2.*mu)/rho );
      vs = math::sqrt( mu/rho );
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);

    } else if (elast_type == TypesElastodynamic::Veloc) {
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
void ElastodynamicModule::
_applyInitialNodeConditions(){

  for (Int32 i = 0, nb = options()->initialNodeCondition().size(); i < nb; ++i) {

    NodeGroup node_group = options()->initialNodeCondition[i]->nodeGroup();
    TypesElastodynamic::eNodeCondition type = options()->initialNodeCondition[i]->type();
    Real3 values = options()->initialNodeCondition[i]->vector();

    // Loop on nodes with this initial condition
    ENUMERATE_NODE(inode, node_group) {
      const Node & node = *inode;
      switch (type) {

      case TypesElastodynamic::Acc:
        m_prev_acc[node] = values;
        break;
      case TypesElastodynamic::Displ:
        m_prev_displ[node] = values;
        break;
      case TypesElastodynamic::Vel:
        m_prev_vel[node] = values;
      case TypesElastodynamic::Force:
        m_force[node] = values;
      case TypesElastodynamic::UnknownCond:
        break;

      }
    }
  }
/*  if (!options()->inputMotion().empty()) {
    // Loop on nodes with this initial condition
    ENUMERATE_NODE (inode, m_input.m_node_group) {
      const Node& node = *inode;
      Real3 val;
      Real time = options()->getStart();
      m_input.m_acc->value(time, val);
      m_prev_acceleration[node] = val.mul(m_input.m_ampli_factors);

      if (m_input.m_is_vel) {
        m_input.m_vel->value(time, val);
        m_prev_velocity[node] = val.mul(m_input.m_ampli_factors);
      }
      if (m_input.m_is_displ) {
        m_input.m_displ->value(time, val);
        m_prev_displacement[node] = val.mul(m_input.m_ampli_factors);
      }
    }*/
  }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_applyInitialCellConditions(){

  for (Integer i = 0, nb = options()->initElastProperties().size(); i < nb; ++i) {

    CellGroup cell_group = options()->initElastProperties[i]->cellGroup();

    // In the future, we will have to find a way to impose different initial
    // properties (stress/strain tensors, densities...) per element from a file
    // (e.g., coming from a previous computation)
    auto rho = options()->initElastProperties[i]->rho();
    Real vp, vs, E, nu, lambda, mu;

    if (elast_type == TypesElastodynamic::YoungNu) {
      E = options()->initElastProperties[i]->young();
      nu = options()->initElastProperties[i]->nu();
      lambda = nu*E/(1. + nu)/(1. - 2.*nu);
      mu = E/2./(1. + nu);
      vp = math::sqrt( (lambda + 2. * mu)/rho );
      vs = math::sqrt( mu/rho );

    } else if (elast_type == TypesElastodynamic::Lame) {
      lambda = options()->initElastProperties[i]->young();
      mu = options()->initElastProperties[i]->nu();
      vp = math::sqrt( (lambda + 2.*mu)/rho );
      vs = math::sqrt( mu/rho );
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);

    } else if (elast_type == TypesElastodynamic::Veloc) {
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

  for (Integer i = 0, nb = options()->initCellCondition().size(); i < nb; ++i) {

    CellGroup cell_group = options()->initCellCondition[i]->cellGroup();
    TypesElastodynamic::eCellCondition type = options()->initCellCondition[i]->type();

    // In the future, we will have to find a way to impose different initial
    // properties (stress/strain tensors, densities...) per element from a file
    // (e.g., coming from a previous computation)
    auto values = options()->initCellCondition[i]->constVolPart();
    if (type == TypesElastodynamic::Stress) {

        // Loop on cells with this initial condition
        ENUMERATE_CELL (icell, cell_group) {
          const Cell& cell = *icell;

          // Initialize the stress tensor for the concerned cell
          m_stress[cell].x = Real3(values.x, 0., 0.);
          m_stress[cell].y = Real3(0., values.y, 0.);
          m_stress[cell].z = Real3(0., 0., values.z);
        }
    } else if (type == TypesElastodynamic::Strain) {

        ENUMERATE_CELL (icell, cell_group) {
          const Cell& cell = *icell;
          // Initialize the strain tensor for the concerned cell
          m_strain[cell].x = Real3(values.x, 0., 0.);
          m_strain[cell].y = Real3(0., values.y, 0.);
          m_strain[cell].z = Real3(0., 0., values.z);
        }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
compute(){

  info() << "Module PASSMO COMPUTE";
  ++linop_nstep_counter;

  // Stop code at exact final time set by user
  auto tf = m_global_final_time();
  auto t = m_global_time();
  auto dt = m_global_deltat();
  auto t0 = options()->getStart();
  dt2 = dt * dt;

  if (t > tf)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

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

  // Apply other Dirichlet/Neumann conditions if any (constant values assumed at present)
  _applyDirichletBoundaryConditions();
  _applyNeumannBoundaryConditions();
  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();

  // Assemble the FEM global operators (LHS matrix/RHS vector b)
  if (NDIM <= 2) {
    _assembleLinearLHS2D();
    _assembleLinearRHS2D();
  }
  else {
    _assembleLinearLHS3D();
    _assembleLinearRHS3D();
  }

  // Solve the linear system AX = B
  _doSolve();

  // Update the nodal variable according to the integration scheme (e.g. Newmark)
  _updateNewmark();

  // Save/Check results
//  _checkResultFile();

  if (t < tf && t + dt > tf) {
    dt = tf - t;
    m_global_deltat = dt;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_updateNewmark(){

  // Updating the nodal accelerations and velocities (after solve) with
  auto dt = m_global_deltat();

  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    auto an = m_prev_acc[node];
    auto vn = m_prev_vel[node];
    auto dn = m_prev_displ[node];
    auto dn1 = m_displ[node];

    if (!is_alfa_method) {
      for (Int32 i = 0; i < NDIM; ++i) {

          auto ba = (bool)m_imposed_acc[node][i];
          auto bv = (bool)m_imposed_vel[node][i];
          auto ui = dn[i] + dt * vn[i] + dt2 * (0.5 - beta) * an[i];
          auto vi = vn[i] + dt * (1. - gamma) * an[i];

        if (!ba)
          m_acc[node][i] = (dn1[i] - ui)/beta/dt2;

        if (!bv)
          m_vel[node][i] = vi + dt*gamma*m_acc[node][i];
      }
    } else {
      // TO DO
    }

    m_prev_acc[node] = m_acc[node];
    m_prev_vel[node] = m_vel[node];
    m_prev_displ[node] = m_displ[node];
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_initBoundaryConditions()
{
  IParallelMng* pm = subDomain()->parallelMng();

  for (const auto& bd : options()->dirichletBoundaryCondition()) {
    FaceGroup face_group = bd->surface();
    if (bd->hasACurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_sacc_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasUCurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_sdispl_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasVCurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_svel_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasFCurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_sforce_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    // Loop on faces of the surface
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;
      Int32 nb_node = face.nbNode();

      // Loop on nodes of the face
      for (Int32 k = 0; k < nb_node; ++k) {
        const Node& node = face.node(k);
        auto coord = m_node_coord[node];
        auto num = node.uniqueId();
        m_imposed_displ[node].x = bd->hasUx();
        m_imposed_displ[node].y = bd->hasUy();
        m_imposed_displ[node].z = bd->hasUz();
        m_imposed_acc[node].x = bd->hasAx();
        m_imposed_acc[node].y = bd->hasAy();
        m_imposed_acc[node].z = bd->hasAz();
        m_imposed_vel[node].x = bd->hasVx();
        m_imposed_vel[node].y = bd->hasVy();
        m_imposed_vel[node].z = bd->hasVz();
        m_imposed_force[node].x = bd->hasFx();
        m_imposed_force[node].y = bd->hasFy();
        m_imposed_force[node].z = bd->hasFz();
      }
    }
  }

  for (const auto& bd : options()->dirichletPointCondition()) {
    NodeGroup nodes = bd->node();

    if (bd->hasACurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_acc_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasUCurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_displ_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasVCurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_vel_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    if (bd->hasFCurve()) {
      String file_name = bd->ACurve();
      if (!file_name.empty()) {
        auto case_table = readFileAsCaseTable(pm, file_name, 3);
        m_force_case_table_list.add(CaseTableInfo{ file_name, case_table });
      }
    }

    // Loop on nodes
    ENUMERATE_NODE (inode, nodes) {
      const Node& node = *inode;
      auto coord = m_node_coord[node];
      auto num = node.uniqueId();
      m_imposed_displ[node].x = bd->hasUx();
      m_imposed_displ[node].y = bd->hasUy();
      m_imposed_displ[node].z = bd->hasUz();
      m_imposed_acc[node].x = bd->hasAx();
      m_imposed_acc[node].y = bd->hasAy();
      m_imposed_acc[node].z = bd->hasAz();
      m_imposed_vel[node].x = bd->hasVx();
      m_imposed_vel[node].y = bd->hasVy();
      m_imposed_vel[node].z = bd->hasVz();
      m_imposed_force[node].x = bd->hasFx();
      m_imposed_force[node].y = bd->hasFy();
      m_imposed_force[node].z = bd->hasFz();
    }
  }

  for (const auto& bs : options()->neumannBoundaryCondition()) {
    FaceGroup face_group = bs->surface();
    String file_name = bs->getCurve();
    if (!file_name.empty()) {
      auto case_table = readFileAsCaseTable(pm, file_name, 3);
      m_traction_case_table_list.add(CaseTableInfo{ file_name, case_table });
    }
  }

  for (const auto& bs : options()->paraxialBoundaryCondition()) {

    FaceGroup face_group = bs->surface();

    // Loop on the faces (=edges in 2D) concerned with the paraxial condition
    // Initializing the local referential per face (just done once) for further use
    ENUMERATE_FACE (iface, face_group) {

      const Face& face = *iface;

      if (face.isSubDomainBoundary() && face.isOwn()) {

        Real3 e1{ 0. }, e2{ 0. }, e3{ 0. };
        DirVectors(face, m_node_coord, NDIM, e1, e2, e3);
        m_e1_boundary[face] = e1;
        m_e2_boundary[face] = e2;
        m_e3_boundary[face] = e3;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_applyDirichletBoundaryConditions(){

  Int32 sac_index{ 0 }, svc_index{ 0 }, suc_index{ 0 }, sfc_index{ 0 };
  for (const auto& bd : options()->dirichletBoundaryCondition()) {
    FaceGroup face_group = bd->surface();

    Real3 acc{};
    bool is_acc_imp{bd->hasACurve() || bd->hasAx() || bd->hasAy() || bd->hasAz()};
    if (bd->hasACurve()) {
      const CaseTableInfo& table_info = m_sacc_case_table_list[sac_index++];
      String file_name = bd->ACurve();
      info() << "Applying acceleration boundary conditions for surface " << face_group.name()
             << " via CaseTable" << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(m_global_time(), acc);
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
             << " via CaseTable" << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(m_global_time(), vel);
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
             << " via CaseTable" << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(m_global_time(), displ);
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
      String file_name = bd->UCurve();
      info() << "Applying force boundary conditions for surface " << face_group.name()
             << " via CaseTable" << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(m_global_time(), force);
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
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;
      Integer nb_node = face.nbNode();

      // Loop on nodes of the face
      for (Integer k = 0; k < nb_node; ++k) {
        const Node& node = face.node(k);
        auto coord = m_node_coord[node];
        auto num = node.uniqueId();

        if (is_acc_imp)
          m_acc[node] = acc;

        if (is_vel_imp)
          m_vel[node] = vel;

        if (is_displ_imp)
          m_displ[node] = displ;

        if (is_force_imp)
          m_force[node] = force;
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
             << " via CaseTable" << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(m_global_time(), acc);
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
             << " via CaseTable" << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(m_global_time(), vel);
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
             << " via CaseTable" << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(m_global_time(), displ);
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
      String file_name = bd->UCurve();
      info() << "Applying force boundary conditions for nodes " << nodes.name()
             << " via CaseTable" << file_name;
      CaseTable* inn = table_info.case_table;

      if (inn != nullptr)
        inn->value(m_global_time(), force);
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
      auto coord = m_node_coord[node];
      auto num = node.uniqueId();

      if (is_acc_imp)
        m_acc[node] = acc;

      if (is_vel_imp)
        m_vel[node] = vel;

      if (is_displ_imp)
        m_displ[node] = displ;

      if (is_force_imp)
        m_force[node] = force;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_applyNeumannBoundaryConditions(){
  Int32 bc_index{ 0 };
  for (const auto& bs : options()->neumannBoundaryCondition()) {
    FaceGroup face_group = bs->surface();
    const CaseTableInfo& case_table_info = m_traction_case_table_list[bc_index];
    ++bc_index;

    Real3 trac{};

    if (bs->curve.isPresent()) {
      String file_name = bs->getCurve();
      info() << "Applying traction boundary conditions for surface " << face_group.name()
             << " via CaseTable" << file_name;
      CaseTable* inn = case_table_info.case_table;

      if (inn != nullptr)
        inn->value(m_global_time(), trac);
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
  // ***TO DO: we may need to add an incident transient wave field for paraxial
  // boundary conditions (e.g., plane wave, etc.), not only an absorbing condition
  // Not implemented yet...
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Computes the Jacobian Matrix of a 3D finite-element at Gauss Point ig
Real3x3 ElastodynamicModule::
_computeJacobian3D(const ItemWithNodes& cell,const Int32& ig, const RealUniqueArray& vec, Real& jacobian) {

  auto	n = cell.nbNode();

  // Jacobian matrix computed at the integration point
  Real3x3	jac;

  for (Int32 inod = 0, indx = 4; inod < n; ++inod) {

    // vector of local derivatives at this integration point, for node inod
    Real3 dPhi {vec[ig + indx + 1], vec[ig + indx + 2], vec[ig + indx + 3]};
    auto coord_nod = m_node_coord[cell.node(inod)];
    jac.x.x += dPhi.x * coord_nod.x;
    jac.x.y += dPhi.x * coord_nod.y;
    jac.x.z += dPhi.x * coord_nod.z;
    jac.y.x += dPhi.y * coord_nod.x;
    jac.y.y += dPhi.y * coord_nod.y;
    jac.y.z += dPhi.y * coord_nod.z;
    jac.z.x += dPhi.z * coord_nod.x;
    jac.z.y += dPhi.z * coord_nod.y;
    jac.z.z += dPhi.z * coord_nod.z;
    indx += 4;
  }
  jacobian = math::matrixDeterminant(jac);
  if (fabs(jacobian) < REL_PREC) {
    ARCANE_FATAL("Cell jacobian is null");
  }
  return jac;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Computes the Jacobian Matrix of a 2D finite-element at Gauss Point ig
Real2x2 ElastodynamicModule::
_computeJacobian2D(const ItemWithNodes& cell,const Int32& ig, const RealUniqueArray& vec, Real& jacobian) {

  auto	n = cell.nbNode();

  // Jacobian matrix computed at the integration point
  Real2x2	jac;

  for (Int32 inod = 0, indx = 4; inod < n; ++inod) {

    // vector of local derivatives at this integration point, for node inod
    Real2 dPhi {vec[ig + indx + 1], vec[ig + indx + 2]};
    auto coord_nod = m_node_coord[cell.node(inod)];
    jac.x.x += dPhi.x * coord_nod.x;
    jac.x.y += dPhi.x * coord_nod.y;
    jac.y.x += dPhi.y * coord_nod.x;
    jac.y.y += dPhi.y * coord_nod.y;
    indx += 4;
  }
  jacobian = jac.x.x * jac.y.y - jac.x.y * jac.y.x;
  if (fabs(jacobian) < REL_PREC) {
    ARCANE_FATAL("Cell jacobian is null");
  }
  return jac;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute stiffness and mass matrix (only in 2D/3D) and cell forces
/*
void ElastodynamicModule::
_computeKMF3D(const Cell& cell,FixedMatrix<24,24>& Ke, FixedMatrix<24,24>& Me, FixedVector<24>& Fe){

  Integer nb_nodes = cell.nbNode();
//  Integer nk{3*nb_nodes};
  Real rho = m_rho(cell);
  Real nu = m_nu(cell);
  Real young = m_young(cell);
//  Integer nb{6};

  // Setting the "B" matrix size for the max number of nodes a lin element can have in 3D (=8 nodes)
  FixedMatrix<3,8> Bmat;
  ElastTensor D(young,nu,3);

  // Loop on the cell Gauss points to compute integrals terms
  Int32 ngauss{ 0 };
  auto vec = cell_fem.getGaussData(cell, integ_order, ngauss);

  for (Int32 igauss = 0, ig = 0; igauss < ngauss; ++igauss, ig += (4+nb_nodes) ) {

    auto wt = vec[ig];
    Real3 pos { vec[ig+1], vec[ig+2], vec[ig+3]};
    auto ijac = _computeInverseJacobian3D(cell,pos);

    //------------------------------------------------------
    // Elementary Derivation Matrix B at current Gauss point
    //------------------------------------------------------
    for (Int32 inod = 0; inod < nb_nodes; ++inod) {
      auto dPhi = cell_fem.getShapeFuncDeriv(cell.type(), inod, pos);
      Bmat(0, inod) += ijac.x.x * dPhi.x + ijac.x.y * dPhi.y + ijac.x.z * dPhi.z;
      Bmat(1, inod) += ijac.y.x * dPhi.x + ijac.y.y * dPhi.y + ijac.y.z * dPhi.z;
      Bmat(2, inod) += ijac.z.x * dPhi.x + ijac.z.y * dPhi.y + ijac.z.z * dPhi.z;
    }

    for (Int32 inod = 0, iig = 4; inod < nb_nodes; ++inod, ++iig) {

      auto rhoPhi_i = wt*rho*vec[ig + iig];

      //----------------------------------------------
      // Elementary Force vector Fe assembly
      //----------------------------------------------
      if (options()->hasBodyf()) {
        //----------------------------------------------
        // Body force terms
        //----------------------------------------------
        Fe(3 * inod + 0) += rhoPhi_i * gravity.x;
        Fe(3 * inod + 1) += rhoPhi_i * gravity.y;
        Fe(3 * inod + 2) += rhoPhi_i * gravity.z;
      }

      //----------------------------------------------
      // Imposed nodal forces
      //----------------------------------------------
      {
        const Node& nodei = cell.node(inod);
        if ((bool)m_imposed_force[nodei].x)
          Fe(3 * inod + 0) += m_force[nodei].x;
        if ((bool)m_imposed_force[nodei].y)
          Fe(3 * inod + 1) += m_force[nodei].y;
        if ((bool)m_imposed_force[nodei].z)
          Fe(3 * inod + 1) += m_force[nodei].z;
      }

      //----------------------------------------------
      // Elementary Mass (Me) Matrix assembly
      //----------------------------------------------
      for (Int32 jnod = 0, jig = 4; jnod < nb_nodes; ++jnod, ++jig) {

        auto Phij = vec[ig + jig];
        auto mij = rhoPhi_i*Phij;

        for (Int32 l = 0; l < 3; ++l){
          int ii = 3 * inod + l;
          int jj = 3 * jnod + l;
          Me(ii,jj) += mij;

          //----------------------------------------------
          // Elementary Stiffness (Ke) Matrix assembly
          //----------------------------------------------
          FixedVector<6> Bii;

          auto ir = ii%3;

          if (!ir) {
            auto i3{ (Int32)(ii / 3) };
            Bii(0) = Bmat(0, i3);
            Bii(1) = 0.;
            Bii(2) = 0.;
            Bii(3) = Bmat(1, i3);
            Bii(4) = Bmat(2, i3);
            Bii(5) = 0.;

          }
          else if (ir == 1) {
            auto i3{ (Int32)((ii - 1) / 3) };
            Bii(0) = 0.;
            Bii(1) = Bmat(1, i3);
            Bii(2) = 0.;
            Bii(3) = Bmat(0, i3);
            Bii(4) = 0.;
            Bii(5) = Bmat(2, i3);

          }
          else if (ir == 2) {
            auto i3{ (Int32)((ii - 2) / 3) };
            Bii(0) = 0.;
            Bii(1) = 0.;
            Bii(2) = Bmat(2, i3);
            Bii(3) = 0.;
            Bii(4) = Bmat(0, i3);
            Bii(5) = Bmat(1, i3);
          }

          for (Int32 ll = 0; ll < 3; ++ll) {

            jj = 3 * jnod + ll;
            FixedVector<6> Bjj;

            auto jr = jj%3;

            if (!jr) {
              auto j3{ (Int32)(jj / 3) };
              Bjj(0) = Bmat(0, j3);
              Bjj(1) = 0.;
              Bjj(2) = 0.;
              Bjj(3) = Bmat(1, j3);
              Bjj(4) = Bmat(2, j3);
              Bjj(5) = 0.;

            }
            else if (jr == 1) {
              auto j3{ (Int32)((jj - 1) / 3) };
              Bjj(0) = 0.;
              Bjj(1) = Bmat(1, j3);
              Bjj(2) = 0.;
              Bjj(3) = Bmat(0, j3);
              Bjj(4) = 0.;
              Bjj(5) = Bmat(2, j3);

            }
            else if (jr == 2) {
              auto j3{ (Int32)((jj - 2) / 3) };
              Bjj(0) = 0.;
              Bjj(1) = 0.;
              Bjj(2) = Bmat(2, j3);
              Bjj(3) = 0.;
              Bjj(4) = Bmat(0, j3);
              Bjj(5) = Bmat(1, j3);
            }

            Ke(ii, jj) += wt * (  Bii(0) * (D(0, 0) * Bjj(0) + D(0, 1) * Bjj(1) + D(0, 2) * Bjj(2)
                                          + D(0, 3) * Bjj(3) + D(0, 4) * Bjj(4) + D(0, 5) * Bjj(5))
                                + Bii(1) * (D(1, 0) * Bjj(0) + D(1, 1) * Bjj(1) + D(1, 2) * Bjj(2)
                                          + D(1, 3) * Bjj(3) + D(1, 4) * Bjj(4) + D(1, 5) * Bjj(5))
                                + Bii(2) * (D(2, 0) * Bjj(0) + D(2, 1) * Bjj(1) + D(2, 2) * Bjj(2)
                                          + D(2, 3) * Bjj(3) + D(2, 4) * Bjj(4) + D(2, 5) * Bjj(5))
                                + Bii(3) * (D(3, 0) * Bjj(0) + D(3, 1) * Bjj(1) + D(3, 2) * Bjj(2)
                                          + D(3, 3) * Bjj(3) + D(3, 4) * Bjj(4) + D(3, 5) * Bjj(5))
                                + Bii(4) * (D(4, 0) * Bjj(0) + D(4, 1) * Bjj(1) + D(4, 2) * Bjj(2)
                                          + D(4, 3) * Bjj(3) + D(4, 4) * Bjj(4) + D(4, 5) * Bjj(5))
                                + Bii(5) * (D(5, 0) * Bjj(0) + D(5, 1) * Bjj(1) + D(5, 2) * Bjj(2)
                                          + D(5, 3) * Bjj(3) + D(5, 4) * Bjj(4) + D(5, 5) * Bjj(5))
                               );
          }
        }
      }
    }
    ig += 4 + nb_nodes;
  }
}
*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute elementary stiffness matrix in 3D at a given Gauss point
void ElastodynamicModule::
_computeK3D(const Cell& cell,const Int32& ig, const RealUniqueArray& vec, const Real3x3& jac, FixedMatrix<60,60>& Ke){

  Int32 nb_nodes = cell.nbNode();
  auto size{3 * nb_nodes};

  // Setting the "B" matrix size for the max number of nodes in 3D:
  // 8 nodes for a lin element/20 nodes for a quadratic one
  FixedMatrix<3,20> Bmat;

  auto lambda = m_lambda(cell);
  auto mu = m_mu(cell);
  auto a{ lambda + 2.*mu };

  for (int i = 0; i <  3; ++i)
    for (int j = 0; j < 20; ++j) {
      Bmat(i, j) = 0.;
    }

  // ! Computes the Inverse Jacobian Matrix of a 3D finite-element
  auto jacobian = math::matrixDeterminant(jac);
  Real3x3 ijac = math::inverseMatrix(jac);

  auto wt = vec[ig] * jacobian;

  //------------------------------------------------------
  // Elementary Derivation Matrix B at current Gauss point
  //------------------------------------------------------
  for (Int32 inod = 0, iig = 4; inod < nb_nodes; ++inod) {
    Real3 dPhi {vec[ig + iig + 1], vec[ig + iig + 2], vec[ig + iig + 3]};
    auto b1 = ijac.x.x * dPhi.x + ijac.x.y * dPhi.y + ijac.x.z * dPhi.z;
    auto b2 = ijac.y.x * dPhi.x + ijac.y.y * dPhi.y + ijac.y.z * dPhi.z;
    auto b3 = ijac.z.x * dPhi.x + ijac.z.y * dPhi.y + ijac.z.z * dPhi.z;
    Bmat(0, inod) = b1;
    Bmat(1, inod) = b2;
    Bmat(2, inod) = b3;
    iig += 4;
  }

  for (Int32 inod = 0; inod < nb_nodes; ++inod) {

    //----------------------------------------------
    // Elementary Stiffness (Ke) Matrix assembly
    //----------------------------------------------
    for (Int32 l = 0; l < 3; ++l){

      auto ii = 3*inod + l;
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

      for (Int32 jj = ii ; jj < size; ++jj) {

        auto ll = jj%3;
        FixedVector<6> Bjj;

        if (!ll){
          auto jnod{ (Int32)(jj / 3) };
          Bjj(0) = Bmat(0, jnod);
          Bjj(1) = 0.;
          Bjj(2) = 0.;
          Bjj(3) = Bmat(1, jnod);
          Bjj(4) = Bmat(2, jnod);
          Bjj(5) = 0.;
        }
        else if (ll == 1){
          auto jnod{ (Int32)((jj-1) / 3) };
          Bjj(0) = 0.;
          Bjj(1) = Bmat(1, jnod);
          Bjj(2) = 0.;
          Bjj(3) = Bmat(0, jnod);
          Bjj(4) = 0.;
          Bjj(5) = Bmat(2, jnod);
        }
        else if (ll == 2){
          auto jnod{ (Int32)((jj-2) / 3) };
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
                      +  Bii(1) * (lambda * Bjj(0) + a * Bjj(1) + lambda * Bjj(2))
                      +  Bii(2) * (lambda * Bjj(0) + lambda * Bjj(1) + a * Bjj(2))
                      +  Bii(3) * (mu * Bjj(3))
                      +  Bii(4) * (mu * Bjj(4))
                      +  Bii(5) * (mu * Bjj(5)) );

        Ke(ii,jj) = kij;
        Ke(jj,ii) = kij;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute elementary mass matrix in 3D at a given Gauss point
void ElastodynamicModule::
_computeM3D(const Cell& cell,const Int32& ig, const RealUniqueArray& vec, const Real& jacobian, FixedMatrix<60,60>& Me){

  Int32 nb_nodes = cell.nbNode();
  auto rho = m_rho(cell);

  auto wt = vec[ig] * jacobian;
  for (Int32 inod = 0, iig = 4; inod < nb_nodes; ++inod) {

    auto rhoPhi_i = wt*rho*vec[ig + iig];

    //----------------------------------------------
    // Elementary Mass (Me) Matrix assembly
    //----------------------------------------------
    for (Int32 jnod = inod, jig = 4*(1+inod) ; jnod < nb_nodes; ++jnod) {

      auto Phi_j = vec[ig + jig];
      auto mij = rhoPhi_i*Phi_j;

      for (Int32 l = 0; l < 3; ++l){
        auto ii = 3*inod + l;
        auto jj = 3*jnod + l;
        Me(ii,jj) = mij;
        Me(jj,ii) = mij;
      }
      jig += 4;
    }
    iig += 4;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute elementary stiffness matrix in 2D at a given Gauss point
void ElastodynamicModule::
_computeK2D(const Cell& cell,const Int32& ig, const RealUniqueArray& vec, const Real2x2& jac, FixedMatrix<18,18>& Ke){

    Int32 nb_nodes = cell.nbNode();
    auto size{2 * nb_nodes};

    // Setting the "B" matrix size for the max number of nodes in 2D:
    // 4 nodes for a lin element/9 nodes for a quadratic one
    FixedMatrix<2,9> Bmat;

    auto lambda = m_lambda(cell);
    auto mu = m_mu(cell);
    auto a{ lambda + 2.*mu };

    for (int i = 0; i <  2; ++i)
      for (int j = 0; j < 9; ++j) {
        Bmat(i, j) = 0.;
      }

    // ! Computes the Inverse Jacobian Matrix of a 2D finite-element
    auto jacobian = jac.x.x * jac.y.y - jac.x.y * jac.y.x;
    Real2x2 ijac;
    ijac.x.x = jac.y.y / jacobian;
    ijac.x.y = -jac.x.y / jacobian;
    ijac.y.x = -jac.y.x / jacobian;
    ijac.y.y = jac.x.x / jacobian;

    auto wt = vec[ig] * jacobian;

    //------------------------------------------------------
    // Elementary Derivation Matrix B at current Gauss point
    //------------------------------------------------------
    for (Int32 inod = 0, iig = 4; inod < nb_nodes; ++inod) {
      Real2 dPhi {vec[ig + iig + 1], vec[ig + iig + 2]};
      auto b1 = ijac.x.x * dPhi.x + ijac.x.y * dPhi.y;
      auto b2 = ijac.y.x * dPhi.x + ijac.y.y * dPhi.y;
      Bmat(0, inod) = b1;
      Bmat(1, inod) = b2;
      iig += 4;
    }

    for (Int32 inod = 0; inod < nb_nodes; ++inod) {

    //----------------------------------------------
    // Elementary Stiffness (Ke) Matrix assembly
    //----------------------------------------------
      for (Int32 l = 0; l < 2; ++l){

        auto ii = 2*inod + l;
        Real3 Bii{};

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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute elementary mass matrix in 2D at a given Gauss point
void ElastodynamicModule::
_computeM2D(const Cell& cell,const Int32& ig, const RealUniqueArray& vec, const Real& jacobian, FixedMatrix<18,18>& Me){

    Int32 nb_nodes = cell.nbNode();
    auto rho = m_rho(cell);

    auto wt = vec[ig] * jacobian;
    for (Int32 inod = 0, iig = 4; inod < nb_nodes; ++inod) {

    auto rhoPhi_i = wt*rho*vec[ig + iig];

    //----------------------------------------------
    // Elementary Mass (Me) Matrix assembly
    //----------------------------------------------
    for (Int32 jnod = inod, jig = 4*(1+inod) ; jnod < nb_nodes; ++jnod) {

      auto Phi_j = vec[ig + jig];
      auto mij = rhoPhi_i*Phi_j;

      for (Int32 l = 0; l < 2; ++l){
        auto ii = 2*inod + l;
        auto jj = 2*jnod + l;
        Me(ii,jj) = mij;
        Me(jj,ii) = mij;
      }
      jig += 4;
    }
    iig += 4;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_computeMFParax3D(const Face& face, const Int32& ig, const RealUniqueArray& vec, const Real& jacobian,
                FixedMatrix<27,27>& Me, FixedVector<27>& Fe,
                const Real& rhocs, const Real& rhocp){

  auto dt = m_global_deltat();
  auto alfa{ gamma / beta / dt };
  auto alfa1{ beta * dt / gamma };
  Real3 ex{ 1., 0., 0. }, ey{ 0., 1., 0. }, ez{ 0., 0., 1. };

  Real3 e1{ m_e1_boundary[face] }, e2{ m_e2_boundary[face] }, e3{ m_e3_boundary[face] };
  Real3x3 Rot({ math::dot(ex, e1), math::dot(ey, e1), math::dot(ez, e1) },
              { math::dot(ex, e2), math::dot(ey, e2), math::dot(ez, e2) },
              { math::dot(ex, e3), math::dot(ey, e3), math::dot(ez, e3) });
  Real3 RhoC{ rhocs, rhocs, rhocp };
  Real3x3 CRot(rhocs * Rot.x, rhocs * Rot.y, rhocp * Rot.z);

  // Loop on the face/edge Gauss points to compute surface integrals terms on the boundary
  Int32 ngauss{ 0 };
  auto wt = vec[ig] * jacobian;

  Int32 nb_nodes = face.nbNode();
  auto size{3 * nb_nodes};

  // Loop on nodes of the face or edge (with no Dirichlet condition)
  for (Int32 inod = 0, iig = 4; inod < nb_nodes; ++inod) {

    auto wtPhi_i = wt*vec[ig + iig];
    Node node1 = face.node(inod);
    auto an = m_prev_acc[node1];
    auto vn = m_prev_vel[node1];
    auto dn = m_prev_displ[node1];

    auto up1 = dn + (1. - alfa1) * vn + (0.5 - alfa1) * dt * an;
    auto alfa_upp1 = alfa * math::multiply(Rot, up1);

    //----------------------------------------------
    // Elementary contribution to force vector
    //----------------------------------------------
    auto ii{3 * inod};
    auto fi0 = - wtPhi_i * rhocs * alfa_upp1.x;
    auto fi1 = - wtPhi_i * rhocs * alfa_upp1.y;
    auto fi2 = - wtPhi_i * rhocp * alfa_upp1.z;
    Fe(ii) += fi0;
    Fe(ii + 1) += fi1;
    Fe(ii + 2) += fi2;

    //----------------------------------------------
    // Elementary contribution to mass matrix
    //----------------------------------------------
    for (Int32 jnod = inod, jig = 4*(1+inod); jnod < nb_nodes; ++jnod) {

      auto Phi_j = vec[ig + jig];

      for (Int32 l = 0; l < 3; ++l) {
        ii = 3 * inod + l;
        auto jj = 3 * jnod + l;
        auto mij = - alfa * wtPhi_i * Phi_j * CRot[l][l];
        Me(ii, jj) = mij;
        Me(jj, ii) = mij;
      }
      jig += 4;
    }
    iig += 4;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_computeMFParax2D(const Face& face, const Int32& ig, const RealUniqueArray& vec, const Real& jacobian,
                  FixedMatrix<6,6>& Me, FixedVector<6>& Fe,
                  const Real& rhocs, const Real& rhocp){

  auto dt = m_global_deltat();
  auto alfa{ (1. - alfaf)*gamma / beta / dt };
  auto alfa1{ beta * dt / gamma };
  Real3 ex{ 1., 0., 0. }, ey{ 0., 1., 0. };

  Real3 e1{ m_e1_boundary[face] }, e2{ m_e2_boundary[face] };
  Real2x2 Rot({ math::dot(ex, e1), math::dot(ey, e1) },
              { math::dot(ex, e2), math::dot(ey, e2) });

  Real2x2 CRot(rhocs * Rot.x, rhocp * Rot.y);

  // Loop on the face/edge Gauss points to compute surface integrals terms on the boundary
  auto wt = vec[ig] * jacobian;

  Int32 nb_nodes = face.nbNode();
  auto size{2 * nb_nodes};

  // Loop on nodes of the face or edge (with no Dirichlet condition)
  for (Int32 inod = 0, iig = 4; inod < nb_nodes; ++inod) {

    auto wtPhi_i = wt*vec[ig + iig];
    Node node1 = face.node(inod);
    auto an = m_prev_acc[node1];
    auto vn = m_prev_vel[node1];
    auto dn = m_prev_displ[node1];

    auto up1 = dn + (1. - alfa1) * vn + (0.5 - alfa1) * dt * an;
    Real2 alfa_upp1;
    alfa_upp1.x = alfa * (Rot.x.x * up1.x + Rot.x.y * up1.y);
    alfa_upp1.y = alfa * (Rot.y.x * up1.x + Rot.y.y * up1.y);

    //----------------------------------------------
    // Elementary contribution to force vector
    //----------------------------------------------
    auto ii{2 * inod};
    auto fi0 = - wtPhi_i * rhocs * alfa_upp1.x;
    auto fi1 = - wtPhi_i * rhocp * alfa_upp1.y;
    Fe(ii) += fi0;
    Fe(ii + 1) += fi1;
    //----------------------------------------------
    // Elementary contribution to mass matrix
    //----------------------------------------------
    for (Int32 jnod = inod, jig = 4*(1+inod); jnod < nb_nodes; ++jnod) {

      auto Phi_j = vec[ig + jig];

      for (Int32 l = 0; l < 2; ++l) {
        ii = 2 * inod + l;
        auto jj = 2 * jnod + l;
        auto mij = - alfa * wtPhi_i * Phi_j * CRot[l][l];
        Me(ii, jj) = mij;
        Me(jj, ii) = mij;
      }
      jig += 4;
    }
    iig += 4;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Assemble the 3D bilinear operator (Left Hand Side A matrix)
void ElastodynamicModule::
_assembleLinearLHS3D()
{
    VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
    rhs_values.fill(0.0);
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    info() << "Assembly of the 3D FEM bilinear (LHS - matrix A) operator ";

    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      auto nb_nodes{ cell.nbNode() };
      // Setting the elementary matrices sizes for the max number of nodes
      // a lin element can have max 8 nodes and a quadratic one, 20
      auto size{3 * nb_nodes};
      FixedMatrix<60,60> Me;
      FixedMatrix<60,60> Ke;

      for (Int32 i = 0; i < size; ++i) {
      for (Int32 j = i; j < size; ++j) {
        Me(i,j) = 0.;
        Me(j,i) = 0.;
        Ke(i,j) = 0.;
        Ke(j,i) = 0.;
      }
      }

      // Loop on the cell Gauss points to compute integrals terms
      Int32 ngauss{ 0 };
      auto vec = cell_fem.getGaussData(cell, integ_order, ngauss);
      auto cm{(1 - alfam)/beta/dt2};
      auto ck{(1 - alfaf)};

      for (Int32 igauss = 0, ig = 0; igauss < ngauss; ++igauss, ig += 4*(1 + nb_nodes)) {

        auto jacobian {0.};
        auto jac = _computeJacobian3D(cell, ig, vec, jacobian);

        // Computing elementary mass matrix at Gauss point ig
        _computeM3D(cell, ig, vec, jacobian, Me);

        // Computing elementary stiffness matrix at Gauss point ig
        _computeK3D(cell, ig, vec, jac, Ke);

        // Considering a simple Newmark scheme here (Generalized-alfa will be done later)
        // Computing Me/beta/dt^2 + Ke
        Int32 n1_index{ 0 };
        for (Node node1 : cell.nodes()) {

          for (Int32 iddl = 0; iddl < 3; ++iddl) {
            DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
            auto ind1 = node1_dofi.asInt32();
            auto ii = 3 * n1_index + iddl;

            // Assemble global bilinear operator (LHS)
            Int32 n2_index{ 0 };
            for (Node node2 : cell.nodes()) {

              for (Int32 jddl = 0; jddl < 3; ++jddl) {
                auto node2_dofj = node_dof.dofId(node2, jddl);
                auto ind2 = node2_dofj.asInt32();
                auto jj = 3 * n2_index + jddl;
                auto mij = Me(ii, jj);
                auto kij = Ke(ii, jj);
                auto aij = cm * mij + ck * kij;

                if (node1.isOwn())
                  m_linear_system.matrixAddValue(node1_dofi, node2_dofj, aij);
              }
              ++n2_index;
            }
          }
          ++n1_index;
        }
      }
    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Assemble the 3D linear operator (Right Hand Side B vector)
void ElastodynamicModule::
_assembleLinearRHS3D(){
    info() << "Assembly of the FEM 3D linear operator (RHS - vector B) ";

    VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
    rhs_values.fill(0.0);
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    auto dt = m_global_deltat();

    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      auto rho = m_rho(cell);
      auto nb_nodes{ cell.nbNode() };

      // Setting the elementary matrices + force vector sizes for the max number of nodes
      // a lin element can have in 2D (=4 nodes)
      auto size{3 * nb_nodes};
      FixedMatrix<60,60> Me;

      for (Int32 i = 0; i < size; ++i) {
        for (Int32 j = i; j < size; ++j) {
          Me(i,j) = 0.;
          Me(j,i) = 0.;
        }
      }

      // Loop on the cell Gauss points to compute integrals terms
      Int32 ngauss{ 0 };
      auto vec = cell_fem.getGaussData(cell, integ_order, ngauss);
      auto cm = (1 - alfam)/beta/dt2;

      for (Int32 igauss = 0, ig = 0; igauss < ngauss; ++igauss, ig += 4*(1 + nb_nodes)) {
        auto jacobian {0.};
        auto jac = _computeJacobian3D(cell, ig, vec, jacobian);

        // Computing elementary mass matrix at Gauss point ig
        _computeM3D(cell, ig, vec, jacobian, Me);

        // Considering a simple Newmark scheme here (Generalized-alfa will be done later)
        // Computing Me/beta/dt^2 + Ke
        Int32 n1_index{ 0 };
        auto iig{4};
        auto wt = vec[ig] * jacobian;

        for (Node node1 : cell.nodes()) {
          for (Int32 iddl = 0; iddl < 3; ++iddl) {
            DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
            auto ind1 = node1_dofi.asInt32();
            auto ii = 3 * n1_index + iddl;

            bool is_node1_dofi_set = (bool)m_imposed_displ[node1][iddl];
            auto rhs_i{ 0. };

            if (node1.isOwn() && !is_node1_dofi_set) {

              /*----------------------------------------------------------
              // Mass contribution to the RHS:
              // (1 - alfm)*Mij/(beta*dt2)*uj_pred - alfm*aj_n
              //----------------------------------------------------------*/
              Int32 n2_index{ 0 };
              for (Node node2 : cell.nodes()) {
                auto an = m_prev_acc[node2][iddl];
                auto vn = m_prev_vel[node2][iddl];
                auto dn = m_prev_displ[node2][iddl];
                auto u_iddl_pred = dn + dt * vn + dt2 * (0.5 - beta) * an;
                auto jj = 3 * n2_index + iddl;
                auto mij = Me(ii, jj);
                rhs_i += mij * (cm * u_iddl_pred - alfam * an);
                ++n2_index;
              }

              /*-------------------------------------------------
              //! Other forces (imposed nodal forces, body forces)
              --------------------------------------------------*/
              {
                if (options()->hasBodyf()) {
                  //----------------------------------------------
                  // Body force terms
                  //----------------------------------------------
                  auto rhoPhi_i = wt*rho*vec[ig + iig];
                  rhs_i += rhoPhi_i * gravity[iddl];
                }

                //----------------------------------------------
                // Imposed nodal forces
                //----------------------------------------------
                if ((bool)m_imposed_force[node1][iddl])
                  rhs_i += m_force[node1][iddl];
              }
              rhs_values[node1_dofi] += rhs_i;
            }
          }
          ++n1_index;
        }
      }
    }

    String dirichletMethod = options()->enforceDirichletMethod();
    info() << "Applying Dirichlet boundary condition via "
           << dirichletMethod << " method ";

    // Looking for Dirichlet boundary nodes & modify linear operators accordingly
    ENUMERATE_ (Node, inode, ownNodes()) {
      auto node = *inode;
      auto num = node.uniqueId().asInt32();

      for (Int32 iddl = 0; iddl < 3; ++iddl) {
        bool is_node_dof_set = (bool)m_imposed_displ[node][iddl];

        if (is_node_dof_set) {
          /*----------------------------------------------------------
              // if Dirichlet node, modify operators (LHS+RHS) allowing to
              // Dirichlet method selected by user
              //----------------------------------------------------------*/
          auto node_dofi = node_dof.dofId(node, iddl);
          auto u_iddl = m_displ[node][iddl];
          if (dirichletMethod == "Penalty") {
            m_linear_system.matrixSetValue(node_dofi, node_dofi, penalty);
            rhs_values[node_dofi] = u_iddl * penalty;
          }
          else if (dirichletMethod == "WeakPenalty") {
            m_linear_system.matrixAddValue(node_dofi, node_dofi, penalty);
            rhs_values[node_dofi] = u_iddl * penalty;
          }
          else if (dirichletMethod == "RowElimination") {
            m_linear_system.eliminateRow(node_dofi, u_iddl);
          }
          else if (dirichletMethod == "RowColumnElimination") {
            m_linear_system.eliminateRowColumn(node_dofi, u_iddl);
          }
        }
      }
    }

    //----------------------------------------------
    // Traction terms assembly
    //----------------------------------------------
    _getTractionContribution(rhs_values);

    //----------------------------------------------
    // Paraxial terms assembly
    //----------------------------------------------
    _getParaxialContribution3D(rhs_values);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Assemble the 2D bilinear operator (Left Hand Side A matrix)
void ElastodynamicModule::
_assembleLinearLHS2D()
{
    VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
    rhs_values.fill(0.0);
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    info() << "Assembly of the FEM 2D bilinear (LHS - matrix A) operator ";

    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      auto nb_nodes{ cell.nbNode() };

      // Setting the 2D elementary matrices sizes for the max number of nodes * 2 dofs per node
      // a lin element can have max 4 nodes and a quadratic element, 9 nodes
      auto size{2*nb_nodes};
      FixedMatrix<18,18> Me;
      FixedMatrix<18,18> Ke;

      for (Int32 i = 0; i < size; ++i) {
        for (Int32 j = i; j < size; ++j) {
          Me(i,j) = 0.;
          Me(j,i) = 0.;
          Ke(i,j) = 0.;
          Ke(j,i) = 0.;
        }
      }

      // Loop on the cell Gauss points to compute integrals terms
      Int32 ngauss{ 0 };
      auto vec = cell_fem.getGaussData(cell, integ_order, ngauss);
      auto cm{(1 - alfam)/beta/dt2};
      auto ck{(1 - alfaf)};

      for (Int32 igauss = 0, ig = 0; igauss < ngauss; ++igauss, ig += 4*(1 + nb_nodes)) {

        auto jacobian {0.};
        auto jac = _computeJacobian2D(cell, ig, vec, jacobian);

        // Computing elementary mass matrix at Gauss point ig
        _computeM2D(cell, ig, vec, jacobian, Me);

        // Computing elementary stiffness matrix at Gauss point ig
        _computeK2D(cell, ig, vec, jac, Ke);

        // Considering a simple Newmark scheme here (Generalized-alfa will be done later)
        // Computing Me/beta/dt^2 + Ke
        Int32 n1_index{ 0 };
        for (Node node1 : cell.nodes()) {

          for (Int32 iddl = 0; iddl < 2; ++iddl) {
            DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
            auto ii = 2 * n1_index + iddl;

            // Assemble global bilinear operator (LHS)
            Int32 n2_index{ 0 };
            for (Node node2 : cell.nodes()) {

              for (Int32 jddl = 0; jddl < 2; ++jddl) {
                auto node2_dofj = node_dof.dofId(node2, jddl);
                auto jj = 2 * n2_index + jddl;
                auto mij = Me(ii, jj);
                auto kij = Ke(ii, jj);
                auto aij = cm * mij + ck * kij;

                if (node1.isOwn())
                  m_linear_system.matrixAddValue(node1_dofi, node2_dofj, aij);
              }
              ++n2_index;
            }
          }
          ++n1_index;
        }
      }
    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Assemble the 2D linear operator (Right Hand Side B vector)
void ElastodynamicModule::
_assembleLinearRHS2D(){
    info() << "Assembly of the FEM 2D linear operator (RHS - vector B) ";

    VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
    rhs_values.fill(0.0);
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    auto dt = m_global_deltat();

    ENUMERATE_ (Cell, icell, allCells()) {
      Cell cell = *icell;
      auto rho = m_rho(cell);
      auto nb_nodes{ cell.nbNode() };

      auto size{2*nb_nodes};
      FixedMatrix<18,18> Me;

      for (Int32 i = 0; i < size; ++i) {
        for (Int32 j = i; j < size; ++j) {
          Me(i,j) = 0.;
          Me(j,i) = 0.;
        }
      }

      // Loop on the cell Gauss points to compute integrals terms
      Int32 ngauss{ 0 };
      auto vec = cell_fem.getGaussData(cell, integ_order, ngauss);
      auto cm = (1 - alfam)/beta/dt2;

      for (Int32 igauss = 0, ig = 0; igauss < ngauss; ++igauss, ig += 4*(1 + nb_nodes)) {
        auto jacobian {0.};
        auto jac = _computeJacobian2D(cell, ig, vec, jacobian);

        // Computing elementary mass matrix at Gauss point ig
        _computeM2D(cell, ig, vec, jacobian, Me);

        // Considering a simple Newmark scheme here (Generalized-alfa will be done later)
        // Computing Me/beta/dt^2 + Ke
        Int32 n1_index{ 0 };
        auto iig{4};
        auto wt = vec[ig] * jacobian;

        for (Node node1 : cell.nodes()) {
          for (Int32 iddl = 0; iddl < 2; ++iddl) {
            DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
            auto ii = 2 * n1_index + iddl;

            bool is_node1_dofi_set = (bool)m_imposed_displ[node1][iddl];
            auto rhs_i{ 0. };

            if (node1.isOwn() && !is_node1_dofi_set) {

              /*----------------------------------------------------------
              // Mass contribution to the RHS:
              // (1 - alfm)*Mij/(beta*dt2)*uj_pred - alfm*aj_n
              //----------------------------------------------------------*/
              Int32 n2_index{ 0 };
              for (Node node2 : cell.nodes()) {
                auto an = m_prev_acc[node2][iddl];
                auto vn = m_prev_vel[node2][iddl];
                auto dn = m_prev_displ[node2][iddl];
                auto u_iddl_pred = dn + dt * vn + dt2 * (0.5 - beta) * an;
                auto jj = 2 * n2_index + iddl;
                auto mij = Me(ii, jj);
                rhs_i += mij * (cm * u_iddl_pred - alfam * an);
                ++n2_index;
              }

              /*-------------------------------------------------
              // Other forces (imposed nodal forces, body forces)
              //-------------------------------------------------*/
              {
                if (options()->hasBodyf()) {
                  //----------------------------------------------
                  // Body force terms
                  //----------------------------------------------
                  auto rhoPhi_i = wt*rho*vec[ig + iig];
                  rhs_i += rhoPhi_i * gravity[iddl];
                }

                //----------------------------------------------
                // Imposed nodal forces
                //----------------------------------------------
                if ((bool)m_imposed_force[node1][iddl])
                  rhs_i += m_force[node1][iddl];
              }
              rhs_values[node1_dofi] += rhs_i;
            }
          }
          ++n1_index;
        }
      }
    }

    String dirichletMethod = options()->enforceDirichletMethod();
    info() << "Applying Dirichlet boundary condition via "
           << dirichletMethod << " method ";

    // Looking for Dirichlet boundary nodes & modify linear operators accordingly
    ENUMERATE_ (Node, inode, ownNodes()) {
      auto node = *inode;

      for (Int32 iddl = 0; iddl < 2; ++iddl) {
        bool is_node_dof_set = (bool)m_imposed_displ[node][iddl];

        if (is_node_dof_set) {
          /*----------------------------------------------------------
            // if Dirichlet node, modify operators (LHS+RHS) allowing to
            // Dirichlet method selected by user
            //----------------------------------------------------------*/
          auto node_dofi = node_dof.dofId(node, iddl);
          auto u_iddl = m_displ[node][iddl];
          if (dirichletMethod == "Penalty") {
            m_linear_system.matrixSetValue(node_dofi, node_dofi, penalty);
            rhs_values[node_dofi] = u_iddl * penalty;
          }
          else if (dirichletMethod == "WeakPenalty") {
            m_linear_system.matrixAddValue(node_dofi, node_dofi, penalty);
            rhs_values[node_dofi] = u_iddl * penalty;
          }
          else if (dirichletMethod == "RowElimination") {
            m_linear_system.eliminateRow(node_dofi, u_iddl);
          }
          else if (dirichletMethod == "RowColumnElimination") {
            m_linear_system.eliminateRowColumn(node_dofi, u_iddl);
          }
        }
      }
    }

    //----------------------------------------------
    // Traction terms assembly
    //----------------------------------------------
    _getTractionContribution(rhs_values);

    //----------------------------------------------
    // Paraxial terms assembly
    //----------------------------------------------
    _getParaxialContribution2D(rhs_values);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_getParaxialContribution3D(Arcane::VariableDoFReal& rhs_values){

    auto dt = m_global_deltat();
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    auto ndim{3};

    for (const auto& bs : options()->paraxialBoundaryCondition()) {

      FaceGroup face_group = bs->surface();
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

        //      ARCANE_FATAL("Paraxial boundary has no elastic properties ");
      }

      info() << "Applying constant paraxial boundary conditions for surface " << face_group.name();

      // Loop on the faces (=edges in 2D) concerned with the paraxial condition
      ENUMERATE_FACE (iface, face_group) {

        const Face& face = *iface;

        if (face.isSubDomainBoundary() && face.isOwn()) {

          if (is_inner) {
            const Cell& cell = face.boundaryCell();
            rho = m_rho[cell];
            cs = m_vs[cell];
            cp = m_vp[cell];
          }

          auto rhocs{ rho * cs };
          auto rhocp{ rho * cp };

          // In 3D, a quadratic face element has max 9 nodes (27 dofs)
          auto nb_nodes{face.nbNode()};
          auto size{ ndim * nb_nodes};
          FixedMatrix<27, 27> Me;
          FixedVector<27> Fe;

          for (Int32 i = 0; i < size; ++i) {
            Fe(i) = 0.;
            for (Int32 j = i; j < size; ++j) {
              Me(i, j) = 0.;
              Me(j, i) = 0.;
            }
          }

          // Loop on the cell Gauss points to compute integrals terms
          Int32 ngauss{ 0 };
          auto vec = cell_fem.getGaussData(face, integ_order, ngauss);

          for (Int32 igauss = 0, ig = 0; igauss < ngauss; ++igauss, ig += 4 * (1 + nb_nodes)) {

            auto jacobian{ 0. };
            auto jac = _computeJacobian3D(face, ig, vec, jacobian);

            _computeMFParax3D(face, ig, vec, jacobian, Me, Fe, rhocs, rhocp);

            // Loop on nodes of the face (with no Dirichlet condition)
            Int32 n1_index{ 0 };
            auto iig{ 4 };
            for (Node node1 : face.nodes()) {

              for (Int32 iddl = 0; iddl < ndim; ++iddl) {

                DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
                auto ii = ndim * n1_index + iddl;

                bool is_node1_dofi_set = (bool)m_imposed_displ[node1][iddl];
                auto rhs_i{ 0. };

                if (node1.isOwn() && !is_node1_dofi_set) {
                  rhs_i += Fe(ii);
                  rhs_values[node1_dofi] += rhs_i;

                  //----------------------------------------------
                  // Elementary contribution to LHS
                  //----------------------------------------------
                  Int32 n2_index{ 0 };
                  for (Node node2 : face.nodes()) {
                    for (Int32 jddl = 0; jddl < ndim; ++jddl) {
                      auto node2_dofj = node_dof.dofId(node2, jddl);
                      auto jj = ndim * n2_index + iddl;
                      auto mij = Me(ii, jj) / beta / dt2;
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
void ElastodynamicModule::
_getParaxialContribution2D(Arcane::VariableDoFReal& rhs_values){

    auto dt = m_global_deltat();
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    for (const auto& bs : options()->paraxialBoundaryCondition()) {

      FaceGroup face_group = bs->surface();
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

        //      ARCANE_FATAL("Paraxial boundary has no elastic properties ");
      }

      info() << "Applying constant paraxial boundary conditions for surface " << face_group.name();

      // Loop on the faces (=edges in 2D) concerned with the paraxial condition
      ENUMERATE_FACE (iface, face_group) {

        const Face& face = *iface;

        if (face.isSubDomainBoundary() && face.isOwn()) {

          if (is_inner) {
            const Cell& cell = face.boundaryCell();
            rho = m_rho[cell];
            cs = m_vs[cell];
            cp = m_vp[cell];
          }

          auto rhocs{ rho * cs };
          auto rhocp{ rho * cp };

          // In 2D, a quadratic edge element has max 3 nodes (6 dofs)
          auto nb_nodes{face.nbNode()};
          auto size{ 2 * nb_nodes};
          FixedMatrix<6, 6> Me;
          FixedVector<6> Fe;

          for (Int32 i = 0; i < size; ++i) {
            Fe(i) = 0.;
            for (Int32 j = i; j < size; ++j) {
              Me(i, j) = 0.;
              Me(j, i) = 0.;
            }
          }

          // Loop on the cell Gauss points to compute integrals terms
          Int32 ngauss{ 0 };
          auto vec = cell_fem.getGaussData(face, integ_order, ngauss);

          for (Int32 igauss = 0, ig = 0; igauss < ngauss; ++igauss, ig += 4 * (1 + nb_nodes)) {

            auto jacobian{ 0. };
            auto jac = _computeJacobian2D(face, ig, vec, jacobian);

            _computeMFParax2D(face, ig, vec, jacobian, Me, Fe, rhocs, rhocp);

            // Loop on nodes of the face (with no Dirichlet condition)
            Int32 n1_index{ 0 };
            auto iig{ 4 };
            for (Node node1 : face.nodes()) {

              for (Int32 iddl = 0; iddl < 2; ++iddl) {

                DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
                auto ii = 2 * n1_index + iddl;

                bool is_node1_dofi_set = (bool)m_imposed_displ[node1][iddl];
                auto rhs_i{ 0. };

                if (node1.isOwn() && !is_node1_dofi_set) {
                  rhs_i += Fe(ii);
                  rhs_values[node1_dofi] += rhs_i;

                  //----------------------------------------------
                  // Elementary contribution to LHS
                  //----------------------------------------------
                  Int32 n2_index{ 0 };
                  for (Node node2 : face.nodes()) {
                    for (Int32 jddl = 0; jddl < 2; ++jddl) {
                      auto node2_dofj = node_dof.dofId(node2, jddl);
                      auto jj = 2 * n2_index + iddl;
                      auto mij = Me(ii, jj) / beta / dt2;
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
void ElastodynamicModule::
_getTractionContribution(Arcane::VariableDoFReal& rhs_values){

    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

    for (const auto& bs : options()->neumannBoundaryCondition()) {
      FaceGroup face_group = bs->surface();

      // Loop on the faces (=edges in 2D) concerned with the traction condition
      ENUMERATE_FACE (j, face_group) {
        const Face& face = *j;

        Real3 trac = m_imposed_traction[face];
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
Real ElastodynamicModule::
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
void ElastodynamicModule::
_doSolve(){
  info() << "Solving Linear system";
  m_linear_system.solve();

  // Re-Apply Dirichlet boundary conditions because the solver has modified the values
  // on all nodes
  _applyDirichletBoundaryConditions(); // ************ CHECK

  {
    VariableDoFReal& dof_d(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;

      auto ux = dof_d[node_dof.dofId(node, 0)];
      auto uy = dof_d[node_dof.dofId(node, 1)];
      auto uz{0.};

      if (NDIM == 3)
          uz = dof_d[node_dof.dofId(node, 2)];

      m_displ[node] = Real3(ux,uy,uz);

      info() << "Node: " << node.localId() << " Ux=" << ux << " Uy=" << uy << " Uz=" << uz;
    }
  }

  m_displ.synchronize();
  m_vel.synchronize();
  m_acc.synchronize();
  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    int p = std::cout.precision();
    std::cout.precision(17);
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      std::cout << "U[" << node.localId() << "][" << node.uniqueId() << "] = "
                << m_displ[node].x << "  " << m_displ[node].y
                << "  " << m_displ[node].z << "\n";
    }
    std::cout.precision(p);
  }

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ARCANE_REGISTER_MODULE_ELASTODYNAMIC(ElastodynamicModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

