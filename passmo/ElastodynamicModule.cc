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

  _initDofs();

  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");

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
  m_global_deltat = options()->getDeltat();
  dt2 = pow(m_global_deltat(),2.);
  m_global_final_time = options()->getFinalTime();
  m_global_time = options()->getStart();
  linop_nstep = options()->getLinopNstep();
  auto nsteps = (int)((m_global_final_time() - m_global_time())/m_global_deltat());
  if (linop_nstep > nsteps) keep_constop = true;

  is_alfa_method = options()->alfa_method();
  if (is_alfa_method) {
    gamma = 0.5 + alfaf - alfam;
    beta = 0.5*pow(0.5 + gamma,2);
  }

  if (options()->getAnalysisType() == TypesElastodynamic::ThreeD)
    NDIM = 3;
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
  //  _initInputMotion();
  _applyInitialNodeConditions();
  _applyInitialCellConditions();
  // This will be useful for nonlinear dynamics only (not used in elastodynamics)
//  _applyInitialCellConditions();

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
    auto type = options()->initElastProperties[i]->type().lower();

    // In the future, we will have to find a way to impose different initial
    // properties (stress/strain tensors, densities...) per element from a file
    // (e.g., coming from a previous computation)
    auto rho = options()->initElastProperties[i]->rho();
    Real vp, vs, E, nu, lambda, mu;

    if (type.contains("young")) {
      E = options()->initElastProperties[i]->young();
      nu = options()->initElastProperties[i]->nu();
      lambda = nu*E/(1. + nu)/(1. - 2.*nu);
      mu = E/2./(1. + nu);
      vp = math::sqrt( (lambda + 2. * mu)/rho );
      vs = math::sqrt( mu/rho );

    } else if (type.contains("lame")) {
      lambda = options()->initElastProperties[i]->young();
      mu = options()->initElastProperties[i]->nu();
      vp = math::sqrt( (lambda + 2.*mu)/rho );
      vs = math::sqrt( mu/rho );
      auto x = lambda/mu;
      nu = x/2./(1. + x);
      E = 2.*mu*(1. + nu);

    } else if (type.contains("vel")) {
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
    bool not_yet_finish = ( m_global_time() < m_global_final_time());
    bool too_much = ( (m_global_time()+m_global_deltat()) > m_global_final_time());

      if (not_yet_finish && too_much ){
        m_global_deltat = m_global_final_time() - m_global_time();
        dt2 = pow(m_global_deltat(),2.);
        subDomain()->timeLoopMng()->stopComputeLoop(true);
      }

      info() << "Time (s) = " << m_global_time()+m_global_deltat();

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

    // Apply the input motion: imposed time history on groups of nodes/faces (acceleration, velocity or displacement)
//    _applyInputMotion();

    // Apply other Dirichlet/Neumann conditions if any (constant values assumed at present)
    _applyBoundaryConditions();

    // Predict the nodal variable according to the integration scheme (e.g. Newmark)
    _predictNewmark();

    info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();

    // Assemble the FEM global operators (LHS matrix/RHS vector b)
    if (NDIM <= 2)
      _assembleLinearGlobal2D();
    else
      _assembleLinearGlobal3D();

    // Solve the linear system AX = B
    _doSolve();

    // Update the nodal variable according to the integration scheme (e.g. Newmark)
    _updateNewmark();

    // Save/Check results
//    _checkResultFile();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_predictNewmark(){

  auto dt = m_global_deltat();

  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    auto an = m_prev_acc[node];
    auto vn = m_prev_vel[node];
    auto dn = m_prev_displ[node];

    for (Int32 i = 0; i < 3; ++i) {

      if (!(bool)m_imposed_displ[node][i])
        m_displ[node][i] = dn[i] + dt * vn[i] + dt2 * (0.5 - beta) * an[i];

      if (!(bool)m_imposed_vel[node][i])
        m_vel[node][i] = vn[i] + dt * (1. - gamma) * an[i];
    }
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

    if (!is_alfa_method) {
      for (Int32 i = 0; i < 3; ++i) {

        if (!(bool)m_imposed_acc[node][i])
          m_acc[node][i] = (m_displ[node][i] - (dn[i] + dt*vn[i] + dt2*(0.5 - beta)*an[i]))/beta/dt2;

        if (!(bool)m_imposed_vel[node][i])
          m_vel[node][i] = vn[i] + dt*( (1. - gamma)*an[i] + gamma*m_acc[node][i] );
      }
    } else {
      // TO DO
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_initBoundaryConditions()
{
  for (const auto& bd : options()->dirichletBoundaryCondition()) {
    FaceGroup face_group = bd->surface();
    NodeGroup node_group = bd->nodeGroup();
    Real value = bd->getConstVal();
    TypesElastodynamic::eBoundaryCondition type = bd->type();

    // Loop on faces of the surface
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;
      Int32 nb_node = face.nbNode();

      // Loop on nodes of the face
      for (Int32 k = 0; k < nb_node; ++k) {
        const Node& node = face.node(k);

        switch (type) {

        case TypesElastodynamic::AX:
          m_imposed_acc[node].x = 1;
          break;
        case TypesElastodynamic::AY:
          m_imposed_acc[node].y = 1;
          break;
        case TypesElastodynamic::AZ:
          m_imposed_acc[node].z = 1;
          break;
        case TypesElastodynamic::UX:
          m_imposed_displ[node].x = 1;
          break;
        case TypesElastodynamic::UY:
          m_imposed_displ[node].y = 1;
          break;
        case TypesElastodynamic::UZ:
          m_imposed_displ[node].z = 1;
          break;
        case TypesElastodynamic::VX:
          m_imposed_vel[node].x = 1;
          break;
        case TypesElastodynamic::VY:
          m_imposed_vel[node].y = 1;
          break;
        case TypesElastodynamic::VZ:
          m_imposed_vel[node].z = 1;
          break;
        case TypesElastodynamic::FX:
          m_imposed_force[node].x = 1;
          break;
        case TypesElastodynamic::FY:
          m_imposed_force[node].y = 1;
          break;
        case TypesElastodynamic::FZ:
          m_imposed_force[node].z = 1;
          break;
        case TypesElastodynamic::Unknown:
          break;
        }
      }
    }

    // Loop on nodes
    ENUMERATE_NODE (inode, node_group) {
      const Node& node = *inode;

      switch (type) {

      case TypesElastodynamic::AX:
        m_imposed_acc[node].x = 1;
        break;
      case TypesElastodynamic::AY:
        m_imposed_acc[node].y = 1;
        break;
      case TypesElastodynamic::AZ:
        m_imposed_acc[node].z = 1;
        break;
      case TypesElastodynamic::UX:
        m_imposed_displ[node].x = 1;
        break;
      case TypesElastodynamic::UY:
        m_imposed_displ[node].y = 1;
        break;
      case TypesElastodynamic::UZ:
        m_imposed_displ[node].z = 1;
        break;
      case TypesElastodynamic::VX:
        m_imposed_vel[node].x = 1;
        break;
      case TypesElastodynamic::VY:
        m_imposed_vel[node].y = 1;
        break;
      case TypesElastodynamic::VZ:
        m_imposed_vel[node].z = 1;
        break;
      case TypesElastodynamic::FX:
        m_imposed_force[node].x = 1;
        break;
      case TypesElastodynamic::FY:
        m_imposed_force[node].y = 1;
        break;
      case TypesElastodynamic::FZ:
        m_imposed_force[node].z = 1;
        break;
      case TypesElastodynamic::Unknown:
        break;
      }
    }
  }

  IParallelMng* pm = subDomain()->parallelMng();

  for (const auto& bs : options()->neumannBoundaryCondition()) {
    FaceGroup face_group = bs->surface();
    String file_name = bs->getCurve();
    if (!file_name.empty()) {
      auto case_table = readFileAsCaseTable(pm, file_name, 3);
      m_traction_case_table_list.add(CaseTableInfo{ file_name, case_table });
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_applyBoundaryConditions(){

  for (const auto& bd : options()->dirichletBoundaryCondition()) {
    FaceGroup face_group = bd->surface();
    NodeGroup node_group = bd->nodeGroup();
    Real value = bd->getConstVal();
    TypesElastodynamic::eBoundaryCondition type = bd->type();

    // Loop on faces of the surface
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;
      Integer nb_node = face.nbNode();

      // Loop on nodes of the face
      for (Integer k = 0; k < nb_node; ++k) {
        const Node& node = face.node(k);

        switch (type) {
        case TypesElastodynamic::AX:
          m_acc[node].x = value;
          break;
        case TypesElastodynamic::AY:
          m_acc[node].y = value;
          break;
        case TypesElastodynamic::AZ:
          m_acc[node].z = value;
          break;
        case TypesElastodynamic::UX:
          m_displ[node].x = value;
          break;
        case TypesElastodynamic::UY:
          m_displ[node].y = value;
          break;
        case TypesElastodynamic::UZ:
          m_displ[node].z = value;
          break;
        case TypesElastodynamic::VX:
          m_vel[node].x = value;
          break;
        case TypesElastodynamic::VY:
          m_vel[node].y = value;
          break;
        case TypesElastodynamic::VZ:
          m_vel[node].z = value;
          break;
        case TypesElastodynamic::FX:
          m_force[node].x = value;
          break;
        case TypesElastodynamic::FY:
          m_force[node].y = value;
          break;
        case TypesElastodynamic::FZ:
          m_force[node].z = value;
          break;
        case TypesElastodynamic::Unknown:
          break;
        }
      }
    }

    // Loop on nodes
    ENUMERATE_NODE (inode, node_group) {
      const Node& node = *inode;

      switch (type) {
      case TypesElastodynamic::AX:
        m_acc[node].x = value;
        break;
      case TypesElastodynamic::AY:
        m_acc[node].y = value;
        break;
      case TypesElastodynamic::AZ:
        m_acc[node].z = value;
        break;
      case TypesElastodynamic::UX:
        m_displ[node].x = value;
        break;
      case TypesElastodynamic::UY:
        m_displ[node].y = value;
        break;
      case TypesElastodynamic::UZ:
        m_displ[node].z = value;
        break;
      case TypesElastodynamic::VX:
        m_vel[node].x = value;
        break;
      case TypesElastodynamic::VY:
        m_vel[node].y = value;
        break;
      case TypesElastodynamic::VZ:
        m_vel[node].z = value;
        break;
      case TypesElastodynamic::FX:
        m_force[node].x = value;
        break;
      case TypesElastodynamic::FY:
        m_force[node].y = value;
        break;
      case TypesElastodynamic::FZ:
        m_force[node].z = value;
        break;
      case TypesElastodynamic::Unknown:
        break;
      }
    }
  }

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
// ! Computes the Inverse Jacobian Matrix of a 3D finite-element
Real3x3 ElastodynamicModule::
_computeInverseJacobian3D(const Cell& cell,const Real3& ref_coord, Real& jacobian) {

  Real3x3 jacmat;
  auto	n = cell.nbNode();

  // Jacobian matrix computed at the integration point
  Real3x3	jac;
  auto cell_type = cell.type();

  for (Int32 inod = 0; inod < n; ++inod) {

    // vector of local derivatives at this integration point, for node inod
    auto dNi = cell_fem.getShapeFuncDeriv(cell_type, inod, ref_coord);
    auto coord_nod = m_node_coord[cell.node(inod)];

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        jac[i][j] += dNi[j] * coord_nod[i];
      }
    }
  }
  jacobian = math::matrixDeterminant(jacmat);
  if (fabs(jacobian) < REL_PREC) {
    ARCANE_FATAL("Cell jacobian is null");
  }

  return math::inverseMatrix(jacmat);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Computes the Inverse Jacobian Matrix of a 2D finite-element
Real2x2 ElastodynamicModule::
_computeInverseJacobian2D(const Cell& cell,const Real3& ref_coord, Real& jacobian) {

  Real2x2 jacmat;
  auto	n = cell.nbNode();

  // Jacobian matrix computed at the integration point
  Real2x2	jac;
  auto cell_type = cell.type();

  for (Int32 inod = 0; inod < n; ++inod) {

    // vector of local derivatives at this integration point, for node inod
    auto dNi = cell_fem.getShapeFuncDeriv(cell_type, inod, ref_coord);
    auto coord_nod = m_node_coord[cell.node(inod)];

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        jac[i][j] += dNi[j] * coord_nod[i];
      }
    }
  }
  jacobian = jacmat.x.x * jacmat.y.y - jacmat.x.y * jacmat.y.x;
  if (fabs(jacobian) < REL_PREC) {
    ARCANE_FATAL("Cell jacobian is null");
  }

  return Real2x2::fromColumns(jacmat.y.y, -jacmat.y.x, -jacmat.x.y, jacmat.x.x) / jacobian;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute stiffness and mass matrix (only in 2D/3D) and cell forces
void ElastodynamicModule::
_computeKMF3D(const Cell& cell,RealUniqueArray2& Ke, RealUniqueArray2& Me, RealUniqueArray& Fe){

  Integer nb_nodes = cell.nbNode();
  Integer nk{3*nb_nodes};
  Real rho = m_rho(cell);
  Real nu = m_nu(cell);
  Real young = m_young(cell);
  Integer nb{6};
  RealUniqueArray2 B(nb,nk);
  ElastTensor D(young,nu);

  for (Int32 i = 0; i < integ_order.m_i; ++i) {
    for (Int32 j = 0; j < integ_order.m_j; ++j) {
      for (Int32 k = 0; k < integ_order.m_k; ++k) {
        Integer3 indices{ i, j, k };
        auto pos = gausspt.getRefPosition(cell,indices);
        Real jacobian;
        auto ijac = _computeInverseJacobian3D(cell,pos,jacobian);
        auto wt = gausspt.getWeight(cell,indices)*jacobian;

        for (Int32 inod = 0; inod < nb_nodes; ++inod) {
          auto dNi = cell_fem.getShapeFuncDeriv(cell.type(), inod, pos);
          B(0, 0) = B(3, 1) = B(5, 2) += ijac.x.x * dNi.x + ijac.x.y * dNi.y + ijac.x.z * dNi.z;
          B(1, 1) = B(3, 0) = B(4, 2) += ijac.y.x * dNi.x + ijac.y.y * dNi.y + ijac.y.z * dNi.z;
          B(2, 2) = B(4, 1) = B(5, 0) += ijac.z.x * dNi.x + ijac.z.y * dNi.y + ijac.z.z * dNi.z;
        }
        RealUniqueArray2 BtDB = bothMultiply(D,B);
        addArray2(Ke,BtDB, wt);

        for (Int32 inod = 0; inod < nb_nodes; ++inod) {

          auto rhoNi = wt*rho*cell_fem.getShapeFuncVal(cell.type(), inod, pos);
          Fe(3*inod + 0) += rhoNi*gravity.x;
          Fe(3*inod + 1) += rhoNi*gravity.y;
          Fe(3*inod + 2) += rhoNi*gravity.z;
          // TO DO: add input motion value : accgt = imposed acc value at current time step
          /*for (Int32 l = 0; l < ndim; ++l){
             Fe(ndim*inod + l) -= rhoNi*accgt;
            */

          for (Int32 jnod = 0; jnod < nb_nodes; ++jnod) {

            auto Nj = cell_fem.getShapeFuncVal(cell.type(), jnod, pos);
            auto mij = rhoNi*Nj;

            for (Int32 l = 0; l < NDIM; ++l){
              int ii = NDIM*inod + l;
              int jj = NDIM*jnod + l;
              Me(ii,jj) += mij;
            }
          }
        }
      }
    }
  }
}

void ElastodynamicModule::
_computeKMF2D(const Cell& cell,RealUniqueArray2& Ke, RealUniqueArray2& Me, RealUniqueArray& Fe){

  Integer nb_nodes = cell.nbNode();
  Integer nk{NDIM*nb_nodes};
  Real rho = m_rho(cell);
  Real nu = m_nu(cell);
  Real young = m_young(cell);
  Integer nb{3}; // 2D axisymmetry is not considered here
  RealUniqueArray2 B(nb,nk);
  ElastTensor D(young,nu);

  for (Int32 i = 0; i < integ_order.m_i; ++i) {
    for (Int32 j = 0; j < integ_order.m_j; ++j) {
        Integer3 indices{ i, j, 0 };
        auto pos = gausspt.getRefPosition(cell,indices);
        Real jacobian;
        auto ijac = _computeInverseJacobian2D(cell,pos,jacobian);
        auto wt = gausspt.getWeight(cell,indices)*jacobian;

        //----------------------------------------------
        // Elementary Derivation Matrix B
        //----------------------------------------------
        for (Int32 inod = 0; inod < nb_nodes; ++inod) {
          auto dNi = cell_fem.getShapeFuncDeriv(cell.type(), inod, pos);
          B(0,0) = B(3,1) += ijac.x.x*dNi.x + ijac.x.y*dNi.y;
          B(1,1) = B(3,0) += ijac.y.x*dNi.x + ijac.y.y*dNi.y;
        }

        //----------------------------------------------
        // Elementary Stiffness Matrix Ke
        //----------------------------------------------
        RealUniqueArray2 BtDB = bothMultiply(D,B);
        addArray2(Ke,BtDB, wt);

        for (Int32 inod = 0; inod < nb_nodes; ++inod) {

          auto rhoNi = wt*rho*cell_fem.getShapeFuncVal(cell.type(), inod, pos);

          //----------------------------------------------
          // Elementary Force vector Fe assembly
          //----------------------------------------------
          if (options()->hasBodyf()) {
            //----------------------------------------------
            // Body force terms
            //----------------------------------------------
            Fe(NDIM * inod + 0) += rhoNi * gravity.x;
            Fe(NDIM * inod + 1) += rhoNi * gravity.y;
          }

          // TO DO: add input motion value : accgt = imposed acc value at current time step
          /*for (Int32 l = 0; l < ndim; ++l){
             Fe(ndim*inod + l) -= rhoNi*accgt;
            */

          //----------------------------------------------
          // Imposed nodal forces
          //----------------------------------------------
          const Node & nodei = cell.node(inod);
          if ((bool)m_imposed_force[nodei].x)
            Fe(NDIM * inod + 0) += m_force[nodei].x;
          if ((bool)m_imposed_force[nodei].y)
            Fe(NDIM * inod + 1) += m_force[nodei].y;

          //----------------------------------------------
          // Elementary Mass Matrix assembly
          //----------------------------------------------
          for (Int32 jnod = 0; jnod < nb_nodes; ++jnod) {

            auto Nj = cell_fem.getShapeFuncVal(cell.type(), jnod, pos);
            auto mij = rhoNi*Nj;

            for (Int32 l = 0; l < NDIM; ++l){
              int ii = NDIM*inod + l;
              int jj = NDIM*jnod + l;
              Me(ii,jj) += mij;
            }
          }
        }
      }
    }
  }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_assembleLinearGlobal3D()
{
  info() << "Assembly of the FEM linear operator (RHS - vector b) ";

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  info() << "Assembly of the FEM bilinear (LHS - matrix A) and linear (RHS - vector B) operators ";

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    auto nb_nodes = cell.nbNode();
    Int32 nk{ 3 * nb_nodes };
    RealUniqueArray2 Ke(nk, nk);
    RealUniqueArray2 Me(nk, nk);
    RealUniqueArray Fe(nk);

    // Computing mass and stiffness matrices for each cell
    _computeKMF3D(cell, Ke, Me, Fe);

    String dirichletMethod = options()->enforceDirichletMethod();

    //----------------------------------------------
    // Weak Methods to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'P' be the penalty term and let 'i' be the set of DOF for which
    //  Dirichlet condition needs to be applied
    //
    //  to apply the Dirichlet condition on 'i'th DOF
    //  - For LHS matrix A, the diag term corresponding to the Dirichlet DOF
    //  "Penalty" => a_{i,i} = 1. * P
    //  "WeakPenalty" => a_{i,i} = a_{i,i} + P
    //
    //  - For RHS vector b the term that corresponds to the Dirichlet DOF
    //  "Penalty" => b_{i} = b_{i} * P
    //  "WeakPenalty" => b_{i} = b_{i} * P
    //----------------------------------------------
    // "RowElimination" method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'i' be the DOF for which  Dirichlet condition 'g_i' needs to be applied
    //
    //  to apply the Dirichlet on 'i'th DOF
    //  - For LHS matrix A, the row terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j
    //           a_{i,j} = 1.  : i==j
    //  - For RHS vector b the terms corresponding to the Dirichlet DOF
    //           b_i = g_i
    //----------------------------------------------
    // "RowColumnElimination" method to enforce Dirichlet BC
    //----------------------------------------------
    //  Let 'I' be the set of DOF for which  Dirichlet condition needs to be applied
    //
    //  to apply the Dirichlet on 'i'th DOF
    //  - For LHS matrix A, the row terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j  for all j
    //           a_{i,j} = 1.  : i==j
    //    also the column terms corresponding to the Dirichlet DOF
    //           a_{i,j} = 0.  : i!=j  for all i
    //----------------------------------------------

    info() << "Applying Dirichlet boundary condition via "
           << dirichletMethod << " method ";

    // Considering a simple Newmark scheme here (Generalized-alfa will be done later)
    // Computing Me/beta/dt^2 + Ke
    Int32 n1_index = 0;

    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;

      for (Int32 iddl = 0; iddl < 3; ++iddl) {
        auto node1_dofi = node_dof.dofId(node1, iddl);
        auto ii = 3 * n1_index + iddl;
        bool is_node1_dofi_set = (bool)m_imposed_displ[node1][iddl];
        Real dimp_iddl = m_displ[node1][iddl];
        auto rhs_i = Me(ii, ii) / beta / dt2 * dimp_iddl + Fe(ii);

        if (node1.isOwn()) {
          if (is_node1_dofi_set) {
            if (dirichletMethod == "Penalty") {
              m_linear_system.matrixSetValue(node1_dofi, node1_dofi, penalty);
              rhs_values[node1_dofi] = rhs_i * penalty;
            }
            else if (dirichletMethod == "WeakPenalty") {
              m_linear_system.matrixAddValue(node1_dofi, node1_dofi, penalty);
              rhs_values[node1_dofi] = rhs_i * penalty;
            }
            else if (dirichletMethod == "RowElimination") {
              m_linear_system.eliminateRow(node1_dofi, dimp_iddl);
            }
            else if (dirichletMethod == "RowColumnElimination") {
              m_linear_system.eliminateRowColumn(node1_dofi, dimp_iddl);
            }
          }
        }

        for (Node node2 : cell.nodes()) {
          for (Int32 jddl = 0; jddl < 3; ++jddl) {
            auto node2_dofj = node_dof.dofId(node2, jddl);
            auto jj = 3 * n2_index + jddl;
            auto aij = Me(ii, jj) / beta / dt2 + Ke(ii, jj);
            if (node1.isOwn() && !is_node1_dofi_set)
              m_linear_system.matrixAddValue(node1_dofi, node2_dofj, aij);
          }
          ++n2_index;
        }
      }
      ++n1_index;
    }
  }

  //----------------------------------------------
  // Traction terms assembly
  //----------------------------------------------
  for (const auto& bs : options()->neumannBoundaryCondition()) {
    FaceGroup face_group = bs->surface();

    // Loop on faces of the surface
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;

      Real3 trac = m_imposed_traction[face];
      auto fac_el = _computeFacLengthOrArea(face);

      // Loop on nodes of the face
      ENUMERATE_NODE (k, face.nodes()){
        const Node& node = *k;

        for (Int32 iddl = 0; iddl < 3; ++iddl)
        if (!(bool)m_imposed_displ[node][iddl] && node.isOwn()) {
          DoFLocalId dof_id = node_dof.dofId(node, iddl);
          rhs_values[dof_id] += trac[iddl] * fac_el;
        }
      }
    }
  }

  //----------------------------------------------
  // Paraxial terms assembly
  //----------------------------------------------
  auto dt = m_global_deltat();

  for (const auto& bs : options()->paraxialBoundaryCondition()) {
    FaceGroup face_group = bs->surface();
    auto rho = bs->getRhopar();
    Real cs,cp;

    if (bs->hasEPar() && bs->hasNuPar()) {

      auto E = bs->getEPar();
      auto nu = bs->getNuPar();
      auto lambda = nu*E/(1. + nu)/(1. - 2.*nu);
      auto mu = E/2./(1. + nu);
      cp = math::sqrt( (lambda + 2. * mu)/rho );
      cs = math::sqrt( mu/rho );

    } else if (bs->hasCp() && bs->hasCs()) {

      cp = bs->getCp();
      cs = bs->getCp();

    } else if (bs->hasLambdaPar() && bs->hasMuPar()) {

      auto mu = bs->getMuPar();
      cp = math::sqrt( (bs->getLambdaPar() + 2. * mu)/rho );
      cs = math::sqrt( mu/rho );

    }
    else {
      info() << "Elastic properties expected for "
             << "Paraxial boundary condition on FaceGroup "
             << face_group.name() << ": \n"
             << "  - (E-par, nu-par) or\n"
             << "  - (lambda-par, mu-par) or\n"
             << "  - (cp, cs)\n";

      ARCANE_FATAL("Paraxial boundary has no elastic properties ");
    }
    auto rhocsdt = rho*cs/dt;
    auto rhocpsdt{rho*cp/dt-rhocsdt};

    info() << "Applying constant paraxial boundary conditions for surface "<< face_group.name();

    // Loop on the faces (=edges in 2D) concerned with the paraxial condition
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;
      auto facint = _computeFacLengthOrArea(face);
      auto VN = FaceNormal(face,m_node_coord);

      // Loop on nodes of the face or edge (with no Dirichlet condition)
      ENUMERATE_NODE (k, face.nodes()){
        const Node& node = *k;
        auto dU{ m_displ[node] - m_prev_displ[node] };

        // Normal displacement (outer edge direction)
        auto dUN{ math::dot(dU, VN) };

        // Tangential displacement (edge direction)
        auto dUT{ dU - dUN * VN };

        for (Int32 iddl = 0; iddl < NDIM; ++iddl)
          if (!(bool)m_imposed_displ[node][iddl] && node.isOwn()) {
            DoFLocalId dof_id = node_dof.dofId(node, iddl);
            rhs_values[dof_id] -= (rhocpsdt * dUN * VN[iddl] + rhocsdt*dUT[iddl]) * facint;
          }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_assembleLinearGlobal2D()
{
  info() << "Assembly of the FEM linear operator in 2D (RHS - vector b) ";

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  info() << "Assembly of the FEM bilinear (LHS - matrix A) and linear (RHS - vector B) operators ";

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    auto nb_nodes = cell.nbNode();
    Integer nk{ NDIM * nb_nodes };
    RealUniqueArray2 Ke(nk, nk);
    RealUniqueArray2 Me(nk, nk);
    RealUniqueArray Fe(nk);

    // Computing mass and stiffness matrices for each cell
    _computeKMF2D(cell, Ke, Me, Fe);

    String dirichletMethod = options()->enforceDirichletMethod();
    info() << "Applying Dirichlet boundary condition via "
           << dirichletMethod << " method ";

    // Considering a simple Newmark scheme here (Generalized-alfa will be done later)
    // Computing Me/beta/dt^2 + Ke
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;

      for (Int32 iddl = 0; iddl < NDIM; ++iddl) {
        DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
        auto ii = NDIM * n1_index + iddl;
        bool is_node1_dofi_set = (bool)m_imposed_displ[node1][iddl];
        Real dimp_iddl = m_displ[node1][iddl];
        auto rhs_i = Me(ii, ii) / beta / dt2 * dimp_iddl + Fe(ii);

        if (node1.isOwn()) {
          if (is_node1_dofi_set) {
            if (dirichletMethod == "Penalty") {
              m_linear_system.matrixSetValue(node1_dofi, node1_dofi, penalty);
              rhs_values[node1_dofi] = rhs_i * penalty;
            }
            else if (dirichletMethod == "WeakPenalty") {
              m_linear_system.matrixAddValue(node1_dofi, node1_dofi, penalty);
              rhs_values[node1_dofi] = rhs_i * penalty;
            }
            else if (dirichletMethod == "RowElimination") {
              m_linear_system.eliminateRow(node1_dofi, dimp_iddl);
            }
            else if (dirichletMethod == "RowColumnElimination") {
              m_linear_system.eliminateRowColumn(node1_dofi, dimp_iddl);
            }
          }
        }

        for (Node node2 : cell.nodes()) {
          for (Int32 jddl = 0; jddl < NDIM; ++jddl) {
            auto node2_dofj = node_dof.dofId(node2, jddl);
            auto jj = NDIM * n2_index + jddl;
            auto aij = Me(ii, jj) / beta / dt2 + Ke(ii, jj);

            if (node1.isOwn() && !is_node1_dofi_set)
              m_linear_system.matrixAddValue(node1_dofi, node2_dofj, aij);
          }
          ++n2_index;
        }
      }
      ++n1_index;
    }
  }

  //----------------------------------------------
  // Traction terms assembly
  //----------------------------------------------
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

        for (Int32 iddl = 0; iddl < NDIM; ++iddl)
          if (!(bool)m_imposed_displ[node][iddl] && node.isOwn()) {
            DoFLocalId dof_id = node_dof.dofId(node, iddl);
            rhs_values[dof_id] += trac[iddl] * facint;
          }
      }
    }
  }

  //----------------------------------------------
  // Paraxial terms assembly
  //----------------------------------------------
  auto dt = m_global_deltat();

  for (const auto& bs : options()->paraxialBoundaryCondition()) {
    FaceGroup face_group = bs->surface();
    auto rho = bs->getRhopar();
    Real cs,cp;
    bool is_inner{false};

    if (bs->hasEPar() && bs->hasNuPar()) {

      auto E = bs->getEPar();
      auto nu = bs->getNuPar();
      auto lambda = nu*E/(1. + nu)/(1. - 2.*nu);
      auto mu = E/2./(1. + nu);
      cp = math::sqrt( (lambda + 2. * mu)/rho );
      cs = math::sqrt( mu/rho );

    } else if (bs->hasCp() && bs->hasCs()) {

      cp = bs->getCp();
      cs = bs->getCp();

    } else if (bs->hasLambdaPar() && bs->hasMuPar()) {

      auto mu = bs->getMuPar();
      cp = math::sqrt( (bs->getLambdaPar() + 2. * mu)/rho );
      cs = math::sqrt( mu/rho );

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


    info() << "Applying constant paraxial boundary conditions for surface "<< face_group.name();

    // Loop on the faces (=edges in 2D) concerned with the paraxial condition
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;

      if (is_inner){
        const Cell& cell = face.boundaryCell();
        rho = m_rho[cell];
        cs = m_vs[cell];
        cp = m_vp[cell];
      }
      auto rhocsdt = rho*cs/dt;
      auto rhocpsdt{rho*cp/dt-rhocsdt};
      auto facint = _computeFacLengthOrArea(face);
      auto VN = EdgeNormal(face.toEdge(),m_node_coord);

      // Loop on nodes of the face or edge (with no Dirichlet condition)
      ENUMERATE_NODE (k, face.nodes()){
        const Node& node = *k;
        auto dU{ m_displ[node] - m_prev_displ[node] };

        // Normal displacement (outer edge direction)
        auto dUN{ math::dot(dU, VN) };

        // Tangential displacement (edge direction)
        auto dUT{ dU - dUN * VN };

        for (Int32 iddl = 0; iddl < NDIM; ++iddl)
          if (!(bool)m_imposed_displ[node][iddl] && node.isOwn()) {
            DoFLocalId dof_id = node_dof.dofId(node, iddl);
            rhs_values[dof_id] -= (rhocpsdt * dUN * VN[iddl] + rhocsdt*dUT[iddl]) * facint;
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

  // Re-Apply boundary conditions because the solver has modified the values
  // on all nodes
  _applyBoundaryConditions(); // ************ CHECK

  {
    VariableDoFReal& dof_d(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      auto ux = dof_d[node_dof.dofId(node, 0)];
      auto uy = dof_d[node_dof.dofId(node, 1)];
      auto uz = dof_d[node_dof.dofId(node, 2)];
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

