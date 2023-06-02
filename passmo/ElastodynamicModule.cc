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
  _initDofs();
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

/*  for (Integer i = 0, nb = options()->initialCellCondition().size(); i < nb; ++i) {

    CellGroup cell_group = options()->initialCellCondition[i]->cellGroup();
    Real3 stress = options()->initialCellCondition[i]->stress();
    Real3 sstress = options()->initialCellCondition[i]->shear_stress();
    Real3 strain = options()->initialCellCondition[i]->strain();
    Real3 sstrain = options()->initialCellCondition[i]->shear_strain();

    // Loop on nodes with this initial condition
    ENUMERATE_CELL(icell, cell_group) {
      const Cell & cell = *icell;

      // Initialize the stress tensor for the concerned cell
      m_cell_stress[cell].x = Real3(stress.x,sstress.x,sstress.y);
      m_cell_stress[cell].y = Real3(sstress.x,stress.y,sstress.z);
      m_cell_stress[cell].z = Real3(sstress.y,sstress.z,stress.z);

      // Initialize the strain tensor for the concerned cell
      m_cell_strain[cell].x = Real3(strain.x,sstrain.x,sstrain.y);
      m_cell_strain[cell].y = Real3(sstrain.x,strain.y,sstrain.z);
      m_cell_strain[cell].z = Real3(sstrain.y,sstrain.z,strain.z);
    }
  }*/
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
      if (m_linear_system.isInitialized() && linop_nstep_counter < linop_nstep){
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
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*void ElastodynamicModule::
_initInputMotion(){

  IParallelMng* pm = subDomain()->parallelMng();
  Int32 ndim{2};
  if (options()->getAnalysisType() == TypesElastodynamic::ThreeD)
    ndim = 3;

  // The following to be potentially removed ------------------------------
  if (options()->inputMotion().size()==1)
    ARCANE_FATAL("Only one input motion is to be defined in the arc file");
  //-----------------------------------------------------------------------

  for (Int32 i = 0, nb = options()->inputMotion().size(); i < nb; ++i) {
    m_input.m_node_group = options()->inputMotion[i]->nodeGroup();

    String fa = options()->inputMotion[i]->accelerationInputFile();
    m_input.m_acc = readFileAsCaseTable(pm, fa, ndim);

    m_input.m_is_vel = options()->inputMotion[i].hasVelocityInputFile();
    if (m_input.m_is_vel ) {
      String fv = options()->inputMotion[i]->velocityInputFile();
      m_input.m_vel = readFileAsCaseTable(pm, fv, ndim);
    }

    m_input.m_is_displ = options()->inputMotion[i].hasDisplacementInputFile();
    if (m_input.m_is_displ) {
      String fu = options()->inputMotion[i]->displacementInputFile();
      m_input.m_displ = readFileAsCaseTable(pm, fu, ndim);
    }
    m_input.m_rigid_base = options()->inputMotion[i]->hasRigidBase();
    m_input.m_ampli_factors = options()->inputMotion[i]->amplificationFactors();
  }
  m_input.m_max_frequency = options()->getMaxFrequency();

  // Loop on nodes
  ENUMERATE_NODE(inode, m_input.m_node_group){
    const Node & node = *inode;
    if (m_input.m_rigid_base) {
      for (Int32 i = 0; i < ndim; ++i) {
        if (m_input.m_component[i]) {
          m_node_has_imposed_acc[node][i] = 1;
          if (m_input.m_is_vel)
            m_node_has_imposed_vel[node][i] = 1;
          if (m_input.m_is_displ)
            m_node_has_imposed_displ[node][i] = 1;
        }
      }
    }
  }*/
  // Print some values
/*{
    Real3 test_value;
    Real param = 1.0e-4;
    table->value(param, test_value);
    tm->info() << "V1 t=" << param << " v=" << test_value;
  }
  {
    Real3 test_value;
    Real param = 1.2e-3;
    table->value(param, test_value);
    tm->info() << "V2 t=" << param << " v=" << test_value;
  }*/
//}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*void ElastodynamicModule::
_applyInputMotion(){

  Int32 ndim{2};
  if (options()->getAnalysisType() == TypesElastodynamic::ThreeD)
    ndim = 3;

  Real time = globalTime();

  // Loop on nodes
  ENUMERATE_NODE(inode, m_input.m_node_group) {
    const Node& node = *inode;
    Real3 values;
    m_input.m_acc->value(time, values);
    for (Int32 i = 0; i < ndim; ++i) {
      if (m_input.m_component[i])
        m_acc[node][i] = values[i];
    }
    if (m_input.m_is_vel) {
      m_input.m_vel->value(time, values);
      for (Int32 i = 0; i < ndim; ++i) {
        if (m_input.m_component[i])
          m_vel[node][i] = values[i];
      }
    }

    if (m_input.m_is_displ) {
      m_input.m_displ->value(time, values);
      for (Int32 i = 0; i < ndim; ++i) {
        if (m_input.m_component[i])
          m_displ[node][i] = values[i];
      }
    }
  }
}*/
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
      auto fac_el = _computeTracFac(face);

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

    // Loop on faces of the surface
    ENUMERATE_FACE (j, face_group) {
      const Face& face = *j;

      Real3 trac = m_imposed_traction[face];
      auto fac_el = _computeTracFac(face);

      // Loop on nodes of the face
      ENUMERATE_NODE (k, face.nodes()){
        const Node& node = *k;

        for (Int32 iddl = 0; iddl < NDIM; ++iddl)
          if (!(bool)m_imposed_displ[node][iddl] && node.isOwn()) {
            DoFLocalId dof_id = node_dof.dofId(node, iddl);
            rhs_values[dof_id] += trac[iddl] * fac_el;
          }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real ElastodynamicModule::
_computeTracFac(const Face& face)
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
  _applyBoundaryConditions();

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

