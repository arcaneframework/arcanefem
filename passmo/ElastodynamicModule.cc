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
  m_linear_system.initialize(subDomain(), m_dofs_on_nodes.dofFamily(), "Solver");
  _initDofs();
  _applyInitialNodeConditions();

  // This will be useful for nonlinear dynamics only (not used in elastodynamics)
//  _applyInitialCellConditions();

  _initBoundaryConditions();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_initDofs(){

  Int32 nb_dofs_per_node = options()->getNb_dofs_per_node();
  m_dofs_on_nodes.initialize(mesh(),nb_dofs_per_node);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_applyInitialNodeConditions(){

  for (Int32 i = 0, nb = options()->initialNodeCondition().size(); i < nb; ++i) {

    NodeGroup node_group = options()->initialNodeCondition[i]->nodeGroup();
    TypesElastodynamic::eNodeCondition type = options()->initialNodeCondition[i]->type();
    Real3 values = options()->initialNodeCondition[i]->values();

    // Loop on nodes with this initial condition
    ENUMERATE_NODE(inode, node_group) {
      const Node & node = *inode;
      switch (type) {

      case TypesElastodynamic::Acceleration:
        m_prev_acceleration[node] = values;
        break;
      case TypesElastodynamic::Displacement:
        m_prev_displacement[node] = values;
        break;
      case TypesElastodynamic::Velocity:
        m_prev_velocity[node] = values;
      case TypesElastodynamic::Force:
        m_force[node] = values;
      case TypesElastodynamic::UnknownCond:
        break;

      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_applyInitialCellConditions(){

  for (Integer i = 0, nb = options()->initialCellCondition().size(); i < nb; ++i) {

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
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
compute(){

    info() << "Module Elastodynamic COMPUTE";

    // Stop code after computations
    if (m_global_iteration() > 0)
        subDomain()->timeLoopMng()->stopComputeLoop(true);

    // Apply the input motion: imposed time history on groups of nodes/faces (acceleration, velocity or displacement)
    _applyInputMotion();

    // Apply other Dirichlet/Neumann conditions if any (constant values assumed at present)
    _applyBoundaryConditions();

    // Predict the nodal variable according to the integration scheme (e.g. Newmark)
    _predictNewmark();

    info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();

    // Assemble the FEM global operators (LHS matrix/RHS vector b)
    _assembleLinearGlobal();

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

  // Prediction of the nodal displacements and velocities (before solve) with
  // Newmark-beta or Generalized-alfa time integration schemes
  Real dt = options()->getDeltat();
  Real dt2 = dt*dt;
  bool is_alfa_method = options()->is_alfa_method();
  Real beta{0.}, gamma{0.};
  Real alfam{0.}, alfaf{0.};

  if (!is_alfa_method) {
    beta = options()->getBeta();
    gamma = options()->getGamma();
  } else {
    alfam = options()->getAlfam();
    alfaf = options()->getAlfaf();
    gamma = 0.5 + alfaf - alfam;
    beta = 0.5*pow(0.5 + gamma,2);
  }

  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    auto an = m_prev_acceleration[node];
    auto vn = m_prev_velocity[node];
    auto dn = m_prev_displacement[node];

    for (Int32 i = 0; i < 3; ++i) {

      if (!(bool)m_node_has_imposed_displ[node][i])
        m_displacement[node][i] = dn[i] + dt * vn[i] + dt2 * (0.5 - beta) * an[i];

      if (!(bool)m_node_has_imposed_vel[node][i])
        m_velocity[node][i] = vn[i] + dt * (1. - gamma) * an[i];
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_updateNewmark(){

  // Updating the nodal accelerations and velocities (after solve) with
  // Newmark-beta or Generalized-alfa time integration schemes
  Real dt = options()->getDeltat();
  Real dt2 = dt*dt;
  bool is_alfa_method = options()->is_alfa_method();
  Real beta{0.}, gamma{0.};
  Real alfam{0.}, alfaf{0.};

  if (!is_alfa_method) {
    beta = options()->getBeta();
    gamma = options()->getGamma();
  } else {
    alfam = options()->getAlfam();
    alfaf = options()->getAlfaf();
    gamma = 0.5 + alfaf - alfam;
    beta = 0.5*pow(0.5 + gamma,2);
  }

  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    auto an = m_prev_acceleration[node];
    auto vn = m_prev_velocity[node];
    auto dn = m_prev_displacement[node];

    if (!is_alfa_method) {
      for (Int32 i = 0; i < 3; ++i) {

        if (!(bool)m_node_has_imposed_acc[node][i])
          m_acceleration[node][i] = (m_displacement[node][i] - (dn[i] + dt*vn[i] + dt2*(0.5 - beta)*an[i]))/beta/dt2;

        if (!(bool)m_node_has_imposed_vel[node][i])
          m_velocity[node][i] = vn[i] + dt*( (1. - gamma)*an[i] + gamma*m_acceleration[node][i] );
      }
    } else {
      // TO DO
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_initBoundaryConditions(){
  for (Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i)
  {
    FaceGroup face_group = options()->boundaryCondition[i]->surface();
    NodeGroup node_group = options()->boundaryCondition[i]->nodeGroup();
    Real value = options()->boundaryCondition[i]->value();
    TypesElastodynamic::eBoundaryCondition type = options()->boundaryCondition[i]->type();

    // Loop on faces of the surface
    ENUMERATE_FACE(j, face_group)
    {
      const Face & face = * j;
      Integer nb_node = face.nbNode();

      // Loop on nodes of the face
      for (Integer k = 0; k < nb_node; ++k)
      {
        const Node & node = face.node(k);

        switch (type)
        {
        case TypesElastodynamic::AccelerationX:
          m_node_has_imposed_acc[node].x = 1;
          break;
        case TypesElastodynamic::AccelerationY:
          m_node_has_imposed_acc[node].y = 1;
          break;
        case TypesElastodynamic::AccelerationZ:
          m_node_has_imposed_acc[node].z = 1;
          break;
        case TypesElastodynamic::DisplacementX:
          m_node_has_imposed_displ[node].x = 1;
          break;
        case TypesElastodynamic::DisplacementY:
          m_node_has_imposed_displ[node].y = 1;
          break;
        case TypesElastodynamic::DisplacementZ:
          m_node_has_imposed_displ[node].z = 1;
          break;
        case TypesElastodynamic::VelocityX:
          m_node_has_imposed_vel[node].x = 1;
          break;
        case TypesElastodynamic::VelocityY:
          m_node_has_imposed_vel[node].y = 1;
          break;
        case TypesElastodynamic::VelocityZ:
          m_node_has_imposed_vel[node].z = 1;
          break;
        case TypesElastodynamic::ForceX:
          m_node_has_imposed_force[node].x = 1;
          break;
        case TypesElastodynamic::ForceY:
          m_node_has_imposed_force[node].y = 1;
          break;
        case TypesElastodynamic::ForceZ:
          m_node_has_imposed_force[node].z = 1;
          break;
        case TypesElastodynamic::Unknown:
          break;
        }
      }
    }

    // Loop on nodes
    ENUMERATE_NODE(inode, node_group)
    {
      const Node & node = *inode;

      switch (type)
      {
      case TypesElastodynamic::AccelerationX:
        m_node_has_imposed_acc[node].x = 1;
        break;
      case TypesElastodynamic::AccelerationY:
        m_node_has_imposed_acc[node].y = 1;
        break;
      case TypesElastodynamic::AccelerationZ:
        m_node_has_imposed_acc[node].z = 1;
        break;
      case TypesElastodynamic::DisplacementX:
        m_node_has_imposed_displ[node].x = 1;
        break;
      case TypesElastodynamic::DisplacementY:
        m_node_has_imposed_displ[node].y = 1;
        break;
      case TypesElastodynamic::DisplacementZ:
        m_node_has_imposed_displ[node].z = 1;
        break;
      case TypesElastodynamic::VelocityX:
        m_node_has_imposed_vel[node].x = 1;
        break;
      case TypesElastodynamic::VelocityY:
        m_node_has_imposed_vel[node].y = 1;
        break;
      case TypesElastodynamic::VelocityZ:
        m_node_has_imposed_vel[node].z = 1;
        break;
      case TypesElastodynamic::ForceX:
        m_node_has_imposed_force[node].x = 1;
        break;
      case TypesElastodynamic::ForceY:
        m_node_has_imposed_force[node].y = 1;
        break;
      case TypesElastodynamic::ForceZ:
        m_node_has_imposed_force[node].z = 1;
        break;
      case TypesElastodynamic::Unknown:
        break;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_applyBoundaryConditions(){

  for (Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i)
  {
    FaceGroup face_group = options()->boundaryCondition[i]->surface();
    NodeGroup node_group = options()->boundaryCondition[i]->nodeGroup();
    Real value = options()->boundaryCondition[i]->value();
    TypesElastodynamic::eBoundaryCondition type = options()->boundaryCondition[i]->type();

    // Loop on faces of the surface
    ENUMERATE_FACE(j, face_group)
    {
      const Face & face = * j;
      Integer nb_node = face.nbNode();

      // Loop on nodes of the face
      for (Integer k = 0; k < nb_node; ++k)
      {
        const Node & node = face.node(k);

        switch (type)
        {
        case TypesElastodynamic::AccelerationX:
          m_acceleration[node].x = value;
          break;
        case TypesElastodynamic::AccelerationY:
          m_acceleration[node].y = value;
          break;
        case TypesElastodynamic::AccelerationZ:
          m_acceleration[node].z = value;
          break;
        case TypesElastodynamic::DisplacementX:
          m_displacement[node].x = value;
          break;
        case TypesElastodynamic::DisplacementY:
          m_displacement[node].y = value;
          break;
        case TypesElastodynamic::DisplacementZ:
          m_displacement[node].z = value;
          break;
        case TypesElastodynamic::VelocityX:
          m_velocity[node].x = value;
          break;
        case TypesElastodynamic::VelocityY:
          m_velocity[node].y = value;
          break;
        case TypesElastodynamic::VelocityZ:
          m_velocity[node].z = value;
          break;
        case TypesElastodynamic::ForceX:
          m_force[node].x = value;
          break;
        case TypesElastodynamic::ForceY:
          m_force[node].y = value;
          break;
        case TypesElastodynamic::ForceZ:
          m_force[node].z = value;
          break;
        case TypesElastodynamic::Unknown:
          break;
        }
      }
    }

    // Loop on nodes
    ENUMERATE_NODE(inode, node_group)
    {
      const Node & node = *inode;

      switch (type)
      {
      case TypesElastodynamic::AccelerationX:
        m_acceleration[node].x = value;
        break;
      case TypesElastodynamic::AccelerationY:
        m_acceleration[node].y = value;
        break;
      case TypesElastodynamic::AccelerationZ:
        m_acceleration[node].z = value;
        break;
      case TypesElastodynamic::DisplacementX:
        m_displacement[node].x = value;
        break;
      case TypesElastodynamic::DisplacementY:
        m_displacement[node].y = value;
        break;
      case TypesElastodynamic::DisplacementZ:
        m_displacement[node].z = value;
        break;
      case TypesElastodynamic::VelocityX:
        m_velocity[node].x = value;
        break;
      case TypesElastodynamic::VelocityY:
        m_velocity[node].y = value;
        break;
      case TypesElastodynamic::VelocityZ:
        m_velocity[node].z = value;
        break;
      case TypesElastodynamic::ForceX:
        m_force[node].x = value;
        break;
      case TypesElastodynamic::ForceY:
        m_force[node].y = value;
        break;
      case TypesElastodynamic::ForceZ:
        m_force[node].z = value;
        break;
      case TypesElastodynamic::Unknown:
        break;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ElastodynamicModule::
_applyInputMotion(){
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Computes the Jacobian Matrix of a finite-element at a Gauss point given by its local coordinates in the element
Real3x3 ElastodynamicModule::
_computeJacobianMatrix(const Cell& cell,const Real3& ref_coord) {
  Integer	n = cell.nbNode();
  CellFEMDispatcher cell_fem(m_node_coord);

  // Jacobian matrix computed at the integration point
  Real3x3	jac;
  Int16 cell_type = cell.type();

  for (Integer inod = 0; inod < n; ++inod) {

    // vector of local derivatives at this integration point, for node inod
    auto dNi = cell_fem.getShapeFuncDeriv(cell_type, inod, ref_coord);
    auto coord_nod = m_node_coord[cell.node(inod)];

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        jac[i][j] += dNi[j] * coord_nod[i];
      }
    }
  }
  return jac;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Computes the Inverse Jacobian Matrix of a finite-element in 2D/3D only
Real3x3 ElastodynamicModule::
_computeInverseJacobianMatrix(const Real3x3& jacmat, const Real& jacobian, const Integer& ndim) {

  if (ndim == 2) {
    auto ijac = Real2x2::fromColumns(jacmat.y.y, -jacmat.y.x, -jacmat.x.y, jacmat.x.x) / jacobian;
    return {Real3(ijac.x.x,ijac.x.y,0.),Real3(ijac.y.x,ijac.y.y,0.),Real3::zero()};
  }

  // ndim = 3
  return math::inverseMatrix(jacmat);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real ElastodynamicModule::
_computeJacobian(const Real3x3& jacmat, const Integer& ndim){

  Real jac{0.};

  if (ndim == 2)
    jac = jacmat.x.x * jacmat.y.y - jacmat.x.y * jacmat.y.x;
  else
    jac = math::matrixDeterminant(jacmat);

  if (fabs(jac) < REL_PREC) {
    ARCANE_FATAL("Cell jacobian is null");
  }
  return jac;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// ! Compute stiffness and mass matrix (only in 2D/3D) and cell forces
void ElastodynamicModule::
_computeKMF(const Cell& cell,RealUniqueArray2& Ke, RealUniqueArray2& Me, RealUniqueArray& Fe){

  Int32 nint1 = options()->gauss_order1();
  Int32 nint2 = options()->gauss_order2();
  Int32 nint3 = options()->gauss_order3();
  Integer3 integ_order(nint1,nint2,nint3);
  Int32 ndim{2};
  if (options()->getAnalysisType() == TypesElastodynamic::ThreeD)
    ndim = 3;
  if (ndim == 2) integ_order[3] = 0;

  Integer nb_nodes = cell.nbNode();
  Integer nk{ndim*nb_nodes};
  Real gravity{-9.81};
  Real density = m_density(cell);
  Real poisson = m_poisson_ratio(cell);
  Real young = m_young_modulus(cell);
  CellFEMDispatcher cell_fem(m_node_coord);
  Integer nb{6};
  if (ndim == 2) nb = 3; // 2D axisymmetry is not considered here
  RealUniqueArray2 B(nb,nk);
  ElastTensor D(young,poisson);

  for (Int32 i = 0; i < nint1; ++i) {
    for (Int32 j = 0; j < nint2; ++j) {
      for (Int32 k = 0; k < nint3; ++k) {
        Integer3 indices{ i, j, k };
        GaussPointDispatcher gausspt(indices, integ_order);
        auto pos = gausspt.getRefPosition(cell);
        auto jac = _computeJacobianMatrix(cell, pos);
        auto jacobian = _computeJacobian(jac,ndim);
        auto wt = gausspt.getWeight(cell)*jacobian;
        auto ijac = _computeInverseJacobianMatrix(jac,jacobian,ndim);

        if (ndim == 2){
          for (Int32 inod = 0; inod < nb_nodes; ++inod) {
            auto dNi = cell_fem.getShapeFuncDeriv(cell.type(), inod, pos);
            B(0,0) = B(3,1) += ijac.x.x*dNi.x + ijac.x.y*dNi.y;
            B(1,1) = B(3,0) += ijac.y.x*dNi.x + ijac.y.y*dNi.y;
          }
        }
        else { //ndim = 3
          for (Int32 inod = 0; inod < nb_nodes; ++inod) {
            auto dNi = cell_fem.getShapeFuncDeriv(cell.type(), inod, pos);
            B(0, 0) = B(3, 1) = B(5, 2) += ijac.x.x * dNi.x + ijac.x.y * dNi.y + ijac.x.z * dNi.z;
            B(1, 1) = B(3, 0) = B(4, 2) += ijac.y.x * dNi.x + ijac.y.y * dNi.y + ijac.y.z * dNi.z;
            B(2, 2) = B(4, 1) = B(5, 0) += ijac.z.x * dNi.x + ijac.z.y * dNi.y + ijac.z.z * dNi.z;
          }
        }
        RealUniqueArray2 BtDB = bothMultiply(D,B);
        addArray2(Ke,BtDB, wt);

        for (Int32 inod = 0; inod < nb_nodes; ++inod) {

          auto rhoNi = wt*density*cell_fem.getShapeFuncVal(cell.type(), inod, pos);
          Fe(ndim*inod + 2) += rhoNi*gravity;
          // TO DO: add input motion value : accgt = imposed acc value at current time step
          /*for (Int32 l = 0; l < ndim; ++l){
             Fe(ndim*inod + l) -= rhoNi*accgt;
            */

          for (Int32 jnod = 0; jnod < nb_nodes; ++jnod) {

            auto Nj = cell_fem.getShapeFuncVal(cell.type(), jnod, pos);
            auto mij = rhoNi*Nj;

            for (Int32 l = 0; l < ndim; ++l){
              int ii = ndim*inod + l;
              int jj = ndim*jnod + l;
              Me(ii,jj) += mij;
            }
          }
        }
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_assembleLinearGlobal()
{
  info() << "Assembly of the FEM linear operator (RHS - vector b) ";
  // Reset is done for elastodynamics only
  m_linear_system.reset();

  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  Int32 ndim{2};
  if (options()->getAnalysisType() == TypesElastodynamic::ThreeD)
    ndim = 3;

  auto gamma = options()->getGamma();
  auto beta = options()->getBeta();
  auto dt = options()->getDeltat();
  auto dt2 = dt*dt;

  info() << "Assembly of the FEM bilinear (LHS - matrix A) and linear (RHS - vector B) operators ";

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    Integer nb_nodes = cell.nbNode();
    Integer nk{ ndim * nb_nodes };
    RealUniqueArray2 Ke(nk, nk);
    RealUniqueArray2 Me(nk, nk);
    RealUniqueArray Fe(nk);

    // Computing mass and stiffness matrices for each cell
    _computeKMF(cell, Ke, Me, Fe);

    // Considering a simple Newmark sheme here (Generalized-alfa will be done later)
    // Computing Me/beta/dt^2 + Ke
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {

      Int32 n2_index = 0;
      for (Int32 iddl = 0; iddl < ndim; ++iddl) {

        DoFLocalId node1_dofi = node_dof.dofId(node1, iddl);
        auto ii = ndim * n1_index + iddl;
        bool is_node1_dofi_set = (bool)m_node_has_imposed_displ[node1][iddl];
        bool is_node1_acci_set = (bool)m_node_has_imposed_acc[node1][iddl];

        if (node1.isOwn()) {
          auto rhs_i = Me(ii, ii) / beta / dt2*m_displacement[node1][iddl] + Fe(ii);
          if (is_node1_dofi_set) {
            //----------------------------------------------
            // penalty method for assembly of Dirichlet BC
            //----------------------------------------------
            rhs_i *= 1.0e40;
          }
          rhs_values[node1_dofi] = rhs_i;
        }

        for (Node node2 : cell.nodes()) {
          for (Int32 jddl = 0; jddl < ndim; ++jddl) {

            DoFLocalId node2_dofj = node_dof.dofId(node2, jddl);
            auto jj = ndim * n2_index + jddl;
            auto aij = Me(ii, jj) / beta / dt2 + Ke(ii, jj);

            if (node1.isOwn()) {
              if (is_node1_dofi_set && ii == jj) {
                //----------------------------------------------
                // penalty method for assembly of Dirichlet BC
                //----------------------------------------------
                aij *= 1.0e40;
              }
              m_linear_system.matrixAddValue(node1_dofi, node2_dofj, aij);
            }
          }
          ++n2_index;
        }
      }
      ++n1_index;
    }
  }
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
    VariableDoFReal& dof_displacement(m_linear_system.solutionVariable());
    auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
    ENUMERATE_ (Node, inode, ownNodes()) {
      Node node = *inode;
      auto ux = dof_displacement[node_dof.dofId(node, 0)];
      auto uy = dof_displacement[node_dof.dofId(node, 1)];
      auto uz = dof_displacement[node_dof.dofId(node, 2)];
      info() << "Node: " << node.localId() << " Ux=" << ux << " Uy=" << uy << " Uz=" << uz;
    }
  }

  m_displacement.synchronize();
  m_velocity.synchronize();
  m_acceleration.synchronize();
  const bool do_print = (allNodes().size() < 200);
  if (do_print) {
    int p = std::cout.precision();
    std::cout.precision(17);
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      std::cout << "U[" << node.localId() << "][" << node.uniqueId() << "] = "
                << m_displacement[node].x << "  " << m_displacement[node].y
                << "  " << m_displacement[node].z << "\n";
    }
    std::cout.precision(p);
  }

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ARCANE_REGISTER_MODULE_ELASTODYNAMIC(ElastodynamicModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

