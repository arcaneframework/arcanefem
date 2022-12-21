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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
Int32 NODE_NDDL = 3;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ElastodynamicModule::ElastodynamicModule(const ModuleBuildInfo& mbi)
        : ArcaneElastodynamicObject(mbi), m_cell_fem_dispatch(m_node_coord) {
    ICaseMng *cm = mbi.subDomain()->caseMng();
    cm->setTreatWarningAsError(true);
    cm->setAllowUnkownRootElelement(false);
}

VersionInfo ElastodynamicModule::versionInfo() const {
    return VersionInfo(1, 0, 0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
startInit(){

    info() << "Module Elastodynamic INIT";

    TypesElastodynamic::eAnalysisType anal = options()->getAnalysisType();
    if (anal < TypesElastodynamic::ThreeD)
        NODE_NDDL = 2;
    else
        NODE_NDDL = 3;

    _initBoundaryConditions();
    _initDofs();

    // TO DO : link to dof_variable to initialize linear system...
//    m_linear_system.initialize(subDomain(), dof_variable);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
compute(){

    info() << "Module Elastodynamic COMPUTE";

    // Stop code after computations
    if (m_global_iteration() > 0)
        subDomain()->timeLoopMng()->stopComputeLoop(true);

    _predictNewmark();


    m_linear_system.reset();
//    m_linear_system.initialize(subDomain(), m_node_temperature);

    info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();
    _doSolve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_initDofs(){

    Integer ndim{2};
    if (options()->getAnalysisType() == TypesElastodynamic::ThreeD)
        ndim = 3;

    int j{0};
    ENUMERATE_NODE(inode, allNodes()) {
        Node node = *inode;
        auto node_dof_ids = static_cast<Integer3>(m_node_dof_ids[node]);

        bool b[3] = {(m_node_has_imposed_ax[node] || m_node_has_imposed_vx[node] || m_node_has_imposed_ux[node]),
                     (m_node_has_imposed_ay[node] || m_node_has_imposed_vy[node] || m_node_has_imposed_uy[node]),
                     (m_node_has_imposed_az[node] || m_node_has_imposed_vz[node] || m_node_has_imposed_uz[node])};

        for (int i = 0; i < ndim; ++i) {
           if (!b[i])
               node_dof_ids[i] = j++;
           else
               node_dof_ids[i] = -1;
        }
    }
    m_nb_neqs = j;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::_predictNewmark(){

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

    Integer ndim{2};
    if (options()->getAnalysisType() == TypesElastodynamic::ThreeD)
        ndim = 3;

    ENUMERATE_NODE(inode, allNodes()){
        Node node = *inode;
        auto node_dof_ids = static_cast<Integer3>(m_node_dof_ids[node]);
        auto an = m_prev_acceleration[node];
        auto vn = m_prev_velocity[node];
        auto dn = m_prev_displacement[node];

        Real3 d,v;
        for (int i = 0; i < ndim; ++i) {
            if (node_dof_ids[i] != -1) {
                d[i] = dn[i] + dt * vn[i] + dt2 * (0.5 - beta) * an[i];
                v[i] = vn[i] + dt * (1. - gamma) * an[i];
            }
        }
        m_displacement[node] = d;
        m_velocity[node] = v;

        // TO DO: predict corresponding dof values...
    }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::_updateNewmark(){
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
    Integer ndim{2};
    if (options()->getAnalysisType() == TypesElastodynamic::ThreeD)
        ndim = 3;

    ENUMERATE_NODE(inode, allNodes()){
        Node node = *inode;
        auto node_dof_ids = static_cast<Integer3>(m_node_dof_ids[node]);
        Real3 an = m_prev_acceleration[node];
        Real3 vn = m_prev_velocity[node];
        Real3 dn = m_prev_displacement[node];
        Real3 d = m_displacement[node];

        Real3 a,v;
        if (!is_alfa_method) {
            for (int i = 0; i < ndim; ++i) {
                if (node_dof_ids[i] != -1) {
                    a[i] = (d[i] - (dn[i] + dt * vn[i] + dt2 * (0.5 - beta) * an[i])) / beta / dt2;
                    v[i] = vn[i] + dt * (1. - gamma) * an[i] + dt * gamma * a[i];
                }
            }
        } else {
            // TO DO
        }
        m_acceleration[node] = a;
        m_velocity[node] = v;
    }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ElastodynamicModule::
_initBoundaryConditions(){

    info() << " Module Elastodynamic INIT BOUNDARY CONDITIONS";
    _applyBoundaryConditions();
}

void ElastodynamicModule::
_applyBoundaryConditions(){
    for (Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i)
    {
        FaceGroup face_group = options()->boundaryCondition[i]->surface();
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
                Real3 & velocity = m_velocity[node];
                Real3 & displacement = m_displacement[node];
                Real3 & acceleration = m_acceleration[node];
                Real3 & force = m_force[node];

                switch (type)
                {
                    case TypesElastodynamic::AccelerationX:
                        acceleration.x = value;
                        break;
                    case TypesElastodynamic::AccelerationY:
                        acceleration.y = value;
                        break;
                    case TypesElastodynamic::AccelerationZ:
                        acceleration.z = value;
                        break;
                    case TypesElastodynamic::DisplacementX:
                        displacement.x = value;
                        break;
                    case TypesElastodynamic::DisplacementY:
                        displacement.y = value;
                        break;
                    case TypesElastodynamic::DisplacementZ:
                        displacement.z = value;
                        break;
                    case TypesElastodynamic::VelocityX:
                        velocity.x = value;
                        break;
                    case TypesElastodynamic::VelocityY:
                        velocity.y = value;
                        break;
                    case TypesElastodynamic::VelocityZ:
                        velocity.z = value;
                        break;
                    case TypesElastodynamic::ForceX:
                        force.x = value;
                        break;
                    case TypesElastodynamic::ForceY:
                        force.y = value;
                        break;
                    case TypesElastodynamic::ForceZ:
                        force.z = value;
                        break;
                    case TypesElastodynamic::Unknown:
                        break;
                }
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ElastodynamicModule::
_initInitialConditions(){

    info() << " Module Elastodynamic INIT INITIAL CONDITIONS";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ElastodynamicModule::
_applyInputMotion(){
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_computeRHS(){
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
_doSolve(){
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_ELASTODYNAMIC(ElastodynamicModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

