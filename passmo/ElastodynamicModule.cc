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
        : ArcaneElastodynamicObject(mbi) {
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

    m_global_deltat = options()->deltat();
    TypesElastodynamic::eAnalysisType anal = options()->getAnalysisType();
    if (anal < TypesElastodynamic::ThreeD)
        NODE_NDDL = 2;
    else
        NODE_NDDL = 3;


//    m_linear_system.initialize(subDomain(), m_node_temperature);
    //m_k_matrix.resize(nb_node, nb_node);
    //m_k_matrix.fill(0.0);

    //m_rhs_vector.resize(nb_node);
    //m_rhs_vector.fill(0.0);

    // # init behavior
    // # init behavior on mesh entities
    // # init BCs
    _initBoundaryConditions();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
compute(){

    info() << "Module Elastodynamic COMPUTE";

    // Stop code after computations
    if (m_global_iteration() > 0)
        subDomain()->timeLoopMng()->stopComputeLoop(true);

    m_linear_system.reset();
//    m_linear_system.initialize(subDomain(), m_node_temperature);

    info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();
    _doSolve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
predict() {

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void ElastodynamicModule::
update() {
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

