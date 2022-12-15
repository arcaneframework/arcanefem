//
// Created by EF on 12/12/22.
//

#ifndef PASSMO_NODEDOFS_H
#define PASSMO_NODEDOFS_H

#include "arcane/mesh/DoFManager.h"
#include "arcane/mesh/DoFFamily.h"
#include "arcane/mesh/ItemConnectivity.h"
#include "arcane/mesh/ItemConnectivityMng.h"
#include "arcane/mesh/GhostLayerFromConnectivityComputer.h"

#include "arcane/IItemConnectivitySynchronizer.h"
#include "arcane/IMesh.h"
#include "arcane/VariableTypedef.h"
#include "arcane/VariableBuildInfo.h"
#include "arcane/ItemEnumerator.h"
#include "arcane/ItemUniqueId.h"
#include "arcane/VariableTypes.h"
#include "arcane/ItemTypes.h"
#include "arcane/mesh/CellFamily.h"
#include "arcane/ItemVector.h"
#include "arcane/Item.h"
#include "arcane/mesh/ParticleFamily.h"
#include "arcane/mesh/NodeFamily.h"
#include "arcane/IParallelMng.h"
#include "arcane/IMeshModifier.h"
#include "arcane/VariableCollection.h"
#include <arcane/IParticleExchanger.h>

#include "INodeDoFs.h"
#include "NodeDoFs_axl.h"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
extern Integer NODE_NDDL;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class NodeDoFsService:
        public ArcaneNodeDoFsObject {

public:

    explicit NodeDoFsService(const ServiceBuildInfo& sbi)
            : ArcaneNodeDoFsObject(sbi)
            , m_connectivity_mng(sbi.subDomain()->traceMng())
            , m_dof_mng(sbi.mesh(),&m_connectivity_mng)
            , m_dof_family_name("DoFFamily")
            , m_dof_on_node_family_name("DoFOnNode")
            , m_dofs_on_node_family_name("DoFsOnNode")
            , m_dofs_multi_on_node_family_name("DoFsMultiOnNode"){}

     virtual ~ NodeDoFsService() {};

public:
    typedef ItemConnectivityT          <Node, DoF> NodeToDoFConnectivity;
    typedef ItemArrayConnectivityT     <Node, DoF> NodeToDoFsConnectivity;
    typedef ItemMultiArrayConnectivityT<Node, DoF> NodeToDoFsMultiConnectivity;

    typedef SharedArray2<Int32> Int32SharedArray2;
    typedef SharedArray2<Int64> Int64SharedArray2;

public:

    Integer addDoF(const Integer& begin_uid) override;
    Integer removeDoF(const Integer& begin_lid) override;
    Integer addDoFs(const Integer& begin_uid, const Integer& size) override;
    Integer removeDoFs(const Integer& begin_lid, const Integer& size) override;
    DoFGroup createDoFGroup(const Integer& begin_lid, const Integer& size) override;
    VariableDoFReal doFVariable() override;
    VariableDoFArrayReal doFVectorVariable(const Integer& size) override;

private:

    DoFManager& dofMng() {return m_dof_mng;}
    mesh::DoFFamily& getNodeDof() { return m_dof_mng.family(m_dof_on_node_family_name); }
    mesh::DoFFamily& getNodeDofs() { return m_dof_mng.family(m_dofs_on_node_family_name); }
    mesh::DoFFamily& getMultiNodeDofs() { return m_dof_mng.family(m_dofs_multi_on_node_family_name); }
    // Connectivity tests
    void _node2DoFsConnectivityRegistered();
    void _node2DoFsMultiConnectivityRegistered();
    void _removeGhost(mesh::DoFFamily& dof_family);
    void _addNodes(Int32Array2View new_nodes_lids, Int64Array2View new_nodes_uids);
    void _removeNodes(Int32ConstArray2View new_nodes_lids,
                      const Integer nb_removed_nodes,
                      Int32SharedArray2& removed_node_lids,
                      Int64SharedArray2& removed_node_uids,
                      Int32SharedArray2& remaining_node_lids,
                      Int64SharedArray2& remaining_node_uids);

private:
    ItemConnectivityMng m_connectivity_mng;
    DoFManager m_dof_mng;
    String m_dof_family_name;
    String m_dof_on_node_family_name;  // One dof per node
    String m_dofs_on_node_family_name; // Several dofs per node
    String m_dofs_multi_on_node_family_name; // Several dofs per node (non constant size)
};

#endif //PASSMO_NODEDOFS_H
