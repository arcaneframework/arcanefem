//
// Created by EF on 12/12/22.
//
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
#include "NodeDoFs_axl.h"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
extern Integer NODE_NDDL;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class NodeDoFs:
        public ArcaneNodeDoFsObject {

public:

    NodeDoFs(const ServiceBuildInfo& sbi)
            : ArcaneNodeDoFsObject(sbi)
            , m_connectivity_mng(sbi.subDomain()->traceMng())
            , m_dof_mng(sbi.mesh(),&m_connectivity_mng)
            , m_dof_family_name("DoFFamily")
            , m_dof_on_node_family_name("DoFOnNode")
            , m_dofs_on_node_family_name("DoFsOnNode")
            , m_dofs_multi_on_node_family_name("DoFsMultiOnNode"){}

public:
    typedef ItemConnectivityT          <Node, DoF> NodeToDoFConnectivity;
    typedef ItemArrayConnectivityT     <Node, DoF> NodeToDoFsConnectivity;
    typedef ItemMultiArrayConnectivityT<Node, DoF> NodeToDoFsMultiConnectivity;

    typedef SharedArray2<Int32> Int32SharedArray2;
    typedef SharedArray2<Int64> Int64SharedArray2;

public:

    DoFManager& dofMng() {return m_dof_mng;}
    mesh::DoFFamily& getNodeDof() { return m_dof_mng.family(m_dof_on_node_family_name); }
    mesh::DoFFamily& getNodeDofs() { return m_dof_mng.family(m_dofs_on_node_family_name); }
    mesh::DoFFamily& getMultiNodeDofs() { return m_dof_mng.family(m_dofs_multi_on_node_family_name); }
    Integer addDoF(const Integer& begin_uid);
    Integer removeDoF(const Integer& begin_lid);
    Integer addDoFs(const Integer& begin_uid, const Integer& size);
    Integer removeDoFs(const Integer& begin_lid, const Integer& size);
    DoFGroup createDoFGroup(const Integer& begin_lid, const Integer& size);
    VariableDoFReal doFVariable();
    VariableDoFArrayReal doFVectorVariable(const Integer& size);
    void doFConnectivity();

private:

    // Connectivity tests
    void _node2DoFConnectivity();
    void _node2DoFConnectivityRegistered();
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Integer NodeDoFs::
addDoF(const Integer& begin_uid)
{
    // Create DoF
    Int32 size = 1;
    Int64UniqueArray dof_uids(size);
    for (Integer i = 0; i < size; ++i) dof_uids[i] = i+begin_uid;

    mesh::DoFFamily& dof_family = getNodeDof();
    info() << "=== add dof to family " << dof_family.name();

    Int32UniqueArray dof_lids(size);
    dof_family.addDoFs(dof_uids,dof_lids);
    dof_family.endUpdate();
    return (begin_uid+size);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Integer NodeDoFs::
addDoFs(const Integer& begin_uid, const Integer& size)
{
    // Create DoFs
    Int64UniqueArray dof_uids(size);
    for (Integer i = 0; i < size; ++i) dof_uids[i] = i+begin_uid;

    mesh::DoFFamily& dof_family = getNodeDofs();
    info() << "=== add dofs to family " << dof_family.name();

    Int32UniqueArray dof_lids(size);
    dof_family.addDoFs(dof_uids,dof_lids);
    dof_family.endUpdate();
    return (begin_uid+size);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Integer
NodeDoFs::
removeDoF(const Integer& begin_lid)
{
    mesh::DoFFamily& dof_family = getNodeDof();
    info() << "=== remove dof from family " << dof_family.name();

    Int32ConstArrayView dof_lids = dof_family.view().localIds();
    dof_family.removeDoFs(dof_lids.subConstView(begin_lid,1));
    dof_family.endUpdate();

    return dof_family.nbItem();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Integer
NodeDoFs::
removeDoFs(const Integer& begin_lid, const Integer& size)
{
    mesh::DoFFamily& dof_family = getNodeDofs();
    info() << "=== remove dofs from family " << dof_family.name();

    Int32ConstArrayView dof_lids = dof_family.view().localIds();
    dof_family.removeDoFs(dof_lids.subConstView(begin_lid,size));
    dof_family.endUpdate();

    // Enumerate DoF
/*    Int64 dof_id;
    ENUMERATE_DOF(idof,dof_family.view())
    {
        dof_id = idof->uniqueId().asInt64();
        info() << "= Dof id : " << dof_id;
    }*/
    return dof_family.nbItem();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

DoFGroup
NodeDoFs::
createDoFGroup(const Integer& begin_lid, const Integer& size)
{
    mesh::DoFFamily& dof_family = getNodeDof();

    // CREATE NEW GROUP; try ArcGeoSim ItemGroupBuilder
    Int32ConstArrayView lids = dof_family.view().localIds().subConstView(begin_lid,size);
    String dof_group_name("DoFGroup");
    return dof_family.createGroup(dof_group_name,lids);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void
NodeDoFs::
_node2DoFsConnectivityRegistered()
{
    info() << "====================================================";
    info() << "== REGISTER A CONNECTIVITY (ARRAY) AND FOLLOW MESH EVOLUTION ";
    info() << "====================================================";

    // Done for node since it's easier to add new items
    mesh::DoFFamily& dofs_on_node_family = getNodeDof();
    Integer nb_dof_per_node = NODE_NDDL;
    // Create the DoFs
    Int64UniqueArray uids(ownNodes().size()*nb_dof_per_node);
    Int64 max_node_uid = mesh::DoFUids::getMaxItemUid(mesh()->nodeFamily());
    Int64 max_dof_uid  = mesh::DoFUids::getMaxItemUid(&dofs_on_node_family);
    Integer j = 0;
    ENUMERATE_NODE(inode,ownNodes()) {
        for (Integer i = 0; i < nb_dof_per_node; ++i)
            uids[j++] = mesh::DoFUids::uid(max_dof_uid,max_node_uid,inode->uniqueId().asInt64(),i);
    }
    Int32UniqueArray lids(uids.size());
    dofs_on_node_family.addDoFs(uids,lids);
    dofs_on_node_family.endUpdate();

    IItemFamily* node_family = mesh()->nodeFamily();

    NodeToDoFsConnectivity node2dofs(node_family,&dofs_on_node_family,nb_dof_per_node,"NodeToDoFs");

    // Create the ghosts
    GhostLayerFromConnectivityComputer ghost_builder(&node2dofs);
    IItemConnectivitySynchronizer* synchronizer = dofMng().connectivityMng()->createSynchronizer(&node2dofs,&ghost_builder);
    synchronizer->synchronize();

    // Save your connectivity
    dofMng().connectivityMng()->registerConnectivity(&node2dofs);

    // Local mesh change: add own and ghost nodes
    Integer nb_subdomain = subDomain()->parallelMng()->commSize();
    Integer nb_new_nodes = 2;
    Int64UniqueArray2 new_nodes_uids(nb_subdomain,nb_new_nodes);
    Int32UniqueArray2 new_nodes_lids(nb_subdomain,nb_new_nodes);
    _addNodes(new_nodes_lids,new_nodes_uids);

    node_family->endUpdate(); // Connectivity is updated in this call
    node_family->computeSynchronizeInfos(); // Connectivity ghosts are updated in this call

    // Follow mesh changes in connectivity
    // Add and remove in FromFamily
    IItemConnectivityMng* connectivity_mng = dofMng().connectivityMng();
    if (! connectivity_mng->isUpToDate(&node2dofs))
    {
        // Handle added nodes : create nb_dof_per_node dofs for each own node added
        Int32ArrayView source_family_added_items_lids;
        Int32ArrayView source_family_removed_items_lids;
        connectivity_mng->getSourceFamilyModifiedItems(&node2dofs,source_family_added_items_lids,source_family_removed_items_lids);
        ItemVector source_family_added_items_own(node_family);
        ENUMERATE_NODE(inode,node_family->view(source_family_added_items_lids)) if (inode->isOwn()) source_family_added_items_own.add(inode.localId());
        // Create new dofs on these new nodes : on the owned node only
        Int64UniqueArray uids(source_family_added_items_own.size()*nb_dof_per_node);
        Integer j = 0;
        Int32SharedArray source_family_lids_in_connectivity(source_family_added_items_own.size()*nb_dof_per_node);
        Int64 max_item_uid= mesh::DoFUids::getMaxItemUid(node_family);
        Int64 max_dof_uid = mesh::DoFUids::getMaxItemUid(&dofs_on_node_family);
        ENUMERATE_NODE(inode,source_family_added_items_own)
        {
            for (Integer i = 0; i < nb_dof_per_node; ++i)
            {
                uids[j] = mesh::DoFUids::uid(max_dof_uid,max_item_uid,inode->uniqueId().asInt64(),i);
                source_family_lids_in_connectivity[j++] = inode.localId(); // Replicate the from item lid each time it is used in a connectivity (needed to use updateConnectivity)
            }
        }
        Int32SharedArray lids(uids.size());
        dofs_on_node_family.addDoFs(uids,lids);
        dofs_on_node_family.endUpdate();
        // Update connectivity
        node2dofs.updateConnectivity(source_family_lids_in_connectivity,lids);
        // Update ghost
        synchronizer->synchronize();
        connectivity_mng->setUpToDate(&node2dofs);
    }

    // Remove the added nodes
    Integer nb_removed_nodes = 1;
    Int32SharedArray2 removed_nodes_lids(nb_subdomain,nb_removed_nodes);
    Int64SharedArray2 removed_nodes_uids(nb_subdomain,nb_removed_nodes);
    Integer nb_remaining_nodes = new_nodes_lids.dim2Size()-nb_removed_nodes;
    Int32SharedArray2 remaining_nodes_lids(nb_subdomain,nb_remaining_nodes);
    Int64SharedArray2 remaining_nodes_uids(nb_subdomain,nb_remaining_nodes);
    _removeNodes(new_nodes_lids,nb_removed_nodes,removed_nodes_lids,removed_nodes_uids,remaining_nodes_lids,remaining_nodes_uids);
    node_family->endUpdate(); // Connectivity and ghosts are updated (since own and ghost dof are removed)
    node_family->computeSynchronizeInfos(); // Not needed by connectivity but needed to have NodeFamily synchronization info up to date

    // Update connectivity : set the removed nodes to Null item lid
    if(!connectivity_mng->isUpToDate(&node2dofs))
    {
        Int32ArrayView source_family_added_items_lids;
        Int32ArrayView source_family_removed_items_lids;
        connectivity_mng->getSourceFamilyModifiedItems(&node2dofs,source_family_added_items_lids,source_family_removed_items_lids);
        // Get dof connected to removed nodes, to remove them (used to test dof family compaction)
        Integer nb_removed_dofs = nb_removed_nodes*nb_dof_per_node;
        Int32UniqueArray removed_dofs;
        removed_dofs.reserve(nb_removed_dofs);
        ItemInternal internal;
        ConnectivityItemVector con(node2dofs);
        for (Integer i = 0; i < nb_removed_nodes; ++i)
        {
            internal.setLocalId(source_family_removed_items_lids[i]);
            Node node(&internal);
            removed_dofs.addRange(node2dofs(node,con).localIds());
        }
        // Prepare data to update connectivity
        Integer nb_connections = source_family_removed_items_lids.size()*nb_dof_per_node;
        Int32SharedArray source_family_removed_items_lids_in_connectivity(nb_connections);
        for (Integer i = 0; i < source_family_removed_items_lids.size(); ++i)
        {
            for(Integer j = 0; j< nb_dof_per_node;++j) source_family_removed_items_lids_in_connectivity[i*nb_dof_per_node+j] = source_family_removed_items_lids[i];
        }
        Int32SharedArray null_item_lids(nb_connections,NULL_ITEM_LOCAL_ID);
        // Update connectivity
        node2dofs.updateConnectivity(source_family_removed_items_lids_in_connectivity,null_item_lids);
        connectivity_mng->setUpToDate(&node2dofs);
        // Unused dof can be removed if desired. Usefull to remove them to test compaction
        dofs_on_node_family.removeDoFs(removed_dofs);
        dofs_on_node_family.endUpdate();
        debug() << "*** REMOVED DOFS " << removed_dofs;

    }

    // Compact Source family
    node_family->compactItems(true);
    debug() << "NODES " << Int32UniqueArray(node_family->view().localIds());
    // update lids
    for (Arcane::Integer rank = 0; rank < nb_subdomain; ++rank)
    {
        node_family->itemsUniqueIdToLocalId(remaining_nodes_lids[rank],remaining_nodes_uids[rank],true);
    }

    // Compact Target family
    debug() << "DOFS " << dofs_on_node_family.view().localIds();
    debug() << "DOF family size " << dofs_on_node_family.nbItem();
    dofs_on_node_family.compactItems(true);
    debug() << "DOF family size " << dofs_on_node_family.nbItem();
    debug() << "DOFS " << dofs_on_node_family.view().localIds();

    dofMng().connectivityMng()->unregisterConnectivity(&node2dofs);

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void
NodeDoFs::
_node2DoFsMultiConnectivityRegistered()
{
    info() << "====================================================";
    info() << "== REGISTER A CONNECTIVITY (MULTI ARRAY) AND FOLLOW MESH EVOLUTION ";
    info() << "====================================================";
    // Done for node since it's easier to add new items
    mesh::DoFFamily& dofs_multi_on_node_family = getMultiNodeDofs();
    IItemFamily* node_family = mesh()->nodeFamily();
    IntegerUniqueArray nb_dof_per_node(node_family->maxLocalId(),0);
    // Create the DoFs
    Int64UniqueArray uids;
    Int64 max_node_uid = mesh::DoFUids::getMaxItemUid(mesh()->nodeFamily());
    Int64 max_dof_uid =  mesh::DoFUids::getMaxItemUid(&dofs_multi_on_node_family);
    ENUMERATE_NODE(inode,ownNodes()) {
        nb_dof_per_node[inode->localId()] = NODE_NDDL; // constant size in initialization
        for (Integer i = 0; i < nb_dof_per_node[inode->localId()]; ++i)
            uids.add(mesh::DoFUids::uid(max_dof_uid,max_node_uid,inode->uniqueId().asInt64(),i));
    }
    Int32UniqueArray lids(uids.size());
    dofs_multi_on_node_family.addDoFs(uids,lids);
    dofs_multi_on_node_family.endUpdate();

    NodeToDoFsMultiConnectivity node2dofs(node_family,&dofs_multi_on_node_family,nb_dof_per_node,"NodeToDoFsMulti");

    // Create the ghosts
    GhostLayerFromConnectivityComputer ghost_builder(&node2dofs);
    IItemConnectivitySynchronizer* synchronizer = dofMng().connectivityMng()->createSynchronizer(&node2dofs,&ghost_builder);
    synchronizer->synchronize();

    // Save your connectivity
    dofMng().connectivityMng()->registerConnectivity(&node2dofs);

    // Local mesh change: add own and ghost nodes
    Integer nb_subdomain = subDomain()->parallelMng()->commSize();
    Integer nb_new_nodes = 2;
    Int64UniqueArray2 new_nodes_uids(nb_subdomain,nb_new_nodes);
    Int32UniqueArray2 new_nodes_lids(nb_subdomain,nb_new_nodes);
    _addNodes(new_nodes_lids,new_nodes_uids);

    node_family->endUpdate(); // Connectivity is updated in this call
    node_family->computeSynchronizeInfos(); // Connectivity ghosts are updated in this call

    // Follow mesh changes in connectivity.
    IItemConnectivityMng* connectivity_mng = dofMng().connectivityMng();
    IntegerUniqueArray nb_dof_per_new_node;

    if (! connectivity_mng->isUpToDate(&node2dofs))
    {
        // Handle added nodes : create a variable number of dofs for each own node added
        Int32ArrayView source_family_added_items_lids;
        Int32ArrayView source_family_removed_items_lids;
        connectivity_mng->getSourceFamilyModifiedItems(&node2dofs,source_family_added_items_lids,source_family_removed_items_lids);
        ItemVector source_family_added_items_own(node_family);
        ENUMERATE_NODE(inode,node_family->view(source_family_added_items_lids)) if (inode->isOwn()) source_family_added_items_own.add(inode.localId());
        IntegerUniqueArray source_family_added_items_own_lids(source_family_added_items_own.viewAsArray());
        nb_dof_per_new_node.resize(source_family_added_items_own_lids.size());
        Integer nb_new_dofs = 0;
        for (Arcane::Integer i = 0; i < nb_dof_per_new_node.size(); ++i)
        {
            nb_dof_per_new_node[i] = i+1;
            nb_new_dofs += nb_dof_per_new_node[i];
        }
        Integer nb_connections  = nb_new_dofs;
        // Create new dofs on these new nodes : on the owned node only
        Int64UniqueArray uids(nb_connections);
        Integer j = 0;
        Int32SharedArray source_family_lids_in_connectivity(nb_connections);
        Int64 max_item_uid= mesh::DoFUids::getMaxItemUid(node_family);
        Int64 max_dof_uid = mesh::DoFUids::getMaxItemUid(&dofs_multi_on_node_family);
        ENUMERATE_NODE(inode,source_family_added_items_own)
        {
            for (Integer i = 0; i < nb_dof_per_new_node[inode.index()]; ++i)
            {
                uids[j] = mesh::DoFUids::uid(max_dof_uid,max_item_uid,inode->uniqueId().asInt64(),i);
                source_family_lids_in_connectivity[j++] = inode.localId(); // Replicate the from item lid each time it used in a connectivity (needed to use updateConnectivity)
            }
        }
        Int32SharedArray lids(uids.size());
        dofs_multi_on_node_family.addDoFs(uids,lids);
        dofs_multi_on_node_family.endUpdate();
        // Update connectivity
        node2dofs.updateConnectivity(source_family_lids_in_connectivity,lids);
        // Update ghost
        synchronizer->synchronize();
        connectivity_mng->setUpToDate(&node2dofs);
    }

    debug() << "NEW NODE LIDS " << new_nodes_lids[mesh()->parallelMng()->commRank()];

    // Remove the added nodes
    Integer nb_removed_nodes = 1;
    Int32SharedArray2 removed_nodes_lids(nb_subdomain,nb_removed_nodes);
    Int64SharedArray2 removed_nodes_uids(nb_subdomain,nb_removed_nodes);
    Integer nb_remaining_nodes = new_nodes_lids.dim2Size()-nb_removed_nodes;
    Int32SharedArray2 remaining_nodes_lids(nb_subdomain,nb_remaining_nodes);
    Int64SharedArray2 remaining_nodes_uids(nb_subdomain,nb_remaining_nodes);
    _removeNodes(new_nodes_lids, nb_removed_nodes,removed_nodes_lids,removed_nodes_uids, remaining_nodes_lids,remaining_nodes_uids);
    debug() << "REMOVED NODES " << removed_nodes_lids[mesh()->parallelMng()->commRank()] << removed_nodes_uids[mesh()->parallelMng()->commRank()];

    node_family->endUpdate(); // Connectivity and ghosts are updated (since own and ghost dof are removed)
    node_family->computeSynchronizeInfos(); // Not needed by connectivity but needed to have NodeFamily synchronization info up to date

    debug() << "NODES " << Int32UniqueArray(node_family->view().localIds());

    // Update connectivity : set the removed nodes to Null item lid
    if(!connectivity_mng->isUpToDate(&node2dofs))
    {
        Int32ArrayView source_family_added_items_lids;
        Int32ArrayView source_family_removed_items_lids;
        connectivity_mng->getSourceFamilyModifiedItems(&node2dofs,source_family_added_items_lids,source_family_removed_items_lids);
        // Prepare data to update connectivity
        Int32SharedArray source_family_removed_items_lids_in_connectivity;
        nb_dof_per_new_node = node2dofs.itemProperty().dim2Sizes();
        Integer nb_removed_dofs = 0;
        for (Integer i = 0; i < source_family_removed_items_lids.size(); ++i)
        {
            Int32 removed_item_lid = source_family_removed_items_lids[i];
            for(Integer j = 0; j< nb_dof_per_new_node[removed_item_lid];++j)
            {
                source_family_removed_items_lids_in_connectivity.add(source_family_removed_items_lids[i]);
                nb_removed_dofs++;
            }
        }
        // Get dof connected to removed nodes, to remove them (used to test dof family compaction)
        Int32UniqueArray removed_dofs;
        removed_dofs.reserve(nb_removed_dofs);
        ItemInternal internal;
        ConnectivityItemVector con(node2dofs);
        for (Integer i = 0; i < nb_removed_nodes; ++i)
        {
            internal.setLocalId(source_family_removed_items_lids[i]);
            Node node(&internal);
            removed_dofs.addRange(node2dofs(node,con).localIds());
        }
        Int32SharedArray null_item_lids(source_family_removed_items_lids_in_connectivity.size(),NULL_ITEM_LOCAL_ID);
        // Update connectivity
        node2dofs.updateConnectivity(source_family_removed_items_lids_in_connectivity,null_item_lids);
        // unused dof can be removed if desired. Usefull to remove them to test compaction
        dofs_multi_on_node_family.removeDoFs(removed_dofs);
        dofs_multi_on_node_family.endUpdate();
        connectivity_mng->setUpToDate(&node2dofs);
    }

    // Compact Source family
    node_family->compactItems(true);
    debug() << "NODES " << node_family->view().localIds();
    // update lids
    for (Arcane::Integer rank = 0; rank < nb_subdomain; ++rank)
    {
        node_family->itemsUniqueIdToLocalId(remaining_nodes_lids[rank],remaining_nodes_uids[rank],true);
    }

    // Compact Target family
    debug() << "DOFS " << dofs_multi_on_node_family.view().localIds();
    debug() << "DOF family size " << dofs_multi_on_node_family.nbItem();
    dofs_multi_on_node_family.compactItems(true);
    debug() << "DOF family size " << dofs_multi_on_node_family.nbItem();
    debug() << "DOFS " << dofs_multi_on_node_family.view().localIds();
    dofMng().connectivityMng()->unregisterConnectivity(&node2dofs);

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void
NodeDoFs::
_removeGhost(mesh::DoFFamily& dof_family)
{
    Int32SharedArray removed_items;
    ENUMERATE_DOF(idof,dof_family.ghostDoFs())
    {
        removed_items.add(idof->localId());
    }
    dof_family.removeDoFs(removed_items);
    dof_family.endUpdate();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void
NodeDoFs::
_addNodes(Int32Array2View new_nodes_lids, Int64Array2View new_nodes_uids)
{
    // Change mesh by removing and adding nodes
    Integer nb_new_nodes = new_nodes_lids.dim2Size();
    Integer nb_subdomain = subDomain()->parallelMng()->commSize();
    Integer local_rank   = subDomain()->parallelMng()->commRank();
    Int64 max_uid = mesh::DoFUids::getMaxItemUid(mesh()->nodeFamily());
    IMeshModifier* mesh_modifier = mesh()->modifier();
    UniqueArray<ItemVectorView> added_items(nb_subdomain);
    // Add node and ghost nodes on each subdomain. Each subdomain has all the nodes of the other as ghosts (just for the demo)
    for (Integer rank = 0; rank < nb_subdomain; ++rank)
    {
        for (Integer i = 0; i < nb_new_nodes;++i) new_nodes_uids[rank][i] = max_uid*(rank+1) +i+1;
        mesh_modifier->addNodes(new_nodes_uids[rank],new_nodes_lids[rank]);
        added_items[rank] = mesh()->nodeFamily()->view(new_nodes_lids[rank]);
        ENUMERATE_NODE(inode,added_items[rank])
        {
            inode->internal()->setOwner(rank,local_rank);
            info() << "== Add item " << inode->localId() << " on rank " << local_rank << " with owner " << rank;
        }
    }

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void
NodeDoFs::
_removeNodes(Int32ConstArray2View new_nodes_lids,
             const Integer nb_removed_nodes,
             Int32SharedArray2& removed_node_lids,
             Int64SharedArray2& removed_node_uids,
             Int32SharedArray2& remaining_node_lids,
             Int64SharedArray2& remaining_node_uids)
{
    ARCANE_ASSERT((nb_removed_nodes <= new_nodes_lids.dim2Size()),("Cannot removed more nodes than available"));
    mesh::NodeFamily* node_family = static_cast<mesh::NodeFamily*> (mesh()->nodeFamily());
    Integer nb_subdomain = subDomain()->parallelMng()->commSize();
    Integer nb_remaining_nodes = remaining_node_lids.dim2Size();
    UniqueArray<ItemVectorView> removed_nodes(nb_subdomain);
    UniqueArray<ItemVectorView> remaining_nodes(nb_subdomain);
    for (Integer rank = 0; rank < nb_subdomain; ++rank)
    {
        removed_node_lids[rank].copy(new_nodes_lids[rank].subConstView(0,nb_removed_nodes));
        removed_nodes[rank] = node_family->view(removed_node_lids[rank]);
        remaining_node_lids[rank].copy(new_nodes_lids[rank].subConstView(nb_removed_nodes,nb_remaining_nodes));
        remaining_nodes[rank] = node_family->view(remaining_node_lids[rank]);
        Int32 i = 0;
        ENUMERATE_NODE(inode,removed_nodes[rank])
        {
            node_family->removeNodeIfNotConnected(inode->internal());
            removed_node_uids[rank][i++] = inode->uniqueId().asInt64();
        }
        i = 0;
        ENUMERATE_NODE(inode,remaining_nodes[rank])
        {
            remaining_node_uids[rank][i++] = inode->uniqueId().asInt64();
        }
    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

VariableDoFReal
NodeDoFs::
doFVariable()
{
    mesh::DoFFamily& dof_family = getNodeDofs();

    // Create and fillDoF Scalar Variable (with internal support)
    info() << "=== CreateDoFScalarVariable : ";
    VariableDoFReal dof_variable(VariableBuildInfo(mesh(),"DoFScalarVariable",m_dofs_on_node_family_name));

    dof_variable.fill(0.);
    ENUMERATE_DOF(idof,dof_family.ownDoFs())
    {
        dof_variable[idof] = (Real)(idof->uniqueId().asInt64());
        info() << "= dof_variable[idof] = " << dof_variable[idof];
    }
    dof_variable.synchronize();
    return dof_variable;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
VariableDoFArrayReal
NodeDoFs::
doFVectorVariable(const Integer& size)
{
    mesh::DoFFamily& dof_family = getNodeDofs();

    // Create and fillDoF Array Variable (with internal support)
    info() << "=== CreateDoFVariable : ";
    VariableDoFArrayReal dof_array_variable(VariableBuildInfo(mesh(),"DoFArrayVariable",m_dofs_on_node_family_name));

    dof_array_variable.resize(size);
    // Initialize with -1
    ENUMERATE_DOF(idof,dof_family.allDoFs()) {dof_array_variable[idof].fill(-1);}

    ENUMERATE_DOF(idof,dof_family.ownDoFs())
    {
        for (Integer i = 0; i < size; ++i)
        {
            dof_array_variable[idof][i] = (Real)((i+1)*idof->uniqueId().asInt64());
            info() << "= dof_array_variable[idof] = " << dof_array_variable[idof][i];
        }
    }
    dof_array_variable.synchronize();
    return dof_array_variable;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_NODEDOFS(NodeDoFs,NodeDoFs);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
