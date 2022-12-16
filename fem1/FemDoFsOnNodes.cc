// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemDoFsOnNodes.cc                                           (C) 2022-2022 */
/*                                                                           */
/* Utilitary classes for FEM.                                                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemDoFsOnNodes.h"

#include <arcane/mesh/DoFFamily.h>
#include "arcane/IIndexedIncrementalItemConnectivityMng.h"
#include "arcane/IIndexedIncrementalItemConnectivity.h"
#include "arcane/IndexedItemConnectivityView.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class FemDoFsOnNodes::Impl
: public TraceAccessor
{
 public:

  Impl(ITraceMng* tm)
  : TraceAccessor(tm)
  {
  }

 public:

  void initialize(IMesh* mesh, Int32 nb_dof_per_node);

 public:

  Ref<IIndexedIncrementalItemConnectivity> m_node_dof_connectivity;
  IItemFamily* m_dof_family = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FemDoFsOnNodes::
FemDoFsOnNodes(ITraceMng* tm)
: m_p(new Impl(tm))
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FemDoFsOnNodes::
~FemDoFsOnNodes()
{
  delete m_p;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemDoFsOnNodes::Impl::
initialize(IMesh* mesh, Int32 nb_dof_per_node)
{
  IItemFamily* dof_family_interface = mesh->findItemFamily(Arcane::IK_DoF, "DoFNodeFamily", true);
  mesh::DoFFamily* dof_family = ARCANE_CHECK_POINTER(dynamic_cast<mesh::DoFFamily*>(dof_family_interface));
  m_dof_family = dof_family_interface;

  // Create the DoFs
  Int64UniqueArray uids(mesh->allNodes().size() * nb_dof_per_node);
  Int64 max_node_uid = mesh::DoFUids::getMaxItemUid(mesh->nodeFamily());
  {
    Integer dof_index = 0;
    ENUMERATE_NODE (inode, mesh->allNodes()) {
      Node node = *inode;
      Int64 node_unique_id = node.uniqueId().asInt64();
      for (Integer i = 0; i < nb_dof_per_node; ++i) {
        uids[dof_index] = node_unique_id * nb_dof_per_node + i;
        ++dof_index;
      }
    }
  }
  info() << "ADD_Dofs list=" << uids;
  Int32UniqueArray dof_lids(uids.size());
  dof_family->addDoFs(uids, dof_lids);
  dof_family->endUpdate();
  info() << "NB_DOF=" << dof_family->allItems().size();

  // Create Node -> DoF connectivity.
  m_node_dof_connectivity = mesh->indexedConnectivityMng()->findOrCreateConnectivity(mesh->nodeFamily(), dof_family, "DoFNode");
  auto* cn = m_node_dof_connectivity->connectivity();
  {
    Integer dof_index = 0;
    ENUMERATE_NODE (inode, mesh->allNodes()) {
      NodeLocalId node = *inode;
      for (Integer i = 0; i < nb_dof_per_node; ++i) {
        cn->addConnectedItem(node, DoFLocalId(dof_lids[dof_index]));
        ++dof_index;
      }
    }
  }
  info() << "End build Dofs";

  IndexedNodeDoFConnectivityView node_dof(m_node_dof_connectivity->view());
  {
    // Set the owners of the DoF.
    IParallelMng* pm = mesh->parallelMng();
    Int32 my_rank = pm->commRank();
    ItemInternalList dofs = m_dof_family->itemsInternal();
    ENUMERATE_ (Node, inode, mesh->allNodes()) {
      Node node = *inode;
      Int32 node_owner = node.owner();
      for (DoFLocalId dof : node_dof.dofs(node)) {
        dofs[dof]->setOwner(node_owner, my_rank);
      }
    }
    dof_family->notifyItemsOwnerChanged();
    dof_family->computeSynchronizeInfos();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemDoFsOnNodes::
initialize(IMesh* mesh, Int32 nb_dof_per_node)
{
  m_p->initialize(mesh, nb_dof_per_node);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

IndexedNodeDoFConnectivityView FemDoFsOnNodes::
nodeDoFConnectivityView() const
{
  return m_p->m_node_dof_connectivity->view();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

IItemFamily* FemDoFsOnNodes::
dofFamily() const
{
  return m_p->m_dof_family;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
