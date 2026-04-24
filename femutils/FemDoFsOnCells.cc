// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemDoFsOnCells.cc                                           (C) 2022-2026 */
/*                                                                           */
/* Utilitary classes for FEM.                                                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemDoFsOnCells.h"
#include <arcane/core/MeshUtils.h>
#include <arcane/mesh/DoFFamily.h>
#include "arcane/IIndexedIncrementalItemConnectivityMng.h"
#include "arcane/IIndexedIncrementalItemConnectivity.h"
#include "arcane/IndexedItemConnectivityView.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class FemDoFsOnCells::Impl
: public TraceAccessor
{
 public:

  Impl(ITraceMng* tm)
  : TraceAccessor(tm)
  {
  }

 public:

  void initialize(IMesh* mesh, Int32 nb_dof_per_cell);

 public:

  Ref<IIndexedIncrementalItemConnectivity> m_cell_dof_connectivity;
  IItemFamily* m_dof_family = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FemDoFsOnCells::
FemDoFsOnCells(ITraceMng* tm)
: m_p(new Impl(tm))
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FemDoFsOnCells::
~FemDoFsOnCells()
{
  delete m_p;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemDoFsOnCells::Impl::
initialize(IMesh* mesh, Int32 nb_dof_per_cell)
{
  IItemFamily* dof_family_interface = mesh->findItemFamily(Arcane::IK_DoF, "DoFCellFamily", true);
  mesh::DoFFamily* dof_family = ARCANE_CHECK_POINTER(dynamic_cast<mesh::DoFFamily*>(dof_family_interface));
  m_dof_family = dof_family_interface;

  // Create the DoFs
  Int64UniqueArray uids(mesh->allCells().size() * nb_dof_per_cell);
  {
    Integer dof_index = 0;
    // Use a mask to make sure the uniqueId() of the dof
    // can not be negative if we multiply the uniqueId().
    const UInt64 uid_mask = (1<<28) - 1;
    ENUMERATE_CELL (icell, mesh->allCells()) {
      Cell cell = *icell;
      Int64 cell_unique_id = cell.uniqueId().asInt64();
      for (Integer i = 0; i < nb_dof_per_cell; ++i) {
        uids[dof_index] = (cell_unique_id & uid_mask) * nb_dof_per_cell + i;
        ++dof_index;
      }
    }
  }
  //info() << "ADD_Dofs list=" << uids;
  Int32UniqueArray dof_lids(uids.size());
  dof_family->addDoFs(uids, dof_lids);
  dof_family->endUpdate();
  info() << "NB_DOF=" << dof_family->allItems().size();

  // Create Cell -> DoF connectivity.
  m_cell_dof_connectivity = mesh->indexedConnectivityMng()->findOrCreateConnectivity(mesh->cellFamily(), dof_family, "DoFCell");
  auto* cn = m_cell_dof_connectivity->connectivity();
  {
    Integer dof_index = 0;
    ENUMERATE_CELL (icell, mesh->allCells()) {
      CellLocalId cell = *icell;
      for (Integer i = 0; i < nb_dof_per_cell; ++i) {
        cn->addConnectedItem(cell, DoFLocalId(dof_lids[dof_index]));
        ++dof_index;
      }
    }
  }
  info() << "End build Dofs";

  IndexedCellDoFConnectivityView cell_dof(m_cell_dof_connectivity->view());
  {
    // Set the owners of the DoF.
    IParallelMng* pm = mesh->parallelMng();
    Int32 my_rank = pm->commRank();
    ItemInternalList dofs = m_dof_family->itemsInternal();
    ENUMERATE_ (Cell, icell, mesh->allCells()) {
      Cell cell = *icell;
      Int32 cell_owner = cell.owner();
      for (DoFLocalId dof : cell_dof.dofs(cell)) {
        dofs[dof]->setOwner(cell_owner, my_rank);
      }
    }
    dof_family->notifyItemsOwnerChanged();
    dof_family->computeSynchronizeInfos();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemDoFsOnCells::
initialize(IMesh* mesh, Int32 nb_dof_per_cell)
{
  m_p->initialize(mesh, nb_dof_per_cell);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

IndexedCellDoFConnectivityView FemDoFsOnCells::
cellDoFConnectivityView() const
{
  return m_p->m_cell_dof_connectivity->view();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

IItemFamily* FemDoFsOnCells::
dofFamily() const
{
  return m_p->m_dof_family;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/