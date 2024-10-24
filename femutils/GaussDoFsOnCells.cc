// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/* GaussDoFsOnCells.cc                                         (C) 2022-2024 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "GaussDoFsOnCells.h"
#include <arcane/mesh/DoFFamily.h>
#include "arcane/IIndexedIncrementalItemConnectivityMng.h"
#include "arcane/IIndexedIncrementalItemConnectivity.h"
#include "arcane/IndexedItemConnectivityView.h"
#include <arcane/VariableTypes.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::mesh;
//using namespace Arcane::FemUtils;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
class GaussDoFsOnCells::Impl
: public TraceAccessor
{
 public:

  explicit Impl(ITraceMng* tm)
  : TraceAccessor(tm)
  {
  }

 public:

  void initialize(IMesh* mesh, Int32 max_nb_gauss_per_cell);

 public:

  Ref<IIndexedIncrementalItemConnectivity> m_cell_gauss_connectivity;
  IItemFamily* m_gauss_family = nullptr; //treated as DoFFamily

  VariableDoFArrayReal* m_gauss_shape = nullptr;
  VariableDoFArrayReal3* m_gauss_shapederiv = nullptr;
  VariableDoFReal3* m_gauss_refpos = nullptr;
  VariableDoFReal3x3* m_gauss_jacobmat = nullptr;
  VariableDoFReal* m_gauss_weight = nullptr;
  VariableDoFReal* m_gauss_jacobian = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

GaussDoFsOnCells::
GaussDoFsOnCells(ITraceMng* tm)
: m_p(new Impl(tm))
{}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

GaussDoFsOnCells::
~GaussDoFsOnCells()
{
  delete m_p->m_gauss_refpos;
  delete m_p->m_gauss_weight;
  delete m_p->m_gauss_jacobian;
  delete m_p->m_gauss_jacobmat;
  delete m_p->m_gauss_shape;
  delete m_p->m_gauss_shapederiv;
  delete m_p;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void GaussDoFsOnCells::Impl::
initialize(IMesh* mesh, Int32 max_nb_gauss_per_cell)
{
  IItemFamily* dof_family_interface = mesh->findItemFamily(Arcane::IK_DoF, "GaussCellFamily", true);
  mesh::DoFFamily* dof_family = ARCANE_CHECK_POINTER(dynamic_cast<mesh::DoFFamily*>(dof_family_interface));
  m_gauss_family = dof_family_interface;

  // Create the Gauss points as Arcane "DoFs" attached to cells
  Int64UniqueArray uids(mesh->allCells().size() * max_nb_gauss_per_cell);
  Int64 max_cell_uid = mesh::DoFUids::getMaxItemUid(mesh->cellFamily());
  {
    Integer gauss_index{ 0 };
    ENUMERATE_CELL (icell, mesh->allCells()) {
      Cell cell = *icell;
      auto cell_type = cell.type();
      Int64 cell_unique_id = cell.uniqueId().asInt64();
      for (Integer i = 0; i < max_nb_gauss_per_cell; ++i) {
        uids[gauss_index++] = cell_unique_id * max_nb_gauss_per_cell + i;
      }
    }
  }

  Integer uidsize = uids.size();
  Int32UniqueArray gauss_lids(uidsize);
  dof_family->addDoFs(uids, gauss_lids);
  dof_family->endUpdate();
  Integer nb_gauss = dof_family->allItems().size();
  info() << "NB_GAUSS=" << nb_gauss << " max_per_cell=" << max_nb_gauss_per_cell;

  // Create Cell -> Gauss (DoF) connectivity.
  m_cell_gauss_connectivity = mesh->indexedConnectivityMng()->findOrCreateConnectivity(mesh->cellFamily(), dof_family, "GaussCell");
  auto* cn = m_cell_gauss_connectivity->connectivity();
  {
    Integer gauss_index{ 0 };
    ENUMERATE_CELL (icell, mesh->allCells()) {
      CellLocalId cell = *icell;
      for (Integer i = 0; i < max_nb_gauss_per_cell; ++i) {
        cn->addConnectedItem(cell, DoFLocalId(gauss_lids[gauss_index++]));
      }
    }
  }
  info() << "End build Gauss points";

  IndexedCellDoFConnectivityView cell_gauss(m_cell_gauss_connectivity->view());
  {
    // Set the owners of the Gauss point (=DoF).
    IParallelMng* pm = mesh->parallelMng();
    Int32 my_rank = pm->commRank();
    ItemInternalList gauss_dofs = m_gauss_family->itemsInternal();
    ENUMERATE_CELL (icell, mesh->allCells()) {
      Cell cell = *icell;
      Int32 cell_owner = cell.owner();
      for (DoFLocalId gauss : cell_gauss.dofs(cell)) {
        gauss_dofs[gauss]->setOwner(cell_owner, my_rank);
      }
    }
    dof_family->notifyItemsOwnerChanged();
    dof_family->computeSynchronizeInfos();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void GaussDoFsOnCells::
initialize(IMesh* mesh, Int32 max_nb_gauss_per_cell)
{
  m_p->initialize(mesh, max_nb_gauss_per_cell);
  m_p->m_gauss_refpos = new VariableDoFReal3(VariableBuildInfo(mesh, "GaussRefPos", "GaussCellFamily"));
  m_p->m_gauss_weight = new VariableDoFReal(VariableBuildInfo(mesh, "GaussWeight", "GaussCellFamily"));
  m_p->m_gauss_jacobian = new VariableDoFReal(VariableBuildInfo(mesh, "GaussJacobian", "GaussCellFamily"));
  m_p->m_gauss_jacobmat = new VariableDoFReal3x3(VariableBuildInfo(mesh, "GaussJacobMat", "GaussCellFamily"));
  m_p->m_gauss_shape = new VariableDoFArrayReal(VariableBuildInfo(mesh, "GaussShape", "GaussCellFamily"));
  m_p->m_gauss_shapederiv = new VariableDoFArrayReal3(VariableBuildInfo(mesh, "GaussShapeDeriv", "GaussCellFamily"));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

VariableDoFReal3& GaussDoFsOnCells::
gaussRefPosition()
{
  return *m_p->m_gauss_refpos;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

VariableDoFReal& GaussDoFsOnCells::
gaussWeight()
{
  return *m_p->m_gauss_weight;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

VariableDoFReal& GaussDoFsOnCells::
gaussJacobian()
{
  return *m_p->m_gauss_jacobian;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

VariableDoFReal3x3& GaussDoFsOnCells::
gaussJacobMat()
{
  return *m_p->m_gauss_jacobmat;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

VariableDoFArrayReal& GaussDoFsOnCells::
gaussShape()
{
  return *m_p->m_gauss_shape;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

VariableDoFArrayReal3& GaussDoFsOnCells::
gaussShapeDeriv()
{
  return *m_p->m_gauss_shapederiv;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

IndexedCellDoFConnectivityView GaussDoFsOnCells::
gaussCellConnectivityView() const
{
  return m_p->m_cell_gauss_connectivity->view();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

IItemFamily* GaussDoFsOnCells::
gaussFamily() const
{
  return m_p->m_gauss_family;
}
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
