// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* GaussDoFsOnCells.h                                          (C) 2022-2024 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef PASSMO_GAUSSDOFSONCELLS_H
#define PASSMO_GAUSSDOFSONCELLS_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemUtils.h"
#include <arcane/ItemTypes.h>
#include <arcane/IndexedItemConnectivityView.h>
#include <arcane/IMesh.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*!
 * \brief Manage one or more DoFs on Nodes.
 *
 * Before using an instance of this class you need to call method
 * initialize().
 *
 * After initialization, you can access DoFs like that:
 *
 * \code
 * FemDoFsOnNodes dofs_on_nodes = ...;
 * auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
 * Node node = ...
 * DoFLocalId dof0 = node_dof.dofId(node, 0); //< First DoF of the node
 * DoFLocalId dof1 = node_dof.dofId(node, 1); //< Second DoF of the node
 * \endcode
 */
class GaussDoFsOnCells
{
  class Impl;

 public:

  GaussDoFsOnCells(Arcane::ITraceMng* tm);
  ~GaussDoFsOnCells();

 public:

  GaussDoFsOnCells(const GaussDoFsOnCells&) = delete;
  GaussDoFsOnCells(GaussDoFsOnCells&&) = delete;
  GaussDoFsOnCells& operator=(GaussDoFsOnCells&&) = delete;
  GaussDoFsOnCells& operator=(const GaussDoFsOnCells&) = delete;

 public:

  /*!
   * \brief Initialize the instance.
   */
  void initialize(Arcane::IMesh* mesh, Arcane::Int32 max_nb_gauss_per_cell);

 public:

  [[nodiscard]] Arcane::IndexedCellDoFConnectivityView gaussCellConnectivityView() const;
  [[nodiscard]] Arcane::IItemFamily* gaussFamily() const;
  Arcane::VariableDoFArrayReal& gaussShape();
  Arcane::VariableDoFArrayReal3& gaussShapeDeriv();
  Arcane::VariableDoFReal3& gaussRefPosition();
  Arcane::VariableDoFReal3x3& gaussJacobMat();
  Arcane::VariableDoFReal& gaussWeight();
  Arcane::VariableDoFReal& gaussJacobian();

 private:

  Impl* m_p = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif //PASSMO_GAUSSDOFSONCELLS_H
