// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemDoFsOnNodes.h                                            (C) 2022-2022 */
/*                                                                           */
/* Manage one or more DoFs on Nodes.                                         */
/*---------------------------------------------------------------------------*/
#ifndef FEMTEST_FEMDOFSONNODES_H
#define FEMTEST_FEMDOFSONNODES_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/ItemTypes.h>
#include <arcane/IndexedItemConnectivityView.h>

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
class FemDoFsOnNodes
{
  class Impl;

 public:

  FemDoFsOnNodes(Arcane::ITraceMng* tm);
  ~FemDoFsOnNodes();

 public:

  FemDoFsOnNodes(const FemDoFsOnNodes&) = delete;
  FemDoFsOnNodes(FemDoFsOnNodes&&) = delete;
  FemDoFsOnNodes& operator=(FemDoFsOnNodes&&) = delete;
  FemDoFsOnNodes& operator=(const FemDoFsOnNodes&) = delete;

 public:

  /*!
   * \brief Initialize the instance.
   */
  void initialize(Arcane::IMesh* mesh, Arcane::Int32 nb_dof_per_node);

 public:

  Arcane::IndexedNodeDoFConnectivityView nodeDoFConnectivityView() const;
  Arcane::IItemFamily* dofFamily() const;

 private:

  Impl* m_p = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
