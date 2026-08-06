// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemDoFsOnCells.h                                            (C) 2022-2026 */
/*                                                                           */
/* Manage one or more DoFs on Cells.                                         */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_FEMDOFSONCELLS_H
#define ARCANEFEM_FEMUTILS_FEMDOFSONCELLS_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/ItemTypes.h>
#include <arcane/IndexedItemConnectivityView.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*!
 * \brief Manage one or more DoFs on Cells.
 *
 * Before using an instance of this class you need to call method
 * initialize().
 *
 * After initialization, you can access DoFs like that:
 *
 * \code
 * FemDoFsOnCells dofs_on_cells = ...;
 * auto cell_dof(m_dofs_on_cells.cellDoFConnectivityView());
 * Cell cell = ...
 * DoFLocalId dof0 = cell_dof.dofId(cell, 0); //< First DoF of the cell
 * DoFLocalId dof1 = cell_dof.dofId(cell, 1); //< Second DoF of the cell
 * \endcode
 */
class FemDoFsOnCells
{
  class Impl;

 public:

  FemDoFsOnCells(Arcane::ITraceMng* tm);
  ~FemDoFsOnCells();

 public:

  FemDoFsOnCells(const FemDoFsOnCells&) = delete;
  FemDoFsOnCells(FemDoFsOnCells&&) = delete;
  FemDoFsOnCells& operator=(FemDoFsOnCells&&) = delete;
  FemDoFsOnCells& operator=(const FemDoFsOnCells&) = delete;

 public:

  /*!
   * \brief Initialize the instance.
   */
  void initialize(Arcane::IMesh* mesh, Arcane::Int32 nb_dof_per_cell);

 public:

  Arcane::IndexedCellDoFConnectivityView cellDoFConnectivityView() const;
  Arcane::IItemFamily* dofFamily() const;

 private:

  Impl* m_p = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif