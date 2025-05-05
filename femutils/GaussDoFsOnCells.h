// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/* GaussDoFsOnCells.h                                          (C) 2022-2025 */
/* Manage one or more Gauss points for FEM integration on Cells              */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef GAUSSDOFSONCELLS_H
#define GAUSSDOFSONCELLS_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemUtils.h"
#include <arcane/ItemTypes.h>
#include <arcane/IndexedItemConnectivityView.h>
#include <arcane/IMesh.h>
#include "MeshTensorVariable.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
/*!
 * \brief Manage one or more Gauss Points on Cells through the Arcane DoFs mechanism.
 *
 * Before using an instance of this class you need to call method
 * initialize().
 */
class GaussDoFsOnCells
{
  class Impl;

 public:

  explicit GaussDoFsOnCells(Arcane::ITraceMng* tm);
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

  Arcane::VariableDoFArrayTensor2& gaussStress();
  Arcane::VariableDoFArrayTensor2& gaussStrain();
  Arcane::VariableDoFArrayTensor2& gaussStrainPlastic();

 private:

  Impl* m_p = nullptr;
};
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif //GAUSSDOFSONCELLS_H
