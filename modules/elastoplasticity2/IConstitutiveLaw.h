// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* IConstitutiveLaw.h                                          (C) 2000-2026 */
/*                                                                           */
/* Interface of a constitutive law.                                          */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_ICONSTITUTIVELAW_H
#define ARCANEFEM_ICONSTITUTIVELAW_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/core/ItemTypes.h>
#include <arccore/common/accelerator/RunQueue.h>

#include <femutils/FemUtilsGlobal.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::ArcaneFem
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

struct ConstitutiveLawInitInfo
{
  Int16 m_nGP = 1;
  bool m_use_gpu = false;
  bool m_hex_quad_mesh = false;
  bool m_use_gpu_functions = false;
  Int16 m_nodes_per_cell = 0;
  RunQueue m_run_queue;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class IConstitutiveLaw
{
 public:

  virtual ~IConstitutiveLaw() = default;

 public:

  virtual String lawName() const =0;
  virtual void initialize(const ConstitutiveLawInitInfo&) = 0;
  virtual void getMaterialProperties() = 0;
  virtual void integrateAndSave() = 0;
  virtual void restoreConvergedState() =0;
  virtual void commitInternalVariables() =0;
  virtual Real getMu() const = 0;
  virtual Real getLambda() const = 0;
  virtual FemUtils::RealMatrix<3, 3> getElasticityMatrix2D() const = 0;

  // TODO: Should not be here
  virtual void setTimeStep(Real dt) =0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::ArcaneFem

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
