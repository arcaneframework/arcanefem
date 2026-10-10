// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ConstitutiveLawBase.h                                       (C) 2000-2026 */
/*                                                                           */
/* Base class of a constitutive law.                                         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_CONSTITUTIVELAWBASE_H
#define ARCANEFEM_CONSTITUTIVELAWBASE_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/String.h>
#include <arcane/core/BasicService.h>

#include <femutils/FemUtils.h>
#include <modules/elastoplasticity2/IConstitutiveLaw.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::ArcaneFem
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class ConstitutiveLawBase
: public BasicService
, public IConstitutiveLaw
{
 public:

  ConstitutiveLawBase(const ServiceBuildInfo& sbi);

 public:

  String lawName() const override { return m_law_name; }
  Real getMu() const override { return mu; }
  Real getLambda() const override { return lambda; }
  FemUtils::RealMatrix<3,3> getElasticityMatrix2D() const override { return m_C_elas_2d; }

  void setTimeStep(Real v) override { dt = v; }

 protected:

  String m_law_name;
  Int16 m_nGP = 1;
  bool m_use_gpu = false;
  bool m_hex_quad_mesh = false;
  bool m_use_gpu_functions = false;
  Int16 m_nodes_per_cell = 0;

  Real dt = 0.0; // This is set by the service user

  Real mu = 0.0;
  Real lambda = 0.0;
  Real H = 0.0;

  VariableNodeReal3 m_node_coord;

  FemUtils::RealMatrix<3, 3> m_C_elas_2d;
  FemUtils::RealMatrix<6, 6> m_C_elas_3d;

 protected:

  void _initialize(const ConstitutiveLawInitInfo& x);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::ArcaneFem

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
