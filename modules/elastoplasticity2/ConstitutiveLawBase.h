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

  ConstitutiveLawBase(const ServiceBuildInfo& sbi)
  : BasicService(sbi)
  {}

 public:

  String lawName() const override { return m_law_name; }
  Real getMu() const override { return mu; }
  Real getLambda() const override { return lambda; }

 protected:

  String m_law_name;
  Int16 m_nGP = 1;
  bool m_use_gpu = false;
  bool m_hex_quad_mesh = false;
  bool m_use_gpu_functions = false;
  Int16 m_nodes_per_cell = 0;

  Real mu = 0.0;
  Real lambda = 0.0;
  Real H = 0.0;

  RealMatrix<3, 3> m_C_elas_2d;
  RealMatrix<6, 6> m_C_elas_3d;

 protected:

  void _initialize(const ConstitutiveLawInitInfo& x)
  {
    m_nGP = x.m_nGP;
    m_use_gpu = x.m_use_gpu;
    m_hex_quad_mesh = x.m_hex_quad_mesh;
    m_use_gpu_functions = x.m_use_gpu_functions;
    m_nodes_per_cell = x.m_nodes_per_cell;
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::ArcaneFem

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
