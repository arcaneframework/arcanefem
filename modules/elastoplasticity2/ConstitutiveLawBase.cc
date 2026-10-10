// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ConstitutiveLawBase.cc                                      (C) 2000-2026 */
/*                                                                           */
/* Base class of a constitutive law.                                         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <modules/elastoplasticity2/ConstitutiveLawBase.h>

#include <arcane/core/ServiceBuildInfo.h>
#include <arcane/core/IMesh.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::ArcaneFem
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ConstitutiveLawBase::
ConstitutiveLawBase(const ServiceBuildInfo& sbi)
: BasicService(sbi)
, m_node_coord(sbi.mesh()->nodesCoordinates())
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void ConstitutiveLawBase::
_initialize(const ConstitutiveLawInitInfo& x)
{
  m_nGP = x.m_nGP;
  m_use_gpu = x.m_use_gpu;
  m_hex_quad_mesh = x.m_hex_quad_mesh;
  m_use_gpu_functions = x.m_use_gpu_functions;
  m_nodes_per_cell = x.m_nodes_per_cell;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::ArcaneFem

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
