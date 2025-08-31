// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------

#include "IArcaneFemBC.h"
#include "FemBoundaryConditions_axl.h"

using namespace Arcane;


class FemBoundaryConditionsService 
: public ArcaneFemBoundaryConditionsObject
{
 public:
  explicit FemBoundaryConditionsService(const ServiceBuildInfo & sbi)
  : ArcaneFemBoundaryConditionsObject(sbi) {}
  
  ConstArrayView<BC::IDirichletPointCondition*> dirichletPointConditions() { return options()->getDirichletPoint(); }
  ConstArrayView<BC::IDirichletBoundaryCondition*> dirichletBoundaryConditions() { return options()->getDirichlet(); }
  ConstArrayView<BC::INeumannBoundaryCondition*> neumannBoundaryConditions() { return options()->getNeumann(); }
  ConstArrayView<BC::ITractionBoundaryCondition*> tractionBoundaryConditions() { return options()->getTraction(); }
  ConstArrayView<BC::IManufacturedSolution*> manufacturedSolutions() { return options()->getManufacturedSolution(); }
  String getHandler() { return options()->getHandler();}
};
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_FEMBOUNDARYCONDITIONS(FemBoundaryConditions, FemBoundaryConditionsService);