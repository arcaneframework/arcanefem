// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
#ifndef IARCANEFEMBC_H
#define IARCANEFEMBC_H

#include <arcane/ItemTypes.h>
#include <arcane/VariableTypedef.h>
#include <arcane/core/IStandardFunction.h>

using namespace Arcane;

namespace BC
{

class IDirichletPointCondition
{
 public:
  virtual Arcane::NodeGroup getNode() =0;
  virtual StringConstArrayView getValue() =0;
  virtual Real getPenalty() =0;
  virtual String getEnforceDirichletMethod() =0;
};

class IDirichletBoundaryCondition
{
 public:
  virtual Arcane::FaceGroup getSurface() =0;
  virtual StringConstArrayView getValue() =0;
  virtual Real getPenalty() =0;
  virtual String getEnforceDirichletMethod() =0;
};

class INeumannBoundaryCondition
{
 public:
  virtual Arcane::FaceGroup getSurface() =0;
  virtual Real getValue() =0;
  virtual Real getValueX() =0;
  virtual Real getValueY() =0;
  virtual Real getValueZ() =0;
  virtual bool hasValue() const =0;
  virtual bool hasValueX() const =0;
  virtual bool hasValueY() const =0;
  virtual bool hasValueZ() const =0;
};

class IManufacturedSolution
{
 public:
  virtual bool getManufacturedSource() =0;
  virtual bool getManufacturedDirichlet() =0;
  virtual Real getPenalty() =0;
  virtual String getEnforceDirichletMethod() =0;
  virtual ICaseFunction* getManufacturedDirichletFunction() =0;
  virtual IStandardFunction* getManufacturedDirichletStandardFunction() =0;
  virtual ICaseFunction* getManufacturedSourceFunction() =0;
  virtual IStandardFunction* getManufacturedSourceStandardFunction() =0;
};

class IArcaneFemBC
{
 public:
  virtual ~IArcaneFemBC() = default; 
  virtual ConstArrayView<BC::IDirichletPointCondition*> dirichletPointConditions() =0;
  virtual ConstArrayView<BC::IDirichletBoundaryCondition*> dirichletBoundaryConditions() =0;
  virtual ConstArrayView<BC::INeumannBoundaryCondition*> neumannBoundaryConditions() =0;
  virtual ConstArrayView<BC::IManufacturedSolution*> manufacturedSolutions() =0;
  virtual String getHandler() =0;
};

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif