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
  virtual String getDirichletInputFile() =0;
  virtual String getEnforceDirichletMethod() =0;
};

class IDirichletBoundaryCondition
{
 public:
  virtual Arcane::FaceGroup getSurface() =0;
  virtual StringConstArrayView getValue() =0;
  virtual Real getPenalty() =0;
  virtual String getDirichletInputFile() =0;
  virtual String getEnforceDirichletMethod() =0;
};

class INeumannBoundaryCondition
{
 public:
  virtual Arcane::FaceGroup getSurface() =0;
  virtual StringConstArrayView getValue() =0;
  virtual String getNeumannInputFile() =0;
};

class ITractionBoundaryCondition
{
 public:
  virtual Arcane::FaceGroup getSurface() =0;
  virtual StringConstArrayView getValue() =0;
  virtual String getTractionInputFile() =0;
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

class IPointCondition
{
public:
  virtual Arcane::NodeGroup getNode() =0;
  virtual StringConstArrayView getValue() =0;
  virtual String getPointInputFile() =0;
  virtual String getPointConditionType() =0;
};

class IParaxialCondition
{
public:
  virtual Arcane::FaceGroup getSurface() =0;
  virtual String getAInputFile() =0;
  virtual String getVInputFile() =0;
  virtual String getUInputFile() =0;
  virtual Integer getInputMotionType() =0;
  virtual Real getCp() =0;
  virtual Real getCs() =0;
  virtual Real getRho() =0;
  virtual Real getTp() =0;
  virtual Real getTs() =0;
  virtual Real getCoef() =0;
  virtual Real getAmplit() =0;
  virtual Real getPhase() =0;
  virtual Integer getOrder() =0;
  virtual Real getNormalAngle() =0;
  virtual Real getInPlaneAngle() =0;
};

class IArcaneFemBC
{
 public:
  virtual ~IArcaneFemBC() = default; 
  virtual ConstArrayView<BC::IPointCondition*> pointConditions() =0;
  virtual ConstArrayView<BC::IDirichletPointCondition*> dirichletPointConditions() =0;
  virtual ConstArrayView<BC::IDirichletBoundaryCondition*> dirichletBoundaryConditions() =0;
  virtual ConstArrayView<BC::INeumannBoundaryCondition*> neumannBoundaryConditions() =0;
  virtual ConstArrayView<BC::ITractionBoundaryCondition*> tractionBoundaryConditions() =0;
  virtual ConstArrayView<BC::IManufacturedSolution*> manufacturedSolutions() =0;
  virtual ConstArrayView<BC::IParaxialCondition*> paraxialBoundaryConditions() =0;
  virtual String getHandler() =0;
};

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif