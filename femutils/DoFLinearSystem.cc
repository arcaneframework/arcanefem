// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoFLinearSystem.cc                                          (C) 2022-2023 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "DoFLinearSystem.h"

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/TraceAccessor.h>
#include <arcane/utils/NumArray.h>
#include <arcane/utils/ITraceMng.h>

#include <arcane/VariableTypes.h>
#include <arcane/IItemFamily.h>
#include <arcane/ISubDomain.h>
#include <arcane/IParallelMng.h>

#include "FemUtils.h"
#include "IDoFLinearSystemFactory.h"
#include "SequentialBasicDoFLinearSystemFactory_axl.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

extern "C++" DoFLinearSystemImpl*
createAlephDoFLinearSystemImpl(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class SequentialDoFLinearSystemImpl
: public TraceAccessor
, public DoFLinearSystemImpl
{
 public:

  SequentialDoFLinearSystemImpl(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name)
  : TraceAccessor(sd->traceMng())
  , m_sub_domain(sd)
  , m_dof_family(dof_family)
  , m_rhs_variable(VariableBuildInfo(dof_family, solver_name + "RHSVariable"))
  , m_dof_variable(VariableBuildInfo(dof_family, solver_name + "SolutionVariable"))
  {}

 public:

  void build()
  {
    Int32 nb_node = m_dof_family->allItems().size();
    m_k_matrix.resize(nb_node, nb_node);
    m_k_matrix.fill(0.0);
    m_rhs_vector.resize(nb_node);
    m_rhs_vector.fill(0.0);
  }

 private:

  void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    m_k_matrix(row, column) += value;
  }

  void solve() override
  {
    Int32 matrix_size = m_k_matrix.extent0();
    Arcane::MatVec::Matrix matrix(matrix_size, matrix_size);
    _convertNumArrayToCSRMatrix(matrix, m_k_matrix.span());
    bool is_verbose = true;
    Arcane::MatVec::Vector vector_b(matrix_size);
    Arcane::MatVec::Vector vector_x(matrix_size);
    {
      auto vector_b_view = vector_b.values();
      auto vector_x_view = vector_x.values();
      for (Int32 i = 0; i < matrix_size; ++i) {
        vector_b_view(i) = m_rhs_vector[i];
        if (is_verbose)
          info() << "VectorB[" << i << "] = " << m_rhs_vector[i];
        vector_x_view(i) = 0.0;
      }
    }

    bool use_direct_solver = matrix_size < 500;
    if (use_direct_solver) {
      info() << "Using direct solver";
      Arcane::MatVec::DirectSolver solver;
      solver.solve(matrix, vector_b, vector_x);
    }
    else {
      info() << "Using internal solver with diagonal preconditioner";
      Real epsilon = 1.0e-15;
      Arcane::MatVec::DiagonalPreconditioner p(matrix);
      Arcane::MatVec::ConjugateGradientSolver solver;
      solver.solve(matrix, vector_b, vector_x, epsilon, &p);
      info() << "End solver nb_iteration=" << solver.nbIteration()
             << " residual_norm=" << solver.residualNorm();
    }

    {
      auto vector_x_view = vector_x.values();
      ENUMERATE_ (DoF, idof, m_dof_family->allItems()) {
        DoF dof = *idof;
        Real v = vector_x_view[dof.localId()];
        m_dof_variable[dof] = v;
        if (is_verbose)
          info() << "VectorX[" << dof.localId() << "] = " << v;
      }
    }
  }

  void setRHSValues(Span<const Real> values) override
  {
    Int32 index = 0;
    DoFGroup own_dofs = m_dof_family->allItems().own();
    ENUMERATE_ (DoF, idof, own_dofs) {
      DoFLocalId dof_id = *idof;
      m_rhs_vector[dof_id] = values[index];
      ++index;
    }
  }

  VariableDoFReal& solutionVariable() override
  {
    return m_dof_variable;
  }

  VariableDoFReal& rhsVariable() override
  {
    return m_rhs_variable;
  }

 private:

  ISubDomain* m_sub_domain = nullptr;
  IItemFamily* m_dof_family = nullptr;
  VariableDoFReal m_rhs_variable;
  VariableDoFReal m_dof_variable;

  NumArray<Real, MDDim2> m_k_matrix;
  //! RHS (Right Hand Side) vector
  NumArray<Real, MDDim1> m_rhs_vector;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Default factory for Linear System.
 *
 * It is only used if there is no service specified in the 'arc' file and the
 * method DoFLinearSystem::setLinearSystemFactory() has not been called.
 */
class DefaultDoFLinearSystemFactory
: public IDoFLinearSystemFactory
{
 public:

  DoFLinearSystemImpl* createInstance(ISubDomain* sd, IItemFamily* dof_family,
                                      const String& solver_name) override
  {
    IParallelMng* pm = sd->parallelMng();
    bool is_parallel = pm->isParallel();
    // If true, we use a dense debug matrix in sequential
    bool use_debug_dense_matrix = true;
    DoFLinearSystemImpl* p = nullptr;
    if (is_parallel || !use_debug_dense_matrix) {
      p = createAlephDoFLinearSystemImpl(sd, dof_family, solver_name);
    }
    else {
      auto* x = new SequentialDoFLinearSystemImpl(sd, dof_family, solver_name);
      x->build();
      p = x;
    }
    return p;
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

DoFLinearSystem::
DoFLinearSystem()
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

DoFLinearSystem::
~DoFLinearSystem()
{
  delete m_p;
  delete m_default_linear_system_factory;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystem::
_checkInit()
{
  if (!m_p)
    ARCANE_FATAL("The instance is not initialized. You need to call initialize() before using this class");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystem::
initialize(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name)
{
  ARCANE_CHECK_POINTER(sd);
  ARCANE_CHECK_POINTER(dof_family);
  ITraceMng* tm = dof_family->traceMng();
  if (m_p)
    ARCANE_FATAL("The instance is already initialized");
  if (!m_linear_system_factory) {
    tm->info() << "Creating default Linear System Factory";
    if (!m_default_linear_system_factory)
      m_default_linear_system_factory = new DefaultDoFLinearSystemFactory();
    m_linear_system_factory = m_default_linear_system_factory;
  }
  m_item_family = dof_family;
  m_p = m_linear_system_factory->createInstance(sd, dof_family, solver_name);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystem::
matrixAddValue(DoFLocalId row, DoFLocalId column, Real value)
{
  _checkInit();
  m_p->matrixAddValue(row, column, value);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystem::
setRHSValues(Span<const Real> values)
{
  _checkInit();
  m_p->setRHSValues(values);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystem::
solve()
{
  _checkInit();

  {
    // For the LinearSystem class we need an array
    // with only the values for the ownNodes().
    // The values of 'rhs_values' should not be updated after
    // this call.
    UniqueArray<Real> rhs_values_for_linear_system;
    VariableDoFReal& rhs_values(rhsVariable());
    ENUMERATE_ (DoF, idof, m_item_family->allItems().own()) {
      rhs_values_for_linear_system.add(rhs_values[idof]);
    }
    setRHSValues(rhs_values_for_linear_system);
  }

  m_p->solve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

VariableDoFReal& DoFLinearSystem::
solutionVariable()
{
  _checkInit();
  return m_p->solutionVariable();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

VariableDoFReal& DoFLinearSystem::
rhsVariable()
{
  _checkInit();
  return m_p->rhsVariable();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystem::
reset()
{
  delete m_p;
  m_p = nullptr;
  m_item_family = nullptr;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class SequentialBasicDoFLinearSystemFactoryService
: public ArcaneSequentialBasicDoFLinearSystemFactoryObject
{
 public:

  SequentialBasicDoFLinearSystemFactoryService(const ServiceBuildInfo& sbi)
  : ArcaneSequentialBasicDoFLinearSystemFactoryObject(sbi)
  {
  }

  DoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name) override
  {
    IParallelMng* pm = sd->parallelMng();
    if (pm->isParallel())
      ARCANE_FATAL("This service is not available in parallel");
    auto* x = new SequentialDoFLinearSystemImpl(sd, dof_family, solver_name);
    x->build();
    return x;
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_SEQUENTIALBASICDOFLINEARSYSTEMFACTORY(SequentialBasicLinearSystem,
                                                              SequentialBasicDoFLinearSystemFactoryService);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
