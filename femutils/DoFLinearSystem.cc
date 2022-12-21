// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoFLinearSystem.cc                                          (C) 2022-2022 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "DoFLinearSystem.h"

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/TraceAccessor.h>
#include <arcane/utils/NumArray.h>

#include <arcane/VariableTypes.h>
#include <arcane/IItemFamily.h>
#include <arcane/ISubDomain.h>
#include <arcane/IParallelMng.h>

#include "FemUtils.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

extern "C++" DoFLinearSystemImpl*
createAlephFemLinearSystem2Impl(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class SequentialFemLinearSystem2Impl
: public TraceAccessor
, public DoFLinearSystemImpl
{
 public:

  SequentialFemLinearSystem2Impl(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name)
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

    Arcane::MatVec::Vector vector_b(matrix_size);
    Arcane::MatVec::Vector vector_x(matrix_size);
    {
      auto vector_b_view = vector_b.values();
      auto vector_x_view = vector_x.values();
      for (Int32 i = 0; i < matrix_size; ++i) {
        vector_b_view(i) = m_rhs_vector[i];
        vector_x_view(i) = 0.0;
      }
    }

    {
      Real epsilon = 1.0e-15;
      Arcane::MatVec::DiagonalPreconditioner p(matrix);
      Arcane::MatVec::ConjugateGradientSolver solver;
      solver.solve(matrix, vector_b, vector_x, epsilon, &p);
    }

    {
      auto vector_x_view = vector_x.values();
      ENUMERATE_ (DoF, idof, m_dof_family->allItems()) {
        DoF dof = *idof;
        m_dof_variable[dof] = vector_x_view[dof.localId()];
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
  if (m_p)
    ARCANE_FATAL("The instance is already initialized");
  m_item_family = dof_family;
  IParallelMng* pm = sd->parallelMng();
  bool is_parallel = pm->isParallel();
  // If true, we use a dense debug matrix in sequential
  bool use_debug_dense_matrix = false;
  if (is_parallel || !use_debug_dense_matrix) {
    m_p = createAlephFemLinearSystem2Impl(sd, dof_family, solver_name);
  }
  else {
    auto* x = new SequentialFemLinearSystem2Impl(sd, dof_family, solver_name);
    x->build();
    m_p = x;
  }
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

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
