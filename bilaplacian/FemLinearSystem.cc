// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemLinearSystem.cc                                          (C) 2022-2022 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemLinearSystem.h"

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

using namespace Arcane;

extern "C++" FemLinearSystemImpl*
createAlephFemLinearSystemImpl(ISubDomain* sd, const Arcane::VariableNodeReal& node_variable);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class SequentialFemLinearSystemImpl
: public TraceAccessor
, public FemLinearSystemImpl
{
 public:

  SequentialFemLinearSystemImpl(ISubDomain* sd, const Arcane::VariableNodeReal& node_variable)
  : TraceAccessor(sd->traceMng())
  , m_sub_domain(sd)
  , m_node_family(node_variable.variable()->itemFamily())
  , m_node_variable(node_variable)
  {}

 public:

  void build()
  {
    Int32 nb_node = m_node_family->allItems().size();
    m_k_matrix.resize(nb_node*2, nb_node*2);
    m_k_matrix.fill(0.0);
    m_rhs_vector.resize(nb_node*2);
    m_rhs_vector.fill(0.0);
  }

 private:

  void matrixAddValue(NodeLocalId row, NodeLocalId column, Real value1, Real value2, Real value3, Real value4) override
  {
//        Real v1  = K_e(2*n1_index  ,   2*n2_index   );
//        Real v2  = K_e(2*n1_index  ,   2*n2_index+1 );
//        Real v3  = K_e(2*n1_index+1,   2*n2_index   );
//        Real v4  = K_e(2*n1_index+1,   2*n2_index+1 );  
    m_k_matrix.span()(row*2  , column*2  )   += value1;
    m_k_matrix.span()(row*2  , column*2+1)   += value2;
    m_k_matrix.span()(row*2+1, column*2  )   += value3;
    m_k_matrix.span()(row*2+1, column*2+1)   += value4;
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
        vector_b_view(i)   = m_rhs_vector[i];
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
      for (Int32 i = 0; i < matrix_size; ++i) {
                     cout << "SOL[" << i <<"] = "<< vector_x_view(i) << endl;   
      }      
      ENUMERATE_ (Node, inode, m_node_family->allItems()) {
        Node node = *inode;
        m_node_variable[node] = vector_x_view[node.localId()*2];

      }
    }
  }

  void setRHSValues(Span<const Real> values) override
  {
    Int32 index = 0;
    NodeGroup own_nodes = m_node_family->allItems().own();
    ENUMERATE_ (Node, inode, own_nodes) {
      NodeLocalId node_id = *inode;
      m_rhs_vector[node_id] = values[index];
      ++index;
    }
  }

 private:

  ISubDomain* m_sub_domain = nullptr;
  IItemFamily* m_node_family = nullptr;
  VariableNodeReal m_node_variable;

  NumArray<Real, MDDim2> m_k_matrix;
  //! RHS (Right Hand Side) vector
  NumArray<Real, MDDim1> m_rhs_vector;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FemLinearSystem::
FemLinearSystem()
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FemLinearSystem::
~FemLinearSystem()
{
  delete m_p;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemLinearSystem::
_checkInit()
{
  if (!m_p)
    ARCANE_FATAL("The instance is not initialized. You need to call initialize() before using this class");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemLinearSystem::
initialize(ISubDomain* sd, const Arcane::VariableNodeReal& node_variable)
{
  ARCANE_CHECK_POINTER(sd);
  if (m_p)
    ARCANE_FATAL("The instance is already initialized");
  IParallelMng* pm = sd->parallelMng();
  bool is_parallel = pm->isParallel();
  if (is_parallel) {
    m_p = createAlephFemLinearSystemImpl(sd, node_variable);
  }
  else {
    auto* x = new SequentialFemLinearSystemImpl(sd, node_variable);
    x->build();
    m_p = x;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemLinearSystem::
matrixAddValue(NodeLocalId row, NodeLocalId column, Real value1, Real value2, Real value3, Real value4)
{
  _checkInit();
  m_p->matrixAddValue(row, column, value1, value2, value3, value4);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemLinearSystem::
setRHSValues(Span<const Real> values)
{
  _checkInit();
  m_p->setRHSValues(values);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemLinearSystem::
solve()
{
  _checkInit();
  m_p->solve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemLinearSystem::
reset()
{
  delete m_p;
  m_p = nullptr;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
