// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoFDoFLinearSystemImplBase.cc                                (C) 2000-2026 */
/*                                                                           */
/* Base class of implementation of IDoFLinearSystemImpl.                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "internal/DoFLinearSystemImplBase.h"

#include <arcane/core/IItemFamily.h>
#include <arcane/core/IParallelMng.h>
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

DoFLinearSystemImplBase::
DoFLinearSystemImplBase(IItemFamily* dof_family, const String& solver_name)
: TraceAccessor(dof_family->traceMng())
, m_dof_family(dof_family)
, m_rhs_variable(VariableBuildInfo(dof_family, solver_name + "RHSVariable"))
, m_solution_variable(VariableBuildInfo(dof_family, solver_name + "SolutionVariable"))
, m_dof_forced_info(VariableBuildInfo(dof_family, solver_name + "DoFForcedInfo"))
, m_dof_forced_value(VariableBuildInfo(dof_family, solver_name + "DoFForcedValue"))
, m_dof_elimination_info(VariableBuildInfo(dof_family, solver_name + "DoFEliminationInfo"))
, m_dof_elimination_value(VariableBuildInfo(dof_family, solver_name + "DoFEliminationValue"))
{}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystemImplBase::
clearValues()
{
  m_dof_forced_info.fill(false);
  m_dof_elimination_info.fill(ELIMINATE_NONE);
  m_dof_elimination_value.fill(0.0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystemImplBase::
setNearNullSpaceVectors(NumArray<Real, MDDim2>& vectors, Int32 block_size)
{
  if (block_size <= 0)
    ARCANE_FATAL("Invalid near-null-space block size '{0}'", block_size);

  m_near_null_space_block_size = block_size;
  if (vectors.extent0() == 0) {
    m_near_null_space_values.resize(0, 0);
    return;
  }

  IItemFamily* dof_family = dofFamily();
  const Int32 nb_dof = dof_family->allItems().size();
  const Int32 nb_mode = vectors.extent0();
  if (vectors.extent1() != nb_dof)
    ARCANE_FATAL("Near-null-space vectors have {0} values per mode, expected {1}", vectors.extent1(), nb_dof);

  m_near_null_space_values.swap(vectors);

  // Gram–Schmidt orthonormalization of the near-null-space vectors.
  // Done once for PETSc and Hypre, which require orthonormal vectors for their AMG preconditioners.
  IParallelMng* pm = dof_family->parallelMng();
  DoFGroup own_dofs = dof_family->allItems().own();
  for (Int32 i = 0; i < nb_mode; ++i) {
    for (Int32 j = 0; j < i; ++j) {
      Real dot = 0.0;
      ENUMERATE_DOF (idof, own_dofs)
        dot += m_near_null_space_values(i, idof.itemLocalId()) * m_near_null_space_values(j, idof.itemLocalId());
      dot = pm->reduce(Parallel::ReduceSum, dot);
      for (Int32 k = 0; k < nb_dof; ++k)
        m_near_null_space_values(i, k) -= dot * m_near_null_space_values(j, k);
    }

    Real norm2 = 0.0;
    ENUMERATE_DOF (idof, own_dofs) {
      Real value = m_near_null_space_values(i, idof.itemLocalId());
      norm2 += value * value;
    }
    norm2 = pm->reduce(Parallel::ReduceSum, norm2);
    if (norm2 <= 1.0e-30)
      ARCANE_FATAL("Near-null-space vector {0} is zero or linearly dependent", i);
    const Real inverse_norm = 1.0 / math::sqrt(norm2);
    for (Int32 k = 0; k < nb_dof; ++k)
      m_near_null_space_values(i, k) *= inverse_norm;
  }

  //info() << "Registered " << nb_mode << " near-null-space vectors (block size=" << block_size << ")";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DoFLinearSystemImplBase::
_applyRowColumnEliminationToRHS(bool is_verbose)
{
  DoFInfoListView item_list_view(m_dof_family);

  auto& dof_elimination_info = getEliminationInfo();
  auto& dof_elimination_value = getEliminationValue();
  auto& rhs_variable = rhsVariable();

  const bool do_print_filling = true;
  for (const auto& rc_value : m_row_column_elimination_map) {
    auto rc = rc_value.first;
    Real matrix_value = rc_value.second;
    DoF dof_row = item_list_view[rc.row_id];
    DoF dof_column = item_list_view[rc.column_id];
    if (dof_row == dof_column)
      continue;
    if (!dof_column.isOwn())
      continue;
    Byte row_elimination_info = dof_elimination_info[dof_row];
    Real elimination_value = dof_elimination_value[dof_row];
    // Subtract the value of RHS vector for current column.
    if (row_elimination_info == ELIMINATE_ROW_COLUMN) {
      Real v = rhs_variable[dof_column];
      rhs_variable[dof_column] = v - matrix_value * elimination_value;
      if (is_verbose)
        info() << "EliminateRowColumn (" << std::setw(4) << rc.row_id
               << "," << std::setw(4) << rc.column_id << ")"
               << " elimination_value=" << std::setw(25) << elimination_value
               << "  old_rhs=" << std::setw(25) << v
               << "  new_rhs=" << std::setw(25) << rhs_variable[dof_column];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
