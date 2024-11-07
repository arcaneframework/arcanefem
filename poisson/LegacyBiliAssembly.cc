// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* LegacyBiliAssembly.hxx                                    (C) 2022-2024   */
/*                                                                           */
/* Methods of the bilinear assembly phase using the DOK data structure       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

void FemModule::
_assembleBilinearOperatorTRIA3()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Timer::Action timer_action(m_time_stats, "AssembleLegacyBilinearOperatorTria3");

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    FixedMatrix<3, 3> K_e;
    {
      K_e = _computeElementMatrixTRIA3(cell); // element stiffness matrix
    }

    {
      Int32 n1_index = 0;
      for (Node node1 : cell.nodes()) {
        Int32 n2_index = 0;
        for (Node node2 : cell.nodes()) {
          Real v = K_e(n1_index, n2_index);
          if (node1.isOwn()) {
            m_linear_system.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
          }
          ++n2_index;
        }
        ++n1_index;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTETRA4()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Timer::Action timer_action(m_time_stats, "AssembleLegacyBilinearOperatorTetra4");

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    auto K_e = _computeElementMatrixTETRA4(cell); // element stiffness matrix
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        Real v = K_e(n1_index, n2_index);
        if (node1.isOwn()) {
          m_linear_system.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/