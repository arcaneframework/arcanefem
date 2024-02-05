// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* LegacyBiliAssembly.hxx                                    (C) 2022-2023   */
/*                                                                           */
/* Methods of the bilinear assembly phase using the legacy data structure    */
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
      //Timer::Action timer_action(m_time_stats, "LegacyComputeElementMatrix");
      K_e = _computeElementMatrixTRIA3(cell); // element stifness matrix
    }

    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]

    {
      //Timer::Action timer_action(m_time_stats, "LegacyAddToGlobalMatrix");
      Int32 n1_index = 0;
      for (Node node1 : cell.nodes()) {
        Int32 n2_index = 0;
        for (Node node2 : cell.nodes()) {
          // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
          Real v = K_e(n1_index, n2_index);
          // m_k_matrix(node1.localId(), node2.localId()) += v;
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
