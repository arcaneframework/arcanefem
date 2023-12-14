// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* LegacyBiliAssembly.hxx                                    (C) 2022-2023   */
/*                                                                           */
/* Methods of the bilinear assembly phase using the legacy data structure    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTRIA3()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Timer::Action timer_action(this->subDomain(), "AssembleLegacyBilinearOperatorTria3");

  std::chrono::_V2::system_clock::time_point lhs_start;
  double compute_average = 0;
  double global_build_average = 0;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using hashmap legacy format\n";
    lhs_start = std::chrono::high_resolution_clock::now();
  }

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;

    std::chrono::_V2::system_clock::time_point compute_El_start;
    if (m_register_time) {
      compute_El_start = std::chrono::high_resolution_clock::now();
    }

    FixedMatrix<3, 3> K_e;
    {
      Timer::Action timer_action(this->subDomain(), "LegacyComputeElementMatrix");
      K_e = _computeElementMatrixTRIA3(cell); // element stifness matrix
    }

    if (m_register_time) {
      auto compute_El_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> compute_duration = compute_El_stop - compute_El_start;
      compute_average += compute_duration.count();
    }

    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]

    std::chrono::_V2::system_clock::time_point global_build_start;
    if (m_register_time) {
      global_build_start = std::chrono::high_resolution_clock::now();
    }

    {
      Timer::Action timer_action(this->subDomain(), "LegacyAddToGlobalMatrix");
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
    if (m_register_time) {
      auto global_build_stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> global_build_duration = global_build_stop - global_build_start;
      global_build_average += global_build_duration.count();
    }
  }

  if (m_register_time) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - lhs_start;
    double lhs_loc_time = duration.count();
    logger << "compute elements average time : " << compute_average / nbCell() << "\n"
           << "compute elements total time : " << compute_average << "\n"
           << "add in global matrix average time : " << global_build_average / nbCell() << "\n"
           << "add in global matrix total time : " << global_build_average << "\n"
           << "lhs total time : " << lhs_loc_time << "\n"
           << "compute element time in lhs : " << compute_average / lhs_loc_time * 100 << "%\n"
           << "add in global matrix time in lhs : " << global_build_average / lhs_loc_time * 100 << "%\n\n"
           << "-------------------------------------------------------------------------------------\n\n";
    lhs_time += lhs_loc_time;
    wbuild << lhs_loc_time << ",";
    timer << lhs_loc_time << ",";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/