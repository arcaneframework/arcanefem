// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrBiliAssembly.hxx                                         (C) 2022-2023 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Initialization of the csr matrix. It only works for p=1 since there is
 * one node per Edge.
 * 
 * 
 */
void FemModule::
_buildMatrixCsr()
{
  //Initialization of the coo matrix;
  //This formula only works in p=1

  /*
  //Create a connection between nodes through the faces
  //Useless here because we only need this information once
  IItemFamily* node_family = mesh()->nodeFamily();
  NodeGroup nodes = node_family->allItems();
  auto idx_cn = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(node_family, node_family, "NodeToNeighbourFaceNodes");
  auto* cn = idx_cn->connectivity();
  ENUMERATE_NODE (node, allNodes()) {
  }
  */

  Int32 nnz = nbFace() * 2 + nbNode();
  m_csr_matrix.initialize(m_dof_family, nnz, nbNode());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  //We iterate through the node, and we do not sort anymore : we assume the nodes ID are sorted, and we will iterate throught the column to avoid making < and > comparison
  ENUMERATE_NODE (inode, allNodes()) {
    Node node = *inode;

    m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

    for (Face face : node.faces()) {
      if (face.nodeId(0) == node.localId())
        m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
      else
        m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleCsrBilinearOperatorTRIA3()
{

  Timer::Action timer_csr_bili(m_time_stats, "AssembleCsrBilinearOperatorTria3");

  double compute_average = 0;
  double global_build_average = 0;
  double build_time = 0;
  std::chrono::_V2::system_clock::time_point lhs_start;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using CPU CSR with NumArray format\n";
    lhs_start = std::chrono::high_resolution_clock::now();
  }
  {
    Timer::Action timer_csr_build(m_time_stats, "CsrBuildMatrix");
    // Build the csr matrix
    _buildMatrixCsr();
  }

  std::chrono::_V2::system_clock::time_point build_stop;
  if (m_register_time) {
    build_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> build_duration = build_stop - lhs_start;
    build_time = build_duration.count();
  }

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    std::chrono::_V2::system_clock::time_point compute_El_start;
    if (m_register_time)
      compute_El_start = std::chrono::high_resolution_clock::now();

    FixedMatrix<3, 3> K_e;
    {
      //Timer::Action timer_csr_compute_add(m_time_stats, "CsrComputeElementMatrixTria3");
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
    //Timer::Action timer_action(m_time_stats, "CsrAddToGlobalMatrix");
    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      Int32 n2_index = 0;
      for (Node node2 : cell.nodes()) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        Real v = K_e(n1_index, n2_index);
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (node1.isOwn()) {
          m_csr_matrix.matrixAddValue(node_dof.dofId(node1, 0), node_dof.dofId(node2, 0), v);
        }
        ++n2_index;
      }
      ++n1_index;
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
    logger << "Building time of the csr matrix :" << build_time << "\n"
           << "Compute Elements average time : " << compute_average / nbCell() << "\n"
           << "Compute Elements total time : " << compute_average << "\n"
           << "Add in global matrix average time : " << global_build_average / nbCell() << "\n"
           << "Add in global matrix total time : " << global_build_average << "\n"
           << "LHS Total time : " << lhs_loc_time << "\n"
           << "Build matrix time in lhs :" << build_time / lhs_loc_time * 100 << "%\n"
           << "Compute element time in lhs : " << compute_average / lhs_loc_time * 100 << "%\n"
           << "Add in global matrix time in lhs : " << global_build_average / lhs_loc_time * 100 << "%\n\n"
           << "-------------------------------------------------------------------------------------\n\n";
    lhs_time += lhs_loc_time;
    wbuild << lhs_loc_time << ",";
    timer << compute_average + global_build_average << ",";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/