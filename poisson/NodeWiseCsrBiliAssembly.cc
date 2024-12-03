// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* NWCSRiliAssembly.hxx                                      (C) 2022-2024   */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/* which avoid to add in the global matrix by iterating through the node.    */
/* It supports GPU Parallelization                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Builds the node-wise CSR matrix for the finite element method.
 *
 * This function initializes and fills the CSR (Compressed Sparse Row) matrix
 * based on the node-wise connectivity of the mesh. It supports different mesh
 * types, specifically "TRIA3" and "TETRA4".
 *
 * For "TRIA3" mesh type, the function iterates over all nodes and sets the
 * coordinates in the CSR matrix based on the node's faces.
 *
 * For "TETRA4" mesh type, the function iterates over all nodes and sets the
 * coordinates in the CSR matrix based on the node's edges.
 *
 * The function also initializes the memory space required for the CSR matrix
 * based on the number of nodes and edges/faces in the mesh.
 */
/*---------------------------------------------------------------------------*/

/* void FemModule::_buildMatrixNodeWiseCsrCPU()
{
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // Compute the number of nnz and initialize the memory space
  Int64 nbnde = nbNode();
  Int64 nedge = options()->meshType == "TETRA4" ? nbEdge() : nbFace();
  Int32 nnz = nedge * 2 + nbnde;
  m_csr_matrix.initialize(m_dof_family, nnz, nbnde);

  if (options()->meshType == "TRIA3") {
    ENUMERATE_NODE (inode, allNodes()) {

      //Since we compute the neighbouring connectivity here, we also fill the csr matrix

      Node node = *inode;

      m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(node, 0));

      for (Face face : node.faces()) {
        if (face.nodeId(0) == node.localId()) {
          //    cn->addConnectedItem(node, face.node(0));
          m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(1), 0));
        }
        else {
          //  cn->addConnectedItem(node, face.node(1));
          m_csr_matrix.setCoordinates(node_dof.dofId(node, 0), node_dof.dofId(face.nodeId(0), 0));
        }
      }
    }
  }
  else if (options()->meshType == "TETRA4") {
    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;
      Int32 node_dof_id = node_dof.dofId(node, 0);
      ItemLocalIdT<DoF> diagonal_entry(node_dof_id);

      m_csr_matrix.setCoordinates(diagonal_entry, diagonal_entry);

      for (Edge edge : node.edges()) {
        if (edge.nodeId(0) == node.localId())
          m_csr_matrix.setCoordinates(diagonal_entry, node_dof.dofId(edge.nodeId(1), 0));
        else
          m_csr_matrix.setCoordinates(diagonal_entry, node_dof.dofId(edge.nodeId(0), 0));
      }
    }
  }
}
*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_buildOffsetsNodeWiseCsr(const SmallSpan<uint>& offsets_smallspan)
{
  Accelerator::RunQueue* queue = acceleratorMng()->defaultQueue();

  // Initialize the neighbors array and shift right by one for CSR format
  NumArray<uint, MDDim1> neighbors(nbNode() + 1);
  neighbors[0] = 0;
  SmallSpan<uint> in_data = neighbors.to1DSmallSpan();

  // Select and execute the appropriate offset update based on mesh type
  if (mesh()->dimension() == 2) { // 2D mesh via node-face connectivity
    UnstructuredMeshConnectivityView connectivity_view(mesh());
    auto node_face_connectivity_view = connectivity_view.nodeFace();

    auto command = makeCommand(queue);
    command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
    {
      in_data[node_id + 1] = node_face_connectivity_view.nbFace(node_id) + 1;
    };
  }
  else { // 3D mesh via node-node connectivity
    auto* connectivity_ptr = m_node_node_via_edge_connectivity.get();
    ARCANE_CHECK_POINTER(connectivity_ptr);
    IndexedNodeNodeConnectivityView node_node_connectivity_view = connectivity_ptr->view();

    auto command = makeCommand(queue);
    command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
    {
      in_data[node_id + 1] = node_node_connectivity_view.nbNode(node_id) + 1;
    };
  }
  queue->barrier();

  // Do the inclusive sum for CSR row array (in_data)
  Accelerator::Scanner<uint> scanner;
  scanner.inclusiveSum(queue, in_data, offsets_smallspan);
}

void FemModule::
_buildMatrixNodeWiseCsr()
{
  Int8 mesh_dim = mesh()->dimension();

  Int32 nb_node = nbNode();
  Int32 nb_non_zero = nb_node + 2 * (mesh_dim == 2 ? nbFace() : m_nb_edge);
  m_csr_matrix.initialize(m_dof_family, nb_non_zero, nb_node, m_queue);

  NumArray<uint, MDDim1> offsets_numarray(nb_node + 1);
  SmallSpan<uint> offsets_smallspan = offsets_numarray.to1DSmallSpan();

  // Compute the array of offsets on Gpu
  _buildOffsetsNodeWiseCsr(offsets_smallspan);

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto out_m_matrix_row = viewOut(command, m_csr_matrix.m_matrix_row);
  auto inout_m_matrix_column = viewInOut(command, m_csr_matrix.m_matrix_column);

  // Select and execute the CSR matrix population based on mesh type
  if (mesh_dim == 2) { // 2D mesh via node-face & face-node connectivity
    UnstructuredMeshConnectivityView connectivity_view(mesh());

    auto node_face_connectivity_view = connectivity_view.nodeFace();
    auto face_node_connectivity_view = connectivity_view.faceNode();

    command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
    {
      // Retrieve the offset from the inclusive sum
      auto offset = offsets_smallspan[node_id];

      // Put the offset into CSR row array
      out_m_matrix_row[node_id] = offset;

      for (auto face_id : node_face_connectivity_view.faceIds(node_id)) {
        auto nodes = face_node_connectivity_view.nodes(face_id);

        // Put the neighbor of the current node into CSR column array
        inout_m_matrix_column[offset] = nodes[0] == node_id ? nodes[1] : nodes[0];

        ++offset;
      }

      inout_m_matrix_column[offset] = node_id; // Self-relation
    };
  }
  else { // 3D mesh via node-node connectivity
    auto* connectivity_ptr = m_node_node_via_edge_connectivity.get();
    ARCANE_CHECK_POINTER(connectivity_ptr);
    IndexedNodeNodeConnectivityView node_node_connectivity_view = connectivity_ptr->view();

    command << RUNCOMMAND_ENUMERATE(Node, node_id, allNodes())
    {
      auto offset = offsets_smallspan[node_id];
      out_m_matrix_row[node_id] = offset;

      for (auto neighbor_idx : node_node_connectivity_view.nodeIds(node_id)) {
        inout_m_matrix_column[offset] = neighbor_idx;
        ++offset;
      }

      inout_m_matrix_column[offset] = node_id;
    };
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_assembleNodeWiseCsrBilinearOperatorTria3()
{
  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_CsrNodeWise");
  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _buildMatrixNodeWiseCsr();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
  Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
  auto in_out_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);

  auto in_node_coord = ax::viewIn(command, m_node_coord);

  UnstructuredMeshConnectivityView m_connectivity_view;
  m_connectivity_view.setMesh(this->mesh());
  auto ncc = m_connectivity_view.nodeCell();
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());

  Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");

  command << RUNCOMMAND_ENUMERATE(Node, inode, allNodes())
  {
    Int32 inode_index = 0;
    for (auto cell : ncc.cells(inode)) {

      // How can I know the right index ?
      // By checking in the global id ?
      // Working currently, but maybe only because p = 1 ?
      auto nodeId1 = cnc.nodeId(cell, 1);
      auto nodeId2 = cnc.nodeId(cell, 2);
      inode_index = (inode == nodeId1) ? 1 : ((inode == nodeId2) ? 2 : 0);

      Real b_matrix[6] = { 0 };
      Real area = _computeCellMatrixGpuTRIA3(cell, cnc, in_node_coord, b_matrix);

      Int32 i = 0;
      for (NodeLocalId node2 : cnc.nodes(cell)) {
        Real x = b_matrix[inode_index * 2] * b_matrix[i * 2] + b_matrix[inode_index * 2 + 1] * b_matrix[i * 2 + 1];

        if (nodes_infos.isOwn(inode)) {
          Int32 row = node_dof.dofId(inode, 0).localId();
          Int32 col = node_dof.dofId(node2, 0).localId();
          Int32 begin = in_row_csr[row];
          Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];
          while (begin < end) {
            if (in_col_csr[begin] == col) {
              in_out_val_csr[begin] += x * area;
              break;
            }
            begin++;
          }
        }
        i++;
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the node-wise CSR bilinear operator for TETRA4 mesh type.
 *
 * This function builds the CSR matrix and computes the bilinear operator
 * for a TETRA4 mesh type. It initializes the CSR matrix, sets up the 
 * necessary views, and performs the computation on the accelerator.
 */
void FemModule::_assembleNodeWiseCsrBilinearOperatorTetra4()
{
  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_CsrNodeWise");
  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _buildMatrixNodeWiseCsr();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
  Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
  auto in_out_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);
  auto in_node_coord = ax::viewIn(command, m_node_coord);

  UnstructuredMeshConnectivityView m_connectivity_view;
  m_connectivity_view.setMesh(this->mesh());
  auto ncc = m_connectivity_view.nodeCell();
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());
  Arcane::ItemGenericInfoListView cells_infos(this->mesh()->cellFamily());

  Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");
  command << RUNCOMMAND_ENUMERATE(Node, inode, allNodes())
  {
    Int32 inode_index = 0;
    for (auto cell : ncc.cells(inode)) {

      if (inode == cnc.nodeId(cell, 1))
        inode_index = 1;
      else if (inode == cnc.nodeId(cell, 2))
        inode_index = 2;
      else if (inode == cnc.nodeId(cell, 3))
        inode_index = 3;
      else
        inode_index = 0;

      Real b_matrix[12] = { 0 };
      Real volume = _computeCellMatrixGpuTETRA4(cell, cnc, in_node_coord, b_matrix);

      Int32 i = 0;
      for (NodeLocalId node2 : cnc.nodes(cell)) {
        Real x = b_matrix[inode_index * 3] * b_matrix[i * 3] + b_matrix[inode_index * 3 + 1] * b_matrix[i * 3 + 1] + b_matrix[inode_index * 3 + 2] * b_matrix[i * 3 + 2];

        if (nodes_infos.isOwn(inode)) {
          Int32 row = node_dof.dofId(inode, 0).localId();
          Int32 col = node_dof.dofId(node2, 0).localId();
          Int32 begin = in_row_csr[row];
          Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];
          while (begin < end) {
            if (in_col_csr[begin] == col) {
              in_out_val_csr[begin] += x * volume;
              break;
            }
            ++begin;
          }
        }
        i++;
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
