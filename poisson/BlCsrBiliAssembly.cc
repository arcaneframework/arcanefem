// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BlcsrBiliAssembly.hxx                                     (C) 2022-2024   */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/* which avoid to add in the global matrix by iterating through the node.    */
/* It supports GPU Parallelization                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"

void FemModule::
_buildMatrixBuildLessCsr()
{

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // Compute the number of nnz and initialize the memory space
  Integer nbnde = nbNode();
  Int64 nedge;

  if (options()->meshType == "TETRA4")
    nedge = nbEdge();
  else if (options()->meshType == "TRIA3")
    nedge = nbFace();
  else
    ARCANE_THROW(NotImplementedException, "");

  Int32 nnz = nedge * 2 + nbnde;
  m_csr_matrix.initialize(m_dof_family, nnz, nbnde, m_queue);

  Integer index = 1;
  m_csr_matrix.m_matrix_row(0) = 0;

  if (options()->meshType == "TETRA4")
    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;
      if (index < nbnde) {
        m_csr_matrix.m_matrix_row(index) = node.nbEdge() + m_csr_matrix.m_matrix_row(index - 1) + 1;
        index++;
      }
    }
  else if (options()->meshType == "TRIA3")
    ENUMERATE_NODE (inode, allNodes()) {
      Node node = *inode;
      if (index < nbnde) {
        m_csr_matrix.m_matrix_row(index) = node.nbFace() + m_csr_matrix.m_matrix_row(index - 1) + 1;
        index++;
      }
    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_buildMatrixGpuBuildLessCsr()
{

  // Compute the number of nnz and initialize the memory space
  Integer nbnde = nbNode();
  Int64 nedge;

  if (options()->meshType == "TETRA4")
    nedge = m_nb_edge;
  else if (options()->meshType == "TRIA3")
    nedge = nbFace();
  else
    ARCANE_THROW(NotImplementedException, "");

  Int32 nnz = nedge * 2 + nbnde;

  NumArray<Int32, MDDim1> tmp_row;
  tmp_row.resize(nbnde);
  m_csr_matrix.initialize(m_dof_family, nnz, nbnde, m_queue);

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);
  auto in_out_tmp_row = ax::viewInOut(command, tmp_row);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  UnstructuredMeshConnectivityView connectivity_view;
  connectivity_view.setMesh(this->mesh());

  if (options()->meshType == "TRIA3") {
    auto nfc = connectivity_view.nodeFace();
    command << RUNCOMMAND_ENUMERATE(Node, inode, allNodes())
    {
      Int64 index = node_dof.dofId(inode, 0).localId();
      in_out_tmp_row[index] = nfc.nbFace(inode) + 1;
    };
  }
  else { // 3D mesh via node-node connectivity
    auto* connectivity_ptr = m_node_node_via_edge_connectivity.get();
    ARCANE_CHECK_POINTER(connectivity_ptr);
    IndexedNodeNodeConnectivityView node_node_connectivity_view = connectivity_ptr->view();
    auto command = makeCommand(queue);
    command << RUNCOMMAND_ENUMERATE(Node, inode, allNodes())
    {
      Int64 index = node_dof.dofId(inode, 0).localId();
      in_out_tmp_row[index] = node_node_connectivity_view.nbNode(inode) + 1;
    };
  }

  ax::Scanner<Int32> scanner;
  scanner.exclusiveSum(queue, tmp_row, m_csr_matrix.m_matrix_row);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE
void FemModule::_addValueToGlobalMatrixGpu(Int32 begin, Int32 end, Int32 col, ax::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> in_out_col_csr, ax::NumArrayView<DataViewGetterSetter<Real>, MDDim1, DefaultLayout> in_out_val_csr, Real x)
{

  // Find the right index in the csr matrix
  for (; begin < end; ++begin) {
    const Int32 currentCol = in_out_col_csr[begin];

    if (currentCol == -1 || currentCol == col) {
      in_out_col_csr[begin] = col;
      in_out_val_csr[begin] += x;
      break;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_assembleBuildLessCsrBilinearOperatorTria3()
{
  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_CsrBuildLess");

  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _buildMatrixGpuBuildLessCsr();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_out_col_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_column);
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
      Int32 row = node_dof.dofId(inode, 0).localId();
      Int32 begin = in_row_csr[row];
      Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];
      for (NodeLocalId node2 : cnc.nodes(cell)) {
        if (nodes_infos.isOwn(inode)) {
          Real x = b_matrix[inode_index * 2] * b_matrix[i * 2] + b_matrix[inode_index * 2 + 1] * b_matrix[i * 2 + 1];
          x = x * area;

          Int32 col = node_dof.dofId(node2, 0).localId();
          _addValueToGlobalMatrixGpu(begin, end, col, in_out_col_csr, in_out_val_csr, x);
        }
        i++;
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_assembleBuildLessCsrBilinearOperatorTetra4()
{
  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_CsrBuildLess");

  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _buildMatrixGpuBuildLessCsr();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_out_col_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_column);
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
      Int32 row = node_dof.dofId(inode, 0).localId();
      Int32 begin = in_row_csr[row];
      Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];
      for (NodeLocalId node2 : cnc.nodes(cell)) {

        if (nodes_infos.isOwn(inode)) {
          Real x = b_matrix[inode_index * 3] * b_matrix[i * 3] + b_matrix[inode_index * 3 + 1] * b_matrix[i * 3 + 1] + b_matrix[inode_index * 3 + 2] * b_matrix[i * 3 + 2];
          x = x * volume;

          Int32 col = node_dof.dofId(node2, 0).localId();
          _addValueToGlobalMatrixGpu(begin, end, col, in_out_col_csr, in_out_val_csr, x);
        }
        i++;
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
