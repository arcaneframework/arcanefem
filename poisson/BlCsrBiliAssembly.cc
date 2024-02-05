// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BlcsrBiliAssembly.hxx                                     (C) 2022-2023   */
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
  Int32 nnz = nbFace() * 2 + nbnde;
  m_csr_matrix.initialize(m_dof_family, nnz, nbnde);

  Integer index = 1;
  m_csr_matrix.m_matrix_row(0) = 0;
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
  Int32 nnz = nbFace() * 2 + nbnde;

  NumArray<Int32, MDDim1> tmp_row;
  tmp_row.resize(nbnde);
  m_csr_matrix.initialize(m_dof_family, nnz, nbnde);

  RunQueue* queue = acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);
  auto in_out_tmp_row = ax::viewInOut(command, tmp_row);
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  UnstructuredMeshConnectivityView connectivity_view;
  connectivity_view.setMesh(this->mesh());
  auto nfc = connectivity_view.nodeFace();

  command << RUNCOMMAND_ENUMERATE(Node, inode, allNodes())
  {
    Int64 index = node_dof.dofId(inode, 0).localId();
    in_out_tmp_row[index] = nfc.nbFace(inode) + 1;
  };
  ax::Scanner<Int32> scanner;
  scanner.exclusiveSum(queue, tmp_row, m_csr_matrix.m_matrix_row);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE
Real FemModule::_computeCellMatrixGpuTETRA4(CellLocalId icell, IndexedCellNodeConnectivityView cnc, ax::VariableNodeReal3InView in_node_coord, Real b_matrix[12])
{
  Real3 m0 = in_node_coord[cnc.nodeId(icell, 0)];
  Real3 m1 = in_node_coord[cnc.nodeId(icell, 1)];
  Real3 m2 = in_node_coord[cnc.nodeId(icell, 2)];
  Real3 m3 = in_node_coord[cnc.nodeId(icell, 3)];

  // Calculate vectors representing edges of the tetrahedron
  Real3 v0 = m1 - m0;
  Real3 v1 = m2 - m0;
  Real3 v2 = m3 - m0;

  // Compute volume using scalar triple product
  Real area = std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2))) / 6.0;

  // Compute gradients of shape functions
  Real3 dPhi0 = Arcane::math::cross(m2 - m1, m1 - m3) ;
  Real3 dPhi1 = Arcane::math::cross(m3 - m0, m0 - m2) ;
  Real3 dPhi2 = Arcane::math::cross(m1 - m0, m0 - m3) ;
  Real3 dPhi3 = Arcane::math::cross(m0 - m1, m1 - m2) ;

  // Construct the B-matrix as a vector
  Real mul = 1.0 / (6.0 * area);
  b_matrix[0] = dPhi0.x * mul;
  b_matrix[1] = dPhi0.y * mul;
  b_matrix[2] = dPhi0.z * mul;

  b_matrix[3] = dPhi1.x * mul;
  b_matrix[4] = dPhi1.y * mul;
  b_matrix[5] = dPhi1.z * mul;

  b_matrix[6] = dPhi2.x * mul;
  b_matrix[7] = dPhi2.y * mul;
  b_matrix[8] = dPhi2.z * mul;

  b_matrix[9]  = dPhi3.x * mul;
  b_matrix[10] = dPhi3.y * mul;
  b_matrix[11] = dPhi3.z * mul;

  /*
  FixedMatrix<3, 4> b_matrix;
  b_matrix(0, 0) = dPhi0.x;
  b_matrix(1, 0) = dPhi0.y;
  b_matrix(2, 0) = dPhi0.z;

  b_matrix(0, 1) = dPhi1.x;
  b_matrix(1, 1) = dPhi1.y;
  b_matrix(2, 1) = dPhi1.z;

  b_matrix(0, 2) = dPhi2.x;
  b_matrix(1, 2) = dPhi2.y;
  b_matrix(2, 2) = dPhi2.z;

  b_matrix(0, 3) = dPhi3.x;
  b_matrix(1, 3) = dPhi3.y;
  b_matrix(2, 3) = dPhi3.z;

  b_matrix.multInPlace(1.0 / (6.0 * volume));

  // Compute the element matrix
  FixedMatrix<4, 4> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(volume);
*/

  return area;
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE
Real FemModule::_computeCellMatrixGpuTRIA3(CellLocalId icell, IndexedCellNodeConnectivityView cnc, ax::VariableNodeReal3InView in_node_coord, Real b_matrix[6])
{
  Real3 m0 = in_node_coord[cnc.nodeId(icell, 0)];
  Real3 m1 = in_node_coord[cnc.nodeId(icell, 1)];
  Real3 m2 = in_node_coord[cnc.nodeId(icell, 2)];

  Real area = 0.5 * ((m1.x - m0.x) * (m2.y - m0.y) - (m2.x - m0.x) * (m1.y - m0.y));//_computeAreaTriangle3Gpu(icell, cnc, in_node_coord);

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  Real mul = 0.5 / area;
  b_matrix[0] = dPhi0.x * mul;
  b_matrix[1] = dPhi0.y * mul;

  b_matrix[2] = dPhi1.x * mul;
  b_matrix[3] = dPhi1.y * mul;

  b_matrix[4] = dPhi2.x * mul;
  b_matrix[5] = dPhi2.y * mul;

  return area;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE
void FemModule::_addValueToGlobalMatrixTria3Gpu(Int32 begin, Int32 end, Int32 col, ax::NumArrayView<DataViewGetterSetter<Int32>, MDDim1, DefaultLayout> in_out_col_csr, ax::NumArrayView<DataViewGetterSetter<Real>, MDDim1, DefaultLayout> in_out_val_csr, Real x)
{

/*
  // Find the right index in the csr matrix
  while (begin < end) {
    if (in_out_col_csr[begin] == -1) {
      in_out_col_csr[begin] = col;
      in_out_val_csr[begin] = x;
      break;
    }
    else if (in_out_col_csr[begin] == col) {
      in_out_val_csr[begin] += x;
      break;
    }
    begin++;
  }
*/

  // Find the right index in the csr matrix
  for (; begin < end; ++begin) {
    const Int32 currentCol = in_out_col_csr[begin];

    if (currentCol == -1 || currentCol == col) {
        in_out_col_csr[begin]  = col;
        in_out_val_csr[begin] += x;
        break;
    }
  }

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_assembleBuildLessCsrBilinearOperatorTria3()
{
  Timer::Action timer_blcsr_bili(m_time_stats, "AssembleBuildLessCsrBilinearOperatorTria3");

  {
    Timer::Action timer_blcsr_build(m_time_stats, "BuildLessCsrBuildMatrixGPU");
    // Build only the row part of the csr matrix on GPU
    // Using scan -> might be improved
    _buildMatrixGpuBuildLessCsr();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();

  // Boucle sur les noeuds déportée sur accélérateur
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


  Timer::Action timer_blcsr_add_compute(m_time_stats, "BuildLessCsrAddAndCompute");

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

/*
      if (inode == cnc.nodeId(cell, 1)) {
        inode_index = 1;
      }
      else if (inode == cnc.nodeId(cell, 2)) {
        inode_index = 2;
      }
      else {
        inode_index = 0;
      }
*/

      Real b_matrix[6] = { 0 };
      Real area = _computeCellMatrixGpuTRIA3(cell, cnc, in_node_coord, b_matrix);

      Int32 i = 0;
      Int32 row = node_dof.dofId(inode, 0).localId();
      Int32 begin = in_row_csr[row];
      Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];
      for (NodeLocalId node2 : cnc.nodes(cell)) {
/*
        Real x = 0.0;
        for (Int32 k = 0; k < 2; k++) {
          x += b_matrix[inode_index * 2 + k] * b_matrix[i * 2 + k];
        }
*/
        Real x = b_matrix[inode_index * 2]     * b_matrix[i * 2] +
                 b_matrix[inode_index * 2 + 1] * b_matrix[i * 2 + 1];

        x = x * area;
        if (nodes_infos.isOwn(inode)) {

          Int32 col = node_dof.dofId(node2, 0).localId();

/*
          if (row == row_csr_size - 1) {
            end = col_csr_size;
          }
          else {
            end = in_row_csr[row + 1];
          }
*/
          _addValueToGlobalMatrixTria3Gpu(begin, end, col, in_out_col_csr, in_out_val_csr, x);
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
  Timer::Action timer_blcsr_bili(m_time_stats, "AssembleBuildLessCsrBilinearOperatorTetra4");

  {
    Timer::Action timer_blcsr_build(m_time_stats, "BuildLessCsrBuildMatrixGPU");
    // Build only the row part of the csr matrix on GPU
    // Using scan -> might be improved
    _buildMatrixGpuBuildLessCsr();
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();

  // Boucle sur les noeuds déportée sur accélérateur
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


  Timer::Action timer_blcsr_add_compute(m_time_stats, "BuildLessCsrAddAndCompute");

  command << RUNCOMMAND_ENUMERATE(Node, inode, allNodes())
  {
    Int32 inode_index = 0;
    for (auto cell : ncc.cells(inode)) {

      // How can I know the right index ?
      // By checking in the global id ?
      // Working currently, but maybe only because p = 1 ?

      if (inode == cnc.nodeId(cell, 1)) {
        inode_index = 1;
      }
      else if (inode == cnc.nodeId(cell, 2)) {
        inode_index = 2;
      }
      else if (inode == cnc.nodeId(cell, 3)) {
        inode_index = 3;
      }
      else {
        inode_index = 0;
      }

      Real b_matrix[12] = { 0 };
      Real area = _computeCellMatrixGpuTETRA4(cell, cnc, in_node_coord, b_matrix);

      Int32 i = 0;
      Int32 row = node_dof.dofId(inode, 0).localId();
      Int32 begin = in_row_csr[row];
      Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];
      for (NodeLocalId node2 : cnc.nodes(cell)) {
/*
        Real x = 0.0;
        for (Int32 k = 0; k < 2; k++) {
          x += b_matrix[inode_index * 2 + k] * b_matrix[i * 2 + k];
        }
*/
        Real x = b_matrix[inode_index * 3]     * b_matrix[i * 3]     +
                 b_matrix[inode_index * 3 + 1] * b_matrix[i * 3 + 1] +
                 b_matrix[inode_index * 3 + 2] * b_matrix[i * 3 + 2];

        x = x * area;
        if (nodes_infos.isOwn(inode)) {

          Int32 col = node_dof.dofId(node2, 0).localId();

/*
          if (row == row_csr_size - 1) {
            end = col_csr_size;
          }
          else {
            end = in_row_csr[row + 1];
          }
*/
          _addValueToGlobalMatrixTria3Gpu(begin, end, col, in_out_col_csr, in_out_val_csr, x);
        }
        i++;
      }
    }
  };
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
