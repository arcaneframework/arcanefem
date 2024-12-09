// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrGpuBiliAssembly.hxx                                      (C) 2022-2024 */
/*                                                                           */
/* Methods of the bilinear assembly phase using the csr data structure       */
/* which handle the parallelization on GPU's or CPU's                        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "ArcaneFemFunctionsGpu.h"
#include "FemModule.h"
#include "BSRMatrix.h"
#include "FemUtils.h"
#include "arccore/base/ArccoreGlobal.h"
#include "arccore/base/NotImplementedException.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::_buildOffsets(const SmallSpan<uint>& offsets_smallspan)
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
_buildMatrixCsrGPU()
{
  Int8 mesh_dim = mesh()->dimension();

  Int32 nb_node = nbNode();
  Int32 nb_non_zero = nb_node + 2 * (mesh_dim == 2 ? nbFace() : m_nb_edge);
  m_csr_matrix.initialize(m_dof_family, nb_non_zero, nb_node, m_queue);

  NumArray<uint, MDDim1> offsets_numarray(nb_node + 1);
  SmallSpan<uint> offsets_smallspan = offsets_numarray.to1DSmallSpan();

  // Compute the array of offsets on Gpu
  _buildOffsets(offsets_smallspan);

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
/**
* @brief Assembles the bilinear operator matrix for the FEM linear system with
 * the CSR sparse matrix format for TRIA3 elements.
 *
 * The method performs the following steps:
 *   1. Builds the CSR matrix using _buildMatrixCsrGPU.
 *   2. For each cell, computes element matrix using _computeElementMatrixTRIA3GPU.
 *   3. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE FixedMatrix<3, 3> computeElementMatrixTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real area = Arcane::FemUtils::Gpu::MeshOperation::computeAreaTria3(cell_lid, cn_cv, in_node_coord);

  Real3 dxU = Arcane::FemUtils::Gpu::MeshOperation::computeGradientXTria3(cell_lid, cn_cv, in_node_coord);
  Real3 dyU = Arcane::FemUtils::Gpu::MeshOperation::computeGradientYTria3(cell_lid, cn_cv, in_node_coord);

  return area * (dxU ^ dxU) + area * (dyU ^ dyU);
}

void FemModule::
_assembleCsrGPUBilinearOperatorTRIA3()
{

  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_Csr_Gpu");

  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    // Build the csr matrix
    _buildMatrixCsrGPU();

    BSRFormat bsr_format(subDomain()->traceMng(), m_queue, *(mesh()), m_dofs_on_nodes);
    bsr_format.initialize(nbFace());
    bsr_format.computeSparsity();

    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(m_queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    bsr_format.assembleBilinear<3>([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3(cell_lid, cn_cv, in_node_coord); });

    for (auto i = 0; i < m_csr_matrix.m_matrix_column.extent0(); ++i)
      m_csr_matrix.m_matrix_column[i] = bsr_format.m_bsr_matrix.columns()[i];
    for (auto i = 0; i < m_csr_matrix.m_matrix_value.extent0(); ++i)
      m_csr_matrix.m_matrix_value[i] = bsr_format.m_bsr_matrix.values()[i];
    for (auto i = 0; i < m_csr_matrix.m_matrix_row.extent0(); ++i)
      m_csr_matrix.m_matrix_row[i] = bsr_format.m_bsr_matrix.rowIndex()[i];

    return;
  }

  RunQueue* queue = acceleratorMng()->defaultQueue();
  // Boucle sur les mailles déportée sur accélérateur
  auto command = makeCommand(queue);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());
  auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
  Int32 row_csr_size = m_csr_matrix.m_matrix_row.dim1Size();
  auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
  Int32 col_csr_size = m_csr_matrix.m_matrix_column.dim1Size();
  auto in_out_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);
  UnstructuredMeshConnectivityView m_connectivity_view;
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_connectivity_view.setMesh(this->mesh());
  auto cnc = m_connectivity_view.cellNode();
  Arcane::ItemGenericInfoListView nodes_infos(this->mesh()->nodeFamily());

  Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");

  command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
  {

    Real K_e[9] = { 0 };

    _computeElementMatrixTRIA3GPU(icell, cnc, in_node_coord, K_e); // element stifness matrix
    //             # assemble elementary matrix into the global one
    //             # elementary terms are positionned into K according
    //             # to the rank of associated node in the mesh.nodes list
    //             for node1 in elem.nodes:
    //                 inode1=elem.nodes.index(node1) # get position of node1 in nodes list
    //                 for node2 in elem.nodes:
    //                     inode2=elem.nodes.index(node2)
    //                     K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
    Int32 n1_index = 0;
    for (NodeLocalId node1 : cnc.nodes(icell)) {
      Int32 n2_index = 0;
      for (NodeLocalId node2 : cnc.nodes(icell)) {
        // K[node1.rank,node2.rank]=K[node1.rank,node2.rank]+K_e[inode1,inode2]
        double v = K_e[n1_index * 3 + n2_index];
        // m_k_matrix(node1.localId(), node2.localId()) += v;
        if (nodes_infos.isOwn(node1)) {

          Int32 row = node_dof.dofId(node1, 0).localId();
          Int32 col = node_dof.dofId(node2, 0).localId();
          Int32 begin = in_row_csr[row];
          Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];

          while (begin < end) {
            if (in_col_csr[begin] == col) {
              // t is necessary to get the right type for the atomicAdd (but that means that we have more operations ?)
              // The Macro is there to avoid compilation error if not in c++ 20
              ax::doAtomic<ax::eAtomicOperation::Add>(in_out_val_csr(begin), v);
              break;
            }
            begin++;
          }
        }
        ++n2_index;
      }
      ++n1_index;
    }
  };
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the bilinear operator matrix for the FEM linear system with
 * the CSR sparse matrix format for TETRA4 elements.
 *
 * The method performs the following steps:
 *   1. Builds the CSR matrix using _buildMatrixCsrGPU.
 *   2. For each cell, computes element matrix using _computeElementMatrixTETRA4GPU.
 *   3. Assembles global matrix by adding contributions from each cell's element 
 *      matrix to the corresponding entries in the global matrix.
 */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE FixedMatrix<4, 4> computeElementMatrixTetra4(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const Accelerator::VariableNodeReal3InView& in_node_coord)
{
  Real volume = Arcane::FemUtils::Gpu::MeshOperation::computeVolumeTetra4(cell_lid, cn_cv, in_node_coord);

  Real4 dxU = Arcane::FemUtils::Gpu::MeshOperation::computeGradientXTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dyU = Arcane::FemUtils::Gpu::MeshOperation::computeGradientYTetra4(cell_lid, cn_cv, in_node_coord);
  Real4 dzU = Arcane::FemUtils::Gpu::MeshOperation::computeGradientZTetra4(cell_lid, cn_cv, in_node_coord);

  return volume * (dxU ^ dxU) + volume * (dyU ^ dyU) + volume * (dzU ^ dzU);
}

void FemModule::
_assembleCsrGPUBilinearOperatorTETRA4()
{
  Timer::Action timer_bili(m_time_stats, "AssembleBilinearOperator_Csr_Gpu");

  {
    Timer::Action timer_build(m_time_stats, "BuildMatrix");
    _buildMatrixCsrGPU();

    BSRFormat bsr_format(subDomain()->traceMng(), m_queue, *(mesh()), m_dofs_on_nodes);
    bsr_format.initialize(m_nb_edge);
    bsr_format.computeSparsity();

    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto command = makeCommand(m_queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    bsr_format.assembleBilinear<4>([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4(cell_lid, cn_cv, in_node_coord); });

    for (auto i = 0; i < m_csr_matrix.m_matrix_column.extent0(); ++i)
      m_csr_matrix.m_matrix_column[i] = bsr_format.m_bsr_matrix.columns()[i];
    for (auto i = 0; i < m_csr_matrix.m_matrix_value.extent0(); ++i)
      m_csr_matrix.m_matrix_value[i] = bsr_format.m_bsr_matrix.values()[i];
    for (auto i = 0; i < m_csr_matrix.m_matrix_row.extent0(); ++i)
      m_csr_matrix.m_matrix_row[i] = bsr_format.m_bsr_matrix.rowIndex()[i];

    return;
  }

  {
    auto command = makeCommand(m_queue);

    Int32 row_csr_size = m_csr_matrix.m_matrix_row.extent0();
    Int32 col_csr_size = m_csr_matrix.m_matrix_column.extent0();

    auto in_row_csr = ax::viewIn(command, m_csr_matrix.m_matrix_row);
    auto in_col_csr = ax::viewIn(command, m_csr_matrix.m_matrix_column);
    auto inout_val_csr = ax::viewInOut(command, m_csr_matrix.m_matrix_value);

    auto in_node_coord = ax::viewIn(command, m_node_coord);

    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cell_node_connectivity_view = m_connectivity_view.cellNode();

    ItemGenericInfoListView nodes_infos(mesh()->nodeFamily());

    Timer::Action timer_add_compute(m_time_stats, "AddAndCompute");
    ax::ProfileRegion ps_region(m_queue, "AddAndComputeBilinearTetra4", 0x00FF7F);
    command << RUNCOMMAND_ENUMERATE(Cell, icell, allCells())
    {

      Real K_e[16]{};
      _computeElementMatrixTETRA4GPU(icell, cell_node_connectivity_view, in_node_coord, K_e);

      Int32 node1_idx_in_cell = 0;
      for (NodeLocalId node1_id : cell_node_connectivity_view.nodes(icell)) {

        Int32 node2_idx_in_cell = 0;
        for (NodeLocalId node2_id : cell_node_connectivity_view.nodes(icell)) {

          if (nodes_infos.isOwn(node1_id)) {
            double v = K_e[node1_idx_in_cell * 4 + node2_idx_in_cell];

            Int32 row = node1_id.localId();
            Int32 col = node2_id.localId();
            Int32 begin = in_row_csr[row];

            Int32 end = (row == row_csr_size - 1) ? col_csr_size : in_row_csr[row + 1];

            while (begin < end) {
              if (in_col_csr[begin] == col) {
                ax::doAtomic<ax::eAtomicOperation::Add>(inout_val_csr(begin), v);
                break;
              }
              begin++;
            }
          }
          ++node2_idx_in_cell;
        }
        ++node1_idx_in_cell;
      }
    };
  }
  // m_csr_matrix.printMatrix("csr_matrix_dump.txt");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
