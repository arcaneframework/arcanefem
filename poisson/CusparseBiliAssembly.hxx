// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CsrGpuiBiliAssembly.hxx                                     (C) 2022-2023 */
/*                                                                           */
/* Methods of the bilinear assembly phase using Cusparse library             */
/* which build the global matrix by merging local ones                       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* USE_CUSPARSE_ADD methods                                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
printCsrMatrix(std::string fileName, cusparseCsr csr, bool is_coo)
{
  ofstream file(fileName);
  file << "size :" << csr.nnz << "\n";
  for (auto i = 0; i < (is_coo ? csr.nnz : nbNode()); i++) {
    file << csr.csrRow[i] << " ";
  }
  file << "\n";
  for (auto i = 0; i < csr.nnz; i++) {
    file << csr.csrCol[i] << " ";
  }
  file << "\n";
  for (auto i = 0; i < csr.nnz; i++) {
    file << csr.csrVal[i] << " ";
  }
  file.close();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_computeCusparseElementMatrix(cusparseCsr& result, cusparseCsr& global, Cell cell, cusparseHandle_t handle, IndexedNodeDoFConnectivityView node_dof,
                              computeTimer& timer)
{

  Timer::Action timer_action(m_time_stats, "ComputeCusparseElementMatrix");

  std::chrono::_V2::system_clock::time_point compute_start;
  if (m_register_time) {
    compute_start = std::chrono::high_resolution_clock::now();
  }
  /*-------------------------------------------------------------------------------------------------------------------------------*/
  //First part : compute element matrix
  // Get coordiantes of the triangle element  TRI3
  //------------------------------------------------
  //                  0 o
  //                   . .
  //                  .   .
  //                 .     .
  //              1 o . . . o 2
  //------------------------------------------------
  Real3 m0 = m_node_coord[cell.nodeId(0)];
  Real3 m1 = m_node_coord[cell.nodeId(1)];
  Real3 m2 = m_node_coord[cell.nodeId(2)];

  Real area = _computeAreaTriangle3(cell); // calculate area

  Real2 dPhi0(m1.y - m2.y, m2.x - m1.x);
  Real2 dPhi1(m2.y - m0.y, m0.x - m2.x);
  Real2 dPhi2(m0.y - m1.y, m1.x - m0.x);

  FixedMatrix<2, 3> b_matrix;
  b_matrix(0, 0) = dPhi0.x;
  b_matrix(0, 1) = dPhi1.x;
  b_matrix(0, 2) = dPhi2.x;

  b_matrix(1, 0) = dPhi0.y;
  b_matrix(1, 1) = dPhi1.y;
  b_matrix(1, 2) = dPhi2.y;

  b_matrix.multInPlace(1.0 / (2.0 * area));

  FixedMatrix<3, 3> int_cdPi_dPj = matrixMultiplication(matrixTranspose(b_matrix), b_matrix);
  int_cdPi_dPj.multInPlace(area);

  std::chrono::_V2::system_clock::time_point convert_coo;
  if (m_register_time) {
    auto compute_el_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> compute_el_duration = compute_el_stop - compute_start;
    timer.compute_el += compute_el_duration.count();
    convert_coo = std::chrono::high_resolution_clock::now();
  }
  /*-------------------------------------------------------------------------------------------------------------------------------*/
  //Second part : putting the matrix in COO format (might want to optimsie that part by doing it earlier) before converting it to csr

  //Must change int_cdPi_dPj in a COO matrix (before converting it to csr);
  void* row_void;
  CHECK_CUDA(cudaMallocManaged(&row_void, 9 * sizeof(Int32)));
  Int32* row_indexes = (Int32*)row_void;
  void* col_void;
  CHECK_CUDA(cudaMallocManaged(&col_void, 9 * sizeof(Int32)));
  Int32* col_indexes = (Int32*)col_void;
  void* vals_void;
  CHECK_CUDA(cudaMallocManaged(&vals_void, 9 * sizeof(float)));
  float* vals = (float*)vals_void;

  cusparseMatDescr_t local_mat;
  CHECK_CUSPARSE(cusparseCreateMatDescr(&local_mat));
  cusparseCsr local;
  local.desc = local_mat;
  local.csrRow = row_indexes;
  local.csrCol = col_indexes;
  local.csrVal = (float*)vals;
  local.nnz = 9;

  int i = 0;
  int j = 0;
  for (NodeLocalId node1 : cell.nodes()) {
    j = 0;
    for (NodeLocalId node2 : cell.nodes()) {
      vals[i * 3 + j] = int_cdPi_dPj(i, j);
      row_indexes[i * 3 + j] = node_dof.dofId(node1, 0);
      col_indexes[i * 3 + j] = node_dof.dofId(node2, 0);
      j++;
    }
    i++;
  }

  std::chrono::_V2::system_clock::time_point sort_coo;
  std::chrono::duration<double> convert_coo_time;
  if (m_register_time) {
    auto convert_coo_stop = std::chrono::high_resolution_clock::now();
    convert_coo_time = convert_coo_stop - convert_coo;
    timer.convert_coo += convert_coo_time.count();
    sort_coo = std::chrono::high_resolution_clock::now();
  }
  //Sorting of the COO values with an insertion sort
  Int32 rj = 0;
  Int32 cj = 0;
  float vj = 0;
  for (i = 1; i < 9; i++) {
    rj = row_indexes[i];
    cj = col_indexes[i];
    vj = vals[i];
    j = i - 1;
    while (j >= 0 && row_indexes[j] > rj) {
      row_indexes[j + 1] = row_indexes[j];
      col_indexes[j + 1] = col_indexes[j];
      vals[j + 1] = vals[j];
      j--;
    }
    row_indexes[j + 1] = rj;
    col_indexes[j + 1] = cj;
    vals[j + 1] = vj;
    Int32 k = j - 1;
    Int32 rk, ck;
    float vk;
    if (j > 0) {
      rk = row_indexes[j];
      ck = col_indexes[j];
      vk = vals[j];
      while (k >= 0 && row_indexes[k] == rk && col_indexes[k] > ck) {
        col_indexes[k + 1] = col_indexes[k];
        vals[k + 1] = vals[k];
        k--;
      }
      col_indexes[k + 1] = ck;
      vals[k + 1] = vk;
    }
  }
  std::chrono::_V2::system_clock::time_point convert_csr;
  if (m_register_time) {
    auto sort_coo_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sort_coo_time = sort_coo_stop - sort_coo;
    std::chrono::duration<double> convert_coo_tot = sort_coo_stop - convert_coo;
    timer.sort_coo += sort_coo_time.count();
    timer.convert_coo_tot += convert_coo_tot.count();
    convert_csr = std::chrono::high_resolution_clock::now();
  }

  //conversion from COO to CSR
  void* csrRowPtr_void;
  CHECK_CUDA(cudaMallocManaged(&csrRowPtr_void, (nbNode() + 1) * sizeof(Int32)));
  Int32* csrRowPtr = (Int32*)csrRowPtr_void;
  CHECK_CUSPARSE(cusparseXcoo2csr(handle, row_indexes, 9, nbNode(), csrRowPtr, CUSPARSE_INDEX_BASE_ZERO));
  local.csrRow = csrRowPtr;
  CHECK_CUDA(cudaFree(row_indexes));

  std::chrono::_V2::system_clock::time_point adding_global;
  if (m_register_time) {
    auto convert_csr_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> convert_csr_time = convert_csr_stop - convert_csr;
    std::chrono::duration<double> convert_tot = convert_csr_stop - convert_coo;
    timer.convert_csr_tot += convert_coo_time.count();
    timer.convert_tot += convert_tot.count();
    adding_global = std::chrono::high_resolution_clock::now();
  }
  /*-------------------------------------------------------------------------------------------------------------------------------*/
  // Third part : adding the local and global, storing result in the res

  //Adding the CSR local matrix to the global one using cusparsecsrgeam2
  //see https://docs.nvidia.com/cuda/cusparse/index.html?highlight=cusparseScsrgeam#cusparse-t-csrgeam2 for the example code
  Int32 baseC,
  nnzC;
  size_t bufferSizeInBytes;
  char* buffer = NULL;
  void* buffer_void = NULL;
  Int32* nnzTotalDevHostPtr = &nnzC;
  CHECK_CUSPARSE(cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST));
  //Int32* csrRowPtrC;
  //CHECK_CUDA(cudaMallocManaged(&csrRowPtrC, nbNode() + 1 * sizeof(Int32)));
  float alpha = 1.0;
  float beta = 1.0;
  Int32 m = nbNode();
  CHECK_CUSPARSE(cusparseScsrgeam2_bufferSizeExt(handle, m, m,
                                                 &alpha,
                                                 global.desc, global.nnz,
                                                 global.csrVal, global.csrRow, global.csrCol,
                                                 &beta,
                                                 local.desc, local.nnz,
                                                 local.csrVal, local.csrRow, local.csrCol,
                                                 result.desc,
                                                 result.csrVal, result.csrRow, result.csrCol, &bufferSizeInBytes));
  CHECK_CUDA(cudaMallocManaged(&buffer_void, bufferSizeInBytes * sizeof(char)));
  buffer = (char*)buffer_void;

  CHECK_CUSPARSE(cusparseXcsrgeam2Nnz(handle, m, m,
                                      local.desc, local.nnz, local.csrRow, local.csrCol,
                                      global.desc, global.nnz, global.csrRow, global.csrCol,
                                      result.desc, result.csrRow, nnzTotalDevHostPtr,
                                      buffer));
  if (NULL != nnzTotalDevHostPtr) {
    nnzC = *nnzTotalDevHostPtr;
  }
  else {
    CHECK_CUDA(cudaMemcpy(&nnzC, result.csrRow + m, sizeof(int), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(&baseC, result.csrRow, sizeof(int), cudaMemcpyDeviceToHost));
    nnzC -= baseC;
  }
  result.nnz = nnzC;
  void* res_col_void;
  void* res_val_void;

  CHECK_CUDA(cudaMallocManaged(&res_col_void, sizeof(Int32) * nnzC));
  CHECK_CUDA(cudaMallocManaged(&res_val_void, sizeof(float) * nnzC));
  result.csrCol = (Int32*)res_col_void;
  result.csrVal = (float*)res_val_void;
  CHECK_CUSPARSE(cusparseScsrgeam2(handle, m, m,
                                   &alpha,
                                   local.desc, local.nnz,
                                   local.csrVal, local.csrRow, local.csrCol,
                                   &beta,
                                   global.desc, global.nnz,
                                   global.csrVal, global.csrRow, global.csrCol,
                                   result.desc,
                                   result.csrVal, result.csrRow, result.csrCol,
                                   buffer));

  CHECK_CUDA(cudaFree(buffer));
  CHECK_CUDA(cudaFree(local.csrVal));
  CHECK_CUDA(cudaFree(local.csrCol));
  CHECK_CUDA(cudaFree(local.csrRow));
  CHECK_CUSPARSE(cusparseDestroyMatDescr(local.desc));

  if (m_register_time) {
    auto adding_global_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> adding_tot = adding_global_stop - adding_global;
    std::chrono::duration<double> compute_tot = adding_global_stop - compute_start;
    timer.add_glob += adding_tot.count();
    timer.compute_tot += compute_tot.count();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * @brief Assemble Bilinear TRIA3 with cusparse help. It only works for p=1 since there is
 * one node per Edge.
 * 
 * 
 */
void FemModule::
_assembleCusparseBilinearOperatorTRIA3()
{

  Timer::Action timer_action(m_time_stats, "AssembleCusparseBilinearOperator");

  computeTimer t = {};
  std::chrono::_V2::system_clock::time_point lhs_s;
  if (m_register_time) {
    logger << "-------------------------------------------------------------------------------------\n"
           << "Using Cusparse for Bilinear assembly\n";
    lhs_s = std::chrono::high_resolution_clock::now();
  }
  //Initialization of the CSR matrix;
  //This formula only works in p=1
  std::chrono::_V2::system_clock::time_point cuda_init_start;
  if (m_register_time) {
    cuda_init_start = std::chrono::high_resolution_clock::now();
  }
  CHECK_CUDA(cudaFree(0));

  double cuda_init = 0;
  std::chrono::_V2::system_clock::time_point cuda_init_stop;
  if (m_register_time) {
    cuda_init_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> cuda_init_time = cuda_init_stop - cuda_init_start;
    cuda_init = cuda_init_time.count();
  }

  Int32 nnz = nbFace() * 2 + nbNode();
  cusparseHandle_t handle;
  CHECK_CUSPARSE(cusparseCreate(&handle));
  //Initialize the global matrix. Everything is in the unified memory
  void* res1_row_void;
  void* res2_row_void;
  CHECK_CUDA(cudaMallocManaged(&res1_row_void, sizeof(Int32) * (nbNode() + 1)));
  CHECK_CUDA(cudaMemset(res1_row_void, 0, sizeof(Int32) * (nbNode() + 1)));
  Int32* res1_row = (Int32*)res1_row_void;
  CHECK_CUDA(cudaMallocManaged(&res2_row_void, sizeof(Int32) * (nbNode() + 1)));
  CHECK_CUDA(cudaMemset(res2_row_void, 0, sizeof(Int32) * (nbNode() + 1)));
  Int32* res2_row = (Int32*)res2_row_void;

  //The number of Node must be changed when p != 1

  //init result matrix
  cusparseCsr res1;
  cusparseCsr res2;

  cusparseMatDescr_t res1_desc;
  cusparseMatDescr_t res2_desc;
  CHECK_CUSPARSE(cusparseCreateMatDescr(&res1_desc));
  CHECK_CUSPARSE(cusparseCreateMatDescr(&res2_desc));

  res1.desc = res1_desc;
  res2.desc = res2_desc;
  res1.csrRow = res1_row;
  res2.csrRow = res2_row;
  res1.csrCol = NULL;
  res2.csrCol = NULL;
  res1.csrVal = NULL;
  res2.csrVal = NULL;

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  Int32 i = 0;

  ENUMERATE_CELL (icell, allCells()) {
    Cell cell = *icell;

    if (i % 2 == 0) {

      CHECK_CUDA(cudaFree(res1.csrCol));
      CHECK_CUDA(cudaFree(res1.csrVal));
      _computeCusparseElementMatrix(res1, res2, cell, handle, node_dof,
                                    t);
    }
    //computation of the local matrix and adding it in the global one
    else {
      CHECK_CUDA(cudaFree(res2.csrCol));
      CHECK_CUDA(cudaFree(res2.csrVal));
      _computeCusparseElementMatrix(res2, res1, cell, handle, node_dof,
                                    t);
    }
    i++;
  }
  /*
  if (nbNode() % 2 == 0)
    printCsrMatrix("csrTest.txt", res1, false);
  else
    printCsrMatrix("csrTest.txt", res2, false);

*/
  CHECK_CUSPARSE(cusparseDestroyMatDescr(res1.desc));
  CHECK_CUSPARSE(cusparseDestroyMatDescr(res2.desc));

  //Must get the results into a linear format
  CHECK_CUDA(cudaFree(res1.csrRow));
  CHECK_CUDA(cudaFree(res1.csrCol));
  CHECK_CUDA(cudaFree(res1.csrVal));
  CHECK_CUDA(cudaFree(res2.csrRow));
  CHECK_CUDA(cudaFree(res2.csrCol));
  CHECK_CUDA(cudaFree(res2.csrVal));
  CHECK_CUSPARSE(cusparseDestroy(handle));
  if (m_register_time) {
    auto lhs_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = lhs_end - cuda_init_stop;
    double lhs_loc_time = duration.count();
    logger << "Average Time to compute element matrix : " << t.compute_el / nbCell() << "\n"
           << "Total Time to compute element matrix : " << t.compute_el << "\n"
           << "Percentage time to compute element matrix : " << t.compute_el / lhs_loc_time * 100 << "%\n"
           << "Average Time to convert to coo : " << t.convert_coo / nbCell() << "\n"
           << "Total Time to convert to coo : " << t.convert_coo << "\n"
           << "Percentage Time to convert to coo : " << t.convert_coo / lhs_loc_time * 100 << "%\n"
           << "Average Time to sort the coo : " << t.sort_coo / nbCell() << "\n"
           << "Total Time to sort the coo : " << t.sort_coo << "\n"
           << "Percentage Time to sort the coo : " << t.sort_coo / lhs_loc_time * 100 << "%\n"
           << "Average Time to convert and sort to coo : " << t.convert_coo_tot / nbCell() << "\n"
           << "Total Time to convert and sort to coo : " << t.convert_coo_tot << "\n"
           << "Percentage Time to convert and sort to coo : " << t.convert_coo_tot / lhs_loc_time * 100 << "%\n"
           << "Average Time to convert to csr : " << t.convert_csr_tot / nbCell() << "\n"
           << "Total Time to convert to csr : " << t.convert_csr_tot << "\n"
           << "Percentage Time to convert to csr : " << t.convert_csr_tot / lhs_loc_time * 100 << "%\n"
           << "Average Time to convert the computed matrix : " << t.convert_tot / nbCell() << "\n"
           << "Total Time to convert the computed matrix : " << t.convert_tot << "\n"
           << "Percentage Time to convert the computed matrix : " << t.convert_tot / lhs_loc_time * 100 << "%\n"
           << "Average Time to add to the global matrix : " << t.add_glob / nbCell() << "\n"
           << "Total Time to add to the global matrix : " << t.add_glob << "\n"
           << "Percentage Time to add to the global matrix : " << t.add_glob / lhs_loc_time * 100 << "%\n"
           << "Average Time to make the computation operation : " << t.compute_tot / nbCell() << "\n"
           << "Total Time to make the computation operation : " << t.compute_tot << "\n"
           << "Percentage Time to make the computation operation : " << t.compute_tot / lhs_loc_time * 100 << "%\n"
           << "Total time for the lhs computation : " << lhs_loc_time << "\n"
           << "Total time of the cuda init : " << cuda_init << "\n"
           << "Total time of lhs with the init : " << cuda_init + lhs_loc_time << "\n"
           << "-------------------------------------------------------------------------------------\n\n";
    wbuild << lhs_loc_time << ",";
    lhs_time += lhs_loc_time;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/