// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* HypreDoFLinearSystem.cc                                     (C) 2022-2024 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "DoFLinearSystem.h"

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/NumArray.h>
#include <arcane/utils/PlatformUtils.h>

#include <arcane/core/VariableTypes.h>
#include <arcane/core/IItemFamily.h>
#include <arcane/core/BasicService.h>
#include <arcane/core/ServiceFactory.h>

#include <arcane/accelerator/core/Runner.h>

#include "FemUtils.h"
#include "IDoFLinearSystemFactory.h"

#include "HypreDoFLinearSystemFactory_axl.h"

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <krylov.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
using namespace Arcane;
namespace
{
  inline void
  check(const char* hypre_func, int error_code)
  {
    if (error_code == 0)
      return;
    char buf[8192];
    HYPRE_DescribeError(error_code, buf);
    cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
         << "\nHYPRE ERROR in function "
         << hypre_func
         << "\nError_code=" << error_code
         << "\nMessage=" << buf
         << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
         << "\n"
         << std::flush;
    throw Exception("HYPRE Check", hypre_func);
  }
  inline void
  hypreCheck(const char* hypre_func, int error_code)
  {
    check(hypre_func, error_code);
    int r = HYPRE_GetError();
    if (r != 0)
      cout << "HYPRE GET ERROR r=" << r
           << " error_code=" << error_code << " func=" << hypre_func << '\n';
  }
} // namespace

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class HypreDoFLinearSystemImpl
: public TraceAccessor
, public DoFLinearSystemImpl
{
 public:

  HypreDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
  : TraceAccessor(dof_family->traceMng())
  , m_dof_family(dof_family)
  , m_rhs_variable(VariableBuildInfo(dof_family, solver_name + "RHSVariable"))
  , m_dof_variable(VariableBuildInfo(dof_family, solver_name + "SolutionVariable"))
  , m_dof_matrix_indexes(VariableBuildInfo(m_dof_family, solver_name + "DoFMatrixIndexes"))
  , m_dof_elimination_info(VariableBuildInfo(m_dof_family, solver_name + "DoFEliminationInfo"))
  , m_dof_elimination_value(VariableBuildInfo(m_dof_family, solver_name + "DoFEliminationValue"))
  {
    info() << "Creating HypreDoFLinearSystemImpl()";
  }

  ~HypreDoFLinearSystemImpl()
  {
  }

 public:

  void build()
  {
  }

 public:

  void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    ARCANE_THROW(NotImplementedException, "");
  }

  void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value) override
  {
    ARCANE_THROW(NotImplementedException, "");
  }

  void eliminateRow(DoFLocalId row, Real value) override
  {
    ARCANE_THROW(NotImplementedException, "");
  }

  void eliminateRowColumn(DoFLocalId row, Real value) override
  {
    ARCANE_THROW(NotImplementedException, "");
  }

  void solve() override;

  VariableDoFReal& solutionVariable() override
  {
    return m_dof_variable;
  }

  VariableDoFReal& rhsVariable() override
  {
    return m_rhs_variable;
  }

  void setSolverCommandLineArguments(const CommandLineArguments& args) override
  {
  }

  void clearValues()
  {
    info() << "Clear values";
    m_csr_view = {};
  }

  void setCSRValues(const CSRFormatView& csr_view) override
  {
    m_csr_view = csr_view;
  }
  bool hasSetCSRValues() const override { return true; }

  void setRunner(Runner* r) override { m_runner = r; }
  Runner* runner() const { return m_runner; }

 private:

  IItemFamily* m_dof_family = nullptr;
  VariableDoFReal m_rhs_variable;
  VariableDoFReal m_dof_variable;
  VariableDoFInt32 m_dof_matrix_indexes;
  VariableDoFByte m_dof_elimination_info;
  VariableDoFReal m_dof_elimination_value;
  Runner* m_runner = nullptr;

  CSRFormatView m_csr_view;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void HypreDoFLinearSystemImpl::
solve()
{
  HYPRE_Init(); /* must be the first HYPRE function call */

  HYPRE_MemoryLocation hypre_memory = HYPRE_MEMORY_HOST;
  HYPRE_ExecutionPolicy hypre_exec_policy = HYPRE_EXEC_HOST;

  bool is_use_device = false;
  if (m_runner) {
    is_use_device = isAcceleratorPolicy(m_runner->executionPolicy());
    info() << "Runner for Hypre=" << m_runner->executionPolicy() << " is_device=" << is_use_device;
  }
  // Si HYPRE n'est pas compilé avec le support GPU, alors on utilise l'hôte.
  // (NOTE: a priori il n'y a pas besoin de faire cela. Si Hypre n'est pas compilé avec
  // support GPU alors HYPRE_MEMORY_DEVICE <=> HYPRE_MEMORY_HOST
  // TODO: détecter la cohérence entre le GPU de Hypre et le notre (.i.e les deux
  // utilisent CUDA ou ROCM)
#ifndef HYPRE_USING_GPU
  if (is_use_device)
    info() << "Hypre is not compiled with GPU support. Using host backend";
#endif

  if (is_use_device) {
    m_runner->setAsCurrentDevice();
    hypre_memory = HYPRE_MEMORY_DEVICE;
    hypre_exec_policy = HYPRE_EXEC_DEVICE;
  }

  hypreCheck("HYPRE_SetMemoryLocation", HYPRE_SetMemoryLocation(hypre_memory));
  /* setup AMG on GPUs */
  hypreCheck("HYPRE_SetExecutionPolicy", HYPRE_SetExecutionPolicy(hypre_exec_policy));

  if (is_use_device) {
    /* use hypre's SpGEMM instead of vendor implementation */
    HYPRE_SetSpGemmUseVendor(false);
    /* use GPU RNG */
    HYPRE_SetUseGpuRand(true);
  }

  /* use hypre's GPU memory pool */
  //HYPRE_SetGPUMemoryPoolSize(bin_growth, min_bin, max_bin, max_bytes);

  /* setup IJ matrix A */
  // TODO: Utiliser le bon communicateur en parallèle.
  MPI_Comm comm = MPI_COMM_WORLD;

  HYPRE_IJMatrix ij_A = nullptr;
  HYPRE_ParCSRMatrix parcsr_A = nullptr;

  const int first_row = 0;
  const int nb_row = m_csr_view.rows().size();
  info() << "NB_ROW=" << nb_row;
  const int last_row = first_row + nb_row - 1;
  {
    int first_col = first_row;
    int last_col = first_col + nb_row - 1;
    info() << "CreateMatrix row=" << first_row << ", " << last_row
           << " col=" << first_col << ", " << last_col;
    HYPRE_IJMatrixCreate(comm, first_row, last_row, first_col, last_col, &ij_A);
  }

  NumArray<Int32, MDDim1> rows_index(nb_row);
  for (Int32 i = 0; i < nb_row; ++i)
    rows_index[i] = i;
  Span<const Int32> rows_index_span = rows_index.to1DSpan();

  int* rows_nb_column_data = const_cast<int*>(m_csr_view.rowsNbColumn().data());

  Real m1 = platform::getRealTime();
  HYPRE_IJMatrixSetObjectType(ij_A, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize_v2(ij_A, hypre_memory);

  /* GPU pointers; efficient in large chunks */
  HYPRE_IJMatrixSetValues(ij_A,
                          nb_row,
                          rows_nb_column_data,
                          rows_index_span.data(),
                          m_csr_view.columns().data(),
                          m_csr_view.values().data());

  HYPRE_IJMatrixAssemble(ij_A);
  HYPRE_IJMatrixGetObject(ij_A, (void**)&parcsr_A);
  Real m2 = platform::getRealTime();
  info() << "Time to create matrix=" << (m2 - m1);

  //HYPRE_IJMatrixPrint(ij_A, "dumpA.txt");

  HYPRE_IJVector ij_vector_b = nullptr;
  HYPRE_ParVector parvector_b = nullptr;
  HYPRE_IJVector ij_vector_x = nullptr;
  HYPRE_ParVector parvector_x = nullptr;

  hypreCheck("IJVectorCreate", HYPRE_IJVectorCreate(comm, first_row, last_row, &ij_vector_b));
  hypreCheck("IJVectorSetObjectType", HYPRE_IJVectorSetObjectType(ij_vector_b, HYPRE_PARCSR));
  HYPRE_IJVectorInitialize_v2(ij_vector_b, hypre_memory);

  hypreCheck("IJVectorCreate", HYPRE_IJVectorCreate(comm, first_row, last_row, &ij_vector_x));
  hypreCheck("IJVectorSetObjectType", HYPRE_IJVectorSetObjectType(ij_vector_x, HYPRE_PARCSR));
  HYPRE_IJVectorInitialize_v2(ij_vector_x, hypre_memory);

  Real v1 = platform::getRealTime();
  hypreCheck("HYPRE_IJVectorSetValues",
             HYPRE_IJVectorSetValues(ij_vector_b, nb_row, rows_index_span.data(),
                                     m_rhs_variable.asArray().data()));

  hypreCheck("HYPRE_IJVectorSetValues",
             HYPRE_IJVectorSetValues(ij_vector_x, nb_row, rows_index_span.data(),
                                     m_dof_variable.asArray().data()));

  hypreCheck("HYPRE_IJVectorAssemble",
             HYPRE_IJVectorAssemble(ij_vector_b));
  HYPRE_IJVectorGetObject(ij_vector_b, (void**)&parvector_b);

  hypreCheck("HYPRE_IJVectorAssemble",
             HYPRE_IJVectorAssemble(ij_vector_x));
  HYPRE_IJVectorGetObject(ij_vector_x, (void**)&parvector_x);
  Real v2 = platform::getRealTime();
  info() << "Time to create vectors=" << (v2 - v1);

  //HYPRE_IJVectorPrint(ij_vector_b, "dumpB.txt");
  //HYPRE_IJVectorPrint(ij_vector_x, "dumpX.txt");

  HYPRE_Solver solver = nullptr;
  HYPRE_Solver precond = nullptr;
  /* setup AMG */
  HYPRE_ParCSRPCGCreate(comm, &solver);

  /* Set some parameters (See Reference Manual for more parameters) */
  HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
  HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
  HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
  HYPRE_PCGSetPrintLevel(solver, 2); /* print solve info */
  HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

  hypreCheck("HYPRE_BoomerAMGCreate", HYPRE_BoomerAMGCreate(&precond));

  HYPRE_BoomerAMGCreate(&precond);
  HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
  HYPRE_BoomerAMGSetCoarsenType(precond, 6);
  HYPRE_BoomerAMGSetOldDefault(precond);
  HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
  HYPRE_BoomerAMGSetNumSweeps(precond, 1);
  HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
  HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

  hypreCheck("HYPRE_ParCSRPCGSetPrecond",
             HYPRE_ParCSRPCGSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond));
  hypreCheck("HYPRE_PCGSetup",
             HYPRE_ParCSRPCGSetup(solver, parcsr_A, parvector_b, parvector_x));
  Real a1 = platform::getRealTime();
  hypreCheck("HYPRE_PCGSolve",
             HYPRE_ParCSRPCGSolve(solver, parcsr_A, parvector_b, parvector_x));
  Real b1 = platform::getRealTime();
  info() << "Time to solve=" << (b1 - a1);

  hypreCheck("HYPRE_IJVectorGetValues",
             HYPRE_IJVectorGetValues(ij_vector_x, nb_row, rows_index_span.data(),
                                     m_dof_variable.asArray().data()));

  HYPRE_Finalize(); /* must be the last HYPRE function call */
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class HypreDoFLinearSystemFactoryService
: public ArcaneHypreDoFLinearSystemFactoryObject
{
 public:

  explicit HypreDoFLinearSystemFactoryService(const ServiceBuildInfo& sbi)
  : ArcaneHypreDoFLinearSystemFactoryObject(sbi)
  {
    info() << "Create HypreDoF";
  }

  DoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name) override
  {
    auto* x = new HypreDoFLinearSystemImpl(dof_family, solver_name);
    x->build();
    return x;
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_HYPREDOFLINEARSYSTEMFACTORY(HypreLinearSystem,
                                                    HypreDoFLinearSystemFactoryService);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
