// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* HypreDoFLinearSystem.cc                                     (C) 2022-2025 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "DoFLinearSystem.h"

#include <arcane/accelerator/RunCommandLoop.h>
#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/PlatformUtils.h>
#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/MemoryUtils.h>
#include <arcane/utils/MemoryView.h>
#include <arcane/utils/ITraceMng.h>
#include <arcane/utils/NumArray.h>

#include <arcane/core/ServiceFactory.h>
#include <arcane/core/VariableTypes.h>
#include <arcane/core/BasicService.h>
#include <arcane/core/VariableUtils.h>
#include <arcane/core/IParallelMng.h>
#include <arcane/core/IItemFamily.h>
#include <arcane/core/ItemPrinter.h>
#include <arcane/core/Timer.h>

#include <arcane/accelerator/VariableViews.h>
#include <arcane/accelerator/core/Runner.h>
#include <arcane/accelerator/core/Memory.h>

#include "IDoFLinearSystemFactory.h"
#include "internal/CsrDoFLinearSystemImpl.h"

#include <petsc.h>
#include <petscpctypes.h>
#include <krylov.h>
#include "PETScDoFLinearSystemFactory_axl.h"

namespace Arcane::FemUtils
{

class PETScDoFLinearSystemImpl
: public CsrDoFLinearSystemImpl
{
 public:

  HypreDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
  : CsrDoFLinearSystemImpl(dof_family, solver_name)
  , m_dof_matrix_numbering(VariableBuildInfo(dof_family, solver_name + "MatrixNumbering"))
  {
    info() << "[Hypre-Info] Creating HypreDoFLinearSystemImpl()";
  }

  ~HypreDoFLinearSystemImpl() override
  {
    info() << "[Hypre-Info] Calling HYPRE_Finalize";
#if HYPRE_RELEASE_NUMBER >= 21500
    HYPRE_Finalize(); // must be the last HYPRE function call //
#endif
  }

 public:

  void build()
  {
#if HYPRE_RELEASE_NUMBER >= 22700
    HYPRE_Init(); // must be the first HYPRE function call //
#endif
  }

 public:

  void solve() override;

  void setSolverCommandLineArguments(const CommandLineArguments& args) override
  {
  }

  void setMaxIter(Int32 v) { m_max_iter = v; }
  void setVerbosityLevel(Int32 v) { m_verbosity = v; }
  void setAmgCoarsener(Int32 v) { m_amg_coarsener = v; }
  void setAmgInterpType(Int32 v) { m_amg_interp_type = v; }
  void setAmgSmoother(Int32 v) { m_amg_smoother = v; }
  void setKrylovDim(Int32 v) { m_krylov_dim = v; }

  void setRelTolerance(Real v) { m_rtol = v; }
  void setAbsTolerance(Real v) { m_atol = v; }
  void setAmgThreshold(Real v) { m_amg_threshold = v; }

  void setSolver(solver v) { m_solver = v; }
  void setPreconditioner(preconditioner v) { m_preconditioner = v; }

 private:

  VariableDoFInt32 m_dof_matrix_numbering;

  NumArray<Int32, MDDim1> m_parallel_columns_index;
  NumArray<Int32, MDDim1> m_parallel_rows_index;
  //! Work array to store values of solution vector in parallel
  NumArray<Real, MDDim1> m_result_work_values;

  Int32 m_first_own_row = -1;
  Int32 m_nb_own_row = -1;
  Int32 m_max_iter = 1000;
  Int32 m_verbosity = 2;
  Int32 m_amg_coarsener = 8;
  Int32 m_amg_interp_type = 6;
  Int32 m_amg_smoother = 6;
  Int32 m_krylov_dim = 2;

  Real m_amg_threshold = 0.25;
  Real m_rtol = 1.0e-7;
  Real m_atol = 0.;

  solver m_solver = solver::CG;
  preconditioner m_preconditioner = preconditioner::AMG;

 private:

  void _computeMatrixNumeration();
};

#include "PETScDoFLinearSystemFactory_axl.h"

ARCANE_REGISTER_SERVICE_PETSCDOFLINEARSYSTEMFACTORY(PETScLinearSystem,
                                                    PetscDoFLinearSystemFactoryService);

}

