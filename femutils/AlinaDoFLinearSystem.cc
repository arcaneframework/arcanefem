// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* AlinaDoFLinearSystem.cc                                     (C) 2000-2026 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/FatalErrorException.h>
#include <arcane/utils/PlatformUtils.h>

#include <arcane/core/VariableTypes.h>
#include <arcane/core/IItemFamily.h>
#include <arcane/core/ItemGroup.h>
#include <arcane/core/IParallelMng.h>

#include <arccore/alina/AlinaLib.h>

#include "FemUtils.h"
#include "internal/DoKDoFLinearSystemImpl.h"
#include "internal/CsrDoFLinearSystemImpl.h"
#include "IDoFLinearSystemFactory.h"
#include "CsrFormatMatrixView.h"

#include "AlinaDoFLinearSystemFactory_axl.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
using namespace Arcane::AlinaLib;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AlinaSolver
: public TraceAccessor
{
  using RowColumn = DoKMatrix::RowColumn;

 public:

  AlinaSolver(IItemFamily* dof_family, const String& solver_name)
  : TraceAccessor(dof_family->traceMng())
  , m_dof_family(dof_family)
  , m_dof_matrix_numbering(VariableBuildInfo(dof_family, solver_name + "MatrixNumbering"))
  {
    info() << "Creating AlinaDoFLinearSystemImpl()";
  }

  ~AlinaSolver() override
  {
  }

 public:

  void build()
  {
    MessagePassing::IMessagePassingMng* mpm = m_dof_family->parallelMng()->messagePassingMng();
    m_solver_parameters = makeRef(new AlinaSolverParameters(mpm,traceMng()));
  }

 public:

  //void applyMatrixTransformation(CsrFormatMatrixView csr_view);
  void solve(CsrFormatMatrixView csr_view, VariableDoFReal& solution_variable, VariableDoFReal& rhs_variable);

 public:

  void setRelTolerance(Real v)
  {
    m_solver_parameters->setSolverRelativeTolerance(v);
  }
  void setAbsTolerance(Real v)
  {
    m_solver_parameters->setSolverAbsoluteTolerance(v);
  }
  void setMaxIteration(Int32 v)
  {
    m_solver_parameters->setSolverMaxIteration(v);
  }

  AlinaSolverParameters* solverParameters() const { return m_solver_parameters.get(); }

 private:

  IItemFamily* m_dof_family = nullptr;
  IParallelMng* m_parallel_mng = m_dof_family->parallelMng();
  VariableDoFInt32 m_dof_matrix_numbering;
  Ref<AlinaSolverParameters> m_solver_parameters;
  Ref<AlinaLib::AlinaSolver> m_solver;
  UniqueArray<double> m_solver_solution;
  UniqueArray<double> m_solver_rhs;
  bool m_do_print_filling = false;

 private:

  void _fillRHSVector(VariableDoFReal& rhs_variable);
  void _fillSolutionVector(VariableDoFReal& rhs_solution);
  void _computeMatrixNumbering();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaSolver::
_fillRHSVector(VariableDoFReal& rhs_variable)
{
  // For the LinearSystem class we need an array
  // with only the values for the ownNodes().
  // The values of 'rhs_values' should not be updated after
  // this call.
  VariableDoFReal& rhs_values(rhs_variable);
  Int32 nb_item = m_dof_family->allItems().own().size();
  m_solver_rhs.resize(nb_item);
  Int32 index = 0;
  ENUMERATE_ (DoF, idof, m_dof_family->allItems().own()) {
    Real v = rhs_values[idof];
    if (m_do_print_filling)
      info() << "SET VECTOR VALUE (" << std::setw(4) << idof.itemLocalId() << ") = " << v;
    m_solver_rhs[index] = rhs_values[idof];
    ++index;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaSolver::
_fillSolutionVector(VariableDoFReal& solution_variable)
{
  VariableDoFReal& solution_values(solution_variable);
  Int32 nb_item = m_dof_family->allItems().own().size();
  m_solver_solution.resize(nb_item);
  Int32 index = 0;
  ENUMERATE_ (DoF, idof, m_dof_family->allItems().own()) {
    Real v = solution_values[idof];
    if (m_do_print_filling)
      info() << "SET SOLUTION VECTOR VALUE (" << std::setw(4) << idof.itemLocalId() << ") = " << v;
    m_solver_solution[index] = solution_values[idof];
    ++index;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaSolver::
_computeMatrixNumbering()
{
  // TODO: This code is similar with the code in Hypre or PETSc.
  // Create a class to handle that.
  IItemFamily* dof_family = m_dof_family;
  IParallelMng* pm = dof_family->parallelMng();
  const bool is_parallel = pm->isParallel();
  const Int32 nb_rank = pm->commSize();
  const Int32 my_rank = pm->commRank();

  DoFGroup all_dofs = dof_family->allItems();
  DoFGroup own_dofs = all_dofs.own();
  const Int32 nb_own_row = own_dofs.size();

  Int32 own_first_index = 0;

  if (is_parallel) {
    // TODO: utiliser un Scan lorsque ce sera disponible dans Arcane
    UniqueArray<Int32> parallel_rows_index(nb_rank, 0);
    pm->allGather(ConstArrayView<Int32>(1, &nb_own_row), parallel_rows_index);
    info() << "ALL_NB_ROW = " << parallel_rows_index;
    for (Int32 i = 0; i < my_rank; ++i)
      own_first_index += parallel_rows_index[i];
  }

  info() << "OwnFirstIndex=" << own_first_index << " NbOwnRow=" << nb_own_row;

  //m_first_own_row = own_first_index;
  //m_nb_own_row = nb_own_row;

  // TODO: Faire avec API accelerateur
  ENUMERATE_DOF (idof, own_dofs) {
    DoF dof = *idof;
    m_dof_matrix_numbering[idof] = own_first_index + idof.index();
    //info() << "Numbering dof_uid=" << dof.uniqueId() << " M=" << m_dof_matrix_numbering[idof];
  }
  info() << " nb_own_row=" << nb_own_row << " nb_item=" << dof_family->nbItem();
  m_dof_matrix_numbering.synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaSolver::
solve(CsrFormatMatrixView csr_view, VariableDoFReal& solution_variable, VariableDoFReal& rhs_variable)
{
  _fillRHSVector(rhs_variable);
  _fillSolutionVector(solution_variable);

  const bool is_parallel = m_parallel_mng->isParallel();
  if (is_parallel)
    _computeMatrixNumbering();

  info() << "Calling AlinaDoFLinearSystemImpl::solve()";

  IItemFamily* dof_family = m_dof_family;
  DoFGroup all_dofs = dof_family->allItems();
  DoFGroup own_dofs = all_dofs.own();
  const Int32 nb_own_dof = own_dofs.size();

  Int32 nb_iteration = 0;
  Real residual_norm = 0.0;

  //! Indexes of the columns (in global numbering)
  NumArray<Int32, MDDim1> m_parallel_columns_index;
  // Indexes of own rows (exluding rows which are not from own items)
  NumArray<Int32, MDDim1> m_parallel_rows_index;
  // Indexes of own rows (exluding rows which are not from own items)
  NumArray<Real, MDDim1> m_parallel_values;

  // csr_view.columns() use matrix coordinates local to sub-domain
  // We need to translate them to global matrix coordinates
  // The original matrix (csr_view) contains values for ghost items
  // so we need to remove them.
  // NOTE: This rebuilding is not needed if the matrix structure and/or
  // values are unchanged.
  SmallSpan<const Int32> rows_index_span = csr_view.rows();
  SmallSpan<const Int32> columns_index_span = csr_view.columns();
  SmallSpan<const Real> values_index_span = csr_view.values();
  if (is_parallel) {
    Int32 nb_orig_row = csr_view.nbRow();
    Int32 nb_orig_value = rows_index_span[nb_orig_row];

    m_parallel_columns_index.resize(nb_orig_value);
    m_parallel_rows_index.resize(csr_view.rows().size());
    m_parallel_values.resize(nb_orig_value);

    Int32 index = 0;
    Int32 parallel_nb_row = 0;
    Int32 parallel_row_index = 0;
    Int32 column_index = 0;
    m_parallel_rows_index[0] = 0;
    ++index;
    ENUMERATE_ (DoF, idof, all_dofs) {
      DoF dof = *idof;
      if (!dof.isOwn())
        continue;
      for (CsrRowColumnIndex rc : csr_view.rowRange(idof.index())) {
        DoFLocalId local_column(csr_view.column(rc));
        m_parallel_columns_index[column_index] = m_dof_matrix_numbering[local_column];
        m_parallel_values[column_index] = csr_view.value(rc);
        ++column_index;
      }
      m_parallel_rows_index[index] = column_index;
      ++index;
    }

    columns_index_span = m_parallel_columns_index.to1DSmallSpan().subSpan(0, column_index);
    values_index_span = m_parallel_values.to1DSmallSpan().subSpan(0, column_index);
    rows_index_span = m_parallel_rows_index.to1DSmallSpan().subSpan(0, index);
  }

  AlinaCSRMatrixView alina_matrix_view(rows_index_span, columns_index_span, values_index_span);

  {
    m_solver_parameters->setString("precond.coarsening.type", "smoothed_aggregation");
    m_solver_parameters->setString("precond.relax.type", "ilu0");
    double t0 = platform::getRealTime();

    m_solver = makeRef(new AlinaLib::AlinaSolver(*m_solver_parameters.get(), alina_matrix_view));

    double t1 = platform::getRealTime();
    info() << "[Alina-Timer] Time to setup = " << (t1 - t0);

    AlinaConvergenceInfo convergence_info = m_solver->solve(m_solver_rhs.view(),
                                                            m_solver_solution.view());
    double t2 = platform::getRealTime();
    info() << "ConvergenceInfo: " << convergence_info.iterations << "  r=" << convergence_info.residual;
    info() << "[Alina-Timer] Time to solve = " << (t2 - t1);
  }

  const bool do_verbose = (nb_own_dof < 200);
  Int32 index = 0;

  ENUMERATE_ (DoF, idof, own_dofs) {
    DoF dof = *idof;

    solution_variable[dof] = m_solver_solution[index];
    if (do_verbose)
      info() << "Node uid=" << dof.uniqueId() << " T=" << solution_variable[dof];
    ++index;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AlinaDoKDoFLinearSystemImpl
: public DoKDoFLinearSystemImpl
{
 public:

  AlinaDoKDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
  : DoKDoFLinearSystemImpl(dof_family, solver_name)
  {
    info() << "Creating AlinaDoFLinearSystemImpl()";
    m_alina_solver = new AlinaSolver(dof_family, solver_name);
  }

  ~AlinaDoKDoFLinearSystemImpl() override
  {
    delete m_alina_solver;
  }

 public:

  void build()
  {
    m_alina_solver->build();
    DoKDoFLinearSystemImpl::clearValues();
  }

 public:

  void applyMatrixTransformation() override
  {
    fillRowColumnEliminationInfos();
    convertToCSRMatrix();
    CsrFormatMatrixView csr_view = getCsrFormatMatrixView();
  }

  void solve() override
  {
    m_alina_solver->solve(getCsrFormatMatrixView(), solutionVariable(), rhsVariable());
  }

  void setSolverCommandLineArguments([[maybe_unused]] const CommandLineArguments& args) override
  {
  }

 public:

  AlinaSolver* underlyingAlinaSolver() const { return m_alina_solver; }

 private:

  AlinaSolver* m_alina_solver = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AlinaCsrDoFLinearSystemImpl
: public CsrDoFLinearSystemImpl
{
 public:

  AlinaCsrDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
  : CsrDoFLinearSystemImpl(dof_family, solver_name)
  {
    info() << "Creating AlinaDoFLinearSystemImpl()";
    m_alina_solver = new AlinaSolver(dof_family, solver_name);
  }

  ~AlinaCsrDoFLinearSystemImpl() override
  {
    delete m_alina_solver;
  }

 public:

  void build()
  {
    m_alina_solver->build();
    CsrDoFLinearSystemImpl::clearValues();
  }

 public:

  void applyMatrixTransformation() override
  {
  }

  void solve() override
  {
    m_alina_solver->solve(getCSRValues(), solutionVariable(), rhsVariable());
  }

  void setSolverCommandLineArguments([[maybe_unused]] const CommandLineArguments& args) override
  {
  }

 public:

  AlinaSolver* underlyingAlinaSolver() const { return m_alina_solver; }

 private:

  AlinaSolver* m_alina_solver = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AlinaDoFLinearSystemFactoryService
: public ArcaneAlinaDoFLinearSystemFactoryObject
{
 public:

  explicit AlinaDoFLinearSystemFactoryService(const ServiceBuildInfo& sbi)
  : ArcaneAlinaDoFLinearSystemFactoryObject(sbi)
  {
  }

  IDoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name) override
  {
    return _createInstance(dof_family, solver_name, eLinearSystemMatrixFormat::Csr);
  }

  IDoFLinearSystemImpl*
  createInstance(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name,
                 eLinearSystemMatrixFormat matrix_format) override
  {
    return _createInstance(dof_family, solver_name, matrix_format);
  }

  IDoFLinearSystemImpl*
  _createInstance(IItemFamily* dof_family, const String& solver_name,
                  eLinearSystemMatrixFormat matrix_format)
  {
    AlinaSolver* x = nullptr;
    IDoFLinearSystemImpl* linear_system = nullptr;
    if (matrix_format == eLinearSystemMatrixFormat::DoK) {
      info() << "Using DoK format Alina linear system";
      auto* v = new AlinaDoKDoFLinearSystemImpl(dof_family, solver_name);
      linear_system = v;
      x = v->underlyingAlinaSolver();
    }
    else if (matrix_format == eLinearSystemMatrixFormat::Csr) {
      info() << "Using Csr format Alina linear system";
      auto* v = new AlinaCsrDoFLinearSystemImpl(dof_family, solver_name);
      linear_system = v;
      x = v->underlyingAlinaSolver();
    }
    else
      ARCANE_FATAL("Unsupported matrix_format '{0}'", static_cast<int>(matrix_format));
    _initializeAlinaSolver(x);
    return linear_system;
  }

  void _initializeAlinaSolver(AlinaSolver* x)
  {
    x->build();
    AlinaSolverParameters* p = x->solverParameters();
    // Setting preconditioner and solver may change other values
    // so they have to be called before others
    p->setSolverType(options()->solver());
    p->setSolverPreconditioner(options()->preconditioner());
    p->setSolverAbsoluteTolerance(options()->atol());
    p->setSolverRelativeTolerance(options()->rtol());
    p->setSolverMaxIteration(options()->maxIter());
    p->setSolverVerbosity(options()->verbosity());
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ALINADOFLINEARSYSTEMFACTORY(AlinaLinearSystem,
                                                    AlinaDoFLinearSystemFactoryService);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
