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
  , m_dof_matrix_indexes(VariableBuildInfo(dof_family, solver_name + "DoFMatrixIndexes"))
  {
    info() << "Creating AlinaDoFLinearSystemImpl()";
  }

  ~AlinaSolver() override
  {
  }

 public:

  void build()
  {
    m_solver_parameters = makeRef(new AlinaParameters());
  }

 public:

  void applyMatrixTransformation(CsrFormatMatrixView csr_view);
  void solve(VariableDoFReal& solution_variable, VariableDoFReal& rhs_variable);

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

  AlinaParameters* solverParameters() const { return m_solver_parameters.get(); }

 private:

  IItemFamily* m_dof_family = nullptr;
  VariableDoFInt32 m_dof_matrix_indexes;
  Ref<AlinaParameters> m_solver_parameters;
  Ref<AlinaSequentialSolver> m_sequential_solver;
  UniqueArray<double> m_solver_solution;
  UniqueArray<double> m_solver_rhs;
  bool m_do_print_filling = false;

 private:

  void _fillRHSVector(VariableDoFReal& rhs_variable);
  void _fillSolutionVector(VariableDoFReal& rhs_solution);
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
applyMatrixTransformation(CsrFormatMatrixView csr_view)
{
  AlinaCSRMatrixView alina_matrix_view(csr_view.rows(), csr_view.columns(), csr_view.values());
  m_sequential_solver = makeRef(new AlinaSequentialSolver(alina_matrix_view, m_solver_parameters.get()));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaSolver::
solve(VariableDoFReal& solution_variable, VariableDoFReal& rhs_variable)
{
  _fillRHSVector(rhs_variable);
  _fillSolutionVector(solution_variable);

  info() << "Calling AlinaDoFLinearSystemImpl::solve()";

  IItemFamily* dof_family = m_dof_family;
  DoFGroup own_dofs = dof_family->allItems().own();
  const Int32 nb_dof = own_dofs.size();

  Int32 nb_iteration = 0;
  Real residual_norm = 0.0;

  AlinaConvergenceInfo convergence_info = m_sequential_solver->solve(m_solver_rhs.view(),
                                                                     m_solver_solution.view());
  info() << "ConvergenceInfo: " << convergence_info.iterations << "  r=" << convergence_info.residual;

  const bool do_verbose = (nb_dof < 200);
  Int32 index = 0;

  ENUMERATE_ (DoF, idof, dof_family->allItems().own()) {
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

class AlinaDoFLinearSystemImpl
: public DoKDoFLinearSystemImpl
{
 public:

  AlinaDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
  : DoKDoFLinearSystemImpl(dof_family, solver_name)
  {
    info() << "Creating AlinaDoFLinearSystemImpl()";
    m_alina_solver = new AlinaSolver(dof_family, solver_name);
  }

  ~AlinaDoFLinearSystemImpl() override
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
    m_alina_solver->applyMatrixTransformation(csr_view);
  }

  void solve() override
  {
    m_alina_solver->solve(solutionVariable(), rhsVariable());
  }

  void setSolverCommandLineArguments([[maybe_unused]] const CommandLineArguments& args) override
  {
  }

 public:

  AlinaParameters* solverParameters() const { return m_alina_solver->solverParameters(); }

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
    IParallelMng* pm = dof_family->parallelMng();
    if (pm->commSize() > 1)
      ARCANE_THROW(NotImplementedException, "Alina solver in parallel is not yet implemented");
    auto* x = new AlinaDoFLinearSystemImpl(dof_family, solver_name);

    x->build();

    AlinaParameters* p = x->solverParameters();

    // Setting preconditioner and solver may change other values
    // so they have to be called before others
    p->setSolverType(options()->solver());
    p->setSolverPreconditioner(options()->preconditioner());

    p->setSolverAbsoluteTolerance(options()->atol());
    p->setSolverRelativeTolerance(options()->rtol());
    p->setSolverMaxIteration(options()->maxIter());
    p->setSolverVerbosity(options()->verbosity());
    return x;
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
