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

class AlinaDoFLinearSystemImpl
: public DoKDoFLinearSystemImpl
{
  using RowColumn = DoKMatrix::RowColumn;

 public:

  AlinaDoFLinearSystemImpl(IItemFamily* dof_family, const String& solver_name)
  : DoKDoFLinearSystemImpl(dof_family, solver_name)
  , m_dof_matrix_indexes(VariableBuildInfo(dof_family, solver_name + "DoFMatrixIndexes"))
  {
    info() << "Creating AlinaDoFLinearSystemImpl()";
  }

  ~AlinaDoFLinearSystemImpl() override
  {
  }

 public:

  void build()
  {
    _computeMatrixInfo();
    m_solver_parameters = makeRef(new AlinaParameters());
    DoKDoFLinearSystemImpl::clearValues();
  }

 private:

  void _computeMatrixInfo();

 public:

  void applyMatrixTransformation() override;
  void solve() override;

  void setSolverCommandLineArguments(const CommandLineArguments& args) override
  {
  }

  void clearValues() override
  {
    info() << "[Alina] Clear values of current solver";
    DoKDoFLinearSystemImpl::clearValues();
    _computeMatrixInfo();
  }

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

  VariableDoFInt32 m_dof_matrix_indexes;
  Ref<AlinaParameters> m_solver_parameters;
  Ref<AlinaSequentialSolver> m_sequential_solver;
  UniqueArray<double> m_solver_solution;
  UniqueArray<double> m_solver_rhs;
  bool m_do_print_filling = false;

 private:

  void _applyMatrixTransformationAndFillAlinaMatrix();
  void _fillRHSVector();
  void _fillSolutionVector();
  void _createRHSAndSolutionVector();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaDoFLinearSystemImpl::
_applyMatrixTransformationAndFillAlinaMatrix()
{
  fillRowColumnEliminationInfos();
  // We provide two ways to fill the Alina matrix.
  // The first one (currently the default) fill the matrix using the DoK Matrix.
  // The second one converts the DoKMatrix to a CSR Matrix then fill the Alina Matrix.
  // The second one will be used if we want to directly call other linear solver
  // like PETSc or Hypre when using a DoK Matrix.
  convertToCSRMatrix();
  CsrFormatMatrixView csr_view = getCsrFormatMatrixView();

  AlinaCSRMatrixView alina_matrix_view(csr_view.rows(), csr_view.columns(), csr_view.values());
  m_sequential_solver = makeRef(new AlinaSequentialSolver(alina_matrix_view, m_solver_parameters.get()));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaDoFLinearSystemImpl::
_createRHSAndSolutionVector()
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaDoFLinearSystemImpl::
_fillRHSVector()
{
  // For the LinearSystem class we need an array
  // with only the values for the ownNodes().
  // The values of 'rhs_values' should not be updated after
  // this call.
  VariableDoFReal& rhs_values(rhsVariable());
  IItemFamily* dof_family = dofFamily();
  Int32 nb_item = dof_family->allItems().own().size();
  m_solver_rhs.resize(nb_item);
  Int32 index = 0;
  ENUMERATE_ (DoF, idof, dof_family->allItems().own()) {
    Real v = rhs_values[idof];
    if (m_do_print_filling)
      info() << "SET VECTOR VALUE (" << std::setw(4) << idof.itemLocalId() << ") = " << v;
    m_solver_rhs[index] = rhs_values[idof];
    ++index;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaDoFLinearSystemImpl::
_fillSolutionVector()
{
  VariableDoFReal& solution_values(solutionVariable());
  IItemFamily* dof_family = dofFamily();
  Int32 nb_item = dof_family->allItems().own().size();
  m_solver_solution.resize(nb_item);
  Int32 index = 0;
  ENUMERATE_ (DoF, idof, dof_family->allItems().own()) {
    Real v = solution_values[idof];
    if (m_do_print_filling)
      info() << "SET SOLUTION VECTOR VALUE (" << std::setw(4) << idof.itemLocalId() << ") = " << v;
    m_solver_solution[index] = solution_values[idof];
    ++index;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaDoFLinearSystemImpl::
_computeMatrixInfo()
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaDoFLinearSystemImpl::
applyMatrixTransformation()
{
  // Matrix transformation
  _applyMatrixTransformationAndFillAlinaMatrix();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AlinaDoFLinearSystemImpl::
solve()
{
  _createRHSAndSolutionVector();
  _fillRHSVector();
  _fillSolutionVector();

  info() << "Calling AlinaDoFLinearSystemImpl::solve()";

  IItemFamily* dof_family = dofFamily();
  DoFGroup own_dofs = dof_family->allItems().own();
  const Int32 nb_dof = own_dofs.size();

  Int32 nb_iteration = 0;
  Real residual_norm = 0.0;

  AlinaConvergenceInfo convergence_info = m_sequential_solver->solve(m_solver_rhs.view(),
                                                                     m_solver_solution.view());
  info() << "ConvergenceInfo: " << convergence_info.iterations << "  r=" << convergence_info.residual;

  const bool do_verbose = (nb_dof < 200);
  Int32 index = 0;

  VariableDoFReal& solution_variable(this->solutionVariable());
  ENUMERATE_ (DoF, idof, dofFamily()->allItems().own()) {
    DoF dof = *idof;

    solution_variable[dof] = m_solver_solution[index];
    if (do_verbose)
      info() << "Node uid=" << dof.uniqueId() << " T=" << solution_variable[dof];
    ++index;
  }
}

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
