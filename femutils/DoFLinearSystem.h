// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoFLinearSystem.h                                           (C) 2022-2024 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b wit DoFs.          */
/*---------------------------------------------------------------------------*/
#ifndef FEMTEST_DOFLINEARSYSTEM_H
#define FEMTEST_DOFLINEARSYSTEM_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/ArrayView.h>
#include <arcane/ItemTypes.h>
#include <arcane/VariableTypedef.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{
class IDoFLinearSystemFactory;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Vue au format CSR pour le solveur linéaire.
 */
class CSRFormatView
{
 public:

  CSRFormatView() = default;
  CSRFormatView(Span<const Int32> rows,
                Span<const Int32> matrix_rows_nb_column,
                Span<const Int32> columns,
                Span<const Real> values)
  : m_matrix_rows(rows)
  , m_matrix_rows_nb_column(matrix_rows_nb_column)
  , m_matrix_columns(columns)
  , m_values(values)
  {
  }

 public:

  Span<const Int32> rows() const { return m_matrix_rows; }
  Span<const Int32> rowsNbColumn() const { return m_matrix_rows_nb_column; }
  Span<const Int32> columns() const { return m_matrix_columns; }
  Span<const Real> values() const { return m_values; }

 private:

  Span<const Int32> m_matrix_rows;
  Span<const Int32> m_matrix_rows_nb_column;
  Span<const Int32> m_matrix_columns;
  Span<const Real> m_values;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \internal.
 * \brief Implementation for NodeLinearSystem.
 */
class DoFLinearSystemImpl
{
 public:

  virtual ~DoFLinearSystemImpl() = default;

 public:

  virtual void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value) = 0;
  virtual void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value) = 0;
  virtual void eliminateRow(DoFLocalId row, Real value) = 0;
  virtual void eliminateRowColumn(DoFLocalId row, Real value) = 0;
  virtual void setEliminationArrays(VariableDoFByte &dof_elimination_info, VariableDoFReal &dof_elimination_value) = 0;
  virtual void solve() = 0;
  virtual VariableDoFReal& solutionVariable() = 0;
  virtual VariableDoFReal& rhsVariable() = 0;
  virtual void setSolverCommandLineArguments(const CommandLineArguments& args) = 0;
  virtual void clearValues() = 0;
  virtual void setCSRValues(const CSRFormatView& csr_view) = 0;
  virtual bool hasSetCSRValues() const = 0;
  virtual void setRunner(Runner* r) =0;
  virtual Runner* runner() const =0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Linear system Matrix A + Vector x + Vector b for Ax=b.
 *
 * Before using an instance of this class you need to call method
 * initialize(). If you want to reuse the same instance for several solving
 * you need to call reset() to destroy the underlying linear system and then
 * you need to call initialize() again.
 *
 * The solve() method solves the current linear system. After this variable
 * returned by the method solutionVariable() will be filled with the values
 * of the solution vector.
 */
class DoFLinearSystem
{
 public:

  DoFLinearSystem();
  ~DoFLinearSystem();
  DoFLinearSystem(const DoFLinearSystem&) = delete;
  DoFLinearSystem(DoFLinearSystem&&) = delete;
  DoFLinearSystem& operator=(DoFLinearSystem&&) = delete;
  DoFLinearSystem& operator=(const DoFLinearSystem&) = delete;

 public:

  /*
   * \brief Initialize the instance.
   */
  void initialize(ISubDomain* sd, IItemFamily* dof_family, const String& solver_name);

  /*
   * \brief Initialize the instance.
   *
   * \a runner may be null.
   */
  void initialize(ISubDomain* sd, Runner* runner, IItemFamily* dof_family, const String& solver_name);

  //! Indicate if method initialize() has been called
  bool isInitialized() const;

  //! Add the value \a value to the (row,column) element of the matrix
  void matrixAddValue(DoFLocalId row, DoFLocalId column, Real value);

  /*!
   * \brief Set the value \a value to the (row,column) element of the matrix.
   *
   * The value is only changed when solve() is called. Any call to matrixAddValue(row,column)
   * done before or after the call to matrixSetValue(row,column) are discarded.
   */
  void matrixSetValue(DoFLocalId row, DoFLocalId column, Real value);

  /*
   * \brief Eliminate the row \a row of the linear system.
   *
   * The elimination is equivalent to the following calls:
   * - matrixSetValue(row,j,0) for j!=row
   * - matrixSetValue(row,row,1.0)
   * - RHS[rc] = value
   *
   * The row is only eliminated when solve() is called. Any call to matrixAddValue(row,...)
   * or matrixSetValue(row,...) are discarded.
   *
   * \note After a row elimination the matrix may no longer be symmetric.
   */
  void eliminateRow(DoFLocalId row, Real value);

  /*
   * \brief Eliminate the row \a rc and column \a rc of the linear system.
   *
   * The elimination is equivalent to the following calls:
   * - matrixSetValue(rc,j,0) for j!=rc
   * - matrixSetValue(i,rc,0) for i!=rc
   * - matrixSetValue(rcw,rc,1.0)
   * - RHS[i] = RHS[i] - A[rc,i] * value for i!=rc
   *
   * The row is only eliminated solve() is called.
   * Any call to matrixAddValue(row,...), matrixSetValue(row,...),
   * matrixAddValue(...,row) or matrixSetValue(...,row) are discarded.
   */
  void eliminateRowColumn(DoFLocalId row, Real value);

  void setEliminationArrays(VariableDoFByte &dof_elimination_info, VariableDoFReal &dof_elimination_value);

  /*!
   * \brief Solve the current linear system.
   */
  void solve();

  /*!
   * \brief Reset the current instance.
   *
   * You have to call initialize() again to re-use the same instance.
   */
  void reset();

  /*!
   * \brief Clear values.
   *
   * After this call, the matrix and the RHS vector will have their
   * values cleared.
   */
  void clearValues();

  /*!
   * \brief Variable containing the solution vector.
   *
   * The values of this variable are only relevant after a call to solve().
   */
  VariableDoFReal& solutionVariable();

  /*!
   * \brief Variable containing the right hand side vector.
   *
   * The values of this variable will be used during the solve() call to
   * fill the right hand side vector.
   */
  VariableDoFReal& rhsVariable();

  //! Set the factory used to create the underlying linear system solver
  void setLinearSystemFactory(IDoFLinearSystemFactory* factory)
  {
    m_linear_system_factory = factory;
  }

  /*
   * \brief Set the arguments for the solver initialisation.
   *
   * This method has to be called before the first call to solve().
   *
   * Currently it only works when the implementation is Aleph and when
   * the underlying solver is PETSc. In this case the arguments in \a args
   * are passed to PetscInitialize().
   *
   * \note The call to PetscInitialize() is only done one time during the
   * simulation so only the first solver using PETSc will proceed
   * to the initialisation.
   */
  void setSolverCommandLineArguments(const CommandLineArguments& args);

  /*!
   * \brief Positionne directement les valeurs de la matrice.
   *
   * Positionne les valeurs de la matrice en considérant le format comme
   * étant au format CSR. Les vues doivent rester valides jusqu'à la
   * résolution du système linéaire (appel à solve()).
   */
  void setCSRValues(const CSRFormatView& csr_view);

  //! Indique si l'implémentation supporte d'utiliser setCSRValue()
  bool hasSetCSRValues() const;

 public:

  IDoFLinearSystemFactory* linearSystemFactory() const
  {
    return m_linear_system_factory;
  }

 private:

  DoFLinearSystemImpl* m_p = nullptr;
  IItemFamily* m_item_family = nullptr;
  IDoFLinearSystemFactory* m_linear_system_factory = nullptr;
  IDoFLinearSystemFactory* m_default_linear_system_factory = nullptr;

 private:

  void _checkInit() const;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
