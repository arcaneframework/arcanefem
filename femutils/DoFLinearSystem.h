// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoFLinearSystem.h                                           (C) 2022-2025 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b wit DoFs.          */
/*---------------------------------------------------------------------------*/
#ifndef ARCANEFEM_FEMUTILS_DOFLINEARSYSTEM_H
#define ARCANEFEM_FEMUTILS_DOFLINEARSYSTEM_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/utils/UtilsTypes.h>

#include <arcane/core/VariableTypedef.h>
#include <arcane/core/ItemTypes.h>

#include "FemUtilsGlobal.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::FemUtils
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Helper class to handle Row elimination.
 *
 * Instances of this class are temporary and are created via a
 * call to DoFLinearSystem::rowEliminationHelper().
 *
 * There is two ways to use this class.
 *
 * 1. use the method addElimination() and specify the row and the value
 * 2. get the variables handling elimination info (getEliminationInfo())
 * and elimination value (getEliminationValue()) and set the corresponding
 * row with ELIMINATION_ROW and the value of the elimination.
 *
 * Only the method 2 is available on accelerator.
*/
class DoFLinearSystemRowEliminationHelper
{
  friend DoFLinearSystem;

 private:

  explicit DoFLinearSystemRowEliminationHelper(DoFLinearSystem* dof_ls);

 public:

  ~DoFLinearSystemRowEliminationHelper();

 public:

  /*
   * \brief Eliminate the row \a row of the linear system.
   *
   * The elimination is equivalent to the following calls:
   * - matrixSetValue(row,j,0) for j!=row
   * - matrixSetValue(row,row,1.0)
   * - RHS[rc] = value
   *
   * The row is only eliminated when solve() is called.
   * Any call to matrixAddValue(row,...)
   * or matrixSetValue(row,...) are discarded.
   *
   * \note After a row elimination the matrix may no longer be symmetric.
   */
  void addElimination(DoFLocalId row, Real value);

  VariableDoFByte& getEliminationInfo();
  VariableDoFReal& getEliminationValue();

 private:

  DoFLinearSystem* m_dof_linear_system = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Helper class to handle RowColumn elimination.
 *
 * Instances of this class are temporary and are created via a
 * call to DoFLinearSystem::rowColumnEliminationHelper().
 *
 * There is two ways to use this class.
 *
 * 1. use the method addElimination() and specify the row and the value
 * 2. get the variables handling elimination info (getEliminationInfo())
 * and elimination value (getEliminationValue()) and set the corresponding
 * row with ELIMINATION_ROW_COLUMN and the value of the elimination.
 *
 * Only the method 2 is available on accelerator.
*/
class DoFLinearSystemRowColumnEliminationHelper
{
  friend DoFLinearSystem;

 private:

  explicit DoFLinearSystemRowColumnEliminationHelper(DoFLinearSystem* dof_ls);

 public:

  ~DoFLinearSystemRowColumnEliminationHelper();

 public:

  /*
   * \brief Eliminate the row \a row of the linear system.
   *
   * The elimination is equivalent to the following calls:
   * - matrixSetValue(row,j,0) for j!=row
   * - matrixSetValue(row,row,1.0)
   * - RHS[rc] = value
   *
   * The row is only eliminated when solve() is called.
   * Any call to matrixAddValue(row,...)
   * or matrixSetValue(row,...) are discarded.
   *
   * \note After a row elimination the matrix may no longer be symmetric.
   */
  void addElimination(DoFLocalId row, Real value);

  VariableDoFByte& getEliminationInfo();
  VariableDoFReal& getEliminationValue();

 private:

  DoFLinearSystem* m_dof_linear_system = nullptr;
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
  friend DoFLinearSystemRowEliminationHelper;
  friend DoFLinearSystemRowColumnEliminationHelper;

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

  /*!
   * \brief Initialize the instance with a specific runner.
   *
   * \a runner may be null.
   */
  void initialize(ISubDomain* sd, Runner* runner, IItemFamily* dof_family, const String& solver_name);

  //! Indicate if method initialize() has been called
  [[nodiscard]] bool isInitialized() const;

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
   * \brief Helper class to eliminate rows in the linear system.
   *
   * The elimination of row \a row is equivalent to the following calls:
   * - matrixSetValue(row,j,0) for j!=row
   * - matrixSetValue(row,row,1.0)
   * - RHS[rc] = value
   *
   * The rows are only eliminated when solve() is called.
   * Any call to matrixAddValue(row,...)
   * or matrixSetValue(row,...) for eliminated rows are discarded.
   *
   * \note After a row elimination the matrix may no longer be symmetric.
   */
  DoFLinearSystemRowEliminationHelper rowEliminationHelper();

  /*
   * \brief Eliminate rows and columns of the linear system.
   *
   * The elimination for a row \a row is equivalent to the following calls:
   * - matrixSetValue(rc,j,0) for j!=rc
   * - matrixSetValue(i,rc,0) for i!=rc
   * - matrixSetValue(rcw,rc,1.0)
   * - RHS[i] = RHS[i] - A[rc,i] * value for i!=rc
   * - RHS[rc] = value
   *
   * The rows are only eliminated solve() when solve() is called.
   * Any call to matrixAddValue(row,...), matrixSetValue(row,...),
   * matrixAddValue(...,row) or matrixSetValue(...,row) for eliminated row
   * are discarded.
   */
  DoFLinearSystemRowColumnEliminationHelper rowColumnEliminationHelper();

  /*!
   * \brief Apply the underlying transformation and Solve the current linear system.
   *
   * This is equivalent to calling the following methods in order:
   *
   * - applyMatrixTransformation()
   * - applyRHSTransformation()
   * - solve()
   */
  void applyLinearSystemTransformationAndSolve();

  /*!
   * \brief Apply Matrix Transformation
   *
   * Apply operation like Row or RowColumn elimination on the matrix.
   * This has to be done one and only time before calling solve().
   * If the values of the matrix are cleared, this method has to be called again.
   */
  void applyMatrixTransformation();

  /*!
   * \brief Apply RHS Transformation
   *
   * Apply operation like Row or RowColumn elimination on the RHS.
   * This method has to be called at least one time calling solve().
   * You may call this method several time if you have manually changed the values
   * of the RHS.
   */
  void applyRHSTransformation();

  /*!
   * \brief Solve the current linear system.
   *
   * You have to call applyMatrixTransformation() and applyRHSTransformation() before
   * calling this method().
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
  [[nodiscard]] bool hasSetCSRValues() const;

 public:

  CSRFormatView& getCSRValues();
  VariableDoFReal& getForcedValue();
  VariableDoFBool& getForcedInfo();

  void setConstantMatrixSparsity(bool val);
  void setConstantMatrixValues(bool val);
  bool isMatrixSparsityConstant();
  bool isMatrixValuesConstant();

  IDoFLinearSystemFactory* linearSystemFactory() const
  {
    return m_linear_system_factory;
  }

 private:

  IDoFLinearSystemImpl* m_p = nullptr;
  IItemFamily* m_item_family = nullptr;
  IDoFLinearSystemFactory* m_linear_system_factory = nullptr;
  IDoFLinearSystemFactory* m_default_linear_system_factory = nullptr;

 private:

  void _checkInit() const;

  // Used by DoFLinearSystemRowEliminationHelper
  void _eliminateRow(DoFLocalId row, Real value);
  // Used by DoFLinearSystemRowColumnEliminationHelper
  void _eliminateRowColumn(DoFLocalId row, Real value);

  // Used by DoFLinearSystemRowEliminationHelper and DoFLinearSystemRowColumnEliminationHelper
  VariableDoFByte& _getEliminationInfo();
  // Used by DoFLinearSystemRowEliminationHelper and DoFLinearSystemRowColumnEliminationHelper
  VariableDoFReal& _getEliminationValue();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemUtils

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
