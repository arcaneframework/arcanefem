// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DoFLinearSystem.h                                           (C) 2022-2022 */
/*                                                                           */
/* Linear system: Matrix A + Vector x + Vector b for Ax=b wit DoFs.          */
/*---------------------------------------------------------------------------*/
#ifndef FEMTEST_DOFLINEARSYSTEM_H
#define FEMTEST_DOFLINEARSYSTEM_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

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
 * \internal.
 * \brief Implementation for NodeLinearSystem.
 */
class DoFLinearSystemImpl
{
 public:

  virtual ~DoFLinearSystemImpl() = default;

 public:

  virtual void matrixAddValue(Arcane::DoFLocalId row, Arcane::DoFLocalId column, Arcane::Real value) = 0;
  virtual void setRHSValues(Arcane::Span<const Arcane::Real> values) = 0;
  virtual void solve() = 0;
  virtual Arcane::VariableDoFReal& solutionVariable() = 0;
  virtual Arcane::VariableDoFReal& rhsVariable() = 0;
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

  /*!
   * \brief Initialize the instance.
   *
   * The variable dof_variable will be filled with the solution value after
   * the call to the method solve().
   */
  void initialize(Arcane::ISubDomain* sd, Arcane::IItemFamily* dof_family, const Arcane::String& solver_name);

  //! Add the value \a value to the (row,column) element of the matrix
  void matrixAddValue(Arcane::DoFLocalId row, Arcane::DoFLocalId column, Arcane::Real value);

  /*!
   * \brief Set the values for vector B.
   *
   * There is one value in values for each own nodes of the current sub-domain.
   */
  void setRHSValues(Span<const Real> values);

  /*!
   * \brief Solve the current linear system.
   */
  void solve();

  /*!
   * \brief Reset the current instance.
   *
   * You have to call initialize() again to re-use the same instance
   */
  void reset();

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
    factory = m_linear_system_factory;
  }

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

  void _checkInit();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
