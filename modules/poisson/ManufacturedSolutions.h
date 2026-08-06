// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ManufacturedSolutions.h                                     (C) 2026      */
/*                                                                           */
/* Built-in manufactured solutions for the Poisson module.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef POISSON_MANUFACTURED_SOLUTIONS_H
#define POISSON_MANUFACTURED_SOLUTIONS_H

#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/String.h>

#include <cmath>

/*---------------------------------------------------------------------------*/

namespace Arcane::FemModulePoissonManufactured
{

/*---------------------------------------------------------------------------*/

class ManufacturedSolutions
{
 public:

  static bool isEnabled(const String& name)
  {
    return !name.empty() && name != "none";
  }

  static Real exact(const String& name, const Real3& x)
  {
    return dirichlet(name, x);
  }

  static Real dirichlet(const String& name, const Real3& x)
  {
    if (name == "sine")
      return _sineExact(x);
    if (name == "poly")
      return _polyExact(x);
    if (name == "quadratic")
      return _quadraticExact(x);
    ARCANE_FATAL("Unknown Poisson manufactured solution '{0}'. Available cases are: sine, poly, quadratic", name);
  }

  static Real source(const String& name, const Real3& x)
  {
    if (name == "sine")
      return _sineSource(x);
    if (name == "poly")
      return _polySource(x);
    if (name == "quadratic")
      return _quadraticSource(x);
    ARCANE_FATAL("Unknown Poisson manufactured solution '{0}'. Available cases are: sine, poly, quadratic", name);
  }

  static Real3 gradient(const String& name, const Real3& x)
  {
    if (name == "sine")
      return _sineGradient(x);
    if (name == "poly")
      return _polyGradient(x);
    if (name == "quadratic")
      return _quadraticGradient(x);
    ARCANE_FATAL("Unknown Poisson manufactured solution '{0}'. Available cases are: sine, poly, quadratic", name);
  }

  static Real neumann(const String& name, const Real3& x, const Real3& normal)
  {
    const Real3 grad = gradient(name, x);
    return grad.x * normal.x + grad.y * normal.y + grad.z * normal.z;
  }

  static Real normalFlux(const String& name, const Real3& x, const Real3& normal)
  {
    return neumann(name, x, normal);
  }

 private:

  static constexpr Real Pi = 3.141592653589793238462643383279502884;

  static Real _sineExact(const Real3& x)
  {
    return std::sin(Pi * x.x) * std::sin(Pi * x.y);
  }

  static Real _sineSource(const Real3& x)
  {
    return 2.0 * Pi * Pi * _sineExact(x);
  }

  static Real3 _sineGradient(const Real3& x)
  {
    return {
      Pi * std::cos(Pi * x.x) * std::sin(Pi * x.y),
      Pi * std::sin(Pi * x.x) * std::cos(Pi * x.y),
      0.0
    };
  }

  static Real _polyExact(const Real3& x)
  {
    return x.x * (1.0 - x.x) * x.y * (1.0 - x.y);
  }

  static Real _polySource(const Real3& x)
  {
    return 2.0 * (x.x * (1.0 - x.x) + x.y * (1.0 - x.y));
  }

  static Real3 _polyGradient(const Real3& x)
  {
    return {
      (1.0 - 2.0 * x.x) * x.y * (1.0 - x.y),
      x.x * (1.0 - x.x) * (1.0 - 2.0 * x.y),
      0.0
    };
  }

  static Real _quadraticExact(const Real3& x)
  {
    return x.x * x.x + x.y * x.y;
  }

  static Real _quadraticSource(const Real3&)
  {
    return -4.0;
  }

  static Real3 _quadraticGradient(const Real3& x)
  {
    return { 2.0 * x.x, 2.0 * x.y, 0.0 };
  }
};

/*---------------------------------------------------------------------------*/

} // namespace Arcane::FemModulePoissonManufactured

/*---------------------------------------------------------------------------*/

using ManufacturedSolutions = Arcane::FemModulePoissonManufactured::ManufacturedSolutions;

/*---------------------------------------------------------------------------*/

#endif
