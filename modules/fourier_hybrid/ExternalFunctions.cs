//----------------------------------------------------------------
// manufactured solution
// u(x,y) = lambda *  ( sin(x) + cos(y) )
//
// for this case the manufactured dirichlet will be
// the same as manufactured solution
//
// u_D(x,y) = lambda *  ( sin(x) + cos(y) )
//
// this is given by function
//    Real manufacturedDirichlet( )
//
// manufactured source term f(x,y) will be d^2(u):
//
// f(x,y) = - lambda * ( sin(x) + cos(y) ) * A/3
//
// this is given by function
//    Real manufacturedSource( )
//----------------------------------------------------------------

using System;
using Arcane;
using Real = System.Double;

namespace FemModuleFourier
{
  public class CaseFunctions
  {
    public
    Real manufacturedDirichlet(Real lambda, Real3 position)
    {
      return lambda * (System.Math.Sin(position.x)  +  System.Math.Cos(position.y));
    }
    public

    Real manufacturedSource(Real alpha, Real3 position)
    {
      return alpha * (System.Math.Sin(position.x)  +  System.Math.Cos(position.y));
    }
  }
}
