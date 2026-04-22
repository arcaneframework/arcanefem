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
//
//----------------------------------------------
// exact solution:
// u(x,y) = ((2^(m+1) - 1) x + 1 )^(1/(m+1)) - 1
//
// this is given by function
//  Real exactSolution(m, (x, y, z))
//
// the nonlinearity term or
// nonconstant coefficient lamda:
//
// lambda(m, u) = ( 1 + u ) ^ m
//
// this is given by function
// Real labmda_funk(m, uk_node)
//----------------------------------------------------------------

//----------------------------------------------------------------

using System;
using Arcane;
using Real = System.Double;

namespace FemModuleFourierNL
{
  public class CaseFunctions
  {
//     public
//     Real manufacturedDirichlet(Real lambda, Real3 position)
//     {
//       return lambda * (System.Math.Sin(position.x)  +  System.Math.Cos(position.y));
//     }
    public

    Real manufacturedSource(Real alpha, Real3 position)
    {
      return alpha * (System.Math.Sin(position.x)  +  System.Math.Cos(position.y));
    }
//     public

//     Real exactSolution(Real m, Real3 position)
//     Real manufacturedDirichlet(Real lambda, Real3 position)
//     {
//       return System.Math.Pow( (System.Math.Pow(2, m + 1) - 1) * position.x + 1, 1 / (m+1) ) - 1;
//     }
//     public
//
//     Real lambdaFunc(Real m , Real uk_node)
//     {
//       return System.Math.Pow(1 + uk_node, m);
    }
  }
}
