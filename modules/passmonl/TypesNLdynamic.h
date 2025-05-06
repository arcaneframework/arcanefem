//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* TypesNLDynamic.h                                                      */
/*                                                                           */
/* Types for NLDynamic module (in axl)                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef PASSMO_TYPESNLDYNAMIC_H
#define PASSMO_TYPESNLDYNAMIC_H

#include <arcane/ItemGroup.h>

struct TypesNLDynamic
{
  enum eAnalysisType {
    PlaneStrain = 0, //!< Plane strain analysis (2D)
    PlaneStress,//!< Plane stress analysis (2D)
    Axi, //!< Axisymmetric analysis (2D)
    ThreeD, //!< 3D analysis
    NoAnalysis //!< Unknown type
  };

  // Type of elastic properties initialized on mesh cells
  enum eElastType {
    YoungNu = 0, //!< Young Modulus & Poisson Ratio
    Lame,//!< Lame coefficients (lamba & mu)
    Veloc, //!< Compression (P) & shear (S) wave velocities
    NoElastPropType //!< Unknown type
  };

  // Type of integration depending on the numerical method used (FEM, ...)
  enum eIntegType {
    FemCell = 0, //!< Looping on cells to compute on integration (Gauss) points (standard FEM method)
    FemPoint,//!< Looping on integration (Gauss) points with pre-defined neighbours (nodes and/or cells) for computations
    NoInteg //!< Unknown type
  };

  // Type of algorithm selected to solve the linear system (AX = B) in case of linear/nonlinear problems
  enum eAlgoType {
    Linear = 0, //!< Linear solving (no iterations)
    NewtonRaphson, //!< Standard Newton-Raphson: A is assembled at each timestep
    ModNewtonRaphson //!< Modified Newton-Raphson: A is assembled every 'linop_nstep' timesteps
  };

  // Type of nonlinear constitutive model
  enum eLawType {
    HOOKE = 0,
    DRUCKP, //!< Drucker-Prager
    MOHRC,//!< Mohr-Coulomb
    UNKNOWN //!< Unknown type
  };
};

#endif //PASSMO_TYPESNLDYNAMIC_H