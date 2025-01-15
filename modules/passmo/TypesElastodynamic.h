//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* TypesElastodynamic.h                                                      */
/*                                                                           */
/* Types for Elastodynamic module (in axl)                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef PASSMO_TYPESELASTODYNAMIC_H
#define PASSMO_TYPESELASTODYNAMIC_H

#include <arcane/ItemGroup.h>

struct TypesElastodynamic
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
};

#endif //PASSMO_TYPESELASTODYNAMIC_H