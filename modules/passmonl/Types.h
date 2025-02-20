//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Types.h                                                                   */
/*                                                                           */
/* Types for NLdynamic module (in axl)                                       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef PASSMONL_TYPES_H
#define PASSMONL_TYPES_H

#include <arcane/ItemGroup.h>

struct Types
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

#endif //PASSMONL_TYPES_H