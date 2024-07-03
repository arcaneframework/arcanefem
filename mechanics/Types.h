//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* TypesElastodynamic.h                                                      */
/*                                                                           */
/* Types for Elastodynamic module (in axl)                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#ifndef PASSMO_TYPESMECHANICS_H
#define PASSMO_TYPESMECHANICS_H

#include <arcane/ItemGroup.h>

struct TypesMechanics
{
  enum eAnalysisType {
    PlaneStrain = 0, //!< Plane strain analysis (2D)
    PlaneStress,//!< Plane stress analysis (2D)
    Axi, //!< Axisymmetric analysis (2D)
    ThreeD, //!< 3D analysis
    NoAnalysis //!< Unknown type
  };

  enum eLawType {
    Elastic = 0, //!< Hooke constitutive model
    DruckerP,//!< Drücker-Prager constitutive model
    Hujeux, //!< Hujeux constitutive model
  };

  // Type of elastic properties initialized on mesh cells
  enum eElastType {
    YoungNu = 0, //!< Young Modulus & Poisson Ratio
    Lame,//!< Lame coefficients (lamba & mu)
    Veloc, //!< Compression (P) & shear (S) wave velocities
    Bulk, //!< Bulk (K) and shear (G) modules
    NoElastPropType //!< Unknown type
  };
};

#endif //PASSMO_TYPESMECHANICS_H