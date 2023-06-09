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
  enum eBoundaryCondition {
    UX, //!< Fixed X-Displacement
    UY, //!< Fixed Y-Displacement
    UZ, //!< Fixed Z-Displacement
    VX, //!< Fixed X-Velocity
    VY, //!< Fixed Y-Velocity
    VZ, //!< Fixed Z-Velocity
    AX, //!< Fixed X-Acceleration
    AY, //!< Fixed Y-Acceleration
    AZ, //!< Fixed Z-Acceleration
    FX, //!< Fixed X-Force
    FY, //!< Fixed Y-Force
    FZ, //!< Fixed Z-Force
    Unknown //!< Unknown Type
  };
  enum eNodeCondition {
    Displ,
    Vel,
    Acc,
    Force,
    UnknownCond
  };
  enum eCellCondition {
    Strain,
    Stress,
    UnknownCellCond
  };
  enum eGaussIntegrationOrder {
    Nint1, //!< Integration order along 1st local cell direction
    Nint2, //!< Integration order along 2nd local cell direction
    Nint3, //!< Integration order along 3rd local cell direction
    NoOrder //!< Unknown Type
  };
  enum eAnalysisType {
    PlaneStrain = 0, //!< Plane strain analysis (2D)
    PlaneStress,//!< Plane stress analysis (2D)
    Axi, //!< Axisymmetric analysis (2D)
    ThreeD, //!< 3D analysis
    NoAnalysis //!< Unknown type
  };
};

#endif //PASSMO_TYPESELASTODYNAMIC_H
