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
    DisplacementX, //!< Fixed X-Displacement
    DisplacementY, //!< Fixed Y-Displacement
    DisplacementZ, //!< Fixed Z-Displacement
    VelocityX, //!< Fixed X-Velocity
    VelocityY, //!< Fixed Y-Velocity
    VelocityZ, //!< Fixed Z-Velocity
    AccelerationX, //!< Fixed X-Acceleration
    AccelerationY, //!< Fixed Y-Acceleration
    AccelerationZ, //!< Fixed Z-Acceleration
    ForceX, //!< Fixed X-Force
    ForceY, //!< Fixed Y-Force
    ForceZ, //!< Fixed Z-Force
    Unknown //!< Unknown Type
  };
  enum eNodeCondition {
    Displacement,
    Velocity,
    Acceleration,
    Force,
    UnknownCond
  };
  enum eGaussIntegrationOrder {
    IntegOrder1, //!< Integration order along 1st local cell direction
    IntegOrder2, //!< Integration order along 2nd local cell direction
    IntegOrder3, //!< Integration order along 3rd local cell direction
    NoOrder //!< Unknown Type
  };
  enum eAnalysisType {
    PlaneStrain = 0, //!< Plane strain analysis (2D)
    PlaneStress,//!< Plane stress analysis (2D)
    AxiSym, //!< Axisymmetric analysis (2D)
    ThreeD, //!< 3D analysis
    NoAnalysis //!< Unknown type
  };
};

#endif //PASSMO_TYPESELASTODYNAMIC_H
