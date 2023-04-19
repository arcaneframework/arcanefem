/*****************************************************************************

         This file is a part of PSD (Parallel Structural Dynamics)
         and ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 2023‑05‑29

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This is a Gmsh .geo file which produces a 2D soildynamics mesh which
     includes four points for applying the point source.

     compile-run: gmsh -2 -format msh41  soil_2d.geo

*******************************************************************************/

//==============================================================================
// ---- Geometry scaling factor ----
//==============================================================================

s  = 50.;       // scaling

//==============================================================================
// ---- mesh size parameters ( denoted by h )----
//==============================================================================

h_top = 0.4;         // top surface
h_bot = 4*h_top;     // bottom surface
h_dc  = h_top;       // for the double-couple


//==============================================================================
// ---- Double couple parameters ----
//==============================================================================

lcDC         =  h_dc  ;    // Force arm length for double couple
DCSouthXCord =  0.09  ;    // South coordinate X
DCSouthYCord = -0.29  ;    // South coordinate Y


//==============================================================================
// ---- Build top points for the external borders----
//==============================================================================

// --- points for Spline ---//
Point(1) = {0*s   ,  0*s   , 0*s, h_top};
Point(2) = {1*s   , -0*s   , 0*s, h_top};
Point(3) = {0.1*s , -0*s   , 0*s, h_top};
Point(4) = {0.2*s ,  0.1*s , 0*s, h_top};
Point(5) = {0.4*s ,  0*s   , 0*s, h_top};
Point(6) = {0.5*s ,  0*s   , 0*s, h_top};
Point(7) = {0.6*s , -0.1*s , 0*s, h_top};
Point(8) = {0.8*s , -0*s   , 0*s, h_top};

// --- points for bottom border ---//
Point(9)  = {-0*s , -0.4*s, -0*s, h_bot};
Point(10) = {1*s  , -0.4*s,  0*s, h_bot};

// --- points for Double-couple ---//

Point(11) = {DCSouthXCord *s            , DCSouthYCord *s           , 0.0 , lcDC}; // South
Point(12) = {DCSouthXCord *s            , (DCSouthYCord)*s + lcDC   , 0.0 , lcDC}; // North
Point(13) = {(DCSouthXCord) *s  - lcDC/2, (DCSouthYCord) *s + lcDC/2, 0.0 , lcDC}; // West
Point(14) = {(DCSouthXCord) *s  + lcDC/2, (DCSouthYCord) *s + lcDC/2, 0.0 , lcDC}; // East

//==============================================================================
// ---- Build the external borders----
//==============================================================================

Spline(1) = {1, 3, 4, 5, 6, 7, 8, 2}; // top

Line(2) = {9, 10};    // bottom
Line(3) = {10, 2};    // left
Line(4) = {1, 9};     // right


//==============================================================================
// ---- Build the surface----
//==============================================================================

Curve Loop(1) = {4, 2, 3, -1};
Plane Surface(1) = {1};

//==============================================================================
// ---- embed double couple points----
//==============================================================================

Point{13, 11, 14, 12} In Surface{1};

//==============================================================================
// ---- mesh guidance for double couple ----
//==============================================================================

Line(5) = {14, 12};
Line(6) = {13, 11};


//==============================================================================
// ---- mesh physical groups ----
//==============================================================================

Physical Curve("surfaceLeft", 7) = {4};
Physical Curve("surfaceTop", 8) = {1};
Physical Curve("surfaceRight", 9) = {3};
Physical Curve("surfaceBottom", 10) = {2};
Physical Curve("fault1", 11) = {6};
Physical Curve("fault2", 12) = {5};
Physical Surface("volume", 13) = {1};
Physical Point("DcNorthPointCord", 14) = {12};
Physical Point("DcSouthPointCord", 15) = {11};
Physical Point("DcEastPointCord", 16) = {14};
Physical Point("DcWestPointCord", 17) = {13};
