/*****************************************************************************

         This file is a part of PSD (Parallel Structural Dynamics)
         and ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 06/July/2023

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This is a Gmsh .geo file which produces a 2D square mesh.

     compile-run: gmsh -2 -format msh41  square_-2pi_to_2pi.geo

*******************************************************************************/

//==============================================================================
// ---- parameters for corner points of the square ----
//==============================================================================

x0 = -2*Pi;
x1 =  2*Pi;
y0 = -2*Pi;
y1 =  2*Pi;

//==============================================================================
// ---- mesh size parameters ( denoted by lc )----
//==============================================================================

lc = 4*Pi/17.0;

//==============================================================================
// ---- corner mesh points of the square ----
//==============================================================================

Point(1) = {x0 , y0 , 0.0 , lc};
Point(2) = {x1 , y0 , 0.0 , lc};
Point(3) = {x1 , y1 , 0.0 , lc};
Point(4) = {x0 , y1 , 0.0 , lc};

//==============================================================================
// ---- edges of the square ----
//==============================================================================

Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {3 , 4};
Line(4) = {4 , 1};

//==============================================================================
// ---- surface of the square ----
//==============================================================================

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};


//==============================================================================
// ---- mesh physical groups ----
//==============================================================================

Physical Surface("volume") = {1};

Physical Curve("left") = {4};
Physical Curve("top") = {3};
Physical Curve("right") = {2};
Physical Curve("bottom") = {1};

Physical Point("botLeft", 6) = {1};
Physical Point("topLeft", 7) = {4};
Physical Point("topRight", 8) = {3};
Physical Point("botRight", 9) = {2};
