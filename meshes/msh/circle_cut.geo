/*****************************************************************************

         This file is a part of ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 08-11-2024

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This is a Gmsh .geo file which produces a 3D sphere with box cut

     compile-run: gmsh -2 -setnumber h 4  circle_cut.geo

*******************************************************************************/


SetFactory("OpenCASCADE");

//==============================================================================
// ---- define parameters for commandline ----
//==============================================================================

DefineConstant[ h = {0.2, Min 0, Max 1e10, Step 1,
                         Name "Parameters/h h"} ];


//==============================================================================
// ---- define geometry ----
//==============================================================================

Point(1) = {0, 0, 0,  h};
Point(2) = {1, 0, 0,  h};
Point(3) = {0, 1, 0,  h};
Point(4) = {-1, 0, 0,  h};
Point(5) = {0, -1, 0,  h};
Circle(1) = {3, 1, 4};
Circle(2) = {4, 1, 5};
Circle(3) = {5, 1, 2};
Line(4) = {2, 1};
Line(5) = {1, 3};
Curve Loop(1) = {5, 1, 2, 3, 4};
Plane Surface(1) = {1};

//==============================================================================
// ---- define physical tags ----
//==============================================================================

Physical Surface("domain", 6) = {1};
Physical Curve("vertical", 7) = {5};
Physical Curve("horizontal", 8) = {4};
Physical Curve("curved", 9) = {1, 2, 3};

//==============================================================================
// ---- msh version/algo imposed ----
//==============================================================================

Mesh.Algorithm = 5;
Mesh.MshFileVersion = 4.1;
