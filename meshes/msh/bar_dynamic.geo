/*****************************************************************************

         This file is a part of PSD (Parallel Structural Dynamics)
         and ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 2023‑02‑07

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This is a Gmsh .geo file which produces a 2D bar mesh for dynamic
     solver.

     compile-run: gmsh -2 bar_dynamic.geo

*******************************************************************************/

lc = 1.0/91.0;
lengthBar=1.0;
heightBar=0.1;

Point(1) = {0.0 , 0.0 , 0.0 , lc};
Point(2) = {lengthBar , 0.0 , 0.0 , lc};
Point(3) = {lengthBar , heightBar , 0.0 , lc};
Point(4) = {0.0 , heightBar , 0.0 , lc};

Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {3 , 4};
Line(4) = {4 , 1};


Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

Physical Surface("volume") = {1};
Physical Curve("surface-left") = {4};
Physical Curve("surface-top") = {3};
Physical Curve("surface-right") = {2};
Physical Curve("surface-bottom") = {1};

//==============================================================================
// ---- msh version imposed ----
//==============================================================================

Mesh.MshFileVersion = 4.1;
Mesh 2;
