/*****************************************************************************

         This file is a part of PSD (Parallel Structural Dynamics)
         and ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Ali Asad
     Email    : ali.asad@cea.fr
     Date     : 21/April/2026

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This is a Gmsh .geo file which produces a 3D cube mesh.

     compile-run: gmsh -3 -setnumber threads 4 unit_cube.geo

*******************************************************************************/

//==============================================================================
// ---- parameters for corner points of the cube ----
//==============================================================================

x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
z0 = 0;
z1 = 1;

//==============================================================================
// ---- mesh size parameters ( denoted by lc )----
//==============================================================================

lc = 1.0/5;

//==============================================================================
// ---- corner mesh points of the cube ----
//==============================================================================

Point(1) = {x0 , y0 , z0 , lc};
Point(2) = {x1 , y0 , z0 , lc};
Point(3) = {x1 , y1 , z0 , lc};
Point(4) = {x0 , y1 , z0 , lc};
Point(5) = {x0 , y0 , z1 , lc};
Point(6) = {x1 , y0 , z1 , lc};
Point(7) = {x1 , y1 , z1 , lc};
Point(8) = {x0 , y1 , z1 , lc};

//==============================================================================
// ---- edges of the square ----
//==============================================================================

Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {3 , 4};
Line(4) = {4 , 1};
Line(5) = {5 , 6};
Line(6) = {6 , 7};
Line(7) = {7 , 8};
Line(8) = {8 , 5};
Line(9) = {5 , 1};
Line(10) = {6 , 2};
Line(11) = {7 , 3};
Line(12) = {8 , 4};


//==============================================================================
// ---- surface of the cube ----
//==============================================================================

Line Loop(13) = {1, 2, 3, 4}; //bottom
Plane Surface(14) = {13};
Line Loop(15) = {5, 6, 7, 8}; //top
Plane Surface(16) = {15};
Line Loop(17) = {4, -9, -8, 12}; //left
Plane Surface(18) = {17};
Line Loop(19) = {2, -11, -6, 10}; //right
Plane Surface(20) = {19};
Line Loop(21) = {3, -12, -7, 11}; //back
Plane Surface(22) = {21};
Line Loop(23) = {1, -10, -5, 9}; //front
Plane Surface(24) = {23};
Surface Loop(25) = {18, 14, 16, 22, 24, 20};
Volume(26) = {25};


//==============================================================================
// ---- mesh physical groups ----
//==============================================================================

Physical Surface("left") = {18};
Physical Surface("right") = {20};
Physical Surface("othersurfaces") = {14, 16, 22, 24};
Physical Volume("volume") = {26};

//==============================================================================
// ---- msh version/algo imposed ----
//==============================================================================

General.NumThreads = threads;

Mesh.Algorithm = 5;
Mesh.Algorithm3D = 10;
Mesh.MshFileVersion = 4.1;

