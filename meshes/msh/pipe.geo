/*****************************************************************************

         This file is a part of ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 27-02-2024

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This is a Gmsh .geo file which produces 3D boolean geometry. The file 
     is a copy of Gmsh's tutorial on boolean only groups  have  been added 
     here. 

     compile-run: gmsh -3 pipe.geo -format msh41

*******************************************************************************/

SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.1;
Mesh.MeshSizeMax = 0.1;
Geometry.NumSubEdges = 100; // nicer display of curve

nturns = DefineNumber[ 1, Min 0.1, Max 1, Step 0.01, Name "Parameters/Turn" ];
npts = 20;
r = 1;
rd = 0.1;
h = 1 * nturns;

For i In {0:npts-1}
  theta = i * 2*Pi*nturns/npts;
  Point(i + 1) = {r * Cos(theta), r * Sin(theta), i * h/npts};
EndFor

Spline(1) = {1:npts};
Wire(1) = {1};

Disk(1) = {1,0,0, rd};

Rectangle(2) = {1+2*rd,-rd,0, 2*rd,2*rd,rd/5};
Rotate {{1, 0, 0}, {0, 0, 0}, Pi/2} { Surface{1,2}; }

Extrude { Surface{1,2}; } Using Wire {1}
Delete{ Surface{1,2}; }

Physical Surface("leftSurfaceCircularPipe", 47) = {5};
Physical Surface("rightSurfaceCircularPipe", 48) = {3};
Physical Surface("leftSurfaceSquarePipe", 49) = {15};
Physical Surface("rightSurfaceSquarePipe", 50) = {6};
Physical Volume("SquarePipeVol", 51) = {1};
Physical Volume("CylindricalPipeVol", 52) = {2};
