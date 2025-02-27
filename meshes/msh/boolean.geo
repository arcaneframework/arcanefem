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

     compile-run: gmsh -3 boolean.geo -format msh41

*******************************************************************************/

SetFactory("OpenCASCADE");

// from http://en.wikipedia.org/wiki/Constructive_solid_geometry

Mesh.Algorithm = 6;
Mesh.MeshSizeMin = 0.4;
Mesh.MeshSizeMax = 0.4;

R = DefineNumber[ 1.4 , Min 0.1, Max 2, Step 0.01,
  Name "Parameters/Box dimension" ];
Rs = DefineNumber[ R*.7 , Min 0.1, Max 2, Step 0.01,
  Name "Parameters/Cylinder radius" ];
Rt = DefineNumber[ R*1.25, Min 0.1, Max 2, Step 0.01,
  Name "Parameters/Sphere radius" ];

Box(1) = {-R,-R,-R, 2*R,2*R,2*R};

Sphere(2) = {0,0,0,Rt};

BooleanIntersection(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

Cylinder(4) = {-2*R,0,0, 4*R,0,0, Rs};
Cylinder(5) = {0,-2*R,0, 0,4*R,0, Rs};
Cylinder(6) = {0,0,-2*R, 0,0,4*R, Rs};

BooleanUnion(7) = { Volume{4}; Delete; }{ Volume{5,6}; Delete; };
BooleanDifference(8) = { Volume{3}; Delete; }{ Volume{7}; Delete; };

Physical Surface("left", 42) = {4};
Physical Surface("right", 43) = {6};
Physical Surface("top", 44) = {7};
Physical Surface("bottom", 45) = {8};
Physical Surface("front", 46) = {5};
Physical Surface("back", 47) = {1};
Physical Volume("vol", 48) = {8};