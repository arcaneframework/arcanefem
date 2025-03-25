/*****************************************************************************

         This file is a part of ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 25-03-2025

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This is a Gmsh .geo file which produces a 3D spheres within sphere cut

     compile-run: gmsh -3  -setnumber threads 4 spheres_in_sphere.geo

*******************************************************************************/
//==============================================================================
// ---- define parameters for commandline ----
//==============================================================================

DefineConstant[ threads = {4, Min 0, Max 128, Step 1,
                         Name "Parameters/threads threads"} ];
                         
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 10, -Pi/2, Pi/2, 2*Pi};
Sphere(2) = {2, 2, 0, 1, -Pi/2, Pi/2, 2*Pi};
Sphere(3) = {-5.6, 5.2, 0, 1, -Pi/2, Pi/2, 2*Pi};
Sphere(4) = {-4, -5.5, 0, 1, -Pi/2, Pi/2, 2*Pi};
Sphere(5) = {0.1, -5.9, 3.3, 2, -Pi/2, Pi/2, 2*Pi};
Sphere(6) = {2.8, 6, 2.3, 2, -Pi/2, Pi/2, 2*Pi};
Sphere(7) = {-2.3, -1.7, -2.3, 2, -Pi/2, Pi/2, 2*Pi};
Sphere(8) = {2.8, 1.3, -4.7, 3, -Pi/2, Pi/2, 2*Pi};

BooleanDifference{ Volume{1}; Delete; }{ Volume{7}; Volume{8}; Volume{2}; Volume{3}; Volume{6}; Volume{5}; Volume{4}; Delete; }

Physical Volume("volume", 29) = {1};

Physical Surface("outerSphere", 30) = {9};
Physical Surface("outerSphere1", 31) = {3, 4, 2};
Physical Surface("outerSphere2", 32) = {7, 5, 6};
Physical Surface("outerSphere3", 33) = {8};

//==============================================================================
// ---- msh version/algo imposed ----
//==============================================================================

General.NumThreads = threads;

Mesh.Algorithm = 5;
Mesh.Algorithm3D = 10;
Mesh.MshFileVersion = 4.1;

