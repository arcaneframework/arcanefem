/*****************************************************************************

         This file is a part of ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 17-07-2024

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This is a Gmsh .geo file which produces a 3D sphere with box cut

     compile-run: gmsh -3  -setnumber threads 4 sphere_cut.geo

*******************************************************************************/

//==============================================================================
// ---- define parameters for commandline ----
//==============================================================================

DefineConstant[ threads = {4, Min 0, Max 128, Step 1,
                         Name "Parameters/threads threads"} ];


//==============================================================================
// ---- define geometry ----
//==============================================================================

SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
Box(2) = {0, -0, 0, 1, 1, 1};


BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }


//==============================================================================
// ---- define physical tags ----
//==============================================================================

//Physical Surface("Cut", 10) = {4, 3, 2};
Physical Surface("verticalXY", 13) = {3};
Physical Surface("verticalYZ", 14) = {4};
Physical Surface("horizontal", 15) = {2};
Physical Surface("curved", 11) = {1};
Physical Volume("domain", 12) = {1};

//==============================================================================
// ---- msh version/algo imposed ----
//==============================================================================

General.NumThreads = threads;

Mesh.Algorithm = 5;
Mesh.Algorithm3D = 10;
Mesh.MshFileVersion = 4.1;

