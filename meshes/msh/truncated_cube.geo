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

     This is a Gmsh .geo file which produces a 3D truncated cube

     compile-run: gmsh -3  -setnumber threads 4 truncated_cube.geo

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
Box(1) = {0, 0, 0, 1, 1, 1};
Box(2) = {0.5, 0.5, 0.5, 1, 1, 1};

BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

//==============================================================================
// ---- define physical tags ----
//==============================================================================

Physical Surface("horizontal", 22) = {8};
Physical Surface("verticalXY", 23) = {9};
Physical Surface("verticalYZ", 24) = {7};
Physical Surface("top", 25) = {4};
Physical Surface("bottom", 26) = {2};
Physical Surface("left", 27) = {1};
Physical Surface("right", 28) = {6};
Physical Surface("front", 29) = {3};
Physical Surface("back", 30) = {5};
Physical Point("center", 31) = {14};
Physical Volume("volume", 32) = {1};

//==============================================================================
// ---- msh version/algo imposed ----
//==============================================================================

General.NumThreads = threads;

Mesh.Algorithm = 5;
Mesh.Algorithm3D = 10;
Mesh.MshFileVersion = 4.1;

