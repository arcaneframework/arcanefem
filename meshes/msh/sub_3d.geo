//-----------------------------------------------------------------------------
//
// Name       : sub.geo
// Author     : Mohd Afeef BADRI
// Date       : 04 / April / 2025
//
// ----------------------------------------------------------------------------
// Comment    : simple submarine problem mesh for ill posed acoustics
//
//
// Usage      : gmsh sub_3d.geo -2 -format msh41 -clmin 1 -clmax 1
//              gmsh sub_3d.geo -2 -format msh41 -setnumber MeshSize 1.0 -setnumber Rfactor 2.0 
//
//-----------------------------------------------------------------------------

SetFactory("OpenCASCADE");

DefineConstant[ meshSize = {1.0, Min 1e-8, Max 100, Step 1, Name "MeshSize"} ];
DefineConstant[ Rfactor  = {1.0, Min 1, Max 100, Step 1, Name "Rfactor"} ];
meshSize /= Rfactor;

//==============================================================================
// Outer Sphere (Water Domain)
//==============================================================================
R_sphere = .1;
Sphere(1) = {0, 0, 0, R_sphere};

//==============================================================================
// Inner Torus (Cavity to be subtracted)
//==============================================================================
Sphere(2) = {0, 0, 0, R_sphere/10};

BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

Physical Surface("outer", 10) = {3};
Physical Surface("inner", 11) = {2};
Physical Volume("water", 12) = {1};

Mesh.MshFileVersion = 4.1;
