//-------- GMSH geo file to create a 2D plate mesh -------//
//
//
// NOTE: Use the -format msh41 format to genrate the mesh

//-----------------------------------------------------------------------------
//
// Name       : plate.fine.geo
// Author     : Mohd Afeef BADRI
// Date       : 19 / April / 2023
//
// ----------------------------------------------------------------------------
// Comment    : simple plate problem mesh for elasticity
//
// Parameters : length  - this is the length of the plate
//              width   - this is the width of the plate
//
// Usage      : gmsh plate.fine.geo -setnumber length 30.0 -2 -format msh41
//
//              gmsh plate.fine.geo -setnumber length 30.0        \
//                -setnumber width 5.0 -2 -format msh41
//
//
//-----------------------------------------------------------------------------

//==============================================================================
// ---- mesh size factor ----
//==============================================================================

h = 0.8;

//==============================================================================
// ---- define parameters for commandline ----
//==============================================================================

DefineConstant[ length = {30.0, Min .0001, Max 1000, Step 1,
                         Name "Parameters/length length"} ];

DefineConstant[ width  = {5.0, Min 10, Max 10000000, Step 1,
                         Name "Parameters/width width"} ];

//length=30.0;
//width=5.0;

Point(1) = {0.0    , 0.0   , 0.0 , h};
Point(2) = {length , 0.0   , 0.0 , h};
Point(3) = {length , width , 0.0 , h};
Point(4) = {0.0    , width , 0.0 , h};

Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {3 , 4};
Line(4) = {4 , 1};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

Physical Surface("volume") = {1};
Physical Curve("left") = {4};
Physical Curve("top") = {3};
Physical Curve("right") = {2};
Physical Curve("bottom") = {1};
Physical Point("botLeft", 6) = {1};
Physical Point("topLeft", 7) = {4};
Physical Point("topRight", 8) = {3};
Physical Point("botRight", 9) = {2};

//==============================================================================
// ---- msh version imposed ----
//==============================================================================

Mesh.MshFileVersion = 4.1;
Mesh 2;
