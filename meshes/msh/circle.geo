//-----------------------------------------------------------------------------
//
// Name       : circle.geo
// Author     : Mohd Afeef BADRI
// Date       : 19 / April / 2023
//
// ----------------------------------------------------------------------------
// Comment    : simple circle problem mesh
//
// Parameters : radius  - this is the radius of the circle
//              meshSize- this is the mesh size of the
//
// Usage      : gmsh circle.geo -setnumber radius .5 -2 -format msh41
//
//
//-----------------------------------------------------------------------------

//==============================================================================
// ---- define parameters for commandline ----
//==============================================================================

DefineConstant[ radius = {0.5, Min .0001, Max 1000, Step 1,
                         Name "Parameters/radius radius"} ];

DefineConstant[ meshSize  = {0.4, Min 10, Max 100, Step 1,
                         Name "Parameters/meshSize meshSize"} ];

//==============================================================================
// ---- mesh size factor ----
//==============================================================================

h = meshSize;

//==============================================================================
// ---- build circle via OCC ----
//==============================================================================

SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, radius, 0, 2*Pi};
Curve Loop(1) = {1};
Plane Surface(1) = {1};
MeshSize {1} = h;

//==============================================================================
// ---- add groups ----
//==============================================================================

Physical Curve("boundary", 2) = {1};
Physical Surface("surface", 3) = {1};
