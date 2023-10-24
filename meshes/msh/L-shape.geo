//-----------------------------------------------------------------------------
//
// Name       : L-shape.geo
// Author     : Mohd Afeef BADRI, Kerian ALLAIRE
// Date       : 24 / October / 2023
//
// ----------------------------------------------------------------------------
// Comment    : Poisson problem L-shape mesh. The boundaries of the mesh
//              are named (boundary).
//
// Parameters : rfactor - this is the mesh refinment factor
//
// Usage      : gmsh L-shape.geo -setnumber rfactor 11 -2 -format msh41
//
// ----------------------------------------------------------------------------


//==============================================================================
// ---- define parameters for commandline ----
//==============================================================================

DefineConstant[ rfactor= {11, Min 1, Max 1000, Step 1,
                         Name "Parameters/rfactor rfactor"} ];

// ----------------------------------------------------------------------------

L   = 1.;   // length
W   = 0.5;  // width

h1 = 1./(rfactor);

//corner points
Point(newp) = {0,   0,   0, h1};
Point(newp) = {L,   0,   0, h1};
Point(newp) = {L,   W,   0, h1};
Point(newp) = {W,   W,   0, h1};
Point(newp) = {W,   L,   0, h1};
Point(newp) = {0,   L,   0, h1};

// lines //
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// surface //
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
Plane Surface(1) = {1};

// tags //
Physical Curve("boundary", 7) = {6, 1, 2, 3, 4, 5};
Physical Surface("surface", 8) = {1};
