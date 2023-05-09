//-----------------------------------------------------------------------------
//
// Name       : box-rod.geo
// Author     : Mohd Afeef BADRI
// Date       : 09 / May / 2023
//
// ----------------------------------------------------------------------------
// Comment    : simple rods in box problem mesh for electrostatics
//
// Parameters : length  - this is the length of the surronding box
//              meshSize- this is the mesh size of the
//
// Usage      : gmsh box-rod.geo -setnumber meshSize 1.0 -2 -format msh41
//
//
//-----------------------------------------------------------------------------

//==============================================================================
// ---- define parameters for commandline ----
//==============================================================================

DefineConstant[ length = {5.0, Min .0001, Max 1000, Step 1,
                         Name "Parameters/length length"} ];

DefineConstant[ meshSize  = {0.5, Min 0.0001, Max 30, Step 1,
                         Name "Parameters/MeshSize MeshSize"} ];

DefineConstant[ rodDistance = {1., Min 1, Max 20, Step 1,
                         Name "Parameters/rodDistance rodDistance"} ];

DefineConstant[ rodWidth = {0.5, Min .1, Max 20, Step 1,
                         Name "Parameters/rodWidth rodWidth"} ];

DefineConstant[ rodHeight = {3., Min .1, Max 20, Step 1,
                         Name "Parameters/rodHeight rodHeight"} ];

//==============================================================================
// ---- mesh size factor ----
//==============================================================================

h = meshSize;

//==============================================================================
// ---- Soil box parameters ----
//==============================================================================

lcBox        = h       ;       // Mesh size parameter
lcRods       = h/2.    ;       // Mesh size parameter
Lbox         = length  ;       // Length of soil box
Hbox         = length  ;       // Width of soil box
Rdistance    = rodDistance ;

//==============================================================================
// ---- Build external box ----
//==============================================================================

Point(1) = {0.0 , 0.0 , 0.0 , lcBox};
Point(2) = {Lbox , 0.0 , 0.0 , lcBox};
Point(3) = {Lbox , Lbox , 0.0 , lcBox};
Point(4) = {0.0 ,  Lbox , 0.0 , lcBox};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


//==============================================================================
// ---- Build internall rods ----
//==============================================================================

// rod 1 //
Point(5) = {(Lbox/2. -  rodDistance/2.)-rodWidth, (length-rodHeight)/2., 0, lcRods};
Point(6) = {Lbox/2. -  rodDistance/2., (length-rodHeight)/2., 0, lcRods};
Point(7) = {Lbox/2. -  rodDistance/2., (length-rodHeight)/2. + rodHeight, 0, lcRods};
Point(8) = {(Lbox/2. -  rodDistance/2.)-rodWidth, (length-rodHeight)/2. + rodHeight, 0, lcRods};

// rod 2 //
Point(9) = {Lbox/2. +  rodDistance/2., (length-rodHeight)/2. + rodHeight, 0, lcRods};
Point(10) = {(Lbox/2. +  rodDistance/2.)+rodWidth, (length-rodHeight)/2. + rodHeight, 0, lcRods};
Point(11) = {(Lbox/2. +  rodDistance/2.)+rodWidth, (length-rodHeight)/2., 0, lcRods};
Point(12) = {Lbox/2. +  rodDistance/2., (length-rodHeight)/2., 0, lcRods};

// rod 1 //
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// rod 2 //
Line(9) = {12, 11};
Line(10) = {11, 10};
Line(11) = {10, 9};
Line(12) = {9, 12};

//==============================================================================
// ---- Build surface ----
//==============================================================================

Curve Loop(1) = {4, 1, 2, 3};
Curve Loop(2) = {8, 5, 6, 7};
Curve Loop(3) = {12, 9, 10, 11};
Plane Surface(1) = {1, 2, 3};

//==============================================================================
// ---- Build groups ----
//==============================================================================

Physical Curve("external", 13) = {1, 2, 3, 4};
Physical Curve("rod1", 14) = {8, 5, 6, 7};
Physical Curve("rod2", 15) = {12, 9, 10, 11};
Physical Surface("air", 16) = {1};
