//-----------------------------------------------------------------------------
//
// Name       : box-rod-circle.geo
// Author     : Mohd Afeef BADRI
// Date       : 09 / May / 2023
//
// ----------------------------------------------------------------------------
// Comment    : simple rod and circle in box problem mesh for electrostatics
//
// Parameters : length  - this is the length of the surronding box
//              meshSize- this is the mesh size of the
//
// Usage      : gmsh box-rod-circle.geo -setnumber meshSize 1.0 -2 -format msh41
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
Lbox         = length  ;       // Length box
Hbox         = length  ;       // Height box
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
// ---- Build internal rods ----
//==============================================================================

// rod 1 //
Point(5) = {(Lbox/2. -  rodDistance/2.)-rodWidth, (length-rodHeight)/2., 0, lcRods};
Point(6) = {Lbox/2. -  rodDistance/2., (length-rodHeight)/2., 0, lcRods};
Point(7) = {Lbox/2. -  rodDistance/2., (length-rodHeight)/2. + rodHeight, 0, lcRods};
Point(8) = {(Lbox/2. -  rodDistance/2.)-rodWidth, (length-rodHeight)/2. + rodHeight, 0, lcRods};

// circle //
Point(9) = {Lbox/2. +  rodDistance/2. + rodWidth/2., (length-rodHeight)/2. + rodHeight/2, 0, lcRods};
Point(10) = {Lbox/2. +  rodDistance/2. + rodWidth*2, (length-rodHeight)/2. + rodHeight/2, 0, lcRods};
Point(11) = {Lbox/2. +  rodDistance/2. - rodWidth, (length-rodHeight)/2. + rodHeight/2, 0, lcRods};

// rod 1 //
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// circle // 
Circle(9) = {11, 9, 10};
Circle(10) = {10, 9, 11};


//==============================================================================
// ---- Build surface ----
//==============================================================================

Curve Loop(1) = {4, 1, 2, 3};
Curve Loop(2) = {8, 5, 6, 7};
Curve Loop(3) = {10, 9};
Plane Surface(1) = {1, 2, 3};



//==============================================================================
// ---- Build groups ----
//==============================================================================

Physical Curve("external", 11) = {4, 1, 2, 3};
Physical Curve("rod1", 12) = {8, 5, 6, 7};
Physical Curve("circle", 13) = {10, 9};
Physical Surface("air", 14) = {1};
