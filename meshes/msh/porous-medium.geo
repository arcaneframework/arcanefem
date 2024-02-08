//-------- GMSH geo file to create a 2D porous-medium mesh -------//
//

//-----------------------------------------------------------------------------
//
// Name       : porous-medium.geo
// Author     : Mohd Afeef BADRI, Kerian ALAIRE
// Date       : 15 / December / 2023
//
// ----------------------------------------------------------------------------
// Comment    : random porous medium geometry with 11 different pores
//
//
// Usage      : gmsh -2 porous-medium.geo -clmin 0.34 -clmax 0.34
//
//              gmsh -2 porous-medium.geo -clmin 0.1 -clmax 0.1
//
//
//-----------------------------------------------------------------------------

//==============================================================================
// ---- mesh size factor ----
//==============================================================================

h = 0.34;

SetFactory("OpenCASCADE");

//==============================================================================
// ---- External boundary ----
//==============================================================================

Point(1) = {0, 0, 0, h};
Point(2) = {10, 0, 0, h};
Point(3) = {10, 10, 0, h};
Point(4) = {0, 10, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


//==============================================================================
// ---- random pores ----
//==============================================================================
Circle(5) = {2.5, 1.5, -0, 1.2, 0, 2*Pi}; 
Circle(6) = {5, 5.5, -0, 1.3, 0, 2*Pi};
Circle(7) = {2.5, 6.5, -0, 0.5, 0, 2*Pi};
Circle(8) = {5.5, 7.5, -0, 0.4, 0, 2*Pi};
Circle(9) = {7.5, 3.5, -0, 0.8, 0, 2*Pi}; 
Circle(10) = {7.9, 8, -0, 0.5, 0, 2*Pi};
Circle(11) = {1.5, 8.5, -0, 1.0, 0, 2*Pi};
Circle(12) = {9., .8, -0, 0.4, 0, 2*Pi};
Circle(13) = {5.2, 1.8, -0, 0.8, 0, 2*Pi};
Circle(14) = {8.8, 6.2, -0, 0.9, 0, 2*Pi};
Circle(15) = {2.0, 4.2, -0, 1., 0, 2*Pi};

Curve Loop(1) = {1,2,3,4};
Curve Loop(2) = {5};
Curve Loop(3) = {6};
Curve Loop(4) = {7};
Curve Loop(5) = {8};
Curve Loop(6) = {9};
Curve Loop(7) = {10};
Curve Loop(8) = {11};
Curve Loop(9) = {12};
Curve Loop(10) = {13};
Curve Loop(11) = {14};
Curve Loop(12) = {15};

Plane Surface(1) = {1, 2, 3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12};

//==============================================================================
// ---- Physical labels ----
//==============================================================================

Physical Curve("wall",33) = {1, 2, 3 ,4};
Physical Curve("pore1",34) = {5};
Physical Curve("pore2",35) = {6};
Physical Curve("pore3",36) = {7};
Physical Curve("pore4",37) = {8};
Physical Curve("pore5",38) = {9};
Physical Curve("pore6",39) = {10};
Physical Curve("pore7",40) = {11};
Physical Curve("pore8",41) = {12};
Physical Curve("pore9",42) = {13};
Physical Curve("pore10",43) = {14};
Physical Curve("pore11",44) = {15};
Physical Surface("medium", 45) = {1};

//==============================================================================
// ---- msh version imposed ----
//==============================================================================

Mesh.MeshSizeFromCurvature = 13;
Mesh.MshFileVersion = 4.1;
Mesh 2;
