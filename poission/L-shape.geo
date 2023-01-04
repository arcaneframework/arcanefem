//-----------------------------------------------------------------------------
// Poission problem L-shape mesh
// The  boundaries  of  the mesh are named ( boundary ).
//-----------------------------------------------------------------------------

L   = 1.;   // length
W   = 0.5;  // width

rfactor = 11; // mesh refinment factor
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
