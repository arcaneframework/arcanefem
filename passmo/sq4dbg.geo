//-----------------------------------------------------------------------------
// Debug passmo mesh
// A small square (4x4) with traction loading on top edge and paraxial boundaries
// elsewhere.
// The  boundaries  of  the square are named (top, right, left, bottom) 
// and inner domain is "surface"
//-----------------------------------------------------------------------------

L   = 4;  // square length
rfactor = 2.; // mesh refinment factor
h = 1./rfactor; 

//corner points
Point(newp) = {0,   0, 0, h};
Point(newp) = {L, 0, 0, h};
Point(newp) = {L, L, 0, h};
Point(newp) = {0, L, 0, h};

//edges
Line(newl) = {1, 2};
Line(newl) = {2, 3};
Line(newl) = {3, 4};
Line(newl) = {4, 1};

Curve Loop(1) = {4, 1, 2, 3};

Plane Surface(1) = {1};

Physical Surface("surface") = {1};

Physical Curve("left") = {4};
Physical Curve("top") = {3};
Physical Curve("right") = {2};
Physical Curve("bottom") = {1};

Physical Point("botLeft", 6) = {1};
Physical Point("topLeft", 7) = {4};
Physical Point("topRight", 8) = {3};
Physical Point("botRight", 9) = {2};

//Structured meshing: 4 nodes equally spaced per line
Transfinite Line {1,3} = 4;
Transfinite Line {2,4} = 4;

// Quad meshing
Transfinite Surface "*";
Recombine Surface "*";



