//-----------------------------------------------------------------------------
// Debug passmo mesh
// Small bar (3 q4) with imposed displ. on right bound & blocked ones on left bound
// The  boundaries  of  the bar are named (top, right, left, bottom) 
// and inner domain is "surface"
//-----------------------------------------------------------------------------

L   = 4;  // bar length
H   = 1;  // bar height
rfactor = 1.; // mesh refinement factor
h = 2./rfactor; 

//corner points
Point(newp) = {0,   0, 0, h};
Point(newp) = {L, 0, 0, h};
Point(newp) = {L, H, 0, h};
Point(newp) = {0, H, 0, h};

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

//Structured meshing: 
// 2 nodes on left/right & 4 nodes equally spaced on top/bottom
Transfinite Line {1,3} = 4;
Transfinite Line {2,4} = 2;

// Quad meshing
Transfinite Surface "*";
Recombine Surface "*";



