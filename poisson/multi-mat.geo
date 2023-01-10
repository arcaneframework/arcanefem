//-----------------------------------------------------------------------------
// Thermal problem multi material mesh
// A small square (5x5) is encapsulated within a larger one (10x10)
// The boundaries of the lager square are nmes (Left,Right,Top,Bot)
// and the two squares are named as (Mat1,Mat2). 
//-----------------------------------------------------------------------------

rfactor = 1;         // mesh refinment factor
h1 = 1./(10*rfactor); // mesh size 

Point(1) = {0, 0, 0, h1};
Point(2) = {10, 0, 0, h1};
Point(3) = {10, 10, 0, h1};
Point(4) = {0, 10, 0, h1};
Point(5) = {2.5, 2.5, 0, h1};
Point(6) = {2.5, 7.5, 0, h1};
Point(7) = {7.5, 7.5, 0, h1};
Point(8) = {7.5, 2.5, 0, h1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 8};
Line(6) = {8, 7};
Line(7) = {7, 6};
Line(8) = {6, 5};

Curve Loop(1) = {4, 1, 2, 3};
Curve Loop(2) = {8, 5, 6, 7};

Plane Surface(1) = {1, 2};
Plane Surface(2) = {2};

Physical Curve("Left", 9) = {4};
Physical Curve("Right", 10) = {2};
Physical Curve("Top", 11) = {3};
Physical Curve("Bot", 12) = {1};
Physical Surface("Mat1", 13) = {1};
Physical Surface("Mat2", 14) = {2};
