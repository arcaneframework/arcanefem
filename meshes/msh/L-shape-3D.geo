//-----------------------------------------------------------------------------
//
// Name       : L-shape-3D.geo
// Author     : Mohd Afeef BADRI
// Date       : 03 / Feb / 2024
//
// ----------------------------------------------------------------------------
// Comment    : Simple problem  L-shape mesh in 3D. The boundaries of the mesh
//              are named "bot", "bc", "other", besides volume is named "vol".
//
// Usage      : gmsh L-shape-3D.geo -3
//
// ----------------------------------------------------------------------------

ls=100.;
baselength=250.0;
basewidth=100.0;
baseheight=500.0;
toplength=500.0;
Lheight=250.0;
Lwidth=250.0;

Point(newp) = {0, 0, 0, ls};
Point(newp) = {baselength, 0, 0, ls};
Point(newp) = {0, basewidth, 0, ls};
Point(newp) = {baselength, basewidth, 0, ls};
Point(newp) = {0, 0, baseheight, ls};
Point(newp) = {0, basewidth, baseheight, ls};
Point(newp) = {toplength, 0, baseheight, ls};
Point(newp) = {toplength, basewidth, baseheight, ls};
Point(newp) = {toplength, 0, Lheight, ls};
Point(newp) = {toplength, basewidth, Lheight, ls};
Point(newp) = {Lwidth, 0, Lheight, ls};
Point(newp) = {Lwidth, basewidth, Lheight, ls};
Point(newp) = {Lwidth+(Lwidth-30.), 0, Lheight, ls};
Point(newp) = {Lwidth+(Lwidth-30.), basewidth, Lheight, ls};


Line(1) = {1, 3};
Line(2) = {3, 4};
Line(3) = {4, 2};
Line(4) = {2, 1};
Line(5) = {1, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};
Line(8) = {6, 8};
Line(9) = {8, 7};
Line(10) = {7, 5};
Line(11) = {8, 10};
Line(12) = {10, 9};
Line(13) = {9, 7};
Line(14) = {9, 13};
Line(15) = {13, 14};
Line(16) = {14, 10};
Line(17) = {13, 11};
Line(18) = {11, 12};
Line(19) = {12, 14};
Line(20) = {2, 11};
Line(21) = {12, 4};

Curve Loop(1) = {5, -10, -13, 14, 17, -20, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {8, 11, -16, -19, 21, -2, -7};
Plane Surface(2) = {2};

Curve Loop(3) = {8, 9, 10, 6};
Plane Surface(3) = {3};

Curve Loop(4) = {6, 7, -1, 5};
Plane Surface(4) = {4};

Curve Loop(5) = {2, 3, 4, 1};
Plane Surface(5) = {5};

Curve Loop(6) = {18, 21, 3, 20};
Plane Surface(6) = {6};

Curve Loop(7) = {19, -15, 17, 18};
Plane Surface(7) = {7};

Curve Loop(8) = {12, 14, 15, 16};
Plane Surface(8) = {8};

Curve Loop(9) = {9, -13, -12, -11};
Plane Surface(9) = {9};

Surface Loop(1) = {1, 4, 3, 2, 9, 8, 7, 6, 5};
Volume(1) = {1};

Physical Volume("vol") = {1};
Physical Surface("bot") = {5};
Physical Surface("bc") = {8};
Physical Surface("other") = {1, 7, 9, 3, 2, 6, 4};

Mesh.MshFileVersion = 4.1;
Mesh 2;
