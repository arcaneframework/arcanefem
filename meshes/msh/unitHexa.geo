//------------------------------------
//  GMSH geo file to create a single 
//  hexahedral mesh element 
//
//  NOTE: mesh will be msh41 format
//------------------------------------

lc = 2.0;

Point(1) = {0.0 , 0.0 , 0.0 , lc};
Point(2) = {1.0 , 0.0 , 0.0 , lc};
Point(3) = {0.0 , 1.0 , 0.0 , lc};
Point(4) = {1.0 , 1.0 , 0.0 , lc};
Point(5) = {0.0 , 0.0 , 1.0 , lc};
Point(6) = {1.0 , 0.0 , 1.0 , lc};
Point(7) = {0.0 , 1.0 , 1.0 , lc};
Point(8) = {1.0 , 1.0 , 1.0 , lc};


Line(1) = {1, 5};
Line(2) = {5, 7};
Line(3) = {7, 3};
Line(4) = {3, 1};
Line(5) = {5, 6};
Line(6) = {6, 8};
Line(7) = {8, 7};
Line(8) = {8, 4};
Line(9) = {4, 3};
Line(10) = {4, 2};
Line(11) = {2, 6};
Line(12) = {2, 1};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

Curve Loop(2) = {4, -12, -10, 9};
Plane Surface(2) = {2};

Curve Loop(3) = {10, 11, 6, 8};
Plane Surface(3) = {3};

Curve Loop(4) = {2, -7, -6, -5};
Plane Surface(4) = {4};

Curve Loop(5) = {12, 1, 5, -11};
Plane Surface(5) = {5};

Curve Loop(6) = {3, -9, -8, 7};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 5, 4, 6, 3};
Volume(1) = {1};

Transfinite Line "*" = 2 Using Bump 0.5;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Physical Surface("left", 13) = {2};
Physical Surface("right", 14) = {4};
Physical Surface("top", 15) = {3};
Physical Surface("bot", 16) = {1};
Physical Surface("side", 17) = {5, 6};

Physical Volume("vol", 18) = {1};

Physical Point("left", 19) = {1, 2, 3, 4};
Physical Point("right", 20) = {6, 5, 7, 8};

Mesh.MshFileVersion = 4.1;
Mesh 2;
