//------------------------------------
//  GMSH geo file to create a single 
//  triangular mesh element 
//
//  NOTE: mesh will be msh41 format
//------------------------------------

lc = 2.0;

Point(1) = {0.0 , 0.0 , 0.0 , lc};
Point(2) = {1.0 , 0.0 , 0.0 , lc};
Point(3) = {0.0 , 1.0 , 0.0 , lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};

Curve Loop(1) = {3, 1, 2};
Plane Surface(1) = {1};

Physical Curve("base", 4) = {1};
Physical Curve("hyp", 5) = {2};
Physical Curve("perp", 6) = {3};
Physical Surface("sur", 7) = {1};

Physical Point("left", 8) = {1, 3};
Physical Point("right", 9) = {2};

Mesh.MshFileVersion = 4.1;
Mesh 2;
