//------------------------------------
//  GMSH geo file to create a single 
//  Quadrangular mesh element 
//
//  NOTE: mesh will be msh41 format
//------------------------------------

lc = 2.0;

Point(1) = {0.0 , 0.0 , 0.0 , lc};
Point(2) = {1.0 , 0.0 , 0.0 , lc};
Point(3) = {0.0 , 1.0 , 0.0 , lc};
Point(4) = {1.0 , 1.0 , 0.0 , lc};


Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

Physical Curve("bot", 5) = {1};
Physical Curve("right", 6) = {2};
Physical Curve("top", 7) = {3};
Physical Curve("left", 8) = {4};
Physical Surface("surface", 9) = {1};
Physical Point("left", 10) = {1, 3};
Physical Point("right", 11) = {2, 4};

Mesh.MshFileVersion = 4.1;
Mesh 2;
