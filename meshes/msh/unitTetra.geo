//------------------------------------
//  GMSH geo file to create a single 
//  tetrahedral mesh element 
//
//  NOTE: mesh will be msh41 format
//------------------------------------

lc = 2.0;

Point(1) = {0.0 , 0.0 , 0.0 , lc};
Point(2) = {1.0 , 0.0 , 0.0 , lc};
Point(3) = {0.0 , 1.0 , 0.0 , lc};
Point(4) = {0.0 , 0.0 , 1.0 , lc};

Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 1};
Line(4) = {4, 2};
Line(5) = {2, 1};
Line(6) = {2, 3};


Curve Loop(1) = {3, 1, 2};
Plane Surface(1) = {1};

Curve Loop(2) = {2, -6, -4};
Plane Surface(2) = {2};

Curve Loop(3) = {3, -5, 6};
Plane Surface(3) = {3};

Curve Loop(4) = {1, 4, 5};
Plane Surface(4) = {4};

Surface Loop(1) = {3, 1, 4, 2};
Volume(1) = {1};

Physical Surface("xy", 7) = {3};
Physical Surface("xz", 8) = {4};
Physical Surface("yz", 9) = {1};
Physical Surface("other", 10) = {2};

Physical Point("left", 11) = {3, 2, 1};
Physical Point("right", 12) = {4};

Physical Volume("vol", 13) = {1};

//==============================================================================
// ---- msh version imposed ----
//==============================================================================

Mesh.MshFileVersion = 4.1;
Mesh 2;
