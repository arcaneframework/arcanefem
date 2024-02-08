//-------- GMSH geo file to create a 2D bar mesh -------//
//
//
// NOTE: Use the -format msh41 format to generate the mesh

lc = 1.0/3.0;
lengthBar=5.0;

Point(1) = {0.0 , 0.0 , 0.0 , lc};
Point(2) = {lengthBar , 0.0 , 0.0 , lc};
Point(3) = {lengthBar , 1.0 , 0.0 , lc};
Point(4) = {0.0 , 1.0 , 0.0 , lc};

Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {3 , 4};
Line(4) = {4 , 1};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

Physical Surface("volume") = {1};

Physical Curve("left") = {4};
Physical Curve("top") = {3};
Physical Curve("right") = {2};
Physical Curve("bottom") = {1};

Physical Point("botLeft", 6) = {1};
Physical Point("topLeft", 7) = {4};
Physical Point("topRight", 8) = {3};
Physical Point("botRight", 9) = {2};

//==============================================================================
// ---- msh version imposed ----
//==============================================================================

Mesh.MshFileVersion = 4.1;
Mesh 2;
