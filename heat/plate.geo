//-------- GMSH geo file to create a 2D plate mesh -------//
//
//
// NOTE: Use the -format msh41 format to genrate the mesh

lc = 2.0;
length=30.0;
width=5.0;

Point(1) = {0.0    , 0.0   , 0.0 , lc};
Point(2) = {length , 0.0   , 0.0 , lc};
Point(3) = {length , width , 0.0 , lc};
Point(4) = {0.0    , width , 0.0 , lc};

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
