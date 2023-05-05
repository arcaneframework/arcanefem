//-----------------------------------------------------------------------------
// Laplace  problem ring mesh
// The  boundaries  of  the mesh are named ( inner, outer ).
//-----------------------------------------------------------------------------

SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
Circle(2) = {0, 0, 0, 1.0, 0, 2*Pi};

Curve Loop(1) = {2};
Curve Loop(2) = {1};

Plane Surface(1) = {1, 2};

Physical Curve("outer", 3) = {2};
Physical Curve("inner", 4) = {1};
Physical Surface("surface", 5) = {1};

MeshSize {1} = 0.1;
MeshSize {2} = 0.1;
