/*****************************************************************************

         This file is a part of ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 11‑03-2025

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This is a Gmsh .geo file which produces a 3D soildynamics mesh which
     includes four points for applying the double couple source.

     compile-run: gmsh -3 -format msh41  cube_double_couple_3d.geo

*******************************************************************************/

dcLen = 0.5;

Point(1) = {0, 0, 0, dcLen};
Point(2) = {0, 0, dcLen, dcLen};
Point(3) = {0, 0, -dcLen, dcLen};
Point(4) = {dcLen, 0, 0, dcLen};
Point(5) = {-dcLen, 0, 0, dcLen};

Line(1) = {4, 1};
Line(2) = {1, 5};
Line(3) = {1, 2};
Line(4) = {3, 1};
Line(5) = {2, 5};
Line(6) = {5, 3};
Line(7) = {2, 4};
Line(8) = {4, 3};

SetFactory("OpenCASCADE");
Box(1) = {-25, -25, -25, 50, 50, 50};

Point{1, 2, 3, 4, 5} In Volume{1};

Physical Surface("leftsur", 21) = {2};
Physical Surface("rightsur", 22) = {1};
Physical Surface("frontsur", 23) = {4};
Physical Surface("backsur", 24) = {3};
Physical Surface("topsur", 25) = {6};
Physical Surface("botsur", 26) = {5};
Physical Volume("vol", 27) = {1};
Physical Point("dcNorth", 28) = {2};
Physical Point("dcSouth", 29) = {3};
Physical Point("dcEast", 31) = {5};
Physical Point("dcWest", 30) = {4};

Mesh.Algorithm = 5;
Mesh.Algorithm3D = 10;
Mesh.MshFileVersion = 4.1;
