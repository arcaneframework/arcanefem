ls=1/100.;
lengthBar=1.0;
heightBar=0.1;
widthBar=0.04;

Point(1) = {lengthBar, 0, 0, ls};
Point(2) = {0, 0, 0, ls};
Point(3) = {0, heightBar, 0, ls};
Point(4) = {0, 0, widthBar, ls};
Point(5) = {0, heightBar, widthBar, ls};
Point(6) = {lengthBar, heightBar, widthBar, ls};
Point(7) = {lengthBar, 0, widthBar, ls};
Point(8) = {lengthBar, heightBar, 0, ls};

Line(1) = {5, 3};
Line(2) = {3, 8};
Line(3) = {8, 6};
Line(4) = {6, 5};
Line(5) = {5, 4};
Line(6) = {4, 2};
Line(7) = {2, 3};
Line(8) = {2, 1};
Line(9) = {1, 8};
Line(10) = {6, 7};
Line(11) = {1, 7};
Line(12) = {7, 4};
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};
Line Loop(15) = {2, -9, -8, 7};
Plane Surface(16) = {15};
Line Loop(17) = {1, -7, -6, -5};
Plane Surface(18) = {17};
Line Loop(19) = {8, 11, 12, 6};
Plane Surface(20) = {19};
Line Loop(21) = {10, -11, 9, 3};
Plane Surface(22) = {21};
Line Loop(23) = {4, 5, -12, -10};
Plane Surface(24) = {23};
Surface Loop(25) = {18, 14, 16, 22, 24, 20};
Volume(26) = {25};

Physical Surface("left") = {18};
Physical Surface("right") = {22};
Physical Surface("top") = {14};
Physical Surface("bot") = {20};
Physical Surface("sides") = {16, 24};
Physical Volume("vol") = {26};

// Hexa meshing
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
Recombine Volume "*";

Mesh.MshFileVersion = 4.1;
Mesh 3;
