//+


h_size = 0.05;
G = 0.1; 
W = G; 
N_rods = 4;
L_rods  = 0.6;
L_box = L_rods + 4*G;
H_box = N_rods * W * 2  + N_rods * W * 2 + 3 * G ; 


//=========================================
// EXTERNAL BOX
//=========================================
Point(1) = {0, 0, 0, h_size};
Point(2) = {L_box, 0, 0, h_size};
Point(3) = {L_box, H_box, 0, h_size};
Point(4) = {0, H_box, 0, h_size};

Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,1};

//=========================================
// Capacitor 1
//=========================================
Point(5) = {G, G, 0, h_size};
Point(6) = {G+L_rods, G, 0, h_size};
Point(7) = {G+L_rods, G+W, 0, h_size};
Point(8) = {G+W, G+W, 0, h_size};
Point(9) = {G+W, G+W+G+W+G, 0, h_size};
Point(10) = {G+W+L_rods-W, G+W+G+W+G, 0, h_size};
Point(11) = {G+W+L_rods-W, G+W+G+W+G+W, 0, h_size};
Point(12) = {G+W, G+W+G+W+G+W, 0, h_size};
Point(13) = {G+W, G+W+G+W+G+W+G+W+G, 0, h_size};
Point(14) = {G+W+L_rods-W, G+W+G+W+G+W+G+W+G, 0, h_size};
Point(15) = {G+W+L_rods-W, G+W+G+W+G+W+G+W+G+W, 0, h_size};
Point(16) = {G+W, G+W+G+W+G+W+G+W+G+W, 0, h_size};
Point(17) = {G+W, G+W+G+W+G+W+G+W+G+W+G+W+G, 0, h_size};
Point(18) = {G+W+L_rods-W, G+W+G+W+G+W+G+W+G+W+G+W+G, 0, h_size};
Point(19) = {G+W+L_rods-W, G+W+G+W+G+W+G+W+G+W+G+W+G+W, 0, h_size};
Point(20) = {G+W, G+W+G+W+G+W+G+W+G+W+G+W+G+W, 0, h_size};
Point(21) = {G+W, G+W+G+W+G+W+G+W+G+W+G+W+G+W+G+W+G, 0, h_size};
Point(22) = {G+W+L_rods-W, G+W+G+W+G+W+G+W+G+W+G+W+G+W+G+W+G, 0, h_size};
Point(23) = {G+W+L_rods-W, G+W+G+W+G+W+G+W+G+W+G+W+G+W+G+W+G+W, 0, h_size};
Point(24) = {G, G+W+G+W+G+W+G+W+G+W+G+W+G+W+G+W+G+W, 0, h_size};

Line(5)  = {5,6};
Line(6)  = {6,7};
Line(7)  = {7,8};
Line(8)  = {8,9};
Line(9)  = {9,10};
Line(10)  = {10,11};
Line(11)  = {11,12};
Line(12)  = {12,13};
Line(13)  = {13,14};
Line(14)  = {14,15};
Line(15)  = {15,16};
Line(16)  = {16,17};
Line(17)  = {17,18};
Line(18)  = {18,19};
Line(19)  = {19,20};
Line(20)  = {20,21};
Line(21)  = {21,22};
Line(22)  = {22,23};
Line(23)  = {23,24};
Line(24)  = {24,5};

//=========================================
// Capacitor 2
//=========================================

Point(25) = {3*G, 3*G, 0, h_size};
Point(26) = {3*G+L_rods, 3*G, 0, h_size};
Point(27) = {3*G+L_rods, 3*G + 7*W + 6*G, 0, h_size};
Point(28) = {3*G, 3*G + 7*W + 6*G, 0, h_size};
Point(29) = {3*G, 3*G + 6*W + 6*G, 0, h_size};
Point(30) = {3*G+L_rods-W, 3*G + 6*W + 6*G, 0, h_size};
Point(31) = {3*G+L_rods-W, 3*G + 5*W + 4*G, 0, h_size};
Point(32) = {3*G, 3*G + 5*W + 4*G, 0, h_size};
Point(33) = {3*G, 3*G + 4*W + 4*G, 0, h_size};
Point(34) = {3*G+L_rods-W, 3*G + 4*W + 4*G, 0, h_size};
Point(35) = {3*G+L_rods-W, 3*G + 3*W + 2*G, 0, h_size};
Point(36) = {3*G, 3*G + 3*W + 2*G, 0, h_size};
Point(37) = {3*G, 3*G + 2*W + 2*G, 0, h_size};
Point(38) = {3*G+L_rods-W, 3*G + 2*W + 2*G, 0, h_size};
Point(39) = {3*G+L_rods-W, 3*G + W, 0, h_size};
Point(40) = {3*G, 3*G + W, 0, h_size};


Line(25)  = {25,26};
Line(26)  = {26,27};
Line(27)  = {27,28};
Line(28)  = {28,29};
Line(29)  = {29,30};
Line(30)  = {30,31};
Line(31)  = {31,32};
Line(32)  = {32,33};
Line(33)  = {33,34};
Line(34)  = {34,35};
Line(35)  = {35,36};
Line(36)  = {36,37};
Line(37)  = {37,38};
Line(38)  = {38,39};
Line(39)  = {39,40};
Line(40)  = {40,25};


//=========================================
// Plane Surface
//=========================================

Curve Loop(1) = {3, 4, 1, 2};
Curve Loop(2) = {23, 24, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
Curve Loop(3) = {27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 25, 26};

Plane Surface(1) = {1, 2, 3};



//=========================================
// Label
//=========================================

Physical Curve("Outer", 41) = {1, 2, 3, 4};

Physical Curve("capacitor1", 42) = {23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 24};

Physical Curve("capacitor2", 43) = {27, 26, 25, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28};

Physical Surface("inside", 44) = {1};
