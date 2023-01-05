// geo file for thermal problem //

L   = 0.1;  // square length
Ccy = L/2;  // circle center height
Cr  = 0.02; // circle radius 

rfactor = 1; // mesh refinment factor
h1 = 1./(40*rfactor); 
h2 = 1./(50*rfactor); 

//corner points
Point(newp) = {0,   0, 0, h1};
Point(newp) = {L, 0, 0, h1};
Point(newp) = {L, L, 0, h1};
Point(newp) = {0, L, 0, h1};

//circle 1
Point(newp) = {Ccy,    0, 0, h2};
Point(newp) = {Ccy-Cr, 0, 0, h2};
Point(newp) = {Ccy+Cr, 0, 0, h2};

//circle 2
Point(newp) = {L, Ccy,    0, h2};
Point(newp) = {L, Ccy-Cr, 0, h2};
Point(newp) = {L, Ccy+Cr, 0, h2};

//circle 2
Point(newp) = {0, Ccy,    0, h2};
Point(newp) = {0, Ccy-Cr, 0, h2};
Point(newp) = {0, Ccy+Cr, 0, h2};

//circle 4
Point(newp) = {Ccy,    L, 0, h2};
Point(newp) = {Ccy-Cr, L, 0, h2};
Point(newp) = {Ccy+Cr, L, 0, h2};


Line(1) = {1, 6};
Line(2) = {7, 2};
Line(3) = {2, 9};
Line(4) = {10, 3};
Line(5) = {3, 16};
Line(6) = {15, 4};
Line(7) = {4, 13};
Line(8) = {12, 1};


Circle(9) = {7, 5, 6};
Circle(10) = {10, 8, 9};
Circle(11) = {15, 14, 16};
Circle(12) = {12, 11, 13};

Curve Loop(1) = {7, -12, 8, 1, -9, 2, 3, -10, 4, 5, -11, 6};

Plane Surface(1) = {1};

Physical Curve("boundary", 13) = {1, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2};
Physical Surface("material", 14) = {1};
