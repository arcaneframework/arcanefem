// geo file for thermal problem //

L   = 0.2;  // square length
Ccy = 0.1;  // circle center height
Cr  = 0.02; // circle radius 

rfactor = 11; // mesh refinment factor
h1 = 1./(40*rfactor); 
h2 = 1./(50*rfactor); 

//corners
Point(newp) = {0,   0, 0, h1};
Point(newp) = {L, 0, 0, h1};
Point(newp) = {L, L, 0, h1};
Point(newp) = {0, L, 0, h1};

//circle
Point(newp) = {0, Ccy, 0, h2};
Point(newp) = {0, Ccy-Cr, 0, h2};
Point(newp) = {0, Ccy+Cr, 0, h2};

//+
Line(newl) = {1, 2};
Line(newl) = {2, 3};
Line(newl) = {3, 4};
Line(newl) = {4, 7};
Circle(newl) = {6, 5, 7};
Line(newl) = {6, 1};

Curve Loop(1) = {6, 1, 2, 3, 4, -5};

Plane Surface(1) = {1};
//+
Physical Surface("Planchere", 7) = {1};
//+
Physical Curve("Gauche", 8) = {6, 4};
//+
Physical Curve("Droite", 9) = {2};
//+
Physical Curve("Bas", 10) = {1};
//+
Physical Curve("Haut", 11) = {3};
//+
Physical Curve("Cercle", 12) = {5};
