//-----------------------------------------------------------------------------
//
// Name       : NACA0012.geo
// Author     : Mohd Afeef BADRI
// Date       : 09 / Feb. / 2023
//
// ----------------------------------------------------------------------------
// Comment    : Aerodynamics Kutta–Joukowski problem mesh
//
// Parameters : O_Grid_Radius  - this is the radius of O_Grid
//              Np_Airfoil     - number of points on top/bottom of airfoil
//              h_Airfoil      - mesh size at the airfoil boundary
//              h_Farfield     - mesh size at the Farfield boundary
//
// Usage      : gmsh NACA0012.geo -setnumber O_Grid_Radius 5 -2 -format msh41
//
//              gmsh NACA0012.geo -setnumber O_Grid_Radius 5        \
//                -setnumber h_Airfoil .01 -setnumber h_Farfield 80 \
//                -2 -format msh41
//
//
// Note       : A NACA 0012 airfoil is meshed with O-grid,  you should have  a
//              mesh with physical boundaries [upperAirfoil, lowerAirfoil] and
//              physical surface [fluid].
//
//-----------------------------------------------------------------------------

SetFactory("OpenCASCADE");

DefineConstant[ O_Grid_Radius = {2, Min .1, Max 100, Step 1,
                         Name "Parameters/O_Grid Radius"} ];

DefineConstant[ Np_Airfoil = {35, Min 35, Max 1000, Step 1,
                         Name "Parameters/Airfoil Points"} ];

DefineConstant[ h_Airfoil = {0.5, Min .0001, Max 1000, Step 1,
                         Name "Parameters/Airfoil MeshSize"} ];

DefineConstant[ h_Farfield = {20, Min 10, Max 10000000, Step 1,
                         Name "Parameters/Farfield MeshSize"} ];

//===============================================================================
// Create upper profile of the airfoil
//===============================================================================

upperAirfoil[] = {};
NpLast  = 0;
NpStart = 1;
For t In {1:Np_Airfoil}
	tt = (t-1)/(Np_Airfoil-1);
	x  = tt;
	y  = 0.17735*Sqrt(tt)-0.075597*tt - 0.212836*(tt*tt) + 0.17363*(tt*tt*tt) - 0.06254*(tt*tt*tt*tt);
	Np = newp; Point(Np) = {x , y, 0, h_Airfoil};
	upperAirfoil[] += {Np};
	NpLast = Np;
EndFor

Spline(1) = {upperAirfoil[]};
Physical Curve("upperAirfoil", 3) = {1};

//===============================================================================
// Create lower profile of the airfoil
//===============================================================================

lowerAirfoil[] = {};
lowerAirfoil[] += {NpLast};
For t In {1:Np_Airfoil-2}
	tt =  (1 - 1/(Np_Airfoil-1)) - ((t-1)/(Np_Airfoil-1));
	x  = tt;
	y  = -( 0.17735*Sqrt(tt)-0.075597*tt- 0.212836*(tt*tt) + 0.17363*(tt*tt*tt) - 0.06254*(tt*tt*tt*tt));
	Np = newp; Point(Np) = {x , y, 0, h_Airfoil};
	lowerAirfoil[] += {Np};
EndFor
lowerAirfoil[] += {NpStart};

Spline(2) = {lowerAirfoil[]};
Physical Curve("lowerAirfoil", 4) = {2};

//===============================================================================
// Create external O grid circle
//===============================================================================

Circle(3) = {0, 0, 0, O_Grid_Radius, 0, 2*Pi};
Physical Curve("FarField", 5) = {3};

//===============================================================================
// Create the internal fluid surface
//===============================================================================
Curve Loop(1) = {3};
Curve Loop(2) = {1, 2};
Plane Surface(1) = {1, 2};
Physical Surface("fluid", 6) = {1};
Mesh.MeshSizeFromCurvature = h_Farfield;
