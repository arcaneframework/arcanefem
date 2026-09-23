/**************************************************************************************
*                                                                                     *
*           This file is a part of PSD project                                        *
*           and ArcaneFEM (finite element tool in Arcane)                             *
*                                                                                     *
*       -------------------------------------------------------------------           *
*                                                                                     *
*       Copyright 2019-2026 CEA/DES                                                   *
*                                                                                     *
*       Licensed under the Apache License, Version 2.0  (the "License");              *
*       you may not use this file except in compliance with the License.              *
*       You may obtain a copy of the License at                                       *
*                                                                                     *
*           http://www.apache.org/licenses/LICENSE-2.0                                *
*                                                                                     *
*       Unless required by applicable law or agreed to in writing, software           *
*       distributed under the License is distributed on an  "AS IS"  BASIS,           *
*       WITHOUT  WARRANTIES  OR  CONDITIONS  OF  ANY  KIND,  either express           *
*       or implied. See  the License  for  the  specific language governing           *
*       permissions and limitations under the License.                                *
*                                                                                     *
*       -------------------------------------------------------------------           *
*                                                                                     *
*                                                                                     *
* Comment: The intended geometry is a quater of a cylinder in 2D it reduces tp quater *
*          of a circular shell. The internal and external radius of the shell is  Ri  *
*          and Re. We use `-format msh2` as a Gmsh argument to output .msh file   in  *
*          in version 2 format.                                                       *
*                                                                                     *
*  Run: gmsh -2 -format msh41 quater_cylinder.geo -bin                                 *
*                                                                                     *
**************************************************************************************/

Ri = 1.;
Re = 1.3;
d = 0.03;

Point(0) = {0, 0, 0, d};
Point(1) = {Ri, 0, 0, d};
Point(2) = {Re, 0, 0, d};
Point(3) = {0, Re, 0, d};
Point(4) = {0, Ri, 0, d};

Line(1) = {1, 2};
Circle(2) = {2, 0, 3};
Line(3) = {3, 4};
Circle(4) = {4, 0, 1};

Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};

Physical Line("bottom", 1) = {1};
Physical Line("outer", 2) = {2};
Physical Line("left", 3) = {3};
Physical Line("inner", 4) = {4};
Physical Surface("domain", 1) = {6};
