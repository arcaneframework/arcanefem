// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* GaussQuadrature.h                                           (C) 2022-2024 */
/*                                                                           */
/* PASSMO : Performant Assessment for Seismic Site Modelling with finite-    */
/* element (FEM) numerical modelling approach                                */
/* Created by : E. Foerster                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef PASSMO_GAUSSQUADRATURE_H_
#define PASSMO_GAUSSQUADRATURE_H_

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "FemUtils.h"
#include <arcane/Assertion.h>
#include "Integer3std.h"

using namespace Arcane;
using namespace Arcane::FemUtils;

extern Real	REL_PREC;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Tables of abscissa & integration weights for Gauss-Legendre quadrature
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// For edges & quadrangles only
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
const Real xh{ 0.5000000000000000 };

// order 1
const Real xgauss1 = 0.;
const Real wgauss1 = 2.;

// order 2
const Real xgauss2[2] = { -1. / sqrt(3.), 1. / sqrt(3.) };
const Real wgauss2[2] = { 1., 1. };

// order 3
const Real xgauss3[3] = { -sqrt(0.6), 0., sqrt(0.6) };
const Real wgauss3[3] = { 5. / 9., 8. / 9., 5. / 9. };

// order 4
const Real xgauss4[4] = {
  -sqrt((3. + 2. * sqrt(1.2)) / 7.), -sqrt((3. - 2. * sqrt(1.2)) / 7.),
  sqrt((3. - 2. * sqrt(1.2)) / 7.), sqrt((3. + 2. * sqrt(1.2)) / 7.)
};

const Real wgauss4[4] = {
  xh - sqrt(5. / 6.) / 6., xh + sqrt(5. / 6.) / 6.,
  xh + sqrt(5. / 6.) / 6., xh - sqrt(5. / 6.) / 6.
};

// order 5
const Real xgauss5[5] = {
  -sqrt(245. + 14. * sqrt(70.)) / 21., -sqrt(245. - 14. * sqrt(70.)) / 21., 0.,
  sqrt(245. - 14. * sqrt(70.)) / 21., sqrt(245. + 14. * sqrt(70.)) / 21.
};

const Real wgauss5[5] = {
  (322. - 13 * sqrt(70.)) / 900., (322. + 13 * sqrt(70.)) / 900., 128. / 225.,
  (322. + 13 * sqrt(70.)) / 900., (322. - 13 * sqrt(70.)) / 900.
};

// order 6
const Real xgauss6[6] = {
  -0.932469514203152, -0.661209386466265,
  -0.238619186083197, 0.238619186083197,
  0.661209386466265, 0.932469514203152
};

const Real wgauss6[6] = {
  0.171324492379170, 0.360761573048139,
  0.467913934572691, 0.467913934572691,
  0.360761573048139, 0.171324492379170
};

// order 7
const Real xgauss7[7] = {
  -0.949107912342759, -0.741531185599394,
  -0.405845151377397, 0., 0.405845151377397,
  0.741531185599394, 0.949107912342759
};

const Real wgauss7[7] = {
  0.129484966168870, 0.279705391489277,
  0.381830050505119, 0.417959183673469, 0.381830050505119,
  0.279705391489277, 0.129484966168870
};

// order 8
const Real xgauss8[8] = {
  -0.960289856497536, -0.796666477413627,
  -0.525532409916329, -0.183434642495650,
  0.183434642495650, 0.525532409916329,
  0.796666477413627, 0.960289856497536
};

const Real wgauss8[8] = {
  0.101228536290376, 0.222381034453374,
  0.313706645877887, 0.362683783378362,
  0.362683783378362, 0.313706645877887,
  0.222381034453374, 0.101228536290376
};

// order 9
const Real xgauss9[9] = {
  -0.968160239507626, -0.836031107326636,
  -0.613371432700590, -0.324253423403809, 0.,
  0.324253423403809, 0.613371432700590,
  0.836031107326636, 0.968160239507626
};

const Real wgauss9[9] = {
  0.081274388361574, 0.180648160694857,
  0.260610696402935, 0.312347077040003, 0.330239355001260,
  0.312347077040003, 0.260610696402935,
  0.180648160694857, 0.081274388361574
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// For triangles only
// integration order max = 4 => number of points max = 7
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Local (reference) coordinates on the 1st edge of the triangle
const Real xg1[4][7] = {
  { 1./3., 0., 0., 0., 0., 0., 0. },
  { xh, 0., xh, 0., 0., 0., 0. },
  { 1./3., 3./5., 1./5., 1./5., 0., 0., 0. },
  { 1./3., (9. - 2. * sqrt(15.)) / 21., (6. + sqrt(15.)) / 21., (6. + sqrt(15.)) / 21.,
    (9. + 2. * sqrt(15.)) / 21., (6. - sqrt(15.)) / 21., (6. - sqrt(15.)) / 21. }
};

// Local (reference) coordinates on the 2nd edge of the triangle
const Real xg2[4][7] = {
  { 1./3., 0., 0., 0., 0., 0., 0. },
  { xh, xh, 0., 0., 0., 0., 0. },
  { 1./3., 1./5., 3./5., 1./5., 0., 0., 0. },
  { 1./3., (6. + sqrt(15.)) / 21., (9. - 2. * sqrt(15.)) / 21., (6. + sqrt(15.)) / 21., (6. - sqrt(15.)) / 21.,
    (9. + 2. * sqrt(15.)) / 21., (6. - sqrt(15.)) / 21. }
};

// Local (reference) coordinates on the 3rd edge of the triangle
/*
const Real xg3[4][7] = {
  { 1./3., 0., 0., 0., 0., 0., 0. },
  { 0., 0.5, 0.5, 0., 0., 0., 0. },
  { 1./3., 0.2, 0.2, 0.6, 0., 0., 0. },
  { 1./3., (6. + sqrt(15.)) / 21., (6. + sqrt(15.)) / 21., (9. - 2. * sqrt(15.)) / 21.,
    (6. - sqrt(15.)) / 21., (6. - sqrt(15.)) / 21., (9. + 2. * sqrt(15.)) / 21. }
};
*/

// Integration weights
const Real wg[4][7] = {
  { xh, 0., 0., 0., 0., 0., 0. },
  { 1./6., 1./6., 1./6., 0., 0., 0., 0. },
  { -27./96., 25./96., 25./96., 25./96., 0., 0., 0. },
  { 0.112500000000000, (155. + sqrt(15.)) / 2400., (155. + sqrt(15.)) / 2400., (155. + sqrt(15.)) / 2400.,
    (155. - sqrt(15.)) / 2400., (155. - sqrt(15.)) / 2400., (155. - sqrt(15.)) / 2400. }
};

//  Correspondence between integration order & number of integration points
const Integer nptg[4] = { 1, 3, 4, 7 };

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// For tetrahedrons from 4 to 10 nodes only
// integration order = 2 => number of integration points = 4
// integration order = 3 => number of integration points = 5
// integration order = 4 => number of integration points = 15
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Local (reference) coordinates along x axis
const Real xit[4] = { 0.58541020, 0.13819660, 0.13819660, 0.13819660 };

// Local (reference) coordinates along y axis
const Real yit[4] = { 0.13819660, 0.58541020, 0.13819660, 0.13819660 };

// Local (reference) coordinates along z axis
const Real zit[4] = { 0.13819660, 0.13819660, 0.58541020, 0.13819660 };

// Integration weight
const Real wgtetra = 1./24.;

const Real a2{(5. - sqrt(5.))/20.}, b2{(5. + 3.*sqrt(5.))/20.};
const Real a1{1./4.}, b3{1./6.}, c3{1./2.};
const Real b41{(7. + sqrt(15.))/34.}, b42{(7. - sqrt(15.))/34.};
const Real c41{(13. - 3.*sqrt(15.))/34.}, c42{(13. + 3.* sqrt(15.))/34.};
const Real d4{(5. - sqrt(15.))/20.}, e4{(5. + sqrt(15.))/20.};

// Local (reference) coordinates along x axis
const Real xtet[4][15] = {
  { a1, 0., 0., 0.,  //order 1
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  { a2, a2, a2, b2,  //order 2
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  { a1, b3, b3, b3, c3,  //order 3
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  { a1, b41, b41, b41, c41, b42, b42, b42, c42, d4, d4, e4, d4, e4, e4}  //order 4
};

// Local (reference) coordinates along y axis
const Real ytet[4][15] = {
  { a1, 0., 0., 0.,  //order 1
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  { a2, a2, b2, a2,  //order 2
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  { a1, b3, b3, c3, b3,  //order 3
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  { a1, b41, b41, c41, b41, b42, b42, c42, b42, d4, e4, d4, e4, d4, e4}  //order 4
};

// Local (reference) coordinates along z axis
const Real ztet[4][15] = {
  { a1, 0., 0., 0.,  //order 1
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  { a2, b2, a2, a2,  //order 2
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  { a1, b3, c3, b3, b3,  //order 3
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
  { a1, b41, c41, b41, b41, b42, c42, b42, b42, e4, d4, d4, e4, e4, d4}  //order 4
};
const Real wgtet1{ 1./6. };
const Real wgtet2{ 1./24. };
const Real wgtet3[2] = { -2./15., 3./40. };
const Real wgtet4[4] = { 8./405., (2665 - 14.*sqrt(15))/226800., (2665 + 4.*sqrt(15))/226800., 5./567.};

//  Correspondence between weight / gauss point indices for order 3
const Integer npwgtet3[5] = { 0, 1, 1, 1 };

//  Correspondence between weight / gauss point indices for order 4
const Integer npwgtet4[15] = { 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3 };

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// For Pyramids from 5 (lin) to 13 nodes (quadratic) only
// integration order = 2 => number of integration points = 5
// integration order = 3 => number of integration points = 6
// integration order = 4 => number of integration points = 27 => not implemented
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Local (reference) coordinates along x axis
const Real ap1{ 0.5702963741068025 };
const Real hp21{ 0.1531754163448146 }, hp22{ 0.6372983346207416 };
const Real hp31{1./6.}, hp32{ 0.08063183038464675 }, hp33{ 0.6098484849057127 };
const Real xpyr[2][6] = {
  { xh, 0., -xh, 0., 0., 0. },  //order 2
  { ap1, 0., -ap1, 0., 0., 0. }  //order 3
};

// Local (reference) coordinates along y axis
const Real ypyr[2][6] = {
  { 0., xh, 0., -xh, 0., 0.},  //order 2
  { 0., ap1, 0., -ap1, 0., 0.} //order 3
};

// Local (reference) coordinates along z axis
const Real zpyr[2][6] = {
  { hp21, hp21, hp21, hp21, hp22, 0.},  //order 2
  { hp31, hp31, hp31, hp31, hp32, hp33 }  //order 3
};

// Integration weight
const Real wgpyr2 { 2./15. };
const Real wgpyr3[3] { 0.1024890634400000, 0.1100000000000000, 0.1467104129066667 };

//  Correspondence between weight / gauss point indices for order 3
const Integer npwgpyr3[6] = { 0, 0, 0, 0, 1, 2 };

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// For Pentaedron (Wedge) from 6 (lin) to 15 nodes (quadratic) only
// integration order = 2 => number of integration points = 6
// integration order = 3 => number of integration points = 8
// integration order = 4 => number of integration points = 21 => not implemented
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Local (reference) coordinates along x axis
const Real apt1{1./sqrt(3.)}, apt2{ 0.577350269189626 };
const Real xpent[2][8] = {
  { -apt1, -apt1, -apt1, apt1, apt1, apt1, 0., 0. },  //order 2
  { -apt2, -apt2, -apt2, -apt2, apt2, apt2, apt2, apt2 }  //order 3
};

// Local (reference) coordinates along y axis
const Real ypent[2][8] = {
  { xh, 0., xh, xh, 0., xh, 0., 0.},  //order 2
  { 1./3., 3./5., 1./5., 1./5., 1./3., 3./5., 1./5., 1./5. }  //order 3
};

// Local (reference) coordinates along z axis
const Real zpent[2][8] = {
  { xh, xh, 0., xh, xh, 0., 0., 0.},  //order 2
  { 1./3., 1./5., 3./5., 1./5., 1./3., 1./5., 3./5., 1./5. }  //order 3
};

// Integration weight
const Real wgpent2 { 1./6.};
const Real wgpent3[2] { -27./96., 25./96. };

//  Correspondence between weight / gauss point indices for order 3
const Integer npwgpent3[8] = { 0, 1, 1, 1, 0, 1, 1, 1 };

const Integer maxnint = 9;

/*---------------------------------------------------------------------------*/
extern Real getRefPosition(const Integer& indx, const Integer& ordre);
extern Real getWeight(const Integer& indx, const Integer& ordre);

extern Real3 LineRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real LineWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 TriRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real TriWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 QuadRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real QuadWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 HexaRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real HexaWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 TetraRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real TetraWeight(const Integer3& indices, const Integer3& ordre);

extern Real3 PentaRefPosition(const Integer3& indices, const Integer3& ordre);
extern Real PentaWeight(const Integer3& indices, const Integer3& ordre);

extern Integer getNbGaussPointsfromOrder(const Int16& /*cell_type*/, const Integer& /*ninteg*/);
extern Real3 getGaussRefPosition(const ItemWithNodes& /*cell*/, const Integer& /*ninteg*/, const Int32& /*rank*/);
extern Real getGaussWeight(const ItemWithNodes& /*cell*/, const Integer& /*ninteg*/, const Int32& /*rank*/);
/*---------------------------------------------------------------------------*/
#endif // PASSMO_GAUSSQUADRATURE_H_
