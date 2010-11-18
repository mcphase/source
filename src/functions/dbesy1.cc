/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselY1(x) 			      dbesy1.cc |
|                                                                               |
| MatPack Library Release 1.0                                                   |
| Copyright (C) 1991,1995 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.						|
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.							|
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.			|
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration	|
| and installation of Matpack.							|
|                                                                               |
\*-----------------------------------------------------------------------------*/

#include "../../include/mpspecfunp.h"

//-----------------------------------------------------------------------------//
//
// double BesselY1 (double x)
//
// BesselY1(x) calculates the double precision Bessel function of the 
// second kind of order for double precision argument x.
// 
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10A1, REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BY1        on the interval  0.          to  1.60000E+01 
//                                        with weighted error   8.65E-33 
//                                         log weighted error  32.06 
//                               significant figures required  32.17 
//                                    decimal places required  32.71 
//
//-----------------------------------------------------------------------------//

double BesselY1 (double x)
{
    static double by1cs[20] = { 
	 0.0320804710061190862932352018628015,
	 1.26270789743350044953431725999727,
	 0.00649996189992317500097490637314144,
	-0.0893616452886050411653144160009712,
	 0.0132508812217570954512375510370043,
	-8.97905911964835237753039508298105e-4,
	 3.64736148795830678242287368165349e-5,
	-1.00137438166600055549075523845295e-6,
	 1.99453965739017397031159372421243e-8,
	-3.02306560180338167284799332520743e-10,
	 3.60987815694781196116252914242474e-12,
	-3.48748829728758242414552947409066e-14,
	 2.78387897155917665813507698517333e-16,
	-1.86787096861948768766825352533333e-18,
	 1.06853153391168259757070336e-20,
	-5.27472195668448228943872e-23,
	 2.27019940315566414370133333333333e-25,
	-8.59539035394523108693333333333333e-28,
	 2.88540437983379456e-30,
	-8.64754113893717333333333333333333e-33 
    };

    const double twodpi = 0.636619772367581343075535053490057,
	         xsml   = sqrt(0.5 * DBL_EPSILON * 4.0),
	         xmin   = exp(MpMax(log(DBL_MIN),-log(DBL_MAX)) + 0.01) * 1.571;
    double y;

    static int nty1, first = 1;
    if (first) {
	nty1 = initds(by1cs, 20, 0.5 * DBL_EPSILON * 0.1);
	first = 0;
    }

    if (x <= 0.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "BesselY1",
		      "x is zero or negative");
	return NAN;
    }

    if (x > 4.0) goto L20;

    if (x < xmin) {
	Matpack.Error(Mat::Overflow, "%s: %s", "BesselY1",
		      "x so small Y1(x) overflows");
	return NAN;
    }

    y = 0.0;
    if (x > xsml) y = x * x;
    return twodpi * log(x * 0.5) * BesselJ1(x) 
	               + (dcsevl(y * 0.125 - 1.0, by1cs, nty1) + 0.5) / x;

  L20:
    double theta, ampl;
    d9b1mp(x, ampl, theta);
    return ampl * sin(theta);
}

//-----------------------------------------------------------------------------//
