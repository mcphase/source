/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselK1(x)                             dbesk1.cc |
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
// double BesselK1 (double x);
//
// BesselK1(x) alculates the double precision modified (hyperbolic) 
// Bessel function of the third kind of order one for double precision 
// argument x.  The argument must be large enough that the result does 
// not overflow and small enough that the result does not underflow. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10B1, REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BK1        on the interval  0.          to  4.00000E+00 
//                                        with weighted error   9.16E-32 
//                                         log weighted error  31.04 
//                               significant figures required  30.61 
//                                    decimal places required  31.64 
//
//-----------------------------------------------------------------------------//

double BesselK1 (double x)
{
    static double bk1cs[16] = { 
	 0.025300227338947770532531120868533,
	-0.35315596077654487566723831691801,
	-0.12261118082265714823479067930042,
	-0.0069757238596398643501812920296083,
	-1.7302889575130520630176507368979e-4,
	-2.4334061415659682349600735030164e-6,
	-2.2133876307347258558315252545126e-8,
	-1.4114883926335277610958330212608e-10,
	-6.6669016941993290060853751264373e-13,
	-2.4274498505193659339263196864853e-15,
	-7.023863479386287597178379712e-18,
	-1.6543275155100994675491029333333e-20,
	-3.2338347459944491991893333333333e-23,
	-5.3312750529265274999466666666666e-26,
	-7.5130407162157226666666666666666e-29,
	-9.1550857176541866666666666666666e-32 
    };
    
    const double xmin = exp(MpMax(log(DBL_MIN), -log(DBL_MAX)) + 0.01),
	         xsml = sqrt(0.5 * DBL_EPSILON * 4.0),
	         xmaxt = -log(DBL_MIN),
	         xmax = xmaxt - xmaxt * 0.5 * log(xmaxt) / (xmaxt + 0.5);

    static int ntk1, first = 1;
    if (first) {
	ntk1 = initds(bk1cs, 16, 0.5 * DBL_EPSILON * 0.1);
	first = 0;
    }

    double y;

    if (x <= 0.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "BesselK1",
		      "x is zero or negative");
	return NAN;
    }

    if (x > 2.0) goto L20;

    if (x < xmin) { 
	Matpack.Error(Mat::Overflow, "%s: %s", "BesselK1",
		      "x so small K1 overflows");
	return NAN;
    }
    
    y = 0.0;
    if (x > xsml) y = x * x;

    return log(x * 0.5) * BesselI1(x) 
	          + (dcsevl( y * 0.5 - 1.0, bk1cs, ntk1) + 0.75) / x;

  L20:
    if (x > xmax) {
	Matpack.Warning(Mat::Underflow, "%s: %s", "BesselK1",
			"x so big K1(x) underflows");
	return 0.0;
    }
    return exp(-x) * BesselExpK1(x);
}

//-----------------------------------------------------------------------------//
