/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselK0(x)                             dbesk0.cc |
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
// double BesselK0 (double x);
//
// BesselK0(x) calculates the double precision modified (hyperbolic) 
// Bessel function of the third kind of order zero for double 
// precision argument 0.  The argument must be greater than zero 
// but not so large that the result underflows. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10B1, REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BK0        on the interval  0.          to  4.00000E+00 
//                                        with weighted error   3.08E-33 
//                                         log weighted error  32.51 
//                               significant figures required  32.05 
//                                    decimal places required  33.11 
//
//-----------------------------------------------------------------------------//

double BesselK0 (double x)
{
    static double bk0cs[16] = { 
	-0.0353273932339027687201140060063153,
	 0.344289899924628486886344927529213,
	 0.0359799365153615016265721303687231,
	 0.00126461541144692592338479508673447,
	 2.28621210311945178608269830297585e-5,
	 2.53479107902614945730790013428354e-7,
	 1.90451637722020885897214059381366e-9,
	 1.03496952576336245851008317853089e-11,
	 4.25981614279108257652445327170133e-14,
	 1.3744654358807508969423832544e-16,
	 3.57089652850837359099688597333333e-19,
	 7.63164366011643737667498666666666e-22,
	 1.36542498844078185908053333333333e-24,
	 2.07527526690666808319999999999999e-27,
	 2.7128142180729856e-30,
	 3.08259388791466666666666666666666e-33 
    };

    const double xsml  = sqrt(0.5 * DBL_EPSILON * 4.0),
	         xmaxt = -log(DBL_MIN),
	         xmax = xmaxt - xmaxt * 0.5 * log(xmaxt) / (xmaxt + 0.5);

    static int ntk0, first = 1;
    if (first) {
	ntk0 = initds(bk0cs, 16, 0.5 * DBL_EPSILON * 0.1);
	first = 0;
    }

    double y;

    if (x <= 0.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "BesselK0",
		      "x is zero or negative");
	return NAN;
    }

    if (x > 2.0) goto L20;
    y = 0.0;
    if (x > xsml) y = x * x;
    return  -log(x*0.5) * BesselI0(x) - 0.25 + dcsevl(y*0.5-1.0, bk0cs, ntk0);
    
  L20:
    if (x > xmax) {
	Matpack.Warning(Mat::Underflow, "%s: %s", "BesselK0",
			"x so big K0(x) underflows");
	return 0.0;
    }
    return exp(-x) * BesselExpK0(x);
}

//-----------------------------------------------------------------------------//
