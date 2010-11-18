/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselI0(x) 			      dbesi0.cc |
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
// double BesselI0 (double x)
//
// BesselI0(x) calculates the double precision modified (hyperbolic)
// Bessel function of the first kind of order zero and double precision 
// argument x.
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10B1, REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BI0        on the interval  0.          to  9.00000E+00 
//                                        with weighted error   9.51E-34 
//                                         log weighted error  33.02 
//                               significant figures required  33.31 
//                                    decimal places required  33.65 
//
//-----------------------------------------------------------------------------//

double BesselI0 (double x)
{
    static double bi0cs[18] = { 
	-0.07660547252839144951081894976243285,
	1.927337953993808269952408750881196,
	0.2282644586920301338937029292330415,
	0.01304891466707290428079334210691888,
	4.344270900816487451378682681026107e-4,
	9.422657686001934663923171744118766e-6,
	1.434006289510691079962091878179957e-7,
	1.613849069661749069915419719994611e-9,
	1.396650044535669699495092708142522e-11,
	9.579451725505445344627523171893333e-14,
	5.333981859862502131015107744e-16,
	2.458716088437470774696785919999999e-18,
	9.535680890248770026944341333333333e-21,
	3.154382039721427336789333333333333e-23,
	9.004564101094637431466666666666666e-26,
	2.240647369123670016e-28,
	4.903034603242837333333333333333333e-31,
	9.508172606122666666666666666666666e-34 
    };

    const double xsml = sqrt(0.5 * DBL_EPSILON * 4.5),
	         xmax = log(DBL_MAX);

    double ret_val;

    static int nti0, first = 1;
    if (first) {
	nti0 = initds(bi0cs, 18, 0.5 * DBL_EPSILON * 0.1);
	first = 0;
    }

    double y = fabs(x);
    if (y > 3.0) goto L20;

    ret_val = 1.0;
    if (y > xsml) ret_val = dcsevl(y * y / 4.5 - 1.0, bi0cs, nti0) + 2.75;
    return ret_val;

  L20:
    if (y > xmax) {
	Matpack.Error(Mat::Overflow, "%s: %s", "BesselI0", 
		      "abs(x) so big I0(x) overflows");
	return NAN;
    }

    return exp(y) * BesselExpI0(x);
}

//-----------------------------------------------------------------------------//
