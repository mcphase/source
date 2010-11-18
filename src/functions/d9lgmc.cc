/*-----------------------------------------------------------------------------*\
| Matpack special functions - d9lgmc(x)                               d9lgmc.cc |
|                                                                               |
| Matpack Library Release 1.0                                                   |
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
// double d9lgmc (double x);
//
// Compute the log gamma correction factor for x >= 10.0 so that 
// log (dgamma(x)) = log(sqrt(2*pi)) + (x-0.5)*log(x) - x + d9lgmc(x) 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C7E, REVISION 900720, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for ALGM       on the interval  0.          to  1.00000E-02 
//                                        with weighted error   1.28E-31 
//
//                                         log weighted error  30.89 
//                               significant figures required  29.81 
//                                    decimal places required  31.48 
//
//-----------------------------------------------------------------------------//

double d9lgmc (double x)
{
    static double algmcs[15] = { 
	 0.1666389480451863247205729650822,
	-1.384948176067563840732986059135e-5,
	 9.810825646924729426157171547487e-9,
	-1.809129475572494194263306266719e-11,
	 6.221098041892605227126015543416e-14,
	-3.399615005417721944303330599666e-16,
	 2.683181998482698748957538846666e-18,
	-2.868042435334643284144622399999e-20,
	 3.962837061046434803679306666666e-22,
	-6.831888753985766870111999999999e-24,
	 1.429227355942498147573333333333e-25,
	-3.547598158101070547199999999999e-27,
	 1.025680058010470912e-28,
	-3.401102254316748799999999999999e-30,
	 1.276642195630062933333333333333e-31 
    };

    const double xbig = 1.0 / sqrt(0.5*DBL_EPSILON),
                 xmax = exp( MpMin(log(DBL_MAX / 12.0), -log(DBL_MIN * 12.0)) );

    double ret_val;

    static int nalgm, first = 1;
    if (first) {
	nalgm = initds(algmcs, 15, 0.5*DBL_EPSILON);
	first = 0;
    }

    if (x < 10.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "d9lgmc",
		      "x must be >= 10"); 
	return NAN;
    }

    if (x >= xmax) goto L20;

    ret_val = 1.0 / (x * 12.0);
    if (x < xbig) 
	ret_val = dcsevl(sqr(10.0 / x) * 2.0 - 1.0, algmcs, nalgm) / x;
    return ret_val;
    
  L20:
    ret_val = 0.0;	
    Matpack.Warning(Mat::Underflow, "%s: %s", "d9lgmc",
		    "x so big d9lgmc(x) underflows");
    
    return ret_val;
}

//-----------------------------------------------------------------------------//
