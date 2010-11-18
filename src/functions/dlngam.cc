/*-----------------------------------------------------------------------------*\
| Matpack special functions - LnGamma(x) Log Gamma Function           dlngam.cc |
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
// double LnGamma (double x);
//
// Calculates the double precision logarithm of the 
// absolute value of the Gamma function for double precision 
// argument x. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C7A, REVISION 900727, originally written by Fullerton W.,(LANL)
// to C++.
//
//-----------------------------------------------------------------------------//


double LnGamma (double x)
{
    const double sq2pil = 0.91893853320467274178032973640562,
	         sqpi2l = 0.225791352644727432363097614947441,
	         pi     = 3.1415926535897932384626433832795,
	         temp   = 1.0 / log(DBL_MAX),
	         xmax   = temp * DBL_MAX,
	         dxrel  = sqrt(DBL_EPSILON);

    double y = fabs(x);
    if (y > 10.0) goto L20;

    // log (abs (Gamma(x)) ) for abs(x) <= 10.0 

    return log(fabs(Gamma(x)));

  L20:
    if (y > xmax) {
	Matpack.Error(Mat::Overflow, "%s: %s", "LnGamma",
		      "abs(x) so big LnGamma(x) overflows");
	return NAN;
    }

    if (x > 0.0)
	return sq2pil + (x - 0.5) * log(x) - x + d9lgmc(y);
    
    double sinpiy = fabs(sin(pi * y));

    if (sinpiy == 0.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "LnGamma",
		      "x is a negative integer");
	return NAN;
    }

    if (fabs((x - Dint(x - 0.5)) / x) < dxrel) 
	Matpack.Warning(Mat::PartialPrecisionLoss, "%s: %s", "LnGamma",
			"answer less than half precision because x too near "
			"negative integer");
    
    return sqpi2l + (x - 0.5) * log(y) - x - log(sinpiy) - d9lgmc(y);
}

//-----------------------------------------------------------------------------//
