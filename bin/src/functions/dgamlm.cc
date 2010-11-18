/*-----------------------------------------------------------------------------*\
| Matpack special functions - dgamlm(xmin,xmax)                       dgamlm.cc |
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
// void dgamlm (double& xmin, double& xmax);
//
// Calculate the minimum and maximum legal bounds for x in Gamma(x). 
// xmin and xmax are not the only bounds, but they are the only non- 
// trivial ones to calculate. 
//
// Output Arguments:
//
// 	xmin	double precision minimum legal value of x in gamma(x).  Any 
//       	smaller value of x might result in underflow. 
// 	xmax	double precision maximum legal value of x in gamma(x).  Any 
//      	larger value of x might cause overflow. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C7A, R2 REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
//-----------------------------------------------------------------------------//

void dgamlm (double& xmin, double& xmax)
{
    const double alnsml = log(DBL_MIN),
	         alnbig = log(DBL_MAX);

    double xold,xln;
    int i;

    xmin = -alnsml;
    for (i = 1; i <= 10; ++i) {
	xold = xmin;
	xln = log(xmin);
	xmin -= xmin * ((xmin + 0.5) * xln - xmin - 0.2258 + alnsml) 
	             / (xmin * xln + 0.5);
	if (fabs(xmin - xold) < 0.005) goto L20;
    }

    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "dgamlm"
		  "unable to find xmin");
    xmin = xmax = NAN;
    return;

  L20:
    xmin = -xmin + 0.01;
    xmax = alnbig;
    for (i = 1; i <= 10; ++i) {
	xold = xmax;
	xln = log(xmax);
	xmax -= xmax * ((xmax - 0.5) * xln - xmax + 0.9189 - alnbig) 
	             / (xmax * xln - 0.5);
	if (fabs(xmax - xold) < 0.005) goto L40;
    }

    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "dgamlm"
		  "unable to find xmax");
    xmin = xmax = NAN;
    return;

  L40:
    xmax += -0.01;
    xmin = MpMax(xmin,-xmax + 1.0);
}

//-----------------------------------------------------------------------------//
