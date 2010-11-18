/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselJ(n,x) 			      dbesjn.cc |
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
// double BesselJ (int n, double x)
//
// BesselJ(n,x) calculates the double precision Bessel function of the 
// first kind of order n for double precision argument x. 
//
// Miller's downward recurrence algorithm is used.
// Implemented by B. M. Gammel, last revision 13.03.1996
// 
//-----------------------------------------------------------------------------//

double BesselJ (int n, double x)
{
    const double accuracy  = 40.0,
	         bignum    = 1.0e10,
	         invbignum = 1.0e-10;

    if (n < 0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "BesselJ",
		      "index n less than 0");
	return NAN;
    } else if (n == 0)
	return BesselJ0(x);
    else if (n == 1)
	return BesselJ1(x);
    else {    
	double ret_val;
	double ax = fabs(x);
	if (ax == 0.0)
	    return 0.0;
	else if (ax > (double) n) {
	    double bj,bjm,bjp;
	    double tox = 2.0 / ax;
	    bjm = BesselJ0(ax);
	    bj  = BesselJ1(ax);
	    for (int j = 1; j < n; j++) {
		bjp = j * tox * bj - bjm;
		bjm = bj;
		bj = bjp;
	    }
	    ret_val = bj;
	} else {
	    double bj,bjm,bjp,sum;
	    double tox = 2.0 / ax;
	    int m = 2 * ((n + (int) sqrt(accuracy*n)) / 2);
	    int jsum = 0;
	    bjp = ret_val = sum = 0.0;
	    bj = 1.0;
	    for (int j = m; j > 0; j--) {
		bjm = j * tox * bj - bjp;
		bjp = bj;
		bj = bjm;
		if (fabs(bj) > bignum) {
		    bj *= invbignum;
		    bjp *= invbignum;
		    ret_val *= invbignum;
		    sum *= invbignum;
		}
		if (jsum) sum += bj;
		jsum = !jsum;
		if (j == n) ret_val = bjp;
	    }
	    sum = 2.0 * sum - bj;
	    ret_val /= sum;
	}
	return  (x < 0.0 && n % 2 == 1) ? -ret_val : ret_val;
    } 
}

//-----------------------------------------------------------------------------//
