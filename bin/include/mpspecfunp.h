/*-----------------------------------------------------------------------------*\
| Matpack special functions private include file                   mpspecfunp.h |
|                                                                               |
| MatPack Library Release 1.0                                                   |
| Copyright (C) 1991-1996 by Berndt M. Gammel                                   |
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

#ifndef _MPSPECFUNP_H_
#define _MPSPECFUNP_H_

// SOME NOTES WHEN CONVERTING FORTRAN ROUTINES TO C++
// r1mach(n) returns:
//      FORTRAN                             C++
// d(1)=2.225073858507201e-308   DBL_MIN = 2.225073858507201e-308
// d(2)=1.797693134862316e+308   DBL_MAX = 1.797693134862316e+308
// d(3)=1.110223024625157e-16    0.5 * DBL_EPSILON
// d(4)=2.220446049250313e-16    DBL_EPSILON= 2.220446049250313e-16
// d(5)=0.3010299956639812       log10(basis)
// double d_int(double*) rounds towards zero

#ifndef _FLOAT_H___
#define DBL_EPSILON  2.220446049250313e-16
#define   DBL_MIN  2.225073858507201e-258
#define   DBL_MAX  1.797693134862316e+258
#endif

//#include "common.h"
#include "vector.h"
//#include <complex.h>
#include<cmath>
#include "mpspecfun.h"	// include also the public prototypes

//-----------------------------------------------------------------------------//
// Since this is the core routine for most special functions it is inlined here!
//-----------------------------------------------------------------------------//
//
// double dcsevl (double x, double* cs, int n);
//
// Evaluate the n-term Chebyshev series cs at x.  Adapted from 
// a method presented in the paper by Broucke referenced below. 
//
// Input Arguments:
//
//	x    value at which the series is to be evaluated. 
//	cs   array of n terms of a Chebyshev series. In evaluating 
//	     cs, only half the first coefficient is summed. 
//	n    number of terms in array cs. 
//
// References:
//
//	R. Broucke, Ten subroutines for the manipulation of Chebyshev series, 
//	Algorithm 446, Communications of the A.C.M. 16, (1973) pp. 254-256. 
//
//	L. Fox and I. B. Parker, Chebyshev Polynomials in 
//      Numerical Analysis, Oxford University Press, 1968,  page 56. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C3A2, REVISION  920501, originally written by Fullerton W., (LANL) 
// to C++.
// 
//-----------------------------------------------------------------------------//

inline double dcsevl (double x, const double* cs, int n)
{
#ifdef DEBUG
    if (n < 1) {   
	Matpack.Error("dcsevl: number of terms <= 0");
	return NAN;
    }
    if (n > 1000){ 
	Matpack.Error("dcsevl: number of terms > 1000");
	return NAN;
    }
    if (fabs(x) > DBL_EPSILON + 1.0) 
	Matpack.Warning("%s: %s", "dcsevl",
			"x outside the interval (-1,+1)");
#endif
    double b0 = 0.0, b1 = 0.0, b2, twox = x * 2;
    for (int i = 1; i <= n; i++) {
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + cs[n - i];
    }

    return (b0 - b2) * 0.5;
} 

//-----------------------------------------------------------------------------//

// private auxilliary functions

// real
int 	initds		(double *os, int nos, double eta);
void	dgamlm		(double& xmin, double& xmax);
double	d9lgmc		(double x);
void	d9b0mp		(double x, double& ampl, double& theta);
void	d9b1mp		(double x, double& ampl, double& theta);
void	d9aimp		(double x, double& ampl, double& theta);

// complex
complex<double> c9lgmc  (const complex<double> &z);

//-----------------------------------------------------------------------------//

#endif

