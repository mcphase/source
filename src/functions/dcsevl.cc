/*-----------------------------------------------------------------------------*\
| Matpack special functions - dcsevl(x,cs,n)                          dcsevl.cc |
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

double dcsevl (double x, const double* cs, int n)
{
    if (n < 1) {   
	Matpack.Error("dcsevl: number of terms <= 0");
	return NAN;
    }
    if (n > 1000){ 
	Matpack.Error("dcsevl: number of terms > 1000");
	return NAN;
    }
    if (fabs(x) > DBL_EPSILON + 1.0) 
	Matpack.Warning(Mat::DomainError, "%s: %s", "dcsevl",
			"x outside the interval (-1,+1)");

    double b0 = 0.0, b1 = 0.0, b2, twox = x * 2.0;
    for (int i = 1; i <= n; i++) {
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + cs[n - i];
    }

    return (b0 - b2) * 0.5;
} 

//-----------------------------------------------------------------------------//

