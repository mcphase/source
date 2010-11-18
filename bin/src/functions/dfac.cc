/*-----------------------------------------------------------------------------*\
| Matpack special functions - Fac(n) the factorial                      dfac.cc |
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
// double Fac (int n);
//
// Fac(n) calculates the double precision factorial for the integer argument n.
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C1, REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
//-----------------------------------------------------------------------------//

double Fac (int n)
{
    static double facn[31] = { 
	1.0,
	1.0,
	2.0,
	6.0,
	24.0,
	120.0,
	720.0,
	5040.0,
	40320.0,
	362880.0,
	3628800.0,
	39916800.0,
	479001600.0,
	6227020800.0,
	87178291200.0,
	1.307674368e12,
	2.0922789888e13,
	3.55687428096e14,
	6.402373705728e15,
	1.21645100408832e17,
	2.43290200817664e18,
	5.109094217170944e19,
	1.12400072777760768e21,
	2.585201673888497664e22,
	6.2044840173323943936e23,
	1.5511210043330985984e25,
	4.03291461126605635584e26,
	1.0888869450418352160768e28,
	3.04888344611713860501504e29,
	8.841761993739701954543616e30,
	2.6525285981219105863630848e32 
    };
    
    const double sq2pil = 0.91893853320467274178032973640562;

    static int nmax = 0;
    static double xmin, xmax;
    
    if (nmax == 0) {
	dgamlm(xmin,xmax);
	nmax = (int) (xmax - 1.0);
    }

    if (n < 0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "Fac",
		      "factorial of negative integer undefined");
	return NAN;
    }

    if (n <= 30) return facn[n];
    if (n > nmax) {
	Matpack.Error(Mat::Overflow, "%s: %s", "Fac",
		      "n so big Fac(n) overflows");
	return NAN;
    }

    double x = (double) (n + 1);
    return exp((x - 0.5) * log(x) - x + sq2pil + d9lgmc(x));
}

//-----------------------------------------------------------------------------//
