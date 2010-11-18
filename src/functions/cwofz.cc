/*-----------------------------------------------------------------------------*\
| Matpack special functions - Faddeeva(z)                              cwofz.cc |
|                                                                               |
| Matpack Library Release 1.0                                                   |
| Copyright (C) 1991-1997 by Berndt M. Gammel                                   |
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
// complex<double> Faddeeva (const complex<double>& z)
//
// Given a complex number z = (x,y), this subroutine computes 
// the value of the Faddeeva function w(z) = exp(-z^2)*erfc(-i*z), 
// where erfc is the complex complementary error function and i means sqrt(-1). 
//
// The accuracy of the algorithm for z in the 1st and 2nd quadrant 
// is 14 significant digits; in the 3rd and 4th it is 13 significant 
// digits outside a circular region with radius 0.126 around a zero 
// of the function. 
//
// All real variables in the program are double precision. 
// The parameter M_2_SQRTPI equals 2/sqrt(pi).  
//
// The routine is not underflow-protected but any variable can be 
// put to 0 upon underflow; 
//
// The routine is overflow-protected: Matpack::Error() is called.
//
// References:
//
// (1) G.P.M. Poppe, C.M.J. Wijers; More Efficient Computation of 
//     the Complex Error-Function, ACM Trans. Math. Software,
//     Vol. 16, no. 1, pp. 47.
// (2) Algorithm 680, collected algorithms from ACM.
// 
// The Fortran source code was translated to C++ by B.M. Gammel
// and added to the MatPack library, 1992.
//
// Last change: B. M. Gammel, 18.03.1996 error handling
// 
//-----------------------------------------------------------------------------//

complex<double> Faddeeva (const complex<double>& z)
{
    // The maximum value of rmaxreal equals the root of the largest number 
    // rmax which can still be implemented on the computer in double precision
    // floating-point arithmetic
    const double rmaxreal = sqrt(DBL_MAX);

    // rmaxexp  = ln(rmax) - ln(2)
    const double rmaxexp  = log(DBL_MAX) - log(2.0);

    // the largest possible argument of a double precision goniometric function
    const double rmaxgoni = 1.0 / DBL_EPSILON;

    double xabs, yabs, daux, qrho, xaux, xsum, ysum, c, h, u, v,
           x, y, xabsq, xquad, yquad, h2, u1, v1, u2, v2, w1,
           rx, ry, sx, sy, tx, ty, qlambda, xi, yi;
    int    i, j, n, nu, np1, kapn, a, b;

    xi = real(z);
    yi = imag(z);
    xabs = fabs(xi);
    yabs = fabs(yi);
    x = xabs / 6.3;
    y = yabs / 4.4;

    // the following statement protects qrho = (x^2 + y^2) against overflow
    if ((xabs > rmaxreal) || (yabs > rmaxreal)) {
        Matpack.Error(Mat::Overflow, "%s: %s", "Faddeeva",
		      "absolute value of argument so large w(z) overflows");
	return NAN;
    }

    qrho = x * x + y * y;
    xabsq = xabs * xabs;
    xquad = xabsq - yabs * yabs;
    yquad = xabs * 2 * yabs;
    
    a = qrho < 0.085264;

    if (a) {
	
        // If (qrho < 0.085264) then the Faddeeva-function is evaluated 
        // using a power-series (Abramowitz/Stegun, equation (7.1.5), p.297). 
        // n is the minimum number of terms needed to obtain the required 
        // accuracy.
	
	qrho = (1 - y * 0.85) * sqrt(qrho);
	n = Nint(qrho * 72 + 6);
	j = (n << 1) + 1;
	xsum = 1.0 / j;
	ysum = 0.0;
	for (i = n; i >= 1; --i) {
	    j -= 2;
	    xaux = (xsum * xquad - ysum * yquad) / i;
	    ysum = (xsum * yquad + ysum * xquad) / i;
	    xsum = xaux + 1.0 / j;
	}
	u1 = (xsum * yabs + ysum * xabs) * -M_2_SQRTPI + 1.0;
	v1 = (xsum * xabs - ysum * yabs) * M_2_SQRTPI;
	daux = exp(-xquad);
	u2 = daux * cos(yquad);
	v2 = -daux * sin(yquad);
	
	u = u1 * u2 - v1 * v2;
	v = u1 * v2 + v1 * u2;
	
    } else {
	
        //  If (qrho > 1.0) then w(z) is evaluated using the Laplace continued 
	//  fraction.  nu is the minimum number of terms needed to obtain the
	//  required accuracy. 
        //  if ((qrho > 0.085264) && (qrho < 1.0)) then w(z) is evaluated
        //  by a truncated Taylor expansion, where the Laplace continued
	//  fraction is used to calculate the derivatives of w(z). 
        //  kapn is the minimum number of terms in the Taylor expansion needed
        //  to obtain the required accuracy. 
        //  nu is the minimum number of terms of the continued fraction needed
        //  to calculate the derivatives with the required accuracy. 
	
	if (qrho > 1.0) {
	    h = 0.0;
	    kapn = 0;
	    qrho = sqrt(qrho);
	    nu = (int) (1442 / (qrho * 26 + 77) + 3);
	} else {
	    qrho = (1 - y) * sqrt(1 - qrho);
	    h = qrho * 1.88;
	    h2 = h * 2;
	    kapn = Nint(qrho * 34 + 7);
	    nu   = Nint(qrho * 26 + 16);
	}
	
	b = h > 0.0;
	
	if (b) qlambda = pow(h2, (double) kapn);
	
	rx = ry = sx = sy = 0.0;
	for (n = nu; n >= 0; --n) {
	    np1 = n + 1;
	    tx = yabs + h + np1 * rx;
	    ty = xabs - np1 * ry;
	    c = 0.5 / (tx * tx + ty * ty);
	    rx = c * tx;
	    ry = c * ty;
	    if (b && (n <= kapn)) {
		tx = qlambda + sx;
		sx = rx * tx - ry * sy;
		sy = ry * tx + rx * sy;
		qlambda /= h2;
	    }
	}
	
	if (h == 0.0) {
	    u = rx * M_2_SQRTPI;
	    v = ry * M_2_SQRTPI;
	} else {
	    u = sx * M_2_SQRTPI;
	    v = sy * M_2_SQRTPI;
	}
	
	if (yabs == 0.0) 
	    u = exp(-(xabs * xabs));
    }
    
    //  evaluation of w(z) in the other quadrants 

    if (yi < 0.0) {

	if (a) {
	    u2 *= 2;
	    v2 *= 2;
	} else {
	    xquad = -xquad;

            // the following statement protects 2*exp(-z**2) against overflow
            if ((yquad > rmaxgoni) || (xquad > rmaxexp)) {
		Matpack.Error(Mat::Overflow, "%s: %s", "Faddeeva",
			      "absolute value of argument so large w(z) overflows");
            	return NAN;
	    }

	    w1 = exp(xquad) * 2;
	    u2 =  w1 * cos(yquad);
	    v2 = -w1 * sin(yquad);
	}

	u = u2 - u;
	v = v2 - v;
	if (xi > 0.0) v = -v;

    } else if (xi < 0.0) 
        v = -v;
   
    return complex<double>(u,v);
}

//-----------------------------------------------------------------------------//
