/*-----------------------------------------------------------------------------*\
| Matpack special functions - Dawson(x) Dawson's integral              ddaws.cc |
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
// double Dawson (double x);
//
// Dawson(x) evaluates Dawson's integral for a double precision real argument x. 
//
//                       2  / x   2 
//                     -x   |    t 
//             F(x) = e     |   e    dt 
//                          | 
//                          / 0 
//
// The main computation uses rational Chebyshev approximations 
// published in Math. Comp. 24, 171-178 (1970) by Cody, Paciorek 
// and Thacher.  This transportable program is patterned after the 
// machine-dependent FUNPACK program DDAW(X), but cannot match that 
// version for efficiency or accuracy.  This version uses rational 
// approximations that are theoretically accurate to about 19 
// significant decimal digits.  The accuracy achieved depends on the 
// arithmetic system, the compiler, the intrinsic functions, and 
// proper selection of the machine-dependent constants. 
//
// Underflow: The program returns 0.0 for |X| > XMAX. 
//
// This is a translation from the Fortran version of a SPECFUN (NETLIB)
// REVISION June 15, 1988 , originally written by W. J. Cody, (Mathematics 
// and Computer Science Division, Argonne National Laboratory, Argonne, 
// IL 60439) to C++. 
//
// Explanation of machine-dependent constants:
//
//   XINF   = largest positive machine number 
//            (ANSI C: DBL_MAX)
//   XMIN   = the smallest positive machine number. 
//            (ANSI C: DBL_MIN)
//   EPS    = smallest positive number such that 1+eps > 1. 
//            Approximately  beta**(-p), where beta is the machine 
//            radix and p is the number of significant base-beta 
//            digits in a floating-point number. 
//            (ANSI C: 0.5*DBL_EPSILON)
//   XMAX   = absolute argument beyond which DAW(X) underflows. 
//            XMAX = min(0.5/xmin, xinf). 
//   XSMALL = absolute argument below DAW(X)  may be represented 
//            by X.  We recommend XSMALL = sqrt(eps). 
//   XLARGE = argument beyond which DAW(X) may be represented by 
//            1/(2x).  We recommend XLARGE = 1/sqrt(eps). 
//
// Approximate values for some important machines are 
//
//                        beta  p     eps     xmin       xinf 
//
//  CDC 7600      (S.P.)    2  48  7.11E-15  3.14E-294  1.26E+322 
//  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465 
//  IEEE (IBM/XT, 
//    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38 
//  IEEE (IBM/XT, 
//    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308 
//  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75 
//  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38 
//                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38 
//   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307 
//
//                         XSMALL     XLARGE     XMAX 
//
//  CDC 7600      (S.P.)  5.96E-08   1.68E+07  1.59E+293 
//  CRAY-1        (S.P.)  5.96E-08   1.68E+07  5.65E+2465 
//  IEEE (IBM/XT, 
//    SUN, etc.)  (S.P.)  2.44E-04   4.10E+03  4.25E+37 
//  IEEE (IBM/XT, 
//    SUN, etc.)  (D.P.)  1.05E-08   9.49E+07  2.24E+307 
//  IBM 3033      (D.P.)  3.73D-09   2.68E+08  7.23E+75 
//  VAX 11/780    (S.P.)  2.44E-04   4.10E+03  1.70E+38 
//                (D.P.)  3.73E-09   2.68E+08  1.70E+38 
//   (G Format)   (D.P.)  1.05E-08   9.49E+07  8.98E+307 
//
//
// These values are not neccessary in ANSI C++ because they are calculated
// at compile time from values given in the standard libraries! 
// 
//-----------------------------------------------------------------------------//


double Dawson (double x)
{
    static double zero = 0.0;

    //----------------------------------------------------------------------
    // Coefficients for R(9,9) approximation for  |x| < 2.5
    //----------------------------------------------------------------------
    static double p1[10] = { 
	-2.6902039878870478241e-12, 4.18572065374337710778e-10,
	-1.34848304455939419963e-8, 9.28264872583444852976e-7,
	-1.23877783329049120592e-5, 4.07205792429155826266e-4,
	-0.00284388121441008500446, 0.0470139022887204722217,
	-0.138868086253931995101,   1.00000000000000000004 
    };
    
    static double q1[10] = { 
	1.71257170854690554214e-10, 1.19266846372297253797e-8,
	4.32287827678631772231e-7,  1.03867633767414421898e-5,
	1.7891096528424624934e-4,   0.00226061077235076703171,
	0.0207422774641447644725,   0.132212955897210128811,
	0.527798580412734677256,    1.0 };

    //----------------------------------------------------------------------
    // Coefficients for R(9,9) approximation in J-fraction form
    // for  x in [2.5, 3.5)
    //----------------------------------------------------------------------
    
    static double p2[10] = { 
	 -1.7095380470085549493,  -37.9258977271042880786,
	 26.1935631268825992835,   12.5808703738951251885,
	-22.7571829525075891337,    4.56604250725163310122,
	 -7.3308008989640287075,   46.5842087940015295573,
	-17.3717177843672791149,    0.500260183622027967838 
    };

    static double q2[9] = { 
	 1.82180093313514478378, 1100.67081034515532891,
	-7.08465686676573000364,  453.642111102577727153,
	40.6209742218935689922,   302.890110610122663923,
       170.641269745236227356,    951.190923960381458747,
	 0.206522691539642105009 
    };

    //----------------------------------------------------------------------
    // Coefficients for R(9,9) approximation in J-fraction form
    // for  x in [3.5, 5.0]
    //----------------------------------------------------------------------

    static double p3[10] = { 
	-4.55169503255094815112,  -18.6647123338493852582,
	-7.36315669126830526754,  -66.8407240337696756838,
	48.450726508149145213,     26.9790586735467649969,
       -33.5044149820592449072,     7.50964459838919612289,
	-1.48432341823343965307,    0.499999810924858824981 
    };

    static double q3[9] = { 
	44.7820908025971749852,    99.8607198039452081913,
	14.0238373126149385228,  3488.17758822286353588,
	-9.18871385293215873406, 1240.18500009917163023,
       -68.8024952504512254535,    -2.3125157538514514307,
	 0.250041492369922381761 
    };

    //----------------------------------------------------------------------
    // Coefficients for R(9,9) approximation in J-fraction form
    // for  |x| > 5.0
    //----------------------------------------------------------------------
 
    static double p4[10] = { 
	-8.11753647558432685797,  -38.404388247745445343,
       -22.3787669028751886675,   -28.8301992467056105854,
	-5.99085540418222002197,  -11.3867365736066102577,
	-6.5282872752698074159,    -4.50002293000355585708,
	-2.50000000088955834952,    0.5000000000000004884 
    };

    static double q4[9] = { 
       269.382300417238816428,     50.4198958742465752861,
        61.1539671480115846173,   208.210246935564547889,
        19.7325365692316183531,   -12.2097010558934838708,
        -6.99732735041547247161,   -2.49999970104184464568,
         0.749999999999027092188 
    };

    const double half   = 0.5,
	         one    = 1.0,
	         six25  = 6.25,
	         one225 = 12.25,
	         two5   = 25.0,
	         // calculate limits (at compile time)
	         xsmall = sqrt(0.5*DBL_EPSILON),
	         xlarge = 1.0/sqrt(0.5*DBL_EPSILON),
	         xmax   = MpMin(0.5/DBL_MIN,DBL_MAX);

    double ret_val, frac, sump, sumq, ax, y, w2;

    ax = fabs(x);
    if (ax > xlarge) 
	if (ax <= xmax) 
	    ret_val = half / x;
	else {
	    Matpack.Warning(Mat::Underflow, "%s: %s", "Dawson",
			    "abs(x) so large Dawson(x) underflows");
	    ret_val = zero;
	}
    else if (ax < xsmall) 
	ret_val = x;
    else {
	y = x * x;
	if (y < six25) {
            // --------------------------------------------------------------
            //  abs(x) < 2.5 
            // --------------------------------------------------------------
	    sump = p1[0];
	    sumq = q1[0];
	    for (int i = 1; i < 10; ++i) {
		sump = sump * y + p1[i];
		sumq = sumq * y + q1[i];
	    }
	    ret_val = x * sump / sumq;
	} else if (y < one225) {
            // --------------------------------------------------------------
            //  2.5 <= abs(x) < 3.5 
            // --------------------------------------------------------------
	    frac = zero;
	    for (int i = 0; i < 9; ++i) 
		frac = q2[i] / (p2[i] + y + frac);
	    ret_val = (p2[9] + frac) / x;
	} else if (y < two5) {
            // --------------------------------------------------------------
            //  3.5 <= abs(x) < 5.0 
            // --------------------------------------------------------------
	    frac = zero;
	    for (int i = 0; i < 9; ++i)
		frac = q3[i] / (p3[i] + y + frac);
	    ret_val = (p3[9] + frac) / x;
	} else {
            // --------------------------------------------------------------
            //  5.0 <= abs(x) .<= xlarge 
            // --------------------------------------------------------------
	    w2 = one / x / x;
	    frac = zero;
	    for (int i = 0; i < 9; ++i)
		frac = q4[i] / (p4[i] + y + frac);
	    frac = p4[9] + frac;
	    ret_val = (half + half * w2 * frac) / x;
	}
    }
    return ret_val;
}

//-----------------------------------------------------------------------------//
