/*-----------------------------------------------------------------------------*\
| Matpack special functions - Airy function                              dai.cc |
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
// double AiryAi (double x);
//
// AiryAi(x) calculates the double precision Airy function for double 
// precision argument x. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10D, REVISION 920618, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for AIF        on the interval -1.00000E+00 to  1.00000E+00 
//                                        with weighted error   8.37E-33 
//                                         log weighted error  32.08 
//                               significant figures required  30.87 
//                                    decimal places required  32.63 
//
// Series for AIG        on the interval -1.00000E+00 to  1.00000E+00 
//                                        with weighted error   7.47E-34 
//                                         log weighted error  33.13 
//                               significant figures required  31.50 
//                                    decimal places required  33.68 
//
//-----------------------------------------------------------------------------//

double AiryAi (double x)
{
    static double aifcs[13] = { 
       -0.037971358496669997496197089469414,
	0.059191888537263638574319728013777,
	9.862928057727997536560389104406e-4,
	6.8488438190765667554854830182412e-6,
	2.5942025962194713019489279081403e-8,
	6.1766127740813750329445749697236e-11,
	1.0092454172466117901429556224601e-13,
	1.2014792511179938141288033225333e-16,
	1.0882945588716991878525295466666e-19,
	7.75137721966848870392384e-23,
	4.4548112037175638391466666666666e-26,
	2.1092845231692343466666666666666e-29,
	8.3701735910741333333333333333333e-33 
    };

    static double aigcs[13] = { 
	0.018152365581161273011556209957864,
	0.021572563166010755534030638819968,
	2.5678356987483249659052428090133e-4,
	1.4265214119792403898829496921721e-6,
	4.5721149200180426070434097558191e-9,
	9.5251708435647098607392278840592e-12,
	1.392563460577139905115042068619e-14,
	1.5070999142762379592306991138666e-17,
	1.2559148312567778822703205333333e-20,
	8.3063073770821340343829333333333e-24,
	4.4657538493718567445333333333333e-27,
	1.9900855034518869333333333333333e-30,
	7.4702885256533333333333333333333e-34 
    };

    const double x3sml = pow(0.5 * DBL_EPSILON, 0.3334),
		 xmaxt = pow(log(DBL_MIN) * -1.5, 0.6667),
                 xmax = xmaxt - xmaxt*log(xmaxt) / (sqrt(xmaxt)*4.0+1.0) - 0.01;

    static int naif, naig, first = 1;
    if (first) {
	naif = initds(aifcs, 13, 0.5 * DBL_EPSILON * 0.1);
	naig = initds(aigcs, 13, 0.5 * DBL_EPSILON * 0.1);
	first = 0;
    }

    double z, theta, xm;

    if (x >= -1.0) goto L20;
    d9aimp(x, xm, theta);
    return xm * cos(theta);

  L20:
    if (x > 1.0) goto L30;
    z = 0.0;
    if (fabs(x) > x3sml) z = x * x * x;
    return dcsevl(z, aifcs, naif) - x * (dcsevl(z, aigcs, naig) + 0.25) + 0.375;

  L30:
    if (x > xmax) goto L40;
    return AiryExpAi(x) * exp(x * -2.0 * sqrt(x) / 3.0);

  L40:	
    Matpack.Warning(Mat::Underflow, "%s: %s", "AiryAi",
		    "x so big Ai(x) underflows");
    return 0.0;
}

//-----------------------------------------------------------------------------//
