/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselJ1(x) 			      dbesj1.cc |
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
// double BesselJ1 (double x)
//
// BesselJ1(x) calculates the double precision Bessel function of the 
// first kind of order one for double precision argument x. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10A1, REVISION 910401 originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BJ1        on the interval  0.          to  1.60000E+01 
//                                        with weighted error   1.16E-33 
//                                         log weighted error  32.93 
//                               significant figures required  32.36 
//                                    decimal places required  33.57 
//
//-----------------------------------------------------------------------------//

double BesselJ1 (double x)
{
    static double bj1cs[19] = { 
	-0.117261415133327865606240574524003,
	-0.253615218307906395623030884554698,
	 0.0501270809844695685053656363203743,
	-0.00463151480962508191842619728789772,
	 2.47996229415914024539124064592364e-4,
	-8.67894868627882584521246435176416e-6,
	 2.14293917143793691502766250991292e-7,
	-3.93609307918317979229322764073061e-9,
	 5.59118231794688004018248059864032e-11,
	-6.3276164046613930247769527401488e-13,
	 5.84099161085724700326945563268266e-15,
	-4.48253381870125819039135059199999e-17,
	 2.90538449262502466306018688e-19,
	-1.61173219784144165412118186666666e-21,
	 7.73947881939274637298346666666666e-24,
	-3.24869378211199841143466666666666e-26,
	 1.2022376772274102272e-28,
	-3.95201221265134933333333333333333e-31,
	 1.16167808226645333333333333333333e-33 
    };

    const double xsml = sqrt(0.5 * DBL_EPSILON * 8.0),
	         xmin = DBL_MIN * 2.0;

    double ret_val;

    static int ntj1, first = 1;
    if (first) {
	ntj1 = initds(bj1cs, 19, 0.5 * DBL_EPSILON * 0.1);
	first = 1;
    }

    double y = fabs(x);
    if (y > 4.0) goto L20;
    
    if (y == 0.0) return 0.0;
    
    if (y <= xmin) 
	Matpack.Warning(Mat::Underflow, "%s: %s", "BesselJ1",
			"abs(x) so small J1(x) underflows");
    
    if (y > xmin) 
	ret_val = x * 0.5;
    
    if (y > xsml)
	ret_val = x * (dcsevl(y * 0.125 * y - 1.0, bj1cs, ntj1) + 0.25);
    
    return ret_val;
    
  L20:
    double ampl, theta;
    d9b1mp(y, ampl, theta);
    return CopySign(ampl, x) * cos(theta);
}

//-----------------------------------------------------------------------------//
