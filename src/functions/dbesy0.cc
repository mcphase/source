/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselY0(x) 			      dbesy0.cc |
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
// double BesselY0 (double x)
//
// BesselY0(x) calculates the double precision Bessel function of the 
// second kind of order zero for double precision argument X. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10A1, REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BY0        on the interval  0.          to  1.60000E+01 
//                                        with weighted error   8.14E-32 
//                                         log weighted error  31.09 
//                               significant figures required  30.31 
//                                    decimal places required  31.73 
//
//-----------------------------------------------------------------------------//

double BesselY0 (double x)
{
    static double by0cs[19] = { 
	-0.01127783939286557321793980546028,
	-0.1283452375604203460480884531838,
	-0.1043788479979424936581762276618,
	 0.02366274918396969540924159264613,
	-0.002090391647700486239196223950342,
	 1.039754539390572520999246576381e-4,
	-3.369747162423972096718775345037e-6,
	 7.729384267670667158521367216371e-8,
	-1.324976772664259591443476068964e-9,
	 1.764823261540452792100389363158e-11,
	-1.881055071580196200602823012069e-13,
	 1.641865485366149502792237185749e-15,
	-1.19565943860460608574599100672e-17,
	 7.377296297440185842494112426666e-20,
	-3.906843476710437330740906666666e-22,
	 1.79550366443615794982912e-24,
	-7.229627125448010478933333333333e-27,
	 2.571727931635168597333333333333e-29,
	-8.141268814163694933333333333333e-32 
    };

    const double twodpi = 0.636619772367581343075535053490057,
	         xsml   = sqrt(0.5 * DBL_EPSILON * 4.0);

    static int nty0, first = 1;
    if (first) {
	nty0 = initds(by0cs, 19, 0.5 * DBL_EPSILON * 0.1);
	first = 0;
    }
    
    if (x <= 0.0) {	
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "BesselY0",
		      "x is zero or negative");
	return NAN;
    }

    if (x > 4.0) {
	double ampl, theta;
	d9b0mp(x, ampl, theta);
	return ampl * sin(theta);
    } else {
	double y = 0.0;
	if (x > xsml) y = x * x;
	return twodpi * log(x*0.5) * BesselJ0(x) + 0.375 
	                  + dcsevl(y * 0.125-1.0, by0cs, nty0);
    }
} 

//-----------------------------------------------------------------------------//
