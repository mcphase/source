/*-----------------------------------------------------------------------------*\
| Matpack special functions - Gamma(x) complete gamma function        dgamma.cc |
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
// double Gamma (double x);
//
// Gamma(x) calculates the double precision complete Gamma function 
// for double precision argument x. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C7A, REVISION 920618, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for GAM        on the interval  0.          to  1.00000E+00 
//                                        with weighted error   5.79E-32 
//                                         log weighted error  31.24 
//                               significant figures required  30.00 
//                                    decimal places required  32.05 
//
//-----------------------------------------------------------------------------//

double Gamma (double x)
{
    static double gamcs[42] = { 
	 0.008571195590989331421920062399942,
	 0.004415381324841006757191315771652,
	 0.05685043681599363378632664588789,
	-.004219835396418560501012500186624,
	 0.001326808181212460220584006796352,
	-1.893024529798880432523947023886e-4,
	 3.606925327441245256578082217225e-5,
	-6.056761904460864218485548290365e-6,
	 1.055829546302283344731823509093e-6,
	-1.811967365542384048291855891166e-7,
	 3.117724964715322277790254593169e-8,
	-5.354219639019687140874081024347e-9,
	 9.19327551985958894688778682594e-10,
	-1.577941280288339761767423273953e-10,
	 2.707980622934954543266540433089e-11,
	-4.646818653825730144081661058933e-12,
	 7.973350192007419656460767175359e-13,
	-1.368078209830916025799499172309e-13,
	 2.347319486563800657233471771688e-14,
	-4.027432614949066932766570534699e-15,
	 6.910051747372100912138336975257e-16,
	-1.185584500221992907052387126192e-16,
	 2.034148542496373955201026051932e-17,
	-3.490054341717405849274012949108e-18,
	 5.987993856485305567135051066026e-19,
	-1.027378057872228074490069778431e-19,
	 1.762702816060529824942759660748e-20,
	-3.024320653735306260958772112042e-21,
	 5.188914660218397839717833550506e-22,
	-8.902770842456576692449251601066e-23,
	 1.527474068493342602274596891306e-23,
	-2.620731256187362900257328332799e-24,
	 4.496464047830538670331046570666e-25,
	-7.714712731336877911703901525333e-26,
	 1.323635453126044036486572714666e-26,
	-2.270999412942928816702313813333e-27,
	 3.896418998003991449320816639999e-28,
	-6.685198115125953327792127999999e-29,
	 1.146998663140024384347613866666e-29,
	-1.967938586345134677295103999999e-30,
	 3.376448816585338090334890666666e-31,
	-5.793070335782135784625493333333e-32 
    };
    
    const double pi     = 3.1415926535897932384626433832795,
	         sq2pil = 0.91893853320467274178032973640562,
	         dxrel  = sqrt(DBL_EPSILON);

    double ret_val;
    static int ngam, first = 1;
    static double xmin, xmax;

    if (first) {
	ngam = initds(gamcs, 42, 0.5 * DBL_EPSILON * 0.1);
	dgamlm(xmin, xmax);
	first = 0;
    }

    int i,n;

    double y = fabs(x);
    if (y > 10.0) goto L50;    

    // compute gamma(x) for -xbnd <= x <= xbnd.  Reduce interval and find 
    // gamma(1+y) for 0.0 <= y < 1.0 first of all. 

    n = (int)x;
    if (x < 0.0) --n;
    y = x - n;
    --n;

    ret_val = dcsevl(y * 2.0 - 1.0, gamcs, ngam) + 0.9375;
    if (n == 0) return ret_val;
    
    if (n > 0) goto L30;

    // compute gamma(x) for x < 1.0 

    n = -n;
    if (x == 0.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "Gamma", 
		      "x is 0");
	return NAN;
    }

    if (x < 0.0 && x + n - 2 == 0.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "Gamma", 
		      "x is a negative integer");
	return NAN;
    }

    if (x < -0.5 && fabs((x - Dint(x - 0.5)) / x) < dxrel) {
	Matpack.Warning(Mat::PartialPrecisionLoss, "%s: %s", "Gamma",
			"answer less than half precision because x"
			" too near negative integer");
    }

    for (i = 1; i <= n; ++i) 
	ret_val /= x + i - 1;
    return ret_val;

    // gamma(x) for x >= 2.0 and x <= 10.0 

  L30:
    for (i = 1; i <= n; ++i) 
	ret_val = (y + i) * ret_val;
    return ret_val;
    
    // gamma(x) for abs(x) > 10.0.  recall y = abs(x).

  L50:
    if (x > xmax) {
	Matpack.Error(Mat::Overflow, "%s: %s", "Gamma",
		      "x so big Gamma(x) overflows");	
	return NAN;
    }

    ret_val = 0.0;

    if (x < xmin) 
	Matpack.Warning(Mat::Underflow, "%s: %s", "Gamma",
			"x so small Gamma(x) underflows");
    
    if (x < xmin) return ret_val;
    
    ret_val = exp((y - 0.5) * log(y) - y + sq2pil + d9lgmc(y));
    if (x > 0.0) return ret_val;
    
    if (fabs((x - Dint(x - 0.5)) / x) < dxrel) 
	Matpack.Warning(Mat::PartialPrecisionLoss, "%s: %s", "Gamma",
			"answer less than half precision because x"
			" too near negative integer");

    double sinpiy = sin(pi * y);
    if (sinpiy == 0.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "Gamma",
		      "x is a negative integer");	
	return NAN;
    }

    return -pi / (y * sinpiy * ret_val);
}

//-----------------------------------------------------------------------------//
