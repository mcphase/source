/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselI1(x) 			      dbesi1.cc |
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
// double BesselI1 (double x);
//
// BesselI1(x) calculates the double precision modified (hyperbolic) 
// Bessel function of the first kind of order one and double precision 
// argument x. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10B1, REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BI1        on the interval  0.          to  9.00000E+00 
//                                        with weighted error   1.44E-32 
//                                         log weighted error  31.84 
//                               significant figures required  31.45 
//                                    decimal places required  32.46 
//
//-----------------------------------------------------------------------------//

double BesselI1 (double x)
{
    static double bi1cs[17] = { 
	-0.0019717132610998597316138503218149,
	 0.40734887667546480608155393652014,
	 0.034838994299959455866245037783787,
	 0.0015453945563001236038598401058489,
	 4.188852109837778412945883200412e-5,
	 7.6490267648362114741959703966069e-7,
	 1.0042493924741178689179808037238e-8,
	 9.9322077919238106481371298054863e-11,
	 7.6638017918447637275200171681349e-13,
	 4.741418923816739498038809194816e-15,
	 2.4041144040745181799863172032e-17,
	 1.0171505007093713649121100799999e-19,
	 3.6450935657866949458491733333333e-22,
	 1.1205749502562039344810666666666e-24,
	 2.9875441934468088832e-27,
	 6.9732310939194709333333333333333e-30,
	 1.43679482206208e-32 
    };
    
    const double xmin = DBL_MIN * 2.0,
	         xsml = sqrt(0.5*DBL_EPSILON * 4.5),
	         xmax = log(DBL_MAX);

    double ret_val;

    static int nti1, first = 1;
    
    if (first) {
	nti1 = initds(bi1cs, 17, 0.5 * DBL_EPSILON * 0.1 );
	first = 0;
    }
    
    double y = fabs(x);
    if (y > 3.0) goto L20;
    
    ret_val = 0.0;
    if (y == 0.0) return ret_val;
    
    if (y <= xmin) 
	Matpack.Warning(Mat::Underflow, "%s: %s", "BesselI1",
			"abs(x) so small I1(x) underflows");
    
    if (y > xmin) ret_val = x * 0.5;
    if (y > xsml) 
	ret_val = x * (dcsevl(y * y / 4.5 - 1.0, bi1cs, nti1) + 0.875);
    return ret_val;
    
  L20:
    if (y > xmax) { 	
	Matpack.Error(Mat::Overflow, "%s: %s", "BesselI1",
		      "abs(x) so big I1(x) overflows");
	return NAN;
    }

    ret_val = exp(y) * BesselExpI1(x);

    return ret_val;
}

//-----------------------------------------------------------------------------//
