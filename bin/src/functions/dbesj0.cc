/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselJ0(x) 			      dbesj0.cc |
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
// double BesselJ0 (double x)
//
// BesselJ0(x) calculates the double precision Bessel function of 
// the first kind of order zero for double precision argument x. 
//
// Series for BJ0        on the interval  0.          to  1.60000E+01 
//                                        with weighted error   4.39E-32 
//                                         log weighted error  31.36 
//                               significant figures required  31.21 
//                                    decimal places required  32.00 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10A1, REVISION 891214 originally written by Fullerton W.,(LANL)
// to C++.
//
//-----------------------------------------------------------------------------//

double BesselJ0 (double x)
{
    static double bj0cs[19] = { 
	 0.10025416196893913701073127264074,
	-0.66522300776440513177678757831124,
	 0.2489837034982813137046046872668,
	-0.033252723170035769653884341503854,
	 0.0023114179304694015462904924117729,
	-9.9112774199508092339048519336549e-5,
	 2.8916708643998808884733903747078e-6,
	-6.1210858663032635057818407481516e-8,
	 9.8386507938567841324768748636415e-10,
	-1.2423551597301765145515897006836e-11,
	 1.2654336302559045797915827210363e-13,
	-1.0619456495287244546914817512959e-15,
	 7.4706210758024567437098915584e-18,
	-4.4697032274412780547627007999999e-20,
	 2.3024281584337436200523093333333e-22,
	-1.0319144794166698148522666666666e-24,
	 4.06081782748733227008e-27,
	-1.4143836005240913919999999999999e-29,
	 4.391090549669888e-32 
    };

    const double xsml = sqrt(0.5 * DBL_EPSILON * 8.0);

    static int ntj0, first = 1;
    if (first) {
	ntj0 = initds(bj0cs, 19, 0.5 * DBL_EPSILON * 0.1);
	first = 0;
    }

    double y = fabs(x);
    if (y > 4.0) goto L20;

    if (y > xsml)
	return dcsevl(y * 0.125 * y - 1.0, bj0cs, ntj0);
    else
	return 1.0;
    
  L20:
    double ampl, theta;
    d9b0mp(y, ampl, theta);
    return ampl * cos(theta);
}

//-----------------------------------------------------------------------------//
