/*-----------------------------------------------------------------------------*\
| Matpack special functions - Erf(x) error function                     derf.cc |
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
// double Erf (double x);
//
// Erf(x) calculates the double precision error function for double 
// precision argument x. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C8A, L5A1E, REVISION 920618, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for erf        on the interval  0.          to  1.00000E+00 
//                                        with weighted error   1.28E-32 
//
//                                         log weighted error  31.89 
//                               significant figures required  31.05 
//                                    decimal places required  32.55 
//
//-----------------------------------------------------------------------------//

double Erf (double x)
{

    static double erfcs[21] = { 
	-0.049046121234691808039984544033376,
	-0.14226120510371364237824741899631,
	 0.010035582187599795575754676712933,
	-5.7687646997674847650827025509167e-4,
	 2.7419931252196061034422160791471e-5,
	-1.1043175507344507604135381295905e-6,
	 3.8488755420345036949961311498174e-8,
	-1.1808582533875466969631751801581e-9,
	 3.2334215826050909646402930953354e-11,
	-7.9910159470045487581607374708595e-13,
	 1.7990725113961455611967245486634e-14,
	-3.7186354878186926382316828209493e-16,
	 7.1035990037142529711689908394666e-18,
	-1.2612455119155225832495424853333e-19,
	 2.0916406941769294369170500266666e-21,
	-3.253973102931407298236416e-23,
	 4.7668672097976748332373333333333e-25,
	-6.5980120782851343155199999999999e-27,
	 8.6550114699637626197333333333333e-29,
	-1.0788925177498064213333333333333e-30,
	 1.2811883993017002666666666666666e-32 
    };

    const double sqrtpi = 1.77245385090551602729816748334115,
	         xbig   = sqrt(-log(sqrtpi * 0.5*DBL_EPSILON )),
	         sqeps  = sqrt(DBL_EPSILON);

    static int nterf, first = 1;
    if (first) {
	nterf = initds(erfcs, 21, 0.5 * DBL_EPSILON * 0.1);
	first = 0;
    }

    double y = fabs(x);
    if (y > 1.0) goto L20;

    // erf(x) = 1.0 - erfc(x)  for  -1.0 <= x <= 1.0 

    if (y <= sqeps) 
	return x * 2.0 * x / sqrtpi;
    else // if (y > sqeps) 
	return x * ( dcsevl(x * 2.0 * x - 1.0, erfcs, nterf) + 1.0 );

    // erf(x) = 1.0 - erfc(x) for abs(x) > 1.0

  L20:
    if (y <= xbig) 
	return CopySign(1.0 - Erfc(y), x);
    else // if (y > xbig)  
	return CopySign(1.0, x);
}

//-----------------------------------------------------------------------------//
