/*-----------------------------------------------------------------------------*\
| Matpack special functions - Cbrt(x) cube root                        dcbrt.cc |
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
// double Cbrt (double x);
//
// Calculates the double precision cube root for double precision argument x. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C2, REVISION 891214, originally written by Fullerton W.,(LANL)
// to C++.
//
//-----------------------------------------------------------------------------//

double Cbrt (double x)
{
    static double cbrt2[5] = { 
	0.62996052494743658238360530363911,
	0.79370052598409973737585281963615,
	1.00000000000000000000000000000000,
	1.25992104989487316476721060727823,
	1.58740105196819947475170563927231 
    };
    
    const int niter = (int)(1.443*log(-0.106*log(0.5*DBL_EPSILON*0.1))+1.0);

    if (x == 0.0) return 0.0;
    
    // split into normalized fraction and a base-2 exponent
    int n;
    double y = frexp(fabs(x),&n);

    int ixpnt = n / 3,
	irem  = n - ixpnt * 3 + 2;
    
    // The approximation below is a generalized Chebyshev series converted 
    // to polynomial form.  The approximation is nearly best in the sense of 
    // relative error with 4.085 digits accuracy.
    float z = (float)y;
    double ret_val = z*(z*(z*0.144586F-0.512653F)+0.928549F)+0.439581F;
    
    // iterate
    for (int iter = 0; iter < niter; ++iter) {
	double cbrtsq = ret_val * ret_val;
	ret_val += (y - ret_val * cbrtsq) / (3*cbrtsq);
    }
    
    // combine again
    // changed 10.03.1996
    // originally: return ldexp(cbrt2[irem] * CopySign(ret_val, x), ixpnt);

    return (x > 0) ? ldexp( cbrt2[irem]*ret_val, ixpnt) 
	           : ldexp(-cbrt2[irem]*ret_val, ixpnt);
}

//-----------------------------------------------------------------------------//
