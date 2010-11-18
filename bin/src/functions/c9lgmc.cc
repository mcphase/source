/*-----------------------------------------------------------------------------*\
| Matpack special functions - c9lgmc() log gamma correction term      c9lgmc.cc |
|                                                                               |
| Last change: May 29, 1997							|
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
// april 1978 edition.  w. fullerton c3, los alamos scientific lab.
//
// Compute the log gamma correction term for large abs(z) when real(z) >= 0.0 
// and for large abs(imag(y)) when real(z) < 0.0.  
// We find c9lgmc so that
//   clog(cgamma(z)) = 0.5*alog(2.*pi) + (z-0.5)*clog(z) - z + c9lgmc(z).
//-----------------------------------------------------------------------------//

complex<double> c9lgmc (const complex<double> &zin)
{
  complex<double> z, z2inv, retval;

  static double bern[11] = {
     0.083333333333333333,
    -0.0027777777777777778,
     0.00079365079365079365,
    -0.00059523809523809524,
     0.00084175084175084175,
    -0.0019175269175269175,
     0.0064102564102564103,
    -0.029550653594771242,
     0.17964437236883057,
    -1.3924322169059011,
    13.402864044168392
  };

  static int nterm = 0;
  static double bound = 0.0, 
                xbig  = 0.0,
                xmax  = 0.0;

  if (nterm == 0) {
    nterm = int(-0.30*log(0.5 * DBL_EPSILON));
    bound = pow(0.1170*nterm*(0.1*0.5*DBL_EPSILON),-1.0/(2.0*nterm-1.0));
    xbig = 1.0/sqrt(0.5 * DBL_EPSILON);
    xmax = exp(MpMin(log(DBL_MAX/12.0),-log(12.*DBL_MIN)));
  }

  z = zin;
  double x = real(z),
         y = imag(z),
         cabsz = abs(z);

  if (x < 0.0 && abs(y) < bound) {    
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "c9lgm",
		  "not valid for negative real(z) and small abs(imag(z))");	
    return complex<double>(NAN,NAN);
  }

  if (cabsz < bound) {
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "c9lgm",
		  "not valid for small cabs(z)");	
    return complex<double>(NAN,NAN);
  }

  if (cabsz >= xmax) {
    Matpack.Warning(Mat::Underflow, "%s: %s", "c9lgm",
		    "z so big c9lgmc underflows");
    return complex<double>(0.0, 0.0);

  } else {

    if (cabsz >= xbig) return 1.0/(12.0*z);
    
    z2inv = 1.0/(z*z);
    retval = complex<double>(0.0,0.0);
    for (int i = 1; i <= nterm; i++)
      retval = bern[nterm - i] + retval * z2inv;
    return  retval / z;
  }
} 

//-----------------------------------------------------------------------------//
