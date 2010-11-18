/*-----------------------------------------------------------------------------*\
| Matpack special functions - LnGamma(z) complex log gamma            clngam.cc |
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
// August 1980 edition.  W. Fullerton c3, Los Alamos Scientific Lab.
// Eventually clngam should make use of c8lgmc for all z except for
// z in the vicinity of 1 and 2.
//-----------------------------------------------------------------------------//

complex<double> LnGamma (const complex<double> &zin)
{
  const double sq2pil = 0.91893853320467274;

  complex<double> retval;

  static double bound = 0.0,
                dxrel = 0.0, 
                rmax  = 0.0;

  if (bound == 0.0) {
    int n = int(-0.30 * log(0.5*DBL_EPSILON));
    bound = 0.1171 * n * pow(0.1 * 0.5 * DBL_EPSILON, -1.0/(2.0*n-1.0));
    dxrel = sqrt(DBL_EPSILON);
    rmax  = DBL_MAX/log(DBL_MAX);
  }
  
  complex<double> z(zin);
  double x = real(zin), 
         y = imag(zin);
  
  double cabsz = abs(z);

  if (cabsz > rmax) {
    Matpack.Error(Mat::Overflow, "%s: %s", "LnGamma",
		  "z so big LnGamma(z) overflows");	 
    return complex<double>(NAN,NAN);
  }

  int n;
  double argsum;
  complex<double> corr(0.0, 0.0);
  complex<double> one(1.0, 0.0);

  if (x >= 0.0 && cabsz  > bound) goto L50;
  if (x <  0.0 && abs(y) > bound) goto L50;
				    
  if (cabsz < bound) goto L20;

  // use the reflection formula for real(z) negative, 
  // abs(z) large, and abs(imag(y)) small.
  if (y > 0.0) z = conj(z);
  corr = exp(-complex<double>(0.0,2.0*M_PI)*z);

  if (real(corr) == 1.0 && imag(corr) == 0.0) {
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "LnGamma", 
		  "z is a negative integer");
    return complex<double>(NAN,NAN);
  }

  retval = sq2pil + 1.0 - complex<double>(0.0,M_PI)*(z-0.5) - LogRel(-corr)
            + (z-0.5) * log(1.0-z) - z - c9lgmc(1.0-z);
  if (y > 0.0) retval = conj(retval);

  if (abs(y) > dxrel) return retval;

  if (0.5 * abs((1.0-corr)*retval/z) < dxrel) {
    Matpack.Warning(Mat::PartialPrecisionLoss, "%s: %s", "LnGamma",
		    "answer less than half precision because z"
		    " too near negative integer");
    return retval;
  }

  // Use the recursion relation for cabs(z) small

 L20:
  if (x >= -0.5 || fabs(y) > dxrel) goto L30;

  if (abs((z-one*((double)((int)(x-0.5)))/x)) < dxrel)    
    Matpack.Warning(Mat::PartialPrecisionLoss, "%s: %s", "LnGamma",
		    "answer less than half precision because z"
		    " too near negative integer"); 
 L30:
  n = (int)(sqrt(bound*bound - y*y) - x + 1.0);
  argsum = 0.0;
  corr = complex<double>(1.0, 0.0);
  for (int i = 1; i <= n; i++) {
    argsum += arg(z);
    corr *= z;
    z += 1.0;
  }

  if (real(corr) == 0.0 && imag(corr) == 0.0) {
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "LnGamma", 
		  "z is a negative integer");
    return complex<double>(NAN,NAN);
  }
  
  corr = -complex<double>(log(abs(corr)), argsum);

  // Use Stirling's approximation for large z

 L50:
  retval = sq2pil + (z-0.5)*log(z) - z + corr + c9lgmc(z);
  return retval;
}

//-----------------------------------------------------------------------------//
