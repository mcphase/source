/*-----------------------------------------------------------------------------*\
| Matpack functions - spherical harmonics Y_lm(theta,phi)          harmonicy.cc |
|                                                                               |
| Last change: Jun 17, 1997							|
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
// Computes the spherical harmonics Y_lm(theta,phi) with l and m 
// integers satisfying -l <= m <= l and arbitrary angles theta and phi. 
//
//-----------------------------------------------------------------------------//

complex<double> SphericalHarmonicY (int l, int m, double theta, double phi)
{
  complex<double> e; 
  double p,f,smphi;
  int sm;
    
  // save original m with sign
  sm = m;

  // make m positive
  if (m < 0) m = -m;  

  // check parameters
  if (m > l) {
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "SphericalHarmonicY",
		  "m aout of range -l <= m <= l");
    return complex<double>(NAN,NAN);
  }

  // normalization factor
  f = sqrt( (2*l+1)/M_4PI * Fac(l-m) / Fac(l+m) );
    
  // associated Legendre polynomial
  p = LegendreP(l,m,cos(theta));

  // phase factor with signed m
  smphi = sm*phi;
  e = complex<double>(cos(smphi),sin(smphi));

  // sign convention ... corrected (m+sm) to (m-sm) M. rotter 11.1.08
  sm = odd( (m-sm)/2 ) ? -1 : 1;

  // total
  return (sm * f * p) * e; 
}

//-----------------------------------------------------------------------------//
