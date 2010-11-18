/*-----------------------------------------------------------------------------*\
| Matpack functions - Chebyshev polynomials T_n (x)               chebysheft.cc |
|                                                                               |
| Last change: June 1, 1997							|
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
// void ChebyshevT (int n, double x, double &y, double &dy, double &d2y);
//
// Computes the value of the Chebyshev polynomial of degree n 
// and its first and second derivatives at a given point 
//
// Input:
//   n  = degree of the polynomial >= 0
//   x  = point in which the computation is performed, -1 <= x <= 1
// Output:
//   y  = value of the polynomial in x 
//   dy = value of the first derivative in x 
//   d2y= value of the second derivative in x 
//
// Note: 
//   This C++ implementation is based on the Fortran function 
//      VACHPO
//   from
//       "Fortran routines for spectral methods"
//   by  Daniele Funaro 
//       Department of Mathematics 
//       University of Pavia 
//       Via Abbiategrasso 209, 27100 Pavia, Italy 
//       e-mails: fun18@ipvian.ian.pv.cnr.it 
//                funaro@dragon.ian.pv.cnr.it 
//-----------------------------------------------------------------------------//

void ChebyshevT (int n, double x, double &y, double &dy, double &d2y)
{
  // check parameters
  if (n < 0 || fabs(x) > 1.0)
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "ChebyshevT",
		  "bad argument");

  double ym, yp, dym, dyp, d2ym, d2yp;

  y = 1.0;
  dy = 0.0;
  d2y = 0.0;
  if (n == 0) return;

  y = x;
  dy = 1.0;
  d2y = 0.0;
  if (n == 1) return;

  yp = 1.0;
  dyp = 0.0;
  d2yp = 0.0;
  for (int k = 2; k <= n; ++k) {
    ym = y;
    y = x * 2.0 * y - yp;
    yp = ym;
    dym = dy;
    dy = x * 2.0 * dy + yp * 2.0 - dyp;
    dyp = dym;
    d2ym = d2y;
    d2y = x * 2.0 * d2y + dyp * 4.0 - d2yp;
    d2yp = d2ym;
  }
}

//-----------------------------------------------------------------------------//
