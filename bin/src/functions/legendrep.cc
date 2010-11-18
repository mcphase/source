/*-----------------------------------------------------------------------------*\
| Matpack functions - Legendre polynomials P_n (x)                 legendrep.cc |
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
// void LegendreP (int n, double x, double &y, double &dy, double &d2y);
//
// Computes the value of the Legendre polynomial of degree n 
// and its first and second derivatives at a given point 
//
// Input:
//   n  = degree of the polynomial  >= 0
//   x  = point in which the computation is performed, -1 <= x <= 1
// Output:
//   y  = value of the polynomial in x 
//   dy = value of the first derivative in x 
//   d2y= value of the second derivative in x 
//
// Note: 
//   This C++ implementation is based on the Fortran function 
//      VALEPO
//   from
//       "Fortran routines for spectral methods"
//   by  Daniele Funaro 
//       Department of Mathematics 
//       University of Pavia 
//       Via Abbiategrasso 209, 27100 Pavia, Italy 
//       e-mails: fun18@ipvian.ian.pv.cnr.it 
//                funaro@dragon.ian.pv.cnr.it 
//-----------------------------------------------------------------------------//

void LegendreP (int n, double x, double &y, double &dy, double &d2y)
{
  // check parameters
  if (n < 0 || fabs(x) > 1.0)
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "LegendreP",
	      "bad argument");

  double c1, c2, c4, ym, yp, dym, dyp, d2ym, d2yp;

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
  for (int i = 2; i <= n; ++i) {
    c1 = (double) i;
    c2 = c1 * 2.0 - 1.0;
    c4 = c1 - 1.0;
    ym = y;
    y = (c2 * x * y - c4 * yp) / c1;
    yp = ym;
    dym = dy;
    dy = (c2 * x * dy - c4 * dyp + c2 * yp) / c1;
    dyp = dym;
    d2ym = d2y;
    d2y = (c2 * x * d2y - c4 * d2yp + c2 * 2.0 * dyp) / c1;
    d2yp = d2ym;
  }
}

//-----------------------------------------------------------------------------//
