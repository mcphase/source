/*-----------------------------------------------------------------------------*\
| Matpack functions - Laguerre polynomials L^(a)_n (x)             laguerrel.cc |
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
// void LaguerreL (int n, double a, double x, 
//                 double &y, double &dy, double &d2y);
//
// Computes the value of the Laguerre polynomial of degree n 
// and its first and second derivatives at a given point 
//
// Input:
//   n  = degree of the polynomial 
//   a  = parameter > -1 
//   x  = point in which the computation is performed, x >= 0
// Output:
//   y  = value of the polynomial in x 
//   dy = value of the first derivative in x 
//   d2y= value of the second derivative in x 
//
// Note: 
//   This C++ implementation is based on the Fortran function 
//       VALAPO
//   from
//       "Fortran routines for spectral methods"
//   by  Daniele Funaro 
//       Department of Mathematics 
//       University of Pavia 
//       Via Abbiategrasso 209, 27100 Pavia, Italy 
//       e-mails: fun18@ipvian.ian.pv.cnr.it 
//                funaro@dragon.ian.pv.cnr.it 
//-----------------------------------------------------------------------------//

void LaguerreL (int n, double a, double x, double &y, double &dy, double &d2y)
{
  // check parameters
  if (n < 0 || a <= -1.0 || x < 0.0)
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "LaguerreL",
		  "bad argument");

  double b1, b2, dk, ym, yp, dym, dyp, d2ym, d2yp;
 
  y = 1.0;
  dy = 0.0;
  d2y = 0.0;
  if (n == 0) return;

  y = a + 1.0 - x;
  dy = -1.0;
  d2y = 0.0;
  if (n == 1) return;

  yp = 1.0;
  dyp = 0.0;
  d2yp = 0.0;
  for (int k = 2; k <= n; ++k) {
    dk = (double) k;
    b1 = (dk * 2.0 + a - 1.0 - x) / dk;
    b2 = (dk + a - 1.0) / dk;
    ym = y;
    y = b1 * y - b2 * yp;
    yp = ym;
    dym = dy;
    dy = b1 * dy - yp / dk - b2 * dyp;
    dyp = dym;
    d2ym = d2y;
    d2y = b1 * d2y - dyp * 2.0 / dk - b2 * d2yp;
    d2yp = d2ym;
  }
}

//-----------------------------------------------------------------------------//
