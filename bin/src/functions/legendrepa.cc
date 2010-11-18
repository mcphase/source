/*-----------------------------------------------------------------------------*\
| Matpack functions - associated Legendre polynomial P_lm(x)      legendrepa.cc |
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
// double LegendreP (int l, int m, double x);
//
// Computes the value of the associated Legendre polynomial P_lm (x) 
// of order l at a given point.
//
// Input:
//   l  = degree of the polynomial  >= 0
//   m  = parameter satisfying 0 <= m <= l, 
//   x  = point in which the computation is performed, range -1 <= x <= 1.
// Returns:
//   value of the polynomial in x
//
//-----------------------------------------------------------------------------//

double LegendreP (int l, int m, double x)
{
  // check parameters
  if (m < 0 || m > l || fabs(x) > 1.0) {
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "LegendreP",
		  "bad argument");
    return NAN;
  }

  double pmm = 1.0;
  if (m > 0) {
    double h = sqrt((1.0-x)*(1.0+x)),
           f = 1.0;
    for (int i = 1; i <= m; i++) {
      pmm *= -f * h;
      f += 2.0; 
    }
  }
  if (l == m)
    return pmm;
  else {
    double pmmp1 = x * (2 * m + 1) * pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      double pll;
      for (int ll = m+2; ll <= l; ll++) {
	pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
	pmm = pmmp1;
	pmmp1 = pll;
      }
      return pll;
    }
  }
}

//-----------------------------------------------------------------------------//
