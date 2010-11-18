/*-----------------------------------------------------------------------------*\
| powi: raise double to integer power                                   powi.cc |
|                                                                               |
| Last change: 09.01.1997                                                       |
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

//#include "common.h"
#include "../../include/vector.h"

//-----------------------------------------------------------------------------//
// double powi (double x, int n)
//   Raise the double x to an integer power n with the minimal number of 
//   multiplications. n can also be negative.
//-----------------------------------------------------------------------------//

double powi (double x, int n)
{
  int m = (n >= 0) ? n : -n;
  double g = 1.0;

iterate:

  if (m & 1) g *= x;		// m is odd
  if ( (m /= 2) ) {		// m/2 is non zero
    x *= x;
    goto iterate;
  }

  return (n >= 0) ? g : 1/g;	// negative powers
}

//-----------------------------------------------------------------------------//
