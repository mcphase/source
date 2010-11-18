/*-----------------------------------------------------------------------------*\
| powii: power with two integer arguments                              powii.cc |
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
// int powii (int x, int n)
//   Raise the integer x to an integer power n with the minimal number of 
//   multiplications. The exponent n must be non-negative.
//-----------------------------------------------------------------------------//

int powii (int x, int n)
{
  if (n < 0) Matpack.Error("powii: exponent must be non-negative");
  int g = 1;

iterate:

  if (n & 1) g *= x;	// n is odd
  if ( (n /= 2) ) {	// n/2 is non zero
    x *= x;
    goto iterate;
  }

  return g;
}

//-----------------------------------------------------------------------------//
