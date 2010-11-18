/*-----------------------------------------------------------------------------*\
| Matpack special functions - Clebsch-Gordan ceofficient             clebsch.cc |
|                                                                               |
| Last change: June 5, 1997							|
|                                                                               |
| Matpack Library Release 1.0                                                   |
| Copyright (C) 1991-1997 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.                                                |
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.                                                    |
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.                       |
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration |
| and installation of Matpack.                                                  |
|                                                                               |
\*-----------------------------------------------------------------------------*/

//#include "matpack.h"
#include "../../include/vector.h"
#include "../../include/mpspecfun.h"

//-----------------------------------------------------------------------------//
//
// double ClebschGordan (double l1, double m1, double l2, double m2, 
//                       double l3, double m3, int &errflag);
//
// Calculate Clebsch-Gordan coefficient using the relation to the
// Wigner 3-j symbol:
//
//                               l1-l2+m3         1/2   ( l1  l2  l3 )
//   (l1 m1 l2 m2 | l3 m3) = (-1)         (2*l3+1)      ( m1  m2 -m3 )  
//
//
// Input Arguments:
// ----------------
//
//   double l1
//   double m1	 
//   double l2
//   double m2	
//   double l3
//   double m3		Parameters in 3j symbol.
//
// Output Arguments:
// -----------------
//
//   int &errflag	Error flag.
//                 	errflag=0  No errors.
//                 	errflag=1  Either l1 < abs(m1) or l1+abs(m1) non-integer.
//                 	errflag=2  abs(l1-l2)<= l3 <= l1+l2 not satisfied.
//                 	errflag=3  l1+l2+l3 not an integer.
//                 	errflag=4  m2max-m2min not an integer.
//                 	errflag=5  m2max less than m2min.
//                 	errflag=6  ndim less than m2max-m2min+1.
//                 	errflag=7  m1+m2-m3 is not zero
//
// Returns:
// --------
//
//   The value of the Clebsch-Gordan coefficient.
//
// References:
// -----------
//  1. See routines in "threejj.cc" and "threejm.cc" for references about
//     the calculation of the Wigner 3-j symbols.
//  2. C++ Implementation for the Matpack C++ Numerics and Graphics Library 
//     by Berndt M. Gammel in June 1997.
//
// Note:
// -----
//  Whenever you have to calculate a series of Clebsch-Gordan coefficients for
//  a range of l-values or m-values you should probably use the 3-j symbol
//  routines. These calculate the 3-j symbols iteratively for a series of
//  l-values or m-values and are therefore much more efficient. Use the relation
//  between Clebsch-Gordan coefficients and Wigner 3-j symbols as given above.
//
//-----------------------------------------------------------------------------//

double ClebschGordan (double l1, double m1, double l2, double m2, 
		      double l3, double m3, int &errflag)
{
  const double err = 0.01;
  double CG = 0.0, m2min, m2max, *cofp;

  // static array for calculation of 3-j symbols
  const int ncof = 100;
  static double cof[ncof];

  // reset error flag
  errflag = 0;
  
  // Check for physical restriction. 
  // All other restrictions are checked by the 3-j symbol routine.
  if ( fabs(m1 + m2 - m3) > err) {
    errflag = 7;
    Matpack.Error(Mat::ArgumentDomain, "%s: %s", "ClebschGordan",
		  "m1 + m2 - m3 is not zero.");
    return 0;
  } 
  
  // calculate minimum storage size needed for ThreeJSymbolM()
  // if the dimension becomes negative the 3-j routine will capture it
  int njm = Nint(min(l2,l3-m1) - max(-l2,-l3-m1) + 1); 
  
  // allocate dynamic memory if necessary
  cofp = (njm > ncof) ? new double[njm] : cof;

  // calculate series of 3-j symbols
  ThreeJSymbolM (l1,l2,l3,m1, m2min,m2max, cofp,njm, errflag);

  // calculated Clebsch-Gordan coefficient
  if (! errflag)
    CG = cofp[Nint(m2-m2min)] * (odd(Nint(l1-l2+m3)) ? -1 : 1) * sqrt(2*l3+1); 

  // free dynamic memory if necessary
  if (njm > ncof) delete [] cofp;

  return CG;
}

//-----------------------------------------------------------------------------//
