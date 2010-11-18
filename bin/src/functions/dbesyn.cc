/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselY(n,x) 			      dbesyn.cc |
|                                                                               |
| MatPack Library Release 1.0                                                   |
| Copyright (C) 1991,1995 by Berndt M. Gammel                                   |
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
// double BesselY (int n, double x)
//
// BesselY(n,x) calculates the double precision Bessel function of the 
// second kind of order n for double precision argument x. 
//
// Implemented by B. M. Gammel, last revision 13.03.1996
// 
//-----------------------------------------------------------------------------//

double BesselY (int n, double x)
{
    if (n < 0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "BesselY",
		      "index n less than 0");
	return NAN;
    } else if (n == 0)
	return BesselY0(x);
    else if (n == 1)
	return BesselY1(x);
    else {
	double tox = 2.0/x,
	       by  = BesselY1(x),
	       bym = BesselY0(x);
	for (int j = 1; j < n; j++) {
	    double byp = j * tox * by - bym;
	    bym = by;
	    by = byp;
	}
	return by;
    }
} 

//-----------------------------------------------------------------------------//
