/*-----------------------------------------------------------------------------*\
| Matpack special functions - initds(os,nos,eta)                      initds.cc |
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
// int initds (double *os, int nos, double eta)
//
// Initialize the orthogonal series, represented by the array os, so 
// that initds is the number of terms needed to insure the error is no 
// larger than eta.  Ordinarily, eta will be chosen to be one-tenth 
// machine precision. 
//
// Input Arguments:
//
//   OS     double precision array of NOS coefficients in an orthogonal 
//          series. 
//   NOS    number of coefficients in OS. 
//   ETA    single precision scalar containing requested accuracy of 
//          series. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C3A2, REVISION 900315, originally written by Fullerton W., (LANL)
// to C++.
//
//-----------------------------------------------------------------------------//

int initds (double *os, int nos, double eta)
{
    if (nos < 1) {
	Matpack.Error(Mat::UnspecifiedError, "%s: %s", "initds",
		      "Number of coefficients is less than 1");
	return 0;
    }

    int i;
    double err = 0.0;
    for (int ii = 1; ii <= nos; ii++) {
	i = nos - ii;
	err += fabs(os[i]);
	if (err > eta) break;
    }
    
    if (i == nos) {
	Matpack.Error(Mat::UnspecifiedError, "%s: %s", "initds",
		      "Chebyshev series too short for specified accuracy");
	return 0;
    }

    return i;
}

//-----------------------------------------------------------------------------//
