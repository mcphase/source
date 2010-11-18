/*-----------------------------------------------------------------------------*\
| hypot() for three double arguments                                 hypot3d.cc |
|                                                                               |
| MatPack Libary Release 1.0                                                    |
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

#include <cmath>

//-----------------------------------------------------------------------------//

double hypot (double x, double y, double z)
//
// The standard hypot() function for three arguments taking care
// of overflows and zerodivides. 
//
// Version for double arguments.
//
{
    double f,g,h;
    double ax=fabs(x), ay=fabs(y), az=fabs(z); // use standard fabs()
    if (ax > ay) 
	if (ax > az) {
	    f = ay/ax; g = az/ax; h = ax;
	} else {
	    f = ax/az; g = ay/az; h = az;
	}
    else
	if (ay > az) {
	    f = ax/ay; g = az/ay; h = ay;
	} else if (az != 0){
	    f = ax/az; g = ay/az; h = az;
	} else 
	    return 0;
    return h*sqrt(1+f*f+g*g);
}

//-----------------------------------------------------------------------------//
