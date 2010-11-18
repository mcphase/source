/*-----------------------------------------------------------------------------*\
| Matpack special functions - LogRel(z) = log(1+z)                    clnrel.cc |
|                                                                               |
| Last change: May 29, 1997							|
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
// LogRel(z) = log(1+z) with relative error accuracy near z = 0.
//
// April 1977 version.  W. Fullerton, c3, Los Alamos Scientific Lab.
//
// let   rho = abs(z)  and
//       r**2 = abs(1+z)**2 = (1+x)**2 + y**2 = 1 + 2*x + rho**2 .
// now if rho is small we may evaluate LogRel(z) accurately by
//       log(1+z) = complex (log(r), arg(1+z))
//                = complex (0.5*log(r**2), arg(1+z))
//                = complex (0.5*LogRel(2*x+rho**2), arg(1+z))
//
//-----------------------------------------------------------------------------//

complex<double> LogRel(const complex<double> &z)
{
  if (abs(1.0 + z) < sqrt(DBL_EPSILON)) 
    Matpack.Warning(Mat::PartialPrecisionLoss, "%s: %s", "LogRel",
		    "answer less than half precision because z too near -1");
  double rho = abs(z);
  if (rho > 0.375) return log(1.0 + z);
  return complex<double>(0.5*LogRel(2.0*real(z)+rho*rho), arg(1.0+z));
}

//-----------------------------------------------------------------------------//
