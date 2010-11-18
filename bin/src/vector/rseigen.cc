/*----------------------------------------------------------------------------*\
| Computes the eigensytem of a real symmetric matrix                rseigen.cc |
| for the double precision matrix class of MatPack                             |
|                                                                              |
| MatPack Libary Release 1.0                                                   |
| Copyright (C) 1990-1995 by Berndt M. Gammel                                  |
|                                                                              |
| Permission to use, copy, modify, and distribute this software and its        |
| documentation for any non-commercial purpose and without fee is hereby       |
| granted, provided that the above copyright notice appear in all copies       |
| and that both that copyright notice and this permission notice appear        |
| in supporting documentation. This software is provided "as is" without       |
| express or implied warranty.                                                 |
|                                                                              |
\*----------------------------------------------------------------------------*/

#include "vector.h"

//----------------------------------------------------------------------------//

void EigenSystemSymmetric (Matrix& z, Vector& d, int sort, int maxiter)
//
//  Driver routine to compute  the eigenvalues  and normalized
//  eigenvectors of  the real symmetric matrix z, given by the
//  lower triangle in z[lo..hi,lo..hi]. The eigenvalues are re-
//  turned in d[lo..hi] in ascending sequence if sort = True,
//  otherwise not ordered for  sort = False. The associated
//  eigenvectors overwrite the given matrix z. The storage re-
//  quirement is n*n + 2*n double.
//  The vector d must already be allocated by the user.
//
//  Reference:
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
// 
{
    // allocate auxilliary vector 
    Vector e( z.Rlo(), z.Rhi() );

    // transform to tridiagonal form 
    Tred(z,d,e);

    // compute eigensystem 
    Imtql(z,d,e,sort,maxiter);
}

//----------------------------------------------------------------------------//

void EigenValuesSymmetric (Matrix& z, Vector& d, int sort, int maxiter)
//
//  Driver routine to compute  the eigenvalues only
//  the real symmetric matrix z, given by the
//  lower triangle in z[lo..hi,lo..hi]. The eigenvalues are re-
//  turned in d[lo..hi] in ascending sequence if sort = True,
//  otherwise not ordered for  sort = False. 
//  The storage requirement is n*n + 2*n double.
//  The vector d must already be allocated by the user.
//
//  Reference:
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
// 
{
    // allocate auxilliary vector 
    Vector e( z.Rlo(), z.Rhi() );

    // transform to tridiagonal form 
    Tred(z,d,e);

    // compute eigensystem 
    Imtql(d,e,sort,maxiter);
}

//----------------------------------------------------------------------------//
