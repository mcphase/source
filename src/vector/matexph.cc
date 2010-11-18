/*----------------------------------------------------------------------------*\
| matrix exponential of a hermitean matrix                          matexph.cc |
|                                                                              |
| MatPack Libary Release 1.0                                                   |
| Copyright (C) 1990-1995 by Berndt M. Gammel                                  |
|                                                                              |
| This contribution to Matpack is from Moritz Hilf, May 9, 1995                |
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


ComplexMatrix MatrixExpHermitean (const ComplexMatrix& A)
//
// Returns the matrix exponential of a hermitean matrix A.
// The argument matrix is not checked for hermiticity !
//
{
    // get the dimension information
    int lo = A.Clo(),
        hi = A.Chi();
    
    // columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi) 
	Matpack.Error(Mat::UnspecifiedError,"MatrixExpHermitean: hermitean matrix must be square"); 
    
    Matrix z(lo,hi,lo,hi), 
           zr(lo,hi,lo,hi),
           zi(lo,hi,lo,hi);

    Vector d(lo,hi);

    // create packed form of hermitean matrix A in matrix z
    for (int i = lo; i <= hi; i++) {	
	for (int j = lo; j < i; j++) {
	    z[j][i] = imag(A[i][j]);
	    z[i][j] = real(A[i][j]);
	}
	z[i][i] = real(A[i][i]);
    }
    
    // diagonalize z
    EigenSystemHermitean(z,d,zr,zi,false,30);

    // combine real and imaginary parts of eigenvectors (zr,zi) in complex matrix
    ComplexMatrix U(zr,zi);
    
    // calculate matrix exponential
    U = U * ComplexDiagonal( Exp(d) ) * U.Hermitean();
    
    return U.Value();    
}













