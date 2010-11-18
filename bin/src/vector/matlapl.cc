/*------------------------------------------------------------------------------*
| Laplacian of a matrix                                              matlapl.cc |
|                                                                               |
| MatPack Libary Release 1.0                                                    |
| Copyright (C) 1990-1994 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to use, copy, modify, and distribute this software and its         |
| documentation for any non-commercial purpose and without fee is hereby        |
| granted, provided that the above copyright notice appear in all copies        |
| and that both that copyright notice and this permission notice appear         |
| in supporting documentation. This software is provided "as is" without        |
| express or implied warranty.                                                  |
|                                                                               |
*------------------------------------------------------------------------------*/

#include "vector.h"


Matrix Laplacian (Matrix& A)
//
// Applies a Laplacian filter to the matrix and returns the result in
// a matrix of the same dimension.
//
// The Laplacian 'del A' of the matrix A is defined to be
//
// del A(i,j) = ( A(i+1,j)+A(i,j+1)+A(i-1,j)+A(i,j-1)-4*A(i,j) ) / h  
//
// with  h = 1. For the elements outside the index range
// the value of the adjacent element within the matrix is 
// assumed.
//
{
    int i,j,m,n,r,s;

    m = A.Clo();
    n = A.Chi();
    r = A.Rlo();
    s = A.Rhi();

    Matrix B(r,s,m,n);

    // calculate Laplacian

    // center
    for (i = r+1; i < s; i++)
      for (j = m+1; j < n; j++)
	B(i,j) = A(i+1,j)+A(i,j+1)+A(i-1,j)+A(i,j-1)-4*A(i,j);

    // sides
    for (i = r+1; i < s; i++) {
	B(i,m) = A(i+1,m)+A(i,m+1)+A(i-1,m)-3*A(i,m);
	B(i,n) = A(i+1,n)+A(i-1,n)+A(i,n-1)-3*A(i,n); 
    }
    for (j = m+1; j < n; j++) {
	B(r,j) = A(r+1,j)+A(r,j+1)+A(r,j-1)-3*A(r,j);
	B(s,j) = A(s,j+1)+A(s-1,j)+A(s,j-1)-3*A(s,j);
    }
    
    // corners
    B(r,m) = A(r+1,m)+A(r,m+1)-2*A(r,m);  
    B(r,n) = A(r+1,j)+A(r,n-1)-2*A(r,n);
    B(s,m) = A(s,m+1)+A(s-1,m)-2*A(s,m);
    B(s,n) = A(s-1,n)+A(s,n-1)-2*A(s,n);
    
    return B.Value();
}

