/*----------------------------------------------------------------------------*\
| special matrices for the matrix classes of MatPack                specmat.cc |
|                                                                              |
| MatPack Libary Release 1.0                                                   |
| Copyright (C) 1990-1994 by Berndt M. Gammel                                  |
|                                                                              |
| Permission to use, copy, modify, and distribute this software and its        |
| documentation for any non-commercial purpose and without fee is hereby       |
| granted, provided that the above copyright notice appear in all copies       |
| and that both that copyright notice and this permission notice appear        |
| in supporting documentation. This software is provided "as is" without       |
| express or implied warranty.                                                 |
|                                                                              |
\*---------------------------------------------------------------------------*/

#include "vector.h"
#include "mprandom.h"

// needed for initializing random generator
#include <ctime>

//-----------------------------------------------------------------------------//

Matrix Vandermonde (const Vector& C)
//
// Returns the Vandermonde matrix whose second last column is given by
// the vector C. The j-th column of a Vandermonde matrix is given C^(n-j).
// The index range of the matrix returned will be that of the given vector.
//
{
    int lo = C.Lo(), hi = C.Hi(); 
    if (hi-lo < 0) return NullMatrix;
    Matrix V(lo,hi,lo,hi);
    for (int i = lo; i <= hi; i++) {
	V[i][hi] = 1;
	for (int k = hi-1; k >= lo; k--)
	  V[i][k] = V[i][k+1]*C[i];
    }
    return V.Value();
}

//-----------------------------------------------------------------------------//

Matrix Hilbert (int N)
//
// Returns the Hilbert matrix, which is the N by N matrix with 
// elements 1/(i+j-1), which is a famous example of a badly conditioned matrix.
// The index range for the matrix returned will be (1,N,1,N).
//
{
    if (N < 1) return NullMatrix;
    Matrix H(1,N,1,N);
    for (int i = 1; i <= N; i++) 
      for (int j = 1; j <= N; j++) 
	H[i][j] = 1.0/(i+j-1);
    return H.Value();
}

//-----------------------------------------------------------------------------//

Matrix InverseHilbert (int N)
//
// Returns the inverse of the N by N Hilbert matrix with elements
// 1/(i+j-1), which is a famous example of a badly conditioned
// matrix.  The result is exact for  N less than about 15.
// The index range for the matrix returned will be (1,N,1,N).
//
{
    if (N < 1) return NullMatrix;
    double p,r;
    Matrix H(1,N,1,N);
    
    p = N;
    for (int i = 1; i <= N; i++) {
	if (i > 1) p = ((N-i+1)*p*(N+i-1))/sqr(i-1);
	r = p*p;
	H[i][i] = r/(2*i-1);
	for (int j = i+1; j <= N; j++) {
	    r = -((N-j+1)*r*(N+j-1))/sqr(j-1);
	    H[i][j] = H[j][i] = r/(i+j-1);
	}
    }
    return H.Value();
}

//-----------------------------------------------------------------------------//

Matrix Wilkinson (int N)
//
// Returns J. H. Wilkinson's eigenvalue test matrix WN+.
// It is a symmetric, tridiagonal matrix with pairs of nearly,
// but not exactly, equal eigenvalues.  
// The index range for the matrix returned will be (1,N,1,N).
// The most frequently used case is Wilkinson(21).
// For example, Wilkinson(7) is
//  
//                3  1  0  0  0  0  0
//                1  2  1  0  0  0  0
//                0  1  1  1  0  0  0
//                0  0  1  0  1  0  0
//                0  0  0  1  1  1  0
//                0  0  0  0  1  2  1
//                0  0  0  0  0  1  3
//
{
    if (N < 2) return NullMatrix;
    Matrix W(1,N,1,N, 0.0);
    for (int i = 1; i < N; i++) {
	W[i][i] = fabs(i - 0.5*(N+1));
	W[i][i+1] = W[i+1][i] = 1;
    }
    W[N][N] = 0.5*(N-1);
    return W.Value();
}

//-----------------------------------------------------------------------------//

Matrix RandomMatrix (int N, long seed)
//
// Returns a N by N matrix with uniformly distributed deviates in the 
// range [-0.5,0.5]. The index range for the matrix returned will be (1,N,1,N).
// If the seed argument is omitted or set to Zero then you will get a
// new random matrix whenever you call this function. 
// If you supply a seed (negative and positive numbers are  
// n o t distinct) you can reproduce the matrix of random numbers.
// The random generator Rnd() is used. 
//
{
  static Ran002 rnd;
  static UniformDistribution U(-0.5,0.5,&rnd);

  if (seed) { // a seed was given
    rnd = Ran002(seed);
    U   = UniformDistribution(-0.5,0.5,&rnd);
  }

  if (N < 1) return NullMatrix;
  Matrix R(1,N,1,N);
  
  double *r = R.Store();
  for (int i = 0; i < R.Elements(); i++) r[i] = U();
    
  return R.Value();
}

//-----------------------------------------------------------------------------//

