/*------------------------------------------------------------------------------*
| implementation of the LU decomposition for the double             cludcomp.cc |
| precision complex matrix class of MatPack.                                    |
|                                                                               |
| MatPack Libary Release 1.0                                                    |
| Copyright (C) 1990-1995 by Berndt M. Gammel                                   |
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

//----------------------------------------------------------------------------//

void Decompose (ComplexMatrix& A, IntVector& Index, int &d)
// 
//  Calculates the LU decomposition of the complex (double) matrix A, overwriting
//  the old contents. The vector containing the permutations is returned in 
//  Index. This vector must be allocated (and deleted) by the user and must
//  have the same index range (offset) like the matrix A. The sign of the 
//  permutations is returnd in d. This will be neccessary if you compute the 
//  determinant.
//
// Optimizations:
//   -- avoid call to index operator that optimizes very badely
//
{
    int i,imax = 0,j,k;
    double big,temp;
    complex<double> sum,dum,*ap;

    const char *SingularMatrix = "Decompose: singular matrix";
    
    // lowest and highest column index of the matrix
    int lo = A.Clo();
    int hi = A.Chi();
    
    // columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || Index.Lo() != lo || Index.Hi() != hi) 
	Matpack.Error(Mat::UnspecifiedError,"Decompose: rows and columns don't have the same range (%d,%d,%d,&d)",
	     A.Rlo(),A.Rhi(),A.Clo(),A.Chi()); 

    // allocate auxilliary vector
    Vector V(lo,hi);
    
    // avoid call to index operator that optimizes very badely
    complex<double> **a = &A[0];
    double *v = &V[0];
    int *idx = &Index[0];

    // sign of permutations
    d = 1; 

    for (i = lo; i <= hi; i++) {
	big = 0.0;
	ap = a[i]; 
	for (j = lo; j <= hi; j++)
	  if ((temp = fabs(real(ap[j])) + fabs(imag(ap[j])) ) > big) big = temp;
	if (big == 0.0) Matpack.Error(Mat::MatrixIsSingular,SingularMatrix);
	v[i] = 1.0 / big;
    }

    for (j = lo; j <= hi; j++) {
	for (i = lo; i < j; i++) {
	    sum = a[i][j];
	    for (k = lo; k < i; k++) 
	      sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
	}

	big = 0.0;
	for (i = j; i <= hi; i++) {
	    sum = a[i][j];
	    for (k = lo; k < j; k++)
	      sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
	    if ( (temp = v[i] * ( fabs(real(sum)) + fabs(imag(sum)) )) >= big) {
		big = temp;
		imax = i;
	    }
	}
	
	// interchange rows if neccessary 
	if (j != imax) {

	    for (k = lo; k <= hi; k++) 
	      MpSwap(a[imax][k],a[j][k]);
	   
	    // sign of determinant changes
	    d = -d;
	    v[imax] = v[j];
	}
	idx[j] = imax;
	
	// instead we could set a[j][j] to a tiny value, say 1e-20, 
	// and continue	the computation
	if (a[j][j] == 0.0) Matpack.Error(Mat::MatrixIsSingular,SingularMatrix);
	
	if (j != hi) {
	    dum = 1.0 / a[j][j];
	    for (i = j+1; i <= hi; i++) a[i][j] *= dum;
	}
    }
}

//----------------------------------------------------------------------------//

void Backsubst (ComplexMatrix& A, IntVector& Index, ComplexVector& B)
// 
//  Do the backsubstitution for a complex right hand b side using
//  the complex LU decomposed matrix A as retuned by Decompose().
//
// Optimizations:
//   -- avoid call to index operator that optimizes very badely
//  
{
    int i,ii,ip,j;
    complex<double> sum;
    
    // lowest and highest column index of the matrix
    int lo = A.Clo(),
        hi = A.Chi();

    // columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || B.Lo() != lo || B.Hi() != hi
	|| Index.Lo() != lo || Index.Hi() != hi) 
      Matpack.Error(Mat::UnspecifiedError,"Backsubst: non conformant matrix (%d,%d,%d,%d) or vector (%d,%d)",
	   A.Rlo(),A.Rhi(),A.Clo(),A.Chi(),B.Lo(),B.Hi()); 

    // avoid call to index operator that optimizes very badely
    complex<double> **a = &A[0], *b = &B[0];
    int *idx = &Index[0];

    ii = lo-1;
    for (i = lo; i <= hi; i++) {
	
	ip = idx[i];
	sum = b[ip];
	b[ip] = b[i];
	if (ii >= lo)
	  for (j = ii; j <= i-1; j++) 
	    sum -= a[i][j]*b[j];
	else if (sum != 0.0) 
	  ii = i;
	b[i] = sum;
    }

    for (i = hi; i >= lo; i--) {
	sum = b[i];
	for (j = i+1; j <= hi; j++) 
	    sum -= a[i][j] * b[j];
	b[i] = sum / a[i][i];
    }
}

//----------------------------------------------------------------------------//

void Backsubst (ComplexMatrix& A, IntVector& Index, ComplexMatrix& B)
// 
//  Do the backsubstitution for many complex right hand sides stored
//  in the columns of the matrix B using the complex LU decomposed 
//  matrix A as returned by Decompose().
//  
// Optimizations:
//   -- avoid call to index operator that optimizes very badely
//
{
    int i,ii,ip,j,k;
    complex<double> sum;
    
    // lowest and highest column index of the matrix
    int lo = A.Clo(),
        hi = A.Chi(),
        bclo = B.Clo(),
        bchi = B.Chi();

    // all columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || B.Rlo() != lo || B.Rhi() != hi
	|| Index.Lo() != lo || Index.Hi() != hi) 
      Matpack.Error(Mat::UnspecifiedError,"Backsubst: non conformant matrices (%d,%d,%d,%d) and (%d,%d,%d,%d)",
	   A.Rlo(),A.Rhi(),A.Clo(),A.Chi(),B.Rlo(),B.Rhi(),B.Clo(),B.Chi()); 

    // avoid call to index operator that optimizes very badely
    complex<double> **a = &A[0], **b = &B[0];
    int *idx = &Index[0];

    for (k = bclo; k <= bchi; k++) {

	ii = lo-1;
	for (i = lo; i <= hi; i++) {
	    
	    ip = idx[i];
	    sum = b[ip][k];
	    b[ip][k] = b[i][k];
	    if (ii >= lo)
	      for (j = ii; j <= i-1; j++) 
		sum -= a[i][j]*b[j][k];
	    else if (sum != 0.0) 
	      ii = i;
	    b[i][k] = sum;
	}
	
	for (i = hi; i >= lo; i--) {
	    sum = b[i][k];
	    for (j = i+1; j <= hi; j++) 
	      sum -= a[i][j] * b[j][k];
	    b[i][k] = sum / a[i][i];
	}
    }
}

//----------------------------------------------------------------------------//

void Improve (ComplexMatrix& A, ComplexMatrix& lu, IntVector& Index, 
	      ComplexVector& b, ComplexVector& x)
//
// iterative improvement step
//
// Optimizations:
//   -- avoid call to index operator that optimizes very badely
//
{
    int j,i;
    complex<double> sdp;

    // lowest and highest column index of the matrix
    int lo = A.Clo();
    int hi = A.Chi();
    
    // dimensional check - columns and rows must have the same range
    if ( A.Rlo()   != lo || A.Rhi()  != hi 
      || lu.Rlo()  != lo || lu.Rhi() != hi 
      || lu.Clo()  != lo || lu.Chi() != hi
      || b.Lo()    != lo || b.Hi()   != hi
      || x.Lo()    != lo || b.Hi()   != hi
      || Index.Lo() != lo || Index.Hi()!= hi)
      Matpack.Error(Mat::UnspecifiedError,"Improve: non conformant matrix or vector"); 
    
    // calculate the residuals
    ComplexVector r(lo,hi);

    // avoid call to index operator that optimizes very badely
    complex<double> **a = &A[0];

    for (i = lo; i <= hi; i++) {
	sdp = -b[i];
	for (j = lo; j <= hi; j++) 
	    sdp += a[i][j] * x[j];
	r[i] = sdp;
    }

    // solve for error term 
    Backsubst(lu,Index,r);

    // subtract it 
    for (i = lo; i <= hi; i++) x[i] -= r[i];
}

//----------------------------------------------------------------------------//

void SolveLinear (ComplexMatrix& A, ComplexVector& b)
//
// Solve the linear equation system A*x=b for x. The result x overwrites
// the vector b. A will be destroyed during the calculation.
//
{
    int sgn;
    int lo = A.Clo();
    int hi = A.Chi();
    
    // all columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || b.Lo() != lo || b.Hi() != hi) 
	Matpack.Error(Mat::UnspecifiedError,"SolveLinear: non conformant matrix (%d,%d,%d,%d) or vector (%d,%d)",
	     A.Rlo(),A.Rhi(),A.Clo(),A.Chi(),b.Lo(),b.Hi()); 
    
    // allocate a permutation table
    IntVector idx(lo,hi);
    
    // solve the equation
    Decompose(A,idx,sgn);   // lu-decomposition: result overwrites A
    Backsubst(A,idx,b);     // backsubstitution: result overwrites b 
}

//----------------------------------------------------------------------------//

void SolveLinear (ComplexMatrix& A, ComplexMatrix& B)
//
// Solve the linear equation system A * X = B with many right hand sides for X.
// X = A.Inverse * B  is returned in B overwriting the previous contents. 
// The matrix A will be destroyed during the calculation.
//
{
    int sgn;
    int lo = A.Clo();
    int hi = A.Chi();
    
    // columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || B.Rlo() != lo || B.Rhi() != hi) 
      Matpack.Error(Mat::UnspecifiedError,"SolveLinear: nonconformant matrices (%d,%d,%d,%d) and (%d,%d,%d,%d)",
	   A.Rlo(),A.Rhi(),A.Clo(),A.Chi(),B.Rlo(),B.Rhi(),B.Clo(),B.Chi()); 

    // allocate a permutation table
    IntVector idx(lo,hi);
    
    // solve the equation
    Decompose(A,idx,sgn);   // lu-decomposition: result overwrites A
    Backsubst(A,idx,B);     // backsubstitution: result overwrites B 
}

//----------------------------------------------------------------------------//
