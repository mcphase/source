/*------------------------------------------------------------------------------*
| implementation of the LU decomposition for the double             dludcomp.cc |
| precision matrix class of MatPack.                                            |
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

//----------------------------------------------------------------------------//

void Decompose (Matrix& A, IntVector& Index, int &d)
// 
//  Calculates the  LU decomposition of the real (double) matrix A, overwriting
//  the old contents. The vector containing the permutations is returned in 
//  Index. This vector must be allocated (and deleted) by the user and must
//  have the same index range like the matrix A. The sign of the 
//  permutations is returnd in d. This will be neccessary if you compute the 
//  determinant.
//
// Optimizations:
//   -- avoid call to index operator that optimizes very badely
//
{
    int i,imax = 0,j,k;
    double big,dum,sum,temp,*ap;

    static const char *sing = "Decompose: matrix is singular";

    // lowest and highest column index of the matrix
    int lo = A.Clo();
    int hi = A.Chi();
    
    // columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || Index.Lo() != lo || Index.Hi() != hi) {
	Matpack.Error(Mat::NonConformant,
		      "Decompose: rows and columns don't have "
		      "the same index range (%d,%d,%d,&d)",
		      A.Rlo(),A.Rhi(),A.Clo(),A.Chi()); 
	return;
    }

    // allocate auxilliary vector
    Vector V(lo,hi);
    
    // avoid call to index operator that optimizes very badely
    double **a = &A[0];
    double *v = &V[0];
    int *idx = &Index[0];

    // sign of permutations
    d = 1; 

    for (i = lo; i <= hi; i++) {
	big = 0.0;
	ap = a[i];
	for (j = lo; j <= hi; j++)
	  if ((temp = fabs(ap[j])) > big) big = temp;
	if (big == 0.0) {
	    Matpack.Error(Mat::MatrixIsSingular,sing);
	    return;
	}
	v[i] = 1.0 / big;
    }

    for (j = lo; j <= hi; j++) {
	for (i = lo; i < j; i++) {
	    sum = a[i][j];
	    for (k = lo; k < i; k++) sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
	}

	big = 0.0;
	for (i = j; i <= hi; i++) {
	    sum = a[i][j];
	    for (k = lo; k < j; k++)
	    sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
	    if ( (dum = v[i] * fabs(sum)) >= big) {
		big = dum;
		imax = i;
	    }
	}
	
	// interchange columns if neccessary 
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
	if (a[j][j] == 0.0) {
	    Matpack.Error(Mat::MatrixIsSingular,sing);
	    return;
	}
	
	if (j != hi) {
	    dum = 1.0 / a[j][j];
	    for (i = j+1; i <= hi; i++) a[i][j] *= dum;
	}
    }
}

//----------------------------------------------------------------------------//

void Backsubst (Matrix& A, IntVector& Index, Vector& B)
// 
//  Do the backsubstitution for a double right hand side with
//  a double LU decomposed matrix as retuned by Decompose().
//  
// Optimizations:
//   -- avoid call to index operator that optimizes very badely
//
{
    int i,ii,ip,j;
    double sum;
    
    // lowest and highest column index of the matrix
    int lo = A.Clo();
    int hi = A.Chi();

    // columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || B.Lo() != lo || B.Hi() != hi
	|| Index.Lo() != lo || Index.Hi() != hi) {
	Matpack.Error(Mat::NonConformant,
		      "Backsubst: non conformant matrix (%d,%d,%d,%d) or "
		      "vector (%d,%d)",
		      A.Rlo(),A.Rhi(),A.Clo(),A.Chi(),B.Lo(),B.Hi()); 
	return;
    }

    // avoid call to index operator that optimizes very badely
    double **a = &A[0];
    double *b = &B[0];
    int *idx = &Index[0];

    ii = lo-1;
    for (i = lo; i <= hi; i++) {
	
	ip = idx[i];
	sum = b[ip];
	b[ip] = b[i];
	if (ii >= lo)
	  for (j = ii; j <= i-1; j++) 
	    sum -= a[i][j]*b[j];
	else if (sum) 
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

void Backsubst (Matrix& A, IntVector& Index, Matrix& B)
// 
//  Do the backsubstitution for many right hand sides stored
//  in the columns of the matrix B using the LU decomposed 
//  matrix A as returned by Decompose().
//  
// Optimizations:
//   -- avoid call to index operator that optimizes very badely
//
{
    int i,ii,ip,j,k;
    double sum;

    // lowest and highest column index of the matrix
    int lo = A.Clo(),
        hi = A.Chi(),
        bclo = B.Clo(),
        bchi = B.Chi();

    // all columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || B.Rlo() != lo || B.Rhi() != hi
	|| Index.Lo() != lo || Index.Hi() != hi) {
	Matpack.Error(Mat::NonConformant,
		      "Backsubst: non conformant matrices (%d,%d,%d,%d)"
		      " and (%d,%d,%d,%d)",
		      A.Rlo(),A.Rhi(),A.Clo(),A.Chi(),
		      B.Rlo(),B.Rhi(),B.Clo(),B.Chi()); 
	return;
    }

    // avoid call to index operator that optimizes very badely
    double **a = &A[0], **b = &B[0];
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
	    else if (sum) 
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

void Improve (Matrix& A, Matrix& lu, IntVector& Index, Vector& b, Vector& x)
//
// iterative improvement step
//
// Optimizations:
//   -- avoid call to index operator that optimizes very badely
//
{
    int j,i;
    double sdp;

    // lowest and highest column index of the matrix
    int lo = A.Clo();
    int hi = A.Chi();
    
    // dimensional check - columns and rows must have the same range
    if ( A.Rlo()   != lo || A.Rhi()  != hi 
      || lu.Rlo()  != lo || lu.Rhi() != hi 
      || lu.Clo()  != lo || lu.Chi() != hi
      || b.Lo()    != lo || b.Hi()   != hi
      || x.Lo()    != lo || b.Hi()   != hi
      || Index.Lo() != lo || Index.Hi()!= hi) {
	Matpack.Error(Mat::NonConformant,
		      "Improve: non conformant matrix or vector"); 
	return;
    }
    
    // calculate the residuals
    Vector r(lo,hi);

    // avoid call to index operator that optimizes very badely
    double  **a = &A[0];

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

void SolveLinear (Matrix& A, Vector& b)
//
// Solve the linear equation system A*x=b for x. The result x overwrites
// the vector b. A will be destroyed during the calculation.
//
{
    int sgn;
    int lo = A.Clo();
    int hi = A.Chi();
    
    // all columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || b.Lo() != lo || b.Hi() != hi) {
	Matpack.Error(Mat::NonConformant,
		      "SolveLinear: non conformant matrix (%d,%d,%d,%d) "
		      "or vector (%d,%d)",
		      A.Rlo(),A.Rhi(),A.Clo(),A.Chi(),b.Lo(),b.Hi()); 
	return;
    }

    // allocate a permutation table
    IntVector idx(lo,hi);
    
    // lu-decomposition: result overwrites A
    // backsubstitution: result overwrites b 
    Decompose(A,idx,sgn);
    if (Matpack.Ok()) Backsubst(A,idx,b); 
}

//----------------------------------------------------------------------------//

void SolveLinear (Matrix& A, Matrix& B)
//
// Solve the linear equation system A * X = B with many right hand sides for X.
// X = A.Inverse() * B  is returned in B overwriting the previous contents. 
// The matrix A will be destroyed during the calculation.
//
{
    int sgn;
    int lo = A.Clo();
    int hi = A.Chi();
    
    // all columns and rows must have the same range
    if (A.Rlo() != lo || A.Rhi() != hi || B.Rlo() != lo || B.Rhi() != hi) {
	Matpack.Error(Mat::NonConformant,
		      "SolveLinear: nonconformant matrices (%d,%d,%d,%d)"
		      " and (%d,%d,%d,%d)",
		      A.Rlo(),A.Rhi(),A.Clo(),A.Chi(),
		      B.Rlo(),B.Rhi(),B.Clo(),B.Chi()); 
	return;
    }
    
    // allocate a permutation table
    IntVector idx(lo,hi);
    
    // lu-decomposition: result overwrites A       
    // backsubstitution: result overwrites B   
    Decompose(A,idx,sgn);
    if (Matpack.Ok()) Backsubst(A,idx,B);
}

//----------------------------------------------------------------------------//
