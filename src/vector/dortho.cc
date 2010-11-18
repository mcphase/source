/*------------------------------------------------------------------------------*
| Ortho() - inplace inversion of a double precision matrix            dortho.cc |
|									        |
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
  

void Ortho (Matrix &A)
//
// Replaces a square double precision Matrix A by its inverse. This is an 
// extremely compact version for inplace matrix inversion.
//
// Ref: Handbook for Automatic Computation, Vol 2,
//      Linear Algebra, Contrib.I/9, Springer 1971
//	Algol routine Ortho2
//
// Optimizations:
//   -- avoid call to index operator that optimizes very badely
//
{
    int g,h,i,j,l,n,nn,n2;
    double s,t;

    // get dimensions
    int lo = A.Rlo();
    int hi = A.Rhi();

    // make shure that the matrix is square
    if (lo != A.Clo() || hi != A.Chi()) 
      Matpack.Error(Mat::UnspecifiedError,"Ortho: non square matrix (%d,%d,%d,%d)",
	   A.Rlo(),A.Rhi(),A.Clo(),A.Chi());

    // allocate intermediate storage
    n = hi-lo+1;
    n2 = (n*(n+1))/2;
    Vector P(lo,hi);
    Vector Q(lo,lo+n2-1);

    // avoid call to index operator that optimizes very badely
    double **a = &A[0], *p = &P[0], *q = &Q[0];
    
    l = lo-1;
    for (i = lo; i <= hi; i++) {
	l++;
        for (s = 0, j = lo; j <= hi; j++) {
	    t = a[j][i];
	    p[j] = t;
      	    s += t * t;
        }
	q[l] = s;
	for (g = i+1; g <= hi; g++) {
	    for (t = 0, j = lo; j <= hi; j++)
		t += p[j] * a[j][g];
      	    l++;
	    q[l] = t;
	    t /= s;
	    for (j = lo; j <= hi; j++)
		a[j][g] -= p[j] * t;
	}
    }
    nn = hi + 2;
    for (i = hi; i >= lo; i--) {
	h = l - i;
	t = q[l];
	for (j = lo; j <= hi; j++) {
	    s = a[j][i];
	    for (g = i+1; g <= hi; g++)
		s -= q[g+h] * a[j][g];
	    a[j][i] = s / t;
	}
	l += i-nn;
    }
    for (i = lo; i <= hi; i++)
	for (j = i+1; j <= hi; j++) {
	    s = a[i][j];
	    a[i][j] = a[j][i];
	    a[j][i] = s;
	}
}
