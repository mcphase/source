/*----------------------------------------------------------------------------*\
| singular value  decomposition for the double precision           dsvdcomp.cc |
| matrix class of MatPack.                                                     |
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
\*----------------------------------------------------------------------------*/

#include "vector.h"
#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif

void SVDecompose (Matrix& A, Vector& W, Matrix& V)
//
// Given a matrix A[l..m][l..n] the routine SVDecompose() computes its singular
// value decomposition A = U * W * V.Transpose(). The matrix U replaces A on 
// output. The diagonal matrix of singular values W is output as a 
// vector w[l..n]. The matrix V (attention: not the transpose V.Transpose()) 
// is output as v[l..n][l..n].
// m must be greater or equal to n; If it is smaller then A should be filled
// up to square with zero rows !
//
// This is the algorithm from:
// Numerical Recipes in C, W. H. Press, B. P. Flannery,
// S. A. Teukolsky, W. T. Vetterling, 1988.
//
// The original algorithm written in Algol 60 is from:
// Golub and Reinsch, Numerische Mathematik 14(1970), pp 403-470.
//
// Bug fixes:
//
// There have been some complaints and bug reports on the following lines
//
//       if (fabs(r[l])+anorm == anorm) {
// and
//       if (fabs(w[nm])+anorm == anorm) break;
// and
//       if (fabs(f)+anorm != anorm) {
//
// which are flawed under certain contitions, for example with a math
// coprocessor which computes intermediate results with higher precision
// or with compilers which optimize these lines !
//
// The problems with compilers, e.g. Microsoft's C 6.00, which optimizes 
// the first line like
//
//        if (fabs(r[l]) == 0.0) {
//
// can be fixed by putting brackets () around the addition. 
//
// A better fix which resolves also the other problem is to revert to the 
// original Golub and Reinsch method which is to specify a convergence 
// limit like
//
//       if (fabs(f) <= eps) {
//
// These bug fixes are now included
// 
// April 18, 1993, B.M.Gammel
//
{
    // index ranges
    int lo = A.Rlo(),
	m  = A.Rhi(),
	n  = V.Rhi();

    if (A.Clo() != lo || W.Lo() != lo || V.Clo() != lo || V.Rlo() != lo ||
	W.Hi() != n || V.Chi() != n || A.Chi() != n)
      Matpack.Error(Mat::UnspecifiedError,"SVDecompose: non conformant index ranges"); 

    if (m < n) 
      Matpack.Error(Mat::UnspecifiedError,"SVDecompose: You must augment A with extra zero rows");

    int flag,i,its,j,jj,k,l = 0,nm = 0;
    double c,f,h,s,x,y,z;

    double anorm = 0.0, 
	       g = 0.0, 
    scale = 0.0;
    
    // allocate an auxilliary vector
    Vector R(lo,n);

    // avoid call to index operator that optimizes very badely
    double **a = &A[0];
    double **v = &V[0];
    double  *w = &W[0];
    double  *r = &R[0];

    // Housholder reduction to bidiagonal form
    for (i = lo; i <= n; i++) {
	l = i+1;
	r[i] = scale*g;
	g = s = scale = 0.0;
	if (i  <=  m) {
	    for (k = i; k <= m; k++) scale += fabs(a[k][i]);
	    if (scale) {
		for (k = i; k <= m; k++) {
		    a[k][i] /= scale;
		    s += a[k][i]*a[k][i];
		}
		f = a[i][i];
		g = -CopySign(sqrt(s),f);
		h = f*g-s;
		a[i][i] = f-g;
		if (i != n) {
		    for (j = l; j <= n; j++) {
			for (s = 0.0, k=i; k <= m; k++) s += a[k][i]*a[k][j];
			f = s/h;
			for (k = i; k <= m; k++) a[k][j] += f*a[k][i];
		    }
		}
		for (k = i; k <= m; k++) a[k][i] *= scale;
	    }
	}
	w[i] = scale*g;
	g = s = scale = 0.0;
	if (i  <=  m && i != n) {
	    for (k = l; k <= n; k++) scale += fabs(a[i][k]);
	    if (scale) {
		for (k = l; k <= n; k++) {
		    a[i][k] /= scale;
		    s += a[i][k]*a[i][k];
		}
		f = a[i][l];
		g = -CopySign(sqrt(s),f);
		h = f*g-s;
		a[i][l] = f-g;
		for (k = l; k <= n; k++) r[k] = a[i][k]/h;
		if (i != m) {
		    for (j = l; j <= m; j++) {
			for (s = 0.0, k = l; k <= n; k++) s += a[j][k]*a[i][k];
			for (k = l; k <= n; k++) a[j][k] += s*r[k];
		    }
		}
		for (k = l; k <= n; k++) a[i][k] *= scale;
	    }
	}
	anorm = MpMax(anorm,(fabs(w[i])+fabs(r[i])));
    }
    
    // Accumulation of right-hand transformations
    for (i = n; i >= lo; i--) {
	if (i < n) {
	    if (g) {
		for (j = l; j <= n; j++)
		  // double division to avoid underflows
		  v[j][i] = (a[i][j]/a[i][l])/g; 

		for (j = l; j <= n; j++) {
		    for (s = 0.0, k = l; k <= n; k++) s += a[i][k]*v[k][j];
		    for (k = l; k <= n; k++) v[k][j] += s*v[k][i];
		}
	    }
	    for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
	}
	v[i][i] = 1.0;
	g = r[i];
	l = i;
    }
    
    // Accumulation of left-hand transformations
    for (i = n; i >= lo; i--) {
	l = i+1;
	g = w[i];
	if (i < n)
	  for (j = l; j <= n; j++) a[i][j] = 0.0;
	if (g) {
	    g = 1.0/g;
	    if (i != n) {
		for (j = l; j <= n; j++) {
		    for (s = 0.0, k = l; k <= m; k++) s += a[k][i]*a[k][j];
		    f = (s/a[i][i])*g;
		    for (k = i; k <= m; k++) a[k][j] += f*a[k][i];
		}
	    }
	    for (j = i; j <= m; j++) a[j][i] *= g;
	} else {
	    for (j = i; j <= m; j++) a[j][i] = 0.0;
	}
	++a[i][i];
    }
    
    // Diagonalization of the bidiagonal form
    for (k = n; k >= lo; k--) {               // loop over singular values
	for (its = 1; its <= 30; its++) {     // loop over iterations
	    flag = 1;
	    for (l = k; l >= lo; l--) {       // test for splitting
		nm = l-1;                     // note that r[l] is always zero
		
// BUG:
//		if ((fabs(r[l])+anorm) == anorm) {
//
// FIX:
		if (fabs(r[l]) <= anorm * DBL_EPSILON) {

		    flag = 0;
		    break;
		}
// BUG:
//		if ((fabs(w[nm])+anorm) == anorm) break;
//
// FIX:
		if (fabs(w[nm]) <= anorm * DBL_EPSILON) break;
		
	    }

	    if (flag) {
		c = 0.0;                      // cancellation of r[l], if l > 1
		s = 1.0;
		for (i = l; i <= k; i++) {
		    f = s*r[i];

// BUG:
//		    if ((fabs(f)+anorm) != anorm) {
//
// FIX:
		    if  (fabs(f) > anorm*DBL_EPSILON) {

			g = w[i];
			h = hypot(f,g);
			w[i] = h;
			h = 1.0/h;
			c = g*h;
			s = (-f*h);
			for (j = lo; j <= m; j++) {
			    y = a[j][nm];
			    z = a[j][i];
			    a[j][nm] = y*c+z*s;
			    a[j][i] = z*c-y*s;
			}
		    }
		}
	    }
	    
	    z = w[k];
	    if (l == k) {             // Convergencs
		if (z < 0.0) {        // make non negative singular value
		    w[k] = -z;
		    for (j = lo; j <= n; j++) v[j][k] = (-v[j][k]);
		}
		break;
	    }
	    
	    if (its == 30) 
	      Matpack.Error(Mat::UnspecifiedError,"SVDecompose: No convergence in 30 iterations");
	    
	    // Shift from botton 2-by-2 minor
	    x = w[l];
	    nm = k-1;
	    y = w[nm];
	    g = r[nm];
	    h = r[k];
	    f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	    g = hypot(f,1.0);
	    f = ((x-z)*(x+z)+h*((y/(f+CopySign(g,f)))-h))/x;

	    // next QR transformation
	    c = s = 1.0;
	    for (j = l; j <= nm; j++) {
		i = j+1;
		g = r[i];
		y = w[i];
		h = s*g;
		g = c*g;
		z = hypot(f,h);
		r[j] = z;
		c = f/z;
		s = h/z;
		f = x*c+g*s;
		g = g*c-x*s;
		h = y*s;
		y = y*c;
		for (jj = lo; jj <= n; jj++) {
		    x = v[jj][j];
		    z = v[jj][i];
		    v[jj][j] = x*c+z*s;
		    v[jj][i] = z*c-x*s;
		}
		z = hypot(f,h);
		w[j] = z;                 // rotation can be arbitrary if z=0
		if (z) {
		    z = 1.0/z;
		    c = f*z;
		    s = h*z;
		}
		f = (c*g)+(s*y);
		x = (c*y)-(s*g);
		for (jj = lo; jj <= m; jj++) {
		    y = a[jj][j];
		    z = a[jj][i];
		    a[jj][j] = y*c+z*s;
		    a[jj][i] = z*c-y*s;
		}
	    }
	    r[l] = 0.0;
	    r[k] = f;
	    w[k] = x;
	}
    }
}


void SVBacksubst (Matrix& u, Vector& w, Matrix& v, Vector& b, Vector& x)
//
// Solves A*X=B for a vector X, where A is given by the arrays u[l..m][l..n],
// w[l..n] and v[l..n][l..n] as returned by SVDecompose(). b[l..n] is the input
// right-hand side. x[l..n] is the output solution vector. No input matrices or
// vectors are destroyed, so the routine can be called with different b's.
//
{
    // index ranges
    int lo = u.Rlo(),
	m  = u.Rhi(),
	n  = v.Rhi();

    if (u.Clo() != lo || w.Lo()  != lo || v.Clo() != lo || 
	v.Rlo() != lo || b.Lo()  != lo || x.Lo()  != lo || 
	w.Hi()  != n  || v.Chi() != n  || u.Chi() != n  ||
	b.Hi()  != n  || x.Hi()  != n)
      Matpack.Error(Mat::UnspecifiedError,"SVBacksubst: non conformant index ranges"); 

    int jj,j,i;
    double s;
    
    // auxilliary vector
    Vector tmp(lo,n);

    for (j = lo; j <= n; j++) {
	s = 0.0;
	if (w[j]) {
	    for (i = lo; i <= m; i++) s += u[i][j]*b[i];
	    s /= w[j];
	}
	tmp[j] = s;
    }
    for (j = lo; j <= n; j++) {
	s = 0.0;
	for (jj = lo; jj <= n; jj++) s += v[j][jj]*tmp[jj];
	x[j] = s;
    }
}


void SVSolveLinear (Matrix& A, Vector& b, Vector& x,
		    double& condition, double threshold)
//
// Solve the linear equation system A * x = b for x, using a singular value
// decomposition. This is a powerful method for dealing with equations
// with matrices which are ill-conditioned or singular - in these cases
// the LU-decomposition (c.f. Decompose(), Backsubst(), SolveLinear()) 
// will fail to give satisfactory results.
// A is given by the square array A[l..n][l..n], the right-hand side
// is given by b[l..n]. The result is returned in vector x[l..n].
// The matrix A will be destroyed during the calculation.
// The condition number (quotient of the largest and the smallest singular
// value) is returned. A threshold for the singular value can be given. 
// (defaults to 0.0 if omitted). That means all singular values smaller 
// than threshold*wmax are set to zero, where wmax is the maximum singular 
// value. You have to experiment with your application to find a good choice.
//
{
    // index ranges
    int lo = A.Rlo(), hi = A.Rhi();

    if (A.Clo() != lo || A.Chi() != hi || 
	b.Lo()  != lo || b.Hi()  != hi ||
	x.Lo()  != lo || x.Hi()  != hi)
      Matpack.Error(Mat::UnspecifiedError,"SVSolveLinear: non conformant index ranges"); 
      
    Vector w(lo,hi);
    Matrix V(lo,hi,lo,hi);
    double wmin,wmax;
    int i;

    SVDecompose(A,w,V);             // svd the square matrix A

    wmax = Max(w);                  // find maximal singular value
    condition = wmax/Min(w);        // compute condition number TODO: ZERODIV !!

    wmin = wmax * threshold;        // get threshold

    for (i = lo; i <= hi; i++)      // zero out small singular values
      if (w(i) < wmin) w(i) = 0.0;

    SVBacksubst(A,w,V,b,x);         // now backsubstitute
}
