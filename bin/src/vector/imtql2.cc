//
// imtql2.c   Computes the eigenvalues of a tridiagonal matrix
// -------    for the double precision matrix class of MatPack.
//            See also imtql.c where the complete eigensystem is
//            computed.
//
// MatPack C++ Libary Release 1.0
// Copyright (C) 1990,1991 by Berndt M. Gammel
//

#include "vector.h"


// auxilliary routine
static void eisort (Vector &d);


void Imtql (Vector& d, Vector& e, int  sort, int maxiter)
//
//  Imtql computes the  eigenvalues of a tridiagonal matrix t using 
//  implizit QL transformations. The  matrix t is given by the
//  diagonal elements d[lo..hi] and the subdiagonal elements e[lo+1..hi].
//  The eigenvalues are overwritten on d[lo..hi].
//  e is used as scratchpad and overwritten during computation.
//  If convergence is not reached within 'maxiter' iterations for every
//  eigenvalue an error will be signaled (default = 30).
//  The total operation count is about 3*n*n*n + 30*n*n.
//  If the sort flag is set True  numerically ascending eigenvalues
//  are produced, otherwise if False the eigenvalues are not ordered
//  (default is True).
//
//  Reference:
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
// 
{
    int i,iter,l,m,ii;
    double b,c,f,g,p,r,s;
    int lo = d.Lo();
    int hi = d.Hi();
    
    // the vectors must be conformant 
    if (lo != e.Lo() || hi != e.Hi())
	Matpack.Error(Mat::UnspecifiedError,"Imtql: vectors are not conformant");

    // zero or one dimensional
    if (hi-lo <= 0) return;

    // shift subdiagonal vector 
    for (i = lo+1; i <= hi; i++) e[i-1] = e[i];
    e[hi] = 0;

    for (l = lo; l <= hi; l++) {

	// reset iteration counter for every eigenvalue 
	iter = 0;

	// iteration
	do {

	    // look for small subdiagonal element 
	    for (m = l; m < hi; m++) {
		f = fabs(d[m]) + fabs(d[m+1]);
		if (fabs(e[m]) + f == f) break;
	    }

	    if (m != l) {

		// eigenvalue didn't converge 
		if (iter++ >= maxiter) { 
		   Matpack.Error(Mat::UnspecifiedError,"Imtql: eigenvalue %d not convered within %d iterations",
			l,maxiter);
		   return;
		}
		
		p = d[l];
		g = (d[l+1] - p) / (2.0 * e[l]);
		r = hypot(1.0,g);
		g = d[m] - p + e[l] / (g + (g < 0.0 ? -fabs(r) : fabs(r)));
		s = c = 1.0;
		p = 0.0;

		for (i = m-1; i >= l; i--) {

		    ii = i + 1;
		    f = s * e[i];
		    b = c * e[i];

		    if (fabs(f) < fabs(g)) {
			s = f / g;
			r = hypot(1.0,s);
			e[ii] = g * r;
			c = 1.0 / r;
			s *= c;
		    } else {
			c = g / f;
			r = hypot(1.0,c);
			e[ii] = f * r;
			s = 1.0 / r;
			c *= s;
		    }

		    g = d[ii] - p;
		    r = (d[i] - g) * s + 2.0 * c * b;
		    p = s * r;
		    d[ii] = g + p;
		    g = c * r - b;

		} // for i //

		d[l] -= p;
		e[l] = g;
		e[m] = 0.0;

	    } // if m != l //

	} while (m != l);

     } // for l //

     // order eigenvalues and eigenvectors
     if (sort) eisort(d);
}


static void eisort (Vector &d)
//
//  sort the vector d[lo..hi] in ascending numerical order
//
{
    int i,j,k;
    double p;

    // lowest and highest column index of the matrix
    int lo = d.Lo();
    int hi = d.Hi();

    for (i = lo; i <= hi; i++) {
	k = i;
	p = d[i];
	for (j = i+1; j <= hi; j++)
	    if (d[j] < p) {
		k = j;
		p = d[j];
	    }
	if (k != i) {
	    d[k] = d[i];
	    d[i] = p;
	}
    }
}
