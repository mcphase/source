//
// tred.c    Reduces a real symmetric matrix to tridiagonal form
// ------    for the double precision matrix class of MatPack.
//
// MatPack C++ Libary Release 1.0
// Copyright (C) 1990,1991 by Berndt M. Gammel
//

#include "vector.h"


void Tred (Matrix& z, Vector& d, Vector& e)
//
//  Tred  reduces the given lower triangle of a real symmetric
//  matrix z[lo..hi,lo..hi] to tridiagonal form using Householder
//  transformations. The diagonal of the result is returned in
//  d[lo..hi],the subdiagonal in e[lo+1..hi] with e[lo] set to 0.
//  The transformation  matrices are accumulated in z (overwriting
//  the original matrix) for a subsequent call of imtql.
//
//  Reference:
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
// 
{
    int i,j,k,l;
    double f,g,h,hh,scale,*zi,*zj,*zk;

    // lowest and highest column index of the matrix
    int lo = z.Clo();
    int hi = z.Chi();
    
    // columns and rows must have the same range
    if (z.Rlo() != lo || z.Rhi() != hi) 
	Matpack.Error(Mat::UnspecifiedError,"Tred: transformation matrix must be square"); 

    // the vectors must be conformant 
    if (d.Lo() != lo || d.Hi() != hi || e.Lo() != lo || e.Hi() != hi)
	Matpack.Error(Mat::UnspecifiedError,"Tred: vectors and matrices are not conformant");

    for (i = hi; i > lo; i--) {

	zi = z[i];
	l = i - 1;
	h = scale = 0.0;

	// scale row 
	for (k = lo; k <= l; k++)
	    scale += fabs(zi[k]);

	if (scale == 0.0 || l == lo)
	    e[i] = zi[l];

	else {

	    for (k = lo; k <= l; k++) {
		zi[k] /= scale;
		h += zi[k] * zi[k];
	    }
	    f = zi[l];
	    g = (f > 0.0) ? -sqrt(h) : sqrt(h);
	    e[i] = scale * g;
	    h -= f * g;
	    zi[l] = f - g;
	    f = 0.0;

	    for (j = lo; j <= l; j++) {
		zj = z[j];
		zj[i] = zi[j] / h;
		g = 0.0;

		// form element of z*U 
		for (k = lo; k <= j; k++)
		    g += zj[k] * zi[k];
		for (k = j+1; k <= l; k++)
		    g += z[k][j] * zi[k];
		e[j] = g / h;
		f += e[j] * zi[j];
	    }

	    hh = f / (h + h);

	    // form reduced z 
	    for (j = lo; j <= l; j++) {
		f = zi[j];
		g = e[j] - hh * f;
		e[j] = g;
		zj = z[j];
		for (k = lo; k <= j; k++)
		    zj[k] -= f * e[k] + g * zi[k];
            }

	} // else //

	d[i] = h;

    } // for i //

    e[lo] = d[lo] = 0.0;

    // accumulation of transformation matrices 
    for (i = lo; i <= hi; i++) {
	zi = z[i];
	l = i - 1;

	if (d[i] != 0.0)
	    for (j = lo; j <= l; j++) {
		g = 0.0;
		for (k = lo; k <= l; k++)
		    g += zi[k] * z[k][j];

		for (k = lo; k <= l; k++) {
		    zk = z[k];
		    zk[j] -= g * zk[i];
		}
	    }

	d[i] = zi[i];
	zi[i] = 1.0;
	for (j = lo; j <= l; j++)
	    zi[j] = z[j][i] = 0.0;
    }
}
