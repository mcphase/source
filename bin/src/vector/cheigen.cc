/*----------------------------------------------------------------------------*\
| Computes the eigensytem of a complex hermitean matrix             cheigen.cc |
| for the complex double precision matrix class of MatPack                     |
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

void EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
			   int sort, int maxiter)
//
//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// 
//  References:
//
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
//
{
    int i;

    // get the dimension information
    int lo = z.Clo(),	
        hi = z.Chi();
    
    // columns and rows must have the same range
    if (z.Rlo() != lo || z.Rhi() != hi) 
	Matpack.Error(Mat::UnspecifiedError,"EigenSystemHermitean: matrix must be square"); 

    // the matrices and vectors must be conformant 
    if ( d.Lo() != lo || d.Hi() != hi
      || zr.Clo() != lo || zr.Chi() != hi || zr.Rlo() != lo || zr.Rhi() != hi
      || zi.Clo() != lo || zi.Chi() != hi || zi.Rlo() != lo || zi.Rhi() != hi)
	Matpack.Error(Mat::UnspecifiedError,"EigenSystemHermitean: vectors and matrices are not conformant");

    // allocate auxilliary vectors e,t1,t2
    Vector e(lo,hi);
    Vector t1(lo,hi);
    Vector t2(lo,hi);

    // zr must be initialized to the identity 
    zr = 0.0;
    for (i = lo; i <= hi; i++) zr[i][i] = 1.0;

    // transform z to tridiagonal form. 
    Chtred(z,d,e,t1,t2);

    // calculate eigensystem of the tridiagonal matrix 
    Imtql(zr,d,e,sort,maxiter);

    // backtransform eigensystem 
    Chtrbk(z,t1,t2,zr,zi);
}

//----------------------------------------------------------------------------//

void EigenValuesHermitean (Matrix& z, Vector& d, int sort, int maxiter)
//
//  Driver routine to compute the  eigenvalues  of a complex hermitian
//  matrix z. The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// 
//  References:
//
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
//
{
    // get the dimension information
    int lo = z.Clo();	
    int hi = z.Chi();
    
    // columns and rows must have the same range
    if (z.Rlo() != lo || z.Rhi() != hi) 
	Matpack.Error(Mat::UnspecifiedError,"EigenValuesHermitean: matrix is not square"); 

    // the matrices and vectors must be conformant 
    if ( d.Lo() != lo || d.Hi() != hi)
	Matpack.Error(Mat::UnspecifiedError,"EigenValuesHermitean: vector and matrix are not conformant");

    // allocate auxilliary vectors e,t1,t2
    Vector e(lo,hi);
    Vector t1(lo,hi);
    Vector t2(lo,hi);

    // transform z to tridiagonal form. 
    Chtred(z,d,e,t1,t2);

    // calculate eigenvalues of the tridiagonal matrix 
    Imtql(d,e,sort,maxiter);
}

//----------------------------------------------------------------------------//

void Chtred (Matrix& z, Vector& d, Vector& e, Vector& t1, Vector& t2)
//
//  Chtred() is the complex analogue to Tred() for a hermitian
//  matrix z, which  is transformed  to real tridiagonal form.
//  The  real parts of the elements of z must be given  in the
//  full lower triangle, the imaginary parts in the transposed
//  positions of the strict upper triangle.The diagonal of the
//  result  is  returned  in  d[lo..hi],  the  subdiagonal  in 
//  e[lo+1..hi] with e[lo] set to 0. After returning z and the 
//  vectors  t1  and  t2  contain  full information  about the
//  unitary transformation applied. 
//
// References:
//
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976
//  c.f. algorithm HTRID3
//  note: variable E2 is not used and array TAU is replaced by
//        two vectors t1 and t2.
//
{
    int i,j,k,l ;
    double f,g,h,fi,gi,hh,si,scale;
    double *zi,*zj,*zk,*zl;

    // lowest and highest column index of the matrix
    int lo = z.Clo();
    int hi = z.Chi();

    // initialize t1 and t2
    t1[hi] = 1.0;
    t2[hi] = 0.0;

    for (i = hi; i >= lo; i--) {

	l = i - 1;
	h = scale = 0.0;
	zi = z[i];

	if (l < lo)
	    e[i] = 0.0;

	else {

	    // scale row 
	    for (k = lo; k <= l; k++)
		scale += fabs(zi[k]) + fabs(z[k][i]);

	    if (scale == 0.0) {
		t1[l] = 1.0;
		t2[l] = e[i] = 0.0;

	    } else {

		for (k = lo; k <= l; k++) {
		    f = zi[k] / scale;
		    g = z[k][i] / scale;
		    h += f * f + g * g;
		    zi[k] = f;
		    z[k][i] = g;
		}

		zl = z[l];
		g = sqrt(h);
		e[i] = scale * g;
		f = hypot(zi[l],zl[i]);

		// form next diagonal element of matrix T
		if (f == 0.0) {
		    t1[l] = -t1[i];
		    si = t2[i];
		    zi[l] = g;

		} else {
		    t1[l] = (zl[i] * t2[i] - zi[l] * t1[i]) / f;
		    si    = (zi[l] * t2[i] + zl[i] * t1[i]) / f;
		    h += f * g;
		    g = 1 + g / f;
		    zi[l] *= g;
		    zl[i] *= g;
		    f = 0.0;
		}

		if (l != lo) {

		    for (j = lo; j <= l; j++) {

			g = gi = 0.0;
			zj = z[j];

			// form element of Z * U 
			for (k = lo; k < j; k++) {
			    zk = z[k];
			    g  +=  zj[k] * zi[k] + zk[j] * zk[i];
			    gi -=  zj[k] * zk[i] - zk[j] * zi[k];
			}

			g  += zj[j] * zi[j];
			gi -= zj[j] * zj[i];

			for (k = j+1; k <= l; k++) {
			    zk = z[k];
			    g  += zk[j] * zi[k] - zj[k] * zk[i];
			    gi -= zk[j] * zk[i] + zj[k] * zi[k];
			}

			// form element of p 
			e[j] = g / h;
			t2[j] = gi / h;
			f += e[j] * zi[j] - t2[j] * zj[i];

		    } // for j //

		    hh = f / (h + h);

		    // form reduced Z //
		    for (j = lo; j <= l; j++) {

			zj = z[j];
			f  = zi[j];
			g  = e[j] - hh * f;
			e[j] = g;
			fi = -zj[i];
			gi = t2[j] - hh * fi;
			t2[j] = -gi;
			zj[j] -= 2 * (f * g + fi * gi);

			for (k = lo; k < j; k++) {
			  zk = z[k];
			  zj[k] -= f*e[k] + g*zi[k] - fi*t2[k] - gi*zk[i];
			  zk[j] -= f*t2[k] + g*zk[i] + fi*e[k] + gi*zi[k];
			}

		    } // for j //

		} // if l != lo //

		for (k = lo; k <= l; k++) {
		    zi[k] *= scale;
		    z[k][i] *= scale;
		}

		t2[l] = -si;

	    } // if scale == 0.0 //

	} // if l < lo //

	d[i] = zi[i];
	zi[i] = scale * sqrt(h);

    } // for i //
} 

//----------------------------------------------------------------------------//

void Chtrbk (Matrix& a, Vector& t1, Vector& t2, Matrix& zr, Matrix& zi)
//
//  Chtrbk() does the backtransformation of the eigenvectors
//  of a  tridiagonal  matrix  using the  unitary transform-
//  ation given by a call to Chtred(),stored in a and t1,t2.
//  The real and  imaginary parts  of the transformed eigen-
//  vectors are returned in zr and zi. 
// 
//  References:
//
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
//  c.f. algorithm HTRIB3
//  note: array TAU is replaced by two vectors t1 and t2.
//
{
    int i,j,k,l;
    double f,g,h,s,si;
    double *ai,*zik,*zrk;

    // lowest and highest column index of the matrix
    int lo = a.Clo();
    int hi = a.Chi();

    for (k = lo; k <= hi; k++) { 
	zik = zi[k];
	zrk = zr[k];
	for (j = lo; j <= hi; j++) {
	    zik[j] = -zrk[j] * t2[k];
	    zrk[j] *= t1[k];
	}
    }

    // recover and apply the Householder matrices 
    for (i = lo+1; i <= hi; i++) {

	l = i - 1;
	ai = a[i];
	h = ai[i];

	if (h != 0.0)

	    for (j = lo; j <= hi; j++) {

		s = si = 0.0;

		for (k = lo; k <= l; k++) {
		    zrk = zr[k];
		    zik = zi[k];
		    f = ai[k];
		    g = a[k][i];
		    s  += f * zrk[j] - g * zik[j];
		    si += f * zik[j] + g * zrk[j];
		}

		s  = (s  / h) / h;
		si = (si / h) / h;

		for (k = lo; k <= l; k++) {
		    zrk = zr[k];
		    zik = zi[k];
		    f = ai[k];
		    g = a[k][i];
		    zrk[j] -=  s  * f + si * g;
		    zik[j] -=  si * f - s  * g;
		}
	   }
    }
}

//----------------------------------------------------------------------------//
