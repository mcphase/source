/*----------------------------------------------------------------------------*\
| Computes the eigensytem of the generalized complex               chgeigen.cc |
| hermitean matrix eigen value problem                                         |
| for the complex double precision matrix class of MatPack.                    |
|                                                                              |
| Originally written by Michel Bockstedte                                      |
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

void EigenSystemHermiteanGeneral (Matrix& a, Matrix& b, Vector& e, 
				  Matrix& zr, Matrix& zi, 
				  int sort, int maxiter)
//
//  Driver routine for the generalized hermitan eigenvalue problem:
//
//                     A * z = e * B * z
//
//  The real part of the  complex  hermitean matrices a[lo..hi,lo..hi] and
//  b[lo..hi,lo..hi]  must be stored in  the lower triangle, the imaginary
//  parts must be stored in the strict upper triangle. The eigenvalues are
//  returned in e[lo..hi] in ascending numerical order if the sort flag is 
//  set to  True,  otherwise not  ordered  for sort = False. The real  and 
//  imaginary parts of the eigenvectors are returned in  the columns of zr
//  and zi. C.f. comments on Chreduce(), Chreback() and the 
//  EigenSystemHermitean()  routine for the complex Hermitean problem 
//  (file cheigen.c)
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
    int lo = a.Clo();	
    int hi = a.Chi();
    
    // columns and rows must have the same range
    if (a.Rlo() != lo || a.Rhi() != hi) 
	Matpack.Error(Mat::UnspecifiedError,"EigenSystem: matrix must be square"); 

    // the matrices and vectors must be conformant 
    if ( e.Lo() != lo || e.Hi() != hi
      || b.Clo()  != lo || b.Chi()  != hi || b.Rlo()  != lo || b.Rhi()  != hi
      || zr.Clo() != lo || zr.Chi() != hi || zr.Rlo() != lo || zr.Rhi() != hi
      || zi.Clo() != lo || zi.Chi() != hi || zi.Rlo() != lo || zi.Rhi() != hi)
	Matpack.Error(Mat::UnspecifiedError,"EigenSystem: vectors and matrices are not conformant");

    // allocate auxilliary matrices
    Matrix l(lo,hi,lo,hi);
    Matrix p(lo,hi,lo,hi);

    // transform to standard eigenvalue problem using Cholesky factorization 
    Chreduce(a,b,l,p);

    // compute the eigensystem of the transformed complex hermitean matrix p
    EigenSystemHermitean(p,e,zr,zi,sort,maxiter);
    
    // transform back to generalized problem
    Chreback(l,zr,zi);
}


//----------------------------------------------------------------------------//


void Chreduce (Matrix& A, Matrix& B, Matrix& L, Matrix& P)
//
//  This routine is the complex  counterpart to  the EISPACK routine reduc1
//  It reduces the generalized hermitian eigenvalue problem to the standard 
//  hermitian eigenvalue problem. 
//
//	  A * x = e * B * x    to     P * y = e * y
//   
//  using the Cholesky factorization. The real part of all matrices must be
//  stored in the lower triangle, the imaginary parts must be stored in the
//  strict upper triangle, obeying:
//
//  	a[i][j] = real(a[i][j]) and a[j][i] = imag(a[j][i]) where i >= j.
//   
//  References:
//
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
//
{
    int i,j,k;
    double Rex,Imx,y = 0;

    // get the dimension information
    int lo = A.Clo();	
    int hi = A.Chi();

    // store diagonal element square roots
    Vector Imp(lo,hi);
    
    // avoid call to index operator that optimizes very badely
    double **a = &A[0],
	   **b = &B[0],
	   **l = &L[0], 
           **p = &P[0],
           *imp= &Imp[0];

    for (i = lo; i <= hi; i++)
	for (j = i; j <= hi; j++) {
	    Rex = b[j][i];
	    Imx = b[i][j];
            if (i == j) Imx = 0.0;
	    for (k = i-1; k >= lo; k--) {
		Rex -= l[i][k] * l[j][k] + l[k][i] * l[k][j];
		Imx -= l[k][i] * l[j][k] - l[i][k] * l[k][j];
	    }
	    if (i == j) {
		if (Rex <= 0.0) 
		    Matpack.Error(Mat::UnspecifiedError,"Chreduce: matrix B is not positiv definite - possible reason: magnetic structure in mcdisp.mf  metastable.\n");
		else { 
		    y = sqrt(Rex);		    
		    l[i][i] = y;
		}
	    } else {
		l[j][i] = Rex / y;
		l[i][j] = -Imx / y;
	    }    
	}
    
    for (i = lo; i <= hi; i++)  {
	y = l[i][i];
	for (j = i; j <= hi; j++) { 
	    Rex = a[j][i];
	    Imx = a[i][j];
	    if (i == j) Imx = 0.0;
	    for (k = i-1; k >= lo; k-- ) {
		Rex -= l[i][k] * p[j][k] - l[k][i] * p[k][j];
		Imx -= l[i][k] * p[k][j] + l[k][i] * p[j][k];
	    }
	    p[j][i] = Rex / y;
	    if (i != j) 
	        p[i][j] = Imx / y;
	    else 
	        imp[i] = Imx / y;
	}
    }
    
    for (j = lo; j <= hi; j++) 
        for (i = j; i <= hi; i++) {
	    Rex = p[i][j];
	    Imx = p[j][i];
	    if (i == j) Imx = imp[i];
	    for (k = i-1; k >= j; k--) {   		  
		if (j != k) {
		    Rex -= p[k][j] * l[i][k] + p[j][k] * l[k][i];
		    Imx -= p[j][k] * l[i][k] - p[k][j] * l[k][i];
	        } else {
		    Rex -= p[k][j] * l[i][k] + imp[j] * l[k][i];
		    Imx -= imp[j] * l[i][k] - p[k][j] * l[k][i];
		}
	    }
	    for (k = j-1 ; k >= lo; k--) {
		Rex -= p[j][k] * l[i][k] - p[k][j] * l[k][i];
		Imx += p[k][j] * l[i][k] + p[j][k] * l[k][i];
	    }	
	    p[i][j] = Rex / l[i][i];      
	    if (i != j) 
	        p[j][i] = Imx / l[i][i];
	    else 
	        imp[i] = Imx / l[i][i];
	}	  
}

//----------------------------------------------------------------------------//

void Chreback (Matrix& LL, Matrix& Re_z, Matrix& Im_z)
//
//  This routine is the complex counterpart to the EISPACK  routine  rebaka
//  Given the  Cholesky factorization of b, it transforms  the eigenvectors
//  of the standard eigenvalue problem to eigenvectors of the corresponding 
//  generalized problem.
//
//  References:
//
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
//
{
    // get the dimension information
    int lo = LL.Clo();	
    int hi = LL.Chi();

    // avoid call to index operator that optimizes very badely
    double **Rez = &Re_z[0],
	   **Imz = &Im_z[0],
           **L   = &LL[0];

    for (int j = lo; j <= hi; j++)
	for (int i = hi; i >= lo; i--) {
	    double Rex = Rez[i][j];
	    double Imx = Imz[i][j];
	    for (int k = i+1; k <= hi; k++) {
		Rex -= L[k][i] * Rez[k][j] - L[i][k] * Imz[k][j];
		Imx -= L[k][i] * Imz[k][j] + L[i][k] * Rez[k][j];
	    } 
	    Rez[i][j] = Rex / L[i][i];
	    Imz[i][j] = Imx / L[i][i];
	}
}

//----------------------------------------------------------------------------//
