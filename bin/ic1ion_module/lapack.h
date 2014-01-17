/* -*-C-*- 

// Copyright (C) 2004 
// Christian Stimming <stimming@tuhh.de>

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2, or (at
// your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307,
// USA.

//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
*/

#ifndef LAPACK_H
#define LAPACK_H

/*  Linkage names between C, C++, and Fortran (platform dependent) */

/** @file
 * @brief Platform-dependent macro definitions
 */

//typedef struct { double r,i; } complexdouble;
struct complexdouble {
   double r, i;
   complexdouble operator=(const double v);
};

#if  defined(RIOS) && !defined(CLAPACK)
# define F77NAME(x) x
#else
#if  defined MWLIBS
# define F77NAME(x) x
#else
# define F77NAME(x) x##_
#endif
#endif

#if defined(SGI) && !defined(SGI_DEC)
# define SGI_DEC
# ifdef __cplusplus
extern "C" {
# endif
	void mkidxname() {}
	void mkdatname() {}
# ifdef __cplusplus
}
# endif
#endif

/* Needed for windows DLLs */
#ifndef DLLIMPORT
#  if defined( __declspec ) | defined ( _MSC_VER )
/*     _MSC_VER checks for Microsoft Visual C++. */
/*      Microsoft Visual C++ 7.1  _MSC_VER = 1310 */
/*      Microsoft Visual C++ 7.0  _MSC_VER = 1300 */
/*      Microsoft Visual C++ 6.0  _MSC_VER = 1200 */
/*      Microsoft Visual C++ 5.0  _MSC_VER = 1100 */
#    if BUILDING_LAPACK_DLL
#      define DLLIMPORT __declspec (dllexport)
#    else     /* Not BUILDING_LAPACK_DLL */
#      define DLLIMPORT __declspec (dllimport)
#    endif    /* Not BUILDING_LAPACK_DLL */
#  else
#    define DLLIMPORT
#  endif    /* __declspec */
#endif  /* DLLIMPORT */

extern "C"
{
  // BLAS Utilities
  void F77NAME(dgemv)(char *trans, int *n, int *m, double *am, double *Jmat, int *incx, double *V,
                      int *nev, double *bm, double *vt, int *incy);
  void F77NAME(dsymv)(char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, 
                      double *beta, double *y, int *incy);
  void F77NAME(dtrmv)(char *uplo, char *trans, char *diag, int *n, double *a, int *lda, double *x, int *incx);
  void F77NAME(dsymm)(char *side, char *uplo, int *m, int *n, double *alpha, double *A, int *lda,
                      double *B, int *ldb, double *beta, double *C, int *ldc);
  void F77NAME(dgemm)(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *A,
                      int *lda, double *B, int *ldb, double *beta, double *c, int *ldc);
  void F77NAME(zhemm)(char *side, char *uplo, int *m, int *n, complexdouble *alpha, complexdouble *A, int *lda,
                      complexdouble *B, int *ldb, complexdouble *beta, complexdouble *C, int *ldc);
  void F77NAME(dcopy)(int *n, double *x, int *incx, double *y, int *incy);
  void F77NAME(zgemv)(char *trans, int *n, int *m, complexdouble *alpha, complexdouble *a,
                      int *lda, complexdouble *x, int *incx, complexdouble *beta,
                      complexdouble *y, int *incy);
  void F77NAME(zhemv)(char *uplo, int *n, complexdouble *alpha, complexdouble *a,
                      int *lda, complexdouble *x, int *incx, complexdouble *beta,
                      complexdouble *y, int *incy);
  void F77NAME(zgemm)(char *transa, char *transb, int *m, int *n, int *k, complexdouble *alpha,
                      complexdouble *A, int *lda, complexdouble *B, int *ldb, complexdouble *beta,
                      complexdouble *C, int *ldc); 
  void F77NAME(daxpy)(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
  void F77NAME(zaxpy)(int *n, complexdouble *za, complexdouble *zx, int *incx, complexdouble *zy, int *incy);
#ifdef _G77 
  void F77NAME(ddot)(double retval, int *n, double *x, int *incx, double *y, int *incy);
  void F77NAME(zdotc)(complexdouble *retval, int *n, complexdouble *zx, int *incx, complexdouble *zy, int* incy);
  void F77NAME(zdotu)(complexdouble *retval, int *n, complexdouble *zx, int *incx, complexdouble *zy, int* incy);
#else
  double F77NAME(ddot)(int *n, double *x, int *incx, double *y, int *incy);
  complexdouble F77NAME(zdotc)(int *n, complexdouble *zx, int *incx, complexdouble *zy, int* incy);
  complexdouble F77NAME(zdotu)(int *n, complexdouble *zx, int *incx, complexdouble *zy, int* incy);
#endif
  double F77NAME(dasum)(int *n, double *zx, int *incx);
  double F77NAME(dzasum)(int *n, complexdouble *zx, int *incx);
  int F77NAME(izamax)(int *n, complexdouble *zx, int *incx);

  // LAPACK Utilities
//int F77NAME(ilaenv)(int *ispec, const char *name, char *opts, 
//                    int *n1, int *n2, int *n3, int *n4);

  // Single precision real routines.

  // Double precision real routines.
  void F77NAME(dsyev)(char *jobz, char *uplo, int *N, double *S, int *lda, 
                      double *eig, double *work, int *lwork, int *info);
  void F77NAME(dgesvd)(char *jobu, char *jobvt, int *m, int *n, double *a, 
                       int *lda, double *sing, double *u, int *ldu, double *vt, 
                       int *ldvt, double *work, int *lwork, int *info);
  void F77NAME(dgeqrf)(int *m, int *n, double *a, int *lda, double *tau,
                       double *work, int *lwork, int *info);
  void F77NAME(dorgqr)(int *m, int *n, int *k, double *a, int *lda, double *tau,
                       double *work, int *lwork, int *info);
  void F77NAME(dsyevr)(char *jobz, char *range, char *uplo, int *n, double *mz, int *lda, 
                      double *vl, double *vu, int *il, int *iu, double *abstol, int *numfnd,
                      double *eigval, double *z, int *ldz, int *isuppz, double *work, 
		      int *lwork, int *iwork, int *liwork, int *info);
  void F77NAME(dstegr)(char *jobz, char *range, int *n, double *d, double *e, // Tridiagonal RRR eigenproblem
                      double *vl, double *vu, int *il, int *iu, double *abstol, int *numfnd, double *eigval,
                      double *z, int *ldz, int *isuppz, double *work,
                      int *lwork, int *iwork, int *liwork, int *info);

  // Double precision complex routines.
  void F77NAME(zheev)(char *jobz, char *uplo, int *N, complexdouble *S, int *lda, 
                      double *eig, complexdouble *work, int *lwork, double *rwork, int *info);
  void F77NAME(zheevr)(char *jobz, char *range, char *uplo, int *n, complexdouble *zm, int *lda, 
                      double *vl, double *vu, int *il, int *iu, double *abstol, int *numfnd, double *eigval, 
                      complexdouble *z, int *ldz, int *isuppz, complexdouble *zwork, 
		      int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

  void F77NAME(zhegvd)(int *itype, char *jobz, char *uplo, int *n, // Divide+Conquer generalised eig
                      complexdouble *za, int *lda, complexdouble *zb, int *ldb, double *e, complexdouble *zwork, 
                      int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

#ifndef NO_ARPACK
  // Double precision real ARPACK routines.
  void F77NAME(dnaupd)(int *ido, char *bmat, int *Hsz, char *whichp, int *nev, double *tol,
                      double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr,
                      double *workd, double *workl, int *lworkl, int *info);
  void F77NAME(dneupd)(int *rvec, char *howmny, int *select, double *dr, double *di, double *z, int *ldz, 
                      double *sigmar, double *sigmai, double *workev, char *bmat, int *n, char *whichp,
                      int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam,
                      int *ipntr, double *workd, double *workl, int *lworkl, int *info);
  void F77NAME(dsaupd)(int *ido, char *bmat, int *Hsz, char *whichp, int *nev, double *tol,
                      double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr,
                      double *workd, double *workl, int *lworkl, int *info);
  void F77NAME(dseupd)(int *rvec, char *howmny, int *select, double *d, double *z, int *ldz, 
                      double *sigma, char *bmat, int *n, char *whichp, int *nev, double *tol,
                      double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr,
                      double *workd, double *workl, int *lworkl, int *info);

  // Complex ARPACK routines
  void F77NAME(znaupd)(int *ido, char *bmat, int *n, char *whichp, int *nev, double *tol, complexdouble *resid, 
                      int *ncv, complexdouble *v, int *ldv, int *iparam, int *ipntr, complexdouble *workd, 
                      complexdouble *workl, int *lworkl, double *rwork, int *info);
  void F77NAME(zneupd)(int *rvec, char *howmny, int *select, complexdouble *d, complexdouble *z, int *ldz,
                      complexdouble *sigma, complexdouble *workev, char *bmat, int *n, char *whichp, int *nev,
                      double *tol, complexdouble *resid, int *ncv, complexdouble *v, int *ldv, int *iparam,
                      int *ipntr, complexdouble *workd, complexdouble *workl, int *lworkl, double *rwork, int *info);
  void F77NAME(zsortc)(char *which, int *apply, int *n, complexdouble *x, complexdouble *y);
#endif

}

#endif
