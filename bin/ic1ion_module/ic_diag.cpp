/* ic_diag.cpp
 *
 * Diagonalises the IC Hamiltonian matrix by various methods
 *
 * Functions:
 *    int ic_diag(sMat<double>&Hic, sMat<double>&iH, complexdouble*V, double*E);// Diagonalises complex hermitian Hic+iH
 *    int ic_diag(int n, complexdouble *zm, complexdouble *z, double *eigval);
 *    int ic_diag(sMat<double>&Hic, double *V, double *E);                      // Diagonalises real symmetric Hic 
 *    int ic_leig(sMat<double>&H,sMat<double>&i,complexdouble*V,double*E,int n);// Finds only the n lowest eigenval/vec
 *    int ic_leig(int n, complexdouble *zm, complexdouble *z, double*E, int iu);
 *    int ic_leig(sMat<double>&Hic, double *V, double *E, int n);               // Finds only the n lowest eigenval/vec
 *    int ic_arpackeig(int n, complexdouble*m,complexdouble*z, double*E,int iu);// Finds only lowest eigval/vecs, use ARPACK
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"

/* Physical constants. Taken from NIST Reference on Constants, Units, and Uncertainty,
//     http://physics.nist.gov/cuu/Constants/
const double mu_BJ = 927.400949e-26;   // J/Tesla - Bohr magneton
const double mu_Be = 5.78838263e-2;    // meV/T - Bohr magneton
const double mu_B  = 0.46686437;       // cm^{-1} / Tesla - Bohr magneton
const double k_B   = 1.3806505e-23;    // J/K - Boltzmann constant
const double Q_e   = 1.60217653e-19;   // C - Charge of electron
const double N_A   = 6.0221415e23;     // Avogadro's number
const double h     = 6.62606896e-34;   // Js - Planck's constant
const double c     = 299792458;        // m/s - speed of light in vacuum */
//#define KB     0.08617343183         // meV/K - Boltzmann constant
//#define KBc    0.69503568043         // cm^{-1)/K - Boltzmann constant

// --------------------------------------------------------------------------------------------------------------- //
// Diagonalises the IC Hamiltonian matrix and returns an array of eigenvectors and eigenvalues
// --------------------------------------------------------------------------------------------------------------- //
int ic_diag(sMat<double> &Hic, sMat<double> &iHic, complexdouble *z, double *eigval)
{
   // Paramaters for ZHEEV
   int lda = Hic.nc(), n = Hic.nr(), info = 0;
   char jobz = 'V', uplo = 'U';
   int lwork = 4*n;
   complexdouble *zwork=0, *zm=0;
   {
      char range = 'A'; double vl,vu; int il,iu,numfnd,ldz=n;
      double abstol = 0.00001; int *isuppz = new int[2*n];
      zm = zmat2f(Hic,iHic); zwork = new complexdouble[lwork];
      int lrwork = 24*n,liwork=10*n;
      double *rwork = new double[lrwork];
      int *iwork = new int[liwork];
      F77NAME(zheevr)(&jobz, &range, &uplo, &n, zm, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, eigval, 
              z, &ldz, isuppz, zwork, &lwork, rwork, &lrwork, iwork, &liwork, &info);
      delete []isuppz; delete []rwork; delete []iwork; free(zm); delete []zwork;
   }
   return info;
}
int ic_diag(int n, complexdouble *zm, complexdouble *z, double *eigval)
{
   // Paramaters for ZHEEV
   int lda = n, info = 0;
   char jobz = 'V', uplo = 'U';
   int lwork = 4*n;
   complexdouble *zwork=0;
   {
      char range = 'A'; double vl,vu; int il,iu,numfnd,ldz=n;
      double abstol = 0.00001; int *isuppz = new int[2*n];
      zwork = new complexdouble[lwork];
      int lrwork = 24*n,liwork=10*n;
      double *rwork = new double[lrwork];
      int *iwork = new int[liwork];
      F77NAME(zheevr)(&jobz, &range, &uplo, &n, zm, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, eigval, 
              z, &ldz, isuppz, zwork, &lwork, rwork, &lrwork, iwork, &liwork, &info);
      delete []isuppz; delete []rwork; delete []iwork; delete []zwork;
   }
   return info;
}
int ic_diag(sMat<double> &Hic, double *m, double *eigval)
{
   // Paramaters for DSYEV
   int lda = Hic.nc(), n = Hic.nr(), info = 0;
   char jobz = 'V', uplo = 'U';
   int lwork = 4*n;
   double *work=0, *mz=0;
   {
      char range = 'A'; double vl,vu; int il=1,iu=n/10,numfnd,ldz=n;
      double abstol = 0.00001; int *isuppz = new int[2*n];
      lwork = 26*n; int liwork=10*n;
      mz = Hic.f_array(); work = new double[lwork];
      int *iwork = new int[liwork];
      F77NAME(dsyevr)(&jobz, &range, &uplo, &n, mz, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, eigval, m, 
                      &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
      delete []isuppz; delete[]iwork; free(mz); delete []work;
   }
   return info;
}
int ic_diag(double *mz, int lda, int n, double *m, double *eigval)
{
   // Paramaters for DSYEV
   int info = 0, lwork = 4*n;
   char jobz = 'V', uplo = 'U';
   double *work=0;
   {
      char range = 'A'; double vl,vu; int il=1,iu=n/10,numfnd,ldz=n;
      double abstol = 0.00001; int *isuppz = new int[2*n];
      lwork = 26*n; int liwork=10*n;
      work = new double[lwork];
      int *iwork = new int[liwork];
      F77NAME(dsyevr)(&jobz, &range, &uplo, &n, mz, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, eigval, m, 
                      &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
      delete []isuppz; delete[]iwork; delete []work;
   }
   return info;
}

// --------------------------------------------------------------------------------------------------------------- //
// Diagonalises the IC Hamiltonian matrix and returns an array of eigenvectors and eigenvalues up to specified index
// --------------------------------------------------------------------------------------------------------------- //
int ic_leig(sMat<double> &Hic, sMat<double> &iHic, complexdouble *z, double *eigval, int iu)
{
   // Paramaters for DSYEV
   int lda = Hic.nc(), n = Hic.nr(), info = 0;
   char jobz = 'V', uplo = 'U';
   int lwork = 4*n;
   complexdouble *zwork=0, *zm=0;
   {
      char range = 'I'; double vl,vu; int il=1,numfnd,ldz=n;
      double abstol = 0.00001; int *isuppz = new int[2*n];
      zm = zmat2f(Hic,iHic); zwork = new complexdouble[lwork];
      int lrwork = 24*n,liwork=10*n;
      double *rwork = new double[lrwork];
      int *iwork = new int[liwork];
      F77NAME(zheevr)(&jobz, &range, &uplo, &n, zm, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, eigval, 
              z, &ldz, isuppz, zwork, &lwork, rwork, &lrwork, iwork, &liwork, &info);
      delete []isuppz; delete []rwork; delete []iwork; free(zm); delete []zwork;
   }
   return info;
}
int ic_leig(int n, complexdouble *zm, complexdouble *z, double *eigval, int iu)
{
   // Paramaters for DSYEV
   int lda = n, info = 0;
   char jobz = 'V', uplo = 'U';
   int lwork = 4*n;
   complexdouble *zwork=0;
   {
      char range = 'I'; double vl,vu; int il=1,numfnd,ldz=n;
      double abstol = 0.00001; int *isuppz = new int[2*n];
      zwork = new complexdouble[lwork];
      int lrwork = 24*n,liwork=10*n;
      double *rwork = new double[lrwork];
      int *iwork = new int[liwork];
      F77NAME(zheevr)(&jobz, &range, &uplo, &n, zm, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, eigval, 
              z, &ldz, isuppz, zwork, &lwork, rwork, &lrwork, iwork, &liwork, &info);
      delete []isuppz; delete []rwork; delete []iwork; delete []zwork;
   }
   return info;
}
int ic_leig(sMat<double> &Hic, double *m, double *eigval, int iu)
{
   // Paramaters for DSYEV
   int lda = Hic.nc(), n = Hic.nr(), info = 0;
   char jobz = 'V', uplo = 'U';
   int lwork = 4*n;
   double *work=0, *mz=0;
   {
      char range = 'I'; double vl,vu; int il=1,numfnd,ldz=n;
      double abstol = 0.00001; int *isuppz = new int[2*n];
      lwork = 26*n; int liwork=10*n;
      mz = Hic.f_array(); work = new double[lwork];
      int *iwork = new int[liwork];
      F77NAME(dsyevr)(&jobz, &range, &uplo, &n, mz, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, eigval, m, 
                      &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
      delete []isuppz; delete[]iwork; free(mz);
   }
   return info;
}

#ifndef NO_ARPACK
// --------------------------------------------------------------------------------------------------------------- //
// Function to calculate some of the lowest energy eigenvectors/values of the Hamiltonian using ARPACK
// --------------------------------------------------------------------------------------------------------------- //
int ic_arpackeig(int n, double *zm, double *z, double *eigval, int nev)
{
   // Doesn't work at present - needs debugging!
   int ido=0,ncv=(2*nev>n?n:2*nev),iparam[]={1,0,10000,1,nev,0,1,0,0,0,0};
   //         IPARAM={ISHIFT,,MAXITER,NB,NCONV,,MODE,NP,NUMOP,NUMOPB,NUMREO
   int lworkl=3*ncv*ncv+6*ncv+1,info=0;  // Must be > 3*NCV**2 + 6*NCV
   char bmat='I';
   double tol=1e-12; 

   char uplo = 'U'; int inc = 1; double alpha=1., beta=0.;

   char whichp[]="SR";       // LM/SM==Largest/Smallest Magnitude; LA/SA==Algebraic; BE=balanced (half small/half large)
   int *ipntr=(int*)malloc(11*sizeof(int));
   double *resid=(double*)malloc(n*sizeof(double)), *v=(double*)malloc(n*ncv*sizeof(double));
   double *workd=(double*)malloc(3*n*sizeof(double)), *workl=(double*)malloc(lworkl*sizeof(double));
// char trans='N'; 
   while(ido!=99)
   {
      F77NAME(dnaupd)(&ido, &bmat, &n, whichp, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
      if(ido==1 || ido==-1)
         F77NAME(dsymv)(&uplo, &n, &alpha, zm, &n, &workd[ipntr[0]], &inc, &beta, &workd[ipntr[1]], &inc);
//       F77NAME(dgemv)(&trans, &n, &n, &alpha, zm, &n, &workd[ipntr[0]], &inc, &beta, &workd[ipntr[1]], &inc);
   }
   int rvec=1,*select=(int*)malloc(ncv*sizeof(int)); char howmny='A';
   double *d=(double*)malloc(nev*sizeof(double)), *di=(double*)malloc(nev*sizeof(double)), *workev=(double*)malloc(3*ncv*sizeof(double));
// F77NAME(dseupd)(&rvec, &howmny, select, d, z, &n, &alpha, &bmat, &n, whichp, &nev, &tol, resid, &ncv, v, &n, iparam, 
//                 ipntr, workd, workl, &lworkl, &info);
   F77NAME(dneupd)(&rvec, &howmny, select, d, di, z, &n, &alpha, &alpha, workev, &bmat, &n, whichp, &nev, &tol, resid, &ncv, v, &n, iparam, 
                   ipntr, workd, workl, &lworkl, &info);
   int i=1,j=2,ii; double elem; std::vector<int> ind(nev,0);
   for(ii=0; ii<nev; ii++) { eigval[ii] = d[ii]; ind[ii] = ii; }
   while(i<nev)
   {
      if(eigval[i-1]<=eigval[i]) { i=j; j++; }
      else { elem = eigval[i-1]; eigval[i-1] = eigval[i]; eigval[i] = elem; ii=ind[i-1]; ind[i-1]=ind[i]; ind[i]=ii; i--; if(i==0) i=1; }
   }
   memcpy(v,z,n*nev*sizeof(double)); 
   for(i=0; i<nev; i++) { memcpy(&z[i*n],&v[ind[i]*n+1],n*sizeof(double)); }
   memset(&eigval[nev],0,(n-nev)*sizeof(double)); memset(&z[nev*n],0,(n-nev)*n*sizeof(double));

   free(ipntr); free(resid); free(v); free(workd); free(workl); free(select); free(d); free(di); free(workev);

   return info;
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to calculate some of the lowest energy eigenvectors/values of the Hamiltonian using ARPACK
// --------------------------------------------------------------------------------------------------------------- //
int ic_arpackeig(int n, complexdouble *zm, complexdouble *z, double *eigval, int nev)
{
   // Doesn't work at present - needs debugging!
   int ido=0,ncv=(2*nev>n?n:2*nev),iparam[]={1,0,10000,1,nev,0,1,0,0,0,0};
   //         IPARAM={ISHIFT,,MAXITER,NB,NCONV,,MODE,NP,NUMOP,NUMOPB,NUMREO
   int lworkl=3*ncv*ncv+8*ncv,info=0;  // Must be > 3*NCV**2 + 5*NCV
   char bmat='I';
   double tol=1e-12; 

   char uplo = 'U'; int inc = 1;
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;

   char whichp[]="SR";       // LM/SM==Largest/Smallest Magnitude; LR/SR==RealPart; LI/SI==ImaginaryPart
   double *rwork=(double*)malloc(ncv*sizeof(double)); int *ipntr=(int*)malloc(14*sizeof(int));
   complexdouble *resid=(complexdouble*)malloc(n*sizeof(complexdouble)), *v=(complexdouble*)malloc(n*ncv*sizeof(complexdouble));
   complexdouble *workd=(complexdouble*)malloc(3*n*sizeof(complexdouble)), *workl=(complexdouble*)malloc(lworkl*sizeof(complexdouble));
   while(ido!=99)
   {
      F77NAME(znaupd)(&ido, &bmat, &n, whichp, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
      if(ido==1 || ido==-1)
         F77NAME(zhemv)(&uplo, &n, &zalpha, zm, &n, &workd[ipntr[0]], &inc, &zbeta, &workd[ipntr[1]], &inc);
   }
   int rvec=1,*select=(int*)malloc(ncv*sizeof(int)); char howmny='A';
   complexdouble *d=(complexdouble*)malloc(nev*sizeof(complexdouble)), *workev=(complexdouble*)malloc(2*ncv*sizeof(complexdouble));
   F77NAME(zneupd)(&rvec, &howmny, select, d, z, &n, &zalpha, workev, &bmat, &n, whichp, &nev, &tol, resid, &ncv, v, &n, iparam, 
                   ipntr, workd, workl, &lworkl, rwork, &info);
   int i=1,j=2,ii; double elem; std::vector<int> ind(nev,0);
   for(ii=0; ii<nev; ii++) { eigval[ii] = d[ii].r; ind[ii] = ii; }
   while(i<nev)
   {
      if(eigval[i-1]<=eigval[i]) { i=j; j++; }
      else { elem = eigval[i-1]; eigval[i-1] = eigval[i]; eigval[i] = elem; ii=ind[i-1]; ind[i-1]=ind[i]; ind[i]=ii; i--; if(i==0) i=1; }
   }
   memcpy(v,z,n*nev*sizeof(complexdouble)); 
   for(i=0; i<nev; i++) { memcpy(&z[i*n],&v[ind[i]*n+1],n*sizeof(complexdouble)); }
   memset(&eigval[nev],0,(n-nev)*sizeof(double)); memset(&z[nev*n],0,(n-nev)*n*sizeof(complexdouble));

   if(info==0) { free(ipntr); free(rwork); free(resid); free(v); free(workd); free(workl); free(select); free(d); free(workev); }

   return info;
}
#endif
