/* truncate.cpp
 *
 * Diagonalises the single ion matrix without mean field or Zeeman terms, generates a new basis from this
 * diagonal matrix and truncates this rotated Hamiltonian to include only lowest lying terms in order to
 * save computation time.
 *
 * Functions:
 *   void truncate_hmltn(&pars, &est, &Hic, &iHic, JHi, JLo)       // Calc. a rotated/truncated Hamiltonian
 *   void truncate_expJ(&pars, &est, &gjmbH, &J, T, lnZ, U, *Jm) // Calc. its expectation values
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008-2011 Duc Le - duc.le@ucl.ac.uk
 */

#include "ic1ion.hpp"
#include "vector.h"          // MatPack vector class
#include <fstream>
#include <ctime>

#ifndef _WINDOWS
#include <unistd.h>
#include <fcntl.h>           // For file control options
#include <sys/mman.h>        // For memory map for truncation routines.
#else
#include <windows.h>
#endif

truncRot g_truncRot;

// --------------------------------------------------------------------------------------------------------------- //
// Truncates the matrix, and stores the single ion part in the estates array
// --------------------------------------------------------------------------------------------------------------- //
void truncate_hmltn(icpars &pars, ComplexMatrix &est, sMat<double> &Hic, sMat<double> &iHic, int JHi, int JLo)
{
   std::cout << "Icalc(): Calculating rotated matrix for truncation." << std::flush;
   clock_t start,end; start = clock();
   int info,Hsz=getdim(pars.n,pars.l);
   complexdouble *Vf; Vf = new complexdouble[Hsz*Hsz]; double *Ef; Ef = new double[Hsz]; 

   // Calculates the eigenvectors and puts it into *est matrix for use by truncate_expJ()
   std::cout << " Starting single ion matrix diagonalisation... " << std::flush;
   info = ic_diag(Hic,iHic,Vf,Ef); if(info!=0) { std::cerr << "truncate_hmltn: Error diagonalising, info==" << info << "\n"; }
   delete[]Ef; 
   for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++) { 
      if(fabs(Vf[ii*Hsz+jj].r)<DBL_EPSILON) Vf[ii*Hsz+jj].r=0.; if(fabs(Vf[ii*Hsz+jj].i)<DBL_EPSILON) Vf[ii*Hsz+jj].i=0.; } 
   std::cout << "Finished.";
   int offset=EST_OFFSET; memcpy(&est[offset][0],Vf,Hsz*Hsz*sizeof(complexdouble)); offset+=Hsz*Hsz;

   // Set up directory to store matrices if the user asks for it.
   char nstr[6]; char filename[255]; char basename[255];
   if(pars.save_matrices) 
   {
      #ifndef _WINDOWS
      struct stat status; int dirstat=0; stat("results/mms/",&status); if(!S_ISDIR(status.st_mode)) dirstat = mkdir("results/mms",0777);
      if(dirstat!=0) { std::cerr << "Icalc(): " << errno << "\n"; exit(EXIT_FAILURE); }
      #else
      DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::Jmat(): Cannot create mms directory\n";
      #endif
      nstr[0] = (pars.l==F?102:100); if(pars.n<10) { nstr[1] = pars.n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = pars.n+38; nstr[3] = 0; }
      strcpy(basename,"results/mms/"); strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;  
   } 
   else 
      strcpy(basename,"/dev/null");
   #define NSTR(K,Q) nstr[1] = K+48; nstr[2] = Q+48; nstr[3] = 0
   #define MSTR(K,Q) nstr[1] = K+48; nstr[2] = 109;  nstr[3] = Q+48; nstr[4] = 0

   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
   sMat<double> zeroes; zeroes.zero(Hsz,Hsz); sMat<double> Upq,Umq; complexdouble *zJmat;
   int cb = (int)(pars.truncate_level*(double)Hsz); complexdouble *Hrot,*zmt; zmt = new complexdouble[Hsz*cb];
   Hrot = new complexdouble[cb*cb]; int memloc=0;
   // Calculates the rotated single ion Hamiltonian
   char notranspose='N',transpose='C',uplo='U',side='L'; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   if(iHic.isempty()) zJmat=zmat2f(Hic,zeroes); else zJmat = zmat2f(Hic,iHic);
   F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,Vf,&Hsz,&zbeta,zmt,&Hsz);
   F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,Vf,&Hsz,zmt,&Hsz,&zbeta,Hrot,&cb); free(zJmat);
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(Hrot[ii*cb+jj].r)<DBL_EPSILON) Hrot[ii*cb+jj].r=0.; if(fabs(Hrot[ii*cb+jj].i)<DBL_EPSILON) Hrot[ii*cb+jj].i=0.; } 
   memcpy(&est[memloc+offset][0],Hrot,cb*cb*sizeof(complexdouble)); memloc+=cb*cb;
   // Calculates the rotated multipolar operators for the mean field terms
   std::cout << " Using " << cb << " levels of " << Hsz << ".\nIcalc(): Starting calculation of rotated mean field operators... " << std::flush;
   icmfmat mfmat(pars.n,pars.l,JHi-JLo+1,pars.save_matrices); double redmat;
   for(int iJ=(JLo-1); iJ<JHi; iJ++)
   {
      if(iJ<6) 
      {
         if(im[iJ]==0) zJmat=zmat2f(mfmat.J[iJ],zeroes); else zJmat = zmat2f(zeroes,mfmat.J[iJ]);
      }
      else
      {
         NSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(pars.n,k[iJ],abs(q[iJ]),pars.l); rmzeros(Upq); mm_gout(Upq,filename); }
         MSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(pars.n,k[iJ],-abs(q[iJ]),pars.l); rmzeros(Umq); mm_gout(Umq,filename); }
         redmat = pow(-1.,(double)abs(pars.l)) * (2*pars.l+1) * threej(2*pars.l,2*k[iJ],2*pars.l,0,0,0);
         #ifdef JIJCONV
         if(pars.B.norm().find("Stevens")!=std::string::npos) redmat*=pars.jijconv[iJ+1];
         #endif
         if(q[iJ]<0) { if((q[iJ]%2)==0) Umq -= Upq; else Umq += Upq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; }
         Umq *= redmat; if(im[iJ]==0) zJmat=zmat2f(Umq,zeroes); else zJmat = zmat2f(zeroes,Umq);
      }
      F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,Vf,&Hsz,&zbeta,zmt,&Hsz);
      F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,Vf,&Hsz,zmt,&Hsz,&zbeta,Hrot,&cb); free(zJmat);
      for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
         if(fabs(Hrot[ii*cb+jj].r)<DBL_EPSILON) Hrot[ii*cb+jj].r=0.; if(fabs(Hrot[ii*cb+jj].i)<DBL_EPSILON) Hrot[ii*cb+jj].i=0.; } 
      memcpy(&est[memloc+offset][0],Hrot,cb*cb*sizeof(complexdouble)); memloc+=cb*cb;
   }
   delete[]Vf; delete[]Hrot; delete[]zmt;
   end = clock(); std::cout << "Done. Time to set up rotated matrices = " << (double)(end-start)/CLOCKS_PER_SEC << "s." << std::endl;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine for opmat() function for cluster module. Returns only req matrix in packed format (Real Upper, i Lower).
// --------------------------------------------------------------------------------------------------------------- //
void truncate_hmltn_packed(icpars &pars, sMat<double> &Mat, sMat<double> &iMat, Matrix &retmat, const char* filename)
{
   int info,Hsz=getdim(pars.n,pars.l);
   int cb = (int)(pars.truncate_level*(double)Hsz); Matrix outmat(1,cb,1,cb); outmat=0;
   complexdouble *Vf; //Vf = new complexdouble[Hsz*Hsz]; 
 //double *Ef; Ef = new double[Hsz]; 

   // Checks whether this sipf file has previously been seen; if not add 
   int iV;
   bool found=false;
   for(iV=0; iV<(int)g_truncRot.sipfs.size(); iV++) { 
      if(strcmp(filename,g_truncRot.sipfs[iV].c_str())==0) { found=true; Vf = g_truncRot.V[iV]; break; }
   }

   // Not previously seen - calculate the rotation matrix and push it into the vector
   if(!found) 
   { 
      g_truncRot.sipfs.push_back(std::string(filename)); 
      sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars); Hic/=MEV2CM; if(!iHic.isempty()) iHic/=MEV2CM;
      double *Ef; Ef = new double[Hsz]; Vf = new complexdouble[Hsz*Hsz]; g_truncRot.V.push_back(Vf);
      info = ic_diag(Hic,iHic,Vf,Ef); if(info!=0) { std::cerr << "truncate_hmltn_packed: Error diagonalising, info==" << info << "\n"; }
      delete[]Ef; 
   }

   // Calculates the rotated single ion Hamiltonian
   sMat<double> zeroes; zeroes.zero(Hsz,Hsz); sMat<double> Upq,Umq; complexdouble *zJmat=zmat2f(Mat,iMat);;
   complexdouble *Hrot,*zmt; zmt = new complexdouble[Hsz*cb]; Hrot = new complexdouble[cb*cb];
   char notranspose='N',transpose='C',uplo='U',side='L'; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,Vf,&Hsz,&zbeta,zmt,&Hsz);
   F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,Vf,&Hsz,zmt,&Hsz,&zbeta,Hrot,&cb); free(zJmat);
   for(int ii=0; ii<cb; ii++) for(int jj=ii; jj<cb; jj++) { 
      if(fabs(Hrot[ii*cb+jj].r)>DBL_EPSILON) outmat(ii+1,jj+1)=Hrot[ii*cb+jj].r; 
      if(fabs(Hrot[ii*cb+jj].i)<DBL_EPSILON) outmat(jj+1,ii+1)=Hrot[ii*cb+jj].i; 
   } 
 //memcpy(&outmat[0][0],Hrot,cb*cb*sizeof(complexdouble)); memloc+=cb*cb;
 //delete[]Vf; 
   delete[]Hrot; delete[]zmt;

   retmat = outmat;
}

// --------------------------------------------------------------------------------------------------------------- //
// Uses the stored eigenvectors of the single ion Hamiltonian to truncate the matrix. Calc. expectation values.
// --------------------------------------------------------------------------------------------------------------- //
void truncate_expJ(icpars &pars, ComplexMatrix &est, Vector &gjmbH, Matrix &J, Vector &T, Vector &lnZ, Vector &U)
{
   int Hsz=getdim(pars.n,pars.l);
   char uplo='U'; complexdouble zme;
   int Esz, incx=1; std::vector<double> E, me, eb;
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;

   int cb = (int)(pars.truncate_level*Hsz), offset = EST_OFFSET+Hsz*Hsz; 
   complexdouble *Hrot; Hrot = new complexdouble[cb*cb]; 
   memcpy(Hrot,&est[offset][0],cb*cb*sizeof(complexdouble)); int szapy=cb*cb; complexdouble a; a.r=1.; a.i=0.;
   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};

   // Calculates the mean field Hamiltonian = H_singleion + sum_i(gjmbH[i]*Operator[i])
   for(int iJ=1; iJ<=(gjmbH.Hi()-gjmbH.Lo()+1); iJ++)
   {
      if (q[iJ]<0) a.r = -gjmbH[iJ+gjmbH.Lo()-1]; else a.r = -gjmbH[iJ+gjmbH.Lo()-1];
      if (fabs(a.r)>DBL_EPSILON) F77NAME(zaxpy)(&szapy,&a,(complexdouble*)&est[iJ*cb*cb+offset][0],&incx,Hrot,&incx);
   }

   // Diagonalises the rotated mean field Hamiltonian
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(Hrot[ii+jj*cb].r)<DBL_EPSILON) Hrot[ii+jj*cb].r=0.; if(fabs(Hrot[ii+jj*cb].i)<DBL_EPSILON) Hrot[ii+jj*cb].i=0.; } 
   iceig VE; VE.calc(cb,Hrot); delete[]Hrot;
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(VE.zV(ii,jj).r)<DBL_EPSILON && fabs(VE.zV(ii,jj).i)<DBL_EPSILON) VE.zV(ii,jj) = 0.; } 

   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<cb; Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T(T.Hi())))<DBL_EPSILON || VE.E(Esz+1)==0) break; }

   // Does initialisations in case we need to recalculate the higher order multipolar matrices (e.g. for spins)
   int memloc=cb*cb, oldJhi=J.Rhi();
   complexdouble *opmat=0, *zmt=0;
   sMat<double> zeroes; zeroes.zero(Hsz,Hsz); sMat<double> Upq,Umq; complexdouble *zJmat;
   char nstr[6]; char filename[255]; char basename[255];
   char notranspose='N',transpose='C',side='L'; 
   if(pars.save_matrices) 
   {
      #ifndef _WINDOWS
      struct stat status; int dirstat=0; stat("results/mms/",&status); if(!S_ISDIR(status.st_mode)) dirstat = mkdir("results/mms",0777);
      if(dirstat!=0) { std::cerr << "ic1ion truncate: Cannot create mms directory\n"; exit(EXIT_FAILURE); }
      #else
      DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::Jmat(): Cannot create mms directory\n";
      #endif
      nstr[0] = (pars.l==F?102:100); if(pars.n<10) { nstr[1] = pars.n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = pars.n+38; nstr[3] = 0; }
      strcpy(basename,"results/mms/"); strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;  
   } 
   else 
      strcpy(basename,"/dev/null");

   // Checks that this time we require expectation values of higher order multipoles even though these were not used in mcphasit
   clock_t start,end; start = clock();
   if((J.Rhi()*cb*cb+Hsz*Hsz)>est.Rows()) 
   {
      std::cerr << "ic1ion truncate: Multipolar operators not precalculated. Calculating now..." << std::flush;
      oldJhi = (est.Rows()-offset)/cb/cb-1;
      zmt = new complexdouble[Hsz*cb];
      opmat = new complexdouble[cb*cb];
   }

   // Calculates the rotated operators for the mean field terms
   complexdouble *zt; Vector Z(1,T.Hi());Z=0.; eb.assign(Esz,0.); U=0;
   for(int iJ=(J.Rlo()-1); iJ<J.Rhi(); iJ++)
   {
      me.assign(Esz,0.);
      zt = (complexdouble*)malloc(cb*sizeof(complexdouble));for(int Ti=1;Ti<=J.Chi();++Ti)J[iJ+1][Ti]=0.; 
      if(iJ>=oldJhi) 
      {
         if(iJ<6) 
         {
            icmfmat mfmat(pars.n,pars.l,6,pars.save_matrices);
            if(im[iJ]==0) zJmat=zmat2f(mfmat.J[iJ],zeroes); else zJmat = zmat2f(zeroes,mfmat.J[iJ]);
         }
         else
         {
            NSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(pars.n,k[iJ],abs(q[iJ]),pars.l); rmzeros(Upq); mm_gout(Upq,filename); }
            MSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(pars.n,k[iJ],-abs(q[iJ]),pars.l); rmzeros(Umq); mm_gout(Umq,filename); }
            double redmat = pow(-1.,(double)abs(pars.l)) * (2*pars.l+1) * threej(2*pars.l,2*k[iJ],2*pars.l,0,0,0);
            #ifdef JIJCONV
            if(pars.B.norm().find("Stevens")!=std::string::npos) redmat*=pars.jijconv[iJ+1];
            #endif
            if(q[iJ]<0) { if((q[iJ]%2)==0) Umq -= Upq; else Umq += Upq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; }
            Umq *= redmat; if(im[iJ]==0) zJmat=zmat2f(Umq,zeroes); else zJmat = zmat2f(zeroes,Umq);
         }
         F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,(complexdouble*)&est[EST_OFFSET][0],&Hsz,&zbeta,zmt,&Hsz);
         F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,(complexdouble*)&est[EST_OFFSET][0],&Hsz,zmt,&Hsz,&zbeta,opmat,&cb); free(zJmat);
         for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
            if(fabs(opmat[ii*cb+jj].r)<DBL_EPSILON) opmat[ii*cb+jj].r=0.; if(fabs(opmat[ii*cb+jj].i)<DBL_EPSILON) opmat[ii*cb+jj].i=0.; } 
      }
      for(int ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         if(iJ>=oldJhi) F77NAME(zhemv)(&uplo, &cb, &zalpha, opmat,                                  &cb, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
         else           F77NAME(zhemv)(&uplo, &cb, &zalpha, (complexdouble*)&est[memloc+offset][0], &cb, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
         #ifdef _G77 
         F77NAME(zdotc)(&zme, &cb, VE.zV(ind_j), &incx, zt, &incx);
         #else
         zme = F77NAME(zdotc)(&cb, VE.zV(ind_j), &incx, zt, &incx);
         #endif
         me[ind_j] = zme.r;for(int Ti=1;Ti<=T.Hi();++Ti){
         if(iJ==(J.Rlo()-1)) { eb[ind_j] = exp(-E[ind_j]/(KB*T(Ti))); Z(Ti)+=eb[ind_j]; U(Ti)+=(E[ind_j]+VE.E(0))*eb[ind_j]; }
         J[iJ+1][Ti]+=me[ind_j]*eb[ind_j];}
      }
      free(zt); for(int Ti=1;Ti<=T.Hi();++Ti){J[iJ+1][Ti]/=Z(Ti); if(iJ==(J.Rlo()-1)) U(Ti)/=Z; }
      memloc+=cb*cb;
   }
   if((J.Rhi()*cb*cb+Hsz*Hsz)>est.Rows()) 
   {
      end = clock(); std::cout << " Done. Elapsed time = " << (double)(end-start)/CLOCKS_PER_SEC << "s." << std::endl;
//    if(!opmat) { delete[]opmat; *opmat=0; } if(!zmt) { delete[]zmt; *zmt=0; }
      delete[]opmat; delete[]zmt;
   }
   for(int Ti=1;Ti<=T.Hi();++Ti)lnZ(Ti) = log(Z(Ti))-VE.E(0)/(KB*T(Ti));
}

// --------------------------------------------------------------------------------------------------------------- //
// Uses the stored eigenvectors of the single ion Hamiltonian to truncate the matrix. Calc. magnetisation density
// --------------------------------------------------------------------------------------------------------------- //
void truncate_spindensity_expJ(icpars &pars, ComplexMatrix &est, Vector &gjmbH, Vector &J, double T, int xyz)
{
   int Hsz=getdim(pars.n,pars.l);
   char uplo='U'; complexdouble zme;
   int Esz, incx=1; std::vector<double> E, me, eb;
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;

   int cb = (int)(pars.truncate_level*Hsz), offset = EST_OFFSET+Hsz*Hsz; 
   complexdouble *Hrot; Hrot = new complexdouble[cb*cb]; 
   memcpy(Hrot,&est[offset][0],cb*cb*sizeof(complexdouble)); int szapy=cb*cb; complexdouble a; a.r=1.; a.i=0.;

   int k[] = {0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

   // Calculates the mean field Hamiltonian = H_singleion + sum_i(gjmbH[i]*Operator[i])
   for(int iJ=1; iJ<=(gjmbH.Hi()-gjmbH.Lo()+1); iJ++)
   {
      if (q[iJ]<0) a.r = -gjmbH[iJ+gjmbH.Lo()-1]; else a.r = -gjmbH[iJ+gjmbH.Lo()-1];
      if (fabs(a.r)>DBL_EPSILON) F77NAME(zaxpy)(&szapy,&a,(complexdouble*)&est[iJ*cb*cb+offset][0],&incx,Hrot,&incx);
   }

   // Diagonalises the rotated mean field Hamiltonian
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(Hrot[ii+jj*cb].r)<DBL_EPSILON) Hrot[ii+jj*cb].r=0.; if(fabs(Hrot[ii+jj*cb].i)<DBL_EPSILON) Hrot[ii+jj*cb].i=0.; } 
   iceig VE; VE.calc(cb,Hrot); delete[]Hrot;
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(VE.zV(ii,jj).r)<DBL_EPSILON && fabs(VE.zV(ii,jj).i)<DBL_EPSILON) VE.zV(ii,jj) = 0.; } 

   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<cb; Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || VE.E(Esz+1)==0) break; }

   complexdouble *opmat=0, *zmt=0;
   zmt = new complexdouble[Hsz*cb];
   opmat = new complexdouble[cb*cb];

   char xyzstr[] = "xyz";
   if(xyz>0) { std::cout << "Calculating the expectation values of the spin density operator S" << xyzstr[xyz-1] << "\n"; }
   else      { std::cout << "Calculating the expectation values of the orbital moment density operator L" << xyzstr[-xyz-1] << "\n"; }

   clock_t start,end; start = clock();
   std::cerr << "ic1ion truncate: Calculating rotated M(Q) matrices, after Balcar and Lovesey..." << std::flush;

   // Calculates the rotated operators for the mean field terms
   complexdouble *zt, *zJmat; double Z=0.; eb.assign(Esz,0.);
   char notranspose='N',transpose='C',side='L'; 
   for(int iJ=0; iJ<49; iJ++)
   {
      me.assign(Esz,0.);
      zt = (complexdouble*)malloc(cb*sizeof(complexdouble)); J[iJ+1]=0.; 

      zJmat = balcar_Mq(xyz,k[iJ],q[iJ],pars.n,pars.l);

      F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,(complexdouble*)&est[EST_OFFSET][0],&Hsz,&zbeta,zmt,&Hsz);
      F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,(complexdouble*)&est[EST_OFFSET][0],&Hsz,zmt,&Hsz,&zbeta,opmat,&cb); free(zJmat);
      for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
         if(fabs(opmat[ii*cb+jj].r)<DBL_EPSILON) opmat[ii*cb+jj].r=0.; if(fabs(opmat[ii*cb+jj].i)<DBL_EPSILON) opmat[ii*cb+jj].i=0.; } 

      for(int ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|M(q)|Vi>
         F77NAME(zhemv)(&uplo, &cb, &zalpha, opmat, &cb, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
         #ifdef _G77 
         F77NAME(zdotc)(&zme, &cb, VE.zV(ind_j), &incx, zt, &incx);
         #else
         zme = F77NAME(zdotc)(&cb, VE.zV(ind_j), &incx, zt, &incx);
         #endif
         me[ind_j] = zme.r;
         if(iJ==0) { eb[ind_j] = exp(-E[ind_j]/(KB*T)); Z+=eb[ind_j]; }
         J[iJ+1]+=me[ind_j]*eb[ind_j];
      }
      free(zt); J[iJ+1]/=Z;
   }

   end = clock(); std::cout << " Done. Elapsed time = " << (double)(end-start)/CLOCKS_PER_SEC << "s." << std::endl;
// if(!opmat) { delete[]opmat; *opmat=0; } if(!zmt) { delete[]zmt; *zmt=0; }
   delete[]opmat; delete[]zmt;
}
