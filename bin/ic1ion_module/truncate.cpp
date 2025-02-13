/* truncate.cpp
 *
 * Diagonalises the single ion matrix without mean field or Zeeman terms, generates a new basis from this
 * diagonal matrix and truncates this rotated Hamiltonian to include only lowest lying terms in order to
 * save computation time.
 *
 * Functions:
 *   void truncate_hmltn(&pars, &Hic, &iHic, JHi, JLo)       // Calc. a rotated/truncated Hamiltonian
 *   void truncate_expJ(&pars, &gjmbH, &J, T, lnZ, U, *Jm) // Calc. its expectation values
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008-2011 Duc Le - duc.le@ucl.ac.uk
 */

#include "ic1ion.hpp"
#include "vector.h"          // MatPack vector class
#include "martin.h"
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
// --------------------------------------------------------------------------------------------------------------- //
void ic1ion_module::truncate_hmltn(icpars &pars,  sMat<double> &Hic, sMat<double> &iHic, int JHi, int JLo)
{
   std::cout << "#Icalc(): Calculating rotated matrix for truncation." << std::flush;
   clock_t start,end; start = clock();
   int info,Hsz=getdim(pars.n,pars.l);
   double *Ef; Ef = new double[Hsz]; 
   //complexdouble * Vf; Vf = new complexdouble[Hsz*Hsz];
   // Calculates the eigenvectors and puts it into Vf ... no into mfmat.T[0] matrix for use by truncate_expJ()
   if(mfmat.T[0]==NULL)mfmat.T[0]=new complexdouble[Hsz*Hsz];
   std::cout << " Starting single ion matrix diagonalisation... " << std::flush;
   info = ic_diag(Hic,iHic,mfmat.T[0],Ef); if(info!=0) { std::cerr << "truncate_hmltn: Error diagonalising, info==" << info << "\n"; }
   delete[]Ef; 
   for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++) { 
      if(fabs(mfmat.T[0][ii*Hsz+jj].r)<DBL_EPSILON) {mfmat.T[0][ii*Hsz+jj].r=0.;} if(fabs(mfmat.T[0][ii*Hsz+jj].i)<DBL_EPSILON) {mfmat.T[0][ii*Hsz+jj].i=0.;} } 
   std::cout << "Finished.";
   
   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   //int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   //int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   //int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
   sMat<double> zeroes; zeroes.zero(Hsz,Hsz); sMat<double> Upq,Umq; complexdouble *zJmat;
   int cb = (int)(pars.truncate_level*(double)Hsz); complexdouble *zmt; zmt = new complexdouble[Hsz*cb];
   if(cb<2){std::cerr <<  "Truncate too strong, " << cb << " states are too few - please increase truncate_level.\n";exit(EXIT_FAILURE);}
   if(mfmat.T[1]==NULL)mfmat.T[1]= new complexdouble[cb*cb]; 
   // Calculates the rotated single ion Hamiltonian
   char notranspose='N',transpose='C',uplo='U',side='L'; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   if(iHic.isempty()) zJmat=zmat2f(Hic,zeroes); else zJmat = zmat2f(Hic,iHic);
   F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,mfmat.T[0],&Hsz,&zbeta,zmt,&Hsz);
   F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,mfmat.T[0],&Hsz,zmt,&Hsz,&zbeta,mfmat.T[1],&cb); free(zJmat);
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(mfmat.T[1][ii*cb+jj].r)<DBL_EPSILON) {mfmat.T[1][ii*cb+jj].r=0.;} if(fabs(mfmat.T[1][ii*cb+jj].i)<DBL_EPSILON){ mfmat.T[1][ii*cb+jj].i=0.;} } 
   // Calculates the rotated multipolar operators for the mean field terms
   std::cout << " Using " << cb << " levels of " << Hsz << ".\n#Icalc(): Starting calculation of rotated mean field operators... " << std::flush;
  // icmfmat mfmat(pars.n,pars.l,JHi-JLo+1,pars.save_matrices);
   //double redmat;
   for(int iJ=(JLo-1); iJ<JHi; iJ++)
   {sMat<double>  Jmat=mfmat.op_generate(iJ);
    if(mfmat.iflag[iJ]==0) zJmat=zmat2f(Jmat,zeroes); else zJmat = zmat2f(zeroes,Jmat);
   
     if(mfmat.T[iJ+2]==NULL)mfmat.T[iJ+2]= new complexdouble[cb*cb]; 
     F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,mfmat.T[0],&Hsz,&zbeta,zmt,&Hsz);
     F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,mfmat.T[0],&Hsz,zmt,&Hsz,&zbeta,mfmat.T[iJ+2],&cb); 
     free(zJmat);
 for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
         if(fabs(mfmat.T[iJ+2][ii*cb+jj].r)<DBL_EPSILON) {mfmat.T[iJ+2][ii*cb+jj].r=0.;} 
         if(fabs(mfmat.T[iJ+2][ii*cb+jj].i)<DBL_EPSILON) {mfmat.T[iJ+2][ii*cb+jj].i=0.;} 
         } 
   mfmat.op_free(iJ);
   }
  delete[]zmt;
   end = clock(); std::cout << "Done. Time to set up rotated matrices = " << (double)(end-start)/CLOCKS_PER_SEC << "s." << std::endl;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine for opmat() function . Returns only req matrix in packed format (Real Upper, i Lower).
// --------------------------------------------------------------------------------------------------------------- //
void ic1ion_module::truncate_hmltn_packed(icpars &pars, sMat<double> &Mat, sMat<double> &iMat, Matrix &retmat, const char* filename)
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
void ic1ion_module::truncate_expJ(icpars &pars,  Vector &gjmbH, Matrix &J, Vector &T, Vector &lnZ, Vector &U)
{
   int Hsz=getdim(pars.n,pars.l);
   char uplo='U'; complexdouble zme;
   int Esz, incx=1; std::vector<double> E, me, eb;
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;

   int cb = (int)(pars.truncate_level*Hsz);//, offset = EST_OFFSET+Hsz*Hsz; 
   if(cb<2){std::cerr <<  "Truncate too strong, " << cb << " states are too few - please increase truncate_level.\n";exit(EXIT_FAILURE);}
   complexdouble *Hrot; Hrot = new complexdouble[cb*cb]; 
   memcpy(Hrot,mfmat.T[1],cb*cb*sizeof(complexdouble)); 
   int szapy=cb*cb; complexdouble a; a.r=1.; a.i=0.;
   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   //int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   //int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   //int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};

   // Calculates the mean field Hamiltonian = H_singleion + sum_i(gjmbH[i]*Operator[i])
   for(int iJ=1; iJ<=(gjmbH.Hi()-gjmbH.Lo()+1); iJ++)
   {
      a.r = -gjmbH[iJ+gjmbH.Lo()-1]; 
      if (fabs(a.r)>DBL_EPSILON) F77NAME(zaxpy)(&szapy,&a,mfmat.T[iJ+1],&incx,Hrot,&incx);
   }

   // Diagonalises the rotated mean field Hamiltonian
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(Hrot[ii+jj*cb].r)<DBL_EPSILON) {Hrot[ii+jj*cb].r=0.;} if(fabs(Hrot[ii+jj*cb].i)<DBL_EPSILON) {Hrot[ii+jj*cb].i=0.;} } 
   iceig VE; VE.calc(cb,Hrot); delete[]Hrot;
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(VE.zV(ii,jj).r)<DBL_EPSILON && fabs(VE.zV(ii,jj).i)<DBL_EPSILON) VE.zV(ii,jj) = 0.; } 

   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<cb; Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T(T.Hi())))<DBL_EPSILON || VE.E(Esz+1)==0) break; }

   // Does initialisations in case we need to recalculate the higher order multipolar matrices (e.g. for spins)
   complexdouble  *zmt=0;//*opmat=0,
   sMat<double> zeroes; zeroes.zero(Hsz,Hsz);
   complexdouble *zJmat;
   char notranspose='N',transpose='C',side='L'; 
   
   // Checks that this time we require expectation values of higher order multipoles even though these were not used in mcphasit
   clock_t start,end; start = clock();
   int oldJhi=mfmat.T.size()-2;
    if(J.Rhi()>oldJhi){zmt = new complexdouble[Hsz*cb];
  std::cerr << "#ic1ion truncate: Multipolar operators not precalculated. Calculating now..." << std::flush;}
   // Calculates the rotated operators for the mean field terms
   complexdouble *zt; Vector Z(1,T.Hi());Z=0.; eb.assign(Esz,0.); U=0;
   for(int iJ=(J.Rlo()-1); iJ<J.Rhi(); iJ++)
   {
      me.assign(Esz,0.);
      zt = (complexdouble*)malloc(cb*sizeof(complexdouble));for(int Ti=1;Ti<=J.Chi();++Ti)J[iJ+1][Ti]=0.; 
      if(iJ>=oldJhi) 
      {sMat<double> Jmat=mfmat.op_generate(iJ);
         if(mfmat.T[iJ+2]==NULL)mfmat.T[iJ+2]= new complexdouble[cb*cb];
         if(mfmat.iflag[iJ]==0) zJmat=zmat2f(Jmat,zeroes); else zJmat = zmat2f(zeroes,Jmat);
       mfmat.op_free(iJ);
         F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,mfmat.T[0],&Hsz,&zbeta,zmt,&Hsz);
         F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,mfmat.T[0],&Hsz,zmt,&Hsz,&zbeta,mfmat.T[iJ+2],&cb); free(zJmat);
         for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
            if(fabs(mfmat.T[iJ+2][ii*cb+jj].r)<DBL_EPSILON) {mfmat.T[iJ+2][ii*cb+jj].r=0.;} if(fabs(mfmat.T[iJ+2][ii*cb+jj].i)<DBL_EPSILON) {mfmat.T[iJ+2][ii*cb+jj].i=0.;} } 
      }
      for(int ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
          // my substitute >>>> I believe this is faster because it does not compute imag part zme.i !
         me[ind_j] = expectation_value(cb,mfmat.T[iJ+2],VE.zV(ind_j));
        for(int Ti=1;Ti<=T.Hi();++Ti){
         if(iJ==(J.Rlo()-1)) { eb[ind_j] = exp(-E[ind_j]/(KB*T(Ti))); Z(Ti)+=eb[ind_j]; U(Ti)+=(E[ind_j]+VE.E(0))*eb[ind_j]; }
         J[iJ+1][Ti]+=me[ind_j]*eb[ind_j];}
      }
      free(zt); for(int Ti=1;Ti<=T.Hi();++Ti){J[iJ+1][Ti]/=Z(Ti); if(iJ==(J.Rlo()-1)) U(Ti)/=Z; }
     }
  // if((J.Rhi()*cb*cb+Hsz*Hsz)>Pst.Rows()) 
   if(J.Rhi()>oldJhi){
      end = clock(); std::cout << " Done. Elapsed time = " << (double)(end-start)/CLOCKS_PER_SEC << "s." << std::endl;
      delete[]zmt;
   }
   for(int Ti=1;Ti<=T.Hi();++Ti)lnZ(Ti) = log(Z(Ti))-VE.E(0)/(KB*T(Ti));
}

// --------------------------------------------------------------------------------------------------------------- //
// Uses the stored eigenvectors of the single ion Hamiltonian to truncate the matrix. Calc. magnetisation density
// --------------------------------------------------------------------------------------------------------------- //
void ic1ion_module::truncate_spindensity_expJ(icpars &pars,  Vector &gjmbH, Vector &J, double T, int xyz)
{
   int Hsz=getdim(pars.n,pars.l);
   char uplo='U'; complexdouble zme;
   int Esz, incx=1; std::vector<double> E, me, eb;
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;

   int cb = (int)(pars.truncate_level*Hsz);
   complexdouble *Hrot; Hrot = new complexdouble[cb*cb]; 
     memcpy(Hrot,mfmat.T[1],cb*cb*sizeof(complexdouble)); 
   int szapy=cb*cb; complexdouble a; a.r=1.; a.i=0.;

   int k[] = {0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

   // Calculates the mean field Hamiltonian = H_singleion + sum_i(gjmbH[i]*Operator[i])
   for(int iJ=1; iJ<=(gjmbH.Hi()-gjmbH.Lo()+1); iJ++)
   {
      if (q[iJ]<0) a.r = -gjmbH[iJ+gjmbH.Lo()-1]; else a.r = -gjmbH[iJ+gjmbH.Lo()-1];
      //if (fabs(a.r)>DBL_EPSILON) F77NAME(zaxpy)(&szapy,&a,(complexdouble*)&Pst[iJ*cb*cb+offset][0],&incx,Hrot,&incx);
       if (fabs(a.r)>DBL_EPSILON) F77NAME(zaxpy)(&szapy,&a,mfmat.T[iJ+1],&incx,Hrot,&incx);
  }

   // Diagonalises the rotated mean field Hamiltonian
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(Hrot[ii+jj*cb].r)<DBL_EPSILON) {Hrot[ii+jj*cb].r=0.;} if(fabs(Hrot[ii+jj*cb].i)<DBL_EPSILON) {Hrot[ii+jj*cb].i=0.;} } 
   iceig VE; VE.calc(cb,Hrot); delete[]Hrot;
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(VE.zV(ii,jj).r)<DBL_EPSILON && fabs(VE.zV(ii,jj).i)<DBL_EPSILON) VE.zV(ii,jj) = 0.; } 

   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<cb; Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || VE.E(Esz+1)==0) break; }

   complexdouble *opmat=0, *zmt=0;
   zmt = new complexdouble[Hsz*cb];
   opmat = new complexdouble[cb*cb];

   char xyzstr[] = "xyz";
   if(xyz>0) { std::cout << "#Calculating the expectation values of the spin density operator S" << xyzstr[xyz-1] << "\n"; }
   else      { std::cout << "#Calculating the expectation values of the orbital moment density operator L" << xyzstr[-xyz-1] << "\n"; }

   clock_t start,end; start = clock();
   std::cerr << "#ic1ion truncate: Calculating rotated M(Q) matrices, after Balcar and Lovesey..." << std::flush;

   // Calculates the rotated operators for the mean field terms
   complexdouble *zt, *zJmat; double Z=0.; eb.assign(Esz,0.);
   char notranspose='N',transpose='C',side='L'; 
   for(int iJ=0; iJ<49; iJ++)
   {
      me.assign(Esz,0.);
      zt = (complexdouble*)malloc(cb*sizeof(complexdouble)); J[iJ+1]=0.; 

      zJmat = balcar_Mq(xyz,k[iJ],q[iJ],pars.n,pars.l);

        F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,mfmat.T[0],&Hsz,&zbeta,zmt,&Hsz);
         F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,mfmat.T[0],&Hsz,zmt,&Hsz,&zbeta,opmat,&cb); free(zJmat);
       for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
         if(fabs(opmat[ii*cb+jj].r)<DBL_EPSILON) {opmat[ii*cb+jj].r=0.;} if(fabs(opmat[ii*cb+jj].i)<DBL_EPSILON) {opmat[ii*cb+jj].i=0.;} } 

      for(int ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|M(q)|Vi>
         // my substitute >>>> I believe this is faster because it does not compute imag part zme.i !
         me[ind_j] = expectation_value(cb,opmat,VE.zV(ind_j)); // defined in martin.c

          if(iJ==0) { eb[ind_j] = exp(-E[ind_j]/(KB*T)); Z+=eb[ind_j]; }
         J[iJ+1]+=me[ind_j]*eb[ind_j];
      }
      free(zt); J[iJ+1]/=Z;
   }

   end = clock(); std::cout << " Done. Elapsed time = " << (double)(end-start)/CLOCKS_PER_SEC << "s." << std::endl;
// if(!opmat) { delete[]opmat; *opmat=0; } if(!zmt) { delete[]zmt; *zmt=0; }
   delete[]opmat; delete[]zmt;
}
