/* spectre.cpp
 *
 * Routines to calculate the IC Hamiltonian using the method in the Spectra/XTAL programs of Hannah Crosswhite.
 * References: https://chmwls.chm.anl.gov/downloads/spectra/about/matrix.html
 *             Carnall, Crosswhite, Crosswhite, Conway, J. Chem. Phys., v64, p3582, 1976
 *
 * Functions:
 * void spectre_hmltn(pars, est, parvalsize)                         // Calculates the IC Hamiltonian
 * vector<double> spectre_expJ(pars,est,parvalsize,gjmbH,Jhi,Jlo,T)  // Calculates the expectation values
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2009 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
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

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the IC Hamiltonian so that's basis states are diagonal in the free-ion part (Hc+Hso)
// --------------------------------------------------------------------------------------------------------------- //
void spectre_hmltn(icpars &pars, ComplexMatrix &est, int parvalsize)
{
   int i,Hsz=getdim(pars.n,pars.l);

   complexdouble *Vcso, *Hrot, *zmt, zalpha, zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0; double *Ecso;
   sMat<double> Hcso = ic_Hcso(pars); Hcso/=MEV2CM; int Hcso_sz = Hcso.nc();          // Hf = full(ic_freehmltn(conf,racah_FtoF_k(F),xi,alpha));
   if (pars.spectrelevels > Hcso.nr()) pars.spectrelevels=Hcso.nr(); 
   Ecso = new double[Hcso_sz]; Vcso = new complexdouble[Hcso_sz*Hcso_sz];
   sMat<double> zeros(Hcso_sz,Hcso_sz); int j = ic_diag(Hcso,zeros,Vcso,Ecso);        // [Vcso,Ecso] = eig(Hf);
   if(j!=0) { std::cerr << "error: spectre\n"; return; }
   for(i=0; i<Hcso_sz; i++) for(j=0; j<Hcso_sz; j++) 
      if(fabs(Vcso[i+j*Hcso_sz].r)<DBL_EPSILON*50) Vcso[i+j*Hcso_sz].r=0.;
   char notranspose = 'N', transpose = 'C'; int cb = pars.spectrelevels; 
   complexdouble *zHcso = zmat2f(Hcso,zeros);
   zmt = new complexdouble[Hsz*cb]; Hrot = new complexdouble[cb*cb];
   F77NAME(zgemm)(&notranspose,&notranspose,&Hcso_sz,&cb,&Hcso_sz,&zalpha,zHcso,      // Hrot = Vcso(:,1:lvl)'*Hcso*Vcso(:,1:lvl);
                     &Hcso_sz,Vcso,&Hcso_sz,&zbeta,zmt,&Hcso_sz);
   F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hcso_sz,&zalpha,Vcso,&Hcso_sz,zmt,
                     &Hcso_sz,&zbeta,Hrot,&cb); free(zHcso); delete[]zmt;
   std::vector<int> J2r(pars.spectrelevels,0), Jz2r, id; 
   int imax,icv=0,mj; fconf conf(pars.n,0,pars.l); int incx = 1;
   std::vector< std::vector<int> > cvSO2CF; std::vector<int> cvEl;
   for(i=0; i<pars.spectrelevels; i++)                                                // for i=1:lvl; 
   {
      // IZAMAX returns fortran-style index (starts at 1)
      imax = F77NAME(izamax)(&Hcso_sz,&Vcso[i*Hcso_sz],&incx); 
      J2r[i] = conf.states[imax-1].J2;                                                //    J(i) = st_SO{find(abs(Vcso(:,i))==max(abs(Vcso(:,i))))}{5}; 
      cvEl.push_back(icv); icv += J2r[i]+1; cvEl.push_back(icv-1);                    //    cvSO2CF(i,1) = icv+1; icv=icv+(2*J(i)+1);
      cvSO2CF.push_back(cvEl); cvEl.clear();                                          //    cvSO2CF(i,2) = icv;
      for(mj=-J2r[i]; mj<=J2r[i]; mj+=2) { Jz2r.push_back(mj); id.push_back(i); }     //    for mj = 1:(2*J(i)+1); id(mj+icv) = i; Jz(mj+icv) = mj-J(i)-1; end
   }                                                                                  // end
   complexdouble *Hrmj = new complexdouble[icv*icv]; 
   convH2H(Hrot,Hrmj,cb,icv,cvSO2CF); delete[]Hrot;                                   // Hrot = full(convH2H(Hf,icv,cvSO2CF));
   // Calculates the CF rotated matrix, and expands to |LSJ,mJ>
   int k,iq,qp,qm; sMat<double> UJ; complexdouble *fUJ, *UJr; 
   double em,ep,epe,eme,p = 1./( pow(-1.,(double)abs(pars.l))*(2.*pars.l+1.) );
   double icfact[] = {0,0,p/threej(2*pars.l,4,2*pars.l,0,0,0),0,p/threej(2*pars.l,8,2*pars.l,0,0,0),0,p/threej(2*pars.l,12,2*pars.l,0,0,0)};
   for(k=2; k<=6; k+=2)
   {
      icfact[k] *= MEV2CM;
      UJ = racah_uJ(pars.n,k,pars.l); fUJ = zmat2f(UJ,zeros);                         // [U4,H,HZ,st_SO] = nshell_rme84(conf);
      UJr = new complexdouble[cb*cb]; zmt = new complexdouble[Hsz*cb];
      F77NAME(zgemm)(&notranspose,&notranspose,&Hcso_sz,&cb,&Hcso_sz,&zalpha,fUJ,     // for l = 1:3; UCOEFF{l} = (RVEC(:,NKPT)'*U4{l}*RVEC(:,NKPT)) ./ icfact(l); end
                       &Hcso_sz,Vcso,&Hcso_sz,&zbeta,zmt,&Hcso_sz);
      F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hcso_sz,&zalpha,Vcso,&Hcso_sz,zmt,
                       &Hcso_sz,&zbeta,UJr,&cb); free(fUJ); delete[]zmt;              // for k = 2:2:6
/*    if(pars.B(k,0)!=0)                                                              //    for iq = find(B{k/2})
      {                                                                               //        Uj = zeros(icv); q = iq-1-k;
         for(i=0; i<icv; i++) for(j=0; j<icv; j++)                                    //        for i = 1:icv; for j = 1:icv
         {
            Hrmj[i+j*icv].r += pow(-1.,(J2r[id[i]]-Jz2r[i])/2.)                       //           Uj(i,j) = (-1)^(J(id(i))-Jz(i)-q)  
                               * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],0,Jz2r[j]) //                       * threej([J(id(i)) k J(id(j)); -Jz(i) q Jz(j)])
                               * UJr[id[i]+id[j]*cb].r * pars.B(k,0)/icfact[k];       //                       * UCOEFF{k/2}(id(i),id(j));
         }                                                                            //        end; end
      }                                                                               //        Hrmj = Hrmj + Uj.*B{k/2}(iq);      % Note have to give -q and +q pars!
*/    for(iq=0; iq<=k; iq++) // if(pars.B(k,iq-k)!=0 || pars.B(k,-(iq-k))!=0)         //    end
      {                                                                               // end
         qm = iq-k; qp = -qm;
         if(qm==0)
         {
            char fh[] = "results/ic1ion.uX0"; fh[16]=k+48; std::fstream FH; FH.open(fh, std::fstream::out); FH.precision(24);
            for(i=0; i<icv; i++) for(j=0; j<icv; j++)
            {
               em = ( pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],2*qm,Jz2r[j]) * UJr[id[i]+id[j]*cb].r ) /icfact[k];
               Hrmj[i+j*icv].r += em * pars.B(k,qp); FH << em << "\n";
            }
            FH.close();
         }
         else
         {
            char fhm[] = "results/ic1ion.uXmX"; fhm[16]=k+48; fhm[18]=qp+48; std::fstream FM; FM.open(fhm, std::fstream::out); FM.precision(24);
            char fhp[] = "results/ic1ion.uXX";  fhp[16]=k+48; fhp[17]=qp+48; std::fstream FP; FP.open(fhp, std::fstream::out); FP.precision(24);
            for(i=0; i<icv; i++) for(j=0; j<icv; j++)
            {
               ep = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],2*qp,Jz2r[j]) * UJr[id[i]+id[j]*cb].r;
               em = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],2*qm,Jz2r[j]) * UJr[id[i]+id[j]*cb].r;
               epe = (ep + em*pow(-1.,qp)) / icfact[k]; eme = (ep - em*pow(-1.,qm)) / icfact[k];
               if(pars.n>(2*pars.l+1))
               {
                  if(fabs(eme)>DBL_EPSILON) { FM << -eme << "\n"; if(fabs(pars.B(k,qm))>DBL_EPSILON) Hrmj[i+j*icv].i -= eme * pars.B(k,qm); } else FM << 0 << "\n";
                  if(fabs(eme)>DBL_EPSILON) { FP << -epe << "\n"; if(fabs(pars.B(k,qp))>DBL_EPSILON) Hrmj[i+j*icv].r -= epe * pars.B(k,qp); } else FP << 0 << "\n";
               }
               else
               {
                  if(fabs(eme)>DBL_EPSILON) { FM <<  eme << "\n"; if(fabs(pars.B(k,qm))>DBL_EPSILON) Hrmj[i+j*icv].i += eme * pars.B(k,qm); } else FM << 0 << "\n";
                  if(fabs(eme)>DBL_EPSILON) { FP <<  epe << "\n"; if(fabs(pars.B(k,qp))>DBL_EPSILON) Hrmj[i+j*icv].r += epe * pars.B(k,qp); } else FP << 0 << "\n";
               }
            }
            FM.close(); FP.close();
         }
      }
   }
   est[0][0] = complex<double> (-0.1,-0.1); est[0][parvalsize/2+1] = complex<double> (icv,0.);
   for(i=1; i<=icv; i++) memcpy(&est[i][1],&Hrmj[(i-1)*icv],icv*sizeof(complexdouble)); delete[]Ecso; delete[]UJr; delete[]Hrmj;
   // Calculates the matrix elements of the magnetic operators L and S
   complexdouble *Lrm,*Srm,*Lrj,*Srj; Lrm = new complexdouble[Hcso_sz*Hcso_sz]; Srm = new complexdouble[Hcso_sz*Hcso_sz]; int S2,L2,S2p,L2p,J2,J2p;
   for(i=0; i<Hcso_sz; i++) for(j=0; j<Hcso_sz; j++)
   {
      S2 = conf.states[i].S2; S2p = conf.states[j].S2; L2 = abs(conf.states[i].L)*2; L2p = abs(conf.states[j].L)*2; J2 = conf.states[i].J2; J2p = conf.states[j].J2;
      Lrm[i+j*Hcso_sz].i = 0.; Srm[i+j*Hcso_sz].i = 0.; Lrm[i+j*Hcso_sz].r = 0.; Srm[i+j*Hcso_sz].r = 0.;
      if(abs(J2-J2p)>2 || S2!=S2p || L2!=L2p || conf.states[i].v!=conf.states[j].v || conf.states[i].U!=conf.states[j].U) continue;
      Lrm[i+j*Hcso_sz].r = pow(-1.,(S2p+L2p+J2+2.)/2.)  * sixj(L2,J2,S2p,J2p,L2p,2) * sqrt((L2p+1.)*(J2+1.)*(J2p+1.)*((L2p/2.)*(L2p/2.+1)));
      Srm[i+j*Hcso_sz].r = pow(-1.,(S2p+L2p+J2p+2.)/2.) * sixj(S2,J2,L2p,J2p,S2p,2) * sqrt((S2p+1.)*(J2+1.)*(J2p+1.)*((S2p/2.)*(S2p/2.+1)));
   }
   zmt = new complexdouble[Hsz*cb]; Lrj = new complexdouble[cb*cb]; Srj = new complexdouble[cb*cb];
   F77NAME(zgemm)(&transpose,&notranspose,&Hcso_sz,&cb,&Hcso_sz,&zalpha,Lrm,&Hcso_sz,Vcso,&Hcso_sz,&zbeta,zmt,&Hcso_sz);
   F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hcso_sz,&zalpha,Vcso,&Hcso_sz,zmt,&Hcso_sz,&zbeta,Lrj,&cb); delete[]Lrm;
   F77NAME(zgemm)(&transpose,&notranspose,&Hcso_sz,&cb,&Hcso_sz,&zalpha,Srm,&Hcso_sz,Vcso,&Hcso_sz,&zbeta,zmt,&Hcso_sz);
   F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hcso_sz,&zalpha,Vcso,&Hcso_sz,zmt,&Hcso_sz,&zbeta,Srj,&cb); delete[]Srm; delete[]zmt;

   double sqrt2 = sqrt(2.),elem;
   char fSx[] = "results/ic1ion.m1"; std::fstream FSx; FSx.open(fSx, std::fstream::out); FSx.precision(24);
   char fLx[] = "results/ic1ion.m2"; std::fstream FLx; FLx.open(fLx, std::fstream::out); FLx.precision(24);
   char fSy[] = "results/ic1ion.m3"; std::fstream FSy; FSy.open(fSy, std::fstream::out); FSy.precision(24);
   char fLy[] = "results/ic1ion.m4"; std::fstream FLy; FLy.open(fLy, std::fstream::out); FLy.precision(24);
   char fSz[] = "results/ic1ion.m5"; std::fstream FSz; FSz.open(fSz, std::fstream::out); FSz.precision(24);
   char fLz[] = "results/ic1ion.m6"; std::fstream FLz; FLz.open(fLz, std::fstream::out); FLz.precision(24);
   for(i=0; i<icv; i++) for(j=0; j<icv; j++)
   {
      ep = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2,J2r[id[j]],-Jz2r[i], 2,Jz2r[j])/sqrt2;
      em = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2,J2r[id[j]],-Jz2r[i],-2,Jz2r[j])/sqrt2;
      elem = (em-ep)*Lrj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FLx << elem << "\n"; else FLx << 0 << "\n";
      elem = (em+ep)*Lrj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FLy << elem << "\n"; else FLy << 0 << "\n";
      elem = (em-ep)*Srj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FSx << elem << "\n"; else FSx << 0 << "\n";
      elem = (em+ep)*Srj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FSy << elem << "\n"; else FSy << 0 << "\n";
      ep = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2,J2r[id[j]],-Jz2r[i], 0,Jz2r[j]);
      elem = ep*Lrj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FLz << elem << "\n"; else FLz << 0 << "\n";
      elem = ep*Srj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FSz << elem << "\n"; else FSz << 0 << "\n";
   }
   FSx.close(); FLx.close(); FSy.close(); FLy.close(); FSz.close(); FLz.close(); delete[]Lrj; delete[]Srj; delete[]Vcso;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the expectation values of the multipolar operators for the spectre-like rotated Hamiltonian
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> spectre_expJ(icpars &pars, ComplexMatrix &est, int parvalsize, Vector &gjmbH, int Jlo, int Jhi, double T)
{
   std::vector<double> J(Jhi+3,0.);
   // Reloads the magnetic operator matrix in the rotated basis
   int i,icv = (int)real(est[0][parvalsize/2+1]); char fn[] = "results/ic1ion.mq"; complexdouble *Jm = new complexdouble[icv*icv]; double elem;
 /*int Jhi = (J.Hi()>6) ? 6 : J.Hi();*/ memset(Jm,0,icv*icv*sizeof(complexdouble)); std::fstream FILEIN; 
   int k[] = {0,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   std::vector<double> mgjmbH(Jhi,0.);
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc();
      for (int iq=Jlo; iq<=Jhi; iq++) mgjmbH[iq] = -gjmbH[iq]*pars.jijconv[iq]; }
   else
   #endif
      for (int iq=Jlo; iq<=Jhi; iq++) mgjmbH[iq] = -gjmbH[iq];
   complexdouble **LS = new complexdouble*[6]; 
   int iq,j; for(iq=0; iq<6; iq++) { LS[iq] = new complexdouble[icv*icv]; memset(LS[iq],0,icv*icv*sizeof(complexdouble)); }
   for(iq=Jlo; iq<=Jhi; iq++)
   {
      if(iq<=6)
      {
         fn[16]=iq+48; FILEIN.open(fn,std::fstream::in); FILEIN.precision(24);  // In ascii, 48=="0"
         if(iq==3||iq==4) for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> LS[iq-1][i+j*icv].i; Jm[i+j*icv].i += LS[iq-1][i+j*icv].i*mgjmbH[iq]; }
         else             for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> LS[iq-1][i+j*icv].r; Jm[i+j*icv].r += LS[iq-1][i+j*icv].r*mgjmbH[iq]; }
         FILEIN.close();
      }
      else
      {
         if(fabs(mgjmbH[iq])>(DBL_EPSILON*1000))
         {
            char fh[] = "results/ic1ion.uXmX"; fh[16]=k[iq]+48; 
            if(q[iq]<0) 
            {
               fh[18]=abs(q[iq])+48; FILEIN.open(fh, std::fstream::out); FILEIN.precision(24);
               for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> elem; Jm[i+j*icv].i += elem*mgjmbH[iq]; } FILEIN.close();
            }
            else 
            { 
               fh[17]=q[iq]+48; fh[18]=0; FILEIN.open(fh, std::fstream::out); FILEIN.precision(24);
               for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> elem; Jm[i+j*icv].r += elem*mgjmbH[iq]; } FILEIN.close();
            }
         }
      }
   }
   complex<double> a(1.,0.); int incx = 1;
   for(i=1; i<=icv; i++) F77NAME(zaxpy)(&icv,(complexdouble*)&a,(complexdouble*)&est[i][1],&incx,&Jm[(i-1)*icv],&incx);
   // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
   iceig VE; VE.calc(icv,Jm);
   // Calculates the moments, energy and partition function.
   std::vector<double> E(icv,0.),eb; int Esz; double U = 0.; double Z = 0.; char uplo = 'U';
   for(Esz=0; Esz<icv; Esz++) { E[Esz] = VE.E(Esz)-VE.E(0); eb.push_back(exp(-E[Esz]/(KB*T))); if(eb[Esz]<DBL_EPSILON) break; }
   complexdouble *zt = new complexdouble[icv]; complexdouble zalpha,zbeta,zme; zalpha.r=1.; zalpha.i=0.; zbeta.r=0.; zbeta.i=0.;
   for(iq=Jlo; iq<=Jhi; iq++) 
   { 
      if((iq>6 && k[iq]%2==1) || (k[iq]>4 && pars.l==D)) continue;
      J[iq] = 0.;
      if(iq>6)
      { 
         memset(Jm,0,icv*icv*sizeof(complexdouble)); char fh[] = "results/ic1ion.uXmX"; fh[16]=k[iq]+48; 
         if(q[iq]<0) 
         {
            fh[18]=abs(q[iq])+48; FILEIN.open(fh, std::fstream::out); FILEIN.precision(24); 
            for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> elem; Jm[i+j*icv].i = elem; } FILEIN.close();
         }
         else
         {
            fh[17]=q[iq]+48; fh[18]=0; FILEIN.open(fh, std::fstream::out); FILEIN.precision(24);
            for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> elem; Jm[i+j*icv].r = elem; } FILEIN.close();
         }
      }
      for(i=0; i<Esz; i++)
      {
         if(iq<=6)
            F77NAME(zhemv)(&uplo, &icv, &zalpha, LS[iq-1], &icv, VE.zV(i), &incx, &zbeta, zt, &incx); 
         else
            F77NAME(zhemv)(&uplo, &icv, &zalpha, Jm, &icv, VE.zV(i), &incx, &zbeta, zt, &incx); 
#ifdef _G77
         F77NAME(zdotc)(&zme, &icv, VE.zV(i), &incx, zt, &incx);
#else
         zme = F77NAME(zdotc)(&icv, VE.zV(i), &incx, zt, &incx);
#endif
         J[iq]+=zme.r*eb[i]; if(iq==Jlo) { Z+=eb[i]; U+=(E[i]+VE.E(0))*eb[i]; }
      }
      J[iq]/=Z; if(fabs(J[iq])<DBL_EPSILON) J[iq]=0.;
      #ifdef JIJCONV
      if(pars.B.norm().find("Stevens")!=std::string::npos) J[iq]*=pars.jijconv[iq];
      #endif
   }
   for(i=0; i<6; i++) delete[]LS[i]; delete[]LS;
   delete[]zt; delete[]Jm; 
   J[iq] = log(Z)-VE.E(0)/(KB*T); J[iq+1] = U/Z;
   return J;
}

/* truncate.cpp
 *
 * Diagonalises the single ion matrix without mean field or Zeeman terms, generates a new basis from this
 * diagonal matrix and truncates this rotated Hamiltonian to include only lowest lying terms in order to
 * save computation time.
 *
 * Functions:
 *   void myPrintMatrix(FILE * file,sMat<double> & M,int d)                        // Prints out full matrix
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 */

 /* F  = [0 68878 50346 32901]./8.066; xi = 751.7./8.066; alpha = [16.23 -566.6 1371]./8.066;
    B  = {zeros(1,5) [sqrt(5/14) 0 0 0 1 0 0 0 sqrt(5/14)].*2365.6./8.066 [0 0 -sqrt(7/2) 0 0 0 1 0 0 0 -sqrt(7/2) 0 0].*10};
    for i=1:3; for j=1:(2*i); B{i}(j) = 0; end; end
    conf='f5'; n=5; NKPT=-1;
    st_SI = racah_states(n,3); id = 0; icv = 0; icv1 = 0; icv2 = 0; for i = 1:length(st_SI); L = racah_lconv(st_SI{i}{2}); 
     S = st_SI{i}{1}; Jmin = abs(L-S); Jmax = L+S; cvSI2SO(i,1) = icv1+1; icv1=icv1+(Jmax-Jmin+1); cvSI2SO(i,2) = icv1; 
     cvSI2CF(i,1) = icv2+1; icv2=icv2+sum(2.*[Jmin:Jmax]+1); cvSI2CF(i,2) = icv2; for J = Jmin:Jmax; id = id + 1; cvSO2CF(id,1) = icv+1; 
     icv=icv+(2*J+1); cvSO2CF(id,2) = icv; end; end
    Hcf = full(ic_hmltn(conf,racah_FtoF_k(F),xi,B,alpha)); Hcf = triu(Hcf)+triu(Hcf,1)'; [Vf,E] = eig(Hcf); E=diag(E)'-min(min(E)); 
    H = full(ic_freehmltn(conf,racah_FtoF_k(F),xi,alpha)); Hfull = full(convH2H(H,icv,cvSO2CF)); Hfull=triu(Hfull)+triu(Hfull,1)'; [Vfree,Efree]=eig(Hfull);
    nn=91; [Vcf,Ecf] = eig(Hcf); Ecf=diag(Ecf); Vcf(find(abs(Vcf)<1e-3))=0; Hrot=Vfree(:,1:nn)'*Hcf*Vfree(:,1:nn); 
    Hrot=triu(Hrot)+triu(Hrot,1)'; [Vcl,Ecl]=eig(Hrot); Vcl(find(abs(Vcl)<1e-3))=0; */

// --------------------------------------------------------------------------------------------------------------- //
// Truncates the matrix, and stores the single ion part in a memorymapped array
// --------------------------------------------------------------------------------------------------------------- //
void truncate_hmltn(icpars &pars, ComplexMatrix &est, sMat<double> &Hic, sMat<double> &iHic, int JHi, int JLo)
{
   std::cout << "mcalc(): Calculating rotated matrix for truncation." << std::flush;
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

   // Set up directory to store matrices if the user asks for it.
   char nstr[6]; char filename[255]; char basename[255];
   if(pars.save_matrices) 
   {
      #ifndef _WINDOWS
      struct stat status; int dirstat=0; stat("results/mms/",&status); if(!S_ISDIR(status.st_mode)) dirstat = mkdir("results/mms",0777);
      if(dirstat!=0) { std::cerr << "mcalc(): " << errno << "\n"; exit(EXIT_FAILURE); }
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
   Hrot = new complexdouble[cb*cb]; int memloc=0;//Hsz*Hsz;
   // Calculates the rotated single ion Hamiltonian
   char notranspose='N',transpose='C',uplo='U',side='L'; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   if(iHic.isempty()) zJmat=zmat2f(Hic,zeroes); else zJmat = zmat2f(Hic,iHic);
   F77NAME(zhemm)(&side,&uplo,&Hsz,&cb,&zalpha,zJmat,&Hsz,Vf,&Hsz,&zbeta,zmt,&Hsz);
   F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hsz,&zalpha,Vf,&Hsz,zmt,&Hsz,&zbeta,Hrot,&cb); free(zJmat);
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(Hrot[ii*cb+jj].r)<DBL_EPSILON) Hrot[ii*cb+jj].r=0.; if(fabs(Hrot[ii*cb+jj].i)<DBL_EPSILON) Hrot[ii*cb+jj].i=0.; } 
   memcpy(&est[memloc+100][0],Hrot,cb*cb*sizeof(complexdouble)); memloc+=cb*cb;
   // Calculates the rotated multipolar operators for the mean field terms
   std::cout << " Using " << cb << " levels of " << Hsz << ".\nmcalc(): Starting calculation of rotated mean field operators... " << std::flush;
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
      memcpy(&est[memloc+100][0],Hrot,cb*cb*sizeof(complexdouble)); memloc+=cb*cb;
   }
   delete[]Vf; delete[]Hrot; delete[]zmt;
   end = clock(); std::cout << "Done. Time to set up rotated matrices = " << (double)(end-start)/CLOCKS_PER_SEC << "s." << std::endl;
}

// --------------------------------------------------------------------------------------------------------------- //
// Uses the memory mapped eigenvectors of the single ion Hamiltonian to truncate the matrix
// --------------------------------------------------------------------------------------------------------------- //
void truncate_expJ(icpars &pars, ComplexMatrix &est, Vector &gjmbH, Vector &J, double T, double *lnZ, double *U, complexdouble *Jm)
{
   int Hsz=getdim(pars.n,pars.l);
 //size_t filesize = Hsz*Hsz*sizeof(complexdouble); 
   char uplo='U'; complexdouble zme;
   int Esz, incx=1; std::vector<double> E, me, eb;
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;

   int cb = (int)(pars.truncate_level*Hsz); complexdouble *Hrot; Hrot = new complexdouble[cb*cb]; 
   memcpy(Hrot,&est[100][0],cb*cb*sizeof(complexdouble)); int szapy=cb*cb; complexdouble a; a.r=1.; a.i=0.;
   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

   // Calculates the mean field Hamiltonian = H_singleion + sum_i(gjmbH[i]*Operator[i])
   for(int iJ=1; iJ<=(gjmbH.Hi()-gjmbH.Lo()+1); iJ++)
   {
      if (q[iJ]<0) a.r = -gjmbH[iJ+gjmbH.Lo()-1]; else a.r = -gjmbH[iJ+gjmbH.Lo()-1];
      if (fabs(a.r)>DBL_EPSILON) F77NAME(zaxpy)(&szapy,&a,(complexdouble*)&est[iJ*cb*cb+100][0],&incx,Hrot,&incx);
   }

   // Diagonalises the rotated mean field Hamiltonian
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(Hrot[ii+jj*cb].r)<DBL_EPSILON) Hrot[ii+jj*cb].r=0.; if(fabs(Hrot[ii+jj*cb].i)<DBL_EPSILON) Hrot[ii+jj*cb].i=0.; } 
   iceig VE; VE.calc(cb,Hrot); delete[]Hrot;
   for(int ii=0; ii<cb; ii++) for(int jj=0; jj<cb; jj++) { 
      if(fabs(VE.zV(ii,jj).r)<DBL_EPSILON && fabs(VE.zV(ii,jj).i)<DBL_EPSILON) VE.zV(ii,jj) = 0.; } 

   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<cb; Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || VE.E(Esz+1)==0) break; }

   // Calculates the rotated operators for the mean field terms
   char nstr[6]; char mapname[255]; char mapbasename[255];
   nstr[0] = (pars.l==F?102:100); if(pars.n<10) { nstr[1] = pars.n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = pars.n+38; nstr[3] = 0; }
 //filesize = cb*cb*sizeof(complexdouble); 
   complexdouble *zt; double Z=0.; eb.assign(Esz,0.); *U=0.;
   int memloc=cb*cb;
   for(int iJ=(J.Lo()-1); iJ<J.Hi(); iJ++)
   {
      me.assign(Esz,0.); strcpy(mapname,mapbasename); 
      if(iJ<6) { nstr[1]=49; nstr[2]=iJ+49; nstr[3]=0; strcat(mapname,nstr); }
      else { if(q[iJ]<0) { MSTR(k[iJ],abs(q[iJ])); strcat(mapname,nstr); } else { NSTR(k[iJ],abs(q[iJ])); strcat(mapname,nstr); } }
      strcat(mapname,".mmap");
      zt = (complexdouble*)malloc(cb*sizeof(complexdouble)); J[iJ+1]=0.; 
      for(int ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         F77NAME(zhemv)(&uplo, &cb, &zalpha, (complexdouble*)&est[memloc+100][0], &cb, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
         #ifdef _G77 
         F77NAME(zdotc)(&zme, &cb, VE.zV(ind_j), &incx, zt, &incx);
         #else
         zme = F77NAME(zdotc)(&cb, VE.zV(ind_j), &incx, zt, &incx);
         #endif
         me[ind_j] = zme.r;
         if(iJ==(J.Lo()-1)) { eb[ind_j] = exp(-E[ind_j]/(KB*T)); Z+=eb[ind_j]; *U+=(E[ind_j]+VE.E(0))*eb[ind_j]; }
         J[iJ+1]+=me[ind_j]*eb[ind_j];
      }
      free(zt); J[iJ+1]/=Z; if(iJ==(J.Lo()-1)) *U/=Z; 
      memloc+=cb*cb;
   }
   *lnZ = log(Z)-VE.E(0)/(KB*T);
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the eigenvalues/vectors of the IC Hamiltonian by a method similar to that used in the XTALS84 programs
//   of Hannah Crosswhite. First the basis states are re-ordered with states having the same J-values following
//   each other consecutively. Then the Hc+Hso interactions are calculated in this basis, yielding a matrix of
//   blocks along the diagonal (Hc and Hso are diagonal in J,mJ). Each J-submatrix is diagonalised individually
//   to save time, and the full eigenvector matrix Vrot (also block-diagonal) is reconstructed from the individual 
//   J-submatrices. This matrix is truncated at some energy [elim] or level [cb] (e.g. only cb columns of the 
//   matrix are used), and used to rotate the full IC Hamiltonian by: V(:,1:cb)'*Hic*V(:,1:cb), yielding a cb*cb
//   matrix which is diagonalised. The eigenvalues are the energies, and the eigenvector in the original 
//   basis is given by (Vr*Vrot')' if Vr are the eigenvectors in the rotated basis.
// --------------------------------------------------------------------------------------------------------------- //
iceig spectre_eig(icpars &pars, double elim)
{
   if(elim>=0) elim *= MEV2CM;

   // Calculates the Coulomb and Spin orbit Hamiltonian matrix in the |alpha,LSJ> basis
   sMat<double> Hcso = ic_Hcso(pars);

   // Finds indices to allow reordering the basis states so states with same J are grouped together
   unsigned int i, j, num_states=Hcso.nr(), ns=getdim(pars.n,pars.l), imax=num_states, cb=ns; 
   int ii=0, jj, imj=0, J2;
   fconf conf(pars.n,0,pars.l);
   std::vector< std::vector<int> > JsortM, JsortMmj; std::vector<int> JsortV, Jind, JsortVmj(ns,0);
   for(i=0; i<26; i++) JsortM.push_back(JsortV); JsortMmj = JsortM;
   for(i=0; i<conf.states.size(); i++) 
   { 
      J2=conf.states[i].J2; JsortM[J2].push_back(i); 
      for(j=0; j<(unsigned int)J2+1; j++) JsortMmj[J2].push_back(imj++); 
   }
   for(i=JsortM.size()-1; i<JsortM.size(); i--) 
      if(!JsortM[i].empty()) 
      {
         for(j=0; j<JsortM[i].size(); j++) { JsortV.push_back(JsortM[i][j]); Jind.push_back(i); }
         for(j=0; j<JsortMmj[i].size(); j++) JsortVmj[JsortMmj[i][j]]=ii++; 
      }
   i=0; while(i<JsortM.size()) { if(JsortM[i].empty()) JsortM.erase(JsortM.begin()+i); else i++; }

   // Split Hc+Hso into J-submatrices
   double **submats = new double*[JsortM.size()];
   for(i=0; i<JsortM.size(); i++) {
      submats[i] = new double[JsortM[i].size()*JsortM[i].size()];
      memset(submats[i],0,JsortM[i].size()*JsortM[i].size()*sizeof(double)); }
   std::vector< std::vector<int> > nz = Hcso.find();
   std::vector<int> i1(num_states,0),i2(num_states,0),i3(num_states,0);
   for(i=0; i<JsortM.size(); i++) for(j=0; j<JsortM[i].size(); j++) { i1[JsortM[i][j]] = i; i2[JsortM[i][j]]=j; }
   for(unsigned int kk=0; kk<nz.size(); kk++) 
   {
      ii = nz[kk][0]; jj = nz[kk][1];
      if(i1[ii]==i1[jj]) submats[i1[ii]][i2[jj]*JsortM[i1[ii]].size()+i2[ii]] = Hcso(ii,jj);
   }

   // Diagonalises the submatrices, and reconstruct the full eigenvector matrix from the eigenvector 
   //   matrices of the submatrices, and expand to the |alpha,LSJmJ> basis
   double *Vunrot = new double[num_states*num_states], *Eunrot = new double[num_states]; 
   memset(Vunrot,0,num_states*num_states*sizeof(double));
   int idx = 0;
   iceig VEs;
   for(i=JsortM.size()-1; i<JsortM.size(); i--)
   {
      int isz = JsortM[i].size();
      VEs.calc(isz,submats[i]);
      for(ii=0; ii<isz; ii++) { 
         for(int jj=0; jj<isz; jj++) Vunrot[(jj+idx)*num_states+ii+idx] = VEs.V(ii,jj); Eunrot[ii+idx] = VEs.E(ii); } 
      idx+=isz; 
   }
   for(i=0; i<JsortM.size(); i++) delete[]submats[i];
   delete[]submats;

   // Finds the indices of the sorted Energies of Hc+Hso, to determine the state cb at which to
   //    truncate the eigenvector matrix used to rotate Hic from the maximum energy limit elim
   std::vector<int> Esort(num_states,0); for(i=0; i<num_states; i++) Esort[i] = i;
   i=1;j=2;
   while(i<num_states)
   {
      if(Eunrot[Esort[i-1]]<=Eunrot[Esort[i]]) { i=j; j++; }
      else { ii=Esort[i-1]; Esort[i-1]=Esort[i]; Esort[i]=ii; i--; if(i==0) i=1; }
   }
   if(elim>=0) for(imax=0; imax<num_states; imax++) if((Eunrot[Esort[imax]]-Eunrot[Esort[0]])>elim) break; imax--;

   // Construct the rotation matrix from the full eigenvector matrix, but putting the columns in Energy-order
   complexdouble *Vrot = new complexdouble[ns*ns];
   unsigned int iy=0, jy=0, id, J2i, J2j;
   for(i=0; i<num_states; i++) 
   { 
      jy=0; J2i = Jind[i];
      for(j=0; j<num_states; j++) 
      { 
         J2j = Jind[Esort[j]];
         if(fabs(Vunrot[Esort[j]*num_states+i])>SMALL) for(id=0; id<J2i+1; id++) Vrot[(iy+id)+(jy+id)*ns].r = Vunrot[Esort[j]*num_states+i];
         jy+=J2j+1;
         if(i==0) if(j==imax) cb=jy;
      } 
      iy+=J2i+1;
   }
   delete[]Vunrot;
   if(elim<0) cb = (unsigned int)elim;

   // Calculates the free-ion + crystal field Hamiltonian in the |alpha,LSJmJ> basis
   sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars);
   sMat<double> zeroes; zeroes.zero(ns,ns);
   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(pars.Bx)>DBL_EPSILON || fabs(pars.By)>DBL_EPSILON || fabs(pars.Bz)>DBL_EPSILON)
   {
      std::vector<double> gjmbH(6,0.);
      if(fabs(pars.Bx)>DBL_EPSILON) { gjmbH[1]=-MUBc*pars.Bx; gjmbH[0]=GS*gjmbH[1]; }
      if(fabs(pars.By)>DBL_EPSILON) { gjmbH[3]=-MUBc*pars.By; gjmbH[2]=GS*gjmbH[3]; }
      if(fabs(pars.Bz)>DBL_EPSILON) { gjmbH[5]=-MUBc*pars.Bz; gjmbH[4]=GS*gjmbH[5]; }
      sMat<double> J,iJ; icmfmat mfmat(pars.n,pars.l,6,pars.save_matrices); mfmat.Jmat(J,iJ,gjmbH,pars.save_matrices); Hic+=J; iHic+=iJ;
   }
   Hic.reorder(JsortVmj); iHic.reorder(JsortVmj);

   // Rotates and truncates Hic by using the eigenvectors of Hc+Hso: Hic_r = V(:,1:cb)'*Hic*V(:,1:cb)
   complexdouble *zJmat,*Hrot,*zmt; zmt = new complexdouble[ns*cb]; Hrot = new complexdouble[cb*cb];
   char notranspose='N',transpose='C',uplo='U',side='L'; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   int ins=(int)ns, icb=(int)cb;
   zJmat = zmat2f(Hic,iHic);
   F77NAME(zhemm)(&side,&uplo,&ins,&icb,&zalpha,zJmat,&ins,Vrot,&ins,&zbeta,zmt,&ins);
   F77NAME(zgemm)(&transpose,&notranspose,&icb,&icb,&ins,&zalpha,Vrot,&ins,zmt,&ins,&zbeta,Hrot,&icb); free(zJmat);

   // Diagonalises the truncated and rotated matrix
   iceig VEr; VEr.calc(icb,Hrot); delete[]Hrot;

   // Transform the eigenvectors of the rotated basis back into the original |alpha,LSJmJ> basis
   memset(zmt,0.,ns*cb*sizeof(complexdouble));
   complexdouble *Vr = new complexdouble[cb*ns*sizeof(complexdouble)]; 
   for(i=0; i<ns; i++) for(j=0; j<cb; j++) Vr[i*cb+j] = VEr.zV(i,j);
   F77NAME(zgemm)(&notranspose,&transpose,&icb,&ins,&icb,&zalpha,Vr,&icb,Vrot,&ins,&zbeta,zmt,&icb);  // (Vr*Vrot')'
   delete[]Vr;
   bool imagfl=false; for(i=0; i<ns; i++) for(j=0; j<cb; j++) if(zmt[i*cb+j].i>SMALL) { imagfl=true; break; }
   iceig VE(ns,!imagfl); double *En = VE.E(); for(i=0; i<cb; i++) En[i] = VEr.E(i); for(i=cb; i<ns; i++) En[i] = -DBL_MAX;
   if(imagfl) 
   {
      complexdouble *zV = VE.zV(0);
      for(i=0; i<ns; i++) for(j=0; j<cb; j++) { zV[j*ns+i].r = zmt[JsortVmj[i]*cb+j].r; zV[j*ns+i].i = -zmt[JsortVmj[i]*cb+j].i; }
   }
   else
   {
      double *V = VE.V(0); for(i=0; i<ns; i++) for(j=0; j<cb; j++) V[j*ns+i] = zmt[JsortVmj[i]*cb+j].r; 
   }

   return VE;
}

