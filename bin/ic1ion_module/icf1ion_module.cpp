/* icf1ion_module.cpp
 *
 * Functions:
 *   void myPrintMatrix(FILE * file,sMat<double> & M,int d)                        // Prints out full matrix
 *   bool checkmat(ComplexMatrix &cmat, complexdouble *fmat,int r, int c)          // Compares Matpack and fortran matrices
 *   void mcalc(Vector &J, double *T, Vector &gjmbH, double *gJ, Vector &ABC,      // Calculates the meanfield moment
 *                 char **sipffile, double *lnZ, double *U, ComplexMatrix &est)    //
 *   int dmcalc(int &tn, double &T, Vector &gjmbH, double &g_J, Vector &ABC,       // Calculates the transition
 *                 char **sipffilename, ComplexMatrix &mat, float &delta,          //   matrix elements
 *                 ComplexMatrix &est)                                             //
 *   void mcalc_parameter_storage_matrix_init(ComplexMatrix *est,Vector &gjmbheff  // initialises parameter storage matrix 
 *                 double *g_J,double &T,Vector &ABC,char **sipffilename)          // for mcalc
 *   void estates(ComplexMatrix *est, Vector &gjmbheff, double *g_J, double &T,    // Calculates the energy and wavefunctions
 *                 Vector &ABC, char **sipffilename)                               //   "estates" matrix
 *   bool get_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l,         // Calculates/loads the Q_q operators
 *                 std::vector<double> &Jvec)                                      //   for beyond dipole calculations
 *   void save_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l,        // Saves the Q_q operators to a file
 *                 std::vector<double> &Jvec)                                      //
 *   void mq(ComplexVector &Mq, double &th, double &ph, double &J0, double &J2,    // Calculates the thermal expectation
 *                 double &J4, double &J6, ComplexMatrix &est)                     //   of the magnetisation density
 *   int dncalc(int &tn, double &th, double &ph, double &J0, double &J2,           // Calculates the transition matrix
 *                 double &J4, double &J6, ComplexMatrix &est, double &T,          //   elements beyond the dipole
 *                 ComplexMatrix &mat)                                             //   approximation.
 *   void spindensity_mcalc(Vector & mom, int & xyz,double *T, Vector &gjmbH,      // Calc. coeffs. of expansion of spindensity 
 *                 double *gJ,Vector &ABC, char **sipffile, ComplexMatrix &est)    //   in terms of Zlm R^2(r) at given T / H_eff
 *   void orbmomdensity_mcalc(Vector & mom,int & xyz, double *T, Vector &gjmbH,    // Calc. coeffs. of expansion of orbital moment density
 *                 double *gJ,Vector &ABC, char **sipffile, ComplexMatrix &est)    //   in terms of Zlm F(r) at given T / H_eff
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008-2010 Duc Le, Martin Rotter
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"
#include "vector.h"          // MatPack vector class
#include <fstream>
#include <ctime>
#define SMALL 1e-6           // must match SMALL in mcdisp.c and ionpars.cpp because it is used to decide wether for small
                             // transition, energy the matrix Mijkl contains wn-wn' or wn/kT

#define MAXNOFCHARINLINE 7024

// --------------------------------------------------------------------------------------------------------------- //
// Looks up the Hund's Rule ground state to determine L and S for a particular configuration p^n, d^n, f^n
// --------------------------------------------------------------------------------------------------------------- //
fstates_t hunds_gs(int n, orbital l=F)  // Defaults to f-electrons
{
   qG2 U; 
   if(n>(2*l+1)) n = (4*l+2)-n;
   if(n<1) { std::cerr << "Invalid number of electrons\n"; exit(1); }

   switch(l)
   {
      case P:
         switch(n)
         {
            case 1: return fstates_t(1,P,"2P");                           // p1    2P
            case 2: return fstates_t(2,P,"3P");                           // p2    3P
            case 3: return fstates_t(3,S,"4S");                           // p3    4S
            default: std::cerr << "Invalid number of electrons\n"; exit(1); 
         }
      case D:
         switch(n)
         {
            case 1: return fstates_t(1,D,1,"2D");                         // d1    2D   1 (10)
            case 2: return fstates_t(2,F,2,"3F");                         // d2    3P   2 (11)
            case 3: return fstates_t(3,F,3,"4F");                         // d3    4P   3 (11)
            case 4: return fstates_t(4,D,4,"5D");                         // d4    5D   4 (10)
            case 5: return fstates_t(5,S,5,"6S");                         // d5    6S   5 (00)
            default: std::cerr << "Invalid number of electrons\n"; exit(1); 
         }
      case F:
         switch(n)
         {
            case 1: U.set(1,0); return fstates_t(1,F,1,U,"2F");           // f1    2F   1 100 10
            case 2: U.set(1,1); return fstates_t(2,H,2,U,"3H");           // f2    3P   2 110 11
            case 3: U.set(2,0); return fstates_t(3,I,3,U,"4I");           // f3    4S   3 111 00
            case 4: U.set(2,0); return fstates_t(4,I,4,U,"5I");           // f4    5S   4 111 00
            case 5: U.set(1,1); return fstates_t(5,H,5,U,"6H");           // f5    6P   5 110 11
            case 6: U.set(1,0); return fstates_t(6,F,6,U,"7F");           // f6    7F   6 100 10
            case 7: U.set(0,0); return fstates_t(7,S,7,U,"8S");           // f7    8S   7 000 00
            default: std::cerr << "Invalid number of electrons\n"; exit(1); 
         }
      default: std::cerr << "Only p-, d- and f-electrons supported\n"; exit(1); 
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Determines number of basis states
// --------------------------------------------------------------------------------------------------------------- //
int icf_getdim(icpars &pars)
{
   fstates_t gs = hunds_gs(pars.n, pars.l);
   int J2min, J2max, ns=0;
   J2min = abs(2*abs(gs.L)-gs.S2); J2max = 2*abs(gs.L)+gs.S2;
   for (int J2=J2min; J2<=J2max; J2+=2) ns+=J2+1; 
   return ns;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the Hamilton matrix for the lowest Coulombic term (manifold of constant L and S)
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> icf_hmltn(sMat<double> &Hcfi, icpars &pars)
{
   int n = pars.n; double xi = pars._xi; orbital e_l = pars.l;
   if (e_l>3||e_l<1) { std::cerr << "Sorry only p-, d- and f-electrons supported at present\n"; exit(0); }
   fstates_t gs = hunds_gs(n, e_l);
   int J2min, J2max, ns=0;
   int S2 = gs.S2, L2 = abs(gs.L)*2;

   // Determines number of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2;
   for (int J2=J2min; J2<=J2max; J2+=2) ns+=J2+1; 
   sMat<double> Hcf(ns,ns); Hcfi.zero(ns,ns);
   // Determines angular momentum quantum numbers of basis states
   int ins=0; std::vector<int> J2(ns), mJ2(ns);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2)
     for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2)
     {
        J2[ins]=iJ2; mJ2[ins]=imJ2; ins++;
     }

   // Determines reduced matrix elements (Stevens operator equivalent factor)
   double rmso = 0., rmU[7]={0.,0.,0.,0.,0.,0.,0.};
   double p = 1./( pow(-1.,(double)abs(e_l))*(2.*e_l+1.) );
   if (n>1)
   {
      fconf confp(n-1,e_l);
      std::vector<cfpls> cfps;
      switch(e_l) {
         case P:  cfps = racah_parents(n,gs.S2,gs.L); break;
         case D:  cfps = racah_parents(n,gs.v,gs.S2,gs.L); break;
         default: cfps = racah_parents(n,gs.v,gs.U,gs.S2,gs.L);  }
      int sz = (int)cfps.size();
      for(int k=0; k<sz; k++)
      {
         int pS2 = confp.states[cfps[k].ind].S2, pL2 = abs(confp.states[cfps[k].ind].L)*2;
         rmso   += racahW(pS2,S2,1,2,1,S2) * racahW(pL2,L2,2*e_l,2,2*e_l,L2) * cfps[k].cfp * cfps[k].cfp; if(n!=(4*e_l+1)) {
         rmU[2] += pow(-1.,(double)abs(pL2+L2)/2.+e_l) * sixj(L2,4,L2,2*e_l,pL2,2*e_l) * cfps[k].cfp * cfps[k].cfp * (L2+1.);
         rmU[4] += pow(-1.,(double)abs(pL2+L2)/2.+e_l) * sixj(L2,8,L2,2*e_l,pL2,2*e_l) * cfps[k].cfp * cfps[k].cfp * (L2+1.);
         rmU[6] += pow(-1.,(double)abs(pL2+L2)/2.+e_l) * sixj(L2,12,L2,2*e_l,pL2,2*e_l)* cfps[k].cfp * cfps[k].cfp * (L2+1.); }
      }
      if(n==(4*e_l+1)) { rmU[2]=1./n; rmU[4]=1./n; rmU[6]=1./n; }
      for(int ik=2; ik<=6; ik+=2) rmU[ik] *= n * threej(2*e_l,2*ik,2*e_l,0,0,0) / p;
   }
   else  // Single electron
   {
   // rmso = sqrt(3/2.)*sqrt(e_l*(e_l+1)*(2*e_l+1));             // s=1/2 substituted into eqn 4-12 of Judd 1963
      rmso = 1./((L2+1.)*(S2+1.));
      rmU[2] = threej(2*e_l,4,2*e_l,0,0,0) / p;                  // See Judd 1963, Eqn 5-13. with U^k=V^k/sqrt(2k+1)
      rmU[4] = threej(2*e_l,8,2*e_l,0,0,0) / p;
      rmU[6] = threej(2*e_l,12,2*e_l,0,0,0)/ p;
   }

   double elp, elm; int k, q, iq;
   for (int i=0; i<ns; i++) 
      for (int j=0; j<ns; j++)
      {
         if(i==j)  // Selection rule for Spin-orbit operator J=J', mJ=mJ'
         {
            Hcf(i,i) += -n*xi * racahW(J2[i],L2,S2,2,S2,L2) * (L2+1)*(S2+1) * sqrt( (9./6)*e_l*(e_l+1)*(2*e_l+1) ) * rmso;
         // Hcf(i,i) = -xi * pow(-1.,(S2+L2+J2[i])/2.) * sixj(S2,S2,2,L2,L2,J2[i]) * rmso;
         }
         for(k=2; k<=6; k+=2) for(iq=0; iq<(2*k+1); iq++)
         {
            q = iq-k; if(fabs(pars.B(k,q))<1e-10) continue;
            if(q==0)
            {
               elp = pow(-1.,(S2-L2-J2[j])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * racahW(L2,J2[i],L2,J2[j],S2,2*k) * rmU[k] 
                        * pow(-1.,(J2[i]+mJ2[i])/2.+k) * wigner(J2[i],J2[j],0-mJ2[i],mJ2[j],2*k,0) / sqrt(2.*k+1.); 
               Hcf(i,j) += elp * pars.B(k,q);
            }
            else
            {
               elp = pow(-1.,(S2-L2-J2[j])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * racahW(L2,J2[i],L2,J2[j],S2,2*k) * rmU[k] * pow(-1.,(J2[i]+mJ2[i])/2.+k);
               elm = elp * wigner(J2[i],J2[j],0-mJ2[i],mJ2[j],2*k,-2*abs(q)) / sqrt(2.*k+1.); 
               elp *=      wigner(J2[i],J2[j],0-mJ2[i],mJ2[j],2*k, 2*abs(q)) / sqrt(2.*k+1.); 
               if(q<0)
                 Hcfi(i,j)+= (elp-elm*pow(-1.,q)) * pars.B(k,q); 
               else
                 Hcf(i,j) += (elp+elm*pow(-1.,q)) * pars.B(k,q); 
            }
         }
      }
   return Hcf;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the Multipolar operator matrices
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> icf_ukq(int n, int k, int q, orbital e_l)
{
   if (e_l>3||e_l<1) { std::cerr << "Sorry only p-, d- and f-electrons supported at present\n"; exit(0); }
   fstates_t gs = hunds_gs(n, e_l);
   int J2min, J2max, ns=0;
   int S2 = gs.S2, L2 = abs(gs.L)*2;

   // Determines number of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2;
   for (int J2=J2min; J2<=J2max; J2+=2) ns+=J2+1; 
   sMat<double> Hcf(ns,ns);
   // Determines angular momentum quantum numbers of basis states
   int ins=0; std::vector<int> J2(ns), mJ2(ns);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2)
     for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2)
     {
        J2[ins]=iJ2; mJ2[ins]=imJ2; ins++;
     }

   // Determines reduced matrix elements (Stevens operator equivalent factor)
   double rmU=0., p = 1./( pow(-1.,(double)abs(e_l))*(2.*e_l+1.) );
   if (n>1 && n<(4*e_l+1))
   {
      fconf confp(n-1,e_l);
      std::vector<cfpls> cfps;
      switch(e_l) {
         case P:  cfps = racah_parents(n,gs.S2,gs.L); break;
         case D:  cfps = racah_parents(n,gs.v,gs.S2,gs.L); break;
         default: cfps = racah_parents(n,gs.v,gs.U,gs.S2,gs.L);  }
      int sz = (int)cfps.size();
      for(int kk=0; kk<sz; kk++)
      {
         int pL2 = abs(confp.states[cfps[kk].ind].L)*2;
         rmU += pow(-1.,(double)abs(pL2+L2)/2.+e_l) * sixj(L2,2*k,L2,2*e_l,pL2,2*e_l) * cfps[kk].cfp * cfps[kk].cfp * (L2+1.);
      }
      rmU *= n * threej(2*e_l,2*k,2*e_l,0,0,0) / p;
   }
   else  // Single electron or single hole
      rmU = threej(2*e_l,2*k,2*e_l,0,0,0) / p;                   // See Judd 1963, Eqn 5-13. with U^k=V^k/sqrt(2k+1)

   double elp, elm;
   if(q==0)
      for (int i=0; i<ns; i++)
         for (int j=0; j<ns; j++)
         {
            Hcf(i,j) = pow(-1.,(S2-L2-J2[j])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * racahW(L2,J2[i],L2,J2[j],S2,2*k) * rmU 
                       * pow(-1.,(J2[i]+mJ2[i])/2.+k) * wigner(J2[i],J2[j],0-mJ2[i],mJ2[j],2*k,0) / sqrt(2.*k+1.); 
         }
   else
      for (int i=0; i<ns; i++)
         for (int j=0; j<ns; j++)
         {
            elp = pow(-1.,(S2-L2-J2[j])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * racahW(L2,J2[i],L2,J2[j],S2,2*k) * rmU * pow(-1.,(J2[i]+mJ2[i])/2.+k);
            elm = elp * wigner(J2[i],J2[j],0-mJ2[i],mJ2[j],2*k,-2*abs(q)) / sqrt(2.*k+1.); 
            elp *=      wigner(J2[i],J2[j],0-mJ2[i],mJ2[j],2*k, 2*abs(q)) / sqrt(2.*k+1.); 
            if(q<0)
              Hcf(i,j) += (elp-elm*pow(-1.,q));
            else
              Hcf(i,j) += (elp+elm*pow(-1.,q));
         }
   return Hcf;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the magnetic moment matrices
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> icf_mumat(int n, int ind, orbital e_l=F)
{
   fstates_t gs = hunds_gs(n, e_l);
   int J2min, J2max, ns=0, L2=2*abs(gs.L), S2=gs.S2;
   J2min = abs(L2-S2); J2max = L2+S2;
   sMat<double> rm(J2min-J2max+1,J2min-J2max+1);
   for (int J2=J2min; J2<=J2max; J2+=2) 
   {
      ns+=J2+1; 
      if (ind%2==0)                                     // Sx, Sy or Sz
         for (int J2p=J2min; J2p<=J2max; J2p+=2) 
            rm(J2-J2min,J2p-J2min) = pow(-1.,(S2+L2+J2p)/2.+1.) * sqrt((S2+1.)*(J2+1.)*(J2p+1.)*(S2/2.)*(S2/2.+1.)) * sixj(S2,J2,L2,J2p,S2,2);
      else                                              // Lx, Ly or Lz
         for (int J2p=J2min; J2p<=J2max; J2p+=2) 
            rm(J2-J2min,J2p-J2min) = pow(-1.,(S2+L2+J2)/2.+1.)  * sqrt((L2+1.)*(J2+1.)*(J2p+1.)*(L2/2.)*(L2/2.+1.)) * sixj(L2,J2,S2,J2p,L2,2);
//    for (int J2p=J2min; J2p<=J2max; J2p+=2) 
//    {
//       Lrm = pow(-1.,(S2+L2+J2)/2.+1.)  * sqrt((L2+1.)*(J2+1.)*(J2p+1.)*(L2/2.)*(L2/2.+1.)) * sixj(L2,J2,S2,J2p,L2,2);
//       Srm = pow(-1.,(S2+L2+J2p)/2.+1.) * sqrt((S2+1.)*(J2+1.)*(J2p+1.)*(S2/2.)*(S2/2.+1.)) * sixj(S2,J2,L2,J2p,S2,2);
//       rm(J2-J2min,J2p-J2min)  = -( Lrm + GS*Srm );
//    }
   }

   sMat<double> mu(ns,ns);
   int ins=0; std::vector<int> J2(ns), mJ2(ns), irm(ns);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2)
     for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2)
     {
        J2[ins]=iJ2; mJ2[ins]=imJ2; irm[ins]=iJ2-J2min; ins++;
     }

   double elm, sqrt2=sqrt(2.);
   if (ind<2)                                           // Sx or Lx
      for (int i=0; i<ns; i++) for(int j=0; j<ns; j++)
      {
         elm = rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],-2,mJ2[j]);
         elm-= rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],2,mJ2[j]);
         if(fabs(elm)>SMALL) mu(i,j)=elm/sqrt2;
       //mu(i,j) = (elm-elp)/sqrt2;
      }
   else if (ind<4)                                      // Sy or Ly
      for (int i=0; i<ns; i++) for(int j=0; j<ns; j++)
      {
         elm = rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],-2,mJ2[j]);
         elm+= rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],2,mJ2[j]);
         if(fabs(elm)>SMALL) mu(i,j)=elm/sqrt2;
       //mu(i,j) = (elm+elp)/sqrt2;
      }
   else                                                 // Sz or Lz
      for (int i=0; i<ns; i++) for(int j=0; j<ns; j++) 
      {
         elm = rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],0,mJ2[j]); 
         if(fabs(elm)>SMALL) mu(i,j)=elm; 
      }
// {
//    mu(i,j) = rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],2*q,mJ2[j]);
// }
   return mu;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to initialise the storage matrix "est" for mcalc
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
          void mcalc_parameter_storage_matrix_init(
                      ComplexMatrix *est, // Output Eigenstates matrix (row 0: real==Eigenvalues;imag==population)
                      Vector &gjmbheff,   // Input  Effective mean fields (meV)
 /* Not Used */       double *g_J,        // Input  Lande g-factor
                      double *T,          // Input  temperature
 /* Not Used */       Vector &ABC,        // Input  Vector of parameters from single ion property file
                      char **sipffilename)// Input  Single ion properties filename
{
   // Parses the input file for parameters
   icpars pars;
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

   // If we just want a blank estates matrix for later use (e.g. in mcalc)
   int nfact = (int)ceil(sqrt(gjmbheff.Elements()+1));
   int Hsz = icf_getdim(pars)*nfact;
   (*est) = ComplexMatrix(0,Hsz,0,Hsz);
   (*est)(0,0) = complex<double> (pars.n, pars.l);
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates thermal expectation values 
// --------------------------------------------------------------------------------------------------------------- //
void icf_expJ(icpars &pars, ComplexMatrix &est, complexdouble *zV, double *vE, double *T, Vector &J, double *lnZ, double *U)
{
   int K[] = {-1,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Hsz=icf_getdim(pars), i, ix, iy, incx=1;
   int nfact = (int)ceil(sqrt(J.Hi()-J.Lo()+2));

   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';
   int Esz, ind_j;
   //std::vector< std::vector<double> > matel;
   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   std::vector<double> E, me, eb; /*matel.clear();*/ E.reserve(Hsz);
   for(Esz=0; Esz<Hsz; Esz++) { E.push_back(vE[Esz]-vE[0]); if(exp(-E[Esz]/(KB**T))<DBL_EPSILON || vE[Esz+1]==0) break; }

   if (*T<0){Esz=(int)(-*T);printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                         printf ("Number   Excitation Energy\n");
     for (ind_j=0;ind_j<Esz;++ind_j) printf ("%i    %4.4g meV\n",ind_j+1,E[ind_j]);
     } // MR 10.9.2010

   for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++)
   {
      int ind = ii+Hsz*jj;
      if(fabs(zV[ind].r*zV[ind].r+zV[ind].i*zV[ind].i)<DBL_EPSILON*100000) 
      {
         zV[ind].r=0.; zV[ind].i=0.;  
      } 
   }

   // For first run calculate also the partition function and internal energy
   *U=0.;
   me.assign(Esz,0.); eb.assign(Esz,0.); double Z=0.; complexdouble *zt=0, zme;
   complexdouble *zJmat = (complexdouble*)malloc(Hsz*Hsz*sizeof(complexdouble));
   zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
   sMat<double> zeros(Hsz,Hsz),Hmat;

   for(int iJ=J.Lo(); iJ<=J.Hi(); iJ++)
   {
      me.assign(Esz,0.); J[iJ]=0.;
      // Using the above reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
      if((iJ>6 && K[iJ]%2==1) || (K[iJ]>4 && pars.l==D)) { /*matel.push_back(me);*/ continue; }
      {
         iy = (iJ-J.Lo()+1)/nfact; ix = (iJ-J.Lo()+1)-iy*nfact;
         for(i=1; i<=Hsz; i++) memcpy(&zJmat[(i-1)*Hsz],&est[i+ix*Hsz][1+iy*Hsz],Hsz*sizeof(complexdouble));

         for(ind_j=0; ind_j<Esz; ind_j++)
         {  // Calculates the matrix elements <Vi|J.H|Vi>
            F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, &zV[ind_j*Hsz], &incx, &zbeta, zt, &incx);
#ifdef _G77 
            F77NAME(zdotc)(&zme, &Hsz, &zV[ind_j*Hsz], &incx, zt, &incx);
#else
            zme = F77NAME(zdotc)(&Hsz, &zV[ind_j*Hsz], &incx, zt, &incx);
#endif
            me[ind_j] = zme.r;
            // For first run calculate also the partition function and internal energy
            if(iJ==J.Lo()) 
            {
//MR 10.9.2010
               if (*T<0)
               { char instr[MAXNOFCHARINLINE];
                 printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
                 if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { fprintf(stderr,"Error reading input\n"); exit(1); }
                 eb[ind_j]=strtod(instr,NULL);
               }
               else
               { eb[ind_j] = exp(-E[ind_j]/(KB**T)); } J[iJ]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; *U+=(E[ind_j]+vE[0])*eb[ind_j];
//MRend 10.9.2010
            }
            else
            // Rest of the runs only calculate the new matrix elements
               J[iJ]+=me[ind_j]*eb[ind_j];
         }
	 if(iJ==J.Lo()) *U/=Z;
         /*matel.push_back(me);*/ J[iJ]/=Z;
      }
      if(fabs(J[iJ])<DBL_EPSILON) J[iJ]=0.;
   } 
   free(zJmat); free(zt);
   *lnZ = log(Z)-vE[0]/(KB**T); 
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the magnetic moment at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void mcalc(Vector &J,          // Output single ion momentum vector <Ja>,<Jb>,<Jc>, etc.
                      double *T,          // Input scalar temperature
                      Vector &gjmbH,      // Input vector of mean fields (meV) 
 /* Not Used */       double *gJ,         // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      double *lnZ,        // Output scalar logarithm of partition function
                      double *U,          // Output scalar internal energy 
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   // Parses the input file for parameters
   icpars pars; 
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

   // Prints out single ion parameters and filename.
// std::cout << "#ic1ion: Single ion parameter file: " << filename << "\n";
// std::cout << "#ic1ion: Read in parameters (in " << pars.e_units << "): F2=" << pars.F[1] << ",F4=" << pars.F[2];
// if(pars.l==F) std::cout << ",F6=" << pars.F[3]; std::cout << ",xi=" << pars.xi << "\n";
// std::cout << "#ic1ion: \t" << pars.B.cfparsout(", ") << "in " << pars.B.units() << "\n";

   // Converts the Jij parameters if necessary
   std::vector<double> vgjmbH(J.Hi()+1,0.); 
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc();
      for(int i=J.Lo(); i<=J.Hi(); i++) vgjmbH[i] = -gjmbH[i]*pars.jijconv[i]; }
   else
   #endif
      for(int i=J.Lo(); i<=J.Hi(); i++) vgjmbH[i] = -gjmbH[i];

   int K[] = {-1,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Q[] = {-1,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {-1,0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};

   // Calculates the IC Hamiltonian matrix
   int nfact = (int)ceil(sqrt(J.Hi()-J.Lo()+2));
   int i,k,q,Hsz=icf_getdim(pars),esz=Hsz*nfact;
   int ix, iy, incx=1;
   complexdouble *H=0;
   bool Hicnotcalc = false;
   std::vector<double> parval; parval.reserve(35);
   parval.push_back(pars.xi); for(k=2; k<=(2*pars.l); k+=2) for(q=-k; q<=k; q++) parval.push_back(pars.B(k,q));
   if(parval.size()%2==1) parval.push_back(0.);

   if((est.Cols()!=(esz+1) || est.Rows()!=(esz+1))) Hicnotcalc = true;
   else if(real(est[0][0])==(double)pars.n && imag(est[0][0])==(double)pars.l)  // Hic previously calculated
   {
      for(i=0; i<(int)(parval.size()/2); i++) if(real(est[0][i+1])!=parval[2*i] || imag(est[0][i+1])!=parval[2*i+1]) 
      { Hicnotcalc = true; break; }
   }
   else Hicnotcalc = true;

   if(Hicnotcalc)
   {
      sMat<double> Hcfi, Hcf = icf_hmltn(Hcfi, pars); Hcf/=MEV2CM; Hcfi/=MEV2CM; H = zmat2f(Hcf,Hcfi);
      if(est.Rhi()!=esz||est.Chi()!=esz) {
         std::cerr << "ERROR module icf1ion - mcalc: Hsz recalculation does not agree with eigenstates matrix dimension\n"; exit(EXIT_FAILURE); }
      else if(esz>(int)parval.size()/2)
      {
         est[0][0] = complex<double> (pars.n,pars.l);
         for(i=0; i<(int)(parval.size()/2); i++) est[0][i+1] = complex<double> (parval[2*i],parval[2*i+1]);
      }
      // Copies Hcf to est
      for(i=1; i<=Hsz; i++) memcpy(&est[i][1],&H[(i-1)*Hsz],Hsz*sizeof(complexdouble)); free(H);
      Hcfi.zero(Hsz,Hsz);
      for(int ind=J.Lo(); ind<=J.Hi(); ind++)
      {
         if(ind<=6)         // Calculates Sx,Lx etc and copies them to est too. 
         {
            Hcf = icf_mumat(pars.n, ind-1, pars.l);
            if(ind==3 || ind==4) H = zmat2f(Hcfi,Hcf); else H = zmat2f(Hcf,Hcfi); 
            iy = (ind-J.Lo()+1)/nfact; ix = (ind-J.Lo()+1)-iy*nfact;
            for(i=1; i<=Hsz; i++) memcpy(&est[i+ix*Hsz][1+iy*Hsz],&H[(i-1)*Hsz],Hsz*sizeof(complexdouble)); free(H);
         }
         else               // Calculates multipolar operator matrices and copies them to est. 
         {
            Hcf = icf_ukq(pars.n,K[ind],Q[ind],pars.l);
            if(im[ind]==1) H = zmat2f(Hcfi,Hcf); else H = zmat2f(Hcf,Hcfi);
            iy = (ind-J.Lo()+1)/nfact; ix = (ind-J.Lo()+1)-iy*nfact;
            for(i=1; i<=Hsz; i++) memcpy(&est[i+ix*Hsz][1+iy*Hsz],&H[(i-1)*Hsz],Hsz*sizeof(complexdouble)); free(H);
         }
      }
   }

   // Calculates the mean field matrix from stored matrices
   H = (complexdouble*)malloc(esz*esz*sizeof(complexdouble));
   for(i=1; i<=Hsz; i++) memcpy(&H[(i-1)*Hsz],&est[i][1],Hsz*sizeof(complexdouble));
   for(int ind=J.Lo(); ind<=J.Hi(); ind++)
   {
      iy = ind/nfact; ix = ind-iy*nfact;
      complex<double> a(vgjmbH[ind],0.);
      for(i=1; i<=Hsz; i++) F77NAME(zaxpy)(&Hsz,(complexdouble*)&a,(complexdouble*)&est[i+ix*Hsz][1+iy*Hsz],&incx,&H[(i-1)*Hsz],&incx);
   }

   // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
   double *vE = new double[Hsz]; complexdouble *zV = new complexdouble[Hsz*Hsz];
   int info = ic_diag(Hsz,H,zV,vE); free(H);
   if(info!=0) { std::cerr << "icf1ion - Error diagonalising, info==" << info << "\n"; delete[]vE; vE=0; delete[]zV; zV=0; exit(EXIT_FAILURE); }

   icf_expJ(pars,est,zV,vE,T,J,lnZ,U); delete[]vE; delete[]zV;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the eigenstates of Hic+effective_mean_fields
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
          void estates(ComplexMatrix *est,// Output Eigenstates matrix (row 0: real==Eigenvalues;imag==population)
                      Vector &gjmbheff,   // Input  Effective mean fields (meV)
 /* Not Used */       double &g_J,        // Input  Lande g-factor
                      double &T,          // Input  temperature
 /* Not Used */       Vector &ABC,        // Input  Vector of parameters from single ion property file
                      char **sipffilename)// Input  Single ion properties filename
{
   clock_t start,end; start = clock();

   // Parses the input file for parameters
   icpars pars;
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

   int K[] = {-1,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Q[] = {-1,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {-1,0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
   int i,j;

   // Converts the Jij parameters if necessary
   std::vector<double> vgjmbH(gjmbheff.Hi()+1,0.);
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc();
      for(i=gjmbheff.Lo(); i<=gjmbheff.Hi(); i++) vgjmbH[i] = -gjmbheff[i]*pars.jijconv[i]; }
   else
   #endif
      for(i=gjmbheff.Lo(); i<=gjmbheff.Hi(); i++) vgjmbH[i] = -gjmbheff[i];

   // Calculates the Hamiltonian matrix
   sMat<double> mat, Hcfi, Hcf; Hcf = icf_hmltn(Hcfi, pars); Hcf/=MEV2CM; Hcfi/=MEV2CM;
   for(int ind=gjmbheff.Lo(); ind<=gjmbheff.Hi(); ind++)
   {
      if(fabs(vgjmbH[ind])<SMALL) continue; 
      if(ind<=6) mat = icf_mumat(pars.n, ind-1, pars.l); else mat = icf_ukq(pars.n,K[ind],Q[ind],pars.l); 
      mat *= vgjmbH[ind];
      if(im[ind]==1) Hcfi += mat; else Hcf += mat;
   }

   // Initialises the output matrix
   int Hsz = Hcf.nr();
   (*est) = ComplexMatrix(0,Hsz,0,Hsz);

   // Diagonalises the Hamiltonian H = Hcf + sum_a(gjmbH_a*Ja)
// iceig VE; if(Hcfi.isempty()) VE.calc(Hcfi); else VE.calc(Hcf,Hcfi);
   double *vE = new double[Hsz], *rV=0; complexdouble *zV=0; int info; 
   if(Hcfi.isempty()) {
      rV = new double[Hsz*Hsz]; info = ic_diag(Hcf,rV,vE); }
   else { 
      zV = new complexdouble[Hsz*Hsz]; info = ic_diag(Hcf,Hcfi,zV,vE); }
   if(info!=0) { std::cerr << "icf1ion - Error diagonalising, info==" << info << "\n"; 
      delete[]vE; vE=0; if(Hcfi.isempty()) delete[]rV; else delete[]zV; exit(EXIT_FAILURE); }

   // Stores the number of electrons and the orbital number in element (0,0)
   (*est)(0,0) = complex<double> (pars.n, pars.l);

   // Puts eigenvectors/values into the est matrix
   for(i=0; i<Hsz; i++) (*est)(0,i+1) = complex<double> (vE[i], exp(-(vE[i]-vE[0])/(KB*T)));   // Row 0

   if(!Hcfi.isempty()) // {
      for(i=1; i<=Hsz; i++) memcpy(&(*est)[i][1],&zV[(i-1)*Hsz],Hsz*sizeof(complexdouble)); 
//    std::cout << "\n\nest==VE.zV = " << checkmat((*est),VE.zV(0),1,1) << endl << endl; }
   else
      for(i=0; i<Hsz; i++){//printf("\n");
         for(j=0; j<Hsz; j++) 
            {(*est)(i+1,j+1) = complex<double> (rV[j+i*Hsz], 0.);
//                         printf("%6.3f %+6.3f i  ",rV[i+j*Hsz],0.0);
//           if(rV[i+j*Hsz]!=(*est)[i+1][j+1]){fprintf(stderr,"compiler problem: bad memory mapping of vectors\n");exit(EXIT_FAILURE);}
            }}

   if(Hcfi.isempty()) delete[]rV; else delete[]zV; delete[]vE;
   end = clock(); std::cerr << "Time to do estates() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";

}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           int dmcalc(int &tn,            // Input transition number; if tn>0, print debug info
                      double &T,          // Input temperature
                      Vector &gjmbH,      // Input vector of mean fields (meV) 
 /* Not Used */       double &g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &mat, // Output M or Q matrix
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ 
   int i,j,k;

   int K[] = {-1,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Q[] = {-1,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {-1,0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};

   int sz = gjmbH.Hi();
   sMat<double> zeroes(est.Rows(),est.Cols()), op, Mab(sz,sz), iMab(sz,sz);
   complexdouble *zJmat=0, *zt=0, zme; zme.r=0; zme.i=0.; 
   std::vector<complexdouble> zij(sz,zme);//, zji(6,zme);
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';
   double Z=0., therm;

   // check if printout should be done and make tn positive
   int pr=1; if (tn<0) { pr=0; tn*=-1; }

   // Copies the already calculated energy levels / wavefunctions from *est
   if(est.Rows()!=est.Cols()) { std::cerr << "dmcalc(): Input rows and columns of eigenstates matrix don't match.\n"; return 0; }
   int Hsz = est.Rows()-1, iJ, jJ, incx = 1;
   j=0; k=0; for(i=0; i<Hsz; ++i) { for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   if(est[0][j+1].real()-est[0][i+1].real()<delta)
   {
      double *en = new double[Hsz]; for(k=0; k<Hsz; k++) en[k] = est[0][k+1].real();
//    iceig VE(Hsz,en,(complexdouble*)&est[1][0],1);
//    for(int ii=1; ii<=Hsz; ii++) memcpy(&zV[(ii-1)*Hsz],&(*est)[ii][1],Hsz*sizeof(complexdouble)); 

      // Parses the input file for parameters
      icpars pars; const char *filename = sipffilename[0];
      ic_parseinput(filename,pars);

      // Calculates the transition matrix elements:
      //    M_ab = <i|Ja|j><j|Jb|i> * (exp(-Ei/kT)-exp(-Ej/kT)) / Z    if delta > small
      //    M_ab = <i|Ja|j><j|Jb|i> * (exp(-Ei/kT)) / kTZ              if delta < small (quasielastic scattering)
      for(iJ=0; iJ<sz; iJ++)
      {
         if(iJ<=6) op = icf_mumat(pars.n, iJ-1, pars.l); else op = icf_ukq(pars.n,K[iJ],Q[iJ],pars.l); 
         if(im[iJ]==1) zJmat=zmat2f(zeroes,op); else zJmat=zmat2f(op,zeroes);
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, (complexdouble*)&est[1][j], &incx, &zbeta, zt, &incx);
         #ifdef _G77 
         F77NAME(zdotc)(&zij[iJ], &Hsz, (complexdouble*)&est[1][i], &incx, zt, &incx);
         #else
         zij[iJ] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[1][i], &incx, zt, &incx);
         #endif
         free(zJmat); free(zt);
      }

      if(i==j) //subtract thermal expectation value from zij=zii
      {
         Vector vJ(gjmbH.Lo(),gjmbH.Hi()); double lnZ=0., U=0.;
         icf_expJ(pars,est,(complexdouble*)&est[1][0],en,&T,vJ,&lnZ,&U);
         for(iJ=0; iJ<sz; iJ++) zij[iJ].r-=vJ[iJ];
      }

      // Calculates the matrix M_ab and iM_ab
      for(iJ=0; iJ<sz; iJ++)
         for(jJ=0; jJ<sz; jJ++)
         {  
            Mab(iJ+1,jJ+1) = (zij[iJ].r*zij[jJ].r + zij[iJ].i*zij[jJ].i);
            iMab(iJ+1,jJ+1) = (-zij[iJ].r*zij[jJ].i + zij[iJ].i*zij[jJ].r);
         }

      delta = en[j]-en[i];
      if(delta<-0.000001) {
         std::cerr << "ERROR module ic1ion - dmcalc: energy gain delta gets negative\n"; exit(EXIT_FAILURE); }
      if(j==i)delta=-SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc
   
      // Calculates the partition function
      for(iJ=0; iJ<Hsz; iJ++) { therm = exp(-(en[iJ]-en[0])/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }
   
      // do some printout if wishes and set correct occupation factor
      if (delta>SMALL)
      {
         therm = exp(-(en[i]-en[0])/(KB*T)) - exp(-(en[j]-en[0])/(KB*T));
         if(pr==1)
         {
            printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
            printf(" |<%i|Ja|%i>|^2=%6.3f\n |<%i|Jb|%i>|^2=%6.3f\n |<%i|Jc|%i>|^2=%6.3f\n",i+1,j+1,Mab(1,1),i+1,j+1,Mab(2,2),i+1,j+1,Mab(3,3));
            printf(" |<%i|Jd|%i>|^2=%6.3f\n |<%i|Je|%i>|^2=%6.3f\n |<%i|Jf|%i>|^2=%6.3f\n",i+1,j+1,Mab(4,4),i+1,j+1,Mab(5,5),i+1,j+1,Mab(6,6));
            printf(" n%i-n%i=%6.3f\n",i,j,therm / Z);
         }
      }
      else
      {
         therm = exp(-(en[i]-en[0])/(KB*T))/(KB*T);    // quasielastic scattering has not wi-wj but wj*epsilon/kT
         if(pr==1)
         {
            printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
            printf(" |<%i|Ja-<Ja>|%i>|^2=%6.3f\n |<%i|Jb-<Jb>|%i>|^2=%6.3f\n |<%i|Jc-<Jc>|%i>|^2=%6.3f\n",i+1,j+1,Mab(1,1),i+1,j+1,Mab(2,2),i+1,j+1,Mab(3,3));
            printf(" |<%i|Jd-<Jd>|%i>|^2=%6.3f\n |<%i|Je-<Je>|%i>|^2=%6.3f\n |<%i|Jf-<Jf>|%i>|^2=%6.3f\n",i+1,j+1,Mab(4,4),i+1,j+1,Mab(5,5),i+1,j+1,Mab(6,6));
            printf(" n%i=%6.3f\n",i,therm/Z);
         }
      }
   
      // multiply matrix Mab by occupation factor
      for(iJ=1; iJ<=sz; iJ++)
         for(jJ=1; jJ<=sz; jJ++)
	    mat(iJ,jJ) = complex<double> ( Mab(iJ,jJ)*therm/Z, iMab(iJ,jJ)*therm/Z );

      delete[]en;
   }

   // determine number of thermally reachable states
   int noft=0;for(i=0;(i<Hsz)&(exp(-(est[0][i+1].real()-est[0][1].real())/(KB*T))>SMALL);++i)noft+=Hsz-i-1;
   return noft;
   //return Hsz*(Hsz-1)/2;
}

/*
// --------------------------------------------------------------------------------------------------------------- //
// Loads a Q_q matrix from file if the file exists and has the same parameters n,l,Jvec
// --------------------------------------------------------------------------------------------------------------- //
bool get_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l, std::vector<double> &Jvec)
{
   int i, j, mn, r, c, sz, ml;
   std::vector<double> mJv(6,0.); 
   char filename[] = "results/mcphas.Qq"; filename[16]=q+120;               // 120==x, 121==y, 122==z
   std::fstream FILEIN; FILEIN.open(filename, std::fstream::in);
   if(FILEIN.fail()==true) return false;
   FILEIN >> mn >> ml; for(i=0; i<6; i++) FILEIN >> mJv[i]; FILEIN >> r >> c;
   if(mn!=n || ml!=(int)l) return false; for(i=0; i<6; i++) if(fabs(mJv[i]-Jvec[i])>1e-4) return false;
   Qq.clear(); sMat<double> emptymat(r,c); double Qt;
   for(i=0; i<6; i++) 
   {
      Qq.push_back(emptymat); FILEIN >> sz; for(j=0; j<sz; j++) {  FILEIN >> r >> c >> Qt; Qq[i](r,c) = Qt; }
   }
   return true;
}

// --------------------------------------------------------------------------------------------------------------- //
// Saves a Q_q matrix to a temporary file in the results/ directory
// --------------------------------------------------------------------------------------------------------------- //
void save_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l, std::vector<double> &Jvec)
{
   int i,j,sz;
   std::vector< std::vector<int> > nz;
   char filename[] = "results/mcphas.Qq"; filename[16]=q+120;               // 120==x, 121==y, 122==z
   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out);
   FILEOUT << n << " " << (int)l << " "; for(i=0; i<6; i++) FILEOUT << Jvec[i] << " "; FILEOUT << "\n";
   FILEOUT << Qq[0].nr() << " " << Qq[0].nc() << " ";
   for(i=0; i<6; i++)
   {
      nz = Qq[i].find(); sz = nz.size(); FILEOUT << sz << "\n";
      for(j=0; j<sz; j++) FILEOUT << nz[j][0] << " " << nz[j][1] << " " << Qq[i](nz[j][0],nz[j][1]) << "\n";
   }
   FILEOUT.close();
}
*/

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the a(K,K') matrix - the routine returns a true if the matrix is real, and false if it is imaginary
//    The matrix is returned in the argument aKK
// --------------------------------------------------------------------------------------------------------------- //
bool icf_loveseyAKK(sMat<double> &aKK, int K, int Kp, int n, orbital l)
{
   if (l>3||l<1) { std::cerr << "Sorry only p-, d- and f-electrons supported at present\n"; exit(0); }

   fstates_t gs = hunds_gs(n, l);
   int J2min, J2max, ns=0, S2 = gs.S2, L2 = abs(gs.L)*2; std::vector<int> J2;

   // Determines number and angular momentum quantum numbers of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2;
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) J2.push_back(iJ2); ns=(int)J2.size();

   // Loads a previously save matrix if it exists
/* char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
   nstr[0] = (l==F?102:100); if(n<10) { nstr[1] = n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 65;   // 65 is ASCII for "A", 100=="d" and 102=="f"
   nstr[1] = K+48; nstr[2] = Kp+48; nstr[3] = 0; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
   aKK = mm_gin(filename); if(!aKK.isempty()) return 1;

   bool df;
   if(l==D) df = true;
   else if(l==F) df = false;
   else {  std::cerr << "racah_mumat(): Only d- and f- configurations are implemented.\n"; return false; }

   int np = (n==1)?n:(n-1); 
   fconf confp(np,l);
   fconf conf(n,l);
   int num_states = (int)conf.states.size(), ns=0;
   int i,j,k,kk,v,vp,isz,jsz,iJ,jJ;
   int j2min,j2max,j2pmin,j2pmax;
   int L2,L2p,S2,S2p,J2,J2p,L2b,S2b;
   std::vector<int> indexJstart,indexJstop;
   std::vector<cfpls> cfpsi,cfpsj;
   double rmLS,sumcfp;

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      indexJstart.push_back(ns+1);
      for(j=j2min; j<=j2max; j+=2) ns++;
      indexJstop.push_back(ns);
   }
*/

   aKK.zero(ns,ns);

   // Selection rules on K.
   if(K<0  || K>(2*abs(l))   || (K%2)!=0)      return 1;        // Eqn 4.2.2
   if(Kp<1 || Kp>=(2*abs(l)) || ((Kp+1)%2)!=0) return 1;        // Eqn 4.2.3
   if(K!=(Kp+1) && K!=(Kp-1))                  return 1;        // Eqn 4.2.4  (first 3j symbol, needs 1+K+K' even)

   // Calculates the matrix elements:
   // Eqn. 3.6.8 (d==delta)
   //
   //                                  K'+1           J'+S+L+L'+l  1/2    2                  1/2
   // A(K,K') = ( <j    > + <j    > ) i      d    (-1)            2    [l]  [L,L',K',K',K,J']
   //               K'-1      K'+1            SS'
   //
   //              \/   ( 1 K K' )  { 1  1  1  }  { J' K' J  } 
   //              /\   ( 0 0 0  )  { K' K  K' }  { L  S  L' }  A(K',K',l)
   //                                            _
   //                     ---     _   _          L                                _
   //              \/   n >   (t{|t) (t|}t') (-1)  { L' K' L }                Lb==L
   //              /\     ---                      { l  Lb l }
   //                      t

   // Where: (Eqn 3.4.9)
   //                                    1/2
   //                 l+1 [ (2l+3)(l+1) ]    ( l K' l+1 )  { 1  K'  K }
   // A(K,K',l) = (-1)    [ ----------- ]    ( 0 0  0   )  { l l+1  l }
   //                     [   (2l+1)    ]
   //
   
   // Calculate the coefficient A(K,K',l)
   double AKKl = pow(-1.,l+1.) * sqrt( (2.*l+3.)*(l+1.)/(2.*l+1.) ) * threej(2*l,2*Kp,2*(l+1),0,0,0) * sixj(2,2*Kp,2*Kp,2*l,2*(l+1),2*l);

   // Calculate the non-state dependent part of the matrix elements
   double redmat = sqrt(2) * (2*l+1)*(2*l+1) * (2*Kp+1.) * sqrt(2*K+1.) * threej(2,2*K,2*Kp,0,0,0) * sixj(2,2,2,2*Kp,2*K,2*Kp) * AKKl;

   // Determines the factor i^(K'-1) - because K' is always odd, matrix is always real
   if((Kp-1)%4==0) redmat = -redmat;

   double rmLS=0.,sumcfp=0.;
   if (n>1)
   {
      fconf confp(n-1,l);
      std::vector<cfpls> cfps;
      switch(l) {
         case P:  cfps = racah_parents(n,gs.S2,gs.L); break;
         case D:  cfps = racah_parents(n,gs.v,gs.S2,gs.L); break;
         default: cfps = racah_parents(n,gs.v,gs.U,gs.S2,gs.L);  }
      int sz = (int)cfps.size(), pL2;
      for(int k=0; k<sz; k++)
      {
         pL2 = abs(confp.states[cfps[k].ind].L)*2;
         sumcfp += cfps[k].cfp*cfps[k].cfp * pow(-1.,pL2/2) * sixj(L2,2*Kp,L2,2*l,pL2,2*l);
      }
      rmLS = (L2+1.) * n * sumcfp * redmat;
   }
   else  // Single electron
   {
      rmLS = -redmat;
   }

   for(int i=0; i<ns; i++) for(int j=0; j<ns; j++)
      aKK(i,j) = pow(-1.,(J2[j]+S2+L2+L2)/2.+l) * sqrt(J2[j]+1.) * sixj(J2[j],2*Kp,J2[i],L2,S2,L2) * rmLS;

   // Calculates the matrix elements at particular |J,M>, |J',M'>
/* for(i=0; i<num_states; i++)
   {
      L2 = abs(conf.states[i].L)*2; S2 = conf.states[i].S2; v = conf.states[i].v;
      if(n!=1) 
      {
         if(df) cfpsi = racah_parents(n,v,S2,conf.states[i].L); 
         else   cfpsi = racah_parents(n,v,conf.states[i].U,S2,conf.states[i].L);
      }

      for(j=0; j<num_states; j++)
      {
         L2p = abs(conf.states[j].L)*2; S2p = conf.states[j].S2; vp = conf.states[j].v;

         if(S2!=S2p) continue;                                  // delta_SS'
         if(abs(L2-L2p)>(2*Kp)) continue;                       // Triangular condition on 6j symbol in 3.6.9

	 if(n==1) 
            rmLS = -redmat;
         else
         {
            if(df) cfpsj = racah_parents(n,vp,S2p,conf.states[j].L); 
            else   cfpsj = racah_parents(n,vp,conf.states[j].U,S2p,conf.states[j].L);

            sumcfp = 0.; isz = (int)cfpsi.size(); jsz = (int)cfpsj.size();
            for(k=0; k<isz; k++)
               for(kk=0; kk<jsz; kk++)
                  if(cfpsi[k].ind==cfpsj[kk].ind) {
                     L2b = 2*abs(confp.states[cfpsi[k].ind].L); S2b = confp.states[cfpsi[k].ind].S2;
                     sumcfp += cfpsi[k].cfp*cfpsj[kk].cfp * pow(-1.,L2b/2) * sixj(L2p,2*Kp,L2,2*l,L2b,2*l); }
            rmLS = sqrt( (L2+1.)*(L2p+1.) ) * n * sumcfp * redmat;
         }

         // Caculate the J-dependent reduced matrix elements 
         j2min = abs(L2-S2); j2max = L2+S2; j2pmin = abs(L2p-S2p); j2pmax = L2p+S2p;
         iJ=indexJstart[i]-1; jJ=indexJstart[j]-1;
         for(J2=j2min; J2<=j2max; J2+=2)
         {
            for(J2p=j2pmin; J2p<=j2pmax; J2p+=2)
            {
               if(abs(J2-J2p)<=(2*Kp))                          // Triangular condition on 6j symbol in 3.6.9
                  aKK(iJ,jJ) = pow(-1.,(J2p+S2+L2+L2p)/2.+l) * sqrt(J2p+1.) * sixj(J2p,2*Kp,J2,L2,S2,L2p) * rmLS;
               jJ++;
            }
            iJ++; jJ=indexJstart[j]-1;
         }
      }
   }
*/
// std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out); FILEOUT.close();
// rmzeros(aKK); mm_gout(aKK,filename); return 1;
   rmzeros(aKK); return 1;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the c(K,K') matrix, returning true if matrix is real, false if imaginary
// --------------------------------------------------------------------------------------------------------------- //
bool icf_loveseyCKK(sMat<double> &cKK, int K, int Kp, int n, orbital l)
{
   if (l>3||l<1) { std::cerr << "Sorry only p-, d- and f-electrons supported at present\n"; exit(0); }

   fstates_t gs = hunds_gs(n, l);
   int J2min, J2max, ns=0, S2 = gs.S2, L2 = abs(gs.L)*2; std::vector<int> J2;

   // Determines number and angular momentum quantum numbers of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2;
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) J2.push_back(iJ2); ns=(int)J2.size();

/* // Loads a previously save matrix if it exists
   char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
   nstr[0] = (l==F?102:100); if(n<10) { nstr[1] = n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 67;   // 67 is ASCII for "C", 100=="d" and 102=="f"
   nstr[1] = K+48; nstr[2] = Kp+48; nstr[3] = 0; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
   cKK = mm_gin(filename); if(!cKK.isempty()) return 1;

   bool df;
   if(l==D) df = true;
   else if(l==F) df = false;
   else { std::cerr << "racah_mumat(): Only d- and f- configurations are implemented.\n"; return false; }
   int np = (n==1)?n:(n-1);
   fconf confp(np,l);
   fconf conf(n,l);
   int num_states = (int)conf.states.size(), ns=0;
   int i,j,k,kk,j2min,j2max,j2pmin,j2pmax;
   int L2,L2p,S2,S2p,J2,J2p;
   std::vector<cfpls> cfpsi,cfpsj;
   std::vector<int> indexJstart,indexJstop,irm;
   sMat<double> mJmat_i(1,1);
   sMat<double> rmJ;
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   double rmLS,sumcfp;
   int v,vp,isz,jsz,L2b,S2b,iJ,jJ; 

   // Eqn. 3.6.11
   //
   //                  K      1/2                            1/2      1/2+S'+L'
   // C(K,K') = <j >  i  (1/2)     [l,l,S,S',L,L',J',K,K',K']     (-1)
   //             K                                                       _ _
   //                              { 1  K  K' }    ---     _   _          S+L
   //              \/   ( l K l )  { S' L' J' }  n >   (t{|t) (t|}t') (-1)    { S  1  S' }  { L  K  L' }
   //              /\   ( 0 0 0 )  { S  L  J  }    ---                        { s  Sb s  }  { l  Lb l  }
   //                                               t

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      indexJstart.push_back(ns+1);
      for(j=j2min; j<=j2max; j+=2) { ns++; irm.push_back(i); }
      indexJstop.push_back(ns);
   }
*/
   cKK.zero(ns,ns);

   // The triangular conditions require that K<=2l and K=even (from 3j), and that (K-1)<=K'<=(K+1) (from 9j)
   if(K%2!=0 || abs(Kp-K)>1) return 1;

   // Calculates part of the prefactor which only depends on K,K',l
   double redmat = sqrt(1./2) * (2*l+1.) * sqrt(2*K+1.) * (2*Kp+1.) * threej(2*l,2*K,2*l,0,0,0); 
   if(fabs(redmat)<DBL_EPSILON) return 1;

   // Determines the factor i^K - because K is always even, matrix is always real
   if(K%4==2) redmat = -redmat;

   double rmLS=0.,sumcfp=0.;
   if (n>1)
   {
      fconf confp(n-1,l);
      std::vector<cfpls> cfps;
      switch(l) {
         case P:  cfps = racah_parents(n,gs.S2,gs.L); break;
         case D:  cfps = racah_parents(n,gs.v,gs.S2,gs.L); break;
         default: cfps = racah_parents(n,gs.v,gs.U,gs.S2,gs.L);  }
      int sz = (int)cfps.size(), pL2, pS2;
      for(int k=0; k<sz; k++)
      {
         pS2 = confp.states[cfps[k].ind].S2; pL2 = abs(confp.states[cfps[k].ind].L)*2;
         sumcfp += cfps[k].cfp*cfps[k].cfp * pow(-1.,(1+S2+L2+pL2+pS2)/2.) * sixj(S2,2,S2,1,pS2,1)*sixj(L2,2*K,L2,2*l,pL2,2*l); 
      }
      rmLS = (L2+1.)*(S2+1.) * n * sumcfp * redmat;
   }
   else  // Single electron
   {
      rmLS = redmat;
   }

   for(int i=0; i<ns; i++) for(int j=0; j<ns; j++)
      if(abs(J2[i]-J2[j])<=(2*Kp))                          // Triangular condition on 9j symbol in 3.6.11
         cKK(i,j) = sqrt(J2[j]+1.) * ninej(2,2*K,2*Kp,S2,L2,J2[j],S2,L2,J2[i]) * rmLS;
/*
   rmJ.zero(num_states,num_states);
   // Calculates the matrix elements at particular |J,M>, |J',M'>
   for(i=0; i<num_states; i++)
   {
      L2 = abs(conf.states[i].L)*2; S2 = conf.states[i].S2; v = conf.states[i].v;
      if(n!=1) 
      {
         if(df) cfpsi = racah_parents(n,v,S2,conf.states[i].L); 
         else   cfpsi = racah_parents(n,v,conf.states[i].U,S2,conf.states[i].L); 
      }

      for(j=0; j<num_states; j++)
      {
         L2p = abs(conf.states[j].L)*2; S2p = conf.states[j].S2; vp = conf.states[j].v;

         if(abs(S2-S2p)>2)     continue;                        // Triangular condition on 6j symbol in 3.6.11
         if(abs(L2-L2p)>(2*K)) continue;                        // Triangular condition on 6j symbol in 3.6.11

	 if(n==1)
	    rmLS = redmat;
	 else
	 {
            if(df) cfpsj = racah_parents(n,vp,S2p,conf.states[j].L); 
            else   cfpsj = racah_parents(n,vp,conf.states[j].U,S2p,conf.states[j].L);

            sumcfp = 0.; isz = (int)cfpsi.size(); jsz = (int)cfpsj.size();
            for(k=0; k<isz; k++)
               for(kk=0; kk<jsz; kk++)
                  if(cfpsi[k].ind==cfpsj[kk].ind) {
                     L2b = 2*abs(confp.states[cfpsi[k].ind].L); S2b = confp.states[cfpsi[k].ind].S2;
                     sumcfp += cfpsi[k].cfp*cfpsj[kk].cfp * pow(-1.,(1+S2p+L2p+L2b+S2b)/2.) 
                                  * sixj(S2,2,S2p,1,S2b,1)*sixj(L2,2*K,L2p,2*l,L2b,2*l); }

            rmLS = sqrt( (S2+1.)*(S2p+1.)*(L2+1.)*(L2p+1.) ) * n * sumcfp * redmat;
         }

         // Caculate the J-dependent reduced matrix elements 
         j2min = abs(L2-S2); j2max = L2+S2; j2pmin = abs(L2p-S2p); j2pmax = L2p+S2p;
         iJ=indexJstart[i]-1; jJ=indexJstart[j]-1;
         for(J2=j2min; J2<=j2max; J2+=2)
         {
            for(J2p=j2pmin; J2p<=j2pmax; J2p+=2)
            {
               if(abs(J2-J2p)<=(2*Kp))                          // Triangular condition on 9j symbol in 3.6.11
                  cKK(iJ,jJ) = sqrt(J2p+1.) * ninej(2,2*K,2*Kp,S2p,L2p,J2p,S2,L2,J2) * rmLS;
               jJ++;
            }
            iJ++; jJ=indexJstart[j]-1;
         }
      }
   }

   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out); FILEOUT.close();
   rmzeros(cKK); mm_gout(cKK,filename); return 1;
*/
   rmzeros(cKK); return 1;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the transition operator Q_q
// --------------------------------------------------------------------------------------------------------------- //
void icf_loveseyQq(std::vector< sMat<double> > &Qq, int q, int n, orbital l, std::vector<double> &Jvec)
{
   if (l>3||l<1) { std::cerr << "Sorry only p-, d- and f-electrons supported at present\n"; exit(0); }

   double theta = Jvec[0], phi = Jvec[1], J[]={Jvec[2],0,Jvec[3],0,Jvec[4],0,Jvec[5]};
   std::string errormsg("lovesey_Qq(): Unable to calculate the A(K,K') or B(K,K') matrix\n");

   fstates_t gs = hunds_gs(n, l);
   int J2min, J2max, ns=0;
   int K,Kp,Q,Qp,i,j;
   sMat<double> Akk,Bkk,ckk,Qmat,QLmat,QSmat;

   bool iA,ic;
   complexdouble Ykq; double Tj,Tj2;

/* int k,minJ2,maxJ2,valJ,valJ_i,valJ_j;
   int imJ_i,imJ_j,J2,J2p,Jz2,Jz2p;
   std::vector< std::vector<int> > nz;

   sMat<double> mJmat_i(1,1);
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   std::vector<int> iJst,iJen;
*/
   // Determines number of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2; for (int J2=J2min; J2<=J2max; J2+=2) ns+=J2+1; 
   int ins=0; std::vector<int> J2(ns), Jz2(ns), irm(ns);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2) {
        J2[ins]=iJ2; Jz2[ins]=imJ2; irm[ins]=iJ2-J2min; ins++; }

   Qmat.zero(ns,ns); Qq.clear(); for(i=0; i<6; i++) Qq.push_back(Qmat);

/* // Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
   //    each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
   valJ = (J2max-J2min)/2;
   for(i=0; i<=valJ; i++) mJmat_row.push_back(mJmat_i);
   for(i=0; i<=valJ; i++) mJmat.push_back(mJmat_row);
*/
   for(Kp=0; Kp<=(2*l+1); Kp++)                                 // K may take even values between [0,2l]; K' odd values [1,(2l+1)]
   {
      if(Kp%2==1)                                               // Calculates the matrices for the first term with A+B
      {
         iA = icf_loveseyAKK(Akk,Kp-1,Kp,n,l); if(iA) Akk*= (J[Kp-1]+J[Kp+1])*((2*Kp+1)/(Kp+1.)); else { std::cerr << errormsg; return; }
         ic = icf_loveseyCKK(ckk,Kp-1,Kp,n,l); if(ic) Bkk = ckk*J[Kp-1];                          else { std::cerr << errormsg; return; }
         ic = icf_loveseyCKK(ckk,Kp+1,Kp,n,l); if(ic) Bkk+= ckk*(J[Kp+1]*sqrt(Kp/(Kp+1.)));       else { std::cerr << errormsg; return; }
         Bkk *= ((Kp+1.)/(2*Kp+1.)) * ( (2*Kp+1)/(Kp+1.) ); K = Kp-1; ckk = Bkk+Akk;
      }
      else                                                      // Calculates the matrices for the second term with just B
      {
         ic = icf_loveseyCKK(ckk,Kp,Kp,n,l); if(ic) ckk *= J[Kp]; else { std::cerr << errormsg; return; } K = Kp;
      }

      for(Q=-K; Q<=K; Q++)
      {
         Qp = -(Q-q); if(Qp<-Kp || Qp>Kp) continue;             // 3j symbol requires Q+Q'-q = 0
 
         Ykq = spherical_harmonics(K,Q,theta,phi); 
         Tj = pow(-1.,K-Kp+q) * sqrt(3) * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q);
         if((fabs(Ykq.r)<DBL_EPSILON && fabs(Ykq.i)<DBL_EPSILON) || fabs(Tj)<DBL_EPSILON) continue;

         Qmat.zero(ns,ns); QLmat.zero(ns,ns); QSmat.zero(ns,ns);

/*       nz = ckk.find();
         for(i=0; i<(int)nz.size(); i++)                        // The matrices above already contain the dependence on |vULSJ>
         {                                                      //    we now use the W-E theorem to add the Jz, Q, dependence
            J2 = conf.states[nz[i][0]].J2;
            J2p = conf.states[nz[i][1]].J2;
            valJ_i = (J2-minJ2)/2; valJ_j = (J2p-minJ2)/2;
            if(mJmat[valJ_i][valJ_j].isempty())
            {
               mJmat[valJ_i][valJ_j].zero(J2+1,J2p+1);
               for(imJ_i=0; imJ_i<=J2; imJ_i++)
               {
                  Jz2 = imJ_i*2-J2;
                  for(imJ_j=0; imJ_j<=J2p; imJ_j++)
                  {
                     Jz2p = imJ_j*2-J2p; Tj2 = threej(2*Kp,J2p,J2,2*Qp,Jz2p,-Jz2);
                     if(fabs(Tj2)>DBL_EPSILON) mJmat[valJ_i][valJ_j](imJ_i,imJ_j) = pow(-1.,Kp+(-J2p+Jz2)/2.) * sqrt(J2+1.) * Tj2;
                  }
               }
            }
            Qmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*ckk(nz[i][0],nz[i][1]));
            if(Kp%2==1)
            {
               QLmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*Akk(nz[i][0],nz[i][1]));
               QSmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*Bkk(nz[i][0],nz[i][1]));
            }
         }
         for(i=0; i<=valJ; i++) for(j=0; j<=valJ; j++) mJmat[i][j].clear();
*/
         for(i=0; i<ns; i++) for(j=0; j<ns; j++)
         {
            Tj2 = threej(2*Kp,J2[j],J2[i],2*Qp,Jz2[j],-Jz2[i]); if(fabs(Tj2)<DBL_EPSILON) continue;
            Tj2 *= pow(-1.,Kp+(-J2[j]+Jz2[i])/2.) * sqrt(J2[i]+1.);
            if(Kp%2==1)
            {
               QLmat(i,j) = Tj2 * Akk(irm[i],irm[j]);
               QSmat(i,j) = Tj2 * Bkk(irm[i],irm[j]);
            }
            Qmat(i,j) = Tj2 * ckk(irm[i],irm[j]);
         }

         Qq[0] += Qmat*(Ykq.r*Tj); Qq[1] += Qmat*(Ykq.i*Tj);
         if(Kp%2==1) { Qq[2] += QSmat*(Ykq.r*Tj); Qq[3] += QSmat*(Ykq.i*Tj);  Qq[4] += QLmat*(Ykq.r*Tj); Qq[5] += QLmat*(Ykq.i*Tj); }
         else        { Qq[2] += Qmat*(Ykq.r*Tj); Qq[3] += Qmat*(Ykq.i*Tj); }
      }
   }

   for(i=0; i<6; i++) { rmzeros(Qq[i]); Qq[i] *= sqrt(4*PI); }

/* 
   // Determines number and angular momentum quantum numbers of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2;
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) J2.push_back(iJ2); ns=(int)J2.size();

   sMat<double> aKK(ns,ns), cKK(ns,ns);

   // Selection rules on K for a(K,K')
   if(K<0  || K>(2*abs(l))   || (K%2)!=0)      return 1;        // Eqn 4.2.2
   if(Kp<1 || Kp>=(2*abs(l)) || ((Kp+1)%2)!=0) return 1;        // Eqn 4.2.3
   if(K!=(Kp+1) && K!=(Kp-1))                  return 1;        // Eqn 4.2.4  (first 3j symbol, needs 1+K+K' even)
   // The triangular conditions require that K<=2l and K=even (from 3j), and that (K-1)<=K'<=(K+1) (from 9j) for c(K,K')
   if(K%2!=0 || abs(Kp-K)>1) return 1;

   // Calculate the non-state dependent part of the matrix elements
   double AKKl = pow(-1.,l+1.) * sqrt( (2.*l+3.)*(l+1.)/(2.*l+1.) ) * threej(2*l,2*Kp,2*(l+1),0,0,0) * sixj(2,2*Kp,2*Kp,2*l,2*(l+1),2*l);
   double Aredmat = sqrt(2) * (2*l+1)*(2*l+1) * (2*Kp+1.) * sqrt(2*K+1.) * threej(2,2*K,2*Kp,0,0,0) * sixj(2,2,2,2*Kp,2*K,2*Kp) * AKKl;
   double Credmat = sqrt(1./2) * (2*l+1.) * sqrt(2*K+1.) * (2*Kp+1.) * threej(2*l,2*K,2*l,0,0,0); 

   if(fabs(Aredmat)>DBL_EPSILON || fabs(Credmat)>DBL_EPSILON)
   {
      // Determines the factor i^(K'-1) - because K' is always odd, matrix is always real
      if((Kp-1)%4==0) Aredmat = -Aredmat;
      // Determines the factor i^K - because K is always even, matrix is always real
      if(K%4==2) Credmat = -Credmat;

      double ArmLS=0., Asumcfp=0., CrmLS=0., Csumcfp=0.;
      if (n>1)
      {
         fconf confp(n-1,l);
         std::vector<cfpls> cfps;
         switch(l) {
            case P:  cfps = racah_parents(n,gs.S2,gs.L); break;
            case D:  cfps = racah_parents(n,gs.v,gs.S2,gs.L); break;
            default: cfps = racah_parents(n,gs.v,gs.U,gs.S2,gs.L);  }
         int sz = (int)cfps.size(), pL2, pS2;
         for(int k=0; k<sz; k++)
         {
            pS2 = confp.states[cfps[k].ind].S2; pL2 = abs(confp.states[cfps[k].ind].L)*2;
            Asumcfp += cfps[k].cfp*cfps[k].cfp * pow(-1.,pL2/2) * sixj(L2,2*Kp,L2,2*l,pL2,2*l);
            Csumcfp += cfps[k].cfp*cfps[k].cfp * pow(-1.,(1+S2+L2+pL2+pS2)/2.) * sixj(S2,2,S2,1,pS2,1)*sixj(L2,2*K,L2,2*l,pL2,2*l); 
         }
         ArmLS = (L2+1.) * n * Asumcfp * Aredmat; CrmLS = (L2+1.)*(S2+1.) * n * sumcfp * redmat;
      }
      else  // Single electron
      {
         ArmLS = -redmat; CrmLS = redmat;
      }

      for(int i=0; i<ns; i++) for(int j=0; j<ns; j++)
      {
         if(abs(J2[i]-J2[j])<=(2*Kp))                          // Triangular condition on 9j symbol in 3.6.11
            cKK(i,j) = sqrt(J2[j]+1.) * ninej(2,2*K,2*Kp,S2,L2,J2[j],S2,L2,J2[i]) * CrmLS; 
         aKK(i,j) = pow(-1.,(J2[j]+S2+L2+L2)/2.+l) * sqrt(J2[j]+1.) * sixj(J2[j],2*Kp,J2[i],L2,S2,L2) * ArmLS;
      }
   }

   std::vector<int> mJ2, irm; J2.clear(); mJ2.reserve(100); irm.reserve(100); J2.reserve(100);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2) {
      J2.push_back(iJ2); mJ2.push_back(imJ2); irm.push_back(iJ2-J2min); }
*/
/*
   int i,j,ns=getdim(n,l);                                      // Number of states = ^{4l+2}C_{n}

   std::string errormsg("lovesey_Qq(): Unable to calculate the A(K,K') or B(K,K') matrix\n");
   double theta = Jvec[0], phi = Jvec[1], J[]={Jvec[2],0,Jvec[3],0,Jvec[4],0,Jvec[5]};
   int K,Kp,Q,Qp;
   sMat<double> Akk,Bkk,ckk,Qmat,QLmat,QSmat;
   Qmat.zero(ns,ns); Qq.clear(); for(i=0; i<6; i++) Qq.push_back(Qmat);

   bool iA,ic;
   complexdouble Ykq; double Tj,Tj2;

   int k,minJ2,maxJ2,valJ,valJ_i,valJ_j;
   int imJ_i,imJ_j,J2,J2p,Jz2,Jz2p;
   std::vector< std::vector<int> > nz;

   sMat<double> mJmat_i(1,1);
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   std::vector<int> iJst,iJen;
   fconf conf(n,0,l); minJ2=99; maxJ2=0; k=0;
   for(i=0; i<(int)conf.states.size(); i++) {
      j = conf.states[i].J2; if(j<minJ2) minJ2 = j; if(j>maxJ2) maxJ2 = j; iJst.push_back(k+1); k+=j+1; iJen.push_back(k); }

   // Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
   //    each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
   valJ = (maxJ2-minJ2)/2;
   for(i=0; i<=valJ; i++)
      mJmat_row.push_back(mJmat_i);
   for(i=0; i<=valJ; i++)
      mJmat.push_back(mJmat_row);

   // The scattering operator is given by: (eqn. 3.7.9 and 3.7.11 for the orbital and spin parts respectively)
   //
   //              ____ ---  {  K'-1 ^    2K'+1                                                  K'                      }
   // <i|Q |j> = \/4*pi >    { Y    (k) [ ------ ] [ A(K'-1,K') + B(K'-1,K') ] ( K'-1 K' 1  ) + Y  B(K',K') ( K' K' 1  ) }
   //     q             ---  {  Q    -     K'+1                                ( Q    Q' -q )    Q          ( Q  Q' -q ) }
   //                  K'QQ'
   //                     -J'+q+M  ( K' J'  J )
   //                 (-1)      *  ( Q' M' -M )                  with |i>==|avUSLJM> and |j>==|a'v'U'S'L'J'M'>
   //
   // The matrices A(K,K') and B(K,K') are:                                                              1/2 
   //                                                                K+2 [                        ( K+1 )            ]
   //  A(K,K') = ( <j    > + <j    > ) a(K,K')           B(K,K+1) = ---- [ <j > c(K,K+1) + <j   > ( --- ) c(K+2,K+1) ]
   //                K'-1      K'+1                                 2K+3 [   K               K+2  ( K+2 )            ]
   //
   //                                                                                                    1/2
   //                                                                K-1 [                        (  K  )            ]
   //  B(K,K) = <j > c(K,K)                              B(K,K-1) = ---- [ <j > c(K,K-1) + <j   > ( --- ) c(K-2,K-1) ]
   //             K                                                 2K-1 [   K               K-2  ( K-1 )            ]
   //
   //
   //  The 3j symbols with the matrices a(K,K') and c(K,K') restrict K to even integers up to 2l, and K' to odd integers
   //  less than 2l. I.e. for f-electrons, K=0,2,4,6 and K'=1,3,5. The summation over K' in the first case however is over
   //  all integers up to 2l, so that in the case of even K', the first term (with the A+B) is zero and in the case of odd
   //  K', the second term (with just the B term) is zero.
   //
   //  The matrices a(K,K') and c(K,K') are solely functions of |i> and |j> and do not depend on the bessel functions
   //  <j_k>, or spherical harmonics Y^K_Q, and are calculated in the functions lovesey_aKK() and loveset_cKK() above.

   for(Kp=0; Kp<=(2*l+1); Kp++)                                 // K may take even values between [0,2l]; K' odd values [1,(2l+1)]
   {
      if(Kp%2==1)                                               // Calculates the matrices for the first term with A+B
      {
         iA = lovesey_aKK(Akk,Kp-1,Kp,n,l); if(iA) Akk*= (J[Kp-1]+J[Kp+1])*((2*Kp+1)/(Kp+1.)); else { std::cerr << errormsg; return; }
         ic = lovesey_cKK(ckk,Kp-1,Kp,n,l); if(ic) Bkk = ckk*J[Kp-1];                          else { std::cerr << errormsg; return; }
         ic = lovesey_cKK(ckk,Kp+1,Kp,n,l); if(ic) Bkk+= ckk*(J[Kp+1]*sqrt(Kp/(Kp+1.)));       else { std::cerr << errormsg; return; }
         Bkk *= ((Kp+1.)/(2*Kp+1.)) * ( (2*Kp+1)/(Kp+1.) ); K = Kp-1; ckk = Bkk+Akk;
      }
      else                                                      // Calculates the matrices for the second term with just B
      {
         ic = lovesey_cKK(ckk,Kp,Kp,n,l); if(ic) ckk *= J[Kp]; else { std::cerr << errormsg; return; } K = Kp;
      }

      for(Q=-K; Q<=K; Q++)
      {
         Qp = -(Q-q); if(Qp<-Kp || Qp>Kp) continue;             // 3j symbol requires Q+Q'-q = 0
 
         Ykq = spherical_harmonics(K,Q,theta,phi); 
         Tj = pow(-1.,K-Kp+q) * sqrt(3) * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q);
         if((fabs(Ykq.r)<DBL_EPSILON && fabs(Ykq.i)<DBL_EPSILON) || fabs(Tj)<DBL_EPSILON) continue;

         Qmat.zero(ns,ns); QLmat.zero(ns,ns); QSmat.zero(ns,ns);

         nz = ckk.find();
         for(i=0; i<(int)nz.size(); i++)                        // The matrices above already contain the dependence on |vULSJ>
         {                                                      //    we now use the W-E theorem to add the Jz, Q, dependence
            J2 = conf.states[nz[i][0]].J2;
            J2p = conf.states[nz[i][1]].J2;
            valJ_i = (J2-minJ2)/2; valJ_j = (J2p-minJ2)/2;
            if(mJmat[valJ_i][valJ_j].isempty())
            {
               mJmat[valJ_i][valJ_j].zero(J2+1,J2p+1);
               for(imJ_i=0; imJ_i<=J2; imJ_i++)
               {
                  Jz2 = imJ_i*2-J2;
                  for(imJ_j=0; imJ_j<=J2p; imJ_j++)
                  {
                     Jz2p = imJ_j*2-J2p; Tj2 = threej(2*Kp,J2p,J2,2*Qp,Jz2p,-Jz2);
                     if(fabs(Tj2)>DBL_EPSILON) mJmat[valJ_i][valJ_j](imJ_i,imJ_j) = pow(-1.,Kp+(-J2p+Jz2)/2.) * sqrt(J2+1.) * Tj2;
                  }
               }
            }
            Qmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*ckk(nz[i][0],nz[i][1]));
            if(Kp%2==1)
            {
               QLmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*Akk(nz[i][0],nz[i][1]));
               QSmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*Bkk(nz[i][0],nz[i][1]));
            }
         }
         for(i=0; i<=valJ; i++) for(j=0; j<=valJ; j++) mJmat[i][j].clear();

         Qq[0] += Qmat*(Ykq.r*Tj); Qq[1] += Qmat*(Ykq.i*Tj);
         if(Kp%2==1) { Qq[2] += QSmat*(Ykq.r*Tj); Qq[3] += QSmat*(Ykq.i*Tj);  Qq[4] += QLmat*(Ykq.r*Tj); Qq[5] += QLmat*(Ykq.i*Tj); }
         else        { Qq[2] += Qmat*(Ykq.r*Tj); Qq[3] += Qmat*(Ykq.i*Tj); }
      }
   }

   for(i=0; i<6; i++) { rmzeros(Qq[i]); Qq[i] *= sqrt(4*PI); } */
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the thermal expectation value of the FT of the magnetisation density -2Q in Bohr magnetons
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
          void mq(ComplexVector &Mq,      // Output expectation values -2[<Q>_{x} <Q>_y <Q>_z]
                  double &th, double &ph, // Input polar and azimuth angles theta and phi
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est)     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
{
   int i,q,n=1,Hsz=est.Cols()-1; orbital l;
   n = (int)est[0][0].real(); i = (int)est[0][0].imag(); 
   switch(i) { case 1: l=P; break; case 2: l=D; break; case 3: l=F; break; default: std::cerr << "Error - only P,D,F-electrons supported.\n"; exit(0); }
   std::vector<double> E,Jvec(6,0.); Jvec[0]=th; Jvec[1]=ph; Jvec[2]=J0; Jvec[3]=J2; Jvec[4]=J4; Jvec[5]=J6;
   std::vector< sMat<double> > Qp, Qm; 
   std::vector< std::vector< sMat<double> > > Qmat; for(i=0; i<3; i++) Qmat.push_back(Qp);
   complexdouble *zQmat, *zt, zme, zalpha, zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0;
   double zMqr,zMqi,Z=0.;
   char trans = 'U'; int incx=1;

   Mq = ComplexVector(1,3);

   icf_loveseyQq(Qm,-1,n,l,Jvec); icf_loveseyQq(Qp,1,n,l,Jvec);
   for(i=0; i<6; i++)  
   {
      Qmat[0].push_back( (Qp[i]-Qm[i]) * (-1/sqrt(2.)) );                    // Qx = -1/sqrt(2) * (Q_{+1} - Q_{-1})
      if(i%2==0) Qmat[1].push_back( (Qp[i+1]+Qm[i+1]) * (-1/sqrt(2.)) );     // real(Qy) = i^2/sqrt(2) * imag(Q_{+1}+Q_{-1})
      else       Qmat[1].push_back( (Qp[i-1]+Qm[i-1]) *  (1/sqrt(2.)) );     // imag(Qy) = i/sqrt(2) * real(Q_{+1}+Q_{-1})
   }
   icf_loveseyQq(Qmat[2],0,n,l,Jvec);

/* if(!get_Qq(Qmat[0],0,n,l,Jvec) || !get_Qq(Qmat[1],1,n,l,Jvec))            // Qmat[0]==Qx, Qmat[1]==Qy, Qmat[2]==Qz
   {
      lovesey_Qq(Qm,-1,n,l,Jvec); lovesey_Qq(Qp,1,n,l,Jvec);
      for(i=0; i<6; i++)  
      {
         Qmat[0].push_back( (Qp[i]-Qm[i]) * (-1/sqrt(2.)) );                 // Qx = -1/sqrt(2) * (Q_{+1} - Q_{-1})
         if(i%2==0) Qmat[1].push_back( (Qp[i+1]+Qm[i+1]) * (-1/sqrt(2.)) );  // real(Qy) = i^2/sqrt(2) * imag(Q_{+1}+Q_{-1})
         else       Qmat[1].push_back( (Qp[i-1]+Qm[i-1]) *  (1/sqrt(2.)) );  // imag(Qy) = i/sqrt(2) * real(Q_{+1}+Q_{-1})
      }
      save_Qq(Qmat[0],0,n,l,Jvec); save_Qq(Qmat[1],1,n,l,Jvec); 
   }

   if(!get_Qq(Qmat[2],2,n,l,Jvec))                                           // Loads the Q_q matrix if prev. saved
   {
      lovesey_Qq(Qmat[2],0,n,l,Jvec);                                        // Calcs. scattering operator matrix Q_q
      save_Qq(Qmat[2],2,n,l,Jvec);
  // myPrintMatrix(stdout,Qmat[2][0],Hsz-1);
   }
*/ for(q=0; q<3; q++)
   {
      zQmat = zmat2f(Qmat[q][0],Qmat[q][1]);
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      zMqr = 0.; zMqi = 0.;
      for(i=1; i<=Hsz; i++)
      {
         F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[i][1], &incx, &zbeta, zt, &incx);
         #ifdef _G77
         F77NAME(zdotc)(&zme, &Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx);
         #else
         zme = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx);
         #endif
//         printf ("%i zme=%g %+g i  Ei=%6.3f ni=%6.3f \n",i,zme.r,zme.i,est[0][i].real(),est[0][i].imag());
         zMqr += (-2.)*zme.r*est[0][i].imag(); zMqi += (-2.)*zme.i*est[0][i].imag(); if(q==0) Z += est[0][i].imag();
      }
      free(zQmat); free(zt); Mq[q+1] = complex<double> (zMqr, zMqi)/Z;
   }
//   printf("MQ=(%g %+g i, %g %+g i,%g %+g i)\n",real(Mq(1)),imag(Mq(1)),real(Mq(2)),imag(Mq(2)),real(Mq(3)),imag(Mq(3)));
}

/*
// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the transition matrix using the scattering operator of Balcar and Lovesey
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           int                            // Returns total number of transitions
               dncalc(int &tn,            // Input transition number |tn|. If tn<0 omit printout. If tn>0 print info.
                  double &th,             // Input zenith angle (with the z==b axis) in radians.
                  double &ph,             // Input azimuth angle (with the x==a axis, to projection in x-y plane).
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est,     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
                  double &T,              // Input temperature (K)
                  ComplexMatrix &mat)     // Output transition matrix, N(alpha,beta) = <-|Qa|+><+|Qb|-> (n- - n+)
/ * 
     Note on Qalpha (Qa or Qb)
       if gJ>0:
        Kartesian components of the scattering operator Qalpha, alpha=1,2,3=a,b,c
        according to Lovesey Neutron Scattering equation 6.87b 
        scattering operator is given in  spherical coordinates Q-1,Q0,Q+1 (introduced
        as described above on input of th and ph) these are related to Qa,Qb,Qc by
        Q1=Qbx(axb)=Qy= i/sqrt(2)(Q+1 + Q-1) 
        Q2=Qb      =Qz= Q0                   
        Q3=Qaxb    =Qx=-1/sqrt(2)(Q+1 - Q-1)
                   
       if gJ=0 
        the orbital and spin contributions have to be given as separate components of Qalpha=1,2,3,4,5,6
        according to Lovesey Neutron Scattering equations 11.55 and 11.71 (the spin part 11.71 has to be
        divided by 2), i.e.
        <-|Q1,3,5|+>=<-|QSa,b,c|+>=
          =<-|sum_i exp(i k ri) s_(a,b,c)|+> /2                   as defined by 11.71 / 2
				   
        <-|Q2,4,6|+>=<-|QLa,b,c|+>=
          =<-|sum_i exp(i k ri) (-(k x grad_i)_(a,b,c)/|k|)|+>     as defined by 11.54 /(-|k|)
* /
{
   // check if printout should be done and make tn positive
   int pr=1; if (tn<0) { pr=0; tn*=-1; }

   int i,iJ,q,n=1,Hsz=est.Cols()-1; orbital l;//=D; find_nl_from_dim(Hsz,*&n,*&l,(complexdouble*)&est[1][1]);
   n = (int)est[0][0].real(); i = (int)est[0][0].imag(); l = (i==2) ? D : F;
   std::vector<double> E,Jvec(6,0.); Jvec[0]=th; Jvec[1]=ph; Jvec[2]=J0; Jvec[3]=J2; Jvec[4]=J4; Jvec[5]=J6;
   std::vector< sMat<double> > Qp1, Qm1; 
   std::vector< std::vector< sMat<double> > > Qq; for(i=0; i<3; i++) Qq.push_back(Qp1);
   complexdouble *zQmat, *zt, zalpha, zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0;
   std::vector<complexdouble> zij(7,zbeta), zji(7,zbeta);
   double Z=0., therm;
   char trans = 'U'; int incx=1;

   // Calculates the scattering operator, Q.
   if(!get_Qq(Qq[0],0,n,l,Jvec) || !get_Qq(Qq[1],1,n,l,Jvec))                // Qq[0]==Qx, Qq[1]==Qy, Qq[2]==Qz
   {
      lovesey_Qq(Qm1,-1,n,l,Jvec); lovesey_Qq(Qp1,1,n,l,Jvec);
      for(i=0; i<6; i++)  
      {
         Qq[0].push_back( (Qp1[i]-Qm1[i]) * (-1.0/sqrt(2.0)) );                  // Qx = -1/sqrt(2) * (Q_{+1} - Q_{-1})
         if(i%2==0) Qq[1].push_back( (Qp1[i+1]+Qm1[i+1]) * (-1.0/sqrt(2.0)) );   // real(Qy) = i^2/sqrt(2) * imag(Q_{+1}+Q_{-1})
         else       Qq[1].push_back( (Qp1[i-1]+Qm1[i-1]) *  (1.0/sqrt(2.0)) );   // imag(Qy) = i/sqrt(2) * real(Q_{+1}+Q_{-1})
      }
      save_Qq(Qq[0],0,n,l,Jvec); save_Qq(Qq[1],1,n,l,Jvec); 
     // myPrintMatrix(stdout,Qq[0][2],Hsz-1);
     // myPrintMatrix(stdout,Qq[0][3],Hsz-1);
     // myPrintMatrix(stdout,Qq[0][4],Hsz-1);
    //  myPrintMatrix(stdout,Qq[0][5],Hsz-1);
   }
   if(!get_Qq(Qq[2],2,n,l,Jvec))                                             // Loads the Q_q matrix if prev. saved
   {
      lovesey_Qq(Qq[2],0,n,l,Jvec);                                          // Calcs. scattering operator matrix Q_q
      save_Qq(Qq[2],2,n,l,Jvec);
   }

   for(i=0; i<Hsz; i++) { 
      therm = exp(-(est[0][i+1].real()-est[0][1].real())/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }

   int a,b,j=0,k=0; for(i=0; i<Hsz; ++i) {for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   ++i;++j; // because in est i and j start from 1...Hsz
 
   for(q=0; q<3; q++)
   {  zQmat = zmat2f(Qq[q][2],Qq[q][3]);     // Spin part
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[j][1], &incx, &zbeta, zt, &incx);
#ifdef _G77
      F77NAME(zdotc)(&zij[2*q+1], &Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx) ;
#else
      zij[2*q+1] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx) ;
#endif
//         int k;for(k=0;k<Hsz;++k)printf("%6.3f %+6.3f i  ",real(est(j,1+k)),imag(est(j,1+k)));
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[i][1], &incx, &zbeta, zt, &incx);
#ifdef _G77
      F77NAME(zdotc)(&zji[2*q+1], &Hsz, (complexdouble*)&est[j][1], &incx, zt, &incx) ;
#else
      zji[2*q+1] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[j][1], &incx, zt, &incx) ;
#endif

      if(i==j)                               //subtract thermal expectation value from zij=zii
      {
         complexdouble expQ; double thexp=0;
         for(iJ=1;iJ<=Hsz;++iJ)    
         {
            therm = exp(-(est[0][iJ].real()-est[0][1].real())/(KB*T)); if(therm<DBL_EPSILON) break; 
            F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[iJ][1], &incx, &zbeta, zt, &incx);
#ifdef _G77
            F77NAME(zdotc)(&expQ, &Hsz, (complexdouble*)&est[iJ][1], &incx, zt, &incx);
#else
            expQ = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[iJ][1], &incx, zt, &incx);
#endif
            thexp += expQ.r * therm / Z;
         }
         zij[2*q+1].r-=thexp;zji[2*q+1].r-=thexp;
      }
      free(zQmat); free(zt);

      zQmat = zmat2f(Qq[q][4],Qq[q][5]);     // orbital part
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[j][1], &incx, &zbeta, zt, &incx);
#ifdef _G77
      F77NAME(zdotc)(&zij[2*q+2], &Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx);
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[i][1], &incx, &zbeta, zt, &incx);
      F77NAME(zdotc)(&zji[2*q+2], &Hsz, (complexdouble*)&est[j][1], &incx, zt, &incx);
#else
      zij[2*q+2] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx);
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[i][1], &incx, &zbeta, zt, &incx);
      zji[2*q+2] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[j][1], &incx, zt, &incx);
#endif
      if(i==j)                               //subtract thermal expectation value from zij=zii
      {
         complexdouble expQ;double thexp=0;
         for(iJ=1;iJ<=Hsz;++iJ)    
         {
            therm = exp(-(est[0][iJ].real()-est[0][1].real())/(KB*T)); if(therm<DBL_EPSILON) break; 
            F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[iJ][1], &incx, &zbeta, zt, &incx);
#ifdef _G77
            F77NAME(zdotc)(&expQ, &Hsz, (complexdouble*)&est[iJ][1], &incx, zt, &incx);
#else
            expQ = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[iJ][1], &incx, zt, &incx);
#endif
            thexp += expQ.r * therm / Z;
         }
         zij[2*q+2].r-=thexp;zji[2*q+2].r-=thexp;
      }
      free(zQmat); free(zt);
   }

   // check if zij are complex conjugate
   for(iJ=1;iJ<=6;++iJ)
      if(fabs(zij[iJ].i+zji[iJ].i)>SMALL) { std::cerr << "ERROR module ic1ion - dncalc: <i|Qalpha|j>not hermitian\n"; exit(EXIT_FAILURE); }

                
   complex<double> im(0,1);
   ComplexVector iQalphaj(1,6);
   
   for(a=1; a<=6; a++){iQalphaj(a) = complex<double> (zij[a].r,zij[a].i);if(a%2==1){iQalphaj(a)*=0.5;}} 
                                                                         // divide spin part by 2
   mat = 0;
   for(a=1; a<=6; a++)
      for(b=1; b<=6; b++)
         mat(a,b) = iQalphaj(a)*conj(iQalphaj(b));

   double delta;
   delta = est[0][j].real()-est[0][i].real();
   if(delta<-0.000001) { std::cerr << "ERROR module ic1ion - dncalc: energy gain delta gets negative\n"; exit(EXIT_FAILURE); }

   if(j==i) delta = -SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc

   // do some printout if wishes and set correct occupation factor
   if (delta>SMALL)
   {
      therm = exp(-(est[0][i].real()-est[0][1].real())/(KB*T)) - exp(-(est[0][j].real()-est[0][1].real())/(KB*T));
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV",i,j,delta);
         printf(" |<%i|Qa|%i>|^2=%6.3f |<%i|Qb|%i>|^2=%6.3f |<%i|Qc|%i>|^2=%6.3f",i,j,real(mat(1,1)),i,j,real(mat(2,2)),i,j,real(mat(3,3)));
         printf(" |<%i|Qd|%i>|^2=%6.3f |<%i|Qe|%i>|^2=%6.3f |<%i|Qf|%i>|^2=%6.3f",i,j,real(mat(4,4)),i,j,real(mat(5,5)),i,j,real(mat(6,6)));
         printf(" n%i-n%i=%6.3f\n",i,j,therm / Z);
      }
   }
   else
   {
      therm = exp(-(est[0][i].real()-est[0][1].real())/(KB*T))/(KB*T);
          // quasielastic scattering has not wi-wj but wj*epsilon/kT
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV",i,j,delta);
         printf(" |<%i|Qa-<Qa>|%i>|^2=%6.3f |<%i|Qb-<Qb>|%i>|^2=%6.3f |<%i|Qc-<Qc>|%i>|^2=%6.3f",i,j,real(mat(1,1)),i,j,real(mat(2,2)),i,j,real(mat(3,3)));
         printf(" |<%i|Qd-<Qd>|%i>|^2=%6.3f |<%i|Qe-<Qe>|%i>|^2=%6.3f |<%i|Qf-<Qf>|%i>|^2=%6.3f",i,j,real(mat(4,4)),i,j,real(mat(5,5)),i,j,real(mat(6,6)));
         printf(" n%i=%6.3f\n",i,therm/Z);
      }
   }
         mat *= therm / Z;

   // determine number of thermally reachable states
   int noft=0;for(i=0;(i<Hsz)&((exp(-(est[0][i+1].real()-est[0][1].real())/(KB*T)))>SMALL);++i)noft+=Hsz-i-1;
   return noft;
//   return Hsz*(Hsz-1)/2;
}


// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of spindensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
void sdod_mcalc(Vector &J,           // Output single ion moments==(expectation values) Zlm R^2(r) at given T, H_eff
                int xyz,             // direction 1,2,3 = x,y,z
                double *T,           // Input scalar temperature
                Vector &gjmbH,       // Input vector of mean fields (meV)
                char **sipffilename, // Single ion properties filename
                ComplexMatrix &est)  // Input/output eigenstate matrix (initialized in parstorage)
{  
   // Parses the input file for parameters
   icpars pars;
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

   // Prints out single ion parameters and filename.
// std::cout << "#ic1ion: Single ion parameter file: " << filename << "\n";
// std::cout << "#ic1ion: Read in parameters (in " << pars.e_units << "): F2=" << pars.F[1] << ",F4=" << pars.F[2];
// if(pars.l==F) std::cout << ",F6=" << pars.F[3]; std::cout << ",xi=" << pars.xi << "\n";
// std::cout << "#ic1ion: \t" << pars.B.cfparsout(", ") << "in " << pars.B.units() << "\n";

   // Calculates the IC Hamiltonian matrix
   int i,k,q,Hsz=getdim(pars.n,pars.l);
   complexdouble *H=0,*Jm=0; 
   bool Hicnotcalc = false;
   std::vector<double> parval; parval.reserve(35);
   if(pars.n==1) { parval.push_back(pars.xi); for(k=2; k<=(2*pars.l); k+=2) for(q=-k; q<=k; q++) parval.push_back(pars.B(k,q)); }
   else {
      for(i=0; i<4; i++) parval.push_back(pars.F[i]); parval.push_back(pars.xi); for(i=0; i<3; i++) parval.push_back(pars.alpha[i]);
      for(k=2; k<=(2*pars.l); k+=2) for(q=-k; q<=k; q++) parval.push_back(pars.B(k,q)); }
   if(parval.size()%2==1) parval.push_back(0.);

   if((est.Cols()!=(Hsz+1) || est.Rows()!=(Hsz+1)) && pars.truncate_level==1) Hicnotcalc = true;
   else if(real(est[0][0])==-0.1 && imag(est[0][0])==-0.1)  // Hic previously calculated
   {
      for(i=0; i<(int)(parval.size()/2); i++) if(real(est[0][i+1])!=parval[2*i] && imag(est[0][i+1])!=parval[2*i+1]) { Hicnotcalc = true; break; }
   }
   else Hicnotcalc = true;
   if(Hicnotcalc)
   {
      if(pars.spectrelevels==-1)
      {
         sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars); Hic/=MEV2CM; iHic/=MEV2CM; H = zmat2f(Hic,iHic);
//       est = ComplexMatrix(0,Hsz,0,Hsz); I comment this out - you should not reinitialize est !!!
         if( (est.Rhi()!=Hsz||est.Chi()!=Hsz) && pars.truncate_level==1)
         {
            std::cerr << "ERROR module ic1ion - mcalc: Hsz recalculation does not agree with eigenstates matrix dimension\n"; exit(EXIT_FAILURE);
         }
         if ((int)(parval.size()/2)>Hsz && pars.truncate_level==1)
         {
//          std::cerr << "ERROR module ic1ion - mcalc: storing Hamiltonian in est matrix not possible - parval size too big\n"; exit(EXIT_FAILURE);
         }
         else
         {
            est[0][0] = complex<double> (-0.1,-0.1);
            for(i=0; i<(int)(parval.size()/2); i++) est[0][i+1] = complex<double> (parval[2*i],parval[2*i+1]);
         }
         if(pars.truncate_level!=1)  // Truncates the matrix, and stores the single ion part in a memorymapped array (Linux only)
            truncate_hmltn(pars, est, Hic, iHic, J.Hi(), J.Lo());
         else
            for(i=1; i<=Hsz; i++) memcpy(&est[i][1],&H[(i-1)*Hsz],Hsz*sizeof(complexdouble));
         free(H);
      }
      else  // Calculates using the Spectre method...                                       // lvl = 10;     % Number of |LSJ> levels to keep
      {
         if(pars.spectrelevels==-2)
         {
            std::cerr << "Trying to determine optimal number of levels for spectre... ";
            std::fstream FH; FH.open("results/ic1ion.spectre", std::fstream::out);
            sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars); Hic/=MEV2CM; iHic/=MEV2CM;
//            std::vector<double> vgjmbH((J.Hi()-J.Lo()+1),0.); for(i=J.Lo(); i<=J.Hi(); i++) vgjmbH[i-J.Lo()] = -gjmbH[i];
//            icmfmat mfmat(pars.n,pars.l,J.Hi()-J.Lo()+1,pars.save_matrices); sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH,pars.save_matrices);
            std::vector<double> vgjmbH((51),0.); for(i=gjmbH.Lo(); i<=gjmbH.Hi()&&i<=51; i++) vgjmbH[i-gjmbH.Lo()] = -gjmbH[i];
            icmfmat mfmat(pars.n,pars.l,51,pars.save_matrices); sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH,pars.save_matrices);
            int cbbest=1,ibest=1; Jmat+=Hic; iJmat+=iHic; Jm = zmat2f(Jmat,iJmat); iceig VE; VE.calc(Hsz,Jm);
            double Uref,dbest=DBL_MAX,dnew=DBL_MAX; std::vector<double> vBest,vNew; // double Unew,Ubest;
            std::vector< std::vector<double> > matel; std::vector<double> vJ = mfmat.expJ(VE,*T,matel,pars.save_matrices); Uref=vJ[J.Hi()-J.Lo()+2];
            FH << "\nRef:\t["; for(i=0; i<(int)vJ.size(); i++) FH << vJ[i] << "\t"; FH << "]\n";
            for(int cb=1; cb<=Hic.nr(); cb++)
            {
               pars.spectrelevels=cb; spectre_hmltn(pars,est,parval.size());
               std::vector<double> vJt = spectre_expJ(pars,est,parval.size(),gjmbH,J.Lo(),J.Hi(),*T);
               FH << "cb=" << cb << "\t["; for(i=J.Lo(); i<(int)vJt.size(); i++) FH << vJt[i] << "\t"; FH << "]\n";
               if(cb==1)
               {
                  dbest=0.; for(i=J.Lo(); i<=J.Hi(); i++) dbest+=fabs(vJ[i-J.Lo()]-vJt[i]); vBest=vJt;
               }
               else
               {
                  vNew = vJt; double dt=0.; for(i=J.Lo(); i<=J.Hi(); i++) dt+=fabs(vJ[i-J.Lo()]-vJt[i]);
                  dnew = dt; if(dnew<dbest) { vBest=vNew; dbest=dnew; cbbest=cb; }
                  if(dnew>dbest) ibest++; if(ibest==5) break;
               }
            }
            FH.close();
            pars.spectrelevels=cbbest; std::cerr << cbbest << " levels seems optimal\n";
         }
         spectre_hmltn(pars,est,parval.size());
         for(i=0; i<(int)(parval.size()/2); i++) est[0][i+1] = complex<double> (parval[2*i],parval[2*i+1]);
      }
   }
   if(pars.spectrelevels==-1)
   {
      if(pars.truncate_level!=1)  // Uses the eigenvectors of the single ion Hamiltonian to truncate the matrix
         { std::cerr << "spin/orbmomdensity using truncate otion not yet implemented... exiting ";exit(EXIT_FAILURE);}
//          truncate_expJ(pars,est,gjmbH,J,*T,lnZ,U,Jm);
      else
      {
         // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
//         icmfmat mfmat(pars.n,pars.l,J.Hi()-J.Lo()+1,pars.save_matrices,pars.density);
//         std::vector<double> vgjmbH((J.Hi()-J.Lo()+1),0.); for(i=J.Lo(); i<=J.Hi(); i++) vgjmbH[i-J.Lo()] = -gjmbH[i];
         icmfmat mfmat(pars.n,pars.l,51,pars.save_matrices,pars.density);
         std::vector<double> vgjmbH(51,0.); for(i=gjmbH.Lo(); i<=gjmbH.Hi()&&i<=51; i++) vgjmbH[i-gjmbH.Lo()] = -gjmbH[i];
         sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH,pars.save_matrices);
         complex<double> a(1.,0.); int incx = 1;
         Jm = zmat2f(Jmat,iJmat); for(i=1; i<=Hsz; i++) F77NAME(zaxpy)(&Hsz,(complexdouble*)&a,(complexdouble*)&est[i][1],&incx,&Jm[(i-1)*Hsz],&incx);

         // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
         iceig VE; if(pars.partial) VE.lcalc(pars,Jm); else if(pars.arnoldi) VE.acalc(pars,Jm); else VE.calc(Hsz,Jm); free(Jm);

         // Calculates the expectation values sum_n{ <n|Ja|n> exp(-En/kT) }
         std::vector< std::vector<double> > matel;
         std::vector<double> vJ = mfmat.spindensity_expJ(VE, xyz,*T,matel,pars.save_matrices);
         for(i=J.Lo(); i<=J.Hi(); i++) J[i] = vJ[i-J.Lo()];//printf("ss=%g\n",J(1));
      }
   }
   else
   {   std::cerr << "spin/orbmomdensity using spectre otion not yet implemented... exiting ";exit(EXIT_FAILURE);
      //std::vector<double> vJ = spectre_expJ(pars,est,parval.size(),gjmbH,J.Lo(),J.Hi(),*T);
      //for(i=J.Lo(); i<=J.Hi(); i++) J[i] = vJ[i];
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of spindensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
void spindensity_mcalc(Vector &J,          // Output single ion moments =expectation values of
                                           // of Zlm R^2(r) at a given temperature T and  effective field H
                      int & xyz,           // direction 1,2,3 = x,y,z
                      double *T,           // Input scalar temperature
                      Vector &gjmbH,       // Input vector of mean fields (meV)
 / * Not Used * /       double *gJ,          // Input Lande g-factor
 / * Not Used * /       Vector &ABC,         // Input vector of parameters from single ion property file
                      char **sipffilename, // Single ion properties filename
                      ComplexMatrix &est)  // Input/output eigenstate matrix (initialized in parstorage)
{
   sdod_mcalc(J,xyz,T,gjmbH,sipffilename,est);
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
void orbmomdensity_mcalc(Vector &J,        // Output single ion moments =expectation values of
                                           // of Zlm R^2(r) at a given temperature T and  effective field H
                      int & xyz,           // direction 1,2,3 = x,y,z
                      double *T,           // Input scalar temperature
                      Vector &gjmbH,       // Input vector of mean fields (meV)
 / * Not Used * /       double *gJ,          // Input Lande g-factor
 / * Not Used * /       Vector &ABC,         // Input vector of parameters from single ion property file
                      char **sipffilename, // Single ion properties filename
                      ComplexMatrix &est)  // Input/output eigenstate matrix (initialized in parstorage)
{
   sdod_mcalc(J,-xyz,T,gjmbH,sipffilename,est);
}
*/
// --------------------------------------------------------------------------------------------------------------- //
// Prints out a header to a specified file
// --------------------------------------------------------------------------------------------------------------- //
void icf_printheader(const char *outfile, icpars &pars)
{
   time_t curtime = time(NULL);
   std::fstream FILEOUT; FILEOUT.open(outfile, std::fstream::out);
   std::string Lstr = Lout(pars.l); strtolower(Lstr);
   FILEOUT << "#{ icf1ionmodule version " << IC1IONMODULE_VERSION << " " << ctime(&curtime);
   if(!pars.ionname.empty()) FILEOUT << "# Ion name: " << pars.ionname << "\n";
   FILEOUT << "# Free ion configuration: " << Lstr << "^" << pars.n << "\n";
   FILEOUT << "# Spin orbit parameter (" << pars.e_units << "): xi=" << pars.xi << "\n";
   FILEOUT << "# Crystal Field parameters normalisation: " << pars.B.norm() << "\n";
   std::string norm=pars.B.norm(); strtolower(norm); if(norm.find("stev")!=std::string::npos)
   {
      FILEOUT << "# Stevens Factors: alpha=" << pars.B.alpha() << ", beta=" << pars.B.beta();
      if(pars.l==F) FILEOUT << ", gamma=" << pars.B.gamma(); FILEOUT << "\n";
   }
   FILEOUT << "# Crystal Field parameters (" << pars.B.units() << "): " << pars.B.cfparsout(", ") << "\n";
   if(fabs(pars.Bx)>DBL_EPSILON || fabs(pars.By)>DBL_EPSILON || fabs(pars.Bz)>DBL_EPSILON)
   {
      FILEOUT << "# With magnetic field: Bx=" << pars.Bx << ", By=" << pars.By << ", Bz=" << pars.Bz << " Tesla.\n";
   }
   FILEOUT.close();
}
 
// --------------------------------------------------------------------------------------------------------------- //
// Outputs the energy and wavefunctions to a specified file
// --------------------------------------------------------------------------------------------------------------- //
void icf_showoutput(const char *filename,                       // Output file name - default "results/mcphas.icr"
                    icpars &pars,                               // Input parameters
                    ComplexMatrix &est)                         // Eigenstates class
{
   fstates_t gs = hunds_gs(pars.n, pars.l);
   char Jlabel[12]; std::string id; 
   int J2min, J2max, ns=0;

   // Determines number of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2;
   for (int J2=J2min; J2<=J2max; J2+=2) ns+=J2+1; 
   std::vector<std::string> statesID(ns);
   // Determines angular momentum quantum numbers of basis states
   int ins=0; std::vector<int> J2(ns), mJ2(ns);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2)
     for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2)
     {
        J2[ins]=iJ2; mJ2[ins]=imJ2; 
        id.assign(gs.id); id.append("_");
        if(iJ2%2==0) sprintf(Jlabel,"%i",iJ2/2);      else sprintf(Jlabel,"%i/2",iJ2);      id.append(Jlabel);
        if(imJ2%2==0)sprintf(Jlabel,",mJ=%i",imJ2/2); else sprintf(Jlabel,",mJ=%i/2",imJ2); id.append(Jlabel);
        statesID[ins].assign(id); ins++; 
     }

   int iE,iV,i=1,j=2;
   std::vector<int> isE,isV(ns,0); isE.reserve(ns);
   int ii; 
   double elem,conv=1.; complexdouble elc;

   if(pars.e_units.find("cm")!=std::string::npos) conv = MEV2CM; 
   else if(pars.e_units.find("K")!=std::string::npos) conv = MEV2K;

// if(iconf==1) icf_conv_basis(pars,VE);

   bool iscomplex=false; if(fabs(pars.By)>SMALL) iscomplex=true;
   for(int K=2; K<=6; K+=2) for(int Q=-K; Q<0; Q++) if(fabs(pars.B(K,Q))>SMALL) { iscomplex=true; break; }

   icf_printheader(filename,pars);
   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out | std::fstream::app); // Opens file for appending
   FILEOUT << "# Energy offset, E0=" << real(est(0,1))*conv << pars.e_units << "\n";
   if(!iscomplex) FILEOUT << "# Energy(" << pars.e_units << ")\tWavefunctions(^{2S+1}L_J,mJ) }\n"; else
   FILEOUT << "# Energy(" << pars.e_units << ")\tAmplitude\t|Amplitude|^2\tWavefunctions(^{2S+1}L_J,mJ) }\n";

   double *V=0; complexdouble *zV=0; if(iscomplex) zV = new complexdouble[ns]; else V = new double[ns];
   for(iE=0; iE<ns; iE++)
   {
      if(real(est(0,iE+1))==0.) if(iE<(ns-1) && real(est(0,iE+2))==0.) break; 
      FILEOUT << (real(est(0,iE+1))-real(est(0,1)))*conv << "\t\t";
      for(ii=0; ii<(int)ns; ii++) isV[ii]=ii;
      i=1; j=2;
      if(iscomplex)
      {
         memcpy(zV,&est[iE+1][1],ns*sizeof(complexdouble));
         while(i<ns)
         {
            if((zV[i-1].r*zV[i-1].r+zV[i-1].i*zV[i-1].i) >= (zV[i].r*zV[i].r+zV[i].i*zV[i].i)) { i=j; j++; }
            else { elc = zV[i-1]; zV[i-1] = zV[i]; zV[i] = elc; ii=isV[i-1]; isV[i-1]=isV[i]; isV[i]=ii; i--; if(i==0) i=1; }
         }
         elc = zV[0];
         FILEOUT << " (" << elc.r; if(elc.i>0) FILEOUT << "+"; else FILEOUT << "-";
         FILEOUT << "i" << elc.i << ")\t\t" << (elc.r*elc.r+elc.i*elc.i) << "\t";
         FILEOUT << "|" << statesID[isV[0]] << ">";
         for(iV=1; iV<pars.num_eigv; iV++)
         {
            elc = zV[iV]; FILEOUT << "\n\t\t+";
            FILEOUT << "(" << elc.r; if(elc.i>0) FILEOUT << "+"; else FILEOUT << "-";
	    FILEOUT << "i" << elc.i << ")\t\t" << (elc.r*elc.r+elc.i*elc.i) << "\t";
            FILEOUT << "|" << statesID[isV[iV]] << ">";
         }
         FILEOUT << "\n";
      }
      else
      {
         for(int ei=0; ei<ns; ei++) V[ei]=real(est(iE+1,ei+1));
         while(i<ns)
         {
            if(fabs(V[i-1])>=fabs(V[i])) { i=j; j++; }
            else { elem = V[i-1]; V[i-1] = V[i]; V[i] = elem; ii=isV[i-1]; isV[i-1]=isV[i]; isV[i]=ii; i--; if(i==0) i=1; }
         }
         FILEOUT << V[0] << "|" << statesID[isV[0]] << ">\t";
         for(iV=1; iV<pars.num_eigv; iV++)
         {
            if(iV%4==0) FILEOUT << "\n\t\t";
            if(V[iV]>0) FILEOUT << "+";
            FILEOUT << V[iV] << "|" << statesID[isV[iV]] << ">\t";
         }
         FILEOUT << "\n";
      }
   }
   if(iscomplex) delete[]zV; else delete[]V;
   FILEOUT.close();
}


// --------------------------------------------------------------------------------------------------------------- //
// Reads in parameters from mcphas.ic or a file specified on the command line.
// Calculates the IC Hamilton matrix; solve its eigensystem and outputs the energies and wavefunctions calculated.
// Also calculates the magnetisation and heat capacity if desired (specify iccalcmag=1,2 and/or iccalcCp=1)
// --------------------------------------------------------------------------------------------------------------- //
int main(int argc, char *argv[])
{
   char *infile, outfile[255], physfile[255]; infile=(char*)malloc(255*sizeof(char));
   icpars pars;
   std::string norm,units;

   if(argc>1) strcpy(infile,  argv[1]); else strcpy(infile,"mcphas.icf");
   if(argc>2) strcpy(outfile, argv[2]); else strcpy(outfile,"results/icf1ion.out");
   if(argc>3) strcpy(physfile,argv[3]); else strcpy(physfile,"results/icf1ion.mag");

   // Gets input parameters and what to calculate from input single-ion parameters file
   const char *filename = &infile[0];
   ic_parseinput(filename,pars);

   // Calculates and diagonalises the Hamiltonian
   Vector gjmbH(1,6,0.); ComplexMatrix est; double T=1.;
   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(pars.Bx)>DBL_EPSILON || fabs(pars.By)>DBL_EPSILON || fabs(pars.Bz)>DBL_EPSILON)
   {
      if(fabs(pars.Bx)>DBL_EPSILON) { gjmbH(2)=MUB*pars.Bx; gjmbH(1)=GS*gjmbH(2); }
      if(fabs(pars.By)>DBL_EPSILON) { gjmbH(4)=MUB*pars.By; gjmbH(3)=GS*gjmbH(4); }
      if(fabs(pars.Bz)>DBL_EPSILON) { gjmbH(6)=MUB*pars.Bz; gjmbH(5)=GS*gjmbH(6); }
   }
   estates(&est,gjmbH,T,T,gjmbH,&infile);
   icf_showoutput(outfile,pars,est);

   // If required calculates some single-ion physical properties and saves to <physfile>
   if(pars.calcphys & PHYSPROP_MAGBIT)
   {
      double lnZ=0., U=0., Tmin, Tstep, Tmax, Hmin, Hstep, Hmax, Hm;
      if(pars.xT==0. && pars.yT==0.) { std::cerr << "icf1ion: Either x- or y-axis must be temperature\n"; exit(-1); }
      if(pars.xT!=0.) { Tmin = pars.xMin; Tstep = pars.xStep; Tmax = pars.xMax; Hmin = pars.yMin; Hstep = pars.yStep; Hmax = pars.yMax; }
                 else { Tmin = pars.yMin; Tstep = pars.yStep; Tmax = pars.yMax; Hmin = pars.xMin; Hstep = pars.xStep; Hmax = pars.xMax; }

      icf_printheader(physfile,pars);
      // Opens file for appending and writes header
      std::fstream FILEOUT; FILEOUT.open(physfile, std::fstream::out | std::fstream::app);
      if(pars.xT!=0.) FILEOUT << "# T(K)\tHa(T)\tHb(T)\tHc(T)\t"; else FILEOUT << "# Ha(T)\tHb(T)\tHc(T)\tT(K)\t";
      if(pars.mag_units==0) FILEOUT << "Magnetisation(uB/atom)\tMa\tMb\tMc\tM_parallel\n";
      else if(pars.mag_units==1) FILEOUT << "Magnetisation(emu/mol)\tMa\tMb\tMc\tM_parallel\n";
      else if(pars.mag_units==2) FILEOUT << "Magnetisation(Am^2/mol)\tMa\tMb\tMc\tM_parallel\n";
      int nT = (int)ceil((Tmax-Tmin)/Tstep), nH = (int)ceil((Hmax-Hmin)/Hstep) + 1;
      std::vector<double> vT(nT,0.); for(int iT=0; iT<nT; iT++) vT[iT] = Tmin+iT*Tstep;
      std::vector<double> mag(nT,0.),ma(nT,0.),mb(nT,0.),mc(nT,0.);

      // Determines the require magnetic fields
      double xnorm = sqrt(pars.xHa*pars.xHa+pars.xHb*pars.xHb+pars.xHc*pars.xHc); if(xnorm==0) xnorm=1.;
      double ynorm = sqrt(pars.yHa*pars.yHa+pars.yHb*pars.yHb+pars.yHc*pars.yHc); if(ynorm==0) ynorm=1.;
      Vector gjmbH0(1,6,0.),J(1,6,0.);
      if(pars.xT==0.) gjmbH0(2)=MUB*pars.xHa/xnorm; else gjmbH0(2)=MUB*pars.yHa/ynorm; gjmbH0(1)=GS*gjmbH0(2);
      if(pars.xT==0.) gjmbH0(4)=MUB*pars.xHb/xnorm; else gjmbH0(4)=MUB*pars.yHb/ynorm; gjmbH0(3)=GS*gjmbH0(4);
      if(pars.xT==0.) gjmbH0(6)=MUB*pars.xHc/xnorm; else gjmbH0(6)=MUB*pars.yHc/ynorm; gjmbH0(5)=GS*gjmbH0(6);
      double convfact=1; if(pars.mag_units==1) convfact = NAMUB*1e3; else if(pars.mag_units==2) convfact = NAMUB;  // 1==cgs, 2==SI
      
      // Reinitialises the estates matrix
      est.Remove(); mcalc_parameter_storage_matrix_init(&est,gjmbH,&T,&T,gjmbH,&infile);

      for(int iH=0; iH<nH; iH++)
      {
         for(int al=1; al<=6; al++) gjmbH(al) = gjmbH0(al)*(Hmin+iH*Hstep);
         for(int iT=0; iT<nT; iT++)
         {
            mcalc(J,&vT[iT],gjmbH,&vT[iT],gjmbH,&infile,&lnZ,&U,est); for(int iJ=1; iJ<=6; iJ++) if(fabs(J(iJ))<SMALL) J(iJ)=0.;
            ma[iT]=GS*J(1)+J(2); mb[iT]=GS*J(3)+J(4); mc[iT]=GS*J(5)+J(6); mag[iT] = sqrt(ma[iT]*ma[iT]+mb[iT]*mb[iT]+mc[iT]*mc[iT]);
         }
         if(pars.xT!=0.)
         {
            Hm = (Hmin+iH*Hstep)/ynorm;
            for(int iT=0; iT<nT; iT++) 
               FILEOUT << vT[iT] << "\t" << Hm*pars.yHa << "\t" << Hm*pars.yHb << "\t" << Hm*pars.yHc << "\t" << mag[iT]*convfact
                       << "   \t" << ma[iT]*convfact << "   \t" << mb[iT]*convfact << "   \t" << mc[iT]*convfact << "   \t"
                       << (ma[iT]*pars.yHa + mb[iT]*pars.yHb + mc[iT]*pars.yHc)*convfact << "\n";
         }
         else
         {
            Hm = (Hmin+iH*Hstep)/xnorm;
            for(int iT=0; iT<nT; iT++) 
               FILEOUT << Hm*pars.xHa << "\t" << Hm*pars.xHb << "\t" << Hm*pars.xHc << "\t" << vT[iT] << "\t" << mag[iT]*convfact
                       << "   \t" << ma[iT]*convfact << "   \t" << mb[iT]*convfact << "   \t" << mc[iT]*convfact << "   \t" 
                       << (ma[iT]*pars.xHa + mb[iT]*pars.xHb + mc[iT]*pars.xHc)*convfact << "\n";
         }
      }
      FILEOUT.close();
   }
/*   
#ifdef _INTEGRAL
 //clock_t start,end; end = clock();
   Vector J(1,6,0.), gmbH(1,6,.0578838263), ABC; gmbH[1]*=2; gmbH[3]*=2; gmbH[5]*=2; 
   double T=2.0,lnZ=0.,U=0.,gJ=0.;
   char *filearray[1]; 
   filearray[0] = infile;
   ComplexMatrix est; mcalc_parameter_storage_matrix_init(&est,gmbH,&gJ,&T,ABC,filearray);
 //ComplexMatrix est; int Hsz=getdim(pars.n,pars.l); est = ComplexMatrix(0,Hsz,0,Hsz);
   end = clock();

   mcalc(J,&T,gmbH,&gJ,ABC,filearray,&lnZ,&U,est);
   start = clock(); std::cerr << "Time to do mcalc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";
   std::cerr << "lnZ = " << lnZ << ", U = " << U << "\n";
   std::cerr << "J[1] = " << J[1] << ", J[2] = " << J[2] << ", J[3] = " << J[3] << ", J[4] = " << J[4] << ", J[5] = " << J[5] << ", J[6] = " << J[6] << "\n";

//for (int it=0; it<1000; it++) {
// end = clock();
// mcalc(J,&T,gmbH,&gJ,ABC,filearray,&lnZ,&U,est);
// start = clock(); std::cerr << "Time to do mcalc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";
// std::cerr << "lnZ = " << lnZ << ", U = " << U << "\n";
// std::cerr << "J[1] = " << J[1] << ", J[2] = " << J[2] << ", J[3] = " << J[3] << ", J[4] = " << J[4] << ", J[5] = " << J[5] << ", J[6] = " << J[6] << "\n";
// if(it%100==0) { std::cerr << it << " "; } } std::cerr << "\n";
   gmbH[1] = 0.; gmbH[2] = 0.; gmbH[3] = 0.; gmbH[4] = 0.; gmbH[5] = 0.; gmbH[6] = 0.; 
   estates(est,gmbH,&gJ,&T,ABC,filearray);
   end = clock(); std::cerr << "Time to do estates() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   
   int imq, tn = 2; float delta=0.; ComplexMatrix mat6(1,6,1,6);
   imq = dmcalc(tn,T,gmbH,gJ,ABC,filearray,mat6,delta,est);
   start = clock(); std::cerr << "Time to calculate dmcalc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";

   ComplexVector Mq;
   double th=PI/4, ph=PI/4, J0=1., J2=1., J4=1., J6=1.;
 //double th=0., ph=0., J0=1., J2=1., J4=0., J6=0.;
   imq = mq(Mq,th,ph,J0,J2,J4,J6,est);
   end = clock(); std::cerr << "Time to calculate mq() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.";
   std::cerr << " Mq = [" << Mq[1].real() << "+" << Mq[1].imag() << "i "
                          << Mq[2].real() << "+" << Mq[2].imag() << "i "
                          << Mq[3].real() << "+" << Mq[3].imag() << "i]\n";

   ComplexMatrix mat;
   imq = dncalc(tn,th,ph,J0,J2,J4,J6,est,T,mat);
   start = clock(); std::cerr << "Time to calculate dncalc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";
#endif
*/
/* int i, Hsz = est.Cols(); //complexdouble *cest = new complexdouble[Hsz*Hsz]; memcpy(cest,&est[0][0],Hsz*Hsz*sizeof(complexdouble));
   std::vector<double> vgjmbH(6,.0578838263); //for(i=0; i<6; i++) vgjmbH[i] = gjmbH[i+1];
   icmfmat mfmat(pars.n,pars.l);
   sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH); complexdouble *zJmat=0; zJmat = zmat2f(Jmat,iJmat);
   complexdouble *zVd = new complexdouble[(Hsz-1)*(Hsz-1)]; double *eigval = new double[Hsz-1];
   start = clock();
   int info = ic_peig(Hsz-1, zJmat, (complexdouble*)&est[0][0], zVd, eigval);
   end = clock(); std::cerr << "Time to do peig() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   iceig VE(Hsz-1,eigval,zVd); strcpy(outfile,"results/mcphas.ic1"); ic_showoutput(outfile,pars,VE);
   free(zJmat); delete[]zVd; delete[]eigval; //delete []cest;
   start = clock();
   Hic/=MEV2CM; iHic/=MEV2CM; 
   VE.calc(Hic,iHic); strcpy(outfile,"results/mcphas.icp"); ic_showoutput(outfile,pars,VE);
   end = clock(); std::cerr << "Time to do VE.calc() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   Hic+=Jmat; iHic+=Jmat;
   VE.calc(Hic,iHic); strcpy(outfile,"results/mcphas.ic2"); ic_showoutput(outfile,pars,VE); */

   return 0;
}
