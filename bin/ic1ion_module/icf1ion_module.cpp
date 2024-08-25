/* icf1ion_module.cpp
 *
 * Functions:
 *   void myPrintMatrix(FILE * file,sMat<double> & M,int d)                                   // Prints out full matrix
 *   bool checkmat(ComplexMatrix &cmat, complexdouble *fmat,int r, int c)                     // Compares Matpack and fortran matrices
 *   void Icalc(Vector &J, double *T,Vector &gjmbHxc,Vector&Hext, double *gJ, Vector &ABC,    // Calculates the meanfield moment
 *                 char **sipffile, double *lnZ, double *U, ComplexMatrix &est)
 *   void IMcalc(Matrix &J, Vector &T, Vector &Hxc,Vector&Hext, double *gJ, Vector &ABC,      // Calculates the meanfield moment
 *                 char **sipffile, Vector &lnZ, Vector &U, ComplexMatrix &est)  
 *   int du1calc(int &tn, double &T, Vector &gjmbHxc,Vector&Hext, double &g_J, Vector &ABC,   // Calculates the transition
 *                 char **sipffilename, ComplexVector &u1, float &delta,int & n, int & nd, ComplexMatrix &est)  //   matrix elements
 *   void mcalc(Vector &J, double *T, Vector &Hxc,Vector&Hext, double *gJ, Vector &ABC,       // Calculates the meanfield moment
 *                 char **sipffile,  ComplexMatrix &est)
 *   void mMcalc(Matrix &J, Vector &T, Vector &Hxc,Vector&Hext, double *gJ, Vector &ABC,       // Calculates the meanfield moment
 *                 char **sipffile,  ComplexMatrix &est)
 *   int dm1(int &tn, double &T, Vector &Hxc,Vector&Hext, double &g_J, Vector &ABC,          // Calculates the transition
 *                 char **sipffilename, ComplexVector & m1, float &delta,                    //   matrix elements of the magnetic moment operator
 *                 ComplexMatrix &est)
 *   void Lcalc(Vector &J, double *T, Vector &Hxc,Vector&Hext, double *gJ, Vector &ABC,       // Calculates the meanfield moment
 *                 char **sipffile,  ComplexMatrix &est)
 *   void LMcalc(Matrix &J, Vector &T, Vector &Hxc,Vector&Hext, double *gJ, Vector &ABC,       // Calculates the meanfield moment
 *                 char **sipffile,  ComplexMatrix &est)
 *   int dL1(int &tn, double &T, Vector &Hxc,Vector&Hext, double &g_J, Vector &ABC,          // Calculates the transition
 *                 char **sipffilename, ComplexVector & L1, float &delta,                    //   matrix elements of the angular momentum  operator
 *                 ComplexMatrix &est)
 *   void Scalc(Vector &J, double *T, Vector &Hxc,Vector&Hext, double *gJ, Vector &ABC,       // Calculates the meanfield moment
 *                 char **sipffile,  ComplexMatrix &est)    //
 *   void SMcalc(Matrix &J, Vector &T, Vector &Hxc,Vector&Hext, double *gJ, Vector &ABC,       // Calculates the meanfield moment
 *                 char **sipffile,  ComplexMatrix &est)    //
 *   int dS1(int &tn, double &T, Vector &Hxc,Vector&Hext, double &g_J, Vector &ABC,          // Calculates the transition
 *                 char **sipffilename, ComplexVector & S1, float &delta,                    //   matrix elements of the spin operator
 *                 ComplexMatrix &est)
 *   void Icalc_parameter_storage_matrix_init(ComplexMatrix *est,Vector &gjmbHxc,Vector&Hext  // initialises parameter storage matrix 
 *                 double *g_J,double &T,Vector &ABC,char **sipffilename)                     // for Icalc
 *   void estates(ComplexMatrix *est, Vector &gjmbHxc,Vector&Hext, double *g_J, double &T,    // Calculates the energy and wavefunctions
 *                 Vector &ABC, char **sipffilename)                                          //   "estates" matrix
 *   bool get_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l,                    // Calculates/loads the Q_q operators
 *                 std::vector<double> &Jvec)                                                 //   for beyond dipole calculations
 *   void save_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l,                   // Saves the Q_q operators to a file
 *                 std::vector<double> &Jvec)                                                 //
 *   void mq(ComplexVector &Mq, double &th, double &ph, double &J0, double &J2,               // Calculates the thermal expectation
 *                 double &J4, double &J6, ComplexMatrix &est)                                //   of the magnetisation density
 *   int dv1calc(int &tn, double &th, double &ph, double &J0, double &J2,                     // Calculates the transition matrix
 *                 double &J4, double &J6, ComplexMatrix &est, double &T,                     //   elements beyond the dipole
 *                 ComplexVecto & v1)                                                         //   approximation.
 *   void spindensity_coeff(Vector & mom, int & xyz,double *T, Vector &gjmbHxc,Vector&Hext,   // Calc. coeffs. of expansion of spindensity 
 *                 double *gJ,Vector &ABC, char **sipffile, ComplexMatrix &est)               //   in terms of Zlm R^2(r) at given T / H_eff
 *   void orbmomdensity_coeff(Vector & mom,int & xyz, double *T, Vector &gjmbHxc,Vector&Hext, // Calc. coeffs. of expansion of orbital moment density
 *                 double *gJ,Vector &ABC, char **sipffile, ComplexMatrix &est)               //   in terms of Zlm F(r) at given T / H_eff
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
      case S:
         switch(n)
         {
            case 1: return fstates_t(1,S,"2S");                           // s1    2S
            default: std::cerr << "Invalid number of electrons\n"; exit(1); 
         }
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
            case 2: return fstates_t(2,F,2,"3F");                         // d2    3F   2 (11)
            case 3: return fstates_t(3,F,3,"4F");                         // d3    4F   3 (11)
            case 4: return fstates_t(4,D,4,"5D");                         // d4    5D   4 (10)
            case 5: return fstates_t(5,S,5,"6S");                         // d5    6S   5 (00)
            default: std::cerr << "Invalid number of electrons\n"; exit(1); 
         }
      case F:
         switch(n)
         {
            case 1: U.set(1,0); return fstates_t(1,F,1,U,"2F");           // f1    2F   1 100 10
            case 2: U.set(1,1); return fstates_t(2,H,2,U,"3H");           // f2    3H   2 110 11
            case 3: U.set(2,0); return fstates_t(3,I,3,U,"4I");           // f3    4I   3 111 20
            case 4: U.set(2,0); return fstates_t(4,I,4,U,"5I");           // f4    5I   4 111 20
            case 5: U.set(1,1); return fstates_t(5,H,5,U,"6H");           // f5    6H   5 110 11
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
// Calculates the spin anisotropy interactions D Sa^2 where a=x,y,z
// --------------------------------------------------------------------------------------------------------------- //
void icf_DS2(sMat<double> &Hcf, icpars &pars)
{
   int n = pars.n; orbital e_l = pars.l;
   fstates_t gs = hunds_gs(n, e_l);
   int J2min, J2max, ns=0, L2=2*abs(gs.L), S2=gs.S2;
   J2min = abs(L2-S2); J2max = L2+S2;
   sMat<double> rm(J2min-J2max+1,J2min-J2max+1);

   // Calculates the reduced matrix element
   for (int J2=J2min; J2<=J2max; J2+=2) 
   {
      ns+=J2+1; 
      for (int J2p=J2min; J2p<=J2max; J2p+=2) 
         rm(J2-J2min,J2p-J2min) = pow(-1.,(S2+L2+J2p)/2.+1.) * sqrt((S2+1.)*(J2+1.)*(J2p+1.)*(S2/2.)*(S2/2.+1.)) * sixj(S2,J2,L2,J2p,S2,2);
   }

   // Determines the states indices.
   int ins=0; std::vector<int> J2(ns), mJ2(ns), irm(ns);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2)
     for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2)
     {
        J2[ins]=iJ2; mJ2[ins]=imJ2; irm[ins]=iJ2-J2min; ins++;
     }

   // Loops through the Cartesian directions, expanding reduced matrix by Wigner-Eckart theorem
   double elp, elm, el0, sqrt2=sqrt(2);
   sMat<double> Sx(ns,ns), Sy(ns,ns), Sz(ns,ns);
   if((fabs(pars.Dx2)+fabs(pars.Dy2))>SMALL)                  // Sx or Sy
   {
      for (int i=0; i<ns; i++) for(int j=0; j<ns; j++)
      {
         elm = rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],-2,mJ2[j]);
         elp = rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],2,mJ2[j]);
         if(fabs(pars.Dx2)>SMALL) { el0=elm-elp; if(fabs(el0)>SMALL) Sx(i,j) = el0/sqrt2; }
         if(fabs(pars.Dy2)>SMALL) { el0=elm+elp; if(fabs(el0)>SMALL) Sy(i,j) = el0/sqrt2; }
      }
      if(fabs(pars.Dx2)>SMALL) { Sx*=Sx; rmzeros(Sx); Hcf += Sx * (pars.Dx2*pars._econv); }
      if(fabs(pars.Dy2)>SMALL) { Sy*=Sy; rmzeros(Sy); Hcf -= Sy * (pars.Dy2*pars._econv); }
   }
   if((fabs(pars.Dz2)>SMALL))                                 // Sz
   {
      for (int i=0; i<ns; i++) for(int j=0; j<ns; j++)
      {
         el0 = rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],0,mJ2[j]); 
         if(fabs(el0)>SMALL) Sz(i,j) = el0;
      }
      Sz*=Sz; rmzeros(Sz); Hcf += Sz * (pars.Dz2*pars._econv); 
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the Hamilton matrix for the lowest Coulombic term (manifold of constant L and S)
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> icf_hmltn(sMat<double> &Hcfi, icpars &pars)
{
   int n = pars.n; double xi = pars._xi; orbital e_l = pars.l;
   if (e_l>3||e_l<0) { std::cerr << "Sorry only s-, p-, d- and f-electrons supported at present\n"; exit(0); }
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
   double p = ( pow(-1.,(double)abs(e_l))*(2.*e_l+1.) );
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
         rmso   += pow(-1.,(double)(pL2+L2)/2.+e_l) * sixj(L2,2,L2,2*e_l,pL2,2*e_l) * 
		   pow(-1.,(double)(pS2+S2+1)/2.)   * sixj(S2,2,S2,1,pS2,1) * cfps[k].cfp * cfps[k].cfp;
         rmU[2] += pow(-1.,(double)(pL2+L2)/2.+e_l) * sixj(L2,4,L2,2*e_l,pL2,2*e_l) * cfps[k].cfp * cfps[k].cfp * (L2+1.);
         rmU[4] += pow(-1.,(double)(pL2+L2)/2.+e_l) * sixj(L2,8,L2,2*e_l,pL2,2*e_l) * cfps[k].cfp * cfps[k].cfp * (L2+1.);
         rmU[6] += pow(-1.,(double)(pL2+L2)/2.+e_l) * sixj(L2,12,L2,2*e_l,pL2,2*e_l)* cfps[k].cfp * cfps[k].cfp * (L2+1.); //}
      }
    //if(n==(4*e_l+1)) { rmU[2]=1./n; rmU[4]=1./n; rmU[6]=1./n; }
      for(int ik=2; ik<=6; ik+=2) rmU[ik] *= n * threej(2*e_l,2*ik,2*e_l,0,0,0) * p;
   }
   else  // Single electron
   {
   // rmso = sqrt(3/2.)*sqrt(e_l*(e_l+1)*(2*e_l+1));             // s=1/2 substituted into eqn 4-12 of Judd 1963
      rmso = 1./((L2+1.)*(S2+1.));
      rmU[2] = threej(2*e_l,4,2*e_l,0,0,0) * p;                  // See Judd 1963, Eqn 5-13. with U^k=V^k/sqrt(2k+1)
      rmU[4] = threej(2*e_l,8,2*e_l,0,0,0) * p;
      rmU[6] = threej(2*e_l,12,2*e_l,0,0,0)* p;
   }

   double elp, elm; int k, q, iq;
   for (int i=0; i<ns; i++) 
      for (int j=0; j<ns; j++)
      {
         if(i==j)  // Selection rule for Spin-orbit operator J=J', mJ=mJ'
         {
            Hcf(i,i) += n*xi * pow(-1.,(S2+L2+J2[i])/2.) * sixj(S2,2,S2,L2,J2[i],L2) * (L2+1)*(S2+1) * sqrt( (9./6)*e_l*(e_l+1)*(2*e_l+1) ) * rmso;
         }
         for(k=2; k<=6; k+=2) for(iq=0; iq<(2*k+1); iq++)
         {
            q = iq-k; if(fabs(pars.B(k,q))<1e-10) continue;
            if(q==0)
            {
               elp = pow(-1.,(S2+L2+J2[i])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * sixj(J2[j],2*k,J2[i],L2,S2,L2) * rmU[k] 
                        * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2*k,J2[j],0-mJ2[i],0,mJ2[j]);
               Hcf(i,j) += elp * pars.B(k,q);
            }
            else
            {
               elp = pow(-1.,(S2+L2+J2[i])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * sixj(J2[j],2*k,J2[i],L2,S2,L2) * rmU[k];
               elm = elp * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2*k,J2[j],0-mJ2[i],-2*abs(q),mJ2[j]);
               elp *=      pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2*k,J2[j],0-mJ2[i], 2*abs(q),mJ2[j]);
               if(q<0)
                 Hcfi(i,j)+= (elm-elp*pow(-1.,q)) * pars.B(k,q); 
               else
                 Hcf(i,j) += (elm+elp*pow(-1.,q)) * pars.B(k,q); 
            }
         }
      }
   // Calculates the spin-anisotropy terms D[xzy]2(S[xyz].S[xyz]) if needed.
   if((fabs(pars.Dx2)+fabs(pars.Dy2)+fabs(pars.Dz2))>SMALL) icf_DS2(Hcf,pars);
   return Hcf;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the Multipolar operator matrices
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> icf_ukq(int n, int k, int q, orbital e_l)
{
   if (e_l>3||e_l<0) { std::cerr << "Sorry only s-, p-, d- and f-electrons supported at present\n"; exit(0); }
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
   double rmU=0., p = pow(-1.,(double)abs(e_l))*(2.*e_l+1.);
   if (n>1) // && n<(4*e_l+1))
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
         rmU += pow(-1.,(double)(pL2+L2)/2.+e_l) * sixj(L2,2*k,L2,2*e_l,pL2,2*e_l) * cfps[kk].cfp * cfps[kk].cfp * (L2+1.);
      }
      rmU *= n * threej(2*e_l,2*k,2*e_l,0,0,0) * p;
   }
   else  // Single electron or single hole
      rmU = threej(2*e_l,2*k,2*e_l,0,0,0) * p;                   // See Judd 1963, Eqn 5-13. with U^k=V^k/sqrt(2k+1)

   double elp, elm;
   if(q==0)
      for (int i=0; i<ns; i++)
         for (int j=0; j<ns; j++)
         {
            Hcf(i,j) = pow(-1.,(S2+L2+J2[i])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * sixj(J2[j],2*k,J2[i],L2,S2,L2) * rmU
                       * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2*k,J2[j],0-mJ2[i],0,mJ2[j]);
         }
   else
      for (int i=0; i<ns; i++)
         for (int j=0; j<ns; j++)
         {
            elp = pow(-1.,(S2+L2+J2[i])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * sixj(J2[j],2*k,J2[i],L2,S2,L2) * rmU;
            elm = elp * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2*k,J2[j],0-mJ2[i],-2*abs(q),mJ2[j]);
            elp *=      pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2*k,J2[j],0-mJ2[i], 2*abs(q),mJ2[j]);
            if(q<0)
              Hcf(i,j) += (elm-elp*pow(-1.,q));
            else
              Hcf(i,j) += (elm+elp*pow(-1.,q));
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
            rm(J2-J2min,J2p-J2min) = pow(-1.,(S2+L2+J2p)/2.) * sqrt((S2+1.)*(J2+1.)*(J2p+1.)*(S2/2.)*(S2/2.+1.)) * sixj(J2p,2,J2,S2,L2,S2);
      else                                              // Lx, Ly or Lz
         for (int J2p=J2min; J2p<=J2max; J2p+=2) 
            rm(J2-J2min,J2p-J2min) = pow(-1.,(S2+L2+J2)/2.)  * sqrt((L2+1.)*(J2+1.)*(J2p+1.)*(L2/2.)*(L2/2.+1.)) * sixj(J2p,2,J2,L2,S2,L2);
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
   else if (ind>1 && ind<4)                             // Sy or Ly
      for (int i=0; i<ns; i++) for(int j=0; j<ns; j++)
      {
         elm = rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],-2,mJ2[j]);
         elm+= rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],2,mJ2[j]);
         if(fabs(elm)>SMALL) mu(i,j)=elm/sqrt2;
       //mu(i,j) = (elm+elp)/sqrt2;
      }
   else if (ind>3 && ind<6)                             // Sz or Lz
      for (int i=0; i<ns; i++) for(int j=0; j<ns; j++) 
      {
         elm = rm(irm[i],irm[j]) * pow(-1.,(J2[i]-mJ2[i])/2.) * threej(J2[i],2,J2[j],-mJ2[i],0,mJ2[j]); 
         if(fabs(elm)>SMALL) mu(i,j)=elm; 
      }
   else { std::cerr << "icf_mumat: Error, index > 5 or index < 0.\n"; exit(EXIT_FAILURE); }

   return mu;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to initialise the storage matrix "est" for Icalc
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
void Icalc_parameter_storage_matrix_init(
                      ComplexMatrix *est,   // Output Eigenstates matrix (row 0: real==Eigenvalues;imag==population)
                      Vector &gjmbHxc,      // Input  vector of exchange fields (meV) 
                      Vector &Hext,         // Input  vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,     // Input  Lande g-factor
                      double * /*T*/,       // Input  temperature
 /* Not Used */       Vector & /*ABC*/,     // Input  Vector of parameters from single ion property file
                      char **sipffilename)  // Input  Single ion properties filename
{
   // Parses the input file for parameters
   icpars pars;
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

   // If we just want a blank estates matrix for later use (e.g. in Icalc)
   int nelm = (gjmbHxc.Elements()<6) ? 6 : gjmbHxc.Elements();
   int nfact = (int)ceil(sqrt(nelm+1));
   int Hsz = icf_getdim(pars)*nfact;
   (*est) = ComplexMatrix(0,Hsz,0,Hsz);
   (*est)(0,0) = complex<double> (pars.n, pars.l);
   (*est)(0,1) = complex<double> (nfact, nelm);
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates thermal expectation values 
// --------------------------------------------------------------------------------------------------------------- //
void icf_expJ(icpars &pars, ComplexMatrix &est, complexdouble *zV, double *vE, double *T, Vector &J, double *lnZ, double *U)
{
   int K[] = {-1,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Q[] = {-1,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {-1,0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
   int Hsz=icf_getdim(pars), i, ix, iy, incx=1;
   int nfact = (int)ceil(sqrt(J.Hi()-J.Lo()+2));

   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';
   int Esz, ind_j;
   //std::vector< std::vector<double> > matel;
   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   std::vector<double> E, me, eb; /*matel.clear();*/ E.reserve(Hsz);
   for(Esz=0; Esz<Hsz; Esz++) { E.push_back(vE[Esz]-vE[0]); if(exp(-E[Esz]/(KB**T))<DBL_EPSILON || vE[Esz+1]==0) break; }

   if (*T<0)
   {
      Esz = (int)(-*T);
      printf ("Temperature T=%g<0: please choose probability distribution for the -T=%i lowest energy states by hand\n",*T,(int)(-*T));
                        printf("Number   Excitation Energy\n");
      for (ind_j=0; ind_j<Esz; ++ind_j) printf ("%i    %4.4g meV\n",ind_j+1,E[ind_j]);
   }  // MR 10.9.2010

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
   sMat<double> Hcf;

   // Checks if have been called from spins, where _par_storage() was set up for just the number of operators 
   //    needed in the MF loop but now we need to calculate expectation values for all multipolar operators
   int oldJhi = J.Hi();
   if(nfact != (int)real(est[0][1])) {
      nfact = (int)real(est[0][1]); oldJhi = (int)imag(est[0][1]); }

   for(int iJ=J.Lo(); iJ<=J.Hi(); iJ++)
   {
      me.assign(Esz,0.); J[iJ]=0.;
      // Using the above reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
      if(iJ>6 && (K[iJ]%2==1 || K[iJ]>2*pars.l)) { /*matel.push_back(me);*/ continue; }
      {
         if(iJ<=oldJhi)
         {
            iy = (iJ-J.Lo()+1)/nfact; ix = (iJ-J.Lo()+1)-iy*nfact;
            for(i=1; i<=Hsz; i++) memcpy(&zJmat[(i-1)*Hsz],&est[i+ix*Hsz][1+iy*Hsz],Hsz*sizeof(complexdouble));
         }
         else
         {
            if(iJ<=6) Hcf = icf_mumat(pars.n, iJ-1, pars.l);     // Calculates Sx,Lx etc
            else      Hcf = icf_ukq(pars.n,K[iJ],Q[iJ],pars.l);  // Calculates multipolar operator matrices 
            std::vector< std::vector<int> > u = Hcf.find();
            memset(zJmat,0.,Hsz*Hsz*sizeof(complexdouble));
            complexdouble x;
            if(im[iJ]==1) for (int j=0; j<(int)u.size(); j++) { x.r=0.; x.i=Hcf(u[j][0],u[j][1]); zJmat[Hcf.nr()*u[j][1]+u[j][0]] = x; }
            else          for (int j=0; j<(int)u.size(); j++) { x.r=Hcf(u[j][0],u[j][1]); x.i=0.; zJmat[Hcf.nr()*u[j][1]+u[j][0]] = x; }
         }

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
               { 
                  char instr[MAXNOFCHARINLINE];
                  printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
                  if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { fprintf(stderr,"Error reading input\n"); exit(1); }
                  eb[ind_j]=strtod(instr,NULL);
               }
               else
               { 
                  eb[ind_j] = exp(-E[ind_j]/(KB**T)); 
               } 
               J[iJ] += me[ind_j]*eb[ind_j]; 
               Z += eb[ind_j]; 
               *U += (E[ind_j]+vE[0])*eb[ind_j];
//MRend 10.9.2010
            }
            else  // Rest of the runs only calculate the new matrix elements
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
           void IMcalc(Matrix &Jret,          // Output single ion momentum vector <Ja>,<Jb>,<Jc>, etc.
                      Vector&T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      Vector &lnZ,        // Output scalar logarithm of partition function
                      Vector &U,          // Output scalar internal energy 
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   // sum exchange field and external field
   Vector gjmbH(1,(Hxc.Hi()<6) ? 6 : Hxc.Hi()); gjmbH=0;
   if(gjmbH.Hi()==Hxc.Hi()) gjmbH=Hxc; else for(int i=1; i<=(gjmbH.Hi()<Hxc.Hi()?gjmbH.Hi():Hxc.Hi()); i++) gjmbH[i]=Hxc[i];
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   double dum; dum=gjmbH(2);gjmbH(2)=gjmbH(4);gjmbH(4)=gjmbH(5);gjmbH(5)=gjmbH(3);gjmbH(3)=dum;
   Vector J(1,(Hxc.Hi()<6) ? 6 : Hxc.Hi()); J=0;
   // --------------------------------------------------------------------
   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(Hext(1))>DBL_EPSILON || fabs(Hext(2))>DBL_EPSILON || fabs(Hext(3))>DBL_EPSILON)
   {
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH(2)+=MUB*Hext(1); gjmbH(1)+=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH(4)+=MUB*Hext(2); gjmbH(3)+=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH(6)+=MUB*Hext(3); gjmbH(5)+=GS*MUB*Hext(3); }
   }
   // Parses the input file for parameters
   icpars pars; 
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

   // Converts the Jij parameters if necessary
   std::vector<double> vgjmbH(J.Hi()+1,0.); 
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc();
      for(int i=J.Lo(); i<=J.Hi(); i++) vgjmbH[i] = -gjmbH[i]*pars.jijconv[i]; }
   else
   #endif
      for(int i=J.Lo(); i<=J.Hi(); i++) vgjmbH[i] = -gjmbH[i];  // Vector of Exchange + External field to be added to matrix in lines 628, 637

   //          0 1 2 3 4 5 6  7  8 91011 12 13 1415161718 19 20 21 222324252627 28 29 30 31 32333435363738 39 40 41 42 43 4445464748495051
   int K[] = {-1,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Q[] = {-1,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {-1,0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};

   // Calculates the IC Hamiltonian matrix
   int nfact = (int)ceil(sqrt(J.Hi()-J.Lo()+2));
   int i,k,q,Hsz=icf_getdim(pars),esz=Hsz*nfact;
   int ix, iy, incx=1;
   complexdouble *H=0, *zM=0;
   bool Hicnotcalc = false;
   std::vector<double> parval; parval.reserve(35);
   parval.push_back(pars.xi); for(k=2; k<=(2*pars.l); k+=2) for(q=-k; q<=k; q++) parval.push_back(pars.B(k,q));
   if(parval.size()%2==1) parval.push_back(0.);

// if((est.Cols()!=(esz+1) || est.Rows()!=(esz+1))) Hicnotcalc = true;
   if(real(est[0][0])==(double)pars.n && imag(est[0][0])==(double)pars.l)  // Hic previously calculated
   {
      for(i=0; i<(int)(parval.size()/2); i++) if(real(est[0][i+2])!=parval[2*i] || imag(est[0][i+2])!=parval[2*i+1]) { 
         Hicnotcalc = true; break; }
   }
   else Hicnotcalc = true;

   int oldJhi = J.Hi();
   if((est.Rhi()!=esz||est.Chi()!=esz))   // Probably from spins: need all multipolar operators, not just those used in MF loop.
   {
      nfact = (int)real(est[0][1]); esz = Hsz*nfact;
      oldJhi = (int)imag(est[0][1]);
   }

   if(Hicnotcalc || pars.l==S)
   {
      sMat<double> Hcfi, Hcf = icf_hmltn(Hcfi, pars); Hcf/=MEV2CM; Hcfi/=MEV2CM; H = zmat2f(Hcf,Hcfi);
      if(est.Rhi()!=esz||est.Chi()!=esz) {
         std::cerr << "ERROR module icf1ion - Icalc: Hsz recalculation does not agree with eigenstates matrix dimension\n"; exit(EXIT_FAILURE); }
      else if(esz>(int)parval.size()/2+1)
      {
         est[0][0] = complex<double> (pars.n,pars.l);
         for(i=0; i<(int)(parval.size()/2); i++) est[0][i+2] = complex<double> (parval[2*i],parval[2*i+1]);
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
      complex<double> a(vgjmbH[ind],0.);
      if(ind<=oldJhi)
      {
         iy = ind/nfact; ix = ind-iy*nfact;
         for(i=1; i<=Hsz; i++) F77NAME(zaxpy)(&Hsz,(complexdouble*)&a,(complexdouble*)&est[i+ix*Hsz][1+iy*Hsz],&incx,&H[(i-1)*Hsz],&incx);
      }
      else
      {
         sMat<double> Hcf, Hcfi(Hsz,Hsz);
         if(ind<=6) {       // Calculates Sx,Lx etc and copies them to est too. 
            Hcf = icf_mumat(pars.n, ind-1, pars.l); if(ind==3 || ind==4) zM = zmat2f(Hcfi,Hcf); else zM = zmat2f(Hcf,Hcfi); }
         else       {       // Calculates multipolar operator matrices and copies them to est. 
            Hcf = icf_ukq(pars.n,K[ind],Q[ind],pars.l); if(im[ind]==1)   zM = zmat2f(Hcfi,Hcf); else zM = zmat2f(Hcf,Hcfi); }
         for(i=1; i<=Hsz; i++) F77NAME(zaxpy)(&Hsz,(complexdouble*)&a,&zM[(i-1)*Hsz],&incx,&H[(i-1)*Hsz],&incx); free(zM);
      }
   }

   // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
   double *vE = new double[Hsz]; complexdouble *zV = new complexdouble[Hsz*Hsz];
   int info = ic_diag(Hsz,H,zV,vE); free(H);
   if(info!=0) { std::cerr << "icf1ion - Error diagonalising, info==" << info << "\n"; delete[]vE; vE=0; delete[]zV; zV=0; exit(EXIT_FAILURE); }


for(int Ti=1;Ti<=T.Hi();++Ti){
         icf_expJ(pars,est,zV,vE,&T(Ti),J,&lnZ(Ti),&U(Ti));
    //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   dum=J(2);J(2)=J(3);J(3)=J(5);J(5)=J(4);J(4)=dum;
   for(i=Jret.Rlo(); i<=Jret.Rhi(); i++) Jret(i,Ti) = J(i);
                          }
         delete[]vE; delete[]zV;
     
     // --------------------------------------------------------------------

}

extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void Icalc(Vector &Jret,          // Output single ion momentum vector <Ja>,<Jb>,<Jc>, etc.
                      double *T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      double *lnZ,        // Output scalar logarithm of partition function
                      double *U,          // Output scalar internal energy 
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{Matrix JM(1,Jret.Hi(),1,1);double * d=NULL;Vector dd;
 for(int i=1;i<=Jret.Hi();++i)JM(i,1)=Jret(i);
 Vector TT(1,1);TT(1)=*T;
 Vector lnZZ(1,1);lnZZ(1)=*lnZ;
 Vector UU(1,1);UU(1)=*U;
 IMcalc(JM,TT,Hxc,Hext,d,dd,sipffilename,lnZZ,UU,est);
 *U=UU(1);*lnZ=lnZZ(1);*T=TT(1);
 for(int i=1;i<=Jret.Hi();++i)Jret(i)=JM(i,1);
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the magnetic moment at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void mcalc(Vector &mom,        // Output magnetic moment (mub)
                      double *T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double *g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   Vector J(1,6);
   double lnZ, U;
   Icalc(J,T,Hxc,Hext,g_J,ABC,sipffilename,&lnZ,&U,est);
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   //mom(1)=GS*J(1)+J(2);
   //mom(2)=GS*J(3)+J(4);
   //mom(3)=GS*J(5)+J(6);
   mom(1)=GS*J(1)+J(4);
   mom(2)=GS*J(2)+J(5);
   mom(3)=GS*J(3)+J(6);
}
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void mMcalc(Matrix &mom,        // Output magnetic moment (mub)
                      Vector &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   Matrix J(1,6,1,T.Hi()); Vector ABC; 
   double gJ=0.;Vector lnZ(1,T.Hi()), U(1,T.Hi());
   IMcalc(J,T,Hxc,Hext,&gJ,ABC,sipffilename,lnZ,U,est);
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   //mom(1)=GS*J(1)+J(2);
   //mom(2)=GS*J(3)+J(4);
   //mom(3)=GS*J(5)+J(6);
   for(int Ti=1;Ti<=T.Hi();++Ti){
   mom(1,Ti)=GS*J(1,Ti)+J(4,Ti);
   mom(2,Ti)=GS*J(2,Ti)+J(5,Ti);
   mom(3,Ti)=GS*J(3,Ti)+J(6,Ti); }
}
// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the <L> at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void Lcalc(Vector &L,          // Output magnetic moment (mub)
                      double *T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double *g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   Vector J(1,6); 
   double lnZ, U;
   Icalc(J,T,Hxc,Hext,g_J,ABC,sipffilename,&lnZ,&U,est);
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   //L(1)=J(2);
   //L(2)=J(4);
   //L(3)=J(6);
   L(1)=J(4);
   L(2)=J(5);
   L(3)=J(6);
}
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void LMcalc(Matrix &L,          // Output magnetic moment (mub)
                      Vector &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   Matrix J(1,6,1,T.Hi());Vector ABC; 
   double gJ=0.; Vector lnZ(1,T.Hi()), U(1,T.Hi());
   IMcalc(J,T,Hxc,Hext,&gJ,ABC,sipffilename,lnZ,U,est);
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   //L(1)=J(2);
   //L(2)=J(4);
   //L(3)=J(6);
   for(int Ti=1;Ti<=T.Hi();++Ti){
   L(1,Ti)=J(4,Ti);
   L(2,Ti)=J(5,Ti);
   L(3,Ti)=J(6,Ti);}
}
// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the <S> at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void Scalc(Vector &S,          // Output magnetic moment (mub)
                      double *T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double *g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   Vector J(1,6); 
   double lnZ, U;
   Icalc(J,T,Hxc,Hext,g_J,ABC,sipffilename,&lnZ,&U,est);
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   //S(1)=J(1);
   //S(2)=J(3);
   //S(3)=J(5);
   S(1)=J(1);
   S(2)=J(2);
   S(3)=J(3);
}
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void SMcalc(Matrix &S,          // Output magnetic moment (mub)
                      Vector &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   Matrix J(1,6,1,T.Hi());Vector  ABC; 
   double gJ=0.;Vector lnZ(1,T.Hi()), U(1,T.Hi());
   IMcalc(J,T,Hxc,Hext,&gJ,ABC,sipffilename,lnZ,U,est);
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   //S(1)=J(1);
   //S(2)=J(3);
   //S(3)=J(5);

for(int Ti=1;Ti<=T.Hi();++Ti){
   S(1,Ti)=J(1,Ti);
   S(2,Ti)=J(2,Ti);
   S(3,Ti)=J(3,Ti);}
}
// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the eigenstates of Hic+effective_mean_fields
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
         void estates(ComplexMatrix *est, // Output Eigenstates matrix (row 0: real==Eigenvalues;imag==population)
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (meV) 
 /* Not Used */       double & /*g_J*/,   // Input  Lande g-factor
                      double &T,          // Input  temperature
 /* Not Used */       Vector & /*ABC*/,   // Input  Vector of parameters from single ion property file
                      char **sipffilename)// Input  Single ion properties filename
{  // sum exchange field and external field
   Vector gjmbH(1,(Hxc.Hi()<6) ? 6 : Hxc.Hi()); gjmbH=0;
   if(gjmbH.Hi()==Hxc.Hi()) gjmbH=Hxc; else for(int i=1; i<=(gjmbH.Hi()<Hxc.Hi()?gjmbH.Hi():Hxc.Hi()); i++) gjmbH[i]=Hxc[i];
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   double dum; dum=gjmbH(2);gjmbH(2)=gjmbH(4);gjmbH(4)=gjmbH(5);gjmbH(5)=gjmbH(3);gjmbH(3)=dum;
   // --------------------------------------------------------------------------

   // Calculates the Zeeman term if magnetic field is not zero 
  if(fabs(Hext(1))>DBL_EPSILON || fabs(Hext(2))>DBL_EPSILON || fabs(Hext(3))>DBL_EPSILON)
   {
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH(2)+=MUB*Hext(1); gjmbH(1)+=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH(4)+=MUB*Hext(2); gjmbH(3)+=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH(6)+=MUB*Hext(3); gjmbH(5)+=GS*MUB*Hext(3); }
   }
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
   std::vector<double> vgjmbH(gjmbH.Hi()+1,0.);
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc();
      for(i=gjmbH.Lo(); i<=gjmbH.Hi(); i++) vgjmbH[i] = -gjmbH[i]*pars.jijconv[i]; }
   else
   #endif
      for(i=gjmbH.Lo(); i<=gjmbH.Hi(); i++) vgjmbH[i] = -gjmbH[i];

   // Calculates the Hamiltonian matrix
   sMat<double> mat, Hcfi, Hcf; Hcf = icf_hmltn(Hcfi, pars); Hcf/=MEV2CM; Hcfi/=MEV2CM;
   for(int ind=gjmbH.Lo(); ind<=gjmbH.Hi(); ind++)
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
      for(i=0; i<Hsz; i++)
      {  
         //printf("\n");
         for(j=0; j<Hsz; j++) 
         {
            (*est)(i+1,j+1) = complex<double> (rV[j+i*Hsz], 0.);
//          printf("%6.3f %+6.3f i  ",rV[i+j*Hsz],0.0);
//          if(rV[i+j*Hsz] != (*est)[i+1][j+1]) { fprintf(stderr,"compiler problem: bad memory mapping of vectors\n"); exit(EXIT_FAILURE); }
         }
      }
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
          int du1calc(int &tn,            // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (meV) 
 /* Not Used */       double & /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexVector &u1ret,  // Output eigenvector u1
                      float &delta,       // Output transition energy
                      int &n, int &nd,
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{  // sum exchange field and external field
   Vector gjmbH(1,(Hxc.Hi()<6) ? 6 : Hxc.Hi()); gjmbH=0;
   if(gjmbH.Hi()==Hxc.Hi()) gjmbH=Hxc; else for(int i=1; i<=(gjmbH.Hi()<Hxc.Hi()?gjmbH.Hi():Hxc.Hi()); i++) gjmbH[i]=Hxc[i];
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   double dum; dum=gjmbH(2);gjmbH(2)=gjmbH(4);gjmbH(4)=gjmbH(5);gjmbH(5)=gjmbH(3);gjmbH(3)=dum;
   ComplexVector u1(1,(Hxc.Hi()<6) ? 6 : Hxc.Hi()); u1=0; u1(1)=u1ret(1);
   // --------------------------------------------------------------------------

   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(Hext(1))>DBL_EPSILON || fabs(Hext(2))>DBL_EPSILON || fabs(Hext(3))>DBL_EPSILON)
   {
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH(2)+=MUB*Hext(1); gjmbH(1)+=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH(4)+=MUB*Hext(2); gjmbH(3)+=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH(6)+=MUB*Hext(3); gjmbH(5)+=GS*MUB*Hext(3); }
   }
  
   int i,j,k;

   int K[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};

   int sz = gjmbH.Hi();
   sMat<double> zeroes(est.Rows()-1,est.Cols()-1), op;
   complexdouble *zJmat=0, *zt=0, zme; zme.r=0; zme.i=0.; 
   std::vector<complexdouble> zij(sz,zme);//, zji(6,zme);
   std::vector<double> u(sz+1),iu(sz+1);
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';
   double Z=0., therm;

   // check if printout should be done and make tn positive
   int pr=0; if (tn<0) { pr=1; tn*=-1; }
   double ninit=u1[1].real();
   double pinit=u1[1].imag();

   // Copies the already calculated energy levels / wavefunctions from *est
   if(est.Rows()!=est.Cols()) { std::cerr << "du1calc(): Input rows and columns of eigenstates matrix don't match.\n"; return 0; }
   int Hsz = est.Rows()-1, iJ, incx = 1;
   j=0; k=0; for(i=0; i<Hsz; ++i) { for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   n=i;nd=j;
   double maxE=delta;
   if((delta=(est[0][j+1].real()-est[0][i+1].real()))<=maxE)
   {
      double *en = new double[Hsz]; for(k=0; k<Hsz; k++) en[k] = est[0][k+1].real();
//    iceig VE(Hsz,en,(complexdouble*)&est[1][0],1);
//    for(int ii=1; ii<=Hsz; ii++) memcpy(&zV[(ii-1)*Hsz],&(*est)[ii][1],Hsz*sizeof(complexdouble)); 

      // Parses the input file for parameters
      icpars pars; const char *filename = sipffilename[0];
      ic_parseinput(filename,pars);

      // Calculates the transition matrix elements:
      //    u1 = <i|Ja|j> * sqrt[(exp(-Ei/kT)-exp(-Ej/kT)) / Z ]   if delta > small
      //    u1 = <i|Ja-<Ja>|j> * sqrt[(exp(-Ei/kT)) / kTZ ]             if delta < small (quasielastic scattering)
      for(iJ=0; iJ<sz; iJ++)
      {
         if(iJ<6) op = icf_mumat(pars.n, iJ, pars.l); else op = icf_ukq(pars.n,K[iJ],Q[iJ],pars.l); 
         if(im[iJ]==1) zJmat=zmat2f(zeroes,op); else zJmat=zmat2f(op,zeroes);
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, (complexdouble*)&est[j+1][1], &incx, &zbeta, zt, &incx);
         #ifdef _G77 
         F77NAME(zdotc)(&zij[iJ], &Hsz, (complexdouble*)&est[i+1][1], &incx, zt, &incx);
         #else
         zij[iJ] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[i+1][1], &incx, zt, &incx);
         #endif
         free(zJmat); free(zt);
      }

      if(i==j && T>0) // subtract thermal expectation value from zij=zii
      {
//       Vector vJ(gjmbH.Lo(),gjmbH.Hi());
//       icf_expJ(pars,est,(complexdouble*)&est[1][0],en,&T,vJ,&lnZ,&U);
//       for(iJ=0; iJ<sz; iJ++) zij[iJ].r-=vJ[iJ];
         std::vector<double> eb, E; E.reserve(Hsz);
         int Esz; for(Esz=0; Esz<Hsz; Esz++) { E.push_back(en[Esz]-en[0]); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || en[Esz+1]==0) break; }
	 eb.assign(Esz,0.);
         for(iJ=0; iJ<sz; iJ++)
         {
            if(iJ<6) op = icf_mumat(pars.n, iJ, pars.l); else op = icf_ukq(pars.n,K[iJ],Q[iJ],pars.l); 
            if(im[iJ]==1) zJmat=zmat2f(zeroes,op); else zJmat=zmat2f(op,zeroes); double Jj=0.;
            zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
            for(int ind_j=0; ind_j<Esz; ind_j++)
            {  // Calculates the matrix elements <Vi|J.H|Vi>
               F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, (complexdouble*)&est[ind_j+1][1], &incx, &zbeta, zt, &incx);
               #ifdef _G77 
               F77NAME(zdotc)(&zme, &Hsz, (complexdouble*)&est[ind_j+1][1], &incx, zt, &incx);
               #else
               zme = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[ind_j+1][1], &incx, zt, &incx);
               #endif
               // For first run calculate also the partition function and internal energy
               if(iJ==0)
               {
//MR 10.9.2010
                  if (T<0) 
                  {
                     char instr[MAXNOFCHARINLINE];
                     printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
                     if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { fprintf(stderr,"Error reading input\n"); exit(1); }
                     eb[ind_j]=strtod(instr,NULL); 
                  }
//MRend 10.9.2010
                  else
                     eb[ind_j] = exp(-E[ind_j]/(KB*T));  
		  Z+=eb[ind_j];
	       }
               Jj+=zme.r*eb[ind_j];
            }
            zij[iJ].r -= Jj/Z; free(zJmat); free(zt);
         }
      }

      if (T<0) T=-T;

      // Calculates the vector u , iu
      for(iJ=0; iJ<sz; iJ++)
      {  
         u[iJ+1] = zij[iJ].r;
         iu[iJ+1] =zij[iJ].i;
      }

      delta = en[j]-en[i];
      if(delta<-0.000001) {
         std::cerr << "ERROR module ic1ion - du1calc: energy gain delta gets negative\n"; exit(EXIT_FAILURE); }
      if(j==i)delta=-SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc
   
      // Calculates the partition function
      Z=0.; for(iJ=0; iJ<Hsz; iJ++) { therm = exp(-(en[iJ]-en[0])/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }
   
      // do some printout if wishes and set correct occupation factor
      if (delta>SMALL)
      {
         therm = exp(-(en[i]-en[0])/(KB*T)) - exp(-(en[j]-en[0])/(KB*T));
         if(pr==1)
         {
            printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
            printf(" |<%i|Ja|%i>|^2=%6.3f\n |<%i|Jb|%i>|^2=%6.3f\n |<%i|Jc|%i>|^2=%6.3f\n",i+1,j+1,u[1]*u[1]+iu[1]*iu[1],i+1,j+1,u[2]*u[2]+iu[2]*iu[2],i+1,j+1,u[3]*u[3]+iu[3]*iu[3]);
            printf(" |<%i|Jd|%i>|^2=%6.3f\n |<%i|Je|%i>|^2=%6.3f\n |<%i|Jf|%i>|^2=%6.3f\n",i+1,j+1,u[4]*u[4]+iu[4]*iu[4],i+1,j+1,u[5]*u[5]+iu[5]*iu[5],i+1,j+1,u[6]*u[6]+iu[6]*iu[6]);
            printf(" n%i-n%i=%6.3f\n",i,j,therm / Z);
         }
      }
      else
      {
         therm = exp(-(en[i]-en[0])/(KB*T))/(KB*T);    // quasielastic scattering has not wi-wj but wj*epsilon/kT
         if(pr==1)
         {
            printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
            printf(" |<%i|Ja-<Ja>|%i>|^2=%6.3f\n |<%i|Jb-<Jb>|%i>|^2=%6.3f\n |<%i|Jc-<Jc>|%i>|^2=%6.3f\n",i+1,j+1,u[1]*u[1]+iu[1]*iu[1],i+1,j+1,u[2]*u[2]+iu[2]*iu[2],i+1,j+1,u[3]*u[3]+iu[3]*iu[3]);
            printf(" |<%i|Jd-<Jd>|%i>|^2=%6.3f\n |<%i|Je-<Je>|%i>|^2=%6.3f\n |<%i|Jf-<Jf>|%i>|^2=%6.3f\n",i+1,j+1,u[4]*u[4]+iu[4]*iu[4],i+1,j+1,u[5]*u[5]+iu[5]*iu[5],i+1,j+1,u[6]*u[6]+iu[6]*iu[6]);
            printf(" n%i=%6.3f\n",i,(KB*T)*therm/Z);
         }
      }
   
      // multiply matrix Mab by occupation factor
      for(iJ=1; iJ<=u1.Hi(); iJ++)
	    u1(iJ) = complex<double> ( u[iJ]*sqrt(therm/Z), iu[iJ]*sqrt(therm/Z) );
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   complex<double>dum1;dum1=u1(2);u1(2)=u1(3);u1(3)=u1(5);u1(5)=u1(4);u1(4)=dum1;
   for(iJ=1; iJ<=u1ret.Hi(); iJ++) u1ret[iJ] = u1[iJ];

      delete[]en;
   }

   if (ninit>Hsz)ninit=Hsz;
   //if (pinit<SMALL)pinit=SMALL;
   double zsum=0,zi,x;
   // determine number of thermally reachable states
   int noft = 0;
   for(i=0; (i<ninit)&((((x=(est[0][i+1].real()-est[0][1].real())/(KB*fabs(T)))<200)? zi=exp(-x):zi=0)>=(pinit*zsum)); ++i)
   {
      noft += Hsz-i; 
      zsum += zi;
   }

//   for(i=0; (i<ninit)&((zi=(exp(-(est[0][i+1].real()-est[0][1].real())/(KB*fabs(T)))))>(pinit*zsum)); ++i) { noft += Hsz-i+1; zsum += zi; }
// int noft=0;for(i=0;(i<Hsz)&(exp(-(est[0][i+1].real()-est[0][1].real())/(KB*fabs(T)))>SMALL);++i)noft+=Hsz-i-1; // removed MR  6.9.2011 to allow for mcdisp options -ninit -pinit   return noft;

   return noft;
   //return Hsz*(Hsz-1)/2;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements of magnetic moment operator
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
              int dm1(int &tn,            // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexVector & m1, // Output m1 vector (1,3)
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ 
   ComplexVector u1(1,6);int n,nd;
   u1(1) = m1(1);
   int nt = du1calc(tn,T,Hxc,Hext,g_J,ABC,sipffilename,u1,delta,n,nd,est);
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
  // m1(1)=GS*u1(1)+u1(2);
  // m1(2)=GS*u1(3)+u1(4);
  // m1(3)=GS*u1(5)+u1(6);
   m1(1)=GS*u1(1)+u1(4);
   m1(2)=GS*u1(2)+u1(5);
   m1(3)=GS*u1(3)+u1(6);
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements of orbital angular momentum operator
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
              int dL1(int &tn,            // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexVector & L1, // Output L1 vector (1,3)
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ 
   ComplexVector u1(1,6);int n,nd;
   u1(1) = L1(1);
   int nt=du1calc(tn,T,Hxc,Hext,g_J,ABC,sipffilename,u1,delta,n,nd,est);
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   //L1(1)=u1(2);
   //L1(2)=u1(4);
   //L1(3)=u1(6);
   L1(1)=u1(4);
   L1(2)=u1(5);
   L1(3)=u1(6);
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements of spin operator
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
              int dS1(int &tn,            // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexVector & S1, // Output S1 vector (1,3)
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ 
   ComplexVector u1(1,6);int n,nd;
   u1(1) = S1(1);
   int nt=du1calc(tn,T,Hxc,Hext,g_J,ABC,sipffilename,u1,delta,n,nd,est);
    //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   //S1(1)=u1(1);
   //S1(2)=u1(3);
   //S1(3)=u1(5);
   S1(1)=u1(1);
   S1(2)=u1(2);
   S1(3)=u1(3);
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the a(K,K') matrix - the routine returns a true if the matrix is real, and false if it is imaginary
//    The matrix is returned in the argument aKK
// --------------------------------------------------------------------------------------------------------------- //
bool icf_loveseyAKK(sMat<double> &aKK, int K, int Kp, int n, orbital l)
{
   if (l>3||l<0) { std::cerr << "Sorry only s-, p-, d- and f-electrons supported at present\n"; exit(0); }

   fstates_t gs = hunds_gs(n, l);
   int J2min, J2max, ns=0, S2 = gs.S2, L2 = abs(gs.L)*2; std::vector<int> J2;

   // Determines number and angular momentum quantum numbers of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2;
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) J2.push_back(iJ2); ns=(int)J2.size();

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

   rmzeros(aKK); return 1;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the c(K,K') matrix, returning true if matrix is real, false if imaginary
// --------------------------------------------------------------------------------------------------------------- //
bool icf_loveseyCKK(sMat<double> &cKK, int K, int Kp, int n, orbital l)
{
   if (l>3||l<0) { std::cerr << "Sorry only s-, p-, d- and f-electrons supported at present\n"; exit(0); }

   fstates_t gs = hunds_gs(n, l);
   int J2min, J2max, ns=0, S2 = gs.S2, L2 = abs(gs.L)*2; std::vector<int> J2;

   // Determines number and angular momentum quantum numbers of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2;
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) J2.push_back(iJ2); ns=(int)J2.size();

   // Eqn. 3.6.11
   //
   //                  K      1/2                            1/2      1/2+S'+L'
   // C(K,K') = <j >  i  (1/2)     [l,l,S,S',L,L',J',K,K',K']     (-1)
   //             K                                                       _ _
   //                              { 1  K  K' }    ---     _   _          S+L
   //              \/   ( l K l )  { S' L' J' }  n >   (t{|t) (t|}t') (-1)    { S  1  S' }  { L  K  L' }
   //              /\   ( 0 0 0 )  { S  L  J  }    ---                        { s  Sb s  }  { l  Lb l  }
   //                                               t

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

   rmzeros(cKK); return 1;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the transition operator Q_q
// --------------------------------------------------------------------------------------------------------------- //
void icf_loveseyQq(std::vector< sMat<double> > &Qq, int q, int n, orbital l, std::vector<double> &Jvec)
{
   if (l>3||l<0) { std::cerr << "Sorry only s-, p-, d- and f-electrons supported at present\n"; exit(0); }

   double theta = Jvec[0], phi = Jvec[1], J[]={Jvec[2],0,Jvec[3],0,Jvec[4],0,Jvec[5]};
   std::string errormsg("lovesey_Qq(): Unable to calculate the A(K,K') or B(K,K') matrix\n");

   fstates_t gs = hunds_gs(n, l);
   int J2min, J2max, ns=0;
   int K,Kp,Q,Qp,i,j;
   sMat<double> Akk,Bkk,ckk,Qmat,QLmat,QSmat;

   bool iA,ic;
   complexdouble Ykq; double Tj,Tj2;

   // Determines number of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2; for (int J2=J2min; J2<=J2max; J2+=2) ns+=J2+1; 
   int ins=0; std::vector<int> J2(ns), Jz2(ns), irm(ns);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2) {
        J2[ins]=iJ2; Jz2[ins]=imJ2; irm[ins]=(iJ2-J2min)/2; ins++; }

   Qmat.zero(ns,ns); Qq.clear(); for(i=0; i<6; i++) Qq.push_back(Qmat);

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

         // The matrices above already contain the dependence on |vULSJ> - we now use the W-E theorem to add the Jz, Q, dependence
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

}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the thermal expectation value of the FT of the magnetisation density -2Q in Bohr magnetons
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
      void mqcalc(ComplexVector &Mq,      // Output expectation values -2[<Q>_{x} <Q>_y <Q>_z]
                  double &th, double &ph, // Input polar and azimuth angles theta and phi
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est)     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
{
   int i,q,n=1,Hsz=est.Cols()-1; orbital l;
   n = (int)est[0][0].real(); i = (int)est[0][0].imag(); 
   switch(i) { 
      case 0: l=S; case 1: l=P; break; case 2: l=D; break; case 3: l=F; break; 
      default: std::cerr << "Error - only S,P,D,F-electrons supported.\n"; exit(0); }
   std::vector<double> E,Jvec(6,0.); Jvec[0]=th; Jvec[1]=ph; Jvec[2]=J0; Jvec[3]=J2; Jvec[4]=J4; Jvec[5]=J6;
   std::vector< sMat<double> > Qp, Qm; 
   std::vector< std::vector< sMat<double> > > Qmat; for(i=0; i<3; i++) Qmat.push_back(Qp);
   complexdouble *zQmat, *zt, zme, zalpha, zbeta; zalpha.r=-1; zalpha.i=0; zbeta.r=0; zbeta.i=0;
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

   for(q=0; q<3; q++)
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
//       printf ("%i zme=%g %+g i  Ei=%6.3f ni=%6.3f \n",i,zme.r,zme.i,est[0][i].real(),est[0][i].imag());
         zMqr += (-2.)*zme.r*est[0][i].imag(); zMqi += (-2.)*zme.i*est[0][i].imag(); if(q==0) Z += est[0][i].imag();
      }
      free(zQmat); free(zt); Mq[q+1] = complex<double> (zMqr, zMqi)/Z;
   }
// printf("MQ=(%g %+g i, %g %+g i,%g %+g i)\n",real(Mq(1)),imag(Mq(1)),real(Mq(2)),imag(Mq(2)),real(Mq(3)),imag(Mq(3)));
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the transition matrix using the scattering operator of Balcar and Lovesey
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
         int                              // Returns total number of transitions
             dmq1(int &tn,                // Input transition number |tn|. If tn>0 omit printout. If tn<0 print info.
                  double &th,             // Input zenith angle (with the z==b axis) in radians.
                  double &ph,             // Input azimuth angle (with the x==a axis, to projection in x-y plane).
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est,     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
                  double &T,              // Input temperature (K)
                  ComplexVector & mq1,    // input mq1(1)= ninit + i pinit   
                                          // Output transition vector, mq1=<-|M(Q)|+> sqrt(n- - n+) in units of MU_B
                  double & maxE)          // input maxE maximal transition energy
/* 
     Note on Qalpha (Qa or Qb)
        Kartesian components of the scattering operator Qalpha, alpha=1,2,3=a,b,c
        according to Lovesey Neutron Scattering equation 6.87b 
        scattering operator is given in  spherical coordinates Q-1,Q0,Q+1 (introduced
        as described above on input of th and ph) these are related to Qa,Qb,Qc by
        Q1=Qbx(axb)=Qy= i/sqrt(2)(Q+1 + Q-1) 
        Q2=Qb      =Qz= Q0                   
        Q3=Qaxb    =Qx=-1/sqrt(2)(Q+1 - Q-1)
                   
       
        the orbital and spin contributions 
        according to Lovesey Neutron Scattering equations 11.55 and 11.71 (the spin part 11.71 has to be
        divided by 2), i.e.
        <-|QSa,b,c|+>=
          =<-|sum_i exp(i k ri) s_(a,b,c)|+> /2                   as defined by 11.71 / 2
				   
        <-|QLa,b,c|+>=
          =<-|sum_i exp(i k ri) (-(k x grad_i)_(a,b,c)/|k|)|+>     as defined by 11.54 /(-|k|)

        mq1=<-|M(Q)|+> sqrt(n- - n+)
        <-|M(Q)|+>=-2<-|Q|+>=-2<-|2 QS + QL|+>
*/
{
   // check if printout should be done and make tn positive
   int pr=0; if (tn<0) { pr=1; tn*=-1; }
   double ninit=mq1[1].real();
   double pinit=mq1[1].imag();

   int i,iJ,q,n=1,Hsz=est.Cols()-1; orbital l;//=D; find_nl_from_dim(Hsz,*&n,*&l,(complexdouble*)&est[1][1]);
   n = (int)est[0][0].real(); i = (int)est[0][0].imag(); l=(orbital)i; //l = (i==2) ? D : F;
   std::vector<double> E,Jvec(6,0.); Jvec[0]=th; Jvec[1]=ph; Jvec[2]=J0; Jvec[3]=J2; Jvec[4]=J4; Jvec[5]=J6;
   std::vector< sMat<double> > Qp, Qm; 
   std::vector< std::vector< sMat<double> > > Qmat; for(i=0; i<3; i++) Qmat.push_back(Qp);
   complexdouble *zQmat, *zt, zalpha, zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0;
   std::vector<complexdouble> zij(7,zbeta), zji(7,zbeta);
   double Z=0., therm;
   char trans = 'U'; int incx=1;

   // Calculates the scattering operator, Q.
   icf_loveseyQq(Qm,-1,n,l,Jvec); icf_loveseyQq(Qp,1,n,l,Jvec);
   for(i=0; i<6; i++)  
   {
      Qmat[0].push_back( (Qp[i]-Qm[i]) * (-1./sqrt(2.)) );                     // Qx = -1/sqrt(2) * (Q_{+1} - Q_{-1})
      if(i%2==0) Qmat[1].push_back( (Qp[i+1]+Qm[i+1]) * (-1./sqrt(2.)) );      // real(Qy) = i^2/sqrt(2) * imag(Q_{+1}+Q_{-1})
      else       Qmat[1].push_back( (Qp[i-1]+Qm[i-1]) *  (1./sqrt(2.)) );      // imag(Qy) = i/sqrt(2) * real(Q_{+1}+Q_{-1})
   }
   icf_loveseyQq(Qmat[2],0,n,l,Jvec);

   for(i=0; i<Hsz; i++) { 
      therm = exp(-(est[0][i+1].real()-est[0][1].real())/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }

   int a,j=0,k=0; for(i=0; i<Hsz; ++i) { for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   ++i;++j; // because in est i and j start from 1...Hsz
 
   for(q=0; q<3; q++)
   {
      zQmat = zmat2f(Qmat[q][2],Qmat[q][3]);    // Spin part
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[j][1], &incx, &zbeta, zt, &incx);
      #ifdef _G77
      F77NAME(zdotc)(&zij[2*q+1], &Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx) ;
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[i][1], &incx, &zbeta, zt, &incx);
      F77NAME(zdotc)(&zji[2*q+1], &Hsz, (complexdouble*)&est[j][1], &incx, zt, &incx) ;
      #else
      zij[2*q+1] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx) ;
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[i][1], &incx, &zbeta, zt, &incx);
      zji[2*q+1] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[j][1], &incx, zt, &incx) ;
      #endif
      free(zQmat); free(zt);

      zQmat = zmat2f(Qmat[q][4],Qmat[q][5]);    // orbital part
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
      if(i==j)                                  // subtract thermal expectation value from zij=zii
      {                                         // MR120120 ... reintroduced
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
      if(fabs(zij[iJ].i+zji[iJ].i)>SMALL) { std::cerr << "ERROR module ic1ion - dmq1: <i|Qalpha|j>not hermitian\n"; exit(EXIT_FAILURE); }

   complex<double> im(0,1);
   ComplexVector iQalphaj(1,6);
   
   for(a=1; a<=6; a++) { iQalphaj(a) = complex<double> (zij[a].r,zij[a].i); if(a%2==1) { iQalphaj(a)*=0.5; } } // divide spin part by 2
   mq1 = 0;
   for(a=1; a<=3; a++) mq1(a) = -2.0*( 2.0*iQalphaj(-1+2*a) + iQalphaj(2*a) );

   double delta;
   delta = est[0][j].real()-est[0][i].real();
   if(delta<-0.000001) { std::cerr << "ERROR module ic1ion - dmq1: energy gain delta gets negative\n"; exit(EXIT_FAILURE); }

   if(j==i) delta = -SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc

   // do some printout if wishes and set correct occupation factor
   if (delta>SMALL)
   {
      therm = exp(-(est[0][i].real()-est[0][1].real())/(KB*T)) - exp(-(est[0][j].real()-est[0][1].real())/(KB*T));
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV",i,j,delta);
         printf(" |<%i|MQa|%i>|^2=%6.3f |<%i|MQb|%i>|^2=%6.3f |<%i|MQc|%i>|^2=%6.3f",i,j,abs(mq1(1))*abs(mq1(1)),i,j,abs(mq1(2))*abs(mq1(2)),i,j,abs(mq1(3))*abs(mq1(3)));
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
         printf(" |<%i|MQa|%i>|^2=%6.3f |<%i|MQb|%i>|^2=%6.3f |<%i|MQc|%i>|^2=%6.3f",i,j,abs(mq1(1))*abs(mq1(1)),i,j,abs(mq1(2))*abs(mq1(2)),i,j,abs(mq1(3))*abs(mq1(3)));
         printf(" n%i=%6.3f\n",i,therm/Z);
      }
   }
   mq1 *= sqrt(therm / Z);

    // determine number of thermally reachable states
   if (ninit>Hsz)ninit=Hsz;
   //if (pinit<SMALL)pinit=SMALL;
   double zsum=0,zi,x;
   int noft=0; 
   for(i=0; (i<ninit)&((((x=(est[0][i+1].real()-est[0][1].real())/(KB*fabs(T)))<200)? zi=exp(-x):zi=0)>=(pinit*zsum)); ++i)
   {
      noft += Hsz-i; 
      zsum += zi;
   }
// removed MR  6.9.2011 to allow for mcdisp options -ninit -pinit
// int noft=0;for(i=0; (i<Hsz)&((exp(-(est[0][i+1].real()-est[0][1].real())/(KB*T)))>SMALL); ++i) noft += Hsz-i-1;
   return noft;
// return Hsz*(Hsz-1)/2;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the coefficient of the spin density operator M^S_q(r), after Balcar, J.Phys.C. v8, p1581, 1975, eqn 10.
// --------------------------------------------------------------------------------------------------------------- //
bool icf_balcarMSq(sMat<double> &MSq, int q, int K, int Q, int n, orbital l)  
{
   if (l>3||l<0) { std::cerr << "Sorry only s-, p-, d- and f-electrons supported at present\n"; exit(0); }

   std::string errormsg("icf_balcarMSq(): Unable to calculate the c(K,K') matrix\n");
   fstates_t gs = hunds_gs(n, l);
   int J2min, J2max, ns=0;

   // Determines number and angular momentum quantum numbers of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2; for (int J2=J2min; J2<=J2max; J2+=2) ns+=J2+1;
   int ins=0; std::vector<int> J2(ns), Jz2(ns), irm(ns);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2) {
        J2[ins]=iJ2; Jz2[ins]=imJ2; irm[ins]=(iJ2-J2min)/2; ins++; }

   sMat<double> ckk,Qmat; MSq.zero(ns,ns);

   // Eqn 10 of Balcar 1975:
   //
   //         S                     -2uB     ---  K     2    ---   [     1/2+M+q-J'+L'+S'     1/2
   // <vSLJM|M (r)|v'S'L'J'M'> =  ---------  >   Y (r) U (r) >     [ (-1)                (3/2)
   //         q                   sqrt(4PI)  ---  Q          ---   [
   //                                        K,Q             K',Q'
   //
   //              \/                                1/2  ( l K l )  { 1  K  K' }
   //              /\    [l,l,S,S',L,L',J,J',K,K',K']     ( 0 0 0 )  { S' L' J' }
   //                                                                { S  L  J  }
   //                                            _ _
   //                     ---     _   _          S+L                                                         ]
   //              \/   n >   (t{|t) (t|}t') (-1)    { S  1  S' }  { L  K  L' }    (  J K' J' ) ( K K'  1 )  ]
   //              /\     ---                        { s  Sb s  }  { l  Lb l  }    ( -M Q' M' ) ( Q Q' -q )  ]
   //                      t
   //
   // In this function we calculate only the sum over the square brackets, leaving the K,Q sum to chrgplt.
   // The first 3j symbol imposes the selection rule that K=even,K<=2l. The 9j symbol means (K-1)<=K'<=(K+1), whilst
   // the final 3j symbol means that Q+Q'=q, that is Q'=q-Q, and that -K'<=Q'<=K', so there is only ever one Q' term
   // in the sum.
   //
   // Since for icf1ion, we only calculated for the lowest LS-multiplet, S==S' and L==L'

   if(K%2!=0 || K>(2*abs(l)) || Q<-K || Q>K || q<-1 || q>1) return 1;

   // Simplifying, using the C(K,K') coeficients previously calculated in icf_loveseyCKK()  (Balcar+Lovesey Book, eqn 3.8.5)
   //
   //          S                    -2uB    ---  K     2    ---   [             M+q-J'      1/2  -K                           ]
   //  <vSLJM|M (r)|v'S'L'J'M'> = --------- >   Y (r) U (r) >     [ C(K,K') (-1)      (3[J])    i   (  J K' J' ) ( K K'  1 )  ]
   //          q                  sqrt(4PI) ---  Q          ---   [                                 ( -M Q' M' ) ( Q Q' -q )  ]
   //                                       K,Q             K',Q'
   //

   int Kp, Qp=q-Q;
   int i,j;
   bool ic;
   double Tj, Tj2, prefact = (K%4==2) ? 1 : -1; prefact*=sqrt(3./PI); // Prefactor is (-2/sqrt(4PI))*sqrt(3)i^-K. Elements in units of uB

   for(Kp=(K-1); Kp<=(K+1); Kp++)
   {
      if(Qp<-Kp || Qp>Kp) continue; Qmat.zero(ns,ns);
      Tj = prefact * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q);
      ic = icf_loveseyCKK(ckk,K,Kp,n,l); if(ic) ckk *= Tj; else { std::cerr << errormsg; return 1; }

      for(i=0; i<ns; i++) for(j=0; j<ns; j++)
      {
         Tj2 = threej(J2[i],2*Kp,J2[j],-Jz2[i],2*Qp,Jz2[j]); if(fabs(Tj2)<DBL_EPSILON) continue;
         Tj2 *= pow(-1.,q+(Jz2[i]-J2[j])/2.) * sqrt(J2[i]+1.);
         Qmat(i,j) = Tj2 * ckk(irm[i],irm[j]);
      }
      MSq += Qmat;
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the coefficient of the orbital magnetic density operator M^L_q(r), after Balcar 1975, eqn 12.
// --------------------------------------------------------------------------------------------------------------- //
bool icf_balcarMLq(sMat<double> &MLq, int q, int K, int Q, int n, orbital l)  
{
   if (l>3||l<0) { std::cerr << "Sorry only s-, p-, d- and f-electrons supported at present\n"; exit(0); }

   std::string errormsg("icf_balcarMLq(): Unable to calculate the a(K,K') matrix\n");
   fstates_t gs = hunds_gs(n, l);
   int J2min, J2max, ns=0;

   // Determines number and angular momentum quantum numbers of basis states
   J2min = abs(2*gs.L-gs.S2); J2max = 2*gs.L+gs.S2; for (int J2=J2min; J2<=J2max; J2+=2) ns+=J2+1;
   int ins=0; std::vector<int> J2(ns), Jz2(ns), irm(ns);
   for (int iJ2=J2min; iJ2<=J2max; iJ2+=2) for (int imJ2=-iJ2; imJ2<=iJ2; imJ2+=2) {
        J2[ins]=iJ2; Jz2[ins]=imJ2; irm[ins]=(iJ2-J2min)/2; ins++; }

   sMat<double> akk,Qmat; MLq.zero(ns,ns);

   // Eqn 12 of Balcar 1975:
   //                                                       infty
   //         L                     -2uB     ---  K     1  /     2     ---   [     q+M+L+L'+S
   // <vSLJM|M (r)|v'S'L'J'M'> =  ---------  >   Y (r) --- | dX U (X)  >     [ (-1)            d
   //         q                   sqrt(4PI)  ---  Q     r  /           ---   [                  SS'
   //                                        K,Q            r          K',Q'
   //
   //              \/                             1/2         1/2  ( l K l )  { l  K' l  }  { K' L' L  }
   //              /\    [l,l,l,L,L',J,J',K,K',K']    [l(l+1)]     ( 0 0 0 )  { K  l  1  }  { S  J  J' }
   //                                                                                      
   //                                            _
   //                     ---     _   _          L                                           ]
   //              \/   n >   (t{|t) (t|}t') (-1)  { L  K  L' }    (  J K' J' ) ( K K'  1 )  ]
   //              /\     ---                      { l  Lb l  }    ( -M Q' M' ) ( Q Q' -q )  ]
   //                      t
   //
   // In this function we calculate only the sum over the square brackets and include the factor -2/sqrt(4pi), leaving the K,Q sum to chrgplt.
   // The first 3j symbol imposes the selection rule that K=even,0<=K<=2l. The 9j symbol means (K-1)<=K'<=(K+1), whilst
   // the final 3j symbol means that Q+Q'=q, that is Q'=q-Q, and that -K'<=Q'<=K', so there is only ever one Q' term
   // in the sum.

   if(K%2!=0 || K>(2*l) || Q<-K || Q>K || q<-1 || q>1) return 1;
   int Kp,Qp=q-Q;

   // Noting that the e(K,K') coefficient is, from Balcar and Lovesey, eqn 3.6.4:
   //
   //                    i^K                            1/2          1/2      S+L+L'+J'  ( l K l ) { 1 K' K } { L  K' L' }
   //  e(K,K') = d     --------  [l,l,l,K,K',K',J',L,L']     [l(l+1)]     (-1)           ( 0 0 0 ) { l l  l } { J' S  J  }
   //             SS'  2 sqrt(3)
   //                                            _
   //                     ---     _   _          L
   //              \/   n >   (t{|t) (t|}t') (-1)  { L  K  L' }
   //              /\     ---                      { l  Lb l  }
   //                      t
   //
   // We can simplify the matrix element above, eqn 3.8.9. Note the 2/r rather than 1/r in the previous equation.
   //
   //          L                    -2uB    ---  K      infty
   //  <vSLJM|M (r)|v'S'L'J'M'> = --------- >   Y (r) 2  /     2    ---   [          -K       1/2    M-J'+q                           ]
   //          q                  sqrt(4PI) ---  Q   --- | dX U(X)  >     [ e(K,K') i  (3[J])    (-1)       (  J K' J' ) ( K  K'  1 ) ]
   //                                       K,Q       r  /          ---   [                                 ( -M Q' M' ) ( Q  Q' -q ) ]
   //                                                    r          K',Q'
   //
   // Finally, we use the relation 4.2.14 of Balcar and Lovesey,  e(K,K') = +/- 0.5(2K'+1) a(K,K')  where + is for K'=K+1, - for K'=K-1

   int i,j;
   bool iA;
   double Tj, Tj2, prefact = (K%4==2) ? 1 : -1; prefact*=sqrt(3./PI); // Prefactor is (-2/sqrt(4PI))*sqrt(3)i^-K. Elements in units of uB

   for(Kp=(K-1); Kp<=(K+1); Kp+=2)       // The 3j symbol in A(K,K') means that K==Kp gives zero... (NB. Does this apply to E(K,K') too?)
   {
      if(Qp<-Kp || Qp>Kp) continue; Qmat.zero(ns,ns);
      if(Kp==(K-1)) Tj = prefact * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q) * -(2*Kp+1);
      else          Tj = prefact * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q) *  (2*Kp+1);
      iA = icf_loveseyAKK(akk,K,Kp,n,l); if(iA) akk *= Tj; else { std::cerr << errormsg; return 1; }
      for(i=0; i<ns; i++) for(j=0; j<ns; j++)
      {
         Tj2 = threej(J2[i],2*Kp,J2[j],-Jz2[i],2*Qp,Jz2[j]); if(fabs(Tj2)<DBL_EPSILON) continue;
         Tj2 *= pow(-1.,q+(Jz2[i]-J2[j])/2.) * sqrt(J2[i]+1.);
         Qmat(i,j) = Tj2 * akk(irm[i],irm[j]);
      }
      MLq += Qmat;
   }

   return 1;
}

complexdouble * icf_balcarMq(int xyz, int K, int Q, int n, orbital l)
{
   int Hsz = getdim(n,l);
   sMat<double> retval_r(Hsz,Hsz),retval_i(Hsz,Hsz),qpp,qmp,qpm,qmm;

   if(xyz==1||xyz==2)
   {
      if(Q!=0)
      {
         icf_balcarMSq(qpp, 1,K, abs(Q),n,l); rmzeros(qpp);
         icf_balcarMSq(qmp, 1,K,-abs(Q),n,l); rmzeros(qmp);
         icf_balcarMSq(qpm,-1,K, abs(Q),n,l); rmzeros(qpm);
         icf_balcarMSq(qmm,-1,K,-abs(Q),n,l); rmzeros(qmm);
	 // sum to coeff of Zlm (neglecting a 1/sqrt(2) factor) 
         if(Q<0) { if(Q%2==0) { qmp -= qpp; qmm -= qpm; } else { qmp += qpp; qmm += qpm; } }
         else    { if(Q%2==0) { qmp += qpp; qmm += qpm; } else { qmp -= qpp; qmm -= qpm; } }
         // add spherical components and multiply by addition factor 1/sqrt(2) which was neglected in the line above
         if(xyz==1) { if(Q<0) retval_i = (qmm-qmp)/(-2.);    else retval_r = (qmm-qmp)/2.; }
         if(xyz==2) { if(Q<0) retval_r = (qmm+qmp)/2.;       else retval_i = (qmm+qmp)/2.; } // changed MR 25.5.2010 // Q<0 signs changed MR 28.5.2010
      }
      else
      {
         icf_balcarMSq(qpp, 1,K,0,n,l); rmzeros(qpp);
         icf_balcarMSq(qpm,-1,K,0,n,l); rmzeros(qpm);
         if(xyz==1) {  retval_r = (qpp-qpm)/(-sqrt(2.)); }
         if(xyz==2) {  retval_i = (qpp+qpm)/sqrt(2.);    } // changed MR 25.5.2010
      }
   }
   else if(xyz==3)
   {
      if(Q!=0)
      {
         icf_balcarMSq(qpp,0,K, abs(Q),n,l); rmzeros(qpp);
         icf_balcarMSq(qmp,0,K,-abs(Q),n,l); rmzeros(qmp);
         if(Q<0) {  if(Q%2==0) qmp -= qpp; else qmp += qpp;  retval_i = qmp/(-sqrt(2.)); }
         else    {  if(Q%2==0) qmp += qpp; else qmp -= qpp;  retval_r = qmp/sqrt(2.);    } // changed by MR 28.5.2010
      }
      else
      {
         icf_balcarMSq(retval_r,0,K,0,n,l); rmzeros(retval_r);
      }
   }
   else if(xyz==-1||xyz==-2)
   {
      if(Q!=0)
      {
         icf_balcarMLq(qpp, 1,K, abs(Q),n,l); rmzeros(qpp);
         icf_balcarMLq(qmp, 1,K,-abs(Q),n,l); rmzeros(qmp);
         icf_balcarMLq(qpm,-1,K, abs(Q),n,l); rmzeros(qpm);
         icf_balcarMLq(qmm,-1,K,-abs(Q),n,l); rmzeros(qmm);
	 // sum to coeff of Zlm (neglecting a 1/sqrt(2) factor) 
         if(Q<0) { if(Q%2==0) { qmp -= qpp; qmm -= qpm; } else { qmp += qpp; qmm += qpm; } }
         else    { if(Q%2==0) { qmp += qpp; qmm += qpm; } else { qmp -= qpp; qmm -= qpm; } }
         // add spherical components and multiply by addition factor 1/sqrt(2)which was neglected in the line above
         if(xyz==-1) { if(Q<0) retval_i = (qmm-qmp)/(-2.);    else retval_r = (qmm-qmp)/2.; }
         if(xyz==-2) { if(Q<0) retval_r = (qmm+qmp)/2.;       else retval_i = (qmm+qmp)/2.; } // changed MR 25.5.2010 // Q<0 signs changed MR 28.5.2010
      }
      else
      {
         icf_balcarMLq(qpp, 1,K,0,n,l); rmzeros(qpp);
         icf_balcarMLq(qpm,-1,K,0,n,l); rmzeros(qpm);
         if(xyz==-1) { retval_r = (qpp-qpm)/(-sqrt(2.)); }
         if(xyz==-2) { retval_i = (qpp+qpm)/sqrt(2.);    } // changed MR 25.5.2010
      }
   }
   else if(xyz==-3)
   {
      if(Q!=0)
      {
         icf_balcarMLq(qpp,0,K, abs(Q),n,l); rmzeros(qpp);
         icf_balcarMLq(qmp,0,K,-abs(Q),n,l); rmzeros(qmp);
         if(Q<0) { if(Q%2==0) qmp -= qpp; else qmp += qpp; retval_i = qmp/(-sqrt(2.)); }
         else    { if(Q%2==0) qmp += qpp; else qmp -= qpp; retval_r = qmp/sqrt(2.);    } // changed by MR 28.5.2010
      }
      else
      {
         icf_balcarMLq(retval_r,0,K,0,n,l); rmzeros(retval_r);
      }
   }
   return zmat2f(retval_r,retval_i);
}

//--------------------------------------------------------------------------------------------------------------
void icf_spindensityexpJ(icpars &pars, complexdouble *zV, double *vE, int xyz, double *T, Vector &J)
{
   int k[] = {0,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int Hsz=icf_getdim(pars), incx=1;

   char xyzstr[] = "xyz";
   if(xyz>0) { std::cout << "Calculating the expectation values of the spin density operator S" << xyzstr[xyz-1] << "\n"; }
   else      { std::cout << "Calculating the expectation values of the orbital moment density operator L" << xyzstr[-xyz-1] << "\n"; }

   std::vector<double> E, ex((J.Hi()>6?J.Hi():6)+2,0.), me, eb;
   int iJ, ind_j, Esz;

   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';

   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<Hsz; Esz++) { E.push_back(vE[Esz]-vE[0]); if(exp(-E[Esz]/(KB**T))<DBL_EPSILON || vE[Esz+1]==0) break; }
   if (*T<0)
   {
      Esz = (int)(-*T);
      printf ("Temperature T=%g<0: please choose probability distribution for the -T=%i lowest energy states by hand\n",*T,(int)(-*T));
                        printf("Number   Excitation Energy\n");
      for (ind_j=0; ind_j<Esz; ++ind_j) printf ("%i    %4.4g meV\n",ind_j+1,E[ind_j]);
   }  // MR 10.9.2010

   for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++)
   {
      int ind = ii+Hsz*jj;
      if(fabs(zV[ind].r*zV[ind].r+zV[ind].i*zV[ind].i)<DBL_EPSILON*100000) 
      {
         zV[ind].r=0.; zV[ind].i=0.;  
      } 
   }

   // For first run calculate also the partition function
   me.assign(Esz,0.); eb.assign(Esz,0.); double Z=0.; complexdouble *zt=0, zme;
   complexdouble *zJmat;// = (complexdouble*)malloc(Hsz*Hsz*sizeof(complexdouble));
   zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
   sMat<double> zeros(Hsz,Hsz),mat;

   for(iJ=J.Lo(); iJ<=J.Hi(); iJ++)
   {
      me.assign(Esz,0.); J[iJ]=0.; ex[iJ]=0;
      // Using the above reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
      zJmat = icf_balcarMq(xyz,k[iJ],q[iJ],pars.n,pars.l); // minus sign stands for orbital density coeff
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));

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
            {
               char instr[MAXNOFCHARINLINE];
               printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
               if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { fprintf(stderr,"Error reading input\n"); exit(1); }
               eb[ind_j]=strtod(instr,NULL);
            }
            else
            {
               eb[ind_j] = exp(-E[ind_j]/(KB**T));
            }
            J[iJ] += me[ind_j]*eb[ind_j];
            Z += eb[ind_j];
//MRend 10.9.2010
         }
         else  // Rest of the runs only calculate the new matrix elements
            J[iJ]+=me[ind_j]*eb[ind_j];
      }
      J[iJ]/=Z;
      if(fabs(J[iJ])<DBL_EPSILON) J[iJ]=0.;
      free(zJmat); free(zt);
   }
  
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of spindensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
void sdod_Icalc(Vector &J,           // Output single ion moments==(expectation values) Zlm R^2(r) at given T, H_eff
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

   int K[] = {-1,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Q[] = {-1,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {-1,0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
   int i;

   // Converts the Jij parameters if necessary
   std::vector<double> vgjmbH(gjmbH.Hi()+1,0.);
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc();
      for(i=gjmbH.Lo(); i<=gjmbH.Hi(); i++) vgjmbH[i] = -gjmbH[i]*pars.jijconv[i]; }
   else
   #endif
      for(i=gjmbH.Lo(); i<=gjmbH.Hi(); i++) vgjmbH[i] = -gjmbH[i];

   // Calculates the IC Hamiltonian matrix
   int nfact = (int)ceil(sqrt(J.Hi()-J.Lo()+2));
   int k,q,Hsz=icf_getdim(pars),esz=Hsz*nfact;
   int ix, iy, incx=1;
   complexdouble *H=0;
   bool Hicnotcalc = false;
   std::vector<double> parval; parval.reserve(35);
   parval.push_back(pars.xi); for(k=2; k<=(2*pars.l); k+=2) for(q=-k; q<=k; q++) parval.push_back(pars.B(k,q));
   if(parval.size()%2==1) parval.push_back(0.);

// if((est.Cols()!=(esz+1) || est.Rows()!=(esz+1))) Hicnotcalc = true;
   if(real(est[0][0])==(double)pars.n && imag(est[0][0])==(double)pars.l)  // Hic previously calculated
   {
      for(i=0; i<(int)(parval.size()/2); i++) if(real(est[0][i+2])!=parval[2*i] || imag(est[0][i+2])!=parval[2*i+1]) { 
         Hicnotcalc = true; break; }
   }
   else Hicnotcalc = true;

   int oldJhi = J.Hi();
   if((est.Rhi()!=esz||est.Chi()!=esz) && !Hicnotcalc) // Called from spins: need all multipolar operators, not just those used in MF loop.
   {
      nfact = (int)real(est[0][1]); esz = Hsz*nfact;
      oldJhi = (int)imag(est[0][1]);
   }

   if(Hicnotcalc || pars.l==S)
   {
      sMat<double> Hcfi, Hcf = icf_hmltn(Hcfi, pars); Hcf/=MEV2CM; Hcfi/=MEV2CM; H = zmat2f(Hcf,Hcfi);
      if(est.Rhi()!=esz||est.Chi()!=esz) {
         std::cerr << "ERROR module icf1ion - Icalc: Hsz recalculation does not agree with eigenstates matrix dimension\n"; exit(EXIT_FAILURE); }
      else if(esz>(int)parval.size()/2+1)
      {
         est[0][0] = complex<double> (pars.n,pars.l);
         for(i=0; i<(int)(parval.size()/2); i++) est[0][i+2] = complex<double> (parval[2*i],parval[2*i+1]);
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
      complex<double> a(vgjmbH[ind],0.);
      if(ind<=oldJhi)
      {
         iy = ind/nfact; ix = ind-iy*nfact;
         for(i=1; i<=Hsz; i++) F77NAME(zaxpy)(&Hsz,(complexdouble*)&a,(complexdouble*)&est[i+ix*Hsz][1+iy*Hsz],&incx,&H[(i-1)*Hsz],&incx);
      }
      else
      {
         sMat<double> Hcf, Hcfi(Hsz,Hsz); complexdouble *zM;
         if(ind<=6) {       // Calculates Sx,Lx etc and copies them to est too. 
            Hcf = icf_mumat(pars.n, ind-1, pars.l); if(ind==3 || ind==4) zM = zmat2f(Hcfi,Hcf); else zM = zmat2f(Hcf,Hcfi); }
         else       {       // Calculates multipolar operator matrices and copies them to est. 
            Hcf = icf_ukq(pars.n,K[ind],Q[ind],pars.l); if(im[ind]==1)   zM = zmat2f(Hcfi,Hcf); else zM = zmat2f(Hcf,Hcfi); }
         for(i=1; i<=Hsz; i++) F77NAME(zaxpy)(&Hsz,(complexdouble*)&a,&zM[(i-1)*Hsz],&incx,&H[(i-1)*Hsz],&incx); free(zM);
      }
   }

   // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
   double *vE = new double[Hsz]; complexdouble *zV = new complexdouble[Hsz*Hsz];
   int info = ic_diag(Hsz,H,zV,vE); free(H);
   if(info!=0) { std::cerr << "icf1ion - Error diagonalising, info==" << info << "\n"; delete[]vE; vE=0; delete[]zV; zV=0; exit(EXIT_FAILURE); }

   icf_spindensityexpJ(pars,zV,vE,xyz,T,J); delete[]vE; delete[]zV;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of chargedensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
void chargedensity_coeff(
                      Vector &mom,         // Output single ion moments == expectation values of
                                           //    of Zlm R^2(r) at a given temperature T and  effective field H
                      double *T,           // Input scalar temperature
                      Vector &Hxc,         // Input vector of exchange fields (meV) 
                      Vector &Hext,        // Input vector of external field (T) 
 /* Not Used */       double *g_J,         // Input Lande g-factor
 /* Not Used */       Vector &ABC,         // Input vector of parameters from single ion property file
                      char **sipffilename, // Single ion properties filename
                      ComplexMatrix &est)  // Input/output eigenstate matrix (initialized in parstorage)
{
   Vector moments(1,51); 
   double lnZ, U;
   Vector Hxce(1,51); 
   Hxce=0;
   for(int i=1; i<=Hxc.Hi(); ++i) { Hxce(i)=Hxc(i); }
   Icalc(moments,T,Hxce,Hext,g_J,ABC,sipffilename,&lnZ,&U,est);

   // Parses the input file for parameters
   icpars pars; 
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

// Definitions in Icalc()
// //          0 1 2 3 4 5 6  7  8 91011 12 13 1415161718 19 20 21 222324252627 28 29 30 31 32333435363738 39 40 41 42 43 4445464748495051
// int K[] = {-1,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
// int Q[] = {-1,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
  
// a(0, 0) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// Indices for spindensity
//             0 not used
//             0 1  2  3 4 5 6  7  8  9 101112131415 16 17 18 19 20 2122232425262728 
// int k[] = {-1,0, 2, 2,2,2,2, 4, 4, 4, 4,4,4,4,4,4, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
// int q[] = {-1,0,-2,-1,0,1,2,-4,-3,-2,-1,0,1,2,3,4,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   mom(1) =  pars.n / sqrt(4.0 * 3.1415); // nofelectrons 
   for(int i=2; i<=6; ++i) mom(i) = moments(5+i) *sqrt((2.0*2+1)/8/PI); mom(4)  *= sqrt(2);
   for(int i=7; i<=15;++i) mom(i) = moments(12+i)*sqrt((2.0*4+1)/8/PI); mom(11) *= sqrt(2);
   for(int i=16;i<=28;++i) mom(i) = moments(23+i)*sqrt((2.0*6+1)/8/PI); mom(22) *= sqrt(2);

// {for(l=2;l<=6;l+=2){for(m=-l;m<=l;++m){if(m!=0){a(l,m)*=sqrt((2.0*l+1)/8/PI);}else{a(l,m)*=sqrt((2.0*l+1)/4/PI);}}}
// } // in case of module ic1ion we just take the prefactors of the Zlm ... ??? what should we take here ???
// MR 23.8.2011: if Tkq are define as in our review then the above should be right
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements of chargedensity coefficient operator
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
int dchargedensity_coeff1(int &tn,        // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexVector & dc1,// Output m1 vector (1,3)
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ 
   ComplexVector u1(1,51);int n,nd;
   u1(1) = dc1(1);
   Vector Hxce(1,51);
   Hxce = 0;
   for(int i=1; i<=Hxc.Hi(); ++i) { Hxce(i)=Hxc(i); }
   int nt = du1calc(tn,T,Hxce,Hext,g_J,ABC,sipffilename,u1,delta,n,nd,est);
   dc1(1)=0;
   for(int i=2; i<=6; ++i) dc1(i) = u1(5+i) *sqrt((2.0*2+1)/8/PI); dc1(4) *=sqrt(2);
   for(int i=7; i<=15;++i) dc1(i) = u1(12+i)*sqrt((2.0*4+1)/8/PI); dc1(11)*=sqrt(2);
   for(int i=16;i<=28;++i) dc1(i) = u1(23+i)*sqrt((2.0*6+1)/8/PI); dc1(22)*=sqrt(2);
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of spindensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
void spindensity_coeff(Vector &J,          // Output single ion moments =expectation values of
                                           //    of Zlm R^2(r) at a given temperature T and  effective field H
                      int & xyz,           // direction 1,2,3 = x,y,z
                      double *T,           // Input scalar temperature
                      Vector &Hxc,         // Input vector of exchange fields (meV) 
                      Vector &Hext,        // Input vector of external field (T) 
 /* Not Used */       double * /*g_J*/,    // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,    // Input vector of parameters from single ion property file
                      char **sipffilename, // Single ion properties filename
                      ComplexMatrix &est)  // Input/output eigenstate matrix (initialized in parstorage)
{  // sum exchange field and external field
   Vector gjmbH(1,Hxc.Hi());
   gjmbH=Hxc;
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   double dum; dum=gjmbH(2);gjmbH(2)=gjmbH(4);gjmbH(4)=gjmbH(5);gjmbH(5)=gjmbH(3);gjmbH(3)=dum;
   // --------------------------------------------------------------------

   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(Hext(1))>DBL_EPSILON || fabs(Hext(2))>DBL_EPSILON || fabs(Hext(3))>DBL_EPSILON)
   {
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH(2)+=MUB*Hext(1); gjmbH(1)+=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH(4)+=MUB*Hext(2); gjmbH(3)+=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH(6)+=MUB*Hext(3); gjmbH(5)+=GS*MUB*Hext(3); }
   }
   sdod_Icalc(J,xyz,T,gjmbH,sipffilename,est);
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
void orbmomdensity_coeff(Vector &J,        // Output single ion moments =expectation values of
                                           //    of Zlm R^2(r) at a given temperature T and  effective field H
                      int & xyz,           // direction 1,2,3 = x,y,z
                      double *T,           // Input scalar temperature
                      Vector &Hxc,         // Input vector of exchange fields (meV) 
                      Vector &Hext,        // Input vector of external field (T) 
 /* Not Used */       double * /*g_J*/,    // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,    // Input vector of parameters from single ion property file
                      char **sipffilename, // Single ion properties filename
                      ComplexMatrix &est)  // Input/output eigenstate matrix (initialized in parstorage)
{  // sum exchange field and external field
   Vector gjmbH(1,Hxc.Hi());
   gjmbH=Hxc;
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   double dum; dum=gjmbH(2);gjmbH(2)=gjmbH(4);gjmbH(4)=gjmbH(5);gjmbH(5)=gjmbH(3);gjmbH(3)=dum;
   // --------------------------------------------------------------------
   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(Hext(1))>DBL_EPSILON || fabs(Hext(2))>DBL_EPSILON || fabs(Hext(3))>DBL_EPSILON)
   {
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH(2)+=MUB*Hext(1); gjmbH(1)+=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH(4)+=MUB*Hext(2); gjmbH(3)+=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH(6)+=MUB*Hext(3); gjmbH(5)+=GS*MUB*Hext(3); }
   }
   sdod_Icalc(J,-xyz,T,gjmbH,sipffilename,est);
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the matrix elements of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
int      sdod_du1calc(int xyz,            // Indicating which of x,y,z direction to calculate; >0==spin <0==orbit
                      int &tn,            // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &gjmbH,      // Input vector of exchange fields + external fields (meV)
                      char **sipffilename,// Single ion properties filename
                      ComplexVector &u1,  // Output Slm1/Llm1 vector (1,49)
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
{
   int i,j,k;

   int K[] = {0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Q[] = {0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
// int im[]= {0, 1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};

   int sz = gjmbH.Hi();
   sMat<double> zeroes(est.Rows()-1,est.Cols()-1), op;
   complexdouble *zJmat=0, *zt=0, zme; zme.r=0; zme.i=0.; 
   std::vector<complexdouble> zij(sz,zme);//, zji(6,zme);
   std::vector<double> u(sz+1),iu(sz+1);
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';
   double Z=0., therm;

   // check if printout should be done and make tn positive
   int pr=0; if (tn<0) { pr=1; tn*=-1; }
   double ninit=u1[1].real();
   double pinit=u1[1].imag();

   // Copies the already calculated energy levels / wavefunctions from *est
   if(est.Rows()!=est.Cols()) { std::cerr << "sdod_du1calc(): Input rows and columns of eigenstates matrix don't match.\n"; return 0; }
   int Hsz = est.Rows()-1, iJ, incx = 1;
   j=0; k=0; for(i=0; i<Hsz; ++i) { for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   if(est[0][j+1].real()-est[0][i+1].real()<delta)
   {
      double *en = new double[Hsz]; for(k=0; k<Hsz; k++) en[k] = est[0][k+1].real();

      // Parses the input file for parameters
      icpars pars; const char *filename = sipffilename[0];
      ic_parseinput(filename,pars);

      // Calculates the transition matrix elements:
      //    u1 = <i|Ja|j> * sqrt[(exp(-Ei/kT)-exp(-Ej/kT)) / Z ]   if delta > small
      //    u1 = <i|Ja-<Ja>|j> * sqrt[(exp(-Ei/kT)) / kTZ ]             if delta < small (quasielastic scattering)
      for(iJ=0; iJ<sz; iJ++)
      {
//       if(iJ<6) op = icf_mumat(pars.n, iJ, pars.l); else op = icf_ukq(pars.n,K[iJ],Q[iJ],pars.l); 
//       if(im[iJ]==1) zJmat=zmat2f(zeroes,op); else zJmat=zmat2f(op,zeroes);
         zJmat = icf_balcarMq(xyz,K[iJ],Q[iJ],pars.n,pars.l); // minus sign stands for orbital density coeff
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, (complexdouble*)&est[j+1][1], &incx, &zbeta, zt, &incx);
         #ifdef _G77 
         F77NAME(zdotc)(&zij[iJ], &Hsz, (complexdouble*)&est[i+1][1], &incx, zt, &incx);
         #else
         zij[iJ] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[i+1][1], &incx, zt, &incx);
         #endif

         if(i==j && T>0) // subtract thermal expectation value from zij=zii
         {
            std::vector<double> eb, E; E.reserve(Hsz);
            int Esz; for(Esz=0; Esz<Hsz; Esz++) { E.push_back(en[Esz]-en[0]); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || en[Esz+1]==0) break; }
            eb.assign(Esz,0.); double Jj=0.;

            for(int ind_j=0; ind_j<Esz; ind_j++)
            {  // Calculates the matrix elements <Vi|J.H|Vi>
               F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, (complexdouble*)&est[ind_j+1][1], &incx, &zbeta, zt, &incx);
               #ifdef _G77 
               F77NAME(zdotc)(&zme, &Hsz, (complexdouble*)&est[ind_j+1][1], &incx, zt, &incx);
               #else
               zme = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[ind_j+1][1], &incx, zt, &incx);
               #endif
               // For first run calculate also the partition function and internal energy
               if(iJ==0)
               {
//MR 10.9.2010
                  if (T<0) 
                  {
                     char instr[MAXNOFCHARINLINE];
                     printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
                     if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { fprintf(stderr,"Error reading input\n"); exit(1); }
                     eb[ind_j]=strtod(instr,NULL); 
                  }
//MRend 10.9.2010
                  else
                     eb[ind_j] = exp(-E[ind_j]/(KB*T));  
		  Z+=eb[ind_j];
	       }
               Jj+=zme.r*eb[ind_j];
            }
            zij[iJ].r -= Jj/Z; 
         }
         free(zJmat); free(zt);
      }

      if (T<0) T=-T;

      // Calculates the vector u , iu
      for(iJ=0; iJ<sz; iJ++)
      {  
         u[iJ+1] = zij[iJ].r;
         iu[iJ+1] =zij[iJ].i;
      }

      delta = en[j]-en[i];
      if(delta<-0.000001) {
         std::cerr << "ERROR module ic1ion - du1calc: energy gain delta gets negative\n"; exit(EXIT_FAILURE); }
      if(j==i)delta=-SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc
   
      // Calculates the partition function
      Z=0.; for(iJ=0; iJ<Hsz; iJ++) { therm = exp(-(en[iJ]-en[0])/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }
   
      // do some printout if wishes and set correct occupation factor
      if (delta>SMALL)
      {
         therm = exp(-(en[i]-en[0])/(KB*T)) - exp(-(en[j]-en[0])/(KB*T));
         if(pr==1)
         {
            printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
            printf(" |<%i|Ja|%i>|^2=%6.3f\n |<%i|Jb|%i>|^2=%6.3f\n |<%i|Jc|%i>|^2=%6.3f\n",i+1,j+1,u[1]*u[1]+iu[1]*iu[1],i+1,j+1,u[2]*u[2]+iu[2]*iu[2],i+1,j+1,u[3]*u[3]+iu[3]*iu[3]);
            printf(" |<%i|Jd|%i>|^2=%6.3f\n |<%i|Je|%i>|^2=%6.3f\n |<%i|Jf|%i>|^2=%6.3f\n",i+1,j+1,u[4]*u[4]+iu[4]*iu[4],i+1,j+1,u[5]*u[5]+iu[5]*iu[5],i+1,j+1,u[6]*u[6]+iu[6]*iu[6]);
            printf(" n%i-n%i=%6.3f\n",i,j,therm / Z);
         }
      }
      else
      {
         therm = exp(-(en[i]-en[0])/(KB*T))/(KB*T);    // quasielastic scattering has not wi-wj but wj*epsilon/kT
         if(pr==1)
         {
            printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
            printf(" |<%i|Ja-<Ja>|%i>|^2=%6.3f\n |<%i|Jb-<Jb>|%i>|^2=%6.3f\n |<%i|Jc-<Jc>|%i>|^2=%6.3f\n",i+1,j+1,u[1]*u[1]+iu[1]*iu[1],i+1,j+1,u[2]*u[2]+iu[2]*iu[2],i+1,j+1,u[3]*u[3]+iu[3]*iu[3]);
            printf(" |<%i|Jd-<Jd>|%i>|^2=%6.3f\n |<%i|Je-<Je>|%i>|^2=%6.3f\n |<%i|Jf-<Jf>|%i>|^2=%6.3f\n",i+1,j+1,u[4]*u[4]+iu[4]*iu[4],i+1,j+1,u[5]*u[5]+iu[5]*iu[5],i+1,j+1,u[6]*u[6]+iu[6]*iu[6]);
            printf(" n%i=%6.3f\n",i,(KB*T)*therm/Z);
         }
      }
   
      // multiply matrix Mab by occupation factor
      for(iJ=1; iJ<=sz; iJ++)
	    u1(iJ) = complex<double> ( u[iJ]*sqrt(therm/Z), iu[iJ]*sqrt(therm/Z) );

      delete[]en;
   }

    // determine number of thermally reachable states
   if (ninit>Hsz)ninit=Hsz;
   //if (pinit<SMALL)pinit=SMALL;
   double zsum=0,zi,x;
   int noft=0; 
   for(i=0; (i<ninit)&((((x=(est[0][i+1].real()-est[0][1].real())/(KB*fabs(T)))<200)? zi=exp(-x):zi=0)>=(pinit*zsum)); ++i)
   {
      noft += Hsz-i; 
      zsum += zi;
   }
   return noft;
}                 

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the matrix elements of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
int dspindensity_coeff1(int &tn,          // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &/*g_J*/,    // Input Lande g-factor
 /* Not Used */       Vector &/*ABC*/,    // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexVector &Slm1,// Output Llm1 vector (1,49)
                      int & xyz,            // Indicating which of x,y,z direction to calculate
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ 
   Vector gjmbH(1,49); gjmbH = 0; for(int i=1; i<=Hxc.Hi(); ++i) { gjmbH(i)=Hxc(i); }
  //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   double dum; dum=gjmbH(2);gjmbH(2)=gjmbH(4);gjmbH(4)=gjmbH(5);gjmbH(5)=gjmbH(3);gjmbH(3)=dum;
   // --------------------------------------------------------------------
   int nt = sdod_du1calc(xyz,tn,T,gjmbH,sipffilename,Slm1,delta,est);
// Slm1(1)=0;
// for(int i=2; i<=6; ++i) Slm1(i) = u1(5+i) *sqrt((2.0*2+1)/8/PI); Slm1(4) *=sqrt(2);
// for(int i=7; i<=15;++i) Slm1(i) = u1(12+i)*sqrt((2.0*4+1)/8/PI); Slm1(11)*=sqrt(2);
// for(int i=16;i<=28;++i) Slm1(i) = u1(23+i)*sqrt((2.0*6+1)/8/PI); Slm1(22)*=sqrt(2);
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the matrix elements of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
int dorbmomdensity_coeff1(int &tn,        // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &/*g_J*/,    // Input Lande g-factor
 /* Not Used */       Vector &/*ABC*/,    // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexVector &Llm1,// Output Llm1 vector (1,49)
                      int & xyz,            // Indicating which of x,y,z direction to calculate
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ 
   Vector gjmbH(1,49); gjmbH = 0; for(int i=1; i<=Hxc.Hi(); ++i) { gjmbH(i)=Hxc(i); }
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   double dum; dum=gjmbH(2);gjmbH(2)=gjmbH(4);gjmbH(4)=gjmbH(5);gjmbH(5)=gjmbH(3);gjmbH(3)=dum;
   // --------------------------------------------------------------------
   int nt = sdod_du1calc(-xyz,tn,T,gjmbH,sipffilename,Llm1,delta,est);
// Llm1(1)=0;
// for(int i=2; i<=6; ++i) Llm1(i) = u1(5+i) *sqrt((2.0*2+1)/8/PI); Llm1(4) *=sqrt(2);
// for(int i=7; i<=15;++i) Llm1(i) = u1(12+i)*sqrt((2.0*4+1)/8/PI); Llm1(11)*=sqrt(2);
// for(int i=16;i<=28;++i) Llm1(i) = u1(23+i)*sqrt((2.0*6+1)/8/PI); Llm1(22)*=sqrt(2);
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// returns operator matrices (n=0 Hamiltonian, n=1,...,nofcomponents: operators of moment components)
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif                                    // on input
int    opmat(int &ni,                      // ni     which operator 0=Hamiltonian, 1,2,3=J1,J2,J3
             char **sipffilename,         // Single ion properties filename
             Vector &Hxc,                 // Hext  vector of external field [meV]
             Vector &Hext,                // Hxc   vector of exchange field [meV]
                                          // on output   
             Matrix &outmat)              // operator matrix of Hamiltonian, I1, I2, I3 depending on n
{  //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   int n=ni;
   if(ni==2){n=3;}
   if(ni==3){n=5;}
   if(ni==4){n=2;}
   if(ni==5){n=4;}
   if(ni==-2){n=-3;}
   if(ni==-3){n=-5;}
   if(ni==-4){n=-2;}
   if(ni==-5){n=-4;}

   // Parses the input file for parameters
   const char *filename = sipffilename[0];
   icpars pars; ic_parseinput(filename,pars);

   //          0 1 2 3 4 5 6  7  8 91011 12 13 1415161718 19 20 21 222324252627 28 29 30 31 32333435363738 39 40 41 42 43 4445464748495051
   int K[] = {-1,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int Q[] = {-1,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {-1,0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};

   if(n==0)                               // return Hamiltonian
   {
      std::vector<double> gjmbH(max(6,Hxc.Hi()),0.); for(int i=1; i<=Hxc.Hi(); i++) gjmbH[i-1]=-Hxc(i);
      //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   double dum; dum=gjmbH[2];gjmbH[2]=gjmbH[4];gjmbH[4]=gjmbH[5];gjmbH[5]=gjmbH[3];gjmbH[3]=dum;
   // --------------------------------------------------------------------
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH[1]-=MUB*Hext(1); gjmbH[0]-=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH[3]-=MUB*Hext(2); gjmbH[2]-=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH[5]-=MUB*Hext(3); gjmbH[4]-=GS*MUB*Hext(3); }
      #ifdef JIJCONV
      if(pars.B.norm().find("Stevens")!=std::string::npos) {
         pars.jijconvcalc();
         for(int i=0; i<Hxc.Hi(); i++) gjmbH[i] *= pars.jijconv[i]; }
      #endif

      // check dimensions of vector
      if(Hxc.Hi()>52) {
         fprintf(stderr,"Error module ic1ion: dimension of exchange field=%i > 52 - check number of columns in file mcphas.j\n",Hxc.Hi()); exit(EXIT_FAILURE); }

      sMat<double> Hcfi, Hcf = icf_hmltn(Hcfi, pars); Hcf/=MEV2CM; Hcfi/=MEV2CM;
      for(int ind=1; ind<=(int)gjmbH.size(); ind++)
      {
         if(ind<=6) {       // Calculates Sx,Lx etc and copies them to est too. 
            if(ind==3 || ind==4) Hcfi += icf_mumat(pars.n, ind-1, pars.l)*gjmbH[ind-1]; else Hcf += icf_mumat(pars.n, ind-1, pars.l)*gjmbH[ind-1]; }
         else {             // Calculates multipolar operator matrices and copies them to est. 
            if(im[ind]==1) Hcfi += icf_ukq(pars.n,K[ind],Q[ind],pars.l)*gjmbH[ind-1];   else Hcf += icf_ukq(pars.n,K[ind],Q[ind],pars.l)*gjmbH[ind-1]; }
      }
      
      sMat<double> tmp = Hcf+Hcfi;
      std::vector< std::vector<int> > u = tmp.findlower();
      Matrix retval(1,tmp.nr(),1,tmp.nc()); retval=0;
      for (int j=0; j<(int)u.size(); j++)
      {
         retval(u[j][0]+1,u[j][1]+1) = Hcf(u[j][0],u[j][1]);
         retval(u[j][1]+1,u[j][0]+1) = Hcfi(u[j][0],u[j][1]);
      }
      outmat = retval; return 0;
   }
   else if(n>0)
   {  
      if(n>51) {
         fprintf(stderr,"Error module ic1ion: operatormatrix index=%i > 51 - check number of columns in file mcphas.j\n",n); exit(EXIT_FAILURE); }

      sMat<double> Jmat; if(n<=6) Jmat = icf_mumat(pars.n, n-1, pars.l); else Jmat = icf_ukq(pars.n,K[n],Q[n],pars.l);
      Matrix retval(1,Jmat.nr(),1,Jmat.nc()); retval=0;
      if(im[n]==1)
      {
         std::vector< std::vector<int> > u = Jmat.findlower();
         for (int j=0; j<(int)u.size(); j++)
            retval(u[j][1]+1,u[j][0]+1) = Jmat(u[j][0],u[j][1]);
      }
      else
      {
         std::vector< std::vector<int> > u = Jmat.findlower();
         for (int j=0; j<(int)u.size(); j++)
            retval(u[j][0]+1,u[j][1]+1) = Jmat(u[j][0],u[j][1]);
      }
      outmat = retval; return 0;
   }
   else
   {  
      if(n<-51) {
         fprintf(stderr,"Error module ic1ion: operatormatrix index=%i < -51 - check number of columns in file mcphas.j\n",n); exit(EXIT_FAILURE); }

      // Calculates the single-ion (without external field) Hamiltonian and operator matrix in LS-basis
      sMat<double> Jmat; if(abs(n)>=6) Jmat = icf_mumat(pars.n, abs(n)-1, pars.l); else Jmat = icf_ukq(pars.n,K[abs(n)],Q[abs(n)],pars.l);
      sMat<double> Hcfi, Hcf = icf_hmltn(Hcfi, pars); Hcf/=MEV2CM; Hcfi/=MEV2CM;
      int Hsz=Hcf.nr(); sMat<double> zeros(Hsz,Hsz);
      Matrix retval(1,Hsz,1,Hsz); retval=0;

      // Diagonalise H_singleion, remove small (==zero) elements, converts operator matrix to dense format to use BLAS matrix multiplication.
      complexdouble *Vf; Vf = new complexdouble[Hsz*Hsz]; double *Ef; Ef = new double[Hsz];
      int info = ic_diag(Hcf,Hcfi,Vf,Ef); if(info!=0) { std::cerr << "icf1ion::opmat(): Error diagonalising, info==" << info << "\n"; }
      delete[]Ef;
      for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++) {
         if(fabs(Vf[ii*Hsz+jj].r)<DBL_EPSILON) Vf[ii*Hsz+jj].r=0.; if(fabs(Vf[ii*Hsz+jj].i)<DBL_EPSILON) Vf[ii*Hsz+jj].i=0.; }
      complexdouble *zJmat; if(im[(int)abs(n)]==1) zJmat=zmat2f(zeros,Jmat); else zJmat=zmat2f(Jmat,zeros);

      // Rotates the operator matrix into basis where H_singleion is diagonal using calculated eigenvectors with BLAS routines.
      complexdouble *Hrot,*zmt; zmt = new complexdouble[Hsz*Hsz]; Hrot = new complexdouble[Hsz*Hsz];
      char notranspose='N',transpose='C',uplo='U',side='L'; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
      F77NAME(zhemm)(&side,&uplo,&Hsz,&Hsz,&zalpha,zJmat,&Hsz,Vf,&Hsz,&zbeta,zmt,&Hsz);
      F77NAME(zgemm)(&transpose,&notranspose,&Hsz,&Hsz,&Hsz,&zalpha,Vf,&Hsz,zmt,&Hsz,&zbeta,Hrot,&Hsz); free(zJmat);

      // Populates the output matrix and frees temporary dense matrices (arrays).
      for(int ii=0; ii<Hsz; ii++) for(int jj=ii; jj<Hsz; jj++) {
         if(fabs(Hrot[ii*Hsz+jj].r)>DBL_EPSILON) retval(ii+1,jj+1)=Hrot[ii*Hsz+jj].r;
         if(fabs(Hrot[ii*Hsz+jj].i)<DBL_EPSILON) retval(jj+1,ii+1)=Hrot[ii*Hsz+jj].i;
      }
      delete[]Vf; delete[]Hrot; delete[]zmt;
     outmat=retval;return 0;
   }
   std::cerr << "icf1ion::opmat - failed to calculate operator matrices\n"; return 1; //exit(EXIT_FAILURE);
}

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
   FILEOUT << "# Spin orbit parameter (" << pars.e_units << "): zeta=" << pars.xi << "\n";
   FILEOUT << "# Crystal Field parameters normalisation: " << pars.B.norm() << "\n";
   std::string norm=pars.B.norm(); strtolower(norm); if(norm.find("stev")!=std::string::npos)
   {
      std::string op; if(pars.B.op_equiv==Lt) op.assign("L"); else op.assign("J");
      FILEOUT << "# Stevens Factors: <" << op << "||alpha||" << op << ">=" << pars.B.alpha();
      FILEOUT <<                  ", <" << op << "||beta||" << op << ">=" << pars.B.beta();
      if(pars.l==F) FILEOUT << ", <" << op << "||gamma||" << op << ">=" << pars.B.gamma(); FILEOUT << "\n";
   }
   FILEOUT << "# Crystal Field parameters (" << pars.B.units() << "): " << pars.B.cfparsout(", ") << "\n";
   if(fabs(pars.Bx)>DBL_EPSILON || fabs(pars.By)>DBL_EPSILON || fabs(pars.Bz)>DBL_EPSILON)
   {
      FILEOUT << "# With magnetic field: Bx=" << pars.Bx << ", By=" << pars.By << ", Bz=" << pars.Bz << " Tesla.\n";
   }
   if(fabs(pars.Dx2)>DBL_EPSILON || fabs(pars.Dy2)>DBL_EPSILON || fabs(pars.Dz2)>DBL_EPSILON)
      FILEOUT << "# With spin anisotropy: Dx2=" << pars.Dx2 << ", Dy2=" << pars.Dy2 << ", Dz2=" << pars.Dz2 << " " << pars.e_units << ".\n";
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
   int ii,nV=pars.num_eigv>ns?ns:pars.num_eigv; 
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
    //if(real(est(0,iE+1))==0.) if(iE<(ns-1) && real(est(0,iE+2))==0.) break; 
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
         for(iV=1; iV<nV; iV++)
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
         for(iV=1; iV<nV; iV++)
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
   Vector gjmbHxc(1,6,0.); ComplexMatrix est; double T=1.;
   Vector Hext(1,3);Hext(1)=pars.Bx;Hext(2)=pars.By;Hext(3)=pars.Bz;
   estates(&est,gjmbHxc,Hext,T,T,gjmbHxc,&infile);
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
      Vector gjmbH0(1,3,0.),J(1,6,0.);
      if(pars.xT==0.) gjmbH0(1)=pars.xHa/xnorm; else gjmbH0(1)=pars.yHa/ynorm; 
      if(pars.xT==0.) gjmbH0(2)=pars.xHb/xnorm; else gjmbH0(2)=pars.yHb/ynorm; 
      if(pars.xT==0.) gjmbH0(3)=pars.xHc/xnorm; else gjmbH0(3)=pars.yHc/ynorm; 
      double convfact=1; if(pars.mag_units==1) convfact = NAMUB*1e3; else if(pars.mag_units==2) convfact = NAMUB;  // 1==cgs, 2==SI
      
      // Reinitialises the estates matrix
      est.Remove(); Icalc_parameter_storage_matrix_init(&est,gjmbHxc,Hext,&T,&T,gjmbHxc,&infile);

      for(int iH=0; iH<nH; iH++)
      {
         for(int al=1; al<=3; al++) Hext(al) = gjmbH0(al)*(Hmin+iH*Hstep);
         for(int iT=0; iT<nT; iT++)
         {
            Icalc(J,&vT[iT],gjmbHxc,Hext,&vT[iT],gjmbHxc,&infile,&lnZ,&U,est); for(int iJ=1; iJ<=6; iJ++) if(fabs(J(iJ))<SMALL) J(iJ)=0.;
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

#ifdef _INTEGRAL
   clock_t start,end; end = clock();
   Vector J(1,6,0.), gmbH(1,6,.0578838263), ABC; gmbH[1]*=2; gmbH[3]*=2; gmbH[5]*=2; 
   double /*T=2.0,*/lnZ=0.,U=0.,gJ=0.; T=2.0;
   char *filearray[1]; 
   filearray[0] = infile;
 /*ComplexMatrix est;*/ est.Remove(); Icalc_parameter_storage_matrix_init(&est,gmbH,&gJ,&T,ABC,filearray);
 //ComplexMatrix est; int Hsz=getdim(pars.n,pars.l); est = ComplexMatrix(0,Hsz,0,Hsz);
   end = clock();

   Icalc(J,&T,gmbH,&gJ,ABC,filearray,&lnZ,&U,est);
   start = clock(); std::cerr << "Time to do Icalc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";
   std::cerr << "lnZ = " << lnZ << ", U = " << U << "\n";
   std::cerr << "J[1] = " << J[1] << ", J[2] = " << J[2] << ", J[3] = " << J[3] << ", J[4] = " << J[4] << ", J[5] = " << J[5] << ", J[6] = " << J[6] << "\n";

//for (int it=0; it<1000; it++) {
// end = clock();
// Icalc(J,&T,gmbH,&gJ,ABC,filearray,&lnZ,&U,est);
// start = clock(); std::cerr << "Time to do Icalc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";
// std::cerr << "lnZ = " << lnZ << ", U = " << U << "\n";
// std::cerr << "J[1] = " << J[1] << ", J[2] = " << J[2] << ", J[3] = " << J[3] << ", J[4] = " << J[4] << ", J[5] = " << J[5] << ", J[6] = " << J[6] << "\n";
// if(it%100==0) { std::cerr << it << " "; } } std::cerr << "\n";
 //gmbH[1] = 0.; gmbH[2] = 0.; gmbH[3] = 0.; gmbH[4] = 0.; gmbH[5] = 0.; gmbH[6] = 0.; 
   est.Remove(); estates(&est,gmbH,gJ,T,ABC,filearray);
   end = clock(); std::cerr << "Time to do estates() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   
   int imq, tn = 2,n,nd; float delta=0.; ComplexMatrix mat6(1,6,1,6);
   imq = du1calc(tn,T,gmbH,gJ,ABC,filearray,mat6,delta,n,nd,est);
   start = clock(); std::cerr << "Time to calculate du1calc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";

   ComplexVector Mq;
   double th=PI/4, ph=PI/4, J0=1., J2=1., J4=1., J6=1.;
 //double th=0., ph=0., J0=1., J2=1., J4=0., J6=0.;
   mq(Mq,th,ph,J0,J2,J4,J6,est);
   end = clock(); std::cerr << "Time to calculate mq() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.";
   std::cerr << " Mq = [" << Mq[1].real() << "+" << Mq[1].imag() << "i "
                          << Mq[2].real() << "+" << Mq[2].imag() << "i "
                          << Mq[3].real() << "+" << Mq[3].imag() << "i]\n";

/* ComplexMatrix mat(1,6,1,6);
   imq = dv1calc(tn,th,ph,J0,J2,J4,J6,est,T,mat);
   start = clock(); std::cerr << "Time to calculate dv1calc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n"; */
#endif

   return 0;
}
