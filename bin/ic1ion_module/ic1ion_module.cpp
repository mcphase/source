/* ic1ion_module.cpp
 *********************************************************
 * Functions:
 *   void myPrintMatrix(FILE * file,sMat<double> & M,int d)                        // Prints out full matrix
 *   bool checkmat(ComplexMatrix &cmat, complexdouble *fmat,int r, int c)          // Compares Matpack and fortran matrices
 *   void Icalc(Vector &J, double *T, Vector &gjmbHxc,Vector&Hext, double *gJ, Vector &ABC,      // Calculates the meanfield moment
 *                 char **sipffile, double *lnZ, double *U, ComplexMatrix &est)    //
 *   void mcalc(Vector &J, double *T, Vector &gjmbHxc,Vector&Hext, double *gJ, Vector &ABC,      // Calculates the meanfield moment
 *                 char **sipffile,  ComplexMatrix &est)    //
 *   void Lcalc(Vector &J, double *T, Vector &gjmbHxc,Vector&Hext, double *gJ, Vector &ABC,      // Calculates the meanfield moment
 *                 char **sipffile,  ComplexMatrix &est)    //
 *   void Scalc(Vector &J, double *T, Vector &gjmbHxc,Vector&Hext, double *gJ, Vector &ABC,      // Calculates the meanfield moment
 *                 char **sipffile,  ComplexMatrix &est)    //
 *   int du1calc(int &tn, double &T, Vector &gjmbHxc,Vector&Hext, double &g_J, Vector &ABC,       // Calculates the transition
 *                 char **sipffilename, ComplexVector & u1, float &delta,          //   matrix elements
 *                 ComplexMatrix &est)                                             //
 *   void Icalc_parameter_storage_matrix_init(ComplexMatrix *est,Vector &gjmbHxc,Vector&Hext // initialises parameter storage matrix 
 *                 double *g_J,double &T,Vector &ABC,char **sipffilename)          // for Icalc
 *   void estates(ComplexMatrix *est, Vector &gjmbHxc,Vector&Hext, double *g_J, double &T,    // Calculates the energy and wavefunctions
 *                 Vector &ABC, char **sipffilename)                               //   "estates" matrix
 *   bool get_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l,         // Calculates/loads the Q_q operators
 *                 std::vector<double> &Jvec)                                      //   for beyond dipole calculations
 *   void save_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l,        // Saves the Q_q operators to a file
 *                 std::vector<double> &Jvec)                                      //
 *   void mq(ComplexVector &Mq, double &th, double &ph, double &J0, double &J2,    // Calculates the thermal expectation
 *                 double &J4, double &J6, ComplexMatrix &est)                     //   of the magnetisation density
 *   int dv1calc(int &tn, double &th, double &ph, double &J0, double &J2,           // Calculates the transition matrix
 *                 double &J4, double &J6, ComplexMatrix &est, double &T,          //   elements beyond the dipole
 *                 ComplexVector & v1)                                             //   approximation.
 *   void chargedensity_coeff(Vector & mom,double *T,Vector &gjmbHxc,Vector&Hext,      // Calc. coeffs. of expansion of chargedensity 
 *                 double *gJ,Vector &ABC, char **sipffile, ComplexMatrix &est)    //   in terms of Zlm R^2(r) at given T / H_eff
 *   void spindensity_coeff(Vector & mom, int & xyz,double *T,Vector &gjmbHxc,Vector&Hext,      // Calc. coeffs. of expansion of spindensity 
 *                 double *gJ,Vector &ABC, char **sipffile, ComplexMatrix &est)    //   in terms of Zlm R^2(r) at given T / H_eff
 *   void orbmomdensity_coeff(Vector & mom,int & xyz, double *T,Vector &gjmbHxc,Vector&Hext,    // Calc. coeffs. of expansion of orbital moment density
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

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in truncate.cpp
// --------------------------------------------------------------------------------------------------------------- //
void truncate_hmltn(icpars &pars, ComplexMatrix &est, sMat<double> &Hic, sMat<double> &iHic, int JHi, int JLo);
void truncate_expJ(icpars &pars, ComplexMatrix &est, Vector &gjmbH, Vector &J, double T, double *lnZ, double *U, complexdouble *Jm);

void myPrintMatrix(FILE * file,sMat<double> & M,int d)
{int i1,j1;
    fprintf (file,"Matrix\n");
   for (i1=0;i1<=d;++i1){
    for (j1=0;j1<=d;++j1) fprintf (file,"%6.3f ",M(i1,j1));
    fprintf (file,"\n");
    }
}    

// --------------------------------------------------------------------------------------------------------------- //
// Checks whether the Matpack matrix is the same as the c-array
// --------------------------------------------------------------------------------------------------------------- //
bool checkmat(ComplexMatrix &cmat, complexdouble *fmat,int r, int c)
{
   int i,j;
   for(i=0; i<(cmat.Rows()-r); i++)
      for(j=0; j<(cmat.Cols()-c); j++)
      { 
       //std::cout << "cmat["<<i+r<<"]["<<j+c<<"]=" << cmat[i+1][j+1] << "\tfmat["<<j<<"*" << cmat.Cols()-c << "+"<<i<<"]=";
       //   std::cout << fmat[j*(cmat.Cols()-c)+i].r << "+" << fmat[j*(cmat.Cols()-c)+i].i << "i\n";
       //std::cout << "cmat["<<j+c<<"]["<<i+r<<"]=" << cmat[j+c][i+r] << "\tfmat["<<j<<"*" << cmat.Cols()-c << "+"<<i<<"]=";
       //   std::cout << fmat[j*(cmat.Cols()-c)+i].r << "+" << fmat[j*(cmat.Cols()-c)+i].i << "i\n";
    //   if(real(cmat[i+r][j+c])!=fmat[j*(cmat.Cols()-c)+i].r || imag(cmat[i+r][j+c])!=fmat[j*(cmat.Cols()-c)+i].i) return false; }
         if(real(cmat[j+c][i+r])!=fmat[j*(cmat.Cols()-c)+i].r || imag(cmat[j+c][i+r])!=fmat[j*(cmat.Cols()-c)+i].i) return false; }
   return true;
}


// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the magnetic moment at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void Icalc(Vector &J,          // Output single ion momentum vector <Ja>,<Jb>,<Jc>, etc.
                      double *T,          // Input scalar temperature
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      double *lnZ,        // Output scalar logarithm of partition function
                      double *U,          // Output scalar internal energy 
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
//#ifndef __linux__
// int thrid = GetCurrentThreadId();
//#else
// pthread_t thrid = pthread_self();
//#endif
  // sum exchange field and external field
  Vector gjmbH(1,gjmbHxc.Hi());
  gjmbH=gjmbHxc;
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

   // Prints out single ion parameters and filename.
// std::cout << "#ic1ion: Single ion parameter file: " << filename << "\n";
// std::cout << "#ic1ion: Read in parameters (in " << pars.e_units << "): F2=" << pars.F[1] << ",F4=" << pars.F[2];
// if(pars.l==F) std::cout << ",F6=" << pars.F[3]; std::cout << ",xi=" << pars.xi << "\n";
// std::cout << "#ic1ion: \t" << pars.B.cfparsout(", ") << "in " << pars.B.units() << "\n";

   // Converts the Jij parameters if necessary
   std::vector<double> vgjmbH((J.Hi()-J.Lo()+1),0.); 
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc();
      for(int i=J.Lo(); i<=J.Hi(); i++) vgjmbH[i-J.Lo()] = -gjmbH[i]*pars.jijconv[i]; }
   else
   #endif
      for(int i=J.Lo(); i<=J.Hi(); i++) vgjmbH[i-J.Lo()] = -gjmbH[i];

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
      for(i=0; i<(int)(parval.size()/2); i++) if(real(est[0][i+1])!=parval[2*i] || imag(est[0][i+1])!=parval[2*i+1]) { Hicnotcalc = true; break; }
   }
   else Hicnotcalc = true;
   if(Hicnotcalc)
   {
      sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars); Hic/=MEV2CM; iHic/=MEV2CM; H = zmat2f(Hic,iHic);
//    est = ComplexMatrix(0,Hsz,0,Hsz); I comment this out - you should not reinitialize est !!!
      if( (est.Rhi()!=Hsz||est.Chi()!=Hsz) && pars.truncate_level==1)
      {
         std::cerr << "ERROR module ic1ion - Icalc: Hsz recalculation does not agree with eigenstates matrix dimension\n"; exit(EXIT_FAILURE);
      }
      if (!((int)(parval.size()/2)>Hsz && pars.truncate_level==1))
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
   if(pars.truncate_level!=1)  // Uses the eigenvectors of the single ion Hamiltonian to truncate the matrix
      truncate_expJ(pars,est,gjmbH,J,*T,lnZ,U,Jm);
   else
   {
      // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
      icmfmat mfmat(pars.n,pars.l,J.Hi()-J.Lo()+1,pars.save_matrices,pars.density);
      #ifdef JIJCONV
      if(pars.B.norm().find("Stevens")!=std::string::npos) mfmat.jijconv.assign(pars.jijconv.begin(),pars.jijconv.end());
      #endif
      sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH,pars.save_matrices); 
      complex<double> a(1.,0.); int incx = 1;
      Jm = zmat2f(Jmat,iJmat); for(i=1; i<=Hsz; i++) F77NAME(zaxpy)(&Hsz,(complexdouble*)&a,(complexdouble*)&est[i][1],&incx,&Jm[(i-1)*Hsz],&incx);

      // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
      iceig VE; if(pars.partial) VE.lcalc(pars,Jm);
      #ifndef NO_ARPACK
      else if(pars.arnoldi) VE.acalc(pars,Jm); 
      #endif
      else VE.calc(Hsz,Jm); free(Jm);

      // Calculates the expectation values sum_n{ <n|Ja|n> exp(-En/kT) }
      std::vector< std::vector<double> > matel; std::vector<double> vJ = mfmat.expJ(VE,*T,matel,pars.save_matrices);
      for(i=J.Lo(); i<=J.Hi(); i++) J[i] = vJ[i-J.Lo()]; 
      *lnZ = vJ[J.Hi()-J.Lo()+1]; *U = vJ[J.Hi()-J.Lo()+2];
   }
}
// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the magnetic moment at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void mcalc(Vector &mom,          // Output magnetic moment (mub)
                      double *T,          // Input scalar temperature
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{Vector J(1,6),ABC; double gJ=0.,lnZ,U;
  Icalc(J,T,gjmbHxc,Hext,&gJ,ABC,sipffilename,&lnZ,&U,est);
 mom(1)=GS*J(1)+J(2);
 mom(2)=GS*J(3)+J(4);
 mom(3)=GS*J(5)+J(6);

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
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{Vector J(1,6),ABC; double gJ=0.,lnZ,U;
  Icalc(J,T,gjmbHxc,Hext,&gJ,ABC,sipffilename,&lnZ,&U,est);
 L(1)=J(2);
 L(2)=J(4);
 L(3)=J(6);

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
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{Vector J(1,6),ABC; double gJ=0.,lnZ,U;
  Icalc(J,T,gjmbHxc,Hext,&gJ,ABC,sipffilename,&lnZ,&U,est);
 S(1)=J(1);
 S(2)=J(3);
 S(3)=J(5);

}
// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           int du1calc(int &tn,            // Input transition number; if tn>0, print debug info
                      double &T,          // Input temperature
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double &/*g_J*/,    // Input Lande g-factor
 /* Not Used */       Vector &/*ABC*/,    // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexVector & u1, // Output u1 vector
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{  // sum exchange field and external field
  Vector gjmbH(1,gjmbHxc.Hi());
  gjmbH=gjmbHxc;
   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(Hext(1))>DBL_EPSILON || fabs(Hext(2))>DBL_EPSILON || fabs(Hext(3))>DBL_EPSILON)
   {
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH(2)+=MUB*Hext(1); gjmbH(1)+=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH(4)+=MUB*Hext(2); gjmbH(3)+=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH(6)+=MUB*Hext(3); gjmbH(5)+=GS*MUB*Hext(3); }
   }
  
   int i,j,k;
 
   // check if printout should be done and make tn positive
   int pr=1; if (tn<0) { pr=0; tn*=-1; }
   double ninit=u1[1].real();
   double pinit=u1[1].imag();
   // Copies the already calculated energy levels / wavefunctions from *est
   if(est.Rows()!=est.Cols()) { std::cerr << "du1calc(): Input rows and columns of eigenstates matrix don't match.\n"; return 0; }
   int Hsz = est.Rows()-1;
   j=0; k=0; for(i=0; i<Hsz; ++i) { for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   if(est[0][j+1].real()-est[0][i+1].real()<delta)
   {
      double *en = new double[Hsz]; for(k=0; k<Hsz; k++) en[k] = est[0][k+1].real();

      // Parses the input file for parameters
      icpars pars; const char *filename = sipffilename[0];
      ic_parseinput(filename,pars);

      // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
      int num_op = gjmbH.Hi()-gjmbH.Lo()+1; icmfmat mfmat(pars.n,pars.l,(num_op>6?num_op:6),pars.save_matrices);

      iceig VE(Hsz,en,(complexdouble*)&est[1][0],1);
 
      // Calculates the transition matrix elements:
      //    u1 = <i|Ja|j> * sqrt[(exp(-Ei/kT)-exp(-Ej/kT)) / Z ]   if delta > small
      //    u1 = <i|Ja-<Ja>|j> * sqrt[(exp(-Ei/kT)) / kTZ ]             if delta < small (quasielastic scattering)
      //    See file icpars.cpp, function mfmat::Mab() to see the actual code to calculate this.
  
      std::vector<double> u((num_op>6?num_op:6)+1), iu((num_op>6?num_op:6)+1);
      mfmat.u1(u,iu,VE,T,i,j,pr,delta,pars.save_matrices);

      for(i=1; i<=(num_op>6?num_op:6); i++)
            u1(i) = complex<double> (u[i], iu[i]);
   }
   // determine number of thermally reachable states
   if (ninit>Hsz)ninit=Hsz;
   if (pinit<SMALL)pinit=SMALL;
   double zsum=0,zi;
   int noft=0;for(i=0;(i<ninit)&((zi=(exp(-(est[0][i+1].real()-est[0][1].real())/(KB*fabs(T)))))>(pinit*zsum));++i){noft+=Hsz-i-1;zsum+=zi;}
//   int noft=0;for(i=0;(i<Hsz)     &((exp(-(est[0][i+1].real()-est[0][1].real())/(KB*T)))>SMALL);++i)noft+=Hsz-i-1; // removed MR  6.9.2011 to allow for mcdisp options -ninit -pinit   return noft;
   return noft;
   //return Hsz*(Hsz-1)/2;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to initialise the storage matrix "est" for Icalc
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
          void Icalc_parameter_storage_matrix_init(
                      ComplexMatrix *est, // Output Eigenstates matrix (row 0: real==Eigenvalues;imag==population)
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,   // Input  Lande g-factor
                      double * /*T*/,     // Input  temperature
 /* Not Used */       Vector & /*ABC*/,   // Input  Vector of parameters from single ion property file
                      char **sipffilename)// Input  Single ion properties filename
{ // sum exchange field and external field
  Vector gjmbH(1,gjmbHxc.Hi());
  gjmbH=gjmbHxc;
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

   // If we just want a blank estates matrix for later use (e.g. in Icalc)
   int Hsz = getdim(pars.n,pars.l); 
   if(pars.truncate_level==1)
   {
      (*est) = ComplexMatrix(0,Hsz,0,Hsz);
      // Stores the number of electrons and the orbital number in element (0,0)
      (*est)(0,0) = complex<double> (pars.n, pars.l);
   }
   else
   {
      int Jlo=gjmbH.Lo(), Jhi=gjmbH.Hi(), cb = (int)(pars.truncate_level*(double)Hsz), matsize=cb*cb;
      for (int ii=Jlo; ii<=Jhi; ii++) matsize += (cb*cb);
      (*est) = ComplexMatrix(0,matsize+100,0,0);
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the eigenstates of Hic+effective_mean_fields
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
          void estates(ComplexMatrix *est,// Output Eigenstates matrix (row 0: real==Eigenvalues;imag==population)
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,   // Input  Lande g-factor
                      double &T,          // Input  temperature
 /* Not Used */       Vector & /*ABC*/,   // Input  Vector of parameters from single ion property file
                      char **sipffilename)// Input  Single ion properties filename
{// sum exchange field and external field
  Vector gjmbH(1,gjmbHxc.Hi());
  gjmbH=gjmbHxc;
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

   // Calculates the IC Hamiltonian matrix
   sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars); int Hsz = Hic.nr();
 
   // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
   int num_op = gjmbH.Hi()-gjmbH.Lo()+1; icmfmat mfmat(pars.n,pars.l,(num_op>6?num_op:6),pars.save_matrices);
   int i,j,gLo=gjmbH.Lo(),gHi=gjmbH.Hi(); std::vector<double> vgjmbH(gHi,0.);
   for(i=gLo; i<=gHi; i++) vgjmbH[i-1] = -gjmbH[i];
   // Converts the Jij parameters if necessary
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc(); mfmat.jijconv.assign(pars.jijconv.begin(),pars.jijconv.end());
      for(i=gLo; i<=gHi; i++) vgjmbH[i-1] *= pars.jijconv[i]; }
   #endif
   sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH,pars.save_matrices); 

   // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
   Hic/=MEV2CM; Hic+=Jmat; if(!iHic.isempty()) iHic/=MEV2CM; if(!iJmat.isempty()) iHic+=iJmat; 
   iceig VE; if(iHic.isempty()) VE.calc(Hic); else VE.calc(Hic,iHic);

   // Initialises the output matrix
   (*est) = ComplexMatrix(0,Hsz,0,Hsz);

   // Stores the number of electrons and the orbital number in element (0,0)
   (*est)(0,0) = complex<double> (pars.n, pars.l);

   // Puts eigenvectors/values into the est matrix
   for(i=0; i<Hsz; i++) (*est)(0,i+1) = complex<double> (VE.E(i), exp(-(VE.E(i)-VE.E(0))/(KB*T)));   // Row 0

   if(VE.iscomplex()) {
      for(i=1; i<=Hsz; i++) memcpy(&(*est)[i][1],VE.zV(i-1),Hsz*sizeof(complexdouble));  }
//    std::cout << "\n\nest==VE.zV = " << checkmat((*est),VE.zV(0),1,1) << endl << endl; }
   else
      for(i=0; i<Hsz; i++){//printf("\n");
         for(j=0; j<Hsz; j++) 
            {(*est)(i+1,j+1) = complex<double> (VE.V(i)[j], 0.);
//                         printf("%6.3f %+6.3f i  ",VE.V(i)[j],0.0);
             if(VE.V(i)[j]!=VE.V(j,i)){fprintf(stderr,"compiler problem: bad memory mapping of vectors\n");exit(EXIT_FAILURE);}
             if(VE.V(i)[j]!=(*est)[i+1][j+1]){fprintf(stderr,"compiler problem: bad memory mapping of vectors\n");exit(EXIT_FAILURE);}
            }}

   end = clock(); std::cerr << "Time to do estates() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";

}

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
   n = (int)est[0][0].real(); i = (int)est[0][0].imag(); l = (orbital)i;
   if(i>3 || i<0) { std::cerr << "ic1ion mq(): Error only s-, p-, d-, and f-electrons supported.\n"; exit(EXIT_FAILURE); }
   std::vector<double> E,Jvec(6,0.); Jvec[0]=th; Jvec[1]=ph; Jvec[2]=J0; Jvec[3]=J2; Jvec[4]=J4; Jvec[5]=J6;
   std::vector< sMat<double> > Qp, Qm; 
   std::vector< std::vector< sMat<double> > > Qmat; for(i=0; i<3; i++) Qmat.push_back(Qp);
   complexdouble *zQmat, *zt, zme, zalpha, zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0;
   double zMqr,zMqi,Z=0.;
   char trans = 'U'; int incx=1;

   Mq = ComplexVector(1,3);

   if(!get_Qq(Qmat[0],0,n,l,Jvec) || !get_Qq(Qmat[1],1,n,l,Jvec))            // Qmat[0]==Qx, Qmat[1]==Qy, Qmat[2]==Qz
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
//         printf ("%i zme=%g %+g i  Ei=%6.3f ni=%6.3f \n",i,zme.r,zme.i,est[0][i].real(),est[0][i].imag());
         zMqr += (-2.)*zme.r*est[0][i].imag(); zMqi += (-2.)*zme.i*est[0][i].imag(); if(q==0) Z += est[0][i].imag();
      }
      free(zQmat); free(zt); Mq[q+1] = complex<double> (zMqr, zMqi)/Z;
   }
//   printf("MQ=(%g %+g i, %g %+g i,%g %+g i)\n",real(Mq(1)),imag(Mq(1)),real(Mq(2)),imag(Mq(2)),real(Mq(3)),imag(Mq(3)));
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the transition matrix using the scattering operator of Balcar and Lovesey
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           int                            // Returns total number of transitions
               dv1calc(int &tn,            // Input transition number |tn|. If tn<0 omit printout. If tn>0 print info.
                  double &th,             // Input zenith angle (with the z==b axis) in radians.
                  double &ph,             // Input azimuth angle (with the x==a axis, to projection in x-y plane).
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est,     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
                  double &T,              // Input temperature (K)
                  ComplexVector & v1)     // Output transition vector, v1(alpha) = <-|Qalpha-<Qalpha>|+> sqrt(n- - n+)
/* 
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
*/
{
   // check if printout should be done and make tn positive
   int pr=1; if (tn<0) { pr=0; tn*=-1; }
   double ninit=v1[1].real();
   double pinit=v1[1].imag();

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

   int a,j=0,k=0; for(i=0; i<Hsz; ++i) {for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
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

//      if(i==j)                               //subtract thermal expectation value from zij=zii
//      {
//         complexdouble expQ; double thexp=0;
//         for(iJ=1;iJ<=Hsz;++iJ)
//         {
//            therm = exp(-(est[0][iJ].real()-est[0][1].real())/(KB*T)); if(therm<DBL_EPSILON) break;
//            F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[iJ][1], &incx, &zbeta, zt, &incx);
#ifdef _G77
//            F77NAME(zdotc)(&expQ, &Hsz, (complexdouble*)&est[iJ][1], &incx, zt, &incx);
#else
//            expQ = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[iJ][1], &incx, zt, &incx);
#endif
//            thexp += expQ.r * therm / Z;
//         }
//         zij[2*q+1].r-=thexp;zji[2*q+1].r-=thexp;
//      }
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
      if(i==j)                             //subtract thermal expectation value from zij=zii
      {                                    //MR120120 ... reintroduced
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
      if(fabs(zij[iJ].i+zji[iJ].i)>SMALL) { std::cerr << "ERROR module ic1ion - dv1calc: <i|Qalpha|j>not hermitian\n"; exit(EXIT_FAILURE); }

                
   complex<double> im(0,1);
   ComplexVector iQalphaj(1,6);
   
   for(a=1; a<=6; a++){iQalphaj(a) = complex<double> (zij[a].r,zij[a].i);if(a%2==1){iQalphaj(a)*=0.5;}} 
                                                                         // divide spin part by 2
   v1 = 0;
   for(a=1; a<=6; a++)
         v1(a) = iQalphaj(a);

   double delta;
   delta = est[0][j].real()-est[0][i].real();
   if(delta<-0.000001) { std::cerr << "ERROR module ic1ion - dv1calc: energy gain delta gets negative\n"; exit(EXIT_FAILURE); }

   if(j==i) delta = -SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc

   // do some printout if wishes and set correct occupation factor
   if (delta>SMALL)
   {
      therm = exp(-(est[0][i].real()-est[0][1].real())/(KB*T)) - exp(-(est[0][j].real()-est[0][1].real())/(KB*T));
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV",i,j,delta);
         printf(" |<%i|Qa|%i>|^2=%6.3f |<%i|Qb|%i>|^2=%6.3f |<%i|Qc|%i>|^2=%6.3f",i,j,abs(v1(1))*abs(v1(1)),i,j,abs(v1(2))*abs(v1(2)),i,j,abs(v1(3))*abs(v1(3)));
         printf(" |<%i|Qd|%i>|^2=%6.3f |<%i|Qe|%i>|^2=%6.3f |<%i|Qf|%i>|^2=%6.3f",i,j,abs(v1(4))*abs(v1(4)),i,j,abs(v1(5))*abs(v1(5)),i,j,abs(v1(6))*abs(v1(6)));
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
         printf(" |<%i|Qa|%i>|^2=%6.3f |<%i|Qb|%i>|^2=%6.3f |<%i|Qc|%i>|^2=%6.3f",i,j,abs(v1(1))*abs(v1(1)),i,j,abs(v1(2))*abs(v1(2)),i,j,abs(v1(3))*abs(v1(3)));
         printf(" |<%i|Qd|%i>|^2=%6.3f |<%i|Qe|%i>|^2=%6.3f |<%i|Qf|%i>|^2=%6.3f",i,j,abs(v1(4))*abs(v1(4)),i,j,abs(v1(5))*abs(v1(5)),i,j,abs(v1(6))*abs(v1(6)));
         printf(" n%i=%6.3f\n",i,therm/Z);
      }
   }
         v1 *= sqrt(therm / Z);

 if (ninit>Hsz)ninit=Hsz;
 if (pinit<SMALL)pinit=SMALL;
   double zsum=0,zi;
   // determine number of thermally reachable states
 int noft=0;for(i=0;(i<ninit)&((zi=(exp(-(est[0][i+1].real()-est[0][1].real())/(KB*fabs(T)))))>(pinit*zsum));++i){noft+=Hsz-i-1;zsum+=zi;}
//   int noft=0;for(i=0;(i<Hsz)     &((exp(-(est[0][i+1].real()-est[0][1].real())/(KB*T)))>SMALL);++i)noft+=Hsz-i-1; // removed MR  6.9.2011 to allow for mcdisp options -ninit -pinit
   return noft;
//   return Hsz*(Hsz-1)/2;
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
      sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars); Hic/=MEV2CM; iHic/=MEV2CM; H = zmat2f(Hic,iHic);
//       est = ComplexMatrix(0,Hsz,0,Hsz); I comment this out - you should not reinitialize est !!!
      if( (est.Rhi()!=Hsz||est.Chi()!=Hsz) && pars.truncate_level==1)
      {
         std::cerr << "ERROR module ic1ion - Icalc: Hsz recalculation does not agree with eigenstates matrix dimension\n"; exit(EXIT_FAILURE);
      }
      if (!((int)(parval.size()/2)>Hsz && pars.truncate_level==1))
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
      iceig VE; if(pars.partial) VE.lcalc(pars,Jm); 
      #ifndef NO_ARPACK
      else if(pars.arnoldi) VE.acalc(pars,Jm); 
      #endif
      else VE.calc(Hsz,Jm); free(Jm);

      // Calculates the expectation values sum_n{ <n|Ja|n> exp(-En/kT) }
      std::vector< std::vector<double> > matel;
      std::vector<double> vJ = mfmat.spindensity_expJ(VE, xyz,*T,matel,pars.save_matrices);
      for(i=J.Lo(); i<=J.Hi(); i++) J[i] = vJ[i-J.Lo()];//printf("ss=%g\n",J(1));
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of chargedensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
void chargedensity_coeff(Vector &mom,          // Output single ion moments =expectation values of
                                           // of Zlm R^2(r) at a given temperature T and  effective field H
                      double *T,           // Input scalar temperature
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,    // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,    // Input vector of parameters from single ion property file
                      char **sipffilename, // Single ion properties filename
                      ComplexMatrix &est)  // Input/output eigenstate matrix (initialized in parstorage)
{Vector moments(1,51),ABC; double gJ=0.,lnZ,U;
 Vector Hxc(1,51);Hxc=0;for(int i=1;i<=gjmbHxc.Hi();++i){Hxc(i)=gjmbHxc(i);}
  Icalc(moments,T,Hxc,Hext,&gJ,ABC,sipffilename,&lnZ,&U,est);

   // Parses the input file for parameters
   icpars pars; 
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);
  
// a(0, 0) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// Indices for spindensity
//            0 not used
//            0 1  2  3 4 5 6  7  8  9 101112131415 16 17 18 19 20 2122232425262728 
//int k[] = {-1,0, 2, 2,2,2,2, 4, 4, 4, 4,4,4,4,4,4, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
//int q[] = {-1,0,-2,-1,0,1,2,-4,-3,-2,-1,0,1,2,3,4,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
mom(1) =  pars.n / sqrt(4.0 * 3.1415); // nofelectrons 
for(int i=2;i<=6;++i)mom(i)=moments(5+i)*sqrt((2.0*2+1)/8/PI);
mom(4)*=sqrt(2);
for(int i=7;i<=15;++i)mom(i)=moments(12+i)*sqrt((2.0*4+1)/8/PI);
mom(11)*=sqrt(2);
for(int i=16;i<=28;++i)mom(i)=moments(23+i)*sqrt((2.0*6+1)/8/PI);
mom(22)*=sqrt(2);

//   {for(l=2;l<=6;l+=2){for(m=-l;m<=l;++m){if(m!=0){a(l,m)*=sqrt((2.0*l+1)/8/PI);}else{a(l,m)*=sqrt((2.0*l+1)/4/PI);}}}
//         } // in case
           // of module ic1ion we just take the prefactors of the Zlm ... ??? what should we take here ???
           // MR 23.8.2011: if Tkq are define as in our review then the above should be right
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
                                           // of Zlm R^2(r) at a given temperature T and  effective field H
                      int & xyz,           // direction 1,2,3 = x,y,z
                      double *T,           // Input scalar temperature
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,    // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,    // Input vector of parameters from single ion property file
                      char **sipffilename, // Single ion properties filename
                      ComplexMatrix &est)  // Input/output eigenstate matrix (initialized in parstorage)
{// sum exchange field and external field
  Vector gjmbH(1,gjmbHxc.Hi());
  gjmbH=gjmbHxc;
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
                                           // of Zlm R^2(r) at a given temperature T and  effective field H
                      int & xyz,           // direction 1,2,3 = x,y,z
                      double *T,           // Input scalar temperature
                      Vector &gjmbHxc,      // Input vector of exchange fields (meV) 
                      Vector &Hext,      // Input vector of external field (meV) 
 /* Not Used */       double * /*g_J*/,    // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,    // Input vector of parameters from single ion property file
                      char **sipffilename, // Single ion properties filename
                      ComplexMatrix &est)  // Input/output eigenstate matrix (initialized in parstorage)
{// sum exchange field and external field
  Vector gjmbH(1,gjmbHxc.Hi());
  gjmbH=gjmbHxc;
   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(Hext(1))>DBL_EPSILON || fabs(Hext(2))>DBL_EPSILON || fabs(Hext(3))>DBL_EPSILON)
   {
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH(2)+=MUB*Hext(1); gjmbH(1)+=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH(4)+=MUB*Hext(2); gjmbH(3)+=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH(6)+=MUB*Hext(3); gjmbH(5)+=GS*MUB*Hext(3); }
   }

   sdod_Icalc(J,-xyz,T,gjmbH,sipffilename,est);
}
