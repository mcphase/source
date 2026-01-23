#include "ic1ion_module.hpp"
#include "truncate.cpp"

#if defined(__linux__) || defined(__APPLE__)
extern "C"
{
ic1ion_module *allocator(const char * filename)
 {return new ic1ion_module(filename);
 }
void deleter(ic1ion_module *ptr)
 {delete ptr;
 }
}
#endif
#ifdef __MINGW32__
extern "C"
{
__declspec (dllexport) ic1ion_module *allocator(const char * filename)
{
return new ic1ion_module(filename);
}
__declspec (dllexport) void deleter(ic1ion_module *ptr)
{
delete ptr;
}
}
#endif

ic1ion_module::ic1ion_module(const char * filename)
{pars=icpars();
 ic_parseinput(filename,pars);
 mfmat=icmfmat(pars.n,pars.l,6,pars.save_matrices);
 Sxmat=icmfmat(pars.n,pars.l,0,pars.save_matrices,1);
 Symat=icmfmat(pars.n,pars.l,0,pars.save_matrices,2);
 Szmat=icmfmat(pars.n,pars.l,0,pars.save_matrices,3);
 Lxmat=icmfmat(pars.n,pars.l,0,pars.save_matrices,-1);
 Lymat=icmfmat(pars.n,pars.l,0,pars.save_matrices,-2);
 Lzmat=icmfmat(pars.n,pars.l,0,pars.save_matrices,-3);
 Hic = ic_hmltn(iHic,pars); Hic/=MEV2CM; iHic/=MEV2CM;
}

ic1ion_module::ic1ion_module(const ic1ion_module & pp)
{printf("copying ic1ion_module\n");
pars=pp.pars;
mfmat=pp.mfmat;
Sxmat=pp.Sxmat;
Symat=pp.Symat;
Szmat=pp.Szmat;
Lxmat=pp.Lxmat;
Lymat=pp.Lymat;
Lzmat=pp.Lzmat;
Hic=pp.Hic;iHic=pp.iHic;
}


void myPrintMatrix(FILE * file,sMat<double> & M,int d)
{ int i1,j1;
   fprintf (file,"Matrix\n");
   for (i1=0;i1<=d;++i1)
   {
      for (j1=0;j1<=d;++j1) fprintf (file,"%6.3f ",M(i1,j1));
      fprintf (file,"\n");
   }
}    

void zmat2pack(sMat<double> &r, sMat<double> &i, Matrix &outmat)
{ sMat<double> tmp = r+i;
   std::vector< std::vector<int> > u = tmp.findlower();
   // Allocates an _r*_c array and initiallises all elements to zero.
   Matrix retval(1,tmp.nr(),1,tmp.nc()); retval=0;
   for (int j=0; j<(int)u.size(); j++)
   {
      retval(u[j][0]+1,u[j][1]+1) = r(u[j][0],u[j][1]);
      retval(u[j][1]+1,u[j][0]+1) = i(u[j][0],u[j][1]);
   }
   // Diagonal elements
   for (int j=0; j<tmp.nr(); j++)
      retval(j+1,j+1) = r(j,j);
   outmat = retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Checks whether the Matpack matrix is the same as the c-array
/* --------------------------------------------------------------------------------------------------------------- //
bool checkmat(ComplexMatrix &cmat, complexdouble *fmat,int r, int c)
{ int i,j;
   for(i=0; i<(cmat.Rows()-r); i++)
      for(j=0; j<(cmat.Cols()-c); j++)
      { 
       //std::cout << "cmat["<<i+r<<"]["<<j+c<<"]=" << cmat[i+1][j+1] << "\tfmat["<<j<<"*" << cmat.Cols()-c << "+"<<i<<"]=";
       //   std::cout << fmat[j*(cmat.Cols()-c)+i].r << "+" << fmat[j*(cmat.Cols()-c)+i].i << "i\n";
       //std::cout << "cmat["<<j+c<<"]["<<i+r<<"]=" << cmat[j+c][i+r] << "\tfmat["<<j<<"*" << cmat.Cols()-c << "+"<<i<<"]=";
       //   std::cout << fmat[j*(cmat.Cols()-c)+i].r << "+" << fmat[j*(cmat.Cols()-c)+i].i << "i\n";
       //if(real(cmat[i+r][j+c])!=fmat[j*(cmat.Cols()-c)+i].r || imag(cmat[i+r][j+c])!=fmat[j*(cmat.Cols()-c)+i].i) return false; }
         if(real(cmat[j+c][i+r])!=fmat[j*(cmat.Cols()-c)+i].r || imag(cmat[j+c][i+r])!=fmat[j*(cmat.Cols()-c)+i].i) return false; }
   return true;
}
*/

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate <Ia><Ib>... at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
bool ic1ion_module::IMcalc(Matrix &Jret,          // Output field of single ion momentum vector <Ja>,<Jb>,<Jc>, etc.
                      Vector&T,          // Input Vector of temperatures
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double & g_J,   // Input Lande g-factor
 /* Not Used */       Vector & ABC,   // Input vector of parameters from single ion property file
                      char *sipffilename,// Single ion properties filename
                      Vector &lnZ,        // Output scalar logarithm of partition function
                      Vector &U)          // Output scalar internal energy 
{ expJ(mfmat,Jret,T,Hxc,Hext,lnZ,U);
return true;
}

bool ic1ion_module::Icalc(Vector &Jret,          // Output single ion momentum vector <Ja>,<Jb>,<Jc>, etc.
                      double &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,   // Input Lande g-factor
 /* Not Used */       Vector & ABC,   // Input vector of parameters from single ion property file
                      char *sipffilename,// Single ion properties filename
                      double &lnZ,        // Output scalar logarithm of partition function
                      double &U)         // Output scalar internal energy 
{ expJ(mfmat,Jret,T,Hxc,Hext,lnZ,U);
  return true;
}

void ic1ion_module::expJ(icmfmat & mfm,     // Operators to be calculated
                       Vector &Jret,          // Output single ion momentum vector <Ja>,<Jb>,<Jc>, etc.
                      double &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
                      double &lnZ,        // Output scalar logarithm of partition function
                      double &U)          // Output scalar internal energy 
{ Matrix JM(1,Jret.Hi(),1,1);
 lnZ = 0.;
 U = 0.;
 for(int i=1;i<=Jret.Hi();++i)JM(i,1)=Jret(i);
 Vector TT(1,1);TT(1)=T;
 Vector lnZZ(1,1);lnZZ(1)=lnZ;
 Vector UU(1,1);UU(1)=U;
 expJ(mfm,JM,TT,Hxc,Hext,lnZZ,UU);
 U=UU(1);lnZ=lnZZ(1);T=TT(1);
 for(int i=1;i<=Jret.Hi();++i)
 {Jret(i)=JM(i,1);//printf("Jret(%i)=%g ",i,Jret(i));
 }
}


void ic1ion_module::expJ(icmfmat & mfmOP,     // Operators to be calculated
                      Matrix &Jret,       // Output single ion momentum vector <Ja>,<Jb>,<Jc>, etc.
                      Vector&T,           // Input Vector of temperatures
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
                      Vector &lnZ,        // Output scalar logarithm of partition function
                      Vector &U)          // Output scalar internal energy 
 { Vector gjmbH=sum_Hxc_Hext(Hxc,Hext);
   Matrix J(1,Jret.Rhi(),1,T.Hi()); // matrix for output to be written to Jret
   for(int Ti=1;Ti<=T.Hi();++Ti)for(int i=1;i<=J.Rhi();++i)J(i,Ti)=0;
   // Converts the Jij parameters if necessary
   std::vector<double> vgjmbH((gjmbH.Hi()-gjmbH.Lo()+1),0.); 
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc();
      for(int i=gjmbH.Lo(); i<=gjmbH.Hi(); i++) vgjmbH[i-gjmbH.Lo()] = -gjmbH[i]*pars.jijconv[i]; }
   else
   #endif
      for(int i=gjmbH.Lo(); i<=gjmbH.Hi(); i++) vgjmbH[i-gjmbH.Lo()] = -gjmbH[i];  // Vector of exchange + external fields to be added to matrix below

   // Calculates the IC Hamiltonian matrix
   int i,Hsz=getdim(pars.n,pars.l);
   complexdouble *Jm=0;
   if(pars.truncate_level!=1)     // Uses the eigenvectors of the single ion Hamiltonian to truncate the matrix
   {// check if truncate is to be used and if Hamiltonian was already calculated 
      if(mfmat.T[0]==NULL){truncate_hmltn(pars,  Hic, iHic, J.Rhi(), J.Rlo());}
    truncate_expJ(pars,gjmbH,J,T,lnZ,U);
   }
   else
   {  // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
      #ifdef JIJCONV
      if(pars.B.norm().find("Stevens")!=std::string::npos) mfmat.jijconv.assign(pars.jijconv.begin(),pars.jijconv.end());
      #endif
      sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH); // add J.H to matrix
      Jm = zmat2f(Jmat,iJmat);  Hic.addto(Jm,false);if(!iHic.isempty())iHic.addto(Jm,true);

      // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
      iceig VE; if(pars.partial) VE.lcalc(pars,Jm);
      #ifndef NO_ARPACK
      else if(pars.arnoldi) VE.acalc(pars,Jm); 
      #endif
      else {VE.calc(Hsz,Jm);} free(Jm);
      // Calculates the expectation values sum_n{ <n|Ja|n> exp(-En/kT) }
        std::vector< std::vector<double> > matel; 
      // get expJ to highest T and matrix elements of eigenstates matel 
      // (for number of low energy states necessary for calculation at Ti=T.Hi() )
      int Ti=T.Hi();

      std::vector<double> vJ =  mfmOP.expJ(VE,T(Ti),matel,J.Rhi());
   
      for(i=J.Rlo(); i<=J.Rhi(); i++) {J(i,T.Hi()) = vJ[i-J.Rlo()];//printf("%g ",J(i,Ti));
                                         } 
      lnZ(T.Hi())=vJ[J.Rhi()-J.Rlo()+1];
      U(T.Hi())=vJ[J.Rhi()-J.Rlo()+2];
      
      vector<double> E; // energy vector
      if(T.Hi()>1)for(int ind_j=0; ind_j<(int)matel[0].size(); ind_j++){E.push_back(VE.E(ind_j)-VE.E(0));}
    // use matel to calculate more quickly the other temperatures
    for(Ti=1;Ti<T.Hi();++Ti)
    {int Esz;std::vector<double> eb;
     Esz=matel[0].size();
     if (T(Ti)<0){Esz=(int)(-T(Ti));printf ("Temperature T=%g<0: please choose probability distribution for the -T=%i lowest energy states by hand\n",T(Ti),(int)(-T(Ti)));
                         printf ("Number   Excitation Energy\n");
     for (int ind_j=0;ind_j<Esz;++ind_j) printf ("%i    %4.4g meV\n",ind_j+1,E[ind_j]);
     } 
     U(Ti)=0;double Z=0;eb.assign(Esz,0.);
     for(int iJ=0; iJ<J.Rhi(); iJ++)
     {J(iJ+1,Ti)=0;
      for(int ind_j=0; ind_j<Esz; ind_j++)
      {if(iJ==0) // for iJ==0 sum up also U and Z
        { if (T(Ti)<0)
         {  char instr[MAXNOFCHARINLINE];
            printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
            if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { printf("Error in input. Exiting\n"); exit(-1); }
            eb[ind_j]=strtod(instr,NULL);
         }
         else
         { 
         eb[ind_j] = exp(-E[ind_j]/(KB*T(Ti)));
         } 
         Z+=eb[ind_j]; 
         U(Ti)+=(VE.E(ind_j))*eb[ind_j];
        } 
        J(iJ+1,Ti)+=matel[iJ][ind_j]*eb[ind_j];
      }       
      J(iJ+1,Ti)/=Z; 
      if(fabs(J(iJ+1,Ti))<DBL_EPSILON) J(iJ+1,Ti)=0.; 
     } // iJ
    lnZ(Ti) = log(Z)-VE.E(0)/(KB*T(Ti)); // set lnZ
    U(Ti)/=Z;
    } //  Ti   
   } // fi truncate

  for(int Ti=1;Ti<=T.Hi();++Ti)
   {
   for(i=Jret.Rlo();i<=Jret.Rhi();++i)Jret(i,Ti)=J(i,Ti);
   }

}



// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the magnetic moment at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
bool ic1ion_module::mcalc(Vector &mom,        // Output magnetic moment (mub)
                      double &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,   // Input Lande g-factor
 /* Not Used */       Vector & ABC,   // Input vector of parameters from single ion property file
                      char *sipffilename)// Single ion properties filename
 {  Vector J(1,6); 
   double  lnZ, U;
   expJ(mfmat,J,T,Hxc,Hext,lnZ,U);
   mom(1)=GS*J(1)+J(4);
   mom(2)=GS*J(2)+J(5);
   mom(3)=GS*J(3)+J(6);
return true;
}

bool ic1ion_module::mMcalc(Matrix &mom,        // Output magnetic moment (mub)
                      Vector &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,   // Input Lande g-factor
 /* Not Used */       Vector & ABC,   // Input vector of parameters from single ion property file
                      char *sipffilename)// Single ion properties filename
{  Matrix J(1,6,1,T.Hi()); 
   Vector lnZ(1,T.Hi()), U(1,T.Hi());
   expJ(mfmat,J,T,Hxc,Hext,lnZ,U);
   for(int Ti=1;Ti<=T.Hi();++Ti){
   mom(1,Ti)=GS*J(1,Ti)+J(4,Ti);
   mom(2,Ti)=GS*J(2,Ti)+J(5,Ti);
   mom(3,Ti)=GS*J(3,Ti)+J(6,Ti); }
 return true;
}
// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the <L> at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
bool ic1ion_module::Lcalc(Vector &L,          // Output magnetic moment (mub)
                      double &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,   // Input Lande g-factor
 /* Not Used */       Vector & ABC,   // Input vector of parameters from single ion property file
                      char *sipffilename)// Single ion properties filename
{  Vector J(1,6); 
   double  lnZ, U;
   expJ(mfmat,J,T,Hxc,Hext,lnZ,U);
   L(1)=J(4);
   L(2)=J(5);
   L(3)=J(6);
return true;
}

bool ic1ion_module::LMcalc(Matrix &L,          // Output magnetic moment (mub)
                      Vector &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,   // Input Lande g-factor
 /* Not Used */       Vector & ABC,   // Input vector of parameters from single ion property file
                      char *sipffilename)// Single ion properties filename
{  Matrix J(1,6,1,T.Hi()); 
  Vector lnZ(1,T.Hi()), U(1,T.Hi());
   expJ(mfmat,J,T,Hxc,Hext,lnZ,U);
   for(int Ti=1;Ti<=T.Hi();++Ti){
   L(1,Ti)=J(4,Ti);
   L(2,Ti)=J(5,Ti);
   L(3,Ti)=J(6,Ti);}
 return true;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the <S> at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
bool ic1ion_module::Scalc(Vector &S,          // Output magnetic moment (mub)
                      double &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,   // Input Lande g-factor
 /* Not Used */       Vector & ABC,   // Input vector of parameters from single ion property file
                      char *sipffilename)// Single ion properties filename
{  Vector J(1,6); 
   double  lnZ, U;
   expJ(mfmat,J,T,Hxc,Hext,lnZ,U);
   S(1)=J(1);
   S(2)=J(2);
   S(3)=J(3);
 return true;
}

bool ic1ion_module::SMcalc(Matrix &S,          // Output magnetic moment (mub)
                      Vector &T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,   // Input Lande g-factor
 /* Not Used */       Vector & ABC,   // Input vector of parameters from single ion property file
                      char *sipffilename)// Single ion properties filename
{  Matrix J(1,6,1,T.Hi()); 
   Vector lnZ(1,T.Hi()), U(1,T.Hi());
   expJ(mfmat,J,T,Hxc,Hext,lnZ,U);
   for(int Ti=1;Ti<=T.Hi();++Ti){
   S(1,Ti)=J(1,Ti);
   S(2,Ti)=J(2,Ti);
   S(3,Ti)=J(3,Ti);}
 return true;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements
// --------------------------------------------------------------------------------------------------------------- //
int ic1ion_module::du1calc(int &tn,            // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double & g_J,    // Input Lande g-factor
 /* Not Used */       Vector & ABC,    // Input vector of parameters from single ion property file
                      char *sipffilename,// Single ion properties filename
                      ComplexVector & u1ret, // Output u1 vector
                      float &delta,       // Input maximal Energy / Output transition energy
                      int &n, int &nd,    // Output state numbers of initial and final state
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{  Vector gjmbH=sum_Hxc_Hext(Hxc,Hext);
   ComplexVector u1(1,gjmbH.Hi()); u1=0; u1(1)=u1ret(1);
   int Hsz = est.Rows()-1;
   int noft=mfmat.u1((complexdouble*)&u1[1],u1.Hi(),T,tn,delta,(complexdouble*)&est[1][0],(complexdouble*)&est[0][1],Hsz,n,nd);
   for(int i=1; i<=u1ret.Hi(); i++) u1ret[i] = u1[i];
   return noft;   
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements of magnetic moment operator
// --------------------------------------------------------------------------------------------------------------- //
int ic1ion_module::dm1(int &tn,            // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char *sipffilename,// Single ion properties filename
                      ComplexVector & m1, // Output m1 vector (1,3)
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{  ComplexVector u1(1,6);int n,nd;
   u1(1) = m1(1);
   int nt = du1calc(tn,T,Hxc,Hext,g_J,ABC,sipffilename,u1,delta,n,nd,est);
   m1(1)=GS*u1(1)+u1(4);
   m1(2)=GS*u1(2)+u1(5);
   m1(3)=GS*u1(3)+u1(6);
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements of orbital angular momentum operator
// --------------------------------------------------------------------------------------------------------------- //
int ic1ion_module::dL1(int &tn,            // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char *sipffilename,// Single ion properties filename
                      ComplexVector & L1, // Output L1 vector (1,3)
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{  ComplexVector u1(1,6);int n,nd;
   u1(1) = L1(1);
   int nt=du1calc(tn,T,Hxc,Hext,g_J,ABC,sipffilename,u1,delta,n,nd,est);
   L1(1)=u1(4);
   L1(2)=u1(5);
   L1(3)=u1(6);
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements of spin operator
// --------------------------------------------------------------------------------------------------------------- //
int ic1ion_module::dS1(int &tn,            // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char *sipffilename,// Single ion properties filename
                      ComplexVector & S1, // Output S1 vector (1,3)
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{  ComplexVector u1(1,6);int n,nd;
   u1(1) = S1(1);
   int nt=du1calc(tn,T,Hxc,Hext,g_J,ABC,sipffilename,u1,delta,n,nd,est);
   S1(1)=u1(1);
   S1(2)=u1(2);
   S1(3)=u1(3);
   return nt;
}



// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the eigenstates of Hic+effective_mean_fields
// --------------------------------------------------------------------------------------------------------------- //
bool ic1ion_module::estates(ComplexMatrix &est, // Output Eigenstates matrix (row 0: real==Eigenvalues;imag==population)
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,   // Input  Lande g-factor
                      double &T,          // Input  temperature
 /* Not Used */       Vector & ABC,   // Input  Vector of parameters from single ion property file
                      char *sipffilename)// Input  Single ion properties filename
{ Vector gjmbH=sum_Hxc_Hext(Hxc,Hext);

   clock_t start,end; start = clock();
   // Calculates the IC Hamiltonian matrix
    int Hsz = Hic.nr();
   // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
    int i,j,gLo=gjmbH.Lo(),gHi=gjmbH.Hi(); std::vector<double> vgjmbH(gHi,0.);
   for(i=gLo; i<=gHi; i++) vgjmbH[i-1] = -gjmbH[i];
   // Converts the Jij parameters if necessary
   #ifdef JIJCONV
   if(pars.B.norm().find("Stevens")!=std::string::npos) {
      pars.jijconvcalc(); mfmat.jijconv.assign(pars.jijconv.begin(),pars.jijconv.end());
      for(i=gLo; i<=gHi; i++) vgjmbH[i-1] *= pars.jijconv[i]; }
   #endif
   sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH); 

   // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
    Jmat+=Hic; if(!iHic.isempty()) iJmat+=iHic; 
   iceig VE; if(iJmat.isempty()) VE.calc(Jmat); else VE.calc(Jmat,iJmat);

   // Initialises the output matrix
   est = ComplexMatrix(0,Hsz,0,Hsz);

   // Stores the number of electrons and the orbital number in element (0,0)
   est(0,0) = complex<double> (pars.n, pars.l);

   // Puts eigenvectors/values into the est matrix
   for(i=0; i<Hsz; i++) est(0,i+1) = complex<double> (VE.E(i), exp(-(VE.E(i)-VE.E(0))/(KB*T)));   // Row 0

   if(VE.iscomplex()) {
      for(i=1; i<=Hsz; i++) memcpy(&est[i][1],VE.zV(i-1),Hsz*sizeof(complexdouble));  
//    std::cout << "\n\nest==VE.zV = " << checkmat((*est),VE.zV(0),1,1) << endl << endl; 
      }
   else
      for(i=0; i<Hsz; i++)
      {  //printf("\n");
         for(j=0; j<Hsz; j++) 
         {
            est(i+1,j+1) = complex<double> (VE.V(i)[j], 0.);
//          printf("%6.3f %+6.3f i  ",VE.V(i)[j],0.0);
            if(VE.V(i)[j]!=VE.V(j,i)){fprintf(stderr,"compiler problem: bad memory mapping of vectors\n");exit(EXIT_FAILURE);}
            if(VE.V(i)[j]!=est[i+1][j+1]){fprintf(stderr,"compiler problem: bad memory mapping of vectors\n");exit(EXIT_FAILURE);}
         }
      }
  end = clock(); std::cerr << "#Time to do estates() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
return true;
}

// --------------------------------------------------------------------------------------------------------------- //
// Loads a Q_q matrix from file if the file exists and has the same parameters n,l,Jvec
// --------------------------------------------------------------------------------------------------------------- //
bool get_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l, std::vector<double> &Jvec)
{  int i, j, mn, r, c, sz, ml;
   std::vector<double> mJv(6,0.); 
   char filename[] = "results/mcphas.Qq"; filename[16]=q+120;               // 120==x, 121==y, 122==z
   std::fstream FILEIN; FILEIN.open(filename, std::fstream::in);
   if(FILEIN.fail()==true) return false;
   FILEIN >> mn >> ml; for(i=0; i<6; i++) FILEIN >> mJv[i]; FILEIN >> r >> c;
   if(mn!=n || ml!=(int)l) {return false;} for(i=0; i<6; i++) if(fabs(mJv[i]-Jvec[i])>1e-4) return false;
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
{  int i,j,sz;
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
bool ic1ion_module::mqcalc(ComplexVector &Mq,      // Output expectation values -2[<Q>_{x} <Q>_y <Q>_z]
                  double &th, double &ph, // Input polar and azimuth angles theta and phi
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est)     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
{  int i,q,n=1,Hsz=est.Cols()-1; orbital l;
   n = (int)est[0][0].real(); i = (int)est[0][0].imag(); l = (orbital)i;
   if(i>3 || i<0) { std::cerr << "ic1ion mqcalc(): Error only s-, p-, d-, and f-electrons supported.\n"; exit(EXIT_FAILURE); }
   std::vector<double> E,Jvec(6,0.); Jvec[0]=th; Jvec[1]=ph; Jvec[2]=J0; Jvec[3]=J2; Jvec[4]=J4; Jvec[5]=J6;
   std::vector< sMat<double> > Qp, Qm; 
   std::vector< std::vector< sMat<double> > > Qmat; for(i=0; i<3; i++) Qmat.push_back(Qp);
   complexdouble  zme;
   double zMqr,zMqi,Z=0.;
   //char trans = 'U'; int incx=1;

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
    //myPrintMatrix(stdout,Qmat[2][0],Hsz-1);
   }
   for(q=0; q<3; q++)
   {
       zMqr = 0.; zMqi = 0.;
      for(i=1; i<=Hsz; i++)
      {if(q==1){  zme.i=Qmat[q][0].MultvxMv((complexdouble*)&est[i][1],false);
                  zme.r=Qmat[q][1].MultvxMv((complexdouble*)&est[i][1],true);
        }else{ zme.r=Qmat[q][0].MultvxMv((complexdouble*)&est[i][1],false);
               zme.i=Qmat[q][1].MultvxMv((complexdouble*)&est[i][1],true);}
      // zme.r=Qmat[q][0].MultvxMv((complexdouble*)&est[i][1],false);
      // zme.i=Qmat[q][1].MultvxMv((complexdouble*)&est[i][1],true);
//       printf ("%i zme=%g %+g i  Ei=%6.3f ni=%6.3f \n",i,zme.r,zme.i,est[0][i].real(),est[0][i].imag());
         zMqr += (-2.)*zme.r*est[0][i].imag(); zMqi += (-2.)*zme.i*est[0][i].imag(); if(q==0) Z += est[0][i].imag();
      }
      Mq[q+1] = complex<double> (zMqr, zMqi)/Z;
   }
// printf("MQ=(%g %+g i, %g %+g i,%g %+g i)\n",real(Mq(1)),imag(Mq(1)),real(Mq(2)),imag(Mq(2)),real(Mq(3)),imag(Mq(3)));
return true;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the transition matrix using the scattering operator of Balcar and Lovesey
// --------------------------------------------------------------------------------------------------------------- //
int ic1ion_module::dmq1(int &tn,                // Input transition number |tn|. If tn>0 omit printout. If tn<0 print info.
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
   n = (int)est[0][0].real(); i = (int)est[0][0].imag(); l = (i==2) ? D : F;
   std::vector<double> E,Jvec(6,0.); Jvec[0]=th; Jvec[1]=ph; Jvec[2]=J0; Jvec[3]=J2; Jvec[4]=J4; Jvec[5]=J6;
   std::vector< sMat<double> > Qp1, Qm1; 
   std::vector< std::vector< sMat<double> > > Qq; for(i=0; i<3; i++) Qq.push_back(Qp1);
   complexdouble z1,z2, zbeta;  zbeta.r=0; zbeta.i=0;
   std::vector<complexdouble> zij(7,zbeta), zji(7,zbeta);
   double Z=0., therm;
   
   // Calculates the scattering operator, Q.
   if(!get_Qq(Qq[0],0,n,l,Jvec) || !get_Qq(Qq[1],1,n,l,Jvec))                  // Qq[0]==Qx, Qq[1]==Qy, Qq[2]==Qz
   {
      lovesey_Qq(Qm1,-1,n,l,Jvec); lovesey_Qq(Qp1,1,n,l,Jvec);
      for(i=0; i<6; i++)  
      {
         Qq[0].push_back( (Qp1[i]-Qm1[i]) * (-1.0/sqrt(2.0)) );                // Qx = -1/sqrt(2) * (Q_{+1} - Q_{-1})
         if(i%2==0) Qq[1].push_back( (Qp1[i+1]+Qm1[i+1]) * (-1.0/sqrt(2.0)) ); // real(Qy) = i^2/sqrt(2) * imag(Q_{+1}+Q_{-1})
         else       Qq[1].push_back( (Qp1[i-1]+Qm1[i-1]) *  (1.0/sqrt(2.0)) ); // imag(Qy) = i/sqrt(2) * real(Q_{+1}+Q_{-1})
      }
      save_Qq(Qq[0],0,n,l,Jvec); save_Qq(Qq[1],1,n,l,Jvec); 
    //for(int iQ=2; iQ<6; iQ++) myPrintMatrix(stdout,Qq[0][iQ],Hsz-1);
   }
   if(!get_Qq(Qq[2],2,n,l,Jvec))                                               // Loads the Q_q matrix if prev. saved
   {
      lovesey_Qq(Qq[2],0,n,l,Jvec);                                            // Calcs. scattering operator matrix Q_q
      save_Qq(Qq[2],2,n,l,Jvec);
   }

   for(i=0; i<Hsz; i++) { 
      therm = exp(-(est[0][i+1].real()-est[0][1].real())/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }

   int a,j=0,k=0; for(i=0; i<Hsz; ++i) {for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   ++i;++j; // because in est i and j start from 1...Hsz
 
   for(q=0; q<3; q++)
   {  z1=Qq[q][2].MultuxMv((complexdouble*)&est[i][1],(complexdouble*)&est[j][1],false);
      z2=Qq[q][3].MultuxMv((complexdouble*)&est[i][1],(complexdouble*)&est[j][1],true);
      zij[2*q+1].r=z1.r+z2.r;
      zij[2*q+1].i=z1.i+z2.i;
      z1=Qq[q][2].MultuxMv((complexdouble*)&est[j][1],(complexdouble*)&est[i][1],false);
      z2=Qq[q][3].MultuxMv((complexdouble*)&est[j][1],(complexdouble*)&est[i][1],true);
      zji[2*q+1].r=z1.r+z2.r;
      zji[2*q+1].i=z1.i+z2.i;
      if(i==j)                               //subtract thermal expectation value from zij=zii
      {                                      //MR120120 ... reintroduced
         complexdouble expQ;double thexp=0;
         for(iJ=1;iJ<=Hsz;++iJ)
         {
            therm = exp(-(est[0][iJ].real()-est[0][1].real())/(KB*T)); if(therm<DBL_EPSILON) break;
            if(q==1)
            expQ.r=Qq[q][3].MultvxMv((complexdouble*)&est[iJ][1],true);
            else
            expQ.r=Qq[q][2].MultvxMv((complexdouble*)&est[iJ][1],false);
         
           thexp += expQ.r * therm / Z;
         }
         zij[2*q+1].r-=thexp;zji[2*q+1].r-=thexp;
      }
      
      z1=Qq[q][4].MultuxMv((complexdouble*)&est[i][1],(complexdouble*)&est[j][1],false);
      z2=Qq[q][5].MultuxMv((complexdouble*)&est[i][1],(complexdouble*)&est[j][1],true);
      zij[2*q+2].r=z1.r+z2.r;
      zij[2*q+2].i=z1.i+z2.i;
      z1=Qq[q][4].MultuxMv((complexdouble*)&est[j][1],(complexdouble*)&est[i][1],false);
      z2=Qq[q][5].MultuxMv((complexdouble*)&est[j][1],(complexdouble*)&est[i][1],true);
      zji[2*q+2].r=z1.r+z2.r;
      zji[2*q+2].i=z1.i+z2.i;
      if(i==j)                               //subtract thermal expectation value from zij=zii
      {                                      //MR120120 ... reintroduced
         complexdouble expQ;double thexp=0;
         for(iJ=1;iJ<=Hsz;++iJ)
         {
            therm = exp(-(est[0][iJ].real()-est[0][1].real())/(KB*T)); if(therm<DBL_EPSILON) break;
            if(q==1)
            expQ.r=Qq[q][5].MultvxMv((complexdouble*)&est[iJ][1],true);
            else
            expQ.r=Qq[q][4].MultvxMv((complexdouble*)&est[iJ][1],false);
            thexp += expQ.r * therm / Z;
         }
         zij[2*q+2].r-=thexp;zji[2*q+2].r-=thexp;
      }
   }

   // check if zij are complex conjugate
   for(iJ=1;iJ<=6;++iJ)
      if(fabs(zij[iJ].i+zji[iJ].i)>SMALL) { std::cerr << "ERROR module ic1ion - dmq1: <i|Qalpha|j>not hermitian\n"; exit(EXIT_FAILURE); }
                
   //complex<double> im(0,1);
   ComplexVector iQalphaj(1,6);
   
   for(a=1; a<=6; a++){iQalphaj(a) = complex<double> (zij[a].r,zij[a].i);if(a%2==1){iQalphaj(a)*=0.5;}} 
                                                                         // divide spin part by 2
   mq1 = 0;
   for(a=1; a<=3; a++)
         mq1(a) =-2.0*(2.0*iQalphaj(-1+2*a)+iQalphaj(2*a));

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
   for(i=0; (i<ninit)&&((((x=(est[0][i+1].real()-est[0][1].real())/(KB*fabs(T)))<200)? zi=exp(-x):zi=0)>=(pinit*zsum)); ++i)
   {
      noft += Hsz-i; 
      zsum += zi;
   }
   return noft;
}


// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of chargedensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
bool ic1ion_module::chargedensity_coeff(
                      Vector &mom,         // Output single ion moments == expectation values of
                                           //    of Zlm R^2(r) at a given temperature T and  effective field H
                      double &T,           // Input scalar temperature
                      Vector &Hxc,         // Input vector of exchange fields (meV) 
                      Vector &Hext,        // Input vector of external field (T) 
 /* Not Used */       double &g_J,    // Input Lande g-factor
 /* Not Used */       Vector & ABC,    // Input vector of parameters from single ion property file
                      char *sipffilename) // Single ion properties filename
{
   Vector moments(1,51); 
   double lnZ, U;
   Vector Hxce(1,51); 
   Hxce=0;
   for(int i=1; i<=Hxc.Hi(); ++i) { Hxce(i)=Hxc(i); }
   expJ(mfmat,moments,T,Hxc,Hext,lnZ,U);
   
// a(0, 0) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// Indices for spindensity
//             0 not used
//             0 1  2  3 4 5 6  7  8  9 101112131415 16 17 18 19 20 2122232425262728 
// int k[] = {-1,0, 2, 2,2,2,2, 4, 4, 4, 4,4,4,4,4,4, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
// int q[] = {-1,0,-2,-1,0,1,2,-4,-3,-2,-1,0,1,2,3,4,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   mom(1) =  pars.n / sqrt(4.0 * 3.1415); // nofelectrons 
   for(int i=2; i<=6; ++i) {mom(i) = moments(5+i) *sqrt((2.0*2+1)/8/PI);} mom(4)  *= sqrt(2);
   for(int i=7; i<=15;++i) {mom(i) = moments(12+i)*sqrt((2.0*4+1)/8/PI);} mom(11) *= sqrt(2);
   for(int i=16;i<=28;++i) {mom(i) = moments(23+i)*sqrt((2.0*6+1)/8/PI);} mom(22) *= sqrt(2);

return true;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements of chargedensity coefficient operator
// --------------------------------------------------------------------------------------------------------------- //
int ic1ion_module::dchargedensity_coeff1(int &tn,        // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double & g_J,    // Input Lande g-factor
 /* Not Used */       Vector &ABC,    // Input vector of parameters from single ion property file
                      char *sipffilename,// Single ion properties filename
                      ComplexVector & dc1,// Output m1 vector (1,3)
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ 
   ComplexVector u1(1,51);int n,nd;
   double gJ=0.;
   u1(1) = dc1(1);
   Vector Hxce(1,51);
   Hxce = 0;
   for(int i=1; i<=Hxc.Hi(); ++i) { Hxce(i)=Hxc(i); }
   int nt = du1calc(tn,T,Hxce,Hext, gJ,ABC,sipffilename,u1,delta,n,nd,est);
   dc1(1)=0;
   for(int i=2; i<=6; ++i) {dc1(i) = u1(5+i) *sqrt((2.0*2+1)/8/PI);} dc1(4) *=sqrt(2);
   for(int i=7; i<=15;++i) {dc1(i) = u1(12+i)*sqrt((2.0*4+1)/8/PI);} dc1(11)*=sqrt(2);
   for(int i=16;i<=28;++i) {dc1(i) = u1(23+i)*sqrt((2.0*6+1)/8/PI);} dc1(22)*=sqrt(2);
   return nt;
}


// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of spindensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
bool ic1ion_module::spindensity_coeff(Vector &J,          // Output single ion moments =expectation values of
                                           //    of Zlm R^2(r) at a given temperature T and  effective field H
                      int & xyz,           // direction 1,2,3 = x,y,z
                      double &T,           // Input scalar temperature
                      Vector &Hxc,         // Input vector of exchange fields (meV) 
                      Vector &Hext,        // Input vector of external field (T) 
 /* Not Used */       double &g_J,    // Input Lande g-factor
 /* Not Used */       Vector & ABC,    // Input vector of parameters from single ion property file
                      char *sipffilename) // Single ion properties filename
{  Vector gjmbH=sum_Hxc_Hext(Hxc,Hext);
   if(pars.truncate_level!=1)     // Uses the eigenvectors of the single ion Hamiltonian to truncate the matrix
     {// check if truncate is to be used and if Hamiltonian was already calculated 
     if(mfmat.T[0]==NULL){    truncate_hmltn(pars,  Hic, iHic, J.Hi(), J.Lo());}
      truncate_spindensity_expJ(pars,gjmbH,J,T,xyz);
     }	
   else

   {double lnZ,U;
   switch(xyz)
   {case 1: expJ(Sxmat,J,T,Hxc,Hext,lnZ,U);break;
    case 2: expJ(Symat,J,T,Hxc,Hext,lnZ,U);break;
    case 3: expJ(Szmat,J,T,Hxc,Hext,lnZ,U);break;
   }
   }
return true;

}
 // sum exchange field and external field
 Vector ic1ion_module::sum_Hxc_Hext(Vector & Hxc,Vector & Hext)
  { Vector gjmbH(1,(Hxc.Hi()<6) ? 6 : Hxc.Hi()); gjmbH=0;
   if(gjmbH.Hi()==Hxc.Hi()) gjmbH=Hxc; 
   else for(int i=1; i<=(gjmbH.Hi()<Hxc.Hi()?gjmbH.Hi():Hxc.Hi()); i++) gjmbH[i]=Hxc[i];
   
   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(Hext(1))>DBL_EPSILON || fabs(Hext(2))>DBL_EPSILON || fabs(Hext(3))>DBL_EPSILON)
   {
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH(4)+=MUB*Hext(1); gjmbH(1)+=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH(5)+=MUB*Hext(2); gjmbH(2)+=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH(6)+=MUB*Hext(3); gjmbH(3)+=GS*MUB*Hext(3); }
   } 
   return gjmbH;
  }
// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the coefficients of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
bool ic1ion_module::orbmomdensity_coeff(Vector &J,        // Output single ion moments =expectation values of
                                           //    of Zlm R^2(r) at a given temperature T and  effective field H
                      int & xyz,           // direction 1,2,3 = x,y,z
                      double &T,           // Input scalar temperature
                      Vector &Hxc,         // Input vector of exchange fields (meV) 
                      Vector &Hext,        // Input vector of external field (T) 
 /* Not Used */       double &g_J,    // Input Lande g-factor
 /* Not Used */       Vector & ABC,    // Input vector of parameters from single ion property file
                      char *sipffilename) // Single ion properties filename
{ Vector gjmbH=sum_Hxc_Hext(Hxc,Hext);

  if(pars.truncate_level!=1)     // Uses the eigenvectors of the single ion Hamiltonian to truncate the matrix
     {// check if truncate is to be used and if Hamiltonian was already calculated 
     if(mfmat.T[0]==NULL){    truncate_hmltn(pars,  Hic, iHic, J.Hi(), J.Lo());}
      truncate_spindensity_expJ(pars,gjmbH,J,T,-xyz);
     }	
   else
   {
 double lnZ,U;
   switch(xyz)
   {case 1: expJ(Lxmat,J,T,Hxc,Hext,lnZ,U);break;
    case 2: expJ(Lymat,J,T,Hxc,Hext,lnZ,U);break;
    case 3: expJ(Lzmat,J,T,Hxc,Hext,lnZ,U);break;
   }
    }
  return true;
}


// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the matrix elements of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
int ic1ion_module::dspindensity_coeff1(int &tn,          // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,    // Input Lande g-factor
 /* Not Used */       Vector &ABC,    // Input vector of parameters from single ion property file
                      char *sipffilename,// Single ion properties filename
                      ComplexVector &Slm1,// Output Llm1 vector (1,49)
                      int & xyz,            // Indicating which of x,y,z direction to calculate
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ int nt=0,n,nd, Hsz = est.Rows()-1;
  switch(xyz)
  {case 1: nt=Sxmat.u1((complexdouble*)&Slm1[1],Slm1.Hi(),T,tn,delta,(complexdouble*)&est[1][0],(complexdouble*)&est[0][1],Hsz,n,nd);break;
   case 2: nt=Symat.u1((complexdouble*)&Slm1[1],Slm1.Hi(),T,tn,delta,(complexdouble*)&est[1][0],(complexdouble*)&est[0][1],Hsz,n,nd);break;
   case 3: nt=Szmat.u1((complexdouble*)&Slm1[1],Slm1.Hi(),T,tn,delta,(complexdouble*)&est[1][0],(complexdouble*)&est[0][1],Hsz,n,nd);break;
   }
   
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the matrix elements of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
// --------------------------------------------------------------------------------------------------------------- //
int ic1ion_module::dorbmomdensity_coeff1(int &tn,        // Input transition number; if tn<0, print debug info
                      double &T,          // Input temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double &g_J,    // Input Lande g-factor
 /* Not Used */       Vector &ABC,    // Input vector of parameters from single ion property file
                      char *sipffilename,// Single ion properties filename
                      ComplexVector &Llm1,// Output Llm1 vector (1,49)
                      int & xyz,          // Indicating which of x,y,z direction to calculate
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{   int nt=0,n,nd,Hsz = est.Rows()-1;
  switch(xyz)
  {case 1: nt=Lxmat.u1((complexdouble*)&Llm1[1],Llm1.Hi(),T,tn,delta,(complexdouble*)&est[1][0],(complexdouble*)&est[0][1],Hsz,n,nd);break;
   case 2: nt=Lymat.u1((complexdouble*)&Llm1[1],Llm1.Hi(),T,tn,delta,(complexdouble*)&est[1][0],(complexdouble*)&est[0][1],Hsz,n,nd);break;
   case 3: nt=Lzmat.u1((complexdouble*)&Llm1[1],Llm1.Hi(),T,tn,delta,(complexdouble*)&est[1][0],(complexdouble*)&est[0][1],Hsz,n,nd);break;
   }
   return nt;
}

// --------------------------------------------------------------------------------------------------------------- //
// returns operator matrices (n=0 Hamiltonian, n=1,...,nofcomponents: operators of moment components)
// --------------------------------------------------------------------------------------------------------------- //
bool ic1ion_module::opmat(int &n,                      // ni     which operator 0=Hamiltonian, 1,2,3=J1,J2,J3
             char *sipffilename,         // Single ion properties filename
             Vector &Hxc,                 // Hext  vector of external field [meV]
             Vector &Hext,                // Hxc   vector of exchange field [meV]
                                          // on output   
             Matrix &outmat)              // operator matrix of Hamiltonian, I1, I2, I3 depending on n
{  
    const char *sipffile = sipffilename;

   int nn = abs(n)-1;

   if(n==0)                               // return Hamiltonian
   {
      std::vector<double> gjmbH(max(6,Hxc.Hi()),0.); for(int i=1; i<=Hxc.Hi(); i++) gjmbH[i-1]=-Hxc(i);
   // --------------------------------------------------------------------
      if(fabs(Hext(1))>DBL_EPSILON) { gjmbH[3]-=MUB*Hext(1); gjmbH[0]-=GS*MUB*Hext(1); }
      if(fabs(Hext(2))>DBL_EPSILON) { gjmbH[4]-=MUB*Hext(2); gjmbH[1]-=GS*MUB*Hext(2); }
      if(fabs(Hext(3))>DBL_EPSILON) { gjmbH[5]-=MUB*Hext(3); gjmbH[2]-=GS*MUB*Hext(3); }
  
      // check dimensions of vector
      if(Hxc.Hi()>51) {
         fprintf(stderr,"Error module ic1ion: dimension of exchange field=%i > 51 - check number of columns in file mcphas.j\n",Hxc.Hi()); exit(EXIT_FAILURE); }

      // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
      #ifdef JIJCONV
      if(pars.B.norm().find("Stevens")!=std::string::npos) mfmat.jijconv.assign(pars.jijconv.begin(),pars.jijconv.end());
      #endif
      sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,gjmbH); 

       Jmat+=Hic; if(!iHic.isempty()) iJmat+=iHic; 

      if(pars.truncate_level!=1)  {       // Truncates the matrix, and packs it into real upper / imag lower triangle format
         truncate_hmltn_packed(pars, Jmat, iJmat, outmat, sipffile); return true; }
      else {
         zmat2pack(Jmat,iJmat,outmat); return true; }
   } 
   else
   {  
      if(nn>50) {
         fprintf(stderr,"Error module ic1ion: operatormatrix index=%i > 51 - check number of columns in file mcphas.j\n",n); exit(EXIT_FAILURE); }

      // Indices n 7-11 are k=2 quadrupoles; 12-18:k=3; 19-27:k=4; 28-38:k=5; 39-51:k=6
      //         nn=abs(n)-1
      int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
    //int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
    //int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
      int Hsz=getdim(pars.n,pars.l); sMat<double> zeroes; zeroes.zero(Hsz,Hsz);

      // Checks if the reduced matrix element is zero, if so, return zero without calculating matrix elements.
      double redmat = pow(-1.,(double)abs(pars.l)) * (2*pars.l+1) * threej(2*pars.l,2*k[nn],2*pars.l,0,0,0);
      if(nn>5 && fabs(redmat)<DBL_EPSILON*100)
      {
         if(pars.truncate_level!=1) { int cb = (int)(pars.truncate_level*(double)Hsz); zeroes.zero(cb,cb); }
         zmat2pack(zeroes,zeroes,outmat); return true;
      }
       
      // Calculates the operator matrices <Sx>, <Lx>, etc.
      sMat<double> Jmat=mfmat.op_generate(nn);
      
         if(pars.truncate_level!=1 || n<0) {
            if(mfmat.iflag[nn]==0) truncate_hmltn_packed(pars,Jmat,zeroes,outmat,sipffile); else truncate_hmltn_packed(pars,zeroes,Jmat,outmat,sipffile); }
         else {
            if(mfmat.iflag[nn]==0) zmat2pack(Jmat,zeroes,outmat);                  else zmat2pack(zeroes,Jmat,outmat); }
         return true;
     
      mfmat.op_free(nn);
   }

   std::cerr << "ic1ion::opmat - failed to calculate operator matrices\n"; return false;
}
