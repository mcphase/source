//{\footnotesize \begin{verbatim}

// module phonon.c
// example c file for dynamically loadable module of program
// mcphas ... this must not contain c++ code, but pure c code 
// which is being compiled with gcc and linked 
// with ld  !! 

#include <cstdio>
#include <cmath>
#include <complex>
#include <vector.h>

#define MU_B 0.05788
#define K_B  0.0862
#define SMALL 1e-10


// this is called directly after loading it into memory from dlopen
void _init(void)
{  fprintf(stdout,"phonon.so: is loaded\n");}

// called just before removing from memory
void _fini(void)
{  fprintf(stdout,"phonon.so: is removed\n");}

//routine Icalc for phonon
#ifdef __MINGW32__
extern "C" __declspec(dllexport) void Icalc(Vector & u0,double * T, Vector &Fxc, Vector & Hext,double * g_J, Vector & MODPAR,char ** sipffilename,
                      double * lnZ,double * U,ComplexMatrix & Icalc_parstorage)
#else
extern "C" void Icalc(Vector & u0,double * T,Vector &Fxc, Vector & Hext,double * g_J, Vector & MODPAR,char ** sipffile,
                      double * lnZ,double * U,ComplexMatrix & Icalc_parstorage)
#endif
{   
    /*on input
    T		temperature[K]
    Hext	vector of external magnetic field [T] [not used]
    Fxc         exchange Force on the oscillator in [meV]
    gJ          Lande factor [not used]
    MODPAR      Einstein oscillator parameters
                MODPAR[1]   mass (m0)  m0=atomic mass unit=1.660539e-27 kg
                MODPAR[2]   Kxx
                MODPAR[3]   Kyy
                MODPAR[4]   Kzz            K is minus second derivative 
                MODPAR[5]   Kxy            of the potential energy with respect
                MODPAR[6]   Kxz             to nuclear displacements [meV]
                MODPAR[7]   Kyz
                MODPAR[8]   umax    maximum displacement  [a0=0.5219 A]
                MODPAR[9]   0       umax restriction in all directions
                            1,2,3   umax restriction in x y z direction only
                            4       umax restriction in x and y direction
                            5       umax restriction in x and z direction
                            6       umax restriction in y and z direction

  on output    
    u0		equilibrium position vector [a0] , a0=0.5219 A
    Z		single ion partition function
    U		single ion energy (meV)
*/


//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
// 			   int sort, int maxiter)

//
//  Driver routine to compute  the eigenvalues  and normalized
//  eigenvectors of  the real symmetric matrix z, given by the
//  lower triangle in z[lo..hi,lo..hi]. The eigenvalues are re-
//  turned in d[lo..hi] in ascending sequence if sort = True,
//  otherwise not ordered for  sort = False. The associated
//  eigenvectors overwrite the given matrix z. The storage re-
//  quirement is n*n + 2*n double.
//  The vector d must already be allocated by the user.
//
//  Reference:
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
//  void EigenSystemSymmetric (Matrix& z, Vector& d, int sort, int maxiter)

Vector F(1,3);
F(1)=Fxc(1);F(2)=Fxc(2);F(3)=Fxc(3);
 if(MODPAR.Hi()<9)
   {fprintf(stderr,"Error loadable module phonon.so:  MODPAR1 .. MODPAR9 all have to be set in sipf file\n");
    exit(EXIT_FAILURE);}

 double m=MODPAR[1]*1.660539e-27; // mass of einstein oscillator in kg
 double a0=0.5219e-10; // Bohr radius in meter
 Matrix K(1,3,1,3);K=0;
 //Matrix Sr(1,3,1,3),Si(1,3,1,3);
 Vector Omega(1,3);
 K(1,1)=MODPAR[2];
 K(2,2)=MODPAR[3];
 K(3,3)=MODPAR[4];
 K(2,1)=MODPAR[5];
 K(3,1)=MODPAR[6];
 K(3,2)=MODPAR[7];
// printf("Icalc phonon Kij= %g %g %g %g %g %g %g %g %g\n",K(1,1),K(2,2),K(3,3),K(2,1),K(3,1),K(3,2),K(1,2),K(1,3),K(2,3));
if(Norm(K)==0){(*U)=0; // last term  to correct energy
(*lnZ)=0; // last term  to correct energy
  u0=0;return;}

 int sort=1,maxiter=1000000;double factor=1e-4;
K*=factor;
// EigenSystemHermitean (K,Omega,Sr,Si,sort,maxiter); // K is destroyed by this
 EigenSystemSymmetric (K,Omega,sort,maxiter); // K is destroyed by this and will contain eigenvectors
Omega/=factor; 
// hbar=1.054572e-34 Js=6582e-16meVs
// 1meV=1.6022e-22 J
// Omegai=m a0^2 (Deltai/hbar)^2
// 
double Delta1,Delta2,Delta3,K_BT,X,Y,Z;

Delta1=sqrt(-Omega(1)*1.6022e-22/m/a0/a0)*6582e-16; // phonon einstein frequencies (meV) 
Delta2=sqrt(-Omega(2)*1.6022e-22/m/a0/a0)*6582e-16;
Delta3=sqrt(-Omega(3)*1.6022e-22/m/a0/a0)*6582e-16;
//printf("Icalc phonon Om= %g %g %g\n",Omega[1],Omega[2],Omega[3]);
//printf("Icalc phonon Delta= %g %g %g\n",Delta1,Delta2,Delta3);
if(isnan((Delta1))){fprintf (stderr, "phonon: Delta1=nan Omega(1)=%g m=%g\n",Omega(1),m);exit (EXIT_FAILURE);}
if(isnan((Delta2))){fprintf (stderr, "phonon: Delta2=nan Omega(2)=%g m=%g\n",Omega(2),m);exit (EXIT_FAILURE);}
if(isnan((Delta3))){fprintf (stderr, "phonon: Delta3=nan Omega(3)=%g m=%g\n",Omega(3),m);exit (EXIT_FAILURE);}


K_BT=(*T)*K_B;
X=exp(-Delta1/K_BT);
Y=exp(-Delta2/K_BT);
Z=exp(-Delta3/K_BT);
// calculate phonon function and partition sum Z
(*lnZ)=-Delta1/2/K_BT-Delta2/2/K_BT-Delta3/2/K_BT-log(1-X)-log(1-Y)-log(1-Z);

// the energy U is sum_i Ei exp(-Ei/(kT))/Z
// with Ei=w0(0.5+i)-1/2 FT.u0 (the last term is added to (*U) at the end of Icalc
(*U)=Delta1*(0.5+X/(1-X))+Delta2*(0.5+Y/(1-Y))+Delta3*(0.5+Z/(1-Z));
//printf("Icalc phonon initial u=%g lnz=%g \n"(*U),(*lnZ));

// u0=-ST Omega^-1 S Fxc= (SrT-iSiT)*Om*(Sr+iSi)*F
// =(SrT-iSiT) (Om Sr F + i Om Si F)= 
// SrT Om Sr F + SiT Om Si F + i (must be zero)
// with Om=-Omega^-1
//Matrix Om(1,3,1,3);Om=0;Om(1,1)=-1/Omega(1);Om(2,2)=-1/Omega(2);Om(3,3)=-1/Omega(3);
Vector uu(1,3);
//printf("Icalc phonon Om= %g %g %g\n",Omega[1],Omega[2],Omega[3]);
 K(1,1)=MODPAR[2];
 K(2,2)=MODPAR[3];
 K(3,3)=MODPAR[4];
 K(2,1)=MODPAR[5];
 K(3,1)=MODPAR[6];
 K(3,2)=MODPAR[7];
 K(1,2)=K(2,1);K(1,3)=K(3,1);K(2,3)=K(3,2);
//printf("Icalc phonon Kij= %g %g %g %g %g %g %g %g %g\n",K(1,1),K(2,2),K(3,3),K(2,1),K(3,1),K(3,2),K(1,2),K(1,3),K(2,3));

//uu=Sr.Transpose()*Om*Sr*F+Si.Transpose()*Om*Si*F;
Matrix Ki(1,3,1,3);Ki=K.Inverse();
//printf("Icalc phonon Kijinverse= %g %g %g %g %g %g %g %g %g\n",Ki(1,1),Ki(2,2),Ki(3,3),Ki(2,1),Ki(3,1),Ki(3,2),Ki(1,2),Ki(1,3),Ki(2,3));

uu=-Ki*F;
//printf("Icalc phonon uu= %g %g %g F= %g %g %g\n",uu(1),uu(2),uu(3),F(1),F(2),F(3));


int i=(int)MODPAR[9];

/*
switch(i)
{case 0: if(uu(1)>MODPAR[8])uu(1)=MODPAR[8];
         if(uu(2)>MODPAR[8])uu(2)=MODPAR[8];
         if(uu(3)>MODPAR[8])uu(3)=MODPAR[8];
         if(uu(1)<-MODPAR[8])uu(1)=-MODPAR[8];
         if(uu(2)<-MODPAR[8])uu(2)=-MODPAR[8];
         if(uu(3)<-MODPAR[8])uu(3)=-MODPAR[8];
         break;
case 1 : if(uu(1)>MODPAR[8])uu(1)=MODPAR[8];if(uu(1)<-MODPAR[8])uu(1)=-MODPAR[8];break;
case 2 : if(uu(2)>MODPAR[8])uu(2)=MODPAR[8];if(uu(2)<-MODPAR[8])uu(2)=-MODPAR[8];break;
case 3 : if(uu(3)>MODPAR[8])uu(3)=MODPAR[8];if(uu(3)<-MODPAR[8])uu(3)=-MODPAR[8];break;
case 4 : if(uu(1)>MODPAR[8])uu(1)=MODPAR[8];if(uu(1)<-MODPAR[8])uu(1)=-MODPAR[8];
         if(uu(2)>MODPAR[8])uu(2)=MODPAR[8];if(uu(2)<-MODPAR[8])uu(2)=-MODPAR[8];break;
case 5 : if(uu(1)>MODPAR[8])uu(1)=MODPAR[8];if(uu(1)<-MODPAR[8])uu(1)=-MODPAR[8];
         if(uu(3)>MODPAR[8])uu(3)=MODPAR[8];if(uu(3)<-MODPAR[8])uu(3)=-MODPAR[8];break;
case 6 : if(uu(3)>MODPAR[8])uu(3)=MODPAR[8];if(uu(3)<-MODPAR[8])uu(3)=-MODPAR[8];
         if(uu(2)>MODPAR[8])uu(2)=MODPAR[8];if(uu(2)<-MODPAR[8])uu(2)=-MODPAR[8];break;
default: break;
}
*/

// try tanh restriction
double sat;sat=MODPAR[8];
switch(i)
{case 0: uu(1)=sat*tanh(uu(1)/sat);
         uu(2)=sat*tanh(uu(2)/sat);
         uu(3)=sat*tanh(uu(3)/sat);
         break;
case 1 : uu(1)=sat*tanh(uu(1)/sat);break;
case 2 : uu(2)=sat*tanh(uu(2)/sat);break;
case 3 : uu(3)=sat*tanh(uu(3)/sat);break;
case 4 : uu(1)=sat*tanh(uu(1)/sat);
         uu(2)=sat*tanh(uu(2)/sat);break;
case 5 : uu(1)=sat*tanh(uu(1)/sat);
         uu(3)=sat*tanh(uu(3)/sat);break;
case 6 : uu(3)=sat*tanh(uu(3)/sat);
         uu(2)=sat*tanh(uu(2)/sat);break;
default: break;
}

(*U)-=0.5*uu*F; // last term  to correct energy
(*lnZ)+=0.5*uu*F/K_BT; // last term  to correct energy
 // to easy convergence of mcphasit the linearity of the einstein Oscillator is damped 
//if(isnan((*lnZ))){fprintf (stderr, "lnzi=nan %g %g %g %g %g %g  %g %g\n",uu(1),uu(2),uu(3),F(1),F(2),F(3),uu*F,K_BT);exit (EXIT_FAILURE);}

  u0=0;
  u0[1] = uu(1);
  u0[2] = uu(2); // should in principle be F/m w0^2, but we set it zero to keep atoms in equilibrium position
  u0[3] = uu(3);
//printf("Icalc phonon u= %g %g %g u=%g lnz=%g \n",u0[1],u0[2],u0[3],(*U),(*lnZ));
return;
}
/**************************************************************************/
// for spins display positions and mcdiff this routine is needed
#ifdef __MINGW32__
extern "C" __declspec(dllexport) void pcalc(Vector & u0,double * T, Vector &Fxc, Vector & Hext,double * g_J, Vector & MODPAR,char ** sipffilename,
                     ComplexMatrix & Icalc_parstorage)
#else
extern "C" void pcalc(Vector & u0,double * T, Vector &Fxc, Vector & Hext,double * g_J, Vector & MODPAR,char ** sipffilename,
                     ComplexMatrix & Icalc_parstorage)
#endif
{double lnZ,U;
 Icalc(u0,T,Fxc,Hext,g_J,MODPAR,sipffilename,&lnZ,&U,Icalc_parstorage);
 double a0=0.5219; // Bohr radius in A
 u0*=a0;
 return;
}

/**************************************************************************/
// for mcdisp this routine is needed
#ifdef __MINGW32__
extern "C" __declspec(dllexport) int du1calc(int & tn,double & T,Vector & Fxc, Vector & Hext,
                       double * g_J,Vector & MODPAR, char ** sipffilename,
                       ComplexVector & u1,float & delta,int & n, int & nd, ComplexMatrix & est)
#else
extern "C" int du1calc(int & tn,double & T,Vector & Fxc, Vector & Hext,
                       double * g_J,Vector & MODPAR, char ** sipffilename,
                       ComplexVector & u1,float & delta,int & n , int & nd, ComplexMatrix & est)
#endif
{ 
  /*on input
     T		temperature[K]
    Hext	vector of external magnetic field [T] [not used]
    Fxc         exchange Force in [meV]
    gJ          Lande factor [not used]
    MODPAR      Einstein oscillator parameters
                MODPAR[1]   mass (m0)  m0=atomic mass unit=1.660539e-27 kg
                MODPAR[2]   Kxx
                MODPAR[3]   Kyy
                MODPAR[4]   Kzz            K is the force matrix [meV]
                MODPAR[5]   Kxy
                MODPAR[6]   Kxz
                MODPAR[7]   Kyz
                MODPAR[8]   umax    maximum displacement  [a0=0.5219 A]
                MODPAR[9]   0       umax restriction in all directions
                            1,2,3   umax restriction in x y z direction only
                            4       umax restriction in x and y direction
                            5       umax restriction in x and z direction
                            6       umax restriction in y and z direction
  on output    
    delta	splittings [meV] 
    u1(i)	transition vector elements ...
*/
int pr;

  pr=0;if (tn<0) {pr=1;tn*=-1;}
  if (T<0){T=-T;} 
 if(MODPAR.Hi()<9)
   {fprintf(stderr,"Error loadable module phonon.so:  MODPAR1 .. MODPAR9 all have to be set in sipf file\n");
    exit(EXIT_FAILURE);}

 double m=MODPAR[1]*1.660539e-27; // mass of einstein oscillator in kg
 double a0=0.5219e-10; // Bohr radius in meter
 Matrix K(1,3,1,3);K=0;
 //Matric Sr(1,3,1,3),Si(1,3,1,3);
 Vector Omega(1,3);
 K(1,1)=MODPAR[2];
 K(2,2)=MODPAR[3];
 K(3,3)=MODPAR[4];
 K(2,1)=MODPAR[5];
 K(3,1)=MODPAR[6];
 K(3,2)=MODPAR[7];
 int sort=1,maxiter=1000000;
 //EigenSystemHermitean (K,Omega,Sr,Si,sort,maxiter);
EigenSystemSymmetric (K,Omega,sort,maxiter);
//printf("du1calc phonon Om= %g %g %g\n",Omega[1],Omega[2],Omega[3]);

 
// hbar=1.054572e-34 Js=6582e-16meVs
double hbar=6582e-16;
// 1meV=1.6022e-22 J
double mev2J=1.6022e-22;
// Omegai=m a0^2 (Deltai/hbar)^2
// 
delta=sqrt(-Omega(tn)*1.6022e-22/m/a0/a0)*6582e-16; // phonon einstein frequencies (meV) 

u1=0;
for(int i=1;i<=3;++i)
{//u1(i)=complex <double> (Sr(i,tn),Si(i,tn));
u1(i)=complex <double> (K(i,tn),0.0);
}
u1*=sqrt(mev2J*hbar*hbar/2/m/a0/a0/delta); // multiply by factor and convert into 1/meV units
if (pr==1) printf ("delta=%4.6g meV\n",delta);

// phonon function has 3 effective transitions
return 3;
}

/**************************************************************************/
// for mcdisp nuclear intensities this routine is needed
#ifdef __MINGW32__
extern "C" __declspec(dllexport) int dP1(int & tn,double & T,Vector & Fxc, Vector & Hext,
                       double * g_J,Vector & MODPAR, char ** sipffilename,
                       ComplexVector & P1,float & maxE,ComplexMatrix & est)
#else
extern "C" int dP1(int & tn,double & T,Vector & Fxc, Vector & Hext,
                       double * g_J,Vector & MODPAR, char ** sipffilename,
                       ComplexVector & P1,float & maxE,ComplexMatrix & est)
#endif
{int n,nd;
 int noft=du1calc(tn,T,Fxc,Hext,g_J,MODPAR,sipffilename,P1,maxE,n,nd,est);
 double a0=0.5219; // Bohr radius in A
 P1*=a0;
 return noft;
}

extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void mcalc(Vector &mom,        // Output magnetic moment (mub)
                      double *T,          // Input scalar temperature
                      Vector &Hxc,        // Input vector of exchange fields (meV) 
                      Vector &Hext,       // Input vector of external field (T) 
 /* Not Used */       double * /*g_J*/,   // Input Lande g-factor
 /* Not Used */       Vector & /*ABC*/,   // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   
   mom(1)=0;
   mom(2)=0;
   mom(3)=0;
}