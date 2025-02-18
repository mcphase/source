//\begin{verbatim}
// example c file for dynamically loadable module of program
// mcphas ... it must not be c++, but pure c compiled with gcc and linked 
// with ld  !! The calculation has been compared to the internal (doublet)
// routine of mcphas
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector.h>


// this is called directly after loading it into memory from dlopen
void _init(void)
{  fprintf(stdout,"bosefermi.so: is loaded\n");}

// called just before removing from memory
void _fini(void)
{  fprintf(stdout,"bosefermi.so: is removed\n");}

//routine Icalc for bose/fermi statistics
#ifdef __MINGW32__
extern "C" __declspec(dllexport) void Icalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      double * lnZ,double * U,ComplexMatrix & est)
#else
extern "C" void Icalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      double * lnZ,double * U,ComplexMatrix & est)
#endif
{   
    /*on input
    T		temperature[K]
    gJmbH	vector of effective field [meV]
    gJ          Lande factor
    ABC         ABC(1) ... spin quantum number S=J
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/
Vector gjmbH(1,gjmbHxc.Hi());
gjmbH=gjmbHxc+(*g_J)*MU_B*Hext;
// check dimensions of vector
if(J.Hi()!=3||gjmbH.Hi()!=3||ABC.Hi()!=3)
   {fprintf(stderr,"Error loadable module bosefermi.so: wrong number of dimensions - check number of columns in file mcphas.j or number of parameters in single ion property file\n");
    exit(EXIT_FAILURE);}
    
double JJ,KBT,i,expp,dd,gmhkt,Jav,gmh,sigma,sum,sum0,sigma0,sigman,eps;
int bose,k,maxloop;
// program bose/fermi function for S=J=ABC(1)

JJ=ABC[1];

if (fabs(JJ-rint(JJ))<0.0001) 
{bose=1;sigma=1/2/(JJ+1);}       //bose case
else {bose=-1;sigma=1/2/JJ;}     //fermi case

KBT=(*T)*KB;
gmh=Norm(gjmbH);
gmhkt=gmh/KBT;

//determination of sigma=exp(mu/kT)

//calculation of sum(sigma) routine I
sum=0.0;
for (i=0;i<=2*JJ+0.00001;++i)
 {dd=i*gmhkt; // dd = Ei/kT
  if (dd>700){expp=0;}else{expp=sigma*exp(-dd);}
  sum+=expp/(1-expp*bose);
 }
sum-=1.0;
sigma0=sigma;sum0=sum;

if (bose==1){sigma=0.5;}else{sigma=10;}

eps=0.1;k=0;maxloop=1000;

//iteration loop for sigma
while(fabs(sum)>0.0000001)
{++k;if (k>maxloop)   {
 fprintf(stderr,"Error loadable module bosefermi.so: T=%g gjmbH=%g -  sigma does not converge after %i loops\n",(*T),gmh,maxloop);
 exit(EXIT_FAILURE);}

 //calculation of sum(sigma) routine II (same as I)
 sum=0.0;
 for (i=0;i<=2*JJ+0.00001;++i)
  {dd=i*gmhkt; // dd = Ei/kT
   if (dd>700){expp=0;}else{expp=sigma*exp(-dd);}
   sum+=expp/(1-expp*bose);
  }
 sum-=1.0;

 sigman=sigma0-eps*sum0*(sigma-sigma0)/(sum-sum0);
 sum0=sum;sigma0=sigma;
 sigma=sigman;
 if (fabs(sum)>fabs(sum0)&&eps>1e-7){eps=eps*0.1;}
// printf("sum=%g sigma=%g\n",sum,sigma);
}

//printf("#sigma=%g mu=%gmeV\n",sigma,KBT*log(sigma));

// now sigma is determined - determine Z, U, <J>
double Z;
Z=1.0;
Jav=0;(*U)=0;
for(i=0;i<=2*JJ+0.000001;++i)
{dd=i*gmhkt; // dd= Ei/kT
 if (dd<-700){expp=0;}else{expp=exp(-dd);}
if (bose==1){Z/=(1-bose*sigma*expp);}
else {Z*=(1-bose*sigma*expp);}

(*U)+=(i-JJ)*gmh/(1/expp/sigma-bose); 
 Jav+=(JJ-i)/(1/expp/sigma-bose);
}
(*lnZ)=log(Z);

  J[1] = Jav*gjmbH(1)/gmh;
  J[2] = Jav*gjmbH(2)/gmh;
  J[3] = Jav*gjmbH(3)/gmh;
//  printf ("Ha=%g Hb=%g Hc=%g ma=%g mb=%g mc=%g \n", H[1], H[2], H[3], m[1], m[2], m[3]);
return;
}
/**************************************************************************/
// for mcdisp this routine is needed
#ifdef __MINGW32__
extern "C" __declspec(dllexport) int dIcalc(int & tn,double & T,Vector & gjmbH,double * g_J,Vector & ABC, char ** sipffile,
                       ComplexMatrix & mat,float & delta,ComplexMatrix & est)
#else
extern "C" int dIcalc(int & tn,double & T,Vector & gjmbH,double * g_J,Vector & ABC, char ** sipffile,
                       ComplexMatrix & mat,float & delta,ComplexMatrix & est)
#endif
{ 
  /*on input
    tn          transition-number - meaningless for kramers doublet, because there is only one transition
    ABC         A,M,Ci...saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    g_J		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splitting of kramers doublet [meV]
    mat(i,j)	<-|Ji|+><+|Jj|-> tanh(delta/2kT)
*/
// NOT IMPLEMENTED
 printf("Error: external module bosefermi.c has not implemented dIcalc function");exit(EXIT_FAILURE);
 return 1; 

}

//\end{verbatim}
