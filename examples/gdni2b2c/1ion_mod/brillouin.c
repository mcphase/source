//{\footnotesize \begin{verbatim}

// module brillouin.c
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
{  fprintf(stdout,"brillouin.so: is loaded\n");}

// called just before removing from memory
void _fini(void)
{  fprintf(stdout,"brillouin.so: is removed\n");}

//routine Icalc for brillouin 
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
if(J.Hi()!=3||gjmbH.Hi()!=3||ABC.Hi()!=1)
   {fprintf(stderr,"Error loadable module brillouin.so: wrong number of dimensions - check number of columns in file mcphas.j or number of parameters in single ion property file\n");
    exit(EXIT_FAILURE);}
    
double JJ,K_BT,XJ,gmhkt,Jav,gmh,Z,X;

// program brillouin function for S=J=ABC(1)
JJ=ABC[1];
K_BT=(*T)*K_B;
gmh=Norm(gjmbH);
gmhkt=gmh/K_BT;
if(JJ*gmhkt>100||gmhkt>100){Jav=JJ;(*lnZ)=JJ*gmhkt;}
 else
 {X=exp(gmhkt);
  XJ=exp(JJ*gmhkt);

//printf("1-X=%g gmh=%g",1-X,gmh);

// calculate Brillouin function and partition sum Z
if (X<=1.000001){Z=2*JJ+1;Jav=0;}
else
{Z=(XJ*X-1/XJ)/(X-1.0);
 Jav=JJ*(XJ*X*X-1/XJ)+(JJ+1)*X*(1.0/XJ-XJ);
 Jav/=(X-1);
 Jav/=(XJ*X-1/XJ);
}
// the above formula is equivalent to the following summing routine:
//for (i=-JJ*2;i<=+0.000001;++i)
//{dd=i*gmhkt;
// if (dd<-700){expp=0;}else{expp=exp(dd);}
// Z += expp; //this is not yet Z, a factor exp(J gJ Heff/kT) is missing
//}
//Jav=0;
//for (i=-JJ*2;i<=+0.000001;++i)
//{dd=i*gmhkt;
// if (dd<-700){expp=0;}else{expp=exp(dd);}
// Jav+=(JJ+i)*expp/Z;
//}
//Z*=exp(JJ*gmhkt); //this is now the correct Z

// calculate magnetic energy U
(*lnZ)=log(Z);
}
(*U)=-gmh*Jav;



if (gmh>0)
{  J[1] = Jav*gjmbH(1)/gmh;
  J[2] = Jav*gjmbH(2)/gmh;
  J[3] = Jav*gjmbH(3)/gmh;
 }
 else
 {J=0;}
//  printf ("Ha=%g Hb=%g Hc=%g ma=%g mb=%g mc=%g \n", H[1], H[2], H[3], m[1], m[2], m[3]);
return;
}

#ifdef __MINGW32__
extern "C" __declspec(dllexport) void mcalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      ComplexMatrix & est)
#else
extern "C" void mcalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      ComplexMatrix & est)
#endif
{double lnZ,U;
 Icalc(J,T,gjmbHxc,Hext,g_J,ABC,sipffile,&lnZ,&U,est);
 double GJ=2.0;
 J*=GJ;
}
/**************************************************************************/
// for mcdisp this routine is needed
#ifdef __MINGW32
extern "C" __declspec(dllexport) int du1calc(int & tn,double & T,Vector & gjmbHxc,Vector & Hext,double * g_J,Vector & ABC, char ** sipffile,
                       ComplexVector & u1,float & delta,int & n, int & nd,ComplexMatrix & est)
#else
extern "C" int du1calc(int & tn,double & T,Vector & gjmbHxc,Vector & Hext,double * g_J,Vector & ABC, char ** sipffile,
                       ComplexVector & u1,float & delta,int & n , int & nd, ComplexMatrix & est)
#endif
{ 
  /*on input
    tn          transition-number
    ABC         A,M,Ci...saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    g_J		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splittings [meV] 
    u1(i)	transition vector elements ...
*/
static Vector J(1,3);
int pr;
Vector gjmbH(1,gjmbHxc.Hi());
gjmbH=gjmbHxc+(*g_J)*MU_B*Hext;
// clalculate thermal expectation values (needed for quasielastic scattering)
//  Icalc(J,&T,gjmbH,g_J,ABC,&lnz,&u);
  pr=1;
  if (tn<0) {pr=0;tn*=-1;}
  if (T<0){T=-T;}
  double JJ,K_BT,XJ,gmhkt,gmh,Z,R,X,sinth,hxxyy,jjkt,corr;
  complex <double> i(0,1),bx,by,bz;

// program brillouin function for S=J=ABC(1)
  JJ=ABC[1];
  K_BT=T*K_B;
  gmh=Norm(gjmbH);
  gmhkt=gmh/K_BT;
  X=exp(gmhkt);
  XJ=exp(JJ*gmhkt);
// calculate Z and R
if (X==1.0){Z=2*JJ+1;R=0;}
else
{if(X>1e50){Z=XJ;R=-2.0*JJ*XJ;}
 else
 {Z=(XJ*X-1/XJ)/(X-1.0);
  R=JJ*(1/XJ-XJ*X*X)+(JJ+1)*X*(XJ-1.0/XJ);
  R/=0.5*(X-1)*(X-1);
 }
}

// calculate coefficients bx,by,bz
 hxxyy=gjmbH(1)*gjmbH(1)+gjmbH(2)*gjmbH(2);
 if (hxxyy/gjmbH(3)/gjmbH(3)>SMALL*SMALL)
 {sinth=sqrt(hxxyy)/gmh;
  bx=-gjmbH(2)+i*gjmbH(1)*gjmbH(3)/gmh;
  bx/=2*gmh*sinth;
  by=gjmbH(1)+i*gjmbH(2)*gjmbH(3)/gmh;
  by/=2*gmh*sinth;
  }
 else
 {sinth=0;by=0.5;
  if(gjmbH(3)>0)
  {bx=0.5*i;}
  else
  {bx=-0.5*i;}
 }
  bz=-i*sinth*0.5;
// -----------------------------------------

if (tn==2) // transition to finite energy
 {delta=gmh; //set delta !!!

 if (delta>SMALL)
  {// now lets calculate mat
  u1(1)=bx*sqrt(-R/Z);
  u1(2)=by*sqrt(-R/Z);
  u1(3)=bz*sqrt(-R/Z);
  } else
  {// quasielastic scattering needs epsilon * nm / KT ....
  jjkt=0.6666667*JJ*(JJ+1)/K_BT;
  u1(1)=bx*sqrt(jjkt);
  u1(2)=by*sqrt(jjkt);
  u1(3)=bz*sqrt(jjkt);
  }
 }
 else
 { delta=-SMALL; // tn=1 ... transition within the same level
   if(X==1.0){jjkt=JJ*(2*JJ*JJ+3*JJ+1)/3/K_BT/(2*JJ+1);}
   else {if(X>1e50)
         {jjkt=-JJ*JJ*K_BT;}
         else 
         {jjkt=(1-2*JJ-2*JJ*JJ)/XJ;
         jjkt+=JJ*JJ/X/XJ;
	 jjkt+=(JJ*JJ+2*JJ+1)*X/XJ;
	 jjkt-=(JJ+1)*(JJ+1)*XJ;
	 jjkt+=(2*JJ*JJ+2*JJ-1)*XJ*X;
	 jjkt-=JJ*JJ*XJ*X*X;
	 jjkt*=X/(1-X)/(1-X);
	 jjkt/=(1/XJ-X*XJ)*K_BT;
         jjkt/=(1/XJ-X*XJ);
         corr=JJ*(XJ*X*X-1/XJ)+(JJ+1)*X*(1/XJ-XJ);
         corr/=(1-X)*(1/XJ-X*XJ);
         jjkt-=corr*corr;
         jjkt/=K_BT;
         }
        }
 // now lets calculate mat
 u1(1)=gjmbH(1)*sqrt(jjkt)/gmh;
 u1(2)=gjmbH(2)*sqrt(jjkt)/gmh;
 u1(3)=gjmbH(3)*sqrt(jjkt)/gmh;
 }
if (pr==1) printf ("delta=%4.6g meV\n",delta);

// brillouin function has 2 effective transitions
return 2;
}
// \end{verbatim}}
