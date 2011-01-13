//{\footnotesize \begin{verbatim}
// example c file for dynamically loadable module of program
// mcphas ... it must not be c++, but pure c compiled with gcc and linked 
// with ld  !! The calculation has been compared to the internal (doublet)
// routine of mcphas
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector.h>

#define MU_B 0.05788
#define K_B  0.0862
#define SMALL 1e-10


// this is called directly after loading it into memory from dlopen
#ifdef __linux__
//extern "C" void _init(void)
#else
//extern "C" __declspec(dllexport)void _init(void)
#endif
//{  fprintf(stdout,"kramer.so: is loaded\n");}

// called just before removing from memory
#ifdef __linux__
//extern "C" void _fini(void)
#else
//extern "C" __declspec(dllexport) void _fini(void)
#endif
//{  fprintf(stdout,"kramer.so: is removed\n");}

//routine mcalc for kramers doublet
#ifdef __MINGW32__
extern "C" __declspec(dllexport) void mcalc(Vector & J,double * T, Vector & gjmbH,double * g_J, Vector & ABC,char ** sipffile,
                      double * lnZ,double * U,ComplexMatrix & est)
#else
extern "C" void mcalc(Vector & J,double * T, Vector & gjmbH,double * g_J, Vector & ABC,char ** sipffile,
                      double * lnZ,double * U,ComplexMatrix & est)
#endif
{   
    /*on input
    T		temperature[K]
    gJmbH	vector of effective field [meV]
    gJ          Lande factor
    ABC         single ion parameter values (A, B, C corresponding to <+|Ja|->,<-|Jb|->,<+|Jc|->/i
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/
// check dimensions of vector
if(J.Hi()!=3||gjmbH.Hi()!=3||ABC.Hi()!=3)
   {fprintf(stderr,"Error loadable module kramer.so: wrong number of dimensions - check number of columns in file mcphas.j or number of parameters in single ion property file\n");
    exit(EXIT_FAILURE);}
    
    
  double alpha, betar, betai, lambdap,lambdap_K_BT, lambdap2, expp, expm, np, nm;
  double nennerp, nennerm, jap, jam, jbp, jbm, jcp, jcm,Z;
  double alpha_lambdap,alphaplambdap,alphaxlambdap;
  alpha = ABC[2] * gjmbH[2];
  betar = -ABC[1] * gjmbH[1];
  betai = -ABC[3] * gjmbH[3];

  lambdap2 = alpha * alpha + betar * betar + betai * betai;
  lambdap = sqrt (lambdap2);
  lambdap_K_BT=lambdap/K_B/(*T);
  if (lambdap_K_BT>700){lambdap_K_BT=700;}
  if (lambdap_K_BT<-700){lambdap_K_BT=-700;}
  expm = exp (lambdap_K_BT);
  expp = 1/expm; //=exp (-lambdap_K_BT);
  Z = expp + expm;
  (*lnZ)=log(Z);
  np = expp / Z;
  nm = expm / Z;
  (*U)=lambdap*(np-nm); // energy
//printf("T=%g expp=%g expm=%g \n",(*T),expp,expm);

//  nennerp = (alpha - lambdap) * (alpha - lambdap) + betar * betar + betai * betai;
//  nennerm = (alpha + lambdap) * (alpha + lambdap) + betar * betar + betai * betai;
    alphaxlambdap=alpha*lambdap;
    alpha_lambdap=alpha-lambdap;
    alphaplambdap=alpha+lambdap;
    nennerp=  2.0*(-alphaxlambdap+lambdap2);    
    nennerm=  2.0*(alphaxlambdap+lambdap2);    

  if (nennerp > SMALL)
    {
      jap = -ABC[1] * 2.0 * betar * (alpha_lambdap) / nennerp;
//      jbp = M * ((alpha_lambdap) * (alpha_lambdap) - (betar * betar + betai * betai)) / nennerp;
      jbp = ABC[2] * (2.0 * alpha*alpha_lambdap) / nennerp;
      jcp = -2.0 * ABC[3] * betai * (alpha_lambdap) / nennerp;
    }
  else
    {
      jap = 0;
      if (alpha * alpha > SMALL)
	{
	  jbp = -copysign (ABC[2], alpha);
	}
      else
	{
	  jbp = 0;
	}
      jcp = 0;
    }

  if (nennerm > SMALL)
    {
      jam = -ABC[1] * 2.0 * betar * (alphaplambdap) / nennerm;
//      jbm = M * ((alpha + lambdap) * (alpha + lambdap) - (betar * betar + betai * betai)) / nennerm;
      jbm = ABC[2] * (2.0 * alpha*alphaplambdap) / nennerm;
      jcm = -2.0 * ABC[3] * betai * (alphaplambdap) / nennerm;
    }
  else
    {
      jam = 0;
      if (alpha * alpha > SMALL)
	{
	  jbm = copysign (ABC[2], alpha);
	}
      else
	{
	  jbm = 0;
	}
      jcm = 0;
    }

  J[1] = np * jap + nm * jam;
  J[2] = np * jbp + nm * jbm;
  J[3] = np * jcp + nm * jcm;
// printf ("np=%g nm=%g jap=%g jbp=%g jcp=%g jam=%g jbm=%g jcm=%g \n",np,nm,jap,jbp,jcp,jam,jbm,jcm);
//  printf ("gjmbHa=%g gjmbHb=%g gjmbHc=%g Ja=%g Jb=%g Jc=%g \n", gjmbH[1], gjmbH[2], gjmbH[3], J[1], J[2], J[3]);

return;
}
/**************************************************************************/
// for mcdisp this routine is needed
#ifdef __MINGW32__
extern "C" __declspec(dllexport) int dmcalc(int & tn,double & T,Vector & gjmbH,double * g_J,Vector & ABC, char ** sipffile,
                       ComplexMatrix & mat,float & delta,ComplexMatrix & est)
#else
extern "C" int dmcalc(int & tn,double & T,Vector & gjmbH,double * g_J,Vector & ABC, char ** sipffile,
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
  double alpha, betar, betai, lambdap,lambdap_K_BT, lambdap2, expp, expm, np, nm;
  double nennerp, nennerm, nenner;
  complex<double> ja,jb,jc,i(0,1),jap,jbp,jcp,jam,jbm,jcm;
  double alpha_lambdap,alphaplambdap,alphaxlambdap;
  double Z,lnz,u;
  static int pr;

  static Vector J(1,3);
  // clalculate thermal expectation values (needed for quasielastic scattering)
  mcalc(J,&T,gjmbH,g_J,ABC,sipffile,&lnz,&u,est);
  pr=1;
  if (tn<0) {pr=0;tn*=-1;}

  alpha = ABC[2] * gjmbH[2];
  betar = -ABC[1] * gjmbH[1];
  betai = -ABC[3] * gjmbH[3];
  lambdap2 = alpha * alpha + betar * betar + betai * betai;
  lambdap = sqrt (lambdap2);
  
  lambdap_K_BT=lambdap/K_B/T;
  if (lambdap_K_BT>700){lambdap_K_BT=700;}
  if (lambdap_K_BT<-700){lambdap_K_BT=-700;}
  expm = exp (lambdap_K_BT);
  expp = 1/expm; //=exp (-lambdap_K_BT);
  Z = expp + expm;
  np = expp / Z;
  nm = expm / Z;

//  nennerp = (alpha - lambdap) * (alpha - lambdap) + betar * betar + betai * betai;
//  nennerm = (alpha + lambdap) * (alpha + lambdap) + betar * betar + betai * betai;
    alphaxlambdap=alpha*lambdap;
    alpha_lambdap=alpha-lambdap;
    alphaplambdap=alpha+lambdap;
    nennerp=  2.0*(-alphaxlambdap+lambdap2);    
    nennerm=  2.0*(alphaxlambdap+lambdap2);    

if (tn==2)
 {delta=2*lambdap; //set delta !!!
    nenner=sqrt(nennerp*nennerm);

  if (nenner > SMALL)
    {
      ja = -ABC[1] * 2.0*(alpha * betar+i * betai * lambdap) / nenner;
      jb = -ABC[2] * 2.0 * (betar*betar+betai*betai) / nenner;
      jc = -ABC[3] * 2.0*(alpha*betai -i *betar*lambdap) / nenner;
    }
  else
    {
      if (alpha > SMALL)
	{ja = ABC[1];  // <-| is the ground state
  	 jb = 0;
         jc = -i*ABC[3];
	}
      else
	{ja = ABC[1];  // <+| is the ground state
  	 jb = 0;
         jc = i*ABC[3]; 	
	}
    }
 if (delta>SMALL)
  {// now lets calculate mat
  mat(1,1)=ja*conj(ja)*(nm-np);
  mat(1,2)=ja*conj(jb)*(nm-np);
  mat(1,3)=ja*conj(jc)*(nm-np);
  mat(2,1)=jb*conj(ja)*(nm-np);
  mat(2,2)=jb*conj(jb)*(nm-np);
  mat(2,3)=jb*conj(jc)*(nm-np);
  mat(3,1)=jc*conj(ja)*(nm-np);
  mat(3,2)=jc*conj(jb)*(nm-np);
  mat(3,3)=jc*conj(jc)*(nm-np);
  } else
  {// quasielastic scattering needs epsilon * nm / KT ....
  mat(1,1)=ja*conj(ja)*nm/K_B/T;
  mat(1,2)=ja*conj(jb)*nm/K_B/T;
  mat(1,3)=ja*conj(jc)*nm/K_B/T;
  mat(2,1)=jb*conj(ja)*nm/K_B/T;
  mat(2,2)=jb*conj(jb)*nm/K_B/T;
  mat(2,3)=jb*conj(jc)*nm/K_B/T;
  mat(3,1)=jc*conj(ja)*nm/K_B/T;
  mat(3,2)=jc*conj(jb)*nm/K_B/T;
  mat(3,3)=jc*conj(jc)*nm/K_B/T;  
  }
 }
 else
 { delta=-SMALL; // transition within the same level
  if (nennerp > SMALL)
    {
      jap = -ABC[1] * 2.0 * betar * (alpha_lambdap) / nennerp;
//      jbp = M * ((alpha_lambdap) * (alpha_lambdap) - (betar * betar + betai * betai)) / nennerp;
      jbp = ABC[2] * (2.0 * alpha*alpha_lambdap) / nennerp;
      jcp = -2.0 * ABC[3] * betai * (alpha_lambdap) / nennerp;
    }
  else
    {
      jap = 0;
      if (alpha * alpha > SMALL)
	{
	  jbp = -copysign (ABC[2], alpha);
	}
      else
	{
	  jbp = -ABC[2];
	}
      jcp = 0;
    }

  if (nennerm > SMALL)
    {
      jam = -ABC[1] * 2.0 * betar * (alphaplambdap) / nennerm;
//      jbm = M * ((alpha + lambdap) * (alpha + lambdap) - (betar * betar + betai * betai)) / nennerm;
      jbm = ABC[2] * (2.0 * alpha*alphaplambdap) / nennerm;
      jcm = -2.0 * ABC[3] * betai * (alphaplambdap) / nennerm;
    }
  else
    {
      jam = 0;
      if (alpha * alpha > SMALL)
	{
	  jbm = copysign (ABC[2], alpha);
	}
      else
	{
	  jbm = ABC[2];
	}
      jcm = 0;
    }
 if (tn==1)
 {// now lets calculate mat
 mat(1,1)=(jam-J(1))*(jam-J(1))*nm/K_B/T;
 mat(1,2)=(jam-J(1))*(jbm-J(2))*nm/K_B/T;
 mat(1,3)=(jam-J(1))*(jcm-J(3))*nm/K_B/T;
 mat(2,1)=(jbm-J(2))*(jam-J(1))*nm/K_B/T;
 mat(2,2)=(jbm-J(2))*(jbm-J(2))*nm/K_B/T;
 mat(2,3)=(jbm-J(2))*(jcm-J(3))*nm/K_B/T;
 mat(3,1)=(jcm-J(3))*(jam-J(1))*nm/K_B/T;
 mat(3,2)=(jcm-J(3))*(jbm-J(2))*nm/K_B/T;
 mat(3,3)=(jcm-J(3))*(jcm-J(3))*nm/K_B/T;
 }else{ // tn = 3 in this case
 // now lets calculate mat
 mat(1,1)=(jap-J(1))*(jap-J(1))*np/K_B/T;
 mat(1,2)=(jap-J(1))*(jbp-J(2))*np/K_B/T;
 mat(1,3)=(jap-J(1))*(jcp-J(3))*np/K_B/T;
 mat(2,1)=(jbp-J(2))*(jap-J(1))*np/K_B/T;
 mat(2,2)=(jbp-J(2))*(jbp-J(2))*np/K_B/T;
 mat(2,3)=(jbp-J(2))*(jcp-J(3))*np/K_B/T;
 mat(3,1)=(jcp-J(3))*(jap-J(1))*np/K_B/T;
 mat(3,2)=(jcp-J(3))*(jbp-J(2))*np/K_B/T;
 mat(3,3)=(jcp-J(3))*(jcp-J(3))*np/K_B/T;
 }
}
if (pr==1) printf ("delta=%4.6g meV\n",delta);
 

return 3;// kramers doublet has always exactly one transition + 2 levels (quasielastic scattering)!
}

//\end{verbatim}}
