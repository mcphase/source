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

//routine Icalc for kramers doublet
#ifdef __MINGW32__
extern "C" __declspec(dllexport) void Icalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      double * lnZ,double * U,ComplexMatrix & Pst)
#else
extern "C" void Icalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      double * lnZ,double * U,ComplexMatrix & Pst)
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
if(J.Hi()!=1||gjmbHxc.Hi()!=1||ABC.Hi()!=1)
   {fprintf(stderr,"Error loadable module kramer.so: wrong number of dimensions - check number of columns in file mcphas.j or number of parameters in single ion property file\n");
    exit(EXIT_FAILURE);}
Vector gjmbH(1,gjmbHxc.Hi());
gjmbH(1)=gjmbHxc(1)+(*g_J)*MU_B*Hext(1);
    
    
  double  betar, lambdap,lambdap_KBT, lambdap2, expp, expm, np, nm;
  double nennerp, nennerm, jap, jam,   Z;
  double alpha_lambdap,alphaplambdap;//,alphaxlambdap;
//  alpha = ABC[2] * gjmbH[2];
  betar = -ABC[1] * gjmbH[1];
//  betai = -ABC[3] * gjmbH[3];

  lambdap2 = betar * betar ;
  lambdap = sqrt (lambdap2);
  lambdap_KBT=lambdap/KB/(*T);
  if (lambdap_KBT>700){lambdap_KBT=700;}
  if (lambdap_KBT<-700){lambdap_KBT=-700;}
  expm = exp (lambdap_KBT);
  expp = 1/expm; //=exp (-lambdap_KBT);
  Z = expp + expm;
  (*lnZ)=log(Z);
  np = expp / Z;
  nm = expm / Z;
  (*U)=lambdap*(np-nm); // energy
//printf("T=%g expp=%g expm=%g \n",(*T),expp,expm);

//  nennerp = (alpha - lambdap) * (alpha - lambdap) + betar * betar + betai * betai;
//  nennerm = (alpha + lambdap) * (alpha + lambdap) + betar * betar + betai * betai;
 //   alphaxlambdap=0;//alpha*lambdap;
    alpha_lambdap=-lambdap;
    alphaplambdap=+lambdap;
    nennerp=  2.0*lambdap2;    
    nennerm=  2.0*lambdap2;    

  if (nennerp > SMALL)
    {
      jap = -ABC[1] * 2.0 * betar * (alpha_lambdap) / nennerp;
//      jbp = M * ((alpha_lambdap) * (alpha_lambdap) - (betar * betar + betai * betai)) / nennerp;
//    jbp = ABC[2] * (2.0 * alpha*alpha_lambdap) / nennerp;
//    jcp = -2.0 * ABC[3] * betai * (alpha_lambdap) / nennerp;
    }
  else
    {
      jap = 0;
/*    if (alpha * alpha > SMALL)
	{
	  jbp = -copysign (ABC[2], alpha);
	}
      else
	{
	  jbp = 0;
	}
      jcp = 0;*/
    }

  if (nennerm > SMALL)
    {
      jam = -ABC[1] * 2.0 * betar * (alphaplambdap) / nennerm;
//      jbm = M * ((alpha + lambdap) * (alpha + lambdap) - (betar * betar + betai * betai)) / nennerm;
     //jbm = ABC[2] * (2.0 * alpha*alphaplambdap) / nennerm;
     //jcm = -2.0 * ABC[3] * betai * (alphaplambdap) / nennerm;
    }
  else
    {
      jam = 0;
     /* if (alpha * alpha > SMALL)
	{
	  jbm = copysign (ABC[2], alpha);
	}
      else
	{
	  jbm = 0;
	}
      jcm = 0;*/
    }

  J[1] = np * jap + nm * jam;
//  J[2] = np * jbp + nm * jbm;
//  J[3] = np * jcp + nm * jcm;
// printf ("np=%g nm=%g jap=%g jbp=%g jcp=%g jam=%g jbm=%g jcm=%g \n",np,nm,jap,jbp,jcp,jam,jbm,jcm);
//  printf ("gjmbHa=%g gjmbHb=%g gjmbHc=%g Ja=%g Jb=%g Jc=%g \n", gjmbH[1], gjmbH[2], gjmbH[3], J[1], J[2], J[3]);
return;
}
//\end{verbatim}

#ifdef __MINGW32__
extern "C" __declspec(dllexport) void mcalc(Vector & m,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      ComplexMatrix & est)
#else
extern "C" void mcalc(Vector & m,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      ComplexMatrix & est)
#endif
{double lnZ,U;
 Vector J(1,1);
 Icalc(J,T,gjmbHxc,Hext,g_J,ABC,sipffile,&U,&lnZ,est);
 m=0;m(1)=J(1);

}