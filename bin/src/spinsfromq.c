/***********************************************************************
 * spinfromq - program to create spinsconfiguration from q vector
 * Author: Martin Rotter
 ***********************************************************************/


#include "../../version"
#include "spincf.hpp"
#include "martin.h"
#include<par.hpp>

/**********************************************************************/
// hauptprogramm
int main (int argc, char **argv)
{ 
printf("#*****************************************************\n");
printf("#* spinfromq - create spin configuration from q vector\n");
printf("#* Author: Martin Rotter %s\n",MCPHASVERSION);
printf("#*****************************************************\n");

// check command line
  if (argc < 6)
    { printf (" program spinsfromq - create spinconfiguration from q vector\n \
                use as: spinsfromq  n1 n2 n3 h k l [h2 k2 l2]\n \
		n1 n2 n3 .... periodicity of supercell\n \
		h k l ....... components of qvector\n");
      exit (1);
    }

  par inputpars("./mcphas.j");
 
  int i,j,n1,n2,n3;
  double lnz,u;

  double T;
  Vector h(1,inputpars.nofcomponents),hext(1,3);
  Vector moment(1,inputpars.nofcomponents);
  Vector qvector (1,3);
  Vector nettom(1,inputpars.nofcomponents*inputpars.nofatoms);
  Vector phi(1,inputpars.nofcomponents*inputpars.nofatoms);
  phi=0;
  Vector momentq0(1,inputpars.nofcomponents*inputpars.nofatoms);
  momentq0=0;
  n1=strtol(argv[1],NULL,10);  
  n2=strtol(argv[2],NULL,10);  
  n3=strtol(argv[3],NULL,10);
  qvector(1)=strtod(argv[4],NULL);
  qvector(2)=strtod(argv[5],NULL);
  qvector(3)=strtod(argv[6],NULL);
     
   T=1;hext=0;nettom=0;
   h=0;for(i=1;i<=inputpars.nofcomponents;++i)h(i)=0.1;
  while(Norm(nettom)<0.1){
  for(i=1;i<=inputpars.nofatoms;++i)
  {(*inputpars.jjj[i]).Icalc(moment,T,h,hext,lnz,u,(*inputpars.jjj[i]).Icalc_parameter_storage_init(h,hext,T));
   for(j=1;j<=inputpars.nofcomponents;++j){nettom(j+(i-1)*inputpars.nofcomponents)=moment(j);
                                          // set phases according to atomic position phi=2*pi*q*r
                                           phi(j+(i-1)*inputpars.nofcomponents)=qvector*(*inputpars.jjj[i]).xyz*2.0*PI;                                        
                                          }
  }
  h*=2;
  }
  printf("# q=%g %g %g:\n",qvector(1),qvector(2),qvector(3));
  // transform hkl to primitive reciprocal lattice
  qvector=qvector*inputpars.rez.Inverse();
  printf("# Miller indices of q vector with respect to primitive reciprocal lattice: (%6.4f %6.4f %6.4f)\n",qvector(1),qvector(2),qvector(3));
  
  
 
  spincf savspin (n1,n2,n3,inputpars.nofcomponents,inputpars.nofatoms);
  spincf savspin1 (n1,n2,n3,inputpars.nofcomponents,inputpars.nofatoms);
  spincf savspin2 (n1,n2,n3,inputpars.nofcomponents,inputpars.nofatoms);
  // n1 n2 n3 .. periodicity of supercell
// nettom .... saturation moment (positive)
// qvector ... wave vector in units of reciprocal lattice
// momentq0 .. ferromagnetic component (between 0 and 1)
// phi ....... phase (for each component)
  savspin.spinfromq (n1,n2,n3,qvector,nettom, momentq0, phi);
 if (argc>7)
  {double s,t;
   printf("# double q option\n");
   s=nettom(1);t=nettom(2); nettom(2)=0;nettom(1)=s;
   savspin1.spinfromq (n1,n2,n3,qvector,nettom, momentq0, phi);
   savspin1.print(stdout);

   nettom=0;nettom(2)=t;
   qvector(1)=strtod(argv[7],NULL);
   qvector(2)=strtod(argv[8],NULL);
   qvector(3)=strtod(argv[9],NULL);
   printf("# q2=%g %g %g:\n",qvector(1),qvector(2),qvector(3));
   qvector=qvector*inputpars.rez.Inverse();
   printf("# Miller indices of q2 with respect to primitive reciprocal lattice: (%6.4f %6.4f %6.4f)\n",qvector(1),qvector(2),qvector(3));
   savspin2.spinfromq (n1,n2,n3,qvector,nettom, momentq0, phi);
   savspin2.print(stdout);
   savspin=savspin1+savspin2;
   printf("# double q structure:\n");
  }
 
  inputpars.savelattice(stdout);
  inputpars.saveatoms(stdout);
  printf("0 0 1 0 0 0 0 %i %i %i\n",n1*n2*n3*inputpars.nofatoms,inputpars.nofatoms,inputpars.nofcomponents);
  savspin.print(stdout);
  return 0;
}


