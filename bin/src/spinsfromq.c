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
    { printf ("# program spinsfromq - create spinconfiguration from q vector\n \
#                use as: spinsfromq [-m f1 f2 f3 ...] n1 n2 n3 h k l [fi][h2 k2 l2 [fi2] [h3 k3 l3 [fi3]]]\n \
#		n1 n2 n3 .... periodicity of supercell\n \
#		h k l [fi]....... components of qvector and optional phase fi (in degree)\n \
#                                 (if fi is given it has to be given for all q vectors)\n \
#                options:\n \
#                -m f1 f2 f3 ...fnofcomponents: multiply spin components by f1 f2 f3 ... fnofcomponents\n \
#      formulas: M(r)=(f1*M1,f2*M2,f3*M3,..., fnofcomponents * Mnofcomponents)*cos(Q.r+fi)    ...for single q\n \
#                M(r)=f1*M1*(1 0 0)*cos(Q1.r+fi)+f2*M2*(0 1 0)*cos(Q2.r+fi2) ...for double q\n \
#                M(r)=f1*M1*(1 0 0)*cos(Q1.r+fi)+f2*M2* (0 1 0)*cos(Q2.r+fi2)+f3*M3*(0 0 1)*cos(Q3.r+fi3) ...for triple q\n \
#      M1,M2,M3 are determined by applying an molecular field with equal components and increasing this \n \
#      field until the total length of the moment exceeds a treshold value. \n \
");
      exit (1);
    }

  par inputpars("./mcphas.j");
 
  int i,j,n1,n2,n3;
  double lnz,u,fi1=0,fi2=0,fi3=0;
  int a;
  double T;
  Vector h(1,inputpars.cs.nofcomponents),hext(1,3);
  Vector moment(1,inputpars.cs.nofcomponents);
  Vector factors(1,inputpars.cs.nofcomponents);
  Vector qvector (1,3);
  Vector nettom(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
  Vector nettom1(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
  Vector nettom2(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
  Vector nettom3(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
  nettom1=0;nettom2=0;nettom3=0;
  Vector phi(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
  phi=0;
  Vector momentq0(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
  momentq0=0;a=0;for(i=1;i<=inputpars.cs.nofcomponents;++i)factors(i)=1.0;
  if(argv[1][0]=='-'){// treat options
                   if(argv[1][1]=='m'){for(i=1;i<=inputpars.cs.nofcomponents;++i)factors(i)=strtod(argv[1+i],NULL);
                                   }
                   else {fprintf(stderr,"Error spinsfromq: option -%c not known\n",argv[1][1]);exit(1);}
                   a=1+inputpars.cs.nofcomponents;
                  }
  n1=strtol(argv[a+1],NULL,10);  
  n2=strtol(argv[a+2],NULL,10);  
  n3=strtol(argv[a+3],NULL,10);
  qvector(1)=strtod(argv[a+4],NULL);
  qvector(2)=strtod(argv[a+5],NULL);
  qvector(3)=strtod(argv[a+6],NULL);
     
   T=1;hext=0;nettom=0;
   h=0;for(i=1;i<=inputpars.cs.nofcomponents;++i)h(i)=0.1;
  while(Norm(nettom)<0.1&&Norm(h)<1e10){
  for(i=1;i<=inputpars.cs.nofatoms;++i)
  {(*inputpars.jjj[i]).Icalc(moment,T,h,hext,lnz,u,(*inputpars.jjj[i]).Icalc_parameter_storage_init(h,hext,T));
   for(j=1;j<=inputpars.cs.nofcomponents;++j){nettom(j+(i-1)*inputpars.cs.nofcomponents)=moment(j);
                                          }
    nettom1(1+(i-1)*inputpars.cs.nofcomponents)=moment(1);
    nettom2(2+(i-1)*inputpars.cs.nofcomponents)=moment(2);
    nettom3(3+(i-1)*inputpars.cs.nofcomponents)=moment(3);                                          
  }
  h*=2; // multiply field h by two until moment netoom gets large enough to have some sizable numbers
  }
  // treat case of zero moment loop ended at high h ...
  if(Norm(h)>=1e10)for(i=1;i<=inputpars.cs.nofcomponents;++i)
           for(j=1;j<=inputpars.cs.nofcomponents;++j)nettom(j+(i-1)*inputpars.cs.nofcomponents)=1;
                                         
    

  printf("#! Moment components: ");
  for(i=1;i<=inputpars.cs.nofatoms;++i)
  {for(j=1;j<=inputpars.cs.nofcomponents;++j){if(i==1)printf("M%i=%6.4f ",j,nettom(j+(i-1)*inputpars.cs.nofcomponents));
                                           nettom(j+(i-1)*inputpars.cs.nofcomponents)*=factors(j);
                                          }
    nettom1(1+(i-1)*inputpars.cs.nofcomponents)*=factors(1);
    nettom2(2+(i-1)*inputpars.cs.nofcomponents)*=factors(2);
    nettom3(3+(i-1)*inputpars.cs.nofcomponents)*=factors(3);                                          
  }

  printf("\n#! n1=%i n2=%i n3=%i nofspins=%i\n",n1,n2,n3,n1*n2*n3*inputpars.cs.nofatoms);
    for(i=1;i<=inputpars.cs.nofatoms;++i)for(j=1;j<=inputpars.cs.nofcomponents;++j){// set phases according to atomic position phi=2*pi*q*r
                                           phi(j+(i-1)*inputpars.cs.nofcomponents)=qvector*(*inputpars.jjj[i]).xyz*2.0*PI;                                        
                                          }
  
  // transform hkl to primitive reciprocal lattice
  
 
  spincf savspin (n1,n2,n3,inputpars.cs.nofcomponents,inputpars.cs.nofatoms);
  spincf savspin1 (n1,n2,n3,inputpars.cs.nofcomponents,inputpars.cs.nofatoms);
  spincf savspin2 (n1,n2,n3,inputpars.cs.nofcomponents,inputpars.cs.nofatoms);
  // n1 n2 n3 .. periodicity of supercell
// nettom .... saturation moment (positive)
// qvector ... wave vector in units of reciprocal lattice
// momentq0 .. ferromagnetic component (between 0 and 1)
// phi ....... phase (for each component)
int fis_given=0;
if ((argc-a-4)%4==0){//fi  is given in command line
                     fi1=strtod(argv[a+7],NULL);
                     phi+=PI*fi1/180;++a;fis_given=1;
                    }
   printf("#! q: hkl=%g %g %g  phase fi=%6.4f deg:\n",qvector(1),qvector(2),qvector(3),fi1);
  qvector=qvector*inputpars.rez.Inverse();
  printf("# Miller indices of q vector with respect to primitive reciprocal lattice: (%6.4f %6.4f %6.4f)\n",qvector(1),qvector(2),qvector(3));

  savspin.spinfromq (n1,n2,n3,qvector,nettom, momentq0, phi);
 if (argc>7+a)
  {savspin1.spinfromq (n1,n2,n3,qvector,nettom1, momentq0, phi);
   //savspin1.print(stdout);

   qvector(1)=strtod(argv[a+7],NULL);
   qvector(2)=strtod(argv[a+8],NULL);
   qvector(3)=strtod(argv[a+9],NULL);
   for(i=1;i<=inputpars.cs.nofatoms;++i)for(j=1;j<=inputpars.cs.nofcomponents;++j){// set phases according to atomic position phi=2*pi*q*r
                                           phi(j+(i-1)*inputpars.cs.nofcomponents)=qvector*(*inputpars.jjj[i]).xyz*2.0*PI;                                        
                                          }
   if(fis_given>0){fi2=strtod(argv[a+10],NULL);
                     phi+=PI*fi2/180;++a;}
   printf("#! q2: hkl2=%g %g %g phase fi2=%6.4f deg:\n",qvector(1),qvector(2),qvector(3),fi2);
   qvector=qvector*inputpars.rez.Inverse();
   printf("# Miller indices of q2 with respect to primitive reciprocal lattice: (%6.4f %6.4f %6.4f)\n",qvector(1),qvector(2),qvector(3));
   savspin2.spinfromq (n1,n2,n3,qvector,nettom2, momentq0, phi);
   //savspin2.print(stdout);
   savspin=savspin1+savspin2;
   if (argc>10+a)
   {savspin1=savspin;
    qvector(1)=strtod(argv[a+10],NULL);
    qvector(2)=strtod(argv[a+11],NULL);
    qvector(3)=strtod(argv[a+12],NULL);
   for(i=1;i<=inputpars.cs.nofatoms;++i)for(j=1;j<=inputpars.cs.nofcomponents;++j){// set phases according to atomic position phi=2*pi*q*r
                                           phi(j+(i-1)*inputpars.cs.nofcomponents)=qvector*(*inputpars.jjj[i]).xyz*2.0*PI;                                        
                                          }
   if(fis_given>0){fi3=strtod(argv[a+13],NULL);
                     phi+=PI*fi3/180;++a;}
 printf("#! q3: hkl3=%g %g %g phase fi3=%6.4f deg:\n",qvector(1),qvector(2),qvector(3),fi3);
   qvector=qvector*inputpars.rez.Inverse();
   printf("# Miller indices of q3 with respect to primitive reciprocal lattice: (%6.4f %6.4f %6.4f)\n",qvector(1),qvector(2),qvector(3));
   savspin2.spinfromq (n1,n2,n3,qvector,nettom3, momentq0, phi);
   //savspin2.print(stdout);
   savspin=savspin1+savspin2;  
    printf("# triple q structure:\n");}
   else
   {printf("# double q structure:\n");}

  }
 
  inputpars.savelattice(stdout);
  inputpars.saveatoms(stdout);
  fprintf(stdout,"#! 0 0   0    0 0 0  %i %i  %i #created by ",savspin.n(),savspin.nofatoms,savspin.nofcomponents);
  for(i=0;i<argc;++i){fprintf(stdout,"%s ",argv[i]);}fprintf(stdout,"\n");
  savspin.print(stdout);
  printf("\n");


  return 0;
}


