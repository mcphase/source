/**************************************************************
 * singleion.c - display single ion momentum  at given htpoint
 * Author: Martin Rotter
 **************************************************************/


#include "../../version"
#include "martin.h"
#include<par.hpp>

/**********************************************************************/
void helpexit()
{ printf (" program single ion  - display single ion expectations values <Ia> <Ib> ... and transition energies at given T and H\n \
                use as: singleion [option] T[K] Hexta[T] Hextb[T] Hextc[T] Hxc1 Hxc2 Hxc3 ... Hxcnofcomponents [meV] \n\n \
                      Hext ..... external field in Tesla \n \
                      Hxc... exchange (molecular) field in meV   \n \
      by default only 5 transition energies are output, start with option -nt n to output n transition energies\n\n");
      exit (1);
}

// hauptprogramm
int main (int argc, char **argv)
{ int i,j,i0,nt;
  double lnz,u;
  float d=1e10;
  double T;

printf("#**************************************************************\n");
printf("# * singleion.c - display single ion momentum  at given htpoint\n");
printf("# * Author: Martin Rotter %s\n",MCPHASVERSION);
printf("# **************************************************************\n");
// check command line
  if (argc < 3)helpexit();

  par inputpars("./mcphas.j");
// check command line again
  if (argc-4>inputpars.nofcomponents+3)helpexit();

  Vector h(1,inputpars.nofcomponents),hext(1,3);
  Vector m(1,inputpars.nofcomponents);
  // transition matrix Mij
  ComplexVector u1(1,inputpars.nofcomponents);
  complex<double> np(1000000,0.01);u1(1)=np; // set (ninit,pinit)=(max number of initial states, minimal occupancy) to be
                        // included into list of transitions
int nmax=5;// max number of transtions to  be output
                             i=1;if(!strcmp(argv[i],"-nt")){nmax=(int)strtod(argv[i+1],NULL);i+=2;}
                             T=strtod(argv[i],NULL);

hext=0;for(j=1;j<=3;++j){++i;if(!strcmp(argv[i],"-nt")){nmax=(int)strtod(argv[i+1],NULL);i+=2;i0+=2;}
                             if(i<argc){hext(j)=strtod(argv[i],NULL);}
                            }i0=i;
h=0;for(j=1;j<=argc-i0-1;++j){++i;if(!strcmp(argv[i],"-nt")){nmax=(int)strtod(argv[i+1],NULL);i+=2;i0+=2;}
                             if(i<argc){h(j)=strtod(argv[i],NULL);}
                            }


printf("#\n#atom-number T[K] ");
for(j=1;j<=3;++j)printf("Hext%c(T) ",'a'-1+j);
for(j=1;j<=inputpars.nofcomponents;++j)printf("gjmbHxc%i(meV) ",j);
for(j=1;j<=inputpars.nofcomponents;++j)printf(" <I%c> ",'a'-1+j);
  printf("transition-energies(meV)...\n");

  for(i=1;i<=inputpars.nofatoms;++i)
{
   (*inputpars.jjj[i]).Icalc(m,T,h,hext,lnz,u,(*inputpars.jjj[i]).Icalc_parameter_storage_init(h,hext,T));
                                                                   // here calculate magnetic moment
  printf("%3i %13g ",i,T); // printout ion number and temperature
  for(j=1;j<=3;++j)printf(" %10g ",hext(j)); // printout external field as requested
  for(j=1;j<=inputpars.nofcomponents;++j)printf(" %10g ",h(j)); // printoutexchangefield as requested
  for(j=1;j<=inputpars.nofcomponents;++j)printf(" %4g ",m(j));  // printout corresponding moments 
 double TT=fabs(T);
 if(nmax>0)
  {(*inputpars.jjj[i]).transitionnumber=-1;
   nt=(*inputpars.jjj[i]).du1calc(TT,h,hext,u1,d,(*inputpars.jjj[i]).eigenstates(h,hext,T));

          // get nt = number of transitions
   for(j=1;j<=nt&&j<=nmax;++j)
   {(*inputpars.jjj[i]).transitionnumber=-j;d=1e10;
   (*inputpars.jjj[i]).du1calc(TT,h,hext,u1,d,(*inputpars.jjj[i]).est);
   printf(" %4g ",d);
   }
   if(nmax<nt){printf("...");}
  }
 printf("\n");
}




}
