/**************************************************************
 * singleion.c - display single ion momentum  at given htpoint
 * Author: Martin Rotter
 **************************************************************/


#include "../../version"
#include "martin.h"
#include<par.hpp>

/**********************************************************************/
void helpexit()
{ printf (" program single ion  - display single ion expectations values <Ja> <Jb> ... and transition energies at given T and H\n \
                use as: singleion T[K] gjmbHa gmbHb gjmbHc[meV] ...\n\n \
             by default only 5 transition energies are output, start with option -nt n to output n transition energies\n\n");
      exit (1);
}

// hauptprogramm
int main (int argc, char **argv)
{ int i,j,i0,nt;
  double lnz,u;
  float d;
  double T;

printf("#**************************************************************\n");
printf("# * singleion.c - display single ion momentum  at given htpoint\n");
printf("# * Author: Martin Rotter %s\n",MCPHASVERSION);
printf("# **************************************************************\n");
// check command line
  if (argc < 3)helpexit();

  par inputpars("./mcphas.j");
// check command line again
  if (argc-4>inputpars.nofcomponents)helpexit();

  Vector h(1,inputpars.nofcomponents);
  Vector m(1,inputpars.nofcomponents);
  // transition matrix Mij
  ComplexMatrix Mijkl(1,inputpars.nofcomponents,1,inputpars.nofcomponents);

int nmax=5;// max number of transtions to  be output
                             i=1;if(!strcmp(argv[i],"-nt")){nmax=(int)strtod(argv[i+1],NULL);i+=2;}
                             T=strtod(argv[i],NULL);i0=i;

h=0;for(j=1;j<=argc-i0-1;++j){++i;if(!strcmp(argv[i],"-nt")){nmax=(int)strtod(argv[i+1],NULL);i+=2;i0+=2;}
                             if(i<argc){h(j)=strtod(argv[i],NULL);}
                            }


printf("#\n#atom-number T[K] ");
for(j=1;j<=inputpars.nofcomponents;++j)printf("gjmbH%c(meV) ",'a'-1+j);
for(j=1;j<=inputpars.nofcomponents;++j)printf(" <J%c> ",'a'-1+j);
  printf("transition-energies(meV)...\n");

  for(i=1;i<=inputpars.nofatoms;++i)
{
   (*inputpars.jjj[i]).mcalc(m,T,h,lnz,u,(*inputpars.jjj[i]).mcalc_parameter_storage_init(h,T));
                                                                   // here calculate magnetic moment
  printf("%3i %13g ",i,T); // printout ion number and temperature
  for(j=1;j<=inputpars.nofcomponents;++j)printf(" %10g ",h(j)); // printout meanfield as requested
  for(j=1;j<=inputpars.nofcomponents;++j)printf(" %4g ",m(j));  // printout corresponding moments 
 double TT=fabs(T);
 if(nt>0)
  {(*inputpars.jjj[i]).transitionnumber=-1;
   nt=(*inputpars.jjj[i]).dmcalc(TT,h,Mijkl,d,(*inputpars.jjj[i]).eigenstates(h,T));
          // get nt = number of transitions
   for(j=1;j<=nt&&j<=nmax;++j)
   {(*inputpars.jjj[i]).transitionnumber=-j;
   (*inputpars.jjj[i]).dmcalc(TT,h,Mijkl,d,(*inputpars.jjj[i]).est);
   printf(" %4g ",d);
   }
   if(nmax<nt){printf("...");}
  }
 printf("\n");
}




}
