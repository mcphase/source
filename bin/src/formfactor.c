/*********************=*************************************************
 *
 * program formfactor used to calculate the neutron formfactor for an ion
 * in dipole approximation - parameters taken from single ion paramater file *.*
 ***********************************************************************/


#include "jjjpar.hpp"
#include "../../version"
#include "martin.h"


/**********************************************************************/
// main program
int main (int argc, char **argv)
{// check command line


  if (argc < 2)
    { printf ("\nProgram to calculate magnetic Formfactor for Neutron Scattering\n\n\
            Usage: formfactor *.sipf\n\n\
            The formfactor is calculated from the coefficients given in the single\n\
            ion property file *.sipf, also the expectation values <jl(Q)> are output\n\
            Note: if a radial wave function is parametrized in the sipf file, then\n\
            this will be used to calculate the formfactor. If not, then the formfactor\n\
            will be calculated from coefficients FF*, which have to be given in the sipf-file.\n\
	\n");
      exit (1);
    }


FILE * fout;

// read create class object ionpars from iontype - sets J, gJ, Stevens factors from the
// routine getpar in cfieldrout.c, thus takes the single ion parameters from
// the same source as the cfield program ...
  jjjpar * jjjps;
jjjps=new jjjpar(0,0,0,argv[1],1);
char s[2];

if ((*jjjps).gJ!=0){strcpy(s,"\0");}else{strcpy(s,"S");}
if((fout=fopen("results/formfactor.out","w"))) //read ion parameters from file
 {fprintf(fout,"# Formfactor for %s  ",argv[1]);
  time_t curtime;
  struct tm * loctime;
  curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
    fprintf(fout,"#   Dipole Approximation for Neutron Magnetic Formfactor:\n");
    fprintf(fout,"#        -Spin Form Factor       FS(Q)=<j0(Q)>\n");
    fprintf(fout,"#        -Angular Form Factor    FL(Q)=<j0(Q)>+<j2(Q)>\n");
    fprintf(fout,"#        -Rare Earth Form Factor F(Q) =<j0(Q)>+<j2(Q)>*(2/gJ-1)\n#\n");
    fprintf(fout,"#--------------------------------------------------------------------------------------\n");

   if((*jjjps).Np(1)!=0)
   {fprintf(fout,"# Formfactor calculated from Radial Wavefunction\n");
     fprintf(fout,"#---------------------------------------------------------------------------------------------------\n");
                    fprintf(fout,"# radial wave function parameters, for transition metal ions the the values are tabulated in\n");
                    fprintf(fout,"# Clementi & Roetti Atomic data and nuclear data tables 14 (1974) 177-478, the radial wave\n");
                    fprintf(fout,"# function is expanded as R(r)=sum_p Cp r^(Np-1) . exp(-XIp r) . (2 XIp)^(Np+0.5) / sqrt((2Np)!)\n");
                    fprintf(fout,"# for rare earth ions see Freeman & Watson PR 127(1962)2058, Sovers J. Phys. Chem. Sol. 28(1966)1073\n");
                    fprintf(fout,"#---------------------------------------------------------------------------------------------------\n");
                    for(int i=(*jjjps).Np.Lo();i<=(*jjjps).Np.Hi();++i){if((*jjjps).Np(i)!=0){fprintf(fout,"#N%i=%i XI%i=%g C%i=%g\n",i,(int)(*jjjps).Np(i),i,(*jjjps).Xip(i),i,(*jjjps).Cp(i));}}                   
    fprintf(fout,"#|Q|(1/A)   F%s(Q)     |F%s(Q)|^2     <j0(Q)>     <j2(Q)>     <j4(Q)>     <j6(Q)>      <j1(Q)>     <j3(Q)>     <j5(Q)> \n",s,s);
   }
   else
   {(*jjjps).magFFout("#",fout);
    fprintf(fout,"#\n#\n");
    fprintf(fout,"#|Q|(1/A)  F%s(Q)    |F%s(Q)|^2     <j0(Q)>     <j2(Q)>     <j4(Q)>     <j6(Q)> \n",s,s);
   }
   for(double Q=0.1;Q<20;Q+=0.1)
   {fprintf(fout,"%10.7f %10.7f %10.7f  %10.7f %10.7f %10.7f %10.7f",Q,(*jjjps).F(Q), (*jjjps).F(Q)*(*jjjps).F(Q), (*jjjps).j0(Q),(*jjjps).j2(Q),(*jjjps).j4(Q),(*jjjps).j6(Q));
        if((*jjjps).Np(1)!=0)
        {fprintf(fout,"%10.7f %10.7f %10.7f", (*jjjps).j1(Q),(*jjjps).j3(Q),(*jjjps).j5(Q));
        }
         fprintf(fout,"\n");
   }
  fclose(fout);
}
else
 {fprintf(stderr,"Error - cannot open results/formfactor.out, maybe directory results does not exist\n");exit(1);}

  fprintf(stderr,"#***********************************************************************\n");
  fprintf(stderr,"#               formfactor written to results/formfactor.out\n");
  fprintf(stderr,"#                         end of program formfactor\n");
  fprintf(stderr,"#***********************************************************************\n");

}

