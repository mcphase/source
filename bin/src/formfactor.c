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
jjjps=new jjjpar(0,0,0,argv[1]);
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
    fprintf(fout,"#        -Rare Earth Form Factor F(Q) =<j0(Q)>+<j2(Q)>*(2/gJ-1)\n\n");
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
   {fprintf(fout,"# Formfactor calculated from Form Factor coefficients - thanks to J Brown\n");
    fprintf(fout,"#   d = 2*pi/Q      \n");
    fprintf(fout,"#   s = 1/2/d = Q/4/pi   \n");
    fprintf(fout,"#   sin(theta) = lambda * s\n");
    fprintf(fout,"#    s2= s*s = Q*Q/16/pi/pi\n");
    fprintf(fout,"#\n");
    fprintf(fout,"#   <j0(Q)>=   FFj0A*EXP(-FFj0a*s2) + FFj0B*EXP(-FFj0b*s2) + FFj0C*EXP(-FFj0c*s2) + FFj0D\n");
    fprintf(fout,"#   <j2(Q)>=s2*(FFj2A*EXP(-FFj2a*s2) + FFj2B*EXP(-FFj2b*s2) + FFj2C*EXP(-FFj2c*s2) + FFj2D\n");
    fprintf(fout,"#   <j4(Q)>=s2*(FFj4A*EXP(-FFj4a*s2) + FFj4B*EXP(-FFj4b*s2) + FFj4C*EXP(-FFj4c*s2) + FFj4D\n");
    fprintf(fout,"#   <j6(Q)>=s2*(FFj6A*EXP(-FFj6a*s2) + FFj6B*EXP(-FFj6b*s2) + FFj6C*EXP(-FFj6c*s2) + FFj6D\n");
    fprintf(fout,"#\n");
    fprintf(fout,"#FFj0A=%+7.4f FFj0a=%+7.4f FFj0B=%+7.4f FFj0b=%+7.4f FFj0C=%+7.4f FFj0c=%+7.4f FFj0D=%+7.4f\n",(*jjjps).magFFj0[1],(*jjjps).magFFj0[2],(*jjjps).magFFj0[3],(*jjjps).magFFj0[4],(*jjjps).magFFj0[5],(*jjjps).magFFj0[6],(*jjjps).magFFj0[7]);
    fprintf(fout,"#FFj2A=%+7.4f FFj2a=%+7.4f FFj2B=%+7.4f FFj2b=%+7.4f FFj2C=%+7.4f FFj2c=%+7.4f FFj2D=%+7.4f\n",(*jjjps).magFFj2[1],(*jjjps).magFFj2[2],(*jjjps).magFFj2[3],(*jjjps).magFFj2[4],(*jjjps).magFFj2[5],(*jjjps).magFFj2[6],(*jjjps).magFFj2[7]);
    fprintf(fout,"#FFj4A=%+7.4f FFj4a=%+7.4f FFj4B=%+7.4f FFj4b=%+7.4f FFj4C=%+7.4f FFj4c=%+7.4f FFj4D=%+7.4f\n",(*jjjps).magFFj4[1],(*jjjps).magFFj4[2],(*jjjps).magFFj4[3],(*jjjps).magFFj4[4],(*jjjps).magFFj4[5],(*jjjps).magFFj4[6],(*jjjps).magFFj4[7]);
    fprintf(fout,"#FFj6A=%+7.4f FFj6a=%+7.4f FFj6B=%+7.4f FFj6b=%+7.4f FFj6C=%+7.4f FFj6c=%+7.4f FFj6D=%+7.4f\n",(*jjjps).magFFj6[1],(*jjjps).magFFj6[2],(*jjjps).magFFj6[3],(*jjjps).magFFj6[4],(*jjjps).magFFj6[5],(*jjjps).magFFj6[6],(*jjjps).magFFj6[7]);
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

