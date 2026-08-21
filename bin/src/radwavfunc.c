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
    { printf ("\nProgram to calculate radial wave function\n\n\
            Usage: radwavfunc *.sipf\n\n\
            The radial wave functrion is calculated from the coefficients given in the single\n\
            ion property file *.sipf, output is stored in radwavfunc.out\n\
	\n");
      exit (1);
    }

FILE * sipf_file;
char instr[MAXNOFCHARINLINE];
//float invalues[100];invalues[0]=99;
// set stevens parameters and landefactor, J and <r^l> of ion
parser ob; int verbose=0;

// read create class object ionpars from iontype - sets J, gJ, Stevens factors from the
// routine getpar in cfieldrout.c, thus takes the single ion parameters from
// the same source as the cfield program ...
 ionpars * iops;
 jjjpar * jjjps;
 char *token;
 if((sipf_file=fopen(argv[1],"r"))) //read ion parameters from file
 { iops=new ionpars(2);
   jjjps=new jjjpar(1,1,1);
   while(feof(sipf_file)==false)
  {if(fgets(instr, MAXNOFCHARINLINE, sipf_file)){// strip /r (dos line feed) from line if necessary
                                      while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
 //                                      printf("%s",instr);
                                    }
   
//   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
        parseline(instr,ob,verbose);
        extract(instr,"IONTYPE",(*iops).iontype,(size_t)MAXNOFCHARINLINE,1);
        
        extract(instr,"N1",(*jjjps).Np(1),ob);extract(instr,"XI1",(*jjjps).Xip(1),ob);extract(instr,"C1",(*jjjps).Cp(1),ob);
        extract(instr,"N2",(*jjjps).Np(2),ob);extract(instr,"XI2",(*jjjps).Xip(2),ob);extract(instr,"C2",(*jjjps).Cp(2),ob);
        extract(instr,"N3",(*jjjps).Np(3),ob);extract(instr,"XI3",(*jjjps).Xip(3),ob);extract(instr,"C3",(*jjjps).Cp(3),ob);
        extract(instr,"N4",(*jjjps).Np(4),ob);extract(instr,"XI4",(*jjjps).Xip(4),ob);extract(instr,"C4",(*jjjps).Cp(4),ob);
        extract(instr,"N5",(*jjjps).Np(5),ob);extract(instr,"XI5",(*jjjps).Xip(5),ob);extract(instr,"C5",(*jjjps).Cp(5),ob);
        extract(instr,"N6",(*jjjps).Np(6),ob);extract(instr,"XI6",(*jjjps).Xip(6),ob);extract(instr,"C6",(*jjjps).Cp(6),ob);
        extract(instr,"N7",(*jjjps).Np(7),ob);extract(instr,"XI7",(*jjjps).Xip(7),ob);extract(instr,"C7",(*jjjps).Cp(7),ob);
        extract(instr,"N8",(*jjjps).Np(8),ob);extract(instr,"XI8",(*jjjps).Xip(8),ob);extract(instr,"C8",(*jjjps).Cp(8),ob);
        extract(instr,"N9",(*jjjps).Np(9),ob);extract(instr,"XI9",(*jjjps).Xip(9),ob);extract(instr,"C9",(*jjjps).Cp(9),ob);

        extract(instr,"ALPHA",(*iops).alpha,ob);
        extract(instr,"BETA",(*iops).beta,ob);
        extract(instr,"GAMMA",(*iops).gamma,ob);

//        extract(instr,"R2",  (*jjjps).r2,ob);
//        extract(instr,"R4",  (*jjjps).r4,ob);
//        extract(instr,"R6",  (*jjjps).r6,ob);
//        }
  }
      if((*jjjps).r2==0){(*jjjps).r2_from_radial_wavefunction();printf("#<r^2> from radial wavefunction  in units of a0^2 a0=0.5292 Angstroem\nR2=%g\n",(*jjjps).r2);}
      if((*jjjps).r4==0){(*jjjps).r4_from_radial_wavefunction();printf("#<r^4> from radial wavefunction  in units of a0^4 a0=0.5292 Angstroem\nR4=%g\n",(*jjjps).r4);}
      if((*jjjps).r6==0){(*jjjps).r6_from_radial_wavefunction();printf("#<r^6> from radial wavefunction in units of a0^6 a0=0.5292 Angstroem\nR6=%g\n",(*jjjps).r6);}
      if((*jjjps).Np(1)!=0)
      {// save radial wavefunction
      (*jjjps).save_radial_wavefunction("results/radwavfunc.out");
      }
      (*iops).r2=(*jjjps).r2;
      (*iops).r4=(*jjjps).r4;
      (*iops).r6=(*jjjps).r6;
  fclose(sipf_file);
 }
 else
 {fprintf (stderr,"Error - cannot read %s\n",argv[1]);exit(1); }

// zero parameters in case initialisation put some values to the parameters ...
(*iops).Blm=0;
(*iops).Llm=0;

  fprintf(stderr,"#***********************************************************************\n");
  fprintf(stderr,"#          radial wave function written to results/radwavfunc.out\n");
  fprintf(stderr,"#                         end of program radwavfunc\n");
  fprintf(stderr,"#***********************************************************************\n");

}

