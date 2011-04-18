/*********************=*************************************************
 *
 * pointc - calculates crystal field parameters by the pointcharge model
 * Reference: Martin Rotter and Ernst Bauer - Crystal field effects 
 * in Rare Earth Compounds, in print
 ***********************************************************************/


#include "jjjpar.hpp"
#include "../../version"
#include "martin.h"


/**********************************************************************/
// main program
int main (int argc, char **argv)
{// check command line


  if (argc < 3)
    { printf ("\nProgram to calculate Crystal field Parameters from Point Charges \n\n\
            Usage: pointc Ce3+ 0.2 4 1 5.3\n\n\
             ... meaning calculate Blms (Stevens Parameters) \n\
                             and   Llms (Wybourne Parameters) \n\
             for one pointcharge of +0.2|e| in distance\n\
             x=4 A y=1 A z=5.3 A from a Ce3+ ion.\n\n\
             Alternative Usage: pointc Ce3+ filename\n\
             ... meaning read several charges+coordinates from file,\n\
             file format: column 1=charge, column 2-4 = x y z coordinate.\n\
             (note,progam makenn creates useful files for this option from \n\
             the crystal structure).\n\
             results are written to stdout (including radial matrix elements\n\
             and Stevens factors)\n\n\
          Note: if an ion is not implemented, it's parameters can be \n\
                entered in a single ion property file and pointc is\n\
                started as \n\
                  pointc file.sipf 0.2 4 1 5.3\n\n\
\n\
                the single ion property file must then contain the \n\
                following information (# denotes comments):\n\
\n\
                  #the name of the ion\n\
                  IONTYPE=Ce3+\n\
                  #stevens parameters(optional,necessary for output of Blm)\n\
                  ALPHA=-0.0571429\n\
                  BETA=0.00634921\n\
                  GAMMA=0\n\
                  # the radial matrix elements RN=<r^N> \n\
                  # in units of a0^N (a0=0.5292 A)\n\
                  R2=1.309\n\
                  R4=3.964\n\
                  R6=23.31\n\
                  # radial wave function parameters, for transition metal ions the the values are tabulated in\n\
                  # Clementi & Roetti Atomic data and nuclear data tables 14 (1974) 177-478, the radial wave\n\
                  # function is expanded as R(r)=sum_p Cp r^(Np-1) . exp(-XIp r) . (2 XIp)^(Np+0.5) / sqrt(2Np!)\n\
                  # for rare earth ions see Freeman & Watson PR 127(1962)2058, Sovers J. Phys. Chem. Sol. 28(1966)1073\n\
                  # e.g. Co2+ is isoelectronic to Fe+, looking at page 422\n\
                  # of Clemente & Roetti the parameters are \n\
                  N1=3 XI1=4.95296 C1=0.36301 \n\
                  N2=3 XI2=12.2963 C2=0.02707 \n\
                  N3=3 XI3=7.03565 C3=0.14777\n\
                  N4=3 XI4=2.74850 C4=0.49771 \n\
                  N5=3 XI5=1.69027 C5=0.11388 \n\
                  # if the above parameters are given the radial wave function\n\
                  # is output to file results\\radwavfun.dat\n\
	\n");
      exit (1);
    }

FILE * table_file;
FILE * sipf_file;
char instr[MAXNOFCHARINLINE];
int n=0;
float invalues[100];invalues[0]=99;
  double q,x,y,z;
// set stevens parameters and landefactor, J and <r^l> of ion


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
                                      printf("%s",instr);
                                    }
   
//   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
        extract(instr,"IONTYPE",(*iops).iontype,(size_t)MAXNOFCHARINLINE);
        
        extract(instr,"N1",(*jjjps).Np(1));extract(instr,"XI1",(*jjjps).Xip(1));extract(instr,"C1",(*jjjps).Cp(1));
        extract(instr,"N2",(*jjjps).Np(2));extract(instr,"XI2",(*jjjps).Xip(2));extract(instr,"C2",(*jjjps).Cp(2));
        extract(instr,"N3",(*jjjps).Np(3));extract(instr,"XI3",(*jjjps).Xip(3));extract(instr,"C3",(*jjjps).Cp(3));
        extract(instr,"N4",(*jjjps).Np(4));extract(instr,"XI4",(*jjjps).Xip(4));extract(instr,"C4",(*jjjps).Cp(4));
        extract(instr,"N5",(*jjjps).Np(5));extract(instr,"XI5",(*jjjps).Xip(5));extract(instr,"C5",(*jjjps).Cp(5));
        extract(instr,"N6",(*jjjps).Np(6));extract(instr,"XI6",(*jjjps).Xip(6));extract(instr,"C6",(*jjjps).Cp(6));
        extract(instr,"N7",(*jjjps).Np(7));extract(instr,"XI7",(*jjjps).Xip(7));extract(instr,"C7",(*jjjps).Cp(7));
        extract(instr,"N8",(*jjjps).Np(8));extract(instr,"XI8",(*jjjps).Xip(8));extract(instr,"C8",(*jjjps).Cp(8));
        extract(instr,"N9",(*jjjps).Np(9));extract(instr,"XI9",(*jjjps).Xip(9));extract(instr,"C9",(*jjjps).Cp(9));

        extract(instr,"ALPHA",(*iops).alpha);
        extract(instr,"BETA",(*iops).beta);
        extract(instr,"GAMMA",(*iops).gamma);

        extract(instr,"R2",  (*jjjps).r2);
        extract(instr,"R4",  (*jjjps).r4);
        extract(instr,"R6",  (*jjjps).r6);
//        }
  }
      if((*jjjps).r2==0){(*jjjps).r2_from_radial_wavefunction();printf("#<r^2>  from radial wavefunction in units of a0^2 a0=0.5292 Angstroem\nR2=%g\n",(*jjjps).r2);}
      if((*jjjps).r4==0){(*jjjps).r4_from_radial_wavefunction();printf("#<r^4>  from radial wavefunction in units of a0^4 a0=0.5292 Angstroem\nR4=%g\n",(*jjjps).r4);}
      if((*jjjps).r6==0){(*jjjps).r6_from_radial_wavefunction();printf("#<r^6>  from radial wavefunction in units of a0^6 a0=0.5292 Angstroem\nR6=%g\n",(*jjjps).r6);}
     
      (*iops).r2=(*jjjps).r2;
      (*iops).r4=(*jjjps).r4;
      (*iops).r6=(*jjjps).r6;
  fclose(sipf_file);
 }
 else
 {printf ("#!MODULE=so1ion\n#<!--mcphase.sipf-->\n");
  printf("#*********************=*************************************************\n");
  printf("# program pointc - crystal field parameters by the pointcharge model\n");
  printf("# Martin Rotter, %s\n",MCPHASVERSION);
  printf("# Reference: Ernst Bauer and Martin Rotter - Crystal field effects \n");
  printf("#            in Rare Earth Compounds, in print\n");
  printf("#***********************************************************************\n");
  iops=new ionpars(argv[1]);  // read ion parameters from internal table  
  printf ("IONTYPE=%s\n",(*iops).iontype);
// printout the information used in pointc to output 
  printf("#J=%4g\n",(*iops).J);
  printf("#Lande Factor gJ\n GJ = %4g\n",(*iops).gJ);
  printf("#Stevens factors\nALPHA=%4g\nBETA=%4g\nGAMMA=%4g\n",(*iops).alpha,(*iops).beta,(*iops).gamma);
  printf("#Expectation values of radial wave function <r^k> in units of a0^k a0=0.5292 Angstroem\n");
  printf("R2=%4g\nR4=%4g\nR6=%4g\n\n",(*iops).r2,(*iops).r4,(*iops).r6);
  printf("# calculation of magnetisation density instead of chargedensity\n");
  printf("# switch calcmagdensity triggers: 0 ... normal mode, 1,2,3 calc <J'i>=gJ/2 (<J1,2,3 * Ji>+<Ji*J1,2,3>)\n");
  printf("#                                                    ... gives magnetisationdensity in a b c dir instead, useful only in chrgplt,charges !\n");
  printf("calcmagdensity=0\n\n");
 }

// zero parameters in case initialisation put some values to the parameters ...
(*iops).Blm=0;
(*iops).Llm=0;

if (argc<5) // read pointcharges from file
{table_file=fopen_errchk(argv[2],"r");
 while(n==0&&feof(table_file)==false)n=inputline(table_file, invalues);
  q=invalues[1];
  x=invalues[2];
  y=invalues[3];
  z=invalues[4]; 
} else 
{ n=4;
  q=strtod(argv[2],NULL);
  x=strtod(argv[3],NULL);
  y=strtod(argv[4],NULL);
  z=strtod(argv[5],NULL);
}

 // print information about pointcharges to file and calculate Blms and Llms
  printf ("\n#pointcharges charge[|e|]  x[A] y[A] z[A]\n");
while(n>0)
{

  printf ("pointcharge= %4g         %4g %4g %4g\n",q,x,y,z);
  Vector B(1,45); B=0;
  Vector gamma(1,45); gamma=0; 

 // calculate Blm's and Llm's
 double r,ct,ct2,st,st2,sfi,cfi;
 r = sqrt(x * x + y * y + z * z);
 ct = z/r;                 //z
 ct2 = ct * ct;      
 st = sqrt(x*x+y*y)/r;
 st2 = st * st;
 if((x*x+y*y)==0){sfi=0;cfi=1;}
 else
 {sfi =  y/sqrt(x*x+y*y);
 cfi =  x/sqrt(x*x+y*y);} 

  int l,m;   
 // cnst is the Zlm constants - put them into the matrix (same constants are used in jjjpar.cpp and ionpars.cpp)
 Matrix cnst(0,6,-6,6); 

 cnst(2,0) = 0.3153962;
 cnst(2,1)=  1.092548;
 cnst(2,2)=  0.5462823;
 cnst(4,0)=  0.1057871;
 cnst(4,1)=  0.6690465;
 cnst(4,2)=  0.4730943;
 cnst(4,3)=  1.77013;
 cnst(4,4)=  0.625845;
 cnst(6,0)=  0.06357014;
 cnst(6,1)=  0.582621;
 cnst(6,2)=  0.4606094;
 cnst(6,3)=  0.921205;
 cnst(6,4)=  0.5045723;
 cnst(6,5)=  2.366619;
 cnst(6,6)=  0.6831942;
 for(l=2;l<=6;l+=2){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);}

 // evaluate the Zlm in order to get gamma_lm
 gamma(1)= cnst(2, -2)  * 2 * st2 * sfi * cfi;
 gamma(2)= cnst(2, -1)  * st * sfi * ct;
 gamma(3)= cnst(2, 0)  * (3 * ct2 - 1);
 gamma(4)= cnst(2, 1)  * st * cfi * ct;
 gamma(5)= cnst(2, 2)  * st2 * (cfi * cfi - sfi * sfi); 

 gamma(13)= cnst(4, -4) * st2 * st2 * 4 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
 gamma(14)= cnst(4, -3) * ct * st * st2 * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
 gamma(15)= cnst(4, -2) * (7 * ct2 - 1) * 2 * st2 * cfi * sfi;
 gamma(16)= cnst(4, -1) * st * sfi * ct * (7 * ct2 - 3);
 gamma(17)= cnst(4, 0) * (35 * ct2 * ct2 - 30 * ct2 + 3);
 gamma(18)= cnst(4, 1)  * st * cfi * ct * (7 * ct2 - 3);
 gamma(19)= cnst(4, 2)  * (7 * ct2 - 1) * st2 * (cfi * cfi - sfi * sfi);
 gamma(20)= cnst(4, 3)  * ct * st * st2 * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
 gamma(21)= cnst(4, 4)  * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi); 

 gamma(33)= cnst(6, -6) * st2 * st2 * st2 * (6 * cfi * cfi * cfi * cfi * cfi * sfi - 20 * cfi * cfi * cfi * sfi * sfi * sfi + 6 * cfi * sfi * sfi * sfi * sfi * sfi);
 gamma(34)= cnst(6, -5) * ct * st * st2 * st2 * (5 * cfi * cfi * cfi * cfi * sfi - 10 * cfi * cfi * sfi * sfi * sfi + sfi * sfi * sfi * sfi * sfi);
 gamma(35)= cnst(6, -4) * (11 * ct2 - 1) * 4 * st2 * st2 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
 gamma(36)= cnst(6, -3) * (11 * ct * ct2 - 3 * ct) * st2 * st * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
 gamma(37)= cnst(6, -2) * 2 * st2 * sfi * cfi * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2);
 gamma(38)= cnst(6, -1) * ct * st * sfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
 gamma(39)= cnst(6, 0)  * (231 * ct2 * ct2 * ct2 - 315 * ct2 * ct2 + 105 * ct2 - 5);
 gamma(40)= cnst(6, 1)  * ct * st * cfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
 gamma(41)= cnst(6, 2)  * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2) * st2 * (cfi * cfi - sfi * sfi);
 gamma(42)= cnst(6, 3)  * (11 * ct * ct2 - 3 * ct) * st2 * st * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
 gamma(43)= cnst(6, 4)  * (11 * ct2 - 1) * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);
 gamma(44)= cnst(6, 5)  * ct * st * st2 * st2 * (cfi * cfi * cfi * cfi * cfi - 10 * cfi * cfi * cfi * sfi * sfi + 5 * cfi * sfi * sfi * sfi * sfi);
 gamma(45)= cnst(6, 6) * st2 * st2 * st2 * (cfi * cfi * cfi * cfi * cfi * cfi - 15 * cfi * cfi * cfi * cfi * sfi * sfi + 15 * cfi * cfi * sfi * sfi * sfi * sfi - sfi * sfi * sfi * sfi * sfi * sfi);


 // this is squaring of the coefficients of Zlm, a technical trick in
 // order to save a multiplication later (good for the Blm)
 for(l=2;l<=6;l+=2){for(m=-l;m<=l;++m)cnst(l,m)*=cnst(l,m);}

 //ro = a(0, 0) / sqrt(4.0 * 3.1415);

 //evaluate th Zlm in order to get Blm
 B(1)= cnst(2, -2)  * 2 * st2 * sfi * cfi;
 B(2)= cnst(2, -1)  * st * sfi * ct;
 B(3)= cnst(2, 0)  * (3 * ct2 - 1);
 B(4)= cnst(2, 1)  * st * cfi * ct;
 B(5)= cnst(2, 2)  * st2 * (cfi * cfi - sfi * sfi);

 B(13)= cnst(4, -4) * st2 * st2 * 4 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
 B(14)= cnst(4, -3) * ct * st * st2 * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
 B(15)= cnst(4, -2) * (7 * ct2 - 1) * 2 * st2 * cfi * sfi;
 B(16)= cnst(4, -1) * st * sfi * ct * (7 * ct2 - 3);
 B(17)= cnst(4, 0) * (35 * ct2 * ct2 - 30 * ct2 + 3);
 B(18)= cnst(4, 1)  * st * cfi * ct * (7 * ct2 - 3);
 B(19)= cnst(4, 2)  * (7 * ct2 - 1) * st2 * (cfi * cfi - sfi * sfi);
 B(20)= cnst(4, 3)  * ct * st * st2 * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
 B(21)= cnst(4, 4)  * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);

 B(33)= cnst(6, -6) * st2 * st2 * st2 * (6 * cfi * cfi * cfi * cfi * cfi * sfi - 20 * cfi * cfi * cfi * sfi * sfi * sfi + 6 * cfi * sfi * sfi * sfi * sfi * sfi);
 B(34)= cnst(6, -5) * ct * st * st2 * st2 * (5 * cfi * cfi * cfi * cfi * sfi - 10 * cfi * cfi * sfi * sfi * sfi + sfi * sfi * sfi * sfi * sfi);
 B(35)= cnst(6, -4) * (11 * ct2 - 1) * 4 * st2 * st2 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
 B(36)= cnst(6, -3) * (11 * ct * ct2 - 3 * ct) * st2 * st * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
 B(37)= cnst(6, -2) * 2 * st2 * sfi * cfi * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2);
 B(38)= cnst(6, -1) * ct * st * sfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
 B(39)= cnst(6, 0)  * (231 * ct2 * ct2 * ct2 - 315 * ct2 * ct2 + 105 * ct2 - 5);
 B(40)= cnst(6, 1)  * ct * st * cfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
 B(41)= cnst(6, 2)  * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2) * st2 * (cfi * cfi - sfi * sfi);
 B(42)= cnst(6, 3)  * (11 * ct * ct2 - 3 * ct) * st2 * st * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
 B(43)= cnst(6, 4)  * (11 * ct2 - 1) * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);
 B(44)= cnst(6, 5)  * ct * st * st2 * st2 * (cfi * cfi * cfi * cfi * cfi - 10 * cfi * cfi * cfi * sfi * sfi + 5 * cfi * sfi * sfi * sfi * sfi);
 B(45)= cnst(6, 6) * st2 * st2 * st2 * (cfi * cfi * cfi * cfi * cfi * cfi - 15 * cfi * cfi * cfi * cfi * sfi * sfi + 15 * cfi * cfi * sfi * sfi * sfi * sfi - sfi * sfi * sfi * sfi * sfi * sfi);


 // now calculation of the coefficients gammaLM  in cgs
 int i;
 double eps0=8.854187817e-12; //units C^2/Nm^2
 double echarge=1.60217646e-19;  // units C
 //gamma2M   
 for (i=1;i<=5;++i){B(i)*=q/r/r/r*4*PI/5; gamma(i)*=q*echarge*1e30/r/r/r/5/eps0;}
 //gamma4M
 for (i=13;i<=21;++i){B(i)*=q/r/r/r/r/r*4*PI/9; gamma(i)*=q*echarge*1e50/r/r/r/r/r/9/eps0;}
 //gamma6M
 for (i=33;i<=45;++i){B(i)*=q/r/r/r/r/r/r/r*4*PI/13; gamma(i)*=q*echarge*1e70/r/r/r/r/r/r/r/13/eps0; }

 // ... gammas are calculated in SI units [N m^(2-2L-1) /C]

 double e,a0,umr,ehv2,ehv4,ehv6;
     e = 4.80325E-10; // elementarladung
 // einheit von r in <r^n> ist Bohrradius^n = a0^n in angstroem^n
      a0 = .5292;//(Angstroem)
 //   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in mJ = 10^4
 //   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in THz =1.509166084e22
 //   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in meV =0.624146e23
 //   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in K =0.72429024e24
     umr = 6.24146E+22;
     ehv2 = a0*a0 * umr;
     ehv4 = a0*a0*a0*a0 * umr;
     ehv6 = a0*a0*a0*a0*a0*a0 * umr;

 double J2meV=1/1.60217646e-22; // 1 millielectron volt = 1.60217646 × 10-22 joules

 // now calculation of the B_LM  and L_LM in meV
 for (i=1;i<=5;++i){(*iops).Blm(i)+=-B(i)*e*e*(*iops).r2*(*iops).alpha*ehv2; // printf("B(%i)=%g sum(B)=%g\n",i,B(i),(*iops).Blm(i));
                   if(i!=3){(*iops).Llm(i)+=-echarge*(*iops).r2*a0*a0*1e-20*gamma(i)*sqrt(5.0/8.0/PI)*J2meV;}  //m<>0
                   else    {(*iops).Llm(i)+=-echarge*(*iops).r2*a0*a0*1e-20*gamma(i)*sqrt(5.0/4.0/PI)*J2meV;}  //m=0
                  }
 for (i=13;i<=21;++i){(*iops).Blm(i)+=-B(i)*e*e*(*iops).r4*(*iops).beta*ehv4; 
                   if(i!=17){(*iops).Llm(i)+=-echarge*(*iops).r4*a0*a0*a0*a0*1e-40*gamma(i)*sqrt(9.0/8.0/PI)*J2meV;}  //m<>0
                   else     {(*iops).Llm(i)+=-echarge*(*iops).r4*a0*a0*a0*a0*1e-40*gamma(i)*sqrt(9.0/4.0/PI)*J2meV;}  //m=0
                    }
 for (i=33;i<=45;++i){(*iops).Blm(i)+=-B(i)*e*e*(*iops).r6*ehv6*(*iops).gamma;
                   if(i!=39){(*iops).Llm(i)+=-echarge*(*iops).r6*a0*a0*a0*a0*a0*a0*1e-60*gamma(i)*J2meV*sqrt(13.0/8.0/PI);}  //m<>0
                   else     {(*iops).Llm(i)+=-echarge*(*iops).r6*a0*a0*a0*a0*a0*a0*1e-60*gamma(i)*J2meV*sqrt(13.0/4.0/PI);}  //m=0
                    }
 n=0;
 if (argc<5)
 { while((n==0)&(feof(table_file)==false))n=inputline(table_file, invalues);
  q=invalues[1];
  x=invalues[2];
  y=invalues[3];
  z=invalues[4]; 
 }
}

if (argc<5){fclose(table_file);}

printf ("# 1 meV = 8.066 cm-1\n# 1 cm-1 = 0.124 meV\n");
printf ("# 1 meV = 11.6 K\n# 1 K = 0.0862 meV\n");


printf("#-------------------------------------------------------\n");
printf("# Crystal Field parameters Blm in Stevens Notation (meV)\n");
printf("#-------------------------------------------------------\n");
(*iops).savBlm(stdout);
printf("#--------------------------------------------------------\n");
printf("# Crystal Field parameters Llm in Wybourne Notation (meV)\n");
printf("#--------------------------------------------------------\n");
(*iops).savLlm(stdout);
  fprintf(stderr,"#***********************************************************************\n");
  fprintf(stderr,"#                         end of program pointc\n");
  fprintf(stderr,"# Reference: Ernst Bauer and Martin Rotter - Magnetism of Complex\n");
  fprintf(stderr,"#            Metallic Alloys: Crystalline Electric Field Effects \n");
  fprintf(stderr,"#            Book Series on Complex Metallic Alloys - Vol. 2, edited\n");
  fprintf(stderr,"#            by Esther Belin-Ferré, World Scientific, 2009\n");
  fprintf(stderr,"#***********************************************************************\n");

}

/* COMMENT COMMENT COMMENT COMMENT COMMENT !!!!!!!!!!!!!!!
//REM q     ...........parameter - reduced charges [|e|]
//'   xNN,yNN,zNN ..... position of charge [A]
//REM******berechnung der b#() aus den reduced charges q..*******
//
//sv20 = 0: sv22 = 0: sv40 = 0: sv42 = 0: sv43 = 0: sv44 = 0: sv60 = 0: sv62 = 0: sv63 = 0: sv64 = 0: sv66 = 0
/
/REM berechnung der Entwicklungskoeffizienten svlm der Zlm
/REM svlmadd
/sv20 = sv20 + FNV20(z, r) / r ^ 3 * 4 * pi# / 5 * q     '= gamma20 in cgs
/sv22 = sv22 + FNV22(x, y, z) / r ^ 3 * 4 * pi# / 5 * q
/sv40 = sv40 + FNV40(z, r) / r ^ 5 * 4 * pi# / 9 * q
/sv42 = sv42 + FNV42(x, y, z) / r ^ 5 * 4 * pi# / 9 * q
/sv43 = sv43 + FNV43(x, y, z) / r ^ 5 * 4 * pi# / 9 * q
/sv44 = sv44 + FNV44(x, y, z) / r ^ 5 * 4 * pi# / 9 * q
/sv60 = sv60 + FNV60(z, r) / r ^ 7 * 4 * pi# / 13 * q
/sv62 = sv62 + FNV62(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
/sv63 = sv63 + FNV63(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
/sv64 = sv64 + FNV64(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
/sv66 = sv66 + FNV66(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
/// sub for calculation of charge density given a radiu R and polar angles teta, fi and expansion coeff. alm
/// some functions 
/DEF FNV20 (z, r) = .25 * SQR(5 / pi#) * (3 * z * z - r * r) / r / r
/  DEF FNV22 (x, y, z) = 1 / 4 * SQR(15 / pi#) * (x * x - y * y) / (x * x + y * y + z * z)
/  DEF FNV40 (z, r) = 3 / 16 * SQR(1 / pi#) * (35 * z * z * z * z - 30 * z * z * r * r + 3 * r * r * r * r) / r / r / r / r
/  DEF FNV42 (x, y, z) = 3 / 8 * SQR(5 / pi#) * (7 * z * z - (x * x + y * y + z * z)) * (x * x - y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
/  DEF FNV43 (x, y, z) = 3 / 8 * SQR(70 / pi#) * z * x * (x * x - 3 * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
/  DEF FNV44 (x, y, z) = 3 / 16 * SQR(35 / pi#) * (x * x * x * x - 6 * x * x * y * y + y * y * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
/  DEF FNV60 (z, r) = 1 / 32 * SQR(13 / pi#) * (231 * z * z * z * z * z * z - 315 * z * z * z * z * r * r + 105 * z * z * r * r * r * r - 5 * r * r * r * r * r * r) / r / r / r / r / r / r
/  DEF FNV62 (x, y, z) = 1 / 64 * SQR(2730 / pi#) * (16 * z * z * z * z - 16 * (x * x + y * y) * z * z + (x * x + y * y) * (x * x + y * y)) * (x * x - y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
/  DEF FNV63 (x, y, z) = 1 / 32 * SQR(2730 / pi#) * z * x * (11 * z * z - 3 * (x * x + y * y + z * z)) * (x * x - 3 * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
/  DEF FNV64 (x, y, z) = 21 / 32 * SQR(13 / 7 / pi#) * (11 * z * z - x * x - y * y - z * z) * (x * x * x * x - 6 * x * x * y * y + y * y * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
/  DEF FNV66 (x, y, z) = 231 / 64 * SQR(26 / 231 / pi#) * (x * x * x * x * x * x - 15 * x * x * x * x * y * y + 15 * x * x * y * y * y * y - y * y * y * y * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
/  DEF FNV66S (x, y, z) = 231 / 32 * SQR(26 / 231 / pi#) * y * x * (3 * x * x - y * y) * (x * x - 3 * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)


/REM umrechnung in Koeffizienten von Olm
/     REM <r^n>-werte in bohrradius^n

/IF site% = 1 THEN r2 = 1.2: r4 = 3.455: r6 = 21.226: REM ce3+
/PRINT #1, "<r^2>="; r2; " a0^2  <r^4>="; r4; " a0^4  <r^6>="; r6; " a0^6    a0=0.5292 Angstroem}"

/     REM stevensfaktoren
/IF site% = 1 THEN alpha = -2 / (5 * 7): beta = 2 / 9 / 5 / 7: gamma = 0: REM ce3+


/  b#(2, 0) = -e ^ 2 * r2 * alpha * sv20 * SQR(5 / 16 / pi#) * ehv2
/  b#(2, 2) = -e ^ 2 * r2 * alpha * sv22 * SQR(15 / 16 / pi#) * ehv2

/  b#(4, 0) = -e ^ 2 * r4 * beta * sv40 * 3 / 16 * SQR(1 / pi#) * ehv4
/  b#(4, 2) = -e ^ 2 * r4 * beta * sv42 * 3 / 8 * SQR(5 / pi#) * ehv4
/  b#(4, 3) = -e ^ 2 * r4 * beta * sv43 * 3 / 8 * SQR(70 / pi#) * ehv4
/  b#(4, 4) = -e ^ 2 * r4 * beta * sv44 * 3 / 16 * SQR(35 / pi#) * ehv4

/  b#(6, 0) = -e ^ 2 * r6 * gamma * sv60 / 32 * SQR(13 / pi#) * ehv6
/  b#(6, 2) = -e ^ 2 * r6 * gamma * sv62 / 64 * SQR(2730 / pi#) * ehv6
/  b#(6, 3) = -e ^ 2 * r6 * gamma * sv63 / 32 * SQR(2730 / pi#) * ehv6
/  b#(6, 4) = -e ^ 2 * r6 * gamma * sv64 * 21 / 32 * SQR(13 / 7 / pi#) * ehv6
/  b#(6, 6) = -e ^ 2 * r6 * gamma * sv66 * 231 / 64 * SQR(26 / 231 / pi#) * ehv6
*/
