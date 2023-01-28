/*********************=*************************************************
 *
 * pointc - calculates crystal field parameters by the pointcharge model
 * Reference: Martin Rotter and Ernst Bauer - Crystal field effects 
 * in Rare Earth Compounds, in print
 ***********************************************************************/


#include "jjjpar.hpp"
#include "../../version"
#include "martin.h"
#define SMALL_DISPLACEMENT  0.001

void calcCEFpar(double & q0,double & q2,double & q4,double & q6,double & x ,double & y, double & z, double & r,Vector & Blm, Vector & Llm,ionpars * iops)
{
 Vector B(0,45); B=0;
 Vector gamma(0,45); gamma=0; 
 double ct,ct2,st,st2,sfi,cfi;
 ct = z/r;                 //z
 ct2 = ct * ct;      
 st = sqrt(x*x+y*y)/r;
 st2 = st * st;
 if((x*x+y*y)==0){sfi=0;cfi=1;}
 else
 {sfi =  y/sqrt(x*x+y*y);
 cfi =  x/sqrt(x*x+y*y);} 

  int l,m;   
 // cnst is the Zlm prefactors (plm) - put them into the matrix (same constants are used in jjjpar.cpp and ionpars.cpp)
 Matrix cnst(0,6,-6,6); set_zlm_constants(cnst);

 // evaluate the Zlm in order to get gamma_lm ... in the following lines Zlm(Omega_i) is evaluated
 gamma(0)= cnst(0, 0);
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
  //... now gamma_lm=Zlm(Omega_i)


 // this is squaring of the coefficients of Zlm, a technical trick in
 // order to save a multiplication later (good for the Blm)
 for(l=0;l<=6;l+=2){for(m=-l;m<=l;++m)cnst(l,m)*=cnst(l,m);}

 //ro = a(0, 0) / sqrt(4.0 * 3.1415);

 //evaluate th Zlm in order to get Blm
 B(0)= cnst(0, 0);
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
// ... now B()=plm*Zlm(Omegai)


 // now calculation of the coefficients gammaLM  in cgs
 int i;
 double eps0=8.854187817e-12; //units C^2/Nm^2
 double echarge=1.60217646e-19;  // units C
 //gamma00   
 B(0)*=q0/r*4*PI; gamma(0)*=q0*echarge*1e10/r/eps0;
 //gamma2M   
 for (i=1;i<=5;++i){B(i)*=q2/r/r/r*4*PI/5; gamma(i)*=q2*echarge*1e30/r/r/r/5/eps0;}
 //gamma4M
 for (i=13;i<=21;++i){B(i)*=q4/r/r/r/r/r*4*PI/9; gamma(i)*=q4*echarge*1e50/r/r/r/r/r/9/eps0;}
 //gamma6M
 for (i=33;i<=45;++i){B(i)*=q6/r/r/r/r/r/r/r*4*PI/13; gamma(i)*=q6*echarge*1e70/r/r/r/r/r/r/r/13/eps0; }

 // ... gammas are calculated gamma()=q*Zlm(Omegai)/[r^(l+1)eps0(2l+1)] in SI units [N m^(2-L-1) /C]
 // ... Blm    are calculated B()=4pi*q*plm*Zlm(Omegai)/[r^(l+1)(2l+1)] in |e|A^(-L-1)

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

 double J2meV=1/1.60217646e-22; // 1 millielectron volt = 1.60217646 � 10-22 joules

 // now calculation of the B_LM  and L_LM in meV
                    Blm(0)=-B(0)*e*e*(*iops).nof_electrons*umr;// printf("B(%i)=%g sum(B)=%g\n",0,B(0),(*iops).Blm(0));                  
                    Llm(0)=-echarge*gamma(0)*sqrt(1.0/4.0/PI)*J2meV;//printf("gamma(%i)=%g Llm=%g\n",0,gamma(0),(*iops).Llm(0));  
 for (i=1;i<=5;++i){Blm(i)=-B(i)*e*e*(*iops).r2*(*iops).alpha*ehv2; // printf("B(%i)=%g sum(B)=%g\n",i,B(i),(*iops).Blm(i));
                   if(i!=3){Llm(i)=-echarge*(*iops).r2*a0*a0*1e-20*gamma(i)*sqrt(5.0/8.0/PI)*J2meV;
                           }  //m<>0
                   else    {Llm(i)=-echarge*(*iops).r2*a0*a0*1e-20*gamma(i)*sqrt(5.0/4.0/PI)*J2meV;
                           }  //m=0
                  }
 for (i=13;i<=21;++i){Blm(i)=-B(i)*e*e*(*iops).r4*(*iops).beta*ehv4; 
                      
                   if(i!=17){Llm(i)=-echarge*(*iops).r4*a0*a0*a0*a0*1e-40*gamma(i)*sqrt(9.0/8.0/PI)*J2meV;
                            }  //m<>0
                   else     {Llm(i)=-echarge*(*iops).r4*a0*a0*a0*a0*1e-40*gamma(i)*sqrt(9.0/4.0/PI)*J2meV;
                            }  //m=0
                    }
 for (i=33;i<=45;++i){Blm(i)=-B(i)*e*e*(*iops).r6*ehv6*(*iops).gamma;
                            
                   if(i!=39){Llm(i)=-echarge*(*iops).r6*a0*a0*a0*a0*a0*a0*1e-60*gamma(i)*J2meV*sqrt(13.0/8.0/PI);
                            }  //m<>0
                   else     {Llm(i)=-echarge*(*iops).r6*a0*a0*a0*a0*a0*a0*1e-60*gamma(i)*J2meV*sqrt(13.0/4.0/PI);
                            }  //m=0
                    } 
// now Llm=-|e|q*<r^l> sqrt((2l+1)/4pi) Zlm(Omegai)/[r^(l+1)eps0(2l+1)] in SI units [N m=J] transformed to meV by J2mEV
//     Ll0=-|e|q*<r^l> sqrt((2l+1)/8pi) ...  for m=0

// now Blm=-e^2 <r^l> theta_l 4pi*q*|plm|*Zlm(Omegai)/[r^(l+1)(2l+1)]
}

/**********************************************************************/
// main program
int main (int argc, char **argv)
{
FILE * table_file;
FILE * conv_file;
FILE * sipf_file;
FILE * dBlm_file;
FILE * dLlm_file;
char instr[MAXNOFCHARINLINE];
char module[MAXNOFCHARINLINE];
int i,n=0,ac=0,acold=-1,batchmode=0,omit_pc=0;
float invalues[100];invalues[0]=99;
  double q2,q4,q6,x,y,z;int do_deriv=0;

// treat options
if(argc>1+ac)
 {while(acold!=ac&&argc>1+ac){acold=ac;
if(strcmp(argv[1+ac],"-d")==0) {ac++;do_deriv=1;dBlm_file=fopen_errchk("results/pointc.dBlm","w");
                                        dLlm_file=fopen_errchk("results/pointc.dLlm","w");} // calculate derivatives of Blm Llm
else if(strcmp(argv[1+ac],"-o")==0) {ac++;omit_pc=1;} 
else if(strcmp(argv[1+ac],"-b")==0) {ac++;batchmode=1;}
                                }
 }

// check command line
  if (argc - ac < 3 && batchmode==0)
    { printf (
"\n Program to calculate Crystal field Parameters from Point Charges \n\n"
" Usage: pointc [options] ionname|sipffile  charge_and_position|file.pos\n\n"
"example 1: pointc  Ce3+ 0.2 4 1 5.3\n"
" ... calculate Blm (Stevens Parameters) and Llm (Wybourne Parameters) for one\n"
"     pointcharge of +0.2|e| in distance x=4 A y=1 A z=5.3 A from a Ce3+ ion.\n"
"example 2: pointc Ce3+ file.pos\n"
" ... read several charges+coordinates from file.pos,file format:\n"
"      column 1=charge, column 2-4 = x y z coordinate. (note,progam makenn\n"
"      creates useful files for this option from the crystal structure).\n"
"example 3: pointc Ce3+ C2.pos 5 6\n"
" ... same as example 2 but reduced charge model,i.e. B2m calculated, with\n"
"     charges in col 1 of C2.pos,B4m and B6m with charges in col 5 and 6, respectively.\n"
"example 4: pointc file.sipf 0.2 4 1 5.3 0.1 0.3\n"
" ... read ion from file.sipf,use 0.2|e| for B2m, 0.1|e| for B4m, 0.3|e| for B6m\n"
" ... the first line of the single ion property file.sipf must be\n"
" #!MODULE=so1ion\n"
" # file.sipf should contain the following information (# denotes comments):\n"
" # the name of the ion\n"
" IONTYPE=Ce3+\n"
" #Stevens parameters (optional, necessary for output of Blm)\n"
" ALPHA=-0.0571429 BETA=0.00634921 GAMMA=0\n"
" # the radial matrix elements RN=<r^N> in units of a0^N (a0=0.5292 A)\n"
" R2=1.309 R4=3.964 R6=23.31\n"
" #optional radial wave function parameters, for transition metal ions the the values\n"
" #are tabulated in Clementi & Roetti Atomic data and nuclear data tables 14 \n"
" #(1974) 177-478, the radial wave function is expanded as\n"
" # R(r)=sum_p Cp r^(Np-1) exp(-XIp r)(2 XIp)^(Np+0.5)/sqrt(2Np!)\n"
" #rare earth:Freeman&Watson PR127(1962)2058,Sovers J.Phys.Chem.Sol.28(1966)1073\n"
" #e.g. Co2+ is isoelectronic to Fe+, looking at page 422\n"
" #of Clemente & Roetti the parameters are \n"
" N1=3 XI1=4.95296 C1=0.36301 \n"
" N2=3 XI2=12.2963 C2=0.02707 \n"
" N3=3 XI3=7.03565 C3=0.14777\n"
" \n"
"OUTPUT:  stdout ... sipf file with CEF pars,radial matrix elements,Stevens factors\n"
"                ... pointcharges and positions (omit with option -o)\n"
"         results\\pointc.out ...contains results of convergence when summing up\n"
"                 contributions of different neighbours one by one...\n"
"         results\\pointc.Blm ... Crystal field parameters Blm in Stevens Notation\n" 
"         results\\pointc.Llm ... Crystal field parameters Llm in Wybourne Notation\n" 
"         results\\pointc.dBlm .. for option -d ... derivatives dBlm/du\n" 
"         results\\pointc.dLlm .. for option -d ... derivatives dLlm/du\n" 
"                 ... derivatives are with respect to u=neighborposition(A)/a0\n"
);
      exit (1);
    } else {if(batchmode==0){fprintf(stderr,"#* pointc 221011 *\n");}}

// set Stevens parameters and Lande factor, J and <r^l> of ion
// read create class object ionpars from iontype - sets J, gJ, Stevens factors from the
// routine getpar in cfieldrout.c, thus takes the single ion parameters from
// the same source as the cfield program ...
 ionpars * iops;
 jjjpar * jjjps;
 char *token;
 int sipf_read_module=-1; // -1 no sipf read,  0 module ic1ion,  4 module so1ion
 if((sipf_file=fopen(argv[1+ac],"r"))) //read ion parameters from file
 { fclose(sipf_file);sipf_file=open_sipf(argv[1+ac],module);fprintf(stderr,"\n");
 if(strstr(module,"so1ion")!=NULL){sipf_read_module=4;
   printf("#!MODULE=so1ion\n");}
 if(strstr(module,"ic1ion")!=NULL){sipf_read_module=0;
   printf("#!MODULE=ic1ion\n");
   iops=new ionpars(2);
   while(feof(sipf_file)==false)
   {if(fgets(instr, MAXNOFCHARINLINE, sipf_file)){// strip /r (dos line feed) from line if necessary
                                      while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
                                     // printf("%s",instr);
                                    }
   
   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
        extract(instr,"IONTYPE",(*iops).iontype,(size_t)MAXNOFCHARINLINE,1);

        extract(instr,"nof_electrons",(*iops).nof_electrons);
        extract(instr,"ALPHA",(*iops).alpha);
        extract(instr,"BETA",(*iops).beta);
        extract(instr,"GAMMA",(*iops).gamma);

        extract(instr,"R2",  (*iops).r2);
        extract(instr,"R4",  (*iops).r4);
        extract(instr,"R6",  (*iops).r6);
 
        }
  }
  }  // fi module is ic1ion
  fclose(sipf_file);
  jjjps=new jjjpar(0,0,0,argv[1+ac],1);
  if((*jjjps).module_type!=sipf_read_module){fprintf(stderr,"ERROR pointc: sipf file %s does not start with '#!MODULE=so1ion' or '#!MODULE=ic1ion' !\n",argv[1+ac]);exit(1);}  
  
  if(sipf_read_module==4){iops=(*jjjps).iops;}

      if((*iops).r2==1e300||(*iops).r2==0){(*jjjps).r2_from_radial_wavefunction();
                        printf("#<r^2>  from radial wavefunction in units of a0^2 a0=0.5292 Angstroem\nR2=%g\n",(*jjjps).r2);
                        (*iops).r2=(*jjjps).r2;}
      if((*iops).r4==1e300||(*iops).r4==0){(*jjjps).r4_from_radial_wavefunction();
                        printf("#<r^4>  from radial wavefunction in units of a0^4 a0=0.5292 Angstroem\nR4=%g\n",(*jjjps).r4);
                        (*iops).r4=(*jjjps).r4;}
      if((*iops).r6==1e300||(*iops).r6==0){(*jjjps).r6_from_radial_wavefunction();
                        printf("#<r^6>  from radial wavefunction in units of a0^6 a0=0.5292 Angstroem\nR6=%g\n",(*jjjps).r6);
                        (*iops).r6=(*jjjps).r6;}
(*iops).Blm=0;
(*iops).Llm=0;
if(sipf_read_module==4){(*jjjps).save_sipf(stdout);}
else
{sipf_file=fopen_errchk(argv[1+ac],"rb");double dummy;
           while(feof(sipf_file)==false){fgets(instr, MAXNOFCHARINLINE,sipf_file);
                      // strip /r (dos line feed) from line if necessary
                      while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
                      if(extract(instr,"R2",dummy)==0){sprintf(instr,"R2=%g\n",(*jjjps).r2);}
                      if(extract(instr,"R4",dummy)==0){sprintf(instr,"R4=%g\n",(*jjjps).r4);}
                      if(extract(instr,"R6",dummy)==0){sprintf(instr,"R6=%g\n",(*jjjps).r6);}

                      // remove Blm and Llm and pointcharges from input file 
                      const char lm[]="B00 B22SB21SB20 B21 B22 B33SB32SB31SB30 B31 B32 B33 B44SB43SB42SB41SB40 B41 B42 B43 B44 B55SB54SB53SB52SB51SB50 B51 B52 B53 B54 B55 B66SB65SB64SB63SB62SB61SB60 B61 B62 B63 B64 B65 B66 ";
                      char lm4[5];
                      for(i=0;i<=45;++i){strncpy(lm4,lm+i*4,4);if(lm4[3]!='S')lm4[3]='\0';
                                        if(extract(instr,lm4,dummy)==0)sprintf(instr,"");
                                        lm4[0]='L';if(extract(instr,lm4,dummy)==0)sprintf(instr,"");
                                        }
                     if(extract(instr,"pointcharge",dummy)==0){sprintf(instr,"");}

                      printf("%s",instr);
                                    }
           fclose(sipf_file);

}
 }
 else if(batchmode==1)  // Batch mode - do not write any files; read sipf,pcfile from stdin, prints to stdout
 {
    char ionname[100];
    batchmode = 1;
    fgets(ionname,sizeof(ionname),stdin);
    iops=new ionpars(ionname);
    (*iops).Blm=0;
    (*iops).Llm=0;
    while (n==0 && feof(stdin)==false) n=inputparline("nof_electrons", stdin, invalues); (*iops).nof_electrons = invalues[1]; n=0;
    while (n==0 && feof(stdin)==false) n=inputparline("R2", stdin, invalues);            (*iops).r2 = invalues[1]; n=0;
    while (n==0 && feof(stdin)==false) n=inputparline("R4", stdin, invalues);            (*iops).r4 = invalues[1]; n=0;
    while (n==0 && feof(stdin)==false) n=inputparline("R6", stdin, invalues);            (*iops).r6 = invalues[1]; n=0;
    while (n==0 && feof(stdin)==false) n=inputline(stdin, invalues);
    q2=invalues[1];q4=q2;q6=q2;
    x=invalues[2];
    y=invalues[3];
    z=invalues[4]; 
    if(n>4)q4=invalues[5];
    if(n>5)q6=invalues[6];
    table_file = stdin; 
 }
 else
 {printf ("#!MODULE=so1ion\n#<!--mcphase.sipf-->\n");
  printf("#*********************=*************************************************\n");
  printf("# program pointc - crystal field parameters by the pointcharge model\n");
  printf("# Martin Rotter, %s\n",MCPHASVERSION);
  printf("# Reference: Ernst Bauer and Martin Rotter - Crystal field effects \n");
  printf("#            in Rare Earth Compounds, in print\n");
  printf("#***********************************************************************\n");
  iops=new ionpars(argv[1+ac]);  // read ion parameters from internal table  
//  printf ("IONTYPE=%s\n",(*iops).iontype);
// printout the information used in pointc to output 
//  printf("#Stevens factors\nALPHA=%4g\nBETA=%4g\nGAMMA=%4g\n",(*iops).alpha,(*iops).beta,(*iops).gamma);
//  printf("#Expectation values of radial wave function <r^k> in units of a0^k a0=0.5292 Angstroem\n");
//  printf("R2=%4g\nR4=%4g\nR6=%4g\n\n",(*iops).r2,(*iops).r4,(*iops).r6);
(*iops).Blm=0;
(*iops).Llm=0;
(*iops).save(stdout);  
 printf("#J=%4g\n",(*iops).J);
 printf("#Lande Factor gJ\n GJ = %4g\n",(*iops).gJ);
 }

if(!batchmode) {
// zero parameters in case initialisation put some values to the parameters ...
conv_file=fopen_errchk("results/pointc.out","w");
fprintf(conv_file,"#c0=c2 c4 c6 (|e|) x y z r (A) B00 L00 B22S L22S B21S L21S B20 L20 B21 L21 B22 L22 B44S L44S ... B66 L66\n");
}
if (do_deriv){fprintf(dBlm_file,"#x y z(A) dB00/dux dB00/duy dB00/duz dB22S/dux dB22S/duy dB22S/duz dB21S/du ... dB20/du... dB22/du... dB66/duz\n");
              fprintf(dLlm_file,"#x y z(A) dL00/dux dL00/duy dL00/duz dL22S/dux dL22S/duy dL22S/duz dL21S/du ... dL20/du... dL22/du... dL66/duz\n");}
Vector dBlm0x(0,45),dLlm0x(0,45);dBlm0x=0;dLlm0x=0;
Vector dBlm0y(0,45),dLlm0y(0,45);dBlm0y=0;dLlm0y=0;
Vector dBlm0z(0,45),dLlm0z(0,45);dBlm0z=0;dLlm0z=0;
int q4col=1,q6col=1;
 
if(!batchmode) {
if (argc<6+ac) // read pointcharges from file
{if(argc>3+ac)q4col=(int)strtod(argv[3+ac],NULL);
 if(argc>4+ac)q6col=(int)strtod(argv[4+ac],NULL);
 table_file=fopen_errchk(argv[2+ac],"r");
 while(n==0&&feof(table_file)==false)n=inputline(table_file, invalues);
  q2=invalues[1];q4=q2;q6=q2;
    x=invalues[2];
    y=invalues[3];
    z=invalues[4]; 
    if(q4col>1&&n>=q4col)q4=invalues[q4col];
    if(q6col>1&&n>=q6col)q6=invalues[q6col];
} else 
{ n=4;
  q2=strtod(argv[2+ac],NULL);q4=q2;q6=q2;
  x=strtod(argv[3+ac],NULL);
  y=strtod(argv[4+ac],NULL);
  z=strtod(argv[5+ac],NULL);
  if(argc>6+ac){n=5;q4=strtod(argv[6+ac],NULL);}
  if(argc>7+ac){n=6;q6=strtod(argv[7+ac],NULL);}
}
}
 // print information about pointcharges to file and calculate Blms and Llms
  if(!omit_pc){if(n==4)printf ("\n#pointcharges charge[|e|]  x[A] y[A] z[A]\n");
               else printf ("\n#pointcharges c2[|e|]  x[A] y[A] z[A] c4[|e|] c6[|e|] c0=c2\n");
              }
while(n>0)
{

 if(!omit_pc){if(n==4)printf ("pointcharge= %4g         %4g %4g %4g\n",q2,x,y,z);
               else printf ("pointcharge= %4g         %4g %4g %4g   %4g  %4g\n",q2,x,y,z,q4,q6);
              }
 // calculate Blm's and Llm's
 double r;
 Vector Blm(0,45),Llm(0,45);
 r = sqrt(x * x + y * y + z * z);
 if(!batchmode)fprintf (conv_file," %4g %4g %4g   %4g %4g %4g  %4g  ",q2,q4,q6,x,y,z,r);
 // calculate Blm Llm for this neighbour
 calcCEFpar(q2,q2,q4,q6,x,y,z,r,Blm,Llm,iops);

 // sum to iops.Blm and iops.Lllm
                    (*iops).Blm(0)+=Blm(0); if(!batchmode) fprintf (conv_file,"%g ",Blm(0));
                    (*iops).Llm(0)+=Llm(0); if(!batchmode) fprintf (conv_file,"%g ",Llm(0));
 for (i=1;i<=5;++i){(*iops).Blm(i)+=Blm(i); if(!batchmode) fprintf (conv_file,"%g ",Blm(i));
                    (*iops).Llm(i)+=Llm(i); if(!batchmode) fprintf (conv_file,"%g ",Llm(i));
                   }
 for (i=13;i<=21;++i){(*iops).Blm(i)+=Blm(i); if(!batchmode) fprintf (conv_file,"%g ",Blm(i));
                    (*iops).Llm(i)+=Llm(i); if(!batchmode) fprintf (conv_file,"%g ",Llm(i));
                   }
 for (i=33;i<=45;++i){(*iops).Blm(i)+=Blm(i); if(!batchmode) fprintf (conv_file,"%g ",Blm(i));
                    (*iops).Llm(i)+=Llm(i); if(!batchmode) fprintf (conv_file,"%g ",Llm(i));
                   }
 if(!batchmode)fprintf(conv_file,"\n");

 if(do_deriv){Vector dBlmx(0,45),dLlmx(0,45);
              Vector dBlmy(0,45),dLlmy(0,45);
              Vector dBlmz(0,45),dLlmz(0,45);
              double x1=x+SMALL_DISPLACEMENT;              
              double y1=y+SMALL_DISPLACEMENT;              
              double z1=z+SMALL_DISPLACEMENT;  
              double a0 = .5292;//(Angstroem)            
               r = sqrt(x1 * x1 + y * y + z * z);calcCEFpar(q2,q2,q4,q6,x1,y,z,r,dBlmx,dLlmx,iops);
               dBlmx-=Blm;dLlmx-=Llm;
               dBlmx*=a0/SMALL_DISPLACEMENT;dLlmx*=a0/SMALL_DISPLACEMENT;
               r = sqrt(x * x + y1 * y1 + z * z);calcCEFpar(q2,q2,q4,q6,x,y1,z,r,dBlmy,dLlmy,iops);
               dBlmy-=Blm;dLlmy-=Llm;
               dBlmy*=a0/SMALL_DISPLACEMENT;dLlmy*=a0/SMALL_DISPLACEMENT;
               r = sqrt(x * x + y * y + z1 * z1);calcCEFpar(q2,q2,q4,q6,x,y,z1,r,dBlmz,dLlmz,iops);
               dBlmz-=Blm;dLlmz-=Llm;
               dBlmz*=a0/SMALL_DISPLACEMENT;dLlmz*=a0/SMALL_DISPLACEMENT;
               dBlm0x-=dBlmx;dBlm0y-=dBlmy;dBlm0z-=dBlmz;
               dLlm0x-=dLlmx;dLlm0y-=dLlmy;dLlm0z-=dLlmz;
              fprintf (dBlm_file,"%4g %4g %4g   ",x,y,z);
              fprintf (dLlm_file,"%4g %4g %4g   ",x,y,z);
              fprintf(dBlm_file,"%4g %4g %4g ",dBlmx(0),dBlmy(0),dBlmz(0));
              fprintf(dLlm_file,"%4g %4g %4g ",dLlmx(0),dLlmy(0),dLlmz(0));
              for(i=1;i<=5;++i){fprintf(dBlm_file,"%4g %4g %4g ",dBlmx(i),dBlmy(i),dBlmz(i));
                                fprintf(dLlm_file,"%4g %4g %4g ",dLlmx(i),dLlmy(i),dLlmz(i));
                                }
              for(i=13;i<=21;++i){fprintf(dBlm_file,"%4g %4g %4g ",dBlmx(i),dBlmy(i),dBlmz(i));
                                fprintf(dLlm_file,"%4g %4g %4g ",dLlmx(i),dLlmy(i),dLlmz(i));
                                }
              for(i=33;i<=45;++i){fprintf(dBlm_file,"%4g %4g %4g ",dBlmx(i),dBlmy(i),dBlmz(i));
                                fprintf(dLlm_file,"%4g %4g %4g ",dLlmx(i),dLlmy(i),dLlmz(i));
                                }
             fprintf(dBlm_file,"\n");
             fprintf(dLlm_file,"\n");

             }

 

 n=0;
 if(!batchmode) {
 if (argc<6+ac)
 { while((n==0)&(feof(table_file)==false))n=inputline(table_file, invalues);
  q2=invalues[1];q4=q2;q6=q2;
    x=invalues[2];
    y=invalues[3];
    z=invalues[4]; 
    if(q4col>1&&n>=q4col)q4=invalues[q4col];
    if(q6col>1&&n>=q6col)q6=invalues[q6col];
 }
 } else {
   n=inputline(table_file, invalues);
    q2=invalues[1];q4=q2;q6=q2;
    x=invalues[2];
    y=invalues[3];
    z=invalues[4]; 
    if(n>4)q4=invalues[5];
    if(n>5)q6=invalues[6];
 }
} // next pointcharge

if(!batchmode) {
if (argc<6+ac){fclose(table_file);}
fclose(conv_file);
}
if(do_deriv){
              fprintf (dBlm_file,"0.0 0.0 0.0   ");
              fprintf (dLlm_file,"0.0 0.0 0.0   ");
              fprintf(dBlm_file,"%4g %4g %4g ",dBlm0x(0),dBlm0y(0),dBlm0z(0));
              fprintf(dLlm_file,"%4g %4g %4g ",dLlm0x(0),dLlm0y(0),dLlm0z(0));
              for(i=1;i<=5;++i){fprintf(dBlm_file,"%4g %4g %4g ",dBlm0x(i),dBlm0y(i),dBlm0z(i));
                                fprintf(dLlm_file,"%4g %4g %4g ",dLlm0x(i),dLlm0y(i),dLlm0z(i));
                                }
              for(i=13;i<=21;++i){fprintf(dBlm_file,"%4g %4g %4g ",dBlm0x(i),dBlm0y(i),dBlm0z(i));
                                  fprintf(dLlm_file,"%4g %4g %4g ",dLlm0x(i),dLlm0y(i),dLlm0z(i));
                                }
              for(i=33;i<=45;++i){fprintf(dBlm_file,"%4g %4g %4g ",dBlm0x(i),dBlm0y(i),dBlm0z(i));
                                  fprintf(dLlm_file,"%4g %4g %4g ",dLlm0x(i),dLlm0y(i),dLlm0z(i));
                                }
             fprintf(dBlm_file,"\n");
             fprintf(dLlm_file,"\n");

             fclose(dBlm_file);fclose(dLlm_file);}
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

 if(batchmode) return 0;

table_file=fopen_errchk("./results/pointc.Blm","w");
fprintf(table_file,"#-------------------------------------------------------\n");
fprintf(table_file,"# Crystal Field parameters Blm in Stevens Notation (meV)\n");
fprintf(table_file,"#-------------------------------------------------------\n");
(*iops).savBlm(table_file);
fclose(table_file);
table_file=fopen_errchk("./results/pointc.Llm","w");
fprintf(table_file,"#--------------------------------------------------------\n");
fprintf(table_file,"# Crystal Field parameters Llm in Wybourne Notation (meV)\n");
fprintf(table_file,"#--------------------------------------------------------\n");
(*iops).savLlm(table_file);
fclose(table_file);

  fprintf(stderr,"#***********************************************************************\n");
  fprintf(stderr,"#                         end of program pointc\n");
  fprintf(stderr,"# Reference: Ernst Bauer and Martin Rotter - Magnetism of Complex\n");
  fprintf(stderr,"#            Metallic Alloys: Crystalline Electric Field Effects \n");
  fprintf(stderr,"#            Book Series on Complex Metallic Alloys - Vol. 2, edited\n");
  fprintf(stderr,"#            by Esther Belin-Ferr�, World Scientific, 2009\n");
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
