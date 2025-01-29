/**************************************************************
 * singleion.c - display single ion momentum  at given htpoint
 * Author: Martin Rotter
 **************************************************************/


#include "../../version"
#include "martin.h"
#include "myev.h"
#include<par.hpp>
#include "trs_io.c"   // for in out of trs file
#define HXCMAXDIM 100
// SMALL value of susceptibility in emu treated as zero
#define SMALL_X  1e-5 
/**********************************************************************/
void helpexit()
{ printf (" program single ion  - calculate single ion expectations values <Ia> <Ib> ... \n" 
          "and transition energies at given T and H\n"
          "   use as: singleion [option] T[K] Hexta[T] Hextb[T] Hextc[T] Hxc1 Hxc2 Hxc3 ... Hxcnofcomponents [meV] \n\n"
          "           T    ..... Temperature in Kelvin \n"
          "           Hext ..... external field in Tesla \n"
          "           Hxc... exchange (molecular) field in meV   \n\n"
          "singleion reads mcphas.j and the singleion parameter files quoted therein\n"
          "and calculatesenergies, eigenstates, expectation values <I> for the given\n"
          "temperature, external magnetic field Hext and exchange field Hxc (the\n"
          "interaction constants given in mcphas.j are ignored).\n\n"
          "for each single ion property file the following files are generated:\n"
          "   results/file.sipf.levels.cef .. energy levels and eigenstates and <I>\n"
          "   results/file.sipf.trs ......... transition energies,matrix elements\n"
          "                                   and (powder) neutron intensities\n"
          "   results/_file.sipf    ......... parameters as read by singleion\n"
          "options: -nt ......... by default only 5 transition energies are output,\n"
          "                       if you want more, start e.g. with \n" 
          "                       option -nt 7 to output 7 transition energies\n"
          "         -pinit 0.1 .. consider only transitions with population of initial state > 0.1\n"
          "         -ninit 3  ... consider only transitions from the 3 lowest eigenstates\n"
          "         -maxE 30  ... consider only transitions with energy lower than 30 meV\n"
          "         -E        ... output to stdout energy of cf levels instead of transition energy\n"
          "         -r ion.sipf . do not read mcphas-j but only the single ion\n"
          "                       parameter file ion.sipf\n"
          "         -M  ......... calculate expectation values and transition matrix\n"
          "                       elements for magnetic moment M (muB)instead of I\n"
          "         -MQ 0 0 1 ...  instead of <I> calculate expectation values, transition matrix elements\n"
          "                       for M(Q=(0 0 1)/A), the Fourier Transform  of magnetic moment density M(r) \n"
          "         -L  ......... calculate expectation values and transition matrix\n"
          "                       elements for orbital momentum L\n" 
          "         -S  ......... calculate expectation values and transition matrix\n"
          "                       elements for spin S\n"
          "         -s 0.3 0.1 .. calculate static magnetic single ion susceptibility X in emu/mol with\n"
          "                       constant offset X0=0.3 emu/mol and molecular field constant lambda=0.1 mol/emu\n"
          "                       1/(X-X0)=(1/Xcf)-lambda. Xcf obtained the same way as option -d 0 0 and\n"
          "                       converting the results to emu/mol by multiplying with factor MU_B MU_B NA/10000=\n"
          "                       = 0.0578838263*0.55848973464 = 0.0323275227902. lambda is only applied if  \n"
          "                       Xcf is diagonal.\n"
          "         -is  0.3 0.1  same as -s but output of the inverse susceptibility Y=X^-1 in mol/emu\n"
          "                       (works only if off diagonal elements of X are zero)\n"
          "         -d 23.4 0.001 calculate dynamic magnetic single ion susceptibility X(omega) and neutron\n"
          "                       scattering function S(Q=0,omega) at hbar omega= E= 23.4 meV \n"
          "                       with imaginary part epsilon=0.001 meV (epsilon>0)\n"
          "         -t  ......... reads .trs files from previous run (possibly modified by user, i.e. removing\n"
          "                       some lines to speed up an approximate calculation of spectra.\n"
          "         -Esteps 10 27 for option -d in addition to initial Energy calculate 10 further\n"
          "                       Energies until 27 meV has been reached\n"
          "         -Tsteps 10 27 in addition to initial temperature calculate 10 further temperatures\n"
          "                       until 27K has been reached\n"
          "         -Hsteps 20 0 0 10 in addition to initial field calculate 20 further external fields\n"
          "                       until (0 0 10) Tesla has been reached\n"
          "         -opmat 2 .... output operator matrix number n=2 to results/*.opmat\n"
          "                Operators in results/output op.mat for different values of n:\n"
          "                n=0                    Hamiltonian\n"
          "                n=1,...,nofcomponents  operator Matrix In in standard basis\n"
          "                n=-1,..,-nofomponents  operator Matrix In for Hamiltonian eigenstates basis\n"
          "                n>nofomponents: all operator Matrices (n=0 to n=nofcomponents) in standard basis\n"
          "                n<-nofomponents: all operator Matrices (n=0 to n=-nofcomponents) in Hamiltonian eigenstates basis\n"
          "         -v       .... verbose, output more information on ongoing calculation\n"
          "Note: for calculating T,H dependencies you can instead of using options -Tsteps or -Hsteps put singleion in a LOOP\n"
          "      and pipe the result into a file\n"
          " ... LOOP linux:   for B in $(seq 0 0.1 14); do singleion 2 $B 0 0 0 0 0; done > results/fielddep.dat\n"
          " ... LOOP linux using perl:\n"
          "perl -e 'for($B=1;$B<14;$B+=0.1){system(\"singleion 2 \".$B.\" 0 0  0 0 0\");}' > results/sus1Tesla.clc \n"
          " ... LOOP for windows using perl:\n"
          "perl -e \"for($B=1;$B<14;$B+=0.1){system('singleion 2 '.$B.' 0 0  0 0 0');}\" > results\\sus1Tesla.clc\n"
          );
      exit (1);
}


void colheader(char observable,int observable_nofcomponents,int nofcomponents,Vector & Q,int elevels,double X0,double lambda)
{int j;
   switch(observable){case 's':printf("# Single ion susceptibility X with 1/(X-X0)=(1/Xcf)-lambda, X0=%g emu/mol, lambda=%g mol/emu.\n# For polycrystal Xpoly=Trace(X)/3\n",X0,lambda);break;
                      case 'i':printf("# Inverse Y of single ion susceptibility X with 1/(X-X0)=(1/Xcf)-lambda: Y=X^(-1), X0=%g emu/mol, lambda=%g mol/emu.\n#For polycrystal Xpoly=Trace(X)/3 \n",X0,lambda);break;
                      default: break;
                     }
   printf("# 1         2     ");for(j=1;j<=3;++j)printf("   %i     ",2+j);
                                   for(j=1;j<=nofcomponents;++j)printf("  %2i      ",j+5);
                                   switch(observable)
                                   {case 'Q': printf("                                 ");
                                              for(j=1;j<=observable_nofcomponents;++j){printf("   %2i          %2i          %2i       %2i    ",5+nofcomponents+(j-1)*4+1,5+nofcomponents+(j-1)*4+2,5+nofcomponents+(j-1)*4+3,5+nofcomponents+(j-1)*4+4);}printf("    ");break;
                                    case 's': printf("   %i         ",5+nofcomponents+1);
                                              for(j=2;j<=observable_nofcomponents;++j){printf("%i  ",5+nofcomponents+j);}printf("    ");break;
                                    case 'i': printf("   %i            ",5+nofcomponents+1);
                                              for(j=2;j<=observable_nofcomponents;++j){printf("%i  ",5+nofcomponents+j);}printf("    ");break;
                                    case 'd': printf(" %i    %i                        %i      %i ",5+nofcomponents+1,5+nofcomponents+2,5+nofcomponents+3,5+nofcomponents+4);
                                              for(j=5;j<=observable_nofcomponents;++j){printf("%i   ",5+nofcomponents+j);}printf("    ");break;
                                    default: for(j=1;j<=observable_nofcomponents;++j){printf("  %i  ",5+nofcomponents+j);}printf("    ");break;
                                   }printf("\n");
 printf("#atom-nr   T[K]   ");for(j=1;j<=3;++j)printf("Hext%c(T) ",'a'-1+j);
                                   for(j=1;j<=nofcomponents;++j)printf("Hxc%i(meV) ",j);
                                   switch(observable)
                                   {case 'Q': printf("Q=(%8.5f %8.5f %8.5f)/A ",Q(1),Q(2),Q(3));
                                              for(j=1;j<=observable_nofcomponents;++j){printf(" |<M%c%c>| real(<M%c%c>) imag(<M%c%c>) <M%c>f(Q) ",observable,'a'-1+j,observable,'a'-1+j,observable,'a'-1+j,'a'-1+j);}printf("(muB)");break;
                                    case 'M': for(j=1;j<=observable_nofcomponents;++j)printf(" <%c%c> ",observable,'a'-1+j);printf("(muB)");break;
                                    case 's': printf("Xpoly(emu/mol) X11 X22 X33 X23 X32 X13 X31 X12 X21(emu/mol) ");break;
                                    case 'i': printf("1/Xpoly(mol/emu) Y11 Y22 Y33 Y23 Y32 Y13 Y31 Y12 Y21(mol/emu) ");break;
                                    case 'd': printf("E(meV) Sdip(Q=0,Omega)(barn/meV) Xpolyr Xpolyi(mb^2/meV) X11r X11i X22r X22i X33r X33i X23r X23i X32r X32i X13r X13i X31r X31i X12r X12i X21r X21i(mb^2/meV) ");break;
                                    default: for(j=1;j<=observable_nofcomponents;++j)printf(" <%c%c> ",observable,'a'-1+j);
                                   }
   if(!elevels){printf("transition-energies(meV)...\n");}else{printf("energy levels(meV)...\n");}
if(observable=='d'){printf("#(*)The unpolarized powder average neutron cross section sigma for each transition\n\
#   is calculated neglecting the formfactor, the Debye Wallerfactor, factor k'/k\n");
}
}

void write_trs_file(jjjpar &jjj,int nmax,double pinit,double ninit,double maxE,double TT,Vector & Hext,Vector & Hxc,Vector & Q,char observable,int i)
        {char filename[MAXNOFCHARINLINE];char * pchr;FILE * fout_trs;
         int nt=0;float d=1e10;
  snprintf(filename,MAXNOFCHARINLINE,"./results/%s.trs",jjj.sipffilename);
// if sipffilename contains path (e.g. "./" or "./../")
// do some substitutions to avoid opening error
 pchr=strstr(filename+10,"/");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"/");}
pchr=strstr(filename+10,"\\");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"\\");}

        fout_trs = fopen_errchk (filename,"w");
        trs_header_out(fout_trs,pinit,ninit,maxE,TT,Hext,observable);

        jjj.maxE=maxE;jjj.pinit=pinit;jjj.ninit=ninit;
        jjj.transitionnumber=0;int tc=0;nt=0;
        if(trs_write_next_line(fout_trs,jjj,nt,1,1,1,i,tc,TT,Hxc,Hext,jjj.eigenstates(Hxc,Hext,TT),d,-1e100,maxE,observable,Q))
        {fprintf(stderr,"Warning singleion: no transition found within energy in range [minE,maxE]=[%g,%g]\n"
                        " please increase energy range in option -maxE \n",0.0,maxE);
        }
        else
        {
         while(tc<nmax&&!trs_write_next_line(fout_trs,jjj,nt,1,1,1,i,tc,TT,Hxc,Hext,
                          jjj.est,d,-1e100,maxE,observable,Q)){if(d>=0)--tc;}
        }
          fclose(fout_trs);
      }
       
void read_trs_file(jjjpar & jjj,ComplexMatrix ** X,double Estart,double dE,int Esteps,int elevels,int & jmin,char * trsstring,double epsilon,
 double T,Vector & Hxc,Vector & Hext,int i,int verbose,char observable,int Ti,double maxE)
{ float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
    char * pchr;    char filename[MAXNOFCHARINLINE];FILE * fin;
Vector qijk(1,3);int qcounter=1;
snprintf(filename,MAXNOFCHARINLINE,"./results/%s.trs",jjj.sipffilename);
// if sipffilename contains path (e.g. "./" or "./../")
// do some substitutions to avoid opening error
 pchr=strstr(filename+10,"/");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"/");}
pchr=strstr(filename+10,"\\");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"\\");}
      if(verbose==1)printf("#");
     // load transitions from file
      fin = fopen_errchk(filename,"rb");
      // clear X matrices (put sign of qcounter negative
      switch(observable){case 'd':
                         case 's':
                         case 'i': jjj.eigenstates(Hxc,Hext,T); // ... this is to recalculate population numbers for new temperature
                                   jjj.chi0(X,Estart, dE,Esteps,0.1,qijk,-qcounter,nn[6],T,Hxc,Hext, jjj.est,1,1,1,i);
                                   break;
                         default: break;
                        }
       int i1=0,j1=0;jmin=0;
       while (feof(fin)==0)
       {if ((i1=inputline(fin,nn))>=6)
       {int tn=(int)nn[5];if(nn[6]>=-SMALL_QUASIELASTIC_ENERGY){++jmin;
                        if(Ti==1&&!elevels){ snprintf(trsstring,MAXNOFCHARINLINE,"%s%4g ",trsstring,nn[6]);}
                                                               }
        // calculate delta(single ion excitation energy), 
        // Malphabeta(transition matrix elements)
//      fprintf(stdout,"#transition %i of ion %i of cryst. unit cell at pos  %i %i %i in mag unit cell:\n",tn,l,i,j,k);
//      if(nn[6]<SMALL_QUASIELASTIC_ENERGY){fprintf(stdout,"#-");}else{fprintf(stdout,"#+");}
         j1=jjj.transitionnumber; // try calculation for transition  j
         jjj.transitionnumber=tn; // try calculation for transition  tn
         if(verbose==1)jjj.transitionnumber=-tn;
         // fill X matrices
      switch(observable){case 'd':
                         case 's':
                         case 'i':jjj.chi0(X,Estart, dE,Esteps,-epsilon,qijk,qcounter,nn[6],T,Hxc,Hext, jjj.est,1,1,1,i);
                                 break;
                         default: break;
                        }
           jjj.transitionnumber=j1; // put back transition number for 1st transition
        }} if(jmin==0){fprintf(stderr,"Warning singleion reading %s: no transition found within energy in range [minE,maxE]=[%g,%g] found\n"
                        " (within first crystallographic unit of magnetic unit cell)\n"
                        " please increase energy range in option -maxE \n",filename,0.0,maxE);
            }
          fclose(fin);
if(verbose==1)printf("\n");
        }


void do_a_sipf(jjjpar & jjj,int nmax,double pinit,double ninit,double maxE,Vector & Hext,Vector & Hxc,
              Vector & Q,char observable,int observable_nofcomponents,int nofcomponents,int i,int elevels,
              double Tstart,int Tsteps,Vector & T,double TT,
              double Estart,int Esteps,Vector & E,double dE,
              Vector & Hstart,int Hsteps,Vector & dH,
              double epsilon,double lambda, double X0,int verbose,double opmat,int no_trs_write)
  { char filename[MAXNOFCHARINLINE],trsstring[MAXNOFCHARINLINE];
    float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
    char  * pchr;int j;
    Matrix I(1,observable_nofcomponents,1,Tsteps);
    Vector lnz(1,Tsteps),u(1,Tsteps);
    // transition matrix Mij
    ComplexVector Mq(1,observable_nofcomponents),u1(1,nofcomponents);
    ComplexMatrix MMq(1,observable_nofcomponents,1,Tsteps);
    FILE * fout, * fout_opmat;

   ComplexMatrix **Xcf;
  
   switch(observable)
                  {case 'd':
                   case 's':
                   case 'i':
                   Xcf=new ComplexMatrix *[Esteps+1]; if(Xcf==NULL)exit(EXIT_FAILURE);
                   for(int Ei=0;Ei<Esteps;++Ei){Xcf[Ei]=new  ComplexMatrix(1,3,1,3);
                                                if(Xcf[Ei]==NULL)exit(EXIT_FAILURE);}
                    break;
                  default: break;
                  }
   jjj.Icalc_parameter_storage_init(Hxc,Hext,Tstart);

 if(nmax>0&&no_trs_write==0)write_trs_file(jjj,nmax,pinit,ninit,maxE,TT,Hext,Hxc,Q,observable,i); // write transition trs files
     
      for(int Hi=0;Hi<=Hsteps;++Hi){Hext=Hstart+(double)Hi*dH;
      switch(observable)
      {case 'L': jjj.Lcalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case 'S': jjj.Scalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case 'M': jjj.mcalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case 'Q': for(int Ti=1;Ti<=Tsteps;++Ti)
                 {Vector II(I.Column(Ti));
                 jjj.mcalc(II,T(Ti),Hxc,Hext,jjj.Icalc_parstorage);
                 I.Column(Ti)=II;
                 jjj.eigenstates(Hxc,Hext,T(Ti));
                 jjj.MQ(Mq, Q);
                 for(int ii=1;ii<=observable_nofcomponents;++ii){MMq(ii,Ti)=Mq(ii);}
                 }
                 break;       
       default: jjj.Icalc(I,T,Hxc,Hext,lnz,u,jjj.Icalc_parstorage);
      }  
     
for(int Ti=1;Ti<=Tsteps;++Ti){
    int jmin=0;  

    if(nmax>0)read_trs_file(jjj,Xcf,Estart,dE,Esteps,elevels,jmin,trsstring,epsilon,T(Ti),Hxc,Hext,i,verbose,observable,Ti,maxE);

for(int Ei=0;Ei<Esteps;++Ei){
       
printf("%3i %8g ",i,T(Ti)); // printout ion number and temperature
      for(j=1;j<=3;++j)printf(" %8g ",Hext(j)); // printout external field as requested
      for(j=1;j<=nofcomponents;++j)printf("%8g ",Hxc(j)); // printoutexchangefield as requested
      complex<double> im(0,1.0);
   complex<double> z(E(Ei+1),epsilon);
   complex<double> bose;double S;
   Matrix iX(1,3,1,3);Matrix X(1,3,1,3);
      switch(observable)
       {case 'Q': for(j=1;j<=observable_nofcomponents;++j)printf("%4g %4g %4g %4g   ",abs(MMq(j,Ti)),real(MMq(j,Ti)),imag(MMq(j,Ti)),I(j,Ti)*jjj.F(Norm(Q)));break;
        case 'd': 
		 bose=1.0/(1.0-exp(-z*(1.0/KB/T(Ti))));
		   S=abs(bose/(im)*Trace((*Xcf[Ei])-(*Xcf[Ei]).Transpose().Conjugate()))*2/3/PI/8.0*3.65/4.0/PI;

		printf("%4g %4g %4g %4g  %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g ",E(Ei+1),S,
		real(Trace((*Xcf[Ei])))/3,imag(Trace((*Xcf[Ei])))/3,
		real((*Xcf[Ei])(1,1)),imag((*Xcf[Ei])(1,1)),
		real((*Xcf[Ei])(2,2)),imag((*Xcf[Ei])(2,2)),
		real((*Xcf[Ei])(3,3)),imag((*Xcf[Ei])(3,3)),
		real((*Xcf[Ei])(2,3)),imag((*Xcf[Ei])(2,3)),
		real((*Xcf[Ei])(3,2)),imag((*Xcf[Ei])(3,2)),
		real((*Xcf[Ei])(1,3)),imag((*Xcf[Ei])(1,3)),
		real((*Xcf[Ei])(3,1)),imag((*Xcf[Ei])(3,1)),
		real((*Xcf[Ei])(1,2)),imag((*Xcf[Ei])(1,2)),
		real((*Xcf[Ei])(2,1)),imag((*Xcf[Ei])(2,1))
		      );break;
       case 's': // transform X from mb^2/meV to emu/mol unit
                 // 1. transform to mb/T using Bohr Magneton MU_B = 0.0578838263 meV/Tesla
                 // 2. transform to emu/mol by multiplying with MU_B NA/10000=0.55848973464 
                 (*Xcf[Ei])*=MU_B*0.55848973464;
                 // 1/(X-X0)=(1/Xcf)-lambda.
                 X=Real((*Xcf[Ei]));
                 if(lambda!=0){
                 if(fabs(X(1,2))<SMALL_X&&fabs(X(1,3))<SMALL_X&&fabs(X(2,3))<SMALL_X)
                 {X(1,1)=X(1,1)/(1-X(1,1)*lambda);
                  X(2,2)=X(2,2)/(1-X(2,2)*lambda);
                  X(3,3)=X(3,3)/(1-X(3,3)*lambda);
                 }else
                 {fprintf(stderr,"#Warning: lambda not applied because X is not diagonal\n");}
                 // apply lambda only if matrix is diagonal and nonzero
                 }
                 X+=X0;
                 printf("%4g %4g %4g %4g  %4g %4g %4g %4g %4g %4g ",
                 Trace(X)/3,
		 X(1,1),
		 X(2,2),
		 X(3,3),
		 X(2,3),
		 X(3,2),
		 X(1,3),
		 X(3,1),
		 X(1,2),
		 X(2,1));
                 break;
       case 'i': // transform X from mb^2/meV to emu/mol unit
                 // 1. transform to mb/T using Bohr Magneton MU_B = 0.0578838263 meV/Tesla
                 // 2. transform to emu/mol by multiplying with MU_B NA/10000=0.55848973464 
                 (*Xcf[Ei])*=MU_B*0.55848973464;
                 // 1/(X-X0)=(1/Xcf)-lambda.
                 X=Real((*Xcf[Ei]));
                 if(fabs(X(1,2))<SMALL_X&&fabs(X(1,3))<SMALL_X&&fabs(X(2,3))<SMALL_X)
                 {X(1,1)=X(1,1)/(1-X(1,1)*lambda);
                  X(2,2)=X(2,2)/(1-X(2,2)*lambda);
                  X(3,3)=X(3,3)/(1-X(3,3)*lambda);
                 X+=X0;
                 printf("%4g %4g %4g %4g ",
                 3/Trace(X),
		 1/X(1,1),
		 1/X(2,2),
		 1/X(3,3));
                 }else
                 {fprintf(stderr,"#Warning: X is not diagonal - no inverse susceptibility calculated\n");}
                 break;
       default: for(j=1;j<=observable_nofcomponents;++j)printf("%4g ",I(j,Ti));  // printout corresponding moments      
       } 

       if(nmax>0)
       {if(Ti==1&&Ei==0){
          if(!elevels){printf("%s",trsstring);if(nmax<jmin){printf(" ...");}}
          else
          {for(j=jjj.est.Clo();j<=jjj.est.Chi();++j){printf("%4g ",real(jjj.est(0,j)));}
          }
                 } // fi Ti==1
       } // fi nmax>0
      printf("\n");
    


      }}} // Ei,Ti,Hi

// create levels.cef file   ******************************************
      snprintf(filename,MAXNOFCHARINLINE,"./results/%s.levels.cef",jjj.sipffilename);
// if sipffilename contains path (e.g. "./" or "./../")
// do some substitutions to avoid opening error
 pchr=strstr(filename+10,"/");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"/");}
pchr=strstr(filename+10,"\\");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"\\");}

      fout=fopen_errchk(filename,"w"); 
     fprintf(fout,"#\n#\n#!d=%i sipffile=%s T= %g K ",jjj.est.Chi(),jjj.sipffilename,TT);
                                   for(j=1;j<=3;++j)fprintf(fout,"Hext%c=%g T ",'a'-1+j,Hext(j));
                                   for(j=1;j<=nofcomponents;++j)fprintf(fout,"Hxc%i=%g meV  ",j,Hxc(j));
                                   switch(observable)
                                   {case 'Q': fprintf(fout,"Q=(%g %g %g)/A ",Q(1),Q(2),Q(3));
                                              for(j=1;j<=observable_nofcomponents;++j){fprintf(fout," M%c%c=%g%+gi ",observable,'a'-1+j,real(MMq(j,1)),imag(MMq(j,1)));}fprintf(fout,"(muB) ");break;
                                    case 'M': for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," %c%c=%g ",observable,'a'-1+j,I(j,1));fprintf(fout,"(muB) ");break;
                                    case 'd': 
                                    case 's': 
                                    case 'i': break; 
                                    default: for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," %c%c=%g ",observable,'a'-1+j,I(j,1));
                                   }
                                   fprintf(fout,"\n");jjj.print_eigenstates(fout);fclose(fout);
 
// continue writing op.mat file   ******************************************  
if(opmat<1e10){
     snprintf(filename,MAXNOFCHARINLINE,"./results/%s.opmat",jjj.sipffilename);
// if sipffilename contains path (e.g. "./" or "./../")
// do some substitutions to avoid opening error
 pchr=strstr(filename+10,"/");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"/");}
pchr=strstr(filename+10,"\\");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"\\");}

      fout_opmat=fopen_errchk(filename,"w"); 

 fprintf(fout_opmat,"#! d=%i  ",jjj.est.Chi());
                    if(opmat>nofcomponents){ for(int opmati=0;opmati<=nofcomponents;++opmati)
                                             {Matrix op(jjj.opmat(opmati,Hxc,Hext));
                                             myPrintComplexMatrix(fout_opmat,op);}

                                           }
                    else
                    { if(opmat<-nofcomponents){

                                                Matrix opp(jjj.opmat(0,Hxc,Hext)); opp=0;
                                               for (int opmati=1;opmati<=jjj.est.Chi();++opmati)
                                               {opp(opmati,opmati)=real(jjj.est(0,opmati));}
                                                myPrintComplexMatrix(fout_opmat,opp);
                                             
                                             for(int opmati=-1;opmati>=-nofcomponents;--opmati)
                                             {Matrix op(jjj.opmat(opmati,Hxc,Hext));
                                             myPrintComplexMatrix(fout_opmat,op);}
                                              }
                     else
                     {Matrix op(jjj.opmat((int)opmat,Hxc,Hext));
                      myPrintComplexMatrix(fout_opmat,op);
                     }
                    }
           fclose(fout_opmat);             
      
     }
 switch(observable){
         case 'd':
         case 's':
         case 'i': for(int Ei=0;Ei<Esteps;++Ei)delete Xcf[Ei];
                   if(Xcf!=NULL)delete []Xcf;break;
         default: break;
                    }
 }// fi i

//***************************************************************************************
//***************************************************************************************


// hauptprogramm
int main (int argc, char **argv)
{ int i,j,do_sipf=0,verbose=0;
   double ninit=100000000,pinit=0,maxE=1e10,opmat=1e10,Estart=0,epsilon=0,dE=0,X0=0,lambda=0;
   int Tsteps=0,Hsteps=0,Esteps=0,elevels=0,no_trs_write=0;
   double Eend=0,Tend=0,Tstart=0;
   Vector Hend(1,3),Hstart(1,3);
  int nofcomponents=0;
  Vector Hext(1,3),Q(1,3),Hxc_in(1,HXCMAXDIM);
  char sipffile[MAXNOFCHARINLINE];

  int nmax=5;// default number of transitions to  be output
  char observable='I'; // default is operators I
printf("#***singleion.c - calculate single ion properties - M. Rotter %s*****\n",MCPHASVERSION);
//***************************************************************************************
// check command line parameters 
//***************************************************************************************
for (i=1;i<argc;++i)
 {if(strncmp(argv[i],"-h",2)==0) {helpexit();}
  else {if(strcmp(argv[i],"-M")==0) observable='M';       
  else {if(strcmp(argv[i],"-MQ")==0){observable='Q';
                                      if(i==argc-1){fprintf(stderr,"Error in command: singleion -MQ needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  Q(1)=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -MQ needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  Q(2)=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -MQ needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  Q(3)=strtod(argv[i+1],NULL);++i;
    			            }      
  else {if(strcmp(argv[i],"-S")==0) {observable='S'; }      
  else {if(strcmp(argv[i],"-d")==0) {observable='d'; 
                                      if(i==argc-1){fprintf(stderr,"Error in command: singleion -d needs arguments E and epsilon\n");exit(EXIT_FAILURE);}
	                                  Estart=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -d needs arguments E and epsilon\n");exit(EXIT_FAILURE);}
	                                  epsilon=strtod(argv[i+1],NULL);++i;
                                    }      
  else {if(strcmp(argv[i],"-s")==0) {observable='s'; 
                                      if(i==argc-1){fprintf(stderr,"Error in command: singleion -s needs arguments X0 and lambda\n");exit(EXIT_FAILURE);}
	                                  X0=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -s needs arguments X0 and lambda\n");exit(EXIT_FAILURE);}
	                                  lambda=strtod(argv[i+1],NULL);++i;
                                    }      
  else {if(strcmp(argv[i],"-is")==0) {observable='i'; 
                                      if(i==argc-1){fprintf(stderr,"Error in command: singleion -s needs arguments X0 and lambda\n");exit(EXIT_FAILURE);}
	                                  X0=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -s needs arguments X0 and lambda\n");exit(EXIT_FAILURE);}
	                                  lambda=strtod(argv[i+1],NULL);++i;
                                    }      
  else {if(strcmp(argv[i],"-L")==0) {observable='L'; }      
  else {if(strcmp(argv[i],"-t")==0) {no_trs_write=1; }      
  else {if(strcmp(argv[i],"-nt")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -nt needs argument\n");exit(EXIT_FAILURE);}
	                                  nmax=(int)strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-pinit")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -pinit needs argument\n");exit(EXIT_FAILURE);}
	                                  pinit=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-ninit")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -ninit needs argument\n");exit(EXIT_FAILURE);}
	                                  ninit=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-maxE")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -maxE needs argument\n");exit(EXIT_FAILURE);}
	                                  maxE=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-E")==0) {elevels=1;}       
  else {if(strcmp(argv[i],"-r")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -r needs argument\n");exit(EXIT_FAILURE);}
	                              do_sipf=1;strcpy(sipffile,argv[i+1]);++i;
    			             }       
  else {if(strcmp(argv[i],"-opmat")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -opmat needs argument\n");exit(EXIT_FAILURE);}
	                              opmat=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-Esteps")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -Esteps needs arguments\n");exit(EXIT_FAILURE);}
	                              Esteps=(int)fabs(strtod(argv[i+1],NULL));++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -Esteps needs 2 arguments\n");exit(EXIT_FAILURE);}
	                              Eend=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-Tsteps")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -Tsteps needs argument(s)\n");exit(EXIT_FAILURE);}
	                              Tsteps=(int)fabs(strtod(argv[i+1],NULL));++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -Tsteps needs 2 arguments\n");exit(EXIT_FAILURE);}
	                              Tend=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-Hsteps")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs arguments\n");exit(EXIT_FAILURE);}
	                              Hsteps=(int)fabs(strtod(argv[i+1],NULL));++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs 4 arguments\n");exit(EXIT_FAILURE);}
	                              Hend(1)=strtod(argv[i+1],NULL);++i;
    			             if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs 4 arguments\n");exit(EXIT_FAILURE);}
	                              Hend(2)=strtod(argv[i+1],NULL);++i;
    			             if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs 4 arguments\n");exit(EXIT_FAILURE);}
	                              Hend(3)=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-v")==0)verbose=1;      
  else{Tstart=strtod(argv[i],NULL);++i; if(!Tsteps){Tend=Tstart;} // now read T
       Hext=0;for(j=1;j<=3;++j){if(i<argc){Hext(j)=strtod(argv[i],NULL);}++i;} // read Hexta Hextb Hextc
       Hxc_in=0;for(j=1;i<argc&&j<HXCMAXDIM;++j){++nofcomponents;Hxc_in(j)=strtod(argv[i],NULL);++i;} //read Hxc1 Hxc2 ... Hxcn
      } // T Hext Hxc
    } // verbose
    } // -Hsteps
    } // -Tsteps
    } // -Esteps
    } // -opmat
    } // -r
    } //-maxE         
    } //-E      
    } //-ninit         
    } //-pinit         
    } //-nt  
    } // -t       
    } // -L  
    } // -is 
    } // -s 
    } // -d 
    } // -S 
   } // -MQ  
   } // -M  
 } // help
  if(argc<2){helpexit();}
  if(nofcomponents==0){fprintf(stdout,"ERROR singleion: please enter exchange field Hxc\n");exit(EXIT_FAILURE);}
  if(epsilon<0){fprintf(stdout,"ERROR singleion option -s: epsilon has to be >0\n");exit(EXIT_FAILURE);}
  if(Esteps!=0&&observable!='d'){fprintf(stdout,"ERROR singleion option Esteps makes only sense for dynamical susceptibility option -d\n");exit(EXIT_FAILURE);}
  Vector Hxc(1,nofcomponents);Hxc=0;for(j=1;j<=nofcomponents;++j)Hxc(j)=Hxc_in(j);

  int observable_nofcomponents;
  switch(observable)
   {case 'M':
    case 'Q':
    case 'S':
    case 'L': observable_nofcomponents=3;break;
    case 's': observable_nofcomponents=10;break; // chipoly 3x3 matrix of static susceptibility
    case 'i': observable_nofcomponents=10;break; // inverse chipoly 3x3 matrix of static susceptibility
    case 'd': observable_nofcomponents=22;
              break; // Energy transfer neutron cross section complex polycrystal susceptibility and 3x3x2
                    // components of complex susceptibility tensor
    default: observable_nofcomponents=nofcomponents; // I
   }
double EE=Estart;
if(Estart>Eend){Estart=Eend;Eend=EE;}
++Esteps;Vector E(1,Esteps); // Esteps= number of Energies to calculate for option -s
if(Esteps>1){dE=(Eend-Estart)/(Esteps-1);}
E(1)=Estart;for(int Ei=1;Ei<Esteps;++Ei){E(Ei+1)=E(Ei)+(Eend-Estart)/(Esteps-1);} //set E's
double TT=Tend;
 if(Tend<Tstart){TT=Tstart;Tstart=Tend;Tend=TT;} 
// always calculate ascending temperatures.
//  TT is used to write trs file - do this for highest temperature (to use ninit and pinit for this)
    
++Tsteps;Vector T(1,Tsteps); // Tsteps= number of temperatures to calculate
T(1)=Tstart;for(int Ti=1;Ti<Tsteps;++Ti){T(Ti+1)=T(Ti)+(Tend-Tstart)/(Tsteps-1);} //set T's
Vector dH(1,3);dH=0;
Hstart=Hext;if(Hsteps){dH=Hend-Hstart;dH*=(1.0/Hsteps);}
//myPrintVector(stdout,dH);printf("%i\n",Hsteps);exit(0);
 
if (!do_sipf)
  {par inputpars("./mcphas.j",verbose);
   inputpars.save_sipfs("./results/_");
   if(nofcomponents!=inputpars.cs.nofcomponents)fprintf(stderr,"#Warning: number of exchange field components read from command line not equal to that in mcphas.j - continuing...\n");
    colheader(observable,observable_nofcomponents,nofcomponents,Q,elevels,X0,lambda);
    
                 

  for(i=1;i<=inputpars.cs.nofatoms;++i)
   { do_a_sipf((*inputpars.jjj[i]),nmax,pinit,ninit,maxE,Hext,Hxc,Q,
              observable,observable_nofcomponents,nofcomponents,i,elevels,
              Tstart,Tsteps,T,TT,
              Estart,Esteps,E,dE,
              Hstart,Hsteps,dH,epsilon,lambda,X0,verbose,opmat,no_trs_write);
   }

  
  } else { // option -r sipffile
   jjjpar jjj(0,0,0,sipffile,nofcomponents,verbose);jjj.save_sipf("./results/_");
   colheader(observable,observable_nofcomponents,nofcomponents,Q,elevels,X0,lambda);

   do_a_sipf(jjj,nmax,pinit,ninit,maxE,Hext,Hxc,Q,
              observable,observable_nofcomponents,nofcomponents,1,elevels,
              Tstart,Tsteps,T,TT,
              Estart,Esteps,E,dE,
              Hstart,Hsteps,dH,epsilon,lambda,X0,verbose,opmat,no_trs_write);
            
fprintf(stderr,"# **********************end of program singleion************************\n");
if(verbose)fprintf(stderr,"# ... you can now use 'cpsingleion' to calculate specific heat,\n"
       "#      entropy etc from results/*.levels.cef\n"
       "# **********************************************************************\n");
  }
}
/*
    if(opmat<1e10){fout_opmat=fopen_errchk("./results/op.mat","w");}
    jjj.Icalc_parameter_storage_init(Hxc,Hext,Tstart);

                 switch(observable)
                  {case 'd':
                   case 's':
                   case 'i':
                   Xcf=new ComplexMatrix ** [2];if(Xcf==NULL)exit(EXIT_FAILURE);
                   Xcf[1]=new ComplexMatrix*[Esteps+1];  
                   for(int Ei=0;Ei<Esteps;++Ei){Xcf[1][Ei]=new ComplexMatrix(1,3,1,3);
                                               if(Xcf[1][Ei]==NULL)exit(EXIT_FAILURE);}
                   break;
                  default: break;
                  }


if(nmax>0)write_trs_file(jjj,nmax,pinit,ninit,maxE,TT,Hext,Hxc,Q,observable,1);

   for(int Hi=0;Hi<=Hsteps;++Hi){Hext=Hstart+(double)Hi*dH;
        switch(observable)
      {case 'L': jjj.Lcalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case 'S': jjj.Scalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case 'M': jjj.mcalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case 'Q': for(int Ti=1;Ti<=Tsteps;++Ti)
                 {Vector II(I.Column(Ti));jjj.mcalc(II,T(Ti),Hxc,Hext,jjj.Icalc_parstorage);I.Column(Ti)=II;
                 jjj.eigenstates(Hxc,Hext,T(Ti));
                 jjj.MQ(Mq, Q);
                 for(int ii=1;ii<=observable_nofcomponents;++ii){MMq(ii,Ti)=Mq(ii);}
                 }
                 break;       
      default: jjj.Icalc(I,T,Hxc,Hext,lnz,u,jjj.Icalc_parstorage);
      }          

for(int Ti=1;Ti<=Tsteps;++Ti){
int jmin=0;

      if(nmax>0)read_trs_file(jjj,X,Estart,dE,Esteps,elevels,jmin,trsstring,epsilon,T(Ti),Hxc,Hext,1,verbose,observable,Ti,maxE);
 

 for(int Ei=0;Ei<Esteps;++Ei){
     printf("%3i %8g ",1,T(Ti)); // printout ion number and temperature
      for(j=1;j<=3;++j)printf(" %10g ",Hext(j)); // printout external field as requested
      for(j=1;j<=nofcomponents;++j)printf("%8g ",Hxc(j)); // printoutexchangefield as requested
   complex<double> im(0,1.0);
   complex<double> z(E(Ei+1),epsilon);
   complex<double> bose;double S;
      switch(observable)
       {case 'Q': for(j=1;j<=observable_nofcomponents;++j)printf("%4g %4g %4g %4g   ",abs(MMq(j,Ti)),real(MMq(j,Ti)),imag(MMq(j,Ti)),I(j,Ti)*jjj.F(Norm(Q)));break;
        case 'd': 
   bose=1.0/(1.0-exp(-z*(1.0/KB/T(Ti))));
   S=abs(bose/(im)*Sum((*Xcf[1][Ei])-(*Xcf[1][Ei]).Transpose().Conjugate()))*2/3/PI/8.0*3.65/4.0/PI;
   
printf("%4g %4g %4g %4g  %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g ",E(Ei+1),S,
real(Trace((*Xcf[1][Ei])))/3,imag(Trace((*Xcf[1][Ei])))/3,
real((*Xcf[1][Ei])(1,1)),imag((*Xcf[1][Ei])(1,1)),
real((*Xcf[1][Ei])(2,2)),imag((*Xcf[1][Ei])(2,2)),
real((*Xcf[1][Ei])(3,3)),imag((*Xcf[1][Ei])(3,3)),
real((*Xcf[1][Ei])(2,3)),imag((*Xcf[1][Ei])(2,3)),
real((*Xcf[1][Ei])(3,2)),imag((*Xcf[1][Ei])(3,2)),
real((*Xcf[1][Ei])(1,3)),imag((*Xcf[1][Ei])(1,3)),
real((*Xcf[1][Ei])(3,1)),imag((*Xcf[1][Ei])(3,1)),
real((*Xcf[1][Ei])(1,2)),imag((*Xcf[1][Ei])(1,2)),
real((*Xcf[1][Ei])(2,1)),imag((*Xcf[1][Ei])(2,1))
      );
                  break;
      default: for(j=1;j<=observable_nofcomponents;++j)printf("%4g ",I(j,Ti));  // printout corresponding moments      
       }

    if(nmax>0)
         {
          if(Ti==1&&Ei==0)
          {if(!elevels){printf("%s",trsstring);if(nmax<jmin){printf(" ...");}}
           else
           {for(j=jjj.est.Clo();j<=jjj.est.Chi();++j){printf("%4g ",real(jjj.est(0,j)));}
           }
          }
         } // nmax>0
    printf("\n"); // fi Ti==1
  TT=T(Ti);
   }}} // Ei Ti Hi

 // create levels.cef file   ******************************************
    snprintf(filename,MAXNOFCHARINLINE,"./results/%s.levels.cef",jjj.sipffilename);
// if sipffilename contains path (e.g. "./" or "./../")
// do some substitutions to avoid opening error
 pchr=strstr(filename+10,"/");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"/");}
pchr=strstr(filename+10,"\\");
 while(pchr!=0){memcpy(pchr,"I",1);pchr=strstr(filename+10,"\\");}

      fout=fopen_errchk(filename,"w");  
    fprintf(fout,"#\n#\n#!d=%i sipffile=%s T= %g K ",jjj.est.Chi(),jjj.sipffilename,TT);
                                   for(j=1;j<=3;++j)fprintf(fout,"Hext%c=%g T ",'a'-1+j,Hext(j));
                                   for(j=1;j<=nofcomponents;++j)fprintf(fout,"Hxc%i=%g meV  ",j,Hxc(j));
                                   switch(observable)
                                   {case 'Q': fprintf(fout,"Q=(%g %g %g)/A ",Q(1),Q(2),Q(3));
                                              for(j=1;j<=observable_nofcomponents;++j){fprintf(fout," M%c%c=%g%+gi ",observable,'a'-1+j,real(MMq(j,1)),imag(MMq(j,1)));}fprintf(fout,"(muB)");break;
                                    case 'M': for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," %c%c=%g ",observable,'a'-1+j,I(j,1));fprintf(fout,"(muB)");break;
                                    case 'd': fprintf(fout,"E(meV) Sdip(Q=0,Omega)(barn/meV) Xpolyr Xpolyi(mb^2/meV) X11r X11i X22r X22i X33r X33i X23r X23i X32r X32i X13r X13i X31r X31i X12r X12i X21r X21i(mb^2/meV) ");break;
                                    default: for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," %c%c=%g ",observable,'a'-1+j,I(j,1));
                                   }
                                   fprintf(fout,"\n");jjj.print_eigenstates(fout);fclose(fout);

// continue writing op.mat file   ******************************************
if(opmat<1e10){fout_opmat=fopen_errchk("results/op.mat","w");
                       fprintf(fout_opmat,"#! d=%i  ",jjj.est.Chi());
if(opmat>nofcomponents){ for(int opmati=0;opmati<=nofcomponents;++opmati)
                                             {Matrix op(jjj.opmat(opmati,Hxc,Hext));
                      myPrintComplexMatrix(fout_opmat,op);}

                        }
                    else
                    { if(opmat<-nofcomponents){Matrix opp(jjj.opmat(0,Hxc,Hext)); opp=0;
                                               for (int opmati=1;opmati<=jjj.est.Chi();++opmati)
                                               {opp(opmati,opmati)=real(jjj.est(0,opmati));}
                                                myPrintComplexMatrix(fout_opmat,opp);
                                             for(int opmati=-1;opmati>=-nofcomponents;--opmati)
                                             {Matrix op(jjj.opmat(opmati,Hxc,Hext));
                                              myPrintComplexMatrix(fout_opmat,op);}
                                              }
                     else
                     {Matrix op(jjj.opmat((int)opmat,Hxc,Hext));
                      myPrintComplexMatrix(fout_opmat,op);
                     }
                    }
fclose(fout_opmat);}
 if(observable=='d'){ for(int Ei=0;Ei<Esteps;++Ei)delete Xcf[1][Ei];
                      if(Xcf[1]!=NULL)delete []Xcf[1];
                         delete []X;
                     }
      
   

 */



