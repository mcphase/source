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

enum Xunit { none = 0 , emu = 1, SI =2 , muBT = 3};

/**********************************************************************/
void helpexit()
{ printf (" program single ion  - calculate single ion expectations values <Ia> <Ib> ... \n" 
          "and transition energies at given T and H\n"
          "   use as: singleion [option] T[K] Hexti[T] Hextj[T] Hextk[T] Hxc1 Hxc2 Hxc3 ... Hxcnofcomponents [meV] \n\n"
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
          "         -U  ......... calculate energy U, ln of partition sum Z, free energy F\n"
          "                       instead of <I>\n"
          "         -t  ......... reads .trs files from previous run (possibly modified by user, i.e. removing\n"
          "                       some lines to speed up an approximate calculation of spectra.)\n"
          "         -Esteps 10 27 for option -d in addition to initial Energy calculate 10 further\n"
          "                       Energies until 27 meV has been reached\n"
          "         -Tsteps 10 27 in addition to initial temperature calculate 10 further temperatures\n"
          "                       until 27K has been reached\n"
          "         -Hsteps 20 0 0 10 in addition to initial field calculate 20 further external fields\n"
          "                       until (0 0 10) Tesla has been reached\n"
          "         -HE ......... in addition to magnetic field also apply electrical field in kV/mm, i.e. there will be \n"
          "                       instead of 3 components  Hexti Hextj Hextk in the command line 6 components \n"
          "                       Hexti Hextj Hextk Eexti Eextj Eextk, similar for Hsteps option there will be 6 components\n"
          "                       (therefore: mind that -HE is given before -Hsteps in the command line)\n"
          "         -opmat 2 .... output operator matrix number n=2 to results/*.opmat\n"
          "                Operators in results/output op.mat for different values of n:\n"
          "                n=0                    Hamiltonian\n"
          "                n=1,...,nofcomponents  operator Matrix In in standard basis\n"
          "                n=-1,..,-nofomponents  operator Matrix In for Hamiltonian eigenstates basis\n"
          "                n>nofomponents: all operator Matrices (n=0 to n=nofcomponents) in standard basis\n"
          "                n<-nofomponents: all operator Matrices (n=0 to n=-nofcomponents) in Hamiltonian eigenstates basis\n"
          "         -v       .... verbose, output more information on ongoing calculation\n\n"
          "  Other Observables:\n"
          "         -M  ...magnetic moment: calculate expectation values and transition matrix\n"
          "                       elements for magnetic moment M (muB)instead of I\n"
          "         -P  ...phonon displacement: calculate phonon displacement in A instead of I\n"
          "         -pel ..electric dipole moment: calculate electrical dipole moment pel in |e|pm\n"
          "                       instead of I\n"
          "         -L  ....orbital momentum: calculate expectation values and transition matrix\n"
          "                       elements for orbital momentum L\n" 
          "         -S  ....spin: calculate expectation values and transition matrix\n"
          "                       elements for spin S\n"
          "         -MQ 0 0 1 ...M(Q),Fourier transform of magnetic moment density:  instead of <I>\n"
          "                       calculate expectation values, transition matrix elements\n"
          "                       for M(Q=(0 0 1)/A), the Fourier Transform  of magnetic moment density M(r) \n"
          "         -sx ....spin density: calculate expectation values and transition matrix\n"
          "                       elements for spindensity coefficients aSx(lm) in expansion \n"
          "                       of spindensity-x-component in Ms(r) = sum_lm aS(l,m) R^2(r) Zlm(Omega)\n"
          "                        E. Balcar J. Phys. C. 8 (1975) 1581\n"
          "         -sy -sz ......for y and z components use option -sy and  -sz\n"
          "         -lx ....orbital moment density: calculate expectation values and transition matrix\n"
          "                       elements for orbital moment density coefficients aLx(lm) in expansion \n"
          "                       of orbital moment density-x-component in Ml(r)=sum_lm  aLx(l,m) F(r) Zlm(Omega)\n"
          "                       with F(r)==1/r int_r^inf R^2(x) dx,   E. Balcar J. Phys. C. 8 (1975) 1581\n"
          "         -ly -lz ......for y and z components use option -ly and  -lz\n\n"
          "   Susceptibilities:\n"
          "         -X[observable] 24 0.1  ... calculate susceptibility chi(z) for z=E + i epsilon with\n"
          "                       E=24 meV and epsilon=0.1 meV. The susceptibility  chi is defined as the\n"
          "                       derivative of I (or any observable) with respect to the field, e.g. for\n"
          "                       the observable magnetic moment -XM the susceptibility is chi=dM/dH\n"
          "                       units: standard unit is (unit of observable)^2/meV, e.g. for -XM\n"
          "                              (muB)^2/meV. you can use the options ... \n"
          "         -emu 0.3 0.1 . to obtain the magnetic susceptibility chi in units of emu/mol with\n"
          "                       constant offset X0=0.3 emu/mol and molecular field constant lambda=0.1 mol/emu\n"
          "                       1/(X-X0)=(1/Xcf)-lambda. Xcf obtained the same way as option -d 0 0 and\n"
          "                       converting the results to emu/mol by multiplying with factor MU_B MU_B NA/10000=\n"
          "                       = 0.0578838263*0.55848973464 = 0.0323275227902. lambda is only applied if  \n"
          "                       Xcf is diagonal.\n"
          "                       e.g. -XM 0 0 -emu 0 0 will calculate the static magnetic susceptibility in emu/mol\n"
          "         -SI ......... to obtain the magnetic susceptibility (-XM) in SI units (Am^2/mol)\n"
          "         -muBT ......... to obtain the magnetic susceptibility (-XM) in  (muB/Tesla)\n"
          "         -iX[observable] 24 0.1  ....same as -X, but output inverse susceptibility\n"
          "                       (works only if off diagonal elements of X are zero)\n\n"

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



void write_trs_file(jjjpar &jjj,int nmax,double pinit,double ninit,double maxE,double TT,Vector & Hext,Vector & Hxc,Vector & Q,ob observable,int i,int HEnofcomp )
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
        char outstring[MAXNOFCHARINLINE];
        snprintf(outstring,MAXNOFCHARINLINE," T= %g K Hi=%g Hj=%g Hk=%g T",TT,Hext(1),Hext(2),Hext(3));
       if(HEnofcomp>5)snprintf(outstring+strlen(outstring),MAXNOFCHARINLINE,"  Ei=%g Ej=%g Ek=%g kV/mm",Hext(4),Hext(5),Hext(6));
        trs_header_out(fout_trs,pinit,ninit,maxE,outstring,observable);

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
 double T,Vector & Hxc,Vector & Hext,int i,int verbose,ob observable,int Ti,double maxE,int calcX)
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
      if(calcX){jjj.eigenstates(Hxc,Hext,T); // ... this is to recalculate population numbers for new temperature
                jjj.chi0(X,Estart, dE,Esteps,0.1,qijk,-qcounter,nn[6],T,Hxc,Hext, jjj.est,1,1,1,i);
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
      if(calcX){jjj.chi0(X,Estart, dE,Esteps,-epsilon,qijk,obint(observable),nn[6],T,Hxc,Hext, jjj.est,1,1,1,i);
                                 }
           jjj.transitionnumber=j1; // put back transition number for 1st transition
        }} if(jmin==0){fprintf(stderr,"Warning singleion reading %s: no transition found within energy in range [minE,maxE]=[%g,%g] found\n"
                        " (within first crystallographic unit of magnetic unit cell)\n"
                        " please increase energy range in option -maxE \n",filename,0.0,maxE);
            }
          fclose(fin);
if(verbose==1)printf("\n");
        }

void colheader(ob observable,int observable_nofcomponents,int nofcomponents,Vector & Q,int elevels,double X0,double lambda,int HEnofcomp,int calcX,Xunit unit,int verbose)
{int j;char str[MAXNOFCHARINLINE];
 if(verbose==1&& calcX)
 {printf("# Single Ion Susceptibility X is defined as:\n"
"#                                                                        \n"
"#                    ----    <l|%s-<%s>|j><j|%s-<%s>|l>                  \n"
"#    X(omega+i eps)= >       --------------------------------   w        \n"
"#                    ----          E - E - omega - i eps         lj       \n"
"#                     l,j           l   j                                \n"
"#                                                                        \n"
"#                                                                        \n"
"#      l->j correspond to the nt transitions listed in the file results/*.trs \n"
"#      (compare  option -nt)                                              \n"
"#                                                                         \n"
"#                            with                                        \n"
"#                                                                      \n"
"#      w  = w -w    and w  =w /kT                                      \n"
"#       lj   l  j        ll  l                                        \n"
"#                                                              \n"
"#                      - E /kT                                  \n"
"#                     e   j                                    \n"
"# w     =       --------------                                 \n"
"#  j            ----    - E /kT                                 \n"
"#               >      e   l                                   \n"
"#               ----                                           \n"
"#                l                                             \n"
"#                                                              \n"
"#  units: the default unit is  (unit of observable)^2/meV      \n"
"#                                                              \n"
"#  unit conversions:   the magnetic susceptibility  (options -XM -iXM)     \n"
"#                      has standard unit muB^2/meV and it can be converted   \n"
"#                      to other units with the options -emu -SI -muBT   \n"
"#                     -muBT:  to muB/Tesla express one Bohrmagneton in  \n"
"#                              muB = 0.0578838263 meV/Tesla      \n"
"#                     -emu:   further, to emu/mol by expressing the other MU_B=9.27e-24Am^2\n"
"#                             the unit Tesla=10000 Oe and multiplying by\n"
"#                             the Avogadro number NA=6.022e23/mol, \n"
"#                             i.e. apply conversion factor MuB NA/10000=0.55848973464 \n"
"#                                                              \n"
"#                     -SI:    or further, to m3/mol by inserting MU_B=9.27e-24 Am^2 and\n"
"#                             multiplying by  NA=6.022e23/mol and setting 1Tesla/mu0=A/m \n"
"#                                                              \n"

,obs[obint(observable)]
,obs[obint(observable)]
,obs[obint(observable)]
,obs[obint(observable)]
);
 }

 if(observable==M&& calcX==1&&unit==emu)
  printf("# Single ion susceptibility X with 1/(X-X0)=(1/Xcf)-lambda, X0=%g emu/mol, lambda=%g mol/emu.\n# For polycrystal Xpoly=Trace(X)/3\n",X0,lambda);
 else if(observable==M&& calcX==-1&&unit==emu)
  printf("# Inverse Y of single ion susceptibility X with 1/(X-X0)=(1/Xcf)-lambda: Y=X^(-1), X0=%g emu/mol, lambda=%g mol/emu.\n#For polycrystal Xpoly=Trace(X)/3 \n",X0,lambda);
 else if(observable==pel&& calcX==1)
  printf("# Single ion susceptibility Xij=dpel_i/deps0Ej  Electric field j=1,2,3   \n# For polycrystal Xpoly=Trace(X)/3\n");
 else if(observable==pel&& calcX==-1)
  printf("# Inverse Y of single ion susceptibility Xij=dpel_i/deps0Ej  Electric field j=1,2,3 : Y=X^(-1)\n#For polycrystal Xpoly=Trace(X)/3 \n");
 else if(observable==P&& calcX==1)
  printf("# Single ion susceptibility Xij=dP_i/deps0Ej  Electric field j=1,2,3   \n# For polycrystal Xpoly=Trace(X)/3\n");
 else if(observable==P&& calcX==-1)
  printf("# Inverse Y of single ion susceptibility Xij=dP_i/deps0Ej  Electric field j=1,2,3 : Y=X^(-1)\n#For polycrystal Xpoly=Trace(X)/3 \n");
 else if(calcX==1)
  printf("# Single ion susceptibility Xij=d%s_i/dmu0Hj  Magnetic field j=1,2,3   \n# For polycrystal Xpoly=Trace(X)/3\n",obs[observable]);
 else if(calcX==-1)
  printf("# Inverse Y of single ion susceptibility Xij=d%s_i/dmu0Hj  Magnetic field j=1,2,3 : Y=X^(-1)\n#For polycrystal Xpoly=Trace(X)/3 \n",obs[observable]);
                      


 str[0]='\0';
 snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"#atom-nr   T[K]   ");for(j=1;j<=3;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"Hext%c(T) ",'i'-1+j);
                              if(HEnofcomp>5)for(j=1;j<=3;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"Eext%c(kV/mm) ",'i'-1+j);
                                   for(j=1;j<=nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"Hxc%i(meV) ",j);
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

                                      
if(calcX==0)
{                                  switch(observable)
                                   {case MQ: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"Q=(%8.5f %8.5f %8.5f)/A ",Q(1),Q(2),Q(3));
                                              for(j=1;j<=observable_nofcomponents;++j){snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," |<M%s%c>| real(<M%s%c>) imag(<M%s%c>) <M%c>f(Q) ",obs[observable],'a'-1+j,obs[observable],'a'-1+j,obs[observable],'a'-1+j,'a'-1+j);}snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(muB)");break;
                                    case M: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"  <%s%i>",obs[observable],j);snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(muB)");break;
                                    case P: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"  <%s%i>",obs[observable],j);snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(A)");break;
                                    case pel: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"  <pel%i>",j);snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(|e|pm)");break;
                                    case U: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," lnZ  F(meV) U(meV) ");break;
                                    case sx: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," <aSx(%i,%i)> ",k[j],q[j]);break;
                                    case sy: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," <aSy(%i,%i)> ",k[j],q[j]);break;
                                    case sz: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," <aSz(%i,%i)> ",k[j],q[j]);break;
                                    case lx: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," <aLx(%i,%i)> ",k[j],q[j]);break;
                                    case ly: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," <aLy(%i,%i)> ",k[j],q[j]);break;
                                    case lz: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," <aLz(%i,%i)> ",k[j],q[j]);break;
                                    default: for(j=1;j<=observable_nofcomponents;++j)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," <%s%i> ",obs[observable],j);
                                   }
} 
else if(calcX==1)
{if(observable==M)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"E(meV) Sdip(Q=0,Omega)(barn/meV) Xpolyr Xpolyi X11r X11i X22r X22i X33r X33i X23r X23i X32r X32i X13r X13i X31r X31i X12r X12i X21r X21i");
 else snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"E(meV) Xpolyr Xpolyi X11r X11i X22r X22i X33r X33i X23r X23i X32r X32i X13r X13i X31r X31i X12r X12i X21r X21i");
 switch(unit){ case emu: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(emu/mol)");break;
               case muBT: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(muB/T)");break;
               case SI: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(m^3/mol)");break;
               default: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"[(%s)^2/meV ion]",obunit[observable]);break;
             }
}
else if(calcX==-1)
{if(observable==M)snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"E(meV) Sdip(Q=0,Omega)(barn/meV) 1/Xpolyr 1/Xpolyi Y11r Y11i Y22r Y22i Y33r Y33i Y23r Y23i Y32r Y32i Y13r Y13i Y31r Y31i Y12r Y12i Y21r Y21i");
 else snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"E(meV) 1/Xpolyr 1/Xpolyi Y11r Y11i Y22r Y22i Y33r Y33i Y23r Y23i Y32r Y32i Y13r Y13i Y31r Y31i Y12r Y12i Y21r Y21i");

switch(unit){ case emu: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(mol/emu)");break;
              case muBT: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(T/muB)");break;
               case SI: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"(mol/m^3)");break;
               default: snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"[meVion/(%s)^2]",obunit[observable]);break;
             }

}

                       if(!elevels){snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," transition-energies(meV)...");}else{snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str)," energy levels(meV)...");}
if(observable==M&&calcX){printf("#(*)The unpolarized powder average neutron cross section sigma for each transition\n\
#   is calculated neglecting the formfactor, the Debye Wallerfactor, factor k'/k\n");
                   }

printf("#");print_col_numbers(stdout,str);

 printf("\n%s\n",str);
}

void do_a_sipf(jjjpar & jjj,int nmax,double pinit,double ninit,double maxE,Vector & Hext,Vector & Hxc,
              Vector & Q,ob observable,int observable_nofcomponents,int nofcomponents,int i,int elevels,
              double Tstart,int Tsteps,Vector & T,double TT,
              double Estart,int Esteps,Vector & E,double dE,
              Vector & Hstart,int Hsteps,Vector & dH,
              double epsilon,double lambda, double X00,int verbose,double opmat,int no_trs_write,int HEnofcomp,
              int calcX,Xunit unit)
  { char filename[MAXNOFCHARINLINE],trsstring[MAXNOFCHARINLINE];
    float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
    char  * pchr;int j;
    Matrix I(1,observable_nofcomponents,1,Tsteps);complex <double> X0 (X00,0);
    Vector lnz(1,Tsteps),u(1,Tsteps);
    // transition matrix Mij
    ComplexVector Mq(1,observable_nofcomponents),u1(1,nofcomponents);
    ComplexMatrix MMq(1,observable_nofcomponents,1,Tsteps);
    FILE * fout, * fout_opmat;

   ComplexMatrix **Xcf;
  
   if(calcX){ Xcf=new ComplexMatrix *[Esteps+1]; if(Xcf==NULL)exit(EXIT_FAILURE);
                   for(int Ei=0;Ei<Esteps;++Ei){Xcf[Ei]=new  ComplexMatrix(1,3,1,3);
                                                if(Xcf[Ei]==NULL)exit(EXIT_FAILURE);}
                  }
   jjj.Icalc_parameter_storage_init(Hxc,Hext,Tstart);

 if(nmax>0&&no_trs_write==0)write_trs_file(jjj,nmax,pinit,ninit,maxE,TT,Hext,Hxc,Q,observable,i,HEnofcomp); // write transition trs files
     for(int Hi=0;Hi<=Hsteps;++Hi){Hext=Hstart+(double)Hi*dH;
      switch(observable)
      {case L: jjj.Lcalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case S: jjj.Scalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case M: jjj.mcalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case P: jjj.pcalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case pel: jjj.pelcalc(I,T,Hxc,Hext,jjj.Icalc_parstorage);break;
       case MQ: for(int Ti=1;Ti<=Tsteps;++Ti)
                 {Vector II(I.Column(Ti));
                 jjj.mcalc(II,T(Ti),Hxc,Hext,jjj.Icalc_parstorage);
                 SetColumn(Ti,I,II);
                 jjj.eigenstates(Hxc,Hext,T(Ti));
                 jjj.MQ(Mq, Q);
                 for(int ii=1;ii<=observable_nofcomponents;++ii){MMq(ii,Ti)=Mq(ii);}
                 }
                 break;       
       case sx: jjj.spindensity_coeff (I,-1,T,Hxc,Hext, jjj.Icalc_parstorage);break;
       case sy: jjj.spindensity_coeff (I,-2,T,Hxc,Hext, jjj.Icalc_parstorage);break;
       case sz: jjj.spindensity_coeff (I,-3,T,Hxc,Hext, jjj.Icalc_parstorage);break;
       case lx: jjj.orbmomdensity_coeff (I,-1,T,Hxc,Hext, jjj.Icalc_parstorage);break;
       case ly: jjj.orbmomdensity_coeff (I,-2,T,Hxc,Hext, jjj.Icalc_parstorage);break;
       case lz: jjj.orbmomdensity_coeff (I,-3,T,Hxc,Hext, jjj.Icalc_parstorage);break;
       default: jjj.Icalc(I,T,Hxc,Hext,lnz,u,jjj.Icalc_parstorage);
      }  

for(int Ti=1;Ti<=Tsteps;++Ti){
    int jmin=0;  

    if(nmax>0)read_trs_file(jjj,Xcf,Estart,dE,Esteps,elevels,jmin,trsstring,epsilon,T(Ti),Hxc,Hext,i,verbose,observable,Ti,maxE,calcX);

for(int Ei=0;Ei<Esteps;++Ei){
       
printf("%3i %8g ",i,T(Ti)); // printout ion number and temperature
      for(j=1;j<=3;++j)printf(" %8g ",Hext(j)); // printout external field as requested
      if(HEnofcomp>5)for(j=4;j<=6;++j)printf(" %8g ",Hext(j)); 
      for(j=1;j<=nofcomponents;++j)printf("%8g ",Hxc(j)); // printoutexchangefield as requested
      complex<double> im(0,1.0);
   complex<double> z(E(Ei+1),epsilon);
   complex<double> bose;double S;
   Matrix iX(1,3,1,3);ComplexMatrix X(1,3,1,3);
        if(calcX){   	 bose=1.0/(1.0-exp(-z*(1.0/KB/T(Ti))));
		  if(observable==M)
                     S=abs(bose/(im)*Trace((*Xcf[Ei])-(*Xcf[Ei]).Transpose().Conjugate()))*2/3/PI/8.0*3.65/4.0/PI;
                if(unit==emu)
                 {// transform X from mb^2/meV to emu/mol unit
                 // 1. transform to mb/T using Bohr Magneton MU_B = 0.0578838263 meV/Tesla
                 // 2. transform to emu/mol by multiplying with MU_B NA/10000=0.55848973464 
                 (*Xcf[Ei])*=MU_B*0.55848973464;
                 // 1/(X-X0)=(1/Xcf)-lambda.
                 X=(*Xcf[Ei]);
                 if(lambda!=0){
                 if(abs(X(1,2))<SMALL_X&&abs(X(1,3))<SMALL_X&&abs(X(2,3))<SMALL_X)
                 {complex <double> one(1.0,0);
                  X(1,1)=X(1,1)/(one-X(1,1)*lambda);
                  X(2,2)=X(2,2)/(one-X(2,2)*lambda);
                  X(3,3)=X(3,3)/(one-X(3,3)*lambda);
                 }else
                 {fprintf(stderr,"#Warning: lambda not applied because X is not diagonal\n");}
                 // apply lambda only if matrix is diagonal and nonzero
                 }
                 (*Xcf[Ei])=X+X0;}
                else if(unit==SI)
                 {// transform X from mb^2/meV to m^3/mol unit
                 // 1. transform to mb/T using Bohr Magneton MU_B = 0.0578838263 meV/Tesla
                 // 2. transform to m3/mol by inserting MU_B=9.27e-24 Am^2 and multiplying by 
                 //     NA=6.022e23/mol 1Tesla=A/m 
                 (*Xcf[Ei])*=MU_B*9.27*6.022e-1;}
                else if(unit==muBT)
                 {// transform X from mb^2/meV to m^3/mol unit
                 // 1. transform to mb/T using Bohr Magneton MU_B = 0.0578838263 meV/Tesla
                 (*Xcf[Ei])*=MU_B;}
                
		printf("%4g ",E(Ei+1));
		if(observable==M)printf("%4g ",S);

                if(calcX==-1) // we want inverse susceptibility 
                {printf("%4g %4g ",3/real(Trace((*Xcf[Ei]))),3/imag(Trace((*Xcf[Ei]))));
                (*Xcf[Ei])=(*Xcf[Ei]).Inverse();
                }
                else
                {printf("%4g %4g ",real(Trace((*Xcf[Ei])))/3,imag(Trace((*Xcf[Ei])))/3);
                }
                printf("%4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g ",
		real((*Xcf[Ei])(1,1)),imag((*Xcf[Ei])(1,1)),
		real((*Xcf[Ei])(2,2)),imag((*Xcf[Ei])(2,2)),
		real((*Xcf[Ei])(3,3)),imag((*Xcf[Ei])(3,3)),
		real((*Xcf[Ei])(2,3)),imag((*Xcf[Ei])(2,3)),
		real((*Xcf[Ei])(3,2)),imag((*Xcf[Ei])(3,2)),
		real((*Xcf[Ei])(1,3)),imag((*Xcf[Ei])(1,3)),
		real((*Xcf[Ei])(3,1)),imag((*Xcf[Ei])(3,1)),
		real((*Xcf[Ei])(1,2)),imag((*Xcf[Ei])(1,2)),
		real((*Xcf[Ei])(2,1)),imag((*Xcf[Ei])(2,1))
		      );}

  else   switch(observable)
       {case U: printf("%4g %4g %4g ",lnz(Ti),-KB*T(Ti)*lnz(Ti),u(Ti));break;
        case MQ: for(j=1;j<=observable_nofcomponents;++j)printf("%4g %4g %4g %4g   ",abs(MMq(j,Ti)),real(MMq(j,Ti)),imag(MMq(j,Ti)),I(j,Ti)*jjj.F(Norm(Q)));break;
 
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
                                   if(HEnofcomp>5)for(j=1;j<=3;++j)fprintf(fout,"Eext%c=%g kV/mm ",'a'-1+j,Hext(j+3)); 
                                   for(j=1;j<=nofcomponents;++j)fprintf(fout,"Hxc%i=%g meV  ",j,Hxc(j));
                                   switch(observable)
                                   {case MQ: fprintf(fout,"Q=(%g %g %g)/A ",Q(1),Q(2),Q(3));
                                              for(j=1;j<=observable_nofcomponents;++j){fprintf(fout," M%s%c=%g%+gi ",obs[observable],'a'-1+j,real(MMq(j,1)),imag(MMq(j,1)));}fprintf(fout,"(muB) ");break;
                                    case M: for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," %s%c=%g ",obs[observable],'a'-1+j,I(j,1));fprintf(fout,"(muB) ");break;
                                    case pel: for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," pel%c=%g ",'a'-1+j,I(j,1));fprintf(fout,"(|e|pm) ");break;
                                    default: for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," %s%c=%g ",obs[observable],'a'-1+j,I(j,1));
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
if(calcX){ for(int Ei=0;Ei<Esteps;++Ei)delete Xcf[Ei];
                   if(Xcf!=NULL)delete []Xcf;
                }
 }// fi i

//***************************************************************************************
//***************************************************************************************


// hauptprogramm
int main (int argc, char **argv)
{ int i,j,do_sipf=0,verbose=0;Xunit unit=none;
   double ninit=100000000,pinit=0,maxE=1e10,opmat=1e10,Estart=0,epsilon=0,dE=0,X0=0,lambda=0;
   int Tsteps=0,Hsteps=0,Esteps=0,elevels=0,no_trs_write=0,calcX=0;
   double Eend=0,Tend=0,Tstart=0;
   Vector Hend(1,HEXT_DIMENSION),Hstart(1,HEXT_DIMENSION);
  int nofcomponents=0;
  Vector Hext(1,HEXT_DIMENSION),Q(1,3),Hxc_in(1,HXCMAXDIM);
  char sipffile[MAXNOFCHARINLINE],cmp[10];
  int HEnofcomp=3;
  int nmax=5;// default number of transitions to  be output
  ob observable=I; // default is operators I
printf("#***singleion.c - calculate single ion properties - M. Rotter %s*****\n",MCPHASVERSION);
//***************************************************************************************
// check command line parameters 
//***************************************************************************************
for (i=1;i<argc;++i)
 {for(int j=0;j<NOFOBS;++j){snprintf(cmp,10,"-%s",obs[j]);
                            if(strcmp(argv[i],cmp)==0){observable=obint(j);++i;}
                            snprintf(cmp,10,"-X%s",obs[j]);
                            if(strcmp(argv[i],cmp)==0){observable=obint(j);++i;calcX=1;if(j>6){fprintf(stderr,"Error singleion: option -X not implemented for observable %s\n" ,obs[j]);exit(EXIT_FAILURE);}
                                      if(i==argc-1){fprintf(stderr,"Error in command: singleion -X needs arguments E and epsilon\n");exit(EXIT_FAILURE);}
	                                  Estart=strtod(argv[i],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -X needs arguments E and epsilon\n");exit(EXIT_FAILURE);}
	                                  epsilon=strtod(argv[i],NULL);++i;
                                                     }
                            snprintf(cmp,10,"-iX%s",obs[j]);
                            if(strcmp(argv[i],cmp)==0){observable=obint(j);++i;calcX=-1;if(j>6){fprintf(stderr,"Error singleion: option -iX not implemented for observable %s\n" ,obs[j]);exit(EXIT_FAILURE);}
                                      if(i==argc-1){fprintf(stderr,"Error in command: singleion -X needs arguments E and epsilon\n");exit(EXIT_FAILURE);}
	                                  Estart=strtod(argv[i],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -X needs arguments E and epsilon\n");exit(EXIT_FAILURE);}
	                                  epsilon=strtod(argv[i],NULL);++i;
                                                     }
                           }
  if(strncmp(argv[i],"-h",2)==0) {helpexit();}
  else if(strcmp(argv[i],"-U")==0) observable=U;       
  else if(strcmp(argv[i],"-MQ")==0){observable=MQ;
                                      if(i==argc-1){fprintf(stderr,"Error in command: singleion -MQ needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  Q(1)=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -MQ needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  Q(2)=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -MQ needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  Q(3)=strtod(argv[i+1],NULL);++i;
    			            }         
  else if(strcmp(argv[i],"-emu")==0) {unit=emu; if(i==argc-1){fprintf(stderr,"Error in command: singleion -emu needs arguments X0 and lambda\n");exit(EXIT_FAILURE);}
	                                  X0=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -emu needs arguments X0 and lambda\n");exit(EXIT_FAILURE);}
	                                  lambda=strtod(argv[i+1],NULL);++i;
                                    }      
  else if(strcmp(argv[i],"-SI")==0) {unit=SI; }      
  else if(strcmp(argv[i],"-muBT")==0) {unit=muBT; }      
  else if(strcmp(argv[i],"-t")==0) {no_trs_write=1; }      
  else if(strcmp(argv[i],"-nt")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -nt needs argument\n");exit(EXIT_FAILURE);}
	                                  nmax=(int)strtod(argv[i+1],NULL);++i;
    			             }       
  else if(strcmp(argv[i],"-pinit")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -pinit needs argument\n");exit(EXIT_FAILURE);}
	                                  pinit=strtod(argv[i+1],NULL);++i;
    			             }       
  else if(strcmp(argv[i],"-ninit")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -ninit needs argument\n");exit(EXIT_FAILURE);}
	                                  ninit=strtod(argv[i+1],NULL);++i;
    			             }       
  else if(strcmp(argv[i],"-maxE")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -maxE needs argument\n");exit(EXIT_FAILURE);}
	                                  maxE=strtod(argv[i+1],NULL);++i;
    			             }       
  else if(strcmp(argv[i],"-E")==0) {elevels=1;}       
  else if(strcmp(argv[i],"-r")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -r needs argument\n");exit(EXIT_FAILURE);}
	                              do_sipf=1;strcpy(sipffile,argv[i+1]);++i;
    			             }       
  else if(strcmp(argv[i],"-opmat")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -opmat needs argument\n");exit(EXIT_FAILURE);}
	                              opmat=strtod(argv[i+1],NULL);++i;
    			             }       
  else if(strcmp(argv[i],"-Esteps")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -Esteps needs arguments\n");exit(EXIT_FAILURE);}
	                              Esteps=(int)fabs(strtod(argv[i+1],NULL));++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -Esteps needs 2 arguments\n");exit(EXIT_FAILURE);}
	                              Eend=strtod(argv[i+1],NULL);++i;
    			             }       
  else if(strcmp(argv[i],"-Tsteps")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -Tsteps needs argument(s)\n");exit(EXIT_FAILURE);}
	                              Tsteps=(int)fabs(strtod(argv[i+1],NULL));++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -Tsteps needs 2 arguments\n");exit(EXIT_FAILURE);}
	                              Tend=strtod(argv[i+1],NULL);++i;
    			             }       
  else if(strcmp(argv[i],"-Hsteps")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs arguments\n");exit(EXIT_FAILURE);}
	                              Hsteps=(int)fabs(strtod(argv[i+1],NULL));++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs 4 arguments\n");exit(EXIT_FAILURE);}
	                              Hend(1)=strtod(argv[i+1],NULL);++i;
    			             if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs 4 arguments\n");exit(EXIT_FAILURE);}
	                              Hend(2)=strtod(argv[i+1],NULL);++i;
    			             if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs 4 arguments\n");exit(EXIT_FAILURE);}
	                              Hend(3)=strtod(argv[i+1],NULL);++i;
                                   if(HEnofcomp>5){
                                      if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs 7 arguments\n");exit(EXIT_FAILURE);}
	                              Hend(4)=strtod(argv[i+1],NULL);++i;
    			             if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs 7 arguments\n");exit(EXIT_FAILURE);}
	                              Hend(5)=strtod(argv[i+1],NULL);++i;
    			             if(i==argc-1){fprintf(stderr,"Error in command: singleion -Hsteps needs 7 arguments\n");exit(EXIT_FAILURE);}
	                              Hend(6)=strtod(argv[i+1],NULL);++i;
                                                  }
  			             }       
  else if(strcmp(argv[i],"-HE")==0)HEnofcomp=6;      
  else if(strcmp(argv[i],"-v")==0)verbose=1;      
  else{Tstart=strtod(argv[i],NULL);++i; if(!Tsteps){Tend=Tstart;} // now read T
       Hext=0;for(j=1;j<=HEnofcomp;++j){if(i<argc){Hext(j)=strtod(argv[i],NULL);}++i;} // read Hexta Hextb Hextc and Ei Ej Ek if requested
       Hxc_in=0;for(j=1;i<argc&&j<HXCMAXDIM;++j){++nofcomponents;Hxc_in(j)=strtod(argv[i],NULL);++i;} //read Hxc1 Hxc2 ... Hxcn
      } // T Hext Hxc
   } // next i

 if(argc<2){helpexit();}
  if(nofcomponents==0){fprintf(stdout,"ERROR singleion: please enter exchange field Hxc\n");exit(EXIT_FAILURE);}
  if(epsilon<0){fprintf(stdout,"ERROR singleion option -X, -iX: epsilon has to be >0\n");exit(EXIT_FAILURE);}
  if(Esteps!=0&&calcX==0){fprintf(stdout,"ERROR singleion option Esteps makes only sense for dynamical susceptibility option -X -iX\n");exit(EXIT_FAILURE);}
  if(unit!=none)if(calcX==0||observable!=M){fprintf(stdout,"ERROR singleion options -emu and -SI only make sense for observable M, i.e. options -XM and -iXM\n");exit(EXIT_FAILURE);}

  Vector Hxc(1,nofcomponents);Hxc=0;for(j=1;j<=nofcomponents;++j)Hxc(j)=Hxc_in(j);

  int observable_nofcomponents;
  switch(observable)
   {case U:
    case M:
    case pel:
    case MQ:
    case S:
    case L: observable_nofcomponents=3;break;
    case sx:
    case sy:
    case sz: observable_nofcomponents=ORBMOMDENS_EV_DIM;break; // orbmom and spindensity has 49 coefficients
    case lx:
    case ly:
    case lz: observable_nofcomponents=SPINDENS_EV_DIM;break; // orbmom and spindensity has 49 coefficients
    default: observable_nofcomponents=nofcomponents; // I
   }

    
// for susceptibility tensor we have nxnx2 (real and imag) + 2 (polycrystal) + Energy  components to calculate
if(calcX){observable_nofcomponents=observable_nofcomponents*observable_nofcomponents*2+2+1;
         if(observable==M) ++observable_nofcomponents;}// in case of magnetic moment also calculate neutron cross section 

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
Vector dH(1,HEXT_DIMENSION);dH=0;
Hstart=Hext;if(Hsteps){dH=Hend-Hstart;dH*=(1.0/Hsteps);}
//myPrintVector(stdout,dH);printf("%i\n",Hsteps);exit(0);
 
if (!do_sipf)
  {par inputpars("./mcphas.j",verbose);
   inputpars.save_sipfs("./results/_");
   if(nofcomponents!=inputpars.cs.nofcomponents)fprintf(stderr,"#Warning: number of exchange field components read from command line not equal to that in mcphas.j - continuing...\n");
    colheader(observable,observable_nofcomponents,nofcomponents,Q,elevels,X0,lambda,HEnofcomp,calcX,unit,verbose);
    
                 

  for(i=1;i<=inputpars.cs.nofatoms;++i)
   { do_a_sipf((*inputpars.jjj[i]),nmax,pinit,ninit,maxE,Hext,Hxc,Q,
              observable,observable_nofcomponents,nofcomponents,i,elevels,
              Tstart,Tsteps,T,TT,
              Estart,Esteps,E,dE,
              Hstart,Hsteps,dH,epsilon,lambda,X0,verbose,opmat,no_trs_write,HEnofcomp,calcX,unit);
   }

  
  } else { // option -r sipffile
   jjjpar jjj(0,0,0,sipffile,nofcomponents,verbose);jjj.save_sipf("./results/_");
   colheader(observable,observable_nofcomponents,nofcomponents,Q,elevels,X0,lambda,HEnofcomp,calcX,unit,verbose);

   do_a_sipf(jjj,nmax,pinit,ninit,maxE,Hext,Hxc,Q,
              observable,observable_nofcomponents,nofcomponents,1,elevels,
              Tstart,Tsteps,T,TT,
              Estart,Esteps,E,dE,
              Hstart,Hsteps,dH,epsilon,lambda,X0,verbose,opmat,no_trs_write,HEnofcomp,calcX,unit);
            
fprintf(stderr,"# **********************end of program singleion************************\n");
if(verbose)fprintf(stderr,"# ... you can now use 'cpsingleion' to calculate specific heat,\n"
       "#      entropy etc from results/*.levels.cef\n"
       "# **********************************************************************\n");
  }
}



