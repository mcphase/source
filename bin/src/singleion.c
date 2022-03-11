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

/**********************************************************************/
void helpexit()
{ printf (" program single ion  - display single ion expectations values <Ia> <Ib> ... \n" 
          "and transition energies at given T and H\n"
          "   use as: singleion [option] T[K] Hexta[T] Hextb[T] Hextc[T] Hxc1 Hxc2 Hxc3 ... Hxcnofcomponents [meV] \n\n"
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
          "         -r ion.sipf . do not read mcphas-j but only the single ion\n"
          "                       parameter file ion.sipf\n"
          "         -M  ......... calculate expectation values and transition matrix\n"
          "                       elements for magnetic moment M (muB)instead of I\n"
          "         -MQ ......... calculate expectation values and transition matrix elements\n"
          "                       for Fourier Transform M(Q) of magnetic moment density M(r)  instead of I\n"
          "         -S  ......... calculate expectation values and transition matrix\n"
          "                       elements for spin S\n"
          "         -L  ......... calculate expectation values and transition matrix\n"
          "                       elements for orbital momentum L\n" 
          "         -opmat 2 .... output operator matrix number 2 to results/op.mat\n"
          "    n=0 Hamiltonian, n=1,...,nofcomponents: operator Matrix In in standard basis\n"
          "                     n=-1,..,-nofomponents: operator Matrix In for Hamiltonian eigenstates basis\n"
          "                     n>nofomponents: all operator Matrices (n=0 to n=nofcomponents) in standard basis\n"
          "                     n<-nofomponents: all operator Matrices (n=0 to n=-nofcomponents) in Hamiltonain eigenstates basis\n"
          "Note: for calculating T or H dependencies you can put single ion in a LOOP\n"
          "      and pipe the result into a file\n"
          "  .... LOOP linux:   for B in $(seq 0 0.1 14); do singleion 2 $B 0 0 0 0 0; done > results/fielddep.dat\n"
          " ...  LOOP linux using perl:\n"
          "perl -e 'for($T=1;$T<90;$T+=20){system(\"singleion \".$T.\" 1 0 0  0 0 0\");}' > results/sus1Tesla.clc \n"
          " ... LOOP for windows using perl:\n"
          "perl -e \"for($T=1;$T<90;$T+=20){system('singleion '.$T.' 1 0 0 0 0 0');}\" > results\\sus1Tesla.clc\n"
          );
      exit (1);
}

// hauptprogramm
int main (int argc, char **argv)
{ int i,j,nt=0,do_sipf=0;
  char sipffile[MAXNOFCHARINLINE];char  filename[MAXNOFCHARINLINE];
  double lnz,u; double ninit=100000000,pinit=0,maxE=1e10,opmat=1e10;
  float d=1e10;int nofcomponents=0;FILE * fout,* fout_trs, * fout_opmat;
  double T;Vector Hext(1,3),Q(1,3),Hxc_in(1,HXCMAXDIM);
  int nmax=5;// default number of transitions to  be output
  char observable='I'; // default is operators I
printf("#**************************************************************\n");
printf("# * singleion.c - calculate single ion properties\n");
printf("# * Author: Martin Rotter %s\n",MCPHASVERSION);
printf("# **************************************************************\n");
//***************************************************************************************
// check command line parameters 
//***************************************************************************************
for (i=1;i<argc;++i)
 {if(strncmp(argv[i],"-h",2)==0) {helpexit();}
  else {if(strcmp(argv[i],"-M")==0) observable='M';       
  else {if(strcmp(argv[i],"-MQ")==0) {observable='Q';
                                      if(i==argc-1){fprintf(stderr,"Error in command: singleion -MQ needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  Q(1)=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -MQ needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  Q(2)=strtod(argv[i+1],NULL);++i;
	                              if(i==argc-1){fprintf(stderr,"Error in command: singleion -MQ needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  Q(3)=strtod(argv[i+1],NULL);++i;
    			             }      
  else {if(strcmp(argv[i],"-S")==0) observable='S';       
  else {if(strcmp(argv[i],"-L")==0) observable='L';       
  else {if(strcmp(argv[i],"-nt")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -nt needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  nmax=(int)strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-pinit")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -pinit needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  pinit=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-ninit")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -ninit needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  ninit=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-maxE")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -maxE needs argument(s)\n");exit(EXIT_FAILURE);}
	                                  maxE=strtod(argv[i+1],NULL);++i;
    			             }       
  else {if(strcmp(argv[i],"-r")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -r needs argument(s)\n");exit(EXIT_FAILURE);}
	                              do_sipf=1;strcpy(sipffile,argv[i+1]);++i;
    			             }       
  else {if(strcmp(argv[i],"-opmat")==0) {if(i==argc-1){fprintf(stderr,"Error in command: singleion -opmat needs argument(s)\n");exit(EXIT_FAILURE);}
	                              opmat=strtod(argv[i+1],NULL);++i;
    			             }       
  else{T=strtod(argv[i],NULL);++i;  // now read T
       Hext=0;for(j=1;j<=3;++j){if(i<argc){Hext(j)=strtod(argv[i],NULL);}++i;} // read Hexta Hextb Hextc
       Hxc_in=0;for(j=1;i<argc&&j<HXCMAXDIM;++j){++nofcomponents;Hxc_in(j)=strtod(argv[i],NULL);++i;} //read Hxc1 Hxc2 ... Hxcn
      } // T Hext Hxc
    } // -opmat
    } // -r
    } //maxE         
    } //ninit         
    } //pinit         
    } //nt         
    } // -L  
    } // -S 
   } // -MQ  
   } // -M  
 } // help
  if(argc<2){helpexit();}
  if(nofcomponents==0){fprintf(stdout,"ERROR singleion: please enter exchange field Hxc\n");exit(EXIT_FAILURE);}
  Vector Hxc(1,nofcomponents);Hxc=0;for(j=1;j<=nofcomponents;++j)Hxc(j)=Hxc_in(j);

  int observable_nofcomponents;
  switch(observable)
   {case 'M':
    case 'Q':
    case 'S':
    case 'L': observable_nofcomponents=3;break;
    default: observable_nofcomponents=nofcomponents; // I
   }
  Vector I(1,observable_nofcomponents);
  // transition matrix Mij
  ComplexVector Mq(1,observable_nofcomponents),u1(1,nofcomponents);
 
double TT=T;

if (!do_sipf)
  {par inputpars("./mcphas.j");
   inputpars.save_sipfs("./results/_");
   if(nofcomponents!=inputpars.nofcomponents)fprintf(stderr,"#Warning: number of exchange field components read from command line not equal to that in mcphas.j - continuing...\n");
    printf("#\n#atom-number T[K] ");for(j=1;j<=3;++j)printf("Hext%c(T) ",'a'-1+j);
                                   for(j=1;j<=nofcomponents;++j)printf("Hxc%i(meV) ",j);
                                   switch(observable)
                                   {case 'Q': printf("Q=(%g %g %g)/A ",Q(1),Q(2),Q(3));
                                              for(j=1;j<=observable_nofcomponents;++j)printf(" |<M%c%c>| real(<M%c%c>) imag(<M%c%c>) <M%c>f(Q) ",observable,'a'-1+j,observable,'a'-1+j,observable,'a'-1+j,'a'-1+j);break;
                                    default: for(j=1;j<=observable_nofcomponents;++j)printf(" <%c%c> ",observable,'a'-1+j);
                                   }
                                   printf("transition-energies(meV)...\n");
   if(opmat<1e10){fout_opmat=fopen_errchk("./results/op.mat","w");}
                   
     for(i=1;i<=inputpars.nofatoms;++i)
     {switch(observable)
      {case 'L': (*inputpars.jjj[i]).Lcalc(I,T,Hxc,Hext,(*inputpars.jjj[i]).Icalc_parameter_storage_init(Hxc,Hext,T));break;
       case 'S': (*inputpars.jjj[i]).Scalc(I,T,Hxc,Hext,(*inputpars.jjj[i]).Icalc_parameter_storage_init(Hxc,Hext,T));break;
       case 'M': (*inputpars.jjj[i]).mcalc(I,T,Hxc,Hext,(*inputpars.jjj[i]).Icalc_parameter_storage_init(Hxc,Hext,T));break;
       case 'Q': (*inputpars.jjj[i]).mcalc(I,T,Hxc,Hext,(*inputpars.jjj[i]).Icalc_parameter_storage_init(Hxc,Hext,T));
                 (*inputpars.jjj[i]).eigenstates(Hxc,Hext,T);(*inputpars.jjj[i]).MQ(Mq, Q);
                 break;       
       default: (*inputpars.jjj[i]).Icalc(I,T,Hxc,Hext,lnz,u,(*inputpars.jjj[i]).Icalc_parameter_storage_init(Hxc,Hext,T));
      }  
              
      sprintf(filename,"./results/%s.levels.cef",(*inputpars.jjj[i]).sipffilename);
      fout=fopen_errchk(filename,"w"); 
     fprintf(fout,"#\n#\n#!d=%i sipffile=%s T= %g K ",(*inputpars.jjj[i]).est.Chi(),(*inputpars.jjj[i]).sipffilename,T);
                                   for(j=1;j<=3;++j)fprintf(fout,"Hext%c=%g T ",'a'-1+j,Hext(j));
                                   for(j=1;j<=nofcomponents;++j)fprintf(fout,"Hxc%i=%g meV  ",j,Hxc(j));
                                   switch(observable)
                                   {case 'Q': fprintf(fout,"Q=(%g %g %g)/A ",Q(1),Q(2),Q(3));
                                              for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," M%c%c=%g%+gi ",observable,'a'-1+j,real(Mq(j)),imag(Mq(j)));break;
                                    default: for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," %c%c=%g ",observable,'a'-1+j,I(j));
                                   }
                                   fprintf(fout,"\n");
      printf("%3i %8g ",i,T); // printout ion number and temperature
      for(j=1;j<=3;++j)printf(" %8g ",Hext(j)); // printout external field as requested
      for(j=1;j<=nofcomponents;++j)printf("%8g ",Hxc(j)); // printoutexchangefield as requested
      switch(observable)
       {case 'Q': for(j=1;j<=observable_nofcomponents;++j)printf("%4g %4g %4g %4g   ",abs(Mq(j)),real(Mq(j)),imag(Mq(j)),I(j)*(*inputpars.jjj[i]).F(Norm(Q)));break;
        default: for(j=1;j<=observable_nofcomponents;++j)printf("%4g ",I(j));  // printout corresponding moments      
       } if(nmax>0)
      { sprintf(filename,"./results/%s.trs",(*inputpars.jjj[i]).sipffilename);
        fout_trs = fopen_errchk (filename,"w");
        trs_header_out(fout_trs,pinit,ninit,maxE,TT,Hext,observable);
 
        (*inputpars.jjj[i]).maxE=maxE;(*inputpars.jjj[i]).pinit=pinit;(*inputpars.jjj[i]).ninit=ninit;
        (*inputpars.jjj[i]).transitionnumber=0;int tc=0;nt=0;
        if(trs_write_next_line(fout_trs,(*inputpars.jjj[i]),nt,1,1,1,1,tc,TT,Hxc,Hext,
                                    (*inputpars.jjj[i]).eigenstates(Hxc,Hext,TT),d,-1e100,maxE,observable,Q))
        {fprintf(stderr,"Warning singleion: no transition found within energy in range [minE,maxE]=[%g,%g] found\n"
                        " (within first crystallographic unit of magnetic unit cell)\n"
                        " please increase energy range in option -maxE \n",0.0,maxE);
        }
        else
        {printf("%4g ",d);
         while(tc<nmax&&!trs_write_next_line(fout_trs,(*inputpars.jjj[i]),nt,1,1,1,1,tc,TT,Hxc,Hext,
                          (*inputpars.jjj[i]).est,d,-1e100,maxE,observable,Q)){printf("%4g ",d);if(d>=0)--tc;}
        }
        (*inputpars.jjj[i]).print_eigenstates(fout);fclose(fout);
        fclose(fout_trs);
        if(nmax<nt){printf("...");}
      }
     if(opmat<1e10){fprintf(fout_opmat,"#! d=%i  ",(*inputpars.jjj[i]).est.Chi());
                    if(opmat>nofcomponents){ for(int opmati=0;opmati<=nofcomponents;++opmati)
                                             {Matrix op((*inputpars.jjj[i]).opmat(opmati,Hxc,Hext));
                                             myPrintComplexMatrix(fout_opmat,op);}

                                           }
                    else
                    { if(opmat<-nofcomponents){

                                                Matrix opp((*inputpars.jjj[i]).opmat(0,Hxc,Hext)); opp=0;
                                               for (int opmati=1;opmati<=(*inputpars.jjj[i]).est.Chi();++opmati)
                                               {opp(opmati,opmati)=real((*inputpars.jjj[i]).est(0,opmati));}
                                                myPrintComplexMatrix(fout_opmat,opp);
                                             
                                             for(int opmati=-1;opmati>=-nofcomponents;--opmati)
                                             {Matrix op((*inputpars.jjj[i]).opmat(opmati,Hxc,Hext));
                                             myPrintComplexMatrix(fout_opmat,op);}
                                              }
                     else
                     {Matrix op((*inputpars.jjj[i]).opmat((int)opmat,Hxc,Hext));
                      myPrintComplexMatrix(fout_opmat,op);
                     }
                    }
                   }     
      printf("\n");
     } if(opmat<1e10)fclose(fout_opmat);
   } else { // option -r sipffile
   jjjpar jjj(0,0,0,sipffile,nofcomponents);jjj.save_sipf("./results/_");
      switch(observable)
      {case 'L': jjj.Lcalc(I,T,Hxc,Hext,jjj.Icalc_parameter_storage_init(Hxc,Hext,T));break;
       case 'S': jjj.Scalc(I,T,Hxc,Hext,jjj.Icalc_parameter_storage_init(Hxc,Hext,T));break;
       case 'M': jjj.mcalc(I,T,Hxc,Hext,jjj.Icalc_parameter_storage_init(Hxc,Hext,T));break;
       case 'Q': jjj.mcalc(I,T,Hxc,Hext,jjj.Icalc_parameter_storage_init(Hxc,Hext,T));
                 jjj.eigenstates(Hxc,Hext,T);jjj.MQ(Mq, Q);
                 break;       
      default: jjj.Icalc(I,T,Hxc,Hext,lnz,u,jjj.Icalc_parameter_storage_init(Hxc,Hext,T));
      }          
   printf("#\n#atom-nr T[K] ");for(j=1;j<=3;++j)printf("Hext%c(T) ",'a'-1+j);
                                   for(j=1;j<=nofcomponents;++j)printf("Hxc%i(meV) ",j);
                                   switch(observable)
                                   {case 'Q': printf("Q=(%g %g %g)/A ",Q(1),Q(2),Q(3));
                                              for(j=1;j<=observable_nofcomponents;++j)printf(" |<M%c%c>| real(<M%c%c>) imag(<M%c%c>) <M%c>f(Q) ",observable,'a'-1+j,observable,'a'-1+j,observable,'a'-1+j,'a'-1+j);break;
                                    default: for(j=1;j<=observable_nofcomponents;++j)printf(" <%c%c> ",observable,'a'-1+j);
                                   }
                                   printf("transition-energies(meV)...\n");
      sprintf(filename,"./results/%s.levels.cef",jjj.sipffilename);
      fout=fopen_errchk(filename,"w");  
    fprintf(fout,"#\n#\n#!d=%i sipffile=%s T= %g K ",jjj.est.Chi(),jjj.sipffilename,T);
                                   for(j=1;j<=3;++j)fprintf(fout,"Hext%c=%g T ",'a'-1+j,Hext(j));
                                   for(j=1;j<=nofcomponents;++j)fprintf(fout,"Hxc%i=%g meV  ",j,Hxc(j));
                                   switch(observable)
                                   {case 'Q': fprintf(fout,"Q=(%g %g %g)/A ",Q(1),Q(2),Q(3));
                                              for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," M%c%c=%g%+gi ",observable,'a'-1+j,real(Mq(j)),imag(Mq(j)));break;
                                    default: for(j=1;j<=observable_nofcomponents;++j)fprintf(fout," %c%c=%g ",observable,'a'-1+j,I(j));
                                   }
                                   fprintf(fout,"\n");

   printf("%3i %8g ",1,T); // printout ion number and temperature
      for(j=1;j<=3;++j)printf(" %10g ",Hext(j)); // printout external field as requested
      for(j=1;j<=nofcomponents;++j)printf("%8g ",Hxc(j)); // printoutexchangefield as requested
      switch(observable)
       {case 'Q': for(j=1;j<=observable_nofcomponents;++j)printf("%4g %4g %4g %4g   ",abs(Mq(j)),real(Mq(j)),imag(Mq(j)),I(j)*jjj.F(Norm(Q)));break;
        default: for(j=1;j<=observable_nofcomponents;++j)printf("%4g ",I(j));  // printout corresponding moments      
       }if(nmax>0)
      { sprintf(filename,"./results/%s.trs",jjj.sipffilename);
        fout_trs = fopen_errchk (filename,"w");
        trs_header_out(fout_trs,pinit,ninit,maxE,TT,Hext,observable);
 
        jjj.maxE=maxE;jjj.pinit=pinit;jjj.ninit=ninit;
        jjj.transitionnumber=0;int tc=0;nt=0;
        if(trs_write_next_line(fout_trs,jjj,nt,1,1,1,1,tc,TT,Hxc,Hext,jjj.eigenstates(Hxc,Hext,TT),d,-1e100,maxE,observable,Q))
        {fprintf(stderr,"Warning singleion: no transition found within energy in range [minE,maxE]=[%g,%g]\n"
                        " please increase energy range in option -maxE \n",0.0,maxE);
        }
        else
        {printf("%4g ",d);
         while(tc<nmax&&!trs_write_next_line(fout_trs,jjj,nt,1,1,1,1,tc,TT,Hxc,Hext,
                          jjj.est,d,-1e100,maxE,observable,Q)){if(d>=0)--tc;printf("%4g ",d);}
        }
        jjj.print_eigenstates(fout);fclose(fout);
        fclose(fout_trs);
        if(nmax<nt){printf("...");}
      } if(opmat<1e10){fout_opmat=fopen_errchk("results/op.mat","w");
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
 
      printf("\n");
   }
printf("# **********************************************************************\n"
       "# end of program singleion\n"
       "# ... you can now use 'cpsingleion' to calculate specific heat,\n"
       "#      entropy etc from results/*.levels.cef\n"
       "# **********************************************************************\n");
}
