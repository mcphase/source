// used in mcdisp.c and singleion.c for in/out of .trs files

void trs_header_out(FILE* fout,double & pinit,double & ninit,double & maxE,double & T,Vector & Hext,char observable)
{time_t curtime;
 struct tm *loctime;
   fprintf(fout, "#output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcdisp.trs-->\n");
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(fout,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"#(*)The unpolarized powder average neutron cross section sigma for each transition \n");
   fprintf(fout,"#   is calculated neglecting the formfactor, the Debye Wallerfactor, factor k'/k as follows:\n");
   fprintf(fout,"#-------------------------------------------------------------- \n");
   fprintf(fout,"#            Transition intensities in barn/sr.                |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"# =                                                            |\n");
   fprintf(fout,"# |                           2                                |\n");
   fprintf(fout,"# |     = const wi  |<i|M |k>|                                 |\n");
   fprintf(fout,"# |                      T                                     |\n");
   fprintf(fout,"# =                                                            |\n");
   fprintf(fout,"#  E -> E                                                      |\n");
   fprintf(fout,"#   i    k                                                     |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                            with                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                      - E /T                                  |\n");
   fprintf(fout,"#                     e   i                                    |\n");
   fprintf(fout,"# wi    = const --------------                                 |\n");
   fprintf(fout,"#               ----    - E /T                                 |\n");
   fprintf(fout,"#               >    n e   i                                   |\n");
   fprintf(fout,"#               ----  i                                        |\n");
   fprintf(fout,"#                i                                             |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                -----                         |\n");
   fprintf(fout,"#                      2     2   \\                        2   |\n");
   fprintf(fout,"#        |<i,r|M |k,s>|   = ---   >     |<i,r|M -<M >|k,s>|    |\n");
   fprintf(fout,"#               T            3   /             u   u           |\n");
   fprintf(fout,"#                                -----                         |\n");
   fprintf(fout,"#                             u = x,y,z                        |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                             and                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                   1       2                  |\n");
   fprintf(fout,"#                  const  =      ( --- r   )                   |\n");
   fprintf(fout,"#                                   2   0                      |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                      -12                     |\n");
   fprintf(fout,"#                  r     = -0.53908* 10    cm                  |\n");
   fprintf(fout,"#                   0                                          |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                 M   =  L  + 2 S  = g  J                      |\n");
   fprintf(fout,"#                                     J                        |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#--------------------------------------------------------------|\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                       1.Sum rule :                           |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#  ----  =       2                                             |\n");
   fprintf(fout,"#  >     |     =--- *g *g *const * J(J+1) *wi                  |\n");
   fprintf(fout,"#  ----  =       3    J  J                                     |\n");
   fprintf(fout,"#   k     E -> E                                               |\n");
   fprintf(fout,"#          i    k                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#--------------------------------------------------------------|\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                       2. sum rule :                          |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#            ----  =            2                              |\n");
   fprintf(fout,"#            >     |         = --- * const*g *g *J(J+1)        |\n");
   fprintf(fout,"#            ----  =            3           J  J               |\n");
   fprintf(fout,"#             k,i   E -> E                                     |\n");
   fprintf(fout,"#                    i    k                                    |\n");
   fprintf(fout,"#-------------------------------------------------------------- \n");
   fprintf(fout,"#! ninit= %g (max number of initial states) -do not modify: needed to count transitions\n",ninit);
   fprintf(fout,"#! pinit= %g (minimum population number of initial states)-do not modify: needed to count transitions\n",pinit);
   fprintf(fout,"#! maxE= %g meV(maximum value of transition energy)-do not modify: needed to count transitions\n",maxE);
   fprintf(fout,"#! T= %g K Ha=%g Hb=%g Hc=%g T\n",T,Hext(1),Hext(2),Hext(3));
   fprintf(fout,"#*********************************************************************\n");
   fprintf (fout, "#i j k ionnr transnr energy |gamma_s|  sigma [barn/sr](*) "
                   "   wnn'|<n|%c1-<%c1>|n'>|^2 wnn'|<n|%c2-<%c2>|n'>|^2 ... "
                   "with wnn'=wn-wn' for n!=n'  and wnn=wn/k_B T \n",observable,observable,observable,observable);
}

//****************************************************************************************
// probes transitions and returns 0 if transition is found, transitionnumber is stored in jjj.transitionnumber
//                        returns 1 if no further transition is found within limit minE maxE
int trs_write_next_line(FILE * fout,jjjpar & jjj,int & nt,int  i,int  j,int  k,int  l,int & tc,double & T,Vector & mf,
                     Vector & Hext,ComplexMatrix & est,float & d,double  minE,double  maxE, char observable)    
    {ComplexVector u1(1,mf.Hi());double gamma;
     ++jjj.transitionnumber;nt=jjj.du1calc(T,mf,Hext,u1,d,est);
    while (minE>d||d>maxE) //only consider transition if it is in interval minE/maxE
     {//first and following  transitions out of energy range ... do not consider them
     //fprintf(stdout," .... transition not stored because out of interval [minE,maxE]=[%g,%g]meV\n",minE,maxE);
     ++jjj.transitionnumber;
     fprintf(stdout,"nt=%i transition number %i: ",nt,jjj.transitionnumber);
     jjj.du1calc(T,mf,Hext,u1,d,est);
     if(jjj.transitionnumber>nt){return 1;}
     }
    if(jjj.transitionnumber>nt){return 1;}
    gamma=Norm2(u1);
     // calculate powder neutron intensities 
     double intensityp=0, intensitym=0; ComplexVector dm1(1,3);
     if(jjj.dm1calc(T,mf,Hext,dm1,est)) // if dm1calc is implemented for this ion
     {intensityp+=Norm2(dm1); // Norm2 ... sum of modulus squared
     intensityp*=0.048434541067;intensitym=intensityp;// prefactor for intensity in barn/sr is 2/3*0.53908*0.53908/4= 0.048434541067
     if (d>SMALL_QUASIELASTIC_ENERGY){if(d/T/KB<20){intensitym=-intensityp/(1-exp(d/T/KB));intensityp/=(1-exp(-d/T/KB));}else{intensitym=0;}}
                                  else{intensityp=intensityp*T*KB;intensitym=intensityp;}
     }
     else
     {intensityp=-1;intensitym=-1;}   
    switch(observable)
     {case 'M': break; // leave it, it was calulated above 
      case 'S': jjj.dS1calc(T,mf,Hext,dm1,est);break;
      case 'L': jjj.dL1calc(T,mf,Hext,dm1,est);break;
     }    
     if(minE<d&&d<maxE)
    { fprintf(fout,"%i %i %i  %i     %i     %9.6g  %9.6g  %10.6g   ",i,j,k,l,jjj.transitionnumber,myround(d),myround(gamma),myround(intensityp));
       switch(observable)
     {case 'I':for(int i=1;i<=u1.Hi();++i)fprintf(fout," %9.6g",real(conj(u1(i))*u1(i)));break;
      default: for(int i=1;i<=dm1.Hi();++i)fprintf(fout," %9.6g",real(conj(dm1(i))*dm1(i)));
      }fprintf(fout,"\n");
     ++tc;}
    if(d>=0&&minE<-d&&-d<maxE) // do not print negative energy transition if d<0 (d<0 means transiton to the same level)
    { fprintf(fout,"%i %i %i  %i     %i     %9.6g  %9.6g  %10.6g   ",i,j,k,l,jjj.transitionnumber,myround(-d),myround(gamma),myround(intensitym));
       switch(observable)
     {case 'I':for(int i=1;i<=u1.Hi();++i)fprintf(fout," %9.6g",real(conj(u1(i))*u1(i)));break;
      default: for(int i=1;i<=dm1.Hi();++i)fprintf(fout," %9.6g",real(conj(dm1(i))*dm1(i)));
      } fprintf(fout,"\n");
    ++tc;}
 return 0;
} //write next line  


