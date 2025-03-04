
// OBSERVABLE SEQUENCE
// used in mcdisp.c and singleion.c for in/out of .trs files

const char * obs []= { "I", "M", "pel", "P",  "L", "S", "sx", "sy", "sz", "lx", "ly", "lz" , "MQ"  } ;
const char * obunit []= { "", "muB", "|e|pm",   "A",  "hbar", "hbar", "hbar/A^3",  "hbar/A^3",  "hbar/A^3",  "hbar/A^3",  "hbar/A^3",  "hbar/A^3" , "muB" } ;


void trs_header_out(FILE* fout,double & pinit,double & ninit,double & maxE,char * out,ob observable)
{time_t curtime;char cc[30];cc[0]=obs[observable][0];cc[1]='\0';
 switch(observable)
{case sx:
 case sy:
 case sz:cc[0]='a';cc[1]='S';cc[2]=observable;cc[3]='\0';break;
 case lx:cc[0]='a';cc[1]='L';cc[2]='x';cc[3]='\0';break;
 case ly:cc[0]='a';cc[1]='L';cc[2]='y';cc[3]='\0';break;
 case lz:cc[0]='a';cc[1]='L';cc[2]='z';cc[3]='\0';break;
 case MQ: cc[0]='M';cc[1]='Q';cc[2]='\0';break;
 default: break;
}
 struct tm *loctime;
   fprintf(fout, "#output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcdisp.trs-->\n");
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(fout,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(fout,"#*********************************************************************\n");
   if(observable!=I&&observable!=MQ)
   {fprintf(fout,"#(*)The unpolarized powder average RAMAN cross section sigma for each transition \n");
   fprintf(fout,"#   is calculated neglecting  the Debye Wallerfactor as follows:\n");
   fprintf(fout,"#-------------------------------------------------------------- \n");
   fprintf(fout,"#           Contributions to X(Omega)                          |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                  =                                           |\n");
   fprintf(fout,"#  X%s(omega)= sum |   PI *  delta(E-(Ei-Ek))                 |\n",obs[obint(observable)]);
   fprintf(fout,"#                  =                                           |\n");
   fprintf(fout,"#                  E -> E                                      |\n");
   fprintf(fout,"#                   i    k                                     |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"# =                                2                           |\n");
   fprintf(fout,"# |     =           wi  |<i|%s|k>|                            |\n",obs[obint(observable)]);
   fprintf(fout,"# =                                                            |\n");
   fprintf(fout,"#  E -> E                                                      |\n");
   fprintf(fout,"#   i    k                                                     |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                            with                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                      - E /T                                  |\n");
   fprintf(fout,"#                     e   i                                    |\n");
   fprintf(fout,"# wi    =       --------------                                 |\n");
   fprintf(fout,"#               ----    - E /T                                 |\n");
   fprintf(fout,"#               >    n e   i                                   |\n");
   fprintf(fout,"#               ----  i                                        |\n");
   fprintf(fout,"#                i                                             |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                               -----                          |\n");
   fprintf(fout,"#                  2       1    \\                        2    |\n");
   fprintf(fout,"#   |<i,r|%s|k,s>|   =    -     >   |<i,r|%s -<%s>|k,s>|    |\n",obs[obint(observable)],obs[obint(observable)],obs[obint(observable)]);
   fprintf(fout,"#                          3    /           u      u           |\n");
   fprintf(fout,"#                              -----                           |\n");
   fprintf(fout,"#                            u = x,y,z                         |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#         2                  2                                 |\n");
   fprintf(fout,"#      %s   unit is (%s)                                   |\n",obs[obint(observable)],obunit[obint(observable)]);
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   }
   else {
   fprintf(fout,"#(*)The unpolarized powder average neutron cross section sigma for each transition \n");
   fprintf(fout,"#   is calculated neglecting the formfactor, the Debye Wallerfactor, factor k'/k as follows:\n");
   fprintf(fout,"#-------------------------------------------------------------- \n");
   fprintf(fout,"#            Transition intensities in barn/sr.                |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                      =                                       |\n");
   fprintf(fout,"#          S(E)= N sum |    delta(E-(Ei-Ek))                   |\n");
   fprintf(fout,"#                      =                                       |\n");
   fprintf(fout,"#                       E -> E                                 |\n");
   fprintf(fout,"#                        i    k                                |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#   with N ... number of ions un the beam                      |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"# =                           2                                |\n");
   fprintf(fout,"# |     = const wi  |<i|M |k>|                                 |\n");
   fprintf(fout,"# =                      T                                     |\n");
   fprintf(fout,"#  E -> E                                                      |\n");
   fprintf(fout,"#   i    k                                                     |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                            with                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                      - E /T                                  |\n");
   fprintf(fout,"#                     e   i                                    |\n");
   fprintf(fout,"# wi    =       --------------                                 |\n");
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
   }
   fprintf(fout,"#! ninit= %g (max number of initial states) -do not modify: needed to count transitions\n",ninit);
   fprintf(fout,"#! pinit= %g (minimum population number of initial states)-do not modify: needed to count transitions\n",pinit);
   fprintf(fout,"#! maxE= %g meV(maximum value of transition energy)-do not modify: needed to count transitions\n",maxE);
   fprintf(fout,"#! %s\n",out);
   
   fprintf(fout,"#1 2 3  4 *** 5 ****** 6 *********** 7 ******* 8 ********************"
                   " 9  10  ********* 11 ******************** 12 *********************\n");
   fprintf(fout,"#i j k ionnr transnr energy(meV) |gamma_s| ");
   if(observable!=I&&observable!=MQ)fprintf(fout,"Tr(X%s)/3[(%s)^2](*)  ",obs[obint(observable)],obunit[obint(observable)]);   
   else               fprintf(fout,"sigma_mag_dip[barn/sr](*) ");
                 fprintf(fout, " n  n'  wnn'|<n|%s1-<%s1>|n'>|^2 wnn'|<n|%s2-<%s2>|n'>|^2 ... "
                   "with wnn'=wn-wn' for n!=n'  and wnn=wn/k_B T \n",cc,cc,cc,cc);
}

//****************************************************************************************
// probes transitions and returns 0 if transition is found, transitionnumber is stored in jjj.transitionnumber
//                        returns 1 if no further transition is found within limit minE maxE
int trs_write_next_line(FILE * fout,jjjpar & jjj,int & nt,int  i,int  j,int  k,int  l,int & tc,double & T,Vector & mf,
                     Vector & Hext,ComplexMatrix & est,float & d,double  minE,double  maxE, ob observable, Vector & Q)    
    {ComplexVector u1(1,mf.Hi());double gamma;int n=0,nd=0;
         if(jjj.transitionnumber>=nt&&nt>0){return 1;}
     ++jjj.transitionnumber;nt=jjj.du1calc(T,mf,Hext,u1,d,n,nd,est);
    while (minE>=d||d>=maxE) //only consider transition if it is in interval minE/maxE
     {//first and following  transitions out of energy range ... do not consider them
     //fprintf(stdout," .... transition not stored because out of interval [minE,maxE]=[%g,%g]meV\n",minE,maxE);
     ++jjj.transitionnumber;
     //fprintf(stdout,"nt=%i transition number %i: ",nt,jjj.transitionnumber);
     if(jjj.transitionnumber>nt){return 1;}
     jjj.du1calc(T,mf,Hext,u1,d,n,nd,est);
     }
   //fprintf(stdout,"nt=%i transition number %i: ",nt,jjj.transitionnumber);
    gamma=Norm2(u1);ComplexVector dm1(1,3);double intensityp=0, intensitym=0; dm1=0;
ComplexVector m1(1,SPINDENS_EV_DIM); m1=0;int ch=0;
    switch(observable)
     {case S: ch=jjj.dS1calc(T,mf,Hext,dm1,est);break;
      case L: ch=jjj.dL1calc(T,mf,Hext,dm1,est);break;
      case MQ: ch=jjj.dMQ1calc(Q,T,dm1,d,est);break;
      case P: ch=jjj.dP1calc(T,mf,Hext,dm1,est);break;
      case pel: ch=jjj.dpel1calc(T,mf,Hext,dm1,est);break;
      case sx: jjj.dspindensity_coeff1(1,T,mf,Hext,m1,est);break;
      case sy: jjj.dspindensity_coeff1(2,T,mf,Hext,m1,est);break;
      case sz: jjj.dspindensity_coeff1(3,T,mf,Hext,m1,est);break;
      case lx: jjj.dorbmomdensity_coeff1(1,T,mf,Hext,m1,est);break;
      case ly: jjj.dorbmomdensity_coeff1(2,T,mf,Hext,m1,est);break;
      case lz: jjj.dorbmomdensity_coeff1(3,T,mf,Hext,m1,est);break;
      default: break;
     }    
  
switch(observable)
     {case sx:
      case sy:
      case sz:
      case lx:
      case ly:
      case lz:
      case I:
      case M: 
      case MQ: 
             // calculate powder neutron intensities 
     if(jjj.dm1calc(T,mf,Hext,dm1,est)) // if dm1calc is implemented for this ion
     {intensityp+=Norm2(dm1); // Norm2 ... sum of modulus squared
     intensityp*=0.048434541067;intensitym=intensityp;// prefactor for intensity in barn/sr is 2/3*0.53908*0.53908/4= 0.048434541067
     if (d>SMALL_QUASIELASTIC_ENERGY){if(d/T/KB<20){intensitym=-intensityp/(1-exp(d/T/KB));intensityp/=(1-exp(-d/T/KB));}else{intensitym=0;}}
                                  else{intensityp=intensityp*T*KB;intensitym=intensityp;}
     }
     else
     {intensityp=-1;intensitym=-1;dm1=0;} 
             break; 
       default: 
     // calculate Trace X
    if(ch) // if dObs1calc is implemented for this ion
     {intensityp+=Norm2(dm1); // Norm2 ... sum of modulus squared
     intensityp*=1.0/3;// no prefactor just Trace/3
     intensitym=intensityp;
     if (d>SMALL_QUASIELASTIC_ENERGY){if(d/T/KB<20){intensitym=-intensityp/(1-exp(d/T/KB));intensityp/=(1-exp(-d/T/KB));}else{intensitym=0;}}
                                  else{intensityp=intensityp*T*KB;intensitym=intensityp;}
     }
     else
     {intensityp=-1;intensitym=-1;dm1=0;} 
   break;
     }
      
   
     if(minE<d&&d<maxE)
    { fprintf(fout,"%i %i %i  %i     %i     %9.6g  %9.6g  %10.6g  %i %i ",i,j,k,l,jjj.transitionnumber,myround(d),myround(gamma),myround(intensityp),n,nd);
       switch(observable)
     {case I:for(int i=1;i<=u1.Hi();++i)fprintf(fout," %9.6g",real(conj(u1(i))*u1(i)));break;
      case sx:
      case sy:
      case sz:
      case lx:
      case ly:
      case lz: for(int i=1;i<=m1.Hi();++i)fprintf(fout," %9.6g",real(conj(m1(i))*m1(i)));
     default: for(int i=1;i<=dm1.Hi();++i)fprintf(fout," %9.6g",real(conj(dm1(i))*dm1(i)));
      }fprintf(fout,"\n");
     ++tc;}
    if(d>=0&&minE<-d&&-d<maxE) // do not print negative energy transition if d<0 (d<0 means transiton to the same level)
    { fprintf(fout,"%i %i %i  %i     %i     %9.6g  %9.6g  %10.6g  %i %i ",i,j,k,l,jjj.transitionnumber,myround(-d),myround(gamma),myround(intensitym),nd,n);
       switch(observable)
     {case I:for(int i=1;i<=u1.Hi();++i)fprintf(fout," %9.6g",real(conj(u1(i))*u1(i)));break;
      case sx:
      case sy:
      case sz:
      case lx:
      case ly:
      case lz: for(int i=1;i<=m1.Hi();++i)fprintf(fout," %9.6g",real(conj(m1(i))*m1(i)));
      default: for(int i=1;i<=dm1.Hi();++i)fprintf(fout," %9.6g",real(conj(dm1(i))*dm1(i)));
      } fprintf(fout,"\n");
    ++tc;}
 return 0;
} //write next line  


