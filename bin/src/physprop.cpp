// methods for class parameters 
#include "physprop.hpp"
#include "../../version"
 // *************************************************************************
 // ************************ physproperties *********************************
 // *************************************************************************


void   sort(float * v,int jmin,int jmax,int * jnew) // sorting function
{// input: v[jmin ... jmax] vector to be sorted
 // jnew [jmin .... jmax]   sort index - after sorting the 
 //                         vector jnew contains the indices
 //                         in a sequence such that v[jnew[jmin ... jmax]] is
 //                         sorted in ascending order 
 int gap,i,j,temp,n;
 n=jmax-jmin+1;
 // initialize jnew
 for(i=jmin;i<=jmax;++i){jnew[i]=i;}
 // if nothing is to be done - return
 if (jmax<=jmin) {return;}
 
 for (gap=n/2;gap>0;gap/=2)
   for (i=gap;i<n;++i)
     for (j=i-gap;j>=0 && v[jnew[jmin+j]]>v[jnew[jmin+j+gap]];j-=gap)
      {temp=jnew[jmin+j];jnew[jmin+j]=jnew[jmin+j+gap]; jnew[jmin+j+gap]=temp;}
      
 return;
}
 
//constructor
physproperties::physproperties (int nofspincorrs,int maxnofhkli,int na,int nm)
{washere=0;
 int i;
 nofspincorr=nofspincorrs;
 nofatoms=na;
 nofcomponents=nm;
 
 m=Vector(1,3);
 H=Vector(1,3);
 jj= new Vector [nofspincorrs+1];for(i=0;i<=nofspincorrs;++i){jj[i]=Vector(1,nofcomponents*nofcomponents*nofatoms);} //  ... number of interaction constants (aa bb cc ab ba ac ca bc cb)
   if (jj == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 hkli= new Vector [maxnofhkli+1];for(i=0;i<=maxnofhkli;++i){hkli[i]=Vector(1,10);}
   if (hkli == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 nofhkls=0;
 sps=spincf(1,1,1,nofatoms,nofcomponents);
 mf=mfcf(1,1,1,nofatoms,nofcomponents);
 fe=0;u=0;
}


//kopier-konstruktor
physproperties::physproperties (const physproperties & p)
{int i;
  x=p.x;y=p.y;
  j=p.j;
  T=p.T;
  H=p.H;
  m=p.m;
  nofhkls=p.nofhkls;
  u=p.u;fe=p.fe;Eel=p.Eel;
  sps=p.sps;
  mf=p.mf;
  washere=p.washere;
 maxnofhkls=p.maxnofhkls;
 nofspincorr=p.nofspincorr;
 nofatoms=p.nofatoms;
 nofcomponents=p.nofcomponents;
  
 jj= new Vector [nofspincorr+1];for(i=0;i<=nofspincorr;++i){jj[i]=Vector(1,nofcomponents*nofcomponents*nofatoms);} //  ... number of interaction constants (aa bb cc ab ba ac ca bc cb)
   if (jj == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 hkli= new Vector [maxnofhkls+1];for(i=0;i<=maxnofhkls;++i){hkli[i]=Vector(1,10);}
   if (hkli == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 for(i=1;i<=nofspincorr;++i)
    {jj[i]=p.jj[i];}
 for(i=1;i<=nofhkls;++i)
    {hkli[i]=p.hkli[i];} 
 
 }


//destruktor
physproperties::~physproperties ()
{//printf("hello destruktor physprop\n");  
delete []jj;delete []hkli;
//printf("hello destruktor physprop\n");  
 
}

void physproperties::update_maxnofhkls(int maxnofhkli)
{delete []hkli;
 maxnofhkls=maxnofhkli;
 int i;
 hkli= new Vector [maxnofhkli+1];for(i=0;i<=maxnofhkli;++i){hkli[i]=Vector(1,10);}
   if (hkli == NULL){fprintf (stderr, "physproperties::update_maxnofhkls - Out of memory\n");exit (EXIT_FAILURE);} 
}


// methode save
double physproperties::save (int verbose, const char * filemode, int htfailed,inipar & ini, par & inputpars,char * prefix)
{ FILE *fout;
  char filename[50];
  time_t curtime;
  struct tm *loctime;  
  int i,j2,l,i1,j1,nmax;
  Vector null(1,nofcomponents*nofatoms);null=0;
  Vector null1(1,3);null1=0;
  double sta;
  sta=0; 
  float nn[200];nn[0]=199;
  float nnerr[200];nnerr[0]=199;
  int ortho=1;
  if (inputpars.cs.alpha()!=90||inputpars.cs.beta()!=90||inputpars.cs.gamma()!=90){ortho=0;}
      Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1; // Trick: put a=b=c=1 to have unit vectors in a b and c directions
                       abc(4)=inputpars.cs.alpha(); abc(5)=inputpars.cs.beta(); abc(6)=inputpars.cs.gamma();
        // ... then use dadadc2ijk and ijk2dadbdc with this unit length "lattice"
   Vector mabc(1,3),Hijk(1,3);
   dadbdc2ijk(Hijk,H,abc);
   ijk2dadbdc(mabc,m,abc);

  printf("saving properties for T=%g K  Ha= %g Hb= %g Hc= %g T\n", T,H(1),H(2),H(3));

//-----------------------------------------------------------------------------------------  
  errno = 0;char outfilename[MAXNOFCHARINLINE];
  strcpy(outfilename,"./results/");strcpy(outfilename+10,prefix);
  strcpy(outfilename+10+strlen(prefix),"mcphas.fum");
  if (verbose==1) printf("saving %s \n",outfilename);
  if (washere==0)
  {fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#{output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.fum-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   fprintf (fout, "#note: - for specific heat calculation use unit conversion 1mev/f.u.=96.48J/mol\n");
   fprintf (fout, "#      - below moments and energies are given per ion - not per formula unit !\n");
   if(ini.doeps){
    fprintf (fout, "#      - strain tensor eps is calculated selfconsistently.\n");
    if(ini.linepscf)
    {
      fprintf (fout, "#      ... however, singleion ion Hamiltonian is always diagonalised with eps=0 (option -linepscf).\n");
    }
    if(ini.linepsjj)
    {
      fprintf (fout, "#      ... however, two ion interaction is always evaluated for eps=0 (option -linepsjj).\n");
    }
   }

   if(ortho==0){fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");}
   fprintf (fout, "#   x    y   T[K] H[T] Ha[T] Hb[T] Hc[T] free energy f[meV/ion] energy u[meV/ion] total moment |m|     ma mb mc m||(projection along H) [mb/ion]");
   if(ortho==0){fprintf (fout, "mi   mj   mk[muB/ion]   Hi  Hj  Hk[T]");}
   if(ini.doeps&&ortho==0){fprintf (fout, " Eel [meV/ion] eps1=epsii eps2=epsjj eps3=epskk epse4=2epsjk eps5=2epsik eps6=2epsij");}
   if(ini.doeps&&ortho!=0){fprintf (fout, " Eel [meV/ion] eps1=epsaa eps2=epsbb eps3=epscc epse4=2epsbc eps5=2epsac eps6=2epsab");}
   fprintf (fout,"}\n");
   fclose(fout);
      }
   if (htfailed!=0){fe=0;u=0;m=0;m[1]=0;m[2]=0;m[3]=0;Eel=0;}
   fout = fopen_errchk (outfilename,"a");
   fprintf (fout, "%4.4g %4.4g  %4.4g %4.4g %4.4g %4.4g %4.4g       %8.8g            %8.8g       %4.4g    %4.4g %4.4g %4.4g    %4.4g",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),myround(fe),myround(u),myround(Norm(m)),myround(mabc[1]),myround(mabc[2]),myround(mabc[3]),myround(m*Hijk/Norm(Hijk)));
   if(ortho==0){fprintf (fout, "    %4.4g %4.4g %4.4g   %4.4g %4.4g %4.4g",myround(m(1)),myround(m(2)),myround(m(3)),Hijk(1),Hijk(2),Hijk(3));}
    if(ini.doeps){fprintf (fout, "  %4.4g  %4.4g %4.4g %4.4g   %4.4g %4.4g %4.4g",myround(Eel),myround(sps.epsilon(1)),myround(sps.epsilon(2)),myround(sps.epsilon(3)),myround(sps.epsilon(4)),myround(sps.epsilon(5)),myround(sps.epsilon(6)));}
   fprintf(fout,"\n");
   fclose(fout);
   strcpy(outfilename,"./results/.");strcpy(outfilename+11,prefix);
  strcpy(outfilename+11+strlen(prefix),"mcphas.fum");
   fout = fopen_errchk (outfilename,"a");
   fprintf (fout, "%4.4g %4.4g  %4.4g %4.4g %4.4g %4.4g %4.4g %8.8g %8.8g  %4.4g %4.4g %4.4g %4.4g %4.4g\n",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),myround(fe),myround(u),myround(Norm(m)),myround(mabc[1]),
            myround(mabc[2]),myround(mabc[3]),myround(m*Hijk/Norm(Hijk)));
   fclose(fout);
   if((fout=fopen("./fit/mcphas.fum","rb"))!=NULL)
    {// some measured data should be fitted
     j2=0;
     if(washere==0){fprintf(stdout,"#Mcphas- calculating standard deviation to ./fit/mcphas.fum - magnetisation ma mb mc (col 11,12,13)\n");}
     while(feof(fout)==0)
     {++j2;
      if ((l=inputline(fout,nn,nnerr))!=0)
      {if(0.0001>(nn[3]-T)*(nn[3]-T)+(nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
         { for(i=8;i<=l;++i){double clc=0;if(ortho==0){
                         switch(i) { case 8: clc=fe;break;
                                     case 9: clc=u;break;
                                     case 10: clc=Norm(m);break;
                                     case 11: clc=mabc(1);break;
                                     case 12: clc=mabc(2);break;
                                     case 13: clc=mabc(3);break;
                                     case 14: clc=m*Hijk/Norm(Hijk);break;
                                     case 15: clc=m[1];break;
                                     case 16: clc=m[2];break;
                                     case 17: clc=m[3];break;                                    

                                     case 21: clc=Eel;break;
                                     case 22: sps.epsilon(1);break;
                                     case 23: sps.epsilon(2);break;
                                     case 24: sps.epsilon(3);break;
                                     case 25: sps.epsilon(4);break;
                                     case 26: sps.epsilon(5);break;
                                     case 27: sps.epsilon(6);break;
                                     default: clc=0;
                                   } } else
{switch(i) { case 8: clc=fe;break;
                                     case 9: clc=u;break;
                                     case 10: clc=Norm(m);break;
                                     case 11: clc=mabc(1);break;
                                     case 12: clc=mabc(2);break;
                                     case 13: clc=mabc(3);break;
                                     case 14: clc=m*Hijk/Norm(Hijk);break;
                                     case 15: clc=Eel;break;
                                     case 16: sps.epsilon(1);break;
                                     case 17: sps.epsilon(2);break;
                                     case 18: sps.epsilon(3);break;
                                     case 19: sps.epsilon(4);break;
                                     case 20: sps.epsilon(5);break;
                                     case 21: sps.epsilon(6);break;
                                     default: clc=0;
                                   }
}

                         if(nnerr[i]>0){if(verbose==1&&clc==0){fprintf(stdout,"sta_mcphas.fum warning: exp value cannot be fitted in line %i column %i in file ./fit/mcphas.fum\n",j2,i);}
//fprintf(stdout,"stacalc_mphas.fum: col %i line %i value %g err %g - calcvalue %g\n",i,j2,nn[i],nnerr[i],clc);
                                       sta+=(nn[i]-clc)*(nn[i]-clc)/nnerr[i]/nnerr[i];}
                        }
            if(verbose==1){fprintf(stdout,"sta_mcphas.fum=%g\n",sta);}
	 }
      }
     }
     fclose(fout);
    }else{errno=0;}
//-----------------------------------------------------------------------------------------  
  errno = 0; Vector totalJ(1,nofcomponents);
  strcpy(outfilename,"./results/");strcpy(outfilename+10,prefix);
  strcpy(outfilename+10+strlen(prefix),"mcphas.xyt");
    if (verbose==1)printf("saving %s\n",outfilename);
  if (washere==0)
  {fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#{output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.xyt-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   fprintf (fout, "#x    y   T[K] H[T] Ha[T] Hb[T] Hc[T] phasnumber-j   period-key supercell-nr1 x nr2 x nr3 ");
           for(i1=1;i1<=nofcomponents;++i1)
	      {//fprintf(fout,"<I%c> ",'a'-1+i1);}
	      fprintf(fout,"<I%i> ",i1);}
	      fprintf(fout,"\n");
   fclose(fout);    
     }
  fout = fopen_errchk (outfilename,"a");
  totalJ=0; 
  if (htfailed!=0){j=0;}else{totalJ=sps.totalJ();}
   if(j<0){sps.wasstable=j;}// if qvector generated structure is stable, then take period key = number of qvector
   fprintf (fout, "%4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g       %ip           %ip                %i x %i x %i ",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),j,sps.wasstable,sps.na(),sps.nb(),sps.nc());
           for(i1=1;i1<=nofcomponents;++i1)
	      {fprintf(fout,"%4.4g ",myround(totalJ(i1)));}
	      fprintf(fout,"\n");
   fclose(fout);
    if((fout=fopen("./fit/mcphas.xyt","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.xyt not implemented\n");}
     fclose(fout);
    }else{errno=0;}

//-----------------------------------------------------------------------------------------  
 nmax=nofspincorr; // look how many spincorrelationfunction we have indeed calculated - the 
                   // user wanted nofspincorr, but maybe it was fewer ...
 for (l=1;l<=inputpars.cs.nofatoms;++l)
      {if(nmax>(*inputpars.jjj[l]).paranz)
       {nmax=(*inputpars.jjj[l]).paranz;
        fprintf(stderr,"Warning: calculation of nofspincorr=%i correlation functions not possible, \n",nofspincorr);
fprintf(stderr,"         because in mcphas.j for atom %i  only %i neighbours are given.\n",l,(*inputpars.jjj[l]).paranz);
       }
      }
       

  // only output nmax correlation functions ...
 for(i=1;i<=nmax;++i){for(l=1;l<=nofatoms;++l){
  errno = 0;
  if (verbose==1)printf("saving mcphas%i.j%i - spinspin corr for sublattice %i neighbour %i\n",l,i,l,i);
  strcpy(outfilename,"./results/");strcpy(outfilename+10,prefix);
  strcpy(outfilename+10+strlen(prefix),"mcphas");
  snprintf(filename,sizeof(filename),"%s%i.j%i",outfilename,l,i);
  if (washere==0)  //printout file header
  {  fout = fopen_errchk (filename,filemode);
   fprintf(fout, "#{output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas*.j*-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   fprintf (fout, "# sublattice %i (da=%g a db=%g b dc=%g c)\n",l,(*inputpars.jjj[l]).xyz(1),(*inputpars.jjj[l]).xyz(2),(*inputpars.jjj[l]).xyz(3));
   fprintf (fout, "# correlation function <JJ(%g %g %g)>\n",myround((*inputpars.jjj[l]).dn[i](1)),myround((*inputpars.jjj[l]).dn[i](2)),myround((*inputpars.jjj[l]).dn[i](3)));
   fprintf (fout, "#x     y     T[K]  H[T]   Ha[T] Hb[T] Hc[T]  ");
           for(i1=1;i1<=(*inputpars.jjj[l]).nofcomponents;++i1)
	      {fprintf(fout,"<J%cJ%c> ",'a'-1+i1,'a'-1+i1);}
           for(i1=1;i1<=(*inputpars.jjj[l]).nofcomponents-1;++i1)
              {for(j1=i1+1;j1<=(*inputpars.jjj[l]).nofcomponents;++j1)
                        {fprintf(fout,"<J%cJ%c> <J%cJ%c> ",'a'-1+i1,'a'-1+j1,'a'-1+j1,'a'-1+i1);
			}
	      }
                          
   fprintf (fout,"}\n");
   fclose(fout);
      }
  fout = fopen_errchk (filename,"a");
  if (htfailed!=0){jj[i](1)=0;jj[i](2)=0;jj[i](3)=0;}
   fprintf (fout, "%4.4g %4.4g   %4.4g %4.4g   %4.4g %4.4g %4.4g     ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
        for(j2=1;j2<=nofcomponents*nofcomponents;++j2)               
            {fprintf (fout, "%4.4g ",myround(jj[i](j2+nofcomponents*nofcomponents*(l-1))));
	    }
   fprintf (fout,"\n");
   fclose(fout);
   snprintf(filename,sizeof(filename),"./fit/mcphas%i.j%i",l,i);
   if((fout=fopen(filename,"rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using %s  not implemented\n",filename);}
     fclose(fout);
    }else{errno=0;}
  }}

//-----------------------------------------------------------------------------------------  
 errno = 0;
  strcpy(outfilename,"./results/");strcpy(outfilename+10,prefix);
  strcpy(outfilename+10+strlen(prefix),"mcphas*.hkl");
  if (verbose==1)printf("saving %s - neutrons and xrays\n",outfilename);
  if (washere==0)
  {//neutrons
   strcpy(outfilename+10+strlen(prefix),"mcphas.hkl");
  fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#{output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.hkl-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   fprintf (fout, "#Neutron Intensity - Mind: only structure+polarizationfactor+formfactor+debeywallerfactor - no lorentzfactor is  taken into account\n");
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  int       h   k   l   int       h   k   l   int ...}\n");
   fclose(fout);
   //xray a component
   if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphasi.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasa.hkl");  }
   fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#{output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphasa.hkl-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
      if(ortho==0){
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - i component\n");
   fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  real(mi(Q)) im(mi(Q))    h   k   l   real(mi(Q)) im(mi(Q))     h   k   l   real(mi(Q)) im(mi(Q))[mu_B/atom] ...}\n");
   }else{
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - a component\n"); 
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  real(ma(Q)) im(ma(Q))      h   k   l    real(ma(Q)) im(ma(Q))      h   k   l  real(ma(Q)) im(ma(Q)) [mu_B/atom]...}\n");
   }
   fclose(fout);
   //xray b component
   if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphasj.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasb.hkl");  }
   fout = fopen_errchk (outfilename,filemode);
    fprintf(fout, "#{output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphasb.hkl-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
      if(ortho==0){
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - j component\n");
   fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  real(mj(Q)) im(mj(Q))       h   k   l   real(mj(Q)) im(mj(Q))       h   k   l   real(mj(Q)) im(mj(Q))[mu_B/atom]...}\n");
   }else{
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - b component\n"); 
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  real(mb(Q)) im(mb(Q))        h   k   l  real(mb(Q)) im(mb(Q))     h   k   l  real(mb(Q)) im(mb(Q)) [mu_B/atom] ...}\n");
   }
   fclose(fout);
   //xray c component
   if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphask.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasc.hkl");  }
   fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#{output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphasc.hkl-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
      if(ortho==0){
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - k component\n");
   fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  real(mk(Q)) im(mk(Q))       h   k   l  real(mk(Q)) im(mk(Q))      h   k   l   real(mk(Q)) im(mk(Q)) [mu_B/atom]...}\n");
   }else{
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - c component\n"); 
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  real(mc(Q)) im(mc(Q))       h   k   l  real(mc(Q)) im(mc(Q))       h   k   l  real(mc(Q)) im(mc(Q))  [mu_B/atom]...}\n");
   }
   fclose(fout);

      }
   int * inew;inew=new int[nofhkls+1];float *intensity;intensity=new float[nofhkls+1];
   if(inew==NULL){fprintf (stderr, "Out of memory for inew\n");exit (EXIT_FAILURE);}
   if(intensity==NULL){fprintf (stderr, "Out of memory for intensity\n");exit (EXIT_FAILURE);}
   //printf("nofhkls=%i\n",nofhkls);
   for (i=1;i<=nofhkls;++i) {intensity[i]=hkli[i](4);}
  if (verbose==1)printf(" .... sorting hkl according to neutron intensities\n");
   sort(intensity,1,nofhkls,inew); // sort according to ascending intensity
  //neutrons
   strcpy(outfilename+10+strlen(prefix),"mcphas.hkl"); 
   fout = fopen_errchk (outfilename,"a");
   if (verbose==1)printf(" .... saving %s\n",outfilename);
   fprintf (fout, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
   for (i=nofhkls;i>=1;--i)
    {if (htfailed!=0){hkli[inew[i]](1)=0;hkli[inew[i]](2)=0;hkli[inew[i]](3)=0;hkli[inew[i]](4)=0;}
    fprintf (fout, "%4.4g %4.4g %4.4g  %4.4g     ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](4)));
    } fprintf(fout,"\n");
   fclose(fout);
  //xray a component
  if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphasi.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasa.hkl");  }
   fout = fopen_errchk (outfilename,"a");
   if (verbose==1)printf(" .... saving %s\n",outfilename);
   fprintf (fout, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
   for (i=nofhkls;i>=1;--i)
    {fprintf (fout, "%4.4g %4.4g %4.4g  %4.4g %4.4g    ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](5)),myround(hkli[inew[i]](6)));
    } fprintf(fout,"\n");
   fclose(fout);
  //xray b component
   if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphasj.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasb.hkl");  }
   fout = fopen_errchk (outfilename,"a");
   if (verbose==1)printf(" .... saving %s\n",outfilename);
   fprintf (fout, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
   for (i=nofhkls;i>=1;--i)
    {fprintf (fout, "%4.4g %4.4g %4.4g  %4.4g %4.4g    ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](7)),myround(hkli[inew[i]](8)));
    } fprintf(fout,"\n");
   fclose(fout);
  //xray c component
   if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphask.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasc.hkl");  }
   fout = fopen_errchk (outfilename,"a");
   if (verbose==1)printf(" .... saving %s\n",outfilename);
   fprintf (fout, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
   for (i=nofhkls;i>=1;--i)
    {fprintf (fout, "%4.4g %4.4g %4.4g  %4.4g %4.4g    ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](9)),myround(hkli[inew[i]](10)));
    } fprintf(fout,"\n");
   fclose(fout);


    if((fout=fopen("./fit/mcphas.hkl","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.hkl not implemented\n");}
     fclose(fout);
    }else{errno=0;}
   delete []inew;delete []intensity;
//-----------------------------------------------------------------------------------------  
 errno = 0;
 strcpy(outfilename+10+strlen(prefix),"mcphas.sps");
   if (verbose==1)printf("saving %s- spinconfiguration\n",outfilename);
  if (washere==0)
  {  fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#{output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.sps-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   // printout the lattice and atomic positions
   strcpy(inputpars.rems[2],"#\n");
   inputpars.savelattice(fout);inputpars.saveatoms(fout);
   fprintf (fout, "#!show_abc_unitcell=1.0\n");
   fprintf (fout, "#!show_primitive_crystal_unitcell=1.0\n");
   fprintf (fout, "#!show_magnetic_unitcell=1.0\n");
   fprintf (fout, "#!show_atoms=1.0\n");
   fprintf (fout, "#!spins_scale_moment=1.0\n");
   fprintf (fout, "#!scale_view_1=1.0 scale_view_2=1.0 scale_view_3=1.0\n");
   if(ortho==0){fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");}
   
   fprintf (fout, "#x y T[K] |H| H[T] Ha[T] Hb[T] Hc[T] nofspins nofatoms(in primitive basis) nofmomentum-components errorcode(0=ok,1=failed) eps1=epsii eps2=epsjj eps3=epskk eps4=2epsjk eps5=2epsik eps6=2epsij\n");
   fprintf (fout, "    #<Ja(1)> <Ja(2)> .... selfconsistent Spinconfiguration  \n");
   fprintf (fout, "    #<Jb(1)> <Jb(2)> .... UNITS:  multiply by gJ to get moment [muB]\n");
   fprintf (fout, "    #<Jc(1)> <Jc(2)> ....}\n");
    fclose(fout);
   }  
  fout = fopen_errchk (outfilename,"a");
   fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i ",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),sps.n()*sps.nofatoms,sps.nofatoms,sps.nofcomponents);
   if (htfailed!=0){fprintf(fout,"1 ");sps.spinfromq(1,1,1,null1,null,null,null);} // failed
    else {fprintf(fout,"0 ");}
   fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g\n",myround(sps.epsilon(1)),myround(sps.epsilon(2)),myround(sps.epsilon(3)),myround(sps.epsilon(4)),myround(sps.epsilon(5)),myround(sps.epsilon(6)));
    sps.print(fout);fprintf(fout,"\n");
   fclose(fout);
    if((fout=fopen("./fit/mcphas.sps","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.sps not implemented\n");}
     fclose(fout);
    }else{errno=0;}
 
//-----------------------------------------------------------------------------------------  
 errno = 0;
  strcpy(outfilename+10+strlen(prefix),"mcphas.mf");
  if (verbose==1)printf("saving %s - mean field configuration\n",outfilename);
  if (washere==0)
  {  fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#{output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.mf-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   // printout the lattice and atomic positions
   strcpy(inputpars.rems[2],"#\n");
   inputpars.savelattice(fout);inputpars.saveatoms(fout);
   fprintf (fout, "#!show_abc_unitcell=1.0\n");
   fprintf (fout, "#!show_primitive_crystal_unitcell=1.0\n");
   fprintf (fout, "#!show_magnetic_unitcell=1.0\n");
   fprintf (fout, "#!show_atoms=1.0\n");
   fprintf (fout, "#!show_chargedensity=1.0\n");
   fprintf (fout, "#!spins_scale_moment=1.0\n");
   fprintf (fout, "#!scale_view_1=1.0 scale_view_2=1.0 scale_view_3=1.0\n");
   if(ortho==0){fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");}
   fprintf (fout, "#x y T[K] |H| H[T] Ha[T] Hb[T] Hc[T] nofspins nofatoms(in primitive basis) nofmeanfield-components errorcode(0=ok,1=failed) eps1=epsii eps2=epsjj eps3=epskk eps4=2epsjk eps5=2epsik eps6=2epsij\n");
   fprintf (fout, "    #mfa(1) mfa(2) .... selfconsistent Mean field configuration \n"); 
   fprintf (fout, "    #mfb(1) mfb(2) .... UNITS: mf(i)=gJ*mu_B*hxc(i)[meV] \n"); 
   fprintf (fout, "    #mfc(1) mfc(2) ....         (i.e. divide by gJ and mu_B=0.05788meV/T to get exchange field[T]}\n");
    fclose(fout);
   }  
     fout = fopen_errchk (outfilename,"a");
fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i ",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),
           myround(H[2]),myround(H[3]),mf.n()*mf.nofatoms,mf.nofatoms,mf.nofcomponents);
   if (htfailed!=0){fprintf(fout,"1 %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g\n",myround(sps.epsilon(1)),myround(sps.epsilon(2)),myround(sps.epsilon(3)),myround(sps.epsilon(4)),myround(sps.epsilon(5)),myround(sps.epsilon(6)));
                    sps.print(fout);fprintf(fout,"\n");}
   else
    {fprintf(fout,"0 %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g\n",myround(sps.epsilon(1)),myround(sps.epsilon(2)),myround(sps.epsilon(3)),myround(sps.epsilon(4)),myround(sps.epsilon(5)),myround(sps.epsilon(6)));
     mf.print(fout);fprintf(fout,"\n");}
   fclose(fout);
    if((fout=fopen("./fit/mcphas.mf","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.mf not implemented\n");}
     fclose(fout);
    }else{errno=0;}
//-----------------------------------------------------------------------------------------  
 
 washere=1;
return sta;
 }

// scroll output files and read physical properties from these if possible,
// on success return 0, otherwise
// return 1
int physproperties::read(int verbose, par & inputpars,char * readprefix)
{ FILE *fin;
  char filename[50];
  int i,j2,l,nmax;
  float nn[200];nn[0]=199;
  int ortho=1,found=0;
  if (inputpars.cs.alpha()!=90||inputpars.cs.beta()!=90||inputpars.cs.gamma()!=90){ortho=0;}
      Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1;
                       abc(4)=inputpars.cs.alpha(); abc(5)=inputpars.cs.beta(); abc(6)=inputpars.cs.gamma();

   Vector mabc(1,3),Hijk(1,3);
   dadbdc2ijk(Hijk,H,abc);
   ijk2dadbdc(mabc,m,abc);

  printf("reading properties for T=%g K  Ha= %g Hb= %g Hc= %g T\n", T,H(1),H(2),H(3));

//-----------------------------------------------------------------------------------------  
  errno = 0;char infilename[MAXNOFCHARINLINE];
  strcpy(infilename,"./results/");strcpy(infilename+10,readprefix);
  strcpy(infilename+10+strlen(readprefix),"mcphas.fum");
  if (verbose==1) printf("reading %s \n",infilename);
   fin = fopen_errchk (infilename,"r");
// here read free energy etc if possible ... otherwise return 1
// check x y T H[1] H[2] H[3] agrees and fe nonzero ?
// then read fe, u, 
// read mabc[1-3] convert to ijk

found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
    //x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
   if(0.0001>(nn[1]-x)*(nn[1]-x)+(nn[2]-y)*(nn[2]-y)+
             (nn[3]-T)*(nn[3]-T)+
             (nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
             found=1;
               } 
   fclose(fin);
   // check if stable, i.e. if free energy nn[8] is nonzero
   if(nn[8]==0) return 1; // 
   fe=nn[8];u=nn[9];mabc[1]=nn[11];mabc[2]=nn[12];mabc[3]=nn[13];
   dadbdc2ijk(m,mabc,abc);
//fprintf (fout, "%4.4g %4.4g  %4.4g %4.4g %4.4g %4.4g %4.4g       %8.8g            %8.8g       %4.4g    %4.4g %4.4g %4.4g    %4.4g",
 //           myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),
//myround(H[2]),myround(H[3]),myround(fe),myround(u),myround(Norm(m)),
//myround(mabc[1]),myround(mabc[2]),myround(mabc[3]),myround(m*Hijk/Norm(Hijk)));
 //  if(ortho==0){fprintf (fout, "    %4.4g %4.4g %4.4g   %4.4g %4.4g %4.4g",myround(m(1)),myround(m(2)),myround(m(3)),Hijk(1),Hijk(2),Hijk(3));}
  
//-----------------------------------------------------------------------------------------  
  errno = 0; Vector totalJ(1,nofcomponents);
  strcpy(infilename,"./results/");strcpy(infilename+10,readprefix);
  strcpy(infilename+10+strlen(readprefix),"mcphas.xyt");
    if (verbose==1)printf("reading %s\n",infilename);
  fin = fopen_errchk (infilename,"r");
// check x y T H[1] H[2] H[3] agrees and fe nonzero ?
//then read j ... not necessary
   //fprintf (fin, "%4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g       %ip           %ip      ",
     //       myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),j,sps.wasstable);
        // then read totalJ ... not necessary if sps is read !!
//   for(i1=1;i1<=nofcomponents;++i1)
//	      {fprintf(fin,"%4.4g ",myround(totalJ(i1)));}
//	      fprintf(fin,"\n");
found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
    //x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
   if(0.0001>(nn[1]-x)*(nn[1]-x)+(nn[2]-y)*(nn[2]-y)+
             (nn[3]-T)*(nn[3]-T)+
             (nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
             found=1;
               } 
   fclose(fin);
   
//-----------------------------------------------------------------------------------------  
 nmax=nofspincorr; // look how many spincorrelationfunction we have indeed calculated - the 
                   // user wanted nofspincorr, but maybe it was fewer ...
 for (l=1;l<=inputpars.cs.nofatoms;++l)
      {if(nmax>(*inputpars.jjj[l]).paranz)
       {nmax=(*inputpars.jjj[l]).paranz;
        fprintf(stderr,"Warning: reading of nofspincorr=%i correlation functions not possible, \n",nofspincorr);
fprintf(stderr,"         because in mcphas.j for atom %i  only %i neighbours are given.\n",l,(*inputpars.jjj[l]).paranz);
       }
      }
       

  // only read nmax correlation functions ...
 for(i=1;i<=nmax;++i){for(l=1;l<=nofatoms;++l){
  errno = 0;
  if (verbose==1)printf("reading mcphas%i.j%i - spinspin corr for sublattice %i neighbour %i\n",l,i,l,i);
  strcpy(infilename,"./results/");strcpy(infilename+10,readprefix);
  strcpy(infilename+10+strlen(readprefix),"mcphas");
  snprintf(filename,sizeof(filename),"%s%i.j%i",infilename,l,i);
  fin = fopen_errchk (filename,"r");
   // check x y T H[1] H[2] H[3] agrees and fe nonzero ?
//then read jj[] 

 // fprintf (fin, "%4.4g %4.4g   %4.4g %4.4g   %4.4g %4.4g %4.4g     ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
  //      for(j2=1;j2<=nofcomponents*nofcomponents;++j2)               
   //         {fprintf (fin, "%4.4g ",myround(jj[i](j2+nofcomponents*nofcomponents*(l-1))));
//	    }
 //  fprintf (fin,"\n");
  found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
    //x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
   if(0.0001>(nn[1]-x)*(nn[1]-x)+(nn[2]-y)*(nn[2]-y)+
             (nn[3]-T)*(nn[3]-T)+
             (nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
             found=1;
               } 
 fclose(fin);
 for(j2=1;j2<=nofcomponents*nofcomponents;++j2)jj[i](j2+nofcomponents*nofcomponents*(l-1))=nn[j2+7];

  }}

//-----------------------------------------------------------------------------------------  
 errno = 0;
  strcpy(infilename,"./results/");strcpy(infilename+10,readprefix);
  strcpy(infilename+10+strlen(readprefix),"mcphas*.hkl");
  if (verbose==1)printf("reading %s - neutrons and xrays\n",infilename);

  //neutrons
   strcpy(infilename+10+strlen(readprefix),"mcphas.hkl"); 
   fin = fopen_errchk (infilename,"r");
   if (verbose==1)printf(" .... reading %s\n",infilename);
   // check x y T H[1] H[2] H[3] agrees and
  //then read hkli[1-4]

    //fprintf (fin, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
//   for (i=nofhkls;i>=1;--i)
//    {fprintf (fin, "%4.4g %4.4g %4.4g  %4.4g     ",myround(hkli[i](1)),myround(hkli[i](2)),
//myround(hkli[i](3)),myround(hkli[i](4)));
//    } fprintf(fin,"\n");
found=0;while(found==0){ i=0;
 while (feof(fin)==0&&0==i)i=inputline(fin,nn);  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
    //x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
   if(0.0001>(nn[1]-x)*(nn[1]-x)+(nn[2]-y)*(nn[2]-y)+
             (nn[3]-T)*(nn[3]-T)+
             (nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
             found=1;nofhkls=(i-7)/4;
               } 
   fclose(fin);
   for (i=nofhkls;i>=1;--i){hkli[i](1)=nn[7+(i-1)*4+1];
                            hkli[i](2)=nn[7+(i-1)*4+2];
                            hkli[i](3)=nn[7+(i-1)*4+3];
                            hkli[i](4)=nn[7+(i-1)*4+4];}
 
  //xray a component
  if(ortho==0){strcpy(infilename+10+strlen(readprefix),"mcphasi.hkl");
                } else {strcpy(infilename+10+strlen(readprefix),"mcphasa.hkl");  }
   fin = fopen_errchk (infilename,"r");
   if (verbose==1)printf(" .... reading %s\n",infilename);
   // check x y T H[1] H[2] H[3] agrees and fe nonzero ?
//then read hkli[5,6]
//    for (i=nofhkls;i>=1;--i)
//    {fprintf (fin, "%4.4g %4.4g %4.4g  %4.4g %4.4g    ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](5)),myround(hkli[inew[i]](6)));
//    } fprintf(fin,"\n");
found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
    //x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
   if(0.0001>(nn[1]-x)*(nn[1]-x)+(nn[2]-y)*(nn[2]-y)+
             (nn[3]-T)*(nn[3]-T)+
             (nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
             found=1;
               } 
   fclose(fin);
  for (i=nofhkls;i>=1;--i){hkli[i](5)=nn[7+(i-1)*4+4];
                           hkli[i](6)=nn[7+(i-1)*4+5];
                            }
//xray b component
   if(ortho==0){strcpy(infilename+10+strlen(readprefix),"mcphasj.hkl");
                } else {strcpy(infilename+10+strlen(readprefix),"mcphasb.hkl");  }
   fin = fopen_errchk (infilename,"r");
   if (verbose==1)printf(" .... reading %s\n",infilename);
    // check x y T H[1] H[2] H[3] agrees and fe nonzero ?
//then read hkli[7,8]
//fprintf (fin, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
//   for (i=nofhkls;i>=1;--i)
//    {fprintf (fin, "%4.4g %4.4g %4.4g  %4.4g %4.4g    ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](7)),myround(hkli[inew[i]](8)));
//    } fprintf(fin,"\n");
found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
    //x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
   if(0.0001>(nn[1]-x)*(nn[1]-x)+(nn[2]-y)*(nn[2]-y)+
             (nn[3]-T)*(nn[3]-T)+
             (nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
             found=1;
               } 
   fclose(fin);
  for (i=nofhkls;i>=1;--i){hkli[i](7)=nn[7+(i-1)*4+4];
                           hkli[i](8)=nn[7+(i-1)*4+5];
                            }
  //xray c component
   if(ortho==0){strcpy(infilename+10+strlen(readprefix),"mcphask.hkl");
                } else {strcpy(infilename+10+strlen(readprefix),"mcphasc.hkl");  }
   fin = fopen_errchk (infilename,"r");
   if (verbose==1)printf(" .... reading %s\n",infilename);
   // check x y T H[1] H[2] H[3] agrees and fe nonzero ?
//then read hkli[9,10]
//fprintf (fin, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
 //  for (i=nofhkls;i>=1;--i)
  //  {fprintf (fin, "%4.4g %4.4g %4.4g  %4.4g %4.4g    ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](9)),myround(hkli[inew[i]](10)));
   // } fprintf(fin,"\n");
found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
    //x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
   if(0.0001>(nn[1]-x)*(nn[1]-x)+(nn[2]-y)*(nn[2]-y)+
             (nn[3]-T)*(nn[3]-T)+
             (nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
             found=1;
               } 
   fclose(fin);
  for (i=nofhkls;i>=1;--i){hkli[i](9)=nn[7+(i-1)*4+4];
                           hkli[i](10)=nn[7+(i-1)*4+5];
                            }

//-----------------------------------------------------------------------------------------  
 errno = 0;
 strcpy(infilename+10+strlen(readprefix),"mcphas.sps");
   if (verbose==1)printf("reasding %s- spinconfiguration\n",infilename);
  fin = fopen_errchk (infilename,"r");
   // check x y T H[1] H[2] H[3] agrees and fe nonzero ?
//then read sps  ... probably sps.read ...??

  // fprintf (fin, " %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i ",
   //         myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),sps.n()*sps.nofatoms,sps.nofatoms,sps.nofcomponents);
//    fprintf(fin,"0 = ok\n");sps.print(fin);fprintf(fin,"\n");
found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
    //x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
   sps.load(fin);
   if(0.0001>(nn[1]-x)*(nn[1]-x)+(nn[2]-y)*(nn[2]-y)+
             (nn[3]-T)*(nn[3]-T)+
             (nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
             found=1;
               } 
   fclose(fin);
  
 
//-----------------------------------------------------------------------------------------  
 errno = 0;
  strcpy(infilename+10+strlen(readprefix),"mcphas.mf");
  if (verbose==1)printf("reading %s - mean field configuration\n",infilename);
     fin = fopen_errchk (infilename,"r");
// check x y T H[1] H[2] H[3] agrees and fe nonzero ?
//then read mf

//fprintf (fin, " %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i ",
 //           myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),
  //         myround(H[2]),myround(H[3]),mf.n()*mf.nofatoms,mf.nofatoms,mf.nofcomponents);
//   fprintf(fin,"0 = ok\n");mf.print(fin);fprintf(fin,"\n");
found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
    //x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
   mf.load(fin);
   if(0.0001>(nn[1]-x)*(nn[1]-x)+(nn[2]-y)*(nn[2]-y)+
             (nn[3]-T)*(nn[3]-T)+
             (nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
             found=1;
               } 
   fclose(fin);
   //-----------------------------------------------------------------------------------------  

return 0; 
 }
