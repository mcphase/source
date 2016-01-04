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
  u=p.u;fe=p.fe;
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
double physproperties::save (int verbose, const char * filemode, int htfailed, par & inputpars,char * prefix)
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
  int ortho=1;
  if (inputpars.alpha!=90||inputpars.beta!=90||inputpars.gamma!=90){ortho=0;}
      Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1;
                       abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;

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
   fprintf (fout, "#note: - for specific heat calculation use unit conversion 1mev/ion=96.48J/mol\n");
   fprintf (fout, "#      - moments and energies are given per ion - not per formula unit !\n");
   if(ortho==0){fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");}
   fprintf (fout, "#   x    y   T[K] H[T] Ha[T] Hb[T] Hc[T] free energy f[meV/ion] energy u[meV/ion] total moment |m|     ma mb mc m||(projection along H) [mb/ion]");
   if(ortho==0){fprintf (fout, "mi   mj   mk[muB/f.u.]   Hi  Hj  Hk[T]");}
   fprintf (fout,"}\n");
   fclose(fout);
      }
   if (htfailed!=0){fe=0;u=0;m=0;m[1]=0;m[2]=0;m[3]=0;}
   fout = fopen_errchk (outfilename,"a");
   fprintf (fout, "%4.4g %4.4g  %4.4g %4.4g %4.4g %4.4g %4.4g       %8.8g            %8.8g       %4.4g    %4.4g %4.4g %4.4g    %4.4g",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),myround(fe),myround(u),myround(Norm(m)),myround(mabc[1]),myround(mabc[2]),myround(mabc[3]),myround(m*Hijk/Norm(Hijk)));
   if(ortho==0){fprintf (fout, "    %4.4g %4.4g %4.4g   %4.4g %4.4g %4.4g",myround(m(1)),myround(m(2)),myround(m(3)),Hijk(1),Hijk(2),Hijk(3));}
   fprintf(fout,"\n");
   fclose(fout);
   fout = fopen_errchk (outfilename,"a");
   fprintf (fout, "%4.4g %4.4g  %4.4g %4.4g %4.4g %4.4g %4.4g %8.8g %8.8g  %4.4g %4.4g %4.4g %4.4g %4.4g\n",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),myround(fe),myround(u),myround(Norm(m)),myround(mabc[1]),
            myround(mabc[2]),myround(mabc[3]),myround(m*Hijk/Norm(Hijk)));
   fclose(fout);
   if((fout=fopen("./fit/mcphas.fum","rb"))!=NULL)
    {// some measured data should be fitted
     if(washere==0){fprintf(stdout,"#Mcphas- calculating standard deviation to ./fit/mcphas.fum - magnetisation ma mb mc (col 11,12,13)\n");}
     while(feof(fout)==0)
     {
      if (inputline(fout,nn)!=0)
      {if(0.0001>(nn[3]-T)*(nn[3]-T)+(nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
         {sta+=(nn[11]-m[1])*(nn[11]-m[1])+(nn[12]-m[2])*(nn[12]-m[2])+(nn[13]-m[3])*(nn[13]-m[3]);
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
   fprintf (fout, "#x    y   T[K] H[T] Ha[T] Hb[T] Hc[T] phasnumber-j   period-key ");
           for(i1=1;i1<=nofcomponents;++i1)
	      {fprintf(fout,"<I%c> ",'a'-1+i1);}
	      fprintf(fout,"\n");
   fclose(fout);    
     }
  fout = fopen_errchk (outfilename,"a");
  totalJ=0; 
  if (htfailed!=0){j=0;}else{totalJ=sps.totalJ();}
   if(j<0){sps.wasstable=j;}
   fprintf (fout, "%4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g       %ip           %ip      ",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),j,sps.wasstable);
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
 for (l=1;l<=inputpars.nofatoms;++l)
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
  sprintf(filename,"%s%i.j%i",outfilename,l,i);
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
   fprintf (fout, "# correlation fuction <JJ(%g %g %g)>\n",myround((*inputpars.jjj[l]).dn[i](1)),myround((*inputpars.jjj[l]).dn[i](2)),myround((*inputpars.jjj[l]).dn[i](3)));
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
   sprintf(filename,"./fit/mcphas%i.j%i",l,i);
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
   fprintf (fout, "#x y T[K] |H| H[T] Ha[T] Hb[T] Hc[T] nofspins nofatoms(in primitive basis) nofmomentum-components\n");
   fprintf (fout, "    #<Ja(1)> <Ja(2)> .... selfconsistent Spinconfiguration  \n");
   fprintf (fout, "    #<Jb(1)> <Jb(2)> .... UNITS:  multiply by gJ to get moment [muB]\n");
   fprintf (fout, "    #<Jc(1)> <Jc(2)> ....}\n");
    fclose(fout);
   }  
  fout = fopen_errchk (outfilename,"a");
   fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i ",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),sps.n()*sps.nofatoms,sps.nofatoms,sps.nofcomponents);
   if (htfailed!=0){fprintf(fout,"1 = failed");sps.spinfromq(1,1,1,null1,null,null,null);}
    fprintf(fout,"0 = ok\n");sps.print(fout);fprintf(fout,"\n");
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
   fprintf (fout, "#x y T[K] |H| H[T] Ha[T] Hb[T] Hc[T] nofspins nofatoms(in primitive basis) nofmeanfield-components\n");
   fprintf (fout, "    #mfa(1) mfa(2) .... selfconsistent Mean field configuration \n"); 
   fprintf (fout, "    #mfb(1) mfb(2) .... UNITS: mf(i)=gJ*mu_B*hxc(i)[meV] \n"); 
   fprintf (fout, "    #mfc(1) mfc(2) ....         (i.e. divide by gJ and mu_B=0.05788meV/T to get exchange field[T]}\n");
    fclose(fout);
   }  
     fout = fopen_errchk (outfilename,"a");
fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i ",
            myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),
           myround(H[2]),myround(H[3]),mf.n()*mf.nofatoms,mf.nofatoms,mf.nofcomponents);
   if (htfailed!=0){fprintf(fout,"1 = failed\n");sps.print(fout);fprintf(fout,"\n");}else
    {fprintf(fout,"0 = ok\n");mf.print(fout);fprintf(fout,"\n");}
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
