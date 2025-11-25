// methods for class parameters 
#include "physprop.hpp"
#include "../../version"
 // *************************************************************************
 // ************************ physproperties *********************************
 // *************************************************************************


void   sort(float * v,int jmin,int & jmax,int * jnew) // sorting function
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
physproperties::physproperties (int & nofspincorrs,int & maxnofhkli,cryststruct & csin)
{washere=0;
 int i;
 nofspincorr=nofspincorrs;
 cs=cryststruct(csin);

 m=Vector(1,3); 
 mabc=Vector(1,3); 
 Pelabc=Vector(1,3); Pel=Vector(1,3);
 Pelabc0=Vector(1,3); Pel0=Vector(1,3);
 H=Vector(1,HEXT_DIMENSION);
 totalJ=Vector(1,cs.nofcomponents);
 cel=Matrix(1,6,1,6);

 jj= new Vector [nofspincorrs+1];for(i=0;i<=nofspincorrs;++i){jj[i]=Vector(1,cs.nofcomponents*cs.nofcomponents*cs.nofatoms);} //  ... number of interaction constants (aa bb cc ab ba ac ca bc cb)
   if (jj == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 hkli= new Vector [maxnofhkli+1];for(i=0;i<=maxnofhkli;++i){hkli[i]=Vector(1,10);}
   if (hkli == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 nofhkls=0;
 sps=spincf(1,1,1,cs.nofatoms,cs.nofcomponents);
 mf=mfcf(1,1,1,cs.nofatoms,cs.nofcomponents);
 fe=0;u=0;
}


//kopier-konstruktor
physproperties::physproperties (const physproperties & p)
{int i;
  x=p.x;y=p.y;
  j=p.j;
  T=p.T;
  m=Vector(1,3); 
  mabc=Vector(1,3); 
  Pelabc=Vector(1,3); Pel=Vector(1,3);
  Pelabc0=Vector(1,3); Pel0=Vector(1,3);
  H=Vector(1,HEXT_DIMENSION);
  cel=Matrix(1,6,1,6);cel=p.cel;
  cv=p.cv;
  H=p.H; 
  Pel=p.Pel;Pelabc=p.Pelabc;
  Pel0=p.Pel0;Pelabc0=p.Pelabc0;
  m=p.m;mabc=p.mabc;
  nofhkls=p.nofhkls;
  u=p.u;fe=p.fe;Eelastic=p.Eelastic;
  washere=p.washere;
 maxnofhkls=p.maxnofhkls;
 nofspincorr=p.nofspincorr;
 cs=cryststruct(p.cs);

  sps=spincf(1,1,1,cs.nofatoms,cs.nofcomponents);
  mf=mfcf(1,1,1,cs.nofatoms,cs.nofcomponents);
  sps=p.sps;
  mf=p.mf;
 totalJ=Vector(1,cs.nofcomponents);
 totalJ=p.totalJ;

 jj= new Vector [nofspincorr+1];for(i=0;i<=nofspincorr;++i){jj[i]=Vector(1,cs.nofcomponents*cs.nofcomponents*cs.nofatoms);} //  ... number of interaction constants (aa bb cc ab ba ac ca bc cb)
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

void physproperties::update_maxnofhkls(int & maxnofhkli)
{delete []hkli;
 maxnofhkls=maxnofhkli;
 int i;
 hkli= new Vector [maxnofhkli+1];for(i=0;i<=maxnofhkli;++i){hkli[i]=Vector(1,10);}
   if (hkli == NULL){fprintf (stderr, "physproperties::update_maxnofhkls - Out of memory\n");exit (EXIT_FAILURE);} 
}





  
   // 1.  puts header for fum file columns >8 into string header
   // 2. sets or reads output column field nn -
   //  if setnn true: for all i>8 up to input nofcols ... if nnerr[i]!=0 -> increase sta according to difference 
   //                 (and finally return sta) ... then ...
   //                 set  nn[i] from saved parameters fe,u,etc. and puts into nofcols the number of output columns
   //                 puts into outstr the numbers  nn[i>8] formatted for output into mcphas.fum
   //  if setnn false: reads nn[8-nofcols]  into parameters fe, u, etc  
   // for fum file
double physproperties::fumcols(float * nn,float * nnerr, int & nofcols,bool setnn,char * header,char * outstr,inipar & ini, int & ortho,par & inputpars,int & verbose)
 {double sta=0;double * ptr;char hs[40];char num[40];
   header[0]='\0';outstr[0]='\0';
    int nofcolsin=0;if(!setnn){nofcolsin=nofcols; Pel=0;Pelabc=0;m=0;mabc=0;fe=0;u=0;sps.epsilon=0;}
     double Nm=Norm(m),mp=m*H(1,3)/Norm(H(1,3));
    nofcols=15;for(int i=8;i<=nofcolsin||i<=nofcols;++i)
    {ptr=NULL;
          if(i<15){switch(i) {       case 8: ptr=&fe;
                                              if(ini.nofrndtries>=0)snprintf(hs,40,"free_energy_f[meV/ion]");
                                              else                  snprintf(hs,40,"energy_u[meV/ion]"); // for Monte Carlo
                                              break;
                                     case 9: ptr=&u;snprintf(hs,40,"energy_u[meV/ion]");break;
                                     case 10: ptr=&Nm;snprintf(hs,40,"total_moment|m|[mb/ion]");break;
                                     case 11: ptr=&mabc[1];snprintf(hs,40,"ma[mb/ion]");break;
                                     case 12: ptr=&mabc[2];snprintf(hs,40,"mb[mb/ion]");break;
                                     case 13: ptr=&mabc[3];snprintf(hs,40,"mc[mb/ion]");break;
                                     case 14: ptr=&mp;snprintf(hs,40,"m||(projection_along_H)[mb/ion]");break;
                                     default: ;
                              }
                   } 
           else   { // i>=15
                   nofcols=14;             
                   if(ortho==0){nofcols=17;
                         switch(i) { 
                                     case 15: ptr=&m[1];snprintf(hs,40,"mi[mb/ion]");break;
                                     case 16: ptr=&m[2];snprintf(hs,40,"mj[mb/ion]");break;
                                     case 17: ptr=&m[3];snprintf(hs,40,"mk[mb/ion]");break;                                    
                                     default: ;
                                   }
                               }
                  if(ini.cv!=0){if(i-nofcols==1)ptr=&cv;snprintf(hs,40,"cV[meV/ionK]");
                                nofcols+=1;
                               }
                  if(ini.doeps){switch(i-nofcols) {
                                     case 1: ptr=&Eelastic;snprintf(hs,40,"Eelastic[meV/ion]");break;
                                     case 2: ptr=&sps.epsilon[1];snprintf(hs,40,"eps1=epsii");break;
                                     case 3: ptr=&sps.epsilon[2];snprintf(hs,40,"eps2=epsjj");break;
                                     case 4: ptr=&sps.epsilon[3];snprintf(hs,40,"eps3=epskk");break;
                                     case 5: ptr=&sps.epsilon[4];snprintf(hs,40,"eps4=2epsjk");break;
                                     case 6: ptr=&sps.epsilon[5];snprintf(hs,40,"eps5=2epsik");break;
                                     case 7: ptr=&sps.epsilon[6];snprintf(hs,40,"eps6=2epsij");break;
                                     default: ;
                                           }
                               nofcols+=7;
                               }
                  if(fabs(inputpars.totalcharge)<SMALLCHARGE)
                            {switch(i-nofcols) { 
                                     case 1: ptr=&Pelabc[1];snprintf(hs,40,"Pela[C/m^2]");break;
                                     case 2: ptr=&Pelabc[2];snprintf(hs,40,"Pelb[C/m^2]");break;
                                     case 3: ptr=&Pelabc[3];snprintf(hs,40,"Pelc[C/m^2]");break;                                    
                                     default: ;
                                          }
                       if(ortho==0){nofcols+=3;
                         switch(i-nofcols) { 
                                     case 1: ptr=&Pel[1];snprintf(hs,40,"Peli[C/m^2]");break;
                                     case 2: ptr=&Pel[2];snprintf(hs,40,"Pelj[C/m^2]");break;
                                     case 3: ptr=&Pel[3];snprintf(hs,40,"Pelk[C/m^2]");break;                                    
                                     default: ;
                                   }
                                 } 
                             nofcols+=3;
                            }
                  
           
                 }
    if(ptr==NULL)
    {if(nnerr[i]>0&&i<=nofcolsin&&verbose==1)
     fprintf(stdout,"sta_mcphas.fum warning: exp value %g cannot be fitted in column %i in file ./fit/mcphas.fum\n",nn[i],i);
    }
    else
    {
    if(setnn)
        {
      if(nnerr[i]>0&&i<=nofcolsin)
            {
//fprintf(stdout,"stacalc_mphas.fum: col %i line %i value %g err %g - calcvalue %g\n",i,j2,nn[i],nnerr[i],clc);
             sta+=(nn[i]-(*ptr))*(nn[i]-(*ptr))/nnerr[i]/nnerr[i];
             }
          nn[i]=(*ptr);
          if(i<=9)snprintf(num,40,"%8.8g ",nn[i]);else snprintf(num,40,"%4.4g ",nn[i]);
          snprintf(outstr+strlen(outstr),MAXNOFCHARINLINE-strlen(outstr),"%*s%s",(int)(strlen(hs)+1-strlen(num)),"",num); 
        }else
        {  (*ptr)=nn[i];
        }
     snprintf(header+strlen(header),MAXNOFCHARINLINE-strlen(header),"%s ",hs);
    } 
    } // next i
if(!setnn)
{// here set mabc m Pelabc P if necessary
 Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1;
                      abc(4)=inputpars.cs.alpha(); abc(5)=inputpars.cs.beta(); abc(6)=inputpars.cs.gamma();
    
   if(Norm(m)<SMALL_FIELD)dadbdc2ijk(m,mabc,abc);
   if(Norm(mabc)<SMALL_FIELD)ijk2dadbdc(mabc,m,abc);
   if(Norm(Pel)<SMALL_FIELD)dadbdc2ijk(Pel,Pelabc,abc);
   if(Norm(Pelabc)<SMALL_FIELD)ijk2dadbdc(Pelabc,Pel,abc);
}

return sta;
}

   // for xyt file
double physproperties::xytcols(float * nn,float * nnerr, int & nofcols,bool setnn,char * header,char * outstr,inipar & ini,int & verbose)
 {double sta=0;double * ptr;int * iptr;char hs[40];char num[40];
   header[0]='\0';outstr[0]='\0';int nofa=sps.na(),nofb=sps.nb(),nofc=sps.nc();
    int nofcolsin=0;if(!setnn){nofcolsin=nofcols;j=0;nofa=0;nofb=0;nofc=0;totalJ=0;}
    nofcols=12+nofcomponents();for(int i=8;i<=nofcolsin||i<=nofcols;++i)
    {ptr=NULL;iptr=NULL;
       switch(i) {       case 8: iptr=&j;snprintf(hs,40,"phasnumber-j");break;
                         case 9: iptr=&sps.wasstable;snprintf(hs,40,"period-key");break;
                         case 10: iptr=&nofa;snprintf(hs,40,"supercell-nr1");break;
                         case 11: iptr=&nofb;snprintf(hs,40,"nr2");break;
                         case 12: iptr=&nofc;snprintf(hs,40,"nr3");break;
                         default: ptr=&totalJ[i-12];snprintf(hs,40,"<I%i>",i-12);break;
                 }

if(setnn){if(nnerr[i]>0&&i<=nofcolsin)
          {if(ptr==NULL&&iptr==NULL)
           {
 if(verbose==1)fprintf(stdout,"sta_mcphas.fum warning: exp value %g cannot be fitted in column %i in file ./fit/mcphas.fum\n",nn[i],i);
           }else{
//fprintf(stdout,"stacalc_mphas.fum: col %i line %i value %g err %g - calcvalue %g\n",i,j2,nn[i],nnerr[i],clc);
            if(iptr==NULL)sta+=(nn[i]-(*ptr))*(nn[i]-(*ptr))/nnerr[i]/nnerr[i];
            else          sta+=(nn[i]-(*iptr))*(nn[i]-(*iptr))/nnerr[i]/nnerr[i];
           }
          }
          if(iptr==NULL)nn[i]=(*ptr);else nn[i]=(*iptr);
          switch(i){case 8: case 9:   snprintf(num,40,"%ip ",(*iptr));break;
                    case 10: case 11: snprintf(num,40,"%i x ",(*iptr));break;
                    case 12:          snprintf(num,40,"%i ",(*iptr));break;
                    default: snprintf(num,40,"%4.4g ",myround(nn[i]));
                   }
          snprintf(outstr+strlen(outstr),MAXNOFCHARINLINE-strlen(outstr),"%*s%s",(int)(strlen(hs)+1-strlen(num)),"",num); 
         }else
         {if(iptr==NULL)(*ptr)=nn[i]; else (*iptr)=nn[i];
         }
     snprintf(header+strlen(header),MAXNOFCHARINLINE-strlen(header),"%s ",hs);
    } // next i
return sta;
}

Vector physproperties::mu0M()// magnetisation mu0*M(Tesla) 
{Vector M(1,3);
  // M=sum m/V_unitcell
  // inputpars.cs.pVol() // V_unitcell in Aengstroem^3
  M=9.274010*1.25663706127*nofatoms()*m/cs.pVol(); // total moment per unit cell (muB)
           // muB=9.274010e-24 J/Tesla
           // mu0=1.25663706127(20)×10−6 N⋅A−2
 return M;
}
Vector physproperties::Pdiveps0()
{Vector P(1,3);
   // Pel in C/m2
    P=Pel/8.8541878188e-12; 
   // epsilon0 = 8.8541878188(14)×10−12 C2⋅kg−1⋅m−3⋅s2
 return P;
}
int physproperties::nofatoms()//
{return cs.nofatoms;}
int physproperties::nofcomponents()//
{return cs.nofcomponents;}

//*********************************************************************************************************
// methode save
double physproperties::save (int & verbose, const char * filemode, int & htfailed,inipar & ini, par & inputpars,char * prefix)
{ FILE *fout;
  char filename[MAXNOFCHARINLINE],str[MAXNOFCHARINLINE],outstr[MAXNOFCHARINLINE];
  time_t curtime;
  struct tm *loctime;  
  int i,j2,l,i1,j1,nmax;
  Vector null(1,nofcomponents()*nofatoms());null=0;
  Vector null1(1,3);null1=0;
  double sta=0;

  float nn[200];nn[0]=199;
  float nnerr[200];nnerr[0]=199;for(int i=1;i<=199;++i)nnerr[i]=0;
  int ortho=1;
  if (inputpars.cs.alpha()!=90||inputpars.cs.beta()!=90||inputpars.cs.gamma()!=90)
   {ortho=0;} 
   Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1;
                       abc(4)=inputpars.cs.alpha(); abc(5)=inputpars.cs.beta(); abc(6)=inputpars.cs.gamma();
   ijk2dadbdc(mabc,m,abc);  // transform m and P to abc coordinate system
   ijk2dadbdc(Pelabc,Pel,abc);
   ijk2dadbdc(Pelabc0,Pel0,abc);

  Vector M(1,3); M=mu0M();// magnetisation mu0*M(Tesla) 
  Vector P(1,3); P=Pdiveps0(); //  Polarisation/epsilon0 (V/m)

  printf("saving properties for ");ini.print_usrdefcols(stdout,x,y,T,H,inputpars.cs.abc,M,P,true);
  if(ortho==0){printf(" Hi=%g Hj=%g Hk=%g ",H(1),H(2),H(3));}
  ini.time_estimate_until_end(x,y);
  printf("\n");

//-----------------------------------mcphas.fum ----------------------------------------------------  
  errno = 0;char outfilename[MAXNOFCHARINLINE];
  int nofcols;
  strcpy(outfilename,"./results/");strcpy(outfilename+10,prefix);
  strcpy(outfilename+10+strlen(prefix),"mcphas.fum");
  if (verbose==1) printf("saving %s \n",outfilename);
  if (htfailed!=0){fe=0;u=0;m=0;m[1]=0;m[2]=0;m[3]=0;Eelastic=0;mabc=0;sps.epsilon=0;Pel=0;Pelabc=0;}
  fumcols(nn,nnerr,nofcols,true,str,outstr,ini,ortho,inputpars,verbose);
  if (washere==0)
  {fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.fum-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   fprintf (fout, "#note: - for specific heat calculation use unit conversion 1mev/f.u.=96.48J/mol\n");
   fprintf (fout, "#      - below moments and energies are given per ion - not per formula unit !\n");
   if(ini.nofrndtries<0)fprintf (fout, "# free energy calculation not (yet) implemented in Monte Carlo run - calculating energy u only\n");
   if(ini.doeps)
   {fprintf (fout, "#      - strain tensor eps is calculated selfconsistently.\n");
    if(ini.linepscf)
    { fprintf (fout, "#      ... however, singleion ion Hamiltonian is always diagonalised with eps=0 (option -linepscf).\n");
    }
    if(ini.linepsjj)
    { fprintf (fout, "#      ... however, two ion interaction is always evaluated for eps=0 (option -linepsjj).\n");
    }
   }
    if(ini.include_cd){ fprintf (fout, "#      ... two ion interaction in mcphas.j extended by classical dipolar interaction according to Bowdens expansion.\n");
    }
   if(fabs(inputpars.totalcharge)<SMALLCHARGE){fprintf(fout,"#! Structural Static Electrical Polarisation  Pel0a=%4.4g  Pel0b=%4.4g  Pel0c=%4.4g C/m^2\n",Pelabc0(1),Pelabc0(2),Pelabc0(3)); }
   if(ortho==0){fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
               if(fabs(inputpars.totalcharge)<SMALLCHARGE)fprintf(fout,"#! Structural Static Electrical Polarisation  Pel0i=%4.4g  Pel0j=%4.4g  Pel0k=%4.4g C/m^2\n",Pel0(1),Pel0(2),Pel0(3));
                }

   ini.print_usrdefcolhead(fout,str);
  fclose(fout);
  }

   fout = fopen_errchk (outfilename,"a");ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
         fprintf(fout,"%s\n",outstr); fclose(fout);

   strcpy(outfilename,"./results/.");strcpy(outfilename+11,prefix); strcpy(outfilename+11+strlen(prefix),"mcphas.fum");
   fout = fopen_errchk (outfilename,"a");ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
   fprintf(fout,"%s\n",outstr); fclose(fout);

   if((fout=fopen("./fit/mcphas.fum","rb"))!=NULL)
    {// some measured data should be fitted
    if(washere==0){fprintf(stdout,"#Mcphas- calculating standard deviation to ./fit/mcphas.fum - magnetisation ma mb mc (col 11,12,13)\n");}
     while(feof(fout)==0)
     {if ((l=inputline(fout,nn,nnerr))!=0)
      {if(ini.checkTH(nn,T,H,inputpars.cs.abc)) // checks if T and H is in accordance with nn  
         {sta=fumcols(nn,nnerr,l,true,str,outstr,ini,ortho,inputpars,verbose);
          if(verbose==1){fprintf(stdout,"sta_mcphas.fum=%g\n",sta);}
	 }
      }
     }
     fclose(fout);
    }else{errno=0;}



//--------------------------------------mcphas.xyt---------------------------------------------------  
  errno = 0;  
  strcpy(outfilename,"./results/");strcpy(outfilename+10,prefix);
  strcpy(outfilename+10+strlen(prefix),"mcphas.xyt");
    if (verbose==1)printf("saving %s\n",outfilename);
  if (htfailed!=0){j=0;totalJ=0;}
  xytcols(nn,nnerr,nofcols,true,str,outstr,ini,verbose);
  if (washere==0)
  {fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.xyt-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
  
   ini.print_usrdefcolhead(fout,str); 
   fclose(fout);    
     }
  fout = fopen_errchk (outfilename,"a");
   if(j<0){sps.wasstable=j;}// if qvector generated structure is stable, then take period key = number of qvector
   ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
    fprintf(fout,"%s\n",outstr);fclose(fout);
    if((fout=fopen("./fit/mcphas.xyt","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.xyt not implemented\n");}

     if(washere==0){fprintf(stdout,"#Mcphas- calculating standard deviation to ./fit/mcphas.xyt - <I>\n");}
     while(feof(fout)==0)
     {if ((l=inputline(fout,nn,nnerr))!=0)
      {if(ini.checkTH(nn,T,H,inputpars.cs.abc)) // checks if T and H is in accordance with nn  
         {double s=xytcols(nn,nnerr,l,true,str,outstr,ini,verbose);
          if(verbose==1){fprintf(stdout,"sta_mcphas.xyt=%g\n",s);}
          sta+=s;
	 }
      }
     }
     fclose(fout);
    
   }else{errno=0;}

//-----------------------------------mcphasj*.j*------------------------------------------------------  
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
 for(i=1;i<=nmax;++i){for(l=1;l<=nofatoms();++l){
  errno = 0;
  if (verbose==1)printf("saving mcphas%i.j%i - spinspin corr for sublattice %i neighbour %i\n",l,i,l,i);
  strcpy(outfilename,"./results/");strcpy(outfilename+10,prefix);
  strcpy(outfilename+10+strlen(prefix),"mcphas");
  snprintf(filename,MAXNOFCHARINLINE,"%s%i.j%i",outfilename,l,i);
  if (htfailed!=0){jj[i](1)=0;jj[i](2)=0;jj[i](3)=0;}
  if (washere==0)  //printout file header
  {  fout = fopen_errchk (filename,filemode);
   fprintf(fout, "#output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas*.j*-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   fprintf (fout, "# sublattice %i (da=%g a db=%g b dc=%g c)\n",l,(*inputpars.jjj[l]).xyz(1),(*inputpars.jjj[l]).xyz(2),(*inputpars.jjj[l]).xyz(3));
   fprintf (fout, "# correlation function <JJ(%g %g %g)>\n",myround((*inputpars.jjj[l]).dn[i](1)),myround((*inputpars.jjj[l]).dn[i](2)),myround((*inputpars.jjj[l]).dn[i](3)));
   str[0]='\0';
           for(i1=1;i1<=(*inputpars.jjj[l]).nofcomponents;++i1)
	      {snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"<J%cJ%c> ",'a'-1+i1,'a'-1+i1);}
           for(i1=1;i1<=(*inputpars.jjj[l]).nofcomponents-1;++i1)
              {for(j1=i1+1;j1<=(*inputpars.jjj[l]).nofcomponents;++j1)
                        {snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str),"<J%cJ%c> <J%cJ%c> ",'a'-1+i1,'a'-1+j1,'a'-1+j1,'a'-1+i1);
			}
	      }
                ini.print_usrdefcolhead(fout,str);   

if(ini.nofrndtries<0){fprintf(fout,"# Monte Carlo calculations of correlation functions not (yet) implemented\n"); }
       
   fclose(fout);
      }
 if(ini.nofrndtries>=0)
 {
  fout = fopen_errchk (filename,"a");
    ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
        for(j2=1;j2<=nofcomponents()*nofcomponents();++j2)               
            {fprintf (fout, "%4.4g ",myround(jj[i](j2+nofcomponents()*nofcomponents()*(l-1))));
	    }
   fprintf (fout,"\n");
   fclose(fout);
   snprintf(filename,MAXNOFCHARINLINE,"./fit/mcphas%i.j%i",l,i);
   if((fout=fopen(filename,"rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using %s  not implemented\n",filename);}
     fclose(fout);
    }else{errno=0;}
  }
  }}

//-----------------------------------mcphas.cel ----------------------------------------------------  
if(ini.cel!=0.0)
 {errno = 0;  
  strcpy(outfilename,"./results/");strcpy(outfilename+10,prefix);
  strcpy(outfilename+10+strlen(prefix),"mcphas.cel");
    if (verbose==1)printf("saving %s\n",outfilename);
    if (htfailed!=0){cel=0;}
  if (washere==0)
  {
   fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.cel-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   fprintf(fout,"# Elastic Constants\n");
      str[0]='\0';
   snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str), "  C11 C22  C33 C44 C55 C66 C12 C13 C14 C15 C16 C23 C24 C25 C26 C34 C35 C36 C45 C46 C56(GPa)");
   ini.print_usrdefcolhead(fout,str);
   fclose(fout);
   }
   strcpy(outfilename+10+strlen(prefix),"mcphas.cel"); 
   fout = fopen_errchk (outfilename,"a");
   if (verbose==1)printf(" .... saving elastic constants in %s\n",outfilename);
   ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
   for (i=1;i<=6;++i)
    { fprintf (fout, " %4.4g ",myround(cel(i,i)));
    }
   for (i=1;i<=5;++i)for(j=i+1;j<=6;++j)
    { fprintf (fout, " %4.4g ",myround(cel(i,j)));
    }
 
    fprintf(fout,"\n");
   fclose(fout);
  }

//-----------------------------------------mcphas*.hkl------------------------------------------------  
 errno = 0;
  strcpy(outfilename,"./results/");strcpy(outfilename+10,prefix);
  strcpy(outfilename+10+strlen(prefix),"mcphas*.hkl");
  if (verbose==1)printf("saving %s - neutrons and xrays\n",outfilename);
  if (washere==0)
  {//neutrons
   strcpy(outfilename+10+strlen(prefix),"mcphas.hkl");
  fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.hkl-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
   fprintf (fout, "#Neutron Intensity - Mind: only structure+polarizationfactor+formfactor+debeywallerfactor - no lorentzfactor is  taken into account\n");
   str[0]='\0';
   snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str), "       h   k   l  int       h   k   l   int       h   k   l   int ");
   ini.print_usrdefcolhead(fout,str);
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
   str[0]='\0';
   snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str), "       h   k   l  real(mi(Q)) im(mi(Q))    h   k   l   real(mi(Q)) im(mi(Q))     h   k   l   real(mi(Q)) im(mi(Q))[mu_B/atom] ...");
   }else{
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - a component\n"); 
   str[0]='\0';
   snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str), "       h   k   l  real(ma(Q)) im(ma(Q))      h   k   l    real(ma(Q)) im(ma(Q))      h   k   l  real(ma(Q)) im(ma(Q)) [mu_B/atom]...");
   }
   ini.print_usrdefcolhead(fout,str); 
   if(ini.nofrndtries<0){fprintf(fout,"# Monte Carlo calculations of hkl not (yet) implemented\n"); }

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
   str[0]='\0';
   snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str), "       h   k   l  real(mj(Q)) im(mj(Q))       h   k   l   real(mj(Q)) im(mj(Q))       h   k   l   real(mj(Q)) im(mj(Q))[mu_B/atom]...}");
   }else{
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - b component\n"); 
   str[0]='\0';
   snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str), "      h   k   l  real(mb(Q)) im(mb(Q))        h   k   l  real(mb(Q)) im(mb(Q))     h   k   l  real(mb(Q)) im(mb(Q))[mu_B/atom] ...}");
   }
   ini.print_usrdefcolhead(fout,str);
  if(ini.nofrndtries<0){fprintf(fout,"# Monte Carlo calculations of hkl not (yet) implemented\n"); }
   fclose(fout);
   //xray c component
   if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphask.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasc.hkl");  }
   fout = fopen_errchk (outfilename,filemode);
   fprintf(fout, "#output file of program %s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphasc.hkl-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
   fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(fout,"#**********************************************************\n");
      if(ortho==0){
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - k component\n");
   fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   str[0]='\0';
   snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str), "       h   k   l  real(mk(Q)) im(mk(Q))       h   k   l  real(mk(Q)) im(mk(Q))      h   k   l   real(mk(Q)) im(mk(Q)) [mu_B/atom]...");
   }else{
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - c component\n"); 
   str[0]='\0';
   snprintf(str+strlen(str),MAXNOFCHARINLINE-strlen(str), "       h   k   l  real(mc(Q)) im(mc(Q))       h   k   l  real(mc(Q)) im(mc(Q))       h   k   l  real(mc(Q)) im(mc(Q))  [mu_B/atom]...");
   }
   ini.print_usrdefcolhead(fout,str);
   if(ini.nofrndtries<0){fprintf(fout,"# Monte Carlo calculations of hkl not (yet) implemented\n"); }
   fclose(fout);

      }
   if(ini.nofrndtries>=0){
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
   ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
   for (i=nofhkls;i>=1;--i)
    {if (htfailed!=0){hkli[inew[i]](1)=0;hkli[inew[i]](2)=0;hkli[inew[i]](3)=0;hkli[inew[i]](4)=0;}
    fprintf (fout, "   %4.4g %4.4g %4.4g  %4.4g  ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](4)));
    } fprintf(fout,"\n");
   fclose(fout);
  //xray a component
  if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphasi.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasa.hkl");  }
   fout = fopen_errchk (outfilename,"a");
   if (verbose==1)printf(" .... saving %s\n",outfilename);
   ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
   for (i=nofhkls;i>=1;--i)
    {fprintf (fout, "   %4.4g %4.4g %4.4g  %4.4g %4.4g  ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](5)),myround(hkli[inew[i]](6)));
    } fprintf(fout,"\n");
   fclose(fout);
  //xray b component
   if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphasj.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasb.hkl");  }
   fout = fopen_errchk (outfilename,"a");
   if (verbose==1)printf(" .... saving %s\n",outfilename);
   ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
   for (i=nofhkls;i>=1;--i)
    {fprintf (fout, "   %4.4g %4.4g %4.4g  %4.4g %4.4g  ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](7)),myround(hkli[inew[i]](8)));
    } fprintf(fout,"\n");
   fclose(fout);
  //xray c component
   if(ortho==0){strcpy(outfilename+10+strlen(prefix),"mcphask.hkl");
                } else {strcpy(outfilename+10+strlen(prefix),"mcphasc.hkl");  }
   fout = fopen_errchk (outfilename,"a");
   if (verbose==1)printf(" .... saving %s\n",outfilename);
   ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
   for (i=nofhkls;i>=1;--i)
    {fprintf (fout, "   %4.4g %4.4g %4.4g  %4.4g %4.4g ",myround(hkli[inew[i]](1)),myround(hkli[inew[i]](2)),myround(hkli[inew[i]](3)),myround(hkli[inew[i]](9)),myround(hkli[inew[i]](10)));
    } fprintf(fout,"\n");
   fclose(fout);


    if((fout=fopen("./fit/mcphas.hkl","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.hkl not implemented\n");}
     fclose(fout);
    }else{errno=0;}
   delete []inew;delete []intensity;
  }
//---------------------------------------mcphas.sps--------------------------------------------------  
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
   if(ortho==0){fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   ini.print_usrdefcolcodes(fout);
   fprintf (fout, " nofspins nofatoms(in primitive basis) nofmeanfield-components errorcode(0=ok,1=failed) eps1=epsii eps2=epsjj eps3=epskk eps4=2epsjk eps5=2epsik eps6=2epsij\n");}
else
  {ini.print_usrdefcolcodes(fout);
   fprintf (fout, " nofspins nofatoms(in primitive basis) nofmeanfield-components errorcode(0=ok,1=failed) eps1=epsaa eps2=epsbb eps3=epscc eps4=2epsbc eps5=2epsac eps6=2epsab\n");}
   fprintf (fout, "    #<I1(atom 1)> <I1(atom 2)> .... selfconsistent Spinconfiguration  \n");
   fprintf (fout, "    #<I2(atom 1)> <I2(atom 2)> .... UNITS:  multiply <I>=<J> by Lande factor g to get moment [muB]\n");
   fprintf (fout, "    #<I3(atom 1)> <I3(atom 2)> ....}\n");
    fclose(fout);
   }  
  fout = fopen_errchk (outfilename,"a");
   ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
   fprintf (fout, " %i %i %i ",
            sps.n()*sps.nofatoms,sps.nofatoms,sps.nofcomponents);
   if (htfailed!=0){fprintf(fout,"1 ");int d1=1;sps.spinfromq(d1,d1,d1,null1,null,null,null);} // failed
    else {fprintf(fout,"0 ");}
   fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g\n",myround(sps.epsilon(1)),myround(sps.epsilon(2)),myround(sps.epsilon(3)),myround(sps.epsilon(4)),myround(sps.epsilon(5)),myround(sps.epsilon(6)));
    sps.print(fout);fprintf(fout,"\n");
   fclose(fout);
    if((fout=fopen("./fit/mcphas.sps","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.sps not implemented\n");}
     fclose(fout);
    }else{errno=0;}
 
//---------------------------------------mcphas.mf--------------------------------------------------  
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
   if(ortho==0){fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   ini.print_usrdefcolcodes(fout);
   fprintf (fout, " nofspins nofatoms(in primitive basis) nofmeanfield-components errorcode(0=ok,1=failed) eps1=epsii eps2=epsjj eps3=epskk eps4=2epsjk eps5=2epsik eps6=2epsij\n");}
else
  {ini.print_usrdefcolcodes(fout);
   fprintf (fout, " nofspins nofatoms(in primitive basis) nofmeanfield-components errorcode(0=ok,1=failed) eps1=epsaa eps2=epsbb eps3=epscc eps4=2epsbc eps5=2epsac eps6=2epsab\n");}
   fprintf (fout, "    #mf1(atom 1) mf1(atom 2) .... selfconsistent Mean field configuration \n"); 
   fprintf (fout, "    #mf2(atom 1) mf2(atom 2) .... UNITS: mf(atom i)=gJ*mu_B*hxc(atom i)[meV] \n"); 
   fprintf (fout, "    #mf3(atom 1) mf3(atom 2) ....         (i.e. divide by gJ and mu_B=0.05788meV/Tesla to get exchange field hxc[Tesla]}\n");
   if(ini.nofrndtries<0){fprintf(fout,"# Monte Carlo calculations of mean fields do not make sense - nothing saved\n"); }

    fclose(fout);
   }  
if(ini.nofrndtries>=0){
     fout = fopen_errchk (outfilename,"a");ini.print_usrdefcols(fout,x,y,T,H,inputpars.cs.abc,M,P,false);
fprintf (fout, " %i %i %i ",
            mf.n()*mf.nofatoms,mf.nofatoms,mf.nofcomponents);
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
                     }
//-----------------------------------------------------------------------------------------  
 
 washere=1;
return sta;
 }

// scroll output files and read physical properties from these if possible,
// on success return 0, otherwise
// return 1
int physproperties::read(int & verbose, par & inputpars,char * readprefix,inipar & ini)
{ FILE *fin;int n;float nnerr[200];nnerr[0]=199;for(int i=1;i<=199;++i)nnerr[i]=0;
  char filename[MAXNOFCHARINLINE],str[MAXNOFCHARINLINE],outstr[MAXNOFCHARINLINE];
  int i,j2,l,nmax;
  float nn[200];nn[0]=199;
  int ortho=1; bool found=0;
  if (inputpars.cs.alpha()!=90||inputpars.cs.beta()!=90||inputpars.cs.gamma()!=90){ortho=0;}
   Vector M(1,3),P(1,3); M=0;P=0; 
 if(verbose==1){ printf("reading properties for ");ini.print_usrdefcols(stdout,x,y,T,H,inputpars.cs.abc,M,P,true);}

//-----------------------------------------mcphas.fum------------------------------------------------  
// here read free energy etc if possible ... otherwise return 1
// check x y T H[1] H[2] H[3] agrees and fe nonzero ?
// then read fe, u, 
// read mabc[1-3] convert to ijk
// NOT read: strain tensor epsilon Pel
  errno = 0;char infilename[MAXNOFCHARINLINE];
  strcpy(infilename,"./results/");strcpy(infilename+10,readprefix);
  strcpy(infilename+10+strlen(readprefix),"mcphas.fum");
  if (verbose==1) printf("reading %s \n",infilename);
   fin = fopen_errchk (infilename,"r");

found=0;
   while(found==0){ 
 while (feof(fin)==0&&0==(n=inputline(fin,nn))){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
   found=ini.checkTH(nn,T,H,inputpars.cs.abc);
         }
   fclose(fin);
   fumcols(nn,nnerr,n,false,str,outstr,ini,ortho,inputpars,verbose); // load parameters 
   // check if stable, i.e. if free energy nn[8] is nonzero
   if(nn[8]==0) return 1; // 

//fprintf (fout, "%4.4g %4.4g  %4.4g %4.4g %4.4g %4.4g %4.4g       %8.8g            %8.8g       %4.4g    %4.4g %4.4g %4.4g    %4.4g",
 //           myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),
//myround(H[2]),myround(H[3]),myround(fe),myround(u),myround(Norm(m)),
//myround(mabc[1]),myround(mabc[2]),myround(mabc[3]),myround(m*Hijk/Norm(Hijk)));
 //  if(ortho==0){fprintf (fout, "    %4.4g %4.4g %4.4g   %4.4g %4.4g %4.4g",myround(m(1)),myround(m(2)),myround(m(3)),Hijk(1),Hijk(2),Hijk(3));}
  
//-----------------------------------------mcphas.xyt------------------------------------------------  
  errno = 0; 
  strcpy(infilename,"./results/");strcpy(infilename+10,readprefix);
  strcpy(infilename+10+strlen(readprefix),"mcphas.xyt");
    if (verbose==1)printf("reading %s\n",infilename);
  fin = fopen_errchk (infilename,"r");
// check x y T H[1] H[2] H[3] agrees and fe nonzero ?
//then read j ... not necessary
   //fprintf (fin, "%4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g       %ip           %ip      ",
     //       myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]),j,sps.wasstable);
        // then read totalJ ... not necessary if sps is read !!
//   for(i1=1;i1<=nofcomponents();++i1)
//	      {fprintf(fin,"%4.4g ",myround(totalJ(i1)));}
//	      fprintf(fin,"\n");
found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
   found=ini.checkTH(nn,T,H,inputpars.cs.abc);
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
 for(i=1;i<=nmax;++i){for(l=1;l<=nofatoms();++l){
  errno = 0;
  if (verbose==1)printf("reading mcphas%i.j%i - spinspin corr for sublattice %i neighbour %i\n",l,i,l,i);
  strcpy(infilename,"./results/");strcpy(infilename+10,readprefix);
  strcpy(infilename+10+strlen(readprefix),"mcphas");
  snprintf(filename,MAXNOFCHARINLINE,"%s%i.j%i",infilename,l,i);
  fin = fopen_errchk (filename,"r");
   // check x y T H[1] H[2] H[3] agrees and fe nonzero ?
//then read jj[] 

 // fprintf (fin, "%4.4g %4.4g   %4.4g %4.4g   %4.4g %4.4g %4.4g     ",myround(x),myround(y),myround(T),myround(Norm(Hijk)),myround(H[1]),myround(H[2]),myround(H[3]));
  //      for(j2=1;j2<=nofcomponents()*nofcomponents();++j2)               
   //         {fprintf (fin, "%4.4g ",myround(jj[i](j2+nofcomponents()*nofcomponents()*(l-1))));
//	    }
 //  fprintf (fin,"\n");
  found=0;while(found==0){ 
 while (feof(fin)==0&&0==inputline(fin,nn)){}  // if yes -> input them
   if (feof(fin)!=0) {fclose(fin);return 1;}
   found=ini.checkTH(nn,T,H,inputpars.cs.abc);
               } 
 fclose(fin);
 for(j2=1;j2<=nofcomponents()*nofcomponents();++j2)jj[i](j2+nofcomponents()*nofcomponents()*(l-1))=nn[j2+7];
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
   found=ini.checkTH(nn,T,H,inputpars.cs.abc);nofhkls=(i-7)/4;
                       } 
   fclose(fin);
update_maxnofhkls(nofhkls);
   for (i=nofhkls;i>=1;--i){
                            hkli[i](1)=nn[7+(i-1)*4+1];
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
     found=ini.checkTH(nn,T,H,inputpars.cs.abc);
               } 
   fclose(fin);
  for (i=nofhkls;i>=1;--i){
                           hkli[i](5)=nn[7+(i-1)*4+4];
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
    found=ini.checkTH(nn,T,H,inputpars.cs.abc);
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
   found=ini.checkTH(nn,T,H,inputpars.cs.abc);
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
    sps.load(fin);
   found=ini.checkTH(nn,T,H,inputpars.cs.abc);
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
   found=ini.checkTH(nn,T,H,inputpars.cs.abc);
               } 
   fclose(fin);
   //-----------------------------------------------------------------------------------------  
if(verbose==1)printf("..ok\n"); else printf("reading ");
return 0; 
 }
