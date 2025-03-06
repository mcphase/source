/*****************************************************************/
// functions for programs charges, spindensity, orbmomdensity momdensity
//                        currdensities
/*****************************************************************/

double check_distance(double * aa,float * numbers)
{double dd;
dd=0.000001;for(int i=1;i<=NOF_USERDEF_MCPHAS_COLS&&i<numbers[0];++i)
       {if(aa[i]!=1e100)dd+=(aa[i]-numbers[i])*(aa[i]-numbers[i]);
       }
      dd=sqrt(dd);
return dd;
}


void extract_xyTHext(char * outstr,double & x,double & y,double & T,Vector & Hext,Vector & abc)
{  Vector Habc(1,3),Eabc(1,3);x=0;y=0;Hext=0;Habc=0;Eabc=0;
  extract(outstr,"x",x); 
  extract(outstr,"y",y); 
  extract(outstr,"T",T); 
  extract(outstr,"Ha",Habc[1]); 
  extract(outstr,"Hb",Habc[2]);
  extract(outstr,"Hc",Habc[3]); 
  extract(outstr,"Ea",Eabc[1]); 
  extract(outstr,"Eb",Eabc[2]);
  extract(outstr,"Ec",Eabc[3]); 
  extract(outstr,"Hi",Hext[1]); 
  extract(outstr,"Hj",Hext[2]);
  extract(outstr,"Hk",Hext[3]); 
  extract(outstr,"Ei",Hext[4]); 
  extract(outstr,"Ej",Hext[5]);
  extract(outstr,"Ek",Hext[6]); 
  extract(outstr,"s1",Hext[7]); 
  extract(outstr,"s2",Hext[8]);
  extract(outstr,"s3",Hext[9]); 
  extract(outstr,"s4",Hext[10]); 
  extract(outstr,"s5",Hext[11]);
  extract(outstr,"s6",Hext[12]); 
crosscheck_H_E(Hext,Habc,Eabc,abc); // put all information to Hext
}

bool check_dd(double tcdd,char * outhklstr,char * outstr)
{if(tcdd>1e-7)
               {fprintf(stderr,"Error program spins - inconsistent output files mcphas.mf and mcdisp.*:  temperature/magnetic/electric field/stress in static configuration\n" 
                               "(from results/mcphas.mf) is different (tcdd=%g) from mcdisp output files results/mcdisp.qes,qem,qeo:\n %s\n"
                               "... probably you need to rerun setup_mcdisp_mf and mcdisp with\n %s\n",tcdd,outhklstr,outstr);exit(EXIT_FAILURE);}
 return true;
} 

// determines how good are xyT Hext h k l E to loaded string parameters  
// returns also tcdd: distance of xyT Hext only (without hklE). returns extracted hkl E
double distance_of_str_to_xyTHext_hklE(char * outstr,double & tcdd,double & x,double & y,double & T,Vector & Hext,
                     double  h,double  k,double l,double  E,Vector & abc)
{ double d,dd,tx=1e10,ty=1e10,tT,tE;Vector tHext(1,Hext.Hi()),thkl(1,3);
  extract_xyTHext(outstr,tx,ty,tT,tHext,abc);
  extract(outstr,"h",thkl(1)); 
  extract(outstr,"k",thkl(2)); 
  extract(outstr,"l",thkl(3));
  extract(outstr,"E",tE);
  tHext-=Hext; 
  d=tHext*tHext;
  if(tx!=1e10){dd=x-tx;d+=dd*dd;}
  if(ty!=1e10){dd=y-ty;d+=dd*dd;}
  dd=T-tT;d+=dd*dd;tcdd=sqrt(d+1e-15);
  dd=h-thkl(1);d+=dd*dd;
  dd=k-thkl(2);d+=dd*dd;
  dd=l-thkl(3);d+=dd*dd;
  dd=E-tE;d+=dd*dd;
  d=sqrt(d+1e-15);
  return d;
}

int check_for_best(FILE *fin_coq,double * aa, spincf & savmf,double & x, double & y, double & T,Vector & Hext,char*outstr,char **out,Vector & abc)
{// load mfconfigurations and check which one is nearest -------------------------------
// returns 0 if ok, 1 if no stable configuration was found
int n,j;
   double dd,delta;
 float numbers[21];numbers[9]=1;numbers[10]=3;
 numbers[0]=20;char instr[MAXNOFCHARINLINE];
 long int pos=0;
 if(aa[0]==0) // indicates we should look for spinconfiguration which fits to T,H,E, ...
{
 for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=NOF_USERDEF_MCPHAS_COLS+1   //error in line reading (8 old format, 9 new format)
		    ;)
    { spincf spins(1,1,1,(int)numbers[9],(int)numbers[10]);
      if(spins.load(fin_coq)==1){
      dd=check_distance(aa,numbers);
      if(n>=11){if((int)numbers[11]!=0)dd=delta+10;} // if mcphase failed do not use this structure
      if(n>=17){for(int ii=1;ii<=6;++ii)spins.epsilon(ii)=numbers[11+ii];}
      if (dd<delta)
       {delta=dd;
                snprintf(outstr,MAXNOFCHARINLINE,"%s=%g %s=%g %s=%g %s=%g %s=%g %s=%g %s=%g n=%g spins nofatoms=%i in primitive basis nofcomponents=%i",
                        out[1],myround(numbers[1]),out[2],myround(numbers[2]),out[3],myround(numbers[3]),out[4],myround(numbers[4]),
                        out[5],myround(numbers[5]),out[6],myround(numbers[6]),out[7],myround(numbers[7]),
                        myround(numbers[8]),(int)numbers[9],(int)numbers[10]);
        savmf=spins;
        extract_xyTHext(outstr,x,y,T,Hext,abc);


       }
      pos=ftell(fin_coq); 
                 fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                 while (instr[strspn(instr," \t")]=='#'&&feof(fin_coq)==0) // pointer to 'ltrimstring' 
                  {pos=ftell(fin_coq);fgets(instr,MAXNOFCHARINLINE,fin_coq);}
     j=fseek(fin_coq,pos,SEEK_SET);
     if (j!=0){fprintf(stderr,"Error loading: wrong  file format\n");exit (EXIT_FAILURE);}
   
    } }if(delta==1000.0)return 1; // no stable structure found
 } 
   else
 {// look for config number -Tin
  for(n=1;n<=(int)aa[0];++n)
  {if(savmf.load(fin_coq)==0){fprintf(stderr,"Error program spins: loading configuration number %i\n",n);exit(1); }
  }snprintf(outstr,MAXNOFCHARINLINE,"n=%i spins nofatoms=%i in primitive basis nofcomponents=%i",savmf.n()*savmf.nofatoms,savmf.nofatoms,savmf.nofcomponents);
 }
 return 0; // ok structure found
}

// inputs file header and returns number of atoms 
int headerinput(FILE * fin_coq,FILE* fout,graphic_parameters & gp,cryststruct & cs,char **out)
{ char instr[MAXNOFCHARINLINE];
 // default out values
strcpy(out[1],"x");
strcpy(out[2],"y");
strcpy(out[3],"T");
strcpy(out[4],"H");
strcpy(out[5],"Ha");
strcpy(out[6],"Hb");
strcpy(out[7],"Hc");

 long int pos=0,j;int n=0;
cs.nofatoms=0;cs.nofcomponents=3;
char *token;cs.abc=0;
  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring'
  { pos=ftell(fin_coq);
   if (pos==-1)
       {fprintf(stderr,"Error: wrong mf/sps/tst file format\n");exit (EXIT_FAILURE);}
   fgets_errchk(instr,MAXNOFCHARINLINE,fin_coq);
   // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
   int i;
   for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
   // strip /r (dos line feed) from line if necessary
    while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

   if (instr[strspn(instr," \t")]=='#'){fprintf(fout,"%s",instr);}
   
   extract(instr,"show_abc_unitcell",gp.show_abc_unitcell);
   extract(instr,"show_primitive_crystal_unitcell",gp.show_primitive_crystal_unitcell);
   extract(instr,"show_magnetic_unitcell",gp.show_magnetic_unitcell);
   extract(instr,"show_atoms",gp.show_atoms);
   extract(instr,"spins_scale_moment",gp.spins_scale_moment);
   extract(instr,"show_chargedensity",gp.show_density);
   extract(instr,"out1",out[1], 20 ,1);
   extract(instr,"out2",out[2], 20 ,1);
   extract(instr,"out3",out[3], 20 ,1);
   extract(instr,"out4",out[4], 20 ,1);
   extract(instr,"out5",out[5], 20 ,1);
   extract(instr,"out6",out[6], 20 ,1);
   extract(instr,"out7",out[7], 20 ,1);
   
   extract(instr,"scale_view_1",gp.scale_view_1);
   extract(instr,"scale_view_2",gp.scale_view_2);
   extract(instr,"scale_view_3",gp.scale_view_3);
   cs.cextract(instr);
   extract(instr,"nofatoms",cs.nofatoms);    extract(instr,"nofcomponents",cs.nofcomponents);
   if ((cs.nofatoms>0)&&((extract(instr,"x",cs.x[n+1])+
                    extract(instr,"y",cs.y[n+1])+
  		       extract(instr,"z",cs.z[n+1])==0)||
		       (extract(instr,"da",cs.x[n+1])+
                   extract(instr,"db",cs.y[n+1])+
		       extract(instr,"dc",cs.z[n+1])==0)))
		  {++n;if(n>cs.nofatoms||cs.nofatoms>cs.maxnofatoms)
                    {fprintf(stderr,"ERROR reading file:maximum number of atoms in unit cell exceeded\n");exit(EXIT_FAILURE);}
                   cs.sipffilenames[n]=new char[MAXNOFCHARINLINE];
                   extract(instr,"sipffilename",cs.sipffilenames[n],(size_t)MAXNOFCHARINLINE,1);
//		   printf("%s\n",cs.sipffilenames[n]);
// HERE take care about atoms sitting at nearly the same position (a nucleus and a mangetic shell
// from makenn -cfph ... and put the positions exactly at the same correct value 
// (makenn had shifted the magnetic charge cloud by 0.01 A along c in order to enable 
// the correct evaluation of cf-phonon interactions by mcphas and mcdisp. here we correct for
// this shift in order to get the right output of charge densities and in spins.out the
// right positions for doing mcdiff !!
                   for(int i=1;i<n;++i){//loop all atoms which have been read
                                        if ( fabs((cs.x[n]-cs.x[i])*cs.abc[1])<0.2 &&
                                             fabs((cs.y[n]-cs.y[i])*cs.abc[2])<0.2 &&
                                             fabs((cs.z[n]-cs.z[i])*cs.abc[3])<0.2 )
                                            { // atom i and n are the same atom therefore check if they
                                              // are displaced along c and move magnetic atom back to nuclear position
                                              if(fabs((cs.x[n]-cs.x[i])*cs.abc[1])>0.001 ||
                                                 fabs((cs.y[n]-cs.y[i])*cs.abc[2])>0.001){fprintf(stderr,"Error spins.c: atoms %i (%g %g %g) and %i (%g %g %g) too close\n",i,cs.x[i],cs.y[i],cs.z[i],n,cs.x[n],cs.y[n],cs.z[n]);exit(EXIT_FAILURE);}
                                              if(cs.z[n]>cs.z[i]){cs.z[n]=cs.z[i];}else{cs.z[i]=cs.z[n];}
                                            }
                                       }
                  }
  }
    j=fseek(fin_coq,pos,SEEK_SET);
    if (j!=0){fprintf(stderr,"Error: wrong mf file format\n");exit (EXIT_FAILURE);}
return n;
}
