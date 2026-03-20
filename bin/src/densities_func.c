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

int find_usrdef_out(const char * var,char **out)
// looks in the set of user defined variables for a specific variable and returns it's index
// if not found returns 0
{int i=0;
 for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col){
 if(strcmp(var,out[col])==0)i=col;
 }
 return i;
}

void extract_xyTHext(char * outstr,double & x,double & y,double & T,Vector & Hext,Vector & abc)
{  Vector Habc(1,3),Eabc(1,3);Hext=0;Habc=0;Eabc=0;
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
double distance_of_str_to_xyTHext_hklE(char * str,double & tcdd,double & x,double & y,double & T,Vector & Hext,
                     double  h,double  k,double l,double  E,Vector & abc)
{ double d,dd,tx=1e10,ty=1e10,tT,tE;Vector tHext(1,Hext.Hi()),thkl(1,3);
  extract_xyTHext(str,tx,ty,tT,tHext,abc);
  extract(str,"h",thkl(1)); 
  extract(str,"k",thkl(2)); 
  extract(str,"l",thkl(3));
  extract(str,"E",tE);
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
// spinconfigurations are (loaded from File handle fin_coq
// aa[0]==0: "nearest" means take the list of doubles stored in field aa and
// look in the first line of a mf configuration and compare numbers in this line
// and remembers in savmf and outstr the best match
// aa[0]>0: "nearest" means load exactly a[0] configurations and return the last in savmf and outstr
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
                snprintf(outstr,MAXNOFCHARINLINE,"%s=%g %s=%g %s=%g %s=%g %s=%g %s=%g %s=%g n=%g spins nr1=%i nr2=%i nr3=%i nofatoms=%i in primitive basis nofcomponents=%i",
                        out[1],myround(numbers[1]),out[2],myround(numbers[2]),out[3],myround(numbers[3]),out[4],myround(numbers[4]),
                        out[5],myround(numbers[5]),out[6],myround(numbers[6]),out[7],myround(numbers[7]),
                        myround(numbers[8]),spins.na(),spins.nb(),spins.nc(),(int)numbers[9],(int)numbers[10]);
        savmf=spins;x=0;y=0;
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
  }snprintf(outstr,MAXNOFCHARINLINE,"n=%i spins nr1=%i nr2=%i nr3=%i nofatoms=%i in primitive basis nofcomponents=%i",savmf.n()*savmf.nofatoms,savmf.na(),savmf.nb(),savmf.nc(),savmf.nofatoms,savmf.nofcomponents);
 }
 return 0; // ok structure found
}

// inputs file header and returns number of atoms 
int headerinput(FILE * fin_coq,FILE* fout,graphic_parameters & gp,cryststruct & cs,char **out,Matrix & N)
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
  extract(instr,"Nii",N[1][1]);
  extract(instr,"Nij",N[1][2]);
  extract(instr,"Nik",N[1][3]);
  extract(instr,"Njj",N[2][2]);
  extract(instr,"Njk",N[2][3]);
  extract(instr,"Nkk",N[3][3]);

   
   extract(instr,"scale_view_1",gp.scale_view_1);
   extract(instr,"scale_view_2",gp.scale_view_2);
   extract(instr,"scale_view_3",gp.scale_view_3);
   cs.cextract(instr);
   extract(instr,"nofatoms",cs.nofatoms);
   extract(instr,"nofcomponents",cs.nofcomponents);
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
N(2,1)=N(1,2);N(3,1)=N(1,3);N(3,2)=N(2,3);
return n;
}


void check_for_best_excitation_and_close(FILE *fin,graphic_parameters & gp,const char * ext,
 spincf & spinconf,spincf & spinconfev_real,spincf & spinconfev_imag,
double & x, double & y, double & T,
Vector & Hext,char*outstr,char*outhklstr,const char * oscill_type,Vector & thkl,Vector & hkl,
char **argv, int & os, Vector & abc, int nofatoms, int doijk,double xx, double yy, double zz,int dim
)
{// loads eigenvectors from handle fin and check which one is nearest -------------------------------
 // to x y T Hext, store in spinconfev_real/spinconfev_imag
 //  and then closes file handle fin
// input file header ------------------------------------------------------------------
char instr[MAXNOFCHARINLINE],dumstr[MAXNOFCHARINLINE];long int pos=0,jj;int i,j,k;
double delta,tcdd,checkdd=1e7,dd;             int extended_eigenvector_dimension;
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong %s file format\n",ext);exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
               // load evs and check which one is nearest -------------------------------   
               extract(instr,"spins_wave_amplitude",gp.spins_wave_amplitude);
               extract(instr,"spins_show_ellipses",gp.spins_show_ellipses);
               extract(instr,"spins_show_oscillation",gp.spins_show_oscillation);
               extract(instr,"phonon_wave_amplitude",gp.phonon_wave_amplitude);
               extract(instr,"phonon_scale_static_displacements",gp.phonon_scale_static_displacements);
               extract(instr,"extended_eigenvector_dimension",extended_eigenvector_dimension);
              }
             if(strcmp("qee/qsd/qod",ext)!=0)extended_eigenvector_dimension=3; // set this to 3 if not density plot
              jj=fseek(fin,pos,SEEK_SET);if (jj!=0){fprintf(stderr,"Error: wrong %s file format\n",ext);exit (EXIT_FAILURE);}
   
               for (delta=1000.0;feof(fin)==0&&fgets(instr,MAXNOFCHARINLINE,fin)!=NULL;)
               { if(fgets(dumstr,MAXNOFCHARINLINE,fin)!=NULL) 
                 {spincf ev_real(spinconf.na(),spinconf.nb(),spinconf.nc(),spinconf.nofatoms,extended_eigenvector_dimension);
                  spincf ev_imag(spinconf.na(),spinconf.nb(),spinconf.nc(),spinconf.nofatoms,extended_eigenvector_dimension);
                 ev_real.load(fin);ev_imag.load(fin);
                 dd=distance_of_str_to_xyTHext_hklE(instr,tcdd,x,y,T,Hext,
                        strtod(argv[5+os],NULL),strtod(argv[6+os],NULL),
                        strtod(argv[7+os],NULL),strtod(argv[8+os],NULL),abc);
                 if (dd<delta)
                 {delta=dd;checkdd=tcdd;hkl=thkl;//E=tE;
                  snprintf(outhklstr,MAXNOFCHARINLINE,"%s ",instr);
          
          if(doijk==3&&strcmp("qee/qsd/qod",ext)==0
             && (argv[1][1]=='s'||argv[1][1]=='o')
             ){
               // moments=xx*momentsx+yy*momentsy+zz*momentsz;
                         for(int ii=1;ii<=nofatoms;++ii)
                   for (i=1;i<=spinconf.na();++i)for(j=1;j<=spinconf.nb();++j)for(k=1;k<=spinconf.nc();++k)
                  for(int nt=1;nt<=dim;++nt){spinconfev_real.m(i,j,k)(nt+dim*(ii-1))=xx*ev_real.m(i,j,k)(nt+3*dim*(ii-1))+yy*ev_real.m(i,j,k)(nt+dim+3*dim*(ii-1))+zz*ev_real.m(i,j,k)(nt+2*dim+3*dim*(ii-1));
                                             spinconfev_imag.m(i,j,k)(nt+dim*(ii-1))=xx*ev_imag.m(i,j,k)(nt+3*dim*(ii-1))+yy*ev_imag.m(i,j,k)(nt+dim+3*dim*(ii-1))+zz*ev_imag.m(i,j,k)(nt+2*dim+3*dim*(ii-1));
                                        }
               } else
               {
                  spinconfev_real=ev_real;
                  spinconfev_imag=ev_imag; 
               }                 
                 }
                 pos=ftell(fin); 
                 fgets(instr,MAXNOFCHARINLINE,fin); 
                 while (instr[strspn(instr," \t")]=='#'&&feof(fin)==0) // pointer to 'ltrimstring' 
                  {pos=ftell(fin);fgets(instr,MAXNOFCHARINLINE,fin);}
                 jj=fseek(fin,pos,SEEK_SET);
               }}
              fclose (fin);
              check_dd(checkdd,outhklstr,outstr); // check if .mf and .ext agree
                fprintf(stdout,"#%s - %s - eigenvector\n",outstr,oscill_type);
              fprintf(stdout,"#real\n");
              spinconfev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              spinconfev_imag.print(stdout);

}
