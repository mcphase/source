/*****************************************************************/
// functions for programs spins
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
                     double  h,double  k,double l,double  E,Vector & abc,Vector & thkl,double &tE)
{ double d,dd,tx=1e10,ty=1e10,tT;Vector tHext(1,Hext.Hi());
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

int check_for_best(FILE *fin_coq,double * aa, spincf & savmf,double & x, double & y, double & T,Vector & Hext,char*outstr,char **out,Vector & abc,int nomorethanthis)
{// load mfconfigurations and check which one is nearest -------------------------------
// returns 0 if ok, 1 if no stable configuration was found
// spinconfigurations are (loaded from File handle fin_coq
// aa[0]==0: "nearest" means take the list of doubles stored in field aa and
// look in the first line of a mf configuration and compare numbers in this line
// and remembers in savmf and outstr the best match
// aa[0]>0: "nearest" means load exactly aa[0] configurations and return the last in savmf and outstr
int n,j;
   double dd,delta;
 float numbers[21];numbers[9]=1;numbers[10]=3;
 numbers[0]=20;char instr[MAXNOFCHARINLINE];
 long int pos=0;int ctr=1;
 if(aa[0]==0) // indicates we should look for spinconfiguration which fits to T,H,E, ...
{
 for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=NOF_USERDEF_MCPHAS_COLS+1   //error in line reading (8 old format, 9 new format)
		    &&ctr<=nomorethanthis;)
    { spincf spins(1,1,1,(int)numbers[9],(int)numbers[10]);
      if(spins.load(fin_coq)==1){++ctr;
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
 {// look for config number aa[0]
  for(n=1;n<=(int)fabs(aa[0]);++n)
  {pos=ftell(fin_coq); // process comment lines until new structure begins
                 fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                 while (instr[strspn(instr," \t")]=='#'&&feof(fin_coq)==0) // pointer to 'ltrimstring' 
                  {pos=ftell(fin_coq);fgets(instr,MAXNOFCHARINLINE,fin_coq);}
     j=fseek(fin_coq,pos,SEEK_SET);
     if (j!=0){fprintf(stderr,"Error loading: wrong  file format\n");exit (EXIT_FAILURE);}
   if(aa[0]<0)inputline(fin_coq,numbers); // input a parameter line if it is a .sps or .mf file
   if(savmf.load(fin_coq)==0){return 1; }
         snprintf(outstr,MAXNOFCHARINLINE,"n=%i spins nr1=%i nr2=%i nr3=%i nofatoms=%i in primitive basis nofcomponents=%i",savmf.n()*savmf.nofatoms,savmf.na(),savmf.nb(),savmf.nc(),savmf.nofatoms,savmf.nofcomponents);
   if(aa[0]<0){snprintf(outstr,MAXNOFCHARINLINE,"%s=%g %s=%g %s=%g %s=%g %s=%g %s=%g %s=%g n=%g spins nr1=%i nr2=%i nr3=%i nofatoms=%i in primitive basis nofcomponents=%i",
                        out[1],myround(numbers[1]),out[2],myround(numbers[2]),out[3],myround(numbers[3]),out[4],myround(numbers[4]),
                        out[5],myround(numbers[5]),out[6],myround(numbers[6]),out[7],myround(numbers[7]),
                        myround(numbers[8]),savmf.na(),savmf.nb(),savmf.nc(),(int)numbers[9],(int)numbers[10]);
              extract_xyTHext(outstr,x,y,T,Hext,abc);}
  }
  
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


bool check_for_best_excitation_and_close(FILE *fin,graphic_parameters & gp,const char * ext,
 spincf & spinconf,spincf & spinconfev_real,spincf & spinconfev_imag,
double x, double  y, double  T,
Vector & Hext,char*outstr,char*outhklstr,const char * oscill_type,Vector & hkl,double &E,
int argc,char **argv, int os, Vector & abc, int nofatoms, int doijk,double xx, double yy, double zz,int dim
)
{// loads eigenvectors from handle fin and check which one is nearest -------------------------------
 // to x y T Hext,
 // 1) hkl store in spinconfev_real/spinconfev_imag,
 // 2) checks if xyTHex agrees with information loaded from handle fin (if T>0)
 // 3) if T<=0: store -Tth degenerate eigenvector in spinconfev  
 //  and then closes file handle fin
 // input:
 // fin           file handle
 // gp            graphic parameters
 // ext           extention (any value, will trigger density projection storage if equal to "qee/qsd/qod")
 // spinconf      spinconfiguration of supercell - only dimension of supercell will be used here (.na .nb .nc .nofatoms) 
  // x,y,T        x,y, Temperature,
  // Hext         external field,
 // outstr        string describing the x y T Hext  --> all as read from mcphas.mf
  // argv(os+5,..,os+7) hkl to look for
 // argv(os+8)    Energy to look for (optional), if not take E and look for it
 // doijk,xx,yy,zz,dim for densities-projection parameters
 // oscill_type   string describing the oscillation type (for printout on screen only)
 // output:
 // spinconfev_real, spinconfev_imag
 //               to store eigenvector
 // hkl           best matching hkl
 // E             best matching energy E
 // outhklstr     information string with best matching T Hext hkl 
 // true          degeneracy found
// input file header ------------------------------------------------------------------
char instr[MAXNOFCHARINLINE],dumstr[MAXNOFCHARINLINE];long int pos=0,jj;int i,j,k;
Vector thkl(1,3);double tE;int deg=0;
double delta,tcdd,checkdd=1e7,dd;             int extended_eigenvector_dimension;
double Edesired;
       if(argc==os+8){Edesired=E;}else{Edesired=strtod(argv[8+os],NULL);}
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
                        strtod(argv[7+os],NULL),Edesired,abc,thkl,tE);
                 if(dd==delta){deg++;}
                 if (dd<delta||(dd==delta&&deg<-T+1))
                 {if(dd<delta)deg=0;
                  delta=dd;checkdd=tcdd;hkl=thkl;E=tE;
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
              if(T>0)check_dd(checkdd,outhklstr,outstr); // check if .mf and .ext agree
               if(argc>os+8){ fprintf(stdout,"#%s - %s - eigenvector\n",outhklstr,oscill_type);
              fprintf(stdout,"#real\n");
              spinconfev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              spinconfev_imag.print(stdout);}
if(deg>-T)return true; else return false;

}


int jvx_create(LoadedInput& LI,
    int argc,
    char** argv,
    spincf & densitycf,
    spincf & spinconf,
    spincf & sc_phonon,
    Vector & h,
    Vector & hh,
    double & x,
    double & y)
{
             Vector hkl(1,3);hkl=0;double E;E=0;
             Vector gjmbHxc(1,3);gjmbHxc=0;
             spincf densityev_real(densitycf*0.0);
             spincf densityev_imag(densitycf*0.0);
             spincf spinconfev_real(spinconf*0.0);
             spincf spinconfev_imag(spinconf*0.0);
             spincf spinconfpev_real(spinconf*0.0);
             spincf spinconfpev_imag(spinconf*0.0);
            // to do jvx output of static structure put zeros into these spinconfigurations

// check sipffilenames and put radius= ... in case single ion module is
//  capable of calculating position 
for(int ii=1;ii<=LI.inputpars->cs.nofatoms;++ii)
{Vector pos(1,3);   
  for(int nt=1;nt<=LI.inputpars->cs.nofcomponents;++nt){h(nt)=hh(nt+LI.inputpars->cs.nofcomponents*(ii-1));}
  if(true==(*LI.inputpars->jjj[ii]).pcalc(pos,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage))
 {double charge;charge=(*LI.inputpars->jjj[ii]).charge;if(charge==0)charge=0.01;
   snprintf(LI.cs.sipffilenames[ii],MAXNOFCHARINLINE,"pointcharge %g |e| radius=%g",
           charge,LI.gp.scale_pointcharges*0.529177*signum(charge)*pow((double)fabs(charge),0.3333));
  // if we can determine from sipffilename the element --> put r g b information
    
  unsigned int r,g,b; char element []  ="E\0";double radius=0.3333;
   element[0]=toupper((*LI.inputpars->jjj[ii]).sipffilename[0]);
   if(isalpha((*LI.inputpars->jjj[ii]).sipffilename[1]))
   element[1]=tolower((*LI.inputpars->jjj[ii]).sipffilename[1]);
  for(int i=0;i<NOFELEMENTS_RGB;++i)
   if(0==strncmp(element,elstr[i],2)){
r=(unsigned int)(unsigned char)elstr[i][2];
g=(unsigned int)(unsigned char)elstr[i][3];
b=(unsigned int)(unsigned char)elstr[i][4];
radius=strtod(elstr[i]+5,NULL);
   snprintf(LI.cs.sipffilenames[ii],MAXNOFCHARINLINE,"pointcharge %g |e| radius=%g %s r=%i g=%i b=%i",
           charge,LI.gp.scale_pointcharges*radius*signum(charge)*pow((double)fabs(charge),0.3333),
           element,r,g,b);
//    printf("%s scale=%g radius=%g charge=%g sig(charge)=%g pow=%g  estr=%s\n",LI.cs.sipffilenames[ii],
//LI.gp.scale_pointcharges,radius,charge,signum(charge),pow((double)fabs(charge),0.3333),elstr[i]+5);

}

// printf("#! atom %i %s displacement u%ix=%g u%iy=%g u%iz=%g A\n",ii,LI.cs.sipffilenames[ii],ii,pos(1),ii,pos(2),ii,pos(3));
  printf("#! atom %i  displacement u%ix=%g u%iy=%g u%iz=%g A\n",ii,ii,pos(1),ii,pos(2),ii,pos(3));
 }
}

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   LI.fin = fopen_errchk ("./results/spins.jvx", "w");
    LI.gp.showprim=0;LI.gp.spins_wave_amplitude=0;
     densitycf.jvx_cd(LI.fin,LI.outstr,LI.cs,LI.gp,0.0,densityev_real,densityev_imag,hkl,LI.T,gjmbHxc,LI.Hext,LI.cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon);
    fclose (LI.fin);

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   LI.fin = fopen_errchk ("./results/spins_prim.jvx", "w");
     LI.gp.showprim=1;
     densitycf.jvx_cd(LI.fin,LI.outstr,LI.cs,LI.gp,0.0,densityev_real,densityev_imag,hkl,LI.T,gjmbHxc,LI.Hext,LI.cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon);
    fclose (LI.fin);

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
// try movie - a spinwave picture  ... including phonons and spindensity changes ...
//***************************************************************************************************************
//***************************************************************************************************************
// check if there is a second 
for(int i=LI.os+1;i<argc;++i){
if(strcmp(argv[i],"-prefix")==0){strcpy(LI.prefix,argv[i+1]); // read alternative LI.prefix
                                   fprintf(stdout,"# prefix for qee qep etc input filenames: %s\n",LI.prefix);
 				   LI.os+=2;}
                       }
if(argc-1==NOF_USERDEF_MCPHAS_COLS+4+LI.os)LI.os+=NOF_USERDEF_MCPHAS_COLS-4;
if(argc-1==6+LI.os)LI.os-=2;  // 
//argc-1==8+LI.os  ... ok 
if (argc-1==8+LI.os){
              // double E;
              char outhklstr[MAXNOFCHARINLINE];
             LI.gp.spins_wave_amplitude=1.0;LI.gp.spins_show_ellipses=1.0;LI.gp.spins_show_oscillation=1.0;
             LI.gp.phonon_wave_amplitude=1.0;LI.gp.phonon_scale_static_displacements=1.0;
//----------------------------------------------------------------------------------------------------------
 
           if(LI.arrow>0){
             strcpy(LI.infilename,"./results/");strcpy(LI.infilename+10,LI.prefix);
             switch(LI.arrow)
             {case 1: strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qes");LI.fin = fopen(LI.infilename, "rb");
                      if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qes", "rb");
                      break;
              case 2: strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qeo");LI.fin = fopen(LI.infilename, "rb");
                      if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qeo", "rb");
                      break;
              case 3: strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qem");LI.fin = fopen(LI.infilename, "rb");
                      if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qem", "rb");
                      break;
              case 4: fprintf(stderr,"mcdisp: output of individual moment oscillation in eigenvector file mcdisp.qemi not yet implemented - thus exiting program spins\n");
                      strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qemi");LI.fin = fopen(LI.infilename, "rb");
                      if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qemi", "rb");
                      break;   
              case 5: strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qpe");LI.fin = fopen(LI.infilename, "rb");
                      if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qpe", "rb");
                      break;                   
             }
             check_for_best_excitation_and_close(LI.fin,LI.gp,"qev",spinconf,spinconfev_real,spinconfev_imag,
                          x,y,LI.T,LI.Hext,LI.outstr,outhklstr,"moment oscillation",hkl,E,argc,argv,LI.os,
                          LI.inputpars->cs.abc,LI.inputpars->cs.nofatoms,LI.doijk, LI.xx,LI.yy,LI.zz,LI.dim);

             }//arrow
     
//----------------------------------------------------------------------------------------------------------
            if(LI.phonon>0){
             strcpy(LI.infilename,"./results/");strcpy(LI.infilename+10,LI.prefix);
             strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qep");LI.fin = fopen(LI.infilename, "rb");
                      if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qep", "rb");
             check_for_best_excitation_and_close(LI.fin,LI.gp,"qep",spinconf,spinconfpev_real,spinconfpev_imag,
                          x,y,LI.T,LI.Hext,LI.outstr,outhklstr,"phonon LI.oscillation",hkl,E,argc,argv,LI.os,
                          LI.inputpars->cs.abc,LI.inputpars->cs.nofatoms,LI.doijk, LI.xx,LI.yy,LI.zz,LI.dim);

             }//phonon
//----------------------------------------------------------------------------------------------------------
            if(LI.density){
             strcpy(LI.infilename,"./results/");strcpy(LI.infilename+10,LI.prefix);
             switch(argv[1][1]) // dimension deLI.finition from jjjpar.hpp
                {case 'c': strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qee");LI.fin = fopen(LI.infilename, "rb");
                           if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qee", "rb");
                           break;
                 case 's': strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qsd");LI.fin = fopen(LI.infilename, "rb");
                           if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qsd", "rb");
                           break;
                 case 'o': strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qod");LI.fin = fopen(LI.infilename, "rb");
                           if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qod", "rb");
                           break;
                 case 'm': fprintf(stderr,"Error spins: magnetic moment density oscillation not yet implemented\n");return false;break;
                             // would have to look into qsd and qod files !!
                 case 'j': strcpy(LI.infilename+10+strlen(LI.prefix),"mcdisp.qod");LI.fin = fopen(LI.infilename, "rb");
                           if(LI.fin==NULL)LI.fin = fopen_errchk ("./results/mcdisp.qod", "rb");
                           break;
                }
             check_for_best_excitation_and_close(LI.fin,LI.gp,"qee/qsd/qod",densitycf,densityev_real,densityev_imag,
          x,y,LI.T,LI.Hext,LI.outstr,outhklstr,"density oscillation",hkl,E,argc,argv,LI.os,
          LI.inputpars->cs.abc,LI.inputpars->cs.nofatoms,LI.doijk, LI.xx,LI.yy,LI.zz,LI.dim);

             }//LI.gp.show_density
//----------------------------------------------------------------------------------------------------------           
              LI.gp.read();// read graphic parameters which are set by user in file results/graphic_parameters.set
                        // in case he wants to overwrite some default settings
              if(LI.arrow==0)LI.gp.spins_scale_moment=0;
              if(LI.density==0)LI.gp.show_density=0;
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              double phase;
             // complex <double> im(0,1);
              for(int i=0;i<LI.gp.nof_frames;++i)
              {phase=2*3.1415*i/(LI.gp.nof_frames-1);
               printf("\n********************************************\n");
               printf(" calculating movie sequence %i(%i)\n",i+1,LI.gp.nof_frames);
               printf("********************************************\n");
               char filename[MAXNOFCHARINLINE];
               snprintf(filename,MAXNOFCHARINLINE,"./results/spins.%02d.jvx",i+1);
               LI.fin = fopen_errchk (filename, "w");LI.gp.showprim=0;
                     densitycf.jvx_cd(LI.fin,outhklstr,LI.cs,LI.gp,
                                  phase,densityev_real,densityev_imag,hkl,LI.T,hh,LI.Hext,LI.cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon,spinconfpev_real,spinconfpev_imag);
               fclose (LI.fin);
               snprintf(filename,MAXNOFCHARINLINE,"./results/spins_prim.%02d.jvx",i+1);
               LI.fin = fopen_errchk (filename, "w");LI.gp.showprim=1;
                     densitycf.jvx_cd(LI.fin,outhklstr,LI.cs,LI.gp,
                                  phase,densityev_real,densityev_imag,hkl,LI.T,hh,LI.Hext,LI.cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon,spinconfpev_real,spinconfpev_imag);
               fclose (LI.fin);
              }
          printf("# %s\n",outhklstr);
          }
fprintf(stderr,"# %s\n",LI.outstr);
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program spins running at \n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/spins.jvx\n");
fprintf(stderr,"# * java JvxAnimator \"results/spins.*.jvx\" \n");
fprintf(stderr,"# ************************************************************************\n");
 return 0;
}


void section4header (FILE * fout)
{
fprintf(fout,"\
#\n\
# %%SECTION 4%% DESCRIPTION OF MAGNETIC UNIT CELL AND LIST OF MAGNETIC ATOMS\n\
#\n\
#\n\
# here follows the description of the magnetic unit cell with respect\n\
# to the primitive crystallographic unit cell:\n\
# 'nr1', 'nr2', 'nr3' ...the crystallographic unit cell has to be taken\n\
#                        nr1 nr2 and nr3 times along r1 r2 and r3,\n\
#                        respectively to get magnetic unit cell\n\
# 'nat' denotes the number of magnetic atoms in magnetic unit cell\n\
#\n\
# It follows a list of nat lines with to describe the magnetic moment configuration\n\
# Notes:\n\
# 'atom-filename' means the single ion property filename of this magnetic atom:\n\
#                 -it must contain the Formfactor Coefficients (e.g. see international tables)\n\
#                                      Lande factor\n\
#                                      Neutron Scattering Length (10^-12 cm) \n\
#                 -it may contain a    Debey Waller Factor\n\
# 'da' 'db' and 'dc' are not used by the program (unless you enter a line #! use_dadbdc=1)\n\
# 'dr1','dr2' and 'dr3' refer to the primitive lattice given below\n\
# 'Ma','Mb','Mc' denote the magnetic moment components in Bohr magnetons\n\
#                in case of non orthogonal lattices instead of Ma Mb Mc the components Mx My Mz\n\
#                have to be given, which refer to an right handed orthogonal coordinate system \n\
#                defined by y||b, z||(a x b) and x normal to y and z\n\
#  <Sa> <Sb> <Sc>  <La> <Lb > <Lc>  (optional) denote the spin and orbital angular momentum components \n\
# 'mf1' 'mf2' 'mf3' (optional line, used to go beyond dipole approx for formfactor)\n\
#                                     denote the corresponding exchange fields in meV\n\
#\n");
}

