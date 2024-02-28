/*****************************************************************/
// functions for programs charges, spindensity, orbmomdensity momdensity
//                        currdensities
/*****************************************************************/

int check_for_best(FILE *fin_coq,double Tin, double hain,double hbin, double hcin, spincf & savmf, double & T,Vector & Hext,char*outstr)
{// load mfconfigurations and check which one is nearest -------------------------------
// returns 0 if ok, 1 if no stable configuration was found
int n;
   double ddT,ddHa,ddHb,ddHc,dd,delta;
 float numbers[20];numbers[9]=1;numbers[10]=3;
 numbers[0]=20;char instr[MAXNOFCHARINLINE];
 long int pos=0;
 if(Tin>=0)
{
 for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)
    { spincf spins(1,1,1,(int)numbers[9],(int)numbers[10]);
      if(spins.load(fin_coq)==1){
     if(Tin==0){ddT=0; // here hain and hbin correspond to x and y in phasediagram
                ddHa=hain-numbers[1];ddHa*=ddHa;
                ddHb=hbin-numbers[2];ddHb*=ddHb;
                ddHc=0;
      } else {
      ddT=Tin-numbers[3];ddT*=ddT;
      ddHa=hain-numbers[5];ddHa*=ddHa;
      ddHb=hbin-numbers[6];ddHb*=ddHb;
      ddHc=hcin-numbers[7];ddHc*=ddHc;
      }
      dd=sqrt(ddT+ddHa+ddHb+ddHc+0.000001);
      if(n>=11){if((int)numbers[11]!=0)dd=delta+10;} // if mcphase failed do not use this structure
      if(n>=17){for(int ii=1;ii<=6;++ii)spins.epsilon(ii)=numbers[11+ii];}
      if (dd<delta)
       {delta=dd;
        snprintf(outstr,MAXNOFCHARINLINE,"x=%g y=%g T=%g Ha=%g Hb=%g Hc=%g n=%g spins nofatoms=%i in primitive basis nofcomponents=%i",myround(numbers[1]),myround(numbers[2]),myround(numbers[3]),myround(numbers[5]),myround(numbers[6]),myround(numbers[7]),myround(numbers[8]),(int)numbers[9],(int)numbers[10]);
        savmf=spins;T=numbers[3];Hext(1)=numbers[5];Hext(2)=numbers[6];Hext(3)=numbers[7];
       }
      pos=ftell(fin_coq); 
                 fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                 while (instr[strspn(instr," \t")]=='#'&&feof(fin_coq)==0) // pointer to 'ltrimstring' 
                  {pos=ftell(fin_coq);fgets(instr,MAXNOFCHARINLINE,fin_coq);}
       fseek(fin_coq,pos,SEEK_SET);
    } }if(delta==1000.0)return 1; // no stable structure found
 } 
   else
 {// look for config number -Tin
  for(n=1;n<=-Tin;++n)
  {if(savmf.load(fin_coq)==0){fprintf(stderr,"Error program spins: loading configuration number %i\n",n);exit(1); }
  }snprintf(outstr,MAXNOFCHARINLINE,"n=%i spins nofatoms=%i in primitive basis nofcomponents=%i",savmf.n()*savmf.nofatoms,savmf.nofatoms,savmf.nofcomponents);
 }
 return 0; // ok structure found
}
   
int headerinput(FILE * fin_coq,FILE* fout,graphic_parameters & gp,cryststruct & cs)
{ char instr[MAXNOFCHARINLINE];
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
   if(cs.abc[1]==0){extract(instr,"a",cs.abc[1]);extract(instr,"b",cs.abc[2]); extract(instr,"c",cs.abc[3]);
                 extract(instr,"alpha",cs.abc[4]);  extract(instr,"beta",cs.abc[5]);extract(instr,"gamma",cs.abc[6]);
   }
   extract(instr,"show_abc_unitcell",gp.show_abc_unitcell);
   extract(instr,"show_primitive_crystal_unitcell",gp.show_primitive_crystal_unitcell);
   extract(instr,"show_magnetic_unitcell",gp.show_magnetic_unitcell);
   extract(instr,"show_atoms",gp.show_atoms);
   extract(instr,"spins_scale_moment",gp.spins_scale_moment);
   extract(instr,"show_chargedensity",gp.show_density);

   extract(instr,"scale_view_1",gp.scale_view_1);
   extract(instr,"scale_view_2",gp.scale_view_2);
   extract(instr,"scale_view_3",gp.scale_view_3);

   extract(instr,"r1x",cs.r[1][1]);extract(instr,"r2x",cs.r[1][2]); extract(instr,"r3x",cs.r[1][3]);
   extract(instr,"r1y",cs.r[2][1]); extract(instr,"r2y",cs.r[2][2]); extract(instr,"r3y",cs.r[2][3]);
   extract(instr,"r1z",cs.r[3][1]); extract(instr,"r2z",cs.r[3][2]); extract(instr,"r3z",cs.r[3][3]);
   extract(instr,"r1a",cs.r[1][1]);extract(instr,"r2a",cs.r[1][2]); extract(instr,"r3a",cs.r[1][3]);
   extract(instr,"r1b",cs.r[2][1]); extract(instr,"r2b",cs.r[2][2]); extract(instr,"r3b",cs.r[2][3]);
   extract(instr,"r1c",cs.r[3][1]); extract(instr,"r2c",cs.r[3][2]); extract(instr,"r3c",cs.r[3][3]);
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
