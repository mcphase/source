/************************************************************************************/
// calculates polar coordinates from Vector X(1..3)
/************************************************************************************/
void jjjpar::getpolar(double x,double y, double z, double & r, double & th, double & ph)
{	 r=sqrt(x*x+y*y+z*z);
         th=acos(z/r);
	 if(sin(th)>=SMALL){
	                    if(x>0) {ph=acos(x/(r*sin(th))-SMALL);}
			    else    {ph=acos(x/(r*sin(th))+SMALL);}
			   }
			 else{ph=0;}
	 if (y<0){ph=2*PI-ph;}
}



/************************************************************************************/
// get parameters from sipf file
/************************************************************************************/
void jjjpar::get_parameters_from_sipfile(char * sipffilename)
{FILE * cf_file;
 int i,j;
 float nn[MAXNOFNUMBERSINLINE];
 nn[0]=MAXNOFNUMBERSINLINE;
  char modulefilename[MAXNOFCHARINLINE];

 char instr[MAXNOFCHARINLINE];
  cf_file = fopen_errchk (sipffilename, "rb");
  fgets_errchk (instr, MAXNOFCHARINLINE, cf_file);
  if(extract(instr,"MODULE=",modulefilename,(size_t)MAXNOFCHARINLINE))
   {if(extract(instr,"#!",modulefilename,(size_t)MAXNOFCHARINLINE))
    {fprintf(stderr,"Error: single ion property file %s does not start with '#!' or 'MODULE='\n",sipffilename);
     exit(EXIT_FAILURE);}
   }
   //ic1ion entered without path ?
      if (strncmp(modulefilename,"ic1ion",6)==0)
      {strcpy(modulefilename,getenv("MCPHASE_DIR"));
    //   strcat(modulefilename,"\\bin\\ic1ion_module\\ic1ion.so");
       strcat(modulefilename,"/bin/ic1ion_module/ic1ion.so");
}
      if (strncmp(modulefilename,"icf1ion",6)==0)
      {strcpy(modulefilename,getenv("MCPHASE_DIR")); strcat(modulefilename,"/bin/ic1ion_module/icf1ion.so"); }

  fprintf (stderr,"#parsing single ion property file: %s - loading module %s\n",sipffilename,modulefilename);

  if(strcmp(modulefilename,"kramer")==0)
    {module_type=1;fprintf (stderr,"[internal]\n");
      ABC=Vector(1,3);i=3;
      nof_electrons=0; // not to be used in module kramer !!
      while(feof(cf_file)==false)
      {fgets(instr, MAXNOFCHARINLINE, cf_file);
       if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
                                           i+=extract(instr,"A",ABC(1))-1;
                                           i+=extract(instr,"B",ABC(2))-1;
                                           i+=extract(instr,"C",ABC(3))-1;
                                          }
      }
      // input all  lines starting with comments
      if(i!=0){fprintf(stderr,"Error reading |<+-|Ja|-+>|,|<+-|Jb|-+>|,|<+-|Jc|+->| from file %s\ncorrect file format is:\n",sipffilename);
              fprintf(stderr,"\nMODULE=kramer\n#comment lines ..\n#matrix elements A=|<+-|Ja|-+>| B=|<+-|Jb|-+>| C=|<+-|Jc|+->|\nA=2 \nB=3 \nC=1\n\n");exit(EXIT_FAILURE);}
      // now we have the numbers corresponding to the vector ABC() in nn[]
      fprintf(stderr," ... kramers doublet with A=<+|Ja|->=%g B=<+-|Jb|+->=+-%g C=<+|Jc|->/i=%g\n",ABC(1),ABC(2),ABC(3));
      est=ComplexMatrix(0,2,1,2); // not used, just initialize to prevent errors
      Icalc_parstorage=ComplexMatrix(0,2,1,2); // not used, just initialize to prevent errors
    }
  else
    {if(strcmp(modulefilename,"brillouin")==0)
     {module_type=3;fprintf (stderr,"[internal]\n");
      ABC=Vector(1,1);i=1;
      nof_electrons=0; // not to be used in module brillouin !!
      while(feof(cf_file)==false)
      {fgets(instr, MAXNOFCHARINLINE, cf_file);
       if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
                                           i+=extract(instr,"J",ABC(1))-1;
                                          }
      }// input all  lines starting with comments
      if(i!=0){fprintf(stderr,"Error reading spin quantum number J=S from file %s\ncorrect file format is:\n",sipffilename);
              fprintf(stderr,"\n#!brillouin\n#comment lines ..\n# Quantum number  J\nJ=3.5\n\n");exit(EXIT_FAILURE);}
      // now we have the numbers corresponding to the vector ABC() in nn[]
      fprintf(stderr," ... Brillouin function with J=S=%g\n",ABC(1));
      est=ComplexMatrix(0,2,1,2);est=0;// not used, just initialize to prevent errors
      Icalc_parstorage=ComplexMatrix(0,2,1,2);Icalc_parstorage=0;// not used, just initialize to prevent errors
     }
     else
     {if(strcmp(modulefilename,"cfield")==0)
     {module_type=2;fprintf (stderr,"#[internal]\n");
      fclose(cf_file);cf_file = fopen_errchk (sipffilename, "rb"); // reopen file
       
      iops=new ionpars(cf_file);

      int dj;dj=(int)(2*J()+1);
      est=ComplexMatrix(0,dj,1,dj);
      Icalc_parstorage=ComplexMatrix(0,dj,1,dj);
      nof_electrons=(*iops).nof_electrons;
      // get 1ion parameters - operator matrices

     }
     else
     {if(strcmp(modulefilename,"so1ion")==0)
     {module_type=4;fprintf (stderr,"#[internal]\n");
      fclose(cf_file);cf_file = fopen_errchk (sipffilename, "rb"); // reopen file
     iops=new ionpars(cf_file);
      nof_electrons=(*iops).nof_electrons;
      int dj;dj=(int)(2*J()+1);
      est=ComplexMatrix(0,dj,1,dj);
      Icalc_parstorage=ComplexMatrix(0,dj,1,dj);
      // get 1ion parameters - operator matrices

     }
     else if (strcmp(modulefilename,"cluster")==0)
     {module_type=5;fprintf (stderr,"#[internal]\n");
      ABC=Vector(1,1);i=1;
      char clusterfilename[MAXNOFCHARINLINE];
      nof_electrons=0; // not to be used in module brillouin !!
      while(feof(cf_file)==false)
      {fgets(instr, MAXNOFCHARINLINE, cf_file);
       if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
                                           i+=extract(instr,"structurefile",clusterfilename,sizeof(clusterfilename))-1;
                                          }
      }// input all  lines starting with comments
      if(i!=0){fprintf(stderr,"Error reading structurefile from file %s\ncorrect file format is:\n",sipffilename);
              fprintf(stderr,"\n#!MODULE=cluster\n#comment lines ..\n# next line contains cluster structure filename\nstructurefile=cluster.j\n\n");exit(EXIT_FAILURE);}
      fprintf(stderr," ... reading cluster structure from %s\n",clusterfilename);
      clusterpars =new par(clusterfilename);
      est=ComplexMatrix(0,2,1,2);est=0;// not used, just initialize to prevent errors
      Icalc_parstorage=ComplexMatrix(0,2,1,2);Icalc_parstorage=0;// not used, just initialize to prevent errors
     }
     else
      {fprintf (stderr,"#[external]\n");
      i=0;nof_electrons=0;
      while(feof(cf_file)==false)
      {fgets(instr, MAXNOFCHARINLINE, cf_file);
       if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
                                           i-=extract(instr,"MODPAR1",nn[1])-1;
                                           i-=extract(instr,"MODPAR2",nn[2])-1;
                                           i-=extract(instr,"MODPAR3",nn[3])-1;
                                           i-=extract(instr,"MODPAR4",nn[4])-1;
                                           i-=extract(instr,"MODPAR5",nn[5])-1;
                                           i-=extract(instr,"MODPAR6",nn[6])-1;
                                           i-=extract(instr,"MODPAR7",nn[7])-1;
                                           i-=extract(instr,"MODPAR8",nn[8])-1;
                                           i-=extract(instr,"MODPAR9",nn[9])-1;
                                              extract(instr,"nof_electrons",nof_electrons);
                                          }
      }
       // input all  lines starting with comments
    //while((i=inputparline ("params",cf_file, nn))==0&&feof(cf_file)==false);
    // now we have the numbers corresponding to vector ABC() in nn[] - these are the module parameters !
    fprintf(stderr,"#parameters: ");
    if(i>0){
             ABC=Vector(1,i);for(j=1;j<=i;++j){ABC(j)=nn[j];fprintf(stderr,"%g ",nn[j]);}
            }else{
             ABC=Vector(1,1);
	    }
    fprintf(stderr,"\n");
module_type=0;
#ifdef __MINGW32__
  handle=LoadLibrary(modulefilename);
  if ((intptr_t)handle<= HINSTANCE_ERROR){fprintf (stderr, "jjjpar::jjjpar - Could not load dynamic library\n");
	       exit (EXIT_FAILURE);
	      }

    I=(void(*)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,double*,double*,ComplexMatrix*))GetProcAddress(handle,"Icalc");
    //*(int **)(&m)=GetProcAddress(handle,"Icalc");
     if (I==NULL) {fprintf (stderr,"jjjpar::jjjpar error %d  module %s loading function Icalc not possible\n",(int)GetLastError(),modulefilename);exit (EXIT_FAILURE);}
    m=(void(*)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"mcalc");
    //*(int **)(&m)=GetProcAddress(handle,"mcalc");
     if (m==NULL) {fprintf (stderr,"jjjpar::jjjpar error %d  module %s loading function mcalc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    L=(void(*)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"Lcalc");
    //*(int **)(&L)=GetProcAddress(handle,"Lcalc");
     if (S==NULL) {fprintf (stderr,"jjjpar::jjjpar error %d  module %s loading function Lcalc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    S=(void(*)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"Scalc");
    //*(int **)(&S)=GetProcAddress(handle,"Scalc");
     if (m==NULL) {fprintf (stderr,"jjjpar::jjjpar error %d  module %s loading function Scalc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    dm=(int(*)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*))GetProcAddress(handle,"du1calc");
    //*(void **)(&dm)=GetProcAddress(handle,"du1calc");
     if (dm==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function du1calc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    mq=(void(*)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*))GetProcAddress(handle,"mq");
    //*(void **)(&mq)=GetProcAddress(handle,"mq");
     if (mq==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function mq not possible - continuing\n",(int)GetLastError(),modulefilename);}
    estates=(void(*)(ComplexMatrix*,Vector*,Vector*,double*,double*,Vector*,char**))GetProcAddress(handle,"estates");
    //*(void **)(&estates)=GetProcAddress(handle,"estates");
     if (estates==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function estates not possible - continuing\n",(int)GetLastError(),modulefilename);
                                est=ComplexMatrix(0,2,1,2);// not used, just initialize to prevent errors
                                est=0;
                               }
    Icalc_parameter_storage=(void(*)(ComplexMatrix*,Vector*,Vector*,double*,double*,Vector*,char**))GetProcAddress(handle,"Icalc_parameter_storage_matrix_init");
    //*(void **)(&Icalc_parameter_storage)=GetProcAddress(handle,"Icalc_parameter_storage_matrix_init");
    if (Icalc_parameter_storage==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %X  module %s loading function Icalc_parameter_storage_matrix_init not possible - continuing\n",(int)GetLastError(),modulefilename);
                                  Icalc_parstorage=ComplexMatrix(0,2,1,2);Icalc_parstorage=0;// not used, just initialize to prevent errors
                                  }

  ddnn=(int(*)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexVector*))GetProcAddress(handle,"dv1calc");
    //*(void **)(&dnn)=GetProcAddress(handle,"dv1calc");
     if (ddnn==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function dv1calc not possible - continuing\n",(int)GetLastError(),modulefilename);}

    sd_m=(void(*)(Vector*,int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"spindensity_coeff");
    //*(void **)(&sd_m)=GetProcAddress(handle,"spindensity_coeff");
    if (sd_m==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function spindensity_coeff not possible - continuing\n",(int)GetLastError(),modulefilename);}

    od_m=(void(*)(Vector*,int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"orbmomdensity_coeff");
    //*(void **)(&od_m)=GetProcAddress(handle,"orbmomdensity_coeff");
    if (od_m==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function orbmomdensity_coeff not possible - continuing\n",(int)GetLastError(),modulefilename);}

    ro_calc=(void(*)(double*,double*,double*,double*,Vector*,double*,Vector*,Vector*,double*,Vector*,char**))GetProcAddress(handle,"ro_calc");
    //*(void **)(&sd_m)=GetProcAddress(handle,"spindensity_coeff");
    if (ro_calc==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function ro_calc not possible - continuing\n",(int)GetLastError(),modulefilename);}

#else
  char * error;
  handle=dlopen (modulefilename,RTLD_NOW | RTLD_GLOBAL);
  if (!handle){fprintf (stderr, "jjjpar::jjjpar - Could not load dynamic library\n");
               if ((error=dlerror())!=NULL)
	         {fprintf (stderr,"%s\n",error);}
	       exit (EXIT_FAILURE);
	      }
 *(void **)(&I)=dlsym(handle,"Icalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s\n",error);exit (EXIT_FAILURE);}
 *(void **)(&m)=dlsym(handle,"mcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);m=NULL;}
 *(void **)(&L)=dlsym(handle,"Lcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);L=NULL;}
 *(void **)(&S)=dlsym(handle,"Scalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);S=NULL;}
  *(void **)(&dm)=dlsym(handle,"du1calc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dm=NULL;}
  *(void **)(&mq)=dlsym(handle,"mq");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);mq=NULL;}
  *(void **)(&estates)=dlsym(handle,"estates");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);estates=NULL;
                                est=ComplexMatrix(0,2,1,2);est=0;// not used, just initialize to prevent errors
                               }
 *(void **)(&Icalc_parameter_storage)=dlsym(handle,"Icalc_parameter_storage_matrix_init");

  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);Icalc_parameter_storage=NULL;
                                Icalc_parstorage=ComplexMatrix(0,2,1,2);Icalc_parstorage=0;// not used, just initialize to prevent errors
                               }

  *(void **)(&ddnn)=dlsym(handle,"dv1calc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);ddnn=NULL;}

  *(void **)(&sd_m)=dlsym(handle,"spindensity_coeff");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);sd_m=NULL;}

  *(void **)(&od_m)=dlsym(handle,"orbmomdensity_coeff");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);od_m=NULL;}

  *(void **)(&ro_calc)=dlsym(handle,"ro_calc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);ro_calc=NULL;}



#endif
     }
    }
   }
  }
  fclose(cf_file);

  magFFj0=Vector(1,7);magFFj0=0;  magFFj0[1]=1;
  magFFj2=Vector(1,7);magFFj2=0;
  magFFj4=Vector(1,7);magFFj4=0;
  magFFj6=Vector(1,7);magFFj6=0;
  Zc=Vector(1,7);Zc=0;

   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;

  DWF=0;  gJ=0;

  cf_file = fopen_errchk (sipffilename, "rb");
  while(feof(cf_file)==false)
  {fgets(instr, MAXNOFCHARINLINE, cf_file);
   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
    extract(instr,"SCATTERINGLENGTHREAL",SLR);
    extract(instr,"SCATTERINGLENGTHIMAG",SLI);
    extract(instr,"GJ",gJ);
    extract(instr,"gJ",gJ);
    // read formfactor if given
    extract(instr,"FFj0A",magFFj0[1]);
    extract(instr,"FFj0a",magFFj0[2]);
    extract(instr,"FFj0B",magFFj0[3]);
    extract(instr,"FFj0b",magFFj0[4]);
    extract(instr,"FFj0C",magFFj0[5]);
    extract(instr,"FFj0c",magFFj0[6]);
    extract(instr,"FFj0D",magFFj0[7]);
    extract(instr,"FFj2A",magFFj2[1]);
    extract(instr,"FFj2a",magFFj2[2]);
    extract(instr,"FFj2B",magFFj2[3]);
    extract(instr,"FFj2b",magFFj2[4]);
    extract(instr,"FFj2C",magFFj2[5]);
    extract(instr,"FFj2c",magFFj2[6]);
    extract(instr,"FFj2D",magFFj2[7]);
    extract(instr,"FFj4A",magFFj4[1]);
    extract(instr,"FFj4a",magFFj4[2]);
    extract(instr,"FFj4B",magFFj4[3]);
    extract(instr,"FFj4b",magFFj4[4]);
    extract(instr,"FFj4C",magFFj4[5]);
    extract(instr,"FFj4c",magFFj4[6]);
    extract(instr,"FFj4D",magFFj4[7]);
    extract(instr,"FFj6A",magFFj6[1]);
    extract(instr,"FFj6a",magFFj6[2]);
    extract(instr,"FFj6B",magFFj6[3]);
    extract(instr,"FFj6b",magFFj6[4]);
    extract(instr,"FFj6C",magFFj6[5]);
    extract(instr,"FFj6c",magFFj6[6]);
    extract(instr,"FFj6D",magFFj6[7]);
   // coefficients of Z(K') according to Lovesey chapter 11.6.1 page 233
    extract(instr,"Z1c0",Zc(1));
    extract(instr,"Z1c2",Zc(2));
    extract(instr,"Z3c2",Zc(3));
    extract(instr,"Z3c4",Zc(4));
    extract(instr,"Z5c4",Zc(5));
    extract(instr,"Z5c6",Zc(6));
    extract(instr,"Z7c6",Zc(7));
   // read debeywallerfactor if given
    extract(instr,"DWF",DWF);
   // read radial wavefunction parameters
        extract(instr,"N1",Np(1));extract(instr,"XI1",Xip(1));extract(instr,"C1",Cp(1));
        extract(instr,"N2",Np(2));extract(instr,"XI2",Xip(2));extract(instr,"C2",Cp(2));
        extract(instr,"N3",Np(3));extract(instr,"XI3",Xip(3));extract(instr,"C3",Cp(3));
        extract(instr,"N4",Np(4));extract(instr,"XI4",Xip(4));extract(instr,"C4",Cp(4));
        extract(instr,"N5",Np(5));extract(instr,"XI5",Xip(5));extract(instr,"C5",Cp(5));
        extract(instr,"N6",Np(6));extract(instr,"XI6",Xip(6));extract(instr,"C6",Cp(6));
        extract(instr,"N7",Np(7));extract(instr,"XI7",Xip(7));extract(instr,"C7",Cp(7));
        extract(instr,"N8",Np(8));extract(instr,"XI8",Xip(8));extract(instr,"C8",Cp(8));
        extract(instr,"N9",Np(9));extract(instr,"XI9",Xip(9));extract(instr,"C9",Cp(9));
  }
 }

 fclose (cf_file);
// check gJ
if(module_type==2&&fabs(gJ-(*iops).gJ)>0.00001)
{fprintf(stderr,"Error internal module cfield : Lande Factor read from %s (gJ=%g) does not conform to internal module value gJ=%g\n",sipffilename,gJ,(*iops).gJ);exit(EXIT_FAILURE);}
if(module_type==4&&fabs(gJ-(*iops).gJ)>0.00001)
{fprintf(stderr,"Error internal module so1ion : Lande Factor read from %s (gJ=%g) does not conform to internal module value gJ=%g\n",sipffilename,gJ,(*iops).gJ);exit(EXIT_FAILURE);}
if (gJ==0){printf("# reading gJ=0 in single ion property file %s -> entering intermediate coupling mode by assigning Ja=Sa Jb=La Jc=Sb Jd=Lb Je=Sc Jf=Lc (S... Spin, L... angular momentum)\n",sipffilename);
           if (module_type==1){fprintf(stderr,"Error internal module kramers: intermediate coupling not supported\n");exit(EXIT_FAILURE);}
           if (module_type==2){fprintf(stderr,"Error internal module cfield : intermediate coupling not supported\n");exit(EXIT_FAILURE);}
           if (module_type==3){fprintf(stderr,"Error internal module brillouin: intermediate coupling not supported\n");exit(EXIT_FAILURE);}
           if (module_type==4){fprintf(stderr,"Error internal module so1ion : intermediate coupling not supported\n");exit(EXIT_FAILURE);}
          }

}

/****************************************************************************/
// get list of indices of exchange parameters
/************************************************************************************/
int jjjpar::get_exchange_indices(char *instr, Matrix *exchangeindices)
{
   bool charflag=false;
   char *tk,*tkp,*instrptr,sep[]=" \t\n",allowedch[]="abcdefghijklmnopqrstuvwxyz";
   int num_indices=0,ii,i,j;

   // Checks if using "JaJb" or "1,2" syntax
   if(strstr(instr,"J")!=NULL) charflag=true; else if(strstr(instr,",")==NULL) {
      fprintf(stderr,"Error in indexexchange: Syntax neither of the form JaJb or 1,2\n"); exit(EXIT_FAILURE); }

   // Moves to start of index list in the string
   instrptr = strstr(instr,"indexexchange")+13; instrptr = strstr(instrptr,"=")+1; instrptr+=strspn(instrptr,sep);
   // Clears all whitespaces at the end of the list
   tkp = strrchr(instr,0)-1; while(tkp>instrptr) { if(tkp[0]!=' '&&tkp[0]!='\t'&&tkp[0]!='\n') break; tkp--; } tkp++;
   // Goes through string finding whitespaces to get number of indices
   tk = strpbrk(instrptr,sep); if(strncmp(sep,tk,1)==0) tk += strspn(tk,sep);
   if(tk!=NULL) { num_indices=1; 
     while(tk!=NULL&&tk<tkp) { tk = strpbrk(tk+1,sep); if(strncmp(sep,tk,1)==0) tk += strspn(tk,sep); num_indices++; } 
   } else return 0; 

   (*exchangeindices) = Matrix(1,num_indices,1,2);

   if(charflag)
   {
      tk = strpbrk(instrptr,"J");
      for(ii=1; ii<=num_indices; ii++)
      {
         tk = strpbrk(tk+1,allowedch); i=(int)tk[0]-96; // 'a'==97 in ASCII
         tk = strpbrk(tk+1,allowedch); j=(int)tk[0]-96;
         (*exchangeindices)(ii,1) = i; (*exchangeindices)(ii,2) = j;
         tk = strpbrk(tk,"J");
      }
   }
   else
   {
      tk = instrptr;
      for(ii=1; ii<=num_indices; ii++)
      {
         i = strtol(tk,&tkp,10); tkp++;
	 j = strtol(tkp,&tk,10); 
         (*exchangeindices)(ii,1) = i; (*exchangeindices)(ii,2) = j;
	 tk = strpbrk(tk+1,"123456789");
      }
   }
   return num_indices;
}


/****************************************************************************/
// function to calculate calculate expectation values <Ialpha> alpha=1...nofcomponents
// from exchange field gjmbHxc [meV] and external field Hext
// this is the heart of the meanfield algorithm an it is necessary to
// keep this routine as efficient as possible
/****************************************************************************/
void jjjpar::Icalc (Vector &mom, double & T, Vector &  gjmbHxc,Vector & Hext ,double & lnZ,double & U,ComplexMatrix & parstorage)
{switch (module_type)
  {case 1: kramer(mom,T,gjmbHxc,Hext,lnZ,U);break;
   case 2:
   case 4: (*iops).cfieldJJ(mom,T,gjmbHxc,Hext,lnZ,U,parstorage);break;
   case 3: brillouin(mom,T,gjmbHxc,Hext,lnZ,U);break;
   case 5: cluster_Icalc(mom,T,gjmbHxc,Hext,lnZ,U);break;
   default: (*I)(&mom,&T,&gjmbHxc,&Hext,&gJ,&ABC,&cffilename,&lnZ,&U,&parstorage);
  }
}

/****************************************************************************/
// this function returns n (the number of transitions in the single ion susceptibility)
// the transition matrix mat first eigenvector u1 corresponding to jjjpar.transitionnumber and delta
// for effective field heff and temperature given on input
/****************************************************************************/
int jjjpar::du1calc(double & T,Vector &  gjmbHxc,Vector & Hext,ComplexVector & u1,float & delta,ComplexMatrix & ests)
{ switch (module_type)
  {case 0: if (dm!=NULL){return (*dm)(&transitionnumber,&T,&gjmbHxc,&Hext,&gJ,&ABC,&cffilename,&u1,&delta,&ests);}
           else return 0;
           break;
   case 1: return kramerdm(transitionnumber,T,gjmbHxc,Hext,u1,delta);break;
   case 2:
   case 4: return (*iops).cfielddm(transitionnumber,T,gjmbHxc,Hext,u1,delta,ests);break;
   case 3: return brillouindm(transitionnumber,T,gjmbHxc,Hext,u1,delta);break;
   case 5: return cluster_dm(transitionnumber,T,gjmbHxc,Hext,u1,delta);break;
   default: return 0;
  }
}


/****************************************************************************/
// returns eigenvalues, boltzmann population and eigenstates matrix parameters of ion
/****************************************************************************/
ComplexMatrix & jjjpar::eigenstates (Vector &  gjmbHxc,Vector & Hext,double & T)
{switch (module_type)
  {case 0:  if(estates!=NULL){(*estates)(&est,&gjmbHxc,&Hext,&gJ,&T,&ABC,&cffilename);}
            return est;break;
   case 2:
   case 4: (*iops).cfeigenstates(&est,gjmbHxc,Hext,T);return est;break;
   default: est=0;return est;
  }
}

/****************************************************************************/
// returns eigenvalues, boltzmann population and eigenstates matrix parameters of ion
/****************************************************************************/
ComplexMatrix & jjjpar::Icalc_parameter_storage_init (Vector &  gjmbHxc,Vector & Hext,double & T)
{switch (module_type)
  {case 0:  if(Icalc_parameter_storage!=NULL){(*Icalc_parameter_storage)(&Icalc_parstorage,&gjmbHxc,&Hext,&gJ,&T,&ABC,&cffilename);}
            return Icalc_parstorage;break;
   case 2:
   case 4: (*iops).cfeigenstates(&Icalc_parstorage,gjmbHxc,Hext,T);return Icalc_parstorage;break;
   default: Icalc_parstorage=0;return Icalc_parstorage;
  }
}
/****************************************************************************/
// returns operator matrices (n=0 Hamiltonian, n=1,...,nofcomponents: operators of moment components)
/****************************************************************************/
Matrix jjjpar::opmat(int n,Vector &  gjmbHxc,Vector & Hext)
{switch (module_type)
  {case 1:  return krameropmat(n,gjmbHxc,Hext);break;
   default: fprintf(stderr,"ERROR operator calculation in module jjjpar - opmat function not defined for module %i\n",module_type);exit(EXIT_FAILURE);
  }
}





