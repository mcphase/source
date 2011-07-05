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
      mcalc_parstorage=ComplexMatrix(0,2,1,2); // not used, just initialize to prevent errors
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
      mcalc_parstorage=ComplexMatrix(0,2,1,2);mcalc_parstorage=0;// not used, just initialize to prevent errors
     }
     else
     {if(strcmp(modulefilename,"cfield")==0)
     {module_type=2;fprintf (stderr,"#[internal]\n");
      fclose(cf_file);cf_file = fopen_errchk (sipffilename, "rb"); // reopen file
       
      iops=new ionpars(cf_file);

      int dj;dj=(int)(2*J()+1);
      est=ComplexMatrix(0,dj,1,dj);
      mcalc_parstorage=ComplexMatrix(0,dj,1,dj);
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
      mcalc_parstorage=ComplexMatrix(0,dj,1,dj);
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
      mcalc_parstorage=ComplexMatrix(0,2,1,2);mcalc_parstorage=0;// not used, just initialize to prevent errors
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

    m=(void(*)(Vector*,double*,Vector*,double*,Vector*,char**,double*,double*,ComplexMatrix*))GetProcAddress(handle,"mcalc");
    //*(int **)(&m)=GetProcAddress(handle,"mcalc");
     if (m==NULL) {fprintf (stderr,"jjjpar::jjjpar error %d  module %s loading function mcalc not possible\n",(int)GetLastError(),modulefilename);exit (EXIT_FAILURE);}
    dm=(int(*)(int*,double*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*))GetProcAddress(handle,"du1calc");
    //*(void **)(&dm)=GetProcAddress(handle,"du1calc");
     if (dm==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function du1calc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    mq=(void(*)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*))GetProcAddress(handle,"mq");
    //*(void **)(&mq)=GetProcAddress(handle,"mq");
     if (mq==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function mq not possible - continuing\n",(int)GetLastError(),modulefilename);}
    estates=(void(*)(ComplexMatrix*,Vector*,double*,double*,Vector*,char**))GetProcAddress(handle,"estates");
    //*(void **)(&estates)=GetProcAddress(handle,"estates");
     if (estates==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function estates not possible - continuing\n",(int)GetLastError(),modulefilename);
                                est=ComplexMatrix(0,2,1,2);// not used, just initialize to prevent errors
                                est=0;
                               }
    mcalc_parameter_storage=(void(*)(ComplexMatrix*,Vector*,double*,double*,Vector*,char**))GetProcAddress(handle,"mcalc_parameter_storage_matrix_init");
    //*(void **)(&mcalc_parameter_storage)=GetProcAddress(handle,"mcalc_parameter_storage_matrix_init");
    if (mcalc_parameter_storage==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %X  module %s loading function mcalc_parameter_storage_matrix_init not possible - continuing\n",(int)GetLastError(),modulefilename);
                                  mcalc_parstorage=ComplexMatrix(0,2,1,2);mcalc_parstorage=0;// not used, just initialize to prevent errors
                                  }

  ddnn=(int(*)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexVector*))GetProcAddress(handle,"dv1calc");
    //*(void **)(&dnn)=GetProcAddress(handle,"dv1calc");
     if (ddnn==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function dv1calc not possible - continuing\n",(int)GetLastError(),modulefilename);}

    sd_m=(void(*)(Vector*,int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"spindensity_mcalc");
    //*(void **)(&sd_m)=GetProcAddress(handle,"spindensity_mcalc");
    if (sd_m==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function spindensity_mcalc not possible - continuing\n",(int)GetLastError(),modulefilename);}

    od_m=(void(*)(Vector*,int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"orbmomdensity_mcalc");
    //*(void **)(&od_m)=GetProcAddress(handle,"orbmomdensity_mcalc");
    if (od_m==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function orbmomdensity_mcalc not possible - continuing\n",(int)GetLastError(),modulefilename);}

    ro_calc=(void(*)(double*,double*,double*,double*,Vector*,double*,Vector*,double*,Vector*,char**))GetProcAddress(handle,"ro_calc");
    //*(void **)(&sd_m)=GetProcAddress(handle,"spindensity_mcalc");
    if (ro_calc==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function ro_calc not possible - continuing\n",(int)GetLastError(),modulefilename);}

#else
  char * error;
  handle=dlopen (modulefilename,RTLD_NOW | RTLD_GLOBAL);
  if (!handle){fprintf (stderr, "jjjpar::jjjpar - Could not load dynamic library\n");
               if ((error=dlerror())!=NULL)
	         {fprintf (stderr,"%s\n",error);}
	       exit (EXIT_FAILURE);
	      }
  //m=(void(*)(Vector*,double*,Vector*,double*,Vector*,char**,double*,double*,ComplexMatrix*))dlsym(handle,"mcalc");
  *(void **)(&m)=dlsym(handle,"mcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s\n",error);exit (EXIT_FAILURE);}
  //dm=(int(*)(int*,double*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*))dlsym(handle,"du1calc");
  *(void **)(&dm)=dlsym(handle,"du1calc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dm=NULL;}
  //mq=(void(*)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*))dlsym(handle,"mq");
  *(void **)(&mq)=dlsym(handle,"mq");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);mq=NULL;}
  //estates=(void(*)(ComplexMatrix*,Vector*,double*,double*,Vector*,char**))dlsym(handle,"estates");
  *(void **)(&estates)=dlsym(handle,"estates");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);estates=NULL;
                                est=ComplexMatrix(0,2,1,2);est=0;// not used, just initialize to prevent errors
                               }
  //mcalc_parameter_storage=(void(*)(ComplexMatrix*,Vector*,double*,double*,Vector*,char**))dlsym(handle,"mcalc_parameter_storage_matrix_init");
  *(void **)(&mcalc_parameter_storage)=dlsym(handle,"mcalc_parameter_storage_matrix_init");

  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);mcalc_parameter_storage=NULL;
                                mcalc_parstorage=ComplexMatrix(0,2,1,2);mcalc_parstorage=0;// not used, just initialize to prevent errors
                               }

  //ddnn=(int(*)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexVector*))dlsym(handle,"dv1calc");
  *(void **)(&ddnn)=dlsym(handle,"dv1calc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);ddnn=NULL;}

  //sd_m=(void(*)(Vector*,int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*))dlsym(handle,"spindensity_mcalc");
  *(void **)(&sd_m)=dlsym(handle,"spindensity_mcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);sd_m=NULL;}

  //od_m=(void(*)(Vector*,int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*))dlsym(handle,"orbmomdensity_mcalc");
  *(void **)(&od_m)=dlsym(handle,"orbmomdensity_mcalc");
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
// function to calculate magnetisation M from effective field H
// this is the heart of the meanfield algorithm an it is necessary to
// keep this routine as efficient as possible
/****************************************************************************/
void jjjpar::mcalc (Vector &mom, double & T, Vector &  gjmbH, double & lnZ,double & U,ComplexMatrix & parstorage)
{switch (module_type)
  {case 1: kramer(mom,T,gjmbH,lnZ,U);break;
   case 2:
   case 4: (*iops).cfieldJJ(mom,T,gjmbH,lnZ,U,parstorage);break;
   case 3: brillouin(mom,T,gjmbH,lnZ,U);break;
   case 5: cluster_mcalc(mom,T,gjmbH,lnZ,U);break;
   default: (*m)(&mom,&T,&gjmbH,&gJ,&ABC,&cffilename,&lnZ,&U,&parstorage);
  }
}

/****************************************************************************/
// this function returns n (the number of transitions in the single ion susceptibility)
// the transition matrix mat first eigenvector u1 corresponding to jjjpar.transitionnumber and delta
// for effective field heff and temperature given on input
/****************************************************************************/
int jjjpar::du1calc(double & T,Vector & gjmbheff,ComplexVector & u1,float & delta,ComplexMatrix & ests)
{ switch (module_type)
  {case 0: if (dm!=NULL){return (*dm)(&transitionnumber,&T,&gjmbheff,&gJ,&ABC,&cffilename,&u1,&delta,&ests);}
           else return 0;
           break;
   case 1: return kramerdm(transitionnumber,T,gjmbheff,u1,delta);break;
   case 2:
   case 4: return (*iops).cfielddm(transitionnumber,T,gjmbheff,u1,delta,ests);break;
   case 3: return brillouindm(transitionnumber,T,gjmbheff,u1,delta);break;
   case 5: return cluster_dm(transitionnumber,T,gjmbheff,u1,delta);break;
   default: return 0;
  }
}


/****************************************************************************/
// returns eigenvalues, boltzmann population and eigenstates matrix parameters of ion
/****************************************************************************/
ComplexMatrix & jjjpar::eigenstates (Vector & gjmbheff,double & T)
{switch (module_type)
  {case 0:  if(estates!=NULL){(*estates)(&est,&gjmbheff,&gJ,&T,&ABC,&cffilename);}
            return est;break;
   case 2:
   case 4: (*iops).cfeigenstates(&est,gjmbheff,T);return est;break;
   default: est=0;return est;
  }
}

/****************************************************************************/
// returns eigenvalues, boltzmann population and eigenstates matrix parameters of ion
/****************************************************************************/
ComplexMatrix & jjjpar::mcalc_parameter_storage_init (Vector & gjmbheff,double & T)
{switch (module_type)
  {case 0:  if(mcalc_parameter_storage!=NULL){(*mcalc_parameter_storage)(&mcalc_parstorage,&gjmbheff,&gJ,&T,&ABC,&cffilename);}
            return mcalc_parstorage;break;
   case 2:
   case 4: (*iops).cfeigenstates(&mcalc_parstorage,gjmbheff,T);return mcalc_parstorage;break;
   default: mcalc_parstorage=0;return mcalc_parstorage;
  }
}
/****************************************************************************/
// returns operator matrices (n=0 Hamiltonian, n=1,...,nofcomponents: operators of moment components)
/****************************************************************************/
Matrix jjjpar::opmat(int n,Vector & gjmbH)
{switch (module_type)
  {case 1:  return krameropmat(n,gjmbH);break;
   default: fprintf(stderr,"ERROR operator calculation in module jjjpar - opmat function not defined for module %i\n",module_type);exit(EXIT_FAILURE);
  }
}


// OBSERVABLES *******************************************************
//1 . NEUTRON SCATTERING OPERATOR  --------------------------------------


/****************************************************************************/
// calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
// according to stored eigenstate matrix est
// input: Qvec ..... Q Vector components 123=xyz=cab
/****************************************************************************/
ComplexVector & jjjpar::MQ(Vector & Qvec)
{double J0,J2,J4,J6;
 double Q,th,ph;
            Q = Norm(Qvec); //dspacing
//            d = 2.0 * PI / Q; s=0.5 / d;
      J0=j0(Q);
      J2=j2(Q);
      J4=j4(Q);
      J6=j6(Q);
            complex<double>dummy;
switch (module_type)
  {case 0:  getpolar(Qvec(2),Qvec(3),Qvec(1),Q,th,ph); // for external module we must provide th and ph with respect
                                                       // to abc coordinate system
            (*mq)(&Mq,&th,&ph,&J0,&J2,&J4,&J6,&est);
             // external module provide Mq(123)=Mq(abc)
             // we must transform this to mcdiff internal xyz||cab coordinate system
            dummy=Mq(3);Mq(3)=Mq(2);Mq(2)=Mq(1);Mq(1)=dummy;
            return Mq;break;
   case 2:  getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph); // internal module cfield does not need transformation
            return (*iops).MQ(th,ph,J0,J2,J4,J6,Zc,est);break;
   case 4:  getpolar(Qvec(2),Qvec(3),Qvec(1),Q,th,ph); // for so1ion we muyst th and ph with respect to abc coordinate system
            Mq=(*iops).MQ(th,ph,J0,J2,J4,J6,Zc,est);
             // so1ion module provide Mq(123)=Mq(abc)
             // we must transform this to mcdiff internal xyz||cab coordinate system
            dummy=Mq(3);Mq(3)=Mq(2);Mq(2)=Mq(1);Mq(1)=dummy;
            return Mq;break;
   default: fprintf(stderr,"ERROR in scattering operator function M(Q) for ion %s \nM(Q) is currently only implemented for internal module cfield and so1ion:\n",cffilename);exit(EXIT_FAILURE);
  }
}

/****************************************************************************/
// returns transition element matrix N(Q) first eigenvector v1 in order to be able to go beyond
//
// dipolar approximation in mcdisp - it requires a call to eigenstates first
//
//on input
//    transitionnumber has to be set correctly to that one which is to be computed
//    sign(transitionnumber)... 1... without printout, -1 with extensive printout
//    est		matrix with eigenstates, eigenvalues [meV], population numbers
//    T                 temperature
//     Q                 components of Qvector in euclidian coordinates 123=abc
//  on output
//    int   	total number of transitions
//    N(i,j)	<-|Q|+><+|Q|-> (n+-n-),  n+,n- population numbers
//    with Q the scattering operator according to Lovesey 11.4, p 222, eq 6.87b
//     (note that  <M(Q)>=-2x<Q>_TH in units of mb)
//    .... occupation number of states (- to + transition chosen according to transitionnumber)
//
/****************************************************************************/
int jjjpar::dv1calc(Vector & Qvec,double & T, ComplexVector & v1,ComplexMatrix & ests)

{double J0,J2,J4,J6;
 double Q,th,ph;
 int i;     complex<double>dummy; // introduced 3.4.10 MR
            Q = Norm(Qvec); //dspacing
      //      d = 2.0 * PI / Q; s=0.5 / d;
      J0=j0(Q);
      J2=j2(Q);
      J4=j4(Q);
      J6=j6(Q);
	 // calculate th and ph (polar angles of Q with respect to xyz of CEF)
 switch (module_type)
  {static int washere=0;

   case 0:if (ddnn!=NULL){getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph);
                          return (*ddnn)(&transitionnumber,&th,&ph,&J0,&J2,&J4,&J6,&ests,&T,&v1);break;}
          else {return 0;}
   case 2:  getpolar(Qvec(3),Qvec(1),Qvec(2),Q,th,ph); // for internal module cfield xyz||cba and we have to give cfielddn polar angles with respect to xyz
            i=(*iops).cfielddn(transitionnumber,th,ph,J0,J2,J4,J6,Zc,ests,T,v1);
            // and we have to switch indices in matrix nat(1..3,1..3) to conform with xyz||cba changed MR 3.4.10
            dummy=v1(1);v1(1)=v1(2);v1(2)=v1(3);v1(3)=dummy; // changed MR 3.4.10
            return i;break;
   case 4:  getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph); // for internal module so1ion xyz||abc and we have to give cfielddn polar angles with respect to xyz
            return (*iops).cfielddn(transitionnumber,th,ph,J0,J2,J4,J6,Zc,ests,T,v1);break;
   default: if(washere==0){fprintf(stderr,"Warning in scattering operator function dv1calc - for ion %s \ngoing beyond dipolar approximation is not implemented\n",cffilename);
                           washere=1;}
            return 0;
  }

}



/************************************************************************************/
//  RETURN TOTAL FORMFACTOR,
//    however if gJ=0 and Q>0 return spin form factor FS(Q)=<j0(Q)>
//            if gJ=0 and Q<0 return angular  form factor FL(Q)=<j0(Q)>+<j2(Q)>
//  D = 2 * pi / Q
//  s = 1 / 2 / D: sintheta = lambda * s
/************************************************************************************/
   double jjjpar::F(double Q)
// {if(gJ==0&&Q>0){return j0(Q);} // in case of intermediate coupling return spin form factor
//  if(gJ==0&&Q<0){return j0(-Q)+j2(-Q);} // in case of intermediate coupling return angular form factor
// return (j0(Q) + j2(Q) * (2 / gJ - 1)); // formfactor F(Q) for rare earth
   // Rewrote to use saved value if Q same as previous call.
   {
     for(int iq=1; iq<(Qsaved(MAXSAVEQ)==-1e16?nsaved:6); iq++) if(fabs(Q-Qsaved(iq))<1e-6) return Fsaved(iq);
     double Fval;
     if(gJ==0&&Q>0) Fval = j0(Q);
     else if(gJ==0&&Q<0) Fval = j0(-Q)+j2(-Q);
     else Fval = (j0(Q) + j2(Q) * (2 / gJ - 1));
     nsaved++; if(nsaved>MAXSAVEQ) nsaved=1;
     Qsaved(nsaved) = Q; Fsaved(nsaved) = Fval;
     return Fval;
   }
   double jjjpar::j0(double Q)
  {double value=0,s; if(fabs(Q)<0.1)return 1.0;
   if(Np(1)!=0){value=jl(0,Q);}// here enter calculation from radial wave function parameters
   else
   {s=Q/4/PI;    value=magFFj0(1)*exp(-magFFj0(2)*s*s)+magFFj0(3)*exp(-magFFj0(4)*s*s)+magFFj0(5)*exp(-magFFj0(6)*s*s)+magFFj0(7);
   }return value;
  }
   double jjjpar::j1(double Q)
  {double value=0;   if(fabs(Q)<0.1)return 0.0;
   if(Np(1)!=0){value=jl(1,Q);}// here enter calculation from radial wave function parameters
   return value;
  }
   double jjjpar::j2(double Q)
  {double value=0,s;  if(fabs(Q)<0.1)return 0.0;
   s=Q/4/PI;
   if(Np(1)!=0){value=jl(2,Q);}// here enter calculation from radial wave function parameters
    else
   { value=magFFj2(1)*exp(-magFFj2(2)*s*s)+magFFj2(3)*exp(-magFFj2(4)*s*s)+magFFj2(5)*exp(-magFFj2(6)*s*s)+magFFj2(7);
    value*=s*s;
   }return value;
  }
   double jjjpar::j3(double Q)
  {double value=0; if(fabs(Q)<0.1)return 0.0;
   if(Np(1)!=0){value=jl(3,Q);}// here enter calculation from radial wave function parameters
   return value;
  }
   double jjjpar::j4(double Q)
  {double value=0,s;  if(fabs(Q)<0.1)return 0.0;
     s=Q/4/PI;
         if(Np(1)!=0){value=jl(4,Q);}// here enter calculation from radial wave function parameters
     else
   {  value=magFFj4(1)*exp(-magFFj4(2)*s*s)+magFFj4(3)*exp(-magFFj4(4)*s*s)+magFFj4(5)*exp(-magFFj4(6)*s*s)+magFFj4(7);
      value*=s*s;
   }return value;
  }
   double jjjpar::j5(double Q)
  {double value=0; if(fabs(Q)<0.1)return 0.0;
      if(Np(1)!=0){value=jl(5,Q);}// here enter calculation from radial wave function parameters
      return value;
  }
   double jjjpar::j6(double Q)
  {double value=0,s;   if(fabs(Q)<0.1)return 0.0;
     s=Q/4/PI;
         if(Np(1)!=0){value=jl(6,Q);}// here enter calculation from radial wave function parameters
     else
   {  value=magFFj6(1)*exp(-magFFj6(2)*s*s)+magFFj6(3)*exp(-magFFj6(4)*s*s)+magFFj6(5)*exp(-magFFj6(6)*s*s)+magFFj6(7);
      value*=s*s;
    } return value;
  }

   double jjjpar::jl(int l,double QA){
    int p,q, pmax=0;
    double Q=QA*0.5292;// convert Q from 1/A into 1/a0
    Vector coeff(1,9);
    for(p=1;p<=9;++p){if(Np(p)!=0){pmax=p;
                                   coeff(p)=Cp(p)*pow(2.0*Xip(p)/Q,Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"Warning: calculation of <j%i(Q)> failed due to Xi%i<=0 - continuing with <j%i(Q)>=0\n",l,p,l);return 0;}
                     }            }
    if(pmax==0){fprintf (stderr,"Warning: calculation of <j%i(Q)> failed - continuing with <j%i(Q)>=0\n",l,l);return 0;}

    double value=0;
    for(p=1;p<=pmax;++p){
    for(q=1;q<=pmax;++q){
    value+=coeff(p)*coeff(q)*tl(l,(int)Np(p)+(int)Np(q),(Xip(p)+Xip(q))/Q);
                        }}
     return value;
   }



   double jjjpar::tl(int l,int N,double x)
     {double value=0;
      switch (l)
       { case 0: value=sn(1,N,x);break;
         case 1: value=sn(2,N,x)-cn(1,N,x);break;
         case 2: value=3*sn(3,N,x)-3*cn(2,N,x)-sn(1,N,x);break;
         case 3: value=cn(1,N,x)-15*cn(3,N,x)-6*sn(2,N,x)+15*sn(4,N,x);break;
         case 4: value=10*cn(2,N,x)-105*cn(4,N,x)+sn(1,N,x)-45*sn(3,N,x)+105*sn(5,N,x);break;
         case 5: value=-cn(1,N,x)+105*cn(3,N,x)-945*cn(5,N,x)+15*sn(2,N,x)-420*sn(4,N,x)+945*sn(6,N,x);break;
         case 6: value=-21*cn(2,N,x)+1260*cn(4,N,x)-10395*cn(6,N,x)-sn(1,N,x)+210*sn(3,N,x)-4725*sn(5,N,x)+10395*sn(7,N,x);break;
        default: fprintf(stderr,"Error function jjjpar:tl - value l=%i not implemented\n",l);exit(1);
       }
     return value;
     }
/*
   double jjjpar::sn(int n,int N,double x)
   {complex <double> c(x,-1.0);
    double value;
    value=(double)factorial(N-n)*imag(pow(c,-N+n-1));
    return value;
   }
   double jjjpar::cn(int n,int N,double x)
   {complex <double> c(x,-1.0);
    double value;
    value=(double)factorial(N-n)*real(pow(c,-N+n-1));
    return value;
   }
*/
// ------------------------------------------------------------------------- //
// Rewrite ::sn() and ::cn() to avoid using complex numbers, to run faster
// ------------------------------------------------------------------------- //
/* Complex powers, from Mathematica: z=x+I y; Do[Print[ComplexExpand[z^ex]], {ex, 1, 8}] (type I as <ESC>ii<Esc>)
 * x+I y
 * x^2+2 I x y-y^2
 * x^3-3 x y^2+I (3 x^2 y-y^3)
 * x^4-6 x^2 y^2+y^4+I (4 x^3 y-4 x y^3)
 * x^5-10 x^3 y^2+5 x y^4+I (5 x^4 y-10 x^2 y^3+y^5)
 * x^6-15 x^4 y^2+15 x^2 y^4-y^6+I (6 x^5 y-20 x^3 y^3+6 x y^5)
 * x^7-21 x^5 y^2+35 x^3 y^4-7 x y^6+I (-7 x^6 y+35 x^4 y^3-21 x^2 y^5+y^7)
 * x^8-28 x^6 y^2+70 x^4 y^4-28 x^2 y^6+y^8+I (-8 x^7 y+56 x^5 y^3-56 x^3 y^5+8 x y^7)
 * when y=-1:                          For powers z^{-p}, take -ve of Im part and div by (1-x^2)^p.
 * x                           +I( -1 )
 * -1+x^2                      +I( -2 x )
 *  -3 x+x^3                   +I(  1-3 x^2 )
 *  1-6 x^2+x^4                +I(  4 x-4 x^3 )
 *  5 x-10 x^3+x^5             +I( -1+10 x^2-5 x^4 )
 *  -1+15 x^2-15 x^4+x^6       +I( -6 x+20 x^3-6 x^5 )
 * -7 x+35 x^3-21 x^5+x^7      +I( 1-21 x^2+35 x^4-7 x^6)
 * 1-28 x^2+70 x^4-28 x^6+x^8  +I( 8 x-56 x^3+56 x^5-8 x^7)
 */
 double jjjpar::sn(int n,int N,double x)    // Need imaginary part
 {
    double denom=1.; if((-N+n-1)<0) denom=pow(1+x*x,-(-N+n-1));
    switch(-N+n-1) {
      case  0: return 0.; break;
      case  1: return (double)factorial(N-n) *  -1.; break;
      case -1: return (double)factorial(N-n) * (1. / denom); break;
      case  2: return (double)factorial(N-n) *  -2*x; break;
      case -2: return (double)factorial(N-n) * ( 2*x / denom); break;
      case  3: return (double)factorial(N-n) *   (1-3*x*x); break;
      case -3: return (double)factorial(N-n) * (-(1-3*x*x) / denom); break;
      case  4: return (double)factorial(N-n) *   4*x * (1-x*x); break;
      case -4: return (double)factorial(N-n) * (-4*x * (1-x*x) / denom); break;
      case  5: return (double)factorial(N-n) *   (-1 + x*x * (10 - 5*x*x)); break;
      case -5: return (double)factorial(N-n) * (-(-1 + x*x * (10 - 5*x*x)) / denom); break;
      case  6: return (double)factorial(N-n) *   x * (-6 + x*x * (20 - 6*x*x)); break;
      case -6: return (double)factorial(N-n) * (-x * (-6 + x*x * (20 - 6*x*x)) / denom);  break;
      case  7: return (double)factorial(N-n) *   ( 1 + x*x * (-21 + x*x * (35 - 6*x*x))); break;
      case -7: return (double)factorial(N-n) * (-( 1 + x*x * (-21 + x*x * (35 - 6*x*x))) / denom); break;
      case  8: return (double)factorial(N-n) *   x * (8 + x*x * (-56 + x*x * (56 - 8*x*x))); break;
      case -8: return (double)factorial(N-n) * (-x * (8 + x*x * (-56 + x*x * (56 - 8*x*x))) / denom); break;
      default: //fprintf(stderr,"jjjpar::sn() Bad power %i\n",-N+n-1); exit(-1);
         complex <double> c(x,-1.0); return factorial((double)(N-n))*imag(pow(c,-N+n-1.));
    }
 }
 double jjjpar::cn(int n,int N,double x)    // Need real part
 {
    double denom=1.; if((-N+n-1)<0) denom=pow(1+x*x,-(-N+n-1));
    switch(-N+n-1) {
      case  0: return 1.; break;
      case  1: return (double)factorial(N-n) *  x; break;
      case -1: return (double)factorial(N-n) * (x / denom); break;
      case  2: return (double)factorial(N-n) *  (-1+x*x); break;
      case -2: return (double)factorial(N-n) * ((-1+x*x) / denom); break;
      case  3: return (double)factorial(N-n) *  x * (-3+x*x); break;
      case -3: return (double)factorial(N-n) * (x * (-3+x*x) / denom); break;
      case  4: return (double)factorial(N-n) *  (1 + x*x * (-6+x*x)); break;
      case -4: return (double)factorial(N-n) * ((1 + x*x * (-6+x*x)) / denom); break;
      case  5: return (double)factorial(N-n) *  x * (5 + x*x *(-10+x*x)); break;
      case -5: return (double)factorial(N-n) * (x * (5 + x*x *(-10+x*x)) / denom); break;
      case  6: return (double)factorial(N-n) *  (-1 + x*x * (15 + x*x * (-15+x*x))); break;
      case -6: return (double)factorial(N-n) * ((-1 + x*x * (15 + x*x * (-15+x*x))) / denom); break;
      case  7: return (double)factorial(N-n) *  x * (-7 + x*x *(35 + x*x * (-21 + x*x))); break;
      case -7: return (double)factorial(N-n) * (x * (-7 + x*x *(35 + x*x * (-21 + x*x))) / denom); break;
      case  8: return (double)factorial(N-n) *  ( 1 + x*x * (-28 + x*x * (70 + x*x * (-28+x*x)))); break;
      case -8: return (double)factorial(N-n) * (( 1 + x*x * (-28 + x*x * (70 + x*x * (-28+x*x)))) / denom); break;
      default: //fprintf(stderr,"jjjpar::cn() Bad power %i\n",-N+n-1); exit(-1);
         complex <double> c(x,-1.0); return factorial((double)(N-n))*imag(pow(c,-N+n-1.));
    }
 }

/************************************************************************************/
//   debyewallerfactor = exp(-2 * DWF *s*s)      (sf ~ exp(-2 DWF sin^2(theta) / lambda^2)=EXP (-W),  (2*DWF=B=8 pi^2 <u^2>)
/************************************************************************************/
   double jjjpar::debyewallerfactor(double & Q)
   {double s;
    s=Q/4/PI;
    return exp(-2*DWF*s*s);
   }

// 2. charge density ----------------------------------------------------------

/************************************************************************************/
// evaluate radial wave function
/************************************************************************************/
double jjjpar::radial_wavefunction(double rr) // rr given in Angstroems, returns R(r) in units of 1/A^1.5
   {//printf("%g ",rr);

    double R=0;int p;double a0=0.5292;
    int ok=0;
    double r=rr/a0;// r is the distance in units of a0
    for(p=1;p<=9;++p){if(Np(p)!=0){ok=1;
                                   R+=exp(-Xip(p)*r)*pow(r,Np(p)-1)*Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"\n\nWarning: calculation of radial wave function R(r=%g) failed due to Xi%i<=0 - continuing with R(r=%g)=0\n\n\n",r,p,r);return 0;}
                     }            }
    // now we have R in units of 1/a0^1.5
    R/=sqrt(a0*a0*a0);
    // now we have R in units of 1/A^1.5
//printf("%g ",R);
    if (ok==1) return R;


//  we have to find the 4f wavefunction R4f(r) for each single ion and the Zlm, cfield has nothing: so we have
//     to take this from chrgplt.bas - a little problem: how do we get the correct R4f(r) ? for a first attempt
//     we could just take the same for all RE.

    static int washere=0;
    if(washere==0){washere=1;fprintf (stderr,"\n\n!! Warning !!: radial wave function parameters not found, will use 4f hydrogen radial wave function\n\n\n");}
double rs;
//k^2 = 11 / 10 * 11 / 9 * 11 / 8 * 11 / 7 * 11 / 6 * 11 / 5 * 11 / 4 * 11 / 3 * 11 / 2 * 11 / 1 * 11
rs = rr * exp(-rr);
R = 280.4 * rs * rs * rs * rs  * exp(-1.5 * rr);
//printf("R4f(%g)=%g\n ",rr,R);
return R;
   }
/************************************************************************************/
   //functions to calculate radial matrix elements <r^n> from radial wave function in units of a0=0.5292 A
/************************************************************************************/
   double jjjpar::rk_from_radial_wavefunction(int k)
   {int p,q, pmax=0;
    Vector coeff(1,9);

    for(p=1;p<=9;++p){if(Np(p)!=0){pmax=p;
                                   coeff(p)=Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"Warning: calculation of <r^%i> failed due to Xi%i<=0 - continuing with <r^%i>=0\n",k,p,k);return 0;}
                     }            }
    if(pmax==0){fprintf (stderr,"Warning: calculation of <r^%i> failed - continuing with <r^%i>=0\n",k,k);return 0;}
    double rk=0;
    for(p=1;p<=pmax;++p){
    for(q=1;q<=pmax;++q){// the following does not work because of large factorials and powers
                        // rk+=coeff(p)*coeff(q)*factorial((int)Np(p)+(int)Np(q)+k)/pow(Xip(p)+Xip(q),Np(p)+Np(q)+k+1);
                        // therefore we have to substitute it by  calculating the product "by hand"
                        double product=1.0;
                        int nnk=(int)Np(p)+(int)Np(q)+k;
                        int i;
                        for(i=1;i<=nnk;++i){product*=i/(Xip(p)+Xip(q));}
                       rk+=coeff(p)*coeff(q)*product/(Xip(p)+Xip(q));
                        }}
     return rk;
   }

   int jjjpar::r2_from_radial_wavefunction() {r2=rk_from_radial_wavefunction(2);if(module_type==2||module_type==4){(*iops).r2=r2;}return true;}
   int jjjpar::r4_from_radial_wavefunction() {r4=rk_from_radial_wavefunction(4);if(module_type==2||module_type==4){(*iops).r4=r4;}return true;}
   int jjjpar::r6_from_radial_wavefunction() {r6=rk_from_radial_wavefunction(6);if(module_type==2||module_type==4){(*iops).r6=r6;}return true;}

void jjjpar::save_radial_wavefunction(const char * filename)
   {double r=0.1;
    FILE * fout;
    if (radial_wavefunction(r)==0){fprintf(stderr,"Warning: save_radial_wavefunction not possible\n");return;}
    fout=fopen_errchk(filename,"w");
    fprintf(fout,"# radial wave function for %s\n",cffilename);
    fprintf(fout,"# the radial wave function is expanded as \n");
    fprintf(fout,"# R(r)=sum_p C_p R_Np,XIp(r)\n");
    fprintf(fout,"# R_Np,XIp(r)=r^(Np-1).exp(-XIp * r).(2 * XIp)^(Np+0.5)/sqrt(2Np!)\n");
    fprintf(fout,"# radial wave function parameters Np XIp Cp values are\n");
    fprintf(fout,"# tabulated in clementi & roetti Atomic data and \n");
    fprintf(fout,"# nuclear data tables 14 (1974) 177-478 for the transition metals\n");
    fprintf(fout,"# for rare earth parameters can be found in Freeman and Watson PR 127 (1962) 2058\n");
    fprintf(fout,"# and O. Sovers, J. Phys. Chem. Solids Vol 28 (1966) 1073\n");
    fprintf(fout,"# the parameters used are (a0=0.5292 A): \n");
    int p;
    for(p=1;p<=9;++p){if(Np(p)!=0){fprintf(fout,"#! N%i=%g XI%i=%g /a0 C%i=%g\n",p,Np(p),p,Xip(p),p,Cp(p));}}
    fprintf(fout,"# r[A]  vs R(r)[1/A^1.5]\n");
    for(r=0.01;r<=10;r*=1.05){fprintf(fout,"%8.8g  %8.8g\n",r,radial_wavefunction(r));}
    fclose(fout);
   }


//***********************************************************************
// sub for calculation of charge density given a radiu R and polar angles teta,
// fi and expansion coefficients alm
//***********************************************************************
double jjjpar::rocalc (double & teta,double & fi,double & R, Vector & moments,double & T, Vector &  gjmbH)
{double ro,rr;
if((module_type==0)&&(ro_calc!=NULL)){(*ro_calc)(&ro,&teta,&fi,&R,&moments,&T,&gjmbH,&gJ,&ABC,&cffilename);return ro;}
if (R>4.0||R<0){ro=0;}else{
 int l,m;
  Vector tetan(1,6);
 int offset=0; // for cfield module
 if(module_type==0){offset=3;} // for external modules this is to treat ic1ion correctly
 if(module_type==2||module_type==4){
                    tetan(2)=(*iops).alpha;// Stevens factors
                    tetan(4)=(*iops).beta;
                    tetan(6)=(*iops).gamma;
// if(calcmagdensity>0){fprintf(stderr,"Problem: calcmagdensity>0 in %s, magnetisation densities in module so1ion and cfield do not work correctly yet, quitting... \n",cffilename);
//     exit(EXIT_FAILURE);}  // here I quit because it is yet unclear if the formulas programmed in are correct
                           // I assume that there is proportionality between
                           //         sum_i(2si+li) Zlm(Omega_i) to
                           // and
                           //               1/2 (J Olm(J) + Olm(J) J)   (see ionpars.cpp, function cfield() )
                           // with coefficients being gJ*tetan(l)
                           // ... this is probably not correct: 1. Wigner eckhardt holds only for spherical tensor operators and not for products of such
                           //                                   2. even if it holds, then the coefficients have to be calculated for spin and orbital contributions and added
                           // todo: check Wigner Eckhardt theorem and if it holds calculate coefficients
                           //
                                    }

Matrix a(0,6,-6,6);
 if(nof_electrons==0){fprintf(stderr,"Error: nof_electrons=0 ... perhaps single ion property file %s does not contain the number of electrons in the shell: 'nof_electrons=...'\n",cffilename);
     exit(EXIT_FAILURE);}
 a(0, 0) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
//if(calcmagdensity>0)a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);

a(2,-2)=moments(offset+4);
a(2,-1)=moments(offset+5);
a(2,0)=moments(offset+6);
a(2,1)=moments(offset+7);
a(2,2)=moments(offset+8);

a(4,-4)=moments(offset+16);
a(4,-3)=moments(offset+17);
a(4,-2)=moments(offset+18);
a(4,-1)=moments(offset+19);
a(4, 0)=moments(offset+20);
a(4, 1)=moments(offset+21);
a(4, 2)=moments(offset+22);
a(4, 3)=moments(offset+23);
a(4, 4)=moments(offset+24);

a(6,-6)=moments(offset+36);
a(6,-5)=moments(offset+37);
a(6,-4)=moments(offset+38);
a(6,-3)=moments(offset+39);
a(6,-2)=moments(offset+40);
a(6,-1)=moments(offset+41);
a(6,-0)=moments(offset+42);
a(6, 1)=moments(offset+43);
a(6, 2)=moments(offset+44);
a(6, 3)=moments(offset+45);
a(6, 4)=moments(offset+46);
a(6, 5)=moments(offset+47);
a(6, 6)=moments(offset+48);


// r given in Angstroems, returns R(r) in units of 1/A^1.5
rr=radial_wavefunction(R);
rr=rr*rr;// then the chargedensity will be in units of 1/A^3
for(l=1;l<=5;l+=2){for(m=-l;m<=l;++m){a(l,m)=0;}}
if(module_type==2||module_type==4){for(l=2;l<=6;l+=2){for(m=-l;m<=l;++m){a(l,m)*=tetan(l)*cnst(l,m);}}
         } // these are prefactors in case of module cfield and so1ion(stevens parameters tetan and zlm prefactors)
else     {for(l=2;l<=6;l+=2){for(m=-l;m<=l;++m){if(m!=0){a(l,m)*=sqrt((2.0*l+1)/8/PI);}else{a(l,m)*=sqrt((2.0*l+1)/4/PI);}}}
         } // in case
           // of module ic1ion we just take the prefactors of the Zlm ... ??? what should we take here ???

 ro=-rr*zlmsum(a,teta,fi); // minus, because electrons are negative
 }
return ro;
}


// sum over different Zlm using the coefficients a(l,m)
double jjjpar::zlmsum(Matrix & a, double & teta, double & fi)
{double ro,ct,ct2,st,st2,sfi,cfi,cf2,sf2;
 ct = cos(teta);                      // z/r
 ct2 = ct * ct;
 st = sin(teta);   // y/r=st sfi
 st2 = st * st;
 sfi = sin(fi);    // x/r=st cfi
 sf2=sfi*sfi;
 cfi = cos(fi);
 cf2=cfi*cfi;
 ro = a(0, 0) * cnst(0,0);

 ro = ro + a(1,-1) * cnst(1,-1) * st *sfi;
 ro = ro + a(1,0) * cnst(1,0) * ct;
 ro = ro + a(1,1) * cnst(1,1) * st *cfi;

 ro = ro + a(2, -2)* cnst(2,-2)  * 2 * st2 * sfi * cfi;
 ro = ro + a(2, -1)* cnst(2,-1)  * st * sfi * ct;
 ro = ro + a(2, 0)* cnst(2,0)  * (3 * ct2 - 1);
 ro = ro + a(2, 1)* cnst(2,1)  * st * cfi * ct;
 ro = ro + a(2, 2)* cnst(2,2)  * st2 * (cfi * cfi - sfi * sfi);

 ro = ro + a(3,-3) * cnst(3,-3) * (3* st2*cf2 - st2*sf2)*st*sfi ;
 ro = ro + a(3,-2) * cnst(3,-2) * 2*sfi*cfi*st2*ct;
 ro = ro + a(3,-1) * cnst(3,-1) * st*sfi*(5*ct2-1);
 ro = ro + a(3,0) * cnst(3,0) * ct *(5*ct2-3);
 ro = ro + a(3,1) * cnst(3,1) * st*cfi*(5*ct2-1);
 ro = ro + a(3,2) * cnst(3,2) * (st2*cf2-st2*sf2)*ct;
 ro = ro + a(3,3) * cnst(3,3) * st*cfi*(st2*cf2-3*st2*sf2);


 ro = ro + a(4, -4)* cnst(4,-4) * st2 * st2 * 4 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
 ro = ro + a(4, -3)* cnst(4,-3) * ct * st * st2 * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
 ro = ro + a(4, -2)* cnst(4,-2) * (7 * ct2 - 1) * 2 * st2 * cfi * sfi;
 ro = ro + a(4, -1)* cnst(4,-1) * st * sfi * ct * (7 * ct2 - 3);
 ro = ro + a(4, 0)* cnst(4,0) * (35 * ct2 * ct2 - 30 * ct2 + 3);
 ro = ro + a(4, 1) * cnst(4,1) * st * cfi * ct * (7 * ct2 - 3);
 ro = ro + a(4, 2)* cnst(4,2)  * (7 * ct2 - 1) * st2 * (cfi * cfi - sfi * sfi);
 ro = ro + a(4, 3)* cnst(4,3)  * ct * st * st2 * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
 ro = ro + a(4, 4) * cnst(4,4) * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);

// x/r=st cfi     y/r=st sfi    ct= z/r

 ro = ro + a(5,-5) * cnst(5,-5) *st*sfi*(5*st2*st2*cf2*cf2-10*st2*cf2*st2*sf2+st2*st2*sf2*sf2);
ro = ro + a(5,-4) * cnst(5,-4) *4*st*sfi*(st2*cf2-st2*sf2)*ct;
ro = ro + a(5,-3) * cnst(5,-3) *st*sfi*(3*st2*cf2-st2*sf2)*(9*ct2-1);
ro = ro + a(5,-2) * cnst(5,-2) *2*st*cfi*st*sfi*ct*(3*ct2-1);
ro = ro + a(5,-1) * cnst(5,-1) *st*sfi*(21*ct2*ct2-14*ct2+1);
ro = ro + a(5,0) * cnst(5,0) *ct*(63*ct2*ct2-70*ct2+15);
ro = ro + a(5,1) * cnst(5,1) *st*cfi*(21*ct2*ct2-14*ct2+1);
ro = ro + a(5,2) * cnst(5,2) *(st2*cf2-st2*sf2)*ct*(3*ct2-1);
ro = ro + a(5,3) * cnst(5,3) *st*cfi*(st2*cf2-3*st2*sf2)*(9*ct2-1);
ro = ro + a(5,4) * cnst(5,4) *(st2*st2*cf2*cf2-6*st2*cf2*st2*sf2+st2*st2*sf2*sf2)*ct;
ro = ro + a(5,5) * cnst(5,5) *st*cfi*(st2*st2*cf2*cf2-10*st2*cf2*st2*sf2+5*st2*st2*sf2*sf2);


 ro = ro + a(6, -6)* cnst(6,-6) * st2 * st2 * st2 * (6 * cfi * cfi * cfi * cfi * cfi * sfi - 20 * cfi * cfi * cfi * sfi * sfi * sfi + 6 * cfi * sfi * sfi * sfi * sfi * sfi);
 ro = ro + a(6, -5)* cnst(6,-5) * ct * st * st2 * st2 * (5 * cfi * cfi * cfi * cfi * sfi - 10 * cfi * cfi * sfi * sfi * sfi + sfi * sfi * sfi * sfi * sfi);
 ro = ro + a(6, -4)* cnst(6,-4) * (11 * ct2 - 1) * 4 * st2 * st2 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
 ro = ro + a(6, -3) * cnst(6,-3)* (11 * ct * ct2 - 3 * ct) * st2 * st * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
 ro = ro + a(6, -2)* cnst(6,-2) * 2 * st2 * sfi * cfi * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2);
 ro = ro + a(6, -1)* cnst(6,-1) * ct * st * sfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
 ro = ro + a(6, 0)* cnst(6,0) * (231 * ct2 * ct2 * ct2 - 315 * ct2 * ct2 + 105 * ct2 - 5);
 ro = ro + a(6, 1)* cnst(6,1)  * ct * st * cfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
 ro = ro + a(6, 2)* cnst(6,2)  * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2) * st2 * (cfi * cfi - sfi * sfi);
 ro = ro + a(6, 3)* cnst(6,3)  * (11 * ct * ct2 - 3 * ct) * st2 * st * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
 ro = ro + a(6, 4)* cnst(6,4)  * (11 * ct2 - 1) * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);
 ro = ro + a(6, 5)* cnst(6,5)  * ct * st * st2 * st2 * (cfi * cfi * cfi * cfi * cfi - 10 * cfi * cfi * cfi * sfi * sfi + 5 * cfi * sfi * sfi * sfi * sfi);
 ro = ro + a(6, 6)* cnst(6,6) * st2 * st2 * st2 * (cfi * cfi * cfi * cfi * cfi * cfi - 15 * cfi * cfi * cfi * cfi * sfi * sfi + 15 * cfi * cfi * sfi * sfi * sfi * sfi - sfi * sfi * sfi * sfi * sfi * sfi);
 return ro;
 }

void jjjpar::set_zlm_constants()
{// cnst is the Zlm constants - put them into the matrix ... (same code is reused in ionpars.cpp)
 cnst= Matrix(0,6,-6,6);

cnst(0,0)=0.28209479177387814347403972578039;

cnst(1,0)=0.48860251190291992158638462283835;
cnst(1,1)=0.48860251190291992158638462283835;

cnst(2,0) = 0.3153962;
cnst(2,1)=  1.092548;
cnst(2,2)=  0.5462823;

cnst(3,0)=0.373176332590115391414395913199;
cnst(3,1)=0.45704579946446573615802069691665;
cnst(3,2)=1.445305721320277027694690077199;
cnst(3,3)=0.5900435899266435103456102775415;


cnst(4,0)=  0.1057871;
cnst(4,1)=  0.6690465;
cnst(4,2)=  0.4730943;
cnst(4,3)=  1.77013;
cnst(4,4)=  0.625845;

cnst(5,0)=0.11695032245342359643971519209028;
cnst(5,1)=0.45294665119569692129844165821715;
cnst(5,2)=2.3967683924866618775009505697816;
cnst(5,3)=0.48923829943525038768400871584296;
cnst(5,4)=2.0756623148810412789957985225952;
cnst(5,5)=0.65638205684017010281411876637614;


cnst(6,0)=  0.06357014;
cnst(6,1)=  0.582621;
cnst(6,2)=  0.4606094;
cnst(6,3)=  0.921205;
cnst(6,4)=  0.5045723;
cnst(6,5)=  2.366619;
cnst(6,6)=  0.6831942;
int l,m;
for(l=1;l<=6;l+=1){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);}
}

// 3. moment density ----------------------------------------------------------

/****************************************************************************/
// function to calculate coefficients of expansion of spindensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
/****************************************************************************/
void jjjpar::spindensity_mcalc (Vector &mom,int xyz, double & T, Vector &  gjmbH, ComplexMatrix & parstorage)
{switch (module_type)
  {case 1: fprintf(stderr,"Problem: magnetisation densities in module kramer are not possible, quitting... \n");
           exit(EXIT_FAILURE);break;
   case 2:
   case 4: fprintf(stderr,"Problem: magnetisation densities in module so1ion and cfield do not work correctly yet, quitting... \n");
           exit(EXIT_FAILURE);break;
// comment on module so1ion/cfield:
//fprintf(stderr,"Problem: calcmagdensity>0 in %s, magnetisation densities in module so1ion and cfield do not work correctly yet, quitting... \n",cffilename);
//     exit(EXIT_FAILURE);}  // here I quit because it is yet unclear if the formulas programmed in are correct
                           // I assume that there is proportionality between
                           //         sum_i(2si+li) Zlm(Omega_i) to
                           // and
                           //               1/2 (J Olm(J) + Olm(J) J)   (see ionpars.cpp, function cfield() )
                           // with coefficients being gJ*tetan(l)
                           // ... this is probably not correct: 1. Wigner eckhardt holds only for spherical tensor operators and not for products of such
                           //                                   2. even if it holds, then the coefficients have to be calculated for spin and orbital contributions and added
                           // todo: check Wigner Eckhardt theorem and if it holds calculate coefficients
                           //
   case 3: fprintf(stderr,"Problem: magnetisation densities in module brillouin are not possible, quitting... \n");
           exit(EXIT_FAILURE);break;
   case 0: (*sd_m)(&mom,&xyz,&T,&gjmbH,&gJ,&ABC,&cffilename,&parstorage);break;
   case 5:fprintf(stderr,"Problem: magnetisation densities are not possible in module cluster, quitting... \n");
           exit(EXIT_FAILURE);break;
   default:fprintf(stderr,"Problem: magnetisation densities are not possible in module, quitting... \n");
           exit(EXIT_FAILURE);break;
  }
}

/****************************************************************************/
// function to calculate coefficients of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
/****************************************************************************/
void jjjpar::orbmomdensity_mcalc (Vector &mom,int  xyz, double & T, Vector &  gjmbH, ComplexMatrix & parstorage)
{switch (module_type)
  {case 1: fprintf(stderr,"Problem: magnetisation densities in module kramer are not possible, quitting... \n");
           exit(EXIT_FAILURE);break;
   case 2:
   case 4: fprintf(stderr,"Problem: magnetisation densities in module so1ion and cfield do not work correctly yet, quitting... \n");
           exit(EXIT_FAILURE);break;
// comment on module so1ion/cfield:
//fprintf(stderr,"Problem: calcmagdensity>0 in %s, magnetisation densities in module so1ion and cfield do not work correctly yet, quitting... \n",cffilename);
//     exit(EXIT_FAILURE);}  // here I quit because it is yet unclear if the formulas programmed in are correct
                           // I assume that there is proportionality between
                           //         sum_i(2si+li) Zlm(Omega_i) to
                           // and
                           //               1/2 (J Olm(J) + Olm(J) J)   (see ionpars.cpp, function cfield() )
                           // with coefficients being gJ*tetan(l)
                           // ... this is probably not correct: 1. Wigner eckhardt holds only for spherical tensor operators and not for products of such
                           //                                   2. even if it holds, then the coefficients have to be calculated for spin and orbital contributions and added
                           // todo: check Wigner Eckhardt theorem and if it holds calculate coefficients
                           //
   case 3: fprintf(stderr,"Problem: magnetisation densities in module brillouin are not possible, quitting... \n");
           exit(EXIT_FAILURE);break;
   case 0: (*od_m)(&mom,&xyz,&T,&gjmbH,&gJ,&ABC,&cffilename,&parstorage);break;
   case 5:fprintf(stderr,"Problem: magnetisation densities are not possible in module cluster, quitting... \n");
           exit(EXIT_FAILURE);break;
   default:fprintf(stderr,"Problem: magnetisation densities are not possible in module, quitting... \n");
           exit(EXIT_FAILURE);break;
  }
}

//***********************************************************************
// sub for calculation of spin density component given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
double jjjpar::spindensity_calc (double & teta,double & fi,double & R, Vector & moments)
{double ro,rr;

if (R>4.0||R<0){ro=0;}else{

 Matrix a(0,6,-6,6); 
//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849 
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
// R given in Angstroems, returns R(r) in units of 1/A^1.5
   rr=radial_wavefunction(R);
   rr=rr*rr;// then the spindensity will be in units of 1/A^3
int i;
 if(moments.Hi()==49)
 {
 for (i=1;i<=49;++i){a(k[i],q[i])=moments(i);}
   ro=rr*zlmsum(a,teta,fi);
 }
 else
 {fprintf(stderr,"Error jjjpar.spindensitycalc: dimension of moments=%i must be 49\n",moments.Hi());
  exit(EXIT_FAILURE);
 }
 }
return ro;
}

//***********************************************************************
// sub for calculation of spin density vector given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
Vector  jjjpar::spindensity_calc (double & teta,double & fi,double & R, Vector & momentsx,Vector & momentsy,Vector & momentsz)
{double rr;
 static Vector mm(1,3);
if (R>4.0||R<0){mm=0;}else{

 Matrix a(0,6,-6,6);
//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
// R given in Angstroems, returns R(r) in units of 1/A^1.5
   rr=radial_wavefunction(R);
   rr=rr*rr;// then the spindensity will be in units of 1/A^3
int i;
  for (i=1;i<=49;++i){a(k[i],q[i])=momentsx(i);}
   mm(1)=rr*zlmsum(a,teta,fi);
  for (i=1;i<=49;++i){a(k[i],q[i])=momentsy(i);}
   mm(2)=rr*zlmsum(a,teta,fi);
  for (i=1;i<=49;++i){a(k[i],q[i])=momentsz(i);}
   mm(3)=rr*zlmsum(a,teta,fi);
 }
return mm;
}

//***********************************************************************
// subs for calculation gradient of spin  density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
 Matrix jjjpar::gradspindensity_calc(double & teta,double & fi,double & R, Vector & momentx, Vector & momenty, Vector & momentz)
{static Matrix grad(1,3,1,3);
 if (R>3.9||R<0){grad=0;}else{
 Vector m0(1,3);Vector m1(1,3);Vector m2(1,3);Vector m3(1,3);
 double d=0.01; // differential in Angstroem
 double teta1,teta2,teta3,fi1,fi2,fi3,R1,R2,R3;
 double ct,st,sf,cf;
 ct = cos(teta); st = sin(teta);   // y/r=st sfi
 sf = sin(fi);   cf = cos(fi);
 // now we use Jacobi Matrix (dr,dth,dfi)=J (dx,dy,dz)
  // dx                      dy             dz
 R1=R+d*st*cf;          R2=R+d*st*sf;       R3=R+d*ct;
 teta1=teta+d*ct*cf/R;  teta2=teta+d*ct*sf/R; teta3=teta-d*st/R;
 if(st>0){fi1=fi-d*sf/st/R;fi2=fi+d*cf/st/R;}else{fi1=fi;fi2=fi;} fi3=fi;

 m0=spindensity_calc(teta,fi,R,momentx,momenty,momentz);
 m1=spindensity_calc(teta1,fi1,R1,momentx,momenty,momentz);
 m2=spindensity_calc(teta2,fi2,R2,momentx,momenty,momentz);
 m3=spindensity_calc(teta3,fi3,R3,momentx,momenty,momentz);
 //            d/dx     d/dy      d/dz
 m1=(m1-m0)/d; m2=(m2-m0)/d; m3=(m3-m0)/d;grad=MatrixfromVectors(m1,m2,m3);
 }
 return grad;
}
//***********************************************************************
// subs for calculation gradient of orbital moment density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
 Matrix jjjpar::gradorbmomdensity_calc(double & teta,double & fi,double & R, Vector & momentx, Vector & momenty, Vector & momentz)
{static Matrix grad(1,3,1,3);
 if (R>3.9||R<0){grad=0;}else{
 Vector m0(1,3);Vector m1(1,3);Vector m2(1,3);Vector m3(1,3);
 double d=0.01; // differential in Angstroem
 double teta1,teta2,teta3,fi1,fi2,fi3,R1,R2,R3;
 double ct,st,sf,cf;
 ct = cos(teta);  st = sin(teta);   // y/r=st sfi
 sf = sin(fi);    cf = cos(fi);
 // now we use Jacobi Matrix (dr,dth,dfi)=J (dx,dy,dz)
  // dx                      dy             dz
 R1=R+d*st*cf;          R2=R+d*st*sf;       R3=R+d*ct;
 teta1=teta+d*ct*cf/R;  teta2=teta+d*ct*sf/R; teta3=teta-d*st/R;
 if(st>0){fi1=fi-d*sf/st/R;fi2=fi+d*cf/st/R;}else{fi1=fi;fi2=fi;} fi3=fi;

 m0=orbmomdensity_calc(teta,fi,R,momentx,momenty,momentz);
 m1=orbmomdensity_calc(teta1,fi1,R1,momentx,momenty,momentz);
 m2=orbmomdensity_calc(teta2,fi2,R2,momentx,momenty,momentz);
 m3=orbmomdensity_calc(teta3,fi3,R3,momentx,momenty,momentz);
 //            d/dx     d/dy      d/dz
 m1=(m1-m0)/d; m2=(m2-m0)/d; m3=(m3-m0)/d;grad=MatrixfromVectors(m1,m2,m3);
 }
 return grad;
}

//***********************************************************************
// subs for calculation gradient of current density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
 Matrix jjjpar::gradcurrdensity_calc(double & teta,double & fi,double & R, Vector & momentx, Vector & momenty, Vector & momentz)
{static Matrix grad(1,3,1,3);
 if (R>3.9||R<0){grad=0;}else{
 Vector m0(1,3);Vector m1(1,3);Vector m2(1,3);Vector m3(1,3);
 double d=0.01; // differential in Angstroem
 double teta1,teta2,teta3,fi1,fi2,fi3,R1,R2,R3;
 double ct,st,sf,cf;
 ct = cos(teta);  st = sin(teta);   // y/r=st sfi
 sf = sin(fi);    cf = cos(fi);
 // now we use Jacobi Matrix (dr,dth,dfi)=J (dx,dy,dz)
  // dx                      dy             dz
 R1=R+d*st*cf;          R2=R+d*st*sf;       R3=R+d*ct;
 teta1=teta+d*ct*cf/R;  teta2=teta+d*ct*sf/R; teta3=teta-d*st/R;
 if(st>0){fi1=fi-d*sf/st/R;fi2=fi+d*cf/st/R;}else{fi1=fi;fi2=fi;} fi3=fi;

 m0=currdensity_calc(teta,fi,R,momentx,momenty,momentz);
 m1=currdensity_calc(teta1,fi1,R1,momentx,momenty,momentz);
 m2=currdensity_calc(teta2,fi2,R2,momentx,momenty,momentz);
 m3=currdensity_calc(teta3,fi3,R3,momentx,momenty,momentz);
 //            d/dx     d/dy      d/dz
 m1=(m1-m0)/d; m2=(m2-m0)/d; m3=(m3-m0)/d;grad=MatrixfromVectors(m1,m2,m3);
 }
 return grad;
}

//***********************************************************************
// sub for calculation of orbital moment density component given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
double jjjpar::orbmomdensity_calc (double & teta,double & fi,double & R, Vector & moments)
{double ro,rr;

if (R>4.0||R<0){ro=0;}else{

Matrix a(0,6,-6,6);
//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

int i;
for (i=1;i<=49;++i){a(k[i],q[i])=moments(i);}

// R given in Angstroems, returns R(r) in units of 1/A^1.5
  rr=Fr(R);
  // then the orbital moment density will be in units of 1/A^3
  ro=rr*zlmsum(a,teta,fi);
 }
return ro;
}

//***********************************************************************
// sub for calculation of orbital moment density vector given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
Vector jjjpar::orbmomdensity_calc (double & teta,double & fi,double & R, Vector & momentsx, Vector & momentsy, Vector & momentsz)
{double rr;
 static Vector mm(1,3);
if (R>4.0||R<0){mm=0;}else{

Matrix a(0,6,-6,6);
//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

int i;
// R given in Angstroems, returns R(r) in units of 1/A^1.5
  rr=Fr(R);
  // then the orbital moment density will be in units of 1/A^3

 for (i=1;i<=49;++i){a(k[i],q[i])=momentsx(i);}
  mm(1)=rr*zlmsum(a,teta,fi);
 for (i=1;i<=49;++i){a(k[i],q[i])=momentsy(i);}
  mm(2)=rr*zlmsum(a,teta,fi);
 for (i=1;i<=49;++i){a(k[i],q[i])=momentsz(i);}
  mm(3)=rr*zlmsum(a,teta,fi);
 }
return mm;
}

//***********************************************************************
// sub for calculation of orbital current density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
Vector jjjpar::currdensity_calc (double & teta,double & fi,double & R, Vector & momentlx, Vector & momently, Vector & momentlz)
{double ro,rr;
 static Vector mm(1,3);
if (R>4.0||R<0){mm=0;}else{
 
 Matrix ax(0,6,-6,6),ay(0,6,-6,6),az(0,6,-6,6);
 Matrix bx(0,6,-6,6),by(0,6,-6,6),bz(0,6,-6,6);
 Matrix dx(0,6,-6,6),dy(0,6,-6,6),dz(0,6,-6,6);
 double ct,st,sf,cf,fp;
 ct = cos(teta);                      // z/r
 st = sin(teta);   // y/r=st sfi
 sf = sin(fi);    // x/r=st cfi
 cf = cos(fi);
 Vector Jr(1,3),Jth(1,3),Jfi(1,3); // Jacobi matrix
 Jr(1)=st*cf; Jr(2)=st*sf; Jr(3)=ct;
 Jth(1)=ct*cf;Jth(2)=ct*sf;Jth(3)=-st;
 if(st>0){Jfi(1)=-sf/st;Jfi(2)=cf/st;Jfi(3)=0;}else{Jfi=0;}
 Vector res(1,3),a(1,3);

//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

int i;
  // R given in Angstroems, returns R(r) in units of 1/A^1.5
   ro=radial_wavefunction(R);
   ro=0.9274e-1*ro*ro/R;// then the spindensity will be in units of 1/A^3  ro=R(R);
   rr=0.9274e-1*Fr(R)/R;
  // then the orbital moment density will be in units of 1/A^3
  // then orbital current density will be in mb/A^4
  // transform to milliAmpere/A^2 by multiplying with 0.9274e-4

for (i=1;i<=49;++i){
ax(k[i],q[i])=momentlx(i);
ay(k[i],q[i])=momently(i);
az(k[i],q[i])=momentlz(i);
}

for (i=1;i<=49;++i){
a(1)=momentlx(i);
a(2)=momently(i);
a(3)=momentlz(i);
xproduct(res,a,Jr);
bx(k[i],q[i])=res(1);
by(k[i],q[i])=res(2);
bz(k[i],q[i])=res(3);
// now blm coefficients are calculated

// we need now the dlm
dx(k[i],q[i])=res(1);
dy(k[i],q[i])=res(2);
dz(k[i],q[i])=res(3);

     // add second term m* Jfi x al-m
if(q[i]!=0){
     a(1)=q[i]*ax(k[i],-q[i]);
     a(2)=q[i]*ay(k[i],-q[i]);
     a(3)=q[i]*az(k[i],-q[i]);
     xproduct(res,Jfi,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
     // add third term cottheta |m| Jth x alm
     if(st>0){
     a(1)=fabs(q[i])*ct/st*ax(k[i],q[i]);
     a(2)=fabs(q[i])*ct/st*ay(k[i],q[i]);
     a(3)=fabs(q[i])*ct/st*az(k[i],q[i]);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
             }
           }
// insert here addition of flm term !!
if(q[i]==-1){
     a(1)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*sf*ax(k[i],0);
     a(2)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*sf*ay(k[i],0);
     a(3)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*sf*az(k[i],0);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
            }
if(q[i]==+1){
     a(1)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*cf*ax(k[i],0);
     a(2)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*cf*ay(k[i],0);
     a(3)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*cf*az(k[i],0);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
            }
if(q[i]>+1){
     fp=-sqrt(k[i]*(double)(k[i]+1)-q[i]*(q[i]-1));
     a(1)=fp*cf*ax(k[i],q[i]-1)-fp*sf*ax(k[i],-q[i]+1);
     a(2)=fp*cf*ay(k[i],q[i]-1)-fp*sf*ay(k[i],-q[i]+1);
     a(3)=fp*cf*az(k[i],q[i]-1)-fp*sf*az(k[i],-q[i]+1);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
            }
if(q[i]<-1){
     fp=-sqrt(k[i]*(double)(k[i]+1)-q[i]*(q[i]+1));
     a(1)=fp*cf*ax(k[i],q[i]+1)+fp*sf*ax(k[i],-q[i]-1);
     a(2)=fp*cf*ay(k[i],q[i]+1)+fp*sf*ay(k[i],-q[i]-1);
     a(3)=fp*cf*az(k[i],q[i]+1)+fp*sf*az(k[i],-q[i]-1);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
            }

}

  mm(1)=ro*zlmsum(bx,teta,fi)+rr*zlmsum(dx,teta,fi);
  mm(2)=ro*zlmsum(by,teta,fi)+rr*zlmsum(dy,teta,fi);
  mm(3)=ro*zlmsum(bz,teta,fi)+rr*zlmsum(dz,teta,fi);

 }
return mm;
}


/************************************************************************************/
// evaluate F(r) for orbital momentum density (see mynotes, balcar 1975)
/************************************************************************************/
double jjjpar::Fr(double rr) // evaluate F(r)=1/r integral_r^inf dx R^2(x)
                      // r in units of Angstroems, F(r) in units of 1/A^3
   {//printf("%g ",rr);

    double F_R=0;int p,pp,i;double a0=0.5292;
    double fp,fpp,sumi;
    int ok=0;
    double r=rr/a0;// r is the distance in units of a0
    for(p=1;p<=9;++p){if(Np(p)!=0){ok=1;
                                   if(Xip(p)<=0){fprintf (stderr,"\n\nWarning: calculation of radial integral F(r=%g) failed due to Xi%i<=0 - continuing with F(r=%g)=0\n\n\n",r,p,r);return 0;}
    fp=exp(-Xip(p)*r)*pow(r,Np(p)-1)*Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));

    for(pp=1;pp<=9;++pp){if(Np(pp)!=0){
     fpp=exp(-Xip(pp)*r)*pow(r,Np(pp)-1)*Cp(pp)*pow(2.0*Xip(pp),Np(pp)+0.5)/sqrt((double)factorial(2*(int)Np(pp)));
     sumi=0;for(i=1;i<=Np(p)+Np(pp)-2;++i){sumi+=pow(r,-i)/factorial((int)Np(p)+(int)Np(pp)-2-i)/pow(Xip(p)+Xip(pp),i+1);}
     sumi*=factorial(Np(p)+Np(pp)-2)/r;
                   F_R+=fp*fpp*sumi;
                          }             }
                     }            }
    // now we have F_R in units of 1/a0^3
    F_R/=a0*a0*a0;
    // now we have F_R in units of 1/A^3
    //printf("%g ",F_R);
    if (ok==1) return F_R;


//  we have to find the 4f wavefunction R4f(r) for each single ion and the Zlm, cfield has nothing: so we have
//     to take this from chrgplt.bas - a little problem: how do we get the correct R4f(r) ? for a first attempt
//     we could just take the same for all RE.

    static int washere=0;
    if(washere==0){washere=1;fprintf (stderr,"\n\n!! Warning !!: radial wave function parameters not found, will use 4f hydrogen radial wave function\n\n\n");}
double rs;
//k^2 = 11 / 10 * 11 / 9 * 11 / 8 * 11 / 7 * 11 / 6 * 11 / 5 * 11 / 4 * 11 / 3 * 11 / 2 * 11 / 1 * 11
rs = rr * exp(-rr);
fp = 280.4 * rs * rs * rs * rs  * exp(-1.5 * rr);

sumi=0;for(i=1;i<=8;++i){sumi+=pow(r,-i)/factorial(8-i)/pow(11.0,i+1);}
sumi*=factorial(8)/r;
F_R+=fp*fp*sumi;

    // now we have F_R in units of 1/a0^3
    F_R/=a0*a0*a0;

return F_R;
   }



