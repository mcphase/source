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
void jjjpar::get_parameters_from_sipfile(char * sipf_filename)
{FILE * cf_file;
 int i,j;
 float nn[MAXNOFNUMBERSINLINE];
 nn[0]=MAXNOFNUMBERSINLINE;
 modulefilename=new char[MAXNOFCHARINLINE];

 char instr[MAXNOFCHARINLINE];
  cf_file = fopen_errchk (sipf_filename, "rb");
  fgets_errchk (instr, MAXNOFCHARINLINE, cf_file);
  if(extract(instr,"MODULE=",modulefilename,(size_t)MAXNOFCHARINLINE))
   {if(extract(instr,"#!",modulefilename,(size_t)MAXNOFCHARINLINE))
    {fprintf(stderr,"Error: single ion property file %s does not start with '#!' or 'MODULE='\n",sipf_filename);
     exit(EXIT_FAILURE);}
   }
   //ic1ion entered without path ?
      if (strncmp(modulefilename,"ic1ion",6)==0)
      {strcpy(modulefilename,getenv("MCPHASE_DIR"));strcat(modulefilename,"/bin/ic1ion_module/ic1ion.so");}
      if (strncmp(modulefilename,"icf1ion",7)==0)
      {strcpy(modulefilename,getenv("MCPHASE_DIR")); strcat(modulefilename,"/bin/ic1ion_module/icf1ion.so"); }
      if (strncmp(modulefilename,"phonon",6)==0)
      {strcpy(modulefilename,getenv("MCPHASE_DIR")); strcat(modulefilename,"/bin/phonon_module/phonon.so"); }

  fprintf (stderr,"#parsing single ion property file: %s - loading module %s\n",sipf_filename,modulefilename);

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
      if(i!=0){fprintf(stderr,"Error reading |<+-|Ja|-+>|,|<+-|Jb|-+>|,|<+-|Jc|+->| from file %s\ncorrect file format is:\n",sipf_filename);
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
      if(i!=0){fprintf(stderr,"Error reading spin quantum number J=S from file %s\ncorrect file format is:\n",sipf_filename);
              fprintf(stderr,"\n#!brillouin\n#comment lines ..\n# Quantum number  J\nJ=3.5\n\n");exit(EXIT_FAILURE);}
      // now we have the numbers corresponding to the vector ABC() in nn[]
      fprintf(stderr," ... Brillouin function with J=S=%g\n",ABC(1));
      est=ComplexMatrix(0,2,1,2);est=0;// not used, just initialize to prevent errors
      Icalc_parstorage=ComplexMatrix(0,2,1,2);Icalc_parstorage=0;// not used, just initialize to prevent errors
     }
     else
     {if(strcmp(modulefilename,"cfield")==0)
     {module_type=2;fprintf (stderr,"#[internal]\n");
      //fclose(cf_file);cf_file = fopen_errchk (sipf_filename, "rb"); // reopen file
       fseek(cf_file,0,SEEK_SET);
      iops=new ionpars(cf_file,sipf_filename);

      int dj;dj=(int)(2*J()+1);
      est=ComplexMatrix(0,dj,1,dj);
      Icalc_parstorage=ComplexMatrix(0,dj,1,dj);
      nof_electrons=(*iops).nof_electrons;
      // get 1ion parameters - operator matrices

     }
     else
     {if(strcmp(modulefilename,"so1ion")==0)
     {module_type=4;fprintf (stderr,"#[internal]\n");
     // fclose(cf_file);cf_file = fopen_errchk (sipf_filename, "rb"); // reopen file
      fseek(cf_file,0,SEEK_SET);
      iops=new ionpars(cf_file,sipf_filename);
      nof_electrons=(*iops).nof_electrons;
      int dj;dj=(int)(2*J()+1);
      est=ComplexMatrix(0,dj,1,dj);
      Icalc_parstorage=ComplexMatrix(0,dj,1,dj);
      // get 1ion parameters - operator matrices

     }
     else if (strcmp(modulefilename,"cluster")==0)
     {module_type=5;fprintf (stderr,"#[internal]\n");
      ABC=Vector(1,1);i=1;
      nof_electrons=0; // not to be used in module cluster !!
      while(feof(cf_file)==false)
      {fgets(instr, MAXNOFCHARINLINE, cf_file);
       i+=extract(instr,"structurefile",clusterfilename,MAXNOFCHARINLINE)-1;
      }// input all  lines starting with comments
      if(i!=0){fprintf(stderr,"Error reading structurefile from file %s\ncorrect file format is:\n",sipf_filename);
              fprintf(stderr,"\n#!MODULE=cluster\n#comment lines ..\n# next line contains cluster structure filename\n"
                            "#!structurefile=cluster.j\n"
                            "\n");exit(EXIT_FAILURE);}
      fprintf(stderr," ... reading cluster structure from %s\n",clusterfilename);
      clusterpars =new par(clusterfilename);
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
    du=(int(*)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*))GetProcAddress(handle,"du1calc");
    //*(void **)(&du)=GetProcAddress(handle,"du1calc");
     if (du==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function du1calc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    
    p=(void(*)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"pcalc");
    //*(int **)(&p)=GetProcAddress(handle,"pcalc");
     if (p==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function pcalc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    dP1=(int(*)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*))GetProcAddress(handle,"dP1");
    //*(void **)(&du)=GetProcAddress(handle,"dP1");
     if (dP1==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function dP1 not possible - continuing\n",(int)GetLastError(),modulefilename);}

    m=(void(*)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"mcalc");
    //*(int **)(&m)=GetProcAddress(handle,"mcalc");
     if (m==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function mcalc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    dm1=(int(*)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*))GetProcAddress(handle,"dm1");
    //*(void **)(&dm1)=GetProcAddress(handle,"dm1");
     if (dm1==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function dm1 not possible - continuing\n",(int)GetLastError(),modulefilename);}

    L=(void(*)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"Lcalc");
    //*(int **)(&L)=GetProcAddress(handle,"Lcalc");
     if (L==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function Lcalc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    dL1=(int(*)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*))GetProcAddress(handle,"dL1");
    //*(void **)(&dL1)=GetProcAddress(handle,"dL1");
     if (dL1==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function dL1 not possible - continuing\n",(int)GetLastError(),modulefilename);}
    S=(void(*)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"Scalc");
    //*(int **)(&S)=GetProcAddress(handle,"Scalc");
     if (S==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function Scalc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    dS1=(int(*)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*))GetProcAddress(handle,"dS1");
    //*(void **)(&dS1)=GetProcAddress(handle,"dS1");
     if (dS1==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function dS1 not possible - continuing\n",(int)GetLastError(),modulefilename);}
    mq=(void(*)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*))GetProcAddress(handle,"mqcalc");
    //*(void **)(&mq)=GetProcAddress(handle,"mqcalc");
     if (mq==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function mqcalc not possible - continuing\n",(int)GetLastError(),modulefilename);}
    ddnn=(int(*)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexVector*,double*))GetProcAddress(handle,"dmq1");
    //*(void **)(&dnn)=GetProcAddress(handle,"dmq1");
     if (ddnn==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function dmq1 not possible - continuing\n",(int)GetLastError(),modulefilename);}
    rixs=(int(*)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexVector*,double*))GetProcAddress(handle,"drixs1");
    //*(void **)(&dnn)=GetProcAddress(handle,"rixs1");
     if (rixs==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function drixs1 not possible - continuing\n",(int)GetLastError(),modulefilename);}

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


    cd_m=(void(*)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"chargedensity_coeff");
    //*(void **)(&cd_m)=GetProcAddress(handle,"chargedensity_coeff");
    if (cd_m==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function chargedensity_coeff not possible - continuing\n",(int)GetLastError(),modulefilename);}
    cd_dm=(int(*)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*))GetProcAddress(handle,"dchargedensity_coeff1");
    //*(void **)(&cd_dm)=GetProcAddress(handle,"dchargedensity_coeff1");
     if (cd_dm==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function dchargedensity_coeff1 not possible - continuing\n",(int)GetLastError(),modulefilename);}


    sd_m=(void(*)(Vector*,int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"spindensity_coeff");
    //*(void **)(&sd_m)=GetProcAddress(handle,"spindensity_coeff");
    if (sd_m==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function spindensity_coeff not possible - continuing\n",(int)GetLastError(),modulefilename);}
    sd_dm=(int(*)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,int*,float*,ComplexMatrix*))GetProcAddress(handle,"dspindensity_coeff1");
    //*(void **)(&sd_dm)=GetProcAddress(handle,"dspindensity_coeff1");
     if (sd_dm==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function dspindensity_coeff1 not possible - continuing\n",(int)GetLastError(),modulefilename);}

    od_m=(void(*)(Vector*,int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*))GetProcAddress(handle,"orbmomdensity_coeff");
    //*(void **)(&od_m)=GetProcAddress(handle,"orbmomdensity_coeff");
    if (od_m==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function orbmomdensity_coeff not possible - continuing\n",(int)GetLastError(),modulefilename);}
    od_dm=(int(*)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,int*,float*,ComplexMatrix*))GetProcAddress(handle,"dorbmomdensity_coeff1");
    //*(void **)(&od_dm)=GetProcAddress(handle,"dorbmomdensity_coeff1");
     if (od_dm==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function dorbmomdensity_coeff1 not possible - continuing\n",(int)GetLastError(),modulefilename);}

    ro_calc=(void(*)(double*,double*,double*,double*,Vector*,double*,Vector*,char**))GetProcAddress(handle,"ro_calc");
    //*(void **)(&ro_calc)=GetProcAddress(handle,"spindensity_coeff");
    if (ro_calc==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function ro_calc not possible - continuing\n",(int)GetLastError(),modulefilename);}

    dyn_opmat=(int(*)(int*,char**,Vector*,Vector*,Matrix*))GetProcAddress(handle,"opmat");
    if (dyn_opmat==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function opmat not possible - continuing\n",(int)GetLastError(),modulefilename);}
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
  *(void **)(&du)=dlsym(handle,"du1calc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);du=NULL;}

 *(void **)(&p)=dlsym(handle,"pcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);p=NULL;}
 *(void **)(&dP1)=dlsym(handle,"dP1");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dP1=NULL;}


 *(void **)(&m)=dlsym(handle,"mcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);m=NULL;}
 *(void **)(&dm1)=dlsym(handle,"dm1");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dm1=NULL;}
 *(void **)(&L)=dlsym(handle,"Lcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);L=NULL;}
 *(void **)(&dL1)=dlsym(handle,"dL1");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dL1=NULL;}
 *(void **)(&S)=dlsym(handle,"Scalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);S=NULL;}
 *(void **)(&dS1)=dlsym(handle,"dS1");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dS1=NULL;}
 
 *(void **)(&mq)=dlsym(handle,"mqcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);mq=NULL;}
  *(void **)(&ddnn)=dlsym(handle,"dmq1");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);ddnn=NULL;}
  *(void **)(&rixs)=dlsym(handle,"drixs1");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);rixs=NULL;}

  *(void **)(&estates)=dlsym(handle,"estates");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);estates=NULL;
                                est=ComplexMatrix(0,2,1,2);est=0;// not used, just initialize to prevent errors
                               }
 *(void **)(&Icalc_parameter_storage)=dlsym(handle,"Icalc_parameter_storage_matrix_init");

  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);Icalc_parameter_storage=NULL;
                                Icalc_parstorage=ComplexMatrix(0,2,1,2);Icalc_parstorage=0;// not used, just initialize to prevent errors
                               }


  *(void **)(&cd_m)=dlsym(handle,"chargedensity_coeff");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);cd_m=NULL;}
  *(void **)(&cd_dm)=dlsym(handle,"dchargedensity_coeff1");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);cd_dm=NULL;}

  *(void **)(&sd_m)=dlsym(handle,"spindensity_coeff");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);sd_m=NULL;}
  *(void **)(&sd_dm)=dlsym(handle,"dspindensity_coeff1");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);sd_dm=NULL;}

  *(void **)(&od_m)=dlsym(handle,"orbmomdensity_coeff");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);od_m=NULL;}
  *(void **)(&od_dm)=dlsym(handle,"dorbmomdensity_coeff1");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);od_dm=NULL;}

  *(void **)(&ro_calc)=dlsym(handle,"ro_calc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);ro_calc=NULL;}

  *(void **)(&dyn_opmat)=dlsym(handle,"opmat");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dyn_opmat=NULL;}

#endif
     }
    }
   }
  }
 // fclose(cf_file);
 fseek(cf_file,0,SEEK_SET);

  magFFj0=Vector(1,MAGFF_NOF_COEFF);magFFj0=0;  magFFj0[1]=1;
  magFFj2=Vector(1,MAGFF_NOF_COEFF);magFFj2=0;
  magFFj4=Vector(1,MAGFF_NOF_COEFF);magFFj4=0;
  magFFj6=Vector(1,MAGFF_NOF_COEFF);magFFj6=0;
  Zc=Vector(1,7);Zc=0;

   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;

  DWF=0;  gJ=0;maxE=1e10;pinit=0;ninit=1e10;

 // cf_file = fopen_errchk (sipf_filename, "rb");
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
    extract(instr,"FFj0d",magFFj0[8]);
    extract(instr,"FFj0E",magFFj0[9]);
    extract(instr,"FFj2A",magFFj2[1]);
    extract(instr,"FFj2a",magFFj2[2]);
    extract(instr,"FFj2B",magFFj2[3]);
    extract(instr,"FFj2b",magFFj2[4]);
    extract(instr,"FFj2C",magFFj2[5]);
    extract(instr,"FFj2c",magFFj2[6]);
    extract(instr,"FFj2D",magFFj2[7]);
    extract(instr,"FFj2d",magFFj2[8]);
    extract(instr,"FFj2E",magFFj2[9]);
    extract(instr,"FFj4A",magFFj4[1]);
    extract(instr,"FFj4a",magFFj4[2]);
    extract(instr,"FFj4B",magFFj4[3]);
    extract(instr,"FFj4b",magFFj4[4]);
    extract(instr,"FFj4C",magFFj4[5]);
    extract(instr,"FFj4c",magFFj4[6]);
    extract(instr,"FFj4D",magFFj4[7]);
    extract(instr,"FFj4d",magFFj4[8]);
    extract(instr,"FFj4E",magFFj4[9]);
    extract(instr,"FFj6A",magFFj6[1]);
    extract(instr,"FFj6a",magFFj6[2]);
    extract(instr,"FFj6B",magFFj6[3]);
    extract(instr,"FFj6b",magFFj6[4]);
    extract(instr,"FFj6C",magFFj6[5]);
    extract(instr,"FFj6c",magFFj6[6]);
    extract(instr,"FFj6D",magFFj6[7]);
    extract(instr,"FFj6d",magFFj6[8]);
    extract(instr,"FFj6E",magFFj6[9]);
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
 if(module_type==5)cluster_ini_Imat();
// check gJ
if(module_type==2&&fabs(gJ-(*iops).gJ)>0.00001)
{fprintf(stderr,"Error internal module cfield : Lande Factor read from %s (gJ=%g) does not conform to internal module value gJ=%g\n",sipf_filename,gJ,(*iops).gJ);exit(EXIT_FAILURE);}
if(module_type==4&&fabs(gJ-(*iops).gJ)>0.00001)
{fprintf(stderr,"Error internal module so1ion : Lande Factor read from %s (gJ=%g) does not conform to internal module value gJ=%g\n",sipf_filename,gJ,(*iops).gJ);exit(EXIT_FAILURE);}
//if (gJ==0){printf("# reading gJ=0 in single ion property file %s -> entering intermediate coupling mode by assigning Ja=Sa Jb=La Jc=Sb Jd=Lb Je=Sc Jf=Lc (S... Spin, L... angular momentum)\n",sipf_filename);
//           if (module_type==1){fprintf(stderr,"Error internal module kramers: intermediate coupling not supported\n");exit(EXIT_FAILURE);}
//           if (module_type==2){fprintf(stderr,"Error internal module cfield : intermediate coupling not supported\n");exit(EXIT_FAILURE);}
//           if (module_type==3){fprintf(stderr,"Error internal module brillouin: intermediate coupling not supported\n");exit(EXIT_FAILURE);}
//           if (module_type==4){fprintf(stderr,"Error internal module so1ion : intermediate coupling not supported\n");exit(EXIT_FAILURE);}
//          }

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
// from exchange field Hxc [meV] and external field Hext
// this is the heart of the meanfield algorithm an it is necessary to
// keep this routine as efficient as possible
/****************************************************************************/
void jjjpar::Icalc (Vector &mom, double & T, Vector &  Hxc,Vector & Hext ,double & lnZ,double & U,ComplexMatrix & parstorage)
{switch (module_type)
  {case 1: kramer(mom,T,Hxc,Hext,lnZ,U);break;
   case 2:
   case 4: (*iops).Icalc(mom,T,Hxc,Hext,lnZ,U,parstorage);break;
   case 3: brillouin(mom,T,Hxc,Hext,lnZ,U);break;
   case 5: cluster_Icalc_mcalc_Micalc (1,mom,T,Hxc,Hext,lnZ,U);break;
   default: (*I)(&mom,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&lnZ,&U,&parstorage);
  }
}

/****************************************************************************/
// this function returns n (the number of transitions in the single ion susceptibility)
// the transition matrix mat first eigenvector u1 corresponding to jjjpar.transitionnumber and delta
// for effective field heff and temperature given on input
/****************************************************************************/
int jjjpar::du1calc(double & T,Vector &  Hxc,Vector & Hext,ComplexVector & u1,float & delta,int & n, int & nd, ComplexMatrix & ests)
{delta=maxE;u1(1)=complex <double> (ninit,pinit);
  switch (module_type)
  {case 0: if (du!=NULL){return (*du)(&transitionnumber,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&u1,&delta,&n,&nd,&ests);}
           else return 0;
           break;
   case 1: return kramerdm(transitionnumber,T,Hxc,Hext,u1,delta,n,nd);break;
   case 2:
   case 4: return (*iops).du1calc(transitionnumber,T,Hxc,Hext,u1,delta,n,nd,ests);break;
   case 3: return brillouindm(transitionnumber,T,Hxc,Hext,u1,delta,n,nd);break;
   case 5: return cluster_dm(1,transitionnumber,T,u1,delta,n,nd,ests);break;
   default: return 0;
  }
}

/****************************************************************************/
// this function calculates series of single ion susceptibility matrices for 
// different energies
   // output:returns 0 on success
   //        the Matrices chi0pointer[1....nofstps] must exist and will be filled with values
   //        ...... the contribution of transition transitionnumber is added to these matrices
   // input: emin est nofstps define energies, eps is the imaginary part of the energy
   //        Q       the Q vector in 1/A
   //        qcounter is a counter telling which q vector in the list is calculated
   //                  sign(qcounter) <0 indicates that chi0c matrices should be cleared
   //        delta ... sign determines if energy gain or loss term is added
/****************************************************************************/
int jjjpar:: chi0(ComplexMatrix ** chi0pointer,double & emin, double  estp, int & nofstps, double & epsilon, Vector & Q, 
                  int qcounter,float & delta,double & T,Vector &  Hxc,Vector & Hext, ComplexMatrix & ests,
                   int i1,int j1,int k1,int l1)
{ // for the moment do nothing module specific but use existing module function to calculate internal
  // well defined chi0
  // ... in future we may then do something more clever by putting here values from a file which is created by external
  // programs such as bfk ... this is triggered by epsilon <0
 
 if(fabs(qcounter)<2)// only do something for first q vector (all others will have the same chi0 [currently not q dependence in chi0]
 {if(qcounter<0){
  if(epsilon>0){for(int i=0;i<nofstps;++i)(*chi0pointer[i])=0; // clear matrices
  }else{// load externally chi0 from bfk0.res type of file
   if(module_type!=4||Hxc.Hi()!=3){fprintf(stderr,"Error mcdisp -r <0 cannot load external chi0: not module so1ion or mf dimension !=3\n");exit(EXIT_FAILURE);}
   printf("running singleion and bfk and loading chi0 from bfk0.res\n");
   // 1. output levels.cef
   char command[MAXNOFCHARINLINE],instr[MAXNOFCHARINLINE];
   sprintf(command,"singleion -r %s %g %g %g %g  %g %g %g",sipffilename,T,Hext(1),Hext(2),Hext(3),Hxc(1),Hxc(2),Hxc(3));
   system(command);
   // 2. create bfk.par
   FILE *file;
   file=fopen_errchk("./results/bfk.par","w");
   fprintf(file,"# Parameter file  bfk.par\n"
                "#\n"
                "#!emin=%g\n"
                "#!emax=%g\n"
                "#!Npoints=%i\n"
                "#!E=50\n"
//                "#!k1= 1. 0. 0. \n"
//                "#!k2= 1. 0. 0.\n"
//                "#!formfactorname=results/formfactor.out\n"
//                "#!scatfilename=hklE.dat\n"
                ,emin,emin+nofstps*estp,nofstps);
   fclose(file);
   // 3. start bfk and create bfk0.res
   // currently coupling constant to conduction electrons  - get it from bfk.par "#!g="
    file=fopen_errchk("./bfk.par","r");
        instr[0]='#';  double g;
      while(instr[strspn(instr," \t")]=='#'&&instr[strspn(instr," \t#")]!='!'){fgets(instr,MAXNOFCHARINLINE,file);}
      extract(instr,"g",g); 
    fclose(file);
   sprintf(command,"bfk %g %g 0 1 results/%s.levels.cef results/bfk.par",g,T,sipffilename);
   printf("command: %s\n",command);
   system(command);
   // 4. read in chi0 from bfk0.res
    file=fopen_errchk("results/bfk0.res","r");
   float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE; 
   while(inputline(file,nn)==0);
    for(int i=0;i<nofstps;++i)
      {if(fabs(nn[1]-(emin+i*estp))>0.01){fprintf(stderr,"Error mcdisp -r reading bfk0.res energies %g %g not consistent\n",nn[1],emin+i*estp);exit(EXIT_FAILURE);}
       inputline(file,nn);(*chi0pointer[i])(1,1)=complex<double>(nn[1],nn[2]);(*chi0pointer[i])(1,2)=complex<double>(nn[3],nn[4]);(*chi0pointer[i])(1,3)=complex<double>(nn[5],nn[6]);
       inputline(file,nn);(*chi0pointer[i])(2,1)=complex<double>(nn[1],nn[2]);(*chi0pointer[i])(2,2)=complex<double>(nn[3],nn[4]);(*chi0pointer[i])(2,3)=complex<double>(nn[5],nn[6]);
       inputline(file,nn);(*chi0pointer[i])(3,1)=complex<double>(nn[1],nn[2]);(*chi0pointer[i])(3,2)=complex<double>(nn[3],nn[4]);(*chi0pointer[i])(3,3)=complex<double>(nn[5],nn[6]);
    //myPrintComplexMatrix(stdout,(*chi0pointer[i]));  
     inputline(file,nn);
      }
   fclose(file);
   }
  } else {
  if(epsilon>0){ // use internal chi0
  ComplexVector u1(1,nofcomponents);float dd; int n,nd;
  ComplexMatrix M(1,nofcomponents,1,nofcomponents);
  du1calc(T,Hxc,Hext,u1,dd,n,nd,ests);  
  complex<double> eps(epsilon*3,0),cc,imag(0,1),d(dd,0);
  if(dd>SMALL_QUASIELASTIC_ENERGY)
  { if(delta<0){ //treat correctly energy gain of neutron
                u1=u1.Conjugate();d=-d;
                M=-u1^u1;
               } else {
                M=u1^u1;
               }
  for(int i=0;i<nofstps;++i){
     complex<double> z(emin+i*estp,epsilon);    
     cc=1.0/(d-z);(*chi0pointer[i])+=cc*M;       
                            } //i
  }else{
     //quasielastic intensity ...  artificially we introduce a splitting epsilon !!! compare Jensen 91 p 158
     // factor 0.5 because every transition is counted as half positive and half negative energy...
   M=u1^u1;
   for(int i=0;i<nofstps;++i){
     complex<double> z(emin+i*estp,epsilon);    
     //  cc=eps/(eps-z);(*chi0pointer[i])+=cc*M;
     cc=0.5*eps/(eps-z);(*chi0pointer[i])+=cc*M;
     cc=0.5*eps/(eps+z);(*chi0pointer[i])+=cc*M.Transpose();
                            } //i
  } 
  }  
 }} //qcounter
 return 0; // success
}

/****************************************************************************/
// initialises matrix est and returns in it eigenvalues, boltzmann population and eigenstates matrix parameters of ion
/****************************************************************************/
ComplexMatrix & jjjpar::eigenstates (Vector &  Hxc,Vector & Hext,double & T)
{switch (module_type)
  {case 0:  if(estates!=NULL){(*estates)(&est,&Hxc,&Hext,&gJ,&T,&ABC,&sipffilename);}
            return est;break;
   case 2:
   case 4: (*iops).cfeigenstates(&est,Hxc,Hext,T);return est;break;
   case 5: cluster_est(&est,Hxc,Hext,T);return est;break;
   default: est=ComplexMatrix(0,2,1,2);est=0;return est;
  }
}

void jjjpar::print_eigenstates(FILE *fout)
{fprintf(fout,"#! Eigenvalues = ");
 Vector ev(Real(est.Row(0))(1,est.Chi())); myPrintVector(fout,ev);
 fprintf(fout,"#Eigenvectors [as colunmns]\n");
 ComplexMatrix es(est(1,est.Rhi(),1,est.Chi()));
 if(strstr(modulefilename,"ic1ion.so")!=NULL) {
    es=es.Transpose(); }
 myPrintComplexMatrix(fout,es);
//----------------------------------------------------------------------------//
// Submatrix extraction 
//----------------------------------------------------------------------------//

//Matrix Matrix::operator () (int rlo, int rhi, int clo, int chi) const
//
// The elements of this matrix within the index range [rlo..rhi,clo..chi] 

}  

/****************************************************************************/
// initialises matrix Icalc_parstorage and returns eigenvalues, boltzmann population and eigenstates matrix parameters of ion
/****************************************************************************/
ComplexMatrix & jjjpar::Icalc_parameter_storage_init (Vector &  Hxc,Vector & Hext,double & T)
{switch (module_type)
  {case 0:  if(Icalc_parameter_storage!=NULL){(*Icalc_parameter_storage)(&Icalc_parstorage,&Hxc,&Hext,&gJ,&T,&ABC,&sipffilename);}
            return Icalc_parstorage;break;
   case 2:
   case 4: (*iops).cfeigenstates(&Icalc_parstorage,Hxc,Hext,T);return Icalc_parstorage;break;
   default: Icalc_parstorage=ComplexMatrix(0,2,1,2);Icalc_parstorage=0;return Icalc_parstorage;
  }
}
/****************************************************************************/
// returns operator matrices (n=0 Hamiltonian, n=1,...,nofcomponents: operators of moment components)
/****************************************************************************/
Matrix jjjpar::opmat(int n,Vector &  Hxc,Vector & Hext)
{
 int retval;
 if(n>=0){
 switch (module_type)
  {
   case 0:  if(opmatM[n]==0) {
               if(dyn_opmat!=NULL) {
                  opmatM[n] = new Matrix;
                  retval=(*dyn_opmat)(&n, &sipffilename, &Hxc, &Hext, opmatM[n]); 
                  if(retval==0) { return *opmatM[n]; }
                  else { 
                     fprintf(stderr,"ERROR operator calculation in module jjjpar - opmat failed in module %i\n",module_type); exit(EXIT_FAILURE); } }
               else { 
                  fprintf(stderr,"ERROR operator calculation in module jjjpar - opmat function not defined for module %i\n",module_type); exit(EXIT_FAILURE); } }
            else {
               if(n==0) {  // Need to recalculate Hamiltonian since have different Hxc, Hext from initialisation.
                 if((*dyn_opmat)(&n, &sipffilename, &Hxc, &Hext, opmatM[n])!=0) { fprintf(stderr,"ERROR operator calculation in module jjjpar - opmat failed in module %i\n",module_type); } 
               } 
               return *opmatM[n];
            }
            break;
   case 1:  return krameropmat(n,Hxc,Hext);break;
   case 4:  return (*iops).opmat(n,Hxc,Hext);break;
   default: fprintf(stderr,"ERROR operator calculation in module jjjpar - opmat function not defined for module %i\n",module_type);exit(EXIT_FAILURE);
  }
 } else //n<0 ... matrix within eigenstates
 { // get eigenstates to matrix
  int i1,j1;
   ComplexMatrix es(est(1,est.Rhi(),1,est.Chi()));
   ComplexMatrix In(1,est.Rhi(),1,est.Chi());
  if(strstr(modulefilename,"ic1ion.so")!=NULL) { // this is because 
//ic1ion stores eigenvector in rows instead of columns
                                                es=es.Transpose(); }
  
  if(module_type==1){// this is because kramer module does not use est matrix at all
                     Matrix Ham(this->opmat(0,Hxc,Hext)); // get Hamiltonian and diagonalise
                     for (i1=In.Rlo();i1<=In.Rhi();++i1){ 
                     for (j1=In.Clo();j1<=In.Chi();++j1) { 
                     if(i1<j1)In(i1,j1)=complex <double> (Ham(j1,i1),-Ham(i1,j1));
                         else In(i1,j1)=complex <double> (Ham(i1,j1),Ham(j1,i1)); 
                     }
                     }
                     Vector E(1,est.Rhi());int sort,maxiter=1000000;
                     myEigenSystemHermitean (In,E,es,sort=1,maxiter);
                    }
  Matrix I(this->opmat(-n,Hxc,Hext)); 
  Matrix mat1(1,est.Rhi(),1,est.Chi());
  ComplexMatrix M(1,est.Rhi(),1,est.Chi());
  for (i1=I.Rlo();i1<=I.Rhi();++i1){ 
    for (j1=I.Clo();j1<=I.Chi();++j1) { 
    if(i1<j1)In(i1,j1)=complex <double> (I(j1,i1),-I(i1,j1));else In(i1,j1)=complex <double> (I(i1,j1),I(j1,i1)); 
    }
    }
    
  M=es.Transpose().Conjugate() * In * es;
  for(i1=M.Rlo();i1<=M.Rhi();++i1){for(j1=M.Clo();j1<=M.Chi();++j1){
    mat1(j1,i1)=imag(M(i1,j1)); 
    mat1(i1,j1)=real(M(i1,j1));
   }}
 //myPrintComplexMatrix(stdout,mat1);
  return mat1;
 }
}





