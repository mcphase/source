/**************************************************************************
 *
 * mcphas - program to calculate static magnetic properties (phase diagram)
 *
 * reference: M. Rotter JMMM 272-276 (2004) 481
 **************************************************************************/

#include<mcphas.h>


int verbose=0;
// for statistics 
int nofmaxloopDIV=0,nofmaxspinchangeDIV=0;
int successrate=0;
int nofcalls=0;

const char * filemode="w";

#include "myev.h"
#include "mcphas_htcalc.c"
#include "mcphas_fecalc.c"
#include "mcphas_physpropcalc.c"

// main program
int main (int argc, char **argv)
{ std::clock_t startcputime = std::clock();
  FILE * fin=NULL; 
  char outfilename[MAXNOFCHARINLINE];
  int im,j,l;
  int nofstapoints=0,noffailedpoints=0;
  int options=1; // this integer indicates how many command strings belong to 
                 //options (=1+number of option-strings)
  float x,y,dumm;
  double z,u;
  double T;
  float nn[20];nn[0]=19;
  double sta=0;
  double stamax=1e33;
  Vector xv(1,3);
  Vector yv(1,3);
  Vector h(1,3);
  
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"*\n");
fprintf(stderr,"* mcphas - program to calculate static magnetic properties (phase diagram)\n");
fprintf(stderr,"*\n");
fprintf(stderr,"* reference: M. Rotter JMMM 272-276 (2004) 481\n");
fprintf(stderr,"**************************************************************************\n\n");
  
// check command line
int errexit=0;char prefix [MAXNOFCHARINLINE];prefix[0]='\0';
  for (im=0;im<=argc-1;++im)
  {if (strcmp(argv[im],"-v")==0) {verbose=1;if (options<im)options=im;}// set verbose mode on
   if (strcmp(argv[im],"-h")==0) errexit=1; // display help message
   if (strcmp(argv[im],"-a")==0) {filemode="a";if (options<im)options=im;} // append output files
   if (strcmp(argv[im],"-stamax")==0&&im+1<=argc-1)
                                 {stamax=strtod (argv[im+1], NULL); // read stamax
                                  if (options<im+1)options=im+1;}
   if (strcmp(argv[im],"-prefix")==0&&im+1<=argc-1)
                                 {strcpy(prefix,argv[im+1]); // read prefix
                                  fprintf(stdout,"#prefix for input/ouput filenames: %s\n",prefix);
 				 if (options<im+1)options=im+1;}
  }
    inipar ini("mcphas.ini",prefix);
    if(errexit==1)ini.errexit();

  if (ini.exit_mcphas!=0)
  {ini.exit_mcphas=0;ini.print();} // if exit was 1 - save parameters and set exit=0
   strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);strcpy(prefix+11+strlen(ini.prefix),"mcphas.ini"); 
   ini.print(prefix);  // copy mcphas.ini to results directory


// as class par load  parameters from file
 strcpy(prefix,ini.prefix);strcpy(prefix+strlen(ini.prefix),"mcphas.j");
  if(fopen(prefix,"rb")==NULL)strcpy(prefix,"mcphas.j");
 if(verbose==1){printf("reading parameters from file %s\n",prefix);}
 par inputpars(prefix); 
// here save single ion property files to results
  strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);inputpars.save_sipfs(prefix); 
  strcpy(prefix+11+strlen(ini.prefix),"mcphas.j");inputpars.save(prefix);

  Vector Imax(1,inputpars.nofatoms*inputpars.nofcomponents);
  Vector Imom(1,inputpars.nofcomponents);
  Vector mmax(1,3*inputpars.nofatoms);
  Vector mmom(1,3);
  Vector h1(1,inputpars.nofcomponents),h1ext(1,3);h1ext=0;
 

//determine saturation momentum (used for scaling the plots, generation of qvectors)
T=1.0;for(l=1;l<=inputpars.nofatoms;++l){h1=0;(*inputpars.jjj[l]).Icalc_parameter_storage_init(h1,h1ext,T); // initialize eigenstate matrix
      for (im=1;im<=inputpars.nofcomponents;++im){h1ext=0;h1=0;h1(im)=20*MU_B; //just put some high field
                            (*inputpars.jjj[l]).Icalc(Imom,T,h1,h1ext,z,u,(*inputpars.jjj[l]).Icalc_parstorage);
                            Imax(inputpars.nofcomponents*(l-1)+im)=Imom(im);
                           // printf("Imax(%i)=%g\n",inputpars.nofcomponents*(l-1)+im,Imax(inputpars.nofcomponents*(l-1)+im));
			   }
      for (im=1;im<=3;++im){h1=0;h1ext=0;h1ext(im)=20; //just put some high field
                          (*inputpars.jjj[l]).mcalc(mmom,T,h1,h1ext,(*inputpars.jjj[l]).Icalc_parstorage);
                            mmax(3*(l-1)+im)=mmom(im);
                           // printf("mmax(%i)=%g\n",3*(l-1)+im,mmax(3*(l-1)+im));
			   }
                                  }

T=0.0;h=0;
// load testspinconfigurations (nooftstspinconfigurations,init-file,sav-file)
    strcpy(prefix,ini.prefix);strcpy(prefix+strlen(ini.prefix),"mcphas.tst");
    if(fopen(prefix,"rb")==NULL)strcpy(prefix,"mcphas.tst");
    strcpy(outfilename,"./results/");strcpy(outfilename+10,ini.prefix);strcpy(outfilename+10+strlen(ini.prefix),"mcphas.phs");
    testspincf testspins (ini.maxnoftestspincf,prefix,outfilename,inputpars.nofatoms,inputpars.nofcomponents);
      
    strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);
    strcpy(prefix+11+strlen(ini.prefix),"mcphas.tst");
    testspins.save(prefix,"w");

    strcpy(outfilename,"./results/");strcpy(outfilename+10,ini.prefix);strcpy(outfilename+10+strlen(ini.prefix),"mcphas.qvc");
    qvectors testqs (ini,inputpars.rez,Imax,outfilename,inputpars.nofatoms,inputpars.nofcomponents,verbose);

// declare variable physprop (typa class physproperties)
   physproperties physprop(ini.nofspincorrs,ini.maxnofhkls,inputpars.nofatoms,inputpars.nofcomponents);
                        
if (argc>options+1){ini.xv=0;ini.yv=0;fin=fopen_errchk (argv[argc-1],"rb");}   //input from file
// loop different H /T points in phase diagram
for (x=ini.xmin;x<=ini.xmax;x+=ini.xstep)
 { //begin initialize display file
   strcpy(prefix,"./results/.");strcpy(prefix+11,ini.prefix);
   strcpy(prefix+11+strlen(ini.prefix),"mcphas.fum");
   FILE * fout;fout = fopen_errchk (prefix,"w");
   fprintf (fout, " %4.4g %6.6g  %4.4g %4.4g %4.4g %4.4g %4.4g %8.8g %8.8g  %4.4g %4.4g %4.4g %4.4g\n",
            0.0,ini.ymin,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-Max(mmax),0.0,0.0);
   fprintf (fout, " %4.4g %6.6g  %4.4g %4.4g %4.4g %4.4g %4.4g %8.8g %8.8g  %4.4g %4.4g %4.4g %4.4g\n",
            0.0,ini.ymax+1e-4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,Max(mmax),0.0,0.0);
   fclose(fout); //end initialize display file
  
  for (y=ini.ymin;y<=ini.ymax;y+=ini.ystep)
  {//correct some of input errors in parameters ini (if there are any)
   if (ini.xstep<0){ini.xstep=-ini.xstep;}
   if (ini.ystep<0){ini.ystep=-ini.ystep;}
   if (ini.xmin>ini.xmax){dumm=ini.xmin;ini.xmin=ini.xmax;ini.xmax=dumm;}
   if (ini.ymin>ini.ymax){dumm=ini.ymin;ini.ymin=ini.ymax;ini.ymax=dumm;}
   
   if(argc>options+1)  //should T-H values be read from file ?
   {while (feof(fin)==0&&0==inputline(fin,nn));  // if yes -> input them
    if (feof(fin)!=0) goto endproper;
    x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
    normalizedadbdc(h,nn[4],inputpars);
   }
   else
   {//if parameters outside specified region then put them into it ...
    xv=ini.xv(1,3);normalizedadbdc(xv,1.0,inputpars);
    yv=ini.yv(1,3);normalizedadbdc(yv,1.0,inputpars);
    if (x<ini.xmin) x=ini.xmin;
    if (y<ini.ymin) y=ini.ymin;    

     T=ini.zero(0)+x*ini.xv(0)+y*ini.yv(0);
     h(1)=ini.zero(1)+x*xv(1)+y*yv(1);
     h(2)=ini.zero(2)+x*xv(2)+y*yv(2);
     h(3)=ini.zero(3)+x*xv(3)+y*yv(3);
   } 
     
      physprop.x=x;physprop.y=y;
      physprop.T=T;
      physprop.H=h;

//calculate physical properties at HT- point
   j=htcalc(physprop.H,T,ini,inputpars,testqs,testspins,physprop);
       switch (j)
       {case 0:
            //save physical properties of HT-point
	    //sta=(sta*nofstapoints+physprop.save (verbose,filemode,j,inputpars))/(nofstapoints+1);
          // 12.3.07 fancy calculation above substituted by normal summing of sta
          sta+=physprop.save (verbose,filemode,j,inputpars,ini.prefix);
   	    ++nofstapoints;
          if (sta>stamax){fprintf(stdout,"#! stamax=%g exceeded - exiting\n",stamax);goto endproper;}
	      break; 
	 case 1: goto endproper;
	      break;
         case 2: //ht calculation leads to no results- save dummy line
	         physprop.save (verbose,filemode,j,inputpars,ini.prefix);
		 sta+=1.0; // increment sta because within manifold of spincf no good solution could be found
     	      ++noffailedpoints;
	      break;	 
	 default:  ;
	}
 if(fin!=NULL){y=ini.ymin-ini.ystep;} // this is to switch off xy loop if xy points are read from file

    }
 }  
endproper:
  testspins.save(filemode);testqs.save(filemode);
   if(argc>options+1) fclose(fin);
   printf("#RESULTS saved in directory ./results/  - files:\n");
   printf("#  %smcphas.fum  - total magnetic moment, energy at different T,H\n",ini.prefix);
   printf("#  %smcphas.sps  - stable configurations at different T,H\n",ini.prefix);
   printf("#  %smcphas.mf   - mean fields at different T,H\n",ini.prefix);
   printf("#  %smcphas.hkl  - strong magnetic satellites, neutron diffraction intensity\n",ini.prefix);
   printf("#  %smcphas*.hkl - strong magnetic satellites, Fourier Comp.of moment in * dir\n",ini.prefix);
   printf("#  %smcphas*.j*  - JJ correlation functions (for exchange magnetostriction)\n",ini.prefix);
   printf("#  %smcphas.xyt  - phasediagram (stable conf.nr, angular and multipolar moments)\n",ini.prefix);
   printf("#!  %smcphas.qvc  - ...corresponding table of all nqvc=%i qvector generated test configs\n",ini.prefix,testqs.nofqs ());
   printf("#!  %smcphas.phs  - ...corresponding table of all ntst=%i configurations (except qvecs)\n",ini.prefix,testspins.n);
   printf("#  _%smcphas.*   - parameters read from input parameter files (.tst,.ini,.j)\n",ini.prefix);
   printf("#  ...         - and a copy of the single ion parameter files used.\n\n");
   double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
   std::cout << "#! Finished in cputime=" << cpu_duration << " seconds [CPU Clock] " << std::endl;
   std::cout << "#!nofHTpoints=" << nofstapoints << " H-T points in phasediagram successfully calculated" << std::endl;
   std::cout << "#!noffailedpoints=" << noffailedpoints << " H-T points in phasediagram failed to converge " << std::endl;
   std::cout << "#!fecalc - free energy calculation was attempted noffecalccalls=" << nofcalls << " times"  << std::endl;
   std::cout << "#!fecalc - free energy calculation was successful at noffecalcsuccess=" << successrate << " times"  << std::endl;
   std::cout << "#!fecalc - free energy diverged maxnofloopsDIV=" << nofmaxloopDIV << " times because maxnofloops was reached" << std::endl;
   std::cout << "#!fecalc - free energy diverged maxspinchangeDIV=" << nofmaxspinchangeDIV << " times because maxspinchange was reached" << std::endl;

   fprintf(stdout,"#! sta=%g\n",sta);
#ifdef _THREADS
std::cout << "#! nofthreads= " << NUM_THREADS << " threads were used in parallel processing " << std::endl;
for (int ithread=0; ithread<NUM_THREADS; ithread++) delete tin[ithread];
#else
std::cout << "# mcphas was compiled without parallel processing option " << std::endl;
#endif

   fprintf(stderr,"**********************************************\n");
   fprintf(stderr,"          End of Program mcphas\n");
   fprintf(stderr," reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(stderr,"**********************************************\n");

return(0);
}

int normalizedadbdc(Vector & dadbdc,double n,par & inputpars)
   {if(Norm(dadbdc)>0.00001){ // normalize Vector dadbdc to length n
    Vector Hijk(1,3);
    Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1;
                     abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;
    dadbdc2ijk(Hijk,dadbdc,abc);
    Hijk*=n/Norm(Hijk);
    ijk2dadbdc(dadbdc,Hijk,abc);
    return true;      }
    else
   {return false;}
   }
