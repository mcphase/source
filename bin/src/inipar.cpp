// methods for class inipar to store parameters for mcphasit 
#include "inipar.hpp"
#include "../../version"
#include <martin.h>

#if defined(__linux__)
#include <sys/sysinfo.h>
#elif defined(__FreeBSD__) || defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <windows.h>
#endif

void getnofthread(int & nofthreads){
  if(nofthreads<1) { // User has not set number of threads in mcphas.ini file
    char* c_nofthreads=getenv("MCPHASE_NOFTHREADS");  // Check if system environment variable set from dos.bat/lin.bat
    if (c_nofthreads)
       nofthreads = atoi(c_nofthreads);
    else {
#if defined(__linux__)                               // System-dependent calls to find number of processors (from GotoBLAS)
       nofthreads = get_nprocs();
#elif defined(__FreeBSD__) || defined(__APPLE__)
       int m[2]; size_t len;
       m[0] = CTL_HW; m[1] = HW_NCPU; len = sizeof(int);
       sysctl(m, 2, &nofthreads, &len, NULL, 0);
#else
       SYSTEM_INFO sysinfo; GetSystemInfo(&sysinfo);
       nofthreads = sysinfo.dwNumberOfProcessors;
#endif
    }
    if(nofthreads<1||nofthreads>255) nofthreads=1;   // All else fails: use only 1 thread
  }
}

 // *************************************************************************
 // ************************ inipar *************************************
 // *************************************************************************
 // class of initial parameters for program mcphas

void inipar::errexit() // type info and error exit 
{     printf (" \n%s \n",MCPHASVERSION);
printf (" use as: mcphas \n or as: mcphas [file]\n");
printf (" [file] ... input file  with sets of x y T H Ha Hb Hc points \n");
printf (" (format as output file mcphas.xyt)\n\n");
printf (" Options: -h     print this help screen\n");
printf ("          -stamax 14  ... end mcphas if standard deviation exceeds 14\n");
printf ("          -a     append output files (do not overwrite) \n");
printf ("          -doeps refine strain epsilon selfconsistently using elastic,magnetoelastic constants \n");
printf ("                 read from mcphas.j and mcphas.djdx mcphas.djdy and mcphas.djdz \n");
printf ("          -linepscf with -doeps use zero strain single ion Hamiltoanian for every mean field iteration\n");
printf ("          -linepsjj with -doeps use zero strain two ion interaction Hamiltoanian for every mean field iteration\n");
printf ("          -prefix 001    try to read files starting with 001, e.g.\n");
printf (" 		    001mcphas.ini, if these exist, otherwise take\n"); 
printf (" 		    standard input files, check if in mcphas.ini there are\n");
printf ("                   parameters such as 001xmin and use those. Output goes to files\n");
printf (" 		    results/001mcphas*.* (option for parallel processes)\n");
printf ("          -read  001    try to read output files of previous calculation starting\n");
printf (" 		     with 001, e.g. results/001mcphas.fum, if these exist, and if \n"); 
printf (" 		    the x-y-T-Ha-Hb-Hc point is found and the calculation was stable \n");
printf ("                   (i.e. free energy in results/001mcphas.fum not zero do not \n");
printf ("                   recalculate it but take results from this previous calculation\n");
printf (" 		     and store those. Option to recalculate nonstable points only.\n");
printf ("          -v     verbose mode: \n");
printf ("                 * more information is printed to stdout, \n");
printf (" 		  * the qvectors file mcphas.qom will contain \n");
printf (" 		    the explicit spinconfigurations\n");
printf (" 		  * ./results/.sps.eps will be updated not only \n");
printf (" 		    when a H-T point has been finished but always \n");
printf (" 		    when a structure with smaller free energy \n");
printf (" 		    has been stabilized\n");
printf (" Note: files which must be in current directory -\n");
printf ("       ./mcphas.ini, ./mcphas.j, directory ./results\n\n");
      exit (EXIT_FAILURE);
} 

void inipar::time_estimate_until_end(double x, double y)
{// estimate time until end 
    int nofpoints=nofstapoints+noffailedpoints+1; // add 1 to avoid zero
    int nofysteps=(int)((ymax-ymin)/ystep); if(nofysteps==0){nofysteps=1;}
    int pointstodo=nofysteps*int((xmax-x)/xstep)+int((ymax-y)/ystep);
    print_time_estimate_until_end(pointstodo/nofpoints);
    //printf("%i  %i HTpoints to do.",nofysteps,pointstodo);

 
}

//load parameters from file
int inipar::load ()
{ FILE *fin_coq;
  char instr[MAXNOFCHARINLINE];
  errno = 0;startcputime= std::clock();
  fin_coq = fopen(savfilename, "rb");
  if (fin_coq==NULL) return 1;
  xv=0;yv=0;xmin=1;xmax=0;ymin=1;ymax=0;xstep=0;ystep=0;zero=0;
  qmin(1)=1;qmin(2)=1;qmin(3)=1;qmax=0;deltaq=0;maxqperiod=0;maxnofspins=0;nofrndtries=0;
  maxnofmfloops=0;maxstamf=0;bigstep=0;maxspinchange=0;nofthreads=0;
  nofspincorrs=0;maxnofhkls=0;maxQ=0;maxnoftestspincf=1000;
  
  while (fgets(instr,MAXNOFCHARINLINE,fin_coq)!=NULL)
  {if(instr[strspn(instr," \t")]!='#'&&instr[strspn(instr," \t")]!='[') // comment lines headed by # or [ are ignored in mcphas.ini
   {extract_with_prefix(instr,prefix,"exit",exit_mcphas);extract_with_prefix(instr,prefix,"pause",pause_mcphas);
    extract_with_prefix(instr,prefix,"displayall",displayall);extract_with_prefix(instr,prefix,"logfevsQ",logfevsQ); 
     
    extract_with_prefix(instr,prefix,"xT",xv[0]);      extract_with_prefix(instr,prefix,"xHa",xv[1]);
    extract_with_prefix(instr,prefix,"xHb",xv[2]);    extract_with_prefix(instr,prefix,"xHc",xv[3]);
    extract_with_prefix(instr,prefix,"xmin",xmin);  extract_with_prefix(instr,prefix,"xmax",xmax);
    extract_with_prefix(instr,prefix,"xstep",xstep);
   
    extract_with_prefix(instr,prefix,"yT",yv[0]);     extract_with_prefix(instr,prefix,"yHa",yv[1]);
    extract_with_prefix(instr,prefix,"yHb",yv[2]);   extract_with_prefix(instr,prefix,"yHc",yv[3]);
    extract_with_prefix(instr,prefix,"ymin",ymin); extract_with_prefix(instr,prefix,"ymax",ymax);
    extract_with_prefix(instr,prefix,"ystep",ystep);   
 
    extract_with_prefix(instr,prefix,"T0",zero(0));
    extract_with_prefix(instr,prefix,"Ha0",zero(1));
    extract_with_prefix(instr,prefix,"Hb0",zero(2));
    extract_with_prefix(instr,prefix,"Hc0",zero(3));
    
   
    extract_with_prefix(instr,prefix,"hmin",qmin[1]); 
    extract_with_prefix(instr,prefix,"kmin",qmin[2]); 
    extract_with_prefix(instr,prefix,"lmin",qmin[3]); 
    extract_with_prefix(instr,prefix,"hmax",qmax[1]); 
    extract_with_prefix(instr,prefix,"kmax",qmax[2]); 
    extract_with_prefix(instr,prefix,"lmax",qmax[3]); 
    extract_with_prefix(instr,prefix,"deltah",deltaq[1]); 
    extract_with_prefix(instr,prefix,"deltak",deltaq[2]); 
    extract_with_prefix(instr,prefix,"deltal",deltaq[3]); 
    extract_with_prefix(instr,prefix,"maxqperiod",maxqperiod);
    extract_with_prefix(instr,prefix,"maxnofspins",maxnofspins);
    extract_with_prefix(instr,prefix,"nofrndtries",nofrndtries);

    extract_with_prefix(instr,prefix,"maxnofmfloops",maxnofmfloops);
    extract_with_prefix(instr,prefix,"maxstamf",maxstamf); 
    extract_with_prefix(instr,prefix,"bigstep",bigstep); 
    extract_with_prefix(instr,prefix,"maxspinchange",maxspinchange); 
    extract_with_prefix(instr,prefix,"maxnoftestspincf",maxnoftestspincf);

    extract_with_prefix(instr,prefix,"nofthreads",nofthreads);

    extract_with_prefix(instr,prefix,"nofspincorrs",nofspincorrs); 
    extract_with_prefix(instr,prefix,"maxnofhkls",maxnofhkls); 
    extract_with_prefix(instr,prefix,"maxQ",maxQ); 
    }
   }
  fclose (fin_coq);

  if (Norm(xv)==0){fprintf(stderr,"ERROR reading xT xHa xHb xHc\n");return 1;}
  if (Norm(yv)==0){fprintf(stderr,"ERROR reading yT yHa yHb yHc\n");return 1;}
  if (xmin>xmax){fprintf(stderr,"ERROR reading xmin xmax\n");return 1;}
  if (ymin>ymax){fprintf(stderr,"ERROR reading ymin ymax\n");return 1;}
  if (xstep==0){fprintf(stderr,"Warning reading xstep: xstep=0\n");}
  if (ystep==0){fprintf(stderr,"Warning reading ystep: ystep=0\n");}

  if(qmin(1)>qmax(1)){fprintf(stderr,"ERROR reading hmin hmax\n");return 1;}
  if(qmin(2)>qmax(2)){fprintf(stderr,"ERROR reading kmin kmax\n");return 1;}
  if(qmin(3)>qmax(3)){fprintf(stderr,"ERROR reading lmin lmax\n");return 1;}
  if(Norm(deltaq)==0){fprintf(stderr,"Warning reading deltah k l: deltah=deltak=deltal=0\n");}
  if(deltaq[1]==0){fprintf(stderr,"ERROR reading deltah=0: deltah must be >0\n");return 1;}
  if(deltaq[2]==0){fprintf(stderr,"ERROR reading deltak=0: deltak must be >0\n");return 1;}
  if(deltaq[3]==0){fprintf(stderr,"ERROR reading deltal=0: deltal must be >0\n");return 1;}
  if(maxqperiod==0){fprintf(stderr,"Warning reading maxqperiod=0\n");}
  if(nofrndtries==0){fprintf(stderr,"Warning reading nofrndtries=0\n");}
  if (maxnofspins==0){maxnofspins=maxqperiod*maxqperiod*maxqperiod;
                      fprintf(stderr,"warning ... reading maxnofspins=0: putting it to %i\n",maxnofspins);}
  if (maxnoftestspincf<1){fprintf(stderr,"ERROR maxnoftestspincf<1 not possible\n");return 1;}
  getnofthread(nofthreads);
  
  if(maxnofmfloops==0){fprintf(stderr,"Error reading maxnofmfloops\n");return 1;}
  if(maxnofmfloops==1){fprintf(stderr,"#! Reading maxnofmfloops=1 - mean fields will be calculated from initial spins and free energy will be evaluated using initial spins\n");}
  if(maxnofmfloops==2){fprintf(stderr,"#! Reading maxnofmfloops=2 - mean fields will be calculated from initial spins, new spins will be calculated  and using these the free energy will be evaluated\n");}
  if(maxstamf==0){fprintf(stderr,"Error reading maxstamf\n");return 1;}
  if(bigstep==0){fprintf(stderr,"Error reading bigstep\n");return 1;}
  if(maxspinchange==0){fprintf(stderr,"Error reading maxspinschange\n");return 1;}

  if(nofspincorrs==0){fprintf(stderr,"Warning reading nofspincorrs=0 - no spin correlation functions will be calculated\n");}
  if(maxnofhkls==0){fprintf(stderr,"Warning reading maxnofhkls=0 - no magnetic neutron reflections  will be calculated\n");}
  if(maxQ==0){fprintf(stderr,"Warning reading maxQ=0:magnetic neutron reflections  will be calculated only in primitive cell of reciprocal lattice\n");}

return 0;
}



void inipar::print () // printout initial parameters to file 
{print(savfilename);}

void inipar::print (const char * filename)
{
 FILE * fout;
// we should print to a file all used configurations
 fout = fopen_errchk (filename,"w");
    fprintf(fout,"# Parameters for meanfield calculation - module %s\n#<!--mcphase.mcphas.ini-->\n",MCPHASVERSION);
    fprintf(fout,"#*********************************************************\n");
    fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
    fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
    fprintf(fout,"#**********************************************************\n"); 
    fprintf(fout,"[MCPHASE RUNTIME CONTROL]\n");
    fprintf(fout,"# to stop program set exit to 1 \
                  \nexit=%i \
                  \n# to hold program set pause to 1 \
                  \npause=%i \
                  \n# to display all structures while iterating set displayall to 1\n# (mind that by using this option mcphas gets very slow) \
                  \ndisplayall=%i \
                  \n# to create a logfile of the propagation versus free energy  set logfevsQ to 1\n# (mind this uses a lot of disc space) \
                  \nlogfevsQ=%i\n\n",exit_mcphas,pause_mcphas,displayall,logfevsQ);

    fprintf(fout,"[XY PHASEDIAGRAM PARAMETERS]\n");
    fprintf(fout,"#xy phasediagram axes - parameters\n");
    fprintf(fout,"# structures are calculated in the xy - phasediagram \
                  \n# the direction of x and y can be chosen:     \
                  \n# vector in (H-T) space corresponding to x axis (xT [K] xHa [T] xHb [T] xHc [T])\n");
    fprintf(fout,"xT=%g\nxHa=%g\nxHb=%g\nxHc=%g\n# range of x\nxmin=%g\nxmax=%g\nxstep=%g\n",
           xv(0), xv(1), xv(2), xv(3), xmin,  xmax,  xstep);
    fprintf(fout,"# vector in (H-T) space corresponding to y axis (yT [K] yHa [T] yHb [T] yHc [T])\n");

    fprintf(fout,"yT=%g\nyHa=%g\nyHb=%g\nyHc=%g\n# range of y\nymin=%g\nymax=%g\nystep=%g\n",
           yv(0), yv(1), yv(2), yv(3), ymin,  ymax,  ystep);
    fprintf(fout,"# offset for phase diagram\n");
    fprintf(fout,"T0=%g\nHa0=%g\nHb0=%g\nHc0=%g\n\n",zero(0),zero(1),zero(2),zero(3));       

    fprintf(fout,"input (xHa xHb xHc) (yHa yHb yHc) and (Ha0 Hb0 Hc0) are vectors\n");
    fprintf(fout,"given in terms of components with respect to unit vectors along the\n");
    fprintf(fout,"Bravais lattice ^a=a/|a|, ^b=b/|b|, ^c=c/|c|.\n");
    fprintf(fout,"For the external magnetic field unit is Tesla.\n");

    fprintf(fout,"[GENERATION OF SPIN CONFIGURATIONS]\n");
    fprintf(fout,"# test q vector (qmin qmax deltaq)\n");
    fprintf(fout,"hmin=%g\nhmax=%g\ndeltah=%g\n",qmin(1),qmax(1),deltaq(1));
    fprintf(fout,"kmin=%g\nkmax=%g\ndeltak=%g\n",qmin(2),qmax(2),deltaq(2));
    fprintf(fout,"lmin=%g\nlmax=%g\ndeltal=%g\n",qmin(3),qmax(3),deltaq(3));
    fprintf(fout,"# maximal periodicity of spinconfigurations generated by q vectors\n");
    fprintf(fout,"maxqperiod=%i\n",maxqperiod);
    fprintf(fout,"# maximal number of spins in spinconfigurations generated by q vectors\n");
    fprintf(fout,"maxnofspins=%i\n",maxnofspins);
    fprintf(fout,"# number of random (Monte Carlo) spin inversions  to try for each initial spinconfiguration\n");
    fprintf(fout,"nofrndtries=%i\n\n",nofrndtries);
    fprintf(fout,"# maximum number of test spin configurations in table\n");
    fprintf(fout,"maxnoftestspincf=%i\n\n",maxnoftestspincf);

    fprintf(fout,"[PARAMETERS FOR SUB FECALC SELFCONSISTENCY PROCESS]\
                  \n# maximum number of selfconsistency loops\n");
    fprintf(fout,"maxnofmfloops=%i\n",maxnofmfloops);
    fprintf(fout,"# standard deviation - limit to end selfconsistency process \
                  \n# standard deviation is defined by ...sta=sqrt(sum_{i=1}^{n} (newmf-old mf)i^2/n) \
		  \n# the meanfield is given by mf=gj mb H [meV] (gj...lande factor, mb... bohr magneton)\n");
    fprintf(fout,"maxstamf=%g\n",maxstamf);
    fprintf(fout,"# mean field step ratio (bigstep=actual step/calculated step<1) to perform actually\n");
    fprintf(fout,"# note: if sta increases - then for 10 iterations set step ratio to smallstep=bigstep/n\n");
    fprintf(fout,"# by default n=5. However, if bigstep>1 then n=integervalue(bigstep) and step ratio=bigstep-n \n");
    fprintf(fout,"bigstep=%g\n",bigstep);

    fprintf(fout,"# sum_{i=1}^{n} abs(actual change of angular momentum <Ji> with respect to \
                  \n# initial  configuration) > maxspinchange will  end selfconsistency process\n");
    fprintf(fout,"maxspinchange=%g\n\n",maxspinchange);

    fprintf(fout,"[OUTPUT OF PHYSICAL PROPERTIES]\n");
    fprintf(fout,"#output of physical properties to compare with experiment\n");
    fprintf(fout,"# 1. For thermal expansion and magnetostriction \
                  \n#  how many spinspin correlation functions  \
                  \n#  should be calculated \n");
    fprintf(fout,"nofspincorrs=%i\n",nofspincorrs);
    fprintf(fout,"# 2. For Neutron Diffraction \
                 \n#  calculation of mxnofhkl strongest reflections\n");
    fprintf(fout," maxnofhkls=%i\n",maxnofhkls);
    fprintf(fout,"#  maximum scattering vector |Q|[1/A] for calculated hkl's\n");
    fprintf(fout," maxQ=%g\n",maxQ);

  fclose(fout);
}

//constructor ... load initial parameters from file
inipar::inipar (const char * file,char * pref)
{ savfilename= new char [strlen(file)+strlen(pref)+1];
  if(pref[0]!='\0')strcpy(savfilename,pref);
  strcpy(savfilename+strlen(pref),file);
  prefix = new char[strlen(pref)+1];
  strcpy(prefix,pref);
  xv=Vector(0,3);yv=Vector(0,3);zero=Vector(0,3);
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
  doeps=0;linepscf=0;linepsjj=0;ipx=NULL;ipy=NULL;ipz=NULL;
  printf("reading file %s\n",savfilename);
  if(load()!=0){if(pref[0]!='\0'){fprintf(stderr,"File %s not found - trying %s\n",savfilename,file);
                strcpy(savfilename,file);}
                if(load()!=0){fprintf(stderr,"# Warning: Cannot load file %s - using default values ! \n",savfilename); 
  // set default values
  xv=0;xv(0)=1;yv=0;yv(3)=1;xmin=1;xmax=1;ymin=0;ymax=0;xstep=1;ystep=1;
  qmin=0;qmax=0;deltaq(1)=0.1;deltaq(2)=0.1;deltaq(3)=0.1;maxqperiod=1;maxnofspins=10;nofrndtries=0;
  maxnofmfloops=100;maxstamf=1e-3;bigstep=1;maxspinchange=100;zero=0;
nofthreads=0;getnofthread(nofthreads);
  nofspincorrs=0;maxnofhkls=5;maxQ=3;maxnoftestspincf=1000;
  nofstapoints=0;
  nofmaxloopDIV=0;nofmaxspinchangeDIV=0;
  successrate=0;
  nofcalls=0;
  noffailedpoints=0;
  print();
                              }
                }
}

//kopier-konstruktor 
inipar::inipar (const inipar & p)
{ savfilename= new char [strlen(p.savfilename)+1];
  strcpy(savfilename,p.savfilename);
  prefix = new char[strlen(p.prefix)+1];
  strcpy(prefix,p.prefix);
  doeps=p.doeps;
  linepscf=p.linepscf;
  linepsjj=p.linepsjj;
  ipx=p.ipx;
  ipy=p.ipy;
  ipz=p.ipz;
  startcputime=p.startcputime;
  nofstapoints=p.nofstapoints;
  nofmaxloopDIV=p.nofmaxloopDIV;
  nofmaxspinchangeDIV=p.nofmaxspinchangeDIV;
  successrate=p.successrate; 
  nofcalls=p.nofcalls;
  noffailedpoints=p.noffailedpoints;
  exit_mcphas=p.exit_mcphas;pause_mcphas=p.pause_mcphas;
  displayall=p.displayall;logfevsQ=p.logfevsQ;
  
  
  xv=Vector(0,3);yv=Vector(0,3);
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
  xv=p.xv;xmin=p.xmin;xmax=p.xmax;xstep=p.xstep;
  yv=p.yv;ymin=p.ymin;ymax=p.ymax;ystep=p.ystep;
  
  qmin=p.qmin;
  qmax=p.qmax;
  deltaq=p.deltaq;  
  maxqperiod=p.maxqperiod;
  maxnofspins=p.maxnofspins;
  nofrndtries=p.nofrndtries;
  maxnoftestspincf=p.maxnoftestspincf;

  maxnofmfloops=p.maxnofmfloops;
  maxstamf=p.maxstamf;
  bigstep=p.bigstep;
  maxspinchange=p.maxspinchange;
  
  nofspincorrs=p.nofspincorrs;
  maxnofhkls=p.maxnofhkls;
  maxQ=p.maxQ;
}

//destruktor
inipar::~inipar ()
{//printf("hello destruktor inipar\n");  
 
delete []savfilename;
delete []prefix;
//printf("hello destruktor inipar\n");  
 }
