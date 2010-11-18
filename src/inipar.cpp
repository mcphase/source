// methods for class inipar 
#include "inipar.hpp"
#include "../../version"
#include <martin.h>




 // *************************************************************************
 // ************************ inipar *************************************
 // *************************************************************************
 // class of initial parameters for program mcphas

void inipar::errexit() // type info and error exit 
{     printf (" \n			%s \n",MCPHASVERSION);
printf (" 		use as: mcphas \n or as: mcphas [file]\n");
printf ("                [file] ... input file  with sets of x y T H Ha Hb Hc points \n");
printf (" 	       (format as output file mcphas.xyt)\n\n");
printf (" 	       Options: -h     print this help screen\n");
printf (" 	                -stamax 14  ... end mcphas if standard deviation exceeds 14\n");
printf (" 	                -v     verbose mode: \n");
printf (" 			          * more information is printed to stdout, \n");
printf (" 			          * the qvectors file mcphas.qom will contain \n");
printf (" 				    the explicit spinconfigurations\n");
printf (" 			          * ./results/.sps.eps will be updated not only \n");
printf (" 				    when a H-T point has been finished but always \n");
printf (" 				    when a structure with smaller free energy \n");
printf (" 				    has been stabilized\n\n");
printf (" 	       Note: files which must be in current directory -\n");
printf (" 	       ./mcphas.ini, ./mcphas.j, directory ./results\n\n");
      exit (EXIT_FAILURE);
} 


//load parameters from file
int inipar::load ()
{ FILE *fin_coq;
  char instr[MAXNOFCHARINLINE];
  errno = 0;
  fin_coq = fopen(savfilename, "rb");
  if (fin_coq==NULL) return 1;
  xv=0;yv=0;xmin=1;xmax=0;ymin=1;ymax=0;xstep=0;ystep=0;
  qmin(1)=1;qmin(2)=1;qmin(3)=1;qmax=0;deltaq=0;maxqperiod=0;maxnofspins=0;nofrndtries=0;
  maxnofmfloops=0;maxstamf=0;bigstep=0;maxspinchange=0;
  nofspincorrs=0;maxnofhkls=0;maxQ=0;
  
  while (fgets(instr,MAXNOFCHARINLINE,fin_coq)!=NULL)
  {if(instr[strspn(instr," \t")]!='#'&&instr[strspn(instr," \t")]!='[') // comment lines headed by # or [ are ignored in mcphas.ini
   {extract(instr,"exit",exit_mcphas);extract(instr,"pause",pause_mcphas);
    extract(instr,"displayall",displayall);extract(instr,"logfevsQ",logfevsQ); 
     
    extract(instr,"xT",xv[0]);      extract(instr,"xHa",xv[1]);
    extract(instr,"xHb",xv[2]);    extract(instr,"xHc",xv[3]);
    extract(instr,"xmin",xmin);  extract(instr,"xmax",xmax);
    extract(instr,"xstep",xstep);
   
    extract(instr,"yT",yv[0]);     extract(instr,"yHa",yv[1]);
    extract(instr,"yHb",yv[2]);   extract(instr,"yHc",yv[3]);
    extract(instr,"ymin",ymin); extract(instr,"ymax",ymax);
    extract(instr,"ystep",ystep);   
 
    extract(instr,"T0",zero(0));
    extract(instr,"Ha0",zero(1));
    extract(instr,"Hb0",zero(2));
    extract(instr,"Hc0",zero(3));
    
   
    extract(instr,"hmin",qmin[1]); 
    extract(instr,"kmin",qmin[2]); 
    extract(instr,"lmin",qmin[3]); 
    extract(instr,"hmax",qmax[1]); 
    extract(instr,"kmax",qmax[2]); 
    extract(instr,"lmax",qmax[3]); 
    extract(instr,"deltah",deltaq[1]); 
    extract(instr,"deltak",deltaq[2]); 
    extract(instr,"deltal",deltaq[3]); 
    extract(instr,"maxqperiod",maxqperiod);
    extract(instr,"maxnofspins",maxnofspins);
    extract(instr,"nofrndtries",nofrndtries);

    extract(instr,"maxnofmfloops",maxnofmfloops);
    extract(instr,"maxstamf",maxstamf); 
    extract(instr,"bigstep",bigstep); 
    extract(instr,"maxspinchange",maxspinchange); 

    extract(instr,"nofspincorrs",nofspincorrs); 
    extract(instr,"maxnofhkls",maxnofhkls); 
    extract(instr,"maxQ",maxQ); 
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
  if(maxqperiod==0){fprintf(stderr,"Warning reading maxqperiod=0\n");}
  if(nofrndtries==0){fprintf(stderr,"Warning reading nofrndtries=0\n");}
  if (maxnofspins==0){maxnofspins=maxqperiod*maxqperiod*maxqperiod;
                      fprintf(stderr,"warning ... reading maxnofspins=0: putting it to %i\n",maxnofspins);}

  if(maxnofmfloops==0){fprintf(stderr,"Error reading maxnofmfloops\n");return 1;}
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


    fprintf(fout,"[PARAMETERS FOR SUB FECALC SELFCONSISTENCY PROCESS]\
                  \n# maximum number of selfconsistency loops\n");
    fprintf(fout,"maxnofmfloops=%i\n",maxnofmfloops);
    fprintf(fout,"# standard deviation - limit to end selfconsistency process \
                  \n# standard deviation is defined by ...sta=sqrt(sum_{i=1}^{n} (newmf-old mf)i^2/n) \
		  \n# the meanfield is given by mf=gj mb H [meV] (gj...lande factor, mb... bohr magneton)\n");
    fprintf(fout,"maxstamf=%g\n",maxstamf);
    fprintf(fout,"# mean field step ratio (=actual step/calculated step) to perform actually\n");
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
inipar::inipar (const char * file)
{ savfilename= new char [strlen(file)+1];
  strcpy(savfilename,file);
  xv=Vector(0,3);yv=Vector(0,3);zero=Vector(0,3);
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
  printf("reading file %s\n",file);
  if(load()!=0){fprintf(stderr,"ERROR loading file %s\n",savfilename);errexit();}
}

//kopier-konstruktor 
inipar::inipar (const inipar & p)
{ savfilename= new char [strlen(p.savfilename)+1];
  strcpy(savfilename,p.savfilename);

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
{delete []savfilename;}
