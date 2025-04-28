/**************************************************************************
 *
 * mcphas - program to calculate static magnetic properties (phase diagram)
 *
 * reference: M. Rotter JMMM 272-276 (2004) 481
 **************************************************************************/

#include<mcphas.h>


int verbose=0;
// for statistics 
int isfull=0;

const char * filemode="w";

#include "myev.h"
#include "mcphas_htcalc.c"
#include "mcphas_fecalc.c"
#include "mcphas_physpropcalc.c"

// main program
int main (int argc, char **argv)
{ int j,l,doeps=0,linepscf=0,linepsjj=0;
  int options=0; // this integer indicates how many command strings belong to 
                 //options 

  int im;
  double stamax=1e33;
  double sta=0;

fprintf(stderr,"**************************************************************************\n");
//fprintf(stderr,"*\n");
fprintf(stderr,"*  %s \n",MCPHASVERSION);
fprintf(stderr,"* program to calculate static magnetic properties (phase diagram)\n");
fprintf(stderr,"* reference: M. Rotter JMMM 272-276 (2004) 481\n");
fprintf(stderr,"**************************************************************************\n\n");
  
// check command line
int errexit=0;char prefix [MAXNOFCHARINLINE];prefix[0]='\0';
              char readprefix [MAXNOFCHARINLINE];readprefix[0]='\0';
  for (im=1;im<argc;++im)  // im=0  argv[0] is command "mcphasit"
  {if (strcmp(argv[im],"-v")==0) {verbose=1;if (options<im)options=im;}// set verbose mode on
   if (strcmp(argv[im],"-h")==0) errexit=1; // display help message
   if (strcmp(argv[im],"-doeps")==0) {doeps=1;if (options<im)options=im;} // do strain epsilon calculation
   if (strcmp(argv[im],"-linepscf")==0) {linepscf=1;if (options<im)options=im;} // do cf strain epsilon calculation linear 
   if (strcmp(argv[im],"-linepsjj")==0) {linepsjj=1;if (options<im)options=im;} // do exchange strain epsilon calculation linear
   if (strcmp(argv[im],"-a")==0) {filemode="a";if (options<im)options=im;} // append output files
   if (strcmp(argv[im],"-stamax")==0&&im+1<=argc-1)
                                 {stamax=strtod (argv[im+1], NULL); // read stamax
                                  if (options<im+1)options=im+1;}
   if (strcmp(argv[im],"-prefix")==0&&im+1<=argc-1)
                                 {strcpy(prefix,argv[im+1]); // read prefix
                                  fprintf(stdout,"#prefix for input/ouput filenames: %s\n",prefix);
 				 if (options<im+1)options=im+1;}
  if (strcmp(argv[im],"-read")==0&&im+1<=argc-1)
                                 {strcpy(readprefix,argv[im+1]); // read prefix
                                  fprintf(stdout,"#reading stable points from mcphas ouput files: results/%s*\n",readprefix);
 				 if (options<im+1)options=im+1;}
  }	
    inipars inip("mcphas.ini",prefix);   
   
    if(errexit==1)(*inip.inis[0]).errexit();
// loop for different prefixes ...
   for(int ninis=0;ninis<inip.nofinis;++ninis)
{
  FILE * fin=NULL; 
  char outfilename[MAXNOFCHARINLINE];
    float x,y,dumm;
  double z,u;
  double T;
  float nn[100];nn[0]=99;
  Vector xv(1,3);
  Vector yv(1,3);
  Vector h(1,HEXT_DIMENSION);
  

   inipar ini((*inip.inis[ninis]));
  if(inip.nofinis>1)printf("# Running McPhase with prefix %s\n",ini.prefix);
   ini.doeps=doeps;ini.linepscf=linepscf;ini.linepsjj=linepsjj;

  if (ini.exit_mcphas!=0)
  {ini.exit_mcphas=0;inip.saveexitzero();} // if exit was 1 - save parameters and set exit=0
  if(strcmp(ini.prefix,readprefix)==0&&prefix[0]!='\0')filemode="a";
   strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);strcpy(prefix+11+strlen(ini.prefix),"mcphas.ini"); 
   ini.print(prefix);  // copy mcphas.ini to results directory with filename prefixmcphas.ini


// as class par load  parameters from file
 strcpy(prefix,ini.prefix);strcpy(prefix+strlen(ini.prefix),"mcphas.j");
fin=fopen(prefix,"rb");if(fin==NULL)strcpy(prefix,"mcphas.j"); else fclose(fin);
 if(verbose==1){printf("reading parameters from file %s\n",prefix);}
 par inputpars(prefix,verbose); 
// here save single ion property files to results
  strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);inputpars.save_sipfs(prefix); 
  strcpy(prefix+11+strlen(ini.prefix),"mcphas.j");inputpars.save(prefix,0);

if (inputpars.cs.alpha()!=90||inputpars.cs.beta()!=90||inputpars.cs.gamma()!=90)
   {ini.defaultcolcode(5,4);ini.defaultcolcode(6,5);ini.defaultcolcode(7,6);} // reset default colcode in ini
  
if(doeps) {
if(ini.nofrndtries<0){fprintf(stderr,"# Error - nofrndtries<0 - Monte Carlo calculations not (yet) possible with strain epsilon.\n");exit(1); }
if(verbose==1&&linepscf){printf("option -linepscf: strain epsilon not used in diagonalisation of single ion Hamiltonian\n");}
// as class par load  parameters derivatives from file
 strcpy(prefix,ini.prefix);strcpy(prefix+strlen(ini.prefix),"mcphas.djdx");
 fin=fopen(prefix,"rb"); if(fin==NULL)strcpy(prefix,"mcphas.djdx"); else  fclose(fin);
  fin=fopen(prefix,"rb");
 if(fin!=NULL){
 if(verbose==1){printf("reading parameters from file %s\n",prefix);}
 ini.ipx= new par(prefix,verbose);  

// here save single ion property files to results
  strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);inputpars.save_sipfs(prefix); 
  strcpy(prefix+11+strlen(ini.prefix),"mcphas.djdx");inputpars.save(prefix,0);
                             
 strcpy(prefix,ini.prefix);strcpy(prefix+strlen(ini.prefix),"mcphas.djdy");
FILE*fin1;fin1=fopen(prefix,"rb");  if(fin1==NULL)strcpy(prefix,"mcphas.djdy"); else fclose(fin1);
fin1=fopen(prefix,"rb");
if(fin1!=NULL){
  if(verbose==1){printf("reading parameters from file %s\n",prefix);}
 ini.ipy= new par(prefix,verbose); 
// here save single ion property files to results
  strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);inputpars.save_sipfs(prefix); 
  strcpy(prefix+11+strlen(ini.prefix),"mcphas.djdy");inputpars.save(prefix,0);
                             
 strcpy(prefix,ini.prefix);strcpy(prefix+strlen(ini.prefix),"mcphas.djdz");
  FILE * fin2; fin2=fopen(prefix,"rb"); if(fin2==NULL)strcpy(prefix,"mcphas.djdz"); else fclose(fin2);
fin2=fopen(prefix,"rb");
if(fin2!=NULL){
  if(verbose==1){printf("reading parameters from file %s\n",prefix);}
 ini.ipz= new par(prefix,verbose);  
// here save single ion property files to results
  strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);inputpars.save_sipfs(prefix); 
  strcpy(prefix+11+strlen(ini.prefix),"mcphas.djdz");inputpars.save(prefix,0);
                             }else {fprintf(stderr,"# Error - could not open mcphas.djdz \n");exit(1); }
 fclose(fin2);
 } else {fprintf(stderr,"# Error - could not open mcphas.djdy \n");exit(1); }
 fclose(fin1);
if(verbose)fprintf(stdout,"# strain due to derivatives of 2ion interactions in mcphas.djdx .djdy .djdz will be calculated\n");
fprintf(stderr,"#Comparing mcphas.djdx .djdy .djdz and mcphas.j ...\n");
  // check if inputpars abc nofatoms  atomic positions sipffilenames nofcomponents  agree with
  // ini.ipx y z
 // check if ini ipx and ipy and ipz have the same paranz for each neighbour
 // operator ~ returns 8 7 6 5 4 3 2 1 0depending on agreement of
 //  8 abc 7 nofatoms 6 atomic positions 5 sipffilenames 4 nofcomponents 3 nofneighbours disagreement
 //  2 neighbour position 1 interaction parmeter disagreement i.e. 0 is perfect match
 if((inputpars!=(*ini.ipx))>3){fprintf(stderr,"# Error - mcphas.djdx does not match mcphas.j in nofcomponents, sipffilenames, atomic positions, nofatoms or lattice \n");exit(1);}
 if((inputpars!=(*ini.ipy))>3){fprintf(stderr,"# Error - mcphas.djdy does not match mcphas.j in nofcomponents, sipffilenames, atomic positions, nofatoms or lattice \n");exit(1);}
 if((inputpars!=(*ini.ipz))>3){fprintf(stderr,"# Error - mcphas.djdz does not match mcphas.j in nofcomponents, sipffilenames, atomic positions, nofatoms or lattice \n");exit(1);}

  if(((*ini.ipx)!=(*ini.ipy))>1){fprintf(stderr,"# Error - mcphas.djdx does not match mcphas.djdy in nofneighbours or neighbour positions\n");exit(1);}
  if(((*ini.ipx)!=(*ini.ipz))>1){fprintf(stderr,"# Error - mcphas.djdx does not match mcphas.djdz in nofneighbours or neighbour positions\n");exit(1);}
 fprintf(stderr,"# ... these are no problems, continuing\n");
 if(verbose==1&&linepsjj){printf("option -linepsj: neglecting strain dependence of two ion interactions when calculating mean fields in mean field loop\n");}
  fclose(fin);
  }        }


  Vector Imax(1,inputpars.cs.nofatoms*inputpars.cs.nofcomponents);
  Vector Imom(1,inputpars.cs.nofcomponents);
  Vector mmax(1,3*inputpars.cs.nofatoms);
  Vector mmom(1,3);
  Vector h1(1,inputpars.cs.nofcomponents),h1ext(1,HEXT_DIMENSION);h1ext=0;
 if(doeps){printf("#Inverting Elastic Constants Matrix\n");
  inputpars.Cel.Inverse();
           }
//determine saturation momentum (used for scaling the plots, generation of qvectors)
if(verbose==1){printf("determine saturation momentum running singleion calculations for different fields");}
T=1.0;for(l=1;l<=inputpars.cs.nofatoms;++l){h1=0;(*inputpars.jjj[l]).Icalc_parameter_storage_init(h1,h1ext,T); // initialize eigenstate matrix
      for (im=1;im<=inputpars.cs.nofcomponents;++im){h1ext=0;h1=0;h1(im)=20*MU_B; //just put some high field
                            (*inputpars.jjj[l]).Icalc(Imom,T,h1,h1ext,z,u,(*inputpars.jjj[l]).Icalc_parstorage);
                            Imax(inputpars.cs.nofcomponents*(l-1)+im)=Imom(im);
                           //if(verbose==1)printf("Imax(%i)=%g\n",inputpars.cs.nofcomponents*(l-1)+im,Imax(inputpars.cs.nofcomponents*(l-1)+im));
			   }
      for (im=1;im<=3;++im){h1=0;h1ext=0;h1ext(im)=20; //just put some high field
                          (*inputpars.jjj[l]).mcalc(mmom,T,h1,h1ext,(*inputpars.jjj[l]).Icalc_parstorage);
                            mmax(3*(l-1)+im)=mmom(im);
                           //if(verbose==1)printf("mmax(%i)=%g ",3*(l-1)+im,mmax(3*(l-1)+im));
			   }      }
if(verbose==1)printf("... done\n");
                            

T=0.0;h=0;
// load testspinconfigurations (nooftstspinconfigurations,init-file,sav-file)
    strcpy(prefix,ini.prefix);strcpy(prefix+strlen(ini.prefix),"mcphas.tst");
    fin=fopen(prefix,"rb");if(fin==NULL)strcpy(prefix,"mcphas.tst");else    fclose(fin);
    strcpy(outfilename,"./results/");strcpy(outfilename+10,ini.prefix);strcpy(outfilename+10+strlen(ini.prefix),"mcphas.phs");
    testspincf testspins (ini.maxnoftestspincf,prefix,outfilename,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
    strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);
    strcpy(prefix+11+strlen(ini.prefix),"mcphas.tst");
    testspins.save(prefix,"w");
    strcpy(outfilename,"./results/");strcpy(outfilename+10,ini.prefix);strcpy(outfilename+10+strlen(ini.prefix),"mcphas.qvc");
    qvectors testqs (ini,inputpars,Imax,outfilename,verbose);
// declare variable physprop (typa class physproperties)
   physproperties physprop(ini.nofspincorrs,ini.maxnofhkls,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);

// *** calculate electrical polarisation (static, from ionic charges and positions, electret matter)               
  if(fabs(inputpars.totalcharge)<SMALLCHARGE)
{Vector rijk(1,3);
for (l=1;l<=inputpars.cs.nofatoms;++l){
// Problem : surface effect !! how to calculate polarisation !????
// try - use wigner seitz cell 
    Vector dr1r2r3(1,3),xyz(1,3);
    dr1r2r3=inputpars.rez*(*inputpars.jjj[l]).xyz; // transform dadbdc to primitive lattice
    for(int n=1;n<=3;++n){
    while(dr1r2r3(n)>0.5)dr1r2r3(n)-=1.0; // shift position into the wigner seitz cell
    while(dr1r2r3(n)<-0.5)dr1r2r3(n)+=1.0;
    if(fabs(dr1r2r3(n)-0.5)<SMALLPOSITIONDEVIATION)dr1r2r3(n)=0.0; // if it is exactly at the edge - split atom and put half at +0.5 and -0.5 and sum --> results
    if(fabs(dr1r2r3(n)+0.5)<SMALLPOSITIONDEVIATION)dr1r2r3(n)=0.0; // --> in the same contribution as putting atom at 0. thus put it at origin for this purpose
    } 
    xyz=inputpars.cs.r*dr1r2r3; // transform back to dadbdc 
     dadbdc2ijk(rijk,xyz, inputpars.cs.abc); 
                      // transforms vector xyz given in terms of abc
                      // to ijk coordinate system (in summing dipole moments it is not necessary to consider
                     //  the subcell vector of the magnetic unit cell - this will cancel in summation anyway)
     physprop.Pel0+=rijk*(*inputpars.jjj[l]).charge; // sum dipolar moments
}physprop.Pel0/=inputpars.cs.pVol(); // divide by number of primitive cells in supercelland Vol of prim unitcell
// thus we have not Pel in units of |e|A/A^3=|e|/A^2. For SI we need it in C/m^2
// |e|=1.602189e-19 Coulomb
// A=1e-10 m
// i.e. |e|/A^2=1.602189e+1=16.02189 C/m^2
 physprop.Pel0*=16.02189;
  if(verbose==1){printf(".. calculating electrical Polarisation\n");}
}else
{physprop.Pel0=0;
 if(verbose==1)printf("...unit cell total charge=%g calculating electrical Polarisation does not make sense\n",inputpars.totalcharge);                              
}
// ***********************************************
if(NUM_THREADS>256){fprintf(stderr,"Error mcphas: too many threads required - change hardcode limit 256 in mcphas_htcalc.c line 69 and recompile\n");exit(EXIT_FAILURE);}
                  for (int ithread=0; ithread<NUM_THREADS; ithread++) 
                    tin[ithread] = new htcalc_input(0,ithread,&inputpars);
 
fin=NULL; if (argc>options+1){ini.xv=0;ini.yv=0;fin=fopen_errchk (argv[argc-1],"rb");}   //input from file
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
   {while (feof(fin)==0&&(0==inputline(fin,nn)||ini.calcTHfromnn(T,h,nn,inputpars.cs)==false)){;}  // if yes -> input them
    if (feof(fin)!=0) goto endproper;
    //x=nn[1];y=nn[2];T=nn[3];h=0;h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];
    // column 5 6 7 are the direction of h given by its components in i j k coordinates
    // --> replaced  above by      -------------^
   // ini.calcTHfromnn(T,h,nn,inputpars.cs): set external field and Temperature given nn as input from file with meaning defined by out1-7 in mcphas.ini
   // returns true if successful  (NormH NormE x y are not used)

   }
   else
   {
   //if parameters outside specified region then put them into it ...
    if (x<ini.xmin) x=ini.xmin;
    if (y<ini.ymin) y=ini.ymin;    
    ini.calcTHfromxy(T,h,x,y,inputpars.cs);
  } 
          
      physprop.x=x;physprop.y=y;
      physprop.T=T;
      physprop.H=h;

if(verbose==1){printf("Ha Hb Hc are components of magnetic field with respect to the Bravais lattice unit vectors ^a=a/|a|  ^b=b/|b| ^c=c/|c|\n");
 if (inputpars.cs.alpha()!=90||inputpars.cs.beta()!=90||inputpars.cs.gamma()!=90)
              {printf("Hi Hj Hk refer to components of magnetic field with respect to Euclidean Coordinates ijk defined by  j||b, k||(a x b) and i normal to k and j\n");}
}

// check if calculation results should and can be read (returns j=0)
int parsread=1;
j=1;if(readprefix[0]!='\0'){j=physprop.read(verbose,inputpars,readprefix,ini);parsread=j;}

// if not (j=1) then calculate physical properties at HT- point
if (j==1){j=htcalc(physprop.H,T,ini,inputpars,testqs,testspins,physprop);}
       switch (j)
       {case 0:
            //save physical properties of HT-point
	    //sta=(sta*ini.nofstapoints+physprop.save (verbose,filemode,j,inputpars))/(ini.nofstapoints+1);
          // 12.3.07 fancy calculation above substituted by normal summing of sta
          if(strcmp(ini.prefix,readprefix)!=0||prefix[0]=='\0'||parsread!=0)
           {sta+=physprop.save (verbose,filemode,j,ini,inputpars,ini.prefix);}
          else
           {ini.print_usrdefcols(stdout,x,y,T,physprop.H,inputpars.cs.abc,true);printf("\n");}
   	    ++ini.nofstapoints;
          if (sta>stamax){fprintf(stdout,"#! stamax=%g exceeded - exiting\n",stamax);goto endproper;}
	      break; 
	 case 1: goto endproper;
	      break;
         case 2: //ht calculation leads to no results- save dummy line
	         physprop.save (verbose,filemode,j,ini,inputpars,ini.prefix);
		 //sta+=1.0; // increment sta because within manifold of spincf no good solution could be found
     	      ++ini.noffailedpoints;
	      break;	 
	 default:  ;
	}
 if(fin!=NULL){   
y=ini.ymin-ini.ystep;} // this is to switch off xy loop if xy points are read from file

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
   double cpu_duration = (std::clock() - ini.startcputime) / (double)CLOCKS_PER_SEC;
   std::cout << "#! Finished in cputime=" << cpu_duration << " seconds [CPU Clock] " << std::endl;
   std::cout << "#!nofHTpoints=" << ini.nofstapoints << " H-T points in phasediagram successfully calculated" << std::endl;
   std::cout << "#!noffailedpoints=" << ini.noffailedpoints << " H-T points in phasediagram failed to converge " << std::endl;
   std::cout << "#!fecalc - free energy calculation was attempted noffecalccalls=" << ini.nofcalls << " times"  << std::endl;
   std::cout << "#!fecalc - free energy calculation was successful at noffecalcsuccess=" << ini.successrate << " times"  << std::endl;
   std::cout << "#!fecalc - free energy diverged maxnofloopsDIV=" << ini.nofmaxloopDIV << " times because maxnofloops was reached" << std::endl;
   std::cout << "#!fecalc - free energy diverged maxspinchangeDIV=" << ini.nofmaxspinchangeDIV << " times because maxspinchange was reached" << std::endl;

if(ini.nofstapoints>0)  { fprintf(stdout,"#! sta=%g\n",(ini.nofstapoints+ini.noffailedpoints)*sta/ini.nofstapoints);}
else { fprintf(stdout,"#! sta=1e10\n");}
#ifdef _THREADS
std::cout << "#! nofthreads= " << NUM_THREADS << " threads were used in parallel processing " << std::endl;
for (int ithread=0; ithread<NUM_THREADS; ithread++) delete tin[ithread];
#else
std::cout << "# mcphas was compiled without parallel processing option " << std::endl;
#endif
if(ini.ipx!=NULL)delete ini.ipx;
if(ini.ipy!=NULL)delete ini.ipy;
if(ini.ipz!=NULL)delete ini.ipz;
}
   fprintf(stderr,"**********************************************\n");
   fprintf(stderr,"          End of Program mcphas\n");
   fprintf(stderr," reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(stderr,"**********************************************\n");

return(0);
}

