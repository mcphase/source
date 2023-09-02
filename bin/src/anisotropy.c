/**************************************************************************
 *
 * anisotropy - program to calculate magnetic anisotropy
 *
 **************************************************************************/

#include<mcphas.h>

int verbose=1;
const char * filemode="w";
// for statistics 
int nofmaxloopDIV=0,nofmaxspinchangeDIV=0;
int successrate=0;
int nofcalls=0;
int isfull=0;

#include "myev.h"
#include "mcphas_htcalc.c"
#include "mcphas_fecalc.c"
#include "mcphas_physpropcalc.c"

// main program
int main (int argc, char **argv)
{  std::clock_t startcputime = std::clock();
  char sipffilename[MAXNOFCHARINLINE];
  int im,l,nofsteps;
  int do_sipffile=0;
  int nofthreads=1;
  double z,u;
  double T,H;
  Vector xv(1,3);
  Vector yv(1,3);
  Vector h(1,3);
  
fprintf(stderr,"#**************************************************************************\n");
fprintf(stderr,"#*\n");
fprintf(stderr,"#* anisotropy - program to calculate magnetic anisotropy \n");
fprintf(stderr,"#*\n");
fprintf(stderr,"#**************************************************************************\n\n");
int poly=0,P=6,doeps=0;//single crystal
// check command line
//T H xn yn zn nofsteps 
  if(argc<4){fprintf(stderr,"ERROR anisotropy: too few parameters\n");exit(EXIT_FAILURE);}
  T=strtod (argv[1], NULL); 
  H=strtod (argv[2], NULL); 
Vector direction(1,3);
if (strcmp(argv[3],"-p")==0){P=4;poly=1;//polycrystal
                    
                            }
else{
  direction(1)=strtod (argv[3], NULL); 
  direction(2)=strtod (argv[4], NULL); 
  direction(3)=strtod (argv[5], NULL); 
}  

nofsteps=(int)strtod (argv[P], NULL); 
  double dtheta=PI/abs(nofsteps)+0.000001;

if (argc>=8){if (strcmp(argv[P+1],"-doeps")==0){doeps=1;}
             if (strcmp(argv[P+1],"-r")==0)
                 {do_sipffile=1;strcpy(sipffilename,argv[P+2]);
                 }
            }
Vector x(1,3),r1(1,3),r2(1,3);x=0;x(1)=1;

if(poly==0){
 direction/=Norm(direction); // normalize
// now get r1 and r2 which are the basis vectors of the plane of rotation
  if (fabs(direction*x)>0.95){x=0;x(2)=1;}
 r1=x - (direction*x)*direction;
 r1/=Norm(r1);xproduct(r2,direction,r1);
 r2/=Norm(r2);
        }

 FILE * fout;fout=fopen_errchk("./results/anisotropy.out","w");
if(poly==0){
fprintf(fout,
"# output file of program: anisotropy @command\n"
"#! displayxtext=azimuth(deg) in plane perpendicular to [%g %g %g] direction\n"
"#! displayytext=M||H(mub)\n"
"#! displaytitle= Anisotropy plot az=0 corresponds to [%g %g %g]\n"
"# phi(deg) theta(deg) T[K] |H|[T] Hx[T] Hy[T] Hz[T] azimuth(deg) |M|[mb] Mx[mb] My[mb] Mz[mb] MparallelH[mb]\n",direction(1),direction(2),direction(3),r1(1),r1(2),r1(3));
}else{

fprintf(fout,
"# output file of program: anisotropy @command\n"
"#! displayxtext=theta(deg)\n"
"#! displayytext=M||H(mub)\n"
"#! displaytitle= Polycrystal Calculation Results - dependence on polar angle theta\n"
"# phi(deg) theta(deg) T[K] |H|[T] Hx[T] Hy[T] Hz[T] azimuth(deg) |M|[mb] Mx[mb] My[mb] Mz[mb] MparallelH[mb]\n");
}

if(do_sipffile){
 jjjpar jjj(0,0,0,sipffilename,argc-9);jjj.save_sipf("./results/_");
 int nofcomponents=argc-9;
 Vector Hxc(1,nofcomponents);
 for(int j=1;j<=nofcomponents;++j)Hxc=(int)strtod (argv[j+2+P], NULL); 
 Hxc=0;h=0;Vector m(1,3);
 jjj.Icalc_parameter_storage_init(Hxc,h,T);
if(poly==0){
 // loop different H 
 for(double az=0;az<2*PI-0.00001;az+=2*PI/nofsteps)
 {h=H*(cos(az)*r1+sin(az)*r2);
  double phi=PI/2;if(h(2)<0){phi=3*PI/2;}
  if(h(1)>0.001&&h(2)>=0){phi=atan(h(2)/h(1));}
  if(h(1)>0.001&&h(2)<0){phi=atan(h(2)/h(1))+2*PI;}
  if(h(1)<-0.001){phi=atan(h(2)/h(1))+PI;}
  double theta=acos(h(3)/H);
  jjj.mcalc(m,T,Hxc,h,jjj.Icalc_parstorage);
            //save physical properties of HT-point
    fprintf(fout,"%6.3f  %6.3f  %6.3f  %6.3f   %6.3f %6.3f %6.3f   %6.3f   %6.3f   %6.3f %6.3f %6.3f %6.3f\n",
           phi*180/PI,theta*180/PI,T,H,h(1),h(2),h(3),az*180/PI,Norm(m),m(1),m(2),m(3),m*h/Norm(h));  
 }      
       } // poly==0
else
 {// loop sphere


 int ct=0; double mpoly=0;
 for (double theta=dtheta;theta<PI;theta+=dtheta){ 
 double  dphi=dtheta*PI/4/sin(theta);
  for (double phi=0;phi<2*PI;phi+=dphi){
 h(1)=H*sin(theta)*cos(phi);
 h(2)=H*sin(theta)*sin(phi);
 h(3)=H*cos(theta);
 jjj.mcalc(m,T,Hxc,h,jjj.Icalc_parstorage);
            //save physical properties of HT-point
    fprintf(fout,"%6.3f  %6.3f  %6.3f  %6.3f   %6.3f %6.3f %6.3f   %6.3f   %6.3f   %6.3f %6.3f %6.3f %6.3f\n",
           phi*180/PI,theta*180/PI,T,H,h(1),h(2),h(3),theta*180/PI,Norm(m),m(1),m(2),m(3),m*h/Norm(h)); 
  ++ct;mpoly+=m*h/Norm(h);
 }} 
 mpoly/=ct;
 fprintf(stdout,"# T(K) H(Tesla) Mpolycrystal (mB) \n %6.3f %6.3f %6.3f \n",T,H,mpoly);
 fprintf(fout,"# T= %6.3f  K Hexternal= %6.3f Tesla Mpolycrystal=%6.3f \n",T,H,mpoly);

 }      // poly==0
}else{  // !do sipf
// as class par load  parameters from file
 if(verbose==1){printf("reading parameters from file mcphas.j\n");}
 char prefix [MAXNOFCHARINLINE];prefix[0]='\0';
 inipar ini("mcphas.ini",prefix);ini.doeps=doeps;
 par inputpars("./mcphas.j",verbose ); inputpars.save("./results/_mcphas.j",0); 
 nofthreads = ini.nofthreads;
  Vector Imax(1,inputpars.cs.nofatoms*inputpars.cs.nofcomponents);
  Vector Imom(1,inputpars.cs.nofcomponents);
  Vector h1(1,inputpars.cs.nofcomponents),h1ext(1,3);h1ext=0; 
  // here save single ion property files to results
  inputpars.save_sipfs("./results/_");
  //determine saturation momentum (used for scaling the plots, generation of qvectors)
  for(l=1;l<=inputpars.cs.nofatoms;++l){h1=0;(*inputpars.jjj[l]).Icalc_parameter_storage_init(h1,h1ext,T); // initialize eigenstate matrix
  for (im=1;im<=inputpars.cs.nofcomponents;++im){h1ext=0;h1=0;h1(im)=20*MU_B; //just put some high field
                            (*inputpars.jjj[l]).Icalc(Imom,T,h1,h1ext,z,u,(*inputpars.jjj[l]).Icalc_parstorage);
                            Imax(inputpars.cs.nofcomponents*(l-1)+im)=Imom(im);
                                              }
                                  }
 // load testspinconfigurations (nooftstspinconfigurations,init-file,sav-file)
   testspincf testspins (ini.maxnoftestspincf,"./mcphas.tst","./results/mcphas.phs",inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
   testspins.save("./results/_mcphas.tst","w");
   qvectors testqs (ini,inputpars,Imax,"./results/mcphas.qvc",verbose);
 // declare variable physprop (typa class physproperties)
   physproperties physprop(ini.nofspincorrs,ini.maxnofhkls,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
   	int nofstapoints=0,noffailedpoints=0,s=0;

if(poly==0){
 // loop different H /T points in phase diagram
 for(double az=0;az<2*PI-0.00001;az+=2*PI/nofsteps)
 {h=H*(cos(az)*r1+sin(az)*r2);
  double phi=PI/2;if(h(2)<0){phi=3*PI/2;}
  if(h(1)>0.001&&h(2)>=0){phi=atan(h(2)/h(1));}
  if(h(1)>0.001&&h(2)<0){phi=atan(h(2)/h(1))+2*PI;}
  if(h(1)<-0.001){phi=atan(h(2)/h(1))+PI;}
  double theta=acos(h(3)/H);
   // set field        
      physprop.T=T;
      physprop.H=h;
 //calculate physical properties at HT- point
   s=htcalc(physprop.H,T,ini,inputpars,testqs,testspins,physprop);
   if(s==1)break;
   //save physical properties of HT-point
   if(s==0)++nofstapoints;
   if(s==2)++noffailedpoints;
    fprintf(fout,"%6.3f  %6.3f  %6.3f  %6.3f   %6.3f %6.3f %6.3f   %6.3f   %6.3f   %6.3f %6.3f %6.3f %6.3f\n",
           phi*180/PI,theta*180/PI,T,H,h(1),h(2),h(3),az*180/PI,Norm(physprop.m),physprop.m(1),physprop.m(2),physprop.m(3),physprop.m*h/Norm(h));  
  } // H/T loop 
 }
else
{int ct=0; double mpoly=0;
 for (double theta=dtheta;theta<PI;theta+=dtheta){ 
 double  dphi=dtheta*PI/4/sin(theta);
  for (double phi=0;phi<2*PI;phi+=dphi){
 h(1)=H*sin(theta)*cos(phi);
 h(2)=H*sin(theta)*sin(phi);
 h(3)=H*cos(theta);
 // set field        
      physprop.T=T;
      physprop.H=h;
 //calculate physical properties at HT- point
   s=htcalc(physprop.H,T,ini,inputpars,testqs,testspins,physprop);
   if(s==1)break;
   //save physical properties of HT-point
   if(s==0)++nofstapoints;
   if(s==2)++noffailedpoints;
    fprintf(fout,"%6.3f  %6.3f  %6.3f  %6.3f   %6.3f %6.3f %6.3f   %6.3f   %6.3f   %6.3f %6.3f %6.3f %6.3f\n",
           phi*180/PI,theta*180/PI,T,H,h(1),h(2),h(3),theta*180/PI,Norm(physprop.m),physprop.m(1),physprop.m(2),physprop.m(3),physprop.m*h/Norm(h));  
   ++ct;mpoly+=physprop.m*h/Norm(h);
 }} 
 mpoly/=ct;
fprintf(stdout,"# T(K) H(Tesla) Mpolycrystal (mB) \n %6.3f %6.3f %6.3f \n",T,H,mpoly);
fprintf(fout,"# T= %6.3f  K Hexternal= %6.3f Tesla Mpolycrystal=%6.3f \n",T,H,mpoly);


}

   std::cout << "#!nofHTpoints=" << nofstapoints << "H-T points  successfully calculated" << std::endl;
   std::cout << "#!noffailedpoints=" << noffailedpoints << "H-T points in phasediagram failed to converge " << std::endl;
   std::cout << "#!fecalc - free energy calculation was attempted noffecalccalls=" << nofcalls << "times"  << std::endl;
   std::cout << "#!fecalc - free energy calculation was successful at noffecalcsuccess=" << successrate << "times"  << std::endl;
   std::cout << "#!fecalc - free energy diverged maxnofloopsDIV=" << nofmaxloopDIV << " times because maxnofloops was reached" << std::endl;
   std::cout << "#!fecalc - free energy diverged maxspinchangeDIV=" << nofmaxspinchangeDIV << " times because maxspinchange was reached" << std::endl;


} // do_sipffile




fclose(fout);
//  testspins.save(filemode);testqs.save(filemode);
   printf("#RESULTS saved in directory ./results/  - files: anisotropy.out\n");
   double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
   std::cout << "#! Finished in cputime=" << cpu_duration << " seconds [CPU Clock] " << std::endl;
   
#ifdef _THREADS
std::cout << "#! nofthreads= " << nofthreads << " threads were used in parallel processing " << std::endl;
for (int ithread=0; ithread<nofthreads; ithread++) delete tin[ithread];
#else
std::cout << "# anisotropy was compiled without parallel processing option " << std::endl;
#endif

 
   fprintf(stderr,"#**********************************************\n");
   fprintf(stderr,"#          End of Program anisotropy\n");
   fprintf(stderr,"#**********************************************\n");

return(0);
}

