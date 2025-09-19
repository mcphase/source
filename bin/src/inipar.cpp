// methods for class inipar to store parameters for mcphasit 
#include "inipar.hpp"
#include "cryststruct.hpp"
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


int normalizedadbdc(Vector & dadbdc,double n,cryststruct & cs)
   {if(Norm(dadbdc)>0.00001){ // normalize Vector dadbdc (da da dc are components with respect to
                              // Bravais lattice a b c) to length n Angstroem
    Vector Hijk(1,3);
    Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1;// !!!! a b c are unit vectors along Bravais Lattice vectors !!!
                     abc(4)=cs.alpha(); abc(5)=cs.beta(); abc(6)=cs.gamma();
    dadbdc2ijk(Hijk,dadbdc,abc);
    Hijk*=n/Norm(Hijk);
    ijk2dadbdc(dadbdc,Hijk,abc);
    return true;      }
    else
   {return false;}
   }
int normalize(Vector & x,double n)
   {
if(Norm(x)>0.00001){ // normalize Vector x to unit length
     x*=n/Norm(x);
    return true;      }
    else
   {return false;}
   }

 // *************************************************************************
 // ************************ inipar *************************************
 // *************************************************************************
 // class of initial parameters for program mcphas

void inipar::errexit() // type info and error exit 
{     printf (" \n%s \n",MCPHASVERSION);
printf (" use as: mcphas [options] [file]\n");
printf (" [file] ... input file  with sets of x y T H Hi Hj Hk points \n");
printf (" (format as output file mcphas.xyt)\n\n");
printf (" Options: -h     print this help screen\n");
printf ("          -stamax 14  ... end mcphas if standard deviation exceeds 14\n");
printf ("          -a     append output files (do not overwrite) \n");
printf ("          -cd     add classical dipole  interaction using Ewald summation, \n");
printf ("                  Bowden J.Phys.C:solid state phys. 14(1981) L827  \n");
printf ("          -doeps refine strain epsilon selfconsistently using elastic,magnetoelastic constants \n");
printf ("                 read from mcphas.j and mcphas.djdx mcphas.djdy and mcphas.djdz and optional djdeps1-6\n");
printf ("          -linepscf with -doeps use zero strain single ion Hamiltonian for every mean field iteration\n");
printf ("          -linepsjj with -doeps use zero strain two ion interaction Hamiltonian for every mean field iteration\n");
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

void inipar::finish_mcphas(int nofqs,int nofspincf)
{
printf("RESULTS saved in directory ./results/  - files:\n");
   printf("#RESULTS saved in directory ./results/  - files:\n");
   printf("#  %smcphas.fum  - total magnetic moment, energy at different T,H\n",prefix);
   printf("#  %smcphas.sps  - stable configurations at different T,H\n",prefix);
   printf("#  %smcphas.mf   - mean fields at different T,H\n",prefix);
   printf("#  %smcphas.hkl  - strong magnetic satellites, neutron diffraction intensity\n",prefix);
   printf("#  %smcphas*.hkl - strong magnetic satellites, Fourier Comp.of moment in * dir\n",prefix);
   printf("#  %smcphas*.j*  - JJ correlation functions (for exchange magnetostriction)\n",prefix);
   printf("#  %smcphas.xyt  - phasediagram (stable conf.nr, angular and multipolar moments)\n",prefix);
   printf("#!  %smcphas.qvc  - ...corresponding table of all nqvc=%i qvector generated test configs\n",prefix,nofqs);
   printf("#!  %smcphas.phs  - ...corresponding table of all ntst=%i configurations (except qvecs)\n",prefix,nofspincf);
   printf("#  _%smcphas.*   - parameters read from input parameter files (.tst,.ini,.j)\n",prefix);
   printf("#  ...         - and a copy of the single ion parameter files used.\n\n");
   double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
   std::cout << "#! Finished in cputime=" << cpu_duration << " seconds [CPU Clock] " << std::endl;
   std::cout << "#!nofHTpoints=" << nofstapoints << " H-T points in phasediagram successfully calculated" << std::endl;
   std::cout << "#!nofreppoints="<< nofreppoints << " points repeated ( nofconvrep=" << nofconvrep << " of which converged after repetition)" << std::endl;
   std::cout << "#!noffailedpoints=" << noffailedpoints << " H-T points in phasediagram failed to converge " << std::endl;
   std::cout << "#!fecalc - free energy calculation was attempted noffecalccalls=" << nofcalls << " times"  << std::endl;
   std::cout << "#!fecalc - free energy calculation was successful at noffecalcsuccess=" << successrate << " times"  << std::endl;
   std::cout << "#!fecalc - free energy diverged maxnofloopsDIV=" << nofmaxloopDIV << " times because maxnofloops was reached" << std::endl;
   std::cout << "#!fecalc - free energy diverged maxspinchangeDIV=" << nofmaxspinchangeDIV << " times because maxspinchange was reached" << std::endl;

if(nofstapoints>0)  { fprintf(stdout,"#! sta=%g\n",(nofstapoints+noffailedpoints)*sta/nofstapoints);}
else { fprintf(stdout,"#! sta=1e10\n");}
#ifdef _THREADS
std::cout << "#! nofthreads= " << nofthreads << " threads were used in parallel processing " << std::endl;
#else
std::cout << "# mcphas was compiled without parallel processing option " << std::endl;
#endif
if(ipx!=NULL)delete ipx;    
if(ipy!=NULL)delete ipy;
if(ipz!=NULL)delete ipz;
if(ipeps1!=NULL)delete ipeps1;    
if(ipeps2!=NULL)delete ipeps2;    
if(ipeps3!=NULL)delete ipeps3;    
if(ipeps4!=NULL)delete ipeps4;    
if(ipeps5!=NULL)delete ipeps5;    
if(ipeps6!=NULL)delete ipeps6;    
}

int usrdefcols[]={7, 1,2,3,4,5,6,7}; // user defined output columns (first number is number of usr def output columns)
                                             // in files mcdisp.qei,qex,qom,dsigma,dsigma.tot
int colcod[]=    {-1,19,20,0,21,1,2,3}; // field to store code for assigning type of data to columns of output,
                                           // set default values here (see list below for different types)
                                           // using the out5 out6  ... commands in mcdisp.par these codes can be modified

#define COLHEADDIM 22	
// different output data for columns 1-7
const char * colhead []= {  "T [K]", //      0                                                 
                            "Ha [T]", //      1                                                  
                            "Hb [T]", //      2                                                  
                            "Hc [T]", //      3 
                            "Hi [T]", //      4                                                 
                            "Hj [T]", //      5                                                 
                            "Hk [T]", //      6
                            "Ea [kV/mm]",  //     7      
                            "Eb [kV/mm]",  //     8      
                            "Ec [kV/mm]",  //     9      
                            "Ei [kV/mm]",  //    10      
                            "Ej [kV/mm]",  //    11      
                            "Ek [kV/mm]",  //    12     
                            "s1 [GPa]",  //    13    
                            "s2 [GPa]",  //    14      
                            "s3 [GPa]",  //    15      
                            "s4 [GPa]",  //    16      
                            "s5 [GPa]",  //    17      
                            "s6 [GPa]",  //    18 
                            "x",  //    19 
                            "y", //     20 
                            "|H| [T]",   //    21
                            "|E| [T]"   //    22
                               };
 bool inipar::defaultcolcode(int col,int colcode) // resets default columns if not set by user (outcolset==true)
{ bool ret=false;                                            // returns true if reset has been successful
 if(!outcolset)for(int i=1;i<=usrdefcols[0];++i)if(usrdefcols[i]==col)
 {colcod[i]=colcode;ret=true;}
return ret;
}
// calculate the value of different output data for user defined columns ...
double inipar::calccolvalue(int i,float & x, float & y,double& T,Vector & Hext,Vector & abc)
{double xx=x,yy=y;
 Vector Habc(1,3);Vector abcu(1,6);abcu=abc;abcu(1)=1;abcu(2)=1;abcu(3)=1;Vector v(1,3);v=Hext(1,3);
          ijk2dadbdc(Habc,v,abcu);
 Vector Eabc(1,3);v=Hext(4,6,-3);
          ijk2dadbdc(Eabc,v,abcu);
Vector Hijk(1,3);Hijk=Hext(1,3);double NormH=Norm(Hijk);
Vector Eijk(1,3);Eijk=Hext(4,6,-3);double NormE=Norm(Eijk);
double ret=(*colvaluepointer(i,xx,yy,T,Hext,Habc,Eabc,NormH,NormE));
return ret;
}

double * inipar::colvaluepointer(int i,double & x, double & y,double& T,Vector & Hext,Vector & Habc,
                 Vector & Eabc,double & NormH, double & NormE)
{      switch (i) {
case 0:  return &T;break;
case 1:  
case 2:  
case 3:  return &Habc(i);
         break;
case 4:  return &Hext(1);break;
case 5:  return &Hext(2);break;
case 6:  return &Hext(3);break;
case 7:  
case 8:  
case 9:  return &Eabc(i-6);
         break;
case 10:  return &Hext(4);break;
case 11:  return &Hext(5);break;
case 12:  return &Hext(6);break;
case 13:  return &Hext(7);break;
case 14:  return &Hext(8);break;
case 15:  return &Hext(9);break;
case 16:  return &Hext(10);break;
case 17:  return &Hext(11);break;
case 18:  return &Hext(12);break;
case 19:  return &x;break;
case 20:  return &y;break;
case 21:  return &NormH ;break;
case 22:  return &NormE;break;
default: fprintf(stderr,"Error mcphas: unknown column code %i\n",i);exit(EXIT_FAILURE);
                    }

return 0;
}



// print user defined column headers
void inipar::print_usrdefcolhead(FILE *fout,char * str)
{fprintf(fout,"#");
 int i;char header [MAXNOFCHARINLINE];header[0]='\0';
 for(i=1;i<=usrdefcols[0];++i)
  {for(int j=0;j<(int)strlen(colhead[colcod[i]]);++j)if(colhead[colcod[i]][j]!=' ')
    snprintf(header+strlen(header),MAXNOFCHARINLINE-strlen(header),"%c",colhead[colcod[i]][j]); 
  snprintf(header+strlen(header),MAXNOFCHARINLINE-strlen(header)," ");
  }
 snprintf(header+strlen(header),MAXNOFCHARINLINE-strlen(header),"%s ",str);
 print_col_numbers(fout,header);
 fprintf(fout,"\n#%s\n",header);
}


// set external field and Temperature
void inipar::calcTHfromxy(double & T,Vector & h,double x, double y,cryststruct & cs)
{    T=zero(0)+x*xv(0)+y*yv(0);
    // this means from input we take the vector xHa xHb xHc, interpret it as fractional 
// coordinates in terms of unit !! vectors along the Bravais lattice vectors
// and normalize this vector to 1 and then multiply it by x
// and then add Ha0 Hb0 Hc0  and store this in physprop.H, further
//  we  use dadbdc2ijk (again with Bravais lattice of unit length)
// to transform Ha Hb Hc to Euclidean ijk coordinates 
// and this will be used in the calculation as external field in Tesla
// -->  input (xHa xHb xHc) and (Ha0 Hb0 Hc0) are vectors 
// given in terms of components with respect to unit vectors along the Bravais lattice a, b, c.
// For the external magnetic field unit is Tesla.
// Therefore a tooltip text will be:
// xHa: Magnetic Field component with respect to Bravais lattice unit vector ^a=a/|a| (normalised to 1 Tesla)
// Ha0: Offset - Magnetic Field component with respect to Bravais lattice unit vector ^a=a/|a| (normalised to 1 Tesla)
    Vector abc(1,6),v(1,3),v1(1,3),v2(1,3),zv(1,3); abc(1)=1; abc(2)=1; abc(3)=1; // trick to get Habc as components along a,b,c
                  abc(4)=cs.alpha(); abc(5)=cs.beta(); abc(6)=cs.gamma();
    v=xv(1,3);normalizedadbdc(v,1.0,cs);// take care that vector xHa xHb Xhc has unit length 1 Tesla
    dadbdc2ijk(v1,v,abc); // transform Habc to ijk coordinates ... this is H
    v=xv(4,6,-3);
    normalize(v,1.0); // take care of Hijk to be of unit length 1 Tesla
    v1=v1+v;

    v=yv(1,3);normalizedadbdc(v,1.0,cs);// take care that vector yHa yHb yHc has unit length 1 Tesla
    dadbdc2ijk(v2,v,abc); // transform Habc to ijk coordinates ... this is H
    v=yv(4,6,-3);normalize(v,1.0); // take care of Hijk to be of unit length 1 Tesla
    v2=v2+v;

    v=zero(1,3);
    dadbdc2ijk(zv,v,abc); // transform Habc to ijk coordinates 
    v=zero(4,6,-3); 
    zv=zv+v;

     h=0;
     h(1)=zv(1)+x*v1(1)+y*v2(1);
     h(2)=zv(2)+x*v1(2)+y*v2(2);
     h(3)=zv(3)+x*v1(3)+y*v2(3);
   
    // similar for the E-field

    v=xv(7,9,-6);normalizedadbdc(v,1.0,cs); // take care that vector xEa xEb XEc has unit length 1 Tesla
    dadbdc2ijk(v1,v,abc); // transform Eabc to ijk coordinates ... this is E
    v=xv(10,12,-9);normalize(v,1.0); // take care of Eijk to be of unit length 1 Tesla
    v1=v1+v;
     
    v=yv(7,9,-6);normalizedadbdc(v,1.0,cs);// take care that vector yEa yEb yEc has unit length 1 Tesla
    dadbdc2ijk(v2,v,abc); // transform Eabc to ijk coordinates ... this is E
    v=yv(10,12,-9);normalize(v,1.0); // take care of Eijk to be of unit length 1 Tesla
    v2=v2+v;
    
    v=zero(7,9,-6);
    dadbdc2ijk(zv,v,abc); // transform Habc to ijk coordinates 
    v=zero(10,12,-9); 
    zv=zv+v;
     h(4)=zv(1)+x*v1(1)+y*v2(1);
     h(5)=zv(2)+x*v1(2)+y*v2(2);
     h(6)=zv(3)+x*v1(3)+y*v2(3);
   
   // now stress tensor
   Vector xs(1,6),ys(1,6);
   xs=xv(13,18,-12);normalize(xs,1.0);
   ys=yv(13,18,-12);normalize(ys,1.0);
   for(int i=1;i<=6;++i){
     h(i+6)=zero(i+12)+x*xs(i)+y*ys(i);}


}

// set external field and Temperature given nn as input from file with meaning defined by out1-7 in mcphas.ini
 // returns true if successful (NormH NormE x y are not used)
bool inipar::calcTHfromnn(double & T,Vector & Hext,float * nn,cryststruct &cs)
{int maxcol=0;for(int i=1;i<=usrdefcols[0];++i)if(usrdefcols[i]>maxcol)maxcol=usrdefcols[i];
 if(nn[0]<maxcol)return false; // array too small, not enough parameters in line
 double xx=0,yy=0,NormH,NormE;T=0;Hext=0;
 Vector Habc(1,3),Eabc(1,3);Habc=0;Eabc=0;
for(int i=1;i<=usrdefcols[0];++i)
{ (*colvaluepointer(colcod[usrdefcols[i]],xx,yy,T,Hext,Habc,Eabc,NormH,NormE))=nn[usrdefcols[i]];
}
// if Habc or Eabc are given - add these to Hext
Vector abc(1,6),v(1,3); abc(1)=1; abc(2)=1; abc(3)=1; 
                  abc(4)=cs.alpha(); abc(5)=cs.beta(); abc(6)=cs.gamma();
dadbdc2ijk(v,Habc,abc); 
Hext(1)+=v(1);
Hext(2)+=v(2);
Hext(3)+=v(3);
dadbdc2ijk(v,Eabc,abc); 
Hext(4)+=v(1);
Hext(5)+=v(2);
Hext(6)+=v(3);

 return true;
}

 // given T and Hext check if in array nn[0-7] the values are in accordance with T and Hext
 // if yes, returns true ... 
bool inipar::checkTH(float * nn,double & T,Vector & Hext,Vector & abc)
{int maxcol=0;for(int i=1;i<=usrdefcols[0];++i)if(usrdefcols[i]>maxcol)maxcol=usrdefcols[i];
 if(nn[0]<maxcol)return false; // array too small
 double d;float x=0,y=0; // do not use x and y
 for(int i=1;i<=usrdefcols[0];++i)
 { // different output data for user defined columns ...
  switch(colcod[i])
  {case 19: case 20:  d=0; break; // do not use x,y
   default: d=calccolvalue(colcod[i],x,y, T,Hext, abc)-nn[i];
  }
//  printf("d=%g i=%i nn=%g |",d,i,nn[i]);
  if(fabs(d)>SMALL_FIELD)return false;
 }
 return true;
}

// print user defined column codes variables out1 -- out7 to fout
void inipar::print_usrdefcolcodes(FILE *fout)
{fprintf(fout,"#!");
 for(int i=1;i<=usrdefcols[0];++i)
 fprintf(fout,"out%i=%s ",usrdefcols[i],colhead[colcod[i]]);
}
// print user defined columns
void inipar::print_usrdefcols(FILE *fout,float & x, float & y,double& T,Vector & Hext,Vector & abc,bool withtext)
{bool c[COLHEADDIM+1];for(int i=0;i<=COLHEADDIM;++i)c[i]=false;
 for(int i=1;i<=usrdefcols[0];++i)
 { double val=calccolvalue(colcod[i],x,y,T,Hext,abc);
   if(withtext)fprintf(fout,"%s=%4.4g ",colhead[colcod[i]],myround(val));
   else fprintf(fout,"%*s%4.4g ",(int)(strlen(colhead[colcod[i]])-8 < 0 ? 0 :strlen(colhead[colcod[i]])-8 ),"",myround(val));
   c[colcod[i]]=true;
 }
if (!c[0]){fprintf(stderr,"#Error: Temperature T not stored  - please change settings out out* in mcphas.ini\n");exit(EXIT_FAILURE); }
for(int i=1;i<=HEXT_DIMENSION;++i)
{if(fabs(Hext(i))>SMALL_FIELD)
 {switch(i)
  {case 1: case 2: case 3:{ bool cc=c[1]|c[2]|c[3]|c[4]|c[5]|c[6]|c[21];
   if(!cc){fprintf(stderr,"#Warning: External Magnetic Field H nonzero but not stored in output files  - please change settings out out* in mcphas.ini\n");exit(EXIT_FAILURE); }
                           } break;
   case 4: case 5: case 6: { bool cc=c[7]|c[8]|c[9]|c[10]|c[11]|c[12]|c[22];
   if(!cc){fprintf(stderr,"#Warning: External Electric Field E nonzero but not stored in output files  - please change settings out out* in mcphas.ini\n");exit(EXIT_FAILURE); }
                           }break;
   default: if(!c[i+6]){fprintf(stderr,"#Warning: stress s%i nonzero but not stored in output files  - please change settings out out* in mcphas.ini\n",i-6);exit(EXIT_FAILURE); }
  }
 }
}

}
 


void inipar::time_estimate_until_end(double x, double y)
{// estimate time until end 
    int nofpoints=nofstapoints+noffailedpoints+1; // add 1 to avoid zero
    int nofysteps=(int)((ymax-ymin)/ystep); if(nofysteps==0){nofysteps=1;}
    int pointstodo=nofysteps*int((xmax-x)/xstep)+int((ymax-y)/ystep);
    print_time_estimate_until_end((double)pointstodo/(double)nofpoints);
    //printf("%i  %i HTpoints to do.",nofysteps,pointstodo);

 
}

//load parameters from file
int inipar::load ()
{int n=-1;char ** lp; lp=NULL;
return load(n,lp);
}

// findnewmatch = true: check if new match (which is not listed in the n prefixes of lofpref) can be found in instr using
//                      char prefix (which may contain a wildcard '*'), if it can be found increase n
//                      and put new match prefix in lofpref and put findnewmatch to false and return
//                      extract_with_prefix for the new match
// findnewmatch = false: return extract_with_prefix for the prefix=lofpref[n]
int inipar::extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,int & var)
{int ret; double val;
ret=extract_match(findnewmatch,  n,lofpref ,instr,pref,  parameter, val);
if(ret==0)var=(int)val;
return ret;
}
int inipar::extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,float & var)
{int ret; double val;
ret=extract_match(findnewmatch,  n,lofpref ,instr,pref,  parameter, val);
if(ret==0)var=(float)val;
return ret;
}

int inipar::extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,double & var)
{
 if(findnewmatch==true)
 {if(pref[0]!='\0')
  {char prefvar[MAXNOFCHARINLINE];char * s,*p;
  snprintf(prefvar,MAXNOFCHARINLINE,"%s%s",pref,parameter);
  s=wstrstr(instr,prefvar);
  if(s!=NULL)
  { // ok there seems to be a new match - see if it is already in the list
   snprintf(prefvar,MAXNOFCHARINLINE,"%s",parameter);
   p=wstrstr(instr,prefvar);
   int i=0;for(char * t=s;t<p;++t){prefvar[i]=*t;++i;}prefvar[i]='\0'; // put the prefix to prefvar
    i=0;bool mm=false;while(i<n&&mm==false){mm=match(prefvar,lofpref[i]);
    ++i;}
    if(i==n&&mm==false){lofpref[i]=new char [strlen(prefvar)+2];snprintf(lofpref[i],strlen(prefvar)+1,"%s",prefvar);
                          snprintf(prefix,MAXNOFCHARINLINE,"%s",prefvar);
                     ++n; findnewmatch=false;printf("increase lofpref %s\n",lofpref[i]);
            }
   }// no new match -> extract without prefix
   else
   {
    return extract(instr,parameter,var);
   }
  }
 else
 {if(n==0){int i=0;lofpref[i]=new char [1];lofpref[i][0]='\0';
                          prefix[0]='\0';
                     ++n; findnewmatch=false;}
 }
 }
 if(n>0)
 {
 if(findnewmatch==true) return extract(instr,parameter,var);  // if we still have to find new match: return without prefix 
  else return extract_with_prefix(instr,lofpref[n-1],parameter,var); // if a new match has been found: return parameters with this new prefix
 }
 else
 {return extract_with_prefix(instr,pref,parameter,var);
 }
}

//load parameters from file
int inipar::load (int & nofinis,char**lofpref)
{ FILE *fin;outcolset=false;
  bool findnewmatch=true;if(nofinis==-1){findnewmatch=false;}

  char instr[MAXNOFCHARINLINE];
  char somestring[MAXNOFCHARINLINE];
  errno = 0;startcputime= std::clock();
  fin = fopen(savfilename, "rb");
  if (fin==NULL) return 1;
  xv=0;yv=0;xmin=1;xmax=0;ymin=1;ymax=0;xstep=0;ystep=0;zero=0;
  qmin(1)=1;qmin(2)=1;qmin(3)=1;qmax=0;deltaq=0;maxqperiod=0;maxnofspins=0;nofrndtries=0;nofMCsteps=0;
  minnr1=0;
  minnr2=0;
  minnr3=0;
  maxnofmfloops=-1;maxstamf=0;bigstep=0;maxspinchange=0;nofthreads=0;repeat=0;
  nofspincorrs=0;maxnofhkls=0;maxQ=0;maxnoftestspincf=1000;
  
  while (fgets(instr,MAXNOFCHARINLINE,fin)!=NULL)
  {if(!(instr[strspn(instr," \t")]=='#'&&instr[strspn(instr," \t#")]!='!')&&instr[strspn(instr," \t")]!='[') // comment lines headed by # or [ are ignored in mcphas.ini
   {extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"exit",exit_mcphas);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"pause",pause_mcphas);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"displayall",displayall);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"logfevsQ",logfevsQ); 
     
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xT",xv[0]);      
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xHa",xv[1]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xHb",xv[2]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xHc",xv[3]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xHi",xv[4]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xHj",xv[5]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xHk",xv[6]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xEa",xv[7]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xEb",xv[8]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xEc",xv[9]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xEi",xv[10]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xEj",xv[11]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xEk",xv[12]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xs1",xv[13]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xs2",xv[14]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xs3",xv[15]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xs4",xv[16]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xs5",xv[17]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xs6",xv[18]);

    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xmin",xmin);  
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xmax",xmax);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"xstep",xstep);
   
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yT",yv[0]);      
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yHa",yv[1]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yHb",yv[2]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yHc",yv[3]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yHi",yv[4]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yHj",yv[5]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yHk",yv[6]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yEa",yv[7]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yEb",yv[8]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yEc",yv[9]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yEi",yv[10]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yEj",yv[11]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"yEk",yv[12]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ys1",yv[13]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ys2",yv[14]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ys3",yv[15]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ys4",yv[16]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ys5",yv[17]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ys6",yv[18]);
    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ymin",ymin); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ymax",ymax);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ystep",ystep);   
 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"T0",zero[0]);      
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Ha0",zero[1]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Hb0",zero[2]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Hc0",zero[3]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Hi0",zero[4]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Hj0",zero[5]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Hk0",zero[6]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Ea0",zero[7]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Eb0",zero[8]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Ec0",zero[9]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Ei0",zero[10]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Ej0",zero[11]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Ek0",zero[12]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"s10",zero[13]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"s20",zero[14]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"s30",zero[15]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"s40",zero[16]);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"s50",zero[17]);    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"s60",zero[18]);

    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hmin",qmin[1]); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"kmin",qmin[2]); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"lmin",qmin[3]); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hmax",qmax[1]); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"kmax",qmax[2]); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"lmax",qmax[3]); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltah",deltaq[1]); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltak",deltaq[2]); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltal",deltaq[3]); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"maxqperiod",maxqperiod);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"minnr1",minnr1);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"minnr2",minnr2);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"minnr3",minnr3);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"maxnofspins",maxnofspins);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"nofrndtries",nofrndtries);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"nofMCsteps",nofMCsteps);
    
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"maxnofmfloops",maxnofmfloops);
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"maxstamf",maxstamf); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"bigstep",bigstep); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"repeat",repeat); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"maxspinchange",maxspinchange); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"maxnoftestspincf",maxnoftestspincf);

    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"nofthreads",nofthreads);

    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"nofspincorrs",nofspincorrs); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"maxnofhkls",maxnofhkls); 
    extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"maxQ",maxQ); 

       for(int j=1;j<=usrdefcols[0];++j) // extract user defined output columns
     {snprintf(somestring,MAXNOFCHARINLINE,"out%i",usrdefcols[j]);
      if(0==extract_match( findnewmatch,nofinis,lofpref,instr,prefix, somestring,colcod[usrdefcols[j]]))outcolset=true;
     }

    }
   }
  fclose (fin);
 for(int i=1;i<=usrdefcols[0];++i){if(colcod[i]>COLHEADDIM)
 {fprintf(stderr,"Error reading mcphas.ini - out%i = %i > %i not possible !\n",i,colcod[i],COLHEADDIM);exit(EXIT_FAILURE);}
 }

//if(nofinis>0)printf("prefix=%s Ha0=%g Hb0=%g Hc0=%g\n",lofpref[nofinis-1],zero(1),zero(2),zero(3));

  if (Norm(xv)==0){fprintf(stderr,"ERROR reading xT xHa xHb xHc=0\n");return 1;}
  if (Norm(yv)==0){fprintf(stderr,"ERROR reading yT yHa yHb yHc=0\n");return 1;}
  if (xmin>xmax){fprintf(stderr,"ERROR reading xmin=%g > xmax=%g \n",xmin,xmax);return 1;}
  if (ymin>ymax){fprintf(stderr,"ERROR reading ymin=%g > ymax=%g \n",ymin,ymax);return 1;}
  if (xstep==0){fprintf(stderr,"Warning reading xstep: xstep=0\n");}
  if (ystep==0){fprintf(stderr,"Warning reading ystep: ystep=0\n");}

  if(qmin(1)>qmax(1)){fprintf(stderr,"ERROR reading hmin=%g >  hmax=%g\n",qmin(1),qmax(1));return 1;}
  if(qmin(2)>qmax(2)){fprintf(stderr,"ERROR reading kmin=%g >  kmax=%g\n",qmin(2),qmax(2));return 1;}
  if(qmin(3)>qmax(3)){fprintf(stderr,"ERROR reading lmin=%g >  lmax=%g\n",qmin(3),qmax(3));return 1;}
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
  
  if(maxnofmfloops<0){fprintf(stderr,"Error reading maxnofmfloops<0\n");return 1;}
  if(nofrndtries<0){fprintf(stderr,"Error nofrndtries<0 \n");return 1;}
  if(nofMCsteps<0){fprintf(stderr,"Error nofMCsteps<0 \n");return 1;}
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
  print(fout);
  fclose(fout);
}

bool inipar::checkpr(FILE* fout,const  char * var,int val,int masterval)
{
if(masterval!=val){fprintf(fout,"%s%s=%i\n",prefix,var,val);return true;}
else {return false;}
}
bool inipar::checkpr(FILE* fout, const char * var,double val,double masterval)
{
if(masterval!=val){fprintf(fout,"%s%s=%g\n",prefix,var,val);return true;}
else {return false;}
}


// prints parameters including prefix - printout only if prefixed parameters is not equal to p.parameter
void inipar::print_with_prefix(FILE * fout, inipar p)
{checkpr(fout, "exit",exit_mcphas,p.exit_mcphas);
 checkpr(fout, "pause",pause_mcphas,p.pause_mcphas);
checkpr(fout, "displayall",displayall,p.displayall);
checkpr(fout, "logfevsQ",logfevsQ,p.logfevsQ);
checkpr(fout, "xT",xv(0),p.xv(0));
checkpr(fout, "xHa",xv(1),p.xv(1));
checkpr(fout, "xHb",xv(2),p.xv(2));
checkpr(fout, "xHc",xv(3),p.xv(3));
checkpr(fout, "xmin",xmin,p.xmin);
checkpr(fout, "xmax",xmax,p.xmax);
checkpr(fout, "xstep",xstep,p.xstep);
if(xv(4)!=0)checkpr(fout, "xHi",xv(4),p.xv(4));
if(xv(5)!=0)checkpr(fout, "xHj",xv(5),p.xv(5));
if(xv(6)!=0)checkpr(fout, "xHk",xv(6),p.xv(6));
if(xv(7)!=0)checkpr(fout, "xEa",xv(7),p.xv(7));
if(xv(8)!=0)checkpr(fout, "xEb",xv(8),p.xv(8));
if(xv(9)!=0)checkpr(fout, "xEc",xv(9),p.xv(9));
if(xv(10)!=0)checkpr(fout, "xEi",xv(10),p.xv(10));
if(xv(11)!=0)checkpr(fout, "xEj",xv(11),p.xv(11));
if(xv(12)!=0)checkpr(fout, "xEk",xv(12),p.xv(12));
if(xv(13)!=0)checkpr(fout, "xs1",xv(13),p.xv(13));
if(xv(14)!=0)checkpr(fout, "xs2",xv(14),p.xv(14));
if(xv(15)!=0)checkpr(fout, "xs3",xv(15),p.xv(15));
if(xv(16)!=0)checkpr(fout, "xs4",xv(16),p.xv(16));
if(xv(17)!=0)checkpr(fout, "xs5",xv(17),p.xv(17));
if(xv(18)!=0)checkpr(fout, "xs6",xv(18),p.xv(18));

checkpr(fout, "yT",yv(0),p.yv(0));
checkpr(fout, "yHa",yv(1),p.yv(1));
checkpr(fout, "yHb",yv(2),p.yv(2));
checkpr(fout, "yHc",yv(3),p.yv(3));
checkpr(fout, "ymin",ymin,p.ymin);
checkpr(fout, "ymax",ymax,p.ymax);
checkpr(fout, "ystep",ystep,p.ystep);
if(yv(4)!=0)checkpr(fout, "yHi",yv(4),p.yv(4));
if(yv(5)!=0)checkpr(fout, "yHj",yv(5),p.yv(5));
if(yv(6)!=0)checkpr(fout, "yHk",yv(6),p.yv(6));
if(yv(7)!=0)checkpr(fout, "yEa",yv(7),p.yv(7));
if(yv(8)!=0)checkpr(fout, "yEb",yv(8),p.yv(8));
if(yv(9)!=0)checkpr(fout, "yEc",yv(9),p.yv(9));
if(yv(10)!=0)checkpr(fout, "yEi",yv(10),p.yv(10));
if(yv(11)!=0)checkpr(fout, "yEj",yv(11),p.yv(11));
if(yv(12)!=0)checkpr(fout, "yEk",yv(12),p.yv(12));
if(yv(13)!=0)checkpr(fout, "ys1",yv(13),p.yv(13));
if(yv(14)!=0)checkpr(fout, "ys2",yv(14),p.yv(14));
if(yv(15)!=0)checkpr(fout, "ys3",yv(15),p.yv(15));
if(yv(16)!=0)checkpr(fout, "ys4",yv(16),p.yv(16));
if(yv(17)!=0)checkpr(fout, "ys5",yv(17),p.yv(17));
if(yv(18)!=0)checkpr(fout, "ys6",yv(18),p.yv(18));


checkpr(fout, "T0",zero(0),p.zero(0));
checkpr(fout, "Ha0",zero(1),p.zero(1));
checkpr(fout, "Hb0",zero(2),p.zero(2));
checkpr(fout, "Hc0",zero(3),p.zero(3));
if(zero(4)!=0)checkpr(fout, "Hi0",zero(4),p.zero(4));
if(zero(5)!=0)checkpr(fout, "Hj0",zero(5),p.zero(5));
if(zero(6)!=0)checkpr(fout, "Hk0",zero(6),p.zero(6));
if(zero(7)!=0)checkpr(fout, "Ea0",zero(7),p.zero(7));
if(zero(8)!=0)checkpr(fout, "Eb0",zero(8),p.zero(8));
if(zero(9)!=0)checkpr(fout, "Ec0",zero(9),p.zero(9));
if(zero(10)!=0)checkpr(fout, "Ei0",zero(10),p.zero(10));
if(zero(11)!=0)checkpr(fout, "Ej0",zero(11),p.zero(11));
if(zero(12)!=0)checkpr(fout, "Ek0",zero(12),p.zero(12));
if(zero(13)!=0)checkpr(fout, "s10",zero(13),p.zero(13));
if(zero(14)!=0)checkpr(fout, "s20",zero(14),p.zero(14));
if(zero(15)!=0)checkpr(fout, "s30",zero(15),p.zero(15));
if(zero(16)!=0)checkpr(fout, "s40",zero(16),p.zero(16));
if(zero(17)!=0)checkpr(fout, "s50",zero(17),p.zero(17));
if(zero(18)!=0)checkpr(fout, "s60",zero(18),p.zero(18));

checkpr(fout, "hmin",qmin(1),p.qmin(1));
checkpr(fout, "hmax",qmax(1),p.qmax(1));
checkpr(fout, "deltah",deltaq(1),p.deltaq(1));
checkpr(fout, "kmin",qmin(2),p.qmin(2));
checkpr(fout, "kmax",qmax(2),p.qmax(2));
checkpr(fout, "deltak",deltaq(2),p.deltaq(2));
checkpr(fout, "lmin",qmin(3),p.qmin(3));
checkpr(fout, "lmax",qmax(3),p.qmax(3));
checkpr(fout, "deltal",deltaq(3),p.deltaq(3));
checkpr(fout, "maxqperiod",maxqperiod,p.maxqperiod);
checkpr(fout, "maxnofspins",maxnofspins,p.maxnofspins);
checkpr(fout, "minnr1",minnr1,p.minnr1);
checkpr(fout, "minnr2",minnr2,p.minnr2);
checkpr(fout, "minnr3",minnr3,p.minnr3);
checkpr(fout, "nofrndtries",nofrndtries,p.nofrndtries);
checkpr(fout, "nofMCsteps",nofMCsteps,p.nofMCsteps);
checkpr(fout, "maxnoftestspincf",maxnoftestspincf,p.maxnoftestspincf);

checkpr(fout, "maxnofmfloops",maxnofmfloops,p.maxnofmfloops);
checkpr(fout, "maxstamf",maxstamf,p.maxstamf);
checkpr(fout, "bigstep",bigstep,p.bigstep);
checkpr(fout, "repeat",repeat,p.repeat);
checkpr(fout, "maxspinchange",maxspinchange,p.maxspinchange);
checkpr(fout, "nofspincorrs",nofspincorrs,p.nofspincorrs);
checkpr(fout, "maxnofhkls",maxnofhkls,p.maxnofhkls);
checkpr(fout, "maxQ",maxQ,p.maxQ);

}



void inipar::print (FILE * fout)
{
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
                  \n# vector in (H-T) space corresponding to x axis (xT [K] xHa [T] xHb [T] xHc [T]) \
                  \n# optional are also xHi[T] xHj[T] xHk[T] (magnetic field in ijk coordinates, \
		  \n# defined by  j||b, k||(a x b) and i normal to k and j ), \
		  \n# electric field xEa[kV/mm] xEb[kV/mm] xEc[kV/mm]  xEi[kV/mm] xEj[kV/mm] xEk[kV/mm] \
		  \n# stress tensor in Voigt notation (1,2,3,4,5,6 = ii jj kk jk ik ij) \
                  \n# xs1[GPa] xs2[GPa] xs3[GPa] xs4[GPa] xs5[GPa] xs6[GPa] \
                  \n");

    fprintf(fout,"xT=%g\nxHa=%g\nxHb=%g\nxHc=%g\n# range of x\nxmin=%g\nxmax=%g\nxstep=%g\n",
           xv(0), xv(1), xv(2), xv(3), xmin,  xmax,  xstep);
    if(xv(4)!=0)fprintf(fout,"xHi=%g\n",xv(4));
    if(xv(5)!=0)fprintf(fout,"xHj=%g\n",xv(5));
    if(xv(6)!=0)fprintf(fout,"xHk=%g\n",xv(6));
    if(xv(7)!=0)fprintf(fout,"xEa=%g\n",xv(7));
    if(xv(8)!=0)fprintf(fout,"xEb=%g\n",xv(8));
    if(xv(9)!=0)fprintf(fout,"xEc=%g\n",xv(9));
    if(xv(10)!=0)fprintf(fout,"xEi=%g\n",xv(10));
    if(xv(11)!=0)fprintf(fout,"xEj=%g\n",xv(11));
    if(xv(12)!=0)fprintf(fout,"xEk=%g\n",xv(12));
    if(xv(13)!=0)fprintf(fout,"xs1=%g\n",xv(13));
    if(xv(14)!=0)fprintf(fout,"xs2=%g\n",xv(14));
    if(xv(15)!=0)fprintf(fout,"xs3=%g\n",xv(15));
    if(xv(16)!=0)fprintf(fout,"xs4=%g\n",xv(16));
    if(xv(17)!=0)fprintf(fout,"xs5=%g\n",xv(17));
    if(xv(18)!=0)fprintf(fout,"xs6=%g\n",xv(18));

    fprintf(fout,"# vector in (H-T) space corresponding to y axis (yT [K] yHa [T] yHb [T] yHc [T])\n");

    fprintf(fout,"yT=%g\nyHa=%g\nyHb=%g\nyHc=%g\n# range of y\nymin=%g\nymax=%g\nystep=%g\n",
           yv(0), yv(1), yv(2), yv(3), ymin,  ymax,  ystep);
    if(yv(4)!=0)fprintf(fout,"yHi=%g\n",yv(4));
    if(yv(5)!=0)fprintf(fout,"yHj=%g\n",yv(5));
    if(yv(6)!=0)fprintf(fout,"yHk=%g\n",yv(6));
    if(yv(7)!=0)fprintf(fout,"yEa=%g\n",yv(7));
    if(yv(8)!=0)fprintf(fout,"yEb=%g\n",yv(8));
    if(yv(9)!=0)fprintf(fout,"yEc=%g\n",yv(9));
    if(yv(10)!=0)fprintf(fout,"yEi=%g\n",yv(10));
    if(yv(11)!=0)fprintf(fout,"yEj=%g\n",yv(11));
    if(yv(12)!=0)fprintf(fout,"yEk=%g\n",yv(12));
    if(yv(13)!=0)fprintf(fout,"ys1=%g\n",yv(13));
    if(yv(14)!=0)fprintf(fout,"ys2=%g\n",yv(14));
    if(yv(15)!=0)fprintf(fout,"ys3=%g\n",yv(15));
    if(yv(16)!=0)fprintf(fout,"ys4=%g\n",yv(16));
    if(yv(17)!=0)fprintf(fout,"ys5=%g\n",yv(17));
    if(yv(18)!=0)fprintf(fout,"ys6=%g\n",yv(18));

    fprintf(fout,"# offset for phase diagram\n");
    fprintf(fout,"T0=%g\nHa0=%g\nHb0=%g\nHc0=%g\n\n",zero(0),zero(1),zero(2),zero(3));       
    
    if(zero(4)!=0)fprintf(fout,"Hi0=%g\n",zero(4));
    if(zero(5)!=0)fprintf(fout,"Hj0=%g\n",zero(5));
    if(zero(6)!=0)fprintf(fout,"Hk0=%g\n",zero(6));
    if(zero(7)!=0)fprintf(fout,"Ea0=%g\n",zero(7));
    if(zero(8)!=0)fprintf(fout,"Eb0=%g\n",zero(8));
    if(zero(9)!=0)fprintf(fout,"Ec0=%g\n",zero(9));
    if(zero(10)!=0)fprintf(fout,"Ei0=%g\n",zero(10));
    if(zero(11)!=0)fprintf(fout,"Ej0=%g\n",zero(11));
    if(zero(12)!=0)fprintf(fout,"Ek0=%g\n",zero(12));
    if(zero(13)!=0)fprintf(fout,"s10=%g\n",zero(13));
    if(zero(14)!=0)fprintf(fout,"s20=%g\n",zero(14));
    if(zero(15)!=0)fprintf(fout,"s30=%g\n",zero(15));
    if(zero(16)!=0)fprintf(fout,"s40=%g\n",zero(16));
    if(zero(17)!=0)fprintf(fout,"s50=%g\n",zero(17));
    if(zero(18)!=0)fprintf(fout,"s60=%g\n",zero(18));

    fprintf(fout,"#input (xHa xHb xHc) (yHa yHb yHc) and (Ha0 Hb0 Hc0) are vectors\n");
    fprintf(fout,"#given in terms of components with respect to unit vectors along the\n");
    fprintf(fout,"#Bravais lattice ^a=a/|a|, ^b=b/|b|, ^c=c/|c|.\n");
    fprintf(fout,"#For the external magnetic field unit is Tesla.\n");
    fprintf(fout,"#For the external electric field unit is kV/mm.\n");
    fprintf(fout,"#For the external stress tensor the unit is GPa.\n\n");

    fprintf(fout,"# out variables to control first columns of output files results/mcphas.*:\n");
    for(int i=1;i<=usrdefcols[0];++i)fprintf(fout,"out%i=%i \n",usrdefcols[i],colcod[i]);
    fprintf(fout,"#     ... in out*=n the numbers n have the following meaning:\n");
    for(int i=0;i<=COLHEADDIM;++i){
    fprintf(fout,"#            %i....%s\n",i,colhead[i]);
                   }
 
    fprintf(fout,"\n[GENERATION OF SPIN CONFIGURATIONS]\n");
    fprintf(fout,"# test q vector (qmin qmax deltaq)\n");
    fprintf(fout,"hmin=%g\nhmax=%g\ndeltah=%g\n",qmin(1),qmax(1),deltaq(1));
    fprintf(fout,"kmin=%g\nkmax=%g\ndeltak=%g\n",qmin(2),qmax(2),deltaq(2));
    fprintf(fout,"lmin=%g\nlmax=%g\ndeltal=%g\n",qmin(3),qmax(3),deltaq(3));
    fprintf(fout,"# maximal periodicity of spinconfigurations generated by q vectors\n");
    fprintf(fout,"maxqperiod=%i\n",maxqperiod);
    fprintf(fout,"# maximal number of spins in spinconfigurations generated by q vectors\n");
    fprintf(fout,"maxnofspins=%i\n",maxnofspins);
    fprintf(fout,"# minimal dimension of supercell along r1 r2 r3 \n");
    fprintf(fout,"minnr1=%i\n",minnr1);
    fprintf(fout,"minnr2=%i\n",minnr2);
    fprintf(fout,"minnr3=%i\n",minnr3);
   fprintf(fout,"# number of random  seed spins to try for each initial spinconfiguration\n");
               
    fprintf(fout,"nofrndtries=%i\n",nofrndtries);
    fprintf(fout,"# number of  Metropolis Monte Carlo steps (when maxnofmfloops\n"
                "# has been reached)\n");
    fprintf(fout,"nofMCsteps=%i\n",nofMCsteps);

    fprintf(fout,"# maximum number of test spin configurations in table\n");
    fprintf(fout,"maxnoftestspincf=%i\n\n",maxnoftestspincf);

    fprintf(fout,"[PARAMETERS FOR SUB FECALC SELFCONSISTENCY PROCESS]\
                  \n# maximum number of selfconsistency loops\n");
    fprintf(fout,"maxnofmfloops=%i\n",maxnofmfloops);
    fprintf(fout,"# standard deviation - limit to end selfconsistency process \
                  \n# standard deviation is defined by ...sta=sqrt(sum_{i=1}^{n} (newmf-old mf)i^2/n) \
		  \n# the meanfield is given by mf=gj mb H [meV] (gj...lande factor, mb... bohr magneton)\n");
    fprintf(fout,"maxstamf=%g\n",maxstamf);
    fprintf(fout,"# mean field step ratio bigstep( = step to perform /calculated step<1) \n");
    fprintf(fout,"# note: if sta increases - then for 10 iterations set step ratio to smallstep=bigstep/n\n");
    fprintf(fout,"# by default n=5. However, if bigstep>1 then n=integervalue(bigstep) and step ratio=bigstep-n \n");
    fprintf(fout,"bigstep=%g\n",bigstep);

    fprintf(fout,"# number of repetitions if all mean field loop fail to stabilise\n");
    fprintf(fout,"# at each repetition either maxstamf or manxofloops or maxspinchange is relaxed\n");
    fprintf(fout,"# (depending on which problems occur most) to allow for more computation time,e.g. \n");
    fprintf(fout,"# repeat 3.4 will allow for 3 repetitions and if maxnofmfloops is too small to converge\n");
    fprintf(fout,"# it will be increased by a factor 1/0.4 for each repetition\n");
    fprintf(fout,"repeat=%g\n",repeat);

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
}

//constructor ... load initial parameters from file
inipar::inipar (const char * file,char * pref,const char * prog)
{  program= new char [strlen(prog)+1];
   strcpy(program,prog);
   savfilename= new char [strlen(file)+strlen(pref)+1];
  if(pref[0]!='\0')strcpy(savfilename,pref);
  strcpy(savfilename+strlen(pref),file);
  prefix = new char[MAXNOFCHARINLINE];
  strcpy(prefix,pref);
  xv=Vector(0,EXTERNAL_PARAMETER_DIMENSION-1);yv=Vector(0,EXTERNAL_PARAMETER_DIMENSION-1);zero=Vector(0,EXTERNAL_PARAMETER_DIMENSION-1);
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
  doeps=0;linepscf=0;linepsjj=0;ipx=NULL;ipy=NULL;ipz=NULL;include_cd=false;
  ipeps1=NULL;ipeps2=NULL;ipeps3=NULL;ipeps4=NULL;ipeps5=NULL;ipeps6=NULL;
  printf("reading file %s\n",savfilename);
  if(load()!=0){if(pref[0]!='\0'){fprintf(stderr,"File %s not found - trying %s\n",savfilename,file);
                strcpy(savfilename,file);}
                if(load()!=0){fprintf(stderr,"# Warning: Cannot load file %s - using default values ! \n",savfilename); 
  // set default values
  xv=0;xv(0)=1;yv=0;yv(3)=1;xmin=1;xmax=1;ymin=0;ymax=0;xstep=1;ystep=1;
  qmin=0;qmax=0;deltaq(1)=0.1;deltaq(2)=0.1;deltaq(3)=0.1;maxqperiod=1;maxnofspins=10;nofrndtries=0;
  minnr1=0;
  minnr2=0;
  minnr3=0;
  maxnofmfloops=100;maxstamf=1e-3;bigstep=1;maxspinchange=100;zero=0;repeat=0;
nofthreads=0;getnofthread(nofthreads);
  nofspincorrs=0;maxnofhkls=5;maxQ=3;maxnoftestspincf=1000;
  nofstapoints=0;
  nofreppoints=0;
  nofconvrep=0;
  nofmaxloopDIV=0;nofmaxspinchangeDIV=0;
  successrate=0;
  nofcalls=0;
  noffailedpoints=0;
  sta=0;
  print();
 // append also the default values for different prefixes 
 // scan Ha Hb Hc Hi Hj Hk chia chib chic 
 FILE * fout = fopen_errchk (savfilename,"a");
  fprintf(fout,
   "# prefix for field scan parallel a\n"
   "Ha_yHa=1\n"
   "Ha_yHb=0\n"
   "Ha_yHc=0\n"
   "Ha_ymin=0\n"
   "Ha_ymax=20\n"
   "Ha_ystep=0.5\n"
   "# prefix for field scan parallel b\n"
   "Hb_yHa=0\n"
   "Hb_yHb=1\n"
   "Hb_yHc=0\n"
   "Hb_ymin=0\n"
   "Hb_ymax=20\n"
   "Hb_ystep=0.5\n"
   "# prefix for field scan parallel c\n"
   "Hc_yHa=0\n"
   "Hc_yHb=0\n"
   "Hc_yHc=1\n"
   "Hc_ymin=0\n"
   "Hc_ymax=20\n"
   "Hc_ystep=0.5\n"
   "# prefix for T scan with field parallel a\n"
   "chia_yHa=1\n"
   "chia_yHb=0\n"
   "chia_yHc=0\n"
   "chia_ymin=1\n"
   "chia_ymax=1\n"
   "chia_ystep=0.5\n"
"# prefix for T scan with field parallel a\n"
   "chib_yHa=0\n"
   "chib_yHb=1\n"
   "chib_yHc=0\n"
   "chib_ymin=1\n"
   "chib_ymax=1\n"
   "chib_ystep=0.5\n"
"# prefix for T scan with field parallel a\n"
   "chic_yHa=0\n"
   "chic_yHb=0\n"
   "chic_yHc=1\n"
   "chic_ymin=1\n"
   "chic_ymax=1\n"
   "chic_ystep=0.5\n"
);
  fclose(fout);

                              }
                }
}

//kopier-konstruktor 
inipar::inipar (const inipar & p)
{ program= new char [strlen(p.program)+1];
  strcpy(program,p.program);
  savfilename= new char [strlen(p.savfilename)+1];
  strcpy(savfilename,p.savfilename);
  prefix = new char[MAXNOFCHARINLINE];
  strcpy(prefix,p.prefix);
  doeps=p.doeps;outcolset=p.outcolset;include_cd=p.include_cd;
  linepscf=p.linepscf;
  linepsjj=p.linepsjj;
  ipx=p.ipx;
  ipy=p.ipy;
  ipz=p.ipz;
  ipeps1=p.ipeps1;
  ipeps2=p.ipeps2;
  ipeps3=p.ipeps3;
  ipeps4=p.ipeps4;
  ipeps5=p.ipeps5;
  ipeps6=p.ipeps6;
  sta=p.sta;
  startcputime=p.startcputime;
  nofstapoints=p.nofstapoints;
  nofreppoints=p.nofreppoints;
  nofconvrep=p.nofconvrep;
  nofmaxloopDIV=p.nofmaxloopDIV;
  nofmaxspinchangeDIV=p.nofmaxspinchangeDIV;
  successrate=p.successrate; 
  nofcalls=p.nofcalls;
  noffailedpoints=p.noffailedpoints;
  exit_mcphas=p.exit_mcphas;pause_mcphas=p.pause_mcphas;
  displayall=p.displayall;logfevsQ=p.logfevsQ;
  
  
  xv=Vector(0,EXTERNAL_PARAMETER_DIMENSION-1);yv=Vector(0,EXTERNAL_PARAMETER_DIMENSION-1);
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
  xv=p.xv;xmin=p.xmin;xmax=p.xmax;xstep=p.xstep;
  yv=p.yv;ymin=p.ymin;ymax=p.ymax;ystep=p.ystep;
  zero=Vector(0,EXTERNAL_PARAMETER_DIMENSION-1);
  zero=p.zero;

  qmin=p.qmin;
  qmax=p.qmax;
  deltaq=p.deltaq;  
  maxqperiod=p.maxqperiod;
  minnr1=p.minnr1;
  minnr2=p.minnr2;
  minnr3=p.minnr3;
  maxnofspins=p.maxnofspins;
  nofrndtries=p.nofrndtries;
  nofMCsteps=p.nofMCsteps;
  maxnoftestspincf=p.maxnoftestspincf;

  maxnofmfloops=p.maxnofmfloops;
  maxstamf=p.maxstamf;
  bigstep=p.bigstep;
  repeat=p.repeat;
  maxspinchange=p.maxspinchange;
  
  nofspincorrs=p.nofspincorrs;
  maxnofhkls=p.maxnofhkls;
  maxQ=p.maxQ;
  nofthreads=p.nofthreads;
}

//destruktor
inipar::~inipar ()
{//printf("hello destruktor inipar\n");  
 
delete []savfilename;
delete []prefix;
//printf("hello destruktor inipar\n");  
 }


//***************************************************************
//constructor ... load initial parameters from file
inipars::inipars (const char * file,char * pref,const char * prog)
{ inis=new inipar*[MAXNOFINIS];
  char * lofprefixes[MAXNOFINIS];
  int nofinisold=-1;nofinis=0;
// here we have to load inis[1...nofinis] with different prefixes matching pref - until no new matching
// prefix is found ...
  while(nofinisold<nofinis&&nofinis<MAXNOFINIS)
  {inis[nofinis]=new inipar(file,pref,prog);
   nofinisold=nofinis;(*inis[nofinis]).load(nofinis,lofprefixes);
  }
 

// remove last inis, because it does not contain a new prefix
if(nofinis>0&&nofinis<MAXNOFINIS){delete inis[nofinis];} 
if(nofinis==0)nofinis=1;
}

void inipars::saveexitzero()
{// read file and put exit=0 and save it again
 char *lines[MAXNOFLINES];
 char instr[MAXNOFCHARINLINE];
 FILE * fout = fopen_errchk ((*inis[0]).savfilename,"r");int i=0;
 while (fgets(instr,MAXNOFCHARINLINE,fout)!=NULL)
 {char * p; p=strstr(instr,"exit=1");
  if(p!=NULL)
    {p[5]='0';
    }
  lines[i]=new char [strlen(instr)+1];strncpy(lines[i],instr,strlen(instr));++i;
 }
 fclose(fout);
 fout=fopen_errchk ((*inis[0]).savfilename,"w");
 for(int ii=0;ii<i;++ii)
 {fprintf(fout,"%s",lines[ii]);
  delete [] lines [ii];
 }
 fclose(fout);

}



//kopier-konstruktor  inipars
inipars::inipars (const inipars & p)
{ nofinis=p.nofinis;
  inis=new inipar*[MAXNOFINIS];
  for(int i=0;i<nofinis;++i)inis[i]=new inipar((*p.inis[i]));
}

//destruktor inipars
inipars::~inipars ()
{//printf("hello destruktor inipars %i\n",nofinis);  
 for(int i=0;i<nofinis;++i)delete  inis[i];
delete []inis;
//printf("hello destruktor inipar\n");  
 }
