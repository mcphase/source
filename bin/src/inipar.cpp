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
printf (" use as: mcphas \n or as: mcphas [file]\n");
printf (" [file] ... input file  with sets of x y T H Hi Hj Hk points \n");
printf (" (format as output file mcphas.xyt)\n\n");
printf (" Options: -h     print this help screen\n");
printf ("          -stamax 14  ... end mcphas if standard deviation exceeds 14\n");
printf ("          -a     append output files (do not overwrite) \n");
printf ("          -doeps refine strain epsilon selfconsistently using elastic,magnetoelastic constants \n");
printf ("                 read from mcphas.j and mcphas.djdx mcphas.djdy and mcphas.djdz \n");
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
                            "Ea [V/m]",  //     7      
                            "Eb [V/m]",  //     8      
                            "Ec [V/m]",  //     9      
                            "Ei [V/m]",  //    10      
                            "Ej [V/m]",  //    11      
                            "Ek [V/m]",  //    12     
                            "s1 [Pa]",  //    13    
                            "s2 [Pa]",  //    14      
                            "s3 [Pa]",  //    15      
                            "s4 [Pa]",  //    16      
                            "s5 [Pa]",  //    17      
                            "s6 [Pa]",  //    18 
                            "x",  //    19 
                            "y", //     20 
                            "|H| [T]",   //    21
                            "|E| [T]"   //    22
                               };
 bool inipar::defaultcolcode(int col,int colcode) // resets default columns if not set by user (outcolset==true)
{                                             // returns true if reset has been successful
 if(!outcolset)for(int i=1;i<=usrdefcols[0];++i)if(usrdefcols[i]==col)
 {colcod[i]=colcode;}
}
// different output data for user defined columns ...
double inipar::setcolvalue(int i,float & x, float & y,double& T,Vector & Hext,Vector & abc)
{      switch (i) {
case 0:  return T;break;
case 1:  
case 2:  
case 3:  {Vector Habc(1,3);Vector abcu(1,6);abcu=abc;abcu(1)=1;abcu(2)=1;abcu(3)=1;Vector v(1,3);v=Hext(1,3);
          ijk2dadbdc(Habc,v,abcu);
         return Habc(i);
         }break;
case 4:  return Hext(1);break;
case 5:  return Hext(2);break;
case 6:  return Hext(3);break;
case 7:  
case 8:  
case 9:{Vector Eabc(1,3);Vector abcu(1,6);abcu=abc;abcu(1)=1;abcu(2)=1;abcu(3)=1;Vector v(1,3);v=Hext(4,6,-3);
          ijk2dadbdc(Eabc,v,abcu);
         return Eabc(i-6);
         }break;
case 10:  return Hext(4);break;
case 11:  return Hext(5);break;
case 12:  return Hext(6);break;
case 13:  return Hext(7);break;
case 14:  return Hext(8);break;
case 15:  return Hext(9);break;
case 16:  return Hext(10);break;
case 17:  return Hext(11);break;
case 18:  return Hext(12);break;
case 19:  return x;break;
case 20:  return y;break;
case 21:  {Vector Hijk(1,3);Hijk=Hext(1,3);return Norm(Hijk);}break;
case 22:  {Vector Eijk(1,3);Eijk=Hext(4,6);return Norm(Eijk);}break;
default: fprintf(stderr,"Error mcphas: unknown column code\n");exit(EXIT_FAILURE);
                    }

return 0;
}



// print user defined column headers
void inipar::print_usrdefcolhead(FILE *fout,char * str)
{fprintf(fout,"#");
 int i;
 for(i=1;i<=usrdefcols[0];++i)fprintf(fout,"%i%*s",i,(int)strlen(colhead[colcod[i]]),"");
 char *t;size_t n;
 for(t=str;t[0]!='\0';t+=n)
 {
 //find first nonspace character
 n= strspn(t," \t"); //printf("%i %i\n",n,t);
 fprintf(fout,"%*s",(int)n,"");
 // next - find until a space occurs
 t+=n;n=strcspn(t," \t");
 if(t[n]!='\0')fprintf(fout,"%2.i%*s",i,(int)(n-2),""); // fill the corresponding space with number and spaces
 ++i;
}
 fprintf(fout,"\n#");
 for(i=1;i<=usrdefcols[0];++i)fprintf(fout,"%s ",colhead[colcod[i]]); 
 fprintf(fout,"%s\n",str);
}


// set external field and Temperature
void inipar::getTH(double & T,Vector & h,double x, double y,cryststruct & cs)
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

    v=yv(1,3);normalizedadbdc(v,1.0,cs);// take care that vector xHa xHb Xhc has unit length 1 Tesla
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

    v=xv(7,9,-6);normalizedadbdc(v,1.0,cs); // take care that vector xHa xHb Xhc has unit length 1 Tesla
    dadbdc2ijk(v1,v,abc); // transform Habc to ijk coordinates ... this is H
    v=xv(10,12,-9);normalize(v,1.0); // take care of Hijk to be of unit length 1 Tesla
    v1=v1+v;
     
    v=yv(7,9,-6);normalizedadbdc(v,1.0,cs);// take care that vector xHa xHb Xhc has unit length 1 Tesla
    dadbdc2ijk(v2,v,abc); // transform Habc to ijk coordinates ... this is H
    v=yv(10,12,-9);normalize(v,1.0); // take care of Hijk to be of unit length 1 Tesla
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
     h(i+6)=zero(i+6)+x*xs(i)+y*ys(i);}


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
   default: d=setcolvalue(colcod[i],x,y, T,Hext, abc)-nn[i];
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
 { double val=setcolvalue(colcod[i],x,y,T,Hext,abc);
   if(withtext)fprintf(fout,"%s=%4.4g ",colhead[colcod[i]],myround(val));
   else fprintf(fout,"%*s%4.4g ",(int)(strlen(colhead[colcod[i]])-8 < 0 ? :0),"",myround(val));
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
   if(!cc){fprintf(stderr,"#Warning: External Magnetic Field H nonzero but not stored in output files  - please change settings out out* in mcphas.ini\n");exit(EXIT_FAILURE); }
                           }break;
   default: if(!c[i]){fprintf(stderr,"#Warning: stress s%i nonzero but not stored in output files  - please change settings out out* in mcphas.ini\n",i-6);exit(EXIT_FAILURE); }
  }
 }
}

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
{ FILE *fin_coq;outcolset=false;
  char instr[MAXNOFCHARINLINE];
  char somestring[MAXNOFCHARINLINE];
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
     
    extract_with_prefix(instr,prefix,"xT",xv[0]);      
    extract_with_prefix(instr,prefix,"xHa",xv[1]);
    extract_with_prefix(instr,prefix,"xHb",xv[2]);    
    extract_with_prefix(instr,prefix,"xHc",xv[3]);
    extract_with_prefix(instr,prefix,"xHi",xv[4]);
    extract_with_prefix(instr,prefix,"xHj",xv[5]);    
    extract_with_prefix(instr,prefix,"xHk",xv[6]);
    extract_with_prefix(instr,prefix,"xEa",xv[7]);
    extract_with_prefix(instr,prefix,"xEb",xv[8]);    
    extract_with_prefix(instr,prefix,"xEc",xv[9]);
    extract_with_prefix(instr,prefix,"xEi",xv[10]);
    extract_with_prefix(instr,prefix,"xEj",xv[11]);    
    extract_with_prefix(instr,prefix,"xEk",xv[12]);
    extract_with_prefix(instr,prefix,"xs1",xv[13]);
    extract_with_prefix(instr,prefix,"xs2",xv[14]);    
    extract_with_prefix(instr,prefix,"xs3",xv[15]);
    extract_with_prefix(instr,prefix,"xs4",xv[16]);
    extract_with_prefix(instr,prefix,"xs5",xv[17]);    
    extract_with_prefix(instr,prefix,"xs6",xv[18]);

    extract_with_prefix(instr,prefix,"xmin",xmin);  extract_with_prefix(instr,prefix,"xmax",xmax);
    extract_with_prefix(instr,prefix,"xstep",xstep);
   
    extract_with_prefix(instr,prefix,"yT",yv[0]);      
    extract_with_prefix(instr,prefix,"yHa",yv[1]);
    extract_with_prefix(instr,prefix,"yHb",yv[2]);    
    extract_with_prefix(instr,prefix,"yHc",yv[3]);
    extract_with_prefix(instr,prefix,"yHi",yv[4]);
    extract_with_prefix(instr,prefix,"yHj",yv[5]);    
    extract_with_prefix(instr,prefix,"yHk",yv[6]);
    extract_with_prefix(instr,prefix,"yEa",yv[7]);
    extract_with_prefix(instr,prefix,"yEb",yv[8]);    
    extract_with_prefix(instr,prefix,"yEc",yv[9]);
    extract_with_prefix(instr,prefix,"yEi",yv[10]);
    extract_with_prefix(instr,prefix,"yEj",yv[11]);    
    extract_with_prefix(instr,prefix,"yEk",yv[12]);
    extract_with_prefix(instr,prefix,"ys1",yv[13]);
    extract_with_prefix(instr,prefix,"ys2",yv[14]);    
    extract_with_prefix(instr,prefix,"ys3",yv[15]);
    extract_with_prefix(instr,prefix,"ys4",yv[16]);
    extract_with_prefix(instr,prefix,"ys5",yv[17]);    
    extract_with_prefix(instr,prefix,"ys6",yv[18]);
    
    extract_with_prefix(instr,prefix,"ymin",ymin); extract_with_prefix(instr,prefix,"ymax",ymax);
    extract_with_prefix(instr,prefix,"ystep",ystep);   
 
    extract_with_prefix(instr,prefix,"T0",zero[0]);      
    extract_with_prefix(instr,prefix,"Ha0",zero[1]);
    extract_with_prefix(instr,prefix,"Hb0",zero[2]);    
    extract_with_prefix(instr,prefix,"Hc0",zero[3]);
    extract_with_prefix(instr,prefix,"Hi0",zero[4]);
    extract_with_prefix(instr,prefix,"Hj0",zero[5]);    
    extract_with_prefix(instr,prefix,"Hk0",zero[6]);
    extract_with_prefix(instr,prefix,"Ea0",zero[7]);
    extract_with_prefix(instr,prefix,"Eb0",zero[8]);    
    extract_with_prefix(instr,prefix,"Ec0",zero[9]);
    extract_with_prefix(instr,prefix,"Ei0",zero[10]);
    extract_with_prefix(instr,prefix,"Ej0",zero[11]);    
    extract_with_prefix(instr,prefix,"Ek0",zero[12]);
    extract_with_prefix(instr,prefix,"s10",zero[13]);
    extract_with_prefix(instr,prefix,"s20",zero[14]);    
    extract_with_prefix(instr,prefix,"s30",zero[15]);
    extract_with_prefix(instr,prefix,"s40",zero[16]);
    extract_with_prefix(instr,prefix,"s50",zero[17]);    
    extract_with_prefix(instr,prefix,"s60",zero[18]);

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

       for(int j=1;j<=usrdefcols[0];++j) // extract user defined output columns
     {snprintf(somestring,MAXNOFCHARINLINE,"out%i",usrdefcols[j]);
      if(0==extract(instr, somestring,colcod[usrdefcols[j]]))outcolset=true;
     }

    }
   }
  fclose (fin_coq);
 for(int i=1;i<=usrdefcols[0];++i){if(colcod[i]>COLHEADDIM)
 {fprintf(stderr,"Error reading mcphas.ini - out%i = %i > %i not possible !\n",i,colcod[i],COLHEADDIM);exit(EXIT_FAILURE);}
 }

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
                  \n# vector in (H-T) space corresponding to x axis (xT [K] xHa [T] xHb [T] xHc [T]) \
                  \n# optional are also xHi[T] xHj[T] xHk[T] (magnetic field in ijk coordinates, \
		  \n# defined by  j||b, k||(a x b) and i normal to k and j ), \
		  \n# electric field xEa[V/m] xEb[V/m] xEc[V/m]  xEi[V/m] xEj[V/m] xEk[V/m] \
		  \n# stress tensor in Voigt notation (1,2,3,4,5,6 = ii jj kk jk ik ij) \
                  \n# xs1[Pa] xs2[Pa] xs3[Pa] xs4[Pa] xs5[Pa] xs6[Pa] \
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
  xv=Vector(0,EXTERNAL_PARAMETER_DIMENSION-1);yv=Vector(0,EXTERNAL_PARAMETER_DIMENSION-1);zero=Vector(0,EXTERNAL_PARAMETER_DIMENSION-1);
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
  doeps=p.doeps;outcolset=p.outcolset;
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
