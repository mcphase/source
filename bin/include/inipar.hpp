//  class inipar ... initial parameters for program mcphasit
//
#ifndef INIPAR
#define INIPAR


#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<martin.h>
#include<vector.h>
#include"cryststruct.hpp"
#include"par.hpp"

#define EXTERNAL_PARAMETER_DIMENSION  HEXT_DIMENSION+7  // dimension of xv, yv zero (see below)



class inipar
{ private:
  bool outcolset; // indicates wether in mcphas.ini user has set some output columns
  
  public:
  bool defaultcolcode(int  col,int  colcode); // resets default columns if not set by user (outcolset==true)
                                             // returns true if reset has been successful
   char * savfilename;
   char * program;
  double sta;
  int doeps,linepscf,linepsjj;
  bool include_cd;
  par * ipx;par * ipy;par * ipz; // storage for two ion interaction parameter derivatives (djdx djdy djdz files)
  par * ipeps1;par * ipeps2;par * ipeps3; // storage for two ion interaction parameter derivatives (djdeps1-6 files)
  par * ipeps4;par * ipeps5;par * ipeps6;

  std::clock_t startcputime;
  int nofstapoints; // number of successful calls to htcalc
  int noffailedpoints; // number of failure calls to htcalc
  int nofmaxloopDIV,nofmaxspinchangeDIV,nofconvrep,nofreppoints;
  int successrate; // number of successful calls to fecalc
  int nofcalls; // number of calls to fecalc

  //MCPHASE RUNTIME CONTROL
  int exit_mcphas,pause_mcphas,displayall,logfevsQ;
  
  // XY PHASEDIAGRAM PARAMETERS
  Vector xv,yv,zero; // xT xHa xHb xHc ,  yT yHa yHb yHc, T0 Ha0 Hb0 Hc0
                     // ... extended (optional) xHi xHj xHk xEa xEb XEc xEi xEj xEk xs1 xs2 xs3 xs4 xs5 xs6
                     // with E1 E2 E3 external field
                     // s1 s2 s3 s4 s5 s6 Voigt components of stress tensor
  float  xmin,  xmax,  xstep;
  float  ymin,  ymax,  ystep;
  
  // GENERATION OF SPINCONFIGURATIONS
  
  // test qvectors to be considered
  Vector qmin,qmax,deltaq; // hmin hmax deltah kmin kmax deltak lmin lmax deltal
  // maximal periodicity for q vector generated structures
  int maxqperiod;
  // maximal number of spins in qvector generated structure
  int maxnofspins; 
  // minimum number of periodicity (for Monte Carlo simulations)
  int minnr1; 
  int minnr2; 
  int minnr3; 
  // number of random seed spins  to try
  // at each configuration
  int nofrndtries;
 // number of random (Monte Carlo) spin inversions  to try for each spins
  int nofMCsteps;
  // maximum number of test spinconfigurations 
  int maxnoftestspincf;

  // Number of threads to use in mcphas
  int nofthreads;
  
  //PARAMETER FOR SUB FECALC - SELFCONSISTENCY PROCESS
  // maximum number of selfconsistency loops
  int maxnofmfloops;
  // sta - limit to end selfconsistency process,
  //standard deviation is defined by ...sta=sqrt(sum_{i=1}^{n} (newmf-old mf)^2/n)  [T]
  float maxstamf;
  // a big step ratio (=step/calculated step) to perform actually
  float bigstep;
  // if a point failes, how often should it be repeated with larger computation time
  float repeat;
  // a small step (=step/calculated step) to perform actually when sta rises
  //float smallstep=0.2;
  //  (<sum abs(actual change of m[mb] with respect to
  // initial  configuration)>) >maxspinchange will  end selfconsistency process
  float maxspinchange;
  char * prefix;

  // OUTPUT OF PHYSICAL PROPERTIES
  // how many spinspin correlation functions 
  // should be calculated
  int nofspincorrs;
  // mximal number of hkls - neutron intensitiest to be calculated
  int maxnofhkls;
  // maximum q[1/A] for hkl's
  double maxQ;

// ***************** intelligent functions **********************
 // set external field and Temperature given x and y
 void calcTHfromxy(double & T,Vector & Hext,double x, double y,cryststruct & cs);

// calculate the value of different output data for user defined column with colcod i...
 double   calccolvalue(int & i,float & x, float & y,double& T,Vector & Hext,Vector & abc);

// return pointer to value of user defined column with colcod i
double * colvaluepointer(int & i,double & x, double & y,double& T,Vector & Hext,Vector & Habc,
                 Vector & Eabc,double & NormH, double & NormE);

 // set external field and Temperature given nn as input from file with meaning defined by out1-7 in mcphas.ini
 // returns true if successful  (NormH NormE x y are not used)
bool calcTHfromnn(double & T,Vector & Hext,float * nn,cryststruct &cs);

 // given T and Hext check if in array nn[0-7] the values are in accordance with T and Hext
 // if yes, returns true ... 
bool checkTH(float * nn,double & T,Vector & Hext,Vector & abc);

  // printout initial parameters to file   
   void print();
   void print (const char * file);
   void print (FILE * fout);
void print_with_prefix(FILE * fout, inipar p);
bool checkpr(FILE* fout,const char * var,int & val,int & masterval);
bool checkpr(FILE* fout,const char * var,double val,double masterval);

// estimeate time until finishing of mcphas
   void time_estimate_until_end(double x, double y);


// output the values into columns with external parameters to file fout
 void print_usrdefcols(FILE *fout,float & x, float & y,double& T,Vector & Hext,Vector & abc,bool withtext);

// print user defined column codes variables out1 -- out7 to fout
 void print_usrdefcolcodes(FILE *fout);

// print column headers for columns with external parameters
 void print_usrdefcolhead(FILE *fout,char * str);

 // exit with error message
   void errexit();
   void finish_mcphas(int  nofqs,int  nofspincf);
  //load parameters from file, returns 1 on error, 0 on success
   int load();
   int load (int & nofinis,char**lofpref);
int extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,float & var);
int extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,double & var);
int extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,int & var);

  inipar (const char * file,char * prefix,const char * prog); //constructor

  inipar (const inipar & p);//kopier-konstruktor

 ~inipar ();//destruktor
};

// define a superclass to store all different prefixes and corresponding parameters
class inipars
{public:
  int    nofinis;
  inipar ** inis;
   void saveexitzero(); // puts exit to zero in input file

  inipars (const char * file,char * prefix, const char * prog); //constructor

  inipars (const inipars & p);//kopier-konstruktor

 ~inipars ();//destruktor
};

#endif
