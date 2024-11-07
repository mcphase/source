//  class inipar ... initial parameters for program mcphas
//
#ifndef INIPAR
#define INIPAR


#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<martin.h>
#include<vector.h>
#include<par.hpp>

class inipar
{ private:
  
  public:
  char * savfilename;
  int doeps,linepscf,linepsjj;\
  par * ipx;par * ipy;par * ipz;

  std::clock_t startcputime;
  int nofstapoints; // number of successful calls to htcalc
  int noffailedpoints; // number of failure calls to htcalc
  int nofmaxloopDIV,nofmaxspinchangeDIV;
  int successrate; // number of successful calls to fecalc
  int nofcalls; // number of calls to fecalc

  //MCPHASE RUNTIME CONTROL
  int exit_mcphas,pause_mcphas,displayall,logfevsQ;
  
  // XY PHASEDIAGRAM PARAMETERS
  Vector xv,yv,zero; // xT xHa xHb xHc ,  yT yHa yHb yHc, T0 Ha0 Hb0 Hc0
  float  xmin,  xmax,  xstep;
  float  ymin,  ymax,  ystep;
  
  // GENERATION OF SPINCONFIGURATIONS
  
  // test qvectors to be considered
  Vector qmin,qmax,deltaq; // hmin hmax deltah kmin kmax deltak lmin lmax deltal
  // maximal periodicity for q vector generated structures
  int maxqperiod;
  // maximal number of spins in qvector generated structure
  int maxnofspins; 
  // number of random spininversions  to try
  // at each configuration
  int nofrndtries;
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
 
  // printout initial parameters to file   
   void print();
   void print (const char * file);
   void time_estimate_until_end(double x, double y);

  // exit with error message
   void errexit();
  //load parameters from file
   int load();
  inipar (const char * file,char * prefix); //constructor
  inipar (const inipar & p);//kopier-konstruktor
 ~inipar ();//destruktor
};

#endif
