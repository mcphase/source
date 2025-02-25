// class physprop for the storage of the physical properties
// of the system at given H and T

#ifndef PHYSPROP
#define PHYSPROP


#include<par.hpp>
#include<inipar.hpp>
#include<martin.h>
#include<spincf.hpp>
#include<mfcf.hpp>

class physproperties
{
  private:
 int washere,nofspincorr;
 
  public:
float x,y; // phasediagramm labels  
int j;  // index of spinstructure
double T; // temperature
Vector m,mabc,H,Pel,Pelabc; // moment and H field and electrical polarisation
double fe;
double u; // free energy and mag energy per ion
double Eel; // elastic energy per ion
int nofatoms;
int nofcomponents;

Vector *jj,*hkli; // spin spin correlation functions
int maxnofhkls,nofhkls;
spincf  sps;
mfcf mf;
   
physproperties (int nofspincorrs,int maxnofhkls,int na, int nm);	//konstruktor
//na number of atoms in basis,nm number of spin components
physproperties (const physproperties & props);	// kopier-konstruktor

~physproperties ();		//destruktor

   // 1.  puts header for fum file columns >8 into string header
   // 2. sets or reads output column field nn -
   //  if setnn true: for all i>8 up to input nofcols ... if nnerr[i]!=0 -> increase sta according to difference 
   //                 (and finally return sta) ... then ...
   //                 set  nn[i] from saved parameters fe,u,etc. and puts into nofcols the number of output columns
   //                 puts into outstr the numbers  nn[i>8] formatted for output into mcphas.fum
   //  if setnn false: reads nn  into parameters fe, u, etc  
 // for fum file
double fumcols(float * nn,float * nnerr, int & nofcols,bool setnn,char * header,char * outstr,inipar & ini,int ortho,par & inputpars,int verbose);
 // for xyt file
double xytcols(float * nn,float * nnerr, int & nofcols,bool setnn,char * header,char * outstr,inipar & ini,int verbose, Vector & totalJ);
 
void update_maxnofhkls(int mxnofhkli);
// save physical properties to output files
double save(int verbose,const char * filemode, int j,inipar & ini,par & inputpars,char * prefix);
// scroll output files and read physical properties from these if possible,
// on success return 0, otherwise
// returns 1
int read(int verbose, par & inputpars,char * readprefix,inipar & ini);
};



#endif
