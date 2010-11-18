//  class inipar ... initial parameters for program mcphas
//
#ifndef INIMCDIS
#define INIMCDIS


#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<martin.h>
#include<vector.h>
#include<mfcf.hpp>

class inimcdis
{ private:
  char * savfilename;
    
  public:
  int * hklfile_start_index;
  char * info;
  double ** hkls;
  int hkllist; 
  int nofatoms; //nofatoms in primitive cryst unit cell
  int nofcomponents; //number of components of mean field (including magnetic, quadrupolar fields ...
  int extended_eigenvector_dimension; // for the creation of eiegnevectors in mcdisp.eev which 
   // may be used to plot charge density waves. extended_eigenvevtor_dimension must not be 
   // larger than allowed by dmcalc, the routine which creates the matrix Malphabeta (see model section)
  double T;
  double Ha;
  double Hb;
  double Hc;
  double emax;
  double emin; // energy boundary for dispersion (used for calc. of sta - see manual)
  double ki;
  double kf; // constant ki/kf
  // qvectors to be considered
  Vector qmin,qmax,deltaq;
  mfcf mf;
   void save(); // save parameters to results/_mcdisp.ini results/_mcdisp.mf
   void errexit();
  inimcdis (const char * file,const char * spinfile); //constructor
  inimcdis (const inimcdis & p);//kopier-konstruktor
 ~inimcdis ();//destruktor
};

#endif

