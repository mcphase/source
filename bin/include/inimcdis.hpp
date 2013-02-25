//  class inipar ... initial parameters for program mcphas
//
#ifndef INIMCDIS
#define INIMCDIS


#include<float.h>
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
  int calculate_magmoment_oscillation; //  creates mcdisp.qem
  int calculate_spinmoment_oscillation; //  creates mcdisp.qes
  int calculate_orbmoment_oscillation; //  creates mcdisp.qeo
  int calculate_chargedensity_oscillation; //  creates mcdisp.qee
  int calculate_spindensity_oscillation; //  creates mcdisp.qsd
  int calculate_orbmomdensity_oscillation; //  creates mcdisp.qod
  int calculate_phonon_oscillation; //  creates mcdisp.qep
  double T;
  Vector Hext;
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

