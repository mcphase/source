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

#define NOFHKLCOLUMNS 7

class inimcdis
{ private:
  int do_jqf;
  char * savfilename;
  Vector qmin,qmax,deltaq;
  void read_hkl_list(FILE * finhkl,double ** hkls,int readqxqyqz,int do_jqfile,Vector & abc);   
  double setcolvalue(int i,Vector & Qvec, double & Qincr);

  public:
  int * hklfile_start_index;
  char * info;
  char * prefix;
  double ** hkls;
  int nofhkls; 
  int nofatoms; //nofatoms in primitive cryst unit cell
  int nofcomponents; //number of components of mean field (including magnetic, quadrupolar fields ...
  int calculate_magmoment_oscillation; //  creates mcdisp.qem
  int calculate_spinmoment_oscillation; //  creates mcdisp.qes
  int calculate_orbmoment_oscillation; //  creates mcdisp.qeo
  int calculate_chargedensity_oscillation; //  creates mcdisp.qee
  int calculate_spindensity_oscillation; //  creates mcdisp.qsd
  int calculate_orbmomdensity_oscillation; //  creates mcdisp.qod
  int calculate_phonon_oscillation; //  creates mcdisp.qep
  int outS;
  int nofthreads;
  double T;
  Vector Hext;
  double emax;
  double emin; // energy boundary for dispersion (used for calc. of sta - see manual)
  double ki;
  double kf; // constant ki/kf
  mfcf mf;
   void save(); // save parameters to results/_mcdisp.ini results/_mcdisp.mf
   void print_usrdefcolhead(FILE *fout);
   void print_usrdefcols(FILE *fout,Vector &Qvec, double & Qincr);
  inimcdis (const char * file,const char * spinfile, char * prefix,int do_jqfile,Vector & abc); //constructor
  inimcdis (const inimcdis & p);//kopier-konstruktor
 ~inimcdis ();//destruktor
};

#endif

