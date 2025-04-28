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
#include<myev.h>


#define NOFHKLCOLUMNS 7

class inimcdis
{ private:
  int do_jqf;
  char * parfile;
  Vector qmin,qmax,deltaq;
  void read_hkl_list(FILE * finhkl,double ** hkls,int readqxqyqz,int do_jqfile,Vector & abc);   
  double setcolvalue(int i,Vector & Qvec, double & Qincr, Vector & qprim,Vector & hkl);
  Vector Eabc,Habc;
  bool outcolset;// indicates wether in mcphas.ini user has set some output columns
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
  int calculate_pel_oscillation; //  creates mcdisp.qpe
  int outS;
  int nofthreads;
  double T;
  Vector Hext;
  double emax;
  double emin; // energy boundary for dispersion (used for calc. of sta - see manual)
  double ki;
  double kf; // constant ki/kf
  mfcf mf;
   void save(); // save parameters to results/_mcdisp.par results/_mcdisp.mf
   void save(const char * filename); // save parameters to file filename
   void print_usrdefcolhead(FILE *fout);
   void print_usrdefcols(FILE *fout,Vector &Qvec, double & Qincr, Vector & qprim,Vector & hkl, bool withtxt=false);
   void mfstring(char *str,size_t t);
   void helpexit();

int extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,float & var);
int extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,double & var);
int extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,int & var);
int extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,char * var,size_t ns,int m);

  inimcdis(const char * file,char * prefix,char * spinfile,
             int & do_jqfile,Vector & abc,
             int & nofcomponents,int & nofatoms);
  int load (char * spinfile, char * prefix,int do_jqfile, Vector & abc,int nofcomp,int nofat); //constructor
  int load (int & nofinis,char**lofpref,char * spinfile, char * prefix,int do_jqfile, Vector & abc,int nofcomp,int nofat); //constructor
  inimcdis (const inimcdis & p);//kopier-konstruktor
 ~inimcdis ();//destruktor
};


// define a superclass to store all different prefixes and corresponding parameters
class inimdpars
{public:
  int    nofinis;
  inimcdis ** inis;
  
  inimdpars (const char * file,char * prefix,char * spinfile,
             int & do_jqfile,Vector & abc,
             int & nofcomponents,int & nofatoms); //constructor

  inimdpars (const inimdpars & p);//kopier-konstruktor

 ~inimdpars ();//destruktor
};
#endif

