//  class inipar ... initial parameters for program mcphas
//
#ifndef INIMCDIFF
#define INIMCDIFF


#include<float.h>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<martin.h>
#include<vector.h>
#include<jjjpar.hpp>

class inimcdiff
{ private:
  public:
  int nofoutputcolumns;
   int * colcod;
   char ** colhead;
   int verbose;
   char * infile;
   char * savfilename;
   char * outfilename;
   char * prefix,*unitcellstr;
   double a=0,b=0,c=0,alpha=0,beta=0,gamma=0;
   double T,lambda,thetamax,ovalltemp;
    int lorenz,n,nat, nofatoms,natmagnetic,use_dadbdc=0;
   Vector P,Pxyz,r1,r2,r3,r1s,r2s,r3s,rez1,rez2,rez3;Matrix eps; 
   float *x1;float*y1;float*z1;
   float *da;float*db;float*dc;
   float *sl1r;float*sl1i;float *dwf1;
    jjjpar ** jjjpars;
Matrix rtoijk,rtoijk_rez; // lattice abc and reciprocal lattice in ijk coordinate system
  int nr1=0,nr2=0,nr3=0; 
  int nofthreads; // currently not used in mcdiff
  Vector H;
   void save(); // save parameters to results/_mcdiff.in
   void save(const char * filename); // save parameters to file filename
   void print_usrdefcolhead(FILE *fout);
   void print_usrdefcols(FILE *fout,float ** out,int i);
  inimcdiff (const char * file, char * prefix,int verbose); //constructor
  inimcdiff (const inimcdiff & p);//kopier-konstruktor
 ~inimcdiff ();//destruktor
};

#endif

