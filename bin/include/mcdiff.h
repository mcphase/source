//********************************************
//libraries, classes and  functions in mcdiff
//********************************************

#if defined  (__linux__) || defined (__APPLE__)
#define MAXNOFREFLECTIONS 30000
#else
#define MAXNOFREFLECTIONS 6000
#endif

#define SMALLPOSITIONDEVIATION 1e-4  
#define SMALLINTENSITY 1e-16
#define NOFOUTPUTCOLUMNS 12
#define MAX_NOF_MF_COMPONENTS 51
#include "../../version"
#include <mpspecfunp.h>
#include <martin.h>
#include <myev.h>
#include <jjjpar.hpp>
#include <complex>
#include<cstddef>
#include<spincf.hpp>
#include<mfcf.hpp>
#include<inimcdiff.hpp>


// get intensity of one reflection
//int getint(jjjpar ** jjjpars,int hi,int ki,int li,float thetamax,Vector rez1,Vector rez2, Vector rez3,
// float scale,double T,float lambda,float ovalltemp,int lorenz,int & n,float & d,
// float & Imag,float & Imagdip,float & inuc,float * outn,complex <double> & mqx,
// complex <double> & mqy,complex <double> & mqz,complex <double> & mqxy,complex <double> & mqxz,
// complex <double> & mqyz,complex <double> & mqx2,complex <double> & mqy2,complex <double> & mqz2,
//  Vector & Pxyz);
int getint(inimcdiff & ini,int hi,int ki,int li,Vector rez1,Vector rez2, Vector rez3,
 float scale,float & d,
 float & Imag,float & Imagdip,float & inuc,float * outn,complex <double> & mqx,
 complex <double> & mqy,complex <double> & mqz,complex <double> & mqxy,complex <double> & mqxz,
 complex <double> & mqyz,complex <double> & mqx2,complex <double> & mqy2,complex <double> & mqz2);

// calculate pattern
//void neutint(jjjpar ** jjjpars,int code,double T,float lambda, float thetamax, float ovalltemp,int lorenz,
//             Vector r1,Vector r2,Vector r3,int & n,int & m,Vector *  hkl,float * D,
//             float * intmag,float * ikern,float ** out,complex <double>*mx,
//             complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,
//             complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2,
//             Vector & Pxyz);
void neutint(inimcdiff & ini,int code,int & m,Vector *  hkl,float * D,
             float * intmag,float * ikern,float ** out,complex <double>*mx,
             complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,
             complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2);

// mcdiff - output of results
//void printheader(jjjpar ** jjjpars,int code,const char * filename,const char* infile,char * unitcell,double T,
//Vector & H,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,int * J,int m,
//float a,float b,float c,Vector & P);
void printheader(inimcdiff & ini,int code,int m);

//void printreflist(jjjpar ** jjjpars,int code,const char * filename,const char* infile,char * unitcell,double T,
//              Vector & H,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,int m,
//              Vector * hkl,float * ikern,float * intmag,float * intmagdip,float * D,float ** out,
//              complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,
//              complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2,
//              float a,float b,float c,Vector & P, Vector & Pxyz);
void printreflist(inimcdiff & ini,int code,int m,
              Vector * hkl,float * ikern,float * intmag,float * D,float ** out,
              complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,
              complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2
              );


// output to mcdiff.sps
//void print_sps(int natmagnetic,float a,float b,float c,float alpha,float beta,float gamma,int nr1,int nr2,int nr3,Vector r1s,Vector r2s,Vector r3s,jjjpar ** jjjpars,double T,Vector H);
void print_sps(inimcdiff & ini);

// output to mcdiff.mf
//void print_mf(mfcf & mfields, int natmagnetic,float a,float b,float c,float alpha,float beta,float gamma,int nr1,int nr2,int nr3,Vector r1s,Vector r2s,Vector r3s,jjjpar ** jjjpars,double T,Vector H);
void print_mf(inimcdiff & ini,mfcf & mfields);