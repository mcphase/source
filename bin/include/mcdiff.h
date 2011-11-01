//libraries, classes and  functions in mcdiff




#define MAXNOFREFLECTIONS 6000
#define SMALL 1e-8
#define SMALLINT 1e-4
#define COLHEADERDIM 26
#include "../../version"
#include <mpspecfunp.h>
#include <martin.h>
#include <myev.h>
#include <jjjpar.hpp>
#include <complex>
#include<cstddef>
#include<spincf.hpp>
#include<mfcf.hpp>

// different output data for columns 10 and 11
const char * colheader []= {"LF          ",
                           "|NSF|[b]    ",
                           "Re(NSF)[b]  ",
                           "Im(NSF)[b]  ",
                           "|MSF|       ",
                           "|MSF.P|     ",
                           "Re(MSF.P)   ",
                           "Im(MSF.P)   ",
                           "|MSFdip|    ",
                           "|MSFdip.P|  ",
                           "Re(MSFdip.P)",
                           "Im(MSFdip.P)",
                           "angl(Q,P)[°]",
                           "i(MSFxMSF*).P",
                           "I+          ",
                           "I-          ",
                           "I+/I-       ",
                           "i(MSFxMSF*)dip.P",
                           "Idip+       ",
                           "Idip-       ",
                           "Idip+/Idip- ",
                           "2*|MSF.P|/sin^2(angl(Q,P)",
                           "2*|MSFdip.P|/sin^2(angl(Q,P)",
                           "2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)",
                           "2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)",
                           "2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)",
                           "2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)"
                           };

double setcoloutput(int i,float & scale, double & ovallt,float & lorentzf,complex <double> & nsf,float & msf2,float & msf2dip, Vector & Pxyz,
                   complex <double> & msfx, complex <double> & msfy, complex <double> & msfz,
                   complex <double> & msfdipx, complex <double> &msfdipy, complex <double> &msfdipz,Vector & Qvec);
// get intensity of one reflection
int getint(jjjpar ** jjjpars,int hi,int ki,int li,float thetamax,Vector rez1,Vector rez2, Vector rez3,float scale,
double T,float lambda,float ovalltemp,int lorenz,int & n,int * J,float & d,float & Theta,float & Imag,float & Imagdip,
float & inuc,float & SF,float & lorentzf,complex <double> & mqx,complex <double> & mqy,complex <double> & mqz,
complex <double> & mqxy,complex <double> & mqxz,complex <double> & mqyz,complex <double> & mqx2,
complex <double> & mqy2,complex <double> & mqz2,int * colcode,Vector & Pxyz);

// calculate pattern
void neutint(jjjpar ** jjjpars,int code,double T,float lambda, float thetamax, float ovalltemp,int lorenz,Vector r1,
Vector r2,Vector r3,int & n,int * J,int & m,Vector *  hkl,float * D,float * theta,float * intmag,float * intmagdip,
float * ikern,float * sf,float * lpg,complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,
complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2,
int * colcode,Vector & Pxyz);

// mcdiff - output of results
void printeln(jjjpar ** jjjpars,int code,const char * filename,const char* infile,char * unitcell,double T,
Vector & H,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,int * J,int m,Vector * hkl,
float * ikern,float * intmag,float * intmagdip,float * D,float * theta,float * sf,float * lpg,complex <double>*mx,
complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,complex <double>*mymz,
complex <double>*mx2,complex <double>*my2,complex <double>*mz2,float a,float b,float c,int * colcode,Vector & P);

// output to mcdiff.sps
void print_sps(const char * filename,int natmagnetic,float a,float b,float c,float alpha,float beta,float gamma,int nr1,int nr2,int nr3,Vector r1s,Vector r2s,Vector r3s,jjjpar ** jjjpars,double T,Vector H);

// output to mcdiff.mf
void print_mf(const char * filename,mfcf & mfields, int natmagnetic,float a,float b,float c,float alpha,float beta,float gamma,int nr1,int nr2,int nr3,Vector r1s,Vector r2s,Vector r3s,jjjpar ** jjjpars,double T,Vector H);
