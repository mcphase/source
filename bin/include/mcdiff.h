//********************************************
//libraries, classes and  functions in mcdiff
//********************************************

#if defined  (__linux__) || defined (__APPLE__)
#define MAXNOFREFLECTIONS 30000
#else
#define MAXNOFREFLECTIONS 6000
#endif

#define SMALL 1e-3
#define SMALLINT 1e-4
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

int usrdefoutcols[]={5,           4, 5, 6,         10,11}; // user defined output columns (first number is number of usr def output columns)
int colcode[]={-1,     -1,-1,-1, 30,31,32,-1,-1,-1, 1, 0}; // field to store code for assigning type of data to columns of output,
                                                           // set default values here (see list below for different types)
                                                           // using the out11 command in mcdiff.in these codes can be modified
#define COLHEADERDIM 32
// different output data for columns 10 and 11
const char * colheader []= {"LF          ", // 0
                           "|NSF|[b]    ",  // 1
                           "Re(NSF)[b]  ",  // 2
                           "Im(NSF)[b]  ",  // 3 
                           "|MSF|       ",  // 4
                           "|MSF.P|     ",  // 5
                           "Re(MSF.P)   ",  // 6
                           "Im(MSF.P)   ",  // 7
                           "|MSFdip|    ",  // 8 
                           "|MSFdip.P|  ",  // 9
                           "Re(MSFdip.P)",  // 10
                           "Im(MSFdip.P)",  // 11
                           "angl(Q,P)[deg]",  // 12 
                           "i(MSFxMSF*).P", // 13
                           "I+          ",  // 14
                           "I-          ",  // 15
                           "I+/I-       ",  // 16
                           "i(MSFxMSF*)dip.P", // 17
                           "Idip+       ",  // 18
                           "Idip-       ",  // 19
                           "Idip+/Idip- ",  // 20
                           "2*|MSF.P|/sin^2(angl(Q,P)", // 21
                           "2*|MSFdip.P|/sin^2(angl(Q,P)", // 22
                           "2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)", // 23
                           "2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)", // 24
                           "2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)", // 25
                           "2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)", // 26
                           "Qx[1/A]     ", // 27
                           "Qy[1/A]     ", // 28
                           "Qz[1/A]     ", // 29                           
                           "d[A]        ", // 30
                           "|Q|[1/A]    ", // 31
                           "2theta      "  // 32
                           };

// different output data for user defined columns ...
double setcoloutput(int i,float & scale, double & ovallt,float & lorentzf,complex <double> & nsf,float & msf2,float & msf2dip, Vector & Pxyz,
                   complex <double> & msfx, complex <double> & msfy, complex <double> & msfz,
                   complex <double> & msfdipx, complex <double> &msfdipy, complex <double> &msfdipz,Vector & Qvec,
                   float & d,float & theta,double & Q)
{double cosw,crossx,crossy,crossz,Ip,Im,sinw2,R;
         // here do some precalculations with formulas common to several options
         switch (i) {case 13: // beyond cases
                     case 14:
                     case 15:
                     case 16:
                     case 23:
                     case 24:crossx=imag(msfy*conj(msfz)-msfz*conj(msfy));
                             crossy=imag(-msfx*conj(msfz)+msfz*conj(msfx));
                             crossz=imag(msfx*conj(msfy)-msfy*conj(msfx));
                            Ip=abs(nsf) * abs(nsf)+msf2 * 3.65 / 4 / PI;
                            Ip-=(crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3))* 3.65 / 4 / PI;
                            Ip+=sqrt(3.65/4/PI)*real(nsf*conj(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3))+conj(nsf)*(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3)));
                            Im=abs(nsf) * abs(nsf)+msf2 * 3.65 / 4 / PI;
                            Im+=(crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3))* 3.65 / 4 / PI;
                            Im-=sqrt(3.65/4/PI)*real(nsf*conj(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3))+conj(nsf)*(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3)));
                             break;
                     case 17:  // dipolar cases
                     case 18:
                     case 19:
                     case 20:
                     case 25:
                     case 26: crossx=imag(msfdipy*conj(msfdipz)-msfdipz*conj(msfdipy));
                              crossy=imag(-msfdipx*conj(msfdipz)+msfdipz*conj(msfdipx));
                              crossz=imag(msfdipx*conj(msfdipy)-msfdipy*conj(msfdipx));
                             Ip=abs(nsf) * abs(nsf)+msf2dip * 3.65 / 4 / PI;
                             Ip-=(crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3))* 3.65 / 4 / PI;
                             Ip+=sqrt(3.65/4/PI)*real(nsf*conj(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3))+conj(nsf)*(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3)));
                             Im=abs(nsf) * abs(nsf)+msf2dip * 3.65 / 4 / PI;
                             Im+=(crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3))* 3.65 / 4 / PI;
                             Im-=sqrt(3.65/4/PI)*real(nsf*conj(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3))+conj(nsf)*(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3)));
                              break;
                    default: break;
                    }

         switch (i) {
case 0:  return lorentzf;break;//   "LF          ",
case 1:  return abs(nsf);break;//    "|NSF|[b]    ",
case 2:  return real(nsf);break;//    "Re(NSF)[b]  ",
case 3:  return imag(nsf);break;//    "Im(NSF)[b]  ",
case 4:  return sqrt(msf2+1e-100);break;//    "|MSF|       ",
case 5:  return abs(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3));break;//    "|MSF.P|     ",
case 6:  return real(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3));break;//    "Re(MSF.P)   ",
case 7:  return imag(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3));break;//   "Im(MSF.P)   ",
case 8:  return sqrt(msf2dip+1e-100);break;//    "|MSFdip|    ",
case 9:  return abs(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3));break;//    "|MSFdip.P|  ",
case 10: return real(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3));break;//    "Re(MSFdip.P)",
case 11: return imag(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3));break;//    "Im(MSFdip.P)"
case 12: cosw=(Pxyz/Norm(Pxyz))*Qvec/Norm(Qvec);return 180.0 / PI * atan(sqrt(1 - cosw * cosw)/cosw); break; // "angl(Q.P)[�]"
case 13: return (crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3));break; //"i(MSFxMSF*).P",
case 14: return Ip* lorentzf * scale * ovallt;
                     //              "I+          ",
case 15: return Im* lorentzf * scale * ovallt;
                     //              "I-          ",
case 16: return Ip/Im;     //              "I+/I-       "
case 17: return (crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3));break; //i(MSFdip x MSFdip*).P
case 18: return Ip* lorentzf * scale * ovallt;
          //Idip+
case 19: return Im* lorentzf * scale * ovallt;
         // Idip-
case 20: return Ip/Im;         // Idip+/Idip-
case 21: cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2.0*abs(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3))/sinw2;break;//|MSF.P|/sin^2(angl(Q,P)
case 22: cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2.0*abs(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3))/sinw2;break;//|MSFdip.P|/sin^2(angl(Q,P)
case 23: R= Ip/Im;     //   I+/I-
         cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))-sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
case 24: R= Ip/Im;     //   I+/I-
         cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))+sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
case 25:  R= Ip/Im;         // Idip+/Idip-
          cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
          return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))-sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
case 26:  R= Ip/Im;         // Idip+/Idip-
          cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
          return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))+sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
case 27: return Qvec(1);// Qx
case 28: return Qvec(2);// Qy
case 29: return Qvec(3);// Qz
case 30: return d;//"d[A]       
case 31: return Q; //|Q|[1/A]   
case 32: return 2*theta;  // 2theta  
default: fprintf(stderr,"Error mcdiff: unknown column code\n");exit(EXIT_FAILURE);   
         }

return 0;
}

// get intensity of one reflection
int getint(jjjpar ** jjjpars,int hi,int ki,int li,float thetamax,Vector rez1,Vector rez2, Vector rez3,
 float scale,double T,float lambda,float ovalltemp,int lorenz,int & n,float & d,
 float & Imag,float & Imagdip,float & inuc,float * outn,complex <double> & mqx,
 complex <double> & mqy,complex <double> & mqz,complex <double> & mqxy,complex <double> & mqxz,
 complex <double> & mqyz,complex <double> & mqx2,complex <double> & mqy2,complex <double> & mqz2,
  Vector & Pxyz);

// calculate pattern
void neutint(jjjpar ** jjjpars,int code,double T,float lambda, float thetamax, float ovalltemp,int lorenz,
             Vector r1,Vector r2,Vector r3,int & n,int & m,Vector *  hkl,float * D,
             float * intmag,float * intmagdip,float * ikern,float ** out,complex <double>*mx,
             complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,
             complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2,
             Vector & Pxyz);

// mcdiff - output of results
void printheader(jjjpar ** jjjpars,int code,const char * filename,const char* infile,char * unitcell,double T,
Vector & H,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,int * J,int m,
float a,float b,float c,Vector & P);

void printreflist(jjjpar ** jjjpars,int code,const char * filename,const char* infile,char * unitcell,double T,
              Vector & H,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,int m,
              Vector * hkl,float * ikern,float * intmag,float * intmagdip,float * D,float ** out,
              complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,
              complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2,
              float a,float b,float c,Vector & P, Vector & Pxyz);


// output to mcdiff.sps
void print_sps(int natmagnetic,float a,float b,float c,float alpha,float beta,float gamma,int nr1,int nr2,int nr3,Vector r1s,Vector r2s,Vector r3s,jjjpar ** jjjpars,double T,Vector H);

// output to mcdiff.mf
void print_mf(mfcf & mfields, int natmagnetic,float a,float b,float c,float alpha,float beta,float gamma,int nr1,int nr2,int nr3,Vector r1s,Vector r2s,Vector r3s,jjjpar ** jjjpars,double T,Vector H);
