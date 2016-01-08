// martins functions
// to get on better in c


#ifndef MARTIN_FOPEN_ERRCHK
#define MARTIN_FOPEN_ERRCHK

#define PI 3.141592654
#define KB 0.0862     // Boltzmanns constant in mev/K
#define MU_B  5.788378E-02 // Bohrmagneton in meV/tesla
#define MAXNOFCHARINLINE 7024

#if __GNUC__ > 2
#include <unistd.h>
#endif

#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<ctime>
#include <stdio.h>
#include<vector.h>
#include "sparsecomplex.hpp"

// extract parameter 'parameter'  from string instr (z.B. "blabla dmin=0.2 blabla") -
// output: var ... value of parameter
// returns 1 on error and 0 if successful
extern  int extract(char * instr,const char * parameter,double & var);

// extract a variable named [parmeter] into var out of a string [instr]
extern   int extract(char * instr,const char * parameter,int & var);
extern   int extract(char * instr,const char * parameter,float & var);
extern   int extract(char * instr,const char * parameter,char * var, size_t n);

// extract a variable which is there also if it is preceded by a prefix
extern   int extract_with_prefix(char * instr,char * prefix, const char * parameter,double & var);
extern   int extract_with_prefix(char * instr,char * prefix, const char * parameter,float & var);
extern   int extract_with_prefix(char * instr,char * prefix, const char * parameter,int & var);
extern   int extract_with_prefix(char * instr,char * prefix, const char * parameter,char * var, size_t n);

// open file: like fopen but with error check 
extern  FILE * fopen_errchk (const char * filename, const char * mode);
// get string instr: like fgets but with error check
extern  char * fgets_errchk (char * instr,int size, FILE * file);


// input line and split into numbers nn[1 ... ]
// on input: nn[0]... size of array nn[] (not changed in routine)
//           fin_coq..filepointer
// returns nof numbers read and
// 0 if end of file or no numbers have been read 
extern  int inputline (FILE * fin_coq, float *nn);

// function to input a line of numbers separated by delimiters
// example:
// parname=3 23 542 23
// returns:0 .... it is a comment line (starting with #) or parameter parname not found
//         n .... number of numbers read
extern int inputparline (const char * parname, FILE * fin_coq, float *nn)
;

// return random number between 0 and z
extern  float rnd(float z);

// return integer of floating number (like basic integer function)
extern  float integer (float s);

//return rounded integer of floating number
extern int cint (float s);

// round to 1e-11 precision
extern double myround(double s);
// round to  precision
extern double myround(double prec,double s);

// return threej symbol 
extern float threej (float AJ1,float AJ2,float AJ3,float AM1,float AM2,float AM3);

//extern int factorial(int number);  // Too slow. Use lookup table, and inline
const int facts[13] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800}; 
// For long int (64-bit numbers), use this: 
// const long int facts[21] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800,
// 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000,
// 355687428096001,  6402373705727994, 121645100408832080, 2432902008176640000};
inline
int factorial(int number) {
   if(number<0||number>13) { fprintf(stderr,"factorial(%i) not handled\n",number); exit(1); }
   return facts[number];
// return (number>13)?0:facts[number];  // Fast but if number>13 gives wrong results!
}
extern double factorial(double number);

template <typename T> int signum(T val) {
    return (T(0) < val) - (val < T(0));
}

#ifndef __linux__
// return rounded integer (not needed any more in MINGW 3.1.3)
//extern double rint(double value);

//extern double copysign(double a,double b);
extern int sleep(int a);
#endif

// some vector functions
void xproduct(Vector & result,Vector a, Vector b);

// calculate transition matrix element of eigenvectr <i|op|j> with eigenvectors
// given as column vectors i and j of Matrix (zr + i zi), the Hermitian operator op is described
// by   matrix op.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].
double matelr (int i,int j,Matrix & zr, Matrix & zi, Matrix & op);
double mateli (int i,int j,Matrix & zr, Matrix & zi, Matrix & op);

// calculate some column i of ComplexMatrix opZ which is defined as the 
// product of two complex matrices op * Z, however Z is given as
// as real matrices  Matrix zr + i Matrix zi, the Hermitian operator op is described
// by   Matrix op.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].
void  opZcol(int i,ComplexMatrix & opZ, Matrix & op,Matrix & zr, Matrix & zi);


// some matrix functions for hermitian matrices in
// real notation: The real parts of the elements must be
 //  stored in the lower triangle of z,the imaginary parts (of the elements
 //  corresponding to the lower triangle) in the positions
 //  of the upper triangle of z[lo..hi,lo..hi].
ComplexMatrix toStandard(Matrix A); // transform to conventional matrix notation
Matrix herm_dirprod(Matrix  A, Matrix  B); // direct product
Matrix herm_prod(Matrix  A, Matrix  B); // normal herm product (AB+BA)/2
double aMb_real(Matrix & M, Matrix & zr,Matrix & zc, int ia, int ib,Matrix*V);            // transition matrix element
double aMb_imag(Matrix & M, Matrix & zr,Matrix & zc, int ia, int ib,Matrix*V);            // <a|M|b>  a,b are columns ia and ib
complex<double> aMb_complex(zsMat<double>&M,Matrix&zr,Matrix&zc,int ia,int ib);  // of zr+izc
double aMb_real(zsMat<double> & M, Matrix & zr,Matrix & zc, int ia, int ib);
double aMb_real(zsMat<double> & M, ComplexMatrix & zc, int ia, int ib);
complex<double> aMb_complex(zsMat<double> & M, ComplexMatrix & zc, int ia, int ib);


Matrix MatrixfromVectors(Vector & v1,Vector & v2,Vector & v3);

// calculate reciprocal lattice rezi from real lattice ri
void rezcalc(Vector r1,Vector  r2,Vector  r3,Vector  rez1,Vector  rez2,Vector  rez3);

void get_abc_in_ijk(Matrix & abc_in_ijk,Vector & abc);
 // get lattice vectors in terms of coordinate system ijk
 //defined by  j||b, k||(a x b) and i normal to k and j

void dadbdc2ijk(Matrix & rijk,Matrix & r, Vector & abc);
// transforms primitive lattice vector matrix r given in terms of abc
 // to ijk coordinate system

void dadbdc2ijk(Vector & rijk,Vector & dadbdc, Vector & abc);
// transforms vector r given in terms of abc
 // to ijk coordinate system

void ijk2dadbdc(Vector & dadbdc,Vector & rijk, Vector & abc);
// transforms vector r given in terms of ijk coordinates
// to da db dc

void hkl2ijk(Vector & qijk,Vector & hkl, Vector & abc);
// transforms Miller indices (in terms of reciprocal lattice abc*)
// to Q vector in ijk coordinate system

void ijk2hkl(Vector & hkl, Vector & qijk,Vector & abc);
// transforms Q vector in ijk coordinate system to Miller indices (in terms of reciprocal lattice abc*)


void nlimits_calc(Vector & nmin, Vector & nmax, double radius, Matrix & a);
// problem: we want to find all lattice vectors Rn=ni*ai which are within a
 // sphere of radius r from the origin (ai = column vectors of matrix a)
 // this routine returns the maximum and minimum values of ni i=1,2,3
 // by probing the corners of a cube



void set_zlm_constants(Matrix & cnst);
#endif

