#include<vector.h>
#include"martin.h"

#ifndef MYEV
#define MYEV


// subs to be able to check and directly diagonalize hermitean
// matrizes, inverse a nearly singular matrix
void myPrintMatrix(FILE * file,ComplexMatrix & M,const char * s=NULL);
void myPrintMatrix(ComplexMatrix & M,const char * s=NULL);
int  myReadMatrix(FILE * file, ComplexMatrix & M);
//  print/read a complex Hermitian matrix z.The real parts of the elements must be/are
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].
void myPrintComplexMatrix(FILE * file,Matrix & M,const char * s=NULL);
void myPrintComplexMatrix(Matrix & M,const char * s=NULL);
int  myReadComplexMatrix(FILE * file,Matrix & M);

void myPrintVector(FILE * file,ComplexVector & M,const char * s=NULL);
void myPrintVector(ComplexVector & M,const char * s=NULL);
void myPrintMatrix(FILE * file,Matrix & M,const char * s=NULL);
void myPrintMatrix(Matrix & M,const char * s=NULL);
void myPrintVector(FILE * file,Vector & M,const char * s=NULL);
void myPrintVector(Vector & M,const char * s=NULL);
void myPrintComplexNumber(FILE * file,complex<double> & M,const char * s=NULL);
void myPrintComplexNumber(complex<double> & M,const char * s=NULL);
void myEigenValuesHermitean (ComplexMatrix & M,Vector & lambda,int & sort,int & maxiter);
void myEigenSystemHermitean (ComplexMatrix & M,Vector & lambda,ComplexMatrix & l,int & sort,int & maxiter);
int  myEigenSystemHermiteanGeneral (ComplexMatrix& a, ComplexMatrix& b, Vector & e, ComplexMatrix & T, int & sort, int & maxiter);

#endif
