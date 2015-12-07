#include<vector.h>
#include"martin.h"

#ifndef MYEV
#define MYEV


// subs to be able to check and directly diagonalize hermitean
// matrizes, inverse a nearly singular matrix
void myPrintComplexMatrix(FILE * file,ComplexMatrix & M);
//  print a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].
void myPrintComplexMatrix(FILE * file,Matrix & M); 
void myPrintComplexVector(FILE * file,ComplexVector & M);
int  myReadComplexMatrix(FILE * file, ComplexMatrix & M);
void myPrintMatrix(FILE * file,Matrix & M);
void myPrintVector(FILE * file,Vector & M);
void myPrintComplexNumber(FILE * file,complex<double> & M);
void myEigenValuesHermitean (ComplexMatrix & M,Vector & lambda,int & sort,int & maxiter);
void myEigenSystemHermitean (ComplexMatrix & M,Vector & lambda,ComplexMatrix & l,int & sort,int & maxiter);
int  myEigenSystemHermiteanGeneral (ComplexMatrix& a, ComplexMatrix& b, Vector & e, ComplexMatrix & T, int & sort, int & maxiter);

#endif
