/*-----------------------------------------------------------------------------*\
| implementation of the double precision                             dmatrix.cc |
| matrix class of MatPack.                                                      |
|                                                                               |
| Last change: Apr 6, 1998							|
|                                                                               |
| Matpack Library Release 1.3                                                   |
| Copyright (C) 1991-1998 by Berndt M. Gammel. All rights reserved              |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.						|
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.							|
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.			|
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration	|
| and installation of Matpack.							|
|                                                                               |
\*-----------------------------------------------------------------------------*/

#include "vector.h"
#include "vecinl.h"

#ifdef _MATPACK_USE_BLAS_
#include "blas.h"
#endif

#include <cstdlib>
#include <cstdarg>

//----------------------------------------------------------------------------//

// error messages for vector and matrix functions
static const char *NonSquareMatrix  = "non square double matrix";
static const char *NonConformMatrix = "non conformant matrix or vector";
static double Zero = 0.0;

//----------------------------------------------------------------------------//
// Notes:
//----------------------------------------------------------------------------//
// All inline functions defined in the following really  m u s t  be
// inlined (choose the appropriate compiler options) - otherwise
// the matrix package will be slow. Do it by hand if your compiler
// is incapable of doing it !!!
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// define a null matrix (empty matrix)
//----------------------------------------------------------------------------//

Matrix NullMatrix;


//----------------------------------------------------------------------------//
// macros to make notation more concise
//----------------------------------------------------------------------------//

#define UT(A) (A.attribute == UpperTriangular)
#define LT(A) (A.attribute == LowerTriangular)
#define GM(A)  (A.attribute == General)


//----------------------------------------------------------------------------//
// LogAndSign class
//----------------------------------------------------------------------------//

void LogAndSign::operator *= (double x)
{
   if (x > Zero) { log_value += log(x); }
   else if (x < Zero) { log_value += log(-x); sign = -sign; }
   else sign = 0;
}

//----------------------------------------------------------------------------//

LogAndSign::LogAndSign (double x)
{
   if (x == Zero) { log_value = Zero; sign = 0; return; }
   else if (x < Zero) { sign = -1; x = -x; }
   else sign = 1;
   log_value = log(x);
}


//----------------------------------------------------------------------------//
// double vector and matrix classes auxilliary inlines
//----------------------------------------------------------------------------//

inline void Matrix::checkdim (const Matrix& A)
//
// Checks if the dimensions of the matries match when adding
//
{
    if (A.cl != cl || A.ch != ch || A.rl != rl || A.rh != rh) 
      Matpack.Error(NonConformMatrix);
}

//----------------------------------------------------------------------------//

void checkdim (const Matrix& A, const Matrix& B)
//
// Checks if the dimensions of the matrices match when adding
//
{
    if (A.cl != B.cl || A.ch != B.ch || A.rl != B.rl || A.rh != B.rh) 
      Matpack.Error(NonConformMatrix);
}

//----------------------------------------------------------------------------//

void checksquare (const Matrix &M)
{
    if (M.rl != M.cl || M.rh != M.ch) Matpack.Error(NonSquareMatrix);
}

//----------------------------------------------------------------------------//

static double** newmat (int nrl, int nrh, int ncl, int nch)
//
// allocate the dynamic part of a matrix on the heap, 
// allocate matrix in one block to faciliate many matrix functions !
//
{
    double* v;
    double **M;
    double** m;
    int rsize = nrh-nrl+1; 
    int csize = nch-ncl+1;

    if ((m = new double*[rsize])) {
	M = m - nrl;
	if ((v = new double[csize*rsize])) {
	    for (v -= ncl; rsize--; v += csize)  *m++ = v;
	    return M;
	}
    }
    Matpack.Error("Matrix(%d,%d,%d,%d): newmat: allocation failure",nrl,nrh,ncl,nch);
    return (double**) 0;
}


//----------------------------------------------------------------------------//
// matrix class constructors and destructors
//----------------------------------------------------------------------------//

Matrix::Matrix (void)
// 
//  Constructor:
//  Define a matrix which is not yet allocated.
//
{
    // reset the dimension information (upper must be less than lower)
    rl = cl = 1;  ncol = nrow = rh = ch = 0;

    // reset the dynamic part 
    M = 0;
    D = 0;
    isempty = true;
}

//----------------------------------------------------------------------------//

Matrix::Matrix (int nrl, int nrh, int ncl, int nch)
// 
//  Constructor:
//  Allocates a double Matrix[nrl..nrh,ncl..nch] 
//
{
    // set the dimension information
    rl = nrl; rh = nrh;
    cl = ncl; ch = nch;
    nrow = nrh-nrl+1;
    ncol = nch-ncl+1;

    if ( nrow > 0 && ncol > 0) {
      // allocate the matrix structure
      M = newmat(nrl,nrh,ncl,nch);
      // set the reference count
      D = new Reference;
      isempty = false;
    } else {
      M = 0; D = 0;
      isempty = true;
    }
}     

//----------------------------------------------------------------------------//

Matrix::Matrix (int nrl, int nrh, int ncl, int nch, double value)
// 
//  Constructor:
//  Allocates a double Matrix[nrl..nrh,ncl..nch] and initializes 
//  all elements with the given value
//
{
    // set the dimension information
    rl = nrl; rh = nrh;
    cl = ncl; ch = nch;
    nrow = nrh-nrl+1;
    ncol = nch-ncl+1;

    if ( nrow > 0 && ncol > 0) {
      // allocate the matrix structure
      M = newmat(nrl,nrh,ncl,nch);
      // set all elements
      copyval(M[nrl]+ncl,value,ncol*nrow);
      // set the reference count
      D = new Reference;
      isempty = false;
    } else {
      M = 0; D = 0;
      isempty = true;
    }
}     

//----------------------------------------------------------------------------//

Matrix::Matrix (int nrl, int nrh, int ncl, int nch, double d0, double d1, ...)
// 
//  Constructor:
//  Allocates a double Matrix[nrl..nrh,ncl..nch] and initializes 
//  all elements with the given values
//
{
    // set the dimension information
    rl = nrl; rh = nrh;
    cl = ncl; ch = nch;
    nrow = nrh-nrl+1;
    ncol = nch-ncl+1;

    if ( nrow > 0 && ncol > 0) {
      // allocate the matrix structure
      M = newmat(nrl,nrh,ncl,nch);
      // set the reference count
      D = new Reference;
      isempty = false;
      // initialize matrix
      int nelem = ncol*nrow;
      double *v = M[nrl]+ncl;
      v[0] = d0;                  // static part
      v[1] = d1;
      va_list ap;                 // variable part
      va_start( ap, d1 );
      for (int i = 2; i < nelem; i++)
	v[i] = va_arg( ap, double );
      va_end( ap );
    } else {
      M = 0; D = 0;
      isempty = true;
    }
}

//----------------------------------------------------------------------------//

Matrix::Matrix (const Matrix& A)
// 
//  Copy constructor:
//  adds a reference to the block, is called automatically when
//  returning a matrix from a function !
//
{
    // set attribute
    attribute = A.attribute;

    // copy the dimension information
    rl = A.rl; rh = A.rh; 
    cl = A.cl; ch = A.ch;
    nrow = A.nrow; 
    ncol = A.ncol;

    // copy link 
    M = A.M;
    D = A.D;
    isempty = A.isempty;
    
    // increase the reference count
    A.addref();

    // copy return value flag
    temporary =  A.temporary;
    form = A.form;
}

//----------------------------------------------------------------------------//

Matrix::~Matrix (void)
//
//  Destructor:
//  Removes a double precision matrix from the heap
//
{
    // decrease the reference count and delete if neccessary
    if (!isempty) {
        if(!D) { isempty = true; } else
	if ( (--(D->count) == 0 && temporary == 0) || D->count < 0) {

	    // free data
	    delete[] (M[rl]+cl);
	    delete[] (M+rl);

	    // free info block
            isempty = true;
	    delete D;
	    D = 0;   // this is crucial !!
	    M = 0;
	}
	temporary = 0;
    }
}

//----------------------------------------------------------------------------//
// resizing 
//----------------------------------------------------------------------------//

void Matrix::Remove (void)
//
// Explicitly removes a matrix from the memory
//
{
    if (!isempty) {
	
	// reset reference count and call destructor
	D->count = 0;
	this->Matrix::~Matrix();

	// set attribute
	attribute = General;

	// reset the dimension information (upper must be less than lower)
	rl = cl = 1; nrow = ncol = rh = ch = 0;
    }
} 

//----------------------------------------------------------------------------//

void Matrix::Resize (int nrl, int nrh, int ncl, int nch)
{
    Matpack.Error("Matrix::Resize: NOT YET IMPLEMENTED");

    double **m;
    int csize = nch-ncl+1,
        rsize = nrh-nrl+1;

    // check for identical resize
    if ((ncl == cl) && (nch == ch) && (nrl == rl) && (nrh == rh)) return;

    // allocate a new matrix
    m = newmat(nrl,nrh,ncl,nch);

    if (D) {  // check if already initialized
	
#ifdef UNFINISHED
	// compute range to copy
	int c_Lo = max(ncl,cl);
	int c_lo = c_Lo-ncl;
	int c_hi = min(nch,ch)-ncl;

	int r_Lo = max(nrl,rl);
	int r_lo = r_Lo-nrl;
	int r_hi = min(nrh,rh)-nrl;
	
	// copy the old contents and pad with zero
	if ((csize = cl-ncl)  > 0)     copyval(v,Zero,csize);
	if ((csize = c_hi-c_lo+1) > 0) copyvec(v+c_lo,V+c_Lo,csize);
	if ((csize = nch-ch)  > 0)     copyval(v+c_hi+1,Zero,csize);
#endif
	
	// remove old matrix
	delete (M[rl]+cl);
	delete (M+rl);
	
	
    } else {  // newly initialized matrix
	   
	// initialize to zero
	copyval( m[nrl]+ncl, Zero, rsize*csize );
	
	// create  the reference block
	D = new Reference;
	temporary = 0;	
    }

    // link block
    M = m;

    // set the dimension information
    rl = nrl; rh = nrh;
    cl = ncl; ch = nch;
    nrow = nrh-nrl+1;
    ncol = nch-ncl+1;
}


//----------------------------------------------------------------------------//
// initializing a matrix from an array
//----------------------------------------------------------------------------//

Matrix& Matrix::operator << (const double* src)
//
// Initializes the matrix from the values given in a static array.
// The matrix will be filled row by row
// The number of elements given in the array must match the number of
// elements in the matrix !
//
{
  copyvec(M[rl]+cl,src,ncol*nrow);
  return *this;
}


//----------------------------------------------------------------------------//
// scalar value assignment operators 
//----------------------------------------------------------------------------//

Matrix& Matrix::operator = (double value)
//
// Set to the matrix value*Id, i.e. the diagonal is set to value,
// all other elements are set to zero.
//
{
  if (ncol != nrow)
    Matpack.Error("Matrix& Matrix::operator=(double): non square matrix\n");

  double **m = M + rl;
  int i;
  copyval( M[rl]+cl, Zero, ncol*nrow );
  for (i = cl; i <= ch; i++) (*m++)[i] = value;
  return *this;
}

//----------------------------------------------------------------------------//

Matrix& Matrix::operator += (double value)
//
// Add the matrix value*Id, i.e. add value to the diagonal elements.
//
{
  if (ncol != nrow)
    Matpack.Error("Matrix& Matrix::operator+=(double): non square matrix\n");

  double **m = M + rl;
  int i;
  for (i = cl; i <= ch; i++) (*m++)[i] += value;
  return *this;
} 

//----------------------------------------------------------------------------//

Matrix& Matrix::operator -= (double value)
//
// Subtract the matrix value*Id, i.e. subtract value from the diagonal.
//
{
  if (ncol != nrow)
    Matpack.Error("Matrix& Matrix::operator-=(double): non square matrix\n");

  double **m = M + rl;
  int i;
  for (i = cl; i <= ch; i++) (*m++)[i] -= value;
  return *this;
} 

//----------------------------------------------------------------------------//

Matrix& Matrix::operator *= (double value)
//
// Multiply with the matrix value*Id, i.e. multiply all elements 
// with value.
//
{
  mulval(M[rl]+cl,value,ncol*nrow);
  return *this;
} 

//----------------------------------------------------------------------------//

Matrix& Matrix::operator /= (double value)
//
// Multiply with the matrix (1/value)*Id, i.e. divide all elements
// by value.
//
{
  divval(M[rl]+cl,value,ncol*nrow);
  return *this;
} 


//----------------------------------------------------------------------------//
// matrix value assignment operators
//----------------------------------------------------------------------------//

Matrix& Matrix::operator = (const Matrix& A)
{  
    if (unbound(A)) {        // right side is returned from a function
			     
	if (this->Empty()) {  // left side is not yet initialized
	    cl = A.cl;
	    ch = A.ch;
	    rl = A.rl;
	    rh = A.rh;
	    nrow = A.nrow;
	    ncol = A.ncol;
	} else {             // free left side with destructor
 	    checkdim(A);     // assign only compatible matrices 
	    this->Matrix::~Matrix();
	}

	// link right side to left side
	M = A.M; 
	D = A.D;
        isempty = A.isempty;
        if (D) D->count = 1;
	((Matrix&)A).temporary = 0;
	((Matrix&)A).D = 0;
	((Matrix&)A).isempty = true;
     
    } else {		// right side is bound to a variable -> copy elements

	if (this->Empty()) {   // allocate
	    cl = A.cl;
	    ch = A.ch;
	    rl = A.rl;
	    rh = A.rh;
	    nrow = A.nrow;
	    ncol = A.ncol;
	    if ( nrow > 0 && ncol > 0) {
	      M = newmat(rl,rh,cl,ch);
	      D = new Reference;
              isempty = false;
	    } else {
	      M = 0;
	      D = 0;
              isempty = true;
	    }
	} else
	    checkdim(A);

	// copy right side to left side in one loop	
	if ( nrow > 0 && ncol > 0)	
	  copyvec( M[rl]+cl, A.M[rl]+cl, ncol*nrow );

        // set reference count	      
        if (D) D->count = 1; 
    }
    
    // copy attribute
    attribute = A.attribute;

    // set return value flag
    temporary = 0;  
    return *this;
}

//----------------------------------------------------------------------------//

Matrix& Matrix::operator += (const Matrix& A)
{
  checkdim(A);

  // set attribute
  attribute = General;

  addvec( M[rl]+cl, A.M[rl]+cl, ncol*nrow );
  if (unbound(A)) ((Matrix&)A).Matrix::~Matrix();
  return *this;
} 

//----------------------------------------------------------------------------//

Matrix& Matrix::operator -= (const Matrix& A)
{
  checkdim(A);

  // set attribute
  attribute = General;

  subvec( M[rl]+cl, A.M[rl]+cl, ncol*nrow );
  if (unbound(A)) ((Matrix&)A).Matrix::~Matrix();
  return *this;
} 

//----------------------------------------------------------------------------//

Matrix& Matrix::operator %= (const Matrix& A)
//
// array multiplication - element by element
//
{
  checkdim(A);

  // set attribute
  attribute = General;

  mulvec( M[rl]+cl, A.M[rl]+cl, ncol*nrow );
  if (unbound(A)) ((Matrix&)A).Matrix::~Matrix();
  return *this;
} 

//----------------------------------------------------------------------------//
// Submatrix extraction 
//----------------------------------------------------------------------------//

Matrix Matrix::operator () (int rlo, int rhi, int clo, int chi) const
//
// The elements of this matrix within the index range [rlo..rhi,clo..chi] 
// are returned in a matrix with the corresponding dimension [rlo..rhi,clo..chi].
//
{
  // check for valid submatrix range
  if (rlo < rl || rhi > rh || clo < cl || chi > ch)
    Matpack.Error("Matrix::operator(): submatrix index out of range (%d,%d,%d,%d)", 
		  rlo,rhi,clo,chi);
  int i;
  Matrix W(rlo,rhi,clo,chi);    
  for (i = rlo; i <= rhi; i++)
    copyvec(W.M[i]+clo,M[i]+clo,W.ncol);
  return W.Value();
}


//----------------------------------------------------------------------------//
//  comparison operators
//----------------------------------------------------------------------------//

int operator == (const Matrix& A, const Matrix& B)
//
// Returns 1, if all corresponding elements of the 
// given matrices are equal, and returns 0 otherwise.
//
{
  int retval;
  checkdim(A,B);
  retval = cmpvec( A.M[A.rl]+A.cl, B.M[B.rl]+B.cl, A.ncol*A.nrow );
  if (unbound(A)) ((Matrix&)A).Matrix::~Matrix();
  if (unbound(B)) ((Matrix&)B).Matrix::~Matrix();
  return retval;
}

//----------------------------------------------------------------------------//

int operator == (const Matrix& A, double value)
//
// Returns 1, if A is equal to value*Id, i.e. if the diagonal is set 
// to value and all other elements are 0, and returns 0 otherwise.
//
{
  checksquare(A);

  int rl = A.rl;
  int size = A.nrow;
  double **a = A.M + rl;
  const double *c;
  int i;
  int retval = 1;
  int diag = 0;
  while (size) {
    c = *a++ + rl;
    i = diag++;
    while (i--) if ( *c++ != Zero ) { retval = 0; goto differ; }
    if ( *c++ != value ) { retval = 0; goto differ; }
    i = size--;
    while (--i) if ( *c++ != Zero ) { retval = 0; goto differ; }
  }
differ:
  if (unbound(A)) ((Matrix&)A).Matrix::~Matrix();
  return retval;
}

//----------------------------------------------------------------------------//

int operator == (double value, const Matrix& A)
//
// Returns 1, if A is equal to value*Id, i.e. if the diagonal is set 
// to value and all other elements are 0, and returns 0 otherwise.
//
{
  checksquare(A);

  int rl = A.rl;
  int size = A.nrow;
  double **a = A.M + rl;
  const double *c;
  int i;
  int retval = 1;
  int diag = 0;
  while (size) {
    c = *a++ + rl;
    i = diag++;
    while (i--) if ( *c++ != Zero ) { retval = 0; goto differ; }
    if ( *c++ != value ) { retval = 0; goto differ; }
    i = size--;
    while (--i) if ( *c++ != Zero ) { retval = 0; goto differ; }
  }
differ:
  if (unbound(A)) ((Matrix&)A).Matrix::~Matrix();
  return retval;
}

//----------------------------------------------------------------------------//

int operator != (const Matrix& A, const Matrix& B)
//
// Returns 1, if A is not equal to B, 0 otherwise.
//
{
  int retval;
  checkdim(A,B);
  retval = ! cmpvec( A.M[A.rl]+A.cl, B.M[B.rl]+B.cl, A.ncol*A.nrow );
  if (unbound(A)) ((Matrix&)A).Matrix::~Matrix();
  if (unbound(B)) ((Matrix&)B).Matrix::~Matrix();
  return retval;
}

//----------------------------------------------------------------------------//

int operator != (const Matrix& A, double value)
//
// Returns 1, if A is not equal to value*Id, 0 otherwise.
//
{
  checksquare(A);

  int rl = A.rl;
  int size = A.nrow;
  double **a = A.M + rl;
  const double *c;
  int i;
  int retval = 0;
  int diag = 0;
  while (size) {
    c = *a++ + rl;
    i = diag++;
    while (i--) if ( *c++ != Zero ) { retval = 1; goto differ; }
    if ( *c++ != value ) { retval = 1; goto differ; }
    i = size--;
    while (--i) if ( *c++ != Zero ) { retval = 1; goto differ; }
  }
differ:
  if (unbound(A)) ((Matrix&)A).Matrix::~Matrix();
  return retval;
}

//----------------------------------------------------------------------------//

int operator != (double value, const Matrix& A)
//
// Returns 1, if A is not equal to value*Id, 0 otherwise.
//
{
  checksquare(A);

  int rl = A.rl;
  int size = A.nrow;
  double **a = A.M + rl;
  const double *c;
  int i;
  int retval = 0;
  int diag = 0;
  while (size) {
    c = *a++ + rl;
    i = diag++;
    while (i--) if ( *c++ != Zero ) { retval = 1; goto differ; }
    if ( *c++ != value ) { retval = 1; goto differ; }
    i = size--;
    while (--i) if ( *c++ != Zero ) { retval = 1; goto differ; }
  }
differ:
  if (unbound(A)) ((Matrix&)A).Matrix::~Matrix();
  return retval;
}

//----------------------------------------------------------------------------//

int Matrix::operator ! () const
//
// Negation operator, equivalent to (A == 0)
// Returns 1, if A is equal to 0, 1 otherwise.
//
{
  int retval = cmpval(M[rl]+cl, Zero, ncol*nrow );
  if (unbound(*this)) ((Matrix*)this)->Matrix::~Matrix();
  return retval;  
}

//----------------------------------------------------------------------------//

Matrix::operator int () const
//
// Make if(A){...}  equivalent to  if(A != 0){...}.
//
{
  int retval = ! cmpval( M[rl]+cl, Zero, ncol*nrow );
  if (unbound(*this)) ((Matrix*)this)->Matrix::~Matrix();
  return retval;
}


//----------------------------------------------------------------------------//
// arithmetic operators
//----------------------------------------------------------------------------//

Matrix operator + (const Matrix& A, double value)
//
// Add matrix value*Id, i.e. add value to the diagonal.
//
{
  Matrix B; 
  B = A; // copies attribute
  B += value;
  return B.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator + (double value, const Matrix& A)
//
// Add matrix value*Id, i.e. add value to the diagonal.
//
{
  Matrix B;
  B = A; // copies attribute
  B += value;
  return B.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator + (const Matrix& A, const Matrix& B)
//
// add matrices
//
{
  checkdim(A,B);
  Matrix C(A.rl,A.rh,A.cl,A.ch);
  addvec(C.Store(),A.Store(),B.Store(),A.Elements());
  return C.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator % (const Matrix& A, const Matrix& B)
//
// array multiplication - element by element
//
{
  checkdim(A,B);
  Matrix C(A.rl,A.rh,A.cl,A.ch);
  mulvec(C.Store(),A.Store(),B.Store(),A.Elements());
  return C.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator - (const Matrix& A)
//
// unary minus sign, i.e. multiply all matrix elements by -1.
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);
  B.attribute = A.attribute;    // copy attribute
  negvec(B.Store(),A.Store(),A.Elements());
  return B.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator - (const Matrix& A, double value)
//
// Subtract matrix value*Id, i.e. subtract value from the diagonal.
//
{
  Matrix B;
  B = A; // copy attribute
  B -= value;
  return B.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator - (double value, const Matrix& A)
//
// Subtract this matrix from the matrix value*Id.
//
{
  if (A.ncol != A.nrow)
    Matpack.Error("Matrix operator-(double,Matrix&): non square matrix\n");

  Matrix B(A.rl,A.rh,A.cl,A.ch);
  B.attribute = A.attribute;                // copy attribute
  negvec(B.Store(),A.Store(),A.Elements()); // copy -A to B
  double **b = B.M + B.rl;
  for (int i = B.cl; i <= B.ch; i++) (*b++)[i] += value;
  return B.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator - (const Matrix& A, const Matrix& B)
//
// subtract matrices element by element
//
{
  checkdim(A,B);
  Matrix C(A.rl,A.rh,A.cl,A.ch);
  subvec(C.Store(),A.Store(),B.Store(),A.Elements());
  return C.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator * (const Matrix& A, double value)
//
// multiply all matrix elements by value
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);
  B.attribute = A.attribute;    // copy attribute
  mulval(B.Store(),A.Store(),value,A.Elements());
  return B.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator * (double value, const Matrix& A)
//
// multiply value by all matrix elements, same as above
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);
  B.attribute = A.attribute;    // copy attribute
  mulval(B.Store(),A.Store(),value,A.Elements());
  return B.Value();
} 

//----------------------------------------------------------------------------//

Matrix operator / (const Matrix& A, double value)
//
// divide all matrix elements by value
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);
  B.attribute = A.attribute;    // copy attribute
  divval(B.Store(),A.Store(),value,A.Elements());
  return B.Value();
}	

//----------------------------------------------------------------------------//

Matrix operator / (double value, const Matrix& A)
//
// divide value by all matrix elements
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);
  B.attribute = A.attribute;   // copy attribute
  divval(B.Store(),value,A.Store(),A.Elements());
  return B.Value();
}	

//----------------------------------------------------------------------------//

Vector operator * (const Vector& V, const Matrix& A)
//
// multiply vector from left to matrix
//
{
  if (A.rl != V.cl || A.rh != V.ch) 
    Matpack.Error("Matrix operator * (Vector& A,Matrix& B): "
		  "non conformant arguments");

  Vector B(A.cl,A.ch);
  double **a = A.M+A.rl;
  double *b  = B.V+B.cl;
  double *v  = V.V+V.cl;
  double *vp,**ap;
  double sum;
  int n,col;
  for (col = A.cl; col <= A.ch; col++) {
    vp = v;
    ap = a;
    sum = 0;
    n = V.ncol;
    while (n--) sum += *vp++ * (*ap++)[col];
    *b++ = sum;
  }
  return B.Value();
} 

//----------------------------------------------------------------------------//

Vector operator * (const Matrix& A, const Vector& B)
//
// multiply a matrix A from the right with a vector B, returning C = A * B.
//
{
  int arl = A.rl, arh = A.rh, acl = A.cl, ncol = A.ncol;

  if (acl != B.cl  || A.ch != B.ch)
    Matpack.Error("Vector operator * (Matrix& A, Vector& B): "
		  "non conformant arguments");

  Vector C(arl,arh);

  // avoid call to index operator that optimizes very badely
  const double * const *a = A.M;
  const double *b = B.Store();
  double *c = C.V;

  for (int i = arl; i <= arh; i++) 
    c[i] = dotvec(a[i]+acl,b,ncol);

  return C.Value();
} 

//----------------------------------------------------------------------------//
// friend functions
//----------------------------------------------------------------------------//

double Trace (const Matrix& d)
//
// Returns the trace of the matrix, that means the sum of 
// the diagonal elements. If the matrix is not square an error
// will result.
//
{
  if (d.ncol != d.nrow)
    Matpack.Error("double Trace(const Matrix&): non square matrix\n");

  int j;
  double **m = d.M+d.rl;
  double trace = 0;
  for (j = d.cl; j <= d.ch; j++) trace += (*m++)[j];  // TODO
  return trace;
} 

//----------------------------------------------------------------------------//
// double Norm (const Matrix& A)
// Returns the norm of the matrix A, which is defined by the 
// largest singular value of A.
//----------------------------------------------------------------------------//

double Norm (const Matrix& A)
{
  if (A.Elements() <= 0) Matpack.Error("double Norm (const Matrix& A) -- empty matrix");
  int lo = A.Rlo(), hi = A.Rhi();
  Matrix A2, V(lo,hi,lo,hi);
  Vector w(lo,hi);
  A2 = (Matrix&) A; // copy of A
  SVDecompose(A2,w,V);
  return Max(w);
}

//----------------------------------------------------------------------------//
// double Norm2 (const Matrix& A)
// Same as Norm() above. 
// Returns the norm of the matrix A, which is defined by the 
// largest singular value of A.
//----------------------------------------------------------------------------//

double Norm2 (const Matrix& A)
{
  if (A.Elements() <= 0) Matpack.Error("double Norm2 (const Matrix& A) -- empty matrix");
  int lo = A.Rlo(), hi = A.Rhi();
  Matrix A2, V(lo,hi,lo,hi);
  Vector w(lo,hi);
  A2 = (Matrix&) A; // copy of A
  SVDecompose(A2,w,V);
  return Max(w);
}

//----------------------------------------------------------------------------//
// double NormFro (const Matrix& A)
// Returns the euclidean norm (Frobenius norm, F-norm) of the matrix.
// This is the square root of the sum of squares of all matrix elements.
//----------------------------------------------------------------------------//

double NormFro (const Matrix& A)
{
  int n = A.Elements();
  if (n <= 0) Matpack.Error("double NormFro (const Matrix& A) -- empty matrix");
  double sum = sqrt(norm2(A.Store(),n));
  return sum;
} 

//----------------------------------------------------------------------------//
// double Norm1 (const Matrix& d)
// Returns the 1-norm of the matrix, that is the maximum of the sum
// of absolute values of a column.
//----------------------------------------------------------------------------//

double Norm1 (const Matrix& d)
{
  if (d.Elements() <= 0) 
    Matpack.Error("double Norm1 (const Matrix& d) -- empty matrix");
  double norm1 = 0;
  for (int j = d.cl; j <= d.ch; j++) {
    int rsize = d.nrow;
    double **m = d.M+d.rl;
    double sum = 0;
    while (rsize--) sum += fabs((*m++)[j]);
    if (sum > norm1) norm1 = sum;
  }
  return norm1;
} 

//----------------------------------------------------------------------------//

#ifdef OLDIE
double NormInf (const Matrix& d)
//
// Returns the infinity-norm of a matrix, that is the maximum of sum
// of absolute values of a row.
//
{
    double *v;
    int j;
    double **m = d.M+d.rl;
    double norminf = 0;
    double sum;
    int rsize = d.nrow;
    while (rsize--) { 
	v = *m++ + d.cl;
	j = d.ncol;
	sum = 0;
	while (j--) sum += fabs(*v++);
	if (sum > norminf) norminf = sum;
    }
    return norminf;
} 
#endif

double NormInf (const Matrix& d)
//
// Returns the infinity-norm of a matrix, that is the maximum of the sum
// of absolute values of a row.
//
{
  double **m = d.M+d.rl;
  double norminf = 0;
  int rows = d.nrow, cols = d.ncol;
  for (int i = 0; i < rows; i++) {
    const double *v = m[i] + d.cl;
    double sum = 0;
    for (int j = 0; j < cols; j++) sum += fabs(v[j]);
    if (sum > norminf) norminf = sum;
  }
  return norminf;
} 

//----------------------------------------------------------------------------//

double Sum (const Matrix &d)
//
// Returns the sum of all matrix elements.
//
{   
  double sum = Zero;
  int n = d.ncol*d.nrow;
  const double* src = d.M[d.rl]+d.cl;
  for (int i = 0; i < n; i++) sum += src[i];
  return sum;
} 

//----------------------------------------------------------------------------//
// double Min (const Matrix &d)
// Returns the smallest element of this matrix.
//----------------------------------------------------------------------------//

double Min (const Matrix &d)
{
  int n = d.ncol*d.nrow-1;
  if (n < 0) Matpack.Error("double Min (const Matrix &d) -- empty matrix");
  const double *src  = d.M[d.rl]+d.cl,
               *sp   = src,
               *minp = src;
  while (n--) {
    sp++; 
    if (*sp < *minp) minp = sp; 
  }
  return *minp;
} 

//----------------------------------------------------------------------------//
// double Min (const Matrix &d, int &i, int &k)
// Returns the smallest element of this matrix together with its index (i,k)
//----------------------------------------------------------------------------//

double Min (const Matrix &d, int &i, int &k)
{
  int n = d.ncol*d.nrow-1;
  if (n < 0) Matpack.Error("const Matrix &d, int &i, int &k -- empty matrix");
  const double *src  = d.M[d.rl]+d.cl,
               *sp   = src,
               *minp = src;
  while (n--) {
    sp++; 
    if (*sp < *minp) minp = sp; 
  }
  // calculate index
  n = int(minp-src);
  i = n / d.ncol + d.rl;
  k = n % d.ncol + d.cl;  
  return *minp;
} 

//----------------------------------------------------------------------------//
// double Max (const Matrix &d)
// Returns the largest element of this matrix.
//----------------------------------------------------------------------------//

double Max (const Matrix &d)
{
  int n = d.ncol*d.nrow-1;
  if (n < 0) Matpack.Error("double Max (const Matrix &d) -- empty matrix");
  const double *src  = d.M[d.rl]+d.cl,
               *sp   = src,
               *maxp = src;
  while (n--) {
    sp++; 
    if (*sp > *maxp) maxp = sp; 
  }
  return *maxp;
} 
//----------------------------------------------------------------------------//
// double Max (const Matrix &d, int &i, int &k)
// Returns the largest element of this matrix together with its index (i,k).
//----------------------------------------------------------------------------//

double Max (const Matrix &d, int &i, int &k)
{
  int n = d.ncol*d.nrow-1;
  if (n < 0) Matpack.Error("double Max (const Matrix &d, int &i, int &k) -- empty matrix");
  const double *src  = d.M[d.rl]+d.cl,
               *sp   = src,
               *maxp = src;
  while (n--) {
    sp++; 
    if (*sp > *maxp) maxp = sp; 
  }
  // calculate index
  n = int(maxp-src);
  i = n / d.ncol + d.rl;
  k = n % d.ncol + d.cl;  
  return *maxp;
} 

//----------------------------------------------------------------------------//

LogAndSign LogDet (const Matrix& M)
//
// Returns the logarithm of the determinant of this matrix. The value
// of the logarithm and the sign are returned as members of the 
// LogAndSign class. Retrieve the 
//
//   logarithmized value  by   retval.LogValue()
//   sign                 by   retval.Sign()
//   value                by   retval.Value()
//
// The LU decomposition is generated and the determinant is computed
// as the product of the diagonal values.  If the matrix is not square 
// or if the matrix  turns out to be singular an error will result.
//
{
  checksquare(M);

  int j;
  Matrix A;
  int sgn;
  LogAndSign det;

  // allocate a permutation table
  IntVector idx(M.rl,M.rh);

  // generate a copy of this matrix
  A = M;

  // compute the LU decomposition
  Decompose(A,idx,sgn);

  // sign from permutations in the decomposition
  det = sgn;

  // product the diagonal elements
  for (j = A.cl; j <= A.ch; j++) det *= A[j][j];

  return det;
}

//----------------------------------------------------------------------------//

double Det (const Matrix& M)
//
// Returns the determinant of this matrix.
//
{
  return LogDet(M).Value();
}


//----------------------------------------------------------------------------//
// elementwise elementary functions
//----------------------------------------------------------------------------//

Matrix Cos (const Matrix& A)
//
// Returns the matrix of the cosines of the elements of A
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);

  // copy attribute
  B.attribute = A.attribute;

  int n = A.nrow * A.ncol;
  const double *a = A.M[A.rl]+A.cl;
  double *b = B.M[B.rl]+B.cl;
  for (int i = 0; i < n; i++) b[i] = cos(a[i]);
  return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Sin (const Matrix& A)
//
// Returns the matrix of the sines of the elements of A
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);

  // copy attribute
  B.attribute = A.attribute;

  int n = A.nrow * A.ncol;
  const double *a = A.M[A.rl]+A.cl;
  double *b = B.M[B.rl]+B.cl;
  for (int i = 0; i < n; i++) b[i] = sin(a[i]);
  return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Cosh (const Matrix& A)
//
// Returns the matrix of the hyperbolic cosines of the elements of A
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);

  // copy attribute
  B.attribute = A.attribute;

  int n = A.nrow * A.ncol;
  const double *a = A.M[A.rl]+A.cl;
  double *b = B.M[B.rl]+B.cl;
  for (int i = 0; i < n; i++) b[i] = cosh(a[i]);
  return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Sinh (const Matrix& A)
//
// Returns the matrix of the hyperbolic sines of the elements of A
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);

  // copy attribute
  B.attribute = A.attribute;

  int n = A.nrow * A.ncol;
  const double *a = A.M[A.rl]+A.cl;
  double *b = B.M[B.rl]+B.cl;
  for (int i = 0; i < n; i++) b[i] = sinh(a[i]);
  return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Exp (const Matrix& A)
//
// Returns the matrix of the exponentials of the elements of A
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);

  // copy attribute
  B.attribute = A.attribute;

  int n = A.nrow * A.ncol;
  const double *a = A.M[A.rl]+A.cl;
  double *b = B.M[B.rl]+B.cl;
  for (int i = 0; i < n; i++) b[i] = exp(a[i]);
  return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Log (const Matrix& A)
//
// Returns the matrix of the natural logarithms of the elements of A
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);

  // copy attribute
  B.attribute = A.attribute;

  int n = A.nrow * A.ncol;
  const double *a = A.M[A.rl]+A.cl;
  double *b = B.M[B.rl]+B.cl;
  for (int i = 0; i < n; i++) b[i] = log(a[i]);
  return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Sqrt (const Matrix& A)
//
// Returns the matrix of the square roots of the elements of A
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);

  // copy attribute
  B.attribute = A.attribute;

  int n = A.nrow * A.ncol;
  const double *a = A.M[A.rl]+A.cl;
  double *b = B.M[B.rl]+B.cl;
  for (int i = 0; i < n; i++) b[i] = sqrt(a[i]);
  return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Sqr (const Matrix& A)
//
// Returns the matrix of the squares of the elements of A
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);

  // copy attribute
  B.attribute = A.attribute;

  int n = A.nrow * A.ncol;
  const double *a = A.M[A.rl]+A.cl;
  double *b = B.M[B.rl]+B.cl;
  for (int i = 0; i < n; i++) b[i] = sqr(a[i]);
  return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Cube (const Matrix& A)
//
// Returns the matrix of the cubes of the elements of A
//
{
  Matrix B(A.rl,A.rh,A.cl,A.ch);

  // copy attribute
  B.attribute = A.attribute;

  int n = A.nrow * A.ncol;
  const double *a = A.M[A.rl]+A.cl;
  double *b = B.M[B.rl]+B.cl;
  for (int i = 0; i < n; i++) b[i] = cube(a[i]);
  return B.Value();
}

//----------------------------------------------------------------------------//
// other functions
//----------------------------------------------------------------------------//

int MatchingIndexRange (const Matrix &a, const Matrix &b)
//
// Return True if matrices have the same index range, False otherwise.
//
{
  return (a.cl == b.cl && a.ch == b.ch && a.rl == b.rl && a.rh == b.rh);
}

//----------------------------------------------------------------------------//

Vector Pack (const Matrix &A)
//
// pack the matrix elements row by row into a vector
//
{
  int n = A.ncol*A.nrow;
  Vector V(A.cl,A.cl+n-1);
  copyvec(V.V+V.cl, A.M[A.rl]+A.cl,n);
  return V.Value();
}

//----------------------------------------------------------------------------//
// member functions
//----------------------------------------------------------------------------//

void Matrix::Set (double value)
//
// Sets all matrix elements to the given value
//
{
  copyval( M[rl]+cl, value, ncol*nrow );
}

//----------------------------------------------------------------------------//

void Matrix::Set (double_mapper fcn)
//
//  Sets all matrix elements by applying the given function
//
{
  int n = ncol*nrow;
  double* src = M[rl]+cl;
  while (n--) {*src = (*fcn)(*src); src++;}  // TODO
}

//----------------------------------------------------------------------------//

void Matrix::Set (double_matrix_index_mapper fcn)
//
//  Sets all matrix elements by applying the given function.
// The function takes the index pair and the value as arguments.
//
{
  int r,c;
  double* src = M[rl]+cl;
 
  for (r = rl; r <= rh; r++) 
    for (c = cl; c <= ch; c++) {
      *src = (*fcn)(r,c,*src); src++;
    }
}

//----------------------------------------------------------------------------//

Matrix& Matrix::ShiftIndex (int row, int col)
//
// Shift the index range of the matrix by the given offsets
// in columns and rows. The matrix will be changed inplace.
//
{
  double **v;
  int n;

  if (this->Empty()) 
    Matpack.Error("Matrix::ShiftIndex: matrix is not allocated");

  // shift pointers
  v = M+rl; 
  n = nrow;
  while (n--) { *v -= col; v++; }
  M -= row;

  // shift dimension information
  rl += row; rh += row;
  cl += col; ch += col;

  return *this;
}

//----------------------------------------------------------------------------//

Vector Matrix::Row (int i) const
//
// Returns the i-th row of this matrix as a vector.
//
{
  if (i < rl || i > rh) 
    Matpack.Error("Matrix::Row(int): index out of range (%d)\n",i); 

  Vector A(cl,ch);
  copyvec(A.V+cl,M[i]+cl,ncol);
  return A.Value();
}

//----------------------------------------------------------------------------//

Vector Matrix::Column (int i) const
//
// Returns the i-th column of this matrix as a vector.
//
{
  if (i < cl || i > ch) 
    Matpack.Error("Matrix::Column(int): index out of range (%d)\n",i); 

  Vector A(rl,rh);
  double **m = M + rl;
  double  *v = A.V + rl;
  int rsize = nrow;
  while (rsize--) *v++ = (*m++)[i];
  return A.Value();
}

//----------------------------------------------------------------------------//

Matrix Matrix::Transpose (void) const
//
// Returns the transpose of this matrix.
//
{
  Matrix A(cl,ch,rl,rh);

  if ( (attribute == UpperTriangular) && (ncol == nrow) )
    A.attribute = LowerTriangular;
  else if ( (attribute == LowerTriangular) && (ncol == nrow) )
    A.attribute = UpperTriangular;
      
  double **m = M + rl;
  double **a = A.M + cl;
  double **ap;
  const double *v;
  int j;
  int i;
  for (i = rl; i <= rh; i++) { 
    v = *m++ + cl;
    ap = a;
    j = ncol;
    while (j--) (*ap++)[i] = *v++;
  }
  return A.Value();
} 

//----------------------------------------------------------------------------//

Matrix Matrix::Apply (double_mapper fcn) const
//
// Applies the given function to all elements
//
{
  Matrix A(rl,rh,cl,ch);

  // copy attribute
  A.attribute = attribute;

  int n = ncol*nrow;
  double* src = M[rl]+cl;
  double* dst = A.M[rl]+cl;
  while (n--) *dst++ = (*fcn)(*src++);
  return A.Value();
}

//----------------------------------------------------------------------------//

Matrix Matrix::Apply (double_matrix_index_mapper fcn) const
//
// Applies the given function to all elements.
// The function takes the index pair and the value as arguments.
//
{
  Matrix A(rl,rh,cl,ch);

  // copy attribute
  A.attribute = attribute;

  int r,c;
  double* src = M[rl]+cl;
  double* dst = A.M[rl]+cl;

  for (r = rl; r <= rh; r++) 
    for (c = cl; c <= ch; c++) 
      *dst++ = (*fcn)(r,c,*src++);
    
  return A.Value();
}

//----------------------------------------------------------------------------//

Matrix Matrix::Inverse (void) const
//
// Returns the inverse of this matrix. If the matrix is not square 
// or if the matrix turns out to be singular an error will result.
//
{
  /*
    // this version need less memory but is two times slower than LU decomp.
    Matrix A2;
    A2 = *this;        // generate a copy of this matrix
    Ortho(A2);         // make inplace matrix inversion
    return A2.Value();
    */

  // inversion by LU decomposition
  Matrix Y(rl,rh,rl,rh), LU; 
  Y = 1;              // make unit matrix (diagonal 1)
  LU = *this;         // copy of this matrix - protect from overwriting
  SolveLinear(LU,Y);  // create inverse in Y
  return Y.Value();
}
 
//----------------------------------------------------------------------------//

int Matrix::Attribute (int newAttr)
// 
// Set new attribute checking for validity. 
// Returns True if successful, False if the attribute can't be set:
// For example you can't set the LowerTrianguar attribute if the 
// matrix is not square !
//
{
  int   status = true,
        isSquare = (ncol == nrow); 

  switch (newAttr) {
    case General:         
      attribute = newAttr;
      break;
    case LowerTriangular: 
      if (isSquare) attribute = newAttr; else status = false;
      break;
    case UpperTriangular:
      if (isSquare) attribute = newAttr; else status = false;
      break;
    default: 
      Matpack.Error("Matrix::Attribute: unknown attribute (%d)",newAttr);
  }
  return status;
}

//----------------------------------------------------------------------------//

Matrix FlipUpDown (const Matrix& A)
//
// Returns the matrix with rows flipped, 
// i.e. up and down is mirrored.
//
{    
  int arl = A.rl, arh = A.rh,
      acl = A.cl, ach = A.ch,
      n = arl + arh,
      k = Nint(0.5*n);

  Matrix B(arl,arh,acl,ach); // attributes changed: result is general matrix

  // avoid call to index operator that optimizes very badely
  double **a = A.M;
  double **b = B.M;

  for (int i = arl; i <= k; i++) 
    for (int j = acl; j <= ach; j++) {
      b[i][j] = a[n-i][j];
      b[n-i][j] = a[i][j];
    }

  return B.Value();
}


//----------------------------------------------------------------------------//

Matrix FlipLeftRight (const Matrix& A)
//
// Returns the matrix with columns flipped, 
// i.e. left and right is mirrored.
//
{    
  int arl = A.rl, arh = A.rh,
      acl = A.cl, ach = A.ch,
      n = acl + ach,
      k = Nint(0.5*n);

  Matrix B(arl,arh,acl,ach); // attributes changed: result is general matrix

  // avoid call to index operator that optimizes very badely
  double **a = A.M;
  double **b = B.M;

  for (int i = arl; i <= arh; i++) 
    for (int j = acl; j <= k; j++) {
      b[i][j] = a[i][n-j];
      b[i][n-j] = a[i][j];
    }

  return B.Value();
}

//----------------------------------------------------------------------------//
