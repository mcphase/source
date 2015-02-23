/*-----------------------------------------------------------------------------*\
| implementation of the double precision                             dvector.cc |
| vector class of MatPack.                                                      |
|                                                                               |
| Last change: Apr 6, 1998							|
|                                                                               |
| Matpack Library Release 1.3                                                   |
| Copyright (C) 1991-1998 by Berndt M. Gammel. All rights reserved              |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its  entirety and its | 
| documentation for  non-commercial purpose and  without fee is hereby granted, | 
| provided that this license information and copyright notice appear unmodified | 
| in all copies. This software is provided 'as is'  without  express or implied | 
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
#include <cstdlib>
#include <cstdarg>

//----------------------------------------------------------------------------//
// prototypes for external routines
//----------------------------------------------------------------------------//

void sort (int n, double*, double *);
void sort (int n, double*, double *, double *);

//----------------------------------------------------------------------------//
// error messages for vector and matrix functions
//----------------------------------------------------------------------------//

static const char *NonConformVector = "non conformant matrix or vector";
static const char *AllocationFail   = "Vector: allocation failure";
static double Zero = 0.0;

//----------------------------------------------------------------------------//
// Notes:
//----------------------------------------------------------------------------//
// All inline functions defined in the following really  m u s t  be
// inlined (choose the appropriate compiler options) - otherwise
// the vector package will be slow. Do it by hand if your compiler
// is incapable of doing it !!!
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
// double vector and matrix classes auxilliary inlines
//----------------------------------------------------------------------------//

inline void Vector::checkdim (const Vector& A)
//
// Checks if the dimensions of the given vector matches this vector
//
{
    if (A.cl != cl || A.ch != ch) Matpack.Error(NonConformVector);
}

//----------------------------------------------------------------------------//

void checkdim (const Vector& A, const Vector& B)
//
// Checks if the dimensions of the given vectors match
//
{
    if (A.cl != B.cl || A.ch != B.ch) Matpack.Error(NonConformVector);
}

//----------------------------------------------------------------------------//
// define a null vector (empty vector)
//----------------------------------------------------------------------------//

Vector NullVector;


//----------------------------------------------------------------------------//
// vector class constructors and destructors
//----------------------------------------------------------------------------//

Vector::Vector (void)
// 
//  Constructor:
//  Define a vector which is not yet allocated
//
{
    // set the dimension information (upper must be less than lower)
    cl = 1; ncol = ch = 0;

    // reset the dynamic part
    V = 0;
    D = 0;
    isempty = true;
    temporary = 0;
    form = MpTextFormat;
}                   

//----------------------------------------------------------------------------//

Vector::Vector (int ncl, int nch)
// 
//  Constructor:
//  Allocates a double Vector[ncl..nch]
//
{
    double* v;

    // set the dimension information
    cl = ncl; ch = nch; ncol = nch-ncl+1;

    // allocate the vector structure
    if ((v = new double[ncol])) 
      V = v - ncl;
    else
      Matpack.Error(AllocationFail);
 
    // set the reference count
    D = new Reference;
    isempty = false;
    temporary = 0;
    form = MpTextFormat;
}                   

//----------------------------------------------------------------------------//

Vector::Vector (int ncl, int nch, double value)
// 
//  Constructor:
//  Allocates a double vector[ncl..nch] and initializes
//  all elements with the given value
//
{
    double* v;

    // set the dimension information
    cl = ncl; ch = nch; ncol = nch-ncl+1;

    // allocate the vector structure
    if ((v = new double[ncol])) {
	V = v - ncl;
	copyval(v,value,ncol);
    } else
	Matpack.Error(AllocationFail);

    // set the reference count
    D = new Reference;
    isempty = false;
    temporary = 0;
    form = MpTextFormat;
}                   

//----------------------------------------------------------------------------//

Vector::Vector (int ncl, int nch, double d0, double d1, ...)
// 
//  Constructor:
//  Allocates a double Vector[ncl..nch] and initializes it by the given elements
//
{
    double* v;

    // set the dimension information
    cl = ncl; ch = nch; ncol = nch-ncl+1;
    
    // allocate the vector structure
    if ((v = new double[ncol])) 
      V = v - ncl;
    else
      Matpack.Error(AllocationFail);
    
    // set the reference count
    D = new Reference;
    isempty = false;
    temporary = 0;
    form = MpTextFormat;
    
    // initialize vector
    v[0] = d0;                  // static part
    v[1] = d1;
    va_list ap;                 // variable part
    va_start( ap, d1 );
    for (int i = 2; i < ncol; i++)
      v[i] = va_arg( ap, double );
    va_end( ap );
}

//----------------------------------------------------------------------------//

Vector::Vector (const Vector &A)
// 
//  Copy constructor:
//  adds a reference to the block, is called automatically when
//  returning a vector from a function !
//
{
    // copy the dimension information
    cl = A.cl;  ch = A.ch;  ncol = A.ncol;

    // copy link to block
    V = A.V;
    D = A.D;
    isempty = A.isempty;
    
    // increase the reference count
    A.addref();
    
    // copy return value flag
    temporary = A.temporary;
    form = A.form;
}                   

//----------------------------------------------------------------------------//

Vector::~Vector (void) 
//
// Destructor:
// Removes a double precision vector from the heap
//
{ 
    // decrease the reference count and delete if neccessary
    if (!isempty) {
        if(!D) { isempty = true; } else
	if ( (--(D->count) == 0 && temporary == 0) || D->count < 0) {
	    delete[] (V + cl);
            isempty = true;
	    delete D;
            isempty = true;
	    D = 0;   // this is crucial !!
	    V = 0;
	}
	temporary = 0;    
    }
}


// ----------------------------- resizing -------------------------- //

void Vector::Remove (void)
//
// Explicitly removes a vector from the memory
//
{
    if (!isempty) { 

	// reset reference count and call destructor
	D->count = 0;
	this->Vector::~Vector();

	// set the dimension information (upper must be less than lower)
	cl = 1;	ncol = ch = 0;
    } 
}

//----------------------------------------------------------------------------//

void Vector::Resize (int ncl, int nch)
//
// Resizes a vector copying the old contents padding the new elements 
// with zero
//
{
    double *v;
    int csize = nch-ncl+1;
    
    // check for identical resize
    if (ncl == cl && nch == ch) return;

    // allocate a new vector 
    if ((v = new double[csize])) {
	
	if (!isempty) {  // check if already initialized
    
            // compute range to copy
            int Lo = MpMax(ncl,cl);
            int lo = Lo-ncl;
            int hi = MpMin(nch,ch)-ncl;
	    
	    // copy the old contents and pad with zero
	    if ((csize = cl-ncl)  > 0) copyval(v,Zero,csize);
	    if ((csize = hi-lo+1) > 0) copyvec(v+lo,V+Lo,csize);
	    if ((csize = nch-ch)  > 0) copyval(v+hi+1,Zero,csize);

	    // remove old vector
	    delete (V+cl);

	} else {  // newly initialized vector
	    
	    // initialize to zero
	    copyval(v,Zero,csize);
	    
	    // create  the reference block
	    D = new Reference;
            isempty = false;
	    temporary = 0;
	}
	
	// link new vector
	V = v - ncl;
	
    } else
	Matpack.Error(AllocationFail);   

    // set the new dimension information
    cl = ncl; ch = nch; ncol = nch-ncl+1;
}

// ------------- initializing a vector from an array --------------- //


Vector& Vector::operator << (const double* src)
//
// Initializes the vector from the values given in a static array
// The number of elements given in the array must match the vector size !
//
{
    copyvec(V+cl,src,ncol);
    return *this;
}

// --------------- scalar value assignment operators --------------- //

Vector& Vector::operator = (double value)
{
    copyval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//

Vector& Vector::operator += (double value)
{
    addval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//

Vector& Vector::operator -= (double value)
{
    subval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//

Vector& Vector::operator *= (double value)
{
    mulval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//

Vector& Vector::operator /= (double value)
{
    divval(V+cl,value,ncol);
    return *this;
}


// --------------- vector value assignment operators --------------- //

Vector& Vector::operator = (const Vector& A)
{
    double *v;
    
    if (unbound(A)) {        // right side is returned from a function

	if (this->Empty()) {  // left side is not yet initialized
	    cl = A.cl;
	    ch = A.ch;
	    ncol = A.ncol;

	} else {             // free left side
	    checkdim(A);     // assign only compatible vectors
	    this->Vector::~Vector();
	}

	// link right side A to left side
	V = A.V; 
	D = A.D;
        isempty = A.isempty;
	D->count = 1;
	((Vector&)A).temporary = 0;
	((Vector&)A).D = 0;
			     
    } else {	        // right side is bound to a variable -> copy elements

	if (this->Empty()) {   // allocate
	    cl = A.cl;
	    ch = A.ch;
	    ncol = ch-cl+1;
			      
	    // allocate the vector structure
	    if ((v = new double[ncol])) {
		V = v - cl;
	    } else
	        Matpack.Error(AllocationFail);

	    D = new Reference;
            isempty = false;
	} else
	    checkdim(A);

	// copy right side to left side
        copyvec(V+cl,A.V+cl,ncol);

        // set reference count	      
        D->count = 1; 
    }

    // set return value flag
    temporary = 0;
    return *this;
}

//----------------------------------------------------------------------------//

Vector& Vector::operator += (const Vector& A)
{
    checkdim(A);
    addvec(V+cl,A.V+cl,ncol);
    if (unbound(A)) ((Vector&)A).Vector::~Vector();
    return *this;
}

//----------------------------------------------------------------------------//

Vector& Vector::operator -= (const Vector& A)
{
    checkdim(A);
    subvec(V+cl,A.V+cl,ncol);
    if (unbound(A)) ((Vector&)A).Vector::~Vector();
    return *this;
}

//----------------------------------------------------------------------------//

Vector& Vector::operator %= (const Vector& A)
//
// array multiplication - element by element
//
{
    checkdim(A);
    mulvec(V+cl,A.V+cl,ncol);
    if (unbound(A)) ((Vector&)A).Vector::~Vector();
    return *this;
}


//----------------------------------------------------------------------------//
// Subvector extraction 
//----------------------------------------------------------------------------//

Vector Vector::operator () (int lo, int hi) const
//
// The elements of this vector within the index range [lo..hi] 
// are returned in a vector with the corresponding dimension [lo..hi].
//
{
    // check for valid subvector range
    if (lo < cl || hi > ch)
      Matpack.Error("Vector::operator(): subvector index out of range (%d,%d)", lo,hi);
    Vector W(lo,hi);    
    copyvec(W.V+lo,V+lo,W.ncol);
    return W.Value();
}


// ---------------------- comparison operators --------------------- //

int operator == (const Vector& A, const Vector& B)
{
    int retval;
    checkdim(A,B);
    retval = cmpvec(A.V+A.cl,B.V+B.cl,A.ncol);
    if (unbound(A)) ((Vector&)A).Vector::~Vector();
    if (unbound(B)) ((Vector&)B).Vector::~Vector();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (const Vector& A, double value)
{
    int retval;
    retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((Vector&)A).Vector::~Vector();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (double value, const Vector& A)
{
    int retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((Vector&)A).Vector::~Vector();
    return retval;
}

//----------------------------------------------------------------------------//

int operator != (const Vector& A, const Vector& B)
{
    int retval;
    checkdim(A,B);
    retval = cmpvec(A.V+A.cl,B.V+B.cl,A.ncol);
    if (unbound(A)) ((Vector&)A).Vector::~Vector();
    if (unbound(B)) ((Vector&)B).Vector::~Vector();
    return (! retval);
}

//----------------------------------------------------------------------------//

int operator != (const Vector& A, double value)
{
    int retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((Vector&)A).Vector::~Vector();
    return (! retval);
}

//----------------------------------------------------------------------------//

int operator != (double value, const Vector& A)
{
    int retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((Vector&)A).Vector::~Vector();
    return (! retval);
}

//----------------------------------------------------------------------------//

int Vector::operator ! () const
//
// Negation operator, equivalent to (A == 0)
// Returns 1, if A is equal to 0, 1 otherwise.
//
{
    int retval = cmpval(V+cl,Zero,ncol);
    if (unbound(*this)) ((Vector*)this)->Vector::~Vector();
    return retval;
}

//----------------------------------------------------------------------------//

Vector::operator int () const
//
// Make if(A){...}  equivalent to  if(A != 0){...}.
//
{
    int retval;
    retval =  ! cmpval(V+cl,Zero,ncol);
    if (unbound(*this)) ((Vector*)this)->Vector::~Vector();
    return retval;
}


// --------------------- arithmetic operators ---------------------- //

Vector operator + (const Vector& A, double value)
{
    Vector V(A.cl,A.ch);
    addval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator + (double value, const Vector& A)
{
    Vector V(A.cl,A.ch);
    addval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator + (const Vector& A, const Vector& B)
{
    checkdim(A,B);
    Vector V(A.cl,A.ch);
    addvec(V.V+A.cl,A.V+A.cl,B.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator % (const Vector& A, const Vector& B)
//
// array multiplication - element by element
//
{
    checkdim(A,B);
    Vector V(A.cl,A.ch);
    mulvec(V.V+A.cl,A.V+A.cl,B.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator - (const Vector& A)
//
// unary minus sign
//
{
    Vector V(A.cl,A.ch);
    negvec(V.V+V.cl,A.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator - (const Vector& A, double value)
{
    Vector V(A.cl,A.ch);
    subval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator - (double value, const Vector& A)
{
    Vector V(A.cl,A.ch);
    subval(V.V+V.cl,value,A.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator - (const Vector& A, const Vector& B)
{
    checkdim(A,B);
    Vector V(A.cl,A.ch);
    subvec(V.V+A.cl,A.V+A.cl,B.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator * (const Vector& A, double value)
{
    Vector V(A.cl,A.ch);
    mulval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator * (double value, const Vector& A)
{
    Vector V(A.cl,A.ch);
    mulval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

double operator * (const Vector& A, const Vector& B)
//
// Computes the scalar product
//
{
    checkdim(A,B);
    double sum = dotvec(A.V+A.cl,B.V+B.cl,A.ncol);
    return sum;
}

//----------------------------------------------------------------------------//

Vector operator / (const Vector& A, double value)
{                  
    Vector V(A.cl,A.ch);
    divval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Vector operator / (double value, const Vector& A)
{                  
    Vector V(A.cl,A.ch);
    divval(V.V+V.cl,value,A.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

Matrix operator ^ (const Vector& A, const Vector& B)
//
// Computes the dyadic product of the column vector A with the
// row vector B.
//
{
    int nrl = A.cl, nrh = A.ch, ncl = B.cl, nch = B.ch;
    int csize = B.ncol;
    int rsize = A.ncol;
    Matrix T(nrl,nrh,ncl,nch);
    double *a = A.V+nrl, aval;
    double *b = B.V+ncl;
    double **t = &T[nrl];
    double *tp,*bp;   
    int n;  
    while (rsize--) {
	bp = b;
	tp = *t++ + ncl; 
	n = csize;
	aval = *a;
	while (n--) *tp++ = *bp++ * aval; 
	a++;
    }
    return T.Value();
}


//----------------------------------------------------------------------------//
// friend functions
//----------------------------------------------------------------------------//

Matrix Diagonal (const Vector& d)
//
// Make a diagonal matrix using the elements of the given vector.
//
{
    int cl = d.cl;
    int ch = d.ch;
    Matrix M(cl,ch,cl,ch, Zero);
    double **m = M.M;
    double  *v = d.V;
    for (int i = cl; i <= ch; i++) m[i][i] = v[i];
    return M.Value();
}

//----------------------------------------------------------------------------//
// double Norm (const Vector &A)
// The euclidean norm (2-norm) of a vector is the square root of the sum of 
// squares of all elements.
//----------------------------------------------------------------------------//

double Norm (const Vector &A)
{
  int n = A.Elements();
  if (n <= 0) Matpack.Error("double Norm (const Vector &A) -- empty vector");
  double sum = sqrt(norm2(A.Store(),n)+1e-300);
  return sum;
}

//----------------------------------------------------------------------------//
// double Norm2 (const Vector &A)
// The euclidean norm (2-norm) of a vector is the square root of the sum of 
// squares of all elements.
//----------------------------------------------------------------------------//

double Norm2 (const Vector &A)
{
  int n = A.Elements();
  if (n <= 0) Matpack.Error("double Norm2 (const Vector &A) -- empty vector");
  double sum = sqrt(norm2(A.Store(),n));
  return sum;
}

//----------------------------------------------------------------------------//
// double Norm1 (const Vector &d)
// The 1-norm of a vector is the sum of absolute values of all elements.
//----------------------------------------------------------------------------//

double Norm1 (const Vector &d)
{
  int n = d.ncol;
  if (n <= 0) Matpack.Error("double Norm1 (const Vector &A) -- empty vector");
  double sum = 0;
  const double* src = d.V+d.cl;
  for (int i = 0; i < n; i++) sum += fabs(src[i]);
  return sum;
}

//----------------------------------------------------------------------------//
// double NormInf (const Vector& d)
// The infinity-norm is the maximum absolute value of all elements.
//----------------------------------------------------------------------------//

double NormInf (const Vector& d)
{
  int n = d.ncol;
  if (n <= 0) Matpack.Error("double NormInf (const Vector& A) -- empty vector");
  const double* src = d.V+d.cl;
  double m = fabs(src[0]);
  for (int i = 1; i < n; i++) {
    double h = fabs(src[i]);
    if (h > m) m = h; 
  }
  return m;
}

//----------------------------------------------------------------------------//
// double Sum (const Vector& d)
// The sum of all elements.
//----------------------------------------------------------------------------//

double Sum (const Vector& d)
{
  int n = d.ncol;
  if (n <= 0) Matpack.Error("double Sum (const Vector& A) -- empty vector");
  double sum = 0;
  const double* src = d.V+d.cl;
  for (int i = 0; i < n; i++) sum += src[i];
  return sum;
}

//----------------------------------------------------------------------------//
// double Min (const Vector& d)
// Returns the smallest element
//----------------------------------------------------------------------------//

double Min (const Vector& d)
{
  int n = d.ncol-1;
  if (n < 0) Matpack.Error("double Min (const Vector& A) -- empty vector");
  const double *src = d.V+d.cl, 
               *sp  = src,
               *sm  = src;
  while (n--) {
    sp++; 
    if (*sp < *sm) sm = sp; 
  }
  return *sm;
}

//----------------------------------------------------------------------------//
// double Min (const Vector& d, int& i)
// Returns the smallest element and its index
//----------------------------------------------------------------------------//

double Min (const Vector& d, int& i)
{
  int n = d.ncol-1;
  if (n < 0) Matpack.Error("double Min (const Vector& A, int& i) -- empty vector");
  const double *src = d.V+d.cl, 
               *sp  = src,
               *sm  = src;
  while (n--) {
    sp++; 
    if (*sp < *sm) sm = sp; 
  }
  i = int(sm-src)+d.cl;
  return *sm;
}

//----------------------------------------------------------------------------//
// double Max (const Vector& d)
// Returns the largest element of the vector
//----------------------------------------------------------------------------//

double Max (const Vector& d)
{
  int n = d.ncol-1;
  if (n < 0) Matpack.Error("double Max (const Vector& A) -- empty vector");
  const double *src = d.V+d.cl, 
               *sp  = src,
               *sm  = src;
  while (n--) {
    sp++; 
    if (*sp > *sm) sm = sp; 
  }
  return *sm;
}

//----------------------------------------------------------------------------//
// double Max (const Vector& d, int& i)
// Returns the largest element of the vector and its index
//----------------------------------------------------------------------------//

double Max (const Vector& d, int& i)
{
  int n = d.ncol-1;
  if (n < 0) Matpack.Error("double Max (const Vector& A, int& i) -- empty vector");
  const double *src = d.V+d.cl, 
               *sp  = src,
               *sm  = src;
  while (n--) {
    sp++; 
    if (*sp > *sm) sm = sp; 
  }
  i = int(sm-src)+d.cl;
  return *sm;
}

//----------------------------------------------------------------------------//
// void Sort (Vector& d)
// Sort vector inplace into ascending numerical order using the system
// qsort() function. The auxilliary comp_double() below is used.
//----------------------------------------------------------------------------//

static int comp_double (const void* d1, const void* d2)
{
  if ( *(double*)d1 < *(double*)d2 ) return -1;
  else if ( *(double*)d1 > *(double*)d2 ) return 1;
  else return 0;
}

void Sort (Vector& d)
{
  unsigned size = d.ncol;
  if (size < 2) return;
  double *data = d.V + d.cl;       // qsort expects 0-offset !
  qsort(data,size,sizeof(double),comp_double);
}

//----------------------------------------------------------------------------//
// void Sort (Vector& d1, Vector& d2)
// Sort vector d1 inplace into ascending numerical order using, 
// while making the corresponding rearrangements in the vector d2.
// Uses the MatPack function: void sort(int,double*,double*);
// d1 and d2 must have an equal number of elements, 
// but need not have the same index range (offset).
//----------------------------------------------------------------------------//

void Sort (Vector& d1, Vector& d2)
{
    unsigned size = d1.ncol;

    if (d2.ncol != int(size)) 
      Matpack.Error("void Sort(Vector& d1, Vector& d2): "
	   "vector sizes %d and %d are incompatible", size, d2.ncol);

    double *data1 = d1.V + d1.cl - 1; // sort expects 1-offset
    double *data2 = d2.V + d2.cl - 1; // sort expects 1-offset
    sort(size,data1,data2);
}

//----------------------------------------------------------------------------//

void Sort (Vector& d1, Vector& d2, Vector& d3)
//
// Sort vector d1 inplace into ascending numerical order using, 
// while making the corresponding rearrangements in the vectors d2 and d3.
// Uses the MatPack function: void sort(int,double*,double*,double*);
// d1, d2, and d3 must have an equal number of elements, 
// but need not have the same index range (offset).
//
{
  int size = d1.Elements();

  // check arguments
  if (d2.Elements() != size || d3.Elements() != size) 
    Matpack.Error("void Sort(Vector& d1, Vector& d2, Vector& d3): "
		  "vector sizes %d, %d, %d are incompatible", 
		  size, d2.Elements(), d3.Elements());

  // nothing to do
  if (size <= 1) return;

  // sort expects 1-offset data vectors
  double *data1 = d1.Store() - 1,
         *data2 = d2.Store() - 1,
         *data3 = d3.Store() - 1;

  sort(size,data1,data2,data3);
}

//-----------------------------------------------------------------------------//

void Reverse (Vector& d)
//
// Reverse the vector inplace
//
{
    unsigned size = d.ncol / 2;  // repl by signed  nt !
    double *d1 = d.V+d.cl;
    double *d2 = d.V+d.ch;
    double t;
    while (size--) {
	t = *d1;
	*d1++ = *d2;
	*d2-- = t;
    }
}

//----------------------------------------------------------------------------//

int MatchingIndexRange (const Vector &a, const Vector &b)
//
// Return True if vectors have the same index range, False otherwise.
//
{
    return (a.cl == b.cl && a.ch == b.ch);
}


//----------------------------------------------------------------------------//
// member functions
//----------------------------------------------------------------------------//

void Vector::Set (double value)
//
// Sets all vector elements to the given value.
// The vector will be changed inplace.
//
{
    copyval(V+cl,value,ncol);    
}

//----------------------------------------------------------------------------//

void Vector::Set (double_mapper fcn)
//
// Sets all vector elements by applying the given function.
// The vector will be changed inplace.
//
{
    int n = ncol;
    double* src = V+cl;
    for (int i = 0; i < n; i++) src[i] = fcn(src[i]);
}

//----------------------------------------------------------------------------//

void Vector::Set (double_vector_index_mapper fcn)
//
// Sets all vector elements by applying the given function.
// The function takes the index and the value as arguments.
// The vector will be changed inplace.
//
{
    double* src = V;
    for (int i = cl; i <= ch; i++) src[i] = fcn(i,src[i]);
}

//----------------------------------------------------------------------------//

Vector& Vector::ShiftIndex (int by)
//
// Shift the index range of the vector by the given offset
// The vector will be changed inplace.
//
{
    if (this->Empty()) 
      Matpack.Error("Vector::ShiftIndex: vector is not allocated");
    
    // shift pointers and dimension information
    V -= by; cl += by; ch += by;

    return *this;
}

//----------------------------------------------------------------------------//

Vector Vector::Apply (double_mapper fcn) const
//
// Applies the given function to all elements
//
{
    Vector A(cl,ch);
    int n = ncol;
    double* src = V+cl;
    double* dst = A.V+cl;
    while (n--) *dst++ = (*fcn)(*src++);
    return A.Value();
}

//----------------------------------------------------------------------------//

Vector Vector::Apply (double_vector_index_mapper fcn) const
//
// Applies the given function to all elements
//
{
    Vector A(cl,ch);

    int c;
    double* src = V+cl;
    double* dst = A.V+cl;

    for (c = cl; c <= ch; c++)
      *dst++ = (*fcn)(c,*src++);

    return A.Value();
}

//----------------------------------------------------------------------------//
// elementwise elementary functions
//----------------------------------------------------------------------------//

Vector Cos (const Vector& A)
//
// Returns the vector of the cosines of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    double* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = cos(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Sin (const Vector& A)
//
// Returns the vector of the sines of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    double* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = sin(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Cosh (const Vector& A)
//
// Returns the vector of the hyperbolic cosines of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    double* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = cosh(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Sinh (const Vector& A)
//
// Returns the vector of the hyperbolic sines of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    double* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = sinh(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Exp (const Vector& A)
//
// Returns the vector of the  exponentials of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    double* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = exp(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Log (const Vector& A)
//
// Returns the vector of the natural logarithms of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    double* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = log(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Sqrt (const Vector& A)
//
// Returns the vector of the square roots of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    double* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = sqrt(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Sqr (const Vector& A)
//
// Returns the vector of the squares of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    double* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = sqr(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Cube (const Vector& A)
//
// Returns the vector of the cubes of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    double* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = cube(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//
