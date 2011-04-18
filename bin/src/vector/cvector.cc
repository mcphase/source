/*-----------------------------------------------------------------------------*\
| implementation of the double precision complex                     cvector.cc |
| vector class of MatPack.                                                      |
|                                                                               |
| Last change: Apr 7, 1998							|
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

//----------------------------------------------------------------------------//
// error messages for vector and matrix functions
//----------------------------------------------------------------------------//

static const char *NonConformVector = "non conformant complex matrix or vector";
static const char *AllocationFail   = "ComplexVector: allocation failure";
static complex<double> Zero(0,0);


//----------------------------------------------------------------------------//
// Notes:
//----------------------------------------------------------------------------//
// All inline functions defined in the following really  m u s t  be
// inlined (choose the appropriate compiler options) - otherwise
// the vector package will be slow. Do it by hand if your compiler
// is incapable of doing it !!!
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
// complex vector and matrix classes auxilliary inlines
//----------------------------------------------------------------------------//

inline void ComplexVector::checkdim (const ComplexVector& A)
//
// Checks if the dimensions of the given vector matches this vector
//
{
    if (A.cl != cl || A.ch != ch) Matpack.Error(NonConformVector);
}

void checkdim (const ComplexVector& A, const ComplexVector& B)
//
// Checks if the dimensions of the given vectors match
//
{
    if (A.cl != B.cl || A.ch != B.ch) Matpack.Error(NonConformVector);
}


//----------------------------------------------------------------------------//
// vector class constructors and destructors
//----------------------------------------------------------------------------//

ComplexVector::ComplexVector (void)
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
    temporary = 0;
    form = MpTextFormat;
}  
                 
//----------------------------------------------------------------------------//

ComplexVector::ComplexVector (int ncl, int nch)
// 
//  Constructor:
//  Allocates a complex vector[ncl..nch]
//
{
    complex<double>* v;

    // set the dimension information
    cl = ncl; ch = nch; ncol = nch-ncl+1;

    // allocate the vector structure
    if ((v = new complex<double>[ncol]))
      V = v - ncl;
    else
      Matpack.Error(AllocationFail);
 
    // set the reference count
    D = new Reference;
    temporary = 0;
    form = MpTextFormat;
}               
    
//----------------------------------------------------------------------------//

ComplexVector::ComplexVector (int ncl, int nch, complex<double> value)
// 
//  Constructor:
//  Allocates a complex vector[ncl..nch] and initializes
//  all elements with the given value
//
{
    complex<double>* v;

    // set the dimension information
    cl = ncl; ch = nch; ncol = nch-ncl+1;

    // allocate the vector structure
    if ((v = new complex<double>[ncol])) {
	V = v - ncl;
	copyval(v,value,ncol);
    } else
	Matpack.Error(AllocationFail);

    // set the reference count
    D = new Reference;
    temporary = 0;
    form = MpTextFormat;
}                   
    
//----------------------------------------------------------------------------//

ComplexVector::ComplexVector (Vector& re)
//
// cast a double vector to a complex vector
//
{
    // set dimension information
    cl = re.Lo(); ch = re.Hi(); ncol = ch-cl+1;

    complex<double>* v;
    double* d;
    int csize = ncol;

    // allocate the vector structure
    if ((v = new complex<double>[csize])) {
	V = v - cl;
	d = &re[cl];
	while (csize--) *v++ = *d++;
    } else
	Matpack.Error(AllocationFail);

    // set the reference count
    D = new Reference;
    temporary = 0;
    form = MpTextFormat;

    if (unbound(re)) re.Vector::~Vector();
}
    
//----------------------------------------------------------------------------//

ComplexVector::ComplexVector (Vector& re, Vector& im)
//
// make a complex vector using a vector of real and a 
// vector of imaginary parts
//
{
    // set dimension information
    cl = re.Lo();  ch = re.Hi(); ncol = ch-cl+1;

    // check for conformant vectors
    if (im.Lo() != cl || im.Hi() != ch) Matpack.Error(NonConformVector);

    complex<double> *v;
    double *d, *e;
    int csize = ncol;

    // allocate the vector structure
    if ((v = new complex<double>[csize])) {
	V = v - cl;
	d = &re[cl];
	e = &im[cl];
	while (csize--) *v++ = complex<double>(*d++,*e++);
    } else
	Matpack.Error(AllocationFail);

    // set the reference count
    D = new Reference;
    temporary = 0;
    form = MpTextFormat;

    if (unbound(re)) re.Vector::~Vector();
    if (unbound(im)) im.Vector::~Vector();
}
    
//----------------------------------------------------------------------------//

ComplexVector::ComplexVector (const ComplexVector &A)
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
    
    // increase the reference count
    A.addref();
    
    // copy return value flag
    temporary = A.temporary;
    form = A.form;
}                   
    
//----------------------------------------------------------------------------//

ComplexVector::~ComplexVector (void) 
//
// Destructor:
// Removes a double precision vector from the heap
//
{ 
    // decrease the reference count and delete if neccessary
    if (D) {
	if ( (--(D->count) == 0 && temporary == 0) || D->count < 0) {
	    delete[] (V + cl);
	    delete D;
	    D = 0;   // this is crucial !!
	    V = 0;
	}
	temporary = 0;
    }
}

//----------------------------------------------------------------------------//
// resizing
//----------------------------------------------------------------------------//

void ComplexVector::Remove (void)
//
// Explicitly removes a vector from the memory
//
{
    if (D) { 

	// reset reference count and call destructor
	D->count = 0;
	this->ComplexVector::~ComplexVector();

	// set the dimension information (upper must be less than lower)
	cl = 1; ncol = ch = 0;
    } 
}

//----------------------------------------------------------------------------//

void ComplexVector::Resize (int ncl, int nch)
//
// Resizes a vector copying the old contents padding the new elements 
// with zero
//
{
    complex<double> *v;
    int csize = nch-ncl+1;
    
    // allocate a new vector 
    if ((v = new complex<double>[csize])) {
	
	if (D) {  // check if already initialized
    
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
	    temporary = 0;
	}
	
	// link new vector
	V = v - ncl;
	
    } else
	Matpack.Error(AllocationFail);   

    // set the new dimension information
    cl = ncl; ch = nch;  ncol = nch-ncl+1;
}


//----------------------------------------------------------------------------//
// initializing a vector from an array 
//----------------------------------------------------------------------------//

ComplexVector& ComplexVector::operator << (const complex<double>* src)
//
// Initializes the vector from the values given in a static array
// The number of elements given in the array must match the vector size !
//
{
    copyvec(V+cl,src,ncol);
    return *this;
}

//----------------------------------------------------------------------------//
// scalar value assignment operators 
//----------------------------------------------------------------------------//

ComplexVector& ComplexVector::operator = (complex<double> value)
{
    copyval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//

ComplexVector& ComplexVector::operator += (complex<double> value)
{
    addval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//

ComplexVector& ComplexVector::operator -= (complex<double> value)
{
    subval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//

ComplexVector& ComplexVector::operator *= (complex<double> value)
{
    mulval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//

ComplexVector& ComplexVector::operator /= (complex<double> value)
{
    divval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//
// vector value assignment operators
//----------------------------------------------------------------------------//

ComplexVector& ComplexVector::operator = (const ComplexVector& A)
{
    complex<double> *v;
    
    if (unbound(A)) {        // right side is returned from a function

	if (this->Empty()) {  // left side is not yet initialized
	    cl = A.cl;
	    ch = A.ch;
	    ncol = A.ncol;

	} else {             // free left side
	    checkdim(A);     // assign only compatible vectors
	    this->ComplexVector::~ComplexVector();
	}

	// link right side to left side
	V = A.V; 
	D = A.D;
	D->count = 1;
	((ComplexVector&)A).temporary = 0;
	((ComplexVector&)A).D = 0;
  
    } else {	       // right side is bound to a variable -> copy elements

	if (this->Empty()) {   // allocate
	    cl = A.cl;
	    ch = A.ch;
	    ncol = ch-cl+1;
			      
	    // allocate the vector structure
	    if ((v = new complex<double>[ncol])) {
		V = v - cl;
	    } else
	      Matpack.Error(AllocationFail);

	    D = new Reference;
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

ComplexVector& ComplexVector::operator += (const ComplexVector& A)
{
    checkdim(A);
    addvec(V+cl,A.V+cl,ncol);
    if (unbound(A)) ((ComplexVector&)A).ComplexVector::~ComplexVector();
    return *this;
}

//----------------------------------------------------------------------------//

ComplexVector& ComplexVector::operator -= (const ComplexVector& A)
{
    checkdim(A);
    subvec(V+cl,A.V+cl,ncol);
    if (unbound(A)) ((ComplexVector&)A).ComplexVector::~ComplexVector();
    return *this;
}

//----------------------------------------------------------------------------//

ComplexVector& ComplexVector::operator %= (const ComplexVector& A)
//
// array multiplication - element by element
//
{
    checkdim(A);
    mulvec(V+cl,A.V+cl,ncol);
    if (unbound(A)) ((ComplexVector&)A).ComplexVector::~ComplexVector();
    return *this;
}


//----------------------------------------------------------------------------//
// Subvector extraction 
//----------------------------------------------------------------------------//

ComplexVector ComplexVector::operator () (int lo, int hi) const
//
// The elements of this vector within the index range [lo..hi] 
// are returned in a vector with the corresponding dimension [lo..hi].
//
{
    // check for valid subvector range
    if (lo < cl || hi > ch)
      Matpack.Error("ComplexVector::operator(): subvector index out of range (%d,%d)", lo,hi);
    ComplexVector W(lo,hi);    
    copyvec(W.V+lo,V+lo,W.ncol);
    return W.Value();
}

  
//----------------------------------------------------------------------------//
// comparison operators 
//----------------------------------------------------------------------------//

int operator == (const ComplexVector& A, const ComplexVector& B)
{
    int retval;
    checkdim(A,B);
    retval = cmpvec(A.V+A.cl,B.V+B.cl,A.ncol);
    if (unbound(A)) ((ComplexVector&)A).ComplexVector::~ComplexVector();
    if (unbound(B)) ((ComplexVector&)B).ComplexVector::~ComplexVector();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (const ComplexVector& A, complex<double> value)
{
    int retval;
    retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((ComplexVector&)A).ComplexVector::~ComplexVector();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (complex<double> value, const ComplexVector& A)
{
    int retval;
    retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((ComplexVector&)A).ComplexVector::~ComplexVector();
    return retval;
}

//----------------------------------------------------------------------------//

int operator != (const ComplexVector& A, const ComplexVector& B)
{
    int retval;
    checkdim(A,B);
    retval = cmpvec(A.V+A.cl,B.V+B.cl,A.ncol);
    if (unbound(A)) ((ComplexVector&)A).ComplexVector::~ComplexVector();
    if (unbound(B)) ((ComplexVector&)B).ComplexVector::~ComplexVector();
    return (! retval);
}

//----------------------------------------------------------------------------//

int operator != (const ComplexVector& A, complex<double> value)
{
    int retval;
    retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((ComplexVector&)A).ComplexVector::~ComplexVector();
    return (! retval);
}

//----------------------------------------------------------------------------//

int operator != (complex<double> value, const ComplexVector& A)
{
    int retval;
    retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((ComplexVector&)A).ComplexVector::~ComplexVector();
    return (! retval);
}

//----------------------------------------------------------------------------//

int ComplexVector::operator ! () const
//
// Negation operator, equivalent to (A == 0)
// Returns 1, if A is equal to 0, 1 otherwise.
//
{
    int retval;
    retval = cmpval(V+cl,Zero,ncol);
    if (unbound(*this)) ((ComplexVector*)this)->ComplexVector::~ComplexVector();
    return retval;
}

//----------------------------------------------------------------------------//

ComplexVector::operator int () const
//
// Make  if(A){...}  equivalent to  if(A != 0){...}.
//
{
    int retval;
    retval =  ! cmpval(V+cl,Zero,ncol);
    if (unbound(*this)) ((ComplexVector*)this)->ComplexVector::~ComplexVector();
    return retval;
}


//----------------------------------------------------------------------------//
// arithmetic operators
//----------------------------------------------------------------------------//

ComplexVector operator + (const ComplexVector& A, complex<double> value)
{
    ComplexVector V(A.cl,A.ch);
    addval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator + (complex<double> value, const ComplexVector& A)
{
    ComplexVector V(A.cl,A.ch);
    addval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator + (const ComplexVector& A, const ComplexVector& B)
{
    checkdim(A,B);
    ComplexVector V(A.cl,A.ch);
    addvec(V.V+A.cl,A.V+A.cl,B.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator % (const ComplexVector& A, const ComplexVector& B)
//
// array multiplication - element by element
//
{
    checkdim(A,B);
    ComplexVector V(A.cl,A.ch);
    mulvec(V.V+A.cl,A.V+A.cl,B.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator - (const ComplexVector& A)
//
// unary minus sign
//
{
    ComplexVector V(A.cl,A.ch);
    negvec(V.V+V.cl,A.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator - (const ComplexVector& A, complex<double> value)
{
    ComplexVector V(A.cl,A.ch);
    subval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator - (complex<double> value, const ComplexVector& A)
{
    ComplexVector V(A.cl,A.ch);
    subval(V.V+V.cl,value,A.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator - (const ComplexVector& A, const ComplexVector& B)
{
    checkdim(A,B);
    ComplexVector V(A.cl,A.ch);
    subvec(V.V+A.cl,A.V+A.cl,B.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator * (const ComplexVector& A, complex<double> value)
{                  
    ComplexVector V(A.cl,A.ch);
    mulval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator * (complex<double> value, const ComplexVector& A)
{
    ComplexVector V(A.cl,A.ch);
    mulval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

complex<double> operator * (const ComplexVector& A, const ComplexVector& B)
//
// Computes the scalar product of two complex vectors
// that is Sum(i){ a_i * conj(b_i) }
//
{
    checkdim(A,B);
    complex<double> sum = dotvec(A.V+A.cl,B.V+B.cl,A.ncol);
    return sum;
}

//----------------------------------------------------------------------------//

ComplexVector operator / (const ComplexVector& A, complex<double> value)
{                  
    ComplexVector V(A.cl,A.ch);
    divval(V.V+V.cl,A.V+A.cl,value,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexVector operator / (complex<double> value, const ComplexVector& A)
{                  
    ComplexVector V(A.cl,A.ch);
    divval(V.V+V.cl,value,A.V+A.cl,A.ncol);
    return V.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix operator ^ (const ComplexVector& A, const ComplexVector& B)
//
// Computes the dyadic product of the column vector A 
// with the row vector B.
//
{
    int nrl = A.cl, nrh = A.ch, ncl = B.cl, nch = B.ch;
    int csize = B.ncol;
    int rsize = A.ncol;
    ComplexMatrix T(nrl,nrh,ncl,nch);
    const complex<double> *a = A.V+nrl,
                          *b = B.V+ncl,
                          *bp; 
    complex<double> **t = &T[nrl], *tp;
    int n;  
    while (rsize--) {
	bp = b;
	tp = *t++ + ncl; 
	n = csize;
	while (n--) *tp++ = conj(*bp++) * *a; 
	a++;
    }
    return T.Value();
}

//----------------------------------------------------------------------------//
// friend functions
//----------------------------------------------------------------------------//

ComplexMatrix Diagonal (const ComplexVector& d)
//
// Make a diagonal matrix using the elements of the given vector.
//
{
    int cl = d.cl;
    int ch = d.ch;
    ComplexMatrix M(cl,ch,cl,ch, Zero);
    complex<double> **m = M.M;
    const complex<double>  *v = d.V;
    for (int i = cl; i <= ch; i++) m[i][i] = v[i];
    return M.Value();
}
    
//----------------------------------------------------------------------------//

ComplexMatrix ComplexDiagonal (const Vector& d)
//
// Make a diagonal matrix using the elements of the given vector.
//
{
  int cl = d.cl;
  int ch = d.ch;
  ComplexMatrix M(cl,ch,cl,ch, Zero);
  complex<double> **m = M.M;
  const double   *v = d.V;
  for (int i = cl; i <= ch; i++) m[i][i] = v[i];
  return M.Value();
}

//----------------------------------------------------------------------------//
// complex<double> Sum (const ComplexVector& c)
// The sum of all elements.
//----------------------------------------------------------------------------//

complex<double> Sum (const ComplexVector& c)
{
  int n = c.ncol;
  if (n <= 0) Matpack.Error("complex<double> Sum (const ComplexVector& c) -- empty vector"); 
  complex<double> sum = 0;
  const complex<double>* src = c.V+c.cl;
  for (int i = 0; i < n; i++) sum += src[i];
  return sum;
}

//----------------------------------------------------------------------------//
// double Norm (const ComplexVector& A)
// Returns the Euclidean norm (2-norm) of the vector, that is the sqrt of the sum of
// modulus squared of all vector elements.
//----------------------------------------------------------------------------//

double Norm (const ComplexVector& A)
{
  int n = A.ncol;
  if (n <= 0) Matpack.Error("double Norm (const ComplexVector& A) -- empty vector");  
  double sum = sqrt(norm2(A.Store(),n));
  return sum;
}

//----------------------------------------------------------------------------//
// double Norm2 (const ComplexVector& A)
// Returns the Euclidean norm (2-norm) of the vector, that is the sum of 
// modulus squared of all vector elements.
//----------------------------------------------------------------------------//

double Norm2 (const ComplexVector& A)
{
  int n = A.ncol;
  if (n <= 0) Matpack.Error("double Norm2 (const ComplexVector& A) -- empty vector");  
  double sum = norm2(A.Store(),n);
  return sum;
}

//----------------------------------------------------------------------------//
// double Norm1 (const ComplexVector& c)
// The 1-norm of a vector is the sum of absolute values of all elements.
//----------------------------------------------------------------------------//

double Norm1 (const ComplexVector& c)
{
  int n = c.ncol;
  if (n <= 0) Matpack.Error("double Norm1 (const ComplexVector& c) -- empty vector");
  double sum = 0;
  const complex<double>* src = c.V+c.cl;
  for (int i = 0; i < n; i++) sum += abs(src[i]);
  return sum;
}

//----------------------------------------------------------------------------//
// double NormInf (const ComplexVector& d)
// The infinity-norm is the maximum absolute value of all elements.
//----------------------------------------------------------------------------//

double NormInf (const ComplexVector& d)
{
  int n = d.ch-d.cl;
  if (n < 0) Matpack.Error("double NormInf (const ComplexVector& d) -- empty vector");
  const complex<double>* src = d.V+d.cl;
  double m = abs(*src++);
  double h;
  while (n--) {
    h = abs(*src++);
    if (h > m) m = h; 
  }
  return m;
}

//----------------------------------------------------------------------------//
// complex<double> Min (const ComplexVector& c)
// Returns the element with the smallest modulus
//----------------------------------------------------------------------------//

complex<double> Min (const ComplexVector& c)
{
  int n = c.ncol-1;
  if (n < 0) Matpack.Error("complex<double> Min (const ComplexVector& c) -- empty vector");
  const complex<double> *src  = c.V+c.cl,
                        *sp   = src,
                        *minp = src;
  double tmp, minnorm = norm(*minp);
  while (n--) {
    sp++; 
    if ((tmp = norm(*sp)) < minnorm) {
      minp = sp; 
      minnorm = tmp;
    }
  }
  return *minp;
}

//----------------------------------------------------------------------------//
// complex<double> Min (const ComplexVector& c, int& i)
// Returns the element with the smallest modulus and its index
//----------------------------------------------------------------------------//

complex<double> Min (const ComplexVector& c, int& i)
{
  int n = c.ncol-1;
  if (n < 0) Matpack.Error("complex<double> Min (const ComplexVector& c, int& i) -- empty vector");
  const complex<double> *src  = c.V+c.cl,
                        *sp   = src,
                        *minp = src;
  double tmp, minnorm = norm(*minp);
  while (n--) {
    sp++; 
    if ((tmp = norm(*sp)) < minnorm) {
      minp = sp; 
      minnorm = tmp;
    }
  }
  i = int(minp-src)+c.cl;
  return *minp;
}

//----------------------------------------------------------------------------//
// complex<double> Max (const ComplexVector& c)
// Returns the element with the largest modulus
//----------------------------------------------------------------------------//

complex<double> Max (const ComplexVector& c)
{
  int n = c.ncol-1;
  if (n < 0) Matpack.Error("complex<double> Max (const ComplexVector& c) -- empty vector");
  const complex<double> *src  = c.V+c.cl,
                        *sp   = src,
                        *maxp = src;
  double tmp, maxnorm = norm(*maxp);
  while (n--) {
    sp++; 
    if ((tmp = norm(*sp)) > maxnorm) {
      maxp = sp; 
      maxnorm = tmp;
    }
  }
  return *maxp;
}

//----------------------------------------------------------------------------//
// complex<double> Max (const ComplexVector& c, int& i)
// Returns the element with the largest modulus and its index
//----------------------------------------------------------------------------//

complex<double> Max (const ComplexVector& c, int& i)
{
  int n = c.ncol-1;
  if (n < 0) Matpack.Error("complex<double> Max (const ComplexVector& c, int& i) -- empty vector");
  const complex<double> *src  = c.V+c.cl,
                        *sp   = src,
                        *maxp = src;
  double tmp, maxnorm = norm(*maxp);
  while (n--) {
    sp++; 
    if ((tmp = norm(*sp)) > maxnorm) {
      maxp = sp; 
      maxnorm = tmp;
    }
  }
  i = int(maxp-src)+c.cl;
  return *maxp;
}

//----------------------------------------------------------------------------//

static int comp_Complex (const void* d1, const void* d2)
{
    double n1 = norm(*(complex<double>*)d1);
    double n2 = norm(*(complex<double>*)d2);
    if (n1 < n2) return -1;
    else if (n1 > n2) return 1;
    else return 0;
}

void Sort (ComplexVector& d)
//
// Sort vector inplace by ascending modulus order using the system
// qsort() function
//
{
    unsigned size = d.ncol;
    complex<double> *data = d.V + d.cl;
    qsort(data,size,sizeof(complex<double>),comp_Complex);
}

//----------------------------------------------------------------------------//

void Reverse (ComplexVector& d)
//
// Reverse the vector inplace
//
{
    unsigned size = d.ncol / 2;
    complex<double> *d1 = d.V+d.cl;
    complex<double> *d2 = d.V+d.ch;
    complex<double> t;
    while (size--) {
	t = *d1;
	*d1++ = *d2;
	*d2-- = t;
    }
}

//----------------------------------------------------------------------------//

int MatchingIndexRange (const ComplexVector &a, const ComplexVector &b)
//
// Return True if vectors have the same index range, False otherwise.
//
{
    return (a.cl == b.cl && a.ch == b.ch);
}
 
//----------------------------------------------------------------------------//
// member functions 
//----------------------------------------------------------------------------//

ComplexVector ComplexVector::Apply (complex_mapper fcn) const
//
// Applies the given function to all elements
//
{
    ComplexVector A(cl,ch);
    int n = ncol;
    complex<double>* src = V+cl;
    complex<double>* dst = A.V+cl;
    for (int i = 0; i < n; i++) dst[i] = fcn(src[i]);
    return A.Value();
}

//----------------------------------------------------------------------------//

ComplexVector ComplexVector::Apply (complex_vector_index_mapper fcn) const
//
// Applies the given function to all elements
// The function takes the index and the value as arguments.
//
{
    ComplexVector A(cl,ch);
    complex<double>* src = V;
    complex<double>* dst = A.V;
    for (int i = cl; i <= ch; i++) dst[i] = fcn(i,src[i]);
    return A.Value();
}

//----------------------------------------------------------------------------//

void ComplexVector::Set (complex<double> value)
//
// Sets all vector elements to the given value.
// The vector will be changed inplace.
//
{
    copyval(V+cl,value,ncol);
}

//----------------------------------------------------------------------------//

void ComplexVector::Set (complex_mapper fcn)
//
// Sets all vector elements by applying the given function.
// The vector will be changed inplace.
//
{
    int n = ncol;
    complex<double>* src = V+cl;
    for (int i = 0; i < n; i++) src[i] = fcn(src[i]);
}
 
//----------------------------------------------------------------------------//

void ComplexVector::Set (complex_vector_index_mapper  fcn)
//
// Sets all vector elements by applying the given function.
// The function takes the index and the value as arguments.
// The vector will be changed inplace.
//
{
    complex<double>* src = V;
    for (int i = cl; i <= ch; i++) src[i] = fcn(i,src[i]);
}
 
//----------------------------------------------------------------------------//
   
ComplexVector& ComplexVector::ShiftIndex (int by)
//
// Shift the index range of the vector by the given offset
// The vector will be changed inplace.
//
{
    if (this->Empty()) 
      Matpack.Error("ComplexVector::ShiftIndex: vector is not allocated");
    
    // shift pointers and dimension information
    V -= by; cl += by; ch += by;

    return *this;
}

//----------------------------------------------------------------------------//
// elementwise elementary functions
//----------------------------------------------------------------------------//

ComplexVector ComplexVector::Conjugate (void) const
//
// Returns the complex conjugate of this vector
//
{
    ComplexVector B(cl,ch);
    int n = ncol;
    const complex<double>* a = V+cl;
    complex<double>* b = B.V+cl;
    for (int i = 0; i < n; i++) b[i] = conj(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Real (const ComplexVector& A)
//
// Returns the vector of the real parts of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = real(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Imag (const ComplexVector& A)
//
// Returns the vector of the imaginary parts of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = imag(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Abs (const ComplexVector& A)
//
// Returns the vector of the moduli of the elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = abs(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector AbsSqr (const ComplexVector& A)
//
// Returns the vector of the modulus squared elements of A
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = norm(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

Vector Arg (const ComplexVector& A)
//
// Returns the vector of the arguments of the complex elements of A
// We don't use the arg() standard function from complex.h because
// it is not save when called with a Zero argument. Use the MatPack
// inline function Angle() instead, which handles the Zero case:
// the argument of complex zero is 0 by definition.
//
{
    Vector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a= A.V+A.cl;
    double* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = Angle(real(a[i]),imag(a[i]));
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexVector Cos (const ComplexVector& A)
//
// Returns the vector of the cosines of the elements of A
//
{
    ComplexVector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    complex<double>* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = cos(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexVector Sin (const ComplexVector& A)
//
// Returns the vector of the sines of the elements of A
//
{
    ComplexVector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    complex<double>* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = sin(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexVector Cosh (const ComplexVector& A)
//
// Returns the vector of the hyperbolic cosines of the elements of A
//
{
    ComplexVector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    complex<double>* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = cosh(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexVector Sinh (const ComplexVector& A)
//
// Returns the vector of the hyperbolic sines of the elements of A
//
{
    ComplexVector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    complex<double>* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = sinh(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexVector Exp (const ComplexVector& A)
//
// Returns the vector of the  exponentials of the elements of A
//
{
    ComplexVector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    complex<double>* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = exp(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexVector Log (const ComplexVector& A)
//
// Returns the vector of the natural logarithms of the elements of A
//
{
    ComplexVector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    complex<double>* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = log(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexVector Sqrt (const ComplexVector& A)
//
// Returns the vector of the square roots of the elements of A
//
{
    ComplexVector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    complex<double>* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = sqrt(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexVector Sqr (const ComplexVector& A)
//
// Returns the vector of the squares of the elements of A
//
{
    ComplexVector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    complex<double>* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = sqr(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexVector Cube (const ComplexVector& A)
//
// Returns the vector of the cubes of the elements of A
//
{
    ComplexVector B(A.cl,A.ch);
    int n = A.ncol;
    const complex<double>* a = A.V+A.cl;
    complex<double>* b = B.V+B.cl;
    for (int i = 0; i < n; i++) b[i] = cube(a[i]);
    return B.Value();
}

//----------------------------------------------------------------------------//
