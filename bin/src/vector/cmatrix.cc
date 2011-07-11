/*-----------------------------------------------------------------------------*\
| implementation of the double precision complex                     cmatrix.cc |
| matrix class of MatPack.                                                      |
|                                                                               |
| Last change: Apr 7, 1998							|
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

//----------------------------------------------------------------------------//
// error messages for vector and matrix functions
//----------------------------------------------------------------------------//

static const char *NonSquareMatrix  = "non square complex matrix";
static const char *NonConformMatrix = "non conformant complex matrix or vector";

const complex<double> Zero(0.0,0.0);

//----------------------------------------------------------------------------//
// Notes:
//----------------------------------------------------------------------------//
// All inline functions defined in the following really  m u s t  be
// inlined (choose the appropriate compiler options) - otherwise
// the matrix package will be slow. Do it by hand if your compiler
// is incapable of doing it !!!
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
// complex vector and matrix classes auxilliary inlines 
//----------------------------------------------------------------------------//

inline void ComplexMatrix::checkdim (const ComplexMatrix& A)
//
// Checks if the dimensions of the matries match when adding
//
{
    if (A.cl != cl || A.ch != ch || A.rl != rl || A.rh != rh) 
      Matpack.Error(NonConformMatrix);
}

//----------------------------------------------------------------------------//

void checkdim (const ComplexMatrix& A, const ComplexMatrix& B)
//
// Checks if the dimensions of the matrices match when adding
//
{
    if (A.cl != B.cl || A.ch != B.ch || A.rl != B.rl || A.rh != B.rh) 
      Matpack.Error(NonConformMatrix);
}

//----------------------------------------------------------------------------//


#ifndef _R_CHECKSQUARE
#define _R_CHECKSQUARE 1
void checksquare (const ComplexMatrix& M)
//
// check that matrix is square and has an identical index range in
// both dimensions.
//
{
    if (M.Rlo() != M.Clo() || M.Rhi() != M.Chi()) Matpack.Error(NonSquareMatrix);
}
#endif

//----------------------------------------------------------------------------//
// allocation primitive 
//----------------------------------------------------------------------------//

static complex<double>** newmat (int nrl, int nrh, int ncl, int nch)
//
// allocate the dynamic part of a matrix on the heap
//
{
    complex<double>* v;
    complex<double> **M;
    complex<double>** m;
    int rsize = nrh-nrl+1; 
    int csize = nch-ncl+1;

    if ((m = new complex<double>*[rsize])) {
	M = m - nrl;
	if ((v = new complex<double>[csize*rsize])) {
	    for (v -= ncl; rsize--; v += csize)  *m++ = v;
	    return M;
	}
    }
    Matpack.Error("ComplexMatrix(%d,%d,%d,%d): newmat: allocation failure",
	 nrl,nrh,ncl,nch);
    return (complex<double>**) 0;
}


//----------------------------------------------------------------------------//
// matrix class constructors and destructors 
//----------------------------------------------------------------------------//

ComplexMatrix::ComplexMatrix (void)
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
}

//----------------------------------------------------------------------------//

ComplexMatrix::ComplexMatrix (int nrl, int nrh, int ncl, int nch)
// 
//  Constructor:
//  Allocates a complex matrix[nrl..nrh,ncl..nch] 
//
{
   // set the dimension information
    rl = nrl; rh = nrh;
    cl = ncl; ch = nch;
    nrow = nrh-nrl+1;
    ncol = nch-ncl+1;

    // allocate the matrix structure
    M = newmat(nrl,nrh,ncl,nch);

    // set the reference count
    D = new Reference;
}     

//----------------------------------------------------------------------------//

ComplexMatrix::ComplexMatrix (int nrl, int nrh, int ncl, int nch, complex<double> value)
// 
//  Constructor:
//  Allocates a complex matrix[nrl..nrh,ncl..nch] and initializes 
//  all elements with the given value
//
{
    // set the dimension information
    rl = nrl; rh = nrh;
    cl = ncl; ch = nch;
    nrow = nrh-nrl+1;
    ncol = nch-ncl+1;

    // allocate the matrix structure
    M = newmat(nrl,nrh,ncl,nch);

    // set all elements
    copyval(M[nrl]+ncl,value,ncol*nrow);

    // set the reference count
    D = new Reference;
}     

//----------------------------------------------------------------------------//

ComplexMatrix::ComplexMatrix (Matrix& re)
//
// cast a double matrix to a complex matrix
//
{
    // set the dimension information
    rl = re.Rlo(); rh = re.Rhi();  
    cl = re.Clo(); ch = re.Chi();
    nrow = rh-rl+1; 
    ncol = ch-cl+1;

    // allocate the matrix structure
    M = newmat(rl,rh,cl,ch);
    
    // copy real matrix
    copyrevec(M[rl]+cl,re[rl]+cl,ncol*nrow);

    // set the reference count
    D = new Reference;

    if (unbound(re)) re.Matrix::~Matrix();
}

//----------------------------------------------------------------------------//

ComplexMatrix::ComplexMatrix (Matrix& re, Matrix& im)
//
// make a complex matrix using a matrix of real and a 
// matrix of imaginary parts
//
{
    // set the dimension information
    rl = re.Rlo(); rh = re.Rhi();   
    cl = re.Clo(); ch = re.Chi();
    nrow = rh-rl+1; 
    ncol = ch-cl+1;

    // allocate the matrix structure
    M = newmat(rl,rh,cl,ch);
    
    // copy real matrix
    copyreimvec(M[rl]+cl,re[rl]+cl,im[rl]+cl,ncol*nrow);

    // set the reference count
    D = new Reference;

    if (unbound(re)) re.Matrix::~Matrix();
    if (unbound(im)) im.Matrix::~Matrix();
}

//----------------------------------------------------------------------------//

ComplexMatrix::ComplexMatrix (const ComplexMatrix& A)
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

    // copy link to block
    M = A.M;
    D = A.D;
    
    // increase the reference count
    A.addref();

    // copy return value flag
    temporary =  A.temporary;
    form = A.form;
}

//----------------------------------------------------------------------------//

ComplexMatrix::~ComplexMatrix (void)
//
//  Destructor:
//  Removes a complex double precision matrix from the heap
//
{
    // decrease the reference count and delete if neccessary
    if (D) {
	if ( (--(D->count) == 0 && temporary == 0) || D->count < 0) {

	    // free data
	    delete[] (M[rl]+cl);
	    delete[] (M+rl);

	    // free info block
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

void ComplexMatrix::Remove (void)
//
// Explicitly removes a matrix from the memory
//
{
    if (D) {
	
	// reset reference count and call destructor
	D->count = 0;
	this->ComplexMatrix::~ComplexMatrix();

	// set attribute
	attribute = General;
	
	// reset the dimension information (upper must be less than lower)
	rl = cl = 1; nrow = ncol = rh = ch = 0;
    }
} 

//----------------------------------------------------------------------------//

void ComplexMatrix::Resize (int rlo, int rhi, int clo, int chi)
{
    (void)rlo; (void)rhi; (void)clo; (void)chi;  // Cast to void to skip -Wunused-parameter
    Matpack.Error("ComplexMatrix::Resize: NOT YET IMPLEMENTED");
}


//----------------------------------------------------------------------------//
// initializing a matrix from an array
//----------------------------------------------------------------------------//

ComplexMatrix& ComplexMatrix::operator << (const complex<double>* src)
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

ComplexMatrix& ComplexMatrix::operator = (complex<double> value)
//
// Set to the matrix value*Id, i.e. the diagonal is set to value,
// all other elements are set to zero.
//
{
    complex<double> **m = M + rl;
    int i;

    if (ncol != nrow&&abs(value)!=0) {
//    Matpack.Error("ComplexMatrix & ComplexMatrix::operator=(complex<double>): non square matrix\n");
      copyval( M[rl]+cl, value, ncol*nrow );  // For non-square matrix, set all element to value
      return *this;
    }

    copyval( M[rl]+cl, Zero, ncol*nrow );

    if (ncol == nrow)
    {for (i = cl; i <= ch; i++) (*m++)[i] = value;  // TODO
    }
    
    return *this;
}

//----------------------------------------------------------------------------//

ComplexMatrix& ComplexMatrix::operator += (complex<double> value)
//
// Add the matrix value*Id, i.e. add value to the diagonal elements.
//
{
    if (ncol != nrow)
      Matpack.Error("ComplexMatrix& ComplexMatrix::operator+=(complex<double>): non square matrix\n");

    complex<double> **m = M + rl;
    int i;
    for (i = cl; i <= ch; i++) (*m++)[i] += value;  // TODO
    return *this;
} 

//----------------------------------------------------------------------------//

ComplexMatrix& ComplexMatrix::operator -= (complex<double> value)
//
// Subtract the matrix value*Id, i.e. subtract value from the diagonal.
//
{
    if (ncol != nrow)
      Matpack.Error("ComplexMatrix& ComplexMatrix::operator-=(complex<double>): non square matrix\n");

    complex<double> **m = M + rl;
    int i;
    for (i = cl; i <= ch; i++) (*m++)[i] -= value;  // TODO
    return *this;
} 

//----------------------------------------------------------------------------//

ComplexMatrix& ComplexMatrix::operator *= (complex<double> value)
//
// Multiply with the matrix value*Id, i.e. multiply all elements 
// with value.
//
{
    mulval(M[rl]+cl,value,ncol*nrow);
    return *this;
} 

//----------------------------------------------------------------------------//

ComplexMatrix& ComplexMatrix::operator /= (complex<double> value)
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

ComplexMatrix& ComplexMatrix::operator = (const ComplexMatrix& A)
{
    if (unbound(A)) {        // right side is returned from a function

	if (this->Empty()) {  // left side is not yet initialized
	    cl = A.cl;
	    ch = A.ch;
	    rl = A.rl;
	    rh = A.rh;
	    nrow = A.nrow;
	    ncol = A.ncol;
	} else {             // free left side
	    checkdim(A);     // assign only compatible matrices 
	    this->ComplexMatrix::~ComplexMatrix();
	}

	// link right side to left side
	M = A.M; 
	D = A.D;
	D->count = 1;
	((ComplexMatrix&)A).temporary = 0;
	((ComplexMatrix&)A).D = 0;
			     
    } else {		  // right side is bound to a variable -> copy elements

	if (this->Empty()) {   // allocate
	    cl = A.cl;
	    ch = A.ch;
	    rl = A.rl;
	    rh = A.rh;
	    nrow = A.nrow;
	    ncol = A.ncol;
	    M = newmat(rl,rh,cl,ch);
	    D = new Reference;
	} else
	    checkdim(A);

	// copy right side to left side in one loop		
	copyvec( M[rl]+cl, A.M[rl]+cl, ncol*nrow );
		
        // set reference count	      
        D->count = 1; 
    }
    
    // copy attribute
    attribute = A.attribute;

    // set return value flag
    temporary = 0;
    return *this;
}

//----------------------------------------------------------------------------//

ComplexMatrix& ComplexMatrix::operator += (const ComplexMatrix& A)
{
    checkdim(A);

    // set attribute
    attribute = General;

    addvec( M[rl]+cl, A.M[rl]+cl, ncol*nrow );
    if (unbound(A)) ((ComplexMatrix&)A).ComplexMatrix::~ComplexMatrix();
    return *this;
} 

//----------------------------------------------------------------------------//

ComplexMatrix& ComplexMatrix::operator -= (const ComplexMatrix& A)
{
    checkdim(A);

    // set attribute
    attribute = General;

    subvec( M[rl]+cl, A.M[rl]+cl, ncol*nrow );
    if (unbound(A)) ((ComplexMatrix&)A).ComplexMatrix::~ComplexMatrix();
    return *this;
} 


//----------------------------------------------------------------------------//

ComplexMatrix& ComplexMatrix::operator %= (const ComplexMatrix& A)
{
    checkdim(A);

    // set attribute
    attribute = General;

    mulvec( M[rl]+cl, A.M[rl]+cl, ncol*nrow );
    if (unbound(A)) ((ComplexMatrix&)A).ComplexMatrix::~ComplexMatrix();
    return *this;
} 


//----------------------------------------------------------------------------//
// Submatrix extraction
//----------------------------------------------------------------------------//

ComplexMatrix ComplexMatrix::operator () (int rlo, int rhi, int clo, int chi) const
//
// The elements of this matrix within the index range [rlo..rhi,clo..chi] 
// are returned in a matrix with the corresponding dimension [rlo..rhi,clo..chi].
//

{
    // check for valid submatrix range
    if (rlo < rl || rhi > rh || clo < cl || chi > ch)
      Matpack.Error("ComplexMatrix::operator(): submatrix index out of range (%d,%d,%d,%d)", 
	   rlo,rhi,clo,chi);
    int i;
    ComplexMatrix W(rlo,rhi,clo,chi);    
    for (i = rlo; i <= rhi; i++)
      copyvec(W.M[i]+clo,M[i]+clo,W.ncol);
    return W.Value();
}


//----------------------------------------------------------------------------//
//  comparison operators 
//----------------------------------------------------------------------------//

int operator == (const ComplexMatrix& A, const ComplexMatrix& B)
//
// Returns 1, if all corresponding elements of the 
// given matrices are equal, and returns 0 otherwise.
//
{
    int retval;
    checkdim(A,B);
    retval = cmpvec( A.M[A.rl]+A.cl, B.M[B.rl]+B.cl, A.ncol*A.nrow );
    if (unbound(A)) ((ComplexMatrix&)A).ComplexMatrix::~ComplexMatrix();
    if (unbound(B)) ((ComplexMatrix&)B).ComplexMatrix::~ComplexMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (const ComplexMatrix& A, complex<double> value)
//
// Returns 1, if A is equal to value*Id, i.e. if the diagonal is set 
// to value and all other elements are 0, and returns 0 otherwise.
//
{
    checksquare(A);

    int rl = A.rl;
    int size = A.nrow;
    complex<double> **a = A.M + rl;
    complex<double> *c;
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
    if (unbound(A)) ((ComplexMatrix&)A).ComplexMatrix::~ComplexMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (complex<double> value, const ComplexMatrix& A)
//
// Returns 1, if A is equal to value*Id, i.e. if the diagonal is set 
// to value and all other elements are 0, and returns 0 otherwise.
//
{
    checksquare(A);

    int rl = A.rl;
    int size = A.nrow;
    complex<double> **a = A.M + rl;
    complex<double> *c;
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
    if (unbound(A)) ((ComplexMatrix&)A).ComplexMatrix::~ComplexMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int operator != (const ComplexMatrix& A, const ComplexMatrix& B)
//
// Returns 1, if A is not equal to B, 0 otherwise.
//
{
    int retval;
    checkdim(A,B);
    retval = ! cmpvec( A.M[A.rl]+A.cl, B.M[B.rl]+B.cl, A.ncol*A.nrow );
    if (unbound(A)) ((ComplexMatrix&)A).ComplexMatrix::~ComplexMatrix();
    if (unbound(B)) ((ComplexMatrix&)B).ComplexMatrix::~ComplexMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int operator != (const ComplexMatrix& A, complex<double> value)
//
// Returns 1, if A is not equal to value*Id, 0 otherwise.
//
{
    checksquare(A);

    int rl = A.rl;
    int size = A.nrow;
    complex<double> **a = A.M + rl;
    complex<double> *c;
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
    if (unbound(A)) ((ComplexMatrix&)A).ComplexMatrix::~ComplexMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int operator != (complex<double> value, const ComplexMatrix& A)
//
// Returns 1, if A is not equal to value*Id, 0 otherwise.
//
{
    checksquare(A);

    int rl = A.rl;
    int size = A.nrow;
    complex<double> **a = A.M + rl;
    complex<double> *c;
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
    if (unbound(A)) ((ComplexMatrix&)A).ComplexMatrix::~ComplexMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int ComplexMatrix::operator ! () const
//
// Negation operator, equivalent to (A == 0)
// Returns 1, if A is equal to 0, 1 otherwise.
//
{
    int retval;
    retval = cmpval(M[rl]+cl, Zero, ncol*nrow );
    if (unbound(*this)) ((ComplexMatrix*)this)->ComplexMatrix::~ComplexMatrix();
    return retval;  
}

//----------------------------------------------------------------------------//

ComplexMatrix::operator int () const
//
// Make if(A){...}  equivalent to  if(A != 0){...}.
//
{
    int retval;
    retval =  ! cmpval( M[rl]+cl, Zero, ncol*nrow );
    if (unbound(*this)) ((ComplexMatrix*)this)->ComplexMatrix::~ComplexMatrix();
    return retval;
}


//----------------------------------------------------------------------------//
// arithmetic operators
//----------------------------------------------------------------------------//

ComplexMatrix operator + (const ComplexMatrix& A, complex<double> value)
//
// Add matrix value*Id, i.e. add value to the diagonal.
//
{
    ComplexMatrix B; 
    B = A; // copies attribute
    B += value;
    return B.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator + (complex<double> value, const ComplexMatrix& A)
//
// Add matrix value*Id, i.e. add value to the diagonal.
//
{
    ComplexMatrix B;
    B = A; // copies attribute
    B += value;
    return B.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator + (const ComplexMatrix& A, const ComplexMatrix& B)
//
// add matrices
//
{
    checkdim(A,B);
    ComplexMatrix C(A.rl,A.rh,A.cl,A.ch);
    addvec(C.Store(),A.Store(),B.Store(),A.Elements());
    return C.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator % (const ComplexMatrix& A, const ComplexMatrix& B)
//
// array multiplication - element by element
//
{
    checkdim(A,B);
    ComplexMatrix C(A.rl,A.rh,A.cl,A.ch);
    mulvec(C.Store(),A.Store(),B.Store(),A.Elements());
    return C.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator - (const ComplexMatrix& A)
//
// Unary minus sign, i.e. multiply all matrix elements by -1.
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);  
    B.attribute = A.attribute;  // copy attribute
    negvec(B.Store(),A.Store(),A.Elements());
    return B.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator - (const ComplexMatrix& A, complex<double> value)
//
// Subtract matrix value*Id, i.e. subtract value from the diagonal.
//
{
    ComplexMatrix B;
    B = A; // copies attribute
    B -= value;
    return B.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator - (complex<double> value, const ComplexMatrix& A)
//
// Subtract this matrix from the matrix value*Id.
//
{
    if (A.ncol != A.nrow)
      Matpack.Error("ComplexMatrix operator-(complex<double>,ComplexMatrix&): non square matrix\n");

    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);
    B.attribute = A.attribute;                // copy attribute
    negvec(B.Store(),A.Store(),A.Elements()); // copy -A to B
    complex<double> **b = B.M + B.rl;
    for (int i = B.cl; i <= B.ch; i++) (*b++)[i] += value;
    return B.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator - (const ComplexMatrix& A, const ComplexMatrix& B)
//
// subtract matrices element by element
//
{
    checkdim(A,B);
    ComplexMatrix C(A.rl,A.rh,A.cl,A.ch);
    subvec(C.Store(),A.Store(),B.Store(),A.Elements());
    return C.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator * (const ComplexMatrix& A, complex<double> value)
//
// multiply all matrix elements by value
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);
    B.attribute = A.attribute;  // copy attribute
    mulval(B.Store(),A.Store(),value,A.Elements());
    return B.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator * (complex<double> value, const ComplexMatrix& A)
//
// multiply all matrix elements by value, same as above
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);
    B.attribute = A.attribute;  // copy attribute
    mulval(B.Store(),A.Store(),value,A.Elements());
    return B.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix operator / (const ComplexMatrix& A, complex<double> value)
//
// divide all matrix elements by value
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);
    B.attribute = A.attribute;  // copy attribute
    divval(B.Store(),A.Store(),value,A.Elements());
    return B.Value();
}	

//----------------------------------------------------------------------------//

ComplexMatrix operator / (complex<double> value, const ComplexMatrix& A)
//
// divide value by all matrix elements
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);
    B.attribute = A.attribute;  // copy attribute
    divval(B.Store(),value,A.Store(),A.Elements());
    return B.Value();
}	

//----------------------------------------------------------------------------//

ComplexVector operator * (const ComplexVector& V, const ComplexMatrix& A)
//
// multiply vector from left to matrix
//
{
    if (A.rl != V.cl || A.rh != V.ch) 
      Matpack.Error("ComplexVector operator * (ComplexVector&, ComplexMatrix&): "
	   "non conformant arguments\n");

    ComplexVector B(A.cl,A.ch);
    complex<double> **a = A.M+A.rl;
    complex<double> *b  = B.V+B.cl;
    complex<double> *v  = V.V+V.cl;
    complex<double> *vp,**ap,sum;
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

ComplexVector operator * (const ComplexMatrix& A, const ComplexVector& B)
//
// multiply a matrix A from the right with a vector B, returning C = A * B.
//
{
    int arl = A.rl, arh = A.rh, acl = A.cl, ncol = A.ncol;

    if (acl != B.cl  || A.ch != B.ch)
      Matpack.Error("Vector operator * (Matrix& A, Vector& B): "
	   "non conformant arguments");

    ComplexVector C(arl,arh);

    // avoid call to index operator that optimizes very badely
    const complex<double> * const *a = A.M;
    const complex<double> *b = B.Store();
    complex<double> *c = C.V;

    for (int i = arl; i <= arh; i++) 
      c[i] = matdotvec(a[i]+acl,b,ncol);

    return C.Value();
} 

//----------------------------------------------------------------------------//
// friend functions
//----------------------------------------------------------------------------//

complex<double> Trace (const ComplexMatrix& d)
//
// Returns the trace of the matrix, that means the sum of 
// the diagonal elements. If the matrix is not square an error
// will result.
//
{
    if (d.ncol != d.nrow)
      Matpack.Error("complex<double> Trace(const ComplexMatrix&): non square matrix\n");

    int j;
    complex<double> **m = d.M+d.rl;
    complex<double> trace = 0;
    for (j = d.cl; j <= d.ch; j++) trace += (*m++)[j];  // TODO
    return trace;
} 

//----------------------------------------------------------------------------//

double Norm (const ComplexMatrix& A)
//
// Returns the norm of the matrix A, which is defined by the 
// largest singular value of A.
//
{
    (void)A;  // Cast to void so gcc doesn't throw -Wunused-parameter
    double sum = 0;
    Matpack.Error("double Norm (const ComplexMatrix& A) NOT YET IMPLEMENTED");
    return sum;
}

//----------------------------------------------------------------------------//

double Norm2 (const ComplexMatrix& A)
//
// Same as Norm() above. 
// Returns the norm of the matrix A, which is defined by the 
// largest singular value of A.
//
{
    (void)A;  // Cast to void so gcc doesn't throw -Wunused-parameter
    double sum = 0;
    Matpack.Error("double Norm2 (const ComplexMatrix& A) NOT YET IMPLEMENTED");
    return sum;
}

//----------------------------------------------------------------------------//
// double NormFro (const ComplexMatrix& A)
// Returns the euclidean norm (Frobenius norm, F-norm) of the matrix, 
// that is the sum of modulus squared of all matrix elements.
//----------------------------------------------------------------------------//

double NormFro (const ComplexMatrix& A)
{
  int n = A.Elements();
  if (n <= 0) Matpack.Error("double NormFro (const ComplexMatrix& A) -- empty matrix");
  double sum = sqrt(norm2(A.Store(),n));
  return sum;
} 

//----------------------------------------------------------------------------//
// double Norm1 (const ComplexMatrix& d)
// Returns the 1-norm of the matrix, that is the maximum of sum
// of absolute values of a column.
//----------------------------------------------------------------------------//

double Norm1 (const ComplexMatrix& d)
{
  int j, rsize;
  complex<double> **m;
  double sum;
  double norm1 = 0;
  for (j = d.cl; j <= d.ch; j++) {
    rsize = d.nrow;
    m = d.M+d.rl;
    sum = 0;
    while (rsize--) sum += abs((*m++)[j]);
    if (sum > norm1) norm1 = sum;
  }
  return norm1;
} 

//----------------------------------------------------------------------------//
// double NormInf (const ComplexMatrix& d)
// Returns the infinity-norm of a matrix, that is the maximum of sum
// of absolute values of a row.
//----------------------------------------------------------------------------//

double NormInf (const ComplexMatrix& d)
{
  const complex<double> *v;
  complex<double> **m = d.M+d.rl;
  double norminf = 0;
  double sum;
  int j;
  int rsize = d.nrow;
  while (rsize--) { 
    v = *m++ + d.cl;
    j = d.ncol;
    sum = 0;
    while (j--) sum += abs(*v++);
    if (sum > norminf) norminf = sum;
  }
  return norminf;
} 

//----------------------------------------------------------------------------//
// complex<double> Sum (const ComplexMatrix& d)
// Returns the sum of all matrix elements.
//----------------------------------------------------------------------------//

complex<double> Sum (const ComplexMatrix& d)
{
  int n = d.ncol*d.nrow;  
  if (n <= 0) Matpack.Error("complex<double> Sum (const ComplexMatrix& d) -- empty matrix");
  complex<double> sum = Zero;
  const complex<double>* src = d.M[d.rl]+d.cl;
  while (n--) sum += *src++;  // TODO
  return sum;
}

//----------------------------------------------------------------------------//
// complex<double> Min (const ComplexMatrix& d)
// Returns the element of this matrix with the smallest norm.
//----------------------------------------------------------------------------//

complex<double> Min (const ComplexMatrix& d)
{
  int n = d.ncol*d.nrow-1;
  if (n < 0) Matpack.Error("complex<double> Min (const ComplexMatrix& d) -- empty matrix");
  const complex<double> *src = d.M[d.rl]+d.cl,
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
// complex<double> Min (const ComplexMatrix& d, int &i, int &k)
// Returns the element of this matrix with the smallest norm
// together with its index (i,k).
//----------------------------------------------------------------------------//

complex<double> Min (const ComplexMatrix& d, int &i, int &k)
{
  int n = d.ncol*d.nrow-1;
  if (n < 0) Matpack.Error("complex<double> Min (const ComplexMatrix& d, int &i, int &k) -- empty matrix");
  const complex<double> *src = d.M[d.rl]+d.cl,
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
  // calculate index
  n = int(minp-src);
  i = n / d.ncol + d.rl;
  k = n % d.ncol + d.cl;  
  return *minp;
} 

//----------------------------------------------------------------------------//
// complex<double> Max (const ComplexMatrix& d)
// Returns the element of this matrix with the largest norm.
//----------------------------------------------------------------------------//

complex<double> Max (const ComplexMatrix& d)
{
  int n = d.ncol*d.nrow-1;
  if (n < 0) Matpack.Error("complex<double> Max (const ComplexMatrix& d) -- empty matrix");
  const complex<double> *src = d.M[d.rl]+d.cl,
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
// complex<double> Max (const ComplexMatrix& d, int &i, int &k)
// Returns the element of this matrix with the largest norm
// together with its index (i,k).
//----------------------------------------------------------------------------//

complex<double> Max (const ComplexMatrix& d, int &i, int &k)
{
  int n = d.ncol*d.nrow-1;
  if (n < 0) Matpack.Error("complex<double> Max (const ComplexMatrix& d, int &i, int &k) -- empty matrix");
  const complex<double> *src = d.M[d.rl]+d.cl,
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
  // calculate index
  n = int(maxp-src);
  i = n / d.ncol + d.rl;
  k = n % d.ncol + d.cl;  
  return *maxp;
} 

//----------------------------------------------------------------------------//

complex<double> Det (const ComplexMatrix& M)
//
// Returns the determinant of this matrix. The LU decomposition is
// generated and the determinant is computed as the product of the
// diagonal values. If the matrix is not square or if the matrix 
// turns out to be singular an error will result.
//
{
    checksquare(M);

    int j;
    ComplexMatrix A;
    int sgn;
    complex<double> det;

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

Matrix Real (const ComplexMatrix& A)
//
// Returns the matrix of the real parts of the elements of A
//
{
    Matrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    double  *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = real(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Imag (const ComplexMatrix& A)
//
// Returns the matrix of the imaginary parts of the elements of A
//
{
    Matrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    double  *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = imag(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Abs (const ComplexMatrix& A)
//
// Returns the matrix of the moduli of the elements of A
//
{
    Matrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    double  *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = abs(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

Matrix AbsSqr (const ComplexMatrix& A)
//
// Returns the matrix of the modulus squared elements of A
//
{
    Matrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    double  *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = norm(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

Matrix Arg (const ComplexMatrix& A)
//
// Returns the matrix of the phases [-Pi..Pi] of the elements of A
// We don't use the arg() standard function from Complex.h because
// it is not save when called with a Zero argument. Use the MatPack
// inline function Angle() instead, which handles the Zero case:
// the argument of Complex zero is 0 by definition.
//
{
    Matrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    double  *b = B.M[B.rl]+B.cl;
    while (n--) { *b++ = Angle(real(*a),imag(*a)); a++; }  // TODO
    return B.Value();
}


//----------------------------------------------------------------------------//
// elementwise elementary functions
//----------------------------------------------------------------------------//

ComplexMatrix Cos (const ComplexMatrix& A)
//
// Returns the matrix of the cosines of the elements of A
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    complex<double> *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = cos(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix Sin (const ComplexMatrix& A)
//
// Returns the matrix of the sines of the elements of A
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    complex<double> *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = sin(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix Cosh (const ComplexMatrix& A)
//
// Returns the matrix of the hyperbolic cosines of the elements of A
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    complex<double> *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = cosh(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix Sinh (const ComplexMatrix& A)
//
// Returns the matrix of the hyperbolic sines of the elements of A
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    complex<double> *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = sinh(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix Exp (const ComplexMatrix& A)
//
// Returns the matrix of the exponentials of the elements of A
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    complex<double> *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = exp(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix Log (const ComplexMatrix& A)
//
// Returns the matrix of the natural logarithms of the elements of A
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    complex<double> *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = log(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix Sqrt (const ComplexMatrix& A)
//
// Returns the matrix of the square roots of the elements of A
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    complex<double> *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = sqrt(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix Sqr (const ComplexMatrix& A)
//
// Returns the matrix of the squares of the elements of A
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    complex<double> *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = sqr(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix Cube (const ComplexMatrix& A)
//
// Returns the matrix of the cubes of the elements of A
//
{
    ComplexMatrix B(A.rl,A.rh,A.cl,A.ch);

    // copy attribute
    B.attribute = A.attribute;

    int n = A.nrow * A.ncol;
    complex<double> *a = A.M[A.rl]+A.cl;
    complex<double> *b = B.M[B.rl]+B.cl;
    while (n--) *b++ = cube(*a++);  // TODO
    return B.Value();
}

//----------------------------------------------------------------------------//
// other functions
//----------------------------------------------------------------------------//


int MatchingIndexRange (const ComplexMatrix &a, const ComplexMatrix &b)
//
// Return True if matrices have the same index range, false otherwise.
//
{
    return (a.cl == b.cl && a.ch == b.ch && a.rl == b.rl && a.rh == b.rh);
}

//----------------------------------------------------------------------------//

ComplexVector Pack (const ComplexMatrix &A)
//
// pack the matrix elements row by row into a vector
//
{
    int n = A.ncol*A.nrow;
    ComplexVector V(A.cl,A.cl+n-1);
    copyvec(V.V+V.cl, A.M[A.rl]+A.cl,n);
    return V.Value();
}

//----------------------------------------------------------------------------//
// member functions 
//----------------------------------------------------------------------------//

void ComplexMatrix::Set (complex<double> value)
//
// Sets all matrix elements to the given value.
// The vector will be changed inplace.
//
{
    copyval( M[rl]+cl, value, ncol*nrow );
}

//----------------------------------------------------------------------------//

void ComplexMatrix::Set (complex_mapper fcn)
//
// Sets all matrix elements by applying the given function.
// The vector will be changed inplace.
//
{
    int n = ncol*nrow;
    complex<double>* src = M[rl]+cl;
    while (n--) {*src = (*fcn)(*src); src++;}  // TODO
}

//----------------------------------------------------------------------------//

void ComplexMatrix::Set (complex_matrix_index_mapper fcn)
//
// Sets all matrix elements by applying the given function.
// The function takes the index pair and the value as arguments.
// The vector will be changed inplace.
//
{
    int r,c;
    complex<double>* src = M[rl]+cl;
 
    for (r = rl; r <= rh; r++) 
      for (c = cl; c <= ch; c++) {
	  *src = (*fcn)(r,c,*src); src++;  // TODO
      }
}

//----------------------------------------------------------------------------//

ComplexMatrix& ComplexMatrix::ShiftIndex (int row, int col)
//
// Shift the index range of the matrix by the given offsets
// in columns and rows. The matrix will be changed inplace.
//
// Dimension integers are not checked for overflow !
//
// last change 22.11.1996
{
    if (Empty()) 
      Matpack.Error("ComplexMatrix::ShiftIndex: matrix is not allocated");

    // shift column pointers
    if (col) {
        complex<double> **v = M+rl; 
	int n = nrow;
	while (n--) { *v -= col; v++; }
	cl += col; ch += col;  // shift dimension information
    }

    // shift row pointer
    if (row) {
        M -= row;
        rl += row; rh += row;  // shift dimension information
    }
    return *this;
}

//----------------------------------------------------------------------------//

ComplexVector ComplexMatrix::Row (int i) const
//
// Returns the i-th row of this matrix as a vector.
//
{
    if (i < rl || i > rh) 
      Matpack.Error("ComplexMatrix::Row(int): index out of range (%d)\n",i); 
    
    ComplexVector A(cl,ch);
    copyvec(A.V+cl,M[i]+cl,ncol);
    return A.Value();
}

//----------------------------------------------------------------------------//

ComplexVector ComplexMatrix::Column (int i) const
//
// Returns the i-th column of this matrix as a vector.
//
{
    if (i < cl || i > ch)
      Matpack.Error("ComplexMatrix::Column(int): index out of range (%d)\n",i); 
    
    ComplexVector A(rl,rh);
    complex<double> **m = M + rl;
    complex<double>  *v = A.V + rl;
    int rsize = nrow;
    while (rsize--) *v++ = (*m++)[i];  // TODO
    return A.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix ComplexMatrix::Transpose (void) const
//
// Returns the transpose of this matrix.
//
{
    ComplexMatrix A(cl,ch,rl,rh);

    if ( (attribute == UpperTriangular) && (ncol == nrow) )
      A.attribute = LowerTriangular;
    else if ( (attribute == LowerTriangular) && (ncol == nrow) )
      A.attribute = UpperTriangular;

    complex<double> **m = M + rl;
    complex<double> **a = A.M + cl;
    complex<double> **ap;
    complex<double> *v;
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

ComplexMatrix ComplexMatrix::Hermitean (void) const
//
// Returns the hermitean conjugate of this matrix,
// that means transposed and conjugate complex.
//
{
    ComplexMatrix A(cl,ch,rl,rh);

    if ( (attribute == UpperTriangular) && (ncol == nrow) )
      A.attribute = LowerTriangular;
    else if ( (attribute == LowerTriangular) && (ncol == nrow) )
      A.attribute = UpperTriangular;

    complex<double> **m = M + rl;
    complex<double> **a = A.M + cl;
    complex<double> **ap;
    complex<double> *v;
    int j;
    int i;
    for (i = rl; i <= rh; i++) { 
	v = *m++ + cl;
	ap = a;
	j = ncol;
	while (j--) (*ap++)[i] = conj(*v++);
    }
    return A.Value();
} 

//----------------------------------------------------------------------------//

ComplexMatrix ComplexMatrix::Conjugate (void) const
//
// Returns the complex conjugate of this matrix
//
{
    ComplexMatrix A(rl,rh,cl,ch);

    // copy attribute
    A.attribute = attribute;

    int n = ncol*nrow;
    complex<double>* src = M[rl]+cl;
    complex<double>* dst = A.M[rl]+cl;
    while (n--) *dst++ = conj(*src++);  // TODO
    return A.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix ComplexMatrix::Apply (complex_mapper fcn) const
//
// Applies the given function to all elements
//
{
    ComplexMatrix A(rl,rh,cl,ch);

    // copy attribute
    A.attribute = attribute;

    int n = ncol*nrow;
    complex<double>* src = M[rl]+cl;
    complex<double>* dst = A.M[rl]+cl;
    while (n--) *dst++ = (*fcn)(*src++);  // TODO
    return A.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix ComplexMatrix::Apply (complex_matrix_index_mapper fcn) const
//
// Applies the given function to all elements.
// The function takes the index pair and the value as arguments.
//
{
    ComplexMatrix A(rl,rh,cl,ch);

    // copy attribute
    A.attribute = attribute;

    int r,c;
    complex<double>* src = M[rl]+cl;
    complex<double>* dst = A.M[rl]+cl;

    for (r = rl; r <= rh; r++) 
      for (c = cl; c <= ch; c++) 
	*dst++ = (*fcn)(r,c,*src++);  // TODO
    
    return A.Value();
}

//----------------------------------------------------------------------------//

ComplexMatrix ComplexMatrix::Inverse (void) const
//
// Returns the inverse of this matrix. If the matrix is not square 
// or if the matrix turns out to be singular an error will result.
//
{
    ComplexMatrix Y(rl,rh,rl,rh), LU; 
    Y = 1;              // make unit matrix (diagonal 1)
    LU = *this;         // copy of this matrix - protect from overwriting
    SolveLinear(LU,Y);  // create inverse in Y
    return Y.Value();
}

//----------------------------------------------------------------------------//

int ComplexMatrix::Attribute (int newAttr)
// 
// Set new attribute checking for validity. 
// Returns True if successful, false if the attribute can't be set:
// For example you can't set the LowerTrianguar attribute if the matrix 
// is not square !
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
	Matpack.Error("ComplexMatrix::Attribute(int): unknown attribute (%d)\n",newAttr);
    }
    return status;
}

//----------------------------------------------------------------------------//

ComplexMatrix FlipUpDown (const ComplexMatrix& A)
//
// Returns the matrix with rows flipped, 
// i.e. up and down is mirrored.
//
{    
    int arl = A.rl, arh = A.rh,
        acl = A.cl, ach = A.ch,
	n = arl + arh,
        k = Nint(0.5*n);

    ComplexMatrix B(arl,arh,acl,ach); 
    // attributes changed: result is general matrix

    // avoid call to index operator that optimizes very badely
    complex<double> **a = A.M, **b = B.M;

    for (int i = arl; i <= k; i++) 
	for (int j = acl; j <= ach; j++) {
	    b[i][j] = a[n-i][j];
	    b[n-i][j] = a[i][j];
	}

    return B.Value();
}


//----------------------------------------------------------------------------//

ComplexMatrix FlipLeftRight (const ComplexMatrix& A)
//
// Returns the matrix with columns flipped, 
// i.e. left and right is mirrored.
//
{    
    int arl = A.rl, arh = A.rh,
        acl = A.cl, ach = A.ch,
	n = acl + ach,
        k = Nint(0.5*n);

    ComplexMatrix B(arl,arh,acl,ach); 
    // attributes changed: result is general matrix

    // avoid call to index operator that optimizes very badely
    complex<double> **a = A.M, **b = B.M;

    for (int i = arl; i <= arh; i++) 
	for (int j = acl; j <= k; j++) {
	    b[i][j] = a[i][n-j];
	    b[i][n-j] = a[i][j];
	}

    return B.Value();
}

//----------------------------------------------------------------------------//

