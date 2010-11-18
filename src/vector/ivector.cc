/*----------------------------------------------------------------------------*\
| implementation of the integer                                     ivector.cc |
| vector class of MatPack.                                                     |
|                                                                              |
| MatPack Libary Release 1.0                                                   |
| Copyright (C) 1990-1994 by Berndt M. Gammel                                  |
|                                                                              |
| Permission to use, copy, modify, and distribute this software and its        |
| documentation for any non-commercial purpose and without fee is hereby       |
| granted, provided that the above copyright notice appear in all copies       |
| and that both that copyright notice and this permission notice appear        |
| in supporting documentation. This software is provided "as is" without       |
| express or implied warranty.                                                 |
|                                                                              |
\*----------------------------------------------------------------------------*/

#include "vector.h"
#include "vecinl.h"
#include <cstdlib>

//----------------------------------------------------------------------------//
// error messages for vector and matrix functions
//----------------------------------------------------------------------------//

static const char *NonConformVector = "non conformant integer vector";
static const char *AllocationFail   = "integer vector allocation failure";
static int Zero = 0;


//----------------------------------------------------------------------------//
// Notes:
//----------------------------------------------------------------------------//
// All inline functions defined in the following really  m u s t  be
// inlined (choose the appropriate compiler options) - otherwise
// the matrix package will be slow. Do it by hand if your compiler
// is incapable of doing it !!!
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
// integer vector and matrix classes auxilliary inlines
//----------------------------------------------------------------------------//

inline void IntVector::checkdim (const IntVector& A)
//
// Checks if the dimensions of the given vector matches this vector
//
{
    if (A.cl != cl || A.ch != ch) Matpack.Error(NonConformVector);
}

void checkdim (const IntVector& A, const IntVector& B)
//
// Checks if the dimensions of the given vectors match
//
{
    if (A.cl != B.cl || A.ch != B.ch) Matpack.Error(NonConformVector);
}


//----------------------------------------------------------------------------//
// integer vector class constructors and destructors
//----------------------------------------------------------------------------//

IntVector::IntVector (void)
// 
//  Constructor:
//  Define an integer vector which is not yet allocated
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

IntVector::IntVector (int ncl, int nch)
// 
//  Constructor:
//  Allocates an integer vector[ncl..nch]
//
{
    int* v;

    // set the dimension information
    cl = ncl; ch = nch; ncol = nch-ncl+1;

    // allocate the vector structure
    if ( (v = new int[ncol]) ) 
      V = v - ncl;
    else
      Matpack.Error(AllocationFail);
 
    // set the reference count
    D = new Reference;
    temporary = 0;
    form = MpTextFormat;
}                   
                  
//----------------------------------------------------------------------------//

IntVector::IntVector (int ncl, int nch, int value)
// 
//  Constructor:
//  Allocates an integer vector[ncl..nch] and initializes
//  all elements with the given value
//
{
    register int* v;

    // set the dimension information
    cl = ncl; ch = nch; ncol = nch-ncl+1;

    // allocate the vector structure
    if ((v = new int[ncol])) {
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

IntVector::IntVector (const IntVector &A)
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

IntVector::~IntVector (void) 
//
// Destructor:
// Removes an integer vector from the heap
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

void IntVector::Remove (void)
//
// Explicitly removes an integer vector from the memory
//
{
    if (D) { 

	// reset reference count and call destructor
	D->count = 0;
	this->IntVector::~IntVector();

	// set the dimension information (upper must be less than lower)
	cl = 1;	ncol = ch = 0;
    } 
}

//----------------------------------------------------------------------------//

void IntVector::Resize (int ncl, int nch)
//
// Resizes a vector copying the old contents padding the new elements 
// with zero
//
{
    int *v;
    int csize = nch-ncl+1;
    
    // check for identical resize
    if (ncl == cl && nch == ch) return;

    // allocate a new vector 
    if ((v = new int[csize])) {
	
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
    cl = ncl; ch = nch; ncol = nch-ncl+1;
}

//----------------------------------------------------------------------------//
// initializing a vector from an array
//----------------------------------------------------------------------------//

IntVector& IntVector::operator << (const int* src)
//
// Initializes the vector from the values given in a static array
// The number of elements given in the array must match the vector size !
//
{
    copyvec(V+cl,src,ncol);
    return *this;
}

//----------------------------------------------------------------------------//
// vector value assignment operators
//----------------------------------------------------------------------------//

IntVector& IntVector::operator = (const IntVector& A)
{
    int *v;
    
    if (unbound(A)) {        // right side is returned from a function

	if (this->Empty()) {  // left side is not yet initialized
	    cl = A.cl;
	    ch = A.ch;
	    ncol = A.ncol;

	} else {             // free left side
	    checkdim(A);     // assign only compatible vectors
	    this->IntVector::~IntVector();
	}

	// link right side to left side
	V = A.V; 
	D = A.D;
	D->count = 1;
	((IntVector&)A).temporary = 0;
	((IntVector&)A).D = 0;
			     
    } else {		 // right side is bound to a variable -> copy elements

	if (this->Empty()) {   // allocate
	    cl = A.cl;
	    ch = A.ch;
	    ncol = ch-cl+1;
			      
	    // allocate the vector structure
	    if ((v = new int[ncol])) {
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
// scalar value assignment operators
//----------------------------------------------------------------------------//

IntVector& IntVector::operator = (int value)
{
    copyval(V+cl,value,ncol);
    return *this;
}

//----------------------------------------------------------------------------//
// Subvector extraction
//----------------------------------------------------------------------------//

IntVector IntVector::operator () (int lo, int hi) const
//
// The elements of this vector within the index range [lo..hi] 
// are returned in a vector with the corresponding dimension [lo..hi].
//
{
    // check for valid subvector range
    if (lo < cl || hi > ch)
      Matpack.Error("IntVector::operator(): subvector index out of range (%d,%d)", lo,hi);
    IntVector W(lo,hi);    
    copyvec(W.V+lo,V+lo,W.ncol);
    return W.Value();
}


//----------------------------------------------------------------------------//
// comparison operators 
//----------------------------------------------------------------------------//

int operator == (const IntVector& A, const IntVector& B)
{
    int retval;
    checkdim(A,B);
    retval = cmpvec(A.V+A.cl,B.V+B.cl,A.ncol);
    if (unbound(A)) ((IntVector&)A).IntVector::~IntVector();
    if (unbound(B)) ((IntVector&)A).IntVector::~IntVector();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (const IntVector& A, int value)
{
    int retval;
    retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((IntVector&)A).IntVector::~IntVector();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (int value, const IntVector& A)
{
    int retval;
    retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((IntVector&)A).IntVector::~IntVector();
    return retval;
}

//----------------------------------------------------------------------------//

int operator != (const IntVector& A, const IntVector& B)
{
    int retval;
    checkdim(A,B);
    retval = cmpvec(A.V+A.cl,B.V+B.cl,A.ncol);
    if (unbound(A)) ((IntVector&)A).IntVector::~IntVector();
    if (unbound(B)) ((IntVector&)A).IntVector::~IntVector();
    return (! retval);
}

//----------------------------------------------------------------------------//

int operator != (const IntVector& A, int value)
{
    int retval;
    retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((IntVector&)A).IntVector::~IntVector();
    return (! retval);
}

//----------------------------------------------------------------------------//

int operator != (int value, const IntVector& A)
{
    int retval;
    retval = cmpval(A.V+A.cl,value,A.ncol);
    if (unbound(A)) ((IntVector&)A).IntVector::~IntVector();
    return (! retval);
}

//----------------------------------------------------------------------------//

int IntVector::operator ! () const
//
// Negation operator, equivalent to (A == 0)
// Returns 1, if A is equal to 0, 1 otherwise.
//
{
    int retval;
    retval = cmpval(V+cl,Zero,ncol);
    if (unbound(*this)) ((IntVector*)this)->IntVector::~IntVector();
    return retval;
}

//----------------------------------------------------------------------------//

IntVector::operator int () const
//
// Make if(A){...}  equivalent to  if(A != 0){...}.
//
{
    int retval;
    retval =  ! cmpval(V+cl,Zero,ncol);
    if (unbound(*this)) ((IntVector*)this)->IntVector::~IntVector();
    return retval;
}


//----------------------------------------------------------------------------//
// friend functions 
//----------------------------------------------------------------------------//

int Min (const IntVector& d)
//
// Returns the smallest element
//
{
    register int n = d.ch-d.cl;
    register int* src = d.V+d.cl;
    int m = *src++;
    while (n--) {
	if (*src < m) m = *src; 
	src++; 
    }
    return m;
}

//----------------------------------------------------------------------------//

int Min (const IntVector& d, int& i)
//
// Returns the smallest element and its index i
//
{
    Matpack.Error("int Min (const IntVector& d, int& i) NOT YET IMPLEMENTED");
    return 0;
}

//----------------------------------------------------------------------------//

int Max (const IntVector& d)
//
// Returns the largest element
//
{
    register int n = d.ch-d.cl;
    register int* src = d.V+d.cl;
    int m = *src++;
    while (n--) {
	if (*src > m) m = *src; 
	src++; 
    }
    return m;
}

//----------------------------------------------------------------------------//

int Max (const IntVector& d, int& i)
//
// Returns the largest element and its index i
//
{
    Matpack.Error("int Max (const IntVector& d, int& i) NOT YET IMPLEMENTED");
    return 0;
}

//----------------------------------------------------------------------------//

static int comp_int (const void* d1, const void* d2)
{
    return *(int*)d1 - *(int*)d2;
}

//----------------------------------------------------------------------------//

void Sort (IntVector& d)
//
// Sort integer vector inplace into ascending numerical order using 
// the system qsort() function
//
{
    unsigned size = d.ncol;
    int *data = d.V + d.cl;
    qsort(data,size,sizeof(int),comp_int);
}

//----------------------------------------------------------------------------//

void Reverse (IntVector& d)
//
// Reverse the vector inplace
//
{
    unsigned size = d.ncol / 2;
    register int *d1 = d.V+d.cl;
    register int *d2 = d.V+d.ch;
    register int t;
    while (size--) {
	t = *d1;
	*d1++ = *d2;
	*d2-- = t;
    }
}

//----------------------------------------------------------------------------//

int MatchingIndexRange (const IntVector &a, const IntVector &b)
//
// Return True if vectors have the same index range, False otherwise.
//
{
    return (a.cl == b.cl && a.ch == b.ch);
}


//----------------------------------------------------------------------------//
// member functions
//----------------------------------------------------------------------------//

void IntVector::Set (int value)
//
// Sets all vector elements to the given value.
// The vector will be changed inplace.
//
{
    copyval(V+cl,value,ncol);    
}

//----------------------------------------------------------------------------//

void IntVector::Set (int_mapper fcn)
//
// Sets all vector elements by applying the given function.
// The vector will be changed inplace.
//
{
    int n = ncol;
    int* src = V+cl;
    for (int i = 0; i < n; i++) src[i] = fcn(src[i]);
}

//----------------------------------------------------------------------------//

void IntVector::Set (int_vector_index_mapper fcn)
//
// Sets all vector elements by applying the given function.
// The function takes the index and the value as arguments.
// The vector will be changed inplace.
//
{
    int* src = V;
    for (int i = cl; i <= ch; i++) src[i] = fcn(i,src[i]);
}

//----------------------------------------------------------------------------//

IntVector& IntVector::ShiftIndex (int by)
//
// Shift the index range of the vector by the given offset
// The vector will be changed inplace.
//
{
    if (this->Empty()) 
      Matpack.Error("IntVector::ShiftIndex: vector is not allocated");
    
    // shift pointers and dimension information
    V -= by; cl += by; ch += by;

    return *this;
}

//----------------------------------------------------------------------------//
