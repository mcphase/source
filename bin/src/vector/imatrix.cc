/*----------------------------------------------------------------------------*\
| implementation of the integer matrix class of Matpack             imatrix.cc |
|                                                                              |
| MatPack Libary Release 1.0                                                   |
| Copyright (C) 1990-1995 by Berndt M. Gammel                                  |
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

static const char *NonSquareMatrix  = "non square matrix";
static const char *NonConformMatrix = "non conformant matrix or vector";
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

inline void IntMatrix::checkdim (const IntMatrix& A)
//
// Checks if the dimensions of the matries match when adding
//
{
    if (A.cl != cl || A.ch != ch || A.rl != rl || A.rh != rh) 
      Matpack.Error(Mat::UnspecifiedError,NonConformMatrix);
}

void checkdim (const IntMatrix& A, const IntMatrix& B)
//
// Checks if the dimensions of the matrices match when adding
//
{
    if (A.cl != B.cl || A.ch != B.ch || A.rl != B.rl || A.rh != B.rh) 
      Matpack.Error(Mat::UnspecifiedError,NonConformMatrix);
}


void checksquare (const IntMatrix &M)
{
    if (M.rl != M.cl || M.rh != M.ch) Matpack.Error(Mat::UnspecifiedError,NonSquareMatrix);
}

//----------------------------------------------------------------------------//

static int** newmat (int nrl, int nrh, int ncl, int nch)
//
// allocate the dynamic part of a matrix on the heap, 
// allocate matrix in one block to faciliate many matrix functions !
//
{
    int* v;
    int **M;
    register int** m;
    register int rsize = nrh-nrl+1; 
    register int csize = nch-ncl+1;

    if ((m = new int*[rsize])) {
	M = m - nrl;
	if ((v = new int[csize*rsize])) {
	    for (v -= ncl; rsize--; v += csize)  *m++ = v;
	    return M;
	}
    }
    Matpack.Error("IntMatrix(%d,%d,%d,%d): newmat: allocation failure",nrl,nrh,ncl,nch);
    return (int**) 0;
}


//----------------------------------------------------------------------------//
// matrix class constructors and destructors
//----------------------------------------------------------------------------//

IntMatrix::IntMatrix (void)
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

IntMatrix::IntMatrix (int nrl, int nrh, int ncl, int nch)
// 
//  Constructor:
//  Allocates an integer matrix[nrl..nrh,ncl..nch] 
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

IntMatrix::IntMatrix (int nrl, int nrh, int ncl, int nch, int value)
// 
//  Constructor:
//  Allocates an integer matrix[nrl..nrh,ncl..nch] and initializes 
//  all elements with the given value
//
{
    // set the dimension information
    rl = nrl; rh = nrh; nrow = nrh-nrl+1;
    cl = ncl; ch = nch; ncol = nch-ncl+1;

    // allocate the matrix structure
    M = newmat(nrl,nrh,ncl,nch);

    // set all elements
    copyval(M[nrl]+ncl,value,ncol*nrow);

    // set the reference count
    D = new Reference;
}     

//----------------------------------------------------------------------------//

IntMatrix::IntMatrix (const IntMatrix& A)
// 
//  Copy constructor:
//  adds a reference to the block, is called automatically when
//  returning a matrix from a function !
//
{
    // copy the dimension information
    rl = A.rl; rh = A.rh; 
    cl = A.cl; ch = A.ch;
    nrow = A.nrow; 
    ncol = A.ncol;

    // copy link 
    M = A.M;
    D = A.D;
    
    // increase the reference count
    A.addref();

    // copy return value flag
    temporary =  A.temporary;
    form = A.form;
}

//----------------------------------------------------------------------------//

IntMatrix::~IntMatrix (void)
//
//  Destructor:
//  Removes an integer matrix from the heap
//
{
    // decrease the reference count and delete if neccessary
    if (D) {
	if ( (--(D->count) == 0 && temporary == 0) || D->count < 0) {

	    // free data
	    delete (M[rl]+cl);
	    delete (M+rl);

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

void IntMatrix::Remove (void)
//
// Explicitly removes a matrix from the memory
//
{
    if (D) {
	
	// reset reference count and call destructor
	D->count = 0;
	this->IntMatrix::~IntMatrix();

	// reset the dimension information (upper must be less than lower)
	rl = cl = 1; nrow = ncol = rh = ch = 0;
    }
} 

//----------------------------------------------------------------------------//

void IntMatrix::Resize (int nrl, int nrh, int ncl, int nch)
{

    Matpack.Error("IntMatrix::Resize: NOT YET IMPLEMENTED");

    int **m;
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

IntMatrix& IntMatrix::operator << (const int* src)
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

IntMatrix& IntMatrix::operator = (int value)
//
// Set to the matrix value*Id, i.e. the diagonal is set to value,
// all other elements are set to zero.
//
{
    if (ncol != nrow)
      Matpack.Error("Matrix& Matrix::operator=(double): non square matrix\n");

    register int **m = M + rl;
    register int i;
    copyval( M[rl]+cl, Zero, ncol*nrow );
    for (i = cl; i <= ch; i++) (*m++)[i] = value;
    return *this;
}


//----------------------------------------------------------------------------//
// matrix value assignment operators
//----------------------------------------------------------------------------//

IntMatrix& IntMatrix::operator = (const IntMatrix& A)
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
	    this->IntMatrix::~IntMatrix();
	}

	// link right side to left side
	M = A.M; 
	D = A.D;
        D->count = 1;
	((IntMatrix&)A).temporary = 0;
	((IntMatrix&)A).D = 0;
     
    } else {		 // right side is bound to a variable -> copy elements

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
    
    // set return value flag
    temporary = 0;  
    return *this;
}


//----------------------------------------------------------------------------//
// Submatrix extraction 
//----------------------------------------------------------------------------//

IntMatrix IntMatrix::operator () (int rlo, int rhi, int clo, int chi) const
//
// The elements of this matrix within the index range [rlo..rhi,clo..chi] 
// are returned in a matrix with the corresponding dimension [rlo..rhi,clo..chi].
//
{
    // check for valid submatrix range
    if (rlo < rl || rhi > rh || clo < cl || chi > ch)
      Matpack.Error("IntMatrix::operator(): submatrix index out of range (%d,%d,%d,%d)", 
	   rlo,rhi,clo,chi);
    int i;
    IntMatrix W(rlo,rhi,clo,chi);    
    for (i = rlo; i <= rhi; i++)
      copyvec(W.M[i]+clo,M[i]+clo,W.ncol);
    return W.Value();
}

//----------------------------------------------------------------------------//
// comparison operators
//----------------------------------------------------------------------------//

int operator == (const IntMatrix& A, const IntMatrix& B)
//
// Returns 1, if all corresponding elements of the 
// given matrices are equal, and returns 0 otherwise.
//
{
    int retval;
    checkdim(A,B);
    retval = cmpvec( A.M[A.rl]+A.cl, B.M[B.rl]+B.cl, A.ncol*A.nrow );
    if (unbound(A)) ((IntMatrix&)A).IntMatrix::~IntMatrix();
    if (unbound(B)) ((IntMatrix&)A).IntMatrix::~IntMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (const IntMatrix& A, int value)
//
// Returns 1, if A is equal to value*Id, i.e. if the diagonal is set 
// to value and all other elements are 0, and returns 0 otherwise.
//
{
    checksquare(A);

    int rl = A.rl;
    int size = A.nrow;
    register int **a = A.M + rl;
    register int *c;
    register int i;
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
    if (unbound(A)) ((IntMatrix&)A).IntMatrix::~IntMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int operator == (int value, const IntMatrix& A)
//
// Returns 1, if A is equal to value*Id, i.e. if the diagonal is set 
// to value and all other elements are 0, and returns 0 otherwise.
//
{
    checksquare(A);

    int rl = A.rl;
    int size = A.nrow;
    register int **a = A.M + rl;
    register int *c;
    register int i;
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
    if (unbound(A)) ((IntMatrix&)A).IntMatrix::~IntMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int operator != (const IntMatrix& A, const IntMatrix& B)
//
// Returns 1, if A is not equal to B, 0 otherwise.
//
{
    int retval;
    checkdim(A,B);
    retval = ! cmpvec( A.M[A.rl]+A.cl, B.M[B.rl]+B.cl, A.ncol*A.nrow );
    if (unbound(A)) ((IntMatrix&)A).IntMatrix::~IntMatrix();
    if (unbound(B)) ((IntMatrix&)B).IntMatrix::~IntMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int operator != (const IntMatrix& A, int value)
//
// Returns 1, if A is not equal to value*Id, 0 otherwise.
//
{
    checksquare(A);

    int rl = A.rl;
    int size = A.nrow;
    register int **a = A.M + rl;
    register int *c;
    register int i;
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
    if (unbound(A)) ((IntMatrix&)A).IntMatrix::~IntMatrix();
    return retval;
}

int operator != (int value, const IntMatrix& A)
//
// Returns 1, if A is not equal to value*Id, 0 otherwise.
//
{
    checksquare(A);

    int rl = A.rl;
    int size = A.nrow;
    register int **a = A.M + rl;
    register int *c;
    register int i;
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
    if (unbound(A)) ((IntMatrix&)A).IntMatrix::~IntMatrix();
    return retval;
}

//----------------------------------------------------------------------------//

int IntMatrix::operator ! () const
//
// Negation operator, equivalent to (A == 0)
// Returns 1, if A is equal to 0, 1 otherwise.
//
{
    int retval;
    retval = cmpval(M[rl]+cl, Zero, ncol*nrow );
    if (unbound(*this)) ((IntMatrix*)this)->IntMatrix::~IntMatrix();
    return retval;  
}

//----------------------------------------------------------------------------//

IntMatrix::operator int () const
//
// Make if(A){...}  equivalent to  if(A != 0){...}.
//
{
    int retval;
    retval =  ! cmpval( M[rl]+cl, Zero, ncol*nrow );
    if (unbound(*this)) ((IntMatrix*)this)->IntMatrix::~IntMatrix();
    return retval;
}


//----------------------------------------------------------------------------//
// friend functions
//----------------------------------------------------------------------------//

int Min (const IntMatrix &d)
//
// Returns the smallest element of this matrix.
//
{
    register int n = d.ncol*d.nrow-1;
    register int* src = d.M[d.rl]+d.cl;
    register int* minp = src++;
    while (n--) {
	if (*src < *minp) minp = src; 
	src++; 
    }
    return *minp;
} 

//----------------------------------------------------------------------------//

int Min (const IntMatrix &d, int& i, int& j)
//
// Returns the smallest element of this matrix and its indices i,j
//
{
    Matpack.Error("int Min (const IntMatrix &d, int& i, int& j) NOT YET IMPLEMENTED");
    return 0;
}

//----------------------------------------------------------------------------//

int Max (const IntMatrix &d)
//
// Returns the largest element of this matrix.
//
{
    register int n = d.ncol*d.nrow-1;
    register int* src = d.M[d.rl]+d.cl;
    register int* maxp = src++;
    while (n--) {
	if (*src > *maxp) maxp = src; 
	src++; 
    }
    return *maxp;
} 

//----------------------------------------------------------------------------//

int Max (const IntMatrix &d, int& i, int& j)
//
// Returns the largest element of this matrix and its indices i,j
//
{
    Matpack.Error("int Max (const IntMatrix &d, int& i, int& j) NOT YET IMPLEMENTED");
    return 0;
}

//----------------------------------------------------------------------------//

int MatchingIndexRange (const IntMatrix &a, const IntMatrix &b)
//
// Return True if matrices have the same index range, False otherwise.
//
{
    return (a.cl == b.cl && a.ch == b.ch && a.rl == b.rl && a.rh == b.rh);
}

//----------------------------------------------------------------------------//

IntVector Pack (const IntMatrix &A)
//
// pack the matrix elements row by row into a vector
//
{
    int n = A.ncol*A.nrow;
    IntVector V(A.cl,A.cl+n-1);
    copyvec(V.V+V.cl, A.M[A.rl]+A.cl,n);
    return V.Value();
}

//----------------------------------------------------------------------------//
// member functions
//----------------------------------------------------------------------------//

void IntMatrix::Set (int value)
//
// Sets all matrix elements to the given value
//
{
    copyval(M[rl]+cl, value, ncol*nrow);
}

//----------------------------------------------------------------------------//

void IntMatrix::Set (int_mapper fcn)
//
//  Sets all matrix elements by applying the given function
//
{
    register int n = ncol*nrow;
    register int* src = M[rl]+cl;
    while (n--) {*src = (*fcn)(*src); src++;}
}

//----------------------------------------------------------------------------//

void IntMatrix::Set (int_matrix_index_mapper fcn)
//
//  Sets all matrix elements by applying the given function.
// The function takes the index pair and the value as arguments.
// The vector will be changed inplace.
//
{
    register int r,c;
    register int* src = M[rl]+cl;
 
    for (r = rl; r <= rh; r++) 
      for (c = cl; c <= ch; c++) {
	  *src = (*fcn)(r,c,*src); src++;
      }
}

//----------------------------------------------------------------------------//

IntMatrix& IntMatrix::ShiftIndex (int row, int col)
//
// Shift the index range of the matrix by the given offsets
// in columns and rows. The matrix will be changed inplace.
//
{
    register int **v;
    register int n;

    if (this->Empty()) 
      Matpack.Error("IntMatrix::ShiftIndex: matrix is not allocated");

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

IntVector IntMatrix::Row (int i) const
//
// Returns the i-th row of this matrix as a vector.
//
{
    if (i < rl || i > rh) 
      Matpack.Error("IntMatrix::Row(int): index out of range (%d)\n",i); 

    IntVector A(cl,ch);
    copyvec(A.V+cl,M[i]+cl,ncol);
    return A.Value();
}

//----------------------------------------------------------------------------//

IntVector IntMatrix::Column (int i) const
//
// Returns the i-th column of this matrix as a vector.
//
{
    if (i < cl || i > ch) 
      Matpack.Error("IntMatrix::Column(int): index out of range (%d)\n",i); 

    IntVector A(rl,rh);
    register int **m = M + rl;
    register int  *v = A.V + rl;
    register int rsize = nrow;
    while (rsize--) *v++ = (*m++)[i];
    return A.Value();
}

//----------------------------------------------------------------------------//

IntMatrix FlipUpDown (const IntMatrix& A)
//
// Returns the matrix with rows flipped, 
// i.e. up and down is mirrored.
//
{    
    int arl = A.rl, arh = A.rh,
        acl = A.cl, ach = A.ch,
	n = arl + arh,
        k = Nint(0.5*n);

    IntMatrix B(arl,arh,acl,ach); 
    // attributes changed: result is general matrix

    // avoid call to index operator that optimizes very badely
    int **a = A.M, **b = B.M;

    for (int i = arl; i <= k; i++) 
	for (int j = acl; j <= ach; j++) {
	    b[i][j] = a[n-i][j];
	    b[n-i][j] = a[i][j];
	}

    return B.Value();
}


//----------------------------------------------------------------------------//

IntMatrix FlipLeftRight (const IntMatrix& A)
//
// Returns the matrix with columns flipped, 
// i.e. left and right is mirrored.
//
{    
    int arl = A.rl, arh = A.rh,
        acl = A.cl, ach = A.ch,
	n = acl + ach,
        k = Nint(0.5*n);

    IntMatrix B(arl,arh,acl,ach); 
    // attributes changed: result is general matrix

    // avoid call to index operator that optimizes very badely
    int **a = A.M, **b = B.M;

    for (int i = arl; i <= arh; i++) 
	for (int j = acl; j <= k; j++) {
	    b[i][j] = a[i][n-j];
	    b[i][n-j] = a[i][j];
	}

    return B.Value();
}

//----------------------------------------------------------------------------//
