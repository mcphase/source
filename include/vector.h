/*-----------------------------------------------------------------------------*\
| vector and matrix classes include file			       vector.h |
|								 	        |
| Last change: Oct 18, 2003							|
|                                                                               |
| Matpack Library Release 1.7.3                                                 |
| Copyright (C) 1991-2003 by Berndt M. Gammel. All rights reserved.             |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.                                                |
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.                                                    |
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.                       |
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration |
| and installation of Matpack.                                                  |
|                                                                               |
\*-----------------------------------------------------------------------------*/


//----------------------------------------------------------------------------//
// Implementation hints:
//
// o  Matrices are represented as a vector of pointers into a vector of 
//    doubles (complex<double>). This improves access times while giving a very
//    elegant way of implementing matrices and vectors with arbitrary
//    offset (not neccessary 0 or 1) without any loss of efficiency.
// o  A reference count scheme is used to prevent garbage on the heap
//    when a matrix or vector variable becomes unbound !
//
//----------------------------------------------------------------------------//


// avoid multiple inclusion
#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "common.h"	// uses common definitions - must be first one

// include STL definitions
#include <complex>	// uses Complex number class
#include <iomanip>	// uses I/O manipulators

using namespace std;

//----------------------------------------------------------------------------//
// classes defined here
//----------------------------------------------------------------------------//
class LogAndSign;
class MatrixBase;
class Vector;          // double precision vector
class Matrix;          // double precision matrix
class ComplexVector;   // complex double precision vector
class ComplexMatrix;   // complex double precision matrix
class IntVector;       // integer vector (non-arithmetic)
class IntMatrix;       // integer matrix  (non-arithmetic)

// attribute values for derived matrices
enum Attributes { General           = 1,
		  UpperTriangular   = 2,
		  LowerTriangular   = 4 };

//----------------------------------------------------------------------------//
// additional inlines for complex numbers
//----------------------------------------------------------------------------//

inline complex<double> log (const complex<double>& b, const complex<double>& x)
    { return log(x)/log(b); }

//----------------------------------------------------------------------------//
// define a class to return the logarithm and sign (+, -, 0) of a number
//----------------------------------------------------------------------------//

class LogAndSign {
  private:
    double log_value;
    int sign;
  public:
    LogAndSign () { log_value=0.0; sign=1; }
    LogAndSign (double);
    void   operator *= (double);
    void   ChangeSign (void) { sign = -sign; }
    double LogValue   (void) const { return log_value; }
    int    Sign       (void) const { return sign; }
    double Value      (void) const { return sign * exp(log_value); }
};

//----------------------------------------------------------------------------//
// define mapper function types
//----------------------------------------------------------------------------//

typedef int     (*int_mapper)     (int&);
typedef double  (*double_mapper)  (double&);
typedef complex<double> (*complex_mapper) (complex<double>&);

typedef int     (*int_vector_index_mapper)     (int,int&);
typedef double  (*double_vector_index_mapper)  (int,double&);
typedef complex<double> (*complex_vector_index_mapper) (int,complex<double>&);
 
typedef int     (*int_matrix_index_mapper)     (int,int,int&);
typedef double  (*double_matrix_index_mapper)  (int,int,double&);
typedef complex<double> (*complex_matrix_index_mapper) (int,int,complex<double>&);

//----------------------------------------------------------------------------//
// Reference block definition
//----------------------------------------------------------------------------//

struct Reference { 
    short count; // reference count
    Reference (void) { count = 1; }
};

//----------------------------------------------------------------------------//
// base class for matrices
//----------------------------------------------------------------------------//

class MatrixBase {   

  protected:
    int rl,rh,cl,ch,ncol,nrow;
    short attribute,temporary,form; 
    struct Reference *D;
    
    MatrixBase (void) { attribute = General; temporary = 0; form = MpTextFormat; }

    void addref (void) const { if (D) ++(D->count); }  
    
  public:
    // member functions
    int Rlo (void) const { return rl; }  
    int Rhi (void) const { return rh; }  
    int Clo (void) const { return cl; }  
    int Chi (void) const { return ch; }
    int Rows (void) const { return nrow; }
    int Cols (void) const { return ncol; }
    int Elements (void) const { return ncol*nrow; }
};

//----------------------------------------------------------------------------//
// double precision version of vector class
//----------------------------------------------------------------------------//

extern Vector NullVector; // used by some library functions

class Vector {

  friend class Matrix;
    
  private:
    double* V;
    int cl,ch,ncol;

    short temporary,form;
    struct Reference *D;

    void addref (void) const { if (D) ++(D->count); }  
    void checkdim (const Vector&);
    friend void checkdim (const Vector&, const Vector&);

  public:

    friend int unbound (const Vector& A) { return (A.temporary == 1); }

    // constructors
    Vector (void);
    Vector (int,int); 
    Vector (int,int,double);
    Vector (int,int,double,double,...);

    // copy constructor
    Vector (const Vector&);

    // destructor
    ~Vector (void);
    
    // indexing
    double& operator[] (int n) { return V[n]; }
    const double& operator[] (int n) const { return V[n]; }

    double& operator() (int n) { 
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (VectorIndexRangeChecking)
	if (n < cl || n > ch)
	  Matpack.Error("Vector: index out of range (%d)",n);	
#endif
	return V[n]; 
    }

    const double& operator() (int n) const { 
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (VectorIndexRangeChecking)
	if (n < cl || n > ch)
	  Matpack.Error("Vector: index out of range (%d)",n);	
#endif
	return V[n]; 
    }

    // create return value
    Vector& Value(void) { temporary = 1; --(D->count); return *this; } 

    // subvector extraction
    Vector  operator () (int,int) const;
    
    // member functions
    int     Lo (void) const { return cl; }  
    int     Hi (void) const { return ch; }

    int     Elements (void) const { return ncol; }

    double* Store (void) { return (D == 0) ? (double*)0 : V+cl; }
    const double* Store (void) const { return (D == 0) ? (double*)0 : V+cl; }

    int     Empty (void) const { return (D == 0); }

    Vector  Apply(double_mapper) const;
    Vector  Apply(double_vector_index_mapper) const;

    void    Remove (void);
    void    Resize (int,int);

    void    Set(double);
    void    Set(double_mapper);  
    void    Set(double_vector_index_mapper);  

    Vector& ShiftIndex(int);

    int     Format (void)  const { return form; }
    Vector& Format (int format) { form = format; return *this; }

    // if(..) and if(!..)
    operator int () const;
    int operator !() const;
    
    // assignment
    Vector& operator =  (double);
    Vector& operator << (const double*);
    Vector& operator += (double);
    Vector& operator -= (double);
    Vector& operator *= (double);
    Vector& operator /= (double);
    Vector& operator =  (const Vector&);
    Vector& operator += (const Vector&);
    Vector& operator -= (const Vector&);
    Vector& operator %= (const Vector&);
    
    // friend operators
    friend int    operator == (double,  const Vector&);
    friend int    operator == (const Vector&, double);
    friend int    operator == (const Vector&, const Vector&);
    friend int    operator != (double, const Vector&);
    friend int    operator != (const Vector&, double);
    friend int    operator != (const Vector&, const Vector&);
    friend Vector operator +  (const Vector&, double);
    friend Vector operator +  (double, const Vector&);
    friend Vector operator +  (const Vector&, const Vector&);
    friend Vector operator -  (const Vector&);
    friend Vector operator -  (const Vector&, double);
    friend Vector operator -  (double, const Vector&);
    friend Vector operator -  (const Vector&, const Vector&);
    friend Vector operator *  (const Vector&, double);
    friend Vector operator *  (double, const Vector&);
    friend Vector operator /  (const Vector&, double);
    friend Vector operator /  (double, const Vector&);
    friend double operator *  (const Vector&, const Vector&);
    friend Vector operator *  (const Vector&, const Matrix&);
    friend Vector operator *  (const Matrix&, const Vector&);
    friend Matrix operator ^  (const Vector&, const Vector&);
    friend Vector operator %  (const Vector&, const Vector&);
  
    friend ostream& operator << (ostream&, const Vector&);
    friend istream& operator >> (istream&, Vector&);

    // friend functions
    friend Matrix Diagonal (const Vector&);
    friend ComplexMatrix ComplexDiagonal (const Vector&);
    friend double Norm     (const Vector&);
    friend double Norm1    (const Vector&);
    friend double Norm2    (const Vector&);
    friend double NormInf  (const Vector&);
    friend double Sum      (const Vector&);
    friend double Min      (const Vector&);
    friend double Max      (const Vector&);
    friend double Min      (const Vector&,int&);
    friend double Max      (const Vector&,int&);
    friend void   Sort     (Vector&);
    friend void   Sort     (Vector&,Vector&);
    friend void   Reverse  (Vector&); 

    friend Vector Pack     (const Matrix&);

    friend Vector Real     (const ComplexVector&);
    friend Vector Imag     (const ComplexVector&);
    friend Vector Abs      (const ComplexVector&);
    friend Vector AbsSqr   (const ComplexVector&);
    friend Vector Arg      (const ComplexVector&);

    friend Vector Cos      (const Vector&);
    friend Vector Sin      (const Vector&);
    friend Vector Cosh     (const Vector&);
    friend Vector Sinh     (const Vector&);
    friend Vector Exp      (const Vector&);
    friend Vector Log      (const Vector&);
    friend Vector Sqrt     (const Vector&);
    friend Vector Sqr      (const Vector&);
    friend Vector Cube     (const Vector&);
    
    friend int    MatchingIndexRange (const Vector&,const Vector&);
};

//----------------------------------------------------------------------------//
// double precision version of matrix class
//----------------------------------------------------------------------------//

extern Matrix NullMatrix; // used by some library functions

class Matrix : public MatrixBase {

  friend class Vector;
  friend class ComplexMatrix;
    
  private:
    double** M;

    void checkdim (const Matrix&);
    friend void checkdim    (const Matrix&, const Matrix&);
    friend void checksquare (const Matrix&);

  public:

    friend int unbound (const Matrix& A) { return (A.temporary == 1); }

    // constructors
    Matrix (void);
    Matrix (int,int,int,int);
    Matrix (int,int,int,int,double);
    Matrix (int,int,int,int,double,double,...);

    // copy constructor
    Matrix (const Matrix&);
    
    // destructor
    ~Matrix (void);
    
    // indexing
    double*& operator[] (int n) { return M[n]; }
    //  const double* const & operator[] (int n) const { return (const double*)(M[n]); } // worked until gcc 3.2
    double*& operator[] (int n) const { return M[n]; }

    double& operator() (int n, int m) {
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (MatrixIndexRangeChecking)
	if (n < rl || n > rh || m < cl || m > ch)
	  Matpack.Error("Matrix: index out of range (%d,%d)", n,m);	
#endif
	return M[n][m];
    }

    const double& operator() (int n, int m) const {
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (MatrixIndexRangeChecking)
	if (n < rl || n > rh || m < cl || m > ch)
	  Matpack.Error("Matrix: index out of range (%d,%d)", n,m);	
#endif
	return M[n][m];
    }
    
    // create return value
    Matrix& Value(void) { temporary = 1; --(D->count); return *this; } 

    // submatrix extraction
    Matrix  operator() (int,int,int,int) const;
    
    // member functions

    double*  Store (void) { return (D == 0) ? (double*)0 : M[rl]+cl; }
    const double*  Store (void) const { return (D == 0) ? (double*)0 : M[rl]+cl; }

    int      Empty (void) const { return (D == 0); }
    Vector   Row (int) const;
    Vector   Column (int) const;
    Matrix   Transpose (void) const;
    Matrix   Inverse (void) const;
    Matrix   Apply(double_mapper) const;
    Matrix   Apply(double_matrix_index_mapper) const;

    int      Attribute (int); // const TODO ??????????
    int      Attribute (void) const { return attribute; }

    int      Format (void)  const { return form; }
    Matrix&  Format (int format) { form = format; return *this; }

    void     Remove (void);
    void     Resize (int,int,int,int);

    void     Set(double);
    void     Set(double_mapper);
    void     Set(double_matrix_index_mapper);

    Matrix&  ShiftIndex(int,int);

    // if(..) and if(!..)
    operator int () const;
    int operator !() const;
    
    // assignment
    Matrix& operator =  (double);
    Matrix& operator <<  (const double*);
    Matrix& operator += (double);
    Matrix& operator -= (double);
    Matrix& operator *= (double);
    Matrix& operator /= (double);
    Matrix& operator =  (const Matrix&);
    Matrix& operator += (const Matrix&);
    Matrix& operator -= (const Matrix&);
    Matrix& operator %= (const Matrix&);
    
    // friend operators
    friend int    operator == (double, const Matrix&);
    friend int    operator == (const Matrix&, double);
    friend int    operator == (const Matrix&, const Matrix&);
    friend int    operator != (double, const Matrix&);
    friend int    operator != (const Matrix&, double);
    friend int    operator != (const Matrix&, const Matrix&);
    friend Matrix operator +  (const Matrix&, double);
    friend Matrix operator +  (double, const Matrix&); 
    friend Matrix operator +  (const Matrix&, const Matrix&);
    friend Matrix operator -  (const Matrix&);
    friend Matrix operator -  (const Matrix&, double);
    friend Matrix operator -  (double, const Matrix&);
    friend Matrix operator -  (const Matrix&, const Matrix&);
    friend Matrix operator *  (const Matrix&, double);
    friend Matrix operator *  (double, const Matrix&);
    friend Matrix operator /  (const Matrix&, double);
    friend Matrix operator /  (double, const Matrix&);
    friend Vector operator *  (const Vector&, const Matrix&);
    friend Vector operator *  (const Matrix&, const Vector&);
    friend Matrix operator *  (const Matrix&, const Matrix&);   

    friend Matrix operator %  (const Matrix&, const Matrix&);   

    friend Matrix operator &  (const Matrix&, const Matrix&);   // test version

    friend ostream& operator << (ostream&, const Matrix&);
    friend istream& operator >> (istream&, Matrix&);

    // friend functions
    friend Matrix     Diagonal (const Vector& d);
    friend double     Trace   (const Matrix&);
    friend double     Norm    (const Matrix&);
    friend double     Norm1   (const Matrix&);
    friend double     Norm2   (const Matrix&);
    friend double     NormFro (const Matrix&);
    friend double     NormInf (const Matrix&);
    friend double     Sum     (const Matrix&);
    friend double     Min     (const Matrix&);
    friend double     Max     (const Matrix&);
    friend double     Min     (const Matrix&,int&,int&);
    friend double     Max     (const Matrix&,int&,int&);
    friend double     Det     (const Matrix&);
    friend LogAndSign LogDet  (const Matrix&);

    friend Vector     Pack    (const Matrix&);

    friend Matrix     Real    (const ComplexMatrix&);
    friend Matrix     Imag    (const ComplexMatrix&);
    friend Matrix     Abs     (const ComplexMatrix&);
    friend Matrix     AbsSqr  (const ComplexMatrix&);
    friend Matrix     Arg     (const ComplexMatrix&);

    friend Matrix     Cos   (const Matrix&);
    friend Matrix     Sin   (const Matrix&);
    friend Matrix     Cosh  (const Matrix&);
    friend Matrix     Sinh  (const Matrix&);
    friend Matrix     Exp   (const Matrix&);
    friend Matrix     Log   (const Matrix&);
    friend Matrix     Sqrt  (const Matrix&);
    friend Matrix     Sqr   (const Matrix&);
    friend Matrix     Cube  (const Matrix&);

    friend Matrix     FlipLeftRight (const Matrix&);
    friend Matrix     FlipUpDown    (const Matrix&);

    friend int        MatchingIndexRange (const Matrix&,const Matrix&);
};

//----------------------------------------------------------------------------//
// complex version of vector class
//----------------------------------------------------------------------------//

class ComplexVector {

  friend class ComplexMatrix;
  
  private:
    complex<double>* V;
    int cl,ch,ncol;
    
    short temporary,form;
    struct Reference *D;
    
    void addref (void) const { if (D) ++(D->count); }  
    void checkdim (const ComplexVector&);
    friend void checkdim (const ComplexVector&, const ComplexVector&);
    
  public:

    friend int unbound (const ComplexVector& A) { return (A.temporary == 1); }

    // constructors
    ComplexVector (void);
    ComplexVector (int,int); 
    ComplexVector (int,int,complex<double>);
    explicit ComplexVector (Vector&);
    ComplexVector (Vector&,Vector&);
    
    // copy constructor
    ComplexVector (const ComplexVector&);
    
    // destructor
    ~ComplexVector (void);
    
    // indexing
    complex<double>& operator[] (int n) { return V[n]; }
    const complex<double>& operator[] (int n) const { return V[n]; }

    complex<double>& operator() (int n) {
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (ComplexVectorIndexRangeChecking)
	if (n < cl || n > ch)
	  Matpack.Error("ComplexVector: index out of range (%d)", n);	
#endif
	return V[n]; 
    }

    const complex<double>& operator() (int n) const {
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (ComplexVectorIndexRangeChecking)
	if (n < cl || n > ch)
	  Matpack.Error("ComplexVector: index out of range (%d)", n);	
#endif
	return V[n]; 
    }
    
    // create return value
    ComplexVector& Value(void) { temporary = 1; --(D->count); return *this; } 

    // subvector extraction
    ComplexVector  operator() (int,int) const;
    
    // member functions
    int           Lo (void) const { return cl; }  
    int           Hi (void) const { return ch; }

    int           Elements (void) const { return ncol; }
    complex<double>* Store (void) { return (D == 0) ? (complex<double>*)0 : V+cl; }
    const complex<double>* Store (void) const { return (D == 0) ? (complex<double>*)0 : V+cl; }

    int           Empty (void) const { return (D == 0); }
    ComplexVector Conjugate (void) const;

    ComplexVector Apply(complex_mapper) const;
    ComplexVector Apply(complex_vector_index_mapper) const;

    void          Remove (void);
    void          Resize (int,int);

    void          Set(complex<double>);
    void          Set(complex_mapper);
    void          Set(complex_vector_index_mapper);

    ComplexVector& ShiftIndex(int);

    int            Format (void)  const { return form; }
    ComplexVector& Format (int format) { form = format; return *this; }

    // if(..) and if(!..)   
    operator int () const;
    int operator !() const;

    // assignment
    ComplexVector& operator =  (complex<double>);
    ComplexVector& operator << (const complex<double>*);
    ComplexVector& operator += (complex<double>);
    ComplexVector& operator -= (complex<double>);
    ComplexVector& operator *= (complex<double>);
    ComplexVector& operator /= (complex<double>);
    ComplexVector& operator =  (const ComplexVector&);
    ComplexVector& operator -= (const ComplexVector&);
    ComplexVector& operator += (const ComplexVector&);
    ComplexVector& operator %= (const ComplexVector&);
    
    // friend operators
    friend int           operator == (complex<double>, const ComplexVector&);
    friend int           operator == (const ComplexVector&, complex<double>);
    friend int           operator == (const ComplexVector&, const ComplexVector&);
    friend int           operator != (complex<double>, const ComplexVector&);
    friend int           operator != (const ComplexVector&, complex<double>);
    friend int           operator != (const ComplexVector&, const ComplexVector&);
    friend ComplexVector operator +  (const ComplexVector&, complex<double>);
    friend ComplexVector operator +  (complex<double>, const ComplexVector&);
    friend ComplexVector operator +  (const ComplexVector&, const ComplexVector&);
    friend ComplexVector operator -  (const ComplexVector&);
    friend ComplexVector operator -  (const ComplexVector&, complex<double>);
    friend ComplexVector operator -  (complex<double>, const ComplexVector&);
    friend ComplexVector operator -  (const ComplexVector&, const ComplexVector&);
    friend ComplexVector operator *  (const ComplexVector&, complex<double>);
    friend ComplexVector operator *  (complex<double>, const ComplexVector&);
    friend ComplexVector operator /  (const ComplexVector&, complex<double>);
    friend ComplexVector operator /  (complex<double>, const ComplexVector&);
    friend complex<double> operator * (const ComplexVector&, const ComplexVector&);
    friend ComplexVector operator *  (const ComplexVector&, const ComplexMatrix&);
    friend ComplexVector operator *  (const ComplexMatrix&, const ComplexVector&);
    friend ComplexMatrix operator ^  (const ComplexVector&, const ComplexVector&);
    friend ComplexVector operator %  (const ComplexVector&, const ComplexVector&);

    friend ostream&      operator << (ostream&, const ComplexVector&);
    friend istream&      operator >> (istream&, ComplexVector&);
       
    // friend functions
    friend ComplexMatrix         Diagonal (const ComplexVector&);
    friend double                Norm     (const ComplexVector&);
    friend double                Norm1    (const ComplexVector&);
    friend double                Norm2    (const ComplexVector&);
    friend double                NormInf  (const ComplexVector&);
    friend complex<double>       Sum      (const ComplexVector&);
    friend complex<double>       Min      (const ComplexVector&);
    friend complex<double>       Max      (const ComplexVector&);
    friend complex<double>       Min      (const ComplexVector&,int&);
    friend complex<double>       Max      (const ComplexVector&,int&);
    friend void         	 Sort     (ComplexVector&);
    friend void  	         Reverse  (ComplexVector&); 

    friend ComplexVector Pack     (const ComplexMatrix&);

    friend Vector        Real     (const ComplexVector&);
    friend Vector        Imag     (const ComplexVector&);
    friend Vector        Abs      (const ComplexVector&);
    friend Vector        AbsSqr   (const ComplexVector&);
    friend Vector        Arg      (const ComplexVector&);

    friend ComplexVector Cos      (const ComplexVector&);
    friend ComplexVector Sin      (const ComplexVector&);
    friend ComplexVector Cosh     (const ComplexVector&);
    friend ComplexVector Sinh     (const ComplexVector&);
    friend ComplexVector Exp      (const ComplexVector&);
    friend ComplexVector Log      (const ComplexVector&);
    friend ComplexVector Sqrt     (const ComplexVector&);
    friend ComplexVector Sqr      (const ComplexVector&);
    friend ComplexVector Cube     (const ComplexVector&);

    friend int MatchingIndexRange (const ComplexVector&,const ComplexVector&);

  protected:

    // added for name inference purposes
    friend inline int           operator == (const ComplexVector& v, double x)
      { return operator == (v,complex<double>(x)); }
    friend inline int           operator == (double x, const ComplexVector& v)
      { return operator == (complex<double>(x),v); }
    friend inline int           operator != (const ComplexVector& v, double x)
      { return operator != (v,complex<double>(x)); }
    friend inline int           operator != (double x, const ComplexVector& v)
      { return operator != (complex<double>(x),v); }
    friend inline ComplexVector operator +  (const ComplexVector& v, double x)
      { return operator + (v,complex<double>(x)); }
    friend inline ComplexVector operator +  (double x, const ComplexVector& v)
      { return operator + (complex<double>(x),v); }
    friend inline ComplexVector operator -  (const ComplexVector& v, double x)
      { return operator - (v,complex<double>(x)); }
    friend inline ComplexVector operator -  (double x, const ComplexVector& v)
      { return operator - (complex<double>(x),v); }
    friend inline ComplexVector operator *  (const ComplexVector& v, double x)
      { return operator * (v,complex<double>(x)); }
    friend inline ComplexVector operator *  (double x, const ComplexVector& v)
      { return operator * (complex<double>(x),v); }
    friend inline ComplexVector operator /  (const ComplexVector& v, double x)
      { return operator / (v,complex<double>(x)); }
    friend inline ComplexVector operator /  (double x, const ComplexVector& v)
      { return operator / (complex<double>(x),v); }
};

//----------------------------------------------------------------------------//
// complex version of matrix class
//----------------------------------------------------------------------------//

class ComplexMatrix : public MatrixBase {

  friend class ComplexVector;
  
  private:
    complex<double>** M;

    void checkdim (const ComplexMatrix&);
    friend void checkdim (const ComplexMatrix&, const ComplexMatrix&);
    friend void checksquare (const ComplexMatrix&);
    
  public:

    friend int unbound (const ComplexMatrix& A) { return (A.temporary == 1); }

    // constructors
    ComplexMatrix (void);
    ComplexMatrix (int,int,int,int);
    ComplexMatrix (int,int,int,int,complex<double>);
    explicit ComplexMatrix (Matrix&);
    ComplexMatrix (Matrix&,Matrix&);
    
    // copy constructor
    ComplexMatrix (const ComplexMatrix&);
    
    // destructor
    ~ComplexMatrix (void);
    
    // indexing
    complex<double>*& operator[] (int n) { return M[n]; }
    // const complex<double>* const & operator [] (int n) const { return (const complex<double>*)(M[n]); } // worked until gcc 3.2
    complex<double>*& operator [] (int n) const { return M[n]; }

    complex<double>& operator() (int n, int m) {
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (ComplexMatrixIndexRangeChecking)
	if (n < rl || n > rh || m < cl || m > ch)
	  Matpack.Error("ComplexMatrix: index out of range (%d,%d)", n,m);	
#endif
	return M[n][m];
    }

    const complex<double>& operator() (int n, int m) const {
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (ComplexMatrixIndexRangeChecking)
	if (n < rl || n > rh || m < cl || m > ch)
	  Matpack.Error("ComplexMatrix: index out of range (%d,%d)", n,m);	
#endif
	return M[n][m];
    }
    
    // create return value
    ComplexMatrix& Value(void) { temporary = 1; --(D->count); return *this; } 

    // submatrix extraction
    ComplexMatrix  operator() (int,int,int,int) const;
    
    // member functions
    complex<double>* Store (void) { return (D==0) ? (complex<double>*)0 : M[rl]+cl; }
    const complex<double>* Store (void) const { return (D==0) ? (complex<double>*)0 : M[rl]+cl; }

    int            Empty (void) const { return (D == 0); }
    ComplexVector  Row (int) const;
    ComplexVector  Column (int) const;
    ComplexMatrix  Transpose (void) const;
    ComplexMatrix  Inverse (void) const;
    ComplexMatrix  Conjugate (void) const;
    ComplexMatrix  Hermitean (void) const;

    ComplexMatrix  Apply(complex_mapper) const;
    ComplexMatrix  Apply(complex_matrix_index_mapper) const;

    int            Attribute (int);
    int            Attribute (void) const { return attribute; }

    void           Remove (void);
    void           Resize (int,int,int,int);

    void           Set(complex<double>);
    void           Set(complex_mapper);
    void           Set(complex_matrix_index_mapper);

    ComplexMatrix& ShiftIndex(int,int);

    int            Format (void)  const { return form; }
    ComplexMatrix& Format (int format) { form = format; return *this; }
    
    // if(..) and if(!..)
    operator int () const;
    int operator !() const;

    // assignment
    ComplexMatrix& operator =  (complex<double>);
    ComplexMatrix& operator << (const complex<double>*);
    ComplexMatrix& operator += (complex<double>);
    ComplexMatrix& operator -= (complex<double>);
    ComplexMatrix& operator *= (complex<double>);
    ComplexMatrix& operator /= (complex<double>);
    ComplexMatrix& operator =  (const ComplexMatrix&);
    ComplexMatrix& operator += (const ComplexMatrix&);
    ComplexMatrix& operator -= (const ComplexMatrix&);
    ComplexMatrix& operator %= (const ComplexMatrix&);
   
    // friend operators
    friend int           operator == (complex<double>, const ComplexMatrix&);
    friend int           operator == (const ComplexMatrix&, complex<double>);
    friend int           operator == (const ComplexMatrix&, const ComplexMatrix&);
    friend int           operator != (complex<double>, const ComplexMatrix&);
    friend int           operator != (const ComplexMatrix&, complex<double>);
    friend int           operator != (const ComplexMatrix&, const ComplexMatrix&);
    friend ComplexMatrix operator +  (const ComplexMatrix&, complex<double>);
    friend ComplexMatrix operator +  (complex<double>, const ComplexMatrix&);
    friend ComplexMatrix operator +  (const ComplexMatrix&, const ComplexMatrix&);
    friend ComplexMatrix operator -  (const ComplexMatrix&);
    friend ComplexMatrix operator -  (const ComplexMatrix&, complex<double>);
    friend ComplexMatrix operator -  (complex<double>, const ComplexMatrix&);
    friend ComplexMatrix operator -  (const ComplexMatrix&, const ComplexMatrix&);
    friend ComplexMatrix operator *  (const ComplexMatrix&, complex<double>);
    friend ComplexMatrix operator *  (complex<double>, const ComplexMatrix&);
    friend ComplexMatrix operator /  (const ComplexMatrix&, complex<double>);
    friend ComplexMatrix operator /  (complex<double>, const ComplexMatrix&);
    friend ComplexVector operator *  (const ComplexVector&, const ComplexMatrix&);
    friend ComplexVector operator *  (const ComplexMatrix&, const ComplexVector&);
    friend ComplexMatrix operator *  (const ComplexMatrix&, const ComplexMatrix&);
    friend ComplexVector operator *  (Vector&, const ComplexMatrix&);
    friend ComplexVector operator *  (const ComplexMatrix&, Vector&);
    friend ComplexMatrix operator %  (const ComplexMatrix&, const ComplexMatrix&);

    // test version
    friend ComplexMatrix operator &  (const ComplexMatrix&, const ComplexMatrix&);

    friend ostream&      operator << (ostream&, const ComplexMatrix&);
    friend istream&      operator >> (istream&, ComplexMatrix&);
     
    // friend functions
    friend ComplexMatrix   Diagonal (const ComplexVector&);
    friend complex<double> Trace   (const ComplexMatrix&);
    friend double  	   Norm    (const ComplexMatrix&);
    friend double  	   Norm1   (const ComplexMatrix&);
    friend double  	   Norm2   (const ComplexMatrix&);
    friend double  	   NormFro (const ComplexMatrix&);
    friend double	   NormInf (const ComplexMatrix&);
    friend complex<double> Sum     (const ComplexMatrix&);
    friend complex<double> Min     (const ComplexMatrix&);
    friend complex<double> Max     (const ComplexMatrix&);
    friend complex<double> Min     (const ComplexMatrix&,int&,int&);
    friend complex<double> Max     (const ComplexMatrix&,int&,int&);
    friend complex<double> Det     (const ComplexMatrix&);

    friend ComplexVector   Pack  (const ComplexMatrix&);
    friend ComplexMatrix   ComplexDiagonal (const Vector&);

    friend Matrix  Real  (const ComplexMatrix&);
    friend Matrix  Imag  (const ComplexMatrix&);
    friend Matrix  Abs   (const ComplexMatrix&);
    friend Matrix  AbsSqr(const ComplexMatrix&);
    friend Matrix  Arg   (const ComplexMatrix&);

    friend ComplexMatrix Cos      (const ComplexMatrix&);
    friend ComplexMatrix Sin      (const ComplexMatrix&);
    friend ComplexMatrix Cosh     (const ComplexMatrix&);
    friend ComplexMatrix Sinh     (const ComplexMatrix&);
    friend ComplexMatrix Exp      (const ComplexMatrix&);
    friend ComplexMatrix Log      (const ComplexMatrix&);
    friend ComplexMatrix Sqrt     (const ComplexMatrix&);
    friend ComplexMatrix Sqr      (const ComplexMatrix&);
    friend ComplexMatrix Cube     (const ComplexMatrix&);
 
    friend ComplexMatrix FlipLeftRight (const ComplexMatrix&);
    friend ComplexMatrix FlipUpDown    (const ComplexMatrix&);

    friend int MatchingIndexRange (const ComplexMatrix&,const ComplexMatrix&);
  
  protected:

    // added for name inference purposes
    friend inline int           operator == (double x, const ComplexMatrix& C)
       { return operator == (complex<double>(x),C); }
    friend inline int           operator == (const ComplexMatrix& C, double x)
       { return operator == (C,complex<double>(x)); }
    friend inline int           operator != (double x, const ComplexMatrix& C)
       { return operator != (complex<double>(x),C); }
    friend inline int           operator != (const ComplexMatrix& C, double x)
       { return operator != (C,complex<double>(x)); }
    friend inline ComplexMatrix operator +  (double x, const ComplexMatrix& C)
       { return operator + (complex<double>(x),C); }
    friend inline ComplexMatrix operator +  (const ComplexMatrix& C, double x)
       { return operator + (C,complex<double>(x)); }
    friend inline ComplexMatrix operator -  (double x, const ComplexMatrix& C)
       { return operator - (complex<double>(x),C); }
    friend inline ComplexMatrix operator -  (const ComplexMatrix& C, double x)
       { return operator - (C,complex<double>(x)); }
    friend inline ComplexMatrix operator *  (double x, const ComplexMatrix& C)
       { return operator * (complex<double>(x),C); }
    friend inline ComplexMatrix operator *  (const ComplexMatrix& C, double x)
       { return operator * (C,complex<double>(x)); }
    friend inline ComplexMatrix operator /  (double x, const ComplexMatrix& C)
       { return operator / (complex<double>(x),C); }
    friend inline ComplexMatrix operator /  (const ComplexMatrix& C, double x)
       { return operator / (C,complex<double>(x)); }
};

//----------------------------------------------------------------------------//
// integer version of vector class
//----------------------------------------------------------------------------//
// This is a basically non-arithmetic vector class. That means no arithmetics
// can be performed with integer vectors. They are only intended for indexing
// purposes, sorting, book keeping, etc... but not for arithmetic combination 
// with double or complex vectors. Only indexing, assignment, comparison,I/O,
// sorting and extrema finding is implemented.
//----------------------------------------------------------------------------//

class IntVector {

  friend class IntMatrix;

  private:
    int* V;
    int cl,ch,ncol;

    short temporary,form;
    struct Reference *D;

    void addref (void) const { if (D) ++(D->count); }  
    void checkdim (const IntVector&);
    friend void checkdim (const IntVector&, const IntVector&);
    friend int unbound (const IntVector& A) { return (A.temporary == 1); }

  public:

    // constructors
    IntVector (void);
    IntVector (int,int); 
    IntVector (int,int,int);

    // copy constructor
    IntVector (const IntVector&);

    // destructor
    ~IntVector (void);
    
    // indexing
    int& operator[] (int n) { return V[n]; }
    const int& operator[] (int n) const { return V[n]; }

    int& operator() (int n) { 
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (IntVectorIndexRangeChecking)
	if (n < cl || n > ch)
	  Matpack.Error("IntVector: index out of range (%d)", n);	
#endif
	return V[n]; 
    }

    const int& operator() (int n) const { 
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (IntVectorIndexRangeChecking)
	if (n < cl || n > ch)
	  Matpack.Error("IntVector: index out of range (%d)", n);	
#endif
	return V[n]; 
    }

    // create return value
    IntVector& Value(void) { temporary = 1; --(D->count); return *this; } 

    // subvector extraction
    IntVector  operator() (int,int) const;
    
    // member functions
    int     Lo (void) const { return cl; }  
    int     Hi (void) const { return ch; }

    int     Elements (void) const { return ncol; }
    int*    Store (void) { return (D == 0) ? (int*)0 : V+cl; }
    const int*    Store (void) const { return (D == 0) ? (int*)0 : V+cl; }

    int     Empty (void) const { return (D == 0); }

    void    Remove (void);
    void    Resize (int,int);

    int     Format (void)  const { return form; }
    IntVector& Format (int format) { form = format; return *this; }
 
    void    Set(int);
    void    Set(int_mapper);
    void    Set(int_vector_index_mapper);

    IntVector& ShiftIndex(int);

    // if(..) and if(!..)
    operator int () const;
    int operator !() const;
    
    // assignment
    IntVector& operator =  (int);
    IntVector& operator << (const int*);
    IntVector& operator =  (const IntVector&);
    
    // friend operators
    friend int    operator == (int, const IntVector&);
    friend int    operator == (const IntVector&, int);
    friend int    operator == (const IntVector&, const IntVector&);
    friend int    operator != (int, const IntVector&);
    friend int    operator != (const IntVector&,int);
    friend int    operator != (const IntVector&, const IntVector&);

    // friend functions
    friend int    Min      (const IntVector&);
    friend int    Max      (const IntVector&);
    friend int    Min      (const IntVector&,int&);
    friend int    Max      (const IntVector&,int&);
    friend void   Sort     (IntVector&);
    friend void   Reverse  (IntVector&); 
    friend IntVector  Pack  (const IntMatrix&);
    friend int    MatchingIndexRange (const IntVector&,const IntVector&);

    // stream i/o functions
    int read  (istream& is);
    int write (ostream& os) const;

    friend ostream& operator << (ostream&, const IntVector&);
    friend istream& operator >> (istream&, IntVector&);
};


//----------------------------------------------------------------------------//
// integer version of matrix class
//----------------------------------------------------------------------------//
// This is a basically non-arithmetic matrix class. That means no arithmetics
// can be performed with integer matrices. They are only intended for indexing
// purposes, sorting, book keeping, etc... but not for arithmetic combination 
// with double or complex matricess. Only indexing, assignment, comparison,I/O,
// and extrema finding is implemented.
//----------------------------------------------------------------------------//

class IntMatrix : public MatrixBase {

  friend class IntVector;
    
  private:
    int** M;
 
    void checkdim (const IntMatrix&);
    friend void checkdim    (const IntMatrix&, const IntMatrix&);
    friend void checksquare (const IntMatrix&);

  public:

    friend int unbound (const IntMatrix& A) { return (A.temporary == 1); }

    // constructors
    IntMatrix (void);
    IntMatrix (int,int,int,int);
    IntMatrix (int,int,int,int,int);

    // copy constructor
    IntMatrix (const IntMatrix&);
    
    // destructor
    ~IntMatrix (void);
    
    // indexing
    int*& operator[] (int n) { return M[n]; }
    // const int* const & operator[] (int n) const { return (const int*)(M[n]); } // worked until gcc 3.2
    int*& operator[] (int n) const { return M[n]; }

    int& operator() (int n, int m) {
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (MatrixIndexRangeChecking)
	if (n < rl || n > rh || m < cl || m > ch)
	  Matpack.Error("IntMatrix: index out of range (%d,%d)", n,m);	
#endif
	return M[n][m];
    }

    const int& operator() (int n, int m) const {
#if defined(DEBUG) || defined (IndexRangeChecking) || defined (MatrixIndexRangeChecking)
	if (n < rl || n > rh || m < cl || m > ch)
	  Matpack.Error("IntMatrix: index out of range (%d,%d)", n,m);	
#endif
	return M[n][m];
    }
    
    // create return value
    IntMatrix& Value(void) { temporary = 1; --(D->count); return *this; } 

    // submatrix extraction
    IntMatrix  operator() (int,int,int,int) const;
    
    // member functions
    int*       Store (void) { return (D == 0) ? (int*)0 : M[rl]+cl; }
    const int* Store (void) const { return (D == 0) ? (int*)0 : M[rl]+cl; }

    int        Empty (void) const { return (D == 0); }
    IntVector  Row (int) const;
    IntVector  Column (int) const;

    int        Format (void)  const { return form; }
    IntMatrix& Format (int format) { form = format; return *this; }

    void       Remove (void);
    void       Resize (int,int,int,int);

    void       Set(int);
    void       Set(int_mapper);
    void       Set(int_matrix_index_mapper);

    IntMatrix& ShiftIndex(int,int);

    // if(..) and if(!..)
    operator int () const;
    int operator !() const;
    
    // assignment
    IntMatrix& operator =  (int);
    IntMatrix& operator << (const int*);
    IntMatrix& operator =  (const IntMatrix&);
    
    // friend operators
    friend int    operator == (int, const IntMatrix&);
    friend int    operator == (const IntMatrix&, int);
    friend int    operator == (const IntMatrix&, const IntMatrix&);
    friend int    operator != (int, const IntMatrix&);
    friend int    operator != (const IntMatrix&, int);
    friend int    operator != (const IntMatrix&, const IntMatrix&);

    friend ostream& operator << (ostream&, const IntMatrix&);
    friend istream& operator >> (istream&, IntMatrix&);

    // friend functions
    friend int        Min   (const IntMatrix&);
    friend int        Max   (const IntMatrix&);
    friend int        Min   (const IntMatrix&,int&,int&);
    friend int        Max   (const IntMatrix&,int&,int&);
    friend IntVector  Pack  (const IntMatrix&);
    friend int        MatchingIndexRange (const IntMatrix&, const IntMatrix&);

    friend IntMatrix  FlipLeftRight (const IntMatrix&);
    friend IntMatrix  FlipUpDown    (const IntMatrix&);
};


//----------------------------------------------------------------------------//
// prototypes
//----------------------------------------------------------------------------//

// Integral Transforms
// Note: vector length must be power of 2
// --------------------------------------
void    FFT             (Vector&,Vector&);  // Fast Fourier Transform
void    InverseFFT      (Vector&,Vector&);
void    FFT             (Vector&);
void    InverseFFT      (Vector&);
void    FHT             (Vector&);          // Fast Hartley Transform

// Linear Algebra
// --------------

// real
void	Decompose	(Matrix&,IntVector&,int&);
void	Backsubst	(Matrix&,IntVector&,Vector&);
void	Backsubst	(Matrix&,IntVector&,Matrix&);
void    SolveLinear     (Matrix&,Vector&);
void    SolveLinear     (Matrix&,Matrix&);
void	Improve		(Matrix&,Matrix&,IntVector&,Vector&,Vector&);
void	Ortho		(Matrix&);
void	SVDecompose	(Matrix&,Vector&,Matrix&);
void	SVBacksubst	(Matrix&,Vector&,Matrix&,Vector&,Vector&);
void    SVSolveLinear   (Matrix&,Vector&,Vector&,double&,double=0.0);

// complex
void	Decompose	(ComplexMatrix&,IntVector&,int&);
void	Backsubst	(ComplexMatrix&,IntVector&,ComplexVector&);
void	Backsubst	(ComplexMatrix&,IntVector&,ComplexMatrix&);
void    SolveLinear     (ComplexMatrix&,ComplexVector&);
void    SolveLinear     (ComplexMatrix&,ComplexMatrix&);
void	Improve		(ComplexMatrix&,ComplexMatrix&,IntVector&,
			 ComplexVector&,ComplexVector&);

// Eigensystems
// ------------

// real symmetric (source/Vector/rseigen.cc, tred.cc, imtql.cc, imtql2.cc)
void	Tred		(Matrix&,Vector&,Vector&);
void	Imtql		(Vector&,Vector&,int=true,int=30);
void	Imtql		(Matrix&,Vector&,Vector&,int=true,int=30);
void	EigenSystemSymmetric (Matrix&,Vector&,int=true,int=30);
void    EigenValuesSymmetric (Matrix&,Vector&,int=true,int=30);

// complex hermitean (source/Vector/cheigen.cc)
void	Chtred		(Matrix&,Vector&,Vector&,Vector&,Vector&);
void	Chtrbk		(Matrix&,Vector&,Vector&,Matrix&,Matrix&);
void	EigenSystemHermitean (Matrix&,Vector&,Matrix&,Matrix&,int=true,int=30);
void    EigenValuesHermitean (Matrix&,Vector&,int=true,int=30);

// generalized complex hermitean (source/Vector/chgeigen.cc)
void    Chreduce        (Matrix&,Matrix&,Matrix&,Matrix&);
void    Chreback        (Matrix&,Matrix&,Matrix&);
void    EigenSystemHermiteanGeneral (Matrix&,Matrix&,Vector&,Matrix&,
				     Matrix&,int=true,int=30);

// More matrix functions (source/Vector/matexph.cc, matsqth.cc)
// ---------------------
ComplexMatrix MatrixExpHermitean (const ComplexMatrix& A);
ComplexMatrix MatrixSqrtHermitean (const ComplexMatrix& A);

// Special matrices
// ----------------
Matrix Vandermonde (const Vector& C);
Matrix Hilbert (int N);
Matrix InverseHilbert (int N);
Matrix Wilkinson (int N);
Matrix RandomMatrix (int N, long seed = 0);

IntMatrix MpMagicSquare (int n);

// Miscellaneous
// -------------
Matrix  Laplacian (Matrix& A);

// q-th quantile
// -------------

void MpQuantile (const float *a, unsigned n, unsigned m, 
		 float q, float &pqq, float &pstdv);

void MpQuantile (const double *a, unsigned n, unsigned m, 
		 double q, double &pqq, double &pstdv);


inline void MpQuantile (const Vector& a, int m, double q, double &qq, double &stdv)
{ 
  MpQuantile (a.Store(), a.Elements(), (unsigned) m, q, qq, stdv);  
}

inline void MpQuantile (const Matrix& a, int m, double q, double &qq, double &stdv)
{
  MpQuantile (a.Store(), a.Elements(), (unsigned) m, q, qq, stdv); 
}


#endif

//----------------------------------------------------------------------------//
