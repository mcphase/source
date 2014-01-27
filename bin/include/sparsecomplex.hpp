/* sparsecomplex.hpp
 *
 * This file implements a sparse complex matrix class as a C++ template.
 *
 * This file is part of the McPhase package.
 *
 * (c) 2013 Duc Le - mducle@snu.ac.kr
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 *
 */

#ifndef ZMATHS_H
#define ZMATHS_H

#include<cstdlib>
#include<cmath>
#include<vector>
#include<iostream>
#include<sstream>
#include<iomanip>
#include<string>
#include<map>
#include<cfloat>       // For definition of EPSILON etc.
#include "vector.h"

// --------------------------------------------------------------------------------------------------------------- //
// Template Class to hold a complex sparse matrix of any type - and also declares a few operations and methods
//    The class stores elements by coordinate indices in a subclass, _ind, consisting of two ints, which serves
//    as the key to a C++ STL map.
//    The member functions and overloaded operators has been made to be as similar to Matlab syntax as possible
// --------------------------------------------------------------------------------------------------------------- //
template <class T> class zsMat;                                         // google  "Explicit Template Specification"
template <class T> std::ostream & operator << (std::ostream & o, const zsMat<T> & m);
template <class T> std::istream & operator >> (std::istream & i, zsMat<T> & m);
template <class T> T atoT(const std::string &s);
// --------------------------------------------------------------------------------------------------------------- //
template <class T> class zsMat {
   private:
     int _r, _c;                                                        // number of rows and columns
     class _ind {                                                       // Defines a class for indexing the map
        public:                                                         //   in which we store the matrix elements.
           int r;                                                       //   The map requires a sort operator which
           int c;                                                       //   we provide by overloading '<'
           _ind(int r_, int c_) : r(r_), c(c_) {};
           ~_ind() {};
           bool operator < (const _ind & i) const { return (c<i.c || (c==i.c && r<i.r)); }    // Column-major
           bool operator <= (const _ind & i) const { return (c<i.c || (c==i.c && r<=i.r)); }
     };
     std::map<_ind,std::complex<T> > _ls;                              // List of elements

   public:
     // Constructors and Destructors //
     zsMat(size_t r=1, size_t c=1): _r(r), _c(c) {};                    // constructs an empty r x c sparse matrix
     ~zsMat() {};
    
     // Other member functions //
     int nnz() { return (int) _ls.size(); };                            // Returns the number of nonzero elements
     int nr() const { return _r; };
     int nc() const { return _c; };
     std::vector<int> size() const;                                     // Returns the size of the matrix: _r and _c
     std::vector< std::vector<int> > find() const;                      // Returns the indices of the nonzero elements
     std::vector< std::vector<int> > findupper() const;                 // Returns the indices of upper triangle elements
     std::vector< std::vector<int> > findlower() const;                 // Returns the indices of lower triangle elements
     std::vector<int> find_col(int i) const;                            // Returns the non-zero indices of column i
     bool issymm() const;                                               // Determines if matrix is symmetric
     bool isherm() const;                                               // Determines if matrix is Hermitian
     bool issquare() const { return _r==_c; };                          // Determines if matrix is square
     bool isempty() const { return _ls.empty(); }                       // Determines if matrix is empty
     zsMat<T> transpose();                                              // Transposes the matrix (without complex conjugation)
     zsMat<T> hermitian();                                              // Hermitian conjugate the matrix
     void del(int r, int c) { _ls.erase(_ind(r,c)); };                  // Deletes an element
     void clear() { _ls.clear(); _r = 0; _c = 0; };                     // Clears the matrix
     void zero() { _ls.clear(); };                                      // Zeros all entries
     void zero(int r, int c) { _ls.clear(); _r = r; _c = c; };          // Zeros all entries, changes shape as well
     void resize(int rs, int rn, int cs, int cn);                       // Resizes a matrix and deletes extra elements
     void reshape(int r, int c) { _r = r; _c = c; };                    // Changes shape without changing entries!
     std::string display_full() const;                                  // Outputs a string of the full matrix
     std::complex<T> *f_array() const;                                  // Returns matrix as a Fortran style 2D array
     std::complex<T> *f_array_tr() const;                               // Returns matrix as a Fortran style 2D array
     void h_array(std::complex<T>* retval) const;                       // Returns matrix as a Fortran style 2D array (pre-allocated)
     T* cp_array() const;                                               // Returns Hermitian matrix as a C-style packed 2D array
     T* fp_array() const;                                               // Returns Hermitian matrix as a Fortran packed 2D array
     Matrix fp_matrix() const;                                          // Returns Hermitian matrix as a Fortran packed 2D array
     void packed2tril();                                                // Rearranges a Hermitian matrix from packed to lower triangular
     void mcol(int c, T val);                                           // Multiplies a given column by a constant
     void mrow(int r, T val);                                           // Multiplies a given row by a constant
     void mcol(int c, std::complex<T> val);                             // Multiplies a given column by a constant
     void mrow(int r, std::complex<T> val);                             // Multiplies a given row by a constant

     // For iterative eigensolvers
     void MultMv(std::complex<T> *v, std::complex<T> *w);               // Calculates the matrix-vector product w = M*v (assume Hermitian)
     void MultMvNH(std::complex<T> *v, std::complex<T> *w);             // Calculates the matrix-vector product w = M*v (non-Hermitian)

     // Overloaded operators
     zsMat<T> operator =  (const zsMat & m);                            // Copy assignment - overwrites previous matrix
     zsMat<T> operator =  (const T val);                                // Constant assignment (only diagonal elements)
     zsMat<T> operator =  (const std::complex<T> val);                  // Constant assignment (only diagonal elements)
     zsMat<T> operator += (const zsMat & m);                            // Add another matrix to current (element-wise)
     zsMat<T> operator += (const T val);                                // Add a constant value to current (element-wise)
     zsMat<T> operator += (const std::complex<T> val);                  // Add a constant value to current (element-wise)
     zsMat<T> operator -= (const zsMat & m);                            // Subtract another matrix from current
     zsMat<T> operator -= (const T val);                                // Subtract a constant value from current 
     zsMat<T> operator -= (const std::complex<T> val);                  // Subtract a constant value from current 
     zsMat<T> operator *= (const T & c);                                // Matrix scalar multiplication
     zsMat<T> operator *= (const std::complex<T> & c);                  // Matrix scalar multiplication
     zsMat<T> operator /= (const T & c);                                // Matrix scalar division
     zsMat<T> operator /= (const std::complex<T> & c);                  // Matrix scalar division
     zsMat<T> operator *= (const zsMat & m);                            // Matrix multiplication
     zsMat<T> operator *= (const std::vector<T> & v);                   // Matrix*Vector (Matrix operating on a vector)
     zsMat<T> operator *= (const std::vector<std::complex<T> > & v);    // Matrix*Vector (Matrix operating on a vector)
     zsMat<T> operator ^= (const zsMat & m);                            // Element-wise matrix multiplication
     zsMat<T> operator %= (const zsMat & m);                            // Element-wise matrix division
     std::complex<T>& operator () (int r, int c);                       // Gets/sets element by subscript index
     std::complex<T>  operator () (int r, int c) const;
     std::complex<T>  operator [] (std::pair<int,int>);                 // Read-only subscript operator

     // Friend function  to provide input/output via <iostream>
     friend std::ostream & operator << <> (std::ostream & o, const zsMat & m);
     friend std::istream & operator >> <> (std::istream & i, zsMat & m);

};  // End of template <class T> class zsMat

// --------------------------------------------------------------------------------------------------------------- //
// Member functions
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::vector<int> zsMat<T>::size() const
{
   std::vector<int> retval;
   retval.push_back(_r);
   retval.push_back(_c);
   return retval;
}

template <class T> std::vector< std::vector<int> > zsMat<T>::find() const
{
   int r = 0;
   typename std::map<_ind,std::complex<T> >::iterator i;                // Help out compiler by telling it should by
   std::vector<int> row(2);                                             //   looking for a type, otherwise get error
   std::vector< std::vector<int> > retval(_ls.size(),row);              //   on g++.
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   for (i=tmp_ls.begin(); i!=tmp_ls.end(); i++)
   {
      row[0] = i->first.r; row[1] = i->first.c;
      retval[r] = row; r++;
   }

   return retval;
}

template <class T> std::vector< std::vector<int> > zsMat<T>::findupper() const
{
   int c,n = 0;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::vector<int> row(2);
   std::vector< std::vector<int> > retval(_ls.size()/2+_r,row);
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   for (c=1; c<=_c; c++)
      for (i=tmp_ls.lower_bound(_ind(1,c)); i!=tmp_ls.lower_bound(_ind(c,c)); i++)
      {
         row[0] = i->first.r; row[1] = i->first.c;
         retval[n++] = row;
      }

   retval.erase(retval.begin()+n,retval.end());
   return retval;
}

template <class T> std::vector< std::vector<int> > zsMat<T>::findlower() const
{
   int c,n = 0;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::vector<int> row(2);
   std::vector< std::vector<int> > retval(_ls.size()/2+_r,row);
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   for (c=1; c<=_c; c++)
      for (i=tmp_ls.lower_bound(_ind(c,c)); i!=tmp_ls.lower_bound(_ind(_r,c)); i++)
      {
         row[0] = i->first.r; row[1] = i->first.c;
         retval[n++] = row;
      }

   retval.erase(retval.begin()+n,retval.end());
   return retval;
}


template <class T> std::vector<int> zsMat<T>::find_col(int c) const
{
   int n = 0;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::vector<int> retval(_c,0);
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   for (i=tmp_ls.lower_bound(_ind(0,c)); i!=tmp_ls.lower_bound(_ind(_r,c)); i++)
      retval[n++] = i->first.c;

   retval.erase(retval.begin()+n,retval.end());
   return retval;
}

template <class T> bool zsMat<T>::issymm() const
{
   typename std::map<_ind,std::complex<T> >::iterator it;
   int i,j;
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   if(_r!=_c)
      return false;

   for (it=tmp_ls.begin(); it!=tmp_ls.end(); it++)
   {
      i = (*it).first.c; j = (*it).first.r;
      if(i>=j)                                                          // Look at upper triangle only, to compare
      {                                                                 //    to lower triangle.
         if(tmp_ls.find(_ind(i,j))->second != (*it).second)
            return false;
      }
   }

   return true;
}

template <class T> bool zsMat<T>::isherm() const
{
   typename std::map<_ind,std::complex<T> >::iterator it;
   int i,j;
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   if(_r!=_c)
      return false;

   for (it=tmp_ls.begin(); it!=tmp_ls.end(); it++)
   {
      i = (*it).first.c; j = (*it).first.r;
      if(i>=j)                                                          // Look at upper triangle only, to compare
      {                                                                 //    to lower triangle.
         if(tmp_ls.find(_ind(i,j))->second != conj((*it).second))
            return false;
      }
   }

   return true;
}

template <class T> zsMat<T> zsMat<T>::transpose()
{
   typename std::map<_ind,std::complex<T> >::iterator it;
   int i,j;
   std::map<_ind,std::complex<T> > tmp_ls;

   for(it=_ls.begin(); it!=_ls.end(); it++)
   {
      i = (*it).first.r; j = (*it).first.c;
      tmp_ls[_ind(j,i)] = (*it).second;
   }
   _ls = tmp_ls;
   i = _c; _c = _r; _r = i;
   return *this;
}

template <class T> zsMat<T> zsMat<T>::hermitian()
{
   typename std::map<_ind,std::complex<T> >::iterator it;
   int i,j;
   std::map<_ind,std::complex<T> > tmp_ls;

   for(it=_ls.begin(); it!=_ls.end(); it++)
   {
      i = (*it).first.r; j = (*it).first.c;
      tmp_ls[_ind(j,i)] = conj((*it).second);
   }
   _ls = tmp_ls;
   i = _c; _c = _r; _r = i;
   return *this;
}

template <class T> void zsMat<T>::resize(int rs, int rn, int cs, int cn)// Resizes a matrix and deletes extra elements
{
   int i;
   typename std::map<_ind,std::complex<T> >::iterator it;

   if(rs==-1) rs = 1; if(rn==-1) rn = _r;
   if(cs==-1) cs = 1; if(cn==-1) cn = _r;
   
 //rs--; cs--; //rn--; cn--;
   _ls.erase(_ls.upper_bound(_ind(rn,cn)),_ls.end());
   if(rs==1 && cs==1)
   {
      for(i=1; i<=cn; i++)
         _ls.erase(_ls.upper_bound(_ind(rn,i)),_ls.lower_bound(_ind(0,i+1)));
      _r = rn; _c = cn;
   }
   else
   {
      _ls.erase(_ls.begin(),_ls.lower_bound(_ind(rs,cs)));
      for(i=cs; i<cn; i++)
      {  
         _ls.erase(_ls.lower_bound(_ind(0,i)),_ls.lower_bound(_ind(rs,i)));
         _ls.erase(_ls.upper_bound(_ind(rn,i)),_ls.lower_bound(_ind(0,i+1)));
         for(it=_ls.lower_bound(_ind(rs,i)); it!=_ls.upper_bound(_ind(rn,i)); i++)
         {
            _ls[_ind(i-rs,it->first.c-cs)] = it->second; _ls.erase(it);
         }
      }
      _r = (rn-rs)+1; _c = (cn-cs)+1;
   }
}

template <class T> std::string zsMat<T>::display_full() const           // Prints the full matrix to a string.
{
   int r,c,mr,mc;
   std::stringstream retval;

   retval << "[";
   for (r=1; r<=(_r-1); r++)
   {
      for (c=1; c<=(_c-1); c++)
      {
         if(_ls.find(_ind(r,c))==_ls.end())
            retval << "0\t";
         else
         {
            mr = _ls.find(_ind(r,c))->first.r;
            mc = _ls.find(_ind(r,c))->first.c;
            if( (mr>_r+1 || mr<0) || (mc>_c+1 || mc<0) )
               retval << "0\t";
            else
               retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << "\t";
         }
      }
      if(_ls.find(_ind(r,c))==_ls.end())
         retval << "0;\n";
      else
      {
         mr = _ls.find(_ind(r,c))->first.r;
         mc = _ls.find(_ind(r,c))->first.c;
         if( (mr>_r+1 || mr<0) || (mc>_c+1 || mc<0) )
            retval << "0;\n";
         else
            retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << ";\n";
      }
   }
   for (c=1; c<=(_c-1); c++)
   {
      if(_ls.find(_ind(r,c))==_ls.end())
         retval << "0\t";
      else
      {
         mr = _ls.find(_ind(r,c))->first.r;
         mc = _ls.find(_ind(r,c))->first.c;
         if( (mr>_r+1 || mr<0) || (mc>_c+1 || mc<0) )
            retval << "0\t";
         else
            retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << "\t";
      }
   }
   if(_ls.find(_ind(r,c))==_ls.end())
      retval << "0];\n";
   else
   {
      mr = _ls.find(_ind(r,c))->first.r;
      mc = _ls.find(_ind(r,c))->first.c;
      if( (mr>_r+1 || mr<0) || (mc>_c+1 || mc<0) )
         retval << "0];\n";
      else
         retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << "];\n";
   }
   return retval.str();
}

template <class T> std::complex<T>* zsMat<T>::f_array() const   // Returns matrix as a Fortran style 2D array
{
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   std::complex<T> *retval;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   // Allocates an _r*_c array and initiallises all elements to zero.
   retval = (std::complex<T>*) calloc(_r*_c,sizeof(std::complex<T>));

   for (i=tmp_ls.begin(); i!=tmp_ls.end(); i++) 
      retval[_r*(i->first.c-1)+(i->first.r-1)] = i->second; 

   return retval;
}
template <class T> void zsMat<T>::h_array(std::complex<T>* retval) const   // Assumes matrix already allocated
{
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::map<_ind,std::complex<T> > tmp_ls = _ls;
   memset(retval,0,_r*_c*sizeof(std::complex<T>));
   // _Assume_(!) the matrix is Hermitian and only loop through the lower triangle.
   for (int c=1; c<=_c; c++)
      for (i=tmp_ls.lower_bound(_ind(c,c)); i!=tmp_ls.lower_bound(_ind(_r+1,c)); i++)
      {
         retval[_r*(i->first.r-1)+(i->first.c-1)] = i->second;
         if(i->first.c!=i->first.r)
            retval[_c*(i->first.c-1)+(i->first.r-1)] = conj(i->second);
      }
}
template <class T> std::complex<T>* zsMat<T>::f_array_tr() const// Returns the transposed matrix as a Fortran style 2D array
{
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   std::complex<T> *retval;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   // Allocates an _r*_c array and initiallises all elements to zero.
   retval = (std::complex<T>*) calloc(_r*_c,sizeof(std::complex<T>));

   for (i=tmp_ls.begin(); i!=tmp_ls.end(); i++) 
      retval[_c*(i->first.r-1)+(i->first.c-1)] = i->second; 

   return retval;
}
template <class T> T* zsMat<T>::fp_array() const                // Returns Hermitian matrix as a Fortran style packed 2D array
{                                                               // with real in upper, imag in lower triangles.
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   T *retval;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   retval = (T*) calloc(_r*_c,sizeof(T));                       // Allocates an _r*_c array and initiallises all elements to zero.

   // _Assume_(!) the matrix is Hermitian and only loop through the lower triangle.
   for (int c=1; c<=_c; c++)
      for (i=tmp_ls.lower_bound(_ind(c,c)); i!=tmp_ls.lower_bound(_ind(_r+1,c)); i++)
      {
         retval[_r*(i->first.r-1)+(i->first.c-1)] = real(i->second);
         if(i->first.c!=i->first.r)
            retval[_c*(i->first.c-1)+(i->first.r-1)] = imag(i->second);
      }

   return retval;
}
template <class T> Matrix zsMat<T>::fp_matrix() const           // Returns Hermitian matrix as a Fortran style packed 2D array
{                                                               // with real in upper, imag in lower triangles.
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   Matrix retval(1,_r,1,_c); retval=0;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   // _Assume_(!) the matrix is Hermitian and only loop through the lower triangle.
   for (int c=1; c<=_c; c++)
      for (i=tmp_ls.lower_bound(_ind(c,c)); i!=tmp_ls.lower_bound(_ind(_r+1,c)); i++)
      {
         retval(i->first.r,i->first.c) = real(i->second);
         if(i->first.c!=i->first.r)
            retval(i->first.c,i->first.r) = imag(i->second);
      }

   return retval;
}
template <class T> T* zsMat<T>::cp_array() const                // Returns Hermitian matrix as a C style packed 2D array
{                                                               // with real in upper, imag in lower triangles.
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   T *retval;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   retval = (T*) calloc(_r*_c,sizeof(T));                       // Allocates an _r*_c array and initiallises all elements to zero.

   // _Assume_(!) the matrix is Hermitian and only loop through the lower triangle.
   for (int c=1; c<=_c; c++)
      for (i=tmp_ls.lower_bound(_ind(c,c)); i!=tmp_ls.lower_bound(_ind(_r+1,c)); i++)
      {
         retval[_r*(i->first.c-1)+(i->first.r-1)] = real(i->second);
         if(i->first.c!=i->first.r)
            retval[_r*(i->first.r-1)+(i->first.c-1)] = imag(i->second);
      }

   return retval;
}

template <class T> void zsMat<T>::packed2tril()                 // Rearranges a Hermitian matrix from packed to lower triangular
{
   int r,c;
   typename std::map<_ind,std::complex<T> >::iterator i;

   for (i=_ls.begin(); i!=_ls.end(); i++)
   {
      if(imag(i->second)!=0) {
         if(real(i->second)==0) _ls.erase(_ind(i->first.r,i->first.c));
         else i->second = complex<T>(real(i->second),0); }
   }

   for (i=_ls.begin(); i!=_ls.end(); i++)
   {
      r = i->first.r; c = i->first.c;
      if(c>r) 
      {
         _ls[_ind(c,r)] += complex<T>(0,real(i->second));
         _ls.erase(_ind(r,c)); 
      }
   }
}

template <class T> void zsMat<T>::mcol(int c, T val)            // Multiplies a particular column by a constant
{
   typename std::map<_ind,std::complex<T> >::iterator it;
   for (it=_ls.lower_bound(_ind(1,c)); it!=_ls.lower_bound(_ind(0,c+1)); it++)
      it->second *= val;
}
template <class T> void zsMat<T>::mrow(int r, T val)            // Multiplies a particular column by a constant
{
   int c;
   typename std::map<_ind,std::complex<T> >::iterator it;
   for(c=1; c<=_r; c++)
   {
      it = _ls.find(_ind(r,c));
      if(it!=_ls.end())
         it->second *= val;
   }
}
template <class T> void zsMat<T>::mcol(int c, std::complex<T> val)
{
   typename std::map<_ind,std::complex<T> >::iterator it;
   for (it=_ls.lower_bound(_ind(1,c)); it!=_ls.lower_bound(_ind(0,c+1)); it++)
      it->second *= val;
}
template <class T> void zsMat<T>::mrow(int r, std::complex<T> val) 
{
   int c;
   typename std::map<_ind,std::complex<T> >::iterator it;
   for(c=1; c<=_r; c++) {
      it = _ls.find(_ind(r,c));
      if(it!=_ls.end())
         it->second *= val; }
}

// --------------------------------------------------------------------------------------------------------------- //
// Matrix-Vector multiplication for iterative methods
// --------------------------------------------------------------------------------------------------------------- //
template <class T> void zsMat<T>::MultMv(std::complex<T> *v, std::complex<T> *w)  // Calculates the matrix-vector product w = M*v
{                                                                                 // Needed by ARPACK
   typename std::map<_ind,std::complex<T> >::iterator i;
   // We have to assume that the size of the vector v is equal to _r
   memset(w,0,_r*sizeof(std::complex<T>));
   for (int c=1; c<=_c; c++)
   {
      // Assume(!) matrix is Hermitian and loops only over the lower triangle
      for (i=_ls.lower_bound(_ind(c,c)); i!=_ls.lower_bound(_ind(_r+1,c)); i++)
      {
         w[i->first.r-1] += (i->second * v[c-1]);
         if(i->first.c!=i->first.r)
            w[c-1] += conj(i->second)*v[i->first.r-1];
      }
   }

}
template <class T> void zsMat<T>::MultMvNH(std::complex<T>*v, std::complex<T>*w)  // Calculates the matrix-vector product w = M*v
{                                                                                 // Needed by ARPACK
   typename std::map<_ind,std::complex<T> >::iterator i;
   // We have to assume that the size of the vector v is equal to _r
   memset(w,0,_r*sizeof(std::complex<T>));
   for (int c=1; c<=_c; c++)
   {
      for (i=_ls.lower_bound(_ind(1,c)); i!=_ls.lower_bound(_ind(0,c+1)); i++)
         w[i->first.r-1] += (i->second * v[c-1]);
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Overloaded operators
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::complex<T>& zsMat<T>::operator () (int r, int c)
{ 
   if(r>_r) _r = r; if(c>_c) _c = c;
   return _ls[_ind(r,c)]; 
}
template <class T> std::complex<T> zsMat<T>::operator () (int r, int c) const
{
   int mr,mc;
   if(_ls.find(_ind(r,c))==_ls.end())
      return std::complex<T> (0,0);
   else
   {
      mr = _ls.find(_ind(r,c))->first.r;
      mc = _ls.find(_ind(r,c))->first.c;
      if( (mr>_r+1 || mr<1) || (mc>_c+1 || mc<1) ) 
         return std::complex<T> (0,0);
      else 
         return _ls.find(_ind(r,c))->second;
   }
}
template <class T> std::complex<T> zsMat<T>::operator [] (std::pair<int,int> rc)
{
   int mr,mc, r=rc.first, c=rc.second;
   if(_ls.find(_ind(r,c))==_ls.end())
      return std::complex<T> (0,0);
   else
   {
      mr = _ls.find(_ind(r,c))->first.r;
      mc = _ls.find(_ind(r,c))->first.c;
      if( (mr>_r+1 || mr<1) || (mc>_c+1 || mc<1) ) 
         return std::complex<T> (0,0);
      else 
         return _ls.find(_ind(r,c))->second;
   }
}

template <class T> zsMat<T> zsMat<T>::operator = (const zsMat<T> & m)
{
   _r = m.nr(); _c = m.nc();
   _ls = m._ls;
   return *this;
}
template <class T> zsMat<T> zsMat<T>::operator = (const T val)          // Assignment of value to diagonal
{
   if(val==0) { _ls.clear(); } else {
   int i_max = (_r>_c)?_c:_r;
   for(int i=1; i<=i_max; i++) {
      _ls[_ind(i,i)] = std::complex<T>(val,0); } }
   return *this;
}
template <class T> zsMat<T> zsMat<T>::operator = (const std::complex<T> val) 
{
   for(int i=1; i<(_r>_c)?_c:_r; i++)
      _ls[_ind(i,i)] = val;
   return *this;
}

template <class T> zsMat<T> zsMat<T>::operator += (const zsMat<T> & m)  // Adds another matrix to current
{
   typename std::map<_ind,std::complex<T> >::const_iterator it;
   if(!m._ls.empty())
      for(it=m._ls.begin(); it!=m._ls.end(); it++)
         _ls[_ind(it->first.r,it->first.c)] += it->second;

   return *this;
}
template <class T> zsMat<T> operator + (const zsMat<T> & m1, const zsMat<T> & m2)
{
   zsMat<T> tmp = m1; tmp += m2; return tmp;                            // Binary addition of two sparse matrices
}

template <class T> zsMat<T> zsMat<T>::operator += (const T val)         // Adds a constant to the current matrix
{
/*
   typename std::map<_ind,std::complex<T> >::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second += val;
*/
   int i_max = (_r>_c)?_c:_r;
   for(int i=1; i<=i_max; i++) {
      _ls[_ind(i,i)] += std::complex<T>(val,0); } 
   return *this;
}
template <class T> zsMat<T> zsMat<T>::operator += (const std::complex<T> val)
{
/*
   typename std::map<_ind,std::complex<T> >::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second += val;
*/
   int i_max = (_r>_c)?_c:_r;
   for(int i=1; i<=i_max; i++) {
      _ls[_ind(i,i)] += val; } 
   return *this;
}
template <class T> zsMat<T> operator + (const zsMat<T> & m, const T val) {
   zsMat<T> tmp = m; tmp += val; return tmp; }
template <class T> zsMat<T> operator + (const T val, const zsMat<T> & m) {
   zsMat<T> tmp = m; tmp += val; return tmp; }
template <class T> zsMat<T> operator + (const zsMat<T> & m, const std::complex<T> val) {
   zsMat<T> tmp = m; tmp += val; return tmp; }
template <class T> zsMat<T> operator + (const std::complex<T> val, const zsMat<T> & m) {
   zsMat<T> tmp = m; tmp += val; return tmp; }

template <class T> zsMat<T> zsMat<T>::operator -= (const zsMat<T> & m)  // Subtracts another matrix from current
{
   typename std::map<_ind,std::complex<T> >::const_iterator it;
   if(!m._ls.empty())
      for(it=m._ls.begin(); it!=m._ls.end(); it++)
      {
         _ls[_ind(it->first.r,it->first.c)] -= it->second;
         if ((_ls.find(_ind(it->first.r,it->first.c))->second) == 0.0)
            _ls.erase(_ind(it->first.r,it->first.c));
      }

   return *this;
}
template <class T> zsMat<T> operator - (const zsMat<T> & m1, const zsMat<T> & m2)
{
   zsMat<T> tmp = m1; tmp -= m2; return tmp;                            // Binary subtraction of two sparse matrices
}

template <class T> zsMat<T> zsMat<T>::operator -= (const T val)         // Subtract a constant from current matrix
{
/*
   typename std::map<_ind,std::complex<T> >::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second -= val;
*/
   int i_max = (_r>_c)?_c:_r;
   for(int i=1; i<=i_max; i++) {
      _ls[_ind(i,i)] -= std::complex<T>(val,0); } 
   return *this;
}
template <class T> zsMat<T> zsMat<T>::operator -= (const std::complex<T> val)
{
/*
   typename std::map<_ind,std::complex<T> >::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second -= val;
*/
   int i_max = (_r>_c)?_c:_r;
   for(int i=1; i<=i_max; i++) {
      _ls[_ind(i,i)] -= val; } 
   return *this;
}
template <class T> zsMat<T> operator - (const zsMat<T> & m, const T val) {
   zsMat<T> tmp = m; tmp -= val; return tmp; }
template <class T> zsMat<T> operator - (const T val, const zsMat<T> & m) {
   zsMat<T> tmp = m; tmp -= val; return tmp; }

template <class T> zsMat<T> zsMat<T>::operator *= (const T & c)         // Matrix scalar multiplication
{
   typename std::map<_ind,std::complex<T> >::iterator i;

   for (i=_ls.begin(); i!=_ls.end(); i++)
   {
      i->second *= c;
   }
   return *this;
}
template <class T> zsMat<T> zsMat<T>::operator *= (const std::complex<T> & c)
{
   typename std::map<_ind,std::complex<T> >::iterator i;
   for (i=_ls.begin(); i!=_ls.end(); i++)
      i->second *= c;
   return *this;
}
template <class T> zsMat<T> operator * (const zsMat<T> & m1, const T & c) {
   zsMat<T> tmp = m1; tmp *= c; return tmp; }                           // Binary matrix scalar multiplication
template <class T> zsMat<T> operator * (const T & c, const zsMat<T> & m1) {
   zsMat<T> tmp = m1; tmp *= c; return tmp; }

template <class T> zsMat<T> zsMat<T>::operator /= (const T & c)         // Matrix scalar division
{
   typename std::map<_ind,std::complex<T> >::iterator i;

   for (i=_ls.begin(); i!=_ls.end(); i++)
   {
      i->second /= c;
   }
   return *this;
}
template <class T> zsMat<T> operator / (const zsMat<T> & m1, const T & c)
{
   zsMat<T> tmp = m1; tmp /= c; return tmp;                             // Binary matrix scalar division
}

template <class T> zsMat<T> zsMat<T>::operator *= (const zsMat<T> & m)  // Matrix multiplication
{
   zsMat<T> tmp;

   if(_c != m._r)
      std::cerr << "Error: zsMat<T>::operator *= (const zsMat<T>): Matrix sizes are incommensurate!\n";
   else
   {
      // Using the traditional method, except ignoring zeros in the sparse structure
      tmp._r = _r; 
      tmp._c = m._c;
      typename std::map<_ind,std::complex<T> >::iterator it;
      std::map<_ind,std::complex<T> > mls = m._ls;
      std::complex<T> elem;
      
      for (int c=1; c<=_c; c++)       // Loops through the columns of the multiplier
         for(int r=1; r<=_r; r++)     // Loops through the rows of the multiplican
         {
            elem = 0.;                // Loops through the rows of the mulitplier, ignoring zeros
            for (it=mls.lower_bound(_ind(c,c)); it!=mls.lower_bound(_ind(_r+1,c)); it++)
            {
                  elem += (_ls.find(_ind(r,it->first.r))->second * it->second);
            }
            if(elem!=0.) tmp(r,c) = elem;
         }
      *this = tmp;
   }
   return *this;
}
template <class T> zsMat<T> operator * (const zsMat<T> & m1, const zsMat<T> & m2)
{
   zsMat<T> tmp = m1;                                                   // Binary matrix multiplication
   tmp *= m2;
   return tmp;
}
template <class T> zsMat<T> zsMat<T>::operator *= (const std::vector<T> & v)
{
   zsMat<T> tmp;                                                        // Binary matrix.vector multiplication
   typename std::map<_ind,std::complex<T> >::iterator i;

   if(_c != (int)v.size())
      std::cerr << "Error: zsMat<T>::operator * (zsMat<T>, std::vector<T>): Matrix and vector sizes incommensurate!\n";
   else
   {
      tmp._r = _r; tmp._c = 1;
      for (int c=1; c<=_c; c++)
      {
         for (i=_ls.lower_bound(_ind(1,c)); i!=_ls.lower_bound(_ind(0,c+1)); i++)
            tmp(c,1) += (i->second * v[c-1]);
      }
      *this = tmp;
   }

   return *this;
}
template <class T> std::vector<T> operator * (const zsMat<T> & m1, const std::vector<T> & v)
{
   std::vector<T> tmp;                                                  // Binary matrix.vector multiplication
   zsMat<T> tmpM = m1;
   int i;

   tmpM *= v; 
   for(i=0; i<tmpM.nr(); i++)
      tmp.push_back(tmpM(i,0));

   return tmp;
}
template <class T> zsMat<T> zsMat<T>::operator *= (const std::vector<std::complex<T> > & v)
{
   zsMat<T> tmp;                                                        // Binary matrix.vector multiplication
   typename std::map<_ind,std::complex<T> >::iterator i;

   if(_c != (int)v.size())
      std::cerr << "Error: zsMat<T>::operator * (zsMat<T>, std::vector<T>): Matrix and vector sizes incommensurate!\n";
   else
   {
      tmp._r = _r; tmp._c = 1;
      for (int c=1; c<=_c; c++)
      {
         for (i=_ls.lower_bound(_ind(1,c)); i!=_ls.lower_bound(_ind(0,c+1)); i++)
            tmp(c,1) += (i->second * v[c-1]);
      }
      *this = tmp;
   }

   return *this;
}
template <class T> std::vector<T> operator * (const zsMat<T> & m1, const std::vector<std::complex<T> > & v)
{
   std::vector<T> tmp;                                                  // Binary matrix.vector multiplication
   zsMat<T> tmpM = m1;
   int i;

   tmpM *= v; 
   for(i=0; i<tmpM.nr(); i++)
      tmp.push_back(tmpM(i,0));

   return tmp;
}

template <class T> zsMat<T> zsMat<T>::operator ^= (const zsMat<T> &m)   // Element-wise matrix multiplication 
{
   if(m.nr()!=_r || m.nc()!=_c) { std::cerr << "zsMat ^ operator: M1 and M2 not same size\n"; return *this; }
   typename std::map<_ind,std::complex<T> >::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second *= m(it->first.r,it->first.c);
   return *this;
}
template <class T> zsMat<T> operator ^ (const zsMat<T> & m1, const zsMat<T> & m2)
{
   zsMat<T> tmp = m1;
   tmp ^= m2;
   return tmp;
}

template <class T> zsMat<T> zsMat<T>::operator %= (const zsMat<T> &m)    // Element-wise matrix division
{
   if(m.nr()!=_r || m.nc()!=_c) { std::cerr << "zsMat % operator: M1 and M2 not same size\n"; return *this; }
   typename std::map<_ind,std::complex<T> >::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second /= m(it->first.r,it->first.c);
   return *this;
}
template <class T> zsMat<T> operator % (const zsMat<T> & m1, const zsMat<T> & m2)
{
   zsMat<T> tmp = m1;
   tmp %= m2;
   return tmp;
}

// --------------------------------------------------------------------------------------------------------------- //
// Friend functions
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::ostream & operator << (std::ostream & o, const zsMat<T> & m)
{
   int i,r,c;
   std::vector< std::vector<int> > nz = m.find();
   int sz = (int)nz.size();
   if(sz != 0)
   {
      for (i=0; i<sz; i++)
      {
         r = nz[i][0]; c = nz[i][1];
         o << "x(" << (r+1) << "," << (c+1) << ")\t=\t" << m(r,c) << ";\n";
      }
   }
   else 
   {
      o << "empty matrix\n";
   }
   
   return o;
}

template <class T> std::istream & operator >> (std::istream & i, zsMat<T> & m)
{
   int r,c;
   T x,y;
   for (r=1; r<=m._r; r++)
      for (c=1; c<=m._c; c++)
      {
         i >> x; i >> y;
         if (x!=0.0)
            m(r,c) = std::complex<T> (x,y);
      }

   return i;
}

// --------------------------------------------------------------------------------------------------------------- //
// Declaration for functions in feast.cpp
// --------------------------------------------------------------------------------------------------------------- //
void feast_zheev(zsMat<double> &H, Vector &en, Matrix &zr, Matrix &zc);

#endif
