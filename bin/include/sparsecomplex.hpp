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
           bool operator < (const _ind & i) const { return (r<i.r || (r==i.r && c<i.c)); }
           bool operator <= (const _ind & i) const { return (r<i.r || (r==i.r && c<=i.c)); }
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
   //std::map<_ind,std::complex<T> > ls() { return _ls; }
     void mcol(int c, T val);                                           // Multiplies a given column by a constant
     void mrow(int r, T val);                                           // Multiplies a given row by a constant
     void mcol(int c, std::complex<T> val);                             // Multiplies a given column by a constant
     void mrow(int r, std::complex<T> val);                             // Multiplies a given row by a constant

     // For TRLAN
     void trlanmv(int*nrow, int*ncol, double*xin, int*ldx, double*yout, int*ldy);

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

/* --------------------------------------------------------------------------------------------------------------- //
   Template functions declared in this file
// --------------------------------------------------------------------------------------------------------------- //
   enum normp { inf = 3, fro = 4 };                                     // Enumeration for the norm function
  
   T atoT(const std::string &s)                                         // Convert a string to a T using istreams
   std::string dispvect(const std::vector<T> &V)                        // Displays a std::vector in row format
   T vmax(const std::vector<T> &V)                                      // Finds the maximum value of a vector
   std::vector<T> mmax(const sMat<T> &M)                                // Finds the maximum values along a column
   std::vector<T> mmax(const sMat<T> &M, std::vector<int> &ind)         //    with indeces
   std::vector<T> msum(const sMat<T> &M)                                // Calcs. sum of each column of a matrix
   T vsum(const std::vector<T> &V)                                      // Calcs. sum of all elements of a vector
   std::vector<T> diag(const sMat<T> &M, int d=0)                       // Returns the diagonal elements as vector
   sMat<T> triu(const sMat<T> &M)                                       // Returns the upper triangle   
   sMat<T> triu(const sMat<T> &M, int d)                                //    from diagonal d
   sMat<T> tril(const sMat<T> &M)                                       // Returns the lower triangle 
   sMat<T> tril(const sMat<T> &M, int d)                                //    from diagonal d
   std::vector<T> all(const sMat<T> &M, int dim)                        // Returns a vector whose element is 1 if all
   std::vector<T> all(sMat<T> & M)                                      //    elements along dimension dim is nonzero
   std::vector<T> vsort(const std::vector<T> &v)                        // Sorts a vector using the gnome sort
   std::vector<T> vsort(const std::vector<T> &v, std::vector<int> &ind) //    with indices
   sMat<T> msort(const sMat<T> &m)                                      // Sorts a matrix along its column
   sMat<T> msort(const sMat<T> &m, sMat<int> &ind, const char *descend) //    with indices, descending if specified
   std::vector<T> setunion(std::vector<T> &a, std::vector<T> &b)        // Calculates the set-theoretic union
   std::vector<T> setdiff(std::vector<T> &a, std::vector<T> &b)         // Calculates the set-theoretic difference
   std::vector<T> setxor(std::vector<T> &a, std::vector<T> &b)          // Calculates the set-theoretic xor
// --------------------------------------------------------------------------------------------------------------- */

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
   int r,n = 0;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::vector<int> row(2);
   std::vector< std::vector<int> > retval(_ls.size()/2+_r,row);
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   for (r=1; r<=_r; r++)
      for (i=tmp_ls.lower_bound(_ind(r,r)); i!=tmp_ls.lower_bound(_ind(r,_c)); i++)
      {
         row[0] = i->first.r; row[1] = i->first.c;
         retval[n++] = row;
      }

   retval.erase(retval.begin()+n,retval.end());
   return retval;
}

template <class T> std::vector<int> zsMat<T>::find_col(int r) const
{
   int n = 0;
   typename std::map<_ind,std::complex<T> >::iterator i;
   std::vector<int> retval(_c,0);
   std::map<_ind,std::complex<T> > tmp_ls = _ls;

   for (i=tmp_ls.lower_bound(_ind(r,0)); i!=tmp_ls.lower_bound(_ind(r,_c)); i++)
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
      for(i=1; i<=rn; i++)
         _ls.erase(_ls.upper_bound(_ind(i,cn)),_ls.lower_bound(_ind(i+1,0)));
      _r = rn; _c = cn;
   }
   else
   {
      _ls.erase(_ls.begin(),_ls.lower_bound(_ind(rs,cs)));
      for(i=rs; i<rn; i++)
      {  
         _ls.erase(_ls.lower_bound(_ind(i,0)),_ls.lower_bound(_ind(i,cs)));
         _ls.erase(_ls.upper_bound(_ind(i,cn)),_ls.lower_bound(_ind(i+1,0)));
         for(it=_ls.lower_bound(_ind(i,cs)); it!=_ls.upper_bound(_ind(i,cn)); i++)
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
   for (int r=1; r<=_r; r++)
      for (i=tmp_ls.lower_bound(_ind(r,1)); i!=tmp_ls.lower_bound(_ind(r,r+1)); i++)
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
   for (int r=1; r<=_r; r++)
      for (i=tmp_ls.lower_bound(_ind(r,1)); i!=tmp_ls.lower_bound(_ind(r,r+1)); i++)
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
   for (int r=1; r<=_r; r++)
      for (i=tmp_ls.lower_bound(_ind(r,1)); i!=tmp_ls.lower_bound(_ind(r,r+1)); i++)
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
   for (int r=1; r<=_r; r++)
      for (i=tmp_ls.lower_bound(_ind(r,1)); i!=tmp_ls.lower_bound(_ind(r,r+1)); i++)
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
   int r;
   typename std::map<_ind,std::complex<T> >::iterator it;
   for(r=1; r<=_r; r++)
   {
      it = _ls.find(_ind(r,c));
      if(it!=_ls.end())
         it->second *= val;
   }
}
template <class T> void zsMat<T>::mrow(int r, T val)            // Multiplies a particular column by a constant
{
   typename std::map<_ind,std::complex<T> >::iterator lb,ub,it;
   lb = _ls.lower_bound(_ind(r,0));
   ub = _ls.upper_bound(_ind(r,_c));
   if (lb->first < ub->first)
      for(it=lb; it->first<ub->first; it++)
         it->second *= val;
}
template <class T> void zsMat<T>::mcol(int c, std::complex<T> val)
{
   int r;
   typename std::map<_ind,std::complex<T> >::iterator it;
   for(r=1; r<=_r; r++) {
      it = _ls.find(_ind(r,c));
      if(it!=_ls.end())
         it->second *= val; }
}
template <class T> void zsMat<T>::mrow(int r, std::complex<T> val) 
{
   typename std::map<_ind,std::complex<T> >::iterator lb,ub,it;
   lb = _ls.lower_bound(_ind(r,0)); ub = _ls.upper_bound(_ind(r,_c));
   if (lb->first < ub->first)
      for(it=lb; it->first<ub->first; it++)
         it->second *= val;
}

// --------------------------------------------------------------------------------------------------------------- //
// Matrix-Vector multiplication for TRLAN
// --------------------------------------------------------------------------------------------------------------- //
template <class T> void zsMat<T>::trlanmv(int*nrow, int*ncol, double*xin, int*ldx, double*yout, int*ldy)
{
   if(*nrow!=_r) { std::cerr << "trlanmv: TRLAN expects " <<  *nrow << "rows but matrix has " << _r << ".\n"; }
   
   typename std::map<_ind,std::complex<T> >::iterator it;

   for(int j=0; j<*ncol; j++)
   {
      for(int i=1; i<=_c; i++)
      {
         yout[(i-1)+(*ldy)*j] = 0.;
         for(it=_ls.lower_bound(_ind(i,0)); it!=_ls.lower_bound(_ind(i,_c)); it++)
            yout[(i-1)+(*ldy)*j] += real(it->second * xin[(it->first.c-1)+(*ldx)*j]);
      }
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
      int r,c;
      typename std::map<_ind,std::complex<T> >::iterator lb,ub,it;
      std::complex<T> elem;

      // Using the traditional method, except ignoring zeros in the sparse structure
      tmp._r = _r;
      tmp._c = m._c;
      for(r=1; r<=_r; r++)
      {
         lb = _ls.lower_bound(_ind(r,0));
         ub = _ls.upper_bound(_ind(r,_c));
         if(lb->first < ub->first)
            for(c=1; c<=m._c; c++)
            {
               elem = 0.0;
               for(it=lb; it->first<ub->first; it++)
                  elem += (it->second * m._ls.find(_ind(it->first.c,c))->second);
               tmp(r,c) = elem;
            }
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
   int i;
   typename std::map<_ind,std::complex<T> >::iterator lb,ub,it;
   T elem;

   if(_c != (int)v.size())
      std::cerr << "Error: zsMat<T>::operator * (zsMat<T>, std::vector<T>): Matrix and vector sizes incommensurate!\n";
   else
   {
      tmp._r = _r; tmp._c = 1;
      for(i=1; i<=_c; i++)
      {
         lb = _ls.lower_bound(_ind(i,0));
         ub = _ls.upper_bound(_ind(i,_c));
         if (lb->first < ub->first)
         {
            elem = 0.0;
            for(it=lb; it->first<ub->first; it++)
               elem += (it->second * v[it->first.c]);
            tmp(i,0) = elem;
         }
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
   int i;
   typename std::map<_ind,std::complex<T> >::iterator lb,ub,it;
   T elem;

   if(_c != (int)v.size())
      std::cerr << "Error: zsMat<T>::operator * (zsMat<T>, std::vector<T>): Matrix and vector sizes incommensurate!\n";
   else
   {
      tmp._r = _r; tmp._c = 1;
      for(i=1; i<=_c; i++)
      {
         lb = _ls.lower_bound(_ind(i,0));
         ub = _ls.upper_bound(_ind(i,_c));
         if (lb->first < ub->first)
         {
            elem = 0.0;
            for(it=lb; it->first<ub->first; it++)
               elem += (it->second * v[it->first.c]);
            tmp(i,0) = elem;
         }
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

#if 0
// --------------------------------------------------------------------------------------------------------------- //
// Function to convert a string to a T for the sMat<T> template class using input streams.
// --------------------------------------------------------------------------------------------------------------- //
template <class T> T atoT(const std::string &s)
{
   std::istringstream iss(s);
   T t;
   iss >> t;
   return t;
}

// --------------------------------------------------------------------------------------------------------------- //
// Displays a std::vector in row format
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::string dispvect(const std::vector<T> &V)
{
   int i;
   int sz = (int)V.size();
   std::stringstream retval;
   retval << "[ ";
   for (i=0; i<sz; i++)
      retval << V[i] << " ";
   retval << "]";

   return retval.str();
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the maximum value of a vector 
// --------------------------------------------------------------------------------------------------------------- //
template <class T> T vmax(const std::vector<T> &V)
{
   int i;
   T retval = 0;
   int n = V.size();

   for(i=0; i<n; i++)
      if(V[i] > retval) retval = V[i];

   return retval;
}
template <class T> T vmax(const std::vector<T> &V, int &ind)  // Returns the index of the maximum value as well.
{
   int i,n = V.size(); T retval = 0;
   for(i=0; i<n; i++)
      if(V[i] > retval) { retval = V[i]; ind = i; }
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the maximum values of a matrix along a column
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::vector<T> mmax(const sMat<T> &M)
{
   std::vector< std::vector<int> > nz = M.find();
   int i,sz = (int)nz.size();
   std::vector<T> retval(M.nc(),0);

   for (i=0; i<sz; i++)
      if( M(nz[i][0],nz[i][1]) > retval[nz[i][1]] ) retval[nz[i][1]] = M(nz[i][0],nz[i][1]);

   return retval;
}
// Returns the index of the maximum value as well.
template <class T> std::vector<T> mmax(const sMat<T> &M, std::vector<int> &ind)
{
   std::vector< std::vector<int> > nz = M.find(); int i,sz = (int)nz.size();
   std::vector<T> retval(M.nc(),0); ind.clear(); ind.resize(M.nc(),0);
   for (i=0; i<sz; i++)
      if( M(nz[i][0],nz[i][1]) > retval[nz[i][1]] ) { retval[nz[i][1]] = M(nz[i][0],nz[i][1]); ind[nz[i][1]] = nz[i][0]; }
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the sum of each column of a matrix
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::vector<T> msum(const sMat<T> &M)
{
   std::vector< std::vector<int> > nz = M.find();
   int i,sz = (int)nz.size();
   std::vector<T> retval(M.nc(),0);

   for (i=0; i<sz; i++)
      retval[nz[i][1]] += M(nz[i][0],nz[i][1]);

   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the sum of all the elements of a vector
// --------------------------------------------------------------------------------------------------------------- //
template <class T> T vsum(const std::vector<T> &V)
{
   int i;
   T retval = 0;
   int n = V.size();

   for(i=0; i<n; i++)
      retval += V[i];

   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Returns the diagonal elements of a matrix as a vector
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::vector<T> diag(const sMat<T> &M, int d=0)
{
   int i;
   int n = M.nc()>M.nr() ? M.nc() : M.nr();
   std::vector<T> r(n-d,0.);
   for(i=0; i<(n-d); i++)
      r[i] = M(i,i+d);
   return r;
}

// --------------------------------------------------------------------------------------------------------------- //
// Returns the upper or lower triangle as a sparse matrix
// --------------------------------------------------------------------------------------------------------------- //
template <class T> sMat<T> triu(const sMat<T> &M, int d)
{
   int n = (int)M.nr(); int m = (int)M.nc();
   sMat<T> r(n,m);
   std::vector< std::vector<int> > nz = M.find();
   int i,sz = (int)nz.size();

   for (i=0; i<sz; i++)
      if(nz[i][1]>=(nz[i][0]+d))
         r(nz[i][0],nz[i][1]) = M(nz[i][0],nz[i][1]);

   return r;
}
template <class T> sMat<T> triu(const sMat<T> &M)
{
   return triu<T>(M,0);
}
template <class T> sMat<T> tril(const sMat<T> &M, int d)
{
   int n = (int)M.nr(); int m = (int)M.nc();
   sMat<double> r(n,m);
   std::vector< std::vector<int> > nz = M.find();
   int i,sz = (int)nz.size();

   for (i=0; i<sz; i++)
      if(nz[i][1]<=(nz[i][0]+d))
         r(nz[i][0],nz[i][1]) = M(nz[i][0],nz[i][1]);

   return r;
}
template <class T> sMat<T> tril(const sMat<T> &M)
{
   return tril<T>(M,0);
}

// --------------------------------------------------------------------------------------------------------------- //
// Returns a vector whose element is 1 if all elements along a row/column is nonzero, and 0 otherwise
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::vector<T> all(const sMat<T> &M, int dim)
{
   int i,j,n,m;
   if(dim==2) { n = (int)M.nr(); m = (int)M.nc(); } else { n = (int)M.nc(); m = (int)M.nr(); }
   std::vector<T> retval(n,1);

   if(dim==2)
   {
      for(i=0; i<n; i++)
         for(j=0; j<m; j++)
            if(M(i,j)==0)
            {
               retval[i] = 0;
               break;
            }
   }
   else
   {
      for(i=0; i<n; i++)
         for(j=0; j<m; j++)
            if(M(j,i)==0)
            {
               retval[i] = 0;
               break;
            }
   }

   return retval;
}
template <class T> std::vector<T> all(sMat<T> & M)
{
   return all(M,1);
}

// --------------------------------------------------------------------------------------------------------------- //
// Sorts using the Gnome sort
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::vector<T> vsort(const std::vector<T> &v)
{
   unsigned int i=1,j=2;
   std::vector<T> r = v;
   T elem;
   while(i<v.size())
   {
      if(r[i-1]<=r[i]) { i=j; j++; }
      else { elem = r[i-1]; r[i-1] = r[i]; r[i] = elem; i--; if(i==0) i=1; }
   }
   return r;
}
template <class T> std::vector<T> vsort(const std::vector<T> &v, std::vector<int> &ind)
{
   unsigned int i=1,j=2;
   int ii;
   std::vector<T> r = v;
   T elem;

   ind.clear(); for(ii=0; ii<(int)v.size(); ii++) ind.push_back(ii);

   while(i<v.size())
   {
      if(r[i-1]<=r[i]) { i=j; j++; }
      else { elem = r[i-1]; r[i-1] = r[i]; r[i] = elem; ii=ind[i-1]; ind[i-1]=ind[i]; ind[i]=ii; i--; if(i==0) i=1; }
   }
   return r;
}
template <class T> sMat<T> msort(const sMat<T> &m)                   // Sorts along each column
{
   unsigned int i,j,col;
   T elem;
   sMat<T> r = m;
   for(col=0; col<(unsigned int)m.nc(); col++)
   {
      i=1; j=2;
      while(i<(unsigned int)m.nr())
      {
         if(r(i-1,col)<=r(i,col)) { i=j; j++; }
         else { elem = r(i-1,col); r(i-1,col) = r(i,col); r(i,col) = elem; i--; if(i==0) i=1; }
      }
   }
   return r;
}
template <class T> sMat<T> msort(const sMat<T> &m, sMat<int> &ind)   // Sorts along each column
{
   unsigned int i,j,col;
   int ielm;
   T elem;
   sMat<T> r = m;
   ind.zero(m.nr(),m.nc());
   for(col=0; col<(unsigned int)m.nc(); col++)
   {
      for(i=0; i<(unsigned int)m.nr(); i++) ind(i,col) = (int)i;
      i=1; j=2;
      while(i<(unsigned int)m.nr())
      {
         if(r(i-1,col)<=r(i,col)) { i=j; j++; }
         else 
         { 
            elem = r(i-1,col); r(i-1,col) = r(i,col); r(i,col) = elem; 
            ielm = ind(i-1,col); ind(i-1,col)=ind(i,col); ind(i,col)=ielm; i--; if(i==0) i=1; 
         }
      }
   }
   return r;
}
template <class T> sMat<T> msort(const sMat<T> &m, sMat<int> &ind, const char *descend) 
{
   unsigned int i,j,col;
   int ielm;
   T elem;
   sMat<T> r = m;
   ind.zero(m.nr(),m.nc());
   for(col=0; col<(unsigned int)m.nc(); col++)
   {
      for(i=0; i<(unsigned int)m.nr(); i++) ind(i,col) = (int)i;
      i=1; j=2;
      while(i<(unsigned int)m.nr())
      {
         if(r(i-1,col)>=r(i,col)) { i=j; j++; }
         else 
         { 
            elem = r(i-1,col); r(i-1,col) = r(i,col); r(i,col) = elem; 
            ielm = ind(i-1,col); ind(i-1,col)=ind(i,col); ind(i,col)=ielm; i--; if(i==0) i=1; 
         }
      }
   }
   return r;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the set-theoretic union, that is the elements that are in both a and b, same as Matlab union(a,b)
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::vector<T> setunion(const std::vector<T> &a, const std::vector<T> &b)
{
   unsigned int i;
   std::vector<T> r = a; r.reserve(a.size()+b.size());

   // Does the union
   for(i=0; i<b.size(); i++) r.push_back(b[i]);

   // Sorts using the Gnome sort
   r = vsort(r);

   // Removes duplicate
   i = 1;
   while(i<r.size())
      if(r[i]==r[i-1]) { r.erase(r.begin()+i); } else { i++; }

   return r;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the set-theoretic difference, that is the elements of a that is not in b
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::vector<T> setdiff(const std::vector<T> &a, const std::vector<T> &b)
{
   unsigned int i,j;
   std::vector<T> r = a;

   // For each element of b, loops over elements of a, and if same, delete element of a... not efficient!
   for(i=0; i<b.size(); i++)
      for(j=0; j<r.size(); j++)
         if(r[j]==b[i]) { r.erase(r.begin()+j); break; }
   
   // Sorts using the Gnome sort
   r = vsort(r);

   return r;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the set-theoretic xor, that is the elements of a and b that is not in the other
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::vector<T> setxor(const std::vector<T> &a, const std::vector<T> &b)
{
   return setunion(setdiff(a,b),setdiff(b,a));
}

// --------------------------------------------------------------------------------------------------------------- //
// Enumeration for the norm function
// --------------------------------------------------------------------------------------------------------------- //
enum normp { one = 1, lsvd = 2, inf = 3, fro = 4 };

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in maths.cpp
// --------------------------------------------------------------------------------------------------------------- //
float sign(float val);                                                    // Determines the sign of a float 
double sign(double val);                                                  // Determines the sign of a double
eigVE<double> eig(const sMat<double>  & M);                               // Diagonalises M by tri-diag and QL fact.
std::vector<double> svd(const sMat<double> &M);                           // Calculates the single value decomposition
std::vector<double> svd(const sMat<double> &M, sMat<double> &V);          // Calculates svd and the left matrix too
std::vector<double> svd(const sMat<double> &M, sMat<double> &U,sMat<double> &V); // Calculates svd and both U/V 
sMat<double> qr(const sMat<double> & M, sMat<double> & Q);                // Calculates the QR decomposition of M
sMat<double> qr(const sMat<double> & M, sMat<double> & Q, int zero);      // Calculates the 'economy' QR decomposition
sMat<double> qr(const sMat<double> & M, int zero);                        // Calculates only R of M = Q*R
sMat<double> mabs(const sMat<double> & M);                                // Calculates the absolute of each element
std::vector<double> vabs(const std::vector<double> & V);                  // Calculates the absolute of each element
double norm(const sMat<double> & M);                                      // Calculates the norm, in similar way to
double norm(const sMat<double> & M, normp p);                             //    the matlab function (see maths.cpp)
sMat<double> randn(int m, int n);                                         // Generates a normal-dist. random matrix
sMat<double> orth(const sMat<double> & M);                                // Generates an orthonormal basis using svd
int rank(const sMat<double> & M, double tol);                             // Estimate # linearly indep. rows/columns
void rmzeros(sMat<double> & M);                                           // Removes entries < eps
std::vector<double> f2vec(double *v, int n);                              // Converts a 1D C-array into a std::vector
sMat<double> f2mat(double *M, int m, int n);                              // Converts a 2D C-array into an sMat
complexdouble* zmat2f(sMat<double> &r, sMat<double> &i);                  // Converts two sMat to complex C-array
complexdouble spherical_harmonics(int k, int q, double th, double phi);   // Calculates the spherical harmonic
#endif

#endif
