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
#include<functional>

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
   //int _nzmax;                                                        // maximum number of entries
     int _m;                                                            // number of rows
     int _n;                                                            // number of columns
     std::vector<int> _p;                                               // column pointers (size n+1)
     std::vector<int> _i;                                               // row indices, size nzmax
     std::vector< std::complex<T> > _x;                                 // numerical values, size nzmax
   //int _nz;                                                           // # of entries in triplet matrix, -1 for compressed-col
     bool _iscsc;

     size_t _nnz;                                                       // A hash of the non-zero pattern (of the vector _p)
     void _genhash();                                                   // Generates a hash of the non-zero pattern [_p].
     std::vector<size_t> _subset;                                       // Vector of hashes of matrices that are superset of self

     // Algorithms from CXSparse. Reference: Direct Methods for Sparse Linear Systems, Timothy A. Davis, SIAM 2006.
     std::vector<int> _wi;                                              // Integer workspace for CSC algorithms
     std::vector< std::complex<T> > _wx;                                // Workspace for CSC algorithms
/*   std::vector<int> _csc_pinv;                                        // Row pivots
     std::vector<int> _csc_sym_pinv;                                    //  inverse row perm. for QR, fill red. perm for Chol
     std::vector<int> _csc_sym_q;                                       //  fill-reducing column permutation for LU and QR
     std::vector<int> _csc_sym_parent;                                  //  elimination tree for Cholesky and QR
     std::vector<int> _csc_sym_cp;                                      //  column pointers for Cholesky, row counts for QR
     std::vector<int> _csc_sym_leftmost;                                //  leftmost[i] = min(find(A(i,:))), for QR
     int _csc_sym_m2, unz, lnz;                                         //  m2=rows for QR, unz,lnz= entries in L,U or V,R
*/
     void _cs_compress();                                               // Converts from triplet to compressed sparse column format
     void _cs_dupl();                                                   // Sums duplicated entries
     void _cs_fkeep (int(*fkeep)(int,int,std::complex<T>,void*),void*); // Keeps entries where function fkeep()==true
     int _cs_scatter(int j, std::complex<T> b, int m, zsMat&C, int nz); // Scatter - construct non-zero column as a dense vector
     zsMat<T> _cs_multiply(const zsMat &B);                             // Calcs. the product C = self*B
     zsMat<T> _cs_add(const zsMat &B, bool sub=false);                  // Calcs. the sum C = self + B ( -B if sub=true)
     zsMat<T> _cs_adds(const zsMat &B, bool sub=false);                 // Calcs. the sum self += B ( -=B if sub=true)
     zsMat<T> _cs_add_pat(const zsMat &B, bool sub=false);              // Calcs. the sum self += B assuming B sparsity == self

   public:
     // Constructors and Destructors //
     zsMat(size_t r=1, size_t c=1): _m(r), _n(c), _iscsc(false) {};     // constructs an empty r x c sparse matrix
     ~zsMat() {};

     // Other member functions //
     int nr() const { return _m; };
     int nc() const { return _n; };
     std::vector< std::vector<int> > find() const;                      // Returns the indices of the nonzero elements
     zsMat<T> transpose();                                              // Transposes the matrix (without complex conjugation)
     zsMat<T> hermitian();                                              // Hermitian conjugate the matrix
     void del(int r, int c);                                            // Deletes an element
     void clear() { _p.clear(); _i.clear(); _x.clear(); _m=0; _n=0; };  // Clears the matrix
     void zero() { _p.clear(); _i.clear(); _x.clear(); }                // Zeros all entries
     void zero(int r, int c) { zero(); _m = r; _n = c; };               // Zeros all entries, changes shape as well
     void rmzeros();                                                    // Removes (numerically) zero entries
     std::complex<T> *f_array() const;                                  // Returns matrix as a Fortran style 2D array
     void f_array(std::complex<T>* retval) const;                       // Returns matrix as a Fortran style 2D array (pre-allocated)
     std::complex<T> *f_array_tr() const;                               // Returns matrix as a Fortran style 2D array
     std::complex<T> *h_array() const;                                  // Returns matrix as a Fortran style 2D array
     void h_array(std::complex<T>* retval) const;                       // Returns matrix as a Fortran style 2D array (pre-allocated)
     T* cp_array() const;                                               // Returns Hermitian matrix as a C-style packed 2D array
     T* fp_array() const;                                               // Returns Hermitian matrix as a Fortran packed 2D array
     Matrix fp_matrix() const;                                          // Returns Hermitian matrix as a Fortran packed 2D array

     // For iterative eigensolvers
     void MultMv(std::complex<T> *v, std::complex<T> *w);               // Calculates the matrix-vector product w = M*v (assume Hermitian)
     void MultMvNH(std::complex<T> *v, std::complex<T> *w);             // Calculates the matrix-vector product w = M*v (non-Hermitian)
     void MultMMH(std::complex<T> *A, std::complex<T> *B, int c);       // Calculates the matrix-matrix product A = M*B (assume Hermitian)
     void MultMM(std::complex<T> *A, std::complex<T> *B, int c);        // Calculates the matrix-matrix product A = M*B (non-Hermitian)

     // Overloaded operators
     zsMat<T> operator =  (const zsMat & m);                            // Copy assignment - overwrites previous matrix
     zsMat<T> operator =  (const T val);                                // Constant assignment (only diagonal elements)
     zsMat<T> operator =  (const std::complex<T> val);                  // Constant assignment (only diagonal elements)
     zsMat<T> operator += (const zsMat & m);                            // Add another matrix to current (element-wise)
     zsMat<T> operator +  (const zsMat & m);                            // Add another matrix to current (element-wise)
     zsMat<T> operator += (const T val);                                // Add a constant value to current (diagonal)
     zsMat<T> operator += (const std::complex<T> val);                  // Add a constant value to current (diagonal)
     zsMat<T> operator -= (const zsMat & m);                            // Subtract another matrix from current
     zsMat<T> operator -  (const zsMat & m);                            // Subtract another matrix from current
     zsMat<T> operator -= (const T val);                                // Subtract a constant value from current 
     zsMat<T> operator -= (const std::complex<T> val);                  // Subtract a constant value from current 
     zsMat<T> operator *= (const T & c);                                // Matrix scalar multiplication
     zsMat<T> operator *= (const std::complex<T> & c);                  // Matrix scalar multiplication
     zsMat<T> operator /= (const T & c);                                // Matrix scalar division
     zsMat<T> operator /= (const std::complex<T> & c);                  // Matrix scalar division
     zsMat<T> operator *= (const zsMat & m);                            // Matrix multiplication
     zsMat<T> operator *  (const zsMat & m);                            // Matrix multiplication
     zsMat<T> operator *= (const std::vector<T> & v);                   // Matrix*Vector (Matrix operating on a vector)
     zsMat<T> operator *= (const std::vector<std::complex<T> > & v);    // Matrix*Vector (Matrix operating on a vector)
     std::complex<T>& operator () (int r, int c);                       // Gets/sets element by subscript index
     std::complex<T>  operator () (int r, int c) const;
     std::complex<T>  operator [] (std::pair<int,int>);                 // Read-only subscript operator

     // Algorithms from CXSparse. Reference: Direct Methods for Sparse Linear Systems, Timothy A. Davis, SIAM 2006.
     zsMat<T> cs_add(zsMat &B, std::complex<T> a, std::complex<T> b);   // Returns the sum C = a*self + b*B, a,b=scalar
     zsMat<T> add_scal(const std::complex<T> v);                        // Returns the sum C = self + v*I, I=identity
     void to_csc() { _cs_compress(); _cs_dupl(); }                      // Converts from triplet to compressed column form
     void to_tri();                                                     // Converts from compressed column to triplet form
     void tridupl();                                                    // Removes duplicated elements in triplet form
     bool is_subset(const zsMat &B);                                    // Determines if sparsity of B is subset of self

     // Friend function  to provide input/output via <iostream>
     friend std::ostream & operator << <> (std::ostream & o, const zsMat & m);
   //friend std::istream & operator >> <> (std::istream & i, zsMat & m);

};  // End of template <class T> class zsMat

// --------------------------------------------------------------------------------------------------------------- //
// Functions copied from CSparse package. Reference: Direct Methods for Sparse Linear Systems, Timothy A. Davis
//    (2006) SIAM, Philadephia. https://www.cise.ufl.edu/research/sparse/CSparse/
// --------------------------------------------------------------------------------------------------------------- //
inline double cs_cumsum(std::vector<int> &p, std::vector<int> &c, int n)
{
   int nz=0; double nz2=0;
   if(p.empty()||c.empty()) return(-1);
   for(int i=0; i<n; i++)
   {
      p[i] = nz; nz += c[i]; nz2 += c[i]; c[i] = p[i];  // sum in nz2 double to avoid int overflow
   }
   p[n] = nz; return nz2;
}
template <class T> void zsMat<T>::_cs_compress()
{
   if(_iscsc) return;
   int nz = (int)_x.size(), p;
   std::vector<int> Cp(_n+1,0), w(_n,0), Ci(nz,0); std::vector< std::complex<T> > Cx(nz,0);
   for(int k=0; k<nz; k++) w[_p[k]]++;                  // Column counts
   cs_cumsum(Cp,w,_n);                                  // Column pointers
   for(int k=0; k<nz; k++)
   {
      Ci[p=w[_p[k]]++] = _i[k];                         // A(i,j) is the pth entry in the compressed matrix
      Cx[p] = _x[k];
   }
   _p = Cp; _i = Ci; _x = Cx; 
   _iscsc = true;
   _genhash();
}
template <class T> void zsMat<T>::_cs_dupl()
{
   if(!_iscsc) return;
   int nz = 0, i,j, q;
   std::vector<int> w(_m,0);
   for (i=0; i<_m; i++) w[i] = -1;                      // Row i not yet seen
   for (j=0; j<_n; j++)
   {
      q = nz;                                           // Column j will start at q
      for (int p=_p[j]; p<_p[j+1]; p++)
      {
         i = _i[p];                                     // A(i,j) is non-zero
         if(w[i]>=q)
            _x[w[i]] += _x[p];                          // A(i,j) is a duplicate
         else
         {
            w[i] = nz;                                  // Record where row i occurs
            _i[nz] = i;                                 // Keep A(i,j)
            _x[nz++] = _x[p];
         }
      }
      _p[j] = q;                                        // Record start of column j
   }
   _p[_n] = nz;                                         // Finalize A
   _i.erase(_i.begin()+nz,_i.end());
   _x.erase(_x.begin()+nz,_x.end());
   _genhash();
}
template <class T> void zsMat<T>::to_tri()              // Converts from compressed sparse column to triplet format
{
   if(!_iscsc) return;
   if(_x.empty()) { _iscsc = false; return; }
   int nz = (int)_x.size(), c=0;
   std::vector<int> Cp(nz,0), Ci(nz,0); std::vector< std::complex<T> > Cx(nz,0);
   for(int j=0; j<_n; j++) for(int i=_p[j]; i<_p[j+1]; i++) {
      Cp[c] = j; Ci[c] = _i[i]; Cx[c++] = _x[i]; }
   _p = Cp; _i = Ci; _x = Cx; _iscsc = false;
}
template <class T> void zsMat<T>::tridupl()             // Sums duplicated entries in triplet format
{
   if(_iscsc) return;
   // Do the gnome sort of the entries
   int i=1,j=2,vp,vi;
   std::complex<T> vx;
   while(i<_p.size())
   {
      if( (_p[i-1]<=_p[i-1])&&(_i[i-1]<=_i[i-1]) ) { i=j; j++; }
      else { 
         vp=_p[i-1]; _p[i-1]=_p[i]; _p[i]=vp; 
         vi=_i[i-1]; _i[i-1]=_i[i]; _i[i]=vi; 
         vx=_x[i-1]; _x[i-1]=_x[i]; _x[i]=vx; i--; if(i==0) i=1; }
   }
   // Runs through sorted list and sums up duplicates.
   j = 0;
   for(i=1; i<_p.size(); i++)
   {
      if(_p[j]==_p[i] && _i[j]==_i[i])
         _x[j]+=_x[i];
      else
      {
         j++; 
         if(i!=j) { 
            _p[j]=_p[i]; _i[j]=_i[i]; _x[j]=_x[i]; }
      }
   }
   _p.erase(_p.begin()+j,_p.end());
   _i.erase(_i.begin()+j,_i.end());
   _x.erase(_x.begin()+j,_x.end());
}
template <class T> void zsMat<T>::_cs_fkeep (int(*fkeep)(int,int,std::complex<T>,void*),void*other)
{                                                       // Keeps an entry if the fkeep function evalates to true 
   int p, nz=0;
   for(int j=0; j<_n; j++)
   {
      p = _p[j];                                        // Get current location of j
      _p[j] = nz;                                       // Record new location of j
      for(; p<_p[j+1]; p++)
      {
         if(fkeep(_i[p],j,_x[p],other)) {
            _x[nz] = _x[p]; _i[nz++] = _i[p]; }         // Keep A(i,j)
      }
   }
   _p[_n] = nz;                                         // Finalize A
   _i.erase(_i.begin()+nz,_i.end());
   _x.erase(_x.begin()+nz,_x.end());
   _genhash();
}
inline int _cs_nonzero(int i, int j, std::complex<double> aij, void *other) {
   return (std::abs(aij)<(DBL_EPSILON*1000)); }
inline int _cs_nonzero(int i, int j, std::complex<float> aij, void *other) {
   return (std::abs(aij)<(FLT_EPSILON*10)); }
template <class T> int _cs_nonzero(int i, int j, std::complex<T> aij, void *other) {
   return (std::abs(aij)==0); }
template <class T> void zsMat<T>::rmzeros()
{
   _cs_fkeep(&_cs_nonzero,NULL);
}
template <class T> int zsMat<T>::_cs_scatter(int j, std::complex<T> beta, int mark, zsMat<T> &c, int nz)
{
   if(!_iscsc) { _cs_compress(); _cs_dupl(); }
   int i;
   for (int p=_p[j]; p<_p[j+1]; p++)
   {
      i = _i[p];                                        // A(i,j) is nonzero
      if (_wi[i]<mark)
      {
         _wi[i] = mark;                                 // i is new entry in column j
         c._i[nz++] = i;                                // add i to pattern of C(:,j)
         _wx[i] = beta * _x[p];                         // wx(i) = beta * A(i,j)
      }
      else _wx[i] += beta * _x[p];                      // i exists in C(:,j) already
   }
   return nz;
}
template <class T> zsMat<T> zsMat<T>::_cs_multiply(const zsMat<T> &B)
{
   if(!_iscsc) { _cs_compress(); _cs_dupl(); }
   int nz=0;
   zsMat<T> C(_m,B._n); 
   C._p.assign(B._n,0); 
   C._i.assign(_x.size()+B._x.size(),0); 
   C._x.assign(_x.size()+B._x.size(),0); 
   C._iscsc=true;
   _wi.assign(_m,0); _wx.assign(_m,0);
   for(int j=0; j<B._n; j++)
   {
      if((nz+_m)>(C._x.size())) { C._i.resize(C._i.size()*2+_m,0); C._x.resize(C._x.size()*2+_m,0); }
      C._p[j] = nz;                                     // column j of C starts here
      for(int p=B._p[j]; p<B._p[j+1]; p++)
      {
         nz = _cs_scatter(B._i[p], B._x[p], j+1, C, nz);
      }
      for(int p=C._p[j]; p<nz; p++) C._x[p] = _wx[C._i[p]];
   }
   C._p[B._n] = nz;                                     // Finalize C
   C._i.erase(C._i.begin()+nz,C._i.end());
   C._x.erase(C._x.begin()+nz,C._x.end());
   C._genhash();
   return C;
}
template <class T> zsMat<T> zsMat<T>::cs_add(zsMat &B, std::complex<T> alpha, std::complex<T> beta)
{
   if(!_iscsc) { _cs_compress(); _cs_dupl(); }
   if(!B._iscsc) { B._cs_compress(); B._cs_dupl(); }
   int nz=0, i;
   zsMat<T> C(_m,B._n); 
   C._p.assign(B._n,0); 
   C._i.assign(_x.size()+B._x.size(),0); 
   C._x.assign(_x.size()+B._x.size(),0); 
   C._iscsc=true;
   _wi.assign(_m,0); _wx.assign(_m,0);
   for(int j=0; j<B._n; j++)
   {
      C._p[j] = nz;                                     // column j of C starts here
      nz = _cs_scatter(j, alpha, j+1, C, nz);           // alpha * A(:,j)
      for(int p=B._p[j]; p<B._p[j+1]; p++) {            // beta * B(:,j)
         i = B._i[p];
         if(_wi[i]<j+1) { 
            _wi[i] = j+1; C._i[nz++] = i; _wx[i] = beta * B._x[p]; }
         else _wx[i] += beta * B._x[p];
      }
      for(int p=C._p[j]; p<nz; p++) C._x[p] = _wx[C._i[p]];
   }
   C._p[B._n] = nz;                                     // Finalize C
   C._i.erase(C._i.begin()+nz,C._i.end());
   C._x.erase(C._x.begin()+nz,C._x.end());
   C._genhash();
   return C;
}
template <class T> zsMat<T> zsMat<T>::_cs_add(const zsMat &B, bool sub)
{
   if(!_iscsc) { _cs_compress(); _cs_dupl(); }
   int nz=0, i;
   zsMat<T> C(_m,B._n); 
   C._p.assign(B._n+1,0); 
   C._i.assign(_x.size()+B._x.size(),0); 
   C._x.assign(_x.size()+B._x.size(),0); 
   C._iscsc=true;
   _wi.assign(_m+1,0); _wx.assign(_m+1,0);
   for(int j=0; j<B._n; j++)
   {
      C._p[j] = nz;                                     // column j of C starts here
      for(int p=_p[j]; p<_p[j+1]; p++) { i = _i[p];     // A(:,j)
         if(_wi[i]<j+1) { _wi[i] = j+1; C._i[nz++] = i; _wx[i] = _x[p]; }
         else _wx[i] += _x[p];
      }
      for(int p=B._p[j]; p<B._p[j+1]; p++) { i=B._i[p]; // B(:,j)
         if(_wi[i]<j+1) { _wi[i] = j+1; C._i[nz++] = i; _wx[i] = (sub) ? -B._x[p] : B._x[p]; } 
         else { if(sub) _wx[i] -= B._x[p]; else _wx[i] += B._x[p]; }
      }
      for(int p=C._p[j]; p<nz; p++) C._x[p] = _wx[C._i[p]];
   }
   C._p[B._n] = nz;                                     // Finalize C
   C._i.erase(C._i.begin()+nz,C._i.end());
   C._x.erase(C._x.begin()+nz,C._x.end());
   C._genhash();
   C._subset.push_back(B._nnz);
   return C;
}
template <class T> zsMat<T> zsMat<T>::_cs_adds(const zsMat &B, bool sub)
{
   if(!_iscsc) { _cs_compress(); _cs_dupl(); }
   int nz=0, i;
   std::vector<int> Cp(B._n+1,0); 
   std::vector<int> Ci(_x.size()+B._x.size(),0); 
   std::vector<std::complex<T> > Cx(_x.size()+B._x.size(),0); 
   _wi.assign(_m+1,0); _wx.assign(_m+1,0);
   for(int j=0; j<B._n; j++)
   {
      Cp[j] = nz;                                       // column j of C starts here
      for(int p=_p[j]; p<_p[j+1]; p++) { i = _i[p];     // A(:,j)
         if(_wi[i]<j+1) { _wi[i] = j+1; Ci[nz++] = i; _wx[i] = _x[p]; }
         else _wx[i] += _x[p];
      }
      for(int p=B._p[j]; p<B._p[j+1]; p++) { i=B._i[p]; // B(:,j)
         if(_wi[i]<j+1) { _wi[i] = j+1; Ci[nz++] = i; _wx[i] = (sub) ? -B._x[p] : B._x[p]; } 
         else { if(sub) _wx[i] -= B._x[p]; else _wx[i] += B._x[p]; }
      }
      for(int p=Cp[j]; p<nz; p++) Cx[p] = _wx[Ci[p]];
   }
   Cp[B._n] = nz; _p = Cp;                              // Finalize C
   Ci.erase(Ci.begin()+nz,Ci.end()); _i = Ci;
   Cx.erase(Cx.begin()+nz,Cx.end()); _x = Cx;
   _genhash();
   _subset.push_back(B._nnz);
   return *this;
}
template <class T> zsMat<T> zsMat<T>::_cs_add_pat(const zsMat &B, bool sub)
{                                                       // Assume non-zero pattern of B is subset of A
   if(!_iscsc) return *this;
   int i;
   _wi.assign(_m+1,0); _wx.assign(_m+1,0);
   for(int j=0; j<B._n; j++)
   {
      i = _p[j];
      for(int p=B._p[j]; p<B._p[j+1]; p++) {            // Scatter B(:,j)
         for(i=_p[j]; i<_p[j+1]; i++) {
            if(_i[i]==B._i[p]) { if(sub) _x[i] -= B._x[p]; else _x[i] += B._x[p]; break; } }
      }
   }
   return *this;
}
template <class T> zsMat<T> zsMat<T>::add_scal(const std::complex<T> v) // Adds a constant to the current matrix
{
   if(_iscsc)
   {
      for(int j=0; j<_n; j++)
      {
         for(int i=_p[j]; i<_p[j+1]; i++)
         {
            if(_i[i]==j) { _x[i] += v; break; }
            else if(_i[i]>j) {
               _i.insert(_i.begin()+i-1,j); _x.insert(_x.begin()+i-1,v);
               for(int m=j+1; m<_n; m++) _p[m]++; 
               break;
            }
         }
      }
      _genhash();
   }
   else
   {
      int sz = (_n>_m) ? _m : _n;
      std::vector<int> isthere(sz,0);
      for(int n=0; n<_x.size(); n++)
      {
         if(_i[n]==_p[n]) { isthere[_i[n]] = 1; _x[n] += v; }
      }
      for(int n=0; n<sz; n++) if(isthere[n]==0) {
         _p.push_back(n); _i.push_back(n); _x.push_back(v); }
   }
   return *this;
}

// --------------------------------------------------------------------------------------------------------------- //
// Member functions
// --------------------------------------------------------------------------------------------------------------- //
template <class T> void zsMat<T>::_genhash()
{
   if(_x.empty()) { _nnz = 0; return; }  // Empty matrix.
   std::hash<int> hash_value;
   _nnz = hash_value(_p[0]+_i[0]);
   // Hash-combiner stolen from Boost: http://www.boost.org/doc/html/hash/reference.html#boost.hash_combine
   for(int j=1; j<_n; j++) for(int i=_p[j]; i<_p[j+1]; i++) {
      _nnz ^= hash_value(_p[j]+_i[i]) + 0x9e3779b9 + (_nnz << 6) + (_nnz >> 2); }
   _subset.clear();
}
template <class T> bool zsMat<T>::is_subset(const zsMat<T> &B)
{
   if(_nnz==B._nnz) return true;
   for(int i=0; i<_subset.size(); i++) { if(B._nnz==_subset[i]) return true; }
   return false;
}
template <class T> std::vector< std::vector<int> > zsMat<T>::find() const
{
   std::vector<int> row(2);
   std::vector< std::vector<int> > retval(_x.size(),row);
   int r=0;

   if(_iscsc)
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++) {
         row[0] = _i[c]+1; row[1] = j+1; retval[r] = row; r++; }
   }
   else
   {
      for(r=0; r<_x.size(); r++) { 
         row[0] = _i[r]+1; row[1] = _p[r]+1; retval[r] = row; }
   }
   return retval;
}

template <class T> zsMat<T> zsMat<T>::transpose()
{
   if(_iscsc)
   {
      int nz = (int)_x.size(), q;
      std::vector<int> Cp(_n+1,0), w(_n,0), Ci(nz,0); std::vector< std::complex<T> > Cx(nz,0);
      for(int k=0; k<nz; k++) w[_i[k]]++;               // Row counts
      cs_cumsum(Cp,w,_m);                               // Row pointers
      for(int j=0; j<_n; j++)
      {  
         for(int p=_p[j]; p<_p[j+1]; p++)
         {
            Ci[q=w[_i[p]]++] = j;                       // Place A(i,j) as C(j,i)
            Cx[q] = _x[p];
         }
      }
      _p = Cp; _i = Ci; _x = Cx;
      _genhash();
   }
   else
   {
      int t;
      for(int i=0; i<_x.size(); i++) { 
         t=_p[i]; _p[i]=_i[i]; _i[i]=t; }
   }
   return *this;
}

template <class T> zsMat<T> zsMat<T>::hermitian()
{
   if(_iscsc)
   {
      int nz = (int)_x.size(), q;
      std::vector<int> Cp(_n+1,0), w(_n,0), Ci(nz,0); std::vector< std::complex<T> > Cx(nz,0);
      for(int k=0; k<nz; k++) w[_i[k]]++;               // Row counts
      cs_cumsum(Cp,w,_m);                               // Row pointers
      for(int j=0; j<_n; j++)
      {  
         for(int p=_p[j]; p<_p[j+1]; p++)
         {
            Ci[q=w[_i[p]]++] = j;                       // Place A(i,j) as C(j,i)
            Cx[q] = conj(_x[p]);
         }
      }
      _p = Cp; _i = Ci; _x = Cx;
      _genhash();
   }
   else
   {
      int t;
      for(int i=0; i<_x.size(); i++) { 
         t=_p[i]; _p[i]=_i[i]; _i[i]=t; }
   }
   return *this;
}

template <class T> void zsMat<T>::del(int r, int c)             // Deletes an element
{
   if(_iscsc)
   {
      for(int i=_p[c]; i<_p[c+1]; i++) if(_i[i]==r) 
      { 
         _i.erase(_i.begin()+i); _x.erase(_x.begin()+i);
         for(int j=c+1; j<=_n; j++) _p[j]--; 
         break;
      }
      _genhash();
   }
   else
   {
      for(int i=0; i<_x.size(); i++) if(_i[i]==_p[i])
      {
         _p.erase(_p.begin()+i); _i.erase(_i.begin()+i); _x.erase(_x.begin()+i);
         break;
      }
   }
}

template <class T> std::complex<T>* zsMat<T>::f_array() const   // Returns matrix as a Fortran style 2D array
{
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   std::complex<T> *retval;
   // Allocates an _r*_c array and initiallises all elements to zero.
   retval = (std::complex<T>*) calloc(_m*_n,sizeof(std::complex<T>));

   if(_iscsc)
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++)
         retval[_m*j+_i[c]] = _x[c];
   }
   else
   {
      for(int i=0; i<_x.size(); i++) retval[_m*_p[i]+_i[i]] += _x[i]; 
   }

   return retval;
}
template <class T> void zsMat<T>::f_array(std::complex<T>*retval) const   // Assumes matrix already allocated
{
   memset(retval,0,_m*_n*sizeof(std::complex<T>));
   if(_iscsc) {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++)
         retval[_m*j+_i[c]] = _x[c]; }
   else {
      for(int i=0; i<_x.size(); i++) retval[_m*_p[i]+_i[i]] += _x[i]; }
}
template <class T> std::complex<T>* zsMat<T>::h_array() const   // Returns matrix as a Fortran style 2D array
{
   std::complex<T> *retval = (std::complex<T>*) calloc(_m*_n,sizeof(std::complex<T>));
   // _Assume_(!) the matrix is Hermitian and only loop through the lower triangle.
   if(_iscsc) 
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++)
      {
         if(_i[c]>=j) {
            retval[_m*_i[c]+j] = conj(_x[c]); if(_i[c]!=j) retval[_n*j+_i[c]] = _x[c];
         }
      }
   }
   else 
   {
      for(int i=0; i<_x.size(); i++) {
         if(_i[i]>=_p[i]) {
            retval[_m*_i[i]+_p[i]] += conj(_x[i]); if(_i[i]!=_p[i]) retval[_n*_p[i]+_i[i]] += _x[i]; } }
   }
   return retval;
}
template <class T> void zsMat<T>::h_array(std::complex<T>* retval) const   // Assumes matrix already allocated
{
   memset(retval,0,_m*_n*sizeof(std::complex<T>));
   if(_iscsc) 
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++)
      {
         if(_i[c]>=j) {
            retval[_m*_i[c]+j] = conj(_x[c]); if(_i[c]!=j) retval[_n*j+_i[c]] = _x[c];
         }
      }
   }
   else 
   {
      for(int i=0; i<_x.size(); i++) {
         // Need to use += here because triplet forms have duplicate elements that should be summed together.
         if(_i[i]>=_p[i]) {
            retval[_m*_i[i]+_p[i]] += conj(_x[i]); if(_i[i]!=_p[i]) retval[_n*_p[i]+_i[i]] += _x[i]; } }
   }
}
template <class T> std::complex<T>* zsMat<T>::f_array_tr() const// Returns the transposed matrix as a Fortran style 2D array
{
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   std::complex<T> *retval = (std::complex<T>*) calloc(_m*_n,sizeof(std::complex<T>));

   if(_iscsc) 
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++)
         retval[_n*_i[c]+j] = _x[c];
   }
   else 
   {
      for(int i=0; i<_x.size(); i++) retval[_n*_i[i]+_p[i]] += _x[i];
   }
   return retval;
}
template <class T> T* zsMat<T>::fp_array() const                // Returns Hermitian matrix as a Fortran style packed 2D array
{                                                               // with real in upper, imag in lower triangles.
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   T *retval = (T*) calloc(_m*_n,sizeof(T));                    // Allocates an _r*_c array and initiallises all elements to zero.

   // _Assume_(!) the matrix is Hermitian and only loop through the lower triangle.
   if(_iscsc) 
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++)
      {
         if(_i[c]>=j) {
            retval[_m*_i[c]+j] = real(_x[c]); if(_i[c]!=j) retval[_n*j+_i[c]] = imag(_x[c]);
         }
      }
   }
   else 
   {
      for(int i=0; i<_x.size(); i++) {
         if(_i[i]>=_p[i]) {
            retval[_m*_i[i]+_p[i]] += real(_x[i]); if(_i[i]!=_p[i]) retval[_n*_p[i]+_i[i]] += imag(_x[i]); } }
   }

   return retval;
}
template <class T> Matrix zsMat<T>::fp_matrix() const           // Returns Hermitian matrix as a Fortran style packed 2D array
{                                                               // with real in upper, imag in lower triangles.
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   Matrix retval(1,_m,1,_n); retval=0;

   // _Assume_(!) the matrix is Hermitian and only loop through the lower triangle.
   if(_iscsc) 
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++)
      {
         if(_i[c]>=j) {
            retval(_i[c]+1,j+1) = real(_x[c]); if(_i[c]!=j) retval(j+1,_i[c]+1) = imag(_x[c]);
         }
      }
   }
   else 
   {
      for(int i=0; i<_x.size(); i++) {
         if(_i[i]>=_p[i]) {
            retval(_i[i]+1,_p[i]+1) += real(_x[i]); if(_i[i]!=_p[i]) retval(_p[i]+1,_i[i]+1) += imag(_x[i]); }
      }
   }

   return retval;
}
template <class T> T* zsMat<T>::cp_array() const                // Returns Hermitian matrix as a C style packed 2D array
{                                                               // with real in upper, imag in lower triangles.
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   T *retval = (T*) calloc(_m*_n,sizeof(T));                    // Allocates an _r*_c array and initiallises all elements to zero.

   // _Assume_(!) the matrix is Hermitian and only loop through the lower triangle.
   if(_iscsc) 
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++)
      {
         if(_i[c]>=j) {
            retval[_m*j+_i[c]] = real(_x[c]); if(_i[c]!=j) retval[_m*_i[c]+j] = imag(_x[c]);
         }
      }
   }
   else 
   {
      for(int i=0; i<_x.size(); i++) {
         if(_i[i]>=_p[i]) {
            retval[_m*_p[i]+_i[i]] += real(_x[i]); if(_i[i]!=_p[i]) retval[_m*_i[i]+_p[i]] += imag(_x[i]); }
      }
   }

   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Matrix-Vector multiplication for iterative methods
// --------------------------------------------------------------------------------------------------------------- //
template <class T> void zsMat<T>::MultMv(std::complex<T> *v, std::complex<T> *w)  // Calculates the matrix-vector product w = M*v
{                                                                                 // Needed by ARPACK
   // We have to assume that the size of the vector v is equal to _r
   memset(w,0,_m*sizeof(std::complex<T>));
   if(_iscsc) 
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++)
      {
         // Assume(!) matrix is Hermitian and loops only over the lower triangle
         if(_i[c]>=j) {
            w[_i[c]] += _x[c] * v[j]; if(_i[c]!=j) w[j] += conj(_x[c]) * v[_i[c]];
         }
      }
   }
   else 
   {
      for(int c=0; c<_x.size(); c++) {
         if(_i[c]<_p[c]) w[_p[c]] += conj(_x[c]) * v[_i[c]];
         if(_i[c]<=_p[c]) w[_i[c]] += _x[c] * v[_p[c]];  }
   }
}
template <class T> void zsMat<T>::MultMvNH(std::complex<T>*v, std::complex<T>*w)  // Calculates the matrix-vector product w = M*v
{                                                                                 // Needed by ARPACK
   // We have to assume that the size of the vector v is equal to _r
   memset(w,0,_m*sizeof(std::complex<T>));
   if(_iscsc) 
   {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++) {
         w[_i[c]] += _x[c] * v[j]; }
   }
   else {
      for(int c=0; c<_x.size(); c++) w[_i[c]] += _x[c] * v[_p[c]]; }
}
template <class T> void zsMat<T>::MultMM(std::complex<T>*A, std::complex<T>*B, int c) 
{                                                                                 // Calc. matrix-matrix product A=M*B
   // Assume A is correct size _r*c, B also correct, _c*c
   memset(A,0,_m*c*sizeof(std::complex<T>));
   if(_iscsc) 
   {
      for(int j=0; j<c; j++)
      {
         for(int k=0; k<_n; k++) for(int n=_p[j]; n<_p[j+1]; n++) {
            A[_m*j+_i[n]] += _x[n] * B[_n*j+k]; }
      }
   }
   else 
   {
      for(int j=0; j<c; j++) {
         for(int n=0; n<_x.size(); n++) A[_m*j+_i[n]] += _x[n] * B[_n*j+_p[n]]; }
   }
}
template <class T> void zsMat<T>::MultMMH(std::complex<T>*A, std::complex<T>*B, int c) 
{                                                                                 // Calc. matrix-matrix product A=M*B (Hermitian)
   // Assume A is correct size _r*c, B also correct, _c*c
   memset(A,0,_m*c*sizeof(std::complex<T>));
   if(_iscsc) 
   {
      for(int k=0; k<c; k++)   // Aik = sum_j Mij Bjk (calculates each column k of A in turn)
      {
         for(int j=0; j<_n; j++)
            for(int n=_p[j]; n<_p[j+1]; n++) {
               if(_i[n]>=j) A[_m*k+_i[n]] += _x[n] * B[_n*k+j]; if(_i[n]!=j) A[_m*k+j] += conj(_x[n]) * B[_n*k+_i[n]];
            }
      }
   }
   else 
   {
      for(int j=0; j<c; j++) {
         for(int n=0; n<_x.size(); n++) {
            if(_i[n]>=_p[n]) {
               A[_m*j+_i[n]] += _x[n] * B[_n*j+_p[n]]; if(_i[n]!=_p[n]) A[_m*j+_p[n]] += conj(_x[n]) * B[_n*j+_i[n]]; }
         }
      }
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Overloaded operators
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::complex<T>& zsMat<T>::operator () (int r, int c)
{ 
   if(r>_m) _m = r; if(c>_n) _n = c;
   r--; c--;
   if(_iscsc)
   {
/*    if((double)r>(double)_m/2.)        // Try to speed up search by looking from either ends
      {
         for(int n=_p[c+1]-1; n>=_p[c]; n--) {
            if(_i[n]<r)
            {
               _i.insert(_i.begin()+n,r); _x.insert(_x.begin()+n,std::complex<T>(0,0));
               for(int m=c+1; m<=_n; m++) _p[m]++; 
               _genhash();
               return _x.at(n);
            }
            else if(_i[n]==r) return _x.at(n);
         }
      }
      else
*/    {
         if(_p[c]==_p[c+1]) {
            _i.insert(_i.begin()+_p[c],r); _x.insert(_x.begin()+_p[c],std::complex<T>(0,0));
            for(int m=c+1; m<=_n; m++) _p[m]++; 
            _genhash();
            return _x.at(_p[c]);
         }
         else
         {
            for(int n=_p[c]; n<_p[c+1]; n++) {
               if(_i[n]>r)
               {
                  _i.insert(_i.begin()+n,r); _x.insert(_x.begin()+n,std::complex<T>(0,0));
                  for(int m=c+1; m<=_n; m++) _p[m]++; 
                  _genhash();
                  return _x.at(n);
               }
               else if(_i[n]==r) return _x.at(n);
            }
            // This should run if the column is empty or the required element is before the current first element
            _i.insert(_i.begin()+_p[c],r); _x.insert(_x.begin()+_p[c],std::complex<T>(0,0));
            for(int m=c+1; m<=_n; m++) _p[m]++; 
            _genhash();
            return _x.at(_p[c]);
         }
      }
   }
   else
   {
      _p.push_back(c); _i.push_back(r); _x.push_back(std::complex<T>(0,0)); 
      return _x.back();
   }
   std::cerr << "zsMat::Error in setting element (" << r << "," << c << ")\n"; exit(-1);
}
template <class T> std::complex<T> zsMat<T>::operator () (int r, int c) const
{
   if (r>_m || c>_n) return std::complex<T> (0,0);
   r--; c--;
   if(_iscsc)
   {
/*    if((double)r>(double)_m/2.)        // Try to speed up search by looking from either ends
      {
         for(int n=_p[c+1]-1; n>=_p[c]; n--) {
            if(_i[n]>r) return std::complex<T> (0,0);
            else if(_i[n]==r) return _x[n];
         }
         return std::complex<T> (0,0);
      }
      else
*/    {
         for(int n=_p[c]; n<_p[c+1]; n++) {
          //if(_i[n]<r) return std::complex<T> (0,0);
          /*else*/ if(_i[n]==r) return _x[n];
         }
         return std::complex<T> (0,0);
      }
   }
   else
   {
      for(int i=0; i<_x.size(); i++)
         if(_p[i]==c && _i[i]==r) return &_x[i];
      return std::complex<T> (0,0);
   }
}
template <class T> std::complex<T> zsMat<T>::operator [] (std::pair<int,int> rc)
{
   int r=rc.first-1, c=rc.second-1;
   if (r>_m || c>_n) return std::complex<T> (0,0);
   if(_iscsc)
   {
/*    if((double)r>(double)_m/2.)        // Try to speed up search by looking from either ends
      {
         for(int n=_p[c+1]-1; n>=_p[c]; n--) {
            if(_i[n]>r) return std::complex<T> (0,0);
            else if(_i[n]==r) return _x[n];
         }
         return std::complex<T> (0,0);
      }
      else
*/    {
         for(int n=_p[c]; n<_p[c+1]; n++) {
          //if(_i[n]<r) return std::complex<T> (0,0);
          /*else*/ if(_i[n]==r) return _x[n];
         }
         return std::complex<T> (0,0);
      }
   }
   else
   {
      std::complex<T> rv(0,0); 
      for(int i=0; i<_x.size(); i++)
         if(_p[i]==c && _i[i]==r) rv += _x[i];
      return rv;
   }
   return std::complex<T> (0,0);
}

template <class T> zsMat<T> zsMat<T>::operator = (const zsMat<T> & m)
{
   _m = m._m; _n = m._n; _p = m._p; _i = m._i; _x = m._x; 
   _iscsc = m._iscsc; _nnz = m._nnz; _subset = m._subset;
   return *this;
}
template <class T> zsMat<T> zsMat<T>::operator = (const T val)          // Assignment of value to diagonal
{
   _p.clear(); _i.clear(); _x.clear();
   if(val!=0)
   {
      int i_max = (_n<_m)?_n:_m;
      _p.assign(_n,0); _i.assign(i_max,0); _x.assign(i_max,0);
      for(int n=0; n<i_max; n++) { 
         _p[n]=n; _i[n]=n; _x[n]=std::complex<T>(val,0); }
      for(int j=i_max; j<_n; j++) _p[j]=i_max;
      if(_iscsc) _genhash();
   }
   return *this;
}
template <class T> zsMat<T> zsMat<T>::operator = (const std::complex<T> val) 
{
   _p.clear(); _i.clear(); _x.clear();
   if(val!=0)
   {
      int i_max = (_n<_m)?_n:_m;
      _p.assign(_n,0); _i.assign(i_max,0); _x.assign(i_max,0);
      for(int n=0; n<i_max; n++) { 
         _p[n]=n; _i[n]=n; _x[n]=val; }
      for(int j=i_max; j<_n; j++) _p[j]=i_max;
      if(_iscsc) _genhash();
   }
   return *this;
}

template <class T> zsMat<T> zsMat<T>::operator += (const zsMat<T> &m)   // Adds another matrix to current
{
   if(_n!=m._n || _m!=m._m) return *this;
   if(_iscsc)
   {
      if(m._iscsc) {
         if(is_subset(m)) return _cs_add_pat(m); else return _cs_adds(m); }
      else { 
         zsMat<T> n(m); n._cs_compress(); n._cs_dupl(); return _cs_adds(n); } 
   }
   else
   {
      _p.insert(_p.end(),m._p.begin(),m._p.end());
      _i.insert(_i.end(),m._i.begin(),m._i.end());
      _x.insert(_x.end(),m._x.begin(),m._x.end());
      return *this;
   }
}
template <class T> zsMat<T> zsMat<T>::operator +  (const zsMat<T> &m) { // Binary addition of two sparse matrices
   if(m._iscsc) return _cs_add(m); else { zsMat<T> n(m); n._cs_compress(); n._cs_dupl(); return _cs_add(n); } }

template <class T> zsMat<T> zsMat<T>::operator += (const T val)         // Adds a constant to the current matrix
{
   return add_scal(std::complex<T>(val,0));
}
template <class T> zsMat<T> zsMat<T>::operator += (const std::complex<T> val)
{
   return add_scal(val);
}
template <class T> zsMat<T> operator + (const zsMat<T> & m, const T val) {
   zsMat<T> tmp = m; tmp += val; return tmp; }
template <class T> zsMat<T> operator + (const T val, const zsMat<T> & m) {
   zsMat<T> tmp = m; tmp += val; return tmp; }
template <class T> zsMat<T> operator + (const zsMat<T> & m, const std::complex<T> val) {
   zsMat<T> tmp = m; tmp += val; return tmp; }
template <class T> zsMat<T> operator + (const std::complex<T> val, const zsMat<T> & m) {
   zsMat<T> tmp = m; tmp += val; return tmp; }

template <class T> zsMat<T> zsMat<T>::operator -= (const zsMat<T> &m)   // Subtracts another matrix from current
{
   if(_n!=m._n || _m!=m._m) return *this;
   if(_iscsc)
   {
      if(m._iscsc) {
         if(is_subset(m)) return _cs_add_pat(m, true); else return _cs_adds(m, true); }
      else { 
         zsMat<T> n(m); n._cs_compress(); n._cs_dupl(); return _cs_adds(n, true); } 
   }
   else
   {  
      int l=_x.size();
      _p.insert(_p.end(),m._p.begin(),m._p.end());
      _i.insert(_i.end(),m._i.begin(),m._i.end());
      _x.insert(_x.end(),m._x.begin(),m._x.end()); 
      for(; l<_x.size(); l++) { _x[l] = -_x[l]; }
      return *this;
   }
}
template <class T> zsMat<T> zsMat<T>::operator -  (const zsMat<T> &m) { // Binary subtraction of two sparse matrices 
   if(m._iscsc) return _cs_add(m, true); else { zsMat<T> n(m); n._cs_compress(); n._cs_dupl(); return _cs_add(n, true); } }

template <class T> zsMat<T> zsMat<T>::operator -= (const T val)         // Subtract a constant from current matrix
{
   return add_scal(std::complex<T>(-val,0));
}
template <class T> zsMat<T> zsMat<T>::operator -= (const std::complex<T> val)
{
   return add_scal(-val);
}
template <class T> zsMat<T> operator - (const zsMat<T> & m, const T val) {
   zsMat<T> tmp = m; tmp -= val; return tmp; }
template <class T> zsMat<T> operator - (const T val, const zsMat<T> & m) {
   zsMat<T> tmp = m; tmp -= val; return tmp; }

template <class T> zsMat<T> zsMat<T>::operator *= (const T & c)         // Matrix scalar multiplication
{
   for(int i=0; i<_x.size(); i++) _x[i] *= c;
   return *this;
}
template <class T> zsMat<T> zsMat<T>::operator *= (const std::complex<T> & c)
{
   for(int i=0; i<_x.size(); i++) _x[i] *= c;
   return *this;
}
template <class T> zsMat<T> operator * (const zsMat<T> & m1, const T & c) {
   zsMat<T> tmp = m1; tmp *= c; return tmp; }                           // Binary matrix scalar multiplication
template <class T> zsMat<T> operator * (const T & c, const zsMat<T> & m1) {
   zsMat<T> tmp = m1; tmp *= c; return tmp; }
template <class T> zsMat<T> operator * (const zsMat<T> & m1, const std::complex<T> & c) {
   zsMat<T> tmp = m1; tmp *= c; return tmp; }                           // Binary matrix scalar multiplication
template <class T> zsMat<T> operator * (const std::complex<T> & c, const zsMat<T> & m1) {
   zsMat<T> tmp = m1; tmp *= c; return tmp; }

template <class T> zsMat<T> zsMat<T>::operator /= (const T & c)         // Matrix scalar division
{
   for(int i=0; i<_x.size(); i++) _x[i] /= c;
   return *this;
}
template <class T> zsMat<T> zsMat<T>::operator /= (const std::complex<T> & c)
{
   for(int i=0; i<_x.size(); i++) _x[i] /= c;
   return *this;
}
template <class T> zsMat<T> operator / (const zsMat<T> & m1, const T & c) {
   zsMat<T> tmp = m1; tmp /= c; return tmp; }                           // Binary matrix scalar division
template <class T> zsMat<T> operator / (const zsMat<T> & m1, const std::complex<T> & c) {
   zsMat<T> tmp = m1; tmp /= c; return tmp; }                           // Binary matrix scalar division

template <class T> zsMat<T> zsMat<T>::operator *= (const zsMat<T> &m) { // Matrix multiplication
   if(m._iscsc) return _cs_multiply(m); else { zsMat<T> n(m); n._cs_compress(); n._cs_dupl(); return _cs_multiply(n); } }
template <class T> zsMat<T> zsMat<T>::operator *  (const zsMat<T> &m) { // Binary matrix multiplication
   if(m._iscsc) return _cs_multiply(m); else { zsMat<T> n(m); n._cs_compress(); n._cs_dupl(); return _cs_multiply(n); } }
template <class T> zsMat<T> zsMat<T>::operator *= (const std::vector<T> & v)
{
   zsMat<T> w(1,v.size());                                              // Binary matrix.vector multiplication
   w._p.assign(2,0); w._p[1]=v.size();
   w._i.assign(v.size(),0); for(int j=0; j<_n; j++) w._i[j]=j;
   w._x.assign(v.size(),0);
   if(_iscsc) {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++) {
         w._x[_i[c]] += _x[c] * v[j]; } _genhash(); }
   else {
      for(int c=0; c<_x.size(); c++) w._x[_i[c]] += _x[c] * v[_p[c]]; }
   return w;
}
template <class T> zsMat<T> zsMat<T>::operator *= (const std::vector<std::complex<T> > & v)
{
   zsMat<T> w(1,v.size());                                              // Binary matrix.vector multiplication
   w._p.assign(2,0); w._p[1]=v.size();
   w._i.assign(v.size(),0); for(int j=0; j<_n; j++) w._i[j]=j;
   w._x.assign(v.size(),0);
   if(_iscsc) {
      for(int j=0; j<_n; j++) for(int c=_p[j]; c<_p[j+1]; c++) {
         w._x[_i[c]] += _x[c] * v[j]; } _genhash(); }
   else {
      for(int c=0; c<_x.size(); c++) w._x[_i[c]] += _x[c] * v[_p[c]]; }
   return w;
}

// --------------------------------------------------------------------------------------------------------------- //
// Friend functions
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::ostream & operator << (std::ostream & o, const zsMat<T> & m)
{
   if(m._iscsc)
   {
      for(int j=0; j<m._n; j++) for(int i=m._p[j]; i<m._p[j+1]; i++) {
         o << "x(" << (m._i[i]+1) << "," << (j+1) << ")\t=\t" << real(m._x[i]) << "+i*" << imag(m._x[i]) << ";\n"; }
   }
   else
   {
      for(int i=0; i<m._x.size(); i++) {
         o << "x(" << (m._i[i]+1) << "," << (m._p[i]+1) << ")\t=\t" << real(m._x[i]) << "+i*" << imag(m._x[i]) << ";\n"; }
   }
   return o;
}

// --------------------------------------------------------------------------------------------------------------- //
// Declaration for functions in feast.cpp
// --------------------------------------------------------------------------------------------------------------- //
//void feast_zheev(zsMat<double> &H, Vector &en, Matrix &zr, Matrix &zc);

#endif
