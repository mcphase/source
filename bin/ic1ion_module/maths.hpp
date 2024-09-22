/* maths.hpp
 *
 * This file implements mathematical and matrix classes with similar syntax to Matlab, and may be used to
 * translate programs from Matlab to C++
 *
 * The template classes sMat<>, and eigVE<> declared here. Declaration of the classes and nonmember functions
 * in this file is found at the start of the file. Declarations for files in maths.cpp is at the end
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 *
 */

#ifndef MATHS_H
#define MATHS_H

#include<cstdlib>
#include<cmath>
#include<vector>
#include<iostream>
#include<sstream>
#include<iomanip>
#include<string>
#include<map>
#include<cfloat>       // For definition of EPSILON etc.
#include "lapack.h"

#define PI 3.1415926535897932384626433832795

// --------------------------------------------------------------------------------------------------------------- //
// Template Class to hold a sparse matrix of any type - and also declares standard matrix algebra methods.
//    The class stores elements by coordinate indices in a subclass, _ind, consisting of two ints, which serves
//    as the key to a C++ STL map.
//    The member functions and overloaded operators has been made to be as similar to Matlab syntax as possible
// --------------------------------------------------------------------------------------------------------------- //
template <class T> class sMat;                                          // google  "Explicit Template Specification"
template <class T> std::ostream & operator << (std::ostream & o, const sMat<T> & m);
template <class T> std::istream & operator >> (std::istream & i, sMat<T> & m);
template <class T> T atoT(const std::string &s);
// --------------------------------------------------------------------------------------------------------------- //
template <class T> class sMat {
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
     std::map<_ind,T> _ls;                                              // List of elements

   public:
     // Constructors and Destructors //
     sMat(size_t r=1, size_t c=1): _r(r), _c(c) {};                     // constructs an empty r x c sparse matrix
     ~sMat() {};
    
     // Other member functions //
     int nnz() { return (int) _ls.size(); };                            // Returns the number of nonzero elements
     int nr() const { return _r; };
     int nc() const { return _c; };
     std::vector<int> size() const;                                     // Returns the size of the matrix: _r and _c
     std::vector< std::vector<int> > find() const;                      // Returns the indices of the nonzero elements
     std::vector< std::vector<int> > findupper() const;                 // Returns the indices of upper triangle elements
     std::vector< std::vector<int> > findlower() const;                 // Returns the indices of lower triangle elements
     bool issymm() const;                                               // Determines if matrix is symmetric
     bool issquare() const { return _r==_c; };                          // Determines if matrix is square
     bool isempty() const { return _ls.empty(); }                       // Determines if matrix is empty
     sMat<T> transpose();                                               // Transposes the matrix
     void mset(const char * m);                                         // Constructs a matrix using Matlab-like syntax
     void pset(int rs, int rn, int cs, int cn, const sMat<T> &m);       // Partially sets elements: x(rs:rn,cs:cn) = m
     sMat<T> setp(int rs, int rn, int cs, int cn);                      // Returns a slice x(rs:rn,cs:cn) of a matrix
     void del(int r, int c) { _ls.erase(_ind(r,c)); };                  // Deletes an element
     void clear() { _ls.clear(); _r = 0; _c = 0; };                     // Clears the matrix
     void zero() { _ls.clear(); };                                      // Zeros all entries
     void zero(int r, int c) { _ls.clear(); _r = r; _c = c; };          // Zeros all entries, changes shape as well
     void resize(int rs, int rn, int cs, int cn);                       // Resizes a matrix and deletes extra elements
     void reshape(int r, int c) { _r = r; _c = c; };                    // Changes shape without changing entries!
     std::string display_full() const;                                  // Outputs a string of the full matrix
     std::string display_nonzero() const;                               // Outputs a string of nonzero elements of the  matrix
     T *f_array() const;                                                // Returns matrix as a Fortran style 2D array
   //std::map<_ind,T> ls() { return _ls; }
     void mcol(int c, T val);                                           // Multiplies a given column by a constant
     void mrow(int r, T val);                                           // Multiplies a given row by a constant

     // Overloaded operators
     sMat<T> operator =  (const sMat & m);                              // Copy assignment - overwrites previous matrix
     sMat<T> operator += (const sMat & m);                              // Add another matrix to current (element-wise)
     sMat<T> operator += (const T val);                                 // Add a constant value to current (element-wise)
     sMat<T> operator -= (const sMat & m);                              // Subtract another matrix from current
     sMat<T> operator -= (const T val);                                 // Subtract a constant value from current 
     sMat<T> operator *= (const T & c);                                 // Matrix scalar multiplication
     sMat<T> operator /= (const T & c);                                 // Matrix scalar division
     sMat<T> operator *= (const sMat & m);                              // Matrix multiplication
     sMat<T> operator *= (const std::vector<T> & v);                    // Matrix*Vector (Matrix operating on a vector)
     sMat<T> operator ^= (const sMat & m);                              // Element-wise matrix multiplication
     sMat<T> operator %= (const sMat & m);                              // Element-wise matrix division
     T& operator () (int r, int c);                                     // Gets/sets element by subscript index
     T  operator () (int r, int c) const;

     // Friend function  to provide input/output via <iostream>
     friend std::ostream & operator << <> (std::ostream & o, const sMat & m);
     friend std::istream & operator >> <> (std::istream & i, sMat & m);

     // Function required by ARPACK++
     void MultMv(T *v, T *w);                                           // Calculates the matrix-vector product w = M*v

};  // End of template <class T> class sMat

/* --------------------------------------------------------------------------------------------------------------- //
   Template functions declared in this file
// --------------------------------------------------------------------------------------------------------------- //
   class eigVE { std::vector<T> E; sMat<T> V; };                        // Class hold eigenval/vectors of a matrix
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
template <class T> std::vector<int> sMat<T>::size() const
{
   std::vector<int> retval;
   retval.push_back(_r);
   retval.push_back(_c);
   return retval;
}

template <class T> std::vector< std::vector<int> > sMat<T>::find() const
{
   int r = 0;
   typename std::map<_ind,T>::iterator i;                               // Help out compiler by telling it should by
   std::vector<int> row(2);                                             //   looking for a type, otherwise get error
   std::vector< std::vector<int> > retval(_ls.size(),row);              //   on g++.
   std::map<_ind,T> tmp_ls = _ls;

   for (i=tmp_ls.begin(); i!=tmp_ls.end(); i++)
   {
      row[0] = i->first.r; row[1] = i->first.c;
      retval[r] = row; r++;
   }

   return retval;
}

template <class T> std::vector< std::vector<int> > sMat<T>::findlower() const
{
   int r,n = 0;
   typename std::map<_ind,T>::iterator i;
   std::vector<int> row(2);
   std::vector< std::vector<int> > retval(_ls.size()/2+_r,row);
   std::map<_ind,T> tmp_ls = _ls;

   for (r=0; r<_r; r++)
      for (i=tmp_ls.lower_bound(_ind(r,0)); i!=tmp_ls.lower_bound(_ind(r,r)); i++)
      {
         row[0] = i->first.r; row[1] = i->first.c;
         retval[n++] = row;
      }

   retval.erase(retval.begin()+n,retval.end());
   return retval;
}

template <class T> std::vector< std::vector<int> > sMat<T>::findupper() const
{
   int r,n = 0;
   typename std::map<_ind,T>::iterator i;
   std::vector<int> row(2);
   std::vector< std::vector<int> > retval(_ls.size()/2+_r,row);
   std::map<_ind,T> tmp_ls = _ls;

   for (r=0; r<_r; r++)
      for (i=tmp_ls.lower_bound(_ind(r,r)); i!=tmp_ls.lower_bound(_ind(r,_c)); i++)
      {
         row[0] = i->first.r; row[1] = i->first.c;
         retval[n++] = row;
      }

   retval.erase(retval.begin()+n,retval.end());
   return retval;
}

template <class T> bool sMat<T>::issymm() const
{
   typename std::map<_ind,T>::iterator it;
   int i,j;
   std::map<_ind,T> tmp_ls = _ls;

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

template <class T> sMat<T> sMat<T>::transpose()
{
   typename std::map<_ind,T>::iterator it;
   int i,j;
   std::map<_ind,T> tmp_ls;

   for(it=_ls.begin(); it!=_ls.end(); it++)
   {
      i = (*it).first.r; j = (*it).first.c;
      tmp_ls[_ind(j,i)] = (*it).second;
   }
   _ls = tmp_ls;
   i = _c; _c = _r; _r = i;
   return *this;
}

template <class T> void sMat<T>::mset(const char *m)                    // Constructs a matrix using Matlab syntax:
{                                                                       //   e.g: m.mset("1 0 0;0 1 0;0 0 1");
   std::string m_str(m);
   std::string entry;
   int c=0,r=0;                                                         // Column and Row indices
   size_t wspc,stmk,enmk,mathop,sqrtop,semicol;                         // Indices in constructor string
   //double el,el2;
   T el,el2;

   _ls.clear(); _c = 1; _r = 1;
   wspc = 0;
   while(wspc!=std::string::npos)
   {
      stmk = m_str.find_first_of(".0123456789",wspc);                   // Each entry in the matrix is delimited by
      if(stmk!=std::string::npos)                                       //   a whitespace character. The loop finds
      {                                                                 //   all instances of white spaces and
         wspc = m_str.find_first_of(" \t;",stmk+1);                     //   parses the substrings formed between.
         enmk = m_str.find_last_of(".0123456789",wspc);                 //   -> Numerical strings become doubles
         entry = m_str.substr(stmk,enmk-stmk+1);                        //   -> Maths operators (+,-,*,/,sqrt)
         mathop = entry.find_first_of("+-*/");                          //      and the numerical strings they 
         if(mathop!=std::string::npos)                                  //      delimit are parsed separately. Then
         {                                                              //      these are recombined to make an entry
            mathop += stmk;                                             //      in our sparse matrix, sMat.
            entry = m_str.substr(stmk,mathop-stmk+1);
            el = atoT<T>(entry);
            entry = m_str.substr(mathop+1,enmk-stmk+1);
            el2 = atoT<T>(entry);
            switch(m_str[mathop])
            {
               case '+': el += el2; break;
               case '-': el -= el2; break;
               case '*': el *= el2; break;
               case '/': el /= el2; break;
            }
         }
         else
         {
            el = atoT<T>(entry);
         }

         stmk = m_str.find_last_of(" \t",wspc-1); 
         if (stmk==std::string::npos)
            stmk = 0;
         entry = m_str.substr(stmk,wspc-stmk+1);
         sqrtop = entry.find("sqrt(");
         if(sqrtop!=std::string::npos)
         {
            if(m_str.find_first_of(")",enmk,wspc-enmk)!=std::string::npos)
            {
               el = sqrt(el);
            }
         }
      }                                                                 // End of parsing.

      if(el!=0.0)
      {
         _ls[_ind(r,c)] = el;
      }
      
      if(wspc>(enmk+1))
      {
         if(!m_str.compare(enmk+1,1,";"))                               // Each row of the matrix is delimited by
         {                                                              //   a semicolon.         
            r++;
            c=0; 
         } 
         else if(wspc!=std::string::npos)                               // To get the "c" count correct, otherwise
         {                                                              //   c is too large by 1 when we exit, so
            c++;                                                        //   _c is also too large by 1.
         }
      }
      else if(!m_str.compare(wspc+1,1,";")) 
      { 
         r++; 
         c=0; 
      } 
      else if(wspc!=std::string::npos) 
      { 
         c++; 
      }
   }                                                                    // End of while loop
   if(++r>_r)
      _r = r;                                                           // Updates the matrix dimensions if needed
   if(++c>_c)
      _c = c;
}

template <class T> void sMat<T>::pset(int rs, int rn, int cs, int cn, const sMat<T> &m)
{
   typename std::map<_ind,T>::iterator it_s,it_n;
   typename std::map<_ind,T>::const_iterator it;

   int r = m.nr(); int c = m.nc();
   if(rs==-1 && rn==-1) { rs = 1; rn = (r>_r)?r:_r; } else if(rn==-1) rn = rs+r-1;
   if(cs==-1 && cn==-1) { cs = 1; cn = (c>_c)?c:_c; } else if(cn==-1) cn = cs+c-1;
   if(rn>_r) _r = rn;
   if(cn>_c) _c = cn;

   rs--; cs--;
   if(!m._ls.empty())
   {
      // Sets elements with rs:rn,cs:cn to zero if any are nonzero
      it_s = _ls.lower_bound(_ind(rs,cs));
      it_n = _ls.lower_bound(_ind(rn,cn));
      while(it_s->first < it_n->first)
      {
         if(it_s->first.r<rn && it_s->first.c<cn && it_s->first.r>=rs && it_s->first.c>=cs)
            _ls.erase(it_s);
         it_s++;
      }

      // Sets new elements rs:rc,cs:cn from the values of the matrix m
      for (it=m._ls.begin(); it!=m._ls.end(); it++)
         _ls[_ind(it->first.r+rs,it->first.c+cs)] = it->second;
   }
}

template <class T> sMat<T> sMat<T>::setp(int rs, int rn, int cs, int cn)// Returns a slice x(rs:rn,cs:cn) of a matrix
{
   int i;
   typename std::map<_ind,T>::iterator it;
   
   if(rs==-1 && rn==-1) { rs = 1; rn = _r; }
   if(cs==-1 && cn==-1) { cs = 1; cn = _c; }

   sMat<T> retval(rn-rs+1,cn-cs+1);

   rs--; cs--; //rn--; cn--;
   if((int)_ls.size() != 0)
   {
      for(i=rs; i<rn; i++)
      {  
         for(it=_ls.lower_bound(_ind(i,cs)); it!=_ls.lower_bound(_ind(i,cn)); it++)
         {
            retval(it->first.r-rs,it->first.c-cs) = it->second;
         }
      }  
   }
   return retval;
}

template <class T> void sMat<T>::resize(int rs, int rn, int cs, int cn) // Resizes a matrix and deletes extra elements
{
   int i;
   typename std::map<_ind,T>::iterator it;

   if(rs==-1) rs = 1; if(rn==-1) rn = _r;
   if(cs==-1) cs = 1; if(cn==-1) cn = _r;
   
   rs--; cs--; //rn--; cn--;
   _ls.erase(_ls.upper_bound(_ind(rn,cn)),_ls.end());
   if(rs==0 && cs==0)
   {
      for(i=0; i<rn; i++)
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
      _r = (rn-rs); _c = (cn-cs);
   }
}

template <class T> std::string sMat<T>::display_full() const            // Prints the full matrix to a string.
{
   int r,c,mr,mc;
   std::stringstream retval;

   retval << "[";
   for (r=0; r<(_r-1); r++)
   {
      for (c=0; c<(_c-1); c++)
      {
         if(_ls.find(_ind(r,c))==_ls.end())
            retval << "0 \t";
         else
         {
            mr = _ls.find(_ind(r,c))->first.r;
            mc = _ls.find(_ind(r,c))->first.c;
            if( (mr>_r || mr<0) || (mc>_c || mc<0) )
               retval << "0 \t";
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
         if( (mr>_r || mr<0) || (mc>_c || mc<0) )
            retval << "0;\n";
         else
            retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << ";\n";
      }
   }
   for (c=0; c<(_c-1); c++)
   {
      if(_ls.find(_ind(r,c))==_ls.end())
         retval << "0 \t";
      else
      {
         mr = _ls.find(_ind(r,c))->first.r;
         mc = _ls.find(_ind(r,c))->first.c;
         if( (mr>_r || mr<0) || (mc>_c || mc<0) )
            retval << "0 \t";
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
      if( (mr>_r || mr<0) || (mc>_c || mc<0) )
         retval << "0];\n";
      else
         retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << "];\n";
   }
//  retval << _r;  retval << " x ";  retval << _c;  retval << "\n";

   return retval.str();
}

template <class T> std::string sMat<T>::display_nonzero() const            // Prints the full matrix to a string.
{
   int r,c,mr,mc;
   std::stringstream retval;

   retval << "[";
   for (r=0; r<(_r-1); r++)
   {
      for (c=0; c<(_c-1); c++)
      {
         if(_ls.find(_ind(r,c))!=_ls.end())
         {
            mr = _ls.find(_ind(r,c))->first.r;
            mc = _ls.find(_ind(r,c))->first.c;
            if( ((mr<=_r || mr>=0) || (mc<=_c || mc>=0))&& fabs(_ls.find(_ind(r,c))->second)>1e-15)
               retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << "\t";
         }
      }
      if(_ls.find(_ind(r,c))!=_ls.end())
      {
         mr = _ls.find(_ind(r,c))->first.r;
         mc = _ls.find(_ind(r,c))->first.c;
         if( ((mr<=_r || mr>=0) || (mc<=_c || mc>=0)) && fabs(_ls.find(_ind(r,c))->second)>1e-15)
             retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << ";\n";
      }
   }
   for (c=0; c<(_c-1); c++)
   {
      if(_ls.find(_ind(r,c))!=_ls.end())
       {
         mr = _ls.find(_ind(r,c))->first.r;
         mc = _ls.find(_ind(r,c))->first.c;
         if( ((mr<=_r || mr>=0) || (mc<=_c || mc>=0))&& fabs(_ls.find(_ind(r,c))->second)>1e-15)
            retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << "\t";
      }
   }
   if(_ls.find(_ind(r,c))!=_ls.end())
       {
      mr = _ls.find(_ind(r,c))->first.r;
      mc = _ls.find(_ind(r,c))->first.c;
      if( ((mr<=_r || mr>=0) || (mc<=_c || mc>=0))&& fabs(_ls.find(_ind(r,c))->second)>1e-15)
                retval << std::setprecision(16) << _ls.find(_ind(r,c))->second << "];\n";
   }
 retval << _r;  retval << " x ";  retval << _c;  retval << "\n";

   return retval.str();
}

template <class T> T* sMat<T>::f_array() const                  // Returns matrix as a Fortran style 2D array
{
   // Fortran 2D arrays are column-major dense contiguous blocks of memory, unlike a 2D C-array which is a 1D array
   // of pointers to other 1D arrays. So effectively what we output is a pointer to a N*M length 1D C-array.
   T *retval;
   typename std::map<_ind,T>::iterator i;
   std::map<_ind,T> tmp_ls = _ls;

   retval = (T*) calloc(_r*_c,sizeof(T));  // Allocates an _r*_c array and initiallises all elements to zero.

   for (i=tmp_ls.begin(); i!=tmp_ls.end(); i++)
      retval[_r*i->first.c+i->first.r] = i->second;

   return retval;
}

template <class T> void sMat<T>::mcol(int c, T val)             // Multiplies a particular column by a constant
{
   int r;
   typename std::map<_ind,T>::iterator it;
   for(r=0; r<_r; r++)
   {
      it = _ls.find(_ind(r,c));
      if(it!=_ls.end())
         it->second *= val;
   }
}
template <class T> void sMat<T>::mrow(int r, T val)             // Multiplies a particular column by a constant
{
   typename std::map<_ind,T>::iterator lb,ub,it;
   lb = _ls.lower_bound(_ind(r,0));
   ub = _ls.upper_bound(_ind(r,_c));
   if (lb->first < ub->first)
      for(it=lb; it->first<ub->first; it++)
         it->second *= val;
}

// --------------------------------------------------------------------------------------------------------------- //
// Overloaded operators
// --------------------------------------------------------------------------------------------------------------- //
template <class T> T& sMat<T>::operator () (int r, int c)
{ 
   if(r>_r) _r = r; if(c>_c) _c = c;
   return _ls[_ind(r,c)]; 
}
template <class T> T sMat<T>::operator () (int r, int c) const
{
   int mr,mc;
   if(_ls.find(_ind(r,c))==_ls.end())
      return 0;
   else
   {
      mr = _ls.find(_ind(r,c))->first.r;
      mc = _ls.find(_ind(r,c))->first.c;
      if( (mr>_r || mr<0) || (mc>_c || mc<0) ) 
         return 0;
      else 
         return _ls.find(_ind(r,c))->second;
   }
}

template <class T> sMat<T> sMat<T>::operator = (const sMat<T> & m)
{
   _r = m.nr(); _c = m.nc();
   _ls = m._ls;
   return *this;
}

template <class T> sMat<T> sMat<T>::operator += (const sMat<T> & m)     // Adds another matrix to current
{
   typename std::map<_ind,T>::const_iterator it;
   if(!m._ls.empty())
      for(it=m._ls.begin(); it!=m._ls.end(); it++)
         _ls[_ind(it->first.r,it->first.c)] += it->second;

   return *this;
}
template <class T> sMat<T> operator + (const sMat<T> & m1, const sMat<T> & m2)
{
   sMat<T> tmp = m1;                                                    // Binary addition of two sparse matrices
   tmp += m2;
   return tmp;
}

template <class T> sMat<T> sMat<T>::operator += (const T val)           // Adds a constant to the current matrix
{
   typename std::map<_ind,T>::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second += val;
   
   return *this;
}
template <class T> sMat<T> operator + (const sMat<T> & m, const T val)
{
   sMat<T> tmp = m;
   tmp += val;
   return tmp;
}

template <class T> sMat<T> sMat<T>::operator -= (const sMat<T> & m)     // Subtracts another matrix from current
{
   typename std::map<_ind,T>::const_iterator it;
   if(!m._ls.empty())
      for(it=m._ls.begin(); it!=m._ls.end(); it++)
      {
         _ls[_ind(it->first.r,it->first.c)] -= it->second;
         if ((_ls.find(_ind(it->first.r,it->first.c))->second) == 0.0)
            _ls.erase(_ind(it->first.r,it->first.c));
      }

   return *this;
}
template <class T> sMat<T> operator - (const sMat<T> & m1, const sMat<T> & m2)
{
   sMat<T> tmp = m1;                                                    // Binary subtraction of two sparse matrices
   tmp -= m2;
   return tmp;
}

template <class T> sMat<T> sMat<T>::operator -= (const T val)           // Subtract a constant from current matrix
{
   typename std::map<_ind,T>::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second -= val;
   
   return *this;
}
template <class T> sMat<T> operator - (const sMat<T> & m, const T val)
{
   sMat<T> tmp = m;
   tmp -= val;
   return tmp;
}

template <class T> sMat<T> sMat<T>::operator *= (const T & c)           // Matrix scalar multiplication
{
   typename std::map<_ind,T>::iterator i;

   for (i=_ls.begin(); i!=_ls.end(); i++)
   {
      i->second *= c;
   }
   return *this;
}
template <class T> sMat<T> operator * (const sMat<T> & m1, const T & c)
{
   sMat<T> tmp = m1;                                                    // Binary matrix scalar multiplication
   tmp *= c;
   return tmp;
}
template <class T> sMat<T> sMat<T>::operator /= (const T & c)           // Matrix scalar division
{
   typename std::map<_ind,T>::iterator i;

   for (i=_ls.begin(); i!=_ls.end(); i++)
   {
      i->second /= c;
   }
   return *this;
}
template <class T> sMat<T> operator / (const sMat<T> & m1, const T & c)
{
   sMat<T> tmp = m1;                                                    // Binary matrix scalar division
   tmp /= c;
   return tmp;
}

template <class T> sMat<T> sMat<T>::operator *= (const sMat<T> & m)     // Matrix multiplication
{
   sMat<T> tmp;

   if(_c != m._r)
      std::cerr << "Error: sMat<T>::operator *= (const sMat<T>): Matrix sizes are incommensurate!\n";
   else
   {
      int r,c;
      typename std::map<_ind,T>::iterator lb,ub,it;
      T elem;

      // Using the traditional method, except ignoring zeros in the sparse structure
      tmp._r = _r;
      tmp._c = m._c;
      for(r=0; r<_r; r++)
      {
         lb = _ls.lower_bound(_ind(r,0));
         ub = _ls.upper_bound(_ind(r,_c));
         if(lb->first < ub->first)
            for(c=0; c<m._c; c++)
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
template <class T> sMat<T> operator * (const sMat<T> & m1, const sMat<T> & m2)
{
   sMat<T> tmp = m1;                                                    // Binary matrix multiplication
   tmp *= m2;
   return tmp;
}
template <class T> sMat<T> sMat<T>::operator *= (const std::vector<T> & v)
{
   sMat<T> tmp;                                                         // Binary matrix.vector multiplication
   int i;
   typename std::map<_ind,T>::iterator lb,ub,it;
   T elem;

   if(_c != (int)v.size())
      std::cerr << "Error: sMat<T>::operator * (sMat<T>, std::vector<T>): Matrix and vector sizes incommensurate!\n";
   else
   {
      tmp._r = _r; tmp._c = 1;
      for(i=0; i<_c; i++)
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
template <class T> void sMat<T>::MultMv(T *v, T *w)                     // Calculates the matrix-vector product w = M*v
{                                                                       // Needed by ARPACK++
   int i;
   typename std::map<_ind,T>::iterator it;
   // We have to assume that the size of the vector v is equal to _c
   for(i=0; i<_c; i++)
   {
      w[i] = 0.0;
      for(it=_ls.lower_bound(_ind(i,0)); it!=_ls.lower_bound(_ind(i,_c)); it++)
         w[i] += (it->second * v[it->first.c]);
   }
}
template <class T> std::vector<T> operator * (const sMat<T> & m1, const std::vector<T> & v)
{
   std::vector<T> tmp;                                                  // Binary matrix.vector multiplication
   sMat<T> tmpM = m1;
   int i;

   tmpM *= v; 
   for(i=0; i<tmpM.nr(); i++)
      tmp.push_back(tmpM(i,0));

   return tmp;
}

template <class T> sMat<T> sMat<T>::operator ^= (const sMat<T> &m)      // Element-wise matrix multiplication 
{
   if(m.nr()!=_r || m.nc()!=_c) { std::cerr << "sMat ^ operator: M1 and M2 not same size\n"; return *this; }
   typename std::map<_ind,T>::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second *= m(it->first.r,it->first.c);
   return *this;
}
template <class T> sMat<T> operator ^ (const sMat<T> & m1, const sMat<T> & m2)
{
   sMat<T> tmp = m1;
   tmp ^= m2;
   return tmp;
}

template <class T> sMat<T> sMat<T>::operator %= (const sMat<T> &m)      // Element-wise matrix division
{
   if(m.nr()!=_r || m.nc()!=_c) { std::cerr << "sMat % operator: M1 and M2 not same size\n"; return *this; }
   typename std::map<_ind,T>::iterator it;
   for (it=_ls.begin(); it!=_ls.end(); it++)
      (*it).second /= m(it->first.r,it->first.c);
   return *this;
}
template <class T> sMat<T> operator % (const sMat<T> & m1, const sMat<T> & m2)
{
   sMat<T> tmp = m1;
   tmp %= m2;
   return tmp;
}

// --------------------------------------------------------------------------------------------------------------- //
// Friend functions
// --------------------------------------------------------------------------------------------------------------- //
template <class T> std::ostream & operator << (std::ostream & o, const sMat<T> & m)
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

template <class T> std::istream & operator >> (std::istream & i, sMat<T> & m)
{
   int r,c;
   T x;
   for (r=0; r<m._r; r++)
      for (c=0; c<m._c; c++)
      {
         i >> x;
         if (x!=0.0)
            m(r,c) = x;
      }

   return i;
}

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
// Class hold eigenvalues and eigenvectors of a diagonalised matrix
// --------------------------------------------------------------------------------------------------------------- //
template <class T> class eigVE {
   public:
   std::vector<T> E;
   sMat<T> V;
};

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
