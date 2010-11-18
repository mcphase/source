/*-----------------------------------------------------------------------------*\
| vector and matrix classes private include file                       vecinl.h |
|                                                                               |
| MatPack Library Release 1.0                                                   |
| Copyright (C) 1991-1996 by Berndt M. Gammel                                   |
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

// avoid multiple inclusion
#ifndef _VECINL_H_
#define _VECINL_H_

#include <complex>
#include <cstring>

//----------------------------------------------------------------------------//
// This include file contains high speed inlined loops for elementary 
// vector operations, such as copying, adding, multiplying, etc..
// Don't change these inline functions before you havn't run extensive tests.
//
// The routines below seem to be the fastest versions at least with the 
//    DEC cxx 1.3  compiler on a DEC Alpha
//    with flags -O2 
// or when using the
//    gcc 2.7.2 compiler on a DEC Alpha and Linux Pentim PC
//    with flags -O3 -finline-functions -fomit-frame-pointer -fforce-addr
//
// Notes:
// 
//   (1) inline void addvec (double* dst, const double* src, int n)
//         { while (n--) *dst++ += *src++; }
//   (2) inline void addvec (double* dst, const double* src, int n)
//         { for (int i = 0; i < n; i++) dst[i] += src[i]; }
//
//   Version (2) is 1.6 times as fast as the pointer arithmetics version (1)
//   with and without optimization. Optimization speeds up both by 1.3.
//
//   (3) inline void mulval1 (double* dst, double value, int n)
//         { while (n--) *dst++ *= value; }
//   (4) inline void mulval2 (double* dst, double value, int n)
//         { for (int i = 0; i < n; i++) dst[i] *= value; }
//
//   In the case of combining vector with constant there is no difference
//   between pointer version (3) and explicit version (4). But if optimization
//   is switched off then (3) is twice as fast as (4).
//
//   Copying should always be done with memcpy(), which is between 2 and
//   10 times as fast as any loop version (depending on the vector dimension).
//
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
// float - single precision
//----------------------------------------------------------------------------//

inline void copyval (float* dst, float value, int n)
{
    while (n--) *dst++ = value;
}

inline void addval (float* dst, float value, int n)
{
    while (n--) *dst++ += value;
}

inline void addval (float* dst, const float* src, float value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] + value;
}

inline void subval (float* dst, float value, int n)
{
    while (n--) *dst++ -= value;
}

inline void subval (float* dst, const float* src, float value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] - value;
}

inline void subval (float* dst, float value, const float* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value - src[i];
}

inline void mulval (float* dst, float value, int n)
{
    while (n--) *dst++ *= value;
}

inline void mulval (float* dst, const float* src, float value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] * value;
}

inline void divval (float* dst, float value, int n)
{
    while (n--) *dst++ /= value;
}

inline void divval (float* dst, const float* src, float value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] / value;
}

inline void divval (float value, float* dst, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value / dst[i];
}

inline void divval (float* dst, float value, const float* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value / src[i];
}

inline int cmpval (const float* dst, float value, int n)
{
    while (n--) if (*dst++ != value) return 0;
    return 1;
}

//----------------------------------------------------------------------------//

inline float dotvec (const float* src1, const float* src2, int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++) sum += src1[i] * src2[i];
    return sum;
}

//----------------------------------------------------------------------------//

inline void copyvec (float* dst, const float* src, int n)
{
    memcpy((void*)dst,(const void *)src, n*sizeof(float));
}

//----------------------------------------------------------------------------//

inline void negvec (float* dst, const float* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = -src[i];
}

inline void addvec (float* dst, const float* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] += src[i];
}

inline void addvec (float* dst, const float* src1, const float* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] + src2[i];
}

inline void subvec (float* dst, const float* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] -= src[i];
}

inline void subvec (float* dst, const float* src1, const float* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] - src2[i];
}

inline void mulvec (float* dst, const float* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] *= src[i];
}

inline void mulvec (float* dst, const float* src1, const float* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] * src2[i];
}

inline int cmpvec (const float* dst, const float* src, int n)
{
    while (n--) if (*dst++ != *src++) return 0;
    return 1;
}

inline float norm2 (const float* src, int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++) sum += src[i] * src[i];
    return sum;
}

//----------------------------------------------------------------------------//
// double - double precision
//----------------------------------------------------------------------------//

inline void copyval (double* dst, double value, int n)
{
    while (n--) *dst++ = value;
}

inline void addval (double* dst, double value, int n)
{
    while (n--) *dst++ += value;
}

inline void addval (double* dst, const double* src, double value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] + value;
}

inline void subval (double* dst, double value, int n)
{
    while (n--) *dst++ -= value;
}

inline void subval (double* dst, const double* src, double value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] - value;
}

inline void subval (double* dst, double value, const double* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value - src[i];
}

inline void mulval (double* dst, double value, int n)
{
    while (n--) *dst++ *= value;
}

inline void mulval (double* dst, const double* src, double value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] * value;
}

inline void divval (double* dst, double value, int n)
{
    while (n--) *dst++ /= value;
}

inline void divval (double* dst, const double* src, double value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] / value;
}

inline void divval (double value, double* dst, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value / dst[i];
}

inline void divval (double* dst, double value, const double* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value / src[i];
}

inline int cmpval (const double* dst, double value, int n)
{
    while (n--) if (*dst++ != value) return 0;
    return 1;
}

//----------------------------------------------------------------------------//

inline double dotvec (const double* src1, const double* src2, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++) sum += src1[i] * src2[i];
    return sum;
}

//----------------------------------------------------------------------------//

inline void copyvec (double* dst, const double* src, int n)
{
    memcpy((void*)dst,(const void *)src, n*sizeof(double));
}

//----------------------------------------------------------------------------//

inline void negvec (double* dst, const double* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = -src[i];
}

inline void addvec (double* dst, const double* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] += src[i];
}

inline void addvec (double* dst, const double* src1, const double* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] + src2[i];
}

inline void subvec (double* dst, const double* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] -= src[i];
}

inline void subvec (double* dst, const double* src1, const double* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] - src2[i];
}

inline void mulvec (double* dst, const double* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] *= src[i];
}

inline void mulvec (double* dst, const double* src1, const double* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] * src2[i];
}

inline int cmpvec (const double* dst, const double* src, int n)
{
    while (n--) if (*dst++ != *src++) return 0;
    return 1;
}

inline double norm2 (const double* src, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++) sum += src[i] * src[i];
    return sum;
}

//----------------------------------------------------------------------------//
// complex<double> - complex double precision
//----------------------------------------------------------------------------//

inline void copyval (complex<double>* dst, const complex<double>& value, int n)
{
    while (n--) *dst++ = value;
}

inline void addval (complex<double>* dst, const complex<double>& value, int n)
{
    while (n--) *dst++ += value;
}

inline void addval (complex<double>* dst, const complex<double>* src, complex<double>& value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] + value;
}

inline void subval (complex<double>* dst, const complex<double>& value, int n)
{
    while (n--) *dst++ -= value;
}

inline void subval (complex<double>* dst, const complex<double>* src, complex<double>& value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] - value;
}

inline void subval (complex<double>* dst, const complex<double>& value, const complex<double>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value - src[i];
}

inline void mulval (complex<double>* dst, complex<double>& value, int n)
{
    while (n--) *dst++ *= value;
}

inline void mulval (complex<double>* dst, const complex<double>* src, const complex<double>& value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] * value;
}

inline void divval (complex<double>* dst, const complex<double>& value, int n)
{
    while (n--) *dst++ /= value;
}

inline void divval (complex<double>* dst, const complex<double>* src, const complex<double>& value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] / value;
}

inline void divval (const complex<double>& value, complex<double>* dst, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value / dst[i];
}

inline void divval (complex<double>* dst, const complex<double>& value, const complex<double>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value / src[i];
}
inline int cmpval (const complex<double>* dst, const complex<double>& value, int n)
{
    while (n--) if (*dst++ != value) return 0;
    return 1;
}

//----------------------------------------------------------------------------//

inline complex<double> matdotvec (const complex<double>* src1, const complex<double>* src2, int n)
{
    complex<double> sum = 0;
    for (int i = 0; i < n; i++) sum += src1[i] * src2[i];
    return sum;
}

inline complex<double> dotvec (const complex<double>* src1, const complex<double>* src2, int n)
{
    complex<double> sum = 0;
    for (int i = 0; i < n; i++) sum += src1[i] * conj(src2[i]);
    return sum;
}

//----------------------------------------------------------------------------//

inline void copyvec (complex<double>* dst, const complex<double>* src, int n)
{
    memcpy((void*)dst,(const void *)src, n*sizeof(complex<double>));
}

inline void copyrevec (complex<double>* dst, const double* re, int n)
{
    for (int i = 0; i < n; i++) dst[i] = re[i];
}

inline void copyreimvec (complex<double>* dst, const double* re, const double* im, int n)
{
    for (int i = 0; i < n; i++) dst[i] = complex<double>(re[i],im[i]);
}

//----------------------------------------------------------------------------//

inline void negvec (complex<double>* dst, const complex<double>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = -src[i];
}

inline void addvec (complex<double>* dst, const complex<double>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] += src[i];
}

inline void addvec (complex<double>* dst, const complex<double>* src1, const complex<double>* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] + src2[i];
}

inline void subvec (complex<double>* dst, const complex<double>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] -= src[i];
}

inline void subvec (complex<double>* dst, const complex<double>* src1, const complex<double>* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] - src2[i];
}

inline void mulvec (complex<double>* dst, const complex<double>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] *= src[i];
}

inline void mulvec (complex<double>* dst, const complex<double>* src1, const complex<double>* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] * src2[i];
}

inline int cmpvec (const complex<double>* dst, const complex<double>* src, int n)
{
    while (n--) if (*dst++ != *src++) return 0;
    return 1;
}

inline double norm2 (const complex<double>* src, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++) sum += norm(src[i]);
    return sum;
}

//----------------------------------------------------------------------------//
// complex<float> - complex single precision
//----------------------------------------------------------------------------//

inline void copyval (complex<float>* dst, const complex<float>& value, int n)
{
    while (n--) *dst++ = value;
}

inline void addval (complex<float>* dst, const complex<float>& value, int n)
{
    while (n--) *dst++ += value;
}

inline void addval (complex<float>* dst, const complex<float>* src, complex<float>& value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] + value;
}

inline void subval (complex<float>* dst, const complex<float>& value, int n)
{
    while (n--) *dst++ -= value;
}

inline void subval (complex<float>* dst, const complex<float>* src, complex<float>& value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] - value;
}

inline void subval (complex<float>* dst, const complex<float>& value, const complex<float>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value - src[i];
}

inline void mulval (complex<float>* dst, complex<float>& value, int n)
{
    while (n--) *dst++ *= value;
}

inline void mulval (complex<float>* dst, const complex<float>* src, const complex<float>& value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] * value;
}

inline void divval (complex<float>* dst, const complex<float>& value, int n)
{
    while (n--) *dst++ /= value;
}

inline void divval (complex<float>* dst, const complex<float>* src, const complex<float>& value, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src[i] / value;
}

inline void divval (const complex<float>& value, complex<float>* dst, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value / dst[i];
}

inline void divval (complex<float>* dst, const complex<float>& value, const complex<float>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = value / src[i];
}
inline int cmpval (const complex<float>* dst, const complex<float>& value, int n)
{
    while (n--) if (*dst++ != value) return 0;
    return 1;
}

//----------------------------------------------------------------------------//

inline complex<float> matdotvec (const complex<float>* src1, const complex<float>* src2, int n)
{
    complex<float> sum = 0;
    for (int i = 0; i < n; i++) sum += src1[i] * src2[i];
    return sum;
}

inline complex<float> dotvec (const complex<float>* src1, const complex<float>* src2, int n)
{
    complex<float> sum = 0;
    for (int i = 0; i < n; i++) sum += src1[i] * conj(src2[i]);
    return sum;
}

//----------------------------------------------------------------------------//

inline void copyvec (complex<float>* dst, const complex<float>* src, int n)
{
    memcpy((void*)dst,(const void *)src, n*sizeof(complex<float>));
}

inline void copyrevec (complex<float>* dst, const double* re, int n)
{
    for (int i = 0; i < n; i++) dst[i] = re[i];
}

inline void copyreimvec (complex<float>* dst, const double* re, const double* im, int n)
{
    for (int i = 0; i < n; i++) dst[i] = complex<float>(re[i],im[i]);
}

//----------------------------------------------------------------------------//

inline void negvec (complex<float>* dst, const complex<float>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] = -src[i];
}

inline void addvec (complex<float>* dst, const complex<float>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] += src[i];
}

inline void addvec (complex<float>* dst, const complex<float>* src1, const complex<float>* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] + src2[i];
}

inline void subvec (complex<float>* dst, const complex<float>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] -= src[i];
}

inline void subvec (complex<float>* dst, const complex<float>* src1, const complex<float>* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] - src2[i];
}

inline void mulvec (complex<float>* dst, const complex<float>* src, int n)
{
    for (int i = 0; i < n; i++) dst[i] *= src[i];
}

inline void mulvec (complex<float>* dst, const complex<float>* src1, const complex<float>* src2, int n)
{
    for (int i = 0; i < n; i++) dst[i] = src1[i] * src2[i];
}

inline int cmpvec (const complex<float>* dst, const complex<float>* src, int n)
{
    while (n--) if (*dst++ != *src++) return 0;
    return 1;
}

inline double norm2 (const complex<float>* src, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++) sum += norm(src[i]);
    return sum;
}

//----------------------------------------------------------------------------//
// integer 
//----------------------------------------------------------------------------//

inline void copyval (int* dst, int value, int n)
{
    while (n--) *dst++ = value;
}

inline void copyvec (int* dst, const int* src, int n)
{
    memcpy((void*)dst,(const void *)src, n*sizeof(int));
}

inline int cmpval (const int* dst, int value, int n)
{
    while (n--) if (*dst++ != value) return 0;
    return 1;
}

inline int cmpvec (const int* dst, const int* src, int n)
{
    while (n--) if (*dst++ != *src++) return 0;
    return 1;
}

//----------------------------------------------------------------------------//
// unsigned char - bytes
//----------------------------------------------------------------------------//

inline void copyval (unsigned char* dst, unsigned char value, int n)
{
    while (n--) *dst++ = value;
}

inline void copyvec (unsigned char* dst, const unsigned char* src, int n)
{
    memcpy((void*)dst,(const void *)src, n*sizeof(unsigned char));
}

inline int cmpval (const unsigned char* dst, unsigned char value, int n)
{
    while (n--) if (*dst++ != value) return 0;
    return 1;
}

inline int cmpvec (const unsigned char* dst, const unsigned char* src, int n)
{
    while (n--) if (*dst++ != *src++) return 0;
    return 1;
}

//----------------------------------------------------------------------------//

#endif
