/*-----------------------------------------------------------------------------*\
| Fast Fourier Transform, Fast Hartley Transform                          fft.c |
|                                                                               |
|                                                                               |
| MatPack Libary Release 1.0                                                    |
| Copyright (C) 1990-1995 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to use, copy, modify, and distribute this software and its         |
| documentation for any non-commercial purpose and without fee is hereby        |
| granted, provided that the above copyright notice appear in all copies        |
| and that both that copyright notice and this permission notice appear         |
| in supporting documentation. This software is provided "as is" without        |
| express or implied warranty.                                                  |
|                                                                               |
\*-----------------------------------------------------------------------------*/

#include "vector.h"

// 
// The following algorithms for the 'best and fastest implementation' of the
// FFT in C are due to a posting in the newsgroups
//
//    sci.math.num-analysis,sci.math,comp.lang.c,comp.lang.c++,comp.dsp
// by 
//    Ron Mayer, 19 Mar 1993, from  ACUSON, Mountain View, CA
//    mayer@acuson.com
//
// This seems to be a piece of really smart code for doing the fast Fourier
// transform FFT by means of a fast Hartley transform FHT. It is superior to
// the code published in the Numerical Recipies and some other published code.
// That is the reason why I decided to include it into the MatPack C++ libary.
// I left the original code unchanged exept for introducing prototypes for all 
// functions according to the ANSI standard. All functions now compile with
// ANSI C and C++. Note that the file "trigtbl.h" is included.
//
// Improvements:
// The FFT routines in their original version competely failed if the vector
// exceeded a certain size. This was due to the limited size of the trig tables.
// I extended the tables to 60 entries with 50 digits precision. Thus I
// completely rewrote them!
//
// Important Note:
// Ron Mayer notes (below) that the Fast Hartley Trasform code below is 
// restriceted by patents (U.S. Patent No. 4,646,256 (1987)). 
// As noted in Computer in Pysics, Vol. 9, No. 4,
// Jul/Aug 1995 pp 373-379 it was placed in the public domain by the
// Board of Trustees of Stanford University in 1994!
//
// Berndt M. Gammel, 1994.
//
//
// FFT and FHT routines  Copyright 1988, 1993 by Ron Mayer
// -------------------------------------------------------
//
//  void fht (double* fz, int n);
//      Does a hartley transform of 'n' points in the array 'fz'.
//
//  void fft (int n, double* real, double* imag);
//      Does a fourier transform of 'n' points of the 'real' and
//      'imag' arrays.
//
//  void ifft (int n, double* real, double* imag);
//      Does an inverse fourier transform of 'n' points of the 'real'
//      and 'imag' arrays.
//
//  void realfft (int n, double* real);
//      Does a real-valued fourier transform of 'n' points of the
//      'real' and 'imag' arrays.  The real part of the transform ends
//      up in the first half of the array and the imaginary part of the
//      transform ends up in the second half of the array.
//
//  void realifft (int n, double* real);
//      The inverse of the realfft() routine above.
// 
//      
// NOTE: This routine uses at least 2 patented algorithms, and may be
//       under the restrictions of a bunch of different organizations.
//       Although I wrote it completely myself; it is kind of a derivative
//       of a routine I once authored and released under the GPL, so it
//       may fall under the free software foundation's restrictions;
//       it was worked on as a Stanford Univ project, so they claim
//       some rights to it; it was further optimized at work here, so
//       I think this company claims parts of it.  The patents are
//       held by R. Bracewell (the FHT algorithm) and O. Buneman (the
//       trig generator), both at Stanford Univ.
//       If it were up to me, I'd say go do whatever you want with it;
//       but it would be polite to give credit to the following people
//       if you use this anywhere:
//           Euler     - probable inventor of the fourier transform.
//           Gauss     - probable inventor of the FFT.
//           Hartley   - probable inventor of the hartley transform.
//           Buneman   - for a really cool trig generator
//           Mayer(me) - for authoring this particular version and
//                       including all the optimizations in one package.
//       Thanks,
//       Ron Mayer; mayer@acuson.com
//
//
// The follwing comments from Ron Mayer came along with the code:
// --------------------------------------------------------------
// 
// As I'm sure you realize, 'best' depends quite a bit on what you
// consider best, and 'fastest' depends largely on the number of points,
// what compiler, and what architecture you're using.  If it's amazingly
// critical to have the best and fastest, you're probably best off
// writing your own to meet your specific application.  Short of
// specialized hardware, if you can use a real valued FFT this would be
// your best optimization (should be faster by a factor of 2).
// 
// It seems that two of the source code sections in your summary used
// 'a duhamel-holman split radix fft'.  Perhaps this is considered a
// 'standard piece of code', but I do know of two different fft's which
// in my opinion qualify as both 'faster' and 'better'. 
// 
// I tried comparing my code to the code by Dave Edelblute that you
// included in your previous posting; but the times for that code was so
// much worse (twice as slow) that it's unlikely that such code would
// have been in a posting labeled best and fastest; so I probably made an
// error in compiling it.  In the nearly identical test programs I
// include in my other posting, they both seem to produce the same
// results; but who knows.  I'll include the code duhamel-holman code I
// tested in the second part of this posting so someone can point out
// where I screwed up.  (If anyone tries it, please let me know if you
// get the same results...)
// 
// The best' published FFT routine I have studied recently is based on
// Singleton's mixed-radix FFT algorithm.  His routine has a number of
// optimizations which include a decent trig function generator, and
// doing a quite optimized radix-4 transforms and only 0 or 1 radix-2
// transforms for any power of 2.  It is also neat in that it works with
// any sized array, not just powers of 2!  I believe the published
// version is in fortran; but someone at work translated it to C, and 
// it seems to be ~25% faster than the duhamel-holman routine you posted
// in your summary (if I compiled it correctly).  I can probably dig up a
// reference to this code next time I dig through all my school stuff if
// anyone really needs it.
// 
// The 'fastest' (for typical computers: single processor, non-vector,
// fast-integer-math, slower-floating-point-math, slow-trig-function) FFT
// routine I have ever seen is one I wrote myself; trying to incorporate
// as many optimizations I could find in various IEEE publications while
// I was at college.  As you can see in the file 'summary.sparc' included
// below, it is nearly twice as fast as the duhamel-holman routine you
// posted in your posting.
// 
// The routine I came up with includes the following optimizations:
// 
//   1) It is a wrapper around a highly optimized FHT (hartley transform)
//      A FHT better localizes memory accesses for large transforms,
//      thereby avoiding paging.  Hartley transforms are also much easier
//      to add optimization tricks too; more than making up for the
//      overhead of converting the hartley transfrom to a fourier
//      transform.  Another advantage is that the transformation from
//          FHT -> real_valued_FFT
//      is faster than the transformation from
//          1/2pointFFT-> real_valued_FFT
//      so my real-valued fft is even better when compared to most
//      published real valued ffts.
// 
//   2) Avoid multiplications by 1 and zero (and sometimes sqrt2).
//      Many published routines such as Numerical Recipes seem to spend
//      a lot of time multiplying by cos(0), cos(pi), etc. and almost all
//      seem to do 1/sqrt_2*x+1/sqrt_2*y instead of 1/sqrt_2*(x+y).
// 
//   3) Faster trig generation.
//      Most algorithms use 1 'sin' library call for each level of the
//      transform; and 2 real multiplications and 2 real additions for
//      each needed trig value within it's loop.
// 
//      I use a stable algorithm to generate each trig value using 1 real
//      multiplication and 1 real addition for each value using a small
//      (log(n)) table of trig values.  The tradeoff is that I require
//      much more integer arithmetic for this calculation, including a
//      (n*log(n)) loop; but for multiples of pi/16 or so, my routine
//      still seems faster.  By taking advantage of the fact that
//      values required for FFTs are for evenly spaced angles, I avoid
//      all calls to slow trig library functions which are unnecessarily
//      complex because they need to work for arbitrary values.
// 
//   4) Generate less trig values
//      I use the identities sin(x)=sin(pi-x)=-sin(pi+x)=-sin(-x),etc. to
//      avoid excessive trig calculations, and sin(2x) = 2*cos(x)*sin(x)
//      to allow simpler trig calculations when accuracy permits.  A more
//      stable than average trig generator mentioned in (3) above allows
//      me to use the unstable sin(2x) = 2*cos(x)*sin(x) hack for every
//      other 'level' in the FFT without the usual loss of accuracy.
// 
//   5) Mixed 2,4-radix inner loop.
//      By doing two levels in the inner loop, I gain all the advantages
//      of a radix-4 algorithm; primarily reducing integer arithmetic and
//      memory access.  This has a great affect on large arrays when
//      paging occurs.
// 
//   6) Unrolling loops and variable storage to optimize register
//      allocation.  I try not to require storing too many values in
//      variables at any one time to ease a compilers register
//      allocation. 
// 
//   7) Remove low levels of the transform out of the loop.  It's
//      significantly faster to do 8 point transforms explicitly; rather
//      than using the general loop.
// 
// One catch to this routine is that at least two of the algorithms used
// by it are patented(!) (just about any FHT is patented by R. Bracewell;
// and the stable trig generator is patented by O. Buneman; both at
// Stanford Univ.)  Who owns the copyright rights to it is also probably
// being debated; since it is a derivative work of a GNU-licensed
// routine, so subject to their restrictions; it was worked on for a
// Stanford project, so they have a claim on it;  and I optimized it
// further working for this company, so they probably claim parts of it.
// 
// Considering Gauss apparently used the equivalent of real valued FFTs
// in 1805; and Euler did fourier transforms it in the mid 1700s; I'm
// amazed that people still want to claim this math.
// 
//
// Here are the test results posted by Ron Mayer
// ---------------------------------------------
//
// This file contains a benchmark results of a number of popular FFT
// algorithms.  The algorithms compared are:
// 
//     FFT-numrec
//         The FFT from numerical recipies, converted to double precision
//     FFT-duhamel
//         A 'duhamel-holman split radix fft' from "electronics letters,
//         jan. 5, 1994", coded by Dave Edelblute, edelblut@cod.nosc.mil
//     FFT-wang
//         Singleton's arbitrary-radix FFT translated to C and coded by
//         John Wang, wang@acuson.com
//     FFT-mayer
//         An original FFT by Ron Mayer (mayer@acuson.com)
//     real-FFT-numrec
//         The real valued FFT from numerical recipies, converted to
//         double precision.
//     real-FFT-mayer
//         An original real valued FFT by Ron Mayer (mayer@acuson.com)
// 
// I compiled each of the programs using gcc 2.0 with the -O4 flag on a
// Sun Sparc 1; and timed (using the "clock()" function in SunOS) a
// number of iterations of forward and reverse transforms of a known data
// set.  At the end of the iterations of forward and reverse transforms I
// compared the data with the original to check for accumulated errors.
// 
// algorithm                  # of       # of     time           errors
//   used                   iterations   points
//
// n=4
// FFT-numrec                (16386       4):   4466488 CPU us ;ssq errors 0.0
// FFT-duhamel               (16386       4):   2016586 CPU us ;ssq errors 0.0
// FFT-wang                  (16386       4):   3299868 CPU us ;ssq errors 0.0
// FFT-mayer                 (16386       4):   1333280 CPU us ;ssq errors 0.0
// real-FFT-numrec           (16386       4):   3133208 CPU us ;ssq errors 0.0
// real-FFT-mayer            (16386       4):    666640 CPU us ;ssq errors 0.0
//
// n=128
// FFT-numrec                (514       128):   3883178 CPU us ;ssq errors 4.1e-21
// FFT-duhamel               (514       128):   6349746 CPU us ;ssq errors 8.6e-22
// FFT-wang                  (514       128):   3866512 CPU us ;ssq errors 1.5e-09
// FFT-mayer                 (514       128):   2999880 CPU us ;ssq errors 6.9e-22
// real-FFT-numrec           (514       128):   2333240 CPU us ;ssq errors 4.1e-21
// real-FFT-mayer            (514       128):   1433276 CPU us ;ssq errors 6.9e-22
//
// n=2048
// FFT-numrec                (34       2048):   5733104 CPU us ;ssq errors 8.6e-19
// FFT-duhamel               (34       2048):   8849646 CPU us ;ssq errors 3.2e-20
// FFT-wang                  (34       2048):   5783102 CPU us ;ssq errors 2.2e-08
// FFT-mayer                 (34       2048):   4649814 CPU us ;ssq errors 9.4e-20
// real-FFT-numrec           (34       2048):   3116542 CPU us ;ssq errors 1.6e-18
// real-FFT-mayer            (34       2048):   2183246 CPU us ;ssq errors 9.4e-20
//
// n=32768
// FFT-numrec                (4       32768):  18732584 CPU us ;ssq errors 1.5e-16
// FFT-duhamel               (4       32768):  22632428 CPU us ;ssq errors 3.7e-18
// FFT-wang                  (4       32768):  16299348 CPU us ;ssq errors 1.1e-06
// FFT-mayer                 (4       32768):  13849446 CPU us ;ssq errors 1.2e-17
// real-FFT-numrec           (4       32768):   9999600 CPU us ;ssq errors 1.9e-16
// real-FFT-mayer            (4       32768):   6716398 CPU us ;ssq errors 1.2e-17
//

// prototypes
void fht (double* fz, int n);
void fft (int n, double* real, double* imag);
void ifft (int n, double* real, double* imag);
void realfft (int n, double* real);
void realifft (int n, double* real);



void FFT (Vector& real, Vector& imag)
//
// Computes the Fast Fourier Transform of N points given by the vectors
// of real and imaginary parts. N must be a power of two. 
//
{
    // verify index range
    if ( ! MatchingIndexRange(real,imag) )
	 Matpack.Error(Mat::UnspecifiedError,"FFT: index range mismatch (%d,%d) and (%d,%d)",
	       real.Lo(),real.Hi(),imag.Lo(),imag.Hi());

    int m, lo = real.Lo(), hi = real.Hi(), n = hi-lo+1;

    // make shure we have a power of two   
    for (m = n; even(m) && m > 0; m /= 2);
    if (m != 1) Matpack.Error(Mat::UnspecifiedError,"FFT: number of elements must be a power of 2 (but is %d)",n);

    // adjust pointers to vectors so that indexing starts from 0
    double *r = &real[lo],
	   *i = &imag[lo];

    // do it
    fft(n,r,i);
}


void InverseFFT (Vector& real, Vector& imag)
//
// Computes the inverse of the Fourier Transform
// given in the vectors of real and imaginary parts. 
// N must be a power of two. 
//
{
    // verify index range
    if ( ! MatchingIndexRange(real,imag) )
	 Matpack.Error(Mat::UnspecifiedError,"InverseFFT: index range mismatch (%d,%d) and (%d,%d)",
	       real.Lo(),real.Hi(),imag.Lo(),imag.Hi());

    int m, lo = real.Lo(), hi = real.Hi(), n = hi-lo+1;

    // make shure we have a power of two   
    for (m = n; even(m) && m > 0; m /= 2);
    if (m != 1) 
	Matpack.Error(Mat::UnspecifiedError,"InverseFFT: number of elements must be a power of 2 (but is %d)",n);

    // adjust pointers to vectors so that indexing starts from 0
    double *r = &real[lo],
	   *i = &imag[lo];

    // do it
    ifft(n,r,i);
}


void FFT (Vector& real)
//
// Computes the Fast Fourier Transform of N points given in the vector
// of real parts. N must be a power of two. 
// The real part of the transform ends up in the first half of the array 
// and the imaginary part of the transform ends up in the second half 
// of the array.
//
{
    int m, lo = real.Lo(), hi = real.Hi(), n = hi-lo+1;

    // make shure we have a power of two   
    for (m = n; even(m) && m > 0; m /= 2);
    if (m != 1) 
	Matpack.Error(Mat::UnspecifiedError,"FFT: number of elements must be a power of 2 (but is %d)", n);

    // adjust pointers to vectors so that indexing starts from 0
    double *r = &real[lo];

    // do it
    realfft(n,r);
}


void InverseFFT (Vector& real)
//
// Counterpart to the real valued FFT above. Computes the inverse of the 
// Fourier Transform in the vector which contains the real parts in the 
// first half of the vector and the imaginary parts in the second half 
// of the vector.  N must be a power of two. 
//
{
    int m, lo = real.Lo(), hi = real.Hi(), n = hi-lo+1;

    // make shure we have a power of two   
    for (m = n; even(m) && m > 0; m /= 2);
    if (m != 1) 
	Matpack.Error(Mat::UnspecifiedError,"InverseFFT: number of elements must be a power of 2 (but is %d)",n);

    // adjust pointers to vectors so that indexing starts from 0
    double *r = &real[lo];

    // do it
    realifft(n,r);
}


void FHT (Vector& f)
//
// Computes the Fast Hartley Transform of N points given in the vector f.
// N must be a power of two. 
//
{
    int m, lo = f.Lo(), hi = f.Hi(), n = hi-lo+1;

    // make shure we have a power of two   
    for (m = n; even(m) && m > 0; m /= 2);
    if (m != 1) 
	Matpack.Error(Mat::UnspecifiedError,"FHT: number of elements must be a power of 2 (but is %d)",n);

    // adjust pointers to vectors so that indexing starts from 0
    double *r = &f[lo];

    // do it
    fht(r,n);
}



//  REAL is usually defined to be double, could also be float,
//  but all routines execpt fht use double ! So don't change it !

#define REAL double
#define GOOD_TRIG       //could also use #define FAST_TRIG
#include "trigtbl.h"

#define SQRT2_2   0.70710678118654752440084436210484
#define SQRT2   2*SQRT2_2


//----------------------------------------------------------------------------//


void fht (REAL* fz, int n)
{
    REAL a,b;
    REAL c1,s1,s2,c2,s3,c3,s4,c4;
    REAL f0,g0,f1,g1,f2,g2,f3,g3;
    int i,k,k1,k2,k3,k4,kx;
    REAL *fi,*fn,*gi;

    TRIG_VARS;

    for (k1=1,k2=0;k1<n;k1++) {
	REAL a;
	for (k=n >> 1; (!((k2^=k)&k)); k >>= 1);
	if (k1>k2) {
	    a=fz[k1];fz[k1]=fz[k2];fz[k2]=a;
	}
    }

    for ( k=0 ; (1 << k) < n ; k++ );

    k &= 1;

    if (k == 0) {

	for (fi=fz,fn=fz+n;fi<fn;fi += 4) {
	    REAL f0,f1,f2,f3;
	    f1     = fi[0]-fi[1];
	    f0     = fi[0]+fi[1];
	    f3     = fi[2]-fi[3];
	    f2     = fi[2]+fi[3];
	    fi[2 ] = (f0-f2);	
	    fi[0 ] = (f0+f2);
	    fi[3 ] = (f1-f3);	
	    fi[1 ] = (f1+f3);
	}

    } else {

	for (fi=fz,fn=fz+n,gi=fi+1;fi<fn;fi += 8,gi += 8) {
	    REAL s1,c1,s2,c2,s3,c3,s4,c4,g0,f0,f1,g1,f2,g2,f3,g3;
	    c1     = fi[0] - gi[0];
	    s1     = fi[0] + gi[0];
	    c2     = fi[2] - gi[2];
	    s2     = fi[2] + gi[2];
	    c3     = fi[4] - gi[4];
	    s3     = fi[4] + gi[4];
	    c4     = fi[6] - gi[6];
	    s4     = fi[6] + gi[6];
	    f1     = (s1 - s2);	
	    f0     = (s1 + s2);
	    g1     = (c1 - c2);	
	    g0     = (c1 + c2);
	    f3     = (s3 - s4);	
	    f2     = (s3 + s4);
	    g3     = SQRT2*c4;		
	    g2     = SQRT2*c3;
	    fi[4 ] = f0 - f2;
	    fi[0 ] = f0 + f2;
	    fi[6 ] = f1 - f3;
	    fi[2 ] = f1 + f3;
	    gi[4 ] = g0 - g2;
	    gi[0 ] = g0 + g2;
	    gi[6 ] = g1 - g3;
	    gi[2 ] = g1 + g3;
	}

    }

    if (n<16) return;

    do {
	REAL s1,c1;
	k  += 2;
	k1  = 1  << k;
	k2  = k1 << 1;
	k4  = k2 << 1;
	k3  = k2 + k1;
	kx  = k1 >> 1;
	fi  = fz;
	gi  = fi + kx;
	fn  = fz + n;

	do {
	    REAL g0,f0,f1,g1,f2,g2,f3,g3;
	    f1      = fi[0 ] - fi[k1];
	    f0      = fi[0 ] + fi[k1];
	    f3      = fi[k2] - fi[k3];
	    f2      = fi[k2] + fi[k3];
	    fi[k2]  = f0	  - f2;
	    fi[0 ]  = f0	  + f2;
	    fi[k3]  = f1	  - f3;
	    fi[k1]  = f1	  + f3;
	    g1      = gi[0 ] - gi[k1];
	    g0      = gi[0 ] + gi[k1];
	    g3      = SQRT2  * gi[k3];
	    g2      = SQRT2  * gi[k2];
	    gi[k2]  = g0	  - g2;
	    gi[0 ]  = g0	  + g2;
	    gi[k3]  = g1	  - g3;
	    gi[k1]  = g1	  + g3;
	    gi     += k4;
	    fi     += k4;
	} while (fi<fn);

	TRIG_INIT(k,c1,s1);

	for (i=1;i<kx;i++) {
	    REAL c2,s2;
	    TRIG_NEXT(k,c1,s1);
	    c2 = c1*c1 - s1*s1;
	    s2 = 2*(c1*s1);
	    fn = fz + n;
	    fi = fz +i;
	    gi = fz +k1-i;

	    do {
		REAL a,b,g0,f0,f1,g1,f2,g2,f3,g3;
		b       = s2*fi[k1] - c2*gi[k1];
		a       = c2*fi[k1] + s2*gi[k1];
		f1      = fi[0 ]    - a;
		f0      = fi[0 ]    + a;
		g1      = gi[0 ]    - b;
		g0      = gi[0 ]    + b;
		b       = s2*fi[k3] - c2*gi[k3];
		a       = c2*fi[k3] + s2*gi[k3];
		f3      = fi[k2]    - a;
		f2      = fi[k2]    + a;
		g3      = gi[k2]    - b;
		g2      = gi[k2]    + b;
		b       = s1*f2     - c1*g3;
		a       = c1*f2     + s1*g3;
		fi[k2]  = f0        - a;
		fi[0 ]  = f0        + a;
		gi[k3]  = g1        - b;
		gi[k1]  = g1        + b;
		b       = c1*g2     - s1*f3;
		a       = s1*g2     + c1*f3;
		gi[k2]  = g0        - a;
		gi[0 ]  = g0        + a;
		fi[k3]  = f1        - b;
		fi[k1]  = f1        + b;
		gi     += k4;
		fi     += k4;
	    } while (fi<fn);
	}

	TRIG_RESET(k,c1,s1);

    } while (k4<n);
}


//----------------------------------------------------------------------------//


void ifft (int n, double* real, double* imag)
{
    fht(real,n);
    fht(imag,n);

    for (int i=1,j=n-1,k=n/2;i<k;i++,j--) {
	double a,b,c,d,q,r,s,t;
	a = real[i]; b = real[j];  q=a+b; r=a-b;
	c = imag[i]; d = imag[j];  s=c+d; t=c-d;
	imag[i] = (s+r)*0.5;  imag[j] = (s-r)*0.5;
	real[i] = (q-t)*0.5;  real[j] = (q+t)*0.5;
    }
}


//----------------------------------------------------------------------------//


void realfft (int n, double* real)
{
    fht(real,n);

    for (int i=1,j=n-1,k=n/2;i<k;i++,j--) {
	double a,b;
	a = real[i];
	b = real[j];
	real[j] = (a-b)*0.5;
	real[i] = (a+b)*0.5;
    }
}


//----------------------------------------------------------------------------//


void fft (int n, double* real, double* imag)
{
    for (int i=1,j=n-1,k=n/2;i<k;i++,j--) {
	double a,b,c,d, q,r,s,t;
	a = real[i]; b = real[j];  q=a+b; r=a-b;
	c = imag[i]; d = imag[j];  s=c+d; t=c-d;
	real[i] = (q+t)*.5; real[j] = (q-t)*.5;
	imag[i] = (s-r)*.5; imag[j] = (s+r)*.5;
    }

    fht(real,n);
    fht(imag,n);
}


//----------------------------------------------------------------------------//


void realifft (int n, double* real)
{
    for (int i=1,j=n-1,k=n/2;i<k;i++,j--) {
	double a,b;
	a = real[i];
	b = real[j];
	real[j] = (a-b);
	real[i] = (a+b);
    }

    fht(real,n);
}


//----------------------------------------------------------------------------//
