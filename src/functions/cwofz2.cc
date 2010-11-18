/*-----------------------------------------------------------------------------*\
| Matpack special functions - Faddeeva_2(z)                           cwofz2.cc |
|                                                                               |
| Matpack Library Release 1.0                                                   |
| Copyright (C) 1991,1995 by Berndt M. Gammel                                   |
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

#include "../../include/mpspecfunp.h"

//-----------------------------------------------------------------------------//
//
// complex<double> Faddeeva_2 (const complex<double>& z)
//
// Given a complex number z = (x,y), this subroutine computes 
// the value of the Faddeeva-function w(z) = exp(-z^2)*erfc(-i*z), 
// where erfc is the complex complementary error function and i 
// means sqrt(-1). 
// The relative precision is 1e-14. Only around zero, for |z| < 1e-4
// a Pade approximation is used, here the precision is
// better than 1e-11. References:
//
//     - Matta, Reichel: Math. Comp., 25 (1971), pp.339-344.
//     - Y.Luke: Mathematical Functions and their Approximations,
//               Acad. Press, N.Y.,San Francisco,London,1975, pp.119-138.
//   
// The routine is not underflow-protected but any variable can be 
// put to 0 upon underflow.
//
// The routine is not overflow-protected.
//
// This algorithm was implemented by B. M. Gammel 1990.
// It is the fastest algorithm known to me. On the average it is about two
// times as fast as algorithm wofz(z) from ACM TOMS 680.
//
//-----------------------------------------------------------------------------//

complex<double> Faddeeva_2 (const complex<double>& z)
{
    // table 1: coefficients for h = 0.5  
    static double n1[12] =
             { 0.25, 1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0,
               20.25, 25.0, 30.25, 36.0 };
    static double e1[12] =
             { 0.7788007830714049,    0.3678794411714423,
               1.053992245618643e-1,  1.831563888873418e-2,
               1.930454136227709e-3,  1.234098040866795e-4,
               4.785117392129009e-6,  1.125351747192591e-7,
               1.605228055185612e-9,  1.388794386496402e-11,
               7.287724095819692e-14, 2.319522830243569e-16 };

    // table 2: coefficients for h = 0.53 
    static double n2[12] =
             { 0.2809, 1.1236, 2.5281, 4.4944, 7.0225, 10.1124,
               13.7641, 17.9776, 22.7529, 28.09, 33.9889, 40.4496 };
    static double e2[12] =
             { 0.7551038420890235,    0.3251072991205958, 
               7.981051630007964e-2,  1.117138143353082e-2,
               0.891593719995219e-3,  4.057331392320188e-5,
               1.052755021528803e-6,  1.557498087816203e-8,
               1.313835773243312e-10, 6.319285885175346e-13,
               1.733038792213266e-15, 2.709954036083074e-18 };
    
   // tables for Pade approximation 
   static double C[7] =
            { 65536.0, -2885792.0, 69973904.0, -791494704.0,
              8962513560.0, -32794651890.0, 175685635125.0 };
   static double D[7] =
            { 192192.0, 8648640.0, 183783600.0, 2329725600.0,
              18332414100.0, 84329104860.0, 175685635125.0 };


    double *n,*e,t,u,r,s,d,f,g,h;
    complex<double> c,d2,v,w,zz;
    int i;
    
    // use Pade approximation 
    s = norm(z);
    if (s < 1e-7) {
	zz = z*z;
	v  = exp(zz);
	c  = C[0];
	d2 = D[0];
	for (i = 1; i <= 6; i++) {
	    c  = c  * zz + C[i];
	    d2 = d2 * zz + D[i];
	}
	w = 1.0 / v + complex<double>(0.0,M_2_SQRTPI) * c/d2 * z * v;
	return w;

    // use trapezoid rule 
    } else {

        // select default table 1 
        n = n1;
        e = e1;
        r = M_1_PI * 0.5;
 
        // if z is too close to a pole select table 2 
	if (fabs(imag(z)) < 0.01 && fabs(real(z)) < 6.01) {
	    h = modf(2*fabs(real(z)),&g);
	    if (h < 0.02 || h > 0.98) {
		n = n2;
		e = e2;
		r = M_1_PI * 0.53;
	    }
	   }
	
	d = (imag(z) - real(z)) * (imag(z) + real(z));
	f = 4 * real(z) * real(z) * imag(z) * imag(z);

	g = h = 0.0;
	for (i = 0; i < 12; i++) {
	    t = d + n[i];
	    u = e[i] / (t * t + f);
	    g += (s + n[i]) * u;
	    h += (s - n[i]) * u;
	}
	u = 1 / s;
	c = r * complex<double>(imag(z) * (u + 2.0 * g),
				real(z) * (u + 2.0 * h) );
	
	if (imag(z) < M_2PI) {
	    s = 2.0 / r;
	    t = s * real(z);
	    u = s * imag(z);
	    s = sin(t);
	    h = cos(t);
	    f = exp(- u) - h;
	    g = 2.0 * exp(d-u) / (s * s + f * f);
	    u = 2.0 * real(z) * imag(z);
	    h = cos(u);
	    t = sin(u);
	    c += g * complex<double>( (h * f - t * s), -(h * s + t * f));
	}
	return c;
    }
}

//-----------------------------------------------------------------------------//
