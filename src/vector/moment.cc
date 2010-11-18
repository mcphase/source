/*------------------------------------------------------------------------*
|  Moments()  returns several statistical moments                moment.c |
|                                                                         |
|  (C) B.M.GAMMEL, 1990                                                   |
*------------------------------------------------------------------------*/

#include "vector.h"

//-----------------------------------------------------------------------------//

void Moments (Vector& data,
	      double &ave, double &adev, double &sdev,
	      double &svar, double &skew, double &kurt)
// 
//   returns the mean 'ave', absolute deviation 'adev',
//   standard deviation 'sdev', variance 'svar', skewness 'skew',
//   and kurtosis 'kurt' of a vector of doubles.
// 
{
    int j;
    double s,p;
    int lo = data.Lo(), 
	hi = data.Hi(),
	n  = hi-lo+1;

    if (hi <= lo) Matpack.Error("Moments: at least 2 elements expected");

    for (s = 0.0, j = lo; j <= hi; j++) s += data[j];
    ave = s / n;
    adev = svar = skew = kurt = 0.0;
    for (j = lo; j <= hi; j++) {
	adev += fabs( s = data[j] - ave );
	svar += (p = s*s);
	skew += (p *= s);
	kurt += (p *= s);
    }
    adev /= n;
    svar /= (n-1);
    sdev = sqrt(svar);
    if (svar) {
	skew /= (n * svar * sdev);
	kurt = kurt / (n * svar * svar) - 3.0;
    } else
	Matpack.Error("Moments: no skew/kurtosis when variance = 0");
}

//-----------------------------------------------------------------------------//
