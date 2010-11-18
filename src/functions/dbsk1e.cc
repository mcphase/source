/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselExpK1(x)                          dbsk1e.cc |
|                                                                               |
| MatPack Library Release 1.0                                                   |
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
// double BesselExpK1 (double x);
//
// BesselExpK1(x) computes the double precision exponentially scaled 
// modified (hyperbolic) Bessel function of the third kind of order 
// one for positive double precision argument x. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10B1, REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BK1        on the interval  0.          to  4.00000E+00 
//                                        with weighted error   9.16E-32 
//                                         log weighted error  31.04 
//                               significant figures required  30.61 
//                                    decimal places required  31.64 
//
// Series for AK1        on the interval  1.25000E-01 to  5.00000E-01 
//                                        with weighted error   3.07E-32 
//                                         log weighted error  31.51 
//                               significant figures required  30.71 
//                                    decimal places required  32.30 
//
// Series for AK12       on the interval  0.          to  1.25000E-01 
//                                        with weighted error   2.41E-32 
//                                         log weighted error  31.62 
//                               significant figures required  30.25 
//                                    decimal places required  32.38 
//
//-----------------------------------------------------------------------------//

double BesselExpK1 (double x)
{
    static double bk1cs[16] = { 
	 0.025300227338947770532531120868533,
	-0.35315596077654487566723831691801,
	-0.12261118082265714823479067930042,
	-0.0069757238596398643501812920296083,
	-1.7302889575130520630176507368979e-4,
	-2.4334061415659682349600735030164e-6,
	-2.2133876307347258558315252545126e-8,
	-1.4114883926335277610958330212608e-10,
	-6.6669016941993290060853751264373e-13,
	-2.4274498505193659339263196864853e-15,
	-7.023863479386287597178379712e-18,
	-1.6543275155100994675491029333333e-20,
	-3.2338347459944491991893333333333e-23,
	-5.3312750529265274999466666666666e-26,
	-7.5130407162157226666666666666666e-29,
	-9.1550857176541866666666666666666e-32 
    };
    
    static double ak1cs[38] = { 
	 0.27443134069738829695257666227266,
	 0.07571989953199367817089237814929,
	-0.0014410515564754061229853116175625,
	 6.6501169551257479394251385477036e-5,
	-4.3699847095201407660580845089167e-6,
	 3.5402774997630526799417139008534e-7,
	-3.3111637792932920208982688245704e-8,
	 3.4459775819010534532311499770992e-9,
	-3.8989323474754271048981937492758e-10,
	 4.7208197504658356400947449339005e-11,
	-6.047835662875356234537359156289e-12,
	 8.1284948748658747888193837985663e-13,
	-1.1386945747147891428923915951042e-13,
	 1.654035840846228232597294820509e-14,
	-2.4809025677068848221516010440533e-15,
	 3.8292378907024096948429227299157e-16,
	-6.0647341040012418187768210377386e-17,
	 9.8324256232648616038194004650666e-18,
	-1.6284168738284380035666620115626e-18,
	 2.7501536496752623718284120337066e-19,
	-4.7289666463953250924281069568e-20,
	 8.2681500028109932722392050346666e-21,
	-1.4681405136624956337193964885333e-21,
	 2.6447639269208245978085894826666e-22,
	-4.82901575648563878979698688e-23,
	 8.9293020743610130180656332799999e-24,
	-1.6708397168972517176997751466666e-24,
	 3.1616456034040694931368618666666e-25,
	-6.0462055312274989106506410666666e-26,
	 1.1678798942042732700718421333333e-26,
	-2.277374158265399623286784e-27,
	 4.4811097300773675795305813333333e-28,
	-8.8932884769020194062336e-29,
	 1.7794680018850275131392e-29,
	-3.5884555967329095821994666666666e-30,
	 7.2906290492694257991679999999999e-31,
	-1.4918449845546227073024e-31,
	 3.0736573872934276300799999999999e-32 
    };

    static double ak12cs[33] = { 
	 0.06379308343739001036600488534102,
	 0.02832887813049720935835030284708,
	-2.475370673905250345414545566732e-4,
	 5.771972451607248820470976625763e-6,
	-2.068939219536548302745533196552e-7,
	 9.739983441381804180309213097887e-9,
	-5.585336140380624984688895511129e-10,
	 3.732996634046185240221212854731e-11,
	-2.825051961023225445135065754928e-12,
	 2.372019002484144173643496955486e-13,
	-2.176677387991753979268301667938e-14,
	 2.157914161616032453939562689706e-15,
	-2.290196930718269275991551338154e-16,
	 2.582885729823274961919939565226e-17,
	-3.07675264126846318762109817344e-18,
	 3.851487721280491597094896844799e-19,
	-5.0447948976415289771172825088e-20,
	 6.888673850418544237018292223999e-21,
	-9.77504154195011830300213248e-22,
	 1.437416218523836461001659733333e-22,
	-2.185059497344347373499733333333e-23,
	 3.4262456218092206316453888e-24,
	-5.531064394246408232501248e-25,
	 9.176601505685995403782826666666e-26,
	-1.562287203618024911448746666666e-26,
	 2.725419375484333132349439999999e-27,
	-4.865674910074827992378026666666e-28,
	 8.879388552723502587357866666666e-29,
	-1.654585918039257548936533333333e-29,
	 3.145111321357848674303999999999e-30,
	-6.092998312193127612416e-31,
	 1.202021939369815834623999999999e-31,
	-2.412930801459408841386666666666e-32 
    };

    const double eta  = 0.5 * DBL_EPSILON * 0.1,
	         xmin = exp(MpMax(log(DBL_MIN),-log(DBL_MAX)) + 0.01),
	         xsml = sqrt(0.5 * DBL_EPSILON * 4.0);

    static int ntk1, ntak1, ntak12, first = 1;
    if (first) {
	ntk1   = initds(bk1cs, 16, eta);
	ntak1  = initds(ak1cs, 38, eta);
	ntak12 = initds(ak12cs,33, eta);
	first = 0;
    }

    double y;

    if (x <= 0.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "BesselExpK1",
		      "x is zero or negative");
	return NAN;
    }

    if (x > 2.0) goto L20;

    if (x < xmin) {
	Matpack.Error(Mat::Overflow, "%s: %s", "BesselExpK1", 
		      "x so small K1(x) overflows");
    	return NAN;
    }

    y = 0.0;
    if (x > xsml) y = x * x;
    return exp(x) * (log(x * 0.5) * BesselI1(x) 
		     + (dcsevl( y * 0.5 - 1.0, bk1cs, ntk1) + 0.75) / x);

  L20:
    if (x <= 8.0) 
	return (dcsevl((16. / x - 5.0)/3.0, ak1cs, ntak1) + 1.25) / sqrt(x);
    else // if (x > 8.0) 
	return (dcsevl(16.0 / x - 1.0, ak12cs, ntak12) + 1.25) / sqrt(x);
}

//-----------------------------------------------------------------------------//
