/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselExpK0(x)                          dbsk0e.cc |
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
// double BesselExpK0 (double x);
//
// BesselExpK0 (x)  computes the double precision exponentially scaled 
// modified (hyperbolic) Bessel function of the third kind of 
// order zero for positive double precision argument X. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10B1, REVISION 900315, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BK0        on the interval  0.          to  4.00000E+00 
//                                        with weighted error   3.08E-33 
//                                         log weighted error  32.51 
//                               significant figures required  32.05 
//                                    decimal places required  33.11 
//
// Series for AK0        on the interval  1.25000E-01 to  5.00000E-01 
//                                        with weighted error   2.85E-32 
//                                         log weighted error  31.54 
//                               significant figures required  30.19 
//                                    decimal places required  32.33 
//
// Series for AK02       on the interval  0.          to  1.25000E-01 
//                                        with weighted error   2.30E-32 
//                                         log weighted error  31.64 
//                               significant figures required  29.68 
//                                    decimal places required  32.40 
//
//-----------------------------------------------------------------------------//

double BesselExpK0 (double x)
{
    static double bk0cs[16] = { 
	-0.0353273932339027687201140060063153,
	 0.344289899924628486886344927529213,
	 0.0359799365153615016265721303687231,
	 0.00126461541144692592338479508673447,
	 2.28621210311945178608269830297585e-5,
	 2.53479107902614945730790013428354e-7,
	 1.90451637722020885897214059381366e-9,
	 1.03496952576336245851008317853089e-11,
	 4.25981614279108257652445327170133e-14,
	 1.3744654358807508969423832544e-16,
	 3.57089652850837359099688597333333e-19,
	 7.63164366011643737667498666666666e-22,
	 1.36542498844078185908053333333333e-24,
	 2.07527526690666808319999999999999e-27,
	 2.7128142180729856e-30,
	 3.08259388791466666666666666666666e-33 
    };

    static double ak0cs[38] = { 
	-0.07643947903327941424082978270088,
	-0.02235652605699819052023095550791,
	 7.734181154693858235300618174047e-4,
	-4.281006688886099464452146435416e-5,
	 3.08170017386297474365001482666e-6,
	-2.639367222009664974067448892723e-7,
	 2.563713036403469206294088265742e-8,
	-2.742705549900201263857211915244e-9,
	 3.169429658097499592080832873403e-10,
	-3.902353286962184141601065717962e-11,
	 5.068040698188575402050092127286e-12,
	-6.889574741007870679541713557984e-13,
	 9.744978497825917691388201336831e-14,
	-1.427332841884548505389855340122e-14,
	 2.156412571021463039558062976527e-15,
	-3.34965425514956277218878205853e-16,
	 5.335260216952911692145280392601e-17,
	-8.693669980890753807639622378837e-18,
	 1.446404347862212227887763442346e-18,
	-2.452889825500129682404678751573e-19,
	 4.2337545262321715728217063424e-20,
	-7.427946526454464195695341294933e-21,
	 1.3231505293926668662779674624e-21,
	-2.390587164739649451335981465599e-22,
	 4.376827585923226140165712554666e-23,
	-8.113700607345118059339011413333e-24,
	 1.521819913832172958310378154666e-24,
	-2.886041941483397770235958613333e-25,
	 5.530620667054717979992610133333e-26,
	-1.070377329249898728591633066666e-26,
	 2.091086893142384300296328533333e-27,
	-4.121713723646203827410261333333e-28,
	 8.19348397112130764013568e-29,
	-1.642000275459297726780757333333e-29,
	 3.316143281480227195890346666666e-30,
	-6.746863644145295941085866666666e-31,
	 1.382429146318424677635413333333e-31,
	-2.851874167359832570811733333333e-32 
    };
    
    static double ak02cs[33] = { 
	-0.01201869826307592239839346212452,
	-0.009174852691025695310652561075713,
	 1.444550931775005821048843878057e-4,
	-4.013614175435709728671021077879e-6,
	 1.567831810852310672590348990333e-7,
	-7.77011043852173771031579975446e-9,
	 4.611182576179717882533130529586e-10,
	-3.158592997860565770526665803309e-11,
	 2.435018039365041127835887814329e-12,
	-2.074331387398347897709853373506e-13,
	 1.925787280589917084742736504693e-14,
	-1.927554805838956103600347182218e-15,
	 2.062198029197818278285237869644e-16,
	-2.341685117579242402603640195071e-17,
	 2.805902810643042246815178828458e-18,
	-3.530507631161807945815482463573e-19,
	 4.645295422935108267424216337066e-20,
	-6.368625941344266473922053461333e-21,
	 9.0695213109865155676223488e-22,
	-1.337974785423690739845005311999e-22,
	 2.03983602185995231552208896e-23,
	-3.207027481367840500060869973333e-24,
	 5.189744413662309963626359466666e-25,
	-8.629501497540572192964607999999e-26,
	 1.4721611831025598552080384e-26,
	-2.573069023867011283812351999999e-27,
	 4.60177408664351658737664e-28,
	-8.411555324201093737130666666666e-29,
	 1.569806306635368939301546666666e-29,
	-2.988226453005757788979199999999e-30,
	 5.796831375216836520618666666666e-31,
	-1.145035994347681332155733333333e-31,
	 2.301266594249682802005333333333e-32 
    };

    const double eta  = 0.5 * DBL_EPSILON * 0.1,
 	         xsml = sqrt(0.5 * DBL_EPSILON * 4.0);

    static int ntk0, ntak0, ntak02, first = 1;
    if (first) {	
	ntk0 = initds(bk0cs, 16, eta);
	ntak0 = initds(ak0cs, 38, eta);
	ntak02 = initds(ak02cs, 33, eta);
	first = 0;
    }

    if (x <= 0.0) {
	Matpack.Error(Mat::ArgumentDomain, "%s: %s", "BesselExpK0",
		      "x is zero or negative");
	return NAN;
    }

    double y;

    if (x > 2.0) goto L20;

    y = 0.0;
    if (x > xsml) y = x * x;
    
    return exp(x) * (-log(x * 0.5) * BesselI0(x) - 0.25 
		     + dcsevl(y * 0.5 - 1.0, bk0cs, ntk0));
    
  L20:
    if (x <= 8.0) 
	return (dcsevl((16.0 / x - 5.0) / 3.0, ak0cs, ntak0) + 1.25) / sqrt(x);
    else // if (x > 8.0) 
	return (dcsevl(16.0 / x - 1.0, ak02cs, ntak02) + 1.25) / sqrt(x);
}

//-----------------------------------------------------------------------------//
