/*-----------------------------------------------------------------------------*\
| Matpack special functions - Airy function                              dbi.cc |
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
// double AiryBi (double x);
//
// AiryBi(x) calculates the double precision Airy function of the 
// second kind for double precision argument x. 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10D, REVISION 920618, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BIF        on the interval -1.00000E+00 to  1.00000E+00 
//                                        with weighted error   1.45E-32 
//                                         log weighted error  31.84 
//                               significant figures required  30.85 
//                                    decimal places required  32.40 
//
// Series for BIG        on the interval -1.00000E+00 to  1.00000E+00 
//                                        with weighted error   1.29E-33 
//                                         log weighted error  32.89 
//                               significant figures required  31.48 
//                                    decimal places required  33.45 
//
// Series for BIF2       on the interval  1.00000E+00 to  8.00000E+00 
//                                        with weighted error   6.08E-32 
//                                         log weighted error  31.22 
//                        approx significant figures required  30.8 
//                                    decimal places required  31.80 
//
// Series for BIG2       on the interval  1.00000E+00 to  8.00000E+00 
//                                        with weighted error   4.91E-33 
//                                         log weighted error  32.31 
//                        approx significant figures required  31.6 
//                                    decimal places required  32.90 
//
//-----------------------------------------------------------------------------//

double AiryBi (double x)
{
    static double bifcs[13] = { 
       -0.016730216471986649483537423928176,
	0.10252335834249445611426362777757,
	0.0017083092507381516539429650242013,
	1.186254546774468117921645921004e-5,
	4.4932907017792133694531887927242e-8,
	1.0698207143387889067567767663628e-10,
	1.7480643399771824706010517628573e-13,
	2.0810231071761711025881891834399e-16,
	1.8849814695665416509927971733333e-19,
	1.3425779173097804625882666666666e-22,
	7.7159593429658887893333333333333e-26,
	3.6533879617478566399999999999999e-29,
	1.4497565927953066666666666666666e-32 
    };
    
    static double bigcs[13] = { 
	0.022466223248574522283468220139024,
	0.037364775453019545441727561666752,
	4.4476218957212285696215294326639e-4,
	2.4708075636329384245494591948882e-6,
	7.9191353395149635134862426285596e-9,
	1.6498079851827779880887872402706e-11,
	2.4119906664835455909247501122841e-14,
	2.6103736236091436985184781269333e-17,
	2.1753082977160323853123792e-20,
	1.4386946400390433219483733333333e-23,
	7.7349125612083468629333333333333e-27,
	3.4469292033849002666666666666666e-30,
	1.2938919273216e-33 
    };

    static double bif2cs[15] = { 
	0.0998457269381604104468284257993,
	0.47862497786300553772211467318231,
	0.025155211960433011771324415436675,
	5.8206938852326456396515697872216e-4,
	7.4997659644377865943861457378217e-6,
	6.1346028703493836681403010356474e-8,
	3.4627538851480632900434268733359e-10,
	1.4288910080270254287770846748931e-12,
	4.49627042983346418950564721792e-15,
	1.1142323065833011708428300106666e-17,
	2.2304791066175002081517866666666e-20,
	3.6815778736393142842922666666666e-23,
	5.0960868449338261333333333333333e-26,
	6.0003386926288554666666666666666e-29,
	6.0827497446570666666666666666666e-32 
    };
    
    static double big2cs[15] = {
	0.033305662145514340465176188111647,
	0.161309215123197067613287532084943,
	0.00631900730961342869121615634921173,
	1.18790456816251736389780192304567e-4,
	1.30453458862002656147116485012843e-6,
	9.37412599553521729546809615508936e-9,
	4.74580188674725153788510169834595e-11,
	1.78310726509481399800065667560946e-13,
	5.1675919278495818037427635664e-16,
	1.19004508386827125129496251733333e-18,
	2.22982880666403517277063466666666e-21,
	3.46551923027689419722666666666666e-24,
	4.53926336320504514133333333333333e-27,
	5.07884996513522346666666666666666e-30,
	4.91020674696533333333333333333333e-33 
    };
    
    const double eta   = 0.5 * DBL_EPSILON * 0.1,
	         x3sml = pow(eta, 0.3333),
	         xmax  = pow(log(DBL_MAX) * 1.5, 0.6666);

    static int nbif, nbig, nbif2, nbig2, first = 1;
    if (first) {
	nbif  = initds(bifcs,  13, eta);
	nbig  = initds(bigcs,  13, eta);
	nbif2 = initds(bif2cs, 15, eta);
	nbig2 = initds(big2cs, 15, eta);
	first = 0;
    }

    double z, theta, xm;

    if (x >= -1.0) goto L20;
    d9aimp(x, xm, theta);
    return xm * sin(theta);

  L20:
    if (x > 1.0) goto L30;
    z = 0.0;
    if (fabs(x) > x3sml) z = x * x * x;
    return dcsevl(z, bifcs, nbif) + 0.625 + x * (dcsevl(z, bigcs, nbig) + 0.4375);

  L30:
    if (x > 2.0) goto L40;
    z = (x * x * x * 2.0 - 9.0) / 7.0;
    return dcsevl(z,bif2cs,nbif2) + 1.125 + x * (dcsevl(z,big2cs,nbig2) + 0.625);

  L40:
    if (x > xmax) { 
	Matpack.Error(Mat::Overflow, "%s: %s", "AiryBi",
		      "x so big that Bi(x) overflows");
	return NAN;
    }

    return AiryExpBi(x) * exp(x * 2.0 * sqrt(x) / 3.0);
}

//-----------------------------------------------------------------------------//
