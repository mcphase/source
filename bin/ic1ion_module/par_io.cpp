/* pario.cpp
 *
 * Parameter input/output routines
 *
 * Functions:
 *   int  getdim(int n, orbital l);                                            // Number of states = ^{4l+2}C_{n}
 *   void getfromionname(std::string &ion, icpars &flags);                     // Gets free ion parameters from tables
 *   void ic_parsecfpars(std::string &n, std::string &v, icpars &p, int l=1);  // Parses CF parameter for k and q
 *   void ic_parseinput(const char *file, icpars &flags);                      // Parses file for 1-ion pars & phys prop.
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2010 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"
#include <fstream>
#include <iomanip>
#include <ctime>

#define Fk2F(F) F[1] *= 225.; F[2] *= 1089.; F[3] *= (184041./25); 

// --------------------------------------------------------------------------------------------------------------- //
// Looks up values of Slater and spin orbit radial integrals from published spectroscopic data
// --------------------------------------------------------------------------------------------------------------- //
void getfromionname(std::string &ionname, icpars &pars)
{
   int n,i; 
   orbital l = (orbital)3;    // Defaults for f-electrons
   std::vector<double> F,a;
   double xi = 0.;
   double B=0.,C=0.; bool flg3d = false, flgBC = false; int S2;
   F.assign(4,0.); a.assign(3,0.);
   pars.ionname.assign(ionname);
   strtolower(ionname);
   ionname.erase(ionname.find("+")+1);
#define IONCMP ionname.compare
   // Trivalent Lanthanides, from Carnall et al. J. Chem. Phys. v90, pp3443, 1989. All values in cm^{-1}, obtained from RE3+:LaCl_3
   //                               F^2            F^4            F^6        spin-orbit              alpha      beta        gamma
        if(IONCMP("ce3+")==0)                                              { xi = 647.3; n = 1; }
   else if(IONCMP("pr3+")==0) { F[1] = 68878.; F[2] = 50347.; F[3] = 32901.; xi = 751.7; n = 2;  a[0] = 16.23; a[1] = -566.6; a[2] = 1371.; }
   else if(IONCMP("nd3+")==0) { F[1] = 73018.; F[2] = 52789.; F[3] = 35757.; xi = 885.3; n = 3;  a[0] = 21.34; a[1] = -593.0; a[2] = 1445.; }
   else if(IONCMP("pm3+")==0) { F[1] = 76400.; F[2] = 54900.; F[3] = 37700.; xi = 1025.; n = 4;  a[0] = 20.50; a[1] = -560.;  a[2] = 1475.; }
   else if(IONCMP("sm3+")==0) { F[1] = 79805.; F[2] = 57175.; F[3] = 40250.; xi = 1176.; n = 5;  a[0] = 20.16; a[1] = -566.9; a[2] = 1500.; }
   else if(IONCMP("eu3+")==0) { F[1] = 83125.; F[2] = 59268.; F[3] = 42560.; xi = 1338.; n = 6;  a[0] = 20.16; a[1] = -566.9; a[2] = 1500.; }
   else if(IONCMP("gd3+")==0) { F[1] = 85669.; F[2] = 60825.; F[3] = 44776.; xi = 1508.; n = 7;  a[0] = 18.92; a[1] = -600.;  a[2] = 1575.; }
   else if(IONCMP("tb3+")==0) { F[1] = 88995.; F[2] = 62919.; F[3] = 47252.; xi = 1707.; n = 8;  a[0] = 18.40; a[1] = -590.9; a[2] = 1650.; }
   else if(IONCMP("dy3+")==0) { F[1] = 91903.; F[2] = 64372.; F[3] = 49386.; xi = 1913.; n = 9;  a[0] = 18.02; a[1] = -633.4; a[2] = 1790.; }
   else if(IONCMP("ho3+")==0) { F[1] = 94564.; F[2] = 66397.; F[3] = 52022.; xi = 2145.; n = 10; a[0] = 17.15; a[1] = -607.9; a[2] = 1800.; }
   else if(IONCMP("er3+")==0) { F[1] = 97483.; F[2] = 67904.; F[3] = 54010.; xi = 2376.; n = 11; a[0] = 17.79; a[1] = -582.1; a[2] = 1800.; }
   else if(IONCMP("tm3+")==0) { F[1] = 100134.;F[2] = 69613.; F[3] = 55975.; xi = 2636.; n = 12; a[0] = 17.26; a[1] = -624.5; a[2] = 1820.; }
   else if(IONCMP("yb3+")==0)                                              { xi = 2928.; n = 13;}
   // Trivalent Actinides, from Carnall et al. J. Chem. Phys. v90, pp3443, 1989. All values in cm^{-1}, from An3+:LaCl_3 and An3+:LaF_3
   else if(IONCMP("u3+")==0)  { F[1] = 39611.; F[2] = 32960.; F[3] = 23084.; xi = 1626.; n = 3;  a[0] = 29.26; a[1] = -824.6; a[2] = 1820.; }
   else if(IONCMP("np3+")==0) { F[1] = 45382.; F[2] = 37242.; F[3] = 25644.; xi = 1937.; n = 4;  a[0] = 31.78; a[1] = -728.0; a[2] = 1820.; }
   else if(IONCMP("pu3+")==0) { F[1] = 48679.; F[2] = 39333.; F[3] = 27647.; xi = 2242.; n = 5;  a[0] = 30.00; a[1] = -678.3; a[2] = 1820.; }
   else if(IONCMP("am3+")==0) { F[1] = 51900.; F[2] = 41600.; F[3] = 29400.; xi = 2564.; n = 6;  a[0] = 26.71; a[1] = -426.6; a[2] = 1820.; }
   else if(IONCMP("cm3+")==0) { F[1] = 55055.; F[2] = 43938.; F[3] = 32876.; xi = 2889.; n = 7;  a[0] = 29.42; a[1] = -362.9; a[2] = 1820.; }
   else if(IONCMP("bk3+")==0) { F[1] = 57697.; F[2] = 45969.; F[3] = 32876.; xi = 3210.; n = 8;  a[0] = 29.56; a[1] = -564.9; a[2] = 1820.; }
   else if(IONCMP("cf3+")==0) { F[1] = 60464.; F[2] = 48026.; F[3] = 34592.; xi = 3572.; n = 9;  a[0] = 27.36; a[1] = -587.5; a[2] = 1820.; }
   else if(IONCMP("es3+")==0) { F[1] = 63174.; F[2] = 50034.; F[3] = 36199.; xi = 3944.; n = 10; a[0] = 30.21; a[1] = -761.0; a[2] = 1820.; }
   else if(IONCMP("fm3+")==0) { F[1] = 65850.; F[2] = 52044.; F[3] = 37756.; xi = 4326.; n = 11; a[0] = 30.;   a[1] = -600.;  a[2] = 1820.; }
   else if(IONCMP("md3+")==0) { F[1] = 68454.; F[2] = 54048.; F[3] = 39283.; xi = 4715.; n = 12; a[0] = 30.;   a[1] = -600.;  a[2] = 1820.; }
   else if(IONCMP("no3+")==0)                                              { xi = 5144.; n = 13;}
   // Trivalent parameters from Sytsma et al., Phys. Rev. B, v52, pp12668, 1995, in cm^{-1}, obtained from An4+:LuPO_4
 //else if(IONCMP("cm3+")==0) { F[1] = 54669.1;F[2] = 44759.8; F[3]= 33021.4; xi= 2867.7;n = 7;  a[0] = 30.27; a[1] = -981.6; a[2] = 749.3; }
 //else if(IONCMP("gd3+")==0) { F[1] = 84075.; F[2] = 61410.8; F[3]= 44425.9; xi= 1494.; n = 7;  a[0] = 18.92; a[1] = -600.;  a[2] = 1575.; }
   // Tetravalent Actinides, from Conway, J. Chem. Phys. v41, pp904, 1964. All values in cm^{-1}, using hydrogenic wavefunctions.
   //                               F_2            F_4            F_6        spin-orbit
   else if(IONCMP("pa4+")==0)                                              {  xi= 1490.; n = 1; }
 //else if(IONCMP("u4+")==0)  { F[1] = 206.;   F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 1870.; n = 2; Fk2F(F); }
 //else if(IONCMP("np4+")==0) { F[1] = 223.8;  F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 2193.; n = 3; Fk2F(F); }
 //else if(IONCMP("pu4+")==0) { F[1] = 242.9;  F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 2429.; n = 3; Fk2F(F); }
   else if(IONCMP("am4+")==0) { F[1] = 282.1;  F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 2821.; n = 3; Fk2F(F); }
   else if(IONCMP("cm4+")==0) { F[1] = 307.0;  F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 3042.; n = 3; Fk2F(F); }
   // Tetravalent Actinides parameters from Poirot et al., Phys. Rev. B, v39, pp6388, 1989, in cm^{-1}, obtained from An4+:ZrSiO_4
   //                               F^2            F^4            F^6        spin-orbit              alpha      beta        gamma
   else if(IONCMP("u4+")==0)  { F[1] = 44258.; F[2] = 40293.; F[3] = 31287.; xi = 1740.; n = 2;  a[0] = 23.; }
   else if(IONCMP("np4+")==0) { F[1] = 47949.; F[2] = 41455.; F[3] = 26528.; xi = 2088.; n = 3;  a[0] = 39.2;  a[1] = -610.;  a[2] = 1200.; }
   else if(IONCMP("pu4+")==0) { F[1] = 49394.; F[2] = 39495.; F[3] = 30684.; xi = 2366.; n = 3;  a[0] = 32.3;  a[1] = -783.;  a[2] = 1200.; }
   // U4+:CsCdBr_3 and other parameters, from Karbowiak et al., Chemical Physics, v308, p135, 2005 (Elsevier)
 //else if(IONCMP("u4+")==0)  { F[1] = 45601.; F[2] = 38622.; F[3] = 28423.; xi = 1718.; n = 2;  a[0] = 29.;   a[1] = -1440.; a[2] = 1532.; }
 //else if(IONCMP("u3+")==0)  { F[1] = 36918.; F[2] = 32942.; F[3] = 19906.; xi = 1601.; n = 3;  a[0] = 27.;   a[1] = -830.;  a[2] = 1093.; }
 //else if(IONCMP("cm3+")==0) { F[1] = 53309.; F[2] = 45993.; F[3] = 31047.; xi = 2800.; n = 7;  a[0] = 30.2;  a[1] = -947.;  a[2] = 910.;  }
 //else if(IONCMP("pr3+")==0) { F[1] = 67459.; F[2] = 49029.; F[3] = 32366.; xi = 741.07;n = 2;  a[0] = 22.98; a[1] = -682.98;a[2] = 1422.; }
  
   // All d-electron parameters from Appendix of AS Chakravarty, Introduction to Magnetic Properties of Solid, Wiley, 1980. Original work also cited
   // 3d ions parameters from JS Grifiths, The Theory of Transition Metal Ions, CUP, 1961 (A and B); lambda from TM Dunn, Trans. Faraday Soc. v57, 1441 (1961)
   //   where 0 is shown, parameter not known... (where 0. is shown parameter is zero!) where theory_lamdba shows "-" parameter is actually theoretical...
   //                                      lambda_experimental lambda_theory-[from M Blume and RE Watson, Proc. R. Soc. Lon. A v271, 565 (1963)]
   else if(IONCMP("sc2+")==0) { B = 0.;    C = 0.;     xi = 79.;  /*xi=86.;*/ n = 1; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("ti")==0)   { B = 560.;  C = 1840.;  xi = 0;    /*xi=   ;*/ n = 4; l=D; flg3d=1; flgBC=1; } //
   else if(IONCMP("ti+")==0)  { B = 682.;  C = 2481.;  xi = 0;    /*xi=   ;*/ n = 3; l=D; flg3d=1; flgBC=1; } //
   else if(IONCMP("ti2+")==0) { B = 718.;  C = 2629.;  xi = 60.;  /*xi=61.;*/ n = 2; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("ti3+")==0) { B = 0.;    C = 0.;     xi = 154.; /*xi=159;*/ n = 1; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("v")==0)    { B = 578.;  C = 2273.;  xi = 0;    /*xi=   ;*/ n = 5; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("v+")==0)   { B = 659.;  C = 2417.;  xi = 0;    /*xi=   ;*/ n = 4; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("v2+")==0)  { B = 766.;  C = 2855.;  xi = 55.;  /*xi=57.;*/ n = 3; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("v3+")==0)  { B = 861.;  C = 4165.;  xi = 106.; /*xi=104;*/ n = 2; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("v4+")==0)  { B = 0.;    C = 0.;     xi = 248.; /*xi=255;*/ n = 1; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("cr")==0)   { B = 790.;  C = 2520.;  xi = 0;    /*xi=   ;*/ n = 6; l=D; flg3d=1; flgBC=1; } //
   else if(IONCMP("cr+")==0)  { B = 710.;  C = 2790.;  xi = 0;    /*xi=   ;*/ n = 5; l=D; flg3d=1; flgBC=1; } //
   else if(IONCMP("cr2+")==0) { B = 830.;  C = 3430.;  xi = 58.;  /*xi=59.;*/ n = 4; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("cr3+")==0) { B = 1030.; C = 3850.;  xi = 91.;  /*xi=91.;*/ n = 3; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("cr4+")==0) { B = 1039.; C = 4238.;  xi = 164.; /*xi=163;*/ n = 2; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("mn")==0)   { B = 720.;  C = 3087.;  xi = 0;    /*xi=   ;*/ n = 7; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("mn+")==0)  { B = 873.;  C = 3130.;  xi = 64.;  /*xi=64.;*/ n = 6; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("mn2+")==0) { B = 960.;  C = 3325.;  xi = 68.57;/*xi=68.57*/n = 5; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("mn3+")==0) { B = 1140.; C = 3675.;  xi = 88.;  /*xi=87.;*/ n = 4; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("mn4+")==0) { F[1]=87044;F[2]=54316; xi = 134.; /*xi=135;*/ n = 3; l=D; flg3d=1;          } // (Expt.) Uylings et al., J. Phys. B. 17 (1984) 4103
   else if(IONCMP("fe")==0)   { B = 806.;  C = 3506.;  xi = 0;    /*xi=   ;*/ n = 8; l=D; flg3d=1; flgBC=1; } //
   else if(IONCMP("fe+")==0)  { B = 869.;  C = 3638.;  xi = 119.; /*xi=115;*/ n = 7; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("fe2+")==0) { B = 1058.; C = 3091.;  xi = 103.; /*xi=114;*/ n = 6; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("fe3+")==0) { F[1]=97130;F[2]=60769; xi = 476.; /*xi=   ;*/ n = 5; l=D; flg3d=1;          } // Havercort Thesis 2p6 3d5
   else if(IONCMP("fe4+")==0) { B = 1144.; C = 4459.;  xi = 129.; /*xi=125;*/ n = 4; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("co")==0)   { B = 798.;  C = 4167.;  xi = 0;    /*xi=   ;*/ n = 9; l=D; flg3d=1; flgBC=1; } //
   else if(IONCMP("co+")==0)  { B = 878.;  C = 3828.;  xi = 228.; /*xi=228;*/ n = 8; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("co2+")==0) { B = 1115.; C = 4366.;  xi = 178.; /*xi=189;*/ n = 7; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("co3+")==0) { B = 1065.; C = 5120.;  xi = 128.6;/*xi=145;*/ n = 6; l=D; flg3d=1; flgBC=1; } // Abragam Bleaney 1970 p 391 for B,C. xi from PRB 67 172401
// else if(IONCMP("ni")==0)   { B = 1025.; C = 4226.;  xi = 0;    /*xi=   ;*/ n =10; l=D; flg3d=1; flgBC=1; } //
   else if(IONCMP("ni+")==0)  { B = 1037.; C = 4314.;  xi = 0;    /*xi=   ;*/ n = 9; l=D; flg3d=1; flgBC=1; } //
   else if(IONCMP("ni2+")==0) { B = 1084.; C = 4831.;  xi = 324.; /*xi=343;*/ n = 8; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("ni3+")==0) { B = 1184.; C = 5105.;  xi = 272.; /*xi= - ;*/ n = 7; l=D; flg3d=1; flgBC=1; } // B,C, from fit to NIST data using Racah II, eqn 84
   else if(IONCMP("ni4+")==0) { F[1]=100185;F[2]=64787;xi = 197.; /*xi= - ;*/ n = 6; l=D; flg3d=1;          } // (Expt.) Uylings et al., J. Phys. B. 17 (1984) 4103
// else if(IONCMP("cu+")==0)  { B = 1216.; C = 4745.;  xi = 0;    /*xi=   ;*/ n =10; l=D; flg3d=1; flgBC=1; } //
   else if(IONCMP("cu2+")==0) { B = 1238.; C = 4659.;  xi = 830.; /*xi=830;*/ n = 9; l=D; flg3d=1; flgBC=1; }
   else if(IONCMP("cu3+")==0) { F[1]=111996;F[2]=69924;xi = 903.; /*xi= - ;*/ n = 8; l=D; flg3d=1;          } // Thesis Havercort Koeln Cu 2p6 3d8
   else if(IONCMP("zn3+")==0) { F[1]=116868;F[2]=72923;xi = 1097.;/*xi=   ;*/ n = 9; l=D; flg3d=1;          } // Havercort Thesis 2p6 3d9
    
   // 4d ions parameters from Richardson, Blackman and Ranschak, J. Chem. Phys. v58, 3010 (1973).
   //   xi from calculations of Blume, Freeman, Watson, Phys. Rev. v134, A320 (1964), or where not calculated from TM Dunn, Trans. Faraday Soc. v57, 1441 (1961)
   else if(IONCMP("y2+")==0)  { B = 0.;    C = 0.;     xi = 312.; /*xi=300;*/ n = 1; l=D; flgBC=1; }
   else if(IONCMP("zr2+")==0) { B = 333.;  C = 3.96*B; xi = 432.; /*xi=425;*/ n = 2; l=D; flgBC=1; }
   else if(IONCMP("zr3+")==0) { B = 0.;    C = 0.;     xi = 507.; /*xi=500;*/ n = 1; l=D; flgBC=1; }
   else if(IONCMP("nb+")==0)  { B = 324.;  C = 3.89*B; xi = 0;    /*xi=   ;*/ n = 4; l=D; flgBC=1; }
   else if(IONCMP("nb2+")==0) { B = 360.;  C = 3.97*B; xi = 560.; /*xi=555;*/ n = 3; l=D; flgBC=1; }
   else if(IONCMP("nb3+")==0) { B = 391.;  C = 4.03*B; xi = 644.; /*xi=670;*/ n = 2; l=D; flgBC=1; }
   else if(IONCMP("nb4+")==0) { B = 0.;    C = 0.;     xi = 750.; /*xi=   ;*/ n = 1; l=D; flgBC=1; }
   else if(IONCMP("mo+")==0)  { B = 355.;  C = 3.91*B; xi = 630.; /*xi=   ;*/ n = 5; l=D; flgBC=1; }
   else if(IONCMP("mo2+")==0) { B = 387.;  C = 3.98*B; xi = 717.; /*xi=695;*/ n = 4; l=D; flgBC=1; }
   else if(IONCMP("mo3+")==0) { B = 416.;  C = 4.03*B; xi = 812.; /*xi=800;*/ n = 3; l=D; flgBC=1; }
   else if(IONCMP("mo4+")==0) { B = 440.;  C = 4.08*B; xi = 950.; /*xi=   ;*/ n = 2; l=D; flgBC=1; }
   else if(IONCMP("mo5+")==0) { B = 0.;    C = 0.;     xi =1030.; /*xi=   ;*/ n = 1; l=D; flgBC=1; }
   else if(IONCMP("tc2+")==0) { B = 414.;  C = 3.99*B; xi = 850.; /*xi=   ;*/ n = 5; l=D; flgBC=1; }
   else if(IONCMP("tc3+")==0) { B = 440.;  C = 4.03*B; xi = 990.; /*xi=   ;*/ n = 4; l=D; flgBC=1; }
   else if(IONCMP("tc4+")==0) { B = 0;     C = 0;      xi =1150.; /*xi=   ;*/ n = 3; l=D; flgBC=1; } //
   else if(IONCMP("ru2+")==0) { B = 436.;  C = 3.99*B; xi =1077.; /*xi=1000*/ n = 6; l=D; flgBC=1; }
   else if(IONCMP("ru3+")==0) { B = 464.;  C = 4.04*B; xi =1197.; /*xi=1180*/ n = 5; l=D; flgBC=1; }
   else if(IONCMP("ru4+")==0) { B = 400;   C = 1270;   xi =1350.; /*xi=   ;*/ n = 4; l=D; flgBC=1; } //
   else if(IONCMP("rh+")==0)  { B = 427.;  C = 3.93*B; xi =0;     /*xi=   ;*/ n = 8; l=D; flgBC=1; } //
   else if(IONCMP("rh2+")==0) { B = 458.;  C = 3.98*B; xi =1664.; /*xi=1640*/ n = 7; l=D; flgBC=1; }
   else if(IONCMP("rh3+")==0) { B = 484.;  C = 4.03*B; xi =1416.; /*xi=1400*/ n = 6; l=D; flgBC=1; }
   else if(IONCMP("rh4+")==0) { B = 0;     C = 0;      xi =1570.; /*xi=   ;*/ n = 5; l=D; flgBC=1; } //
   else if(IONCMP("pd+")==0)  { B = 451.;  C = 3.94*B; xi =0;     /*xi=   ;*/ n = 9; l=D; flgBC=1; } //
   else if(IONCMP("pd2+")==0) { B = 480.;  C = 3.99*B; xi =1529.; /*xi=1600*/ n = 8; l=D; flgBC=1; }
   else if(IONCMP("pd3+")==0) { B = 506.;  C = 4.03*B; xi =1529.; /*xi=1600*/ n = 7; l=D; flgBC=1; }
   else if(IONCMP("ag2+")==0) { B = 502.;  C = 3.99*B; xi =1794.; /*xi=1840*/ n = 9; l=D; flgBC=1; }
   else if(IONCMP("ag3+")==0) { B = 528.;  C = 4.03*B; xi =1940.; /*xi=1930*/ n = 8; l=D; flgBC=1; }
   // 5d ions parameters from G Burns, J. Chem. Phys. v41, 1521 (1964) B,C only.
   else if(IONCMP("re2+")==0) { F[1]=30870;F[2]=22050; xi =0.;    /*xi=   ;*/ n = 5; l=D; }
   
   else { std::cerr << "getfromionname(): Name of ion " << ionname << " not recognised.\n"; return; }

   if(flg3d) // 3d ion - "xi" is actually lambda parameter -> |lambda| = xi/2S
   {
      S2 = (n<=(2*l+1)) ? n : ((4*l+2)-n);    // Finds 2S (maximise S from Hund's Rules)
      xi *= (double)S2;
   }
   if(flgBC) // Need to convert B and C parameters to F^2, F^4 slater integrals
   {
      // Eqn 77 of Racah II, A = F_0-49F_4 = F^0-F^4/9; B = F_2-5F_4 = (9F^2-5F^4)/441; C = 35F_4 = 5F^4/63;
      F[2] = (63./5)*C; F[1] = (441.*B+5*F[2])/9.;
   }

   pars.n = n; pars.l = l; pars._F = F; pars._xi = xi; pars._alpha = a;
   if(pars.e_units.find("cm")!=std::string::npos || pars.e_units.find("wave")!=std::string::npos) {}
   // converts from cm^{-1} to meV
   else if(pars.e_units.find("meV")!=std::string::npos) { for(i=0;i<4;i++) F[i]/=MEV2CM; xi/=MEV2CM; for(i=0;i<3;i++) a[i]/=MEV2CM; } 
   // converts from cm^{-1} to K
   else if(pars.e_units.find("K")!=std::string::npos)   { for(i=0;i<4;i++) F[i]*=CM2K;   xi*=CM2K;   for(i=0;i<3;i++) a[i]*=CM2K; }   
   else std::cerr << "getfromionname(): Energy units " << pars.e_units << " not recognised. Accepted units are cm^{-1}, meV, K.\n";

   pars.F = F; pars.xi = xi; pars.alpha = a; pars.B.find_rk(ionname); pars.B.calc_stevfact(n,l); pars.B.convback();
}

void ic_parsecfpars(std::string &varname, std::string &varval, icpars &pars, int length)
{
   std::string tmpstr = varname.substr(length,1);
   int k = atoi(tmpstr.c_str()), q;
   if((int)varname.length()==(length+3))
   {
      if(varname.compare(length+1,1,"m")==0 || varname.compare(length+1,1,"-")==0) tmpstr = varname.substr(length+2,1);
      else if(varname.compare(length+2,1,"s")==0) tmpstr = varname.substr(length+1,1);
      else { std::cerr << "ic_parseinput: CF parameter name " << varname << " not recognised.\n"; return; }
      q = -atoi(tmpstr.c_str());
   }
   else { tmpstr = varname.substr(length+1,1); q = atoi(tmpstr.c_str()); }
   tmpstr = varname.substr(0,length);
   pars.B.assign(tmpstr,k,q,atof(varval.c_str()));
}

// --------------------------------------------------------------------------------------------------------------- //
// Parses input file for IC parameters
// --------------------------------------------------------------------------------------------------------------- //
void ic_parseinput(const char *filename, icpars &pars)
{
   std::string strline,varname,varval,tmpstr,issstr;
   size_t equal_sym_pos,first_nonspace_pos,first_num_pos=std::string::npos,last_num_pos=std::string::npos;
   std::istringstream iss;
   int k;

   std::fstream FILEIN; FILEIN.open(filename, std::fstream::in);
   if(FILEIN.fail()==true) { std::cerr << "ic_parseinput(): Cannot open file " << filename << "\n"; return; }

   // If the parameters file is from cfield, assume parameters are in meV, in Stevens normalisation, with theta_k=1.
   getline(FILEIN,strline);
   if(strline.compare(0,8,"#!cfield")==0)
   {
      tmpstr.assign("meV"); conv_e_units(pars,tmpstr); tmpstr.assign("B"); pars.B.conv(tmpstr);
   }

   while(!FILEIN.eof())
   {
      strline.clear(); getline(FILEIN,strline);
      if(strline.find("#",0,1)!=std::string::npos)              // Line is a comment, ignores it
         continue; 

      first_nonspace_pos = strline.find_first_not_of(" \t");
      equal_sym_pos = strline.find("=");
    
      varname.clear(); varval.clear();
      if(first_nonspace_pos==std::string::npos || equal_sym_pos==std::string::npos) 
         varname.assign(strline);
      else
      {
         iss.clear();
         varname = strline.substr(first_nonspace_pos,equal_sym_pos-first_nonspace_pos);
         varname = varname.substr(0,varname.find_first_of(" \t"));
         strtolower(varname);
         varval = strline.substr(equal_sym_pos+1);
         size_t subfirst = varval.find_first_not_of(" \t"); 
         size_t sublast  = varval.find_last_not_of(" \t")+1; 
         if(subfirst!=std::string::npos)
         {
            varval = varval.substr(subfirst,sublast);
            first_num_pos = varval.find_first_of("-0123456789.EDed");
            last_num_pos = varval.find_first_not_of("-0123456789.EDed",first_num_pos);
            if(first_num_pos!=std::string::npos)
               iss.str(varval.substr(first_num_pos,last_num_pos));
            else
               iss.str(varval);
         }
         else
            iss.str(varval);
      }

      if(varname.find("iontype")!=std::string::npos)
         getfromionname(varval,pars);                           // Looks up values of Fk and xi from tables
      else if(varname.find("unit")!=std::string::npos)
         conv_e_units(pars,varval);
      else if(varname.find("conf")!=std::string::npos)
      {  tmpstr = varval.substr(0,1); pars.l = Lin(tmpstr); tmpstr = varval.substr(1,2); pars.n = atoi(tmpstr.c_str()); }
      else if(varname.compare("n")==0)
      {  iss >> pars.n; pars.B.calc_stevfact(pars.n,pars.l); }
      else if(varname.compare("l")==0)
      {  tmpstr = varval.substr(0,1); pars.l = Lin(tmpstr); pars.B.calc_stevfact(pars.n,pars.l); }
      else if(varname.compare("f")==0)                          // Slater integrals expressed as four floats
         for(k=0; k<4; k++)
         {
            iss >> pars.F[k];
            first_num_pos = varval.find_first_of("0123456789.-EeDd",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.-EeDd",first_num_pos);
            iss.str(varval.substr(first_num_pos,last_num_pos));
         }
      else if(varname.compare("f0")==0) { iss >> pars.F[0]; pars._F[0]=pars.F[0]*pars._econv; }
      else if(varname.compare("f2")==0) { iss >> pars.F[1]; pars._F[1]=pars.F[1]*pars._econv; }
      else if(varname.compare("f4")==0) { iss >> pars.F[2]; pars._F[2]=pars.F[2]*pars._econv; }
      else if(varname.compare("f6")==0) { iss >> pars.F[3]; pars._F[3]=pars.F[3]*pars._econv; }
      else if(varname.compare("xi")==0) { iss >> pars.xi; pars._xi = pars.xi*pars._econv; }
      else if(varname.compare("zeta")==0) { iss >> pars.xi; pars._xi = pars.xi*pars._econv; }
      else if(varname.compare("alpha")==0)
      {
         if(varval.find(",")!=std::string::npos)
            for(k=0; k<3; k++)
            {
               iss >> pars.alpha[k]; pars._alpha[k] = pars.alpha[k]*pars._econv;
               first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
               last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
               iss.str(varval.substr(first_num_pos,last_num_pos));
            }
         else                              { iss >> pars.alpha[0]; pars._alpha[0]=pars.alpha[0]*pars._econv; }
      }
      else if(varname.compare("beta")==0)  { iss >> pars.alpha[1]; pars._alpha[1]=pars.alpha[1]*pars._econv; }
      else if(varname.compare("gamma")==0) { iss >> pars.alpha[2]; pars._alpha[2]=pars.alpha[2]*pars._econv; }

      else if(varname.find_first_of("awbvld")==0 && varname.find_first_of("0123456789")==1) 
         ic_parsecfpars(varname, varval, pars);
      else if(varname.compare(0,2,"ar")==0 && varname.find_first_of("0123456789")==2)
         ic_parsecfpars(varname, varval, pars, 2);
      else if(varname.find("density")!=std::string::npos)
         pars.density = varval;
//    else if(varname.find("observable")!=std::string::npos)
//       pars.observable = varval;
      else if(varname.find("eigenvectors")!=std::string::npos)
         iss >> pars.num_eigv;
      else if(varname.find("basis")!=std::string::npos)
         pars.basis = varval;
      else if(varname.find("calc")!=std::string::npos)
      {
         if(varname.find("mag")!=std::string::npos)  // Physical properties calculation flags.
            pars.calcphys += PHYSPROP_MAGBIT; 
         else if(varname.find("cp")!=std::string::npos || varname.find("heat")!=std::string::npos)
            pars.calcphys += PHYSPROP_CP_BIT; 
         else if(varname.find("susc")!=std::string::npos) {
            if(varname.find("inv")!=std::string::npos) 
               pars.calcphys += PHYSPROP_INVBIT;
            else 
               pars.calcphys += PHYSPROP_SUSBIT; }
      }
      else if(varname.find("partial_standalone")!=std::string::npos)
         pars.partial_standalone = true;
      else if(varname.find("partial")!=std::string::npos)
         pars.partial = true;
      #ifndef NO_ARPACK
      else if(varname.find("arnoldi_standalone")!=std::string::npos)
         pars.arnoldi_standalone = true;
      else if(varname.find("arnoldi")!=std::string::npos)
         pars.arnoldi = true;
      #endif
      else if(varname.find("save_matrices")!=std::string::npos)
         pars.save_matrices = true;
      else if(varname.find("truncate_matrix")!=std::string::npos)
      {  
         iss >> pars.truncate_level; if(pars.truncate_level<=0 || pars.truncate_level>=1) pars.truncate_level=0.5;
      }  
      else if(varname.find("use_J_operator_equivalent")!=std::string::npos)
      {  pars.B.op_equiv = Jt; pars.B.calc_stevfact(pars.n,pars.l); pars.B.convback(); }
      else if(varname.find("use_L_operator_equivalent")!=std::string::npos)
      {  pars.B.op_equiv = Lt; pars.B.calc_stevfact(pars.n,pars.l); pars.B.convback(); }
      else if(varname.find("dx2")!=std::string::npos) iss >> pars.Dx2;
      else if(varname.find("dy2")!=std::string::npos) iss >> pars.Dy2;
      else if(varname.find("dz2")!=std::string::npos) iss >> pars.Dz2;
      else if(varname.find("emu")!=std::string::npos)
         pars.mag_units = 1;
      else if(varname.find("simag")!=std::string::npos)
         pars.mag_units = 2;
      else if(varname.compare("xval")==0)
      {
         tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xT = atof(tmpstr.c_str());
         first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
         last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
         if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
         {
            tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xHa = atof(tmpstr.c_str());
            first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
            if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
            {
               tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xHb = atof(tmpstr.c_str());
               first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
               last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
               if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
               {
                  tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xHc = atof(tmpstr.c_str());
               }
            }
         }
      }
      else if(varname.compare("xrange")==0)
      {
         tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xMin = atof(tmpstr.c_str());
         first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
         last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	 if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
         {
            tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xStep = atof(tmpstr.c_str());
            first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	    if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
            {
               tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xMax = atof(tmpstr.c_str());
            }
         }
      }
      else if(varname.compare("yval")==0)
      {
         tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yT = atof(tmpstr.c_str());
         first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
         last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	 if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
         {
            tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yHa = atof(tmpstr.c_str());
            first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
            if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
            {
               tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yHb = atof(tmpstr.c_str());
               first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
               last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
               if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
               {
                  tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yHc = atof(tmpstr.c_str());
               }
            }
         }
      }
      else if(varname.compare("yrange")==0)
      {
         tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yMin = atof(tmpstr.c_str());
         first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
         last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	 if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
         {
            tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yStep = atof(tmpstr.c_str());
            first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	    if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
            {
               tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yMax = atof(tmpstr.c_str());
            }
         }
      }
      else if(varname.compare("xt")==0)    { tmpstr=iss.str(); pars.xT    = atof(tmpstr.c_str()); }
      else if(varname.compare("xha")==0)   { tmpstr=iss.str(); pars.xHa   = atof(tmpstr.c_str()); }
      else if(varname.compare("xhb")==0)   { tmpstr=iss.str(); pars.xHb   = atof(tmpstr.c_str()); }
      else if(varname.compare("xhc")==0)   { tmpstr=iss.str(); pars.xHc   = atof(tmpstr.c_str()); }
      else if(varname.compare("yt")==0)    { tmpstr=iss.str(); pars.yT    = atof(tmpstr.c_str()); }
      else if(varname.compare("yha")==0)   { tmpstr=iss.str(); pars.yHa   = atof(tmpstr.c_str()); }
      else if(varname.compare("yhb")==0)   { tmpstr=iss.str(); pars.yHb   = atof(tmpstr.c_str()); }
      else if(varname.compare("yhc")==0)   { tmpstr=iss.str(); pars.yHc   = atof(tmpstr.c_str()); }
      else if(varname.compare("xmin")==0)  { tmpstr=iss.str(); pars.xMin  = atof(tmpstr.c_str()); }
      else if(varname.compare("xstep")==0) { tmpstr=iss.str(); pars.xStep = atof(tmpstr.c_str()); }
      else if(varname.compare("xmax")==0)  { tmpstr=iss.str(); pars.xMax  = atof(tmpstr.c_str()); }
      else if(varname.compare("ymin")==0)  { tmpstr=iss.str(); pars.yMin  = atof(tmpstr.c_str()); }
      else if(varname.compare("ystep")==0) { tmpstr=iss.str(); pars.yStep = atof(tmpstr.c_str()); }
      else if(varname.compare("ymax")==0)  { tmpstr=iss.str(); pars.yMax  = atof(tmpstr.c_str()); }
      else if(varname.compare("bx")==0)    { tmpstr=iss.str(); pars.Bx    = atof(tmpstr.c_str()); }
      else if(varname.compare("by")==0)    { tmpstr=iss.str(); pars.By    = atof(tmpstr.c_str()); }
      else if(varname.compare("bz")==0)    { tmpstr=iss.str(); pars.Bz    = atof(tmpstr.c_str()); }
   }  // while(!FILEIN.eof())
   FILEIN.close();
}

