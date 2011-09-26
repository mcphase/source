/* coulomb.cpp
 *
 * Calculates the Coulomb (electrostatic) interaction operator matrices, after the method of Racah
 *
 * Functions:
 *   double racah_xwu(qR7 W, qG2 U, qG2 Up);                                    // Looks up values of x(W|UU')
 *   double racah_chi(orbital L, orbital Lp, qG2 U, qG2 Up);                    // Looks up values of (U|chi(L)|U')
 *   int racah_e2sign(int S2, int v);                                           // Looks up sign (+/-) of e2 (Racah 4, eq 73)
 *   double racah_e2prod(qR7 W, qG2 U, qG2 Up, orbital L, orbital Lp);          // Calculates sum of x(W|UU')*(U|chi(L)|U')
 *   sMat<double> racah_e2(int n);                                              // Calculates the matrix elements of e2
 *   double racah_yfn(int n, int v, int S2, qG2 U, int vp, qG2 Up);             // Looks up values of y(f^n,vSU|v'S'U')
 *   double racah_phi(qG2 U, qG2 Up, orbital Lp, orbital L);                    // Looks up values of (U|phi(L)|U')
 *   double racah_g(qG2 U, bool R5flag=false);                                  // Calcs. eigenvalue of Casimir's op for G2 or R5
 *   double racah_g(qR7 W);                                                     // Calcs. eigenvalue of Casimir's op for R7
 *   sMat<double> racah_e3(int n);                                              // Calculates the matrix elements of e3
 *   sMat<double> racah_emat(int n, double E0, double E1, double E2, double E3);//Calculates the Coulomb interaction matrix
 *   sMat<double> racah_emat(int n, double F0, double F2, double F4);           // Calculates the Coulomb matrix for d-electrons
 *   std::vector<double> racah_FtoE(std::vector<double> F);                     // Converts from F_k to E
 *   std::vector<double> racah_EtoF(std::vector<double> E);                     // Converts from E to F_k
 *   std::vector<double> racah_FtoF_k(std::vector<double> F);                   // Converts from F^k to F_k
 *   std::vector<double> racah_F_ktoF(std::vector<double> F_k);                 // Converts from F_k to F^k
 *   sMat<double> racah_ci(int n, double alpha, double beta, double gamma);     // Calcs. conf. interaction matrix for f-elec.
 *   sMat<double> racah_ci(int n, double alpha, double beta);                   // Calcs. conf. interaction matrix for d-elec.
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 *
 * Refs: Racah I   - G. Racah, Phys. Rev. 61, 186 (1942)
 *       Racah II  - G. Racah, Phys. Rev. 62, 438 (1942)
 *       Racah III - G. Racah, Phys. Rev. 63, 367 (1943)
 *       Racah IV  - G. Racah, Phys. Rev. 76, 1352 (1949)
 *       Nielson and Koster, Spectroscopic Coefficients for the $p^n$, $d^n$, and $f^n$ Configurations, MIT Press, NY (1963)
 *       Rajnak and Wybourne, Phys. Rev., 132, 280 (1963)
 */

#include "ic1ion.hpp"

// --------------------------------------------------------------------------------------------------------------- //
// Function to look up the value of the coefficients x(W,UU'), after Racah IV. 
// --------------------------------------------------------------------------------------------------------------- //
double racah_xwu(qR7 W, qG2 U, qG2 Up)
{
   double retval = 0;
   int id;

   if(W.isequal("200"))
   {
      if(U.isequal("20") && Up.isequal("20"))
         retval = 2;
   }
   else if (W.isequal("210"))
   {
      // Racah IV, Table VII, W = (210)
      // U =          (11)            (20)            (21)       U'
      double t[] = {    0,              0,  12*sqrt(455.),   // (11)
                        0,          -6./7,  6*sqrt(66.)/7,   // (20)
            12*sqrt(455.),  6*sqrt(66.)/7,              0 }; // (21)   x((210),(21)(21)) = [3/7, 0]

      id = ( 2*(Up.u1-1)+Up.u2-1 )*3 + ( 2*(U.u1-1)+U.u2-1 ); 
      retval = t[id];
   }
   else if (W.isequal("211"))
   {
      // Racah IV, Table VIII, W = (211)
      // U =          (10)            (11)            (20)            (21)            (30)       U'
      double t[] = {    0,              0,              0,              0, -20*sqrt(143.),   // (10)
                        0,              0,              0,  10*sqrt(182.),            10.,   // (11)
                        0,              0,          -8./7,  4*sqrt(33.)/7,     4*sqrt(3.),   // (20)
                        0,  10*sqrt(182.),  4*sqrt(33.)/7,              0,             2.,   // (21) [4/7 3]
           -20*sqrt(143.),            10.,     4*sqrt(3.),             2.,             2. }; // (30)
 
      id = ( 2*(Up.u1-1)+Up.u2 )*5 + ( 2*(U.u1-1)+U.u2 ); 
      retval = t[id];
   }
   else if (W.isequal("220"))
   {
      // Racah Table IX, W = (220)
      // U =          (20)           (21)           (22)        U'
      double t[] = {3./14, 3*sqrt(55.)/7, -3*sqrt(5./28),   // (20)
            3*sqrt(55.)/7,             0,     3./sqrt(7),   // (21)  x((220),(21)(21)) = [6/7 3]
          -3.*sqrt(5./28),    3./sqrt(7),           3./2};  // (22)

      id = ( 2*(Up.u1-1)+Up.u2-2 )*3 + ( 2*(U.u1-1)+U.u2-2 ); 
      retval = t[id];
   }
   else if (W.isequal("221"))
   {
      // Racah IV, Table X, W = (221)
      // U =          (10)             (11)            (20)             (21)             (30)            (31)       U'
      double t[] = {    0,               0,              0,               0,     5*sqrt(143.),  -15.*sqrt(429),   // (10)
                        0,               0,              0,14*sqrt(910./11),      2*sqrt(10.),  2.*sqrt(39)/11,   // (11)
                        0,               0,           2./7,  -10*sqrt(6.)/7,         sqrt(3.),    9*sqrt(3./7),   // (20)
                        0,14*sqrt(910./11), -10*sqrt(6.)/7,               0,    5*sqrt(2./11),   3.*sqrt(2)/11,   // (21) [-1/7 12/11]
             5*sqrt(143.),     2*sqrt(10.),       sqrt(3.),   5*sqrt(2./11),            -1./2,   3./2/sqrt(11),   // (30)
           -15*sqrt(429.),  2*sqrt(39.)/11,   9*sqrt(3./7),   3*sqrt(2.)/11,    3./2/sqrt(11),           1./22};  // (31)

      id = ( 2*(Up.u1-1)+Up.u2 )*6 + ( 2*(U.u1-1)+U.u2 ); 
      retval = t[id];
   }
   else if (W.isequal("222"))
   {
      // Racah IV, Table XI, W = (222)
      // U =          (00)            (10)            (20)            (30)            (40)       U'
      double t[] = {    0,              0,              0,              0, -30.*sqrt(143),   // (00)
                        0,              0,              0, -3.*sqrt(1430),  9.*sqrt(1430),   // (10)
                        0,              0,          6./11,-3*sqrt(42./11),  9.*sqrt(2)/11,   // (20)
                        0, -3*sqrt(1430.),-3*sqrt(42./11),            -3.,    1./sqrt(11),   // (30)
           -30*sqrt(143.),  9*sqrt(1430.),  9.*sqrt(2)/11,     1/sqrt(11),          3./11};  // (40)

      id = Up.u1*5 + U.u1;
      retval = t[id];
   }
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to look up the value of (U|chi(L)|U') after Racah IV.
// --------------------------------------------------------------------------------------------------------------- //
double racah_chi(orbital Lp, orbital L, qG2 U, qG2 Up)
{
   double retval = 0;
   int id;

   if(abs(L)!=abs(Lp)) { return retval; }     // Numerical values of L must equal - but can have L and Lp different

   if(Up.isequal("31"))      // Use table VIb
   {
      if(L>0) { id = (L-1) + ( L>3 ? 1 : 0 ) + ( L>5 ? 1 : 0 ) + ( L>6 ? 1 : 0 ) + ( L>7 ? 1 : 0 ); }
      else    { id = (L==-3 ? 3 : 0) + (L!=-3 ? -2*L-4 : 0 ); }
      //    (10)            (11)            (20)            (21)             (30)            (31)               (31)'   L Lnum id L-1 diff
      //   ------------------------------------------------------------------------------------------------------------------------------
      double t[] = {0,11.*sqrt(330),           0,              0,   76.*sqrt(143),         -6644.,                 0, // P   1  0   0  0
              0,               0,   -8.*sqrt(78),-60*sqrt(39./7),               0,          4792.,                 0, // D   2  1   1  0
              0,               0,              0,  -312.*sqrt(5),   -48.*sqrt(39),          4420.,    336.*sqrt(143), // F   3  2   2  0
             1.,               0,              0,  12.*sqrt(715),   -98.*sqrt(33),          -902.,    336.*sqrt(143), // F' -3  3
              0,               0,    5.*sqrt(65),   2024/sqrt(7),  20.*sqrt(1001),         -2684.,                 0, // G   4  4   3  1
              0,    11.*sqrt(85),             0,31*sqrt(1309./3),  -20.*sqrt(374),         -2024.,   -48.*sqrt(6545), // H   5  5   4  1
              0,   -25.*sqrt(77),              0, 103*sqrt(5./3),   -44.*sqrt(70),          2680.,   -48.*sqrt(6545), // H' -5  6
              0,               0,   10.*sqrt(21),              0,   -57.*sqrt(33),      -12661./5, -3366.*sqrt(34)/5, // I   6  7   5  2
              0,               0,              0,              0,  18.*sqrt(1122),       17336./5, -3366.*sqrt(34)/5, // I' -6  8
              0,               0,              0,-52*sqrt(323./23),-494*sqrt(19./23), 123506./23,144*sqrt(21318.)/23, // K   7  9   6  3
              0,               0,              0,-336*sqrt(66./23),73*sqrt(1122./23), -85096./23,144*sqrt(21318.)/23, // K' -7 10
              0,               0,              0,  -24*sqrt(190),               0,         -4712.,                 0, // L   8 11   7  4
              0,               0,              0,              0,  -21.*sqrt(385),          -473.,                 0, // M   9 12   8  4
              0,               0,              0,              0,               0,          1672.,                 0, // N  10 13   9  4
              0,               0,              0,              0,               0,           220.,                 0};// O  11 14  10  4

      if(U.isequal("10")) { retval = t[id*7]; }
      else if(U.isequal("11")) { retval = t[id*7+1]; }
      else if(U.isequal("20")) { retval = t[id*7+2]; }
      else if(U.isequal("21")) { retval = t[id*7+3]; }
      else if(U.isequal("30")) { retval = t[id*7+4]; }
      else if(U.isequal("31")) { 
         if(L==Lp) { retval = t[id*7+5]; }
         else      { retval = t[id*7+6]; } }
   }
   else if(Up.isequal("40")) // Use table VIc
   {
#define CD(c,d) case c: id = d; break
      switch (L) { CD(0,0); CD(2,1); CD(3,2); CD(4,3); CD(-4,4); CD(5,5); CD(6,6); 
                   CD(-6,7); CD(7,8); CD(8,9); CD(-8,10); CD(9,11); CD(10,12); CD(12,13); default: return retval; }
      //    (00)            (10)             (20)            (30)            (40)              (40)'   L Lnum id
      //    ----------------------------------------------------------------------------------------------------
      double t[] = {1.,       0,                0,               0,        -1408.,                 0, // S   0   0
              0,              0,    -88.*sqrt(13),               0,          -44.,                 0, // D   2   1
              0,             1.,                0,    90.*sqrt(11),         1078.,                 0, // F   3   2
              0,              0, 53*sqrt(715./27), -16.*sqrt(1001),     -16720./9, -34.*sqrt(2618)/9, // G   4   3
              0,              0,7*sqrt(15470./27),   64.*sqrt(442),      10942./9, -34.*sqrt(2618)/9, // G' -4   4
              0,              0,                0,  -72.*sqrt(462),         -704.,                 0, // H   5   5
              0,              0,34*sqrt(1045./31),-9*sqrt(21945./31),   -2453./31,60.*sqrt(74613)/31, // I   6   6
              0,              0,-12*sqrt(1785./31),756*sqrt(85./31),    36088./31,60.*sqrt(74613)/31, // I' -6   7
              0,              0,                0,   -84.*sqrt(33),         -132.,                 0, // K   7   8
              0,              0,                0,               0,     -4268./31,924.*sqrt(1995)/31, // L   8   9
              0,              0,                0,               0,     11770./31,924.*sqrt(1995)/31, // L' -8  10
              0,              0,                0,   -99.*sqrt(15),        -1067.,                 0, // M   9  11
              0,              0,                0,               0,          528.,                 0, // N  10  12
              0,              0,                0,               0,           22.,                 0};// Q  12  13
      
      if(U.isequal("40")) { if(L==Lp) { retval = t[id*6+4]; } else { retval = t[id*6+5]; } }
      else { retval = t[id*6+U.u1]; }
   }
   else                      // Use table VIa 
   {
      //      S        P        D       F        G       H        I        K        L        M       N   // U    Up
      //    ---------------------------------------------------------------------------------------------------
      double t[] = {0 ,0,    143.,      0,   -130.,      0,     35.,       0,       0,       0,      0,  // 20   20 
              0,       0,       0,      0,       0,     1.,       0,       0,       0,       0,      0,  // 11   21 
              0,       0,-39.*sqrt(2),  0, 4.*sqrt(65),  0,       0,       0,       0,       0,      0,  // 20   21 
              0,       0,    377.,   455.,   -561.,    49.,       0,   -315.,    245.,       0,      0,  // 21   21 
              0,       0,     13.,   -65.,     55.,   -75.,       0,    133.,    -75.,       0,      0,  // 21   21 
              0,       0,       0,     1.,       0,      0,       0,       0,       0,       0,      0,  // 10   30 
              0, -13.*sqrt(11), 0,      0,       0,  sqrt(39),    0,       0,       0,       0,      0,  // 11   30 
              0,       0,       0,      0,-13.*sqrt(5),  0,     30.,       0,       0,       0,      0,  // 20   30 
              0,       0,       0,12.*sqrt(195),8.*sqrt(143),11.*sqrt(42), 0,-4.*sqrt(17),0, 0,      0,  // 21   30 
              0,    -52.,       0,    38.,    -52.,    88.,     25.,    -94.,       0,     25.,      0,  // 30   30 
              0,       0, 3.*sqrt(429), 0,-38.*sqrt(65), 0, 21.*sqrt(85),  0,       0,       0,      0,  // 20   22  
              0,       0, 45.*sqrt(78), 0, 12.*sqrt(11),-12.*sqrt(546), 0, 0,-8.*sqrt(665),  0,      0,  // 21   22 
           260.,       0,    -25.,      0,     94.,   104.,   -181.,       0,    -36.,       0,    40.}; // 22   22 

      if(L>=0) {
      if(Up.isequal("20") && U.isequal("20")) { retval = t[L]; }
      else if(Up.isequal("21")) {
         if(U.isequal("11")) { retval = t[L+11]; }
	 else if(U.isequal("20")) { retval = t[L+22]; }
	 else if(U.isequal("21")) { retval = 0; } }
      else if(Up.isequal("30")) {
         if(U.isequal("10")) { retval = t[L+55]; }
	 else if(U.isequal("11")) { retval = t[L+66]; }
	 else if(U.isequal("20")) { retval = t[L+77]; }
	 else if(U.isequal("21")) { retval = t[L+88]; }
	 else if(U.isequal("30")) { retval = t[L+99]; } }
      else if(Up.isequal("22")) {
         if(U.isequal("20")) { retval = t[L+110]; }
	 else if(U.isequal("21")) { retval = t[L+121]; }
	 else if(U.isequal("22")) { retval = t[L+132]; }
      }}
   }
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Looks up the sign of the matrix element of e2, after Racah IV, eqn 73 and table I.
// --------------------------------------------------------------------------------------------------------------- //
int racah_e2sign(int S2, int v)
{
   int s = 0;
   switch(S2)
   {
      case 0: if(v==0 || v==2 || v==4 || v==6) { s = 1; } break;
      case 1: if(v==1 || v==3 || v==5) { s = 1; } else if(v==7) { s = -1; } break;
      case 2: if(v==2 || v==4) { s = 1; } else if (v==6) { s = -1; } break;
      case 3: if(v==3) { s = 1; } else if (v==5 || v==7) { s = -1; } break;
      case 4: if(v==4 || v==6) { s = -1; } break;
      case 5: if(v==5 || v==7) { s = -1; } break;
      case 6: if(v==6) { s = -1; } break;
      case 7: if(v==7) { s = -1; } break;
   }
   return s;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the sum of the values of x(W,UU')*(U|chi(L)|U') [there are two values for U=U'=(21)]
// --------------------------------------------------------------------------------------------------------------- //
double racah_e2prod(qR7 W, qG2 U, qG2 Up, orbital L, orbital Lp)
{
   double x1=0,x2=0,c1,c2;
   double retval = 0;
   if(U.isequal("21") && Up.isequal("21"))
   {
      if(W.isequal("210")) { x1 = (3./7); x2 = 0.; }
      else if(W.isequal("211")) { x1 = (4./7);  x2 = 3.; }
      else if(W.isequal("220")) { x1 = (-6./7); x2 = -3.; }
      else if(W.isequal("221")) { x1 = (-1./7); x2 = (12./11); }
      switch(L)
      {
         case 0: c1= 0;    c2= 0;    break;
	 case 1: c1= 0;    c2= 0;    break;
	 case 2: c1= 377.; c2= 13.;  break;
	 case 3: c1= 455.; c2=-65.;  break;
	 case 4: c1=-561.; c2= 55.;  break;
	 case 5: c1= 49.;  c2=-75.;  break;
	 case 6: c1= 0;    c2= 0;    break;
	 case 7: c1=-315.; c2= 133.; break;
	 case 8: c1= 245.; c2=-75.;  break;
	 case 9: c1= 0;    c2= 0;    break;
	 case 10: c1=0;    c2= 0;    break;
	 default: c1=0;    c2= 0;
      }
      retval = x1*c1 + x2*c2;
   }
   else
   {
      retval = racah_xwu(W,U,Up)*racah_chi(L,Lp,U,Up);
   }
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the matrix elements of the coulomb operator e2, after Racah IV, section 6.2
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_e2(int n)
{
   int i,j,sign,S2,v;
   fconf conf(n);
   int num_states = (int)conf.states.size();
   sMat<double> e2(num_states,num_states);
   double elem;
   
   for(j=0; j<num_states; j++)
   {
      S2 = conf.states[j].S2; v = conf.states[j].v;
      sign = racah_e2sign(S2,v);
      for(i=0; i<=j; i++)
      {
         if(abs(conf.states[i].L)==abs(conf.states[j].L) && S2==conf.states[i].S2 && v==conf.states[i].v)
	 {
	    elem = racah_e2prod(conf.states[i].W, conf.states[i].U, conf.states[j].U, conf.states[i].L, conf.states[j].L);
	    e2(i,j) = (double)sign * elem;
	    if(i!=j) { e2(j,i) = (double)sign * elem; }
	 }
      }
   }
   return e2;
}

// --------------------------------------------------------------------------------------------------------------- //
// Looks up the value of the coefficient y(f^n, vSU, v'SU') from tables in Racah IV.
// --------------------------------------------------------------------------------------------------------------- //
double racah_yfn(int n, int v, int S2, qG2 U, int vp, qG2 Up)
{
   double y = 0;
   int id;
   
   switch(n)
   {
      case 2: if(S2==0 && vp==2 && v==2 && U.isequal("20") && Up.isequal("20")) { y = 2; } break;
      case 3: if(S2==1 && vp==3)
              {  // Table XV from Racah IV.
	         // Up = (11)         (20)               (21)       U
                 double t[] = {0,       0,       -6*sqrt(22.),  // (10)
                    2.,                 0,                  0,  // (11)
                    0,              10./7,      2*sqrt(66.)/7,  // (20)
                    0,      2*sqrt(66.)/7,               2./7}; // (21)
		 id = ( 2*(U.u1-1)+U.u2 )*3 + ( 2*(Up.u1-1)+Up.u2-1 );
		 y = t[id];
	      } break;
      case 4: if(vp==4)
              {  //  Table XVI and XVII from Racah IV.
		 //  Sp =     1       1       1          1        1           0            0             0         
	         //  Up =   (10)    (11)    (20)       (21)     (30)        (20)         (21)          (22)     //  v S  U 
                 double t[] = {0,     0,      0, -12*sqrt(33./5), 0,          0,           0,            0,     //  2 1 (10)
                              0,   6./5,      0,         0,      6.,          0,           0,            0,     //  2 1 (11)
                              0,      0,      0, 8*sqrt(11./15),  0,          0,           0,            0,     //  4 1 (10)
                              0, 29./15,      0,         0,   -1./3,          0,           0,            0,     //  4 1 (11)
                              0,      0,   6./7, -8*sqrt(11./147), 4./sqrt(3),0,           0,            0,     //  4 1 (20)
                 8*sqrt(11./15),      0,-8*sqrt(11./147),-2./21, -4./3,       0,           0,            0,     //  4 1 (21)
                              0,  -1./3, 4./sqrt(3), -4./3,    1./3,          0,           0,            0,     //  4 1 (30)
                              0,      0,      0,         0,       0,          0,           0, -12.*sqrt(22),    //  0 0 (00)
                              0,      0,      0,         0,       0, 3*sqrt(3./175), -4*sqrt(33./35), -sqrt(3./5),//2 0 (20)
                              0,      0,      0,         0,       0,   221./140, 8*sqrt(11./245), -sqrt(7./80), //  4 0 (20)
                              0,      0,      0,         0,       0,8*sqrt(11./245),    2./7,            0,     //  4 0 (21)
                              0,      0,      0,         0,       0,  -sqrt(7./80),        0,          1./4};   //  4 0 (22)
	         if(S2==2) { if(v==2) { id = U.u2*8 + 2*(Up.u1-1)+Up.u2; y = t[id]; }
		             else if(v==4) { id = (2*(U.u1-1)+U.u2+2)*8 + 2*(Up.u1-1)+Up.u2; y = t[id]; } }
	         else if(S2==0) { if(v==0 && U.isequal("00")) { y = t[Up.u2+61]; } 
		             else if(v==2 && U.isequal("20")) { y = t[Up.u2+69]; }
		             else if(v==4) { id = (U.u2+9)*8 + Up.u2+5; y = t[id]; } }
	      } break;
      case 5: if(vp==5)
              {  
	         if(S2==3)
	         {  //  Table XVIII from Racah IV. vp=5, Sp=3/2
                    //  Up =   (10)  (11)        (20)          (21)           (30)     v 2S  U
                    double t[] = {0,  0,           0,            0,             0,  // 3  3 (00)
                                0,    0,           0,  9.*sqrt(11),             0,  // 3  3 (10)
                                0,    0,  3./sqrt(7),  sqrt(33./7),  -2.*sqrt(21),  // 3  3 (20)
                                0,    0,           0, -sqrt(55./3),             0,  // 5  3 (10)
                                0, -1./3,          0,            0,         -5./3,  // 5  3 (11)
                                0,    0,        5./7, 5*sqrt(11./147), 2./sqrt(3),  // 5  3 (20)
                     -sqrt(55./3),    0, 5*sqrt(11./147),   -4./21,         -2./3,  // 5  3 (21)
                                0, -5./3, 2./sqrt(3),        -2./3,         -1./3}; // 5  3 (30)
                    if(v==3) { id = U.u1*5 + 2*(Up.u1-1)+Up.u2; y = t[id]; }
		    else if(v==5) { id = (2*(U.u1-1)+U.u2+3)*5 + 2*(Up.u1-1)+Up.u2; y = t[id]; }
                 }
		 else if(S2==1) 
                 {  //   Table XIX from Racah IV. vp=5, Sp=1/2
                    //   Up =  (10)  (11)       (20)             (21)            (30)            (31)    v 2S  U
                    double t[] = {0,  0,          0,      36./sqrt(5),             0,   -36.*sqrt(2), // 1  2 (10)
                                0, 3./sqrt(2),    0,                0,  3.*sqrt(5)/2,   -sqrt(39./8), // 3  2 (11)
                                0,    0,       3./7,   -11.*sqrt(6)/7,   -4.*sqrt(3),              0, // 3  2 (20)
                   3*sqrt(33./10),    0, 3*sqrt(33./98), 3./7/sqrt(11), -3./2/sqrt(2), 3./2/sqrt(22), // 3  2 (21)
                                0,    0,          0,     43./sqrt(30),             0,     4.*sqrt(3), // 5  2 (10)
                                0, -5./6,         0,                0,-5*sqrt(5./72),  -sqrt(13./48), // 5  2 (11)
                                0,    0,      11./7,   -11./7/sqrt(6),    4./sqrt(3),              0, // 5  2 (20)
                     43./sqrt(30),    0, 11./7/sqrt(6),      25./231, 29./6/sqrt(22),  1./22/sqrt(2), // 5  2 (21)
                                0, -5*sqrt(5./72), 4./sqrt(3), 29./6/sqrt(22), -1./12, 1./4/sqrt(11), // 5  2 (30)
                       4.*sqrt(3), -sqrt(13./48), 0,   1./22/sqrt(2), 1./4/sqrt(11),          1./44}; // 5  2 (31)
		    if(v==1) { y = t[2*(Up.u1-1)+Up.u2]; }
		    else if(v==3) { id = (2*(U.u1-1)+U.u2)*6 + 2*(Up.u1-1)+Up.u2; y = t[id]; }
		    else if(v==5) { id = (2*(U.u1-1)+U.u2+4)*6 + 2*(Up.u1-1)+Up.u2; y = t[id]; }
		 }
	      } break;
      case 6: if(vp==6)
              {
	         if(S2==4)
		 {  //   Table XX from Racah IV. vp=6, Sp=2, v=4, S=2
		    //   Up =    (11)          (20)             (21)      U
                    double t[] = {0,             0,               0,  // (00)
                                  0,             0,    -6.*sqrt(11),  // (10)
                                  0, -2*sqrt(2./7),   2*sqrt(33./7)}; // (20)
		    if(v==4) { id = U.u1*3 + 2*(Up.u1-1)+Up.u2-1; y = t[id]; }
                 }
	         if(S2==2)
		 {  //   Table XXI from Racah IV. vp=6, Sp=1
		    //   Up =    (10)      (11)           (20)             (21)           (30)           (31)     v S  U
                    double t[] = {0,         0,             0,  -48*sqrt(2./5),             0,          -36.,  // 2 1 (10)
                                0,  sqrt(6./5),             0,               0,      sqrt(3.),3*sqrt(13./10),  // 2 1 (11)
                                0,           0,             0,    46./sqrt(15),             0,   -8.*sqrt(6),  // 4 1 (10)
                                0, 11./3/sqrt(5),           0,               0,-19./3/sqrt(2),  sqrt(13./60),  // 4 1 (11)
                                0,           0, -6.*sqrt(2)/7,  -22./7/sqrt(3),  8*sqrt(2./3),             0,  // 4 1 (20)
                    -sqrt(110./3),           0, sqrt(22./147),-16./21/sqrt(11),  5./3/sqrt(2),   1./sqrt(22),  // 4 1 (21)
                                0, -sqrt(5.)/3,  4*sqrt(2./3),   4./3/sqrt(11),  1./3/sqrt(2),  -1./sqrt(22)}; // 4 1 (30)
		    if(v==2) { id = U.u2*6 + 2*(Up.u1-1)+Up.u2; y = t[id]; }
		    else if(v==4) { id = (2*(U.u1-1)+U.u2+2)*6 + 2*(Up.u1-1)+Up.u2; y = t[id]; }
                 }
	         if(S2==0)
		 {  //   Table XXII from Racah IV. vp=6, Sp=0
		    //   Up =    (00)       (10)            (20)           (30)           (40)     v S  U
                    double t[] = {0,          0,    6./sqrt(55), 2*sqrt(42./5), 6*sqrt(2./55),  // 2 0 (20)
                                 0,           0, -61./sqrt(770),  8*sqrt(3./5), -6./sqrt(385),  // 4 0 (20)
                                 0, 3.*sqrt(22),     sqrt(2./7),     -sqrt(3.),    1./sqrt(7),  // 4 0 (21)
                    -4*sqrt(33./5),           0,   -1./sqrt(22),             0,   2./sqrt(11)}; // 4 0 (22)
		    if(v==2) { y = t[Up.u1]; }
		    else if(v==4) { id = (U.u2+1)*5 + Up.u1; y = t[id]; }
                 }
	      } break;
      case 7: if(vp==7)
              {
	         if(S2==3)
		 {  //   Table XXIII from Racah IV. vp=7, Sp=3/2, v=3, S=3/2
		    //   Up =   (20)           (21)            (22)     U
                    double t[] = {0,             0, -12.*sqrt(11),  // (00)
                                  0,   6.*sqrt(33),             0,  // (10)
                        -sqrt(5./7), 2*sqrt(11./7),           -1.}; // (20)
	            if(v==3) { y = t[U.u1*3 + Up.u2]; }
		 }
		 else if(S2==1)
		 {  //   Table XXIV from Racah IV. vp=7, Sp=1/2, v=3, S=1/2
		    //   Up =   (00)        (10)           (20)         (30)           (40)      U
                    double t[] = {0,          0,             0, 2.*sqrt(10),             0,  // (11)
                                  0,          0, -16./sqrt(77), -2.*sqrt(6), 6*sqrt(2./77),  // (20)
                                  0, -sqrt(66.),    sqrt(6./7),          1.,    sqrt(3./7)}; // (21)
		    if(v==3) { y = t[(2*(U.u1-1)+U.u2-1)*5 + Up.u1]; }
		 }
	      } break;
   }
   return y;
}

double racah_phi(qG2 U, qG2 Up, orbital Lp, orbital L)
{
   double retval=0;
   int id;

   if(Up.isequal("31"))
   {  //  Table XIVb from Racah IV. U' = (31)
      //  U = (10)   (11)         (21)         (30)              (31)               (31)'      L  Lnum id
      double t[]={0,sqrt(330.),     0,    17*sqrt(143.),         209.,                 0,  //  P   1   0
              0,      0,   12*sqrt(273.),        0,             -200.,                 0,  //  D   2   1
              1.,     0,     -36*sqrt(5.),  -16*sqrt(39.),       624.,    -80*sqrt(143.),  //  F   3   2
              0,      0,    -3*sqrt(715.),   24*sqrt(33.),      -616.,    -80*sqrt(143.),  //  F' -3   3
              0,      0,      11*sqrt(7.),    4*sqrt(1001.),     836.,                 0,  //  G   4   4
              0,  sqrt(85.),-2*sqrt(1309./3), sqrt(187./2),  -1353./2,  -5*sqrt(6545.)/2,  //  H   5   5
              0,  sqrt(77.), -74*sqrt(5./3), 31*sqrt(35./2),   703./2,  -5*sqrt(6545.)/2,  //  H' -5   6
              0,      0,            0,     30*sqrt(33.),     -2662./5,   528*sqrt(34.)/5,  //  I   6   7
              0,      0,            0,           0,            -88./5,   528*sqrt(34.)/5,  //  I' -6   8
              0,      0,-28*sqrt(323./23),   4*sqrt(437.),   6652./23,96*sqrt(21318.)/23,  //  K   7   9
              0,      0,  42*sqrt(66./23),       0,         -5456./23,96*sqrt(21318.)/23,  //  K' -7  10
              0,      0,    -6*sqrt(190.),       0,             -464.,                 0,  //  L   8  11
              0,      0,            0,    -6*sqrt(385.),         814.,                 0,  //  M   9  12
              0,      0,            0,           0,             -616.,                 0,  //  N  10  13
              0,      0,            0,           0,              352.,                 0}; //  O  11  14
      switch (L) { CD(1,0); CD(2,1); CD(3,2); CD(-3,3); CD(4,4); CD(5,5); CD(-5,6); CD(6,7); 
                   CD(-6,8); CD(7,9); CD(-7,10); CD(8,11); CD(9,12); CD(10,13); CD(11,14); default: return retval; }
      if(U.isequal("31")) { if(L==Lp) { retval = t[id*6+4]; } else { retval = t[id*6+5]; } }
      else if(U.isequal("30")) { retval = t[id*6+3]; } 
      else if(U.isequal("21")) { retval = t[id*6+2]; }
      else if(U.isequal("11")) { retval = t[id*6+1]; }
      else if(U.isequal("10")) { retval = t[id*6]; }
   }
   else if(Up.isequal("40"))
   {  //  Table XIVc from Racah IV. Up=(40)
      //  U =       (20)               (21)               (22)       L Lnum id
      double t[] = {   0,                 0,     2*sqrt(2145.),  //  S   0   0
            11*sqrt(13.),      -6*sqrt(26.),       9*sqrt(33.),  //  D   2   1
                       0,      3*sqrt(455.),                 0,  //  F   3   2
        -4*sqrt(715./27), -131*sqrt(11./27),   -4*sqrt(11./27),  //  G   4   3
         sqrt(15470./27),  17*sqrt(238./27), -17*sqrt(238./27),  //  G' -4   4
                       0,     -12*sqrt(21.),      3*sqrt(286.),  //  H   5   5
        7*sqrt(1045./31),                 0,  3*sqrt(3553./31),  //  I   6   6
        3*sqrt(1785./31),                 0,   75*sqrt(21./31),  //  I' -6   7
                       0,     -2*sqrt(119.),                 0,  //  K   7   8
                       0,  22*sqrt(105./31),   4*sqrt(627./31),  //  L   8   9
                       0,  -84*sqrt(19./31),  12*sqrt(385./31),  //  L' -8  10
                       0,  -84*sqrt(19./31),       sqrt(2530.)}; //  N  10  11
      switch (L) { CD(0,0); CD(2,1); CD(3,2); CD(4,3); CD(-4,4); CD(5,5); CD(6,6); 
                   CD(-6,7); CD(7,8); CD(8,9); CD(-8,10); CD(10,11); default: return retval; }
      if(U.isequal("20")) { retval = t[id*3]; }
      else if(U.isequal("21")) { retval = t[id*3+1]; }
      else if(U.isequal("22")) { retval = t[id*3+2]; }
   }
   else if(L>=0)
   {  //  Table XIVa from Racah IV.
      //  L =     S    P       D        F         G        H        I        K       L     M     N        U   U'
      double t[]={0, -11.,     0,       0,        0,      3.,       0,       0,      0,    0,    0,  // (11) (11)
              0,      0,    -11.,       0,      -4.,       0,      7.,       0,      0,    0,    0,  // (20) (20)
              // Added following check with Wybourne tables in J. Chem. Phys. v31, p340, 1961
              0,      0, -6*sqrt(2.),   0, -sqrt(65.),     0,       0,       0,      0,    0,    0,  // (21) (20)
	      // Original table resumes
              0,      0,      0,       1.,        0,       0,       0,       0,      0,    0,    0,  // (10) (21)
              0,      0, 6*sqrt(2.),    0,  sqrt(65.),     0,       0,       0,      0,    0,    0,  // (20) (21)
              0,      0,    -57.,     63.,      55.,   -105.,       0,    -14.,    42.,    0,    0,  // (21) (21)
              0, sqrt(11.),   0,        0,        0,   sqrt(39.),   0,       0,      0,    0,    0,  // (11) (30)
              0,      0,      0,        0,  2*sqrt(5.),    0,      3.,       0,      0,    0,    0,  // (20) (30)
              0,      0,      0, sqrt(195.),-sqrt(143.),-2*sqrt(42.),0,-4*sqrt(17.), 0,    0,    0,  // (21) (30)
              0,     83.,     0,     -72.,      20.,    -15.,     42.,    -28.,      0,   6.,    0,  // (30) (30)
             1.,      0,      0,        0,        0,       0,       0,       0,      0,    0,    0,  // (00) (22)
              0,      0,3*sqrt(429.),   0, 4*sqrt(65.),    0, 3*sqrt(85.),   0,      0,    0,    0,  // (20) (22)
           144.,      0,     69.,       0,     -148.,    72.,     39.,       0,   -96.,    0,  56.}; // (22) (22)
      if(U.isequal("11") && Up.isequal("11")) { retval = t[L]; }
      else if(Up.isequal("20")) 
      { 
         if (U.isequal("20")) { retval = t[L+11]; } 
	 else if(U.isequal("21")) { retval = t[L+22]; }
      }
      else if(Up.isequal("10") && U.isequal("21")) { retval = t[L+33]; }
      else if(Up.isequal("21")) 
      { 
         if (U.isequal("10")) { retval = t[L+33]; } 
	 else if(U.isequal("20")) { retval = t[L+44]; }
	 else if(U.isequal("21")) { retval = t[L+55]; }
      }
      else if(Up.isequal("30")) 
      { 
         if (U.isequal("11")) { retval = t[L+66]; }
	 else if(U.isequal("20")) { retval = t[L+77]; }
	 else if(U.isequal("21")) { retval = t[L+88]; }
	 else if(U.isequal("30")) { retval = t[L+99]; }
      }
      else if(Up.isequal("22")) 
      { 
         if (U.isequal("00")) { retval = t[L+110]; } 
         else if(U.isequal("20")) { retval = t[L+121]; }
	 else if(U.isequal("22")) { retval = t[L+132]; }
      }
   }
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to calculate the eigenvalues of Casimir's operators for the group G2, after Racah IV
// --------------------------------------------------------------------------------------------------------------- //
double racah_g(qG2 U, bool R5flag)
{
   if(R5flag)
      return (U.u1*(U.u1+3.) + U.u2*(U.u2+1.))/6.;                        // From Racah IV, eqn 19, or Judd, 5-50
   else
      return (U.u1*U.u1 + U.u1*U.u2 + U.u2*U.u2 + 5.*U.u1 + 4.*U.u2)/12.; // Racah IV, eqn. 28
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to calculate the eigenvalues of Casimir's operators for the group R7, after Judd, "Operator Techniques"
// --------------------------------------------------------------------------------------------------------------- //
double racah_g(qR7 W)
{
   return (W.w1*(W.w1+5.) + W.w2*(W.w2+3.) + W.w3*(W.w3+1.))/10.;         // From Racah IV, eqn 19, or Judd, 5-50
}


// --------------------------------------------------------------------------------------------------------------- //
// Calculates the matrix elements of the coulomb operator e3, after Racah IV.
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_e3(int n)
{
   int i,j,S2,v;
   orbital L;
   qG2 U;
   qR7 W;
   fconf conf(n);
   int num_states = (int)conf.states.size();
   sMat<double> e3(num_states,num_states);
   double e;
   
   for(j=0; j<num_states; j++)
   {
      S2 = conf.states[j].S2; v = conf.states[j].v; L = conf.states[j].L; W = conf.states[j].W; U = conf.states[j].U;
      for(i=0; i<=j; i++)
      {
	 // e = 0;
         if(i==j && S2==n)                                              // States of maximum spin, diagonal elements
	    e3(i,i) = -3. * (abs(L)*(abs(L)+1.)/2 - 12.*racah_g(U));    //  Eqn 81, Racah IV.
	 //{e = -3. * (abs(L)*(abs(L)+1.)/2 - 12.*racah_g(U));    //  Eqn 81, Racah IV.
	 // e3(i,i) = e; }
	 else if(abs(L)==abs(conf.states[i].L) && S2==conf.states[i].S2)
	 {
	    if(v==conf.states[i].v)
	    {
	       if((n==6 && v==6) || (n==7 && v==7)) { e = 0.; }
	       else if(n==(v+2) || n==(v+4))
	       {
	          if(v==S2) e = -2. * (abs(L)*(abs(L)+1.)/2 - 12.*racah_g(U));           // Eqn 81
	          else      e = racah_yfn(v,v,S2,conf.states[i].U,v,U)
		               * racah_phi(conf.states[i].U,U,conf.states[i].L,L);
                  e *= ( n==(v+2) ? ((1.-v)/(7.-v)) : (-4./(7.-v)) );                    // Eqn 82 and 84.
	       }
	       else                                                                      // Eqn 87 ---\/
	          e = racah_yfn(n,v,S2,conf.states[i].U,v,U) * racah_phi(conf.states[i].U,U,conf.states[i].L,L);
            }
	    else if(n==5 && conf.states[i].v==1 && v==3 && conf.states[i].S2==1)         // Eqn 87
	       e = sqrt(2./5) * racah_yfn(3,1,1,conf.states[i].U,3,U) * racah_phi(conf.states[i].U,U,conf.states[i].L,L);
	    else if(n==6 && conf.states[i].v==0 && v==4 && conf.states[i].S2==0)         // Eqn 86a
	       e = sqrt(9./5) * racah_yfn(4,0,0,conf.states[i].U,4,U) * racah_phi(conf.states[i].U,U,conf.states[i].L,L);
	    else if(n==6 && conf.states[i].v==2 && v==4)                                 // Eqn 86c
	       e = sqrt(1./6) * racah_yfn(4,2,S2,conf.states[i].U,4,U) * racah_phi(conf.states[i].U,U,conf.states[i].L,L);
	    else if(n==7 && conf.states[i].v==1 && v==5 && conf.states[i].S2==1)         // Eqn 86d
	       e = sqrt(3./2) * racah_yfn(5,1,1,conf.states[i].U,5,U) * racah_phi(conf.states[i].U,U,conf.states[i].L,L);
	    else                                                                         // Eqn 87 
	       e = racah_yfn(n,conf.states[i].v,S2,conf.states[i].U,v,U) * racah_phi(conf.states[i].U,U,conf.states[i].L,L);
	    if(i!=j && e!=0) { e3(i,j) = e; e3(j,i) = e; }
	    if(i==j)         { e3(i,i) = e - (abs(L)*(abs(L)+1.)/2 - 12.*racah_g(U)); }  // Eqn 87 and 78
	 }
      }
   }
   return e3;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the coulomb interaction matrix elements, for f-electrons using Racah operators, after Racah IV.
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_emat(int n, double E0, double E1, double E2, double E3)
{
   int nn = n;
   sMat<double> emat;
   if (n>7) nn=14-n; if(nn<1) { std::cerr << "racah_emat: number of f-electrons n > 14 or < 1\n"; return emat; }
   if (nn==1) { emat(0,0) = E0+E1+E2+E3; return emat; }
   sMat<double> e2 = racah_e2(nn);
   sMat<double> e3 = racah_e3(nn);
   emat.zero(e2.nr(),e2.nc());
   fconf conf(nn);
   int i,v;
   double S;

   for(i=0; i<e2.nr(); i++)
   {
      v = conf.states[i].v; S = conf.states[i].S2/2.;
      emat(i,i) = (n*(n-1.)/2)*E0 + ( (9./2)*(n-v) + v*(v+2.)/4 - S*(S+1.) )*E1;
   }
   e2 *= E2;
   e3 *= E3;
   return emat+e2+e3;
}
// --------------------------------------------------------------------------------------------------------------- //
// Looks up the coulomb interaction matrix elements for d-electrons from Tables from Neilson and Koster, 1963.
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_emat(int n, double F0, double F2, double F4)
{
   int nn = n;
   sMat<double> e;
   if (n>5) nn=10-n; if(nn<1) { std::cerr << "racah_emat: number of d-electrons n > 10 or < 1\n"; return e; }

   switch(nn) {
      case 1: e(0,0) = 0; break; // F0+F2+F4; break;
      case 2:
         e.zero(5,5); 
         e(0,0) = F0+F2/7.-4*F4/21.;          e(1,1) = F0-8*F2/49.-F4/49.;          e(2,2) = F0+2*F2/7.+2*F4/7.;
         e(3,3) = F0-3*F2/49.+4*F4/49.;       e(4,4) = F0+4*F2/49.+F4/441.;
         break;
      case 3:
         e.zero(8,8);
         e(0,0) = 3*F0-F4/3.;                 e(1,1) = 3*F0-(15*F2+8*F4)/49.;       e(2,2) = 3*F0-6*F2/49.-4*F4/147.; 
         e(3,3) = 3*F0+(F2+F4)/7.;            e(4,4) = 3*F0+3*F2/49.-19*F4/147.;    e(5,5) = 3*F0+9*F2/49.-29*F4/147.;
         e(6,6) = 3*F0-11*F2/49.+13*F4/441.;  e(7,7) = 3*F0-6*F2/49.-4*F4/147.;
         e(3,4) = (3*F2/49.-5*F4/147.)*sqrt(21.);  e(4,3) = e(3,4); 
         break;
      case 4:
         e.zero(16,16);
         e(0,0) = 6*F0-3*(F2+F4)/7.;          e(1,1) = 6*F0-F2/7.-2*F4/63.;         e(2,2) = 6*F0-3*F2/49.-139*F4/441.;
         e(3,3) = 6*F0-5*F2/49.-43*F4/147.;   e(4,4) = 6*F0-2*F2/49.-13*F4/147.;    e(5,5) = 6*F0-8*F2/49.-38*F4/147.;
         e(6,6) = 6*F0-12*F2/49.-94*F4/441.;  e(7,7) = 6*F0-17*F2/49.-23*F4/147.;   e(8,8) = 6*F0+2*(F2+F4)/7.;
         e(9,9) = 6*F0+6*F2/49.-38*F4/147.;   e(10,10)=6*F0+(15*F2-6*F4)/49.;       e(11,11)=6*F0+(3*F2-11*F4)/49.;
         e(12,12)=6*F0-4*F4/21.;              e(13,13)=6*F0-6*F2/49.+17*F4/147.;    e(14,14)=6*F0-4*F2/49.-64*F4/441.;
         e(15,15)=6*F0-(15*F2+F4)/49.;
         e(1,2) = (4*F2/49.-20*F4/441.)*sqrt(14.);  e(2,1) = e(1,2);
         e(4,5) = 12*F2/49.-20*F4/147.;             e(5,4) = e(4,5);
         e(8,9) = (6*F2/49.-10*F4/147.)*sqrt(21.);  e(9,8) = e(8,9);
         e(10,11)=(12*F2/49.-20*F4/147.)*sqrt(2.);  e(11,10)=e(10,11);
         e(13,14)=(4*F2/49.-20*F4/441.)*sqrt(11.);  e(14,13)=e(13,14);
         break;
      case 5:
         e.zero(16,16);
         e(0,0) = 10*F0-5*(F2+F4)/7.;         e(1,1) = 10*F0-4*F2/7.-5*F4/21.;      e(2,2) = 10*F0-(18*F2+25*F4)/49.;
         e(3,3) = 10*F0-(13*F2+20*F4)/49.;    e(4,4) = 10*F0-25*F2/49.-190*F4/441.; e(5,5) = 10*F0-3*F2/49.-65*F4/147.;
         e(6,6) = 10*F0+20*F2/49.-80*F4/147.; e(7,7) = 10*F0;                       e(8,8) = 10*F0-4*F2/49.-40*F4/147.;
         e(9,9) = 10*F0-(6*F2+20*F4)/49.;     e(10,10)=10*F0-25*F2/49.-5*F4/147.;   e(11,11)=10*F0-9*F2/49.-55*F4/147.;
         e(12,12)=10*F0+3*F2/49.-155*F4/441.; e(13,13)=10*F0-13*F2/49.-145*F4/441.; e(14,14)=10*F0-22*F2/49.-10*F4/147.;
         e(15,15)=10*F0-(24*F2+10*F4)/49.;
         e(7,9) = (-6*F2/49.+10*F4/147.)*sqrt(14.); e(9,7)=e(7,9);
         break;
      default: return e;
   }
   return e;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the coulomb interaction matrix elements for p-electrons from a formula of van Vleck PR 45, 405 (1936)
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_emat(int n, double F0, double F2)
{
   int nn = n;
   sMat<double> e;
   if (n>3) nn=6-n;  if(nn<1) { std::cerr << "racah_emat: number of p-electrons n > 6 or < 1\n"; return e; }

   int Hsz[]={1,3,3}, L[]={1,1,0,2,0,1,2},                 // p1:2P - p2: 3P,1S,1D - p3: 4S,2P,2D
                     S2[]={1,2,0,0,3,1,1}, st=0;

   F2/=25; nn--;  // Converts from Slater nomalisation to Condon and Shortley normalisation

   for (int ii=0; ii<nn; ii++) st+=Hsz[ii];
   // Uses formula 38 of van Vleck
   for (int ii=0; ii<Hsz[nn]; ii++) e(ii,ii) = (n*(n+1.)/2.)*F0 + ((-5*n*n+20*n-3*L[ii+st]*(L[ii+st]+1)-12*(S2[ii+st]/2.)*((S2[ii+st]/2.)+1))/2.)*F2;

   return e;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the conversions between Coulomb interaction parameters for Slater integrals and Racah's e_k operators
//    NB These conversion factors are only valid for f-electrons!!! 
//    For d-electrons to convert from Condon/Shortley F_k to Slater F^k: F^2=F_2*49 and F^4=F_4*441
//    For p-electrons to convert from Condon/Shortley F_k to Slater F^k: F^2=F_2*25
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> racah_FtoE(std::vector<double> F)      // Converts from F_k to E
{
   std::vector<double> E(4,0);
   E[0] = ( F[0] - 10*F[1] -  33*F[2] -  286*F[3] );
   E[1] = (        70*F[1] + 231*F[2] + 2002*F[3] )/9.;
   E[2] = (           F[1] -   3*F[2] +    7*F[3] )/9.;
   E[3] = (         5*F[1] +   6*F[2] -   91*F[3] )/3.;
   return E;
}
std::vector<double> racah_EtoF(std::vector<double> E)      // Converts from E to F_k
{
   std::vector<double> F(4,0);
   F[0] = ( E[0] + 9*E[1]/7. );
   F[1] = (          E[1] + 143*E[2] + 11*E[3] )/42.;
   F[2] = (          E[1] - 130*E[2] +  4*E[3] )/77.;
   F[3] = (          E[1] +  35*E[2] -  7*E[3] )/462.;
   return E;
}
std::vector<double> racah_FtoF_k(std::vector<double> F)    // Converts from F^k to F_k
{
   std::vector<double> F_k(4,0);
   F_k[0] = F[0];
   F_k[1] = F[1]/225.;
   F_k[2] = F[2]/1089.;
   F_k[3] = F[3]/(184041./25);
   return F_k;
}
std::vector<double> racah_F_ktoF(std::vector<double> F_k)  // Converts from F_k to F^k
{
   std::vector<double> F(4,0);
   F[0] = F_k[0]; 
   F[1] = F_k[1]*225.; 
   F[2] = F_k[2]*1089.; 
   F[3] = F_k[3]*(184041./25);
   return F;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the linear configuration interaction Hamilton matrix, after Rajnak and Wybourne, PR, v132, 280 (1963)
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_ci(int n, double alpha, double beta, double gamma)    // For f-electrons
{
   int nn = n;
   sMat<double> ci;
   if (n>7) nn=14-n; if(nn<1) { std::cerr << "racah_ci: number of f-electrons n > 14 or < 1\n"; return ci; }
   if (nn==1) { ci(0,0) = alpha*3.*(3.+1.) + beta/2. + gamma*6./10.; return ci; } // f1 has only 1 state: 2F [100][10]
   fconf conf(nn);
   ci.zero(conf.states.size(),conf.states.size());
   int i,L;

   for(i=0; i<ci.nr(); i++)
   {
      L = abs(conf.states[i].L);                                         // See Eqn 19 and similar...
      ci(i,i) = alpha*L*(L+1.) + beta*racah_g(conf.states[i].U) + gamma*racah_g(conf.states[i].W);
   }

   return ci;
}

sMat<double> racah_ci(int n, double alpha, double beta)                  // For d-electrons
{
   int nn = n;
   sMat<double> ci;
   if (n>5) nn=10-n; if(nn<1) { std::cerr << "racah_ci number of d-electrons n > 10 or < 1\n"; return ci; }
   if (nn==1) { ci(0,0) = alpha*2.*(2.+1.) + beta*2./3.; return ci; }    // d1 has only 1 state: 2D [10]
   fconf conf(nn,D);
   ci.zero(conf.states.size(),conf.states.size());
   int i,L;

   for(i=0; i<ci.nr(); i++)
   {
      L = abs(conf.states[i].L);
      ci(i,i) = alpha*L*(L+1.) + beta*racah_g(conf.states[i].U,true);    // Eqn 22 of Rajnak and Wybourne
   }

   return ci;
}

sMat<double> racah_ci(int n, double alpha)                               // For p-electrons
{
   int nn = n;
   sMat<double> ci;
   if (n>3) nn=6-n; if(nn<1) { std::cerr << "racah_ci number of p-electrons n > 6 or < 1\n"; return ci; }
   nn--;
   int Hsz[]={1,3,3}, L[]={1,1,0,2,0,1,2}, st=0;                         // p1:2P - p2: 3P,1S,1D - p3: 4S,2P,2D
   for (int i=1; i<nn; i++) st+=Hsz[i-1];
   for (int i=0; i<Hsz[nn]; i++) ci(i,i)=alpha*L[i+st]*(L[i+st]+1.);     // Eqn 21 of Rajnak and Wybourne

   return ci;
}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile: g++ coulomb.cpp states.cpp; ./a.out 2
// --------------------------------------------------------------------------------------------------------------- //
/*int main(int argc, char *argv[])
{
   int n;
   if(argc>1) { n = atoi(argv[1]); } else { n = 2; }
   //qR7 W(2,1,1); qG2 U(2,1); qG2 Up(2,1);

   //std::vector<numdenom> xwu = racah_xwu(W,U,Up);
   //std::cout << xwu[0] << ", " << xwu[1] << "\n";
   //Up.set(1,1); std::vector<numdenom> xwu2 = racah_xwu(W,U,Up);
   //xwu = xwu2; std::cout << xwu[0] << ", " << xwu[1] << "\n";
   
   //qG2 U(0,0); std::cout << U.isequal("00") << "\n";
   //qG2 Up(2,2);
   //std::cout << racah_phi(U,Up,S,S) << "\n";
   //std::vector< std::vector<double> > e2 = racah_e2(2);
   sMat<double> e2 = racah_e2(n);
   std::cout << "\ne2 =\n"; std::cout << e2;
   //std::cout << e2.display_full();
   sMat<double> e3 = racah_e3(n);
   std::cout << "\ne3 =\n"; std::cout << e3;
   sMat<double> emat = racah_emat(n,1.,1.,1.,1.);
   std::cout << emat;

   return 0;
}*/
