/* icpars.cpp
 *
 * Holds classes which manipulate the CF and free ion parameters, and calculates the moment.
 *
 * Classes (in icpars.hpp)
 *   cfpars  - Contains the CF parameters. Converts between different normalisations and energy units.
 *   icpars  - Contains the free ion parameters. Loads defaults for known ions. Converts between energy units.
 *   iceig   - Contains and calculates eigenvectors/values by various methods (full,partial,arnoldi)
 *   icmfmat - Contains the L and S matrices and calculates the magnetic moment
 *
 * Functions:
 *   void strtolower(std::string &instring);                           // Converts a string to lower case
 *   void conv_e_units(icpars &flags, std::string &newunits);          // Converts 1-ion pars to diff. energy units
 *   std::vector<double> stev_thetak(int n, orbital l);                // Calculates the Stevens factors.
 *   std::vector<double> rk_int(std::string &ionname);                 // Looks up value of radial integrals
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"
#include <cctype>                  // For std::tolower
#include <fstream>

#define SMALL 1e-6   // must match SMALL in mcdisp.c and ionpars.cpp because it is used to decide wether for small
		     // transition, energy the matrix Mijkl contains wn-wn' or wn/kT
#define MAXNOFCHARINLINE 144
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the number of allowed states from the number of electrons and l
// --------------------------------------------------------------------------------------------------------------- //
int getdim(int n, orbital l)                                    // Number of states = ^{4l+2}C_{n}
{
   //      4l+2       N        N!          (N-k+1).(N-k+2)...N
   // ns =     C   ==  C  = ---------  =  ---------------------   with N=4l+2, k=n
   //           n       k    k!(N-k)!         1.2.3...k
   //
   int i,j=1,nn=(n>(2*abs(l)+1))?(4*abs(l)+2-n):n;              // Use n<2l+1 equivalents only to avoid overflow
   long int ns=1;
   for(i=(4*abs(l)+2-nn+1); i<=(4*abs(l)+2); i++) ns*=i;        // Computes (4l+2-n+1).(4l+2-n+1)...(4l+2)
   for(i=nn; i>1; i--) j*=i; ns/=j;                             // Computes n.(n-1)...1
   return ns;
/* // Wikipedia Algorithm
   int k=4*abs(l)+2;
   if (n > k)   return 0;
   if (n > k/2) n = k-n;                                        // Take advantage of symmetry
   long double accum = 1;
   for (int i = 1; i <= n; i++)
      accum *= ( (k-n+i) / i );
   return (int) (accum + 0.5);                                  // avoid rounding error */
}

// --------------------------------------------------------------------------------------------------------------- //
// Looks up conversion factor for Wybourne to Stevens...
// --------------------------------------------------------------------------------------------------------------- //
//double wy2stev(int i)
//{
//   double lll[] = {0,0,0,0,0,0,sqrt(6.)/2., -sqrt(6.), 1./2., -sqrt(6.), sqrt(6.)/2., 0,0,0,0,0,0,0, /*k=4*/ sqrt(70.)/8., -sqrt(35.)/2., 
//                 sqrt(10.)/4., -sqrt(5.)/2., 1./8., -sqrt(5.)/2., sqrt(10.)/4., -sqrt(35.)/2., sqrt(70.)/8., 0,0,0,0,0,0,0,0,0,0,0,
//         /*k=6*/ sqrt(231.)/16., -3*sqrt(77.)/8., 3*sqrt(14.)/16., -sqrt(105.)/8., sqrt(105.)/16., -sqrt(42.)/8., 
//                 1./16., -sqrt(42.)/8., sqrt(105.)/16., -sqrt(105.)/8., 3*sqrt(14.)/16., -3*sqrt(77.)/8., sqrt(231.)/16.};
//   return lll[i];
//}
//
//double clm(int i)
//{
//   double clm[]={sqrt(3./4./PI),sqrt(3./4./PI),sqrt(3./4./PI),sqrt(3./4./PI),sqrt(3./4./PI),sqrt(3./4./PI),
//                 sqrt(15./PI)/4.,sqrt(15./PI)/2.,sqrt(5./PI)/4.,sqrt(15./PI)/2.,sqrt(15./PI)/4.,
//                 sqrt(35./32./PI),sqrt(105./16./PI),sqrt(21./32./PI),sqrt(7./16./PI),sqrt(21./32./PI),sqrt(105./16./PI),sqrt(35./32./PI),
//                 (3./16)*sqrt(35./PI),(3./8)*sqrt(70./PI),(3./8)*sqrt(5./PI),(3./4)*sqrt(5./2./PI),(3./16)*sqrt(1./PI),(3./4)*sqrt(5./2./PI),(3./8)*sqrt(5./PI),
//                    (3./8)*sqrt(70./PI),(3./16)*sqrt(35./PI),
//                 sqrt(693./512./PI), sqrt(3465./256./PI), sqrt(385./512./PI), sqrt(1155./64./PI), sqrt(165./256./PI), sqrt(11./256./PI),
//                 sqrt(165./256./PI), sqrt(1155./64./PI), sqrt(385./512./PI), sqrt(3465./256./PI), sqrt(693./512./PI),
//                 (231./64.)*sqrt(26./231./PI), sqrt(9009./512./PI), (21./32.)*sqrt(13./7./PI), (1./32.)*sqrt(2730./PI), (1./64.)*sqrt(2730./PI), (1./8.)*sqrt(273./4./PI),
//                 (1./32)*sqrt(13./PI),
//                 (1./8.)*sqrt(273./4./PI), (1./64.)*sqrt(2730./PI), (1./32.)*sqrt(2730./PI), (21./32.)*sqrt(13./7./PI), sqrt(9009./512./PI), (231./64.)*sqrt(26./231./PI)};
//   return clm[i];
//}                

// --------------------------------------------------------------------------------------------------------------- //
// Converts a C++ string to lower case
// --------------------------------------------------------------------------------------------------------------- //
void strtolower(std::string &instring)
{
   int i,strlength = instring.length();
   for(i=0; i<strlength; i++)
      instring[i] = std::tolower(instring[i]);
}
// --------------------------------------------------------------------------------------------------------------- //
void str2upper(std::string &instring)
{
   int i,strlength = instring.length();
   for(i=0; i<strlength; i++)
      instring[i] = std::toupper(instring[i]);
}

// --------------------------------------------------------------------------------------------------------------- //
// Converts the single ion parameters to a different energy unit. Supported units: cm^{-1}, meV, K
// --------------------------------------------------------------------------------------------------------------- //
void conv_e_units(icpars &pars, std::string &newunit)
{
   int k;
   pars.B.conv_e_units(newunit);
#define LOOP(X) for(k=0; k<X; k++)
   if(pars.e_units.find("cm")!=std::string::npos || pars.e_units.find("wave")!=std::string::npos) 
   {
      if(newunit.find("cm")!=std::string::npos || newunit.find("wave")!=std::string::npos) { pars._econv=1.; }
      else if(newunit.find("meV")!=std::string::npos)  // Convert cm^{-1} to meV
      {  
         LOOP(4) pars.F[k]/=MEV2CM; pars.xi/=MEV2CM; pars.e_units="meV"; LOOP(3) pars.alpha[k]/=MEV2CM;     pars._econv=MEV2CM;
      }
      else if(newunit.find("K")!=std::string::npos)    // Convert cm^{-1} to K
      {  
         LOOP(4) pars.F[k]*=CM2K;   pars.xi*=CM2K;   pars.e_units="Kelvin"; LOOP(3) pars.alpha[k]*=CM2K;    pars._econv=1/CM2K;
      }
      else std::cerr << "conv_e_units(): Energy units " << newunit << " not recognised. Accepted units are cm^{-1}, meV, K.\n";
   }
   else if(pars.e_units.find("meV")!=std::string::npos) 
   { 
      if(newunit.find("cm")!=std::string::npos || newunit.find("wave")!=std::string::npos) // Convert meV to cm^{-1}
      {  
         LOOP(4) pars.F[k]*=MEV2CM; pars.xi*=MEV2CM; pars.e_units="cm^{-1}"; LOOP(3) pars.alpha[k]*=MEV2CM; pars._econv=1.;
      }
      else if(newunit.find("meV")!=std::string::npos) { pars._econv=MEV2CM; }
      else if(newunit.find("K")!=std::string::npos)    // Convert meV to K
      {  
         LOOP(4) pars.F[k]*=MEV2K;  pars.xi*=MEV2K;  pars.e_units="Kelvin"; LOOP(3) pars.alpha[k]*=MEV2K;   pars._econv=1/CM2K;
      }
      else std::cerr << "conv_e_units(): Energy units " << newunit << " not recognised. Accepted units are cm^{-1}, meV, K.\n";
   }
   else if(pars.e_units.find("K")!=std::string::npos) 
   { 
      if(newunit.find("cm")!=std::string::npos || newunit.find("wave")!=std::string::npos) // Convert K to cm^{-1}
      {  
         LOOP(4) pars.F[k]/=CM2K;   pars.xi/=CM2K;   pars.e_units="cm^{-1}"; LOOP(3) pars.alpha[k]/=CM2K;   pars._econv=1.;
      }
      else if(newunit.find("meV")!=std::string::npos)  // Convert K to meV
      {  
         LOOP(4) pars.F[k]/=MEV2K;  pars.xi/=MEV2K;  pars.e_units="meV"; LOOP(3) pars.alpha[k]/=MEV2K;      pars._econv=MEV2CM;
      }
      else if(newunit.find("K")!=std::string::npos) { pars._econv=1/CM2K; }
      else std::cerr << "conv_e_units(): Energy units " << newunit << " not recognised. Accepted units are cm^{-1}, meV, K.\n";
   }
   else 
      std::cerr << "conv_e_units(): Energy units " << pars.e_units << " not recognised. Accepted units: cm^{-1}, meV, K.\n";
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the Stevens operator equivalent factors after Elliot, Judd and Runciman.
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> stev_thetak(int n, orbital l)
{
   std::vector<double> theta_k(3); 
   if(n==1) { theta_k[0] = -2/35.; theta_k[1] = 2/315.; theta_k[2] = 0.; return theta_k; }
   theta_k[0] = -8.*sqrt(7./15.); theta_k[1] = 16.*sqrt(14./11.); theta_k[2] = -640.*sqrt(7./429.);

   orbital L;
   int S2,J2;
   int ml, Li = 0, lm = abs(l), ln = (n<=(2*l+1)) ? (abs(l)-n) : (abs(l)-n+(2*l+1)); 
   fconf conf(n,l);
   fconf confp(n-1,l);
   std::vector<cfpls> cfps;
   double sumcfp,noncfpprod;
   int k,is,ic;
   bool df = false; if(l==D) df = true;

   // Calculates L,S,J from Hund's rules
   S2 = (n<=(2*l+1)) ? n : ((4*l+2)-n);                       // Maximise S.
   for(ml=lm; ml>ln; ml--) Li+=ml;                            // Maximise L.
   J2 = (n<=(2*l+1)) ? abs(2*Li-S2) : (2*Li+S2);              // J=|L-S| for <half-full, J=L+S for >half-full shell.
   L = (orbital)Li;

   // Finds the seniority and other quantum numbers from the list of states
   for(is=0; is<(int)conf.states.size(); is++) if(conf.states[is].S2==S2 && conf.states[is].L==L) break;

   noncfpprod = pow(-1.,-(double)l-Li) * (2.*Li+1.);
   if(df) cfps = racah_parents(n,conf.states[is].v,S2,L); else cfps = racah_parents(n,conf.states[is].v,conf.states[is].U,S2,L);

   // We use the formulae of Elliot, Judd and Runciman, Proc. R. Soc. Lon. A, v240, pp509, 1957 to calculate theta_k
   for(k=0; k<3; k++)
   {
      sumcfp = 0.;
      for(ic=0; ic<(int)cfps.size(); ic++) 
         sumcfp += racahW(l*2,Li*2,l*2,Li*2,abs(confp.states[cfps[ic].ind].L)*2,(k+1)*4) * cfps[ic].cfp*cfps[ic].cfp 
	           * pow(-1.,(double)abs(confp.states[cfps[ic].ind].L)+2.*(k+1.)) * noncfpprod;
      theta_k[k] *= sqrt(factorial(J2-2*(k+1))/factorial(J2+2*(k+1)+1)) * n * sumcfp * pow(-1.,2.*(k+1.)+Li+(S2+J2)/2.) 
                    * (J2+1) * sixj(Li*2,J2,S2,J2,Li*2,4*(k+1));
   }
   return theta_k;
}
// --------------------------------------------------------------------------------------------------------------- //
// Looks up value of radial integrals
// --------------------------------------------------------------------------------------------------------------- //
#define IONCMP ionname.compare
std::vector<double> rk_int(std::string &ionname)
{
   std::vector<double> rk(3);
   ionname.assign(ionname); strtolower(ionname);
   // Values taken from program cfield, by Peter Fabi, FZ Juelich, file theta.c
        if(IONCMP("ce3+")==0) { rk[0] = 1.309;  rk[1] = 3.964;  rk[2] = 23.31;  }  
   else if(IONCMP("pr3+")==0) { rk[0] = 1.1963; rk[1] = 3.3335; rk[2] = 18.353; }  /* U. Walter Diss.         */
   else if(IONCMP("nd3+")==0) { rk[0] = 1.114;  rk[1] = 2.910;  rk[2] = 15.03;  }  
   else if(IONCMP("pm3+")==0) { rk[0] = 1.0353; rk[1] = 2.5390; rk[2] = 12.546; }  /*          -"-            */
   else if(IONCMP("sm3+")==0) { rk[0] = 0.9743; rk[1] = 2.260;  rk[2] = 10.55;  }  
   else if(IONCMP("eu3+")==0) { rk[0] = 0.9175; rk[1] = 2.020;  rk[2] = 9.039;  }  
   else if(IONCMP("gd3+")==0) { rk[0] = 0.8671; rk[1] = 1.820;  rk[2] = 7.831;  }  
   else if(IONCMP("tb3+")==0) { rk[0] = 0.8220; rk[1] = 1.651;  rk[2] = 6.852;  }  
   else if(IONCMP("dy3+")==0) { rk[0] = 0.7814; rk[1] = 1.505;  rk[2] = 6.048;  }  
   else if(IONCMP("ho3+")==0) { rk[0] = 0.7446; rk[1] = 1.379;  rk[2] = 5.379;  }  
   else if(IONCMP("er3+")==0) { rk[0] = 0.7111; rk[1] = 1.270;  rk[2] = 4.816;  }  
   else if(IONCMP("tm3+")==0) { rk[0] = 0.6804; rk[1] = 1.174;  rk[2] = 4.340;  }  
   else if(IONCMP("yb3+")==0) { rk[0] = 0.6522; rk[1] = 1.089;  rk[2] = 3.932;  }  
   else if(IONCMP("u4+")==0)  { rk[0] = 2.042;  rk[1] = 7.632;  rk[2] = 47.774; }  /* Freeman et al. PRB 13 (1976) 1168 */
   else if(IONCMP("u3+")==0)  { rk[0] = 2.346;  rk[1] = 10.906; rk[2] = 90.544; }  /* Freeman et al. PRB 13 (1976) 1168 */
   else if(IONCMP("np4+")==0) { rk[0] = 1.884;  rk[1] = 6.504;  rk[2] = 37.80;  }  /* Lewis et al. J. Chem Phys. 53 (1970) 809*/
   else if(IONCMP("nd2+")==0) { rk[0] = 1.392;  rk[1] = 5.344;  rk[2] = 45.450; }  
   else if(IONCMP("sm2+")==0) { rk[0] = 1.197;  rk[1] = 3.861;  rk[2] = 28.560; }  
   else if(IONCMP("eu2+")==0) { rk[0] = 1.098;  rk[1] = 3.368;  rk[2] = 23.580; }  
   else if(IONCMP("gd2+")==0) { rk[0] = 1.028;  rk[1] = 2.975;  rk[2] = 19.850; }  
   else if(IONCMP("tb2+")==0) { rk[0] = 0.968;  rk[1] = 2.655;  rk[2] = 16.980; }  
   else if(IONCMP("dy2+")==0) { rk[0] = 0.913;  rk[1] = 2.391;  rk[2] = 14.730; }  
   else if(IONCMP("ho2+")==0) { rk[0] = 0.866;  rk[1] = 2.169;  rk[2] = 12.920; }  
   else if(IONCMP("er2+")==0) { rk[0] = 0.824;  rk[1] = 1.979;  rk[2] = 11.450; }  
   else if(IONCMP("tm2+")==0) { rk[0] = 0.785;  rk[1] = 1.819;  rk[2] = 10.240; }  
   else { rk[0] = 1.; rk[1] = 1.; rk[2] = 1.; }
   
   return rk;
}

// --------------------------------------------------------------------------------------------------------------- //
// Constructor function for the icpars class 
// --------------------------------------------------------------------------------------------------------------- //
icpars::icpars()
{
   alpha.assign(3,0.); F.assign(4,0.); xi = 0.; _econv = 1.; _alpha.assign(3,0.); _F.assign(4,0.); 
   #ifdef JIJCONV
   jijconv.assign(52,1.); _jijconvalreadycalc = false;
   #endif
   n = 1; l = (orbital)3; e_units.assign("cm^{-1}"); calcphys = 0; mag_units = 0;
   xT=0.; xHa=0.; xHb=0.; xHc=0.; xMin=0.; xStep=0.; xMax=0.;
   yT=0.; yHa=0.; yHb=0.; yHc=0.; yMin=0.; yStep=0.; yMax=0.;
   Bx=0.; By=0.;  Bz=0.; basis.assign("JmJ"); save_matrices = false;
   perturb = false; partial = false; arnoldi = false; spectrelevels = -1; truncate_level = 1; num_eigv = 4;
   partial_standalone = false; arnoldi_standalone = false;
}
// --------------------------------------------------------------------------------------------------------------- //
// Methods functions for the class icpars
// --------------------------------------------------------------------------------------------------------------- //
#ifdef JIJCONV
void icpars::jijconvcalc()
{
   if(_jijconvalreadycalc) return;
   int k[] = {0,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   double pna[]={0,0,0,0,0,0,0,sqrt(15/PI)/4,sqrt(15/PI)/2,sqrt(5/PI)/4,sqrt(15/PI)/2,sqrt(15/PI)/4,
      sqrt(35./32/PI),sqrt(105./16/PI),sqrt(21./32/PI),sqrt(7./16/PI),sqrt(21./32/PI),sqrt(105./16/PI),sqrt(35./32/PI),
      3*sqrt(35/PI)/16,3*sqrt(70/PI)/8,3*sqrt(5/PI)/8,3*sqrt(5./2/PI)/4,3./16/sqrt(PI),3*sqrt(5./2/PI)/4,3*sqrt(5/PI)/8,3*sqrt(70/PI)/8,3*sqrt(35/PI)/16,
      sqrt(693./512/PI),sqrt(3465./256/PI),sqrt(385./512/PI),sqrt(1155./64/PI),sqrt(165./256/PI),sqrt(11./256/PI),sqrt(165./256/PI),sqrt(1155./64/PI),
         sqrt(385./512/PI),sqrt(3465./256/PI),sqrt(693./512/PI),
      231*sqrt(26./231/PI)/64,sqrt(9009./512/PI),21*sqrt(13./7/PI)/32,sqrt(2730/PI)/32,sqrt(2730/PI)/64,sqrt(273./4/PI)/8,sqrt(13/PI)/32,sqrt(273./4/PI)/8,
         sqrt(2730/PI)/64,sqrt(2730/PI)/32,21*sqrt(13./7/PI)/32,sqrt(9009./512/PI),231*sqrt(26./231/PI)/64};
   for(int iq=7; iq<=51; iq++)
   {
      if(q[iq]==0) jijconv[iq]=sqrt((2*k[iq]+1.)/4./PI)/pna[iq]; else jijconv[iq]=sqrt((2*k[iq]+1.)/8./PI)/pna[iq];
      if(k[iq]==2) jijconv[iq]/=B.alpha(); if(k[iq]==4) jijconv[iq]/=B.beta(); if(k[iq]==6) jijconv[iq]/=B.gamma();
   }
   _jijconvalreadycalc = true;
}
#endif
// --------------------------------------------------------------------------------------------------------------- //
// Overloaded operators for icpars:: class
// --------------------------------------------------------------------------------------------------------------- //
bool icpars::operator ==(icpars c) const
{
   if(c.n!=n) return false; if(c.l!=l) return false;
   if(fabs(c._xi-_xi)>1e-5) return false;
   if(c.B!=B) return false;
   int i; 
   for(i=0; i<4; i++) if(fabs(c._F[i]-_F[i])>1e-5) return false;
   for(i=0; i<3; i++) if(fabs(c._alpha[i]-_alpha[i])>1e-5) return false;
   return true;
}
bool icpars::operator !=(icpars c) const
{
   if(c.n!=n) return true; if(c.l!=l) return true;
   if(fabs(c._xi-_xi)>1e-5) return true;
   if(c.B!=B) return true;
   int i; 
   for(i=0; i<4; i++) if(fabs(c._F[i]-_F[i])>1e-5) return true;
   for(i=0; i<3; i++) if(fabs(c._alpha[i]-_alpha[i])>1e-5) return true;
   return false;
}

// --------------------------------------------------------------------------------------------------------------- //
// Constructor functions for the cfpars class 
// --------------------------------------------------------------------------------------------------------------- //
cfpars::cfpars()                                              // Blank constructor
{
   int i; for(i=0; i<27; i++) { _Bi[i] = 0.; _Bo[i] = 0.; }
   _minustype = false; _withiontype = false;
   _normalisation.assign("Wybourne"); _cfname.assign("L"); 
   _stevfact.assign(3,1.); _istevfact.assign(3,1.); _rk.assign(3,1.);
   _units.assign("cm^{-1}");
}
cfpars::cfpars(std::string &ionname, int n, orbital l)        // Constructor function with stevens factor and <r^k>
{
   int i;
   for(i=0; i<27; i++) { _Bi[i] = 0.; _Bo[i] = 0.; }
   _stevfact = stev_thetak(n,l); _rk = rk_int(ionname);
   for(i=0; i<3; i++) if(_stevfact[i]!=0) _istevfact[i] = 1/_stevfact[i]; else _istevfact[i] = 0.;
   _minustype = false; if(_rk[0]!=1.) _withiontype = true;
   _normalisation.assign("Stevens"); _cfname.assign("B");
   _units.assign("cm^{-1}");
}

// --------------------------------------------------------------------------------------------------------------- //
// Properties functions for the class cfpars
// --------------------------------------------------------------------------------------------------------------- //
void cfpars::assign(std::string &S, int &k, int &q, double v) // Assign a particular parameter of type S
{
   int i;
   double val = v;
   double l[] = {sqrt(6.)/2., -sqrt(6.), 1./2., -sqrt(6.), sqrt(6.)/2., /*k=4*/ sqrt(70.)/8., -sqrt(35.)/2., 
                 sqrt(10.)/4., -sqrt(5.)/2., 1./8., -sqrt(5.)/2., sqrt(10.)/4., -sqrt(35.)/2., sqrt(70.)/8.,
         /*k=6*/ sqrt(231.)/16., -3*sqrt(77.)/8., 3*sqrt(14.)/16., -sqrt(105.)/8., sqrt(105.)/16., -sqrt(42.)/8., 
                 1./16., -sqrt(42.)/8., sqrt(105.)/16., -sqrt(105.)/8., 3*sqrt(14.)/16., -3*sqrt(77.)/8., sqrt(231.)/16.};
   str2upper(S);

   if(k==2) i = 2+q; else if(k==4) i = 5+(4+q); else if(k==6) i = 14+(6+q); 
   else { std::cerr << "cfpars::assign() invalid rank k=" << k << " must be 2,4,6.\n"; return; }

#define MTP _minustype =
   if(_units.find("meV")!=std::string::npos) val *= MEV2CM;   // Converts parameters to 1/cm for internal use.
   else if(_units.find("K")!=std::string::npos) val /= CM2K;

   if(S.compare("A")==0)
   {
      if(_cfname.compare("A")!=0) conv(S); _Bo[i] = v; _Bi[i] = val*_rk[k/2-1]/l[i]; if(q<0) _Bi[i] = -_Bi[i]; MTP 0;
   } 
   else if(S.compare("W")==0)
   {
      if(_cfname.compare("W")!=0) conv(S); _Bo[i] = v; _Bi[i] = val*_rk[k/2-1]/l[i]; if(q!=0) _Bi[i] /= 2; MTP 1;
   }
   else if(S.compare("B")==0)
   {
      if(_cfname.compare("B")!=0) conv(S); _Bo[i] = v; _Bi[i] = val*_istevfact[k/2-1]/l[i]; if(q<0) _Bi[i] = -_Bi[i]; MTP 0;
   }
   else if(S.compare("V")==0)
   {
      if(_cfname.compare("V")!=0) conv(S); _Bo[i] = v; _Bi[i] = val*_istevfact[k/2-1]/l[i]; if(q!=0) _Bi[i] /= 2; MTP 1;
   }
   else if(S.compare("L")==0 || S.compare("D")==0)
   {
      if(_cfname.compare("L")!=0 && _cfname.compare("D")!=0) conv(S); _cfname.assign(S); _Bo[i] = v; _Bi[i] = val; MTP 1;
   }
   else if(S.compare("AR")==0)
   {
      if(_cfname.compare("AR")!=0) conv(S); _Bo[i] = v; _Bi[i] = val/l[i]; MTP 0;
   }
   else { std::cerr << "cfpars::assign() parameter type " << S << " not recognised. Must be either A,W,B,V,L,D,AR\n"; }
}
// --------------------------------------------------------------------------------------------------------------- //
void cfpars::calc_stevfact(int n, orbital l)
{
   _stevfact = stev_thetak(n,l);
   int i; for(i=0; i<3; i++) if(_stevfact[i]!=0) _istevfact[i] = 1/_stevfact[i]; else _istevfact[i] = 0.;
}
// --------------------------------------------------------------------------------------------------------------- //
void cfpars::find_rk(std::string &ionname)
{
   _rk = rk_int(ionname); if(_rk[0]!=1. && _rk[1]!=1. && _rk[2]!=1.) _withiontype = true;
}
// --------------------------------------------------------------------------------------------------------------- //
double cfpars::get(int k, int q)
{
   double val=0.;
   if(k==2) val = _Bo[2+q]; else if(k==4) val = _Bo[5+(4+q)]; else if(k==6) val = _Bo[14+(6+q)];
   return val;
}
// --------------------------------------------------------------------------------------------------------------- //

// --------------------------------------------------------------------------------------------------------------- //
// Methods functions for the class cfpars
// --------------------------------------------------------------------------------------------------------------- //
std::string cfpars::cfparsout(const char *delimiter)          // Prints the cf parameters as a delimited string
{
   int k,q,ik=0;
   std::stringstream retval;

   for(k=0; k<3; k++)
   {
      for(q=0; q<(4*(k+1)+1); q++)
         if(fabs(_Bo[ik+q])>DBL_EPSILON)
         {
            retval << _cfname << (k+1)*2;
            if(q<(2*(k+1)))
               if(_minustype) retval << "-" << (2*(k+1)-q) << "="; else retval << (2*(k+1)-q) << "S=";
            else retval << -(2*(k+1)-q) << "=";
            retval << _Bo[ik+q] << delimiter;
         }
      ik += q;
   }
   return retval.str();
}
// --------------------------------------------------------------------------------------------------------------- //
void cfpars::conv_B_norm(std::string &newnorm_)               // Converts cf parameters between Stevens and Wybourne
{
   int q;
   double lambda_2q[] = {sqrt(6.)/2., -sqrt(6.), 1./2., -sqrt(6.), sqrt(6.)/2.};
   double lambda_4q[] = {sqrt(70.)/8., -sqrt(35.)/2., sqrt(10.)/4., -sqrt(5.)/2., 1./8., -sqrt(5.)/2., sqrt(10.)/4., 
                         -sqrt(35.)/2., sqrt(70.)/8.};
   double lambda_6q[] = {sqrt(231.)/16., -3*sqrt(77.)/8., 3*sqrt(14.)/16., -sqrt(105.)/8., sqrt(105.)/16., -sqrt(42.)/8., 
                         1./16., -sqrt(42.)/8., sqrt(105.)/16., -sqrt(105.)/8., 3*sqrt(14.)/16., -3*sqrt(77.)/8., sqrt(231.)/16.};
   std::string oldnorm; oldnorm.assign(_normalisation); strtolower(oldnorm);
   std::string newnorm; newnorm.assign(newnorm_); strtolower(newnorm);
   if(oldnorm.find("stev")!=std::string::npos && newnorm.find("wy")!=std::string::npos)
   {
      for(q=0; q<5; q++)  _Bo[q]    /= lambda_2q[q];
      for(q=0; q<9; q++)  _Bo[5+q]  /= lambda_4q[q];
      for(q=0; q<13; q++) _Bo[14+q] /= lambda_6q[q];
      _normalisation.assign("Wybourne");
   }
   else if(oldnorm.find("wy")!=std::string::npos && newnorm.find("stev")!=std::string::npos)
   {
      for(q=0; q<5; q++)  _Bo[q]    *= lambda_2q[q];
      for(q=0; q<9; q++)  _Bo[5+q]  *= lambda_4q[q];
      for(q=0; q<13; q++) _Bo[14+q] *= lambda_6q[q];
      _normalisation.assign("Stevens");
   }
}
// --------------------------------------------------------------------------------------------------------------- //
void cfpars::conv_e_units(std::string &units)                 // Converts parameters to different energy units
{
   int k;
   if(units.find("cm")!=std::string::npos || units.find("wave")!=std::string::npos) 
   {
      if(_units.find("cm")!=std::string::npos || _units.find("wave")!=std::string::npos) {}
      else if(_units.find("meV")!=std::string::npos)          // Convert meV cm^{-1}
         for(k=0;k<27;k++) _Bo[k]*=MEV2CM; 
      else if(_units.find("K")!=std::string::npos)            // Convert K to cm^{-1}
         for(k=0;k<27;k++) _Bo[k]/=CM2K;
      else std::cerr << "conv_e_units(): Energy units " << _units << " not recognised. Accepted units are cm^{-1}, meV, K.\n";
      _units.assign("cm^{-1}");
   }
   else if(units.find("meV")!=std::string::npos) 
   { 
      if(_units.find("cm")!=std::string::npos || _units.find("wave")!=std::string::npos)
         for(k=0;k<27;k++) _Bo[k]/=MEV2CM;                    // Convert cm^{-1} to meV
      else if(_units.find("meV")!=std::string::npos) {}
      else if(_units.find("K")!=std::string::npos)            // Convert K to meV
         for(k=0;k<27;k++) _Bo[k]/=MEV2K;
      else std::cerr << "conv_e_units(): Energy units " << _units << " not recognised. Accepted units are cm^{-1}, meV, K.\n";
      _units.assign("meV");
   }
   else if(units.find("K")!=std::string::npos) 
   { 
      if(_units.find("cm")!=std::string::npos || _units.find("wave")!=std::string::npos)
         for(k=0;k<27;k++) _Bo[k]*=CM2K;                      // Convert cm^{-1} to K
      else if(_units.find("meV")!=std::string::npos)          // Convert meV to K
         for(k=0;k<27;k++) _Bo[k]*=MEV2K;
      else if(_units.find("K")!=std::string::npos) {}
      else std::cerr << "conv_e_units(): Energy units " << _units << " not recognised. Accepted units are cm^{-1}, meV, K.\n";
      _units.assign("K");
   }
   else 
      std::cerr << "conv_e_units(): Energy units " << units << " not recognised. Accepted units: cm^{-1}, meV, K.\n";
}
// --------------------------------------------------------------------------------------------------------------- //
void cfpars::conv(std::string &newcfname)                     // Converts parameters to a new type (A,W,B,V,L,D,AR)
{
   int k,q,i;
   double l[] = {sqrt(6.)/2., -sqrt(6.), 1./2., -sqrt(6.), sqrt(6.)/2., /*k=4*/ sqrt(70.)/8., -sqrt(35.)/2., 
                 sqrt(10.)/4., -sqrt(5.)/2., 1./8., -sqrt(5.)/2., sqrt(10.)/4., -sqrt(35.)/2., sqrt(70.)/8.,
         /*k=6*/ sqrt(231.)/16., -3*sqrt(77.)/8., 3*sqrt(14.)/16., -sqrt(105.)/8., sqrt(105.)/16., -sqrt(42.)/8., 
                 1./16., -sqrt(42.)/8., sqrt(105.)/16., -sqrt(105.)/8., 3*sqrt(14.)/16., -3*sqrt(77.)/8., sqrt(231.)/16.};
   double half[] = {.5,.5,1.,.5,.5, .5,.5,.5,.5,1.,.5,.5,.5,.5, .5,.5,.5,.5,.5,.5,1.,.5,.5,.5,.5,.5,.5};
   double imin[] = {-1.,-1.,1.,1.,1.,-1.,-1.,-1.,-1.,1.,1.,1.,1.,1.,-1.,-1.,-1.,-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.};
   strtoupper(newcfname);
#define CFCMP newcfname.compare
#define CFLOOP(ARG) i=0; for(k=0;k<3;k++) for(q=0; q<(4*(k+1)+1); q++) { ARG i++; } 
#define normstev _normalisation.assign("Stevens")
#define normwy _normalisation.assign("Wybourne")
        if(CFCMP("A")==0) { CFLOOP( _Bo[i]= _Bi[i]*l[i]/_rk[k]*imin[i]; _cfname.assign("A"); normstev; ) }
   else if(CFCMP("W")==0) { CFLOOP( _Bo[i]= _Bi[i]*l[i]/_rk[k]*half[i]; _cfname.assign("W"); normstev; ) }
   else if(CFCMP("B")==0) { CFLOOP( _Bo[i]= _Bi[i]*l[i]*_stevfact[k]*imin[i]; _cfname.assign("B"); normstev; ) }
   else if(CFCMP("V")==0) { CFLOOP( _Bo[i]= _Bi[i]*l[i]*_stevfact[k]*half[i]; _cfname.assign("V"); normstev; ) }
   else if(CFCMP("L")==0 || CFCMP("D")==0) { for(i=0;i<27;i++) _Bo[i] = _Bi[i]; _cfname.assign(newcfname); normwy; }
   else if(CFCMP("AR")==0) { CFLOOP( _Bo[i]= _Bi[i]*l[i]; _cfname.assign("AR"); normstev; ) }
   else std::cerr << "cfpars::conv(): CF parameter type " << newcfname << "kq not recognised. Accepted: A,W,B,V,L,D,AR\n";

   if(_units.find("K")!=std::string::npos)   for(i=0; i<27; i++) _Bo[i]*=CM2K;
   if(_units.find("meV")!=std::string::npos) for(i=0; i<27; i++) _Bo[i]/=MEV2CM;
}
// --------------------------------------------------------------------------------------------------------------- //

// --------------------------------------------------------------------------------------------------------------- //
// Overloaded operators for class cfpars
// --------------------------------------------------------------------------------------------------------------- //
double cfpars::operator ()(int k, int q) const                // Operator to access internal parameters Bi
{
   if(k==2) return _Bi[2+q]; else if(k==4) return _Bi[5+(4+q)]; else if(k==6) return _Bi[14+(6+q)]; else return 0.;

}
bool cfpars::operator ==(cfpars c) const
{
   int i;
   for(i=0; i<27; i++) if(fabs(c._Bi[i]-_Bi[i])>1e-5) return false;
   return true;
}
bool cfpars::operator !=(cfpars c) const
{
   int i;
   for(i=0; i<27; i++) if(fabs(c._Bi[i]-_Bi[i])>1e-5) return true;
   return false;
}
// --------------------------------------------------------------------------------------------------------------- //

