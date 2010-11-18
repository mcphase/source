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
   n = 1; l = (orbital)3; e_units.assign("cm^{-1}"); calcphys = 0; mag_units = 0;
   xT=0.; xHa=0.; xHb=0.; xHc=0.; xMin=0.; xStep=0.; xMax=0.;
   yT=0.; yHa=0.; yHb=0.; yHc=0.; yMin=0.; yStep=0.; yMax=0.;
   Bx=0.; By=0.;  Bz=0.; basis.assign("JmJ"); save_matrices = false;
   perturb = false; partial = false; arnoldi = false; spectrelevels = -1; truncate_level = 1; num_eigv = 4;
   partial_standalone = false; arnoldi_standalone = false;
}
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


// --------------------------------------------------------------------------------------------------------------- //
// Constructors for class iceig::
// --------------------------------------------------------------------------------------------------------------- //
iceig::iceig(sMat<double>&H)
{
   _Hsz = H.nc(); _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz]; _zV = 0;
   int info = ic_diag(H,_V,_E); 
   if(info!=0) { std::cerr << "iceig(H) - Error diagonalising, info==" << info << "\n"; }
}
iceig::iceig(sMat<double>&H, sMat<double>&iH)
{
   _Hsz = H.nc(); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; _V = 0;
   int info = ic_diag(H,iH,_zV,_E); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; }
}
iceig::iceig(int Hsz, double *E, double *V)
{
   _Hsz = Hsz; _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz]; _zV = 0;
   memcpy(_E,E,Hsz*sizeof(double)); memcpy(_V,V,Hsz*Hsz*sizeof(double));
}
iceig::iceig(int Hsz, double *E, complexdouble *zV)
{
   _Hsz = Hsz; _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; _V = 0;
   memcpy(_E,E,Hsz*sizeof(double)); memcpy(_zV,zV,Hsz*Hsz*sizeof(complexdouble));
}
iceig::iceig(int Hsz, double *E, complexdouble *zV, int step)
{
   _Hsz = Hsz; _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; _V = 0;
   memcpy(_E,E,Hsz*sizeof(double)); 
   int i; for(i=0; i<Hsz; i++) memcpy(&_zV[i*Hsz],&zV[i*Hsz+i+step],Hsz*sizeof(complexdouble));
}
// --------------------------------------------------------------------------------------------------------------- //
// Destructor
// --------------------------------------------------------------------------------------------------------------- //
iceig::~iceig()
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
}
// --------------------------------------------------------------------------------------------------------------- //
// Methods for iceig:: class
// --------------------------------------------------------------------------------------------------------------- //
void iceig::calc(sMat<double>&H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz];
   int info = ic_diag(H,_V,_E); 
   if(info!=0) { std::cerr << "iceig(H) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_V; _V=0; }
}
void iceig::calc(sMat<double>&H, sMat<double>&iH)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz];
   int info = ic_diag(H,iH,_zV,_E); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
}
void iceig::calc(int Hsz, complexdouble *H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = Hsz; _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz];
   int info = ic_diag(Hsz,H,_zV,_E); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
}
void iceig::lcalc(icpars &pars, sMat<double>&H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_V,0,_Hsz*_Hsz*sizeof(double));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm
   }
   if(nev>_Hsz) nev=_Hsz;
   int info = ic_leig(H,_V,_E,nev); 
   if(info!=0) { std::cerr << "iceig(H) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_V; _V=0; }
}
void iceig::lcalc(icpars &pars, sMat<double>&H, sMat<double>&iH)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_zV,0,_Hsz*_Hsz*sizeof(complexdouble));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm 
   }
   if(nev>_Hsz) nev=_Hsz;
   int info = ic_leig(H,iH,_zV,_E,nev); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
}
void iceig::lcalc(icpars &pars, complexdouble *H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = getdim(pars.n,pars.l); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz];
   memset(_E,0,_Hsz*sizeof(double)); memset(_zV,0,_Hsz*_Hsz*sizeof(complexdouble));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm 
   }
   if(nev>_Hsz) nev=_Hsz;
   int info = ic_leig(_Hsz,H,_zV,_E,nev); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
}
void iceig::pcalc(icpars &pars, complexdouble *zV, sMat<double> &J, sMat<double> &iJ)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = getdim(pars.n,pars.l); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz];
   memset(_E,0,_Hsz*sizeof(double)); memset(_zV,0,_Hsz*_Hsz*sizeof(complexdouble));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm 
   }
   complexdouble *zJmat = zmat2f(J,iJ); int info = ic_peig(_Hsz, zJmat, zV, _zV, _E, nev); free(zJmat);
   if(info!=0) { std::cerr << "iceig::pcalc() - Error\n"; }
}
void iceig::acalc(icpars &pars, sMat<double>&H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_V,0,_Hsz*_Hsz*sizeof(double));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm
   }
   if(nev>=_Hsz)   // We want all eigenvalues - better not to use the Arnoldi method
   {
      int info = ic_diag(H,_V,_E); 
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_V; _V=0; }
   }
   else
   {
      double *dH = H.f_array(); int info = ic_arpackeig(_Hsz,dH,_V,_E,nev); free(dH);
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_V; _V=0; }
   }
}
void iceig::acalc(icpars &pars, sMat<double>&H, sMat<double>&iH)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_zV,0,_Hsz*_Hsz*sizeof(complexdouble));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm
   }
   if(nev>=_Hsz)   // We want all eigenvalues - better not to use the Arnoldi method
   {
      int info = ic_diag(H,iH,_zV,_E); 
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
   }
   else
   {
      complexdouble *zH = zmat2f(H,iH); int info = ic_arpackeig(_Hsz,zH,_zV,_E,nev); free(zH);
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
   }
}
void iceig::acalc(icpars &pars, complexdouble *H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = getdim(pars.n,pars.l); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_zV,0,_Hsz*_Hsz*sizeof(complexdouble));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm 
   }
   if(nev>=_Hsz)   // We want all eigenvalues - better not to use the Arnoldi method
   {
      int info = ic_diag(_Hsz,H,_zV,_E); 
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; exit(EXIT_FAILURE); }
   }
   else
   {
      int info = ic_arpackeig(_Hsz,H,_zV,_E,nev); 
      if(info!=0) { std::cerr << "iceig::acalc() - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; exit(EXIT_FAILURE); }
   }
}
std::string iceig::strout() 
{
   int i,j;
   std::stringstream ss;
   for(i=0; i<_Hsz; i++) ss << _E[i]; ss << "\n";
   for(i=0; i<_Hsz; i++) { for(j=0; j<_Hsz; j++) ss << _V[j*_Hsz+i]; ss << "\n"; }
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------- //
// Constructor for class icmfmat::
// --------------------------------------------------------------------------------------------------------------- //
icmfmat::icmfmat()
{ 
   sMat<double> t; J.assign(6,t); 
   iflag.assign(6,0); iflag[2]=1; iflag[3]=1;
   _n = 1; _l = S; _num_op = 1;
}
icmfmat::icmfmat(int n, orbital l, int num_op, bool save_matrices, std::string density)
{
   _n = n; _l = l; _num_op = num_op; _density = density;
   sMat<double> t; J.assign(6,t); 
   iflag.assign(num_op>6?num_op:6,0); iflag[2]=1; iflag[3]=1;
   // Determines the filename strings for where the moment operator matrices are stored if previously calculated
   int nn = n; if(n>(2*l+1)) nn = (4*l+2)-n;
   char nstr[6]; char basename[255]; char Lfilestr[255], Sfilestr[255]; strcpy(basename,"results/mms/");
   if(save_matrices) {
   #ifndef _WINDOWS
   struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
      if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::(): Can't create mms dir, " << strerror(errno) << "\n";
   #else
   DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::(): Cannot create mms directory\n";
   #endif
   nstr[0] = (l==3?102:100); if(n<10) { nstr[1] = n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 76;   // 76 is ASCII for "L", 85=="U", 100=="d" and 102=="f"
   } else { strcpy(basename,"nodir/"); }
   nstr[1] = 49; // 49=="1"

   // Calculates the L and S operator matrix for the each directions
   sMat<double> Sp1, Sm1, Lp1, Lm1; 
   nstr[2]=120; nstr[3]=0; strcpy(Lfilestr,basename); strcat(Lfilestr,nstr); strcat(Lfilestr,".mm");   
   nstr[0]=83;             strcpy(Sfilestr,basename); strcat(Sfilestr,nstr); strcat(Sfilestr,".mm");
   J[1] = mm_gin(Lfilestr); J[0] = mm_gin(Sfilestr); 
   nstr[2]=121; nstr[3]=0; strcpy(Sfilestr,basename); strcat(Sfilestr,nstr); strcat(Sfilestr,".mm");   
   nstr[0]=76;             strcpy(Lfilestr,basename); strcat(Lfilestr,nstr); strcat(Lfilestr,".mm");
   J[3] = mm_gin(Lfilestr); J[2] = mm_gin(Sfilestr); 
   if(J[0].isempty() || J[1].isempty() || J[2].isempty() || J[3].isempty())
   { 
      racah_mumat(n,1,Lp1,Sp1,l); rmzeros(Sp1); rmzeros(Lp1);
      racah_mumat(n,-1,Lm1,Sm1,l); rmzeros(Sm1); rmzeros(Lm1);
      J[0] = (Sm1-Sp1)/sqrt(2); J[2] = (Sm1+Sp1)/sqrt(2); Sm1.clear(); Sp1.clear();   // Sx and Sy
      J[1] = (Lm1-Lp1)/sqrt(2); J[3] = (Lm1+Lp1)/sqrt(2); Lm1.clear(); Lp1.clear();   // Lx ans Ly
      mm_gout(J[2],Sfilestr); mm_gout(J[3],Lfilestr);
      nstr[2]=120; nstr[3]=0; strcpy(Lfilestr,basename); strcat(Lfilestr,nstr); strcat(Lfilestr,".mm");   
      nstr[0]=83;             strcpy(Sfilestr,basename); strcat(Sfilestr,nstr); strcat(Sfilestr,".mm");
      mm_gout(J[0],Sfilestr); mm_gout(J[1],Lfilestr);
   }

   nstr[0]=76; nstr[2]=122; nstr[3]=0; strcpy(Lfilestr,basename); strcat(Lfilestr,nstr); strcat(Lfilestr,".mm");
   nstr[0]=83;                         strcpy(Sfilestr,basename); strcat(Sfilestr,nstr); strcat(Sfilestr,".mm");
   J[4] = mm_gin(Sfilestr); J[5] = mm_gin(Lfilestr);                               // Sz and Lz
   if(J[4].isempty() || J[5].isempty()) { 
      racah_mumat(n,0,J[5],J[4],l); rmzeros(J[4]); rmzeros(J[5]); mm_gout(J[4],Sfilestr); mm_gout(J[5],Lfilestr); }

/* // Checks the moment operator matrices against those given by Chan and Lam.
// int ii,jj=0; sMat<double> mu; double g_s = 2.0023193043622; // electronic g-factor
// chanlam_mumat(n,1,mu,l); for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    if(fabs(-mu(ii,jj)-J[1](ii,jj)-g_s*J[0](ii,jj))>10*DBL_EPSILON) { std::cerr << "icmfmat: Magnetic moment operator x does not agree.\n"; break; }
//    if(ii==mu.nr() && jj==mu.nc()) std::cerr << "icmfmat: Magnetic moment operator x agrees.\n";
// chanlam_mumat(n,2,mu,l); for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    if(fabs(-mu(ii,jj)-J[3](ii,jj)-g_s*J[2](ii,jj))>10*DBL_EPSILON) { std::cerr << "icmfmat: Magnetic moment operator y does not agree.\n"; break; }
//    if(ii==mu.nr() && jj==mu.nc()) std::cerr << "icmfmat: Magnetic moment operator y agrees.\n";
// chanlam_mumat(n,3,mu,l); for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    if(fabs(-mu(ii,jj)-J[5](ii,jj)-g_s*J[4](ii,jj))>10*DBL_EPSILON) { std::cerr << "icmfmat: Magnetic moment operator z does not agree.\n"; break; }
//    if(ii==mu.nr() && jj==mu.nc()) std::cerr << "icmfmat: Magnetic moment operator z agrees.\n";
   double sumcheck;
   chanlam_mumat(n,1,mu,l); sumcheck = 0.; for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    std::cout << -mu(ii,jj) << "\t" << J[1](ii,jj)+g_s*J[0](ii,jj)  << "\t" << fabs(-mu(ii,jj)-J[1](ii,jj)-g_s*J[0](ii,jj)) << "\n";
      sumcheck += fabs(-mu(ii,jj)-J[1](ii,jj)-g_s*J[0](ii,jj)); std::cout << "Moment Matrix Check: sum(-mu_x(ChanLam) - (Lx+gSx)) = " << sumcheck << "\n";
   chanlam_mumat(n,2,mu,l); sumcheck = 0.; for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    std::cout << -mu(ii,jj) << "\t" << J[3](ii,jj)+g_s*J[2](ii,jj)  << "\t" << fabs(-mu(ii,jj)-J[3](ii,jj)-g_s*J[2](ii,jj)) << "\n";
      sumcheck += fabs(-mu(ii,jj)-J[3](ii,jj)-g_s*J[2](ii,jj)); std::cout << "Moment Matrix Check: sum(-mu_y(ChanLam) - (Ly+gSy)) = " << sumcheck << "\n";
   chanlam_mumat(n,3,mu,l); sumcheck = 0.; for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
      sumcheck += fabs(-mu(ii,jj)-J[5](ii,jj)-g_s*J[4](ii,jj)); std::cout << "Moment Matrix Check: sum(-mu_z(ChanLam) - (Lz+gSz)) = " << sumcheck << "\n"; */
}
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the mean field matrix sum_i (H_i*J_i)
// --------------------------------------------------------------------------------------------------------------- //
void icmfmat::Jmat(sMat<double>&Jmat, sMat<double>&iJmat, std::vector<double>&gjmbH, bool save_matrices)
{
   int i; Jmat.zero(J[0].nr(),J[0].nc()); iJmat.zero(J[0].nr(),J[0].nc()); 
   if(_num_op<(int)gjmbH.size()) { iflag.resize(_num_op,0); _num_op = (int)gjmbH.size(); }
   for(i=0; i<(_num_op>6?6:_num_op); i++)
      if(fabs(gjmbH[i])>DBL_EPSILON*100) { if(iflag[i]==1) iJmat += J[i]*gjmbH[i]; else Jmat += J[i]*gjmbH[i]; }
   // Higher order than dipole operators needed
   if(_num_op>6)
   {
      char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
      if(save_matrices) {
      #ifndef _WINDOWS
      struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
         if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::Jmat(): Can't create mms dir, " << strerror(errno) << "\n";
      #else
      DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
         if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::Jmat(): Cannot create mms directory\n";
      #endif
      nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
      strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
      } else { strcpy(basename,"nodir/"); }
#define NSTR(K,Q) nstr[1] = K+48; nstr[2] = Q+48; nstr[3] = 0
#define MSTR(K,Q) nstr[1] = K+48; nstr[2] = 109;  nstr[3] = Q+48; nstr[4] = 0
      // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
      int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
      int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
    //int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
      sMat<double> Upq,Umq; double redmat; int n = _n; //if(n>(2*_l+1)) n = 4*_l+2-n; 

      for(i=6; i<_num_op; i++)
      {
         if(q[i]<0) iflag[i]=1; 
         if (fabs(gjmbH[i])>DBL_EPSILON) 
         {
            redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[i],2*_l,0,0,0);// * wy2stev(i);
            if(k[i]%2==1) continue;   // Using the above reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
            if(k[i]>4 && _l==D) continue;
            NSTR(k[i],abs(q[i])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[i],abs(q[i]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
            if(q[i]==0) { Jmat += Upq * (gjmbH[i]*redmat); continue; }
            MSTR(k[i],abs(q[i])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[i],-abs(q[i]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
            if(q[i]<0) { 
            // if((q[i]%2)==0) iJmat += (Upq - Umq) * (gjmbH[i]*redmat); else iJmat += (Upq + Umq) * (gjmbH[i]*redmat); } changed MR 15.12.09
               if((q[i]%2)==0) iJmat -= (Upq - Umq) * (gjmbH[i]*redmat); else iJmat -= (Upq + Umq) * (gjmbH[i]*redmat); }
            else {
               if((q[i]%2)==0)  Jmat += (Upq + Umq) * (gjmbH[i]*redmat); else  Jmat += (Upq - Umq) * (gjmbH[i]*redmat); } 
         } 
      }
   }
}
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the expectation values <V|J|V>exp(-beta*T) given a set of eigenstates
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> icmfmat::expJ(iceig &VE, double T, std::vector< std::vector<double> > &matel, bool save_matrices)
{
   double *vt=0, Z=0., U=0.; complexdouble *zt=0, zme;
   std::vector<double> E, ex((_num_op>6?_num_op:6)+2,0.), me, eb; matel.clear();
   int iJ, ind_j, Esz, Hsz=VE.Hsz(), incx=1; 
   if(Hsz!=J[0].nr()) { std::cerr << "icmfmat::expJ() - Hamiltonian matrix size not same as mean field operator!\n"; return E; }
   sMat<double> zeroes; zeroes.zero(J[0].nr(),J[0].nc());
   double alpha = 1, beta = 0; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';
   // Checks that the eigenvalues are orthonormal
/* char transa='C', transb='N'; double summm=0.;
   if(VE.iscomplex())
   {
      complexdouble *zmm = (complexdouble*)malloc(Hsz*Hsz*sizeof(complexdouble)); 
      complexdouble *vet = (complexdouble*)malloc(Hsz*Hsz*sizeof(complexdouble)); memcpy(vet,VE.zV(0),Hsz*Hsz*sizeof(complexdouble));
      F77NAME(zgemm)(&transa, &transb, &Hsz, &Hsz, &Hsz, &zalpha, vet, &Hsz, VE.zV(0), &Hsz, &zbeta, zmm, &Hsz);
      for(int ii=0; ii<Hsz; ii++) { zmm[ii*Hsz+ii].r-=1.; summm += F77NAME(dzasum)(&Hsz, &zmm[ii*Hsz], &incx); if(VE.E(ii+1)==0) break; }
      std::cout << "#ic1ion: Sum(V^TV-I) = " << summm << "\n";
      free(zmm); free(vet);
   }
   else
   {
      double *dmm = (double*)malloc(Hsz*Hsz*sizeof(double)); 
      double *vet = (double*)malloc(Hsz*Hsz*sizeof(double)); memcpy(vet,VE.V(0),Hsz*Hsz*sizeof(double));
      F77NAME(dgemm)(&transa, &transb, &Hsz, &Hsz, &Hsz, &alpha, vet, &Hsz, VE.V(0), &Hsz, &beta, dmm, &Hsz);
      for(int ii=0; ii<Hsz; ii++) { dmm[ii*Hsz+ii]-=1.; summm += F77NAME(dasum)(&Hsz, &dmm[ii*Hsz], &incx); if(VE.E(ii+1)==0) break; }
      std::cout << "#ic1ion: Sum(V^TV-I) = " << summm << "\n";
      free(dmm); free(vet);
   }*/

   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<J[0].nr(); Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || VE.E(Esz+1)==0) break; }

   if (T<0){Esz=(int)(-T);printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                         printf ("Number   Excitation Energy\n");
     for (ind_j=0;ind_j<Esz;++ind_j) printf ("%i    %4.4g meV\n",ind_j+1,E[ind_j]);
     } // MR 10.9.2010

   for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++) 
      if(fabs(VE.zV(ii,jj).r*VE.zV(ii,jj).r+VE.zV(ii,jj).i*VE.zV(ii,jj).i)<DBL_EPSILON*100000) 
      {
         VE.zV(ii,jj).r=0.; VE.zV(ii,jj).i=0.;  
      }  
   // For first run calculate also the partition function and internal energy
   me.assign(Esz,0.); eb.assign(Esz,0.); Z=0.;
   if(!VE.iscomplex()) 
   {
      double *fJmat=J[0].f_array(); vt = (double*)malloc(Hsz*sizeof(double)); 
      for(ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, VE.V(ind_j), &incx, &beta, vt, &incx);
#ifdef _G77 
         F77NAME(ddot)(me[ind_j],&Hsz, VE.V(ind_j), &incx, vt, &incx);
#else
         me[ind_j] = F77NAME(ddot)(&Hsz, VE.V(ind_j), &incx, vt, &incx);
#endif
//MR 10.9.2010
     if (T<0)
     { char instr[MAXNOFCHARINLINE];
      printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
       fgets(instr, MAXNOFCHARINLINE, stdin);
       eb[ind_j]=strtod(instr,NULL);
     }
      else
        { eb[ind_j] = exp(-E[ind_j]/(KB*T));} ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MRend 10.9.2010                                                                   !!!!    -----------------!!!!
      }
      free(fJmat); free(vt); matel.push_back(me); ex[0]/=Z; U/=Z;
   }
   else
   {
      complexdouble *zJmat;
      zeroes.zero(J[0].nr(),J[0].nc()); if(iflag[0]==0) zJmat=zmat2f(J[0],zeroes); else zJmat = zmat2f(zeroes,J[0]);
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      for(ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
#ifdef _G77 
         F77NAME(zdotc)(&zme, &Hsz, VE.zV(ind_j), &incx, zt, &incx);
#else
         zme = F77NAME(zdotc)(&Hsz, VE.zV(ind_j), &incx, zt, &incx);
#endif
         me[ind_j] = zme.r;
//MR 10.9.2010
     if (T<0)
     { char instr[MAXNOFCHARINLINE];
      printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
       fgets(instr, MAXNOFCHARINLINE, stdin);
       eb[ind_j]=strtod(instr,NULL);
     }
      else
        { eb[ind_j] = exp(-E[ind_j]/(KB*T));} ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MRend 10.9.2010
      }
      free(zJmat); free(zt); matel.push_back(me); ex[0]/=Z; U/=Z;
   }

   char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
   if(save_matrices) {
   #ifndef _WINDOWS
   struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
      if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::expJ(): Can't create mms dir, " << strerror(errno) << "\n";
   #else
   DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::expJ(): Cannot create mms directory\n";
   #endif
   nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
   } else { strcpy(basename,"nodir/"); }
   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   sMat<double> Upq,Umq; double redmat; int n = _n; //if(n>(2*_l+1)) n = 4*_l+2-n; 

   if(!_density.empty()) 
   { 
      std::cout << "Calculating the expectation of the moment density operator ";
      if(_density.find("l")!=std::string::npos || _density.find("s")!=std::string::npos) std::cout << _density << "\n"; 
      else
      {
         if(_density.find("1")!=std::string::npos) std::cout << "Sx\n";
         if(_density.find("2")!=std::string::npos) std::cout << "Lx\n";
         if(_density.find("3")!=std::string::npos) std::cout << "Sy\n";
         if(_density.find("4")!=std::string::npos) std::cout << "Ly\n";
         if(_density.find("5")!=std::string::npos) std::cout << "Sz\n";
         if(_density.find("6")!=std::string::npos) std::cout << "Lz\n";
      }
   }

   // Rest of the runs only calculate the new matrix elements
   for(iJ=1; iJ<(_num_op>6?_num_op:6); iJ++)
   {
      me.assign(Esz,0.);
      // Using the above reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
      if((iJ>6 && k[iJ]%2==1) || (k[iJ]>4 && _l==D)) { matel.push_back(me); continue; }
      if(!VE.iscomplex() && _density.empty())
      {
         if(iflag[iJ]==0) {
            double *fJmat; vt = (double*)malloc(Hsz*sizeof(double)); 
            if(iJ<6) 
               fJmat=J[iJ].f_array(); 
            else 
            {
               NSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
               Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[iJ],abs(q[iJ]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
               MSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
               Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[iJ],-abs(q[iJ]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
               redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[iJ],2*_l,0,0,0);// * wy2stev(iJ);
//             if(q[iJ]<0) { if((q[iJ]%2)==0) Upq -= Umq; else Upq += Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } changed MR 15.12.09
               if(q[iJ]<0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; }
               Upq *= redmat; fJmat = Upq.f_array();
            }
            for(ind_j=0; ind_j<Esz; ind_j++)
            {  // Calculates the matrix elements <Vi|J.H|Vi>
               F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, VE.V(ind_j), &incx, &beta, vt, &incx);
#ifdef _G77 
               F77NAME(ddot)(me[ind_j],&Hsz, VE.V(ind_j), &incx, vt, &incx);
#else
               me[ind_j] = F77NAME(ddot)(&Hsz, VE.V(ind_j), &incx, vt, &incx);
#endif
               ex[iJ]+=me[ind_j]*eb[ind_j];
            }
            free(fJmat); free(vt); matel.push_back(me); ex[iJ]/=Z; 
         } 
         else { me.assign(Esz,0.); matel.push_back(me); }
      }
      else
      {
         complexdouble *zJmat; zeroes.zero(J[0].nr(),J[0].nc());
         if(iJ<6) 
         {
            if(iflag[iJ]==0) zJmat=zmat2f(J[iJ],zeroes); else zJmat = zmat2f(zeroes,J[iJ]);
         }
         else 
         {
         //   if(!_density.empty()) { zJmat = balcar_Mq(_density,k[iJ],q[iJ],_n,_l); } else {
            NSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[iJ],abs(q[iJ]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
            MSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[iJ],-abs(q[iJ]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
            redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[iJ],2*_l,0,0,0);// * wy2stev(iJ);
//          if(q[iJ]<0) { if((q[iJ]%2)==0) Upq -= Umq; else Upq += Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } changed MR 15.12.09
            if(q[iJ]<0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; }
            Upq *= redmat; if(iflag[iJ]==0) zJmat=zmat2f(Upq,zeroes); else zJmat = zmat2f(zeroes,Upq);
//      }
         }
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         for(ind_j=0; ind_j<Esz; ind_j++)
         {  // Calculates the matrix elements <Vi|J.H|Vi>
            F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
#ifdef _G77 
            F77NAME(zdotc)(&zme, &Hsz, VE.zV(ind_j), &incx, zt, &incx);
#else
            zme = F77NAME(zdotc)(&Hsz, VE.zV(ind_j), &incx, zt, &incx);
#endif
            me[ind_j] = zme.r;
            ex[iJ]+=me[ind_j]*eb[ind_j];
         }
         free(zJmat); free(zt); matel.push_back(me); ex[iJ]/=Z;
      }
      if(fabs(ex[iJ])<DBL_EPSILON) ex[iJ]=0.; 
   }
   ex[iJ] = log(Z)-VE.E(0)/(KB*T); ex[iJ+1] = U;
   return ex;
}
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the matrix M_ab=<i|Ja|j><j|Jb|i>{exp(-beta_i*T)-exp(-beta_j*T)} for some state i,j
// --------------------------------------------------------------------------------------------------------------- //
void icmfmat::Mab(sMat<double>&Mab, sMat<double>&iMab, iceig&VE, double T, int i, int j,int pr,float & delta, bool save_matrices)
{
   double *vt=0, Z=0., therm; complexdouble *zt=0, zme; zme.r=0; zme.i=0.;
   int sz = (_num_op>6?_num_op:6);
   std::vector<double> mij(sz,0.);//, mji(6,0.);
   std::vector<complexdouble> zij(sz,zme);//, zji(6,zme);
   Mab.zero(sz,sz); iMab.zero(sz,sz);
   int iJ, jJ, Hsz=VE.Hsz(), incx=1; 
   if(Hsz!=J[0].nr()) { std::cerr << "icmfmat::Mab() - Hamiltonian matrix size not same as mean field operator!\n"; return; }
   sMat<double> zeroes; zeroes.zero(J[0].nr(),J[0].nc());
   double alpha = 1, beta = 0; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   complexdouble *zJmat=0;
   char uplo = 'U';

   char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
   if(save_matrices) {
   #ifndef _WINDOWS
   struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
      if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::Mab(): Can't create mms dir, " << strerror(errno) << "\n";
   #else
   DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::Mab(): Cannot create mms directory\n";
   #endif
   nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
   } else { strcpy(basename,"nodir/"); }
   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
                 
   sMat<double> Upq,Umq; double redmat; int n = _n; if(n>(2*_l+1)) n = 4*_l+2-n; 

   // Calculates the matrix elements: <i|Ja|j> and <j|Ja|i> for each of the six Ja's
   for(iJ=0; iJ<sz; iJ++)
   {
//    if(k[iJ]%2==1) { if(VE.iscomplex()) { zij[iJ].r=0.; zij[iJ].i=0.; } else mij[iJ]=0.; continue; }
      if(iJ>=6)
      {
         NSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[iJ],abs(q[iJ]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
         MSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[iJ],-abs(q[iJ]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
         redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[iJ],2*_l,0,0,0);
//       if(q[iJ]<0) { if((q[iJ]%2)==0) Upq -= Umq; else Upq += Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } changed MR 15.12.09
         if(q[iJ]<0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; }
         Upq *= redmat;
      }

      if(!VE.iscomplex() && iflag[iJ]==0)
      {
         vt = (double*)malloc(Hsz*sizeof(double)); 
         double *fJmat; if(iJ>=6) fJmat=Upq.f_array(); else fJmat=J[iJ].f_array();
         F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, VE.V(j), &incx, &beta, vt, &incx);
#ifdef _G77 
         F77NAME(ddot)(mij[iJ], &Hsz, VE.V(i), &incx, vt, &incx); zij[iJ].r = mij[iJ];
#else
         mij[iJ] = F77NAME(ddot)(&Hsz, VE.V(i), &incx, vt, &incx); zij[iJ].r = mij[iJ];
#endif
         free(fJmat); free(vt);
      } 
      else
      {
         zeroes.zero(J[0].nr(),J[0].nc());
         if(iJ>=6) { if(iflag[iJ]==0) zJmat=zmat2f(Upq,zeroes);   else zJmat = zmat2f(zeroes,Upq); }
         else      { if(iflag[iJ]==0) zJmat=zmat2f(J[iJ],zeroes); else zJmat = zmat2f(zeroes,J[iJ]); }
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(j), &incx, &zbeta, zt, &incx);
#ifdef _G77 
         F77NAME(zdotc)(&zij[iJ], &Hsz, VE.zV(i), &incx, zt, &incx);
#else
         zij[iJ] = F77NAME(zdotc)(&Hsz, VE.zV(i), &incx, zt, &incx);
#endif
//       int k;for(k=0;k<Hsz;++k)printf("%6.3f %+6.3f i  ",VE.zV(j)[k].r,VE.zV(j)[k].i);
         free(zJmat); free(zt);
      }
   }

   // // Calculates the matrix M_ab
   // for(iJ=0; iJ<6; iJ++)
   //    for(jJ=0; jJ<6; jJ++)
   //       Mab(i+1,j+1) = (mij[iJ]*mji[jJ]) * ( exp(-(VE.E(i)-VE.E(0))/(KB*T)) - exp(-(VE.E(j)-VE.E(0))/(KB*T)) );


   if(i==j) {//subtract thermal expectation value from zij=zii
            std::vector< std::vector<double> > matel;
            std::vector<double> vJ = expJ(VE,T,matel,save_matrices);
            for(iJ=0; iJ<sz; iJ++)zij[iJ].r-=vJ[iJ];
            }

   // Calculates the matrix M_ab and iM_ab
   for(iJ=0; iJ<sz; iJ++)
      for(jJ=0; jJ<sz; jJ++)
      {  
         Mab(iJ+1,jJ+1) = (zij[iJ].r*zij[jJ].r + zij[iJ].i*zij[jJ].i);
         iMab(iJ+1,jJ+1) = (-zij[iJ].r*zij[jJ].i + zij[iJ].i*zij[jJ].r);
      }

   delta = VE.E(j)-VE.E(i);
   if(delta<-0.000001)
   {
      std::cerr << "ERROR module ic1ion - dmcalc: energy gain delta gets negative\n"; 
      exit(EXIT_FAILURE);
   }
   if(j==i)delta=-SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc

   // Calculates the partition function
   for(iJ=0; iJ<Hsz; iJ++) { therm = exp(-(VE.E(iJ)-VE.E(0))/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }

   // do some printout if wishes and set correct occupation factor
   if (delta>SMALL)
   {
      therm = exp(-(VE.E(i)-VE.E(0))/(KB*T)) - exp(-(VE.E(j)-VE.E(0))/(KB*T));
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
         printf(" |<%i|Ja|%i>|^2=%6.3f\n |<%i|Jb|%i>|^2=%6.3f\n |<%i|Jc|%i>|^2=%6.3f\n",i+1,j+1,Mab(1,1),i+1,j+1,Mab(2,2),i+1,j+1,Mab(3,3));
         printf(" |<%i|Jd|%i>|^2=%6.3f\n |<%i|Je|%i>|^2=%6.3f\n |<%i|Jf|%i>|^2=%6.3f\n",i+1,j+1,Mab(4,4),i+1,j+1,Mab(5,5),i+1,j+1,Mab(6,6));
         printf(" n%i-n%i=%6.3f\n",i,j,therm / Z);
      }
   }
   else
   {
      therm = exp(-(VE.E(i)-VE.E(0))/(KB*T))/(KB*T);    // quasielastic scattering has not wi-wj but wj*epsilon/kT
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
         printf(" |<%i|Ja-<Ja>|%i>|^2=%6.3f\n |<%i|Jb-<Jb>|%i>|^2=%6.3f\n |<%i|Jc-<Jc>|%i>|^2=%6.3f\n",i+1,j+1,Mab(1,1),i+1,j+1,Mab(2,2),i+1,j+1,Mab(3,3));
         printf(" |<%i|Jd-<Jd>|%i>|^2=%6.3f\n |<%i|Je-<Je>|%i>|^2=%6.3f\n |<%i|Jf-<Jf>|%i>|^2=%6.3f\n",i+1,j+1,Mab(4,4),i+1,j+1,Mab(5,5),i+1,j+1,Mab(6,6));
         printf(" n%i=%6.3f\n",i,(KB*T)*therm/Z);
      }
   }

   // multiply matrix Mab by occupation factor
   for(iJ=0; iJ<sz; iJ++)
      for(jJ=0; jJ<sz; jJ++) { Mab(iJ+1,jJ+1) *= therm/Z; iMab(iJ+1,jJ+1) *= therm/Z; }

}

//--------------------------------------------------------------------------------------------------------------
std::vector<double> icmfmat::orbmomdensity_expJ(iceig &VE,int xyz, double T, std::vector< std::vector<double> > &matel, bool save_matrices)
{
   return spindensity_expJ(VE,-xyz, T, matel, save_matrices);
}

//--------------------------------------------------------------------------------------------------------------
std::vector<double> icmfmat::spindensity_expJ(iceig &VE,int xyz, double T, std::vector< std::vector<double> > &matel, bool save_matrices)
{
   double *vt=0, Z=0., U=0.; complexdouble *zt=0, zme;
   std::vector<double> E, ex((_num_op>6?_num_op:6)+2,0.), me, eb; matel.clear();
   int iJ, ind_j, Esz, Hsz=VE.Hsz(), incx=1;
   if(Hsz!=J[0].nr()) { std::cerr << "icmfmat::expJ() - Hamiltonian matrix size not same as mean field operator!\n"; return E; }
   sMat<double> zeroes; zeroes.zero(J[0].nr(),J[0].nc());
   double alpha = 1, beta = 0; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';
   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<J[0].nr(); Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || VE.E(Esz+1)==0) break; }
   if (T<0){Esz=(int)(-T);printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                         printf ("Number   Excitation Energy\n");
     for (ind_j=0;ind_j<Esz;++ind_j) printf ("%i    %4.4g meV\n",ind_j+1,E[ind_j]);
     } // MR 17.9.2010
   for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++)
      if(fabs(VE.zV(ii,jj).r*VE.zV(ii,jj).r+VE.zV(ii,jj).i*VE.zV(ii,jj).i)<DBL_EPSILON*100000)
      {
         VE.zV(ii,jj).r=0.; VE.zV(ii,jj).i=0.;
      }

   // For first run calculate also the partition function and internal energy
   me.assign(Esz,0.); eb.assign(Esz,0.); Z=0.;
   if(!VE.iscomplex())
   {
      double *fJmat=J[0].f_array(); vt = (double*)malloc(Hsz*sizeof(double));
      for(ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, VE.V(ind_j), &incx, &beta, vt, &incx);
#ifdef _G77
         F77NAME(ddot)(me[ind_j],&Hsz, VE.V(ind_j), &incx, vt, &incx);
#else
         me[ind_j] = F77NAME(ddot)(&Hsz, VE.V(ind_j), &incx, vt, &incx);
#endif
 //MR 17.9.2010
     if (T<0)
     { char instr[MAXNOFCHARINLINE];
      printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
       fgets(instr, MAXNOFCHARINLINE, stdin);
       eb[ind_j]=strtod(instr,NULL);
     }
      else
        { eb[ind_j] = exp(-E[ind_j]/(KB*T));} ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MRend 17.9.2010
     }
      free(fJmat); free(vt); matel.push_back(me); ex[0]/=Z; U/=Z;
   }
   else
   {
      complexdouble *zJmat;
      zeroes.zero(J[0].nr(),J[0].nc()); if(iflag[0]==0) zJmat=zmat2f(J[0],zeroes); else zJmat = zmat2f(zeroes,J[0]);
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      for(ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
#ifdef _G77
         F77NAME(zdotc)(&zme, &Hsz, VE.zV(ind_j), &incx, zt, &incx);
#else
         zme = F77NAME(zdotc)(&Hsz, VE.zV(ind_j), &incx, zt, &incx);
#endif
         me[ind_j] = zme.r;
//         eb[ind_j] = exp(-E[ind_j]/(KB*T)); ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MR 17.9.2010
     if (T<0)
     { char instr[MAXNOFCHARINLINE];
      printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
       fgets(instr, MAXNOFCHARINLINE, stdin);
       eb[ind_j]=strtod(instr,NULL);
     }
      else
        { eb[ind_j] = exp(-E[ind_j]/(KB*T));} ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MRend 17.9.2010                                                                   !!!!    -----------------!!!!
      }
      free(zJmat); free(zt); matel.push_back(me); ex[0]/=Z; U/=Z;
   }

   char nstr[6]; char basename[255]; strcpy(basename,"results/mms/");
   if(save_matrices) {
   #ifndef _WINDOWS
   struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
      if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::expJ(): Can't create mms dir, " << strerror(errno) << "\n";
   #else
   DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY))
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::expJ(): Cannot create mms directory\n";
   #endif
   nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
   } else { strcpy(basename,"nodir/"); }
   int k[] = {0,1, 1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   sMat<double> Upq,Umq;  //if(n>(2*_l+1)) n = 4*_l+2-n;
   if(xyz>0)
   {std::cout << "Calculating the expectation values of the spin density operator\n";}
   else
   {std::cout << "Calculating the expectation values of the orbital moment density operator \n";}

   // Rest of the runs only calculate the new matrix elements
   for(iJ=0; iJ<(_num_op>6?_num_op:6); iJ++)
   {ex[iJ]=0;
      me.assign(Esz,0.);
      // Using the above reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
      if((k[iJ]%2==1) || (k[iJ]>4 && _l==D)) { matel.push_back(me); continue; }
         complexdouble *zJmat; zeroes.zero(J[0].nr(),J[0].nc());
         
         zJmat = balcar_Mq(xyz,k[iJ],q[iJ],_n,_l); // minus sign stands for orbital density coeff
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         for(ind_j=0; ind_j<Esz; ind_j++)
         {  // Calculates the matrix elements <Vi|J.H|Vi>
            F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
#ifdef _G77
            F77NAME(zdotc)(&zme, &Hsz, VE.zV(ind_j), &incx, zt, &incx);
#else
            zme = F77NAME(zdotc)(&Hsz, VE.zV(ind_j), &incx, zt, &incx);
#endif
            me[ind_j] = zme.r;
            ex[iJ]+=me[ind_j]*eb[ind_j];
         }
         free(zJmat); free(zt); matel.push_back(me); ex[iJ]/=Z;
      
      if(fabs(ex[iJ])<DBL_EPSILON) ex[iJ]=0.;
   }
   //ex[iJ] = log(Z)-VE.E(0)/(KB*T); ex[iJ+1] = U;
   return ex;
}
