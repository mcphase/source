#include "martin.h"
#include "singleion_module.hpp"


class phonon_module : public singleion_module
{private:  
  double charge;
public:
phonon_module(const char * filename="");
phonon_module(const phonon_module & pp);
~phonon_module();
//void Icalc() override;
 bool Icalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext,double & gJ,Vector & ABC,char * sipffilename ,double & lnZ,double & U)
           override;
// bool IMcalc(Matrix &mom, Vector & T, Vector &  Hxc,Vector & Hext,double & gJ,Vector & ABC,char * sipffilename ,Vector & lnZ,Vector & U)
//           override;
 int du1calc(int &transitionnumber,double &T,Vector&Hxc,Vector&Hext,double&gJ,Vector&ABC,char * sipffilename,ComplexVector&u1,float&delta,int&n,int&nd,ComplexMatrix&ests)
           override;
// bool estates(ComplexMatrix & est, Vector & Hxc,Vector & Hext, double & g_J, double &T,       // Calculates the energy and wavefunctions
//                Vector &ABC, char *sipffilename)
//           override;
 //bool opmat(int &ni,                      // ni     which operator 0=Hamiltonian, 1,2,3=J1,J2,J3
 //            char *sipffilename,         // Single ion properties filename
 //            Vector &Hxc,                 // Hext  vector of external field [meV]
 //            Vector &Hext,                // Hxc   vector of exchange field [meV]
 //                                         // on output   
 //            Matrix &outmat)       
 //          override;
 bool pcalc(Vector & u0,double & T, Vector &Fxc, Vector & Hext,double & g_J, Vector & MODPAR,char * sipffilename) 
           override;
 int dP1(int & tn,double & T,Vector & Fxc, Vector & Hext,
                       double & g_J,Vector & MODPAR, char * sipffilename,
                       ComplexVector & P1,float & maxE,ComplexMatrix & est)
           override;
bool pelcalc(Vector & u0,double & T, Vector &Fxc, Vector & Hext,double & g_J, Vector & MODPAR,char * sipffilename) 
           override;
 int dpel1(int & tn,double & T,Vector & Fxc, Vector & Hext,
                       double & g_J,Vector & MODPAR, char * sipffilename,
                       ComplexVector & pel1,float & maxE,ComplexMatrix & est)
           override;
 bool mcalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext,double & gJ,Vector & ABC,char * sipffilename)
           override;

/* 
  bool mMcalc(Matrix &mom, Vector & T, Vector &  Hxc,Vector & Hext,double & gJ,Vector & ABC,char * sipffilename)
           override;
 int dm1(int &transitionnumber,double &T,Vector&Hxc,Vector&Hext,double&gJ,Vector&ABC,char * sipffilename,ComplexVector&u1,float&delta,ComplexMatrix&ests)override;
 bool Lcalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext,double & gJ,Vector & ABC,char * sipffilename)
           override;
 bool LMcalc(Matrix &mom, Vector & T, Vector &  Hxc,Vector & Hext,double & gJ,Vector & ABC,char * sipffilename)
           override;
 int dL1(int &transitionnumber,double &T,Vector&Hxc,Vector&Hext,double&gJ,Vector&ABC,char * sipffilename,ComplexVector&u1,float&delta,ComplexMatrix&ests)override;
 bool Scalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext,double & gJ,Vector & ABC,char * sipffilename)
           override;
 bool SMcalc(Matrix &mom, Vector & T, Vector &  Hxc,Vector & Hext,double & gJ,Vector & ABC,char * sipffilename)
           override;
 int dS1(int &transitionnumber,double &T,Vector&Hxc,Vector&Hext,double&gJ,Vector&ABC,char * sipffilename,ComplexVector&u1,float&delta,ComplexMatrix&ests)override;
 bool mqcalc(ComplexVector &Mq,      // Output expectation values -2[<Q>_{x} <Q>_y <Q>_z]
                  double &th, double &ph, // Input polar and azimuth angles theta and phi
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est)     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
                   override;
        int                              // Returns total number of transitions
             dmq1(int &tn,                // Input transition number |tn|. If tn>0 omit printout. If tn<0 print info.
                  double &th,             // Input zenith angle (with the z==b axis) in radians.
                  double &ph,             // Input azimuth angle (with the x==a axis, to projection in x-y plane).
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est,     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
                  double &T,              // Input temperature (K)
                  ComplexVector & mq1,    // input mq1(1)= ninit + i pinit   
                                          // Output transition vector, mq1=<-|M(Q)|+> sqrt(n- - n+) in units of MU_B
                  double & maxE)          // input maxE maximal transition energy
                       override;
 bool chargedensity_coeff(
                      Vector &mom,         // Output single ion moments == expectation values of
                                           //    of Zlm R^2(r) at a given temperature T and  effective field H
                      double &T,           // Input scalar temperature
                      Vector &Hxc,         // Input vector of exchange fields (meV) 
                      Vector &Hext,        // Input vector of external field (T) 
                      double & g_J,    // Input Lande g-factor
                      Vector & ABC,    // Input vector of parameters from single ion property file
                      char *sipffilename, // Single ion properties filename
                      ComplexMatrix &est)  // Input/output eigenstate matrix 
                    override;
 //bool ro_calc(double & ro,double & teta,double & fi,double & R,Vector & moments,double & gJ,Vector & ABC,char * sipffilename) 
 //                  override;
 bool spindensity_coeff(Vector & mom, int & xyz,double &T,Vector &Hxc,Vector&Hext,       // Calc. coeffs. of expansion of spindensity 
                double &gJ,Vector &ABC, char *sipffile, ComplexMatrix &est)
                  override;
 int dchargedensity_coeff1(int &tn, double &T, Vector &Hxc,Vector&Hext, double &g_J,     // Calculates the transition
                  Vector &ABC, char *sipffilename, ComplexVector & dc1, float &delta,      //   matrix elements of the chargedensity coefficients
                  ComplexMatrix &est)
                 override;
 int dspindensity_coeff1(int &tn, double &T, Vector &Hxc,Vector&Hext, double &g_J,     // Calculates the transition
                  Vector &ABC, char *sipffilename, ComplexVector & dc1,int & xyz, float &delta,      //   matrix elements of the chargedensity coefficients
                  ComplexMatrix &est)
                  override;
 bool orbmomdensity_coeff(Vector & mom,int & xyz, double  & T,Vector &Hxc,Vector&Hext,     // Calc. coeffs. of expansion of orbital moment density
                  double &gJ,Vector &ABC, char *sipffile, ComplexMatrix &est)              //   in terms of Zlm F(r) at given T / H_eff
                  override;        
 int dorbmomdensity_coeff1(int &tn, double &T, Vector &Hxc,Vector&Hext, double &g_J,     // Calculates the transition
                  Vector &ABC, char *sipffilename, ComplexVector & dc1,int & xyz, float &delta,      //   matrix elements of the chargedensity coefficients
                  ComplexMatrix &est)
                  override;

*/
 //int drixs1(int & transitionnumber,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,ComplexMatrix & ests,double & T,ComplexVector & drixs,double & delta) 
 //                 override;
};
