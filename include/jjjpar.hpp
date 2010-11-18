// jjjpar is a class to store all parameters associated
// with a specific ion, for example CEF parameters
// and exchange parameters 
// it is used by many programs in the package
// moreover, it loads also the user defined single ion module functions (linux only)
#ifndef JJJPAR
#define JJJPAR

typedef void fnc_t();

#ifdef __linux__
#include<dlfcn.h>
#else
#include <windows.h>
#endif

#include<martin.h>
#include<ionpars.hpp>
#include<myev.h>
#include<stdlib.h>


class par;
class jjjpar
{
public:
//******************************************************88
// basic parameters
  char * cffilename; // single ion parameter filename
  char * modulefilename; // module name
    double J(); // returns total angular momentum if possible
   Vector &  tetan(); //returns stevens parameters if possible

  Vector xyz,mom; // atom position, moment
  int paranz;   // number of exchange parameters
  int nofcomponents; // number of moments (components of moment vector)
  double gJ;
  Matrix *jij; // exchange constants 
  Vector *dn; // exchange - coordinates of neighbors 
  int *sublattice; // sublattice of neighbours
  int diagonalexchange;  // switch 1=exchange is diagonal, 0=exchange is not diagonal
   void increase_nofcomponents(int n); // increase nofcomponents by n
   void add(jjjpar & b, Vector & abc); // add parameters b to this
   void addpars (int number, jjjpar & addjjj); // enlarge the set of parameters by
                                                        // inserting a new exchange parameters addjjj
							// into field at position number



   void save (FILE *file); // to save the parameters to a filehandle
   void saveatom (FILE *file); // to save the atom coordinates and properties to a filehandle
   void save_sipf(const char * path); //save single ion parameter file filename to path*


   jjjpar (FILE * fin, int nofcomp); //konstruktor with filehandle of mcphas.j file
   jjjpar (double x, double y, double z,char * sipffile); // constructor with filename of single ion parameter file
               // constructor with positions scattering length dwf
   jjjpar(double x,double y,double z, double slr,double sli, double dwf);
   jjjpar (int n=1,int diag=0,int nofmom=3); // konstructor without file
   jjjpar (const jjjpar & jjjpars);	// kopier-konstruktor
   
  ~jjjpar ();		//destruktor
  



// BASIC SIPF MODULE FUNCTIONS    *************************************************

  // integer to tell which module is loaded 0 - external, 1 - kramer, 2- cfield, 3 - brillouin
  int module_type;
  Matrix cnst;// cnst is the Zlm constants - put them into the matrix
   int nof_electrons; // no of electrons in d or f shell
private:
  Vector ABC;   // storage for single ion module paramters
  void getpolar(double x,double y, double z, double & r, double & th, double & ph);// calculates polar coordinates from Vector X(1..3)
  void get_parameters_from_sipfile(char * sipffilename); // function to read single ion parameter files


public:


   // subroutine to calculate momentum <J> from effective field gjmbH [meV]
   void  mcalc (Vector &mom, double & T, Vector &  gjmbH, double & lnZ,double & U,ComplexMatrix & ests);

   // returns transition element matrix M  and transition energy delta (to calculate chi0 in mcdisp,see manual)
   int  dmcalc (double & T,Vector &  gjmbheff, ComplexMatrix & mat,float & delta,ComplexMatrix & ests);
   int transitionnumber; // the transition associated with the ion (important if there are more in the single ion spectrum)

   ComplexMatrix est; // eigenstates
   ComplexMatrix mcalc_parstorage; // paramter storage for mcalc
   // returns eigenvalues and eigenstates matrix parameters of ion (if possible)
   ComplexMatrix & eigenstates (Vector & gjmbheff, double & T);
   // initialisis parameter storage for mcalc parameters (if possible)
   ComplexMatrix & mcalc_parameter_storage_init (Vector & gjmbheff,double & T);
   // returns operator matrices (n=0 Hamiltonian, n=1,...,nofcomponents: operators of moment components)
   Matrix opmat(int n,Vector & gjmbH);

private:
  // external module functions, intern_mcalc=0
 
  void (*m)(Vector*,double*,Vector*,double*,Vector*,char**,double*,double*,ComplexMatrix*);
  int  (*dm)(int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*,float*,ComplexMatrix*);

  void (*estates)(ComplexMatrix*,Vector*,double*,double*,Vector*,char**);
  void (*mcalc_parameter_storage)(ComplexMatrix*,Vector*,double*,double*,Vector*,char**);


public:
// OBSERVABLES *******************************************************
//1 . NEUTRON SCATTERING OPERATOR  --------------------------------------
   // calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
   // according to stored eigenstate matrix est, requires a call to eigenstates first
   ComplexVector & MQ(Vector & Qvec);
   ComplexVector Mq;

  // returns transition element matrix N(Q) in order to be able to go beyond
   // dipolar approximation in mcdisp - it requires a call to eigenstates first
   int dncalc(Vector & Qvec, double & T, ComplexMatrix & nat, ComplexMatrix & ests);

private :
  void (*mq)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*);
  int  (*ddnn)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexMatrix*);


public:
  double SLR,SLI; // scattering length
  double DWF; // DebeyWallerFactor [A^2]
  Vector magFFj0; // magnetic formfactor numbers
  Vector magFFj2; // magnetic formfactor numbers
  Vector magFFj4; // magnetic formfactor numbers
  Vector magFFj6; // magnetic formfactor numbers
  Vector Zc;      // Z-factors from Lovesey table 11.1 for Z(K) calc (needed to go beyond dipole approx)
//  D = 2 * pi / Q
//  s = 1 / 2 / D: sintheta = lambda * s
// 'magnetic formfactors
//  j0 = ff(1) * EXP(-ff(2) * s * s) + ff(3) * EXP(-ff(4) * s * s)
//  j0 = j0 + ff(5) * EXP(-ff(6) * s * s) + ff(7)
//  j2 = ff(8) * s * s * EXP(-ff(9) * s * s) + ff(10) * s * s * EXP(-ff(11) * s * s)
//  j2 = j2 + ff(12) * s * s * EXP(-ff(13) * s * s) + s * s * ff(14)
//  F = (j0 + j2 * (2 / gJ - 1))  formfactor F(Q)
//  RETURN TOTAL FORMFACTOR,
//    however if gJ=0 and Q>0 return spin form factor FS(Q)=<j0(Q)>
//            if gJ=0 and Q<0 return angular  form factor FL(Q)=<j0(Q)>+<j2(Q)>
   double F(double Q);

//   debyewallerfactor = EXP(-2 * DWF *s*s)
   double debyewallerfactor(double & Q);

// 2. charge density ----------------------------------------------------------

   Vector Np,Xip,Cp; // radial wave function parameters
   // evaluate radial wave function // r given in Angstroems, returns R(r) in units of 1/A^1.5
   double radial_wavefunction(double r);
   void save_radial_wavefunction(const char * filename);

   //functions to calculate radial matrix elements <r^n> from radial wave function
   int r2_from_radial_wavefunction();
   int r4_from_radial_wavefunction();
   int r6_from_radial_wavefunction();

   double r2;
   double r4;  // radial wave function exp values
   double r6;

   // calculation of chargedensity
   double rocalc (double & teta,double & fi,double & R, Vector & moments,double & T, Vector &  gjmbH);
   void (*ro_calc)(double*,double*,double*,double*,Vector*,double*,Vector*,double*,Vector*,char**);

// 3. moment density ----------------------------------------------------------

/****************************************************************************/
// function to calculate coefficients of expansion of spindensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
/****************************************************************************/
void spindensity_mcalc (Vector &mom,int xyz, double & T, Vector &  gjmbH, ComplexMatrix & parstorage);

/****************************************************************************/
// function to calculate coefficients of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
/****************************************************************************/
void orbmomdensity_mcalc (Vector &mom,int xyz, double & T, Vector &  gjmbH, ComplexMatrix & parstorage);

//***********************************************************************
// sub for calculation of spin density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
double spindensity_calc (double & teta,double & fi,double & R, Vector & moments);
Vector spindensity_calc (double & teta,double & fi,double & R, Vector & momentsx, Vector & momentsy, Vector & momentsz);

   double Fr(double r); // evaluate F(r)=1/r integral_r^inf dx R^2(x)
                        // r in units of Angstroems, F(r) in units of 1/A^3

//***********************************************************************
// sub for calculation of orbital moment density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
double orbmomdensity_calc (double & teta,double & fi,double & R, Vector & moments);
Vector orbmomdensity_calc (double & teta,double & fi,double & R, Vector & momentsx, Vector & momentsy, Vector & momentsz);
Vector currdensity_calc (double & teta,double & fi,double & R, Vector & momentlx, Vector & momently, Vector & momentlz);

//***********************************************************************
// subs for calculation gradient of spin and orbital moment density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
 Matrix gradspindensity_calc(double & teta,double & fi,double & R, Vector & momentsx, Vector & momentsy, Vector & momentsz);
 Matrix gradorbmomdensity_calc(double & teta,double & fi,double & R, Vector & momentlx, Vector & momently, Vector & momentlz);
 Matrix gradcurrdensity_calc(double & teta,double & fi,double & R, Vector & momentlx, Vector & momently, Vector & momentlz);

private:
#ifdef __linux__
  void *handle;
#else
//  HANDLE handle;
  HINSTANCE__* handle;
#endif
  
  void  (*sd_m)(Vector*,int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*);
  void  (*od_m)(Vector*,int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*);

  double rk_from_radial_wavefunction(int k); // needed for public radial wave function <r^n> calculation
   void set_zlm_constants();
   // sum over different Zlm using the coefficients a(l,m)
   double zlmsum(Matrix & a, double & teta, double & fi);


  // kramers internal module functions, module_type=1
  void kramer (Vector &mom,double & T,Vector & H, double & Z,double & U);
  int  kramerdm (int & tn,double & T,Vector &  heff, ComplexMatrix & mat,float & delta);
  Matrix krameropmat (int & n ,Vector & H);

  // realisation of class iops - cfield internal module functions, intern_mcalc=2
  // the class iops calls for some functionality the program cfield (e.g. for
  // getting stevens factors and other parameters, for the matrices Olm etc.)
  ionpars * iops;

  // brillouin internal module functions,module_type=3
  void brillouin (Vector &mom, double & T,Vector & H, double & Z,double & U);
  int  brillouindm (int & tn,double & T,Vector &  heff, ComplexMatrix & mat,float & delta);

  // cluster internal module functions, module_type=5
  void cluster_mcalc (Vector &mom,double & T,Vector & H, double & Z,double & U);
  int  cluster_dm (int & tn,double & T,Vector &  heff, ComplexMatrix & mat,float & delta);
  par * clusterpars;


};

#include<par.hpp>
#endif
