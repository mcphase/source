// jjjpar is a class to store all parameters associated
// with a specific ion, for example CEF parameters
// and exchange parameters 
// it is used by many programs in the package
// moreover, it loads also the user defined single ion module functions (linux only)
#ifndef JJJPAR
#define JJJPAR

typedef void fnc_t();

#ifdef __MINGW32__
#include <windows.h>
#else
#include<dlfcn.h>
#endif

#include<martin.h>
#include<ionpars.hpp>
#include<perlparse.h>
#include<myev.h>
#include<stdlib.h>

#include "sparsecomplex.hpp"  // For sparse matrices class for cluster module

#define MAXSAVEQ 5   // Number of Q vector values to save in calculation of F(Q)
                     //  so as to not repeat calculations.

#define MAGMOM_EV_DIM 3  // eigenvector dimension for moment oscillation
#define SPIN_EV_DIM 3  // eigenvector dimension for spin oscillation
#define ORBMOM_EV_DIM 3  // eigenvector dimension for orbmom oscillation
#define CHARGEDENS_EV_DIM 28  // eigenvector dimension for chargedensity oscillation
#define SPINDENS_EV_DIM 49  // eigenvector dimension for spindensity oscillation
#define ORBMOMDENS_EV_DIM 49  // eigenvector dimension for orbmomdensity oscillation
#define PHONON_EV_DIM 3 // phonon eigenvector dimension ?? check if this is sensible ??

// Class to store work matrices for iterative eigensolvers (e.g. ARPACK, FEAST) - since these routines are called 
// often each MF iterations and deleting/allocating them many times seems to give memory errors.
class iterwork {
  public:
     complexdouble *zwork;
     double *dwork;
     int *iwork, zsize, dsize, isize;
     iterwork(): zsize(0), dsize(0), isize(0) {}; 
     iterwork(int lzwork, int ldwork, int liwork);
    ~iterwork(); 
     void realloc_z(int lzwork);
     void realloc_d(int ldwork);
     void realloc_i(int liwork);
};

class par;
class jjjpar
{
public:
// ********************************************************************************
//                               basic parameters
// ********************************************************************************
  char * sipffilename; // single ion parameter filename
  char * modulefilename; // module name
  char * clusterfilename;
    double J(); // returns total angular momentum if possible
   Vector &  tetan(); //returns stevens parameters if possible
  Matrix * G; // magnetoelastic coupling constants for the ion (in meV)
  Vector xyz,mom; // atom position, moment
  int paranz;   // number of exchange parameters
  int nofcomponents; // number of moments (components of moment vector)
  double gJ,ninit,pinit,maxE;
  Matrix *jij; // exchange constants 
  Vector *dn; // exchange - coordinates of neighbors da db dc
  Vector *dr; // exchange - coordinates of neighbours in Angstrom in Euclidean ijk system
  int *sublattice; // sublattice of neighbours
  int diagonalexchange;  // switch 1=exchange is diagonal, 0=exchange is not diagonal
   void increase_nofcomponents(int n); // increase nofcomponents by n
   void decrease_nofcomponents(int n); // decrease nofcomponents by n
   void add(jjjpar & b, Vector & abc); // add parameters b to this
   void addpars (int number, jjjpar & addjjj); // enlarge the set of parameters by
                                                        // inserting a new exchange parameters addjjj
							// into field at position number
   void scalepars (double scalefactor); // multiply all exchange parameters with scale factor
   void delpar (int number); // remove a neighbour from list

   void save (FILE *file,int noindexchange); // to save the parameters to a filehandle
   void saveatom (FILE *file); // to save the atom coordinates and properties to a filehandle
   void saveG(FILE * file); // save coupling constants G
   void save_sipf(const char * path); //save single ion parameter file filename to path*
   void save_sipf(FILE *file); //save single ion parameter file filename to path*


   jjjpar (FILE * fin, int nofcomp,int verbose=0); //constructor with filehandle of mcphas.j file
   jjjpar (double x, double y, double z,char * sipffile,int n,int verbose=0); // constructor with filename of single ion parameter file
               // constructor with positions scattering length dwf
   jjjpar(double x,double y,double z, double slr,double sli, double dwf);
   jjjpar (int n=1,int diag=0,int nofmom=3); // constructor without file
   jjjpar (const jjjpar & jjjpars);	// copy-constructor
   
  ~jjjpar ();		//destruktor
  



// ********************************************************************************
//                          BASIC SIPF MODULE FUNCTIONS    
// ********************************************************************************

  // integer to tell which module is loaded 0 - external, 1 - kramer, 2- cfield, 3 - brillouin
  int module_type;
  Matrix cnst;// cnst is the Zlm constants - put them into the matrix
   int nof_electrons; // no of electrons in d or f shell
private:
  Vector ABC;   // storage for single ion module paramters
  void getpolar(double x,double y, double z, double & r, double & th, double & ph);// calculates polar coordinates from Vector X(1..3)
  void get_parameters_from_sipfile(char * sipffilename,int verbose); // function to read single ion parameter files
  int  get_exchange_indices(char *instr, Matrix *exchangeindices,const char * ie);

public:
   // subroutine to calculate expectation values <Ialpha> alpha=1...nofcomponents
   // from exchange field Hxc [meV] and external field Hext
   void  Icalc (Vector &mom, double & T, Vector &  Hxc,Vector & Hext, double & lnZ,double & U,ComplexMatrix & parstorage);
   void  Icalc (Matrix &mom, Vector & T, Vector &  Hxc,Vector & Hext, Vector & lnZ,Vector & U,ComplexMatrix & parstorage);

   // returns transition element matrix M  and transition energy delta (to calculate chi0 in mcdisp,see manual)
   int  du1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & u1,float & delta,int & n, int & nd, ComplexMatrix & ests);
   int transitionnumber; // the transition associated with the ion (important if there are more in the single ion spectrum)

   // calculates series of single ion susceptibility matrices for different energies 
   // output:returns 0 on success
   //        the Matrices chi0pointer[1....nofstps] must exist and will be filled with values
   //        ...... the contribution of transition transitionnumber is added to these matrices
   // input: emin est nofstps define energies, eps is the imaginary part of the energy
   //        Q       the Q vector in 1/A
   //        qcounter is a counter telling which q vector in the list is calculated
   //                  sign(qcounter) <0 indicates that chi0c matrices should be cleared and nothing calculated
   //        delta ... sign determines if energy gain or loss term is added
   int chi0(ComplexMatrix ** chi0pointer,double & emin, double estp, int & nofstps, double & eps,Vector & Q, 
            int qcounter,float & delta, double & T,Vector &  Hxc,Vector & Hext, ComplexMatrix & ests,int i1,int j1,int k1,int l1);

   ComplexMatrix est; // eigenstates
   ComplexMatrix Icalc_parstorage; // paramter storage for Icalc
   // returns eigenvalues and eigenstates matrix parameters of ion (if possible)
   ComplexMatrix & eigenstates (Vector &  Hxc,Vector & Hext, double & T);
   void print_eigenstates(FILE *fout);
   // initialisis parameter storage for Icalc parameters (if possible)
   ComplexMatrix & Icalc_parameter_storage_init (Vector &  Hxc,Vector & Hext,double & T);
   // returns operator matrices (n=0 Hamiltonian, n=1,...,nofcomponents: operators of moment components)
   Matrix opmat(int n,Vector &  Hxc,Vector & Hext);

private:
  // external module functions, intern_Icalc=0
 
  void (*I)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,double*,double*,ComplexMatrix*);
  void (*IM)(Matrix*,Vector*,Vector*,Vector*,double*,Vector*,char**,Vector*,Vector*,ComplexMatrix*);
  int  (*du)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,int*,int*,ComplexMatrix*);

  void (*estates)(ComplexMatrix*,Vector*,Vector*,double*,double*,Vector*,char**);
  void (*Icalc_parameter_storage)(ComplexMatrix*,Vector*,Vector*,double*,double*,Vector*,char**);

  int (*dyn_opmat)(int*,char**,Vector*,Vector*,Matrix*);
  Matrix *opmatM[52];

public:
// ********************************************************************************
//                                       OBSERVABLES 
// ********************************************************************************
//0. PHONON displacement
int pcalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext,ComplexMatrix & ests);
int  dP1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & dP1,ComplexMatrix & ests);
private:
void (*p)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
int  (*dP1)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*);

public:
//1. MAGNETIC MOMENT
   // returns magnetic moment
   int mcalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext,ComplexMatrix & parstorage);
   int mcalc(Matrix &mom, Vector & T, Vector &  Hxc,Vector & Hext,ComplexMatrix & parstorage);
   int micalc(Vector &momi,  double & T, Vector &  Hxc,Vector & Hext,ComplexMatrix & parstorage);
   int Lcalc(Vector &L, double & T, Vector &  Hxc,Vector & Hext,ComplexMatrix & parstorage);
   int Lcalc(Matrix &L, Vector & T, Vector &  Hxc,Vector & Hext,ComplexMatrix & parstorage);
   int Scalc(Vector &S, double & T, Vector &  Hxc,Vector & Hext,ComplexMatrix & parstorage);
   int Scalc(Matrix &S, Vector & T, Vector &  Hxc,Vector & Hext,ComplexMatrix & parstorage);
   int  dm1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & dm1,ComplexMatrix & ests);
   int  dmi1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & dmi1,ComplexMatrix & ests);
   int  dL1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & dL1,ComplexMatrix & ests);
   int  dS1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & dS1,ComplexMatrix & ests);

private:  // handle for mcalc in loadable modules
  void (*m)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
  void (*mM)(Matrix*,Vector*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
  void (*L)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
  void (*LM)(Matrix*,Vector*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
  void (*S)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
  void (*SM)(Matrix*,Vector*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
  int  (*dm1)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*);
  int  (*dL1)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*);
  int  (*dS1)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*);

//2 . NEUTRON SCATTERING OPERATOR  --------------------------------------
   // calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
   // according to stored eigenstate matrix est, requires a call to eigenstates first
public:
  int MQ(ComplexVector & Mq,Vector & Qvec);
  
  // returns transition element matrix N(Q) in order to be able to go beyond
   // dipolar approximation in mcdisp - it requires a call to eigenstates first
   int dMQ1calc(Vector & Qvec, double & T, ComplexVector & dMQ1,float & delta, ComplexMatrix & ests);

private :
  void (*mq)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*);
  int  (*ddnn)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexVector*,double*);

public:
  double SLR,SLI; // scattering length
  double DWF; // DebeyWallerFactor [A^2]
  int FF_type; // use by program mcdisp, mcdiff to store which formfactor this ion is
void FFinfo(FILE * fout); // formfactor information print to fout, for mcdiff and mcdisp
                         // info about formactor is printed according to settings of FFtype
                         // FF_type has to be set by mcdiff / mcdisp correctly before calling
                         // this function
  int checkFFcoeffnonzero(int l);
  void magFFout(const char * linestart ,FILE * fout);
private:
  Vector magFFj0; // magnetic formfactor numbers
  Vector magFFj2; // magnetic formfactor numbers
  Vector magFFj4; // magnetic formfactor numbers
  Vector magFFj6; // magnetic formfactor numbers
public:
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
   double j0(double Q);
   double j1(double Q);
   double j2(double Q);
   double j3(double Q);
   double j4(double Q);
   double j5(double Q);
   double j6(double Q);
  int jl_lmax; // initialized to 6,will be lowered if Np+Np<l+2 at 
               //calculation of jjjpar::F(Q) from radial wave function, used for printout mcdiff.out mdisp*.*

private:
   double jl(int l,double Q);
   long double tl(int l,int N,long double x);
   long double sn(int n,int N,long double x);
   long double cn(int n,int N,long double x);
   double Fsaved[MAXSAVEQ+1],Qsaved[MAXSAVEQ+1]; int nsaved;
   double DBWsaved[MAXSAVEQ+1],DBWQsaved[MAXSAVEQ+1]; int DBWnsaved;

public:
//   debyewallerfactor = EXP(-2 * DWF *s*s)
   double debyewallerfactor(double & Q);

//  DENSITIES ----------------------------------------------------------
// ... radial wave functions
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

   double charge; // charge in units of |e|
   int magnetic;

private:
  double rk_from_radial_wavefunction(int k); // needed for public radial wave function <r^n> calculation
   // sum over different Zlm using the coefficients a(l,m)
   double zlmsum(Matrix & a, double & teta, double & fi);

// 3. charge density ----------------------------------------------------------
public:
   // calculation of chargedensity
   // function to calculate coefficients of expansion of chargedensity in terms
   // of Zlm R^2(r) at a given temperature T and  effective field H
   int chargedensity_coeff (Vector &mom, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage);
   int dchargedensity_coeff1 (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & dchargedensity_coeff1,ComplexMatrix & ests);

   double chargedensity_calc (double & teta,double & fi,double & R, Vector & moments);
private:
   void (*ro_calc)(double*,double*,double*,double*,Vector*,double*,Vector*,char**);
   void  (*cd_m)(Vector*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
   int   (*cd_dm)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,float*,ComplexMatrix*);

public:
// 4. spindensities ----------------------------------------------------------
// function to calculate coefficients of expansion of spindensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
int spindensity_coeff (Vector &mom,int xyz, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage);
int dspindensity_coeff1 (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & dspindensity_coeff1,ComplexMatrix & ests);
// sub for calculation of spin density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
double spindensity_calc (double & teta,double & fi,double & R, Vector & moments);
Vector spindensity_calc (double & teta,double & fi,double & R, Vector & momentsx, Vector & momentsy, Vector & momentsz);
   double Fr(double r); // evaluate F(r)=1/r integral_r^inf dx R^2(x)
                        // r in units of Angstroems, F(r) in units of 1/A^3
// subs for calculation gradient of spin and orbital moment density given a radius R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
 Matrix gradspindensity_calc(double & teta,double & fi,double & R, Vector & momentsx, Vector & momentsy, Vector & momentsz);
private:
  void  (*sd_m)(Vector*,int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
   int   (*sd_dm)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,int*,float*,ComplexMatrix*);

// 5. orbmomdensities ----------------------------------------------------------
public:
// function to calculate coefficients of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
int orbmomdensity_coeff (Vector &mom,int xyz, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage);
int dorbmomdensity_coeff1 (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & dorbmomdensity_coeff1,ComplexMatrix & ests);
// sub for calculation of orbital moment density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
double orbmomdensity_calc (double & teta,double & fi,double & R, Vector & moments);
Vector orbmomdensity_calc (double & teta,double & fi,double & R, Vector & momentsx, Vector & momentsy, Vector & momentsz);
// subs for calculation gradient of spin and orbital moment density given a radius R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
 Matrix gradorbmomdensity_calc(double & teta,double & fi,double & R, Vector & momentlx, Vector & momently, Vector & momentlz);
private:
  void  (*od_m)(Vector*,int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexMatrix*);
  int   (*od_dm)(int*,double*,Vector*,Vector*,double*,Vector*,char**,ComplexVector*,int*,float*,ComplexMatrix*);

// 6. currentdensities ----------------------------------------------------------
public:
Vector currdensity_calc (double & teta,double & fi,double & R, Vector & momentlx, Vector & momently, Vector & momentlz);
// subs for calculation gradient of spin and orbital moment density given a radius R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
Matrix gradcurrdensity_calc(double & teta,double & fi,double & R, Vector & momentlx, Vector & momently, Vector & momentlz);

//7 . Resonant X-ray Scattering  --------------------------------------
   // calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
   // according to stored eigenstate matrix est, requires a call to eigenstates first
public:
  //int RMXS(ComplexVector & Mq,Vector & Qvec); // resonant magnetic xray scattering - maybe to be needed in mcdiff ?
  
  // returns transition elements of RIXS operator - it requires a call to eigenstates first
   int drixs1calc(Vector & Qvec, double & T, ComplexVector & drixs1, ComplexMatrix & ests);

private :
  //void (*rmx)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*);
  int  (*rixs)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexVector*,double*);
            //(transitionnumber,th,ph,J0,J2,J4,J6,ests,T,drixs,maxE)
private:
//#ifdef __linux__
//  void *handle;
//#else
////  HANDLE handle;
//  HINSTANCE__* handle;
//#endif
#ifdef __MINGW32__
  HINSTANCE__* handle;
#else
void *handle;
#endif
  


// ********************************************************************************
//                                INTERNAL MODULE FUNCTIONS 
// ********************************************************************************

  // kramers internal module functions, module_type=1
  void kramer (Vector &mom,double & T,Vector &  Hxc,Vector & Hext, double & Z,double & U);
  int  kramerdm (int & tn,double & T,Vector &  Hxc,Vector & Hext, ComplexVector & u1,float & delta,int & n, int & nd);
  Matrix krameropmat (int & n ,Vector &  Hxc,Vector & Hext);

  // realisation of class iops - cfield internal module functions, intern_Icalc=2
  // the class iops calls for some functionality the program cfield (e.g. for
  // getting stevens factors and other parameters, for the matrices Olm etc.)
public:
  ionpars * iops;
  par * clusterpars;
private:
  // brillouin internal module functions,module_type=3
  void brillouin (Vector &mom, double & T,Vector &  Hxc,Vector & Hext, double & Z,double & U);
  int  brillouindm (int & tn,double & T,Vector &  Hxc,Vector & Hext, ComplexVector & u1,float & delta,int & n, int & nd);

  // cluster internal module functions, module_type=5
  void cluster_Icalc_mcalc_Micalc (int code,Vector &mom,double & T,Vector &  Hxc,Vector & Hext, double & Z,double & U);
  void cluster_Micalc (Vector &mom,ComplexMatrix & ests);
  int  cluster_dm (int code,int & tn,double & T, ComplexVector & u1,float & delta,int & n, int & nd,ComplexMatrix & ests);
  void cluster_est(ComplexMatrix * est,Vector &Hxc,Vector &Hext,double & T);
  void cluster_calcH_and_diagonalize(Vector & En,ComplexMatrix &zc,Vector & Hxc,Vector & Hext);
  void cluster_ini_Imat();
  void cluster_Iaa(zsMat<double> *Iai, int a, int i);
  zsMat<double> ** Ia; zsMat<double> ** cluster_M; 
  zsMat<double> *clusterH; Vector *oldHext; bool justinit; // Added to cache Hamiltonian matrices between iterations when Hext=same
  int * dnn; int dim; bool useperl;
  iterwork *workspace;
  double truncate, feast, arpack; int fdim; bool is1sttrunc, oldeig; complexdouble *zm;

};

#include<par.hpp>
#endif
