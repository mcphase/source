#ifndef IONPARS
#define IONPARS

#include <vector.h>
#include <cstdio>
#include <mpspecfunp.h>
#include <mcdisp.h>

#define IONPARS_MAXNOFCOMPONENTS 54
// standard operator sequence I1,....,I54
//                1  2  3  4    5    6   7   8   9    10        48  49   50   51   52   53   54
// module so1ion: Jx Jy Jz O22S O21S O20 O21 O22 O33S O32S .... O66 Jx^2 Jy^2 Jz^2 Jx^4 Jy^4 Jz^4


#define NOF_OLM_MATRICES 51
#define NOF_RIXS_MATRICES 9

#define SMALL_PROBABILITY 1e-6
#define SMALL_DEVIATION 1e-6
// ionpars: class for so1ion / cfield modules to 
//          - read single ion parameterfile for so1ion /cfield module
//          - load and store matrices for internal module so1ion /cfield
//          - diagonalize singleion problem and calculate obsevables and transition matrix elements

class ionpars  
{private:  
   // calculates scattering operator 
   void MQM(ComplexMatrix & MQXM,ComplexMatrix & MQYM,ComplexMatrix & MQZM, double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc);
   void setup_and_solve_Hamiltonian(Vector &  gjmbHxc,Vector & Hext,Vector & En,Matrix & zr,Matrix & zi,int sort);
   void calculate_Z_wn(Vector & En,double & T,double & Z,Vector & wn);
   void calculate_Z_wn(Vector & En,double & T,double & Zs,double & lnZs,Vector & wn);
   int noft(ComplexMatrix & est,double & T,double & pinit,double & ninit);
   void getijdelta_from_transitionnumber(int & i,int & j,float & delta,int & dj,int & tn,int & pr,ComplexMatrix &ests);
   complex<double> observable1(int & i,int & j,float & delta,Matrix & zr,Matrix & zi,
                         double & T,ComplexMatrix&ests,int & pr,const char *optype,Matrix & O);
  void popnr_diff(ComplexVector& dMQ,int &i,int &j,ComplexMatrix & est,float & delta,double & T,int &pr,const char * n);
 public:
   int so1ion; // switch to identify the coordinate system orientation with respect to the axes abc
              // 0 ...  cfield xyz||cab
              // 1 .... so1ion xyz||abc
   char * iontype; // description string
   int nof_electrons; // nof 4f electrons
   double J;// momentum quantum number 
   double gJ; // Lande factor 
   double alpha;
   double beta;  // stevens factors
   double gamma;
   double r2;
   double r4;  // radial wave function exp values
   double r6;

   complex<double> sigma0,sigma1,sigma2; // optical conductivity for RIXS operator

   Matrix Ja; Matrix Jb; Matrix Jc; Matrix Hcf;
   Matrix cnst;

   Matrix **Olm; // array of matrices
   ComplexMatrix ** Ri; // for RIXS

   Matrix ** In;
   Vector Blm; // Cf parameters  
   Vector Llm; // Cf parameters  

   // functions needed to calculate thermal expectation value of observables 
   void cfeigenstates (ComplexMatrix *est, Vector &  gjmbHxc,Vector & Hext, double & T);
   void Icalc(Vector & JJ,double & T, Vector &  gjmbHxc,Vector & Hext, double & lnZs, double & U, ComplexMatrix & ests);

   Matrix opmat (int & n ,Vector &  Hxc,Vector & Hext);
   // on input
   // n		which operator 0=Hamiltonian, 1,2,3=J1,J2,J3
   // Hext,Hxc	vector of external and exchange field [meV]
   // on output    
   // operator matrix of Hamiltonian, I1, I2, I3 depending on n

   void Jcalc(Vector & JJ,double & T, Vector &  gjmbHxc,Vector & Hext, ComplexMatrix & ests);
   // charge density coefficients
   void chargedensity_coeffcalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage);
   // calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
   // according to stored eigenstate matrix est
   // calculates the scattering operator given the polar angles th, ph (with respect to the CEF coordinate 
   // system xyz and the <jl(qr)> and the eigenstate matrix with eigenstates and thermal population numbers
   ComplexVector & MQ(double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc,ComplexMatrix & est);

    // and transition matrix elements
   int  du1calc (int & tn,double & T,Vector &  gjmbHxc,Vector & Hext, ComplexVector & u1,float & delta,int & n, int & nd,ComplexMatrix & ests);
   int  dJ1calc (int & tn,double & T,Vector &  gjmbHxc,Vector & Hext, ComplexVector & u1,float & delta,ComplexMatrix & ests);
   int dchargedensity_coeff1calc(int & tn,double & T,Vector &  gjmbHxc,Vector & Hext, ComplexVector & chargedensity_coeff1,float & delta,ComplexMatrix & ests);   
   int dMQ1(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexVector & v1);
   // for RIXS
   int cfielddrixs1(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexVector & drixs);
   
   void savBlm(FILE * file); // saving Blm to file 
   void savLlm(FILE * file); // saving Blm to file 
   void save(FILE * file); // save ion parameters to file 

   ionpars(int dimj);
   ionpars(FILE * cf_file, char * cffilename);
   ionpars (char * iontype); // constructor from iontype (mind:no matrices filled with values !)
   ~ionpars();
   ionpars(const ionpars & p);
};
#endif

