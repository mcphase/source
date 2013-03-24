#ifndef IONPARS
#define IONPARS

#include <vector.h>
#include <cstdio>
#include <mpspecfunp.h>

#define IONPARS_MAXNOFCOMPONENTS 51
// standard operator sequence I1,....,I51
// module so1ion: Jx Jy Jz O22S O21S O20 O21 O22 O33S O32S .... O66 Jx^2 Jy^2 Jz^2

// ionpars: class for so1ion / cfield modules to 
//          - read single ion parameterfile for so1ion /cfield module
//          - load and store matrices for internal module so1ion /cfield
//          - diagonalize singleion problem and calculate obsevables and transition matrix elements

class ionpars  
{private:  
   // calculates scattering operator 
   void MQM(ComplexMatrix & MQXM,ComplexMatrix & MQYM,ComplexMatrix & MQZM, double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc);
  
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
   ComplexMatrix Jaa;
   ComplexMatrix Jbb;
   ComplexMatrix Jcc;
   Matrix cnst;

   Matrix **Olm; ComplexMatrix **OOlm; // array of matrices
   ComplexMatrix ** Ri; // for RIXS

   ComplexMatrix ** In;
   Vector Blm; // Cf parameters  
   Vector Llm; // Cf parameters  

   // functions needed to calculate thermal expectation value of moment  
   Vector & cfield (double & T,Vector &  gjmbHxc,Vector & Hext, double & Z,double & U, ComplexMatrix & ests);
   void cfieldJJ(Vector & JJ,double & T, Vector &  gjmbHxc,Vector & Hext, double & lnZs, double & U, ComplexMatrix & ests);
   void cfeigenstates (ComplexMatrix *est, Vector &  gjmbHxc,Vector & Hext, double & T);
   // and transition matrix elements
   int  cfielddm (int & tn,double & T,Vector &  gjmbHxc,Vector & Hext, ComplexVector & u1,float & delta,ComplexMatrix & ests);
   int cfielddn(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexVector & v1);
   // for RIXS
   int cfielddrixs1(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexVector & drixs);
   // charge density coefficients
   void chargedensity_coeffcalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage);
   int dchargedensity_coeff1calc(int & tn,double & T,Vector &  gjmbHxc,Vector & Hext, ComplexVector & chargedensity_coeff1,float & delta,ComplexMatrix & ests);
   // calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
   // according to stored eigenstate matrix est
   // calculates the scattering operator given the polar angles th, ph (with respect to the CEF coordinate 
   // system xyz and the <jl(qr)> and the eigenstate matrix with eigenstates and thermal population numbers
   ComplexVector & MQ(double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc,ComplexMatrix & est);


   void savBlm(FILE * file); // saving Blm to file 
   void savLlm(FILE * file); // saving Blm to file 
   void save(FILE * file); // save ion parameters to file 


   ionpars(int dimj);
   ionpars(FILE * cf_file);
   ionpars (char * iontype); // constructor from iontype (mind:no matrices filled with values !)
   ~ionpars();
   ionpars(const ionpars & p);
};
#endif

