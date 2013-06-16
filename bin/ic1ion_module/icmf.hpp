/* icpars.hpp
 * 
 * Header file for the parameters class icpars and cfpars
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#ifndef ICMF_H

// --------------------------------------------------------------------------------------------------------------- //
// Defines a class to hold the eigenvalues and eigenvectors
// --------------------------------------------------------------------------------------------------------------- //
class iceig
{
   private:
      int _Hsz;
      double *_E;
      double *_V;
      complexdouble *_zV;

   public:
      iceig() { _Hsz=0; _E=0; _V=0; _zV=0; }       // Blank constructor
      ~iceig();                                    // Destructor
      iceig(int Hsz, bool isreal=true);            // Constructs a zero eigenstates object, allocating memory
      iceig(sMat<double>&H);                       // Constructs the eigenstates of H
      iceig(sMat<double>&H, sMat<double>&iH);      // Constructs the eigenstates of H+iH
      iceig(int Hsz, double *E, double *V);        // Constructor for known eigenvalues/vectors
      iceig(int Hsz, double *E, complexdouble *V); // Constructor for known eigenvalues/complex eigenvectors
      iceig(int Hsz,double*E,complexdouble*V,int s);// Constructor for known eigenvalues/complex eigenvectors
      iceig(const iceig &p);                       // Copy constructor
      iceig &operator = (const iceig &p);          // Copy assignment - overwrites previous eigenvalues/vectors
      double E(int i) { return _E[i]; }
      double *E() { return &_E[0]; }
      double *V(int i) { return &_V[i*_Hsz]; }
      double V(int i,int j) {return _V[j*_Hsz+i];}
      complexdouble *zV(int i) { return &_zV[i*_Hsz]; }
      complexdouble zV(int i,int j){return _zV[j*_Hsz+i];}
      void assign(int Hsz,double*E,complexdouble*V,int s);
      void calc(sMat<double>&H);                   // Calculates the eigenstates of H
      void calc(sMat<double>&H, sMat<double>&iH);  // Calculates the eigenstates for a matrix H+iH
      void calc(int Hsz, complexdouble *H);
      void calc(int Hsz, double *H);
      void lcalc(icpars &pars, sMat<double> &H);   // Calculates a partial set of eigenstates of H
      void lcalc(icpars &pars, sMat<double> &H, sMat<double> &iH); 
      void lcalc(icpars &pars, complexdouble *H);  // Calculates a partial set from a fortran-style matrix
      #ifndef NO_ARPACK
      void acalc(icpars &pars, complexdouble *H);  // Calculates a partial set using the Arnoldi method (ARPACK)
      void acalc(icpars &pars, sMat<double> &J);
      void acalc(icpars &pars, sMat<double> &J, sMat<double> &iJ);
      #endif
      bool iscomplex() { return _zV!=0; }          // Determines of eigenvalues are complex
      int Hsz() { return _Hsz; }
      std::string strout();                        // Prints the eigenstates out as a matrix
};

// --------------------------------------------------------------------------------------------------------------- //
// Defines a class to hold the matrices Lx/Sx etc. and to calculate their expectation values
// --------------------------------------------------------------------------------------------------------------- //
class icmfmat
{
   private:
      int _n;
      orbital _l;
      int _num_op;
      std::string _density;                        // Flag to output expectation values of spin/orbital density operator.

   public:
      std::vector<sMat<double> > J;                // A vector of the matrices [Sx Lx Sy Ly Sz Lz]
      std::vector<int> iflag;                      // Vector to determine if matrix is imaginery

      icmfmat();                                   // Blank constructor
      icmfmat(int n, orbital l, int num_op,        // Constructor for l^n configuration
        bool save_matrices, std::string density="");
      void Jmat(sMat<double>&J, sMat<double>&iJ,   // Calculates the mean field matrix sum_i (H_i*J_i)
        std::vector<double>&gjmbH, bool save_matrices);
      std::vector<double> expJ(iceig&VE, double T, // Calculates the expectation values <V|J|V>exp(-beta*T)
        std::vector<std::vector<double> >&matel,   //   matel is an m*n matrix of the elements <n|Jm|n>
        bool save_matrices);
      std::vector<double> spindensity_expJ(iceig&VE,//Calculates the expectation values
        int xyz,double T,                          //   <V|spindensitycoeff_of_Zlm|V>exp(-beta*T)
        std::vector<std::vector<double> >&matel,   
        bool save_matrices);
      std::vector<double> orbmomdensity_expJ(      // Calculates the expectation values
        iceig&VE,int xyz, double T,                //   <V|orbmomdensitycoeff_of_Zlm|V>exp(-beta*T)
        std::vector<std::vector<double> >&matel,
        bool save_matrices);
      void u1(std::vector<double> &u1,             // Calculates the vector u1 = <i|Ja-<Ja>|j>
        std::vector<double>&iu1, iceig&V, double T,// * sqrt{exp(-beta_i*T)-exp(-beta_j*T)}
        int i, int j, int p, float &d,
        bool save_matrices);
      void dod_u1(int xyz, std::vector<double>&u1, // Calculates the vector u1 = <i|M(q)-<M(q)>|j>
        std::vector<double>&iu1, iceig&V, double T,// * sqrt{exp(-beta_i*T)-exp(-beta_j*T)}
        int i, int j, int p, float &d, 
        bool save_matrices);
      #ifdef JIJCONV
      std::vector<double> jijconv;                 // Conversion from Stevens/Wybourne norm of Jij pars
      #endif
};

#endif
