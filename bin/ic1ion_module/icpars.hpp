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

#ifndef ICPARS_H
#define ICPARS_H
//#define JIJCONV

#define PHYSPROP_MAGBIT 1                  // Defined bit values for the flags.calcphys bit mask
#define PHYSPROP_SUSBIT 2
#define PHYSPROP_INVBIT 4
#define PHYSPROP_CP_BIT 8

#define MEV2CM 8.06554486                  // Conversion factor: 1meV = MEV2CM * 1cm^{-1}    [ == (Q_e/1000)/(hc*100) ] 
#define MEV2K  11.6045047                  // Conversion factor: 1meV = MEV2K * 1K           [ == (Q_e/1000)/k_B      ] 
#define CM2K   1.43877505                  // Conversion factor: 1cm^{-1} = CM2K * 1K        [ == (hc*100)/k_B        ]

#define F424    0.13805                    // Ratios of F_4/F_2 slater integrals for 4f hydrogenic wavefunctions
#define F624    0.015018                   // Ratios of F_6/F_2 slater integrals for 4f hydrogenic wavefunctions
#define F425    0.14218                    // Ratios of F_4/F_2 slater integrals for 5f hydrogenic wavefunctions
#define F625    0.016104                   // Ratios of F_6/F_2 slater integrals for 5f hydrogenic wavefunctions

#define KBc     0.69503568043              // cm^{-1)/K - Boltzmann constant
#define MUBc    0.46686437                 // cm^{-1} / Tesla - Bohr magneton
#define KB      0.08617343183              // meV/K - Boltzmann constant
#define MUB     0.0578838263               // meV/T - Bohr magneton
#define NAMUB   5.5849397                  // N_A*mu_B - J/T/mol - product of Bohr magneton and Avogadro's number

#define GS      2.0023193043622            // The electron gyromagnetic ratio

// --------------------------------------------------------------------------------------------------------------- //
// Defines a class to hold the values of the crystal field parameters and methods to convert and print them.
//   Note that the parameters are stored internally in both Wybourne normalisation (type Lkq) and in whichevery
//   format they were specified in the input file mcphas.ic with automatic conversion if necessary.
// --------------------------------------------------------------------------------------------------------------- //
class cfpars
{
   private:
      double _Bi[27];                        // The crystal field parameters for internal calculations (type Lkq).
      double _Bo[27];                        // The crystal field paramaters for external output
      bool _minustype;                       // How -|q| parameters should be printed: true==Bk-q; false==BkqS
      bool _withiontype;                     // Flag to show if constructed with ion name (has stevfact and <r^k>)
      std::vector<double> _stevfact;         // The Stevens factors \theta_k for this ion (if any).
      std::vector<double> _istevfact;        // The inverse Stevens factors \theta_k for this ion (if any).
      std::vector<double> _rk;               // The radial integrals <r^k> looked up from tables
      std::string _cfname;                   // Name of CF parameter: [Akq,Vkq],[Bkq,Wkq],[Lkq,Dkq],ARkq
      std::string _normalisation;            // Normalisation for crystal field parameters (accepted: Stevens, Wybourne)
      std::string _units;                    // The energy units of the (external) parameters [internal in 1/cm]
   public:
      // Properties
      void assign(std::string &S, int &k, int &q, double val); // Assigns a value val to a parameter of type S
      std::string cfname() { return _cfname; }
      std::string norm() { return _normalisation; }
      std::string units() { return _units; }
      void calc_stevfact(int n, orbital l);  // (Re)calculates the stevens factor for l^n
      void find_rk(std::string &ionname);    // (Re)lookup the radial integral <r^k> for a particular ion
      double get(int k, int q);              // Returns the value of the external parameter k,q
      double alpha(){ return _stevfact[0]; }
      double beta() { return _stevfact[1]; }
      double gamma(){ return _stevfact[2]; }

      // Methods
      std::string cfparsout(const char*de);  // Converts nonzero CF pars to string with delimiter de 
      void conv_B_norm(std::string &nrm);    // Converts parameters between Stevens and Wybourne normalisations
      void conv_e_units(std::string &units); // Converts parameters to different units
      void conv(std::string &newcfname);     // Converts parameters to a new type (A,V,B,W,L,D,AR)
      void convback();                       // Converts internal parameters back from external.
      bool check();                          // Checks that internal and external parameters agree.

      // Overloaded operators
      double operator()(int k, int q) const; // Operator to access internal parameters Bi
      bool operator==(cfpars c) const;       // Operator to determine if parameters are the same
      bool operator!=(cfpars c) const;       // Operator to determine if parameters are not the same

      // Constructors
      cfpars();                              // Blank constructor
      cfpars(std::string&i,int n,orbital l); // Constructor for a particular ion, assigning a Stevens factor, and <r^k>.
};

// Defines a class to hold the values needed to for input to the module
// --------------------------------------------------------------------------------------------------------------- //
class icpars
{
   private:
      std::vector<double> _F;                // Internal coulomb parameters in cm^{-1} for calculations
      double _xi;                            // Internal spin-orbit parameters for calculations
      std::vector<double> _alpha;            // Internal CI parameters in cm^{-1} for calculations
      double _econv;                         // Energy conversion factor fgor e_units to cm^{-1}
      // Needs to access _F, _xi, _alpha
      friend sMat<double> ic_hmltn(sMat<double> &H_cfi, icpars &pars);
      friend sMat<double> icf_hmltn(sMat<double> &Hcfi, icpars &pars);
      friend sMat<double> ic_Hcso(icpars &pars);
      friend void getfromionname(std::string &ionname, icpars &pars);
      friend void conv_e_units(icpars &pars, std::string &newunit);
      friend void ic_parseinput(const char *filename, icpars &pars);
    //friend int ic_peig(icpars &pars, double *Vd, complexdouble *zVd, double *eigval);
      #ifdef JIJCONV
      bool _jijconvalreadycalc;              // Flag to show that the conversion factor for Jij already calculated
      #endif
   public:
      orbital l;                             // Orbital angular momentum of electrons, defaults to f-electrons (l=3)
      int n;                                 // Number of electrons in lowest configuration.
      int calcphys;                          // Determines whether to output physical properties and which ones
      int mag_units;                         // Determines unit to output magnetisation: uB/atom, emu/mole, m^3/mol
      std::string e_units;                   // Energy units of parameters and output (accepted: meV, K, cm^{-1})
      std::string ionname;                   // Name of ion if the user specifies it
      std::vector<double> F;                 // The eletrostatic radial integral parameters
      double xi;                             // The spin-orbit radial integral parameter
      std::vector<double> alpha;             // The linear configuration interaction parameters
      cfpars B;                              // The crystal field parameters
      bool perturb;                          // Flag to indicate if perturbation routine should be used to calculate M
      bool partial;                          // Flag to indicate if partial diagonalisation should be used
      bool arnoldi;                          // Flag to indicate if the Arnoldi method should be use to diagonalise H
      bool partial_standalone;               // Flag to indicate if partial diag. should be used for ic1ion.out
      bool arnoldi_standalone;               // Flag to indicate if the Arnoldi method should be use for ic1ion.out
      bool save_matrices;                    // Flag to indicate if matrices for CF, etc. should be saved and reloaded
    //bool bflag;                            // Flag to show if a field norm= or nostevfact is given
      int spectrelevels;                     // If using the spectre method, number of LS levels to keep. -1 means all
      double truncate_level;                 // Fraction of matrix to keep, for matrix truncation.
      int num_eigv;                          // Number of eigenvectors to print in output
      std::string density;                   // Flag to output expectation values of spin/orbital density operator.
    //std::string observable;                // Flag to tell mq() and dv1calc() which observable matrix to calculate
      std::string basis;                     // Name of basis to output eigenvectors, supported: "JmJ" and "mSmL"
      double Bx,By,Bz;                       // For magnetic field for Zeeman term
      double xT,xHa,xHb,xHc;                 // The vector in (H-T) phase space to calculate the x-axis of phase diag
      double xMin,xStep,xMax;                // Start, step and end of x-axis in the phase diagram
      double yT,yHa,yHb,yHc;                 // The vector in (H-T) phase space to calculate the y-axis of phase diag
      double yMin,yStep,yMax;                // Start, step and end of y-axis in the phase diagram
      #ifdef JIJCONV
      std::vector<double> jijconv;           // Conversion factor for Jij coupling parameters from Stevens/Wybourne norm.
      void jijconvcalc();                    // Calculates the conversion factors above.
      #endif

      bool operator==(icpars c) const;       // Operator to determine if parameters are the same
      bool operator!=(icpars c) const;       // Operator to determine if parameters are not the same

      icpars();                              // Blank constructor
};

#endif
