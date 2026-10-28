
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
#define KB 0.08617343183   // Boltzmanns constant in mev/K
#define MU_B  0.0578838263 // Bohrmagneton in meV/tesla
#define ECHARGE 1.602189e-19   // elementary charge |e| in Coulomb
#define MU0    4*PI*1e-7    // mu0 in Vs/Am  
#define A0     5.29177210903e-11  // Bohr radius in meter
#define HBAR   6582.119569e-16   // Planck constant hbar in meVs
#define MEV2J  1.6021766339999e-22 // 1mev=1.6021766339999e-22J
#define M0     1.660539e-27       //m0=atomic mass unit=1.660539e-27 kg

#define MAXNOFCHARINLINE 7024 
#define MAXNOFLINES 7024  // max nof lines in mcphas.ini
#define SMALL 1e-6  // (meV) regulates if energy is treated as degenerate or not after a diagonalisation
                    // used in many modules: singleion_module.hpp mcdisp mcphas ... change with caution !!!
                    // ! this is a central switch !
                   // must match SMALL in mcdisp.c and ionpars.cpp because it is used to decide wether for small
                   // transition, energy the matrix Mijkl contains wn-wn' or wn/kT
                   // also for module kramer - to trigger numerical limited calculation 
                    // also for adding jjpar sets to see what is difference in position or what is equal
                    // also for checking jjj parameters if values are equal 
#define SMALLENERGYOFFSET_PHONON 0.0 // a small offset in calculating the single ion excitations
                    // gap delta in the phonon module  to avoid numerical proiblems at q=0
#define SMALLPOSITIONDEVIATION 1e-4  // used in mcdiff mcphas to see if a lattice atom position matches or not
#define HUGE_EXP   200   // huge value of E/kT in exp(E/kT)
#define MAXNOFINIS 100  // maximum number of different prefixes in mcphas.ini mcdisp.ini

#define HEXT_DIMENSION  12  // dimension of external field Hext (Hi Hj Hk Ei Ej Ek s1 s2 s3 s4 s5 s6)
#define NOF_USERDEF_MCPHAS_COLS 7  // number of user defined columns in mcphas.ini -> mcphas.* outpus files
#define SMALL_FIELD 1e-6     // nonzero value of external field E or H
