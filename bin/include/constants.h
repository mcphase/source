
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
#define KB 0.08617343183   // Boltzmanns constant in mev/K
#define MU_B  0.0578838263 // Bohrmagneton in meV/tesla
#define ECHARGE 1.602189e-19   // elementary charge |e| in Coulomb
#define MU0    4*PI*1e-7    // mu0 in Vs/Am  
   

#define MAXNOFCHARINLINE 7024
#define SMALL 1e-6  // (meV) regulates if energy is treated as degenerate or not after a diagonalisation
                    // used in many modules: singleion_module.hpp mcdisp mcphas ... change with caution !!!
                    // ! this is a central switch !
                   // must match SMALL in mcdisp.c and ionpars.cpp because it is used to decide wether for small
                   // transition, energy the matrix Mijkl contains wn-wn' or wn/kT
                   // also for module kramer - to trigger numerical limited calculation 
                    // also for adding jjpar sets to see what is difference in position or what is equal
                    // also for checking jjj parameters if values are equal 

#define SMALLPOSITIONDEVIATION 1e-4  // used in mcdiff mcphas to see if a lattice atom position matches or not
#define HUGE_EXP   200   // huge value of E/kT in exp(E/kT)