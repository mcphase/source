//libraries, classes and  functions in mcdisp


#define SMALL_NORM 1e-10 // small norm2 of gamma - if norm2(gamma) less than that, it will not be normalised
#define SMALL_XPROD_FOR_PARALLEL_VECTORS   1e-10  // small value of xproduct to decide if vectors are parallel
#define SMALL_GAMMA   1e-10  // small value of transition matrix element gamma
#define SMALL_QUASIELASTIC_ENERGY 1e-6    // deviation from single ion gap delta to take energy into account as not being equal to
                      // delta and therefore being included into output 
		      // transitions of single ions less then SMALL_QUASIELASTIC_ENERGY have in Mijkl wn/kT instead of wn-wn'
		      // SMALL_QUASIELASTIC_ENERGY is also used in jjjpar.cpp and ionpars.hpp !!!!
                      // in ionpars.hpp it is used to decide wether for SMALL_QUASIELASTIC_ENERGY transition
 		     // energy the matrix Mijkl contains wn-wn' or wn/kT
#define SMALLINT 1e-4 // small intensity treshhold in barn/f.u. - only peaks larger than this threshold will
                      // be included in sta_int calculation
#define ANTIPEAK_CUTOFF 100 // if antipeaks are fitted,  sta contribution never increases above this cutoff
#define SMALLEDIF 1e-5 // small difference in calculation of transition energy
                       // used to give error if recalculation of mcdisp.trs
		       // energies gives different results than file





