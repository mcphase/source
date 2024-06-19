//libraries, classes and  functions in mcphas
#include "par.hpp"
#include "qvectors.hpp"
#include "physprop.hpp"
#include "testspin.hpp"
#include "spincf.hpp"
#include "mfcf.hpp"
#include "inipar.hpp"
#include <complex>
#include "martin.h"

 // normalize Vector dadbdc to length n
int normalizedadbdc(Vector & dadbdc,double n,par & inputpars);

// routines for mcphas for calculation of magnetic phases
//htcalc.h
#include <cerrno>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector.h>
     //calculation of physical properties at given HT point

/*check if parameter file mcphas.ini has been changed*/
void checkini(testspincf & testspins,qvectors & testqs,inipar & ini);


/* calculate magnetic structure at a given HT- point
  on input:
    T	Temperature[K]
    H	Vector of External Magnetic Field [T]
    inputpars	Input parameters (exchange constants etc...)
    testqs	Set of propagation vectors to be tested
    testspins	Set of Spinconfigurations to be tested
  on return:
    physprops	physical properties at (HT) point (i.e. magnetic structure
		neutron intensities, thermal expansion ...)
 // returns 0 if successfull
 // returns 1 if too maxnofspinconfigurations is exceeded
 // returns 2 if no spinconfiguration has been found at ht point
 */
int htcalc(Vector H,double T,inipar & ini, par & inputpars,qvectors & testqs,
                     testspincf & testspins,physproperties & physprops);

// calculate physical properties for stabilized HT point
void physpropclc(Vector H,double T,
                               spincf & sps,mfcf & mf,physproperties & physprops,inipar & ini,par & inputpars);

// here the free energy is calculated for a given (initial) spinconfiguration
// using the meanfield algorithm
/*on input:
    T		Temperature[K]
    Hex		Vector of external magnetic field [T]
    inputpars	exchange and other parameters
    sps		initial spinconfiguration
    testspins	all other testspinconfigurations
  on return
    returns free energy[meV]
    sps		selfconsistently stabilized spinconfiguration (may be different
		from initial spinconfiguration)
    u		mangetic energy[meV]

 */
double fecalc(double & U, double & Eel,int & r,double & sc,Vector H,double T,inipar & ini,par & inputpars,
            spincf & sps,mfcf & mf,testspincf & testspins,qvectors & testqs);



// this sub checks if a spinconfiguration has already been added to
// table testspins and adds it if necessary
int checkspincf(int j,spincf & sps,qvectors & testqs, Vector & nettom,
		    Vector & momentq0, Vector & phi,
                     testspincf & testspins,physproperties & physprops,inipar & ini);


// subroutine to calculate magnetisation M from effective field H
// this is the heart of the meanfield algorithm an it is necessary to
// keep this routine as efficient as possible
// at the moment we do only groundstate doublet
/*on input
    A,M,Ci	saturation moment/gJ[muB] of groundstate doublet in a.b.c direction
    g_J		lande factor
    T		temperature[K]
    H		vector of mean field
  on output
    m		single ion magnetic moment vector [muB]
    Z		single ion partition function
    U		single ion magnetic energy
*/
//void Icalc (float & A, float & M, float & C,float & gJ,double & T,
//                          Vector & H, Vector &  m, float & Z,float & U);







