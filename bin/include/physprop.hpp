// class physprop for the storage of the physical properties
// of the system at given H and T

#ifndef PHYSPROP
#define PHYSPROP


#include<par.hpp>
#include<martin.h>
#include<spincf.hpp>
#include<mfcf.hpp>

class physproperties
{
  private:
 int washere,nofspincorr;
    
  public:
float x,y; // phasediagramm labels  
int j;  // index of spinstructure
double T; // temperature
Vector m,H; // moment and H feld
double fe;
double u; // free energy and mag energy per ion
int nofatoms;
int nofcomponents;

Vector *jj,*hkli; // spin spin correlation functions
int maxnofhkls,nofhkls;
spincf  sps;
mfcf mf;
   
physproperties (int nofspincorrs,int maxnofhkls,int na, int nm);	//konstruktor
//na number of atoms in basis,nm number of spin components
physproperties (const physproperties & props);	// kopier-konstruktor

~physproperties ();		//destruktor

void update_maxnofhkls(int mxnofhkli);
double save(int verbose,const char * filemode, int j,par & inputpars,char * prefix);

};



#endif
