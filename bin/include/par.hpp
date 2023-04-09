// class of cf and exchange parameters for a given 
// crystal: corresponds to paramters given in mcphas.j
#ifndef PAR
#define PAR


#include<martin.h>
#include "jjjpar.hpp"

#define MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL 4000
#define SMALL 1e-6  // for module kramer - to trigger numerical limited calculation 
                    // for adding jjpar sets to see what is difference in position or what is equal
                    // for checking jjj parameters if values are equal 

class par
{ 
  private:
   
   
  public:
  
  char *rems[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL];
  jjjpar **jjj; // pointer of field of exchange parameter sets
  //Vector gJ;
     
  //lattice
  float a,b,c,alpha,beta,gamma;
  Vector abc;
  Matrix r,rez;
  int nofatoms;
  int nofcomponents;

   //elastic constants in Voigt notation in meV/atom, i.e. multiply by nofatoms to get 
   // elastic energy for a crystallographic primitive unit cell as described by
   // mcphas.j
   Matrix Cel,CelInv;
   

  //jjjpar 
   
   par (const char *filejjj,int verbose=0);	//konstruktor
   par (float ai,float bi,float ci,float alphai,float betai,float gammai,int nofcompi);
   par (const par & pars);	// kopier-konstruktor
   
~par ();		//destruktor

int newatom(jjjpar * p); //creates new atom from an existing and returns its index
int delatom(int n); //removes atom n and returns new nofatoms
void reduce_unitcell();//checks every atom in the unit cell and removes
                       // any atom, which is connected to another by a lattice vector
void add(par & b); // add exchange parameters
void scale(double scalefactor); // scale all interaction parameters by scalefactor
void save(FILE * fout,int noindexchange); // save lattice, atoms and exchange parameters to file
void save(const char * filename,int noindexchange); // save lattice, atoms and exchange parameters to file
void savelattice(FILE *fout);// save lattice to file
void saveatoms(FILE *fout);// save atom positions and properties  to file
void set_nofcomponents (int n); //sets the number of components in the interaction vector
void increase_nofcomponents (int n); //increases the number of components in the interaction vector
void decrease_nofcomponents (int n); //decreases the number of components in the interaction vector

void save_sipfs(const char *path);   //save single ion parameter files filename to path*

// operator!= returns 8 7 6 5 4 3 2 1 0depending on agreement of
 //  8 abc 7 nofatoms 6 atomic positions 5 sipffilenames 4 nofcomponents 3 nofneighbours disagreement
 //  2 neighbour position 1 interaction parmeter disagreement i.e. 0 is perfect match
 int operator!=(par & op2); // match

};
#endif
