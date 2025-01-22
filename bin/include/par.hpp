// class of cf and exchange parameters for a given 
// crystal: corresponds to paramters given in mcphas.j
#ifndef PAR
#define PAR


#include<martin.h>
#include "jjjpar.hpp"
#include "cryststruct.hpp"

#define MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL 4000
#define SMALL 1e-6  // for module kramer - to trigger numerical limited calculation 
                    // for adding jjpar sets to see what is difference in position or what is equal
                    // for checking jjj parameters if values are equal 
#define SMALLCHARGE 1e-6  // for checking if totalcharge is small and outpu polarisation P

class par
{ 
  private:
   
   
  public:
  
  char *rems[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL];
  jjjpar **jjj; // pointer of field of exchange parameter sets
  //Vector gJ;
     
  //lattice
  Matrix rez;
  cryststruct cs;

   //elastic constants in Voigt notation in meV/atom, i.e. multiply by nofatoms to get 
   // elastic energy for a crystallographic primitive unit cell as described by
   // mcphas.j
   Matrix Cel,CelInv;

   // sum of charges 
   double totalcharge;   

  //jjjpar 
   
   par (const char *filejjj,int verbose=0);	//konstruktor
   par (Vector abc,int nofcompi); // simple contructor: needed in clusterize
   par (const par & pars);	// kopier-konstruktor
   
~par ();		//destruktor

int newatom(jjjpar * p); //creates new atom from an existing and returns its index
int delatom(int n, Matrix & distribute,int verbose); //removes atom n and returns new nofatoms
// if n<0 then atom number |n| is removed and also all interactions of other atoms
// with this atom are removed from the interaction table 
// if n>0 interactions with the other atoms are kept and transferred to 
// a group of atoms (numbers given in column 1 of distribute) with 
// coefficients given in col 2 of distribute. Only interactions
// with atoms given in column 1 of distribute are removed completely.
// Attention: sublattice index is not changed by this function !
//            -- thus after running it sublattice[s] refers still to
//            original numbering with all atoms in the parameters set !

void reduce_unitcell(int verbose);//checks every atom in the unit cell and removes
                       // any atom, which is connected to another by a lattice vector
void add(par & b); // add exchange parameters
void scale(double scalefactor); // scale all interaction parameters by scalefactor
void save(FILE * fout,int noindexchange); // save lattice, atoms and exchange parameters to file
void save(const char * filename,int noindexchange); // save lattice, atoms and exchange parameters to file
void savelattice(FILE *fout);// save lattice to file
void saveatoms(FILE *fout);// save atom positions and properties  to file
void save_sipfs(const char *path);   //save single ion parameter files filename to path*
void save_mcdiff_in (const char * program); // save structure in mcdiff.in program is program name calling this
void set_nofcomponents (int n); //sets the number of components in the interaction vector
void increase_nofcomponents (int n); //increases the number of components in the interaction vector
void decrease_nofcomponents (int n); //decreases the number of components in the interaction vector

// operator!= returns 8 7 6 5 4 3 2 1 0depending on agreement of
 //  8 abc 7 nofatoms 6 atomic positions 5 sipffilenames 4 nofcomponents 3 nofneighbours disagreement
 //  2 neighbour position 1 interaction parmeter disagreement i.e. 0 is perfect match
 int operator!=(par & op2); // match

};
#endif
