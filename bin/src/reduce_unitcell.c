/***********************************************************************
 *
 * reduce_unitcell.c - program to reduce unit cell by removing 
 *
 ***********************************************************************/


#define MAXNOFATOMS 100
#include "par.hpp"
#include "martin.h"

/**********************************************************************/
// hauptprogramm
int main (int argc, char **argv)
{ 
// check command line
  if (argc <= 1)
    { printf (" program reduce_unitcell, output is written to stdout\n \
                use as: reduce_unitcell mcphas.j\n\n \
                This program checks every atom in the unit cell in file mcphas.j and removes\n \
                any atom, which is connected to another by a lattice vector.\n \
                \n");
      exit (1);
    }

 par a(argv[1]);

 a.reduce_unitcell();  

 a.save(stdout);
}


