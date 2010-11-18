/***********************************************************************
 *
 * addj.c - program to add *.j files 
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
  if (argc <= 2)
    { printf (" program addj - add exchange parameter file2 to file1, output is written to stdout\n \
                use as: addj file1.j file2.j\n");
      exit (1);
    }

 par a(argv[1]);
 par b(argv[2]);
  
 a.add(b);

 a.save(stdout);
}


