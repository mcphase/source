/***********************************************************************
 *
 * coq2jjj.c - program to bring mcphas.coq file to new format of exchange 
 *                    parameter input file
 *
 ***********************************************************************/

#include "paramete.hpp"
#include <martin.h>
/**********************************************************************/
/**********************************************************************/
// main program
int main (int argc, char **argv)
{
 
// check command line
  if (argc < 1)
    { printf ("use as: coq2jjj *.coq - parameters from file *.coq are reformatted to stdout");
      exit (1);
    }

  
// load  parameters from file  (as object of class parameters)
   parameters inputpars (argv[1]);

// save parameters
   inputpars.savjjj(stdout);
 

  return 0;
}


