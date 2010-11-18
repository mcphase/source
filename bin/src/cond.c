/***********************************************************************
 *
 * cond.c - program to modify parameters of *.coq file to obey 
 *                    special condition
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
    { printf ("use as: cond *.coq  - parameters from file *.coq are subject to constraints");
      exit (1);
    }

  
// load  parameters from file  (as object of class parameters)
   parameters inputpars (argv[1]);

//--------------------------------------------------------------
// here the condition for  the parameters is stated ( - CeCu2 ...)
 int i;
 float sum=0;
 float mult;
 const float kb=0.0862; 
 const float Tn=3.5;
 for (i=1;i<=inputpars.paranz;++i)
 {// determine multiplicity
  mult=1;
  if (inputpars.dn[i](1)!=0){mult*=2.0;}
  if (inputpars.dn[i](2)!=0){mult*=2.0;}
  if (integer((inputpars.dn[i](1)+inputpars.dn[i](2))/1)!=
                  (inputpars.dn[i](1)+inputpars.dn[i](2))/1){mult*=-1;}
  sum+=mult*inputpars.jij[i](3);
 }
 inputpars.jij[9](3)+=kb*Tn/inputpars.C/inputpars.C-sum;
//-----------------------------------------------------------------------


// save parameters
   inputpars.save(stdout);
 

  return 0;
}


