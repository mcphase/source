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
  if (argc <= 1)
    { printf (" program addj - add exchange parameter file2 to file1, output is written to stdout\n \
                use as: addj [options ]file1.j [file2.j]\n \
                Options: -nofcomponents 23 fixes the nofcomponents to 23 by \n \
                        reducing (removing entries) or increasing (by filling with zeroes) \n \
                        the exchange parameter tables\n \
                        -ni  forces output without indexchange \n \
                        -s  0.2 scales all interactions in file1 by 0.2 before adding \n \
                        -v  verbose \n \
                If file2.j is not given, a copy of the input file is saved \n \
               \n");
      exit (1);
    } else {fprintf(stderr, "#* addj 221011 *\n");}
 int ow=1; int n=0,noindexchange=0;double scale=1.0;int verbose=0;
 while(argv[ow][0]=='-'){
 if(strcmp(argv[ow],"-nofcomponents")==0){ow+=1;
 // option setting nofcomponents
 n=(int)strtod(argv[ow],NULL);
 if(n<1){fprintf(stderr,"Error program add option nofcomponents=%i is less than 1\n",n);exit(1);}
                                        }
 if(strcmp(argv[ow],"-s")==0){ow+=1;scale=strtod(argv[ow],NULL);}

 if(strcmp(argv[ow],"-ni")==0){noindexchange=1;}
 if(strcmp(argv[ow],"-ni")==0){verbose=1;}
 ++ow;}

 par a(argv[ow],verbose);a.scale(scale);

 if(n>0){a.set_nofcomponents(n);}  

 if(argc-ow>1){par b(argv[1+ow],verbose);
 if(n>0){b.set_nofcomponents(n);}  
               a.add(b);}

 a.save(stdout,noindexchange);
}


