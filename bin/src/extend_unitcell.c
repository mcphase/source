/***********************************************************************
 *
 * extend_unitcell.c - program to extend unit cell  
 *
 ***********************************************************************/


#include "par.hpp"

/**********************************************************************/
// hauptprogramm
/**********************************************************************/
int main (int argc, char **argv)
{ 
// check command line
  if (argc <= 1)
    { printf (" program extend_unitcell, program to extend crystallographic unit cell n1 n2 n3 \n \
                times in r1  r2 r3  direction, output is written to stdout\n \
                use as: extend_unitcell [options] n1 n2 n3 mcphas.j\n\n \
                Options: -nofcomponents 23 fixes the nofcomponents to 23 by \n \
                        reducing (removing entries) or increasing (by filling with zeroes) \n \
                        the exchange parameter tables\n \
                        -i  forces output with indexchange \n \
                        -ni  forces output without indexchange \n \
		                \n");
      exit (1);
    } else { fprintf (stderr,"#* extend_unitcell 260403 *\n"); }

int ow=1; int n=0,noindexchange=0;

while(argv[ow][0]=='-'){
 if(strcmp(argv[ow],"-nofcomponents")==0){ow+=1;
 // option setting nofcomponents
 n=(int)strtod(argv[ow],NULL);
 if(n<1){fprintf(stderr,"Error program add option nofcomponents=%i is less than 1\n",n);exit(1);}
                                        }
 if(strcmp(argv[ow],"-ni")==0){noindexchange=1;}
 if(strcmp(argv[ow],"-i")==0){noindexchange=-1;}
 ++ow;}


int n1=atoi(argv[ow]),n2=atoi(argv[ow+1]),n3=atoi(argv[ow+2]);
par a(argv[ow+3]);


if(n>0){a.set_nofcomponents(n);
        //if(verbose){fprintf(stderr,"Setting nofcomponents=%i\n",n);}
        }

 a.extend_unitcell(n1,n2,n3);  

 a.sort(); // sort output parameters according to ascending distance
 a.save(stdout,noindexchange);
 
fprintf(stderr,"# end of extend_unitcell \n");
}


