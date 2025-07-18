/***********************************************************************
 *
 * addj.c - program to add *.j files 
 *
 ***********************************************************************/


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
                Options: -nofcomponents 23  ...fixes the nofcomponents to 23 by \n \
                        reducing (removing entries) or increasing (by filling with zeroes) \n \
                        the exchange parameter tables\n \
                        -ni                 ... forces output without indexchange \n \
                        -s  0.2             ... scales all interactions (including \n \
                                                magnetoelastic G) in file1 by 0.2 before adding \n \
                        -rmcomp 9 15        ... removes components 8 to 15 from result before output\n \
                        -pi 5 7  1 2 1 3 print to stdout only interaction tensor (rows 1-2, columns 1-3) of ion 5 with neighbour 7 as matrix \n \
                        -v  verbose \n \
                If file2.j is not given, a copy of the input file is saved \n \
               \n");
      exit (1);
    } else {fprintf(stderr, "#* addj 250712 *\n");}
 int ow=1; int n=0,noindexchange=0,rml=0,rmh;double scale=1.0;int verbose=0;
 int pa=0,pi,prl,prh,pcl,pch;
 while(argv[ow][0]=='-'){
 if(strcmp(argv[ow],"-nofcomponents")==0){ow+=1;
 // option setting nofcomponents
 n=(int)strtod(argv[ow],NULL);
 if(n<1){fprintf(stderr,"Error program add option nofcomponents=%i is less than 1\n",n);exit(1);}
                                        }
 if(strcmp(argv[ow],"-s")==0){ow+=1;scale=strtod(argv[ow],NULL);}
 if(strcmp(argv[ow],"-rmcomp")==0){ow+=1;rml=(int)strtol(argv[ow], (char **)NULL, 10);
                                   ow+=1;rmh=(int)strtol(argv[ow], (char **)NULL, 10);
                             }
 if(strcmp(argv[ow],"-ni")==0){noindexchange=1;}
 if(strcmp(argv[ow],"-v")==0){verbose=1;}
 if(strcmp(argv[ow],"-pi")==0){
                              ow+=1;pa=(int)strtol(argv[ow], (char **)NULL, 10);
                             ow+=1;pi=(int)strtol(argv[ow], (char **)NULL, 10);
                             ow+=1;prl=(int)strtol(argv[ow], (char **)NULL, 10);
                             ow+=1;prh=(int)strtol(argv[ow], (char **)NULL, 10);
                             ow+=1;pcl=(int)strtol(argv[ow], (char **)NULL, 10);
                             ow+=1;pch=(int)strtol(argv[ow], (char **)NULL, 10);
                              }

 ++ow;}

 par a(argv[ow],verbose);a.scale(scale);

 if(n>0){a.set_nofcomponents(n);}  

 if(argc-ow>1){par b(argv[1+ow],verbose);
 if(n>0){b.set_nofcomponents(n);}  
               a.add(b);}

  if(rml>0){a.remove_components(rml,rmh,verbose);}
if(pa) a.print_interaction(stdout,pa,pi,prl,prh,pcl,pch);
else
 a.save(stdout,noindexchange);
}


