/***********************************************************************
 *
 * reduce_unitcell.c - program to reduce unit cell by removing 
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
    { printf (" program reduce_unitcell, output is written to stdout\n \
                use as: reduce_unitcell [option] mcphas.j\n\n \
                This program checks every atom in the unit cell in file mcphas.j and removes\n \
                any atom, which is connected to another by a lattice vector.\n \
                a list of superfluous sipf file is stored in reduce_unitcell_sipf.del \n \
                Options: -nofcomponents 23 fixes the nofcomponents to 23 by \n \
                        reducing (removing entries) or increasing (by filling with zeroes) \n \
                        the exchange parameter tables\n \
                        -ni  forces output without indexchange \n \
			-delatoms 1,2,5,7    instead of removing atoms connected by a lattice vector \n \
			remove atoms number 1,2,5 and 7 from the list and also all interactions with those \n \
                \n");
      exit (1);
    } else { fprintf (stderr,"#* reduce_unitcell 221011 *\n"); }

int ow=1; int n=0,noindexchange=0;
char * token;
float ns[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL+1];ns[0]=MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL;
 while(argv[ow][0]=='-'){
 if(strcmp(argv[ow],"-nofcomponents")==0){ow+=1;
 // option setting nofcomponents
 n=(int)strtod(argv[ow],NULL);
 if(n<1){fprintf(stderr,"Error program add option nofcomponents=%i is less than 1\n",n);exit(1);}
                                        }
if(strcmp(argv[ow],"-delatoms")==0){ow+=1;
// substitute all commma with spaces
  while ((token=strchr(argv[ow],','))!=NULL){*token=' ';}
// store atoms to be deleted in ns
  ns[0]=splitstring(argv[ow],ns);
}

 if(strcmp(argv[ow],"-ni")==0){noindexchange=1;}
 ++ow;}

 par a(argv[ow]);

if(ns[0]==MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL){
 a.reduce_unitcell();  
}else{
for(n=ns[0];n>0;--n)
  {a.delatom(-(int)ns[n]);
  }
}
 a.save(stdout,noindexchange);
 a.save_mcdiff_in("reduce_unitcell");

fprintf(stderr,"# end of reduce_unitcell - created mcdiff.in and list of redundant sipf files\n");
fprintf(stderr,"# in file reduce_unitcell_sipf.del, to delete these files use:\n");
fprintf(stderr,"# perl -l -n -e \"unlink\" reduce_unitcell_sipf.del\n");
}


