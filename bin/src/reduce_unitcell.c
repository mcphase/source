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
			-delatoms 1:3+-0.75:4+-0.25,2,5,7  removes atoms number 1,2,5,7 for atom 1 the interactions\n \
			 are transferred to 75%% to atom 3 and 25%% to atom 4 and interactions to atoms 3 and 4 are\n \
			 removed, for 2 5 7 all the interactions with other atoms are removed\n \
			-v  verbose mode\n \
                \n");
      exit (1);
    } else { fprintf (stderr,"#* reduce_unitcell 221011 *\n"); }

int ow=1,i=0; int n=0,noindexchange=0,verbose=0;
char * token;
char *substr[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL+1];
float ns[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL+1];ns[0]=MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL;
float nscoeff[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL+1];nscoeff[0]=MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL;

 while(argv[ow][0]=='-'){
 if(strcmp(argv[ow],"-nofcomponents")==0){ow+=1;
 // option setting nofcomponents
 n=(int)strtod(argv[ow],NULL);
 if(n<1){fprintf(stderr,"Error program add option nofcomponents=%i is less than 1\n",n);exit(1);}
                                        }
if(strcmp(argv[ow],"-delatoms")==0){ow+=1;
// substitute all commma with spaces
  while ((token=strrchr(argv[ow],','))!=NULL){++i;substr[i]=token+1;*token='\0';}
++i;substr[i]=argv[ow];

}

 if(strcmp(argv[ow],"-ni")==0){noindexchange=1;}
 if(strcmp(argv[ow],"-v")==0){verbose=1;}
 ++ow;}

 par a(argv[ow]);

if(n>0){a.set_nofcomponents(n);if(verbose){fprintf(stderr,"Setting nofcomponents=%i\n",n);}}


if(i==0){
 a.reduce_unitcell(verbose);  
}else{
for(n=1;n<=i;++n)
  {
 while ((token=strchr(substr[n],':'))!=NULL){*token=' ';}
//fprintf(stderr,"string is %s\n",substr[n]);
 int nn=splitstring(substr[n],ns,nscoeff);
 int dim=1;if(nn>1){dim=nn-1;}
   Matrix dis(1,dim,1,2);int an=ns[1];
 for(int ii=2;ii<=nn;++ii){
// do not use next line, but take care of already deleted atoms from previous loops (renumbering of atoms) ...
//dis(ii-1,1)=ns[ii];
int nnew=ns[ii];
for(int ic=1;ic<n;++ic)
 {int ndone=(int)strtod(substr[ic],NULL);if(ndone<ns[ii])--nnew;
 }

dis(ii-1,1)=nnew;
dis(ii-1,2)=nscoeff[ii];
                          }
   if(n<i){int anp1=strtod (substr[n+1], NULL);if((int)an<=(int)anp1){fprintf(stderr,"Error program reduce_unitcell option -delatoms atom numbers have to be given in ascending order\n");exit(1);}}
if(nn>1){   a.delatom((int)an,dis,verbose);
}else{ a.delatom(-(int)an,dis,verbose);}
  }
}
 a.save(stdout,noindexchange);
 a.save_mcdiff_in("reduce_unitcell");

fprintf(stderr,"# end of reduce_unitcell - created mcdiff.in and list of redundant sipf files\n");
fprintf(stderr,"# in file reduce_unitcell_sipf.del, to delete these files use:\n");
fprintf(stderr,"# perl -l -n -e \"unlink\" reduce_unitcell_sipf.del\n");
}


