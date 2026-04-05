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
                options for creating quantum dots, quantum chains and quantum planes: \n \
                        -r1    remove interactions beyond the extended unit cell in r1 direction\n \
                        -r2    remove interactions beyond the extended unit cell in r2 direction\n \
                        -r3    remove interactions beyond the extended unit cell in r3 direction\n \
		                \n");
      exit (1);
    } else { fprintf (stderr,"#* extend_unitcell 260403 *\n"); }

int ow=1; int n=0,noindexchange=0;bool r1=false,r2=false,r3=false;

while(argv[ow][0]=='-'){
 if(strcmp(argv[ow],"-nofcomponents")==0){ow+=1;
 // option setting nofcomponents
 n=(int)strtod(argv[ow],NULL);
 if(n<1){fprintf(stderr,"Error program add option nofcomponents=%i is less than 1\n",n);exit(1);}
                                        }
 if(strcmp(argv[ow],"-ni")==0){noindexchange=1;}
 if(strcmp(argv[ow],"-i")==0){noindexchange=-1;}
 if(strcmp(argv[ow],"-r1")==0){r1=true;}
 if(strcmp(argv[ow],"-r2")==0){r2=true;}
 if(strcmp(argv[ow],"-r3")==0){r3=true;}
 ++ow;}


int n1=atoi(argv[ow]),n2=atoi(argv[ow+1]),n3=atoi(argv[ow+2]);
par a(argv[ow+3]);


if(n>0){a.set_nofcomponents(n);
        //if(verbose){fprintf(stderr,"Setting nofcomponents=%i\n",n);}
        }

 a.extend_unitcell(n1,n2,n3);  Vector nnr123(1,3);
// treat options r1 r2 r3 to remove interactions at unit cell boundary

if(r1||r2||r3)
for(int n=1;n<=a.cs.nofatoms;++n)
 for(int nn=(*a.jjj[n]).paranz;nn>=1;--nn) // count down because delpar will change numbering for parameters > nn
 {nnr123=a.rez*((*a.jjj[n]).xyz+(*a.jjj[n]).dn[nn]); 
//myPrintVector((*a.jjj[n]).dn[nn]);
  // nnr123:  neighbour nn of atom n -  coordinates with respect to primitive lattice
  if(r1&&(nnr123(1)<0||nnr123(1)>1))(*a.jjj[n]).delpar(nn);
  else if(r2&&(nnr123(2)<0||nnr123(2)>1))(*a.jjj[n]).delpar(nn);
  else if(r3&&(nnr123(3)<0||nnr123(3)>1))(*a.jjj[n]).delpar(nn);
 }
 
 a.sort(); // sort output parameters according to ascending distance
 a.save(stdout,noindexchange);
 
fprintf(stderr,"# end of extend_unitcell \n");
}


