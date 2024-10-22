/***********************************************************************
 *
 * mcdiff - program to calculate neutron and magnetic xray diffraction
 *
 * reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R
 ***********************************************************************/
#include <mcdiff.h>
#include <inimcdiff.hpp>
#include "mcdiff_intcalc.c"
#include "mcdiff_output.c"


void erxit() // type info and error exit 
{ printf (" \n \
                use as: mcdiff [options][hkllistfilename [-z]]\n \
		- for format of input file mcdiff.in see mcphase manual\n \
                - optional an hkl list can be given in file hkllistfilename-in\n \
                  this case the program computes reflections in this list.\n \
                  if it is a 3column list,the azimuth dependence of\n \
                  magnetic scattering is calculated.If neutron intensities\n \
                  are listed in column 4 the program computes\n \
                  rpvalue=100 * sum_i |Icalc_i-Iexp_i|/sum_i |Iobs_i|\n \
                  and does not output the azimuth dependence \n \
                  if in addition experimental errors are given in column 5,\n \
                  chisquared=1/N *sum_i (Icalc_i - Iexp_i)^2/err_i^2 is\n \
                  calculated also.\n \
                - options: -h  ... prints this help message\n \
                           -v  ... verbose - output more information\n \
                           -prefix 001   ... prefix for parameters to be read from mcdiff.in\n \
                                   and used for creation of output file (useful for running in parallel\n \
                                   calculations for different zones: e.g. put in mcdiff.in instead of\n \
                                   #!lambda= several statements #!001lambda= ... #!002lambda=\n \
                                   and start several jobs of mcdisp with -prefix 001, -prefix 002 \n \
                                   simultaneously creating output fiels results/001mcdiff.out \n \
                                   results/002mcdiff.out. \n \
                                   Note: if there exists a file 001mcdiff.in then mcdiff will read\n \
                                   input from this file instead from mcdiff.in when using this option.\n \
                           -z  ... output zero intensity if (hkl) in hkllist does\n \
                                   not correspond to reciprocal lattice of supercell\n \
                                   (default: move hkl to nearest rec latticepoint)\n \
                - results are saved in results/mcdiff.out\n");
        exit (1);
}

// hauptprogramm
int main (int argc, char **argv)
{ 
 // test to test threej function
/* if (argc>6) {printf ("cint(%g)=%i\n", strtod(argv[6],NULL),cint(strtod(argv[6],NULL)));
		// test matpack routine  
                    int n,ndim=20; double thrcof[20];int errflag;
		    double min,max;
		     ThreeJSymbolM	(strtod(argv[1],NULL),
	             strtod(argv[2],NULL),
	             strtod(argv[3],NULL),
	             strtod(argv[4],NULL),
	             min,max, thrcof, ndim, 
			 errflag);
		     n=(int)(strtod(argv[5],NULL)-min);	 

              printf ("threej symbol=%g=%g\n",
              threej(strtod(argv[1],NULL),
	             strtod(argv[2],NULL),
	             strtod(argv[3],NULL),
	             strtod(argv[4],NULL),
	             strtod(argv[5],NULL),
	             strtod(argv[6],NULL)
	             ),thrcof[n]);
		     return 0;}
*/
// test spherical harmonics
/*if (argc>3) {
  int l,m; double theta,phi;
  l=(int)strtod(argv[1],NULL);
  m=(int)strtod(argv[2],NULL);
  theta=strtod(argv[3],NULL);
  phi=strtod(argv[4],NULL);
  
  printf("Y%i%i(%g,%g)=%g %+g i\n",l,m,theta,phi,real(SphericalHarmonicY (l,m,theta,phi)),imag(SphericalHarmonicY (l,m,theta,phi)));
 return 0;}
*/
  std::clock_t startcputime = std::clock();
  FILE * fin;
  char hklfile [MAXNOFCHARINLINE]; //input hkl list file name
  hklfile[0]='\0';
  char prefix [MAXNOFCHARINLINE];prefix[0]='\0';
  int verbose=0;
fprintf(stderr,"***********************************************************************\n");
fprintf(stderr,"*\n");
fprintf(stderr,"* mcdiff - program to calculate neutron and magnetic xray diffraction\n");
fprintf(stderr,"*\n");
fprintf(stderr,"* reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
fprintf(stderr,"***********************************************************************\n");
   int zeronotmatchinghkl=0;
  // check command line
  for (int i=1;i<=argc-1;++i){
   if(strcmp(argv[i],"-h")==0) {erxit();} // print help end exit
                          else {if(strcmp(argv[i],"-prefix")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -prefix needs argument(s)\n");exit(EXIT_FAILURE);}
  		                                  strcpy(prefix,argv[i+1]);++i;
 						  fprintf(stdout,"#prefix for reading parameters from mcdiff.in and for ouput filenames: %s\n",prefix);
 					                         }
                           else {if(strcmp(argv[i],"-v")==0) {verbose=1;}
                           
                            else {strcpy(hklfile,argv[i]);
                                 if (argc>i+1){if(strcmp(argv[i+1],"-z")==0){zeronotmatchinghkl=1;}}
                                 } // hklfile
                                } // -v
                               } // prefix
                              } // for
                                    

inimcdiff ini("mcdiff.in",prefix,verbose);

float * out[ ini.nofoutputcolumns+1];
for(int i=1;i<= ini.nofoutputcolumns;++i)
 {out[i]=new float [MAXNOFREFLECTIONS+1];
  if(out[i]==NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
 }
float totint[MAXNOFREFLECTIONS+1];
Vector * hkl;float * D;
complex <double> * mx;
complex <double> * my;
complex <double> * mz;
complex <double> * mxmy;
complex <double> * mxmz;
complex <double> * mymz;
complex <double> * mx2;
complex <double> * my2;
complex <double> * mz2;

hkl = new Vector[MAXNOFREFLECTIONS+1];
 for(int i=0;i<=MAXNOFREFLECTIONS;++i){hkl[i]=Vector(1,3);}
 D= new float[MAXNOFREFLECTIONS+1];
 
if(ini.colcod[0]<0) // for magnetic xray scattering we need some more storage ...
{mx= new complex <double> [MAXNOFREFLECTIONS+1];
 my= new complex <double> [MAXNOFREFLECTIONS+1];
 mz= new complex <double> [MAXNOFREFLECTIONS+1];
 mxmy= new complex <double> [MAXNOFREFLECTIONS+1];
 mxmz= new complex <double> [MAXNOFREFLECTIONS+1];
 mymz= new complex <double> [MAXNOFREFLECTIONS+1];
 mx2= new complex <double> [MAXNOFREFLECTIONS+1];
 my2= new complex <double> [MAXNOFREFLECTIONS+1];
 mz2= new complex <double> [MAXNOFREFLECTIONS+1];
}

Vector hhkkll(1,3);
int code=0;int m=0;

// if hkllist is given, read the file and put hkls to hkl[i], m is number of reflections to be considered
if (hklfile[0]!='\0'){int nr;
      float nn[20];nn[0]=19;
     // open hkllist file
     fprintf(stdout,"reading hkl list from file %s\n",hklfile);
       fin = fopen_errchk (hklfile, "rb");
       while(feof(fin)==false){nr=inputline(fin,nn);
                               if(nr>2)
                               {hhkkll(1)=nn[1];hhkkll(2)=nn[2];hhkkll(3)=nn[3];++m;
                                code=1;  // set >0 to indicate that hkl list is given  on input                            
                               // transformieren der millerindizes auf magnetische einheitszelle(re-use nn[1-3] to store noninteger values
                                  nn[1]=hhkkll*(ini.rtoijk.Inverse()*ini.r1);
                                  nn[2]=hhkkll*(ini.rtoijk.Inverse()*ini.r2);
                                  nn[3]=hhkkll*(ini.rtoijk.Inverse()*ini.r3);
                                  hkl[m](1)=rint(nn[1]);// round to integer reciprocal lattice point
                                  hkl[m](2)=rint(nn[2]);
                                  hkl[m](3)=rint(nn[3]);
                                 // check if magnetic reflection is indeed on magnetic reciprocal lattice
                              if(fabs(nn[1]-hkl[m](1))>SMALLPOSITIONDEVIATION||fabs(nn[2]-hkl[m](2))>SMALLPOSITIONDEVIATION||fabs(nn[3]-hkl[m](3))>SMALLPOSITIONDEVIATION)
                                {fprintf(stderr,"Warning mcdiff - reading hkl=(%g %g %g): calculation impossible, because this corresponds to ", hhkkll(1),hhkkll(2),hhkkll(3));
                                 fprintf(stderr,"non integer supercell reciprocal lattice point (%g %g %g)", nn[1], nn[2], nn[3]);
                                 if(zeronotmatchinghkl==1)
                                 {hkl[m](1)=nn[1];// do not round to integer reciprocal lattice point
                                  hkl[m](2)=nn[2];
                                  hkl[m](3)=nn[3];
                                  fprintf(stderr,"- will output zero intensity\n");
                                 }
                                 else
                                 {// transform integer back to hkl
                                   hhkkll=hkl[m](1)*ini.rez1+hkl[m](2)*ini.rez2+hkl[m](3)*ini.rez3;
                                   hhkkll/=2.0*PI;
                                   nn[1]=hhkkll*ini.rtoijk.Column(1);
                                   nn[2]=hhkkll*ini.rtoijk.Column(2);
                                   nn[3]=hhkkll*ini.rtoijk.Column(3);
                                 fprintf(stderr,"\n - will calculate hkl=(%g %g %g) instead !\n\n", nn[1], nn[2], nn[3]);
                                 }
                                }
                                if(nr>3){mx[m]=complex <double> (nn[4],0);code=2;}// intensities given                                
                                if(nr>4){my[m]=complex <double> (nn[5],0);code=3;}// errors given                                
                              }}
       fclose(fin);      
           }

printheader(ini,code,m);

neutint(ini,code,m,hkl,D,totint, out,mx,my,mz,mxmy,mxmz,mymz,mx2,my2,mz2);

printreflist(ini,code,m,hkl,totint,D, out,mx,my,mz,mxmy,mxmz,mymz,mx2,my2,mz2);

double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
std::cout << "#! Finished in cputime=" << cpu_duration << " seconds [CPU Clock] " << std::endl;
std::cout << "#! nofhkls=" << m << " different q vectors generated " << std::endl;

fprintf (stderr,"...results written to %s\n",ini.outfilename);
fprintf (stderr,"***********************************************************\n");
fprintf (stderr,"                   End of Program mcdiff\n");
fprintf (stderr,"reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
fprintf (stderr,"***********************************************************\n");
if(ini.colcod[0]<0){
  delete []mx;delete []my;delete []mz;delete []mxmy;
  delete []mxmz;delete []mymz;delete []mx2;delete []my2;delete []mz2;}
delete []hkl;delete []D;
  for(int i=1;i<= ini.nofoutputcolumns;++i){delete []out[i]; } 
 return 0;
}


