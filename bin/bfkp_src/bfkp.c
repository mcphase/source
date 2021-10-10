/**************************************************************
 * bfkq.c - calculate neutron spectra with CF-Phonon interaction
 * Author: Martin Rotter
 **************************************************************/


#include "../../version"
#include "martin.h"
#include "myev.h"

#define K_B  0.0862
/**********************************************************************/
void helpexit()
{ printf (" program bfkp \n "
          " use as:  bfkp T Emin Emax deltaE eta \n"
          "   T    .... Temprature (K)\n"
          "   Emin .... minimal Energy (meV)\n "
          "   Emax .... maximal Energy (meV)\n "
          "   deltaE .... Energy stepwidth (meV)\n "
          "   eta .... imaginary small number for linewidth (meV)\n "   
          "   required: ./results/op.mat for Operator matrices \n  "
          "             ./phonon.int for Frequencies and Intensities of Phonon\n"
          "             ./g_alpha.s for cf-phonon coupling constnats\n"
          "  Output to stdout: nuclear coherent inelastic spectrum (equ. 68 ff) \n"
          
          );
      exit (1);
}

// hauptprogramm
int main (int argc, char **argv)
{int i,imax,k,d=0,n,m,smax=0; double T,Emin, Emax, deltaE, eta;
 FILE * fin_opmat;
 ComplexMatrix *opmatM[52]; double omegaphon[1000],intensity[1000];
printf("#**************************************************************\n");
printf("# * bfkq.c - calculate CF-Phonon Interaction neutron spectra\n");
printf("# * Author: Martin Rotter %s\n",MCPHASVERSION);
printf("# **************************************************************\n");
for (i=1;i<argc;++i)
 {if(strncmp(argv[i],"-h",2)==0) {helpexit();}
 else{T=strtod(argv[i],NULL);++i; if(T<=0){fprintf(stderr,"Error: T<=0\n");exit(EXIT_FAILURE);}
                                 // now read T
      Emin=strtod(argv[i],NULL);++i;  // now read Emin
      Emax=strtod(argv[i],NULL);++i;  // now read Emax
      deltaE=strtod(argv[i],NULL);++i;   // now read deltaE
      eta=strtod(argv[i],NULL);++i;  // now read eta
     }
 }
// reading results/op.mat
// read dimension
char instr[MAXNOFCHARINLINE];
 fin_opmat = fopen_errchk ("./results/op.mat", "rb");
  fgets_errchk (instr, MAXNOFCHARINLINE, fin_opmat);
  extract(instr,"d",d);
fclose(fin_opmat);
if(d<1){  fprintf (stderr,"#reading ./results/op.mat - Error: dimension d=%i\n",d);}
 fin_opmat = fopen_errchk ("./results/op.mat", "rb");

for(imax=0;feof(fin_opmat)==false;++imax)
{opmatM[imax] = new ComplexMatrix(1,d,1,d);
// fprintf(stderr,"%i",myReadComplexMatrix(fin_opmat, (*opmatM[imax])));
myReadComplexMatrix(fin_opmat, (*opmatM[imax]));
}
--imax;
fprintf(stdout,"# Ok - Read Hamiltonian and i=%i Olm Operator Matrices with dimension d=%i from ./results/op.mat \n",imax,d);
fclose(fin_opmat);

/* commented out because we will use phonon.int  ... no eigenvectors, but intensities as input
// read now the phonon Frequencies and Eigenvectors from  ./phonon.ev
 float nn[1000];nn[0]=1000;
fin_opmat = fopen_errchk ("./phonon.ev", "rb");
while(feof(fin_opmat)==false)
if((n=inputline(fin_opmat,nn))!=0)
{ // store phonon information 
 // In the first column the frequency ɷ(k, j) of the mode is given. 
 // Columns from 2 to N + 1 contain components of the real parts of 
 // the polarization vectors Ree(k, j), ordered Ax, Ay, Az, Bx, By, etc.
 // for Cartesian components x, y, z and particles A, B, . . ., and olumns 
 // from N + 2 to 2N + 1 contains components of the imaginary parts of 
  // the polarization vectors Ime(k, j) ordered in the same way.
  N=(int)((n-1)/6); // number of atoms
  omegaphon[smax]=nn[1]; // energy of mode
  eq[smax]= new ComplexMatrix(1,3,1,N);
  for(i=1;i<=3;++i)for(alp=1;alp<=N;++alp){
  (*eq[smax])(i,alp)=complex <double> (nn[3*(alp-1)+i+1],nn[3*N+3*(alp-1)+i+1]);
  }
 // myPrintComplexMatrix(stdout,(*eq[smax]));
  ++smax;
}
--smax;
fclose(fin_opmat);
fprintf(stdout,"# Read smax=%i phonons ./phonon.ev \n",smax+1);
*/

 float nn[1000];nn[0]=1000;
fin_opmat = fopen_errchk ("./phonon.int", "rb");
while(feof(fin_opmat)==false)
if((n=inputline(fin_opmat,nn))!=0)
{ // store phonon information 
 // In the first column the frequency ɷ(k, j) of the mode is given. 
 // Column  2 contains the corresponding intensity of the mode.
  omegaphon[smax]=nn[1]; // energy of mode
  intensity[smax]=nn[2]; // intensity of mode
  ++smax;
}
--smax;
fclose(fin_opmat);
fprintf(stdout,"# Read smax=%i phonons ./phonon.int \n",smax+1);


// read coupling
float galphas[1000]; int alpha[1000],s[1000],Ng=-1;
fin_opmat = fopen_errchk ("./g_alpha.s", "rb");
while(feof(fin_opmat)==false)
if((n=inputline(fin_opmat,nn))>=3)
{ // g_alpha(s)
 ++Ng;if(Ng>1000){fprintf (stderr,"#Error:too many lines in ./g_alpha.s\n");exit (EXIT_FAILURE);}
 alpha[Ng]=nn[1];if(alpha[Ng]>imax){fprintf (stderr,"#Error:too many lines in ./g_alpha.s\n");exit (EXIT_FAILURE);}
 s[Ng]=nn[2];if(s[Ng]>smax+1){fprintf (stderr,"#Error: s = %i > smax = %i in ./g_alpha.s\n",s[Ng],smax+1);exit (EXIT_FAILURE);}

 galphas[Ng]=nn[3]; 
  
}
fclose(fin_opmat);
fprintf(stdout,"# Read Ng+1=%i CF-Phonon coupling constants g_alpha(s) from ./g_alpha.s \n",Ng+1);

// do the calculations
double OOcef,rr,quot,En,Em,beta,Z,omega;
complex <double> sum,Mqs,E0,Dqs,omegeta;
beta=1/(K_B*T);
Z=0;E0=(*opmatM[0])(1,1);
for(n=1;n<=d;++n)
{(*opmatM[0])(n,n)-=E0;
 En=real((*opmatM[0])(n,n));
Z+=exp(-En*beta);
}
double p[1000];
for(n=1;n<=d;++n){En=real((*opmatM[0])(n,n));p[n]=exp(-En*beta)/Z;}
fprintf (stdout, "# omega[meV]   Iphon=dyn structfact for coherent neutron scattering x Imag(Dqs) [1/meV] for s=1,2,3,...\n");
double dsigma=0;

for(omega=Emin;omega<=Emax;omega+=deltaE)
{fprintf (stdout,"%+9.6f ",omega);
 dsigma=0;
  for(int S=0;S<=smax;++S)
{// evaluate equation (80/81)
 Mqs=0;
 for(k=0;k<=Ng;++k)if(S==s[k]-1)
 {
 OOcef=0;
 for(n=1;n<=d;++n)for(m=1;m<=d;++m)
 {rr=abs((*opmatM[alpha[k]])(n,m));
 En=real((*opmatM[0])(n,n));
 Em=real((*opmatM[0])(m,m));
 if(fabs(En-Em)>0.0000001){quot=(p[m]-p[n])/(En-Em);}
 else {quot=p[m]*beta;quot=0;} // quasielastic cf transitions excluded !!!
 OOcef+=rr*rr*quot;
 OOcef/=beta;
 sum=complex <double> (En-Em-omega,-eta);
 sum=OOcef/sum;
 Mqs+=2.0*beta*omegaphon[s[k]-1]*fabs(galphas[k])*fabs(galphas[k])*sum;
 }}
// equation (71)
omegeta=complex <double> (omega,eta);
Dqs=-2.0*omegaphon[S]/(omegeta*omegeta-omegaphon[S]*omegaphon[S]+omegeta*Mqs);

// evaluate (68)
// not possible because only one Mass ... only print Imag (Dqs)
// ... yet we can use the intensity of phonons from input phonon.int to evaluate (68)
// note: delta(x)=1/pi Im (1/(x-i eps) ) = 1/pi  eps/(eps^2+x^2) 
// https://en.wikipedia.org/wiki/Dirac_delta_function
//
// note already from dimension it is clear, that Dqs corresponds to the 
// sum of delta functions - it follows that:
dsigma+=imag(Dqs)*intensity[S];  
 
}
fprintf (stdout,"%+9.6f ",dsigma);
//fprintf (stdout,"%+9.6f + i %+9.6f ",real(Mqs),imag(Mqs));
fprintf(stdout,"\n");

}


// free memory
for(d=0;d<=imax;++d)delete opmatM[d];
//for(d=0;d<=smax;++d) delete eq[d];
}