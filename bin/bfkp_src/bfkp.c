/**************************************************************
 * bfkq.c - calculate neutron spectra with CF-Phonon interaction
 * Author: Martin Rotter
 **************************************************************/


#include "../../version"
#include "martin.h"
#include "myev.h"

#define K_B  0.0862
#define PI   3.141592654
#define GMAX 10000
/**********************************************************************/
void helpexit()
{ printf (" program bfkp \n "
          " use as:  bfkp T Emin Emax deltaE eta \n"
          "   T    .... Temperature (K)\n"
          "   Emin .... minimal Energy (meV)\n "
          "   Emax .... maximal Energy (meV)\n "
          "   deltaE .... Energy stepwidth (meV)\n "
          "   eta .... imaginary small number for linewidth (meV)\n "   
          "   required: ./results/op.mat for Operator matrices \n  "
          "             ./g_alpha.s for cf-phonon coupling constants\n"
          "  T>0         ./phonon.int for Frequencies and Intensities of Phonon\n"
          "  Output to stdout: nuclear coherent inelastic spectrum (equ. 68 ff) \n"
          "  T<0  \n "       
          "  Output to stdout: magnetic (crystal-field) spectrum (equ. 18 ff), with formfactor F(Q)=1 \n"
          
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
if(argc<4) {helpexit();}
for (i=1;i<argc;++i)
 {if(strncmp(argv[i],"-h",2)==0) {helpexit();}
 else{T=strtod(argv[i],NULL);++i; 
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


// read coupling
float nn[1000];nn[0]=1000;
float galphas[GMAX],galphasr[GMAX],galphasi[GMAX],omegaks[GMAX]; 
int alpha[GMAX],s[GMAX],Ng=-1,N;
long int pos=0,posold=1;
fin_opmat = fopen_errchk ("./g_alpha.s", "rb");
while(feof(fin_opmat)==false){
if((n=inputline(fin_opmat,nn))>=3)
{ // g_alpha(s)
 ++Ng;if(Ng>GMAX-1){fprintf (stderr,"#Error:too many lines in ./g_alpha.s\n");exit (EXIT_FAILURE);}
 alpha[Ng]=nn[1];if(alpha[Ng]>imax){fprintf (stderr,"#Error: alpha=%i too large  in ./g_alpha.s - file opmat does not contains so many operators\n",alpha[Ng]);exit (EXIT_FAILURE);}
 s[Ng]=nn[2]; galphas[Ng]=nn[3]; 
 galphasr[Ng]=nn[4]; 
 galphasi[Ng]=nn[5]; 
 omegaks[Ng]=nn[9]; posold=pos;
 pos=ftell(fin_opmat);
}else{if(pos!=posold){posold=pos;
      fseek(fin_opmat,pos,SEEK_SET);
      fgets_errchk(instr,MAXNOFCHARINLINE,fin_opmat);
      instr[0]=' ';
      extract(instr,"N",N);
   } 
}

}
fclose(fin_opmat);
fprintf(stdout,"# Read Ng+1=%i CF-Phonon coupling constants g_alpha(s) from ./g_alpha.s number of k points N=%i \n",Ng+1,N);

double OOcef,rr,quot,En,Em,beta,Z,omega;
complex <double> sum,Mqs,E0,Dqs,omegeta;
beta=fabs(1/(K_B*T));
Z=0;E0=(*opmatM[0])(1,1);
for(n=1;n<=d;++n)
{(*opmatM[0])(n,n)-=E0;
 En=real((*opmatM[0])(n,n));
Z+=exp(-En*beta);
}
double p[1000],gamma;
for(n=1;n<=d;++n){En=real((*opmatM[0])(n,n));p[n]=exp(-En*beta)/Z;}

if(T<=0){//**********************************************************************************
// do the calculations of magnetic crystal field cross section 
//**********************************************************************************
T=-T;
fprintf (stdout, "# Note: Imagpolycrystal [barn/srmeV]= (r0 gJ/2)^2 1/pi  (1-exp(-hbar omega/kT))^(-1) 2/3Trace{Im(chi(omega))} \n "
                 "#        with  r0 =  -0.53908 * 10^-12 cm  , 1 barn=10^-24 cm \n"
                 "# Here follows output:\n"
                 "# omega[meV]  vs   (1-exp(-hbar omega/kT))^(-1) 2/3Trace{Im(chi(omega))}  [1/meV]\n");

// sort transitions into groups belonging to the same transition energy
double fr,dsigma,omegamu[GMAX];int mu,nu,mumax=0,ok,alp,bet;
ComplexMatrix chi(1,3,1,3),OM(1,3,1,3);
ComplexMatrix *P[GMAX];
for(n=1;n<=d;++n)for(m=1;m<=d;++m)
 {En=real((*opmatM[0])(n,n));
  Em=real((*opmatM[0])(m,m));
  ok=0;
  for(mu=1;mu<=mumax;++mu)if(fabs(omegamu[mu]-En+Em)<0.0001){ok=1;
            if (fabs(omegamu[mu])>0.0001){ fr=(p[m]-p[n])/(beta*omegamu[mu]);}
            else {fr=p[m];} // approximate fr for quasielastic line
            for(alp=1;alp<=3;++alp)for(bet=1;bet<=3;++bet){(*P[mu])(alp,bet)+=fr*((*opmatM[alp])(n,m))*((*opmatM[bet])(m,n));}           
             gamma=real((*P[mu])(1,1)+(*P[mu])(2,2)+(*P[mu])(3,3));
            if(fabs(omegamu[mu])<0.0001){gamma*=beta;}else{gamma*=beta*omegamu[mu];}
            fprintf(stdout,"# %+9.6f %+9.6f \n",omegamu[mu],gamma);
             }
  if(ok==0){++mumax;if(mumax>GMAX-1){fprintf (stderr,"#Error:too many transitions - increase GMAX and recompile\n");exit (EXIT_FAILURE);}
            omegamu[mumax]=En-Em;
            P[mumax]= new ComplexMatrix(1,3,1,3);
            if (fabs(omegamu[mumax])>0.0001){ fr=(p[m]-p[n])/(beta*omegamu[mumax]);}
            else {fr=p[m];} // approximate fr for quasielastic line
            for(alp=1;alp<=3;++alp)for(bet=1;bet<=3;++bet){(*P[mumax])(alp,bet)=fr*(*opmatM[alp])(n,m)*(*opmatM[bet])(m,n);}
            gamma=real((*P[mumax])(1,1)+(*P[mumax])(2,2)+(*P[mumax])(3,3));
            if(fabs(omegamu[mumax])<0.0001){gamma*=beta;}else{gamma*=beta*omegamu[mumax];}
            fprintf(stdout,"# %+9.6f %+9.6f \n",omegamu[mumax],gamma);

            }
 }


for(omega=Emin;omega<=Emax;omega+=deltaE)
{fprintf (stdout,"%+9.6f ",omega);
 chi=0;//omegeta=complex <double> (omega,eta);
 for(nu=1;nu<=mumax;++nu){ // sum  in (27) is sufficient over nu (because OM is diagonal, i.e. prop delta_munu)
 OM=0;
for(alp=1;alp<=3;++alp){
 
 OM(alp,alp)=complex <double> (-omega+omegamu[nu],-eta);
           

  // here should be added to O the interation term M(omega)/P equation (28)
  //  ....
    }
 
   chi+=(*P[nu])*OM.Inverse(); //equ (27)
   } 

dsigma=0.66666*imag(chi(1,1)+chi(2,2)+chi(3,3));  //equ(18)
if(fabs(omega*beta)>0.001){dsigma*=(beta*omega)/(1-exp(-beta*omega));}

fprintf (stdout,"%+9.6f ",dsigma);
fprintf(stdout,"\n");

}

// print transition energies and intensities
fprintf(stdout,"# Transition Intensities without CF-Phonon Interaction\n");
fprintf(stdout,"# E[meV]  gamma  I[barn/sr]\n");
 for(nu=1;nu<=mumax;++nu){
gamma=real((*P[nu])(1,1)+(*P[nu])(2,2)+(*P[nu])(3,3));
dsigma=0.66666*gamma;
if(fabs(omegamu[nu])<0.0001){gamma*=beta;}else{gamma*=beta*omegamu[nu];}
if(fabs(omegamu[nu]*beta)>0.001){dsigma*=(beta*omegamu[nu])/(1-exp(-beta*omegamu[nu]));}
fprintf(stdout,"# %+9.6f %+9.6f %+9.6f \n",omegamu[nu],gamma, 0.8571*0.8571*0.54*0.54/4*dsigma);
}

//myPrintComplexMatrix(stdout,(*P[nu]));

// free memory
for(mu=0;mu<=mumax;++mu)delete P[mu];

} else { //**********************************************************************************
// do the calculations of nuclear coherent phonon cross section 
//**********************************************************************************

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
for(k=0;k<=Ng;++k)if(s[k]>smax+1){fprintf (stderr,"#Error: s = %i > smax = %i in ./g_alpha.s\n",s[k],smax+1);exit (EXIT_FAILURE);}



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
 if(fabs(En-Em)>0.00001){quot=(p[m]-p[n])/(En-Em);}
 else {quot=p[m]*beta;} 
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
// note already from dimension it is clear, that Dqs/PI corresponds to the 
// sum of delta functions - it follows that:
dsigma+=imag(Dqs)*intensity[S]/PI;  
 
}
fprintf (stdout,"%+9.6f ",dsigma);
//fprintf (stdout,"%+9.6f + i %+9.6f ",real(Mqs),imag(Mqs));
fprintf(stdout,"\n");

}
}

// free memory
for(d=0;d<=imax;++d)delete opmatM[d];
//for(d=0;d<=smax;++d) delete eq[d];
}