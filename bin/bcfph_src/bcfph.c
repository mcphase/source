/**************************************************************
 * bfkq.c - calculate neutron spectra with CF-Phonon interaction
 * Author: Martin Rotter
 **************************************************************/


#include "../../version"
#include "martin.h"
#include "myev.h"


#if defined(__linux__)
#include <sys/sysinfo.h>
#elif defined(__FreeBSD__) || defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <windows.h>
#endif


#define K_B  0.0862
#define PI   3.141592654
#define TRMAX 10000 /* maximum number of crystal field transitions */
/**********************************************************************/
void helpexit()
{ printf (" program bcfph \n "
          " use as:  bcfph T Emin Emax deltaE eta \n"
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


#ifdef _THREADS
#if defined  (__linux__) || defined (__APPLE__)
#include <pthread.h>
#define MUTEX_LOCK     pthread_mutex_lock
#define MUTEX_UNLOCK   pthread_mutex_unlock
#define MUTEX_TYPE     pthread_mutex_t
#define MUTEX_INIT(mm)  pthread_mutex_init (&mm, NULL)
#define EVENT_TYPE     pthread_cond_t
#define EVENT_INIT(e)  pthread_cond_init (&e, NULL)    
#define EVENT_SIG(e)   pthread_cond_signal (&e)
#define THRLC_TYPE     pthread_key_t
#define THRLC_INIT(k)  pthread_key_create(&k, dataDestructor)
#define THRLC_FREE(k)  pthread_key_delete(k)
#define THRLC_SET(k,v) pthread_setspecific (k,v)
#define THRLC_GET(v)   pthread_getspecific (v)
#define THRLC_GET_FAIL NULL
void dataDestructor(void *data) { }
#else
#include <windows.h>
#define MUTEX_LOCK     EnterCriticalSection
#define MUTEX_UNLOCK   LeaveCriticalSection
#define MUTEX_TYPE     CRITICAL_SECTION
#define MUTEX_INIT(mm)  InitializeCriticalSection (&mm)
#define EVENT_TYPE     HANDLE
#define EVENT_INIT(e)  e = CreateEvent (NULL, TRUE, FALSE, NULL)
#define EVENT_SIG(e)   SetEvent(e)
#define THRLC_TYPE     DWORD
#define THRLC_INIT(k)  k = TlsAlloc()
#define THRLC_FREE(k)  TlsFree(k)
#define THRLC_SET(k,v) TlsSetValue (k,v)
#define THRLC_GET(v)   TlsGetValue (v)
#define THRLC_GET_FAIL 0
#endif
#define NUM_THREADS nofthreads
// ----------------------------------------------------------------------------------- //
// Declares a struct to store all the information needed for each omcalc iteration
// ----------------------------------------------------------------------------------- //
typedef struct{
  double eta;
   double * omegamu;
   int  mumax;
   int  d;
   ComplexMatrix ** opmatM;
   int  N;
   int  imax; 
   int Ng;
   float * omegaks;
   double  beta;
   double * p;
   float * galphasr; 
   float * galphasi; 
   int * alpha;
   ComplexMatrix **P;
   int thread_id;
} omcalc_thread_data;
class omcalc_input { public:
   int j; 
   int thread_id;
   int level;
   double dsigma;
   double omega;
   omcalc_input(int _j, int _tid, double _dsigma, int _level, double _omega) 
   { 
      thread_id = _tid; j = _j; dsigma = _dsigma; level= _level; omega= _omega;
   }
   ~omcalc_input(){}
};
// ----------------------------------------------------------------------------------- //
// Declares these variables global, so all threads can see them
// ----------------------------------------------------------------------------------- //
omcalc_thread_data thrdat;
omcalc_input *tin[256];  // Max number of threads - hard coded because global variable.
MUTEX_TYPE mutex_loop;
MUTEX_TYPE mutex_tests;
MUTEX_TYPE mutex_min;
MUTEX_TYPE mutex_index;
EVENT_TYPE checkfinish;
THRLC_TYPE threadSpecificKey;

#endif // def _THREADS


#ifdef _THREADS
#if defined  (__linux__) || defined (__APPLE__)
void *omcalc(void *input)
#else
DWORD WINAPI omcalc(void *input)
#endif
#else
double omcalc(double & dsigma, double & omega, double & eta,double * omegamu,int & mumax,int & d,
           ComplexMatrix ** opmatM,int& N,int & imax, int&Ng,float * omegaks,double & beta,double * p,
           float * galphasr, float * galphasi, int * alpha,ComplexMatrix **P)
#endif
{


#ifdef _THREADS
   omcalc_input *myinput; myinput = (omcalc_input *)input;
   int thread_id = myinput->thread_id;
   double omega = myinput->omega;
   #define eta thrdat.eta
   #define omegamu thrdat.omegamu
   #define mumax thrdat.mumax
   #define d thrdat.d
   #define opmatM thrdat.opmatM
   #define N thrdat.N
   #define imax thrdat.imax
   #define Ng thrdat.Ng
   #define omegaks thrdat.omegaks
   #define beta thrdat.beta
   #define p thrdat.p
   #define galphasr thrdat.galphasr
   #define galphasi thrdat.galphasi
   #define alpha thrdat.alpha
   #define P thrdat.P
   double dsigma;
#endif


 ComplexMatrix chi(1,3,1,3),OM(1,3,1,3);
 complex <double> M,gbetaks,bmudnunm[52],A,B,factor; 
 complex <double> omegeta;
 double nks,omeganm,omegamud,Ems,bo,omegaksold=0;
 double En,Em;
 chi=0;omegeta=complex <double> (omega,eta);
 for(int nu=1;nu<=mumax;++nu){ // sum  in (27) is sufficient over nu (because OM is diagonal, i.e. prop delta_munu)
 OM=0;
for(int alp=1;alp<=3;++alp){
 
 OM(alp,alp)=complex <double> (-omega+omegamu[nu],-eta);
  // evaluate (47) to get M(omega) ---------------------------------------------------------------
  M=0;

     for(int n=1;n<=d;++n)for(int m=1;m<=d;++m){ // sum over nm
           En=real((*opmatM[0])(n,n));
           Em=real((*opmatM[0])(m,m));
           omeganm=En-Em;
           omegamud=omeganm-omegamu[nu]; 
                    for(int bet=1;bet<=imax;++bet){ //++++ calulate b^betalp_mudnunm (42)(39) +++++++++++++++++++++++++++++
                      bmudnunm[bet]=0;
                     for(int ms=1;ms<=d;++ms){Ems=real((*opmatM[0])(ms,ms));
                                          if(fabs(Ems-Em-omegamu[nu])<0.00001&&fabs(En-Ems-omegamud)<0.00001)bmudnunm[bet]+=(*opmatM[bet])(n,ms)*(*opmatM[alp])(ms,m);
                                          if(fabs(En-Ems-omegamu[nu])<0.00001&&fabs(Ems-En-omegamud)<0.00001)bmudnunm[bet]-=(*opmatM[bet])(ms,m)*(*opmatM[alp])(n,ms);
                                         }
                                        }
                                       //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // sum all gbeta(s) in (48) ||||||||||||||||||||||||||||||||
                omegaksold=0;
                    for(int i=1;i<=Ng;++i)if((bo=fabs(beta*omegaks[i]))>0){ // sum over "ksbeta"=i where beta=alpha[i], s=s[i], k not needed, but omegaks[i]
                          
                if(fabs(omegaks[i]-omegaksold)>0.00001){ // calculate factor only if omegaks[i] changed since last time
                          if(bo>0.001)nks=1/(exp(bo)-1);
                          else nks=1/bo;
                
                A= complex <double> (omeganm-omegaks[i],-eta); // to avoid divergencies also add finite width to A and B
                A=(nks*p[m]-(1+nks)*p[n])/beta/A;
                B=complex <double> (omeganm+omegaks[i],-eta);
                B=((1+nks)*p[m]-nks*p[n])/beta/B;
                factor=(A/(omeganm-omegaks[i]-omegeta)+B/(omeganm-omegaks[i]-omegeta));
                omegaksold=omegaks[i];
                                                                }
                gbetaks=complex <double> (galphasr[i],galphasi[i]);
                M+=abs(bmudnunm[alpha[i]]*gbetaks)*factor;
                 //fprintf(stderr,"%g %g %g %g %g \n",A,B,abs(bmudnunm*gbetaks),real(M),imag(M));
                //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                                                                      }
                    }
  M/=N;
  // here should be subtracte from O the interation term M(omega)/P equation (28) ----------------
//fprintf(stderr,"%g %g %g %g  %g %g \n",real(M),imag(M),real(OM(alp,alp)),imag(OM(alp,alp)),real((*P[nu])(alp,alp)),imag((*P[nu])(alp,alp))); 
 if(abs(M)>0.000001)OM(alp,alp)=OM(alp,alp)-M/(*P[nu])(alp,alp);  
    }
// myPrintComplexMatrix(stderr,OM);
 
   chi+=(*P[nu])*OM.Inverse(); //equ (27)
   } 
dsigma=0.66666*imag(chi(1,1)+chi(2,2)+chi(3,3));  //equ(18)
if(fabs(omega*beta)>0.001){dsigma*=(beta*omega)/(1-exp(-beta*omega));}
//fprintf(stdout,"# omega=%g dsigma=%g\n",omega,dsigma);

#ifdef _THREADS
myinput->dsigma=dsigma;
#undef eta
#undef omegamu
#undef mumax
#undef d
#undef opmatM
#undef N
#undef imax
#undef Ng
#undef omegaks
#undef beta
#undef p
#undef galphasr
#undef galphasi
#undef alpha
#undef P
MUTEX_LOCK(&mutex_loop);
thrdat.thread_id = thread_id;
EVENT_SIG(checkfinish);
MUTEX_UNLOCK(&mutex_loop);
#if defined  (__linux__) || defined (__APPLE__)
pthread_exit(NULL);
#else
return 0;
#endif
#else
return dsigma;	
#endif
}

// hauptprogramm
int main (int argc, char **argv)
{
#ifdef _THREADS
int nofthreads=1;
#if defined(__linux__)                               // System-dependent calls to find number of processors (from GotoBLAS)
       nofthreads = get_nprocs();
#elif defined(__FreeBSD__) || defined(__APPLE__)
       int mm[2]; size_t len;
       mm[0] = CTL_HW; mm[1] = HW_NCPU; len = sizeof(int);
       sysctl(mm, 2, &nofthreads, &len, NULL, 0);
#else
       SYSTEM_INFO sysinfo; GetSystemInfo(&sysinfo);
       nofthreads = sysinfo.dwNumberOfProcessors;
#endif
if(nofthreads<1||nofthreads>255) nofthreads=1; 
int num_threads_started=-1;
#endif



int i,imax,k,d=0,n,m,smax=0; double T,Emin, Emax, deltaE, eta;
 FILE * fin_opmat;
 ComplexMatrix *opmatM[52]; double omegaphon[1000],intensity[1000];
printf("#**************************************************************\n");
printf("# * bfkq.c - calculate CF-Phonon Interaction neutron spectra\n");
printf("# * Author: Martin Rotter %s\n",MCPHASVERSION);
printf("# **************************************************************\n");
printf("# Command: bcfph ");
for (i=1;i<argc;++i)printf(" %s",argv[i]);
printf("\n");
if(argc<4) {helpexit();}
for (i=1;i<argc;++i)
 {
 if(strncmp(argv[i],"-h",2)==0) {helpexit();}
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
float * galphas,*galphasr,*galphasi,*omegaks;
int Ng=-1,N=0,GMAX=1; 
fprintf(stdout,"# reading coupling from g_alpha.s ...\n");
fin_opmat = fopen_errchk ("./g_alpha.s", "rb");
while(feof(fin_opmat)==false){fgets(instr,MAXNOFCHARINLINE,fin_opmat);++GMAX;}
fclose(fin_opmat);
fprintf(stdout,"# ...found %i lines - reserving memory for storage of coupling parameters\n",GMAX); 
galphas = new float[GMAX];
galphasr = new float[GMAX];
galphasi = new float[GMAX];
omegaks = new float[GMAX];
 
int * s; s=new int[GMAX];
int * alpha; alpha=new int[GMAX];
long int pos=0,posold=1;
fprintf(stdout,"# storing coupling paramters\n");
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
      fgets(instr,MAXNOFCHARINLINE,fin_opmat);
      instr[0]=' ';
      extract(instr,"N",N);
                     }
}

}
fclose(fin_opmat);
if(N==0)N=1;
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
fprintf (stdout, "# Calculating CF transition energies and matrix elements  gamma (see manual, output of singleion *.trs) \n"
                 "# E[meV]   vs  gamma \n");

// sort transitions into groups belonging to the same transition energy
double fr,dsigma,omegamu[TRMAX];int mu,nu,mumax=0,ok,alp,bet;
ComplexMatrix *P[TRMAX];
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
  if(ok==0){++mumax;if(mumax>TRMAX-1){fprintf (stderr,"#Error:too many transitions - increase TRMAX and recompile\n");exit (EXIT_FAILURE);}
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
fprintf (stdout, "# Note: Imagpolycrystal [barn/srmeV]= (r0 gJ/2)^2 1/pi  (1-exp(-hbar omega/kT))^(-1) 2/3Trace{Im(chi(omega))} \n"
                 "#        with  r0 =  -0.53908 * 10^-12 cm  , 1 barn=10^-24 cm \n"
                 "# Here follows output:\n"
                 "# omega[meV]  vs   (1-exp(-hbar omega/kT))^(-1) 2/3Trace{Im(chi(omega))}  [1/meV]\n");

#ifdef _THREADS  // Populates the thread data structure
   thrdat.eta=eta;
   thrdat.omegamu=omegamu;
   thrdat.mumax=mumax;
   thrdat.d=d;
   thrdat.opmatM=opmatM;
   thrdat.N=N;
   thrdat.imax=imax;
   thrdat.Ng=Ng;
   thrdat.omegaks=omegaks;
   thrdat.beta=beta;
   thrdat.p=p;
   thrdat.galphasr=galphasr;
   thrdat.galphasi=galphasi;
   thrdat.alpha=alpha;
   thrdat.P=P;
   // Initialises mutual exclusions and threads
   MUTEX_INIT(mutex_loop);
   MUTEX_INIT(mutex_index);
   EVENT_INIT(checkfinish);
   #if defined  (__linux__) || defined (__APPLE__)
   pthread_t threads[NUM_THREADS]; int rc; void *status;
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   #else
   HANDLE threads[NUM_THREADS];
   DWORD tid[NUM_THREADS], dwError;
   int retval;
   #endif
   int ithread;
   omcalc_input *tin[NUM_THREADS];     
   for (ithread=0; ithread<NUM_THREADS; ithread++) 
   { tin[ithread] = new omcalc_input(ithread,1,0.0,0,0);    
   } 
    ithread=0; num_threads_started=0; int oldi=-1;
#endif
                     int ct,ctmax=(int)(Emax-Emin)/deltaE;
double * om; om = new double [ctmax+1];
double * ds; ds = new double [ctmax+1]; 

#ifdef _THREADS  
                  for (ct=0;ct<=ctmax;ct+=NUM_THREADS)
#else
                  
                  for (ct=0;ct<=ctmax;++ct)
#endif
                  {
#ifdef _THREADS  
                     oldi=ct;
                     // Runs threads until all are running - but wait until they are completed in order before printing output.
                     //  This is to ensure that the output is exactly the same as for a single thread (otherwise it would be out of order).
                     for(int th=0; th<NUM_THREADS; th++)
                     { 
                        ct=oldi+th; if(ct>ctmax) break;
                       //printf("calling thread %i ",num_threads_started);
                       tin[num_threads_started]->dsigma=ds[ct]; 
                       tin[num_threads_started]->level = ct;
                       tin[num_threads_started]->omega = Emin+ct*deltaE;

                        #if defined  (__linux__) || defined (__APPLE__)
                        rc = pthread_create(&threads[num_threads_started], &attr, omcalc, (void *) tin[num_threads_started]);
                        if(rc) { printf("Error return code %i from omcalc thread %i\n",rc,num_threads_started+1); exit(EXIT_FAILURE); }
                        #else
                        threads[num_threads_started] = CreateThread(NULL, 0, omcalc, (void *) tin[num_threads_started], 0, &tid[num_threads_started]);
                        if(threads[num_threads_started]==NULL) { dwError=GetLastError(); printf("Error code %lu from omcalc thread %i\n",dwError,num_threads_started+1); exit(EXIT_FAILURE); }
                        #endif  
                        ++num_threads_started;
                       }
                     #if defined  (__linux__) || defined (__APPLE__)
                     for(int th=0; th<num_threads_started; th++)rc = pthread_join(threads[th], &status);
                     #else
                     if(num_threads_started>0){retval=WaitForMultipleObjects(num_threads_started,threads,TRUE,INFINITE);
                     if(retval<WAIT_OBJECT_0||retval>WAIT_OBJECT_0+num_threads_started-1){printf("Error waitformultipleobjects=%i num_threads_started=%i\n",retval,num_threads_started); exit(EXIT_FAILURE); }
                      for(int th=0; th<num_threads_started; th++)CloseHandle(threads[th]);}         
                     #endif
                     num_threads_started=0; 
                       ithread=-1;
                     for(int th=0; th<NUM_THREADS; th++)
                     {
                         ct=oldi+th; if(ct>ctmax) break;
                        ++ithread;
                       ds[tin[ithread]->level] = tin[ithread]->dsigma; 
                       om[tin[ithread]->level] = tin[ithread]->omega; 
                       printf("%i %i %g %g\n",ithread,tin[ithread]->level,om[tin[ithread]->level],ds[tin[ithread]->level]);
                            
#else
                     omega=Emin+ct*deltaE;om[ct]=omega;

 omcalc( dsigma,  omega, eta,omegamu,mumax,d,opmatM,N,imax,Ng,omegaks,beta,p,galphasr,galphasi,alpha,P);

ds[ct]=dsigma;
#endif


#ifdef _THREADS
                     }
                   i=oldi;
#endif
		   }

// print results
for (ct=0;ct<ctmax;++ct){ fprintf (stdout,"%+9.6f %+9.6f\n",om[ct],ds[ct]);}


// print transition energies and intensities
fprintf(stdout,"# CF Transition Matrix Elements\n");
fprintf(stdout,"# E[meV]  gamma\n");
 for(nu=1;nu<=mumax;++nu){
gamma=real((*P[nu])(1,1)+(*P[nu])(2,2)+(*P[nu])(3,3));
dsigma=0.66666*gamma;
if(fabs(omegamu[nu])<0.0001){gamma*=beta;}else{gamma*=beta*omegamu[nu];}
if(fabs(omegamu[nu]*beta)>0.001){dsigma*=(beta*omegamu[nu])/(1-exp(-beta*omegamu[nu]));}
fprintf(stdout,"# %+9.6f %+9.6f\n",omegamu[nu],gamma);

} 
#ifdef _THREADS
std::cout << "#! nofthreads= " << NUM_THREADS << " threads were used in parallel processing " << std::endl;
for (int ithread=0; ithread<NUM_THREADS; ithread++) delete tin[ithread];
#else
std::cout << "# bcfph was compiled without parallel processing option " << std::endl;
#endif
//myPrintComplexMatrix(stdout,(*P[nu]));

// free memory
for(mu=0;mu<=mumax;++mu)delete P[mu];
delete om;
delete ds;
} else { //**********************************************************************************
// do the calculations of nuclear coherent phonon cross section 
//**********************************************************************************
if(N>1){fprintf(stderr,"Error bcfph: N=%i k points found in g_alpha.s - please use g_alpha.s with only one k point\n"
                       " - you can generate this by qep2bcfph.pl or makegalphas.pl\n",N);exit(EXIT_FAILURE);}
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
delete alpha;
delete galphas;
delete galphasr;
delete galphasi;
delete omegaks;

for(d=0;d<=imax;++d)delete opmatM[d];
//for(d=0;d<=smax;++d) delete eq[d];
}