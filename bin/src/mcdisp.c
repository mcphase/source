/***********************************************************************
 *
 * mcdisp - program to calculate the dispersion of magnetic excitations
 *
 * reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751
 *            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400
 ***********************************************************************/
  
 
#include <mcdisp.h> 
#include <ctime>

#include "par.hpp"
#include "mfcf.hpp"
#include "mdcf.hpp"
#include "jq.hpp"
#include "inimcdis.hpp"

#include <complex>
#include <martin.h>

#include "../../version"
#include "myev.c"

void errexit() // type info and error exit 
{     printf (" \n %s \n",MCDISPVERSION);
    printf ("use as: mcdisp\n"); 
    printf (" or as: mcdisp [options] [file]\n");
    printf ("  [file] ... input file with mean field set (default mcdisp.mf)\n");
    printf ("Options:\n");
    printf (" -jq           ... calculate J(Q) (Fourier transform of 2ion coupling) store in mcdisp.jq largest evalue and eigenvector\n");
    printf ("                   if energies are given for hkls in mcdisp.par, output file mcdisp_scaled.jq contains scaled parameters\n");
    printf ("                   such that energy of first hkl set corresponds to highest eigenvalue of J(Q)\n");
    printf (" -jqe          ... calculate J(Q) (Fourier transform of 2ion coupling) store in mcdisp.jq all eigenvalues \n");
    printf (" -max n        ... restrict single ion susceptibility to n lowest\n");
    printf ("                   lying transitions starting from the ground state\n");
    printf (" -minE E       ... an energy range may be given by minE and maxE: only\n");
    printf (" -maxE E           single ion transitions within this energy range will \n");
    printf ("                   be considered\n");
    printf (" -r            ... refine energies\n");
    printf (" -x            ... calculate resonant inelastic x-ray intensities (maximized with respect to azimuth) instead of neutron intensities\n");
    printf (" -xa   stp     ... calculate resonant inelastic x-ray intensities with complete azimuth dependence for each reflection (stp in deg)\n");
    printf (" -xaf  az      ... calculate resonant inelastic x-ray intensities at specified azimuth (deg) for each reflection\n");
    printf (" -d            ... calculate intensities in dipole approximation only\n");
    printf (" -v            ... verbose\n");
    printf (" -a            ... do not overwrite output files in results - append results\n");
    printf (" -A            ... do not overwrite output files - compare hkl's to be calculated with existing list in mcdisp.qom and\n");
    printf ("                   continue calculation at last matching q vector\n");
    printf (" -c            ... only create single ion transition file ./results/mcdisp.trs and exit\n");
    printf (" -t            ... read single ion transition file ./results/mcdisp.trs (do not create it)\n");
    printf (" -ninit n      ... maximum number n of (low energy) initial states (single ion transitions)\n");
    printf ("                   (not functional with all single ion modules)\n");
    printf (" -pinit p      ... minimum populationnumber p of initial state (single ion transitions)\n");
    printf ("                   in order to be considered (not functional with all single ion modules)\n");
    printf (" -prefix 001   ... prefix for parameters to be read from mcdisp.par and used for creation of output files\n"
            "                   (useful for running in parallel calculations for different zones: e.g. put in\n"
            "                   mcdisp.par instead of #!hklline= several statements #!001hklline= ... #!002hklline=\n"
            "                   and start several jobs of mcdisp with -prefix 001, -prefix 002 simultaneously, afterwards merge\n"
            "                   output files, e.g. *mcdisp.qei  with appendfile)\n");
    printf (" -ignore_non_hermitian_matrix_error   ... ignores error when energies get complex due to unphysical mf groundstate\n");
    printf ("\n");
    printf ("Note: files which must be in current directory -\n");
    printf ("      ./mcdisp.par, ./mcphas.j, directory ./results\n");
      exit (EXIT_FAILURE);
} 

int index_s(int i,int j,int k,int l, int t, const mdcf & md, const inimcdis & ini)
{int s=0,i1,j1,k1;
// calculates the index of the Matrix A given 
// ijk ... index of crystallographic unit cell in magnetic unit cell
// l   ... number of atom in crystallographic cell
// t   ... transitionnumber
 for(i1=1;i1<i;++i1){
 for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 s+=md.baseindex_max(i1,j1,k1);
 }}}
 i1=i;
 for(j1=1;j1<j;++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 s+=md.baseindex_max(i1,j1,k1);
 }}
 j1=j;
 for(k1=1;k1<k;++k1){
 s+=md.baseindex_max(i1,j1,k1);
 }
 s+=md.baseindex(i,j,k,l,t);
 return s;
}



#ifdef _THREADS
#define _THREADS_JSSS 1
#define _THREADSREFINE 1
// ----------------------------------------------------------------------------------- //
// Defines to ease interchange between linux and windows thread codes...
// ----------------------------------------------------------------------------------- //
#if defined  (__linux__) || defined (__APPLE__)
#include <pthread.h>
#define MUTEX_LOCK    pthread_mutex_lock
#define MUTEX_UNLOCK  pthread_mutex_unlock
#define MUTEX_TYPE    pthread_mutex_t
#define MUTEX_INIT(m) pthread_mutex_init (&m, NULL)
#define EVENT_TYPE    pthread_cond_t
#define EVENT_INIT(e) pthread_cond_init (&e, NULL)    
#define EVENT_SIG(e)  pthread_cond_signal (&e)
#else
#include <windows.h>
#define MUTEX_LOCK    EnterCriticalSection
#define MUTEX_UNLOCK  LeaveCriticalSection
#define MUTEX_TYPE    CRITICAL_SECTION
#define MUTEX_INIT(m) InitializeCriticalSection (&m)
#define EVENT_TYPE    HANDLE
#define EVENT_INIT(e) e = CreateEvent (NULL, TRUE, FALSE, NULL)
#define EVENT_SIG(e)  SetEvent(e)
#endif
#define NUM_THREADS ini.nofthreads
 
// ----------------------------------------------------------------------------------- //
// Declares a struct to store all the information needed for each disp_calc iteration
// ----------------------------------------------------------------------------------- //
typedef struct{
   ComplexMatrix **chi, **chibey, **chiPhon;
   ComplexMatrix **Echargedensity;mfcf  **qee_real, **qee_imag;
   ComplexMatrix **Espindensity;mfcf  **qsd_real, **qsd_imag;
   ComplexMatrix **Eorbmomdensity;mfcf  **qod_real, **qod_imag;
   ComplexMatrix **Ephonon;mfcf  **qep_real, **qep_imag;
   ComplexMatrix **Emagmom;mfcf  **qem_real, **qem_imag;
   ComplexMatrix **Espin;mfcf  **qes_real, **qes_imag;
   ComplexMatrix **Eorbmom;mfcf  **qel_real, **qel_imag;
   ComplexMatrix **ch;
   Matrix **pol;
   ComplexMatrix **Tau;
   Vector hkl, q;  jq **J;
   inimcdis **ini;
   par **inputpars;
   mdcf **md;
   int thread_id;
} intcalcapr_thread_data;
class intcalcapr_input { public:
   int thread_id;
   int dimA, level, do_verbose,calc_rixs,do_phonon;
   double  En,intensity, intensitybey, intensityP, QQ;
   double epsilon; int iE,Estp;
   intcalcapr_input(int _dimA, int _tid, int _level, int _doverb, int _calcrixs,int _do_phonon, double _En)
   { 
      thread_id = _tid; dimA = _dimA; level = _level; do_verbose = _doverb;calc_rixs= _calcrixs, do_phonon=_do_phonon; En = _En;
   }
};
// ----------------------------------------------------------------------------------- //
// Declares these variables global, so all threads can see them
// ----------------------------------------------------------------------------------- //
intcalcapr_thread_data thrdat;
MUTEX_TYPE mutex_loop;
MUTEX_TYPE mutex_index;
//MUTEX_TYPE mutex_Jlock_write0;
//MUTEX_TYPE mutex_Jlock_write1;
EVENT_TYPE JQfree;
EVENT_TYPE checkfinish;

// ----------------------------------------------------------------------------------- //
// Routine to calculate the matrix product and sum for threaded calculation of jsss
// ----------------------------------------------------------------------------------- //
#endif // if _THREADS

#include "mcdisp_intcalc.c"
#include "mcdisp_output.c"
#include "trs_io.c"   // for in out of trs file
 
#ifdef _THREADS_JSSS
#define inputpars (*thrdat.inputpars[thread_id])
#define ini (*thrdat.ini[thread_id])
#define md (*thrdat.md[thread_id])
#define J (*thrdat.J[thread_id])
#define q thrdat.q
#if defined  (__linux__) || defined (__APPLE__)
void *jsss_mult(void *input)
#else
DWORD WINAPI jsss_mult(void *input)
#endif
#else
void jsss_mult(int ll, long int &nofneighbours, Vector q,  par &inputpars, inimcdis &ini, jq &J, mdcf &md)
#endif
{
#ifdef _THREADS_JSSS
    intcalcapr_input *myinput; myinput = (intcalcapr_input *)input;
    int nofneighbours=myinput->dimA, ll=myinput->level;
    int thread_id = myinput->thread_id;
#endif
    complex<double> ipi(0,2*3.1415926535), expqd;
    int i,j,k,i1,j1,k1,s,ss,sl,tl,tll,m,n;
    int l;
    double REexpqd, IMexpqd, jjval; int jsi,jsj;
    for(l=1;l<=(*inputpars.jjj[ll]).paranz;++l) 
    {
         int sd=(*inputpars.jjj[ll]).sublattice[l];
         Vector d(1,3),d_rint(1,3),xyz(1,3),xyz_rint(1,3);// some vector
         Vector ij(1,3);
         xyz=(*inputpars.jjj[ll]).xyz+(*inputpars.jjj[ll]).dn[l]-(*inputpars.jjj[sd]).xyz; // line added 17.6.09 to remove rounding bug in PCSMO calculation
         d=inputpars.rez*(const Vector&)xyz;
         for (i=1;i<=3;++i)d_rint(i)=rint(d(i)); // rint d for loop below to determine crystallographic unit ss ...

         xyz=(*inputpars.jjj[ll]).dn[l];
         d=inputpars.rez*(const Vector&)xyz;// set d to distance for later use to determine phase factor in J(Q) ...

         expqd = exp(ipi*(q*d)); REexpqd = real(expqd); IMexpqd = imag(expqd);

//	  if (do_verbose==1) {printf("#adding neighbor %i (%6.3f %6.3f %6.3f) of atom %i (%6.3f %6.3f %6.3f)- it contributes to J(s,s'):\n",l,xyz(1),xyz(2),xyz(3),ll,(*inputpars.jjj[ll]).xyz[1],(*inputpars.jjj[ll]).xyz[2],(*inputpars.jjj[ll]).xyz[3]);
//                              } 
   //2. in order to sum up we must take into account that the magnetic unit cell is
   //   larger than the crystallographic one - the ll-l neighbor interaction contributes
   //   to many different components of Js,ss(q) ... note s,ss runs over all the atoms
   //   in the magnetic supercell
         for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
         s=md.in(i1,j1,k1); 

         //calc ss (check in which crystallographic unit ss of the magnetic cell the neighbour l-ll lies)	 
         i=(int)(i1+d_rint(1)-1); // calculate 
	 j=(int)(j1+d_rint(2)-1);
	 k=(int)(k1+d_rint(3)-1);
	 if (i>=0) ij(1)=integer(1.0*i/ini.mf.na())*ini.mf.na();
	 else      ij(1)=(integer(1.0*(i+1)/ini.mf.na())-1)*ini.mf.na();
	 if (j>=0) ij(2)=integer(1.0*j/ini.mf.nb())*ini.mf.nb();
	 else      ij(2)=(integer(1.0*(j+1)/ini.mf.nb())-1)*ini.mf.nb();
	 if (k>=0) ij(3)=integer(1.0*k/ini.mf.nc())*ini.mf.nc();
	 else      ij(3)=(integer(1.0*(k+1)/ini.mf.nc())-1)*ini.mf.nc();
//	 if (do_verbose==1) {printf("#ijk=%i %i %i  ij()=%6.3f %6.3f %6.3f ",i,j,k,ij(1),ij(2),ij(3));}
         i=i-(int)ij(1)+1;
	 j=j-(int)ij(2)+1;
	 k=k-(int)ij(3)+1;
	 ss=md.in(i,j,k);
//          if (do_verbose==1) {printf("#s=%i %i %i  s'=%i %i %i\n",i,j,k,i1,j1,k1);}
          // sum up 

//         mdl - Changed 110710 - To speed up computation by calculating exp(+2i.Pi.Q.d) real and imag parts separately, 
//                                and put into J.mati(s,ss) directly without using intermediate jsss matrix.
           complex<double> **jsss = J.mati(s,ss).M;

//         ComplexMatrix jsss(1,ini.nofcomponents*md.baseindex_max(i1,j1,k1),1,ini.nofcomponents*md.baseindex_max(i,j,k));
//         jsss=0;
        
	  sl=(*inputpars.jjj[ll]).sublattice[l]; // the whole loop has also to be done 
                                                 // for all the other transitions of sublattice sl

//#ifdef _THREADS_JSSS   ... not needed becaues each thread has different ll !! thus writes to a different 
// region of memory (jsi)
//	  MUTEX_LOCK(&mutex_Jlock_write1); 
//           while(Jlock[si]){pthread_cond_wait(&JQfree, &mutex_Jlock_write);}
//           Jlock[si]=1;
//          MUTEX_UNLOCK(&mutex_Jlock_write1);
// #endif
           
          // therefore calculate offset of the set of transitions
          for(tl=1;tl<=md.noft(i1,j1,k1,ll);++tl){ jsi = ini.nofcomponents*(md.baseindex(i1,j1,k1,ll,tl)-1);
	  for(tll=1;tll<=md.noft(i,j,k,sl);++tll){ jsj = ini.nofcomponents*(md.baseindex(i,j,k,sl,tll)-1);

	     for(m=1;m<=ini.nofcomponents;++m){for(n=1;n<=ini.nofcomponents;++n){ //this should also be ok for nofcomponents > 3 !!! (components 1-3 denote the magnetic moment)
         jjval = (*inputpars.jjj[ll]).jij[l](m,n);  
         jsss[jsi+m][jsj+n] += complex<double>(jjval*REexpqd, jjval*IMexpqd);
//         jsss(jsi+m,jsj+n) += complex<double>(jjval*REexpqd, jjval*IMexpqd);
                                              }                                 } // but orbitons should be treated correctly by extending 3 to n !!

	                                         }} 
//#ifdef _THREADS_JSSS
//         MUTEX_LOCK(&mutex_Jlock_write0);
//          Jlock[si]=0;
//         MUTEX_UNLOCK(&mutex_Jlock_write0);
//          pthread_cond_signal(&JQfree);
//#endif
//       J.mati(s,ss)+=jsss;

          ++nofneighbours; // count neighbours summed up
	 }}}
   }
#ifdef _THREADS_JSSS
   myinput->dimA=nofneighbours;
#if defined  (__linux__) || defined (__APPLE__)
    pthread_exit(NULL);
#else
return true;
#endif
#endif

}
#ifdef _THREADS_JSSS
#undef inputpars
#undef ini
#undef md
#undef J
#undef q
#endif


void sortE(Vector & d,ComplexMatrix & z)
{       int i,j,k;
    double p;
    complex <double> p1;
    // lowest and highest column index of the matrix
    int lo = z.Clo();
    int hi = z.Chi();

    for (i = lo; i <= hi; i++) {
	k = i;
	p = d(i);
	for (j = i+1; j <= hi; j++)
	    if (d(j) < p) {
		k = j;
		p = d(j);
	    }
	if (k != i) {
	    d(k) = d(i);
	    d(i) = p;
	    for (j = lo; j <= hi; j++) {
		p1 = z(j,i);
		z(j,i) = z(j,k);
		z(j,k) = p1;
	    }
	}
    }
}
void sortEc(ComplexVector & d,ComplexMatrix & z)
{       int i,j,k;
    complex <double> p;
    complex <double> p1;
    // lowest and highest column index of the matrix
    int lo = z.Clo();
    int hi = z.Chi();

    for (i = lo; i <= hi; i++) {
	k = i;
	p = d(i);
	for (j = i+1; j <= hi; j++)
	    if (real(d(j)) < real(p)) {  // Sorts by the real part of d
		k = j;
		p = d(j);
	    }
	if (k != i) {
	    d(k) = d(i);
	    d(i) = p;
	    for (j = lo; j <= hi; j++) {
		p1 = z(j,i);
		z(j,i) = z(j,k);
		z(j,k) = p1;
	    }
	}
    }
}

// rotate chi(1..3,1..3) from xyz to uvw coordinates
void rottouvw(ComplexMatrix & chi,inimcdis & ini,Vector & abc,int & counter)
{Vector hkl(1,3),u(1,3),v(1,3),w(1,3),q1(1,3),q2(1,3);q1=0;q2=0;
 static Vector wold(1,3);
 ComplexMatrix M(1,3,1,3);
 hkl(1)=ini.hkls[counter][1];
 hkl(2)=ini.hkls[counter][2];
 hkl(3)=ini.hkls[counter][3];
 hkl2ijk(u,hkl,abc);
 int i=counter-1;
 if(i==0){i++;}
 hkl(1)=ini.hkls[i][1];
 hkl(2)=ini.hkls[i][2];
 hkl(3)=ini.hkls[i][3];
 hkl2ijk(q1,hkl,abc);
 while(fabs(fabs(q1*q2)-Norm(q1)*Norm(q2))<SMALL_XPROD_FOR_PARALLEL_VECTORS&&i<ini.nofhkls){
 hkl(1)=ini.hkls[i+1][1];
 hkl(2)=ini.hkls[i+1][2];
 hkl(3)=ini.hkls[i+1][3];
 hkl2ijk(q2,hkl,abc);
     i++;}
//printf("A:q1*q2=%g q1=(%g %g %g) q2=(%g %g %g) i=%i counter=%i\n",q1*q2,q1(1),q1(2),q1(3),q2(1),q2(2),q2(3),i,counter);
 while(fabs(fabs(q1*q2)-Norm(q1)*Norm(q2))<SMALL_XPROD_FOR_PARALLEL_VECTORS&&i>1){i--;
 hkl(1)=ini.hkls[i][1];
 hkl(2)=ini.hkls[i][2];
 hkl(3)=ini.hkls[i][3];
 hkl2ijk(q1,hkl,abc);
     }
//printf("B:q1*q2=%g q1=(%g %g %g) q2=(%g %g %g) i=%i counter=%i\n",q1*q2,q1(1),q1(2),q1(3),q2(1),q2(2),q2(3),i,counter);
 xproduct(w,q1,q2);
 if(Norm(w)<SMALL_XPROD_FOR_PARALLEL_VECTORS){fprintf(stderr,"Error mcdisp: for option outS=3,4 more than 1 linear independent hkl set has to be given in order to determine scattering plane\n");exit(EXIT_FAILURE);}
 xproduct(v,w,u);
 // normalize
 u=u/Norm(u);
 v=v/Norm(v);
 w=w/Norm(w);
 if(Norm(w-wold)>SMALL_NORM){printf("#for this and the following q vectors the vector w (perp to scattering plane) is w=(wx,wy,wz)=(%8.4f %8.4f %8.4f) with x||(b x (a x b)),y||b,z||(a x b)\n",w(1),w(2),w(3));}
 wold=w;
 // now u v w is determined in terms of xyz coordinates
 ComplexMatrix rot(1,3,1,3);
 for(int i=1;i<=3;++i){rot(i,1)=u(i);rot(i,2)=v(i);rot(i,3)=w(i);}
 M=rot.Transpose()*chi(1,3,1,3)*rot;
 for(int i=1;i<=3;++i)for(int j=1;j<=3;++j){chi(i,j)=M(i,j);}
 
}

// *******************************************************************************************
// procedure to calculate the dispersion
void dispcalc(inimcdis & ini,par & inputpars,int calc_rixs,int do_phonon, int do_gobeyond,
              int do_Erefine,int do_jqfile,int do_createtrs,int do_readtrs, int do_verbose,int do_ignore_non_hermitian_matrix_error,
              int maxlevels,double minE,double maxE,double ninit,double pinit,double epsilon, const char * filemode)
{ int i,j,k,l,ll,s,ss,i1,i2,j1,j2,k1,k2,l1,l2,t1,t2,b,bb,m,tn;
 std::clock_t lastcputime = std::clock();
  FILE * fin;
  FILE * fout;
  FILE * foutqom=NULL;
  FILE * foutqei=NULL;
  FILE * foutqep=NULL;
  FILE * foutqee=NULL;
  FILE * foutqem=NULL;
  FILE * foutqes=NULL;
  FILE * foutqel=NULL;
  FILE * foutqsd=NULL;
  FILE * foutqod=NULL;
  FILE * fout1=NULL;
  FILE * foutds=NULL;
  FILE * foutdstot=NULL;
  FILE * foutds1=NULL;
  FILE * jqfile=NULL;
  float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
  double E;
  char filename[MAXNOFCHARINLINE];
  double sta=0,sta_int=0,sta_without_antipeaks=0,sta_int_without_antipeaks=0;
  double sta_without_weights=0,sta_int_without_weights=0,sta_without_antipeaks_weights=0,sta_int_without_antipeaks_weights=0;
  double jqsta=-1.0e10;  double jqsta_int=0;double jqsta_int_scaled=0; double jqsta_scaled=-1.0e10;double scalefactor=1;
  double jq0=0,jqmax=-1e10,hmax,kmax,lmax;
  Vector hkl(1,3),q(1,3),qold(1,3),qijk(1,3);                 
  Vector mf(1,ini.nofcomponents);
  int jmin;
  IntVector noftransitions(1,inputpars.cs.nofatoms); // vector to remember how many transitions are on each atom
  //int offset[inputpars.cs.nofatoms+1]; // vector to remember where higher  transitions are stored
                                    // (as "separate ions on the same unit cell position")
  mf=0;
   int sort=0;int maxiter=1000000;
  time_t curtime;
  struct tm *loctime;
  float d;
  double gamman;
  Vector gamma(1,ini.nofcomponents);
  complex<double> imaginary(0,1);
  ComplexVector u1(1,ini.nofcomponents);
  // transition matrix Mij
  ComplexMatrix Mijkl(1,ini.nofcomponents,1,ini.nofcomponents);
  // transformation matrix Uij
  ComplexMatrix Uijkl(1,ini.nofcomponents,1,ini.nofcomponents);

  ComplexVector chargedensity_coeff1(1,CHARGEDENS_EV_DIM);
  ComplexVector spindensity_coeff1(1,3*SPINDENS_EV_DIM);
  ComplexVector orbmomdensity_coeff1(1,3*ORBMOMDENS_EV_DIM);
  ComplexVector magmom_coeff1(1,MAGMOM_EV_DIM);
  ComplexVector spin_coeff1(1,SPIN_EV_DIM);
  ComplexVector orbmom_coeff1(1,ORBMOM_EV_DIM);
  ComplexVector phonon_coeff1(1,PHONON_EV_DIM);

  //calculate single ion properties of every atom in magnetic unit cell
  int nofEstps=0;if(do_Erefine)nofEstps=(int)((ini.emax-ini.emin)/(fabs(epsilon)/2)+1);
  mdcf md(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),inputpars.cs.nofatoms,ini.nofcomponents,nofEstps,do_Erefine);

  
 if (do_readtrs==0)
 {
 // ********************************************** write mcdisp.trs *******************************************************
 snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.trs",ini.prefix);printf("# saving  %s\n",filename);
  fout = fopen_errchk (filename,"w");
   trs_header_out(fout,pinit,ninit,maxE,ini.T,ini.Hext,'I');
  for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  for(l=1;l<=inputpars.cs.nofatoms;++l){
   if(do_verbose==1)fprintf(stdout,"trying du1calc for ion %i in crystallographic unit cell %i %i %i:\n",l,i,j,k);
    for(ll=1;ll<=ini.nofcomponents;++ll)
     {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} //mf ... mean field vector of atom s in first 
                                                            //crystallographic unit of magnetic unit cell
      md.est_ini(i,j,k,l,(*inputpars.jjj[l]).eigenstates(mf,ini.Hext,ini.T)); 
      (*inputpars.jjj[l]).transitionnumber=0;
      (*inputpars.jjj[l]).maxE=maxE;(*inputpars.jjj[l]).pinit=pinit;(*inputpars.jjj[l]).ninit=ninit;
     noftransitions(l)=0;int noft=0; 
     if(trs_write_next_line(fout,(*inputpars.jjj[l]),noft,i,j,k,l,noftransitions(l),ini.T,mf,ini.Hext,
                    md.est(i,j,k,l),d,minE,maxE,'I',q))
       {fprintf(stderr,"ERROR mcdisp.par: no transition found within energy in range [minE,maxE]=[%g,%g] found\n"
                        " (within first crystallographic unit of magnetic unit cell)\n"
                        " please increase energy range in option -maxE and -minE\n",minE,maxE);
                        exit(EXIT_FAILURE);}
    
     jmin=(*inputpars.jjj[l]).transitionnumber;    // store number of first valid transition
     
// now do  other transitions of the same ion:
   int idummy=0;  
   while(noftransitions(l)<maxlevels&&
         !trs_write_next_line(fout,(*inputpars.jjj[l]),idummy,i,j,k,l,
                              noftransitions(l),ini.T,mf,ini.Hext,md.est(i,j,k,l),d,minE,maxE,'I',q));
 
   (*inputpars.jjj[l]).transitionnumber=jmin; // put back transition number for 1st transition
  }}}}
   fclose(fout);
 // **********************************************end write mcdisp.trs *******************************************************
 } // do_readtrs==0

  if (do_createtrs==1){fprintf(stdout,"single ion transition file ./results/mcdisp.trs created - please comment transitions which should not enter the calculation and restart with option -t\n");exit(0);}
  snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.trs",ini.prefix);
  printf("\n#reading %s\n\n",filename);
// read transitions to be considered from file
 for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
    if((fin = fopen(filename,"rb"))==NULL){snprintf(filename,MAXNOFCHARINLINE,"./results/mcdisp.trs");
                  fin = fopen_errchk(filename,"rb");printf("\n#... not possible, therefore reading %s\n\n",filename);}
noftransitions=0;
 int nparread=0;double Tr,Har,Hbr,Hcr;
 char instr[MAXNOFCHARINLINE];
 while(fgets(instr,MAXNOFCHARINLINE,fin)!=NULL&&nparread<6)
 {nparread+=1-extract(instr,"ninit",ninit);
  nparread+=1-extract(instr,"pinit",pinit);
  nparread+=1-extract(instr,"maxE",maxE);
  nparread+=1-extract(instr,"T",Tr);
  nparread+=1-extract(instr,"Ha",Har);
  nparread+=1-extract(instr,"Hb",Hbr);
  nparread+=1-extract(instr,"Hc",Hcr);
 }

 if (Tr!=ini.T||Har!=ini.Hext(1)||Hbr!=ini.Hext(2)||Hcr!=ini.Hext(3)||nparread!=7){fprintf(stderr,"ERROR: reading mcdisp.trs one of the parameters not set or not in line with mcdisp.mf mcdisp.par: ninit pinit maxE T Ha Hb Hc ! \n");exit(EXIT_FAILURE);}
  while (feof(fin)==0)
  {if ((i1=inputline(fin,nn))>=5)
   {if(i==(int)nn[1]&&j==(int)nn[2]&&k==(int)nn[3])
    {l=(int)nn[4];++noftransitions(l);}
  }} 
     fclose(fin);
  //just read dimensions of matrices in md
    int mqdim=3;if(calc_rixs)mqdim=9;
    md.set_noftransitions(i,j,k,noftransitions,mqdim);
  // for later use:
    for(l=1;l<=inputpars.cs.nofatoms;++l){ //save eigenstates of ions (if possible)
      for(ll=1;ll<=ini.nofcomponents;++ll)
       {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} //mf ... mean field vector of atom s in first 
                                                              //crystallographic unit of magnetic unit cell
       if(do_readtrs!=0)md.est_ini(i,j,k,l,(*inputpars.jjj[l]).eigenstates(mf,ini.Hext,ini.T)); // initialize ests if not already done above
       (*inputpars.jjj[l]).ninit=ninit; // set the constants with values read from mcdisp.trs for each ion so calls to du1calc
       (*inputpars.jjj[l]).pinit=pinit; // dm1calc etc have the same transition number scheme
       (*inputpars.jjj[l]).maxE=maxE;

       }

    md.U(i,j,k)=0; // initialize transformation matrix U
    if(do_Erefine==1)md.M(i,j,k)=0; // initialize matrix M
    md.sqrt_gamma(i,j,k)=0; // and sqrt(gamma^s) matrix sqrt_gamma
 }}}
 
// determine the dimension of the dynamical matrix Ass' s,s'=1....dimA
int dimA=0;
for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 dimA+=md.baseindex_max(i1,j1,k1);
 }}}

// matrix E^s_alpha' used to store the coefficients for extending the eigenvector (see manual)
ComplexMatrix Echargedensity(1,dimA,1,CHARGEDENS_EV_DIM);Echargedensity=0;
ComplexMatrix Espindensity(1,dimA,1,3*SPINDENS_EV_DIM);Espindensity=0;
ComplexMatrix Eorbmomdensity(1,dimA,1,3*ORBMOMDENS_EV_DIM);Eorbmomdensity=0;
ComplexMatrix Ephonon(1,dimA,1,PHONON_EV_DIM);Ephonon=0;
ComplexMatrix Emagmom(1,dimA,1,MAGMOM_EV_DIM);Emagmom=0;
ComplexMatrix Espin(1,dimA,1,SPIN_EV_DIM);Espin=0;
ComplexMatrix Eorbmom(1,dimA,1,ORBMOM_EV_DIM);Eorbmom=0;

  snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.trs",ini.prefix);  
 for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  for(l=1;l<=inputpars.cs.nofatoms;++l){
  if((fin = fopen(filename,"rb"))==NULL){snprintf(filename,MAXNOFCHARINLINE,"./results/mcdisp.trs");
                  fin = fopen_errchk(filename,"rb");}
  jmin=0;
  while (feof(fin)==0)
  {if ((i1=inputline(fin,nn))>=5)
   {if(i==(int)nn[1]&&j==(int)nn[2]&&k==(int)nn[3]&&l==(int)nn[4])
    {tn=(int)nn[5];++jmin;  
    // calculate delta(single ion excitation energy), 
    // Malphabeta(transition matrix elements)

      // do calculation for atom s=(ijkl)
      for(ll=1;ll<=ini.nofcomponents;++ll)
       {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);    //mf ... mean field vector of atom s
        }
     int  n,nd;
      if (do_verbose==1){fprintf(stdout,"#transition %i of ion %i of cryst. unit cell at pos  %i %i %i in mag unit cell:\n",tn,l,i,j,k);
                         if(nn[6]<SMALL_QUASIELASTIC_ENERGY){fprintf(stdout,"#-");}else{fprintf(stdout,"#+");} }     
        j1=(*inputpars.jjj[l]).transitionnumber; // try calculation for transition  j
        if (do_verbose==1){(*inputpars.jjj[l]).transitionnumber=-tn; // try calculation for transition  tn with printout
                          }else{(*inputpars.jjj[l]).transitionnumber=tn;} // no printout
        (*inputpars.jjj[l]).du1calc(ini.T,mf,ini.Hext,u1,d,n,nd,md.est(i,j,k,l));
        if(do_Erefine)Mijkl = u1^u1;
        gamman=Norm2(u1);if(gamman==0)gamman=SMALL_GAMMA;
        u1/=sqrt(gamman);
       if(fabs((fabs(d)-fabs(nn[6]))/(fabs(nn[6])+1.0))>SMALLEDIF)
        {fprintf(stderr,"ERROR mcdisp: reading mcdisp.trs with transition energy delta %g meV different from internal calculation %g meV\n",nn[6],d);	 
         exit(EXIT_FAILURE);}
         // treat correctly case for neutron energy loss
	 if (nn[6]<0) // if transition energy is less than zero do a conjugation of the matrix U
	 {  for(int ii=u1.Lo();ii<=u1.Hi();++ii)u1(ii)=conj(u1(ii));
	 }
       j1=md.baseindex(i,j,k,l,jmin); 
       md.delta(i,j,k)(j1)=nn[6]; // set delta


//----------------------------------OBSERVABLES -------------------------------------------------
if (do_verbose==1){ fprintf(stdout,"# ... recalculate now M(s=%i %i %i %i) with eigenvector dimension for observable\n",i,j,k,l);}
//-----------------------------------------------------------------------------------
if(ini.calculate_chargedensity_oscillation){
   if((*inputpars.jjj[l]).dchargedensity_coeff1(ini.T,mf,ini.Hext,chargedensity_coeff1,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,CHARGEDENS_EV_DIM,chargedensity_coeff1,md,nn,ini,Echargedensity);}
if(ini.calculate_spindensity_oscillation){
   if((*inputpars.jjj[l]).dspindensity_coeff1(ini.T,mf,ini.Hext,spindensity_coeff1,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,3*SPINDENS_EV_DIM,spindensity_coeff1,md,nn,ini,Espindensity);}
if(ini.calculate_orbmomdensity_oscillation){
   if((*inputpars.jjj[l]).dorbmomdensity_coeff1(ini.T,mf,ini.Hext,orbmomdensity_coeff1,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,3*ORBMOMDENS_EV_DIM,orbmomdensity_coeff1,md,nn,ini,Eorbmomdensity);}
if(ini.calculate_phonon_oscillation){
   if((*inputpars.jjj[l]).dP1calc(ini.T,mf,ini.Hext,phonon_coeff1,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,PHONON_EV_DIM,phonon_coeff1,md,nn,ini,Ephonon);}
if(ini.calculate_magmoment_oscillation){
   if((*inputpars.jjj[l]).dm1calc(ini.T,mf,ini.Hext,magmom_coeff1,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,MAGMOM_EV_DIM,magmom_coeff1,md,nn,ini,Emagmom);}
if(ini.calculate_spinmoment_oscillation){
   if((*inputpars.jjj[l]).dS1calc(ini.T,mf,ini.Hext,spin_coeff1,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,SPIN_EV_DIM,spin_coeff1,md,nn,ini,Espin);}
if(ini.calculate_orbmoment_oscillation){
   if((*inputpars.jjj[l]).dL1calc(ini.T,mf,ini.Hext,orbmom_coeff1,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,ORBMOM_EV_DIM,orbmom_coeff1,md,nn,ini,Eorbmom);}
//----------------------------------------------------------------------------------------------

                           // mind in manual the 1st dimension alpha=1 corresponds
			   // to the nth dimension here, because myEigensystmHermitean
			   // sorts the eigenvalues according to ascending order !!!
          if (nn[6]>SMALL_QUASIELASTIC_ENERGY)
			    {md.sqrt_gamma(i,j,k)(j1)=sqrt(gamman);// gamma(ini.nofcomponents)=sqr(gamma^s)
                            }
			    else if (nn[6]<-SMALL_QUASIELASTIC_ENERGY)
                            {md.sqrt_gamma(i,j,k)(j1)=imaginary*sqrt(gamman);// gamma(ini.nofcomponents)=sqr(gamma^s)
                            }
 			    else
			    { //quasielastic line needs gamma=SMALL_QUASIELASTIC_ENERGY .... because Mijkl and therefore gamma have been set to 
			      // wn/kT instead of wn-wn'=SMALL_QUASIELASTIC_ENERGY*wn/kT (in jjjpar.cpp -mdcalc routines)
			      //set fix delta but keep sign
			          if (nn[6]>0){md.delta(i,j,k)(j1)=SMALL_QUASIELASTIC_ENERGY;
  			     md.sqrt_gamma(i,j,k)(j1)=sqrt(SMALL_QUASIELASTIC_ENERGY*gamman);
                                              }
				  else        {md.delta(i,j,k)(j1)=-SMALL_QUASIELASTIC_ENERGY;
                             md.sqrt_gamma(i,j,k)(j1)=imaginary*sqrt(SMALL_QUASIELASTIC_ENERGY*gamman);
			                      }
			    }
        (* inputpars.jjj[l]).transitionnumber=j1; // put back transition number for 1st transition
        for(m=1;m<=ini.nofcomponents;++m){
        md.U(i,j,k)(ini.nofcomponents*(j1-1)+m,j1)=u1(m);
        if(do_Erefine==1)for(n=1;n<=ini.nofcomponents;++n)md.M(i,j,k)(ini.nofcomponents*(j1-1)+m,ini.nofcomponents*(j1-1)+n)=Mijkl(m,n);
        }
if (do_verbose==1){
                  fprintf(stdout,"#Ualpha1 (s=%i %i %i):\n",i,j,k);
                  myPrintComplexVector(stdout,u1); 
                 }

    }}}
    fclose(fin);

  }}}}

//**************************************************************************************
//initialize output files
//************************************************************************************* 
// initialize file with jq matrix
if (do_jqfile)
{  snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.jq",ini.prefix);printf("#saving %s\n",filename);
 if(strcmp(filemode,"A")==0)filemode="a";
 jqfile = fopen_errchk (filename,filemode);
 writeheader(inputpars,jqfile); printf("#saving mcdisp.jq\n");
   fprintf(jqfile,"#!<--mcphas.mcdisp.dsigma.jq-->\n");
   fprintf (jqfile, "#Fourier Transform of 2 Ion Interaction - sta is calculated by comparing the larges eigenvalue\n# to that of the first q vector of the calculation");
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),jqfile);
  if (do_verbose==1){   fprintf (jqfile, "#q=(hkl)\n #spin s() - spin s'()\n #3x3 matrix jss'(q) real im .... [meV]\n");}
  else {if(do_jqfile==1)fprintf(jqfile,"#h  vs  k  vs  l  vs Qincr[1/A] vs largest eigenvalue of J(hkl) matrix (meV) vs components of corresponding eigenvector re im re im re im re im\n");
        if(do_jqfile==2)fprintf(jqfile,"#h  vs  k  vs  l  vs Qincr[1/A] vs eigenvalues of J(hkl) matrix (meV) \n");       
       }
}

// ************************************************************************************************
//MAIN LOOP - do calculation of excitation energy for every Q vector     
// ************************************************************************************************
#ifdef _THREADS
   // Initialises mutual exclusions and threads
   MUTEX_INIT(mutex_loop);
   MUTEX_INIT(mutex_index);
//   MUTEX_INIT(mutex_Jlock_write0);
//   MUTEX_INIT(mutex_Jlock_write1);
   EVENT_INIT(checkfinish);
   #if defined  (__linux__) || defined (__APPLE__)
   pthread_t threads[NUM_THREADS]; int rc; void *status;
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   #else
   HANDLE threads[NUM_THREADS];
   DWORD tid[NUM_THREADS], dwError;
   long unsigned int retval;
   #endif
   int ithread;
   intcalcapr_input *tin[NUM_THREADS];     
   thrdat.ini = new inimcdis*[NUM_THREADS]; thrdat.inputpars = new par*[NUM_THREADS];
   thrdat.md  = new mdcf*[NUM_THREADS];     thrdat.J         = new jq*[NUM_THREADS];
   for (ithread=0; ithread<NUM_THREADS; ithread++) 
   { thrdat.ini[ithread] = new inimcdis(ini); 
     thrdat.inputpars[ithread] = new par(inputpars);
   } 
#endif
int counter,firstcounter=1;qijk=0;double qincr=-1;
if(strcmp(filemode,"A")==0){// check if some q values have already been calculated in a previous run !
                  // ... if hkl values match do not recalculate ... set firstcounter accordingly
                  filemode="a";float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE; 
                 snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qom",ini.prefix);foutqom = fopen(filename,"r");
                 if(foutqom!=NULL){// read file and compare hkls
                                   while (feof(foutqom)==0&&firstcounter<=ini.nofhkls)
                                   {if (inputline(foutqom,nn)>=7)
                                      {if(fabs(ini.hkls[firstcounter][1]-nn[5])<SMALL_HKL_DIFF&&
                                          fabs(ini.hkls[firstcounter][2]-nn[6])<SMALL_HKL_DIFF&&
                                          fabs(ini.hkls[firstcounter][3]-nn[7])<SMALL_HKL_DIFF) 
                                          {++firstcounter;}
                                         else{break;}
                                      }
                                   }
                                   fclose(foutqom);
                                  }               
                 }
if(firstcounter>ini.nofhkls){printf("# mcdisp: all hkl already calculated in previous run - nothing to do - exiting\n");exit(0);}
for(counter=firstcounter;counter<=ini.nofhkls;++counter){
		     hkl(1)=ini.hkls[counter][1];
		     hkl(2)=ini.hkls[counter][2];
		     hkl(3)=ini.hkls[counter][3];

 // transform hkl to primitive lattice
 q=inputpars.cs.r.Transpose()*hkl;

fprintf(stdout,"#q=(%g,%g,%g)",hkl(1),hkl(2),hkl(3));print_time_estimate_until_end((double)(ini.nofhkls-counter)/(counter-firstcounter+1));
fprintf(stdout,"\n");
 if(do_verbose==1){fprintf(stdout,"#Setting up J(q) matrix .... \n");}
 // calculate J(q)
 jq J(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),md.nofcomponents,md);
 jq Jl(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),1,md);

// int signa,signb,signc,sa,sb,sc,
 long int nofneighbours=0;
 // initialize Js,ss(Q)=0 (see manual for description of this matrix)
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
   s=J.in(i1,j1,k1);
   for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
   ss=J.in(i2,j2,k2); 
           J.mati(s,ss)= 0;// set (ini.nofcomponents*nofatoms) x (ini.nofcomponents*nofatoms) matrix Js,ss(q)=0
           Jl.mati(s,ss)= 0;// set Js,ss(q)=0 
   }}}
 }}}

#ifdef _THREADS_JSSS
       thrdat.q = q; thrdat.thread_id = -1;       
   for (ithread=0; ithread<NUM_THREADS; ithread++) 
   {
      tin[ithread] = new intcalcapr_input(dimA,ithread,1,do_verbose,calc_rixs,do_phonon,0.); 
      thrdat.J[ithread] = &J;
     // thrdat.J[ithread] = new jq(J);
      tin[ithread]->dimA=0; 
      thrdat.md[ithread] = new mdcf(md,0);      
   } 
   int thrcount=0, ithread=0;
#endif
#ifdef _THREADS
int num_threads_started=-1;
#endif
 // calculate Js,ss(Q) summing up contributions from the l=1-paranz parameters
   
   for(ll=1;ll<=inputpars.cs.nofatoms;++ll)
   { //sum up l.th neighbour interaction of crystallographic atom ll
     // 1. transform dn(l) to primitive lattice and round it to integer value
  #ifndef _THREADS_JSSS
      jsss_mult(ll,nofneighbours,q,inputpars,ini,J,md);
  #else
      thrcount++;
      tin[ithread]->level=ll;
      #if defined  (__linux__) || defined (__APPLE__)
      rc = pthread_create(&threads[ithread], &attr, jsss_mult, (void *) tin[ithread]);
      if(rc) { printf("Error return code %i from jsss thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
      #else
      threads[ithread] = CreateThread(NULL, 0, jsss_mult, (void *) tin[ithread], 0, &tid[ithread]);
      if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %lu from jsss thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
      #endif
      num_threads_started = ithread+1;
      if(thrcount%NUM_THREADS==0)
      {
         #if defined  (__linux__) || defined (__APPLE__)
         for(int th=0; th<NUM_THREADS; th++)
            rc = pthread_join(threads[th], &status);
         #else
         retval=WaitForMultipleObjects(NUM_THREADS,threads,TRUE,INFINITE);
         if(retval<WAIT_OBJECT_0||retval>WAIT_OBJECT_0+NUM_THREADS-1){printf("Error waitformultipleobjects jsss\n"); exit(EXIT_FAILURE); }
         for(int th=0; th<NUM_THREADS; th++)CloseHandle(threads[th]);
         #endif
         ithread=0;
      }
      else ithread++;
#endif
   }
#ifdef _THREADS_JSSS
    #if defined  (__linux__) || defined (__APPLE__)
    for(int th=0; th<ithread; th++)
       rc = pthread_join(threads[th], &status);
    #else
    if(ithread>0){retval=WaitForMultipleObjects(ithread,threads,TRUE,INFINITE);
    if(retval<WAIT_OBJECT_0||retval>WAIT_OBJECT_0+ithread-1){printf("Error waitformultipleobjects jsssend\n"); exit(EXIT_FAILURE); }
    for(int th=0; th<ithread; th++)CloseHandle(threads[th]);}
    #endif

    //for(int th=0; th<NUM_THREADS; th++) 
   // {
   //    nofneighbours += tin[th]->dimA;
   //    for(int i1=1;i1<=ini.mf.na();++i1) for(int j1=1;j1<=ini.mf.nb();++j1) for(int k1=1;k1<=ini.mf.nc();++k1)
   //       for(int i2=1;i2<=ini.mf.na();++i2) for(int j2=1;j2<=ini.mf.nb();++j2) for(int k2=1;k2<=ini.mf.nc();++k2)
   //          J.mat(i1,j1,k1,i2,j2,k2)+=(*thrdat.J[th]).mat(i1,j1,k1,i2,j2,k2); 
   // }

    for (ithread=0; ithread<NUM_THREADS; ithread++) {
       //delete thrdat.J[ithread];
       delete thrdat.md[ithread];delete tin[ithread]; }
    
#endif


if (do_jqfile){qold=qijk;hkl2ijk(qijk,hkl, inputpars.cs.abc);  if(qincr==-1){qincr=0;qold=qijk;}qincr+=Norm(qijk-qold);
              writehklblocknumber(jqfile,ini,counter);
         
                  if (do_verbose==1){fprintf (jqfile, "#q=(%g, %g, %g) ",hkl(1),hkl(2),hkl(3));
                                     fprintf(jqfile,"nofneighbours= %li\n",nofneighbours);
                                    }
                  else
                  {fprintf (jqfile, "%g  %g  %g %g  ",hkl(1),hkl(2),hkl(3),qincr);}
                  }


if(do_verbose==1){fprintf(stdout,"#Transform J(q) matrix  with U...\n");}

// transform J(s,ss) (with md.U) and multiply 
// (with eigenvalues sqrt(gamma) [here md.sqrt_gamma]
// of matrix Malphabeta ) matrix J to L [here Jl.mati]  ... compare manual
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 s=ini.mf.in(i1,j1,k1);
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  ss=ini.mf.in(i2,j2,k2);
  Jl.mati(s,ss)=md.U(i1,j1,k1).Hermitean()*J.mati(s,ss)*md.U(i2,j2,k2);
   for(int kk=1;kk<=Jl.mati(s,ss).Rhi();++kk)for(int jj=1;jj<=Jl.mati(s,ss).Chi();++jj)
   Jl.mati(s,ss)(kk,jj)*=md.sqrt_gamma(i1,j1,k1)(kk)*conj(md.sqrt_gamma(i2,j2,k2)(jj));

//if (do_verbose==1){
//                  fprintf(stdout,"#J(s=%i%i%i,s''=%i%i%i)=\n",i1,j1,k1,i2,j2,k2);
//                  myPrintComplexMatrix(stdout,J.mati(s,ss)); 
//                  fprintf(stdout,"#sqr(gamma_s=%i%i%i)=\n",i1,j1,k1);
//                  myPrintComplexVector(stdout,md.sqrt_gamma(i1,j1,k1));
//                  fprintf(stdout,"#U(s=%i%i%i)=\n",i1,j1,k1);
//                  myPrintComplexMatrix(stdout,md.U(i1,j1,k1));
//                  fprintf(stdout,"#sqr(gamma_s=%i%i%i)=\n",i2,j2,k2);
//                  myPrintComplexVector(stdout,md.sqrt_gamma(i2,j2,k2));
//                  fprintf(stdout,"#U(s=%i%i%i)=\n",i2,j2,k2);
//                  myPrintComplexMatrix(stdout,md.U(i2,j2,k2));
//                  fprintf(stdout,"#sqr(gamma_s) U(s)T* J(s=%i%i%i,s''=%i%i%i) U(s'') sqr(gamma_s'')*=\n",i1,j1,k1,i2,j2,k2);
//                  myPrintComplexMatrix(stdout,Jl.mati(s,ss)); 
//                 }
  }}}
 }}}


// calculate Ac
if(do_verbose==1){fprintf(stdout,"#calculating matrix A\n");}
// Ac  is the matrix A which is given in manual chapter 9.2.1 (eq (30) ff) 
// -- diagonalization gives omega_r and Tau
   ComplexMatrix Ac(1,dimA,1,dimA);
   ComplexMatrix Lambda(1,dimA,1,dimA);
   ComplexMatrix Tau(1,dimA,1,dimA);
   ComplexMatrix J_Q(1,ini.nofcomponents*dimA,1,ini.nofcomponents*dimA);
   Ac=0;J_Q=0;Lambda=0;
   for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){

    for(l1=1;l1<=inputpars.cs.nofatoms;++l1){
     for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
      s=index_s(i1,j1,k1,l1,t1,md,ini);
      b=md.baseindex(i1,j1,k1,l1,t1);

                       // this is standard DMD as described in the review rotter et al JPcondMat 2012
                       if(md.delta(i1,j1,k1)(b)<0){Lambda(s,s)=-1;}else{Lambda(s,s)=+1;}
                       Ac(s,s)=md.delta(i1,j1,k1)(b)*Lambda(s,s);
                     
      if(do_verbose==1){fprintf(stdout,"#i=%i j=%i k=%i atomnr=%i trans=%i ... s=%i ",i1,j1,k1,l1,t1,s);
                        fprintf(stdout,"#lambda(%i,%i)xdelta(%i)=%g + i %g\n",s,s,s,real(Ac(s,s)),imag(Ac(s,s)));
                       }
      }} 
   for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
    for(l1=1;l1<=inputpars.cs.nofatoms;++l1){ 
     for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
    for(l2=1;l2<=inputpars.cs.nofatoms;++l2){ 
     for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2){
      s=index_s(i1,j1,k1,l1,t1,md,ini);
      ss=index_s(i2,j2,k2,l2,t2,md,ini);
      b=md.baseindex(i1,j1,k1,l1,t1);
      bb=md.baseindex(i2,j2,k2,l2,t2);
//     Ac(s,ss)-=Jl.mati(Jl.in(i1,j1,k1),Jl.in(i2,j2,k2))(ini.nofcomponents*b,ini.nofcomponents*bb);
       Ac(s,ss)-=Jl.mati(Jl.in(i1,j1,k1),Jl.in(i2,j2,k2))(b,bb);
                                             //nofcomponents^th dimension corresponds to 1st in manual 
					     // and it is only necessary to take into 
					     // acount this dimension!!
       for(i=1;i<=ini.nofcomponents;++i){for(j=1;j<=ini.nofcomponents;++j){
        J_Q(ini.nofcomponents*(s-1)+i,ini.nofcomponents*(ss-1)+j)+=J.mati(J.in(i1,j1,k1),J.in(i2,j2,k2))(ini.nofcomponents*(b-1)+i,ini.nofcomponents*(bb-1)+j);
       }}

     }}}}

  }}}
 }}}

// printout Fouriertransform of matrix jq
if (do_jqfile){
       if (do_verbose==1)
       {//fprintf (jqfile, "#spin (%i*r1 %i*r2 %i*r3) - spin (%i*r1 %i*r2 %i*r3)\n",i1,j1,k1,i2,j2,k2);
         myPrintComplexMatrix(jqfile,J_Q); 
       }

	// diagonalize JQ to get eigenvalues (biggest corresponds to Tn) !!!
         Vector Tn(1,ini.nofcomponents*ini.mf.n()*inputpars.cs.nofatoms);
         ComplexMatrix eigenvectors(1,ini.nofcomponents*ini.mf.n()*inputpars.cs.nofatoms,1,ini.nofcomponents*ini.mf.n()*inputpars.cs.nofatoms);
         myEigenSystemHermitean (J_Q,Tn,eigenvectors,sort=1,maxiter);
         i2=ini.nofcomponents*ini.mf.n()*inputpars.cs.nofatoms;
       if(do_verbose==1)
       {fprintf(jqfile,"#eigenvalues(highest corresponds to Tn, predicted magstructure)\n");
         myPrintVector(jqfile,Tn); 
        fprintf(jqfile,"#eigenvectors(moment direction):\n");
         myPrintComplexMatrix(jqfile,eigenvectors); 
       }
       else
       {if(do_jqfile==1){// print largest eigenvalue and eigenvector
        fprintf(jqfile," %g ",Tn(i2));
        for (i1=1;i1<=i2;++i1)
        {fprintf(jqfile," %6.3g ",real(eigenvectors(i1,i2)));
         fprintf(jqfile," %6.3g ",imag(eigenvectors(i1,i2)));
        }
                         }
        if(do_jqfile==2){// print all eigenvalues
        for (i1=i2;i1>0;--i1)
        {fprintf(jqfile," %g ",Tn(i1));}
                        }
        fprintf(jqfile,"\n");
       }
       // if we are calculating initial q-vektor i.e. jqsta=-1e10 ....
       if (jqsta<-0.9e10){jq0=Tn(i2);jqsta=-0.1e10;jqsta_scaled=jqsta;scalefactor=1;jqmax=jq0;hmax=hkl(1);kmax=hkl(2);lmax=hkl(3);
                          // here the first hkl vectors eigenvalue has been determined
                          // ... look in table, if there is a value, to which is might be scaled
                          if(ini.hkls[counter][0]>3&&jq0!=0.0){scalefactor=ini.hkls[counter][4]/jq0;}
                         }
       // ... on subsequent runs 
       else           {if(Tn(i2)>jqmax){jqmax=Tn(i2);hmax=hkl(1);kmax=hkl(2);lmax=hkl(3);}
                       if(Tn(i2)>jq0)
                         {if(jqsta<0){jqsta=0;jqsta_scaled=0;}
                          jqsta+=(Tn(i2)-jq0)*(Tn(i2)-jq0);jqsta_scaled+=(Tn(i2)-jq0)*(Tn(i2)-jq0)*scalefactor*scalefactor;}
                       else{if((jqsta<0)&(Tn(i2)-jq0>jqsta)){jqsta=Tn(i2)-jq0;jqsta_scaled=(Tn(i2)-jq0)*scalefactor;}}
                       
                      }
                      double test;
                     for(int kk=NOFHKLCOLUMNS;kk<=ini.hkls[counter][0];kk+=NOFHKLCOLUMNS-3){
                      for (j1=1;j1<=ini.hkls[counter][kk-NOFHKLCOLUMNS+7];++j1)
	              {test=fabs(Tn(i2+1-j1)-ini.hkls[counter][kk-NOFHKLCOLUMNS+j1+3]);
                       jqsta_int+=test*test;
                       test=fabs(Tn(i2+1-j1)*scalefactor-ini.hkls[counter][kk-NOFHKLCOLUMNS+j1+3]);
                       jqsta_int_scaled+=test*test;
                       //fprintf(stdout,"%i %g test=%g\n",counter,ini.hkls[counter][0],test);                       
                      }
    	                                                                  }
 }
 else
 {// no jqfile but excitations to be calculated
 if(do_verbose==1){fprintf(stdout,"#diagonalizing %ix%i matrix A, A=\n",dimA,dimA);
                           myPrintComplexMatrix(stdout,Ac); 
                           myPrintComplexMatrix(stdout,Lambda); 
                   }

   // diagonalize Ac to get energies  and eigenvectors !!!
   Vector En(1,dimA);
   ComplexVector Enc(1,dimA);  // For complex eigenvalues in case the non-symmetric eigensolver is used.
   Vector ints(1,dimA);
   Vector intsbey(1,dimA);
   Vector intsP(1,dimA);
//   myEigenValuesHermitean (Ac,En,sort,maxiter);
  
//   myEigenSystemHermitean (Ac,En,Tau,sort,maxiter);
//    myPrintVector(stdout,En);
   int eigrval = myEigenSystemHermiteanGeneral (Lambda,Ac,En,Tau,sort=0,maxiter);
   bool notposdef = false;
   // myEigenSystemHermiteanGeneral return values modified (MDL 141019)
   if(eigrval==0) 
   {
     En=1.0/En;
     sortE(En,Tau);
        //  Tau=Tau.Conjugate();
  	// conjugate inserted 31.10.05, because when calculating simple AF - I noticed
	// that the eigensystemhgermitean returns eigenvectors as column vectors, but
	// the components need to be complex conjugated
         // conjugate removed again MR 11.4.2011 because now done correctly in myev.c
   }
   else if(eigrval>0) {  // 0==sucess. +1==Lambda not hermitian. +2==Ac not hermitian, +3==Lambda AND Ac not hermitian
      if(eigrval>1) {
         fprintf(stderr,"# Dynamical Matrix Ac not hermitian. Check exchange parameter file mcphas.j is consistent\n");
         fprintf(stderr,"#   Note that each interaction between pairs of ions must match - e.g. Interaction between\n");
         fprintf(stderr,"#   Atom 1 with Neighbour 2 (which is atom 2) must equal interaction between Atom 2 with Neighbour 1\n");
         fprintf(stderr,"#   Press q to quit, or any other key to ignore this error.\n"); 
         if(do_ignore_non_hermitian_matrix_error==0) {
            if(getchar()=='q') exit(1); 
         } else {
            // Skips q-point, but make sure qincr is correct.
            fprintf(stderr,"# Skipping this q-point.\n");
            if(qincr!=-1) { // In order to keep the q-increments the same - since we're missing a point here.
               qold=qijk; hkl2ijk(qijk,hkl, inputpars.cs.abc); qincr+=Norm(qijk-qold);
               ini.print_usrdefcols(foutqei,qijk,qincr,q);
               fprintf (foutqei, "%4.4g %4.4g %4.4g  %4.4g %4.4g           ",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)), myround(Norm(qijk)),0.);
               fprintf(foutqei, "-1    -1   -1\n");
            }
            continue;
         }
      }
      if(eigrval%2==1) fprintf(stderr,"# Warning: Trace matrix Lambda is not diagonal\n");
      notposdef = true;
                 }
   else if(eigrval<0 && eigrval>=-dimA) { // -dimA<eigrval<0 means algorithm did not converge.
      fprintf(stderr,"# The Hermitian-definite eigensolver failed to converge. Trying again with a non-symmetric eigensolver.\n");
      fprintf(stderr,"#   This can generate complex/imaginary eigenvalues, which will be outputed as comments in mcdisp.qei\n");
      fprintf(stderr,"#   without intensities. Real eigenvalues will be treated normally.\n");
      notposdef = true;
   }
   else if(eigrval<-dimA) {               // eigrval<-dimA means Ac is not positive definite;
      fprintf(stderr,"# At Q=(%g,%g,%g), the dynamical matrix Ac(Q) is not positive definite, so real eigenvalues are not guaranteed.\n",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
      fprintf(stderr,"#   possible reason: magnetic structure in mcdisp.mf is metastable, leading to soft modes at this Q-vector which could\n");
      fprintf(stderr,"#   be a possible ordering wavevector for the actual stable structure. It could also be due to rounding errors.\n");
      fprintf(stderr,"#   McDisp will now use a non-symmetric eigensolver. Real eigenvalues will be treated as normal. Complex/imaginary eigenvalues\n");
      fprintf(stderr,"#   will be outputed to mcdisp.qei in a commented out line, but the corresponding intensities will not be calculated.\n");
      notposdef = true;
   }
   if(notposdef) {
      eigrval = myEigenSystemGeneral(Lambda,Ac,Enc,Tau);
      if(eigrval!=0) { 
         fprintf(stderr,"# The non-symmetric eigensolver failed. This Q point will be skipped.\n");
         if(qincr!=-1) { // In order to keep the q-increments the same - since we're missing a point here.
            qold=qijk; hkl2ijk(qijk,hkl, inputpars.cs.abc); qincr+=Norm(qijk-qold);
            ini.print_usrdefcols(foutqei,qijk,qincr,q);
            fprintf (foutqei, "%4.4g %4.4g %4.4g  %4.4g %4.4g           ",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)), myround(Norm(qijk)),0.);
            fprintf(foutqei, "-1    -1   -1\n");
         }
         continue;
      }
      Enc=1./Enc;
      sortEc(Enc,Tau);    // Sorts by the real part of the eigenvalues
      for(i=1; i<=dimA; i++) En(i) = (imag(Enc(i))==0) ? real(Enc(i)) : -DBL_MAX;  // Sets -DBL_MAX as flag that eigenvalue is complex
   }

 if(do_verbose==1){   ComplexMatrix test(1,dimA,1,dimA);
   // check normalisation of eigenvectors -------------------- only do this in verbose mode MR 5.6.2013
   bool notnorm = false;
   test=Tau.Hermitean()*Ac*Tau;
   if(notposdef) {    // For the non-symmetric solver, only the eigenvectors corresponding to real eigenvalues have been normalised.
     for(i=1; i<=dimA; i++) {
       if(En(i)!=-DBL_MAX) { if(fabs(fabs(real(test(i,i)))-1)>SMALL_QUASIELASTIC_ENERGY) { notnorm = true; break; } }
     }
   } else { 
     ComplexMatrix unit(1,dimA,1,dimA);unit=1;
     if( NormFro(unit-test)>SMALL_QUASIELASTIC_ENERGY) notnorm = true;
   }
   if(notnorm) {
    myPrintComplexMatrix(stdout,test); 
 fprintf(stderr,"Error: eigenvectors t not correctly normalised\n"); 
 fprintf(stderr,"   Press q to quit, or any other key to ignore this error.\n"); if(getchar()=='q') exit(1); }}
//-------------------------------------------------------
 if(do_verbose==1){ fprintf(stdout,"#eigenvectors (matrix Tau):\n");
                    myPrintComplexMatrix(stdout,Tau); 
               fprintf(stdout,"#saving the following eigenvalues (meV) to mcdisp.qom:\n");
   for (i=1;i<=dimA;++i){fprintf(stdout, " %4.4g",En(i));}
                         }        


   // calculate and printout intensities [the energies have already
   // been printed out above, so any refinement of energies during intcalc
   // is not included in the output file]
  double QQ; 
#ifndef _THREADS  
             mfcf qee_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,CHARGEDENS_EV_DIM);
             mfcf qee_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,CHARGEDENS_EV_DIM);
             mfcf qsd_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,3*SPINDENS_EV_DIM);
             mfcf qsd_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,3*SPINDENS_EV_DIM);
             mfcf qod_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,3*ORBMOMDENS_EV_DIM);
             mfcf qod_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,3*ORBMOMDENS_EV_DIM);
             mfcf qep_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,PHONON_EV_DIM);
             mfcf qep_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,PHONON_EV_DIM);
             mfcf qem_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,MAGMOM_EV_DIM);
             mfcf qem_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,MAGMOM_EV_DIM);
             mfcf qes_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPIN_EV_DIM);
             mfcf qes_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPIN_EV_DIM);
             mfcf qel_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOM_EV_DIM);
             mfcf qel_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOM_EV_DIM);
#endif
  double DMDtotint=0,DMDtotintbey=0;
  ComplexMatrix chitot(1,3,1,3),chitotbey(1,3,1,3); chitot=0;chitotbey=0;
  if(do_verbose==1){fprintf(stdout,"\n#calculating  intensities approximately ...\n");}
  intcalc_ini(ini,inputpars,md,do_Erefine,epsilon,do_verbose,do_gobeyond,calc_rixs,do_phonon,hkl,counter);
  qold=qijk;hkl2ijk(qijk,hkl, inputpars.cs.abc);QQ=Norm(qijk);

  if(qincr==-1){qincr=0;qold=qijk;
              // for the first q vector in the loop we have to initialize files ...
              snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qom",ini.prefix);foutqom = fopen_errchk (filename,filemode);printf("#saving %s\n",filename);
              if(calc_rixs){snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qex",ini.prefix);foutqei = fopen_errchk (filename,filemode);printf("#saving %s\n",filename);}
                     else {snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qei",ini.prefix);foutqei = fopen_errchk (filename,filemode);printf("#saving %s\n",filename);
                           snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.dsigma.tot",ini.prefix);foutdstot = fopen_errchk (filename,filemode);printf("#saving %s\n",filename);
                           if(do_Erefine==1){
                           snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.dsigma",ini.prefix);foutds = fopen_errchk (filename,filemode);printf("#saving %s\n",filename);                    
                                            }
                           }
               writeheaders(foutqom,foutqei,foutdstot,foutds,inputpars,ini,calc_rixs,do_Erefine);                  
               //------------observables-----------------------------------
               if(ini.calculate_chargedensity_oscillation){snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qee",ini.prefix);foutqee=evfileinit(filemode,filename,inputpars,"qee",CHARGEDENS_EV_DIM);}
               if(ini.calculate_spindensity_oscillation)  {snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qsd",ini.prefix);foutqsd=evfileinit(filemode,filename,inputpars,"qsd",3*SPINDENS_EV_DIM);}
               if(ini.calculate_orbmomdensity_oscillation){snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qod",ini.prefix);foutqod=evfileinit(filemode,filename,inputpars,"qod",3*ORBMOMDENS_EV_DIM);}
               if(ini.calculate_phonon_oscillation)       {snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qep",ini.prefix);foutqep=evfileinit(filemode,filename,inputpars,"qep",PHONON_EV_DIM);}
               if(ini.calculate_magmoment_oscillation)    {snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qem",ini.prefix);foutqem=evfileinit(filemode,filename,inputpars,"qem",MAGMOM_EV_DIM);}
               if(ini.calculate_spinmoment_oscillation)   {snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qes",ini.prefix);foutqes=evfileinit(filemode,filename,inputpars,"qes",SPIN_EV_DIM);}
               if(ini.calculate_orbmoment_oscillation)    {snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qel",ini.prefix);foutqel=evfileinit(filemode,filename,inputpars,"qel",ORBMOM_EV_DIM);}
               //-----------------------------------------------------------
               lastcputime=std::clock();
              } else
             {// close and reopen files to prevent data loss if process ends 
               if ((std::clock() - lastcputime) / (double)CLOCKS_PER_SEC >60)
              {lastcputime=std::clock();
              fclose(foutqom);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qom",ini.prefix);foutqom = fopen_errchk (filename,"a");
              fclose(foutqei);
              if(calc_rixs){snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qex",ini.prefix);foutqei = fopen_errchk (filename,"a");}
                     else {snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qei",ini.prefix);foutqei = fopen_errchk (filename,"a");
                           fclose(foutdstot);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.dsigma.tot",ini.prefix);foutdstot = fopen_errchk (filename,"a");
                          if(do_Erefine==1){
                           fclose(foutds);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.dsigma",ini.prefix);foutds = fopen_errchk (filename,"a");               
                                            }
                           }
                                
               //------------observables-----------------------------------
               if(ini.calculate_chargedensity_oscillation){fclose(foutqee);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qee",ini.prefix);foutqee=fopen_errchk (filename,"a");}
               if(ini.calculate_spindensity_oscillation)  {fclose(foutqsd);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qsd",ini.prefix);foutqsd=fopen_errchk (filename,"a");}
               if(ini.calculate_orbmomdensity_oscillation){fclose(foutqod);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qod",ini.prefix);foutqod=fopen_errchk (filename,"a");}
               if(ini.calculate_phonon_oscillation)       {fclose(foutqep);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qep",ini.prefix);foutqep=fopen_errchk (filename,"a");}
               if(ini.calculate_magmoment_oscillation)    {fclose(foutqem);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qem",ini.prefix);foutqem=fopen_errchk (filename,"a");}
               if(ini.calculate_spinmoment_oscillation)   {fclose(foutqes);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qes",ini.prefix);foutqes=fopen_errchk (filename,"a");}
               if(ini.calculate_orbmoment_oscillation)    {fclose(foutqel);snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.qel",ini.prefix);foutqel=fopen_errchk (filename,"a");}
               //-----------------------------------------------------------
               }

             }
         qincr+=Norm(qijk-qold); 
         writehklblocknumber(foutqom,foutqei,foutdstot,foutds,foutqee,foutqsd,foutqod,foutqep,foutqem,foutqes,foutqel,
                             ini,calc_rixs,do_Erefine,counter);
                  ini.print_usrdefcols(foutqom,qijk,qincr,q);
                  fprintf (foutqom, "%4.4g %4.4g  %4.4g ",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
                  for (i=1;i<=dimA;++i)fprintf (foutqom, " %4.4g ",myround(En(i)));
                  fprintf (foutqom, " > ");

                  int dim=3;
                  dim=(int)((ini.hkls[counter][0]-3)/4);
      Vector dd(1,dim);dd=0;dd+=100000.0;
      Vector dd_int(1,dim);dd_int=0;  dd_int+=100000.0;
      Vector dd1(1,dim);dd1=0;  dd1+=100000.0;
      Vector dd1_int(1,dim);dd1_int=0;  dd1_int+=100000.0;
      Vector dd_without_antipeaks(1,dim);dd_without_antipeaks=0;dd_without_antipeaks+=100000.0;
      Vector dd_int_without_antipeaks(1,dim);  dd_int_without_antipeaks=0; dd_int_without_antipeaks+=100000.0;
      Vector dd_without_weights(1,dim);dd_without_weights=0;dd_without_weights+=100000.0;
      Vector dd_int_without_weights(1,dim); dd_int_without_weights=0; dd_int_without_weights+=100000.0;
      Vector dd_without_antipeaks_weights(1,dim);dd_without_antipeaks_weights=0;dd_without_antipeaks_weights+=100000.0;
      Vector dd_int_without_antipeaks_weights(1,dim); dd_int_without_antipeaks_weights=0; dd_int_without_antipeaks_weights+=100000.0;
#ifndef _THREADS
                     int dimchi=3,dimchibey=3;if(calc_rixs){dimchi=9;dimchibey=1;}
                     ComplexMatrix chiPhon(1,1,1,1);
                     ComplexMatrix chi(1,dimchi,1,dimchi);
                     ComplexMatrix chibey(1,dimchibey,1,dimchibey);
                     Matrix pol(1,3,1,3);                     
#else               // Populates the thread data structure
                  thrdat.qee_real = new mfcf*[NUM_THREADS];          thrdat.qee_imag = new mfcf*[NUM_THREADS];
                  thrdat.qsd_real = new mfcf*[NUM_THREADS];          thrdat.qsd_imag = new mfcf*[NUM_THREADS];
                  thrdat.qod_real = new mfcf*[NUM_THREADS];          thrdat.qod_imag = new mfcf*[NUM_THREADS];
                  thrdat.qep_real = new mfcf*[NUM_THREADS];          thrdat.qep_imag = new mfcf*[NUM_THREADS];
                  thrdat.qem_real = new mfcf*[NUM_THREADS];          thrdat.qem_imag = new mfcf*[NUM_THREADS];
                  thrdat.qes_real = new mfcf*[NUM_THREADS];          thrdat.qes_imag = new mfcf*[NUM_THREADS];
                  thrdat.qel_real = new mfcf*[NUM_THREADS];          thrdat.qel_imag = new mfcf*[NUM_THREADS];
                  thrdat.chi      = new ComplexMatrix*[NUM_THREADS]; thrdat.chibey   = new ComplexMatrix*[NUM_THREADS];
                  thrdat.chiPhon     = new ComplexMatrix*[NUM_THREADS]; 
                  thrdat.pol      = new Matrix*[NUM_THREADS];        
                  thrdat.Echargedensity       = new ComplexMatrix*[NUM_THREADS]; 
                  thrdat.Espindensity       = new ComplexMatrix*[NUM_THREADS]; 
                  thrdat.Eorbmomdensity       = new ComplexMatrix*[NUM_THREADS]; 
                  thrdat.Ephonon       = new ComplexMatrix*[NUM_THREADS]; 
                  thrdat.Emagmom       = new ComplexMatrix*[NUM_THREADS]; 
                  thrdat.Espin       = new ComplexMatrix*[NUM_THREADS]; 
                  thrdat.Eorbmom       = new ComplexMatrix*[NUM_THREADS]; 
                  thrdat.Tau      = new ComplexMatrix*[NUM_THREADS]; 
                  thrdat.hkl = hkl; thrdat.thread_id = -1;
                  for (ithread=0; ithread<NUM_THREADS; ithread++) 
                  {
                     tin[ithread] = new intcalcapr_input(dimA,ithread,1,do_verbose,calc_rixs,do_phonon,En);
                     int dimchi=3,dimchibey=3;if(calc_rixs){dimchi=9;dimchibey=1;}
                     thrdat.chiPhon[ithread] = new ComplexMatrix(1,1,1,1);
                     thrdat.chi[ithread] = new ComplexMatrix(1,dimchi,1,dimchi);
                     thrdat.chibey[ithread] = new ComplexMatrix(1,dimchibey,1,dimchibey);
                     thrdat.pol[ithread] = new Matrix(1,3,1,3);
                     thrdat.md[ithread] = new mdcf(md,1);

                     if(ini.calculate_chargedensity_oscillation){thrdat.qee_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,CHARGEDENS_EV_DIM);
                     thrdat.qee_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,CHARGEDENS_EV_DIM);
                     thrdat.Echargedensity[ithread] = new ComplexMatrix(1,dimA,1,CHARGEDENS_EV_DIM); *thrdat.Echargedensity[ithread]=Echargedensity;}
                     if(ini.calculate_spindensity_oscillation){thrdat.qsd_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,3*SPINDENS_EV_DIM);
                     thrdat.qsd_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,3*SPINDENS_EV_DIM);
                     thrdat.Espindensity[ithread] = new ComplexMatrix(1,dimA,1,3*SPINDENS_EV_DIM); *thrdat.Espindensity[ithread]=Espindensity;}
                     if(ini.calculate_orbmomdensity_oscillation){thrdat.qod_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,3*ORBMOMDENS_EV_DIM);
                     thrdat.qod_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,3*ORBMOMDENS_EV_DIM);
                     thrdat.Eorbmomdensity[ithread] = new ComplexMatrix(1,dimA,1,3*ORBMOMDENS_EV_DIM); *thrdat.Eorbmomdensity[ithread]=Eorbmomdensity;}
                     if(ini.calculate_phonon_oscillation){thrdat.qep_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,PHONON_EV_DIM);
                     thrdat.qep_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,PHONON_EV_DIM);
                     thrdat.Ephonon[ithread] = new ComplexMatrix(1,dimA,1,PHONON_EV_DIM); *thrdat.Ephonon[ithread]=Ephonon;}
                     if(ini.calculate_magmoment_oscillation){thrdat.qem_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,MAGMOM_EV_DIM);
                     thrdat.qem_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,MAGMOM_EV_DIM);
                     thrdat.Emagmom[ithread] = new ComplexMatrix(1,dimA,1,MAGMOM_EV_DIM); *thrdat.Emagmom[ithread]=Emagmom;}
                     if(ini.calculate_spinmoment_oscillation){thrdat.qes_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPIN_EV_DIM);
                     thrdat.qes_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPIN_EV_DIM);
                     thrdat.Espin[ithread] = new ComplexMatrix(1,dimA,1,SPIN_EV_DIM); *thrdat.Espin[ithread]=Espin;}
                     if(ini.calculate_orbmoment_oscillation){thrdat.qel_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOM_EV_DIM);
                     thrdat.qel_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOM_EV_DIM);
                     thrdat.Eorbmom[ithread] = new ComplexMatrix(1,dimA,1,ORBMOM_EV_DIM); *thrdat.Eorbmom[ithread]=Eorbmom;}

//                     thrdat.Tau[ithread] = new ComplexMatrix(1,dimA,1,dimA); *thrdat.Tau[ithread]=Tau;
                  thrdat.Tau[ithread]=&Tau; // try to save memory by not copying Tau for every thread !
                  }
                  ithread=0; num_threads_started=0; int oldi=-1;// Vector vQQ(1,dimA); removed MR 14.1.2013
#endif
                      

#ifdef _THREADS  
                  for (i=1;i<=dimA;i+=NUM_THREADS)
#else
                  for (i=1;i<=dimA;++i)
#endif
                  {
#ifdef _THREADS  
                     oldi=i;
                     // Runs threads until all are running - but wait until they are completed in order before printing output.
                     //  This is to ensure that the output is exactly the same as for a single thread (otherwise it would be out of order).
                     for(int th=0; th<NUM_THREADS; th++)
                     { 
                        i=oldi+th; if(i>dimA) break;
                        if(do_gobeyond==0) intsbey(i)=-1.1; else intsbey(i)=+1.1;
                      if (En(i)!=-DBL_MAX && // Matrix Ac is not +ve definite, and this eigenvalue is not real - don't do intensity
                          En(i)<=ini.emax&&En(i)>=ini.emin) // only do intensity calculation if within energy range
                      {if(do_verbose)printf("calling thread %i ",num_threads_started);
                       tin[num_threads_started]->En=En(i); tin[num_threads_started]->intensitybey=intsbey(i);tin[num_threads_started]->intensityP=intsP(i); 
                       tin[num_threads_started]->level = i;
                        #if defined  (__linux__) || defined (__APPLE__)
                        rc = pthread_create(&threads[num_threads_started], &attr, intcalc_approx, (void *) tin[num_threads_started]);
                        if(rc) { printf("Error return code %i from intcalc thread %i\n",rc,num_threads_started+1); exit(EXIT_FAILURE); }
                        #else
                        threads[num_threads_started] = CreateThread(NULL, 0, intcalc_approx, (void *) tin[num_threads_started], 0, &tid[num_threads_started]);
                        if(threads[num_threads_started]==NULL) { dwError=GetLastError(); printf("Error code %lu from intcalc thread %i\n",dwError,num_threads_started+1); exit(EXIT_FAILURE); }
                        #endif  
                        ++num_threads_started;
                      }
                      else {ints(i)=-1;intsbey(i)=-1;intsP(i)=-1;}
                      }
                     #if defined  (__linux__) || defined (__APPLE__)
                     for(int th=0; th<num_threads_started; th++)rc = pthread_join(threads[th], &status);
                     #else
                     if(num_threads_started>0){retval=WaitForMultipleObjects(num_threads_started,threads,TRUE,INFINITE);
                     if(retval<WAIT_OBJECT_0||retval>WAIT_OBJECT_0+num_threads_started-1){printf("Error waitformultipleobjects=%li num_threads_started=%i\n",retval,num_threads_started); exit(EXIT_FAILURE); }
                      for(int th=0; th<num_threads_started; th++)CloseHandle(threads[th]);}         
                     #endif
                     num_threads_started=0; 
                     #define chi      (*thrdat.chi[ithread])
                     #define chibey   (*thrdat.chibey[ithread])
                     #define qee_real (*thrdat.qee_real[ithread])
                     #define qee_imag (*thrdat.qee_imag[ithread])
                     #define qsd_real (*thrdat.qsd_real[ithread])
                     #define qsd_imag (*thrdat.qsd_imag[ithread])
                     #define qod_real (*thrdat.qod_real[ithread])
                     #define qod_imag (*thrdat.qod_imag[ithread])
                     #define qep_real (*thrdat.qep_real[ithread])
                     #define qep_imag (*thrdat.qep_imag[ithread])
                     #define qem_real (*thrdat.qem_real[ithread])
                     #define qem_imag (*thrdat.qem_imag[ithread])
                     #define qes_real (*thrdat.qes_real[ithread])
                     #define qes_imag (*thrdat.qes_imag[ithread])
                     #define qel_real (*thrdat.qel_real[ithread])
                     #define qel_imag (*thrdat.qel_imag[ithread])
                     ithread=-1;
                     for(int th=0; th<NUM_THREADS; th++)
                     {
                         i=oldi+th; if(i>dimA) break;
                        if (En(i)!=-DBL_MAX)   // Matrix Ac is not +ve definite, and this eigenvalue is not real - don't do intensity
                        if (En(i)<=ini.emax&&En(i)>=ini.emin) // only do intensity calculation if within energy range
                      {++ithread;
                       ints(tin[ithread]->level) = tin[ithread]->intensity; 
                       intsbey(tin[ithread]->level) = tin[ithread]->intensitybey;
                       intsP(tin[ithread]->level) = tin[ithread]->intensityP;
                       //printf("%i %g",ithread,ints(tin[ithread]->level));
                      }       
#else
                     if(do_gobeyond==0){intsbey(i)=-1.1;}else{intsbey(i)=+1.1;}
                     if (En(i)!=-DBL_MAX && // Matrix Ac is not +ve definite, and this eigenvalue is not real - don't do intensity
                         En(i)<=ini.emax&&En(i)>=ini.emin) // only do intensity calculation if within energy range
                     {
                     ints(i)=intcalc_approx(chi,chibey,chiPhon,pol,intsbey(i),intsP(i),
                                            qee_real,qee_imag,Echargedensity,
                                            qsd_real,qsd_imag,Espindensity,
                                            qod_real,qod_imag,Eorbmomdensity,
                                            qep_real,qep_imag,Ephonon,
                                            qem_real,qem_imag,Emagmom,
                                            qes_real,qes_imag,Espin,
                                            qel_real,qel_imag,Eorbmom,
                                            dimA,Tau,i,En(i),ini,inputpars,hkl,md,do_verbose,calc_rixs,do_phonon,QQ);
                     }
                     else
                     {ints(i)=-1;intsbey(i)=-1;intsP(i)=-1;
                     }
#endif
                      //printout rectangular function to .mcdisp.
	             if(calc_rixs){
                          // determine Isp Ipp Ips Isp  from chi(9x9 matrix)
                         // here calculate azimuth dependence of I and maximize ....
                         // ...
                                   double Iss=-1,azss=0,Isp=-1,azsp=0,Ips=-1,azps=0,Ipp=-1,azpp=0;
                                   double Irr=-1,azrr=0,Irl=-1,azrl=0,Ilr=-1,azlr=0,Ill=-1,azll=0;
                                   double Isst=0,Ispt=0,Ipst=0,Ippt=0;
                                   double Irrt=0,Irlt=0,Ilrt=0,Illt=0;
                                   ComplexVector eis(1,3),eos(1,3),eip(1,3),eop(1,3);
                                   ComplexVector eir(1,3),eor(1,3),eil(1,3),eol(1,3);
                       double daz=PI/90;if (calc_rixs==2)daz=epsilon*PI/180;
                       double azmin=0.0,azmax=2*PI;if (calc_rixs==3){azmin=epsilon*PI/180;azmax=azmin;}
                       if (En(i)!=-DBL_MAX && // Matrix Ac is not +ve definite, and this eigenvalue is not real - don't do intensity
                           En(i)<=ini.emax&&En(i)>=ini.emin){
                         for(double azimuth=azmin;azimuth<=azmax&&ints(i)>-1;azimuth+=daz)                             
                              { calc_eps(eis,eip,eir,eil,eos,eop,eor,eol,ini,azimuth,qijk,hkl, inputpars.cs.abc,QQ,En(i));
                                // eis,p and eos,p are polarisation vectors for sigma/pi plarisation in terms of
                                // eir,l and eor,l are polarisation vectors for righ/left circular plarisation in terms of
                                // the ijk coordinate system ijk form an euclidian righthanded 
                                //coordinate system j||b, k||(a x b) and i normal to j and k,
                                // where abc denot the crystal lattice vectors as defined in mcphas.j
                                Isst=calc_irix(eis,eos,chi);if(Isst>Iss){Iss=Isst;azss=azimuth*180/PI;}
                                Ispt=calc_irix(eis,eop,chi);if(Ispt>Isp){Isp=Ispt;azsp=azimuth*180/PI;}
                                Ipst=calc_irix(eip,eos,chi);if(Ipst>Ips){Ips=Ipst;azps=azimuth*180/PI;}
                                Ippt=calc_irix(eip,eop,chi);if(Ippt>Ipp){Ipp=Ippt;azpp=azimuth*180/PI;}
                                Irrt=calc_irix(eir,eor,chi);if(Irrt>Irr){Irr=Irrt;azrr=azimuth*180/PI;}
                                Irlt=calc_irix(eir,eol,chi);if(Irlt>Irl){Irl=Irlt;azrl=azimuth*180/PI;}
                                Ilrt=calc_irix(eil,eor,chi);if(Ilrt>Ilr){Ilr=Ilrt;azlr=azimuth*180/PI;}
                                Illt=calc_irix(eil,eol,chi);if(Illt>Ill){Ill=Illt;azll=azimuth*180/PI;}
                               if(calc_rixs==2){ini.print_usrdefcols(foutqei,qijk,qincr,q);
                                                fprintf (foutqei, "%4.4g %4.4g %4.4g  %4.4g %4.4g           ",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)),
                                                         myround(QQ),myround(En(i)));
                                                fprintf (foutqei, " %5.4E %3.0f    %5.4E %3.0f    %5.4E %3.0f    %5.4E %3.0f",myround(1e-8,Isst),myround(1e-8,azimuth*180/PI),myround(1e-8,Ispt),myround(1e-8,azimuth*180/PI),myround(1e-8,Ipst),myround(1e-8,azimuth*180/PI),myround(1e-8,Ippt),myround(1e-8,azimuth*180/PI));
                                                fprintf (foutqei, " %5.4E %3.0f    %5.4E %3.0f    %5.4E %3.0f    %5.4E %3.0f",myround(1e-8,Irrt),myround(1e-8,azimuth*180/PI),myround(1e-8,Irlt),myround(1e-8,azimuth*180/PI),myround(1e-8,Ilrt),myround(1e-8,azimuth*180/PI),myround(1e-8,Illt),myround(1e-8,azimuth*180/PI));
                                                fprintf (foutqei, "\n");
                                                }
                              }}
                              if (En(i)==-DBL_MAX) {
                              fprintf (foutqei, "#| "); 
                              ini.print_usrdefcols(foutqei,qijk,qincr,q);
                              fprintf (foutqei, "%4.4g %4.4g %4.4g  %4.4g %4.4g%si%-4.4g    ",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)),
                                                myround(QQ),myround(real(Enc(i))),(imag(Enc(i))<0)?"-":"+",myround(fabs(imag(Enc(i)))));
                              } else {
                              ini.print_usrdefcols(foutqei,qijk,qincr,q);
                              fprintf (foutqei, "%4.4g %4.4g %4.4g  %4.4g %4.4g           ",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)),
                                                myround(QQ),myround(En(i)));
                              }
                                if (En(i)!=-DBL_MAX&&En(i)<=ini.emax&&En(i)>=ini.emin){
                                 fprintf (foutqei, " %5.4E %3.0f    %5.4E %3.0f    %5.4E %3.0f    %5.4E %3.0f",myround(1e-8,Iss),myround(1e-8,azss),myround(1e-8,Isp),myround(1e-8,azsp),myround(1e-8,Ips),myround(1e-8,azps),myround(1e-8,Ipp),myround(1e-8,azpp));
                                 fprintf (foutqei, " %5.4E %3.0f    %5.4E %3.0f    %5.4E %3.0f    %5.4E %3.0f",myround(1e-8,Irr),myround(1e-8,azrr),myround(1e-8,Irl),myround(1e-8,azrl),myround(1e-8,Ilr),myround(1e-8,azlr),myround(1e-8,Ill),myround(1e-8,azll));
                                               } else {fprintf(foutqei, " -1 0  -1 0  -1 0  -1 0   -1 0  -1 0  -1 0  -1 0 ");}
                                   fprintf (foutqei, "\n");
                                                
                     }else{ 
                     if(En(i)!=-DBL_MAX) {
                    double test; // add to sta distance to nearest measured peak squared
 	              for (j1=1;4*j1<=ini.hkls[counter][0]-3;++j1)
	              {if ((test=fabs(En(i)-ini.hkls[counter][4*j1]))<dd1(j1)){dd1(j1)=test;double weight=ini.hkls[counter][4*j1+1];
                                                                               if(weight>0){dd(j1)=sqrt(weight)*test;  // weight>0
                                                                                            dd_without_antipeaks(j1)=sqrt(weight)*test;
                                                                                            dd_without_weights(j1)=test;
                                                                                            dd_without_antipeaks_weights(j1)=test;}
                                                                               if(weight==0){dd(j1)=0.0;  // weight=0
                                                                                             dd_without_antipeaks(j1)=0;
                                                                                             dd_without_weights(j1)=test;
                                                                                             dd_without_antipeaks_weights(j1)=test;}
                                                                               if(weight<0){if(fabs(test)<1/sqrt(ANTIPEAK_CUTOFF))test=1/sqrt(ANTIPEAK_CUTOFF);// prevents division by zero - antipeak cutoff
                                                                                            dd(j1)=sqrt(-weight)/test;  // weight<0
                                                                                            dd_without_antipeaks(j1)=0;
                                                                                            dd_without_weights(j1)=1/test;
                                                                                            dd_without_antipeaks_weights(j1)=0;}
                                                                               }
                       double inten;if((inten=ini.hkls[counter][4*j1+2])==0.0)inten=SMALLINT;
                       if ((test=fabs(En(i)-ini.hkls[counter][4*j1]))<dd1_int(j1)&&ints(i)+intsP(i)>inten){dd1_int(j1)=test;
                                                                               double weight=ini.hkls[counter][4*j1+1];
                                                                               if(weight>0){dd_int(j1)=sqrt(weight)*test;  // weight>0
                                                                                            dd_int_without_antipeaks(j1)=sqrt(weight)*test;
                                                                                            dd_int_without_weights(j1)=test;
                                                                                            dd_int_without_antipeaks_weights(j1)=test;}
                                                                               if(weight==0){dd_int(j1)=0.0;  // weight=0
                                                                                             dd_int_without_antipeaks(j1)=0;
                                                                                             dd_int_without_weights(j1)=test;
                                                                                             dd_int_without_antipeaks_weights(j1)=test;}
                                                                               if(weight<0){if(fabs(test)<1/sqrt(ANTIPEAK_CUTOFF))test=1/sqrt(ANTIPEAK_CUTOFF);// prevents division by zero
                                                                                            dd_int(j1)=sqrt(-weight)/test;  // weight<0
                                                                                            dd_int_without_antipeaks(j1)=0;
                                                                                            dd_int_without_weights(j1)=1/test;
                                                                                            dd_int_without_antipeaks_weights(j1)=0;}
                                                                               }
                      }

                     //if(intsbey(i)<0)intsbey(i)=-1.2;
                      fprintf (foutqom, " %6.6g",myround(intsbey(i)));
                      ini.print_usrdefcols(foutqei,qijk,qincr,q);
                      fprintf (foutqei, "%6.6g %6.6g %6.6g  %6.6g %6.6g  %6.6g  %6.6g %6.6g  ",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)),
                                         myround(QQ),myround(En(i)),myround(1e-8,ints(i)),myround(1e-8,intsbey(i)),myround(1e-8,intsP(i)));
	             } else {
                      fprintf (foutqei, "#| "); 
                      ini.print_usrdefcols(foutqei,qijk,qincr,q);
                      fprintf (foutqei, "%6.6g %6.6g %6.6g  %6.6g %6.6g%si%-6.6g  -1     -1    -1     ",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)),
                                         myround(QQ),myround(real(Enc(i))),(imag(Enc(i))<0)?"-":"+",myround(fabs(imag(Enc(i)))));
                     } 
                  if (En(i)!=-DBL_MAX&&En(i)<=ini.emax&&En(i)>=ini.emin){   
                       DMDtotint+=ints(i);DMDtotintbey+=intsbey(i);
                       chitot=chitot+chi;chitotbey=chitotbey+chibey;                    
                       switch(ini.outS)
                         {case 0: break;
                          case 1: for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutqei," %4.4g %4.4g ",real(chi(i1,j1)),imag(chi(i1,j1)));break;
                          case 2: for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutqei," %4.4g %4.4g ",real(chibey(i1,j1)),imag(chibey(i1,j1)));break;     
                          case 3: rottouvw(chi,ini,inputpars.cs.abc,counter);for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutqei," %4.4g %4.4g ",real(chi(i1,j1)),imag(chi(i1,j1)));break;
                          case 4: rottouvw(chibey,ini,inputpars.cs.abc,counter);for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutqei," %4.4g %4.4g ",real(chibey(i1,j1)),imag(chibey(i1,j1)));break;     
                          case 5: for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutqei," %4.4g %4.4g ",real(chi(i1,j1)),imag(chi(i1,j1)));break;
                          case 6: rottouvw(chi,ini,inputpars.cs.abc,counter);for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutqei," %4.4g %4.4g ",real(chi(i1,j1)),imag(chi(i1,j1)));break;
                         } }
                       fprintf(foutqei,"\n");
                       if(do_verbose==1){fprintf(stdout, "#level %i IdipFF= %4.4g Ibeyonddip=%4.4g Iphonon=%4.4g\n",i,ints(i),intsbey(i),intsP(i));}
                      }
                     // printout eigenvectors only if evaluated during intensity calculation...
                  if (En(i)!=-DBL_MAX&&En(i)<=ini.emax&&En(i)>=ini.emin){ 

if(ini.calculate_chargedensity_oscillation)print_ev(foutqee,i,ini,hkl,QQ,En,ints,intsbey,intsP,qee_real,qee_imag,q);
if(ini.calculate_spindensity_oscillation)print_ev(foutqsd,i,ini,hkl,QQ,En,ints,intsbey,intsP,qsd_real,qsd_imag,q);
if(ini.calculate_orbmomdensity_oscillation)print_ev(foutqod,i,ini,hkl,QQ,En,ints,intsbey,intsP,qod_real,qod_imag,q);
if(ini.calculate_phonon_oscillation)print_ev(foutqep,i,ini,hkl,QQ,En,ints,intsbey,intsP,qep_real,qep_imag,q);
if(ini.calculate_magmoment_oscillation)print_ev(foutqem,i,ini,hkl,QQ,En,ints,intsbey,intsP,qem_real,qem_imag,q);
if(ini.calculate_spinmoment_oscillation)print_ev(foutqes,i,ini,hkl,QQ,En,ints,intsbey,intsP,qes_real,qes_imag,q);
if(ini.calculate_orbmoment_oscillation)print_ev(foutqel,i,ini,hkl,QQ,En,ints,intsbey,intsP,qel_real,qel_imag,q);
                          }
#ifdef _THREADS
                     }
                   i=oldi;
#endif
		   }
#ifdef _THREADS
                  #undef qee_real
                  #undef qee_imag
                  #undef qsd_real
                  #undef qsd_imag
                  #undef qod_real
                  #undef qod_imag
                  #undef qep_real
                  #undef qep_imag
                  #undef qem_real
                  #undef qem_imag
                  #undef qes_real
                  #undef qes_imag
                  #undef qel_real
                  #undef qel_imag
                  #undef chi
                  #undef chibey
                  for (ithread=0; ithread<NUM_THREADS; ithread++) 
                  {
                     delete thrdat.chi[ithread]; delete thrdat.chibey[ithread];delete thrdat.chiPhon[ithread]; 
                     delete thrdat.pol[ithread]; delete thrdat.md[ithread]; 
                     if(ini.calculate_chargedensity_oscillation){delete thrdat.qee_real[ithread]; delete thrdat.qee_imag[ithread];delete thrdat.Echargedensity[ithread]; }
                     if(ini.calculate_spindensity_oscillation){delete thrdat.qsd_real[ithread]; delete thrdat.qsd_imag[ithread];delete thrdat.Espindensity[ithread]; }
                     if(ini.calculate_orbmomdensity_oscillation){delete thrdat.qod_real[ithread]; delete thrdat.qod_imag[ithread];delete thrdat.Eorbmomdensity[ithread];} 
                     if(ini.calculate_phonon_oscillation){delete thrdat.qep_real[ithread]; delete thrdat.qep_imag[ithread];delete thrdat.Ephonon[ithread]; }
                     if(ini.calculate_magmoment_oscillation){delete thrdat.qem_real[ithread]; delete thrdat.qem_imag[ithread];delete thrdat.Emagmom[ithread];} 
                     if(ini.calculate_spinmoment_oscillation){delete thrdat.qes_real[ithread]; delete thrdat.qes_imag[ithread];delete thrdat.Espin[ithread]; }
                     if(ini.calculate_orbmoment_oscillation){delete thrdat.qel_real[ithread]; delete thrdat.qel_imag[ithread];delete thrdat.Eorbmom[ithread]; }
                    // delete thrdat.Tau[ithread];
                     delete tin[ithread]; 
                  }
                  delete[] thrdat.Echargedensity;  
                  delete[] thrdat.Espindensity;  
                  delete[] thrdat.Eorbmomdensity;  
                  delete[] thrdat.Ephonon;  
                  delete[] thrdat.Emagmom;  
                  delete[] thrdat.Espin;  
                  delete[] thrdat.Eorbmom;  
                  delete[] thrdat.Tau;
                  delete[] thrdat.chi; delete[] thrdat.chibey;delete[] thrdat.chiPhon;
                  delete[] thrdat.pol; 
                  delete[] thrdat.qee_real; delete[] thrdat.qee_imag; 
                  delete[] thrdat.qsd_real; delete[] thrdat.qsd_imag; 
                  delete[] thrdat.qod_real; delete[] thrdat.qod_imag; 
                  delete[] thrdat.qep_real; delete[] thrdat.qep_imag; 
                  delete[] thrdat.qem_real; delete[] thrdat.qem_imag; 
                  delete[] thrdat.qes_real; delete[] thrdat.qes_imag; 
                  delete[] thrdat.qel_real; delete[] thrdat.qel_imag; 
#endif
if(!calc_rixs){ini.print_usrdefcols(foutdstot,qijk,qincr,q);
               fprintf (foutdstot, "%4.4g %4.4g  %4.4g %4.4g %4.4g",hkl(1),hkl(2),hkl(3),DMDtotint,DMDtotintbey);
               switch(ini.outS)
                         {case 0: break;
                          case 1: for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutdstot," %4.4g %4.4g ",real(chitot(i1,j1)),imag(chitot(i1,j1)));break;
                          case 2: for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutdstot," %4.4g %4.4g ",real(chitotbey(i1,j1)),imag(chitotbey(i1,j1)));break;     
                          case 3: rottouvw(chitot,ini,inputpars.cs.abc,counter);for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutdstot," %4.4g %4.4g ",real(chitot(i1,j1)),imag(chitot(i1,j1)));break;
                          case 4: rottouvw(chitotbey,ini,inputpars.cs.abc,counter);for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutdstot," %4.4g %4.4g ",real(chitotbey(i1,j1)),imag(chitotbey(i1,j1)));break;     
                          case 5: for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutdstot," %4.4g %4.4g ",real(chitot(i1,j1)),imag(chitot(i1,j1)));break;
                          case 6: rottouvw(chitot,ini,inputpars.cs.abc,counter);for(i1=1;i1<=3;++i1)for(j1=1;j1<=3;++j1) fprintf(foutdstot," %4.4g %4.4g ",real(chitot(i1,j1)),imag(chitot(i1,j1)));break;
                         } 
                      
    sta+=dd*dd;sta_int+=dd_int*dd_int;
    sta_without_antipeaks+=dd_without_antipeaks*dd_without_antipeaks;
    sta_int_without_antipeaks+=dd_int_without_antipeaks*dd_int_without_antipeaks;
    sta_without_weights+=dd_without_weights*dd_without_weights;
    sta_int_without_weights+=dd_int_without_weights*dd_int_without_weights;
    sta_without_antipeaks_weights+=dd_without_antipeaks_weights*dd_without_antipeaks_weights;
    sta_int_without_antipeaks_weights+=dd_int_without_antipeaks_weights*dd_int_without_antipeaks_weights;
              }
              //initialize output file for display
            snprintf(filename,MAXNOFCHARINLINE,"./results/.%smcdisp.qom",ini.prefix);fout1 = fopen_errchk (filename,"w");
            fprintf (fout1, "#%s ",MCDISPVERSION);
            curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout1);
            fprintf (fout1, "#displayytext=I(barns/meV/sr/f.u.)\n");
            fprintf (fout1, "#displayxtext=E(meV)\n");
            fprintf (fout1, "#displaytitle=(%4.4f %4.4f %4.4f) blue: DMD_Dipapprox red: DMD_exact green: Minv_Dipapprox\n",hkl(1),hkl(2),hkl(3));
            //fprintf (fout1,"#Ha[T] Hb[T] Hc[T] T[K] h k l  energies[meV] intensities(dip approx for FF) [barn/meV/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)\n");
		     if (do_Erefine==0) epsilon=(Max(En)-Min(En)+0.001)/100;
		    // if (epsilon<=0) epsilon=0.1;
                  for (i=1;i<=dimA;++i) if(En(i)!=-DBL_MAX)
		    { 
		     if (ints(i)>SMALLINT)  // draw triangles to show calculated intensity
		      {for (E=0;E<=ints(i)/fabs(epsilon);E+=ints(i)/2/fabs(epsilon)/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	                fprintf (fout1, " %4.4g %4.4g %4.4g 0\n",myround(En(i)-fabs(epsilon)+E*fabs(epsilon)*fabs(epsilon)/ints(i)),myround(E),myround(En(i)));
		       }
		       for (E=ints(i)/fabs(epsilon);E>=0;E-=ints(i)/2/fabs(epsilon)/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	                fprintf (fout1, " %4.4g %4.4g %4.4g 0\n",myround(En(i)+fabs(epsilon)-E*fabs(epsilon)*fabs(epsilon)/ints(i)),myround(E),myround(En(i)));
		       }
                       fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	               fprintf (fout1, " %4.4g 0 %4.4g 0\n",myround(En(i)+fabs(epsilon)),myround(En(i)));
		      }
		     if (intsbey(i)>SMALLINT)  // draw triangles to show calculated intensity
		      {for (E=0;E<=intsbey(i)/fabs(epsilon);E+=intsbey(i)/2/fabs(epsilon)/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	                fprintf (fout1, " %4.4g 0 %4.4g %4.4g \n",myround(En(i)),myround(En(i)-fabs(epsilon)+E*fabs(epsilon)*fabs(epsilon)/intsbey(i)),E);
		       }
		       for (E=intsbey(i)/fabs(epsilon);E>=0;E-=intsbey(i)/2/fabs(epsilon)/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	                fprintf (fout1, " %4.4g 0 %4.4g %4.4g\n",myround(En(i)),myround(En(i)+fabs(epsilon)-E*fabs(epsilon)*fabs(epsilon)/intsbey(i)),E);
		       }
                       fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	               fprintf (fout1, " %4.4g 0 %4.4g 0 \n",myround(En(i)),myround(En(i)+fabs(epsilon)));
		      }
		    }
	  fclose(fout1);   
              if(do_verbose==1){fprintf(stdout, "\n");}

//*********************************************************************		    
   // do refinement of energies by output of scattering cross section vs enrgy transfer if required
  if (do_Erefine==1&&!calc_rixs){double totint=0;
                if(do_verbose==1){fprintf(stdout, "#refining calculation with exact calculation of energy dependence of scattering cross section\n");}
          snprintf(filename,MAXNOFCHARINLINE,"./results/.%smcdisp.dsigma",ini.prefix);foutds1 = fopen_errchk (filename,"w");
          fprintf (foutds1, "#{%s ",MCDISPVERSION);
          curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutds1);
          fprintf (foutds1, "#Scattering Cross Section \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energy[meV] dsigma/dOmegadE'[barn/mev/sr/f.u.] (dipolar approx for FF) chixxr chixxi  chixyr chixyi chixzr chixri chiyxr chiyxi chiyyr chiyyi chiyzr chiyzi chizxr chizxi chizyr chizyi chizzr chizzi (1/meV/f.u.) f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
          ComplexMatrix ch(1,md.nofcomponents,1,md.nofcomponents);
#ifdef _THREADSREFINE
          ComplexMatrix *chpointer[nofEstps+1];
          for(int Estp=0;Estp<nofEstps;Estp++){chpointer[Estp]=new ComplexMatrix(1,md.nofcomponents,1,md.nofcomponents);}
          thrdat.ch = new ComplexMatrix*[NUM_THREADS];                  
          thrdat.hkl = hkl; thrdat.q = q; thrdat.thread_id = -1;
          for (ithread=0; ithread<NUM_THREADS; ithread++) 
          {  tin[ithread] = new intcalcapr_input(dimA,ithread,1,do_verbose,calc_rixs,do_phonon,0.); 
             tin[ithread]->epsilon=fabs(epsilon); thrdat.J[ithread] = new jq(J);thrdat.md[ithread] = new mdcf(md,1);
             thrdat.ch[ithread] = new ComplexMatrix(1,md.nofcomponents,1,md.nofcomponents);
          }
          ithread=0; int oldEstp=0.;
          Vector vIntensity(1,(int)((ini.emax-ini.emin)/(fabs(epsilon)/2)+1)); int iE=1;
#endif
	     double intensity;
#ifdef _THREADSREFINE
                     tin[ithread]->Estp=0; tin[ithread]->iE=iE;
                     #if defined  (__linux__) || defined (__APPLE__)
                     rc = pthread_create(&threads[ithread], &attr, intcalc_Erefine, (void *) tin[ithread]); rc = pthread_join(threads[ithread], &status); 
                     if(rc) { printf("Error return code %i from erefine joining thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
                     #else
                     threads[ithread] = CreateThread(NULL, 0, intcalc_Erefine, (void *) tin[ithread], 0, &tid[ithread]);
                     if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %lu from erefine thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
                     if(WaitForSingleObject(threads[ithread],INFINITE)==0xFFFFFFFF) { printf("Error in waiting for erefine thread %i to end\n",ithread+1); exit(EXIT_FAILURE); }
                     CloseHandle(threads[ithread]);
                     #endif
                     intensity=tin[ithread]->intensity;
                     ch=(*thrdat.ch[ithread]);
#else
//		     intensity=intcalc(ch,dimA,ini.emin,ini,inputpars,J,q,hkl,md,do_verbose,fabs(epsilon));   
		     intensity=intcalc_Erefine(ch,0,ini,inputpars,J,q,hkl,md,do_verbose,fabs(epsilon));   
#endif
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	             fprintf (foutds1, " %4.4g %4.4g   ",ini.emin,myround(intensity));
                     for (int ii=1;ii<=ch.Rhi();++ii)for (int jj=1;jj<=ch.Chi();++jj)fprintf(foutds1,"%4.4g %4.4g  ",myround(real(ch(ii,jj))),myround(imag(ch(ii,jj))));
                     fprintf(foutds1,"\n");

#ifdef _THREADSREFINE
                     tin[ithread]->Estp=nofEstps-1; tin[ithread]->iE=iE;
                     #if defined  (__linux__) || defined (__APPLE__)
                     rc = pthread_create(&threads[ithread], &attr, intcalc_Erefine, (void *) tin[ithread]); rc = pthread_join(threads[ithread], &status); 
                     if(rc) { printf("Error return code %i from erefine thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
                     #else
                     threads[ithread] = CreateThread(NULL, 0, intcalc_Erefine, (void *) tin[ithread], 0, &tid[ithread]);
                     if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %lu from erefine thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
                     if(WaitForSingleObject(threads[ithread],INFINITE)==0xFFFFFFFF) { printf("Error in waiting for erefine thread %i to end\n",ithread+1); exit(EXIT_FAILURE); }
                     CloseHandle(threads[ithread]); 
                    #endif
                     intensity=tin[ithread]->intensity;
                     ch=(*thrdat.ch[ithread]);
#else
//		     intensity=intcalc(ch,dimA,ini.emax,ini,inputpars,J,q,hkl,md,do_verbose,fabs(epsilon));   
		     intensity=intcalc_Erefine(ch,nofEstps-1,ini,inputpars,J,q,hkl,md,do_verbose,fabs(epsilon));   
#endif
                    fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	             fprintf (foutds1, " %4.4g %4.4g ",ini.emax,myround(intensity));
                     for (int ii=1;ii<=ch.Rhi();++ii)for (int jj=1;jj<=ch.Chi();++jj)fprintf(foutds1,"%4.4g %4.4g  ",myround(real(ch(ii,jj))),myround(imag(ch(ii,jj))));
                     fprintf(foutds1,"\n");
	  fclose(foutds1);
#ifdef _THREADSREFINE
	  for(int Estp=0;Estp<nofEstps;Estp+=NUM_THREADS)
#else
	  for(int Estp=0;Estp<nofEstps;++Estp)
#endif
	   {E=ini.emin+Estp*fabs(epsilon)/2;
#ifndef _THREADSREFINE
//		    intensity=intcalc(ch,dimA,E,ini,inputpars,J,q,hkl,md,do_verbose,fabs(epsilon));   
		    intensity=intcalc_Erefine(ch,Estp,ini,inputpars,J,q,hkl,md,do_verbose,fabs(epsilon));   
		     totint+=intensity*fabs(epsilon)/2;
          snprintf(filename,MAXNOFCHARINLINE,"./results/.%smcdisp.dsigma",ini.prefix);foutds1 = fopen_errchk (filename,"a");
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	             fprintf (foutds1, " %4.4g %4.4g ",myround(E),myround(intensity));
                     for (int ii=1;ii<=ch.Rhi();++ii)for (int jj=1;jj<=ch.Chi();++jj)fprintf(foutds1,"%4.4g %4.4g  ",myround(real(ch(ii,jj))),myround(imag(ch(ii,jj))));
                     fprintf(foutds1,"\n");
          fclose(foutds1);	   
#else
                   oldEstp=Estp;
                     for(ithread=0; ithread<NUM_THREADS; ithread++)
                     {
                        Estp=oldEstp+ithread; if(Estp>=nofEstps) break;
                        tin[ithread]->Estp=Estp; tin[ithread]->iE=iE++;
                        #if defined  (__linux__) || defined (__APPLE__)
                        rc = pthread_create(&threads[ithread], &attr, intcalc_Erefine, (void *) tin[ithread]);
                        if(rc) { printf("Error return code %i from erefine thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
                        #else
                        threads[ithread] = CreateThread(NULL, 0, intcalc_Erefine, (void *) tin[ithread], 0, &tid[ithread]);
                        if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %lu from erefine thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
                        #endif
                        num_threads_started = ithread+1;
                     }
                     #if defined  (__linux__) || defined (__APPLE__)
                     for(int th=0; th<num_threads_started; th++)
                        rc = pthread_join(threads[th], &status);
                     #else
                     if(num_threads_started>0){retval=WaitForMultipleObjects(num_threads_started,threads,TRUE,INFINITE);
                     if(retval<WAIT_OBJECT_0||retval>WAIT_OBJECT_0+num_threads_started-1){printf("Error waitformultipleobjects erfine\n"); exit(EXIT_FAILURE); }
                     for(int th=0; th<num_threads_started; th++)CloseHandle(threads[th]);}
                        #endif
                     snprintf(filename,MAXNOFCHARINLINE,"./results/.%smcdisp.dsigma",ini.prefix);foutds1 = fopen_errchk (filename,"a");
                      for(ithread=0; ithread<NUM_THREADS; ithread++) {vIntensity(tin[ithread]->iE) = tin[ithread]->intensity;
                                                                     (*chpointer[tin[ithread]->iE-1])  =(*thrdat.ch[ithread]); 
                      fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	              fprintf (foutds1, " %4.4g %4.4g\n ",myround(ini.emin+tin[ithread]->Estp*fabs(epsilon)/2),myround(tin[ithread]->intensity));
                                                                      }
                     fclose(foutds1);	   
                     Estp=oldEstp;
	   }
	  iE=1; for(int Estp=0;Estp<nofEstps;Estp++)
	   {E=ini.emin+Estp*fabs(epsilon)/2;
                     ch=(*chpointer[iE-1]);
                     intensity = vIntensity(iE++); totint+=intensity*fabs(epsilon)/2;
#endif
                     ini.print_usrdefcols(foutds,qijk,qincr,q);
                     fprintf (foutds, "%4.4g %4.4g  %4.4g ",myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	             fprintf (foutds, " %4.4g %4.4g ",myround(E),myround(intensity));
                     for (int ii=1;ii<=ch.Rhi();++ii)for (int jj=1;jj<=ch.Chi();++jj)fprintf(foutds,"%4.4g %4.4g  ",myround(real(ch(ii,jj))),myround(imag(ch(ii,jj))));
                     fprintf(foutds,"\n");
	   }
#ifdef _THREADSREFINE
          for (ithread=0; ithread<NUM_THREADS; ithread++) 
          {
             delete thrdat.md[ithread];delete thrdat.J[ithread]; delete tin[ithread];
             delete thrdat.ch[ithread];
          }
          delete[] thrdat.ch; 
          for(int Estp=0;Estp<nofEstps;Estp++)delete chpointer[Estp];
#endif
	             fprintf (foutdstot, " %4.4g ",totint);
                     }  // do_Erefine
//*********************************************************************		    

   if(!calc_rixs)fprintf (foutdstot, "\n");              
   fprintf (foutqom, "\n");
   } // do jqfile
} // next hkl
#ifdef _THREADS
   for (ithread=0; ithread<NUM_THREADS; ithread++) 
   {  delete thrdat.ini[ithread]; 
      delete thrdat.inputpars[ithread]; 
   }
   delete[] thrdat.inputpars; delete[] thrdat.md; delete[] thrdat.ini;delete[] thrdat.J;
#endif

                                     

    if (do_jqfile) 
     {fprintf(jqfile,"#!the largest eigenvalue of J(q) is jqmax=%g meV at hmax=%g kmax=%g lmax=%g \n",jqmax,hmax,kmax,lmax);
     fprintf(jqfile,"#!for the first q vector in the list jq0=%g meV at h0=%g k0=%g l0=%g \n",jq0,ini.hkls[firstcounter][1],ini.hkls[firstcounter][2],ini.hkls[firstcounter][3]);
      fprintf(jqfile,"#it follows the standard deviation sta defined as:\n");
      fprintf(jqfile,"#A)the sum of squared differences between the highest eigenvalue\n");
      fprintf(jqfile,"#of a q vector and that of the first q-vector in the list in mcdisp.par.\n");
      fprintf(jqfile,"#only those eigenvalues are taken into account in the sum, which are larger\n");
      fprintf(jqfile,"#than that of the first q-vector in the list in mcdisp.par - this is usefule\n");
      fprintf(jqfile,"#for obtaining an exchange interaction with maximum at the first q-vector\n");
      fprintf(jqfile,"#in the list in mcdisp.par\n");
      fprintf(jqfile,"# ... if the first q vector has the largest eigenvalue, then sta is negative and contains the\n");
      fprintf(jqfile,"#distance to the closest eigenvalue\n");
      fprintf(jqfile,"#!sta=%g\n",jqsta);
      fprintf(jqfile,"#B)another standard deviation is given below: calculated as squared sum of differences between\n");
      fprintf(jqfile,"#the highest eigenvalue of J(Q) and energies in column 4 of mcdisp.par, if column 5 and 6  \n");
      fprintf(jqfile,"#in mcdisp.par contain values, then these are compared to the other eigenvalues of J(Q)\n");
      fprintf(jqfile,"#!sta4=%g\n",jqsta_int);
      fprintf(jqfile,"#C)another standard deviation is given below: calculated as squared distance of \n");
      fprintf(jqfile,"#the q-vectors q0 and qmax, i.e. (hmax-h0)^2+(kmax-k0)^2+(lmax-l0)^2 \n");
double staq=(hmax-ini.hkls[firstcounter][1])*(hmax-ini.hkls[firstcounter][1])+(kmax-ini.hkls[firstcounter][2])*(kmax-ini.hkls[firstcounter][2])+(lmax-ini.hkls[firstcounter][3])*(lmax-ini.hkls[firstcounter][3]);
      fprintf(jqfile,"#!staq=%g\n",staq);
      fprintf(jqfile,"#D)other standard deviations are given below: products of  sta / sta4 and staq \n");
      fprintf(jqfile,"#!staxstaq=%g\n",jqsta*staq);
      fprintf(jqfile,"#!sta4xstaq=%g\n",jqsta_int*staq);
      fprintf(jqfile,"#\n#if in mcdisp.par after the hkl of the first q vector a number is given, then\n");
      fprintf(jqfile,"#scaled interaction parameters are computed such that their highest eigenvalue\n");
      fprintf(jqfile,"#agrees with this number. Scaled parameters are saved in results/mcdisp_scaled.j\n");
      fprintf(jqfile,"#and here follow standard deviations computed with the scaled parameters\n");
      if(scalefactor!=1.0){fprintf(jqfile,"#!sta_scaled=%g\n",jqsta_scaled);
                           fprintf(jqfile,"#!sta4_scaled=%g\n",jqsta_int_scaled);
                           fprintf(jqfile,"#!sta_scaledxstaq=%g\n",jqsta_scaled*staq);
                           fprintf(jqfile,"#!sta4_scaledxstaq=%g\n",jqsta_int_scaled*staq);
                           fprintf(jqfile,"#!scalefactor=%g\n",scalefactor);
                           }
       fclose(jqfile);
      if(scalefactor!=1.0){inputpars.scale(scalefactor);
                           snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp_scaled.j",ini.prefix);
                          printf("# saving  %s\n",filename);
                          jqfile = fopen_errchk (filename,"w");
                           inputpars.save(jqfile,0);
                           fclose(jqfile);
                           }
     }
    else
     {
      if(!calc_rixs){staout(foutqom,sta,sta_int,sta_without_antipeaks,sta_int_without_antipeaks,sta_without_weights,sta_int_without_weights,sta_without_antipeaks_weights,sta_int_without_antipeaks_weights);
      staout(foutqei,sta,sta_int,sta_without_antipeaks,sta_int_without_antipeaks,sta_without_weights,sta_int_without_weights,sta_without_antipeaks_weights,sta_int_without_antipeaks_weights);
      staout(stdout,sta,sta_int,sta_without_antipeaks,sta_int_without_antipeaks,sta_without_weights,sta_int_without_weights,sta_without_antipeaks_weights,sta_int_without_antipeaks_weights);
      if (do_Erefine==1){fclose(foutds);}
      fclose(foutdstot);
                    }
      fclose(foutqom);
      fclose(foutqei);

                        if(ini.calculate_chargedensity_oscillation)fclose(foutqee);
                        if(ini.calculate_spindensity_oscillation)fclose(foutqsd);
                        if(ini.calculate_orbmomdensity_oscillation)fclose(foutqod);
                        if(ini.calculate_phonon_oscillation)fclose(foutqep);
                        if(ini.calculate_magmoment_oscillation)fclose(foutqem);
                        if(ini.calculate_spinmoment_oscillation)fclose(foutqes);
                        if(ini.calculate_orbmoment_oscillation)fclose(foutqel);
                        
     } 
}

//*************************************************************************************************
// main program
int main (int argc, char **argv)
{std::clock_t startcputime = std::clock();
 int i,do_Erefine=0,do_jqfile=0,do_verbose=0,maxlevels=10000000,do_createtrs=0;
 int do_ignore_non_hermitian_matrix_error=0;
 int do_readtrs=0,calc_beyond=1,calc_rixs=0;
 char spinfile [MAXNOFCHARINLINE]; //default spin-configuration-input file
  snprintf(spinfile,MAXNOFCHARINLINE,"mcdisp.mf");
 const char * filemode="w";
 char prefix [MAXNOFCHARINLINE];prefix[0]='\0';
 double epsilon=0.05; //imaginary part of omega to avoid divergence
 double minE=-100000.0,maxE=+100000.0,pinit=SMALL_PROBABILITY,ninit=1e10;
 fprintf(stderr,"#***********************************************************************\n");
 fprintf(stderr,"#*\n");
 fprintf(stderr,"#* mcdisp - program to calculate the dispersion of magnetic excitations\n");
 fprintf(stderr,"#*\n");
 fprintf(stderr,"#* reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
 fprintf(stderr,"#*            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
 fprintf(stderr,"#***********************************************************************\n\n");


//***************************************************************************************
// check command line parameters 
//***************************************************************************************
for (i=1;i<=argc-1;++i){
   if(strcmp(argv[i],"-r")==0) {do_Erefine=1; if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -r needs argument epsilon\n");exit(EXIT_FAILURE);}
		                                                epsilon=strtod(argv[i+1],NULL);++i;
							        fprintf(stdout,"#epsilon= %g\n",epsilon);
				     }		
     else {if(strcmp(argv[i],"-xaf")==0) {calc_rixs=3;calc_beyond=0; // rixs with azimuth fixed
                                                  if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -xaf needs argument(s)\n");exit(EXIT_FAILURE);}
		                                  epsilon=strtod(argv[i+1],NULL);++i; // use epsilon to convey azimuth
						  fprintf(stdout,"#maximum number of single ion excitations taken into account (starting with lowest energy): %i\n",maxlevels);
					         }
      else {if(strcmp(argv[i],"-xa")==0) {calc_rixs=2;calc_beyond=0; // rixs with azimuth dependence in steps
                                                  if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -xa needs argument(s)\n");exit(EXIT_FAILURE);}
		                                  epsilon=strtod(argv[i+1],NULL);++i; // use epsilon to convey azimuth step
						  fprintf(stdout,"#maximum number of single ion excitations taken into account (starting with lowest energy): %i\n",maxlevels);
					         }
       else {if(strcmp(argv[i],"-x")==0) {calc_rixs=1;calc_beyond=0;}  // rixs without azimuth dep .. Irixs max only
        else {if(strcmp(argv[i],"-d")==0) {calc_beyond=0;}
         else {if(strcmp(argv[i],"-jq")==0) {do_jqfile=1;minE=SMALL_QUASIELASTIC_ENERGY;maxlevels=1;}
          else {if(strcmp(argv[i],"-jqe")==0) {do_jqfile=2;minE=SMALL_QUASIELASTIC_ENERGY;maxlevels=1;}
           else {if(strcmp(argv[i],"-t")==0) do_readtrs=1;       
            else {if(strcmp(argv[i],"-c")==0) do_createtrs=1;       
             else {if(strcmp(argv[i],"-A")==0) filemode="A";       
              else {if(strcmp(argv[i],"-a")==0) filemode="a";       
               else {if(strcmp(argv[i],"-v")==0||strcmp(argv[i],"-verbose")==0) do_verbose=1;       
                else {if(strcmp(argv[i],"-max")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -max needs argument(s)\n");exit(EXIT_FAILURE);}
  		                                  maxlevels=(int)strtod(argv[i+1],NULL);++i;  
						  fprintf(stdout,"#maximum number of single ion excitations taken into account (starting with lowest energy): %i\n",maxlevels);
  					         }       
                 else {if(strcmp(argv[i],"-maxE")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -maxE needs argument(s)\n");exit(EXIT_FAILURE);}
		                                  maxE=strtod(argv[i+1],NULL);++i;
 						  fprintf(stdout,"#maximum Energy of single ion excitations taken into account: %g\n",maxE);
  					         }       
                  else {if(strcmp(argv[i],"-minE")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -minE needs argument(s)\n");exit(EXIT_FAILURE);}
 		                                  minE=strtod(argv[i+1],NULL);++i;
						  fprintf(stdout,"#minimum Energy of single ion excitations taken into account: %g\n",minE);
					         }
                   else {if(strcmp(argv[i],"-ninit")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -ninit needs argument(s)\n");exit(EXIT_FAILURE);}
 		                                  ninit=strtod(argv[i+1],NULL);++i;
						  fprintf(stdout,"#maximum number of lowest lying initial states to be taken into account in single ion excitations: %g\n",ninit);
					         }
                    else {if(strcmp(argv[i],"-pinit")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -pinit needs argument(s)\n");exit(EXIT_FAILURE);}
 		                                  pinit=strtod(argv[i+1],NULL);++i;
						  fprintf(stdout,"#minimum population of initial state for single ion excitations to be taken into account: %g\n",pinit);
					         }
                     else {if(strcmp(argv[i],"-prefix")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -prefix needs argument(s)\n");exit(EXIT_FAILURE);}
  		                                  strcpy(prefix,argv[i+1]);snprintf(spinfile,MAXNOFCHARINLINE,"%smcdisp.mf",prefix);++i;
 						  fprintf(stdout,"#prefix for reading parameters from mcdisp.par and for ouput filenames: %s\n",prefix);
 					         }
                      else {if(strcmp(argv[i],"-ignore_non_hermitian_matrix_error")==0) {do_ignore_non_hermitian_matrix_error=1;
 						  fprintf(stdout,"#ignoring not positive definite matrices\n");
 					         }
                       else {if(strncmp(argv[i],"-h",2)==0) {errexit();}
              	       else{strcpy(spinfile,argv[i]);}
                          } // help
                         } // do_ignore_non_hermitian_matrix_error
                        } // prefi
	 	       } // pinit
	 	      } // ninit
	 	     } // minE
	 	    } //max          
	 	   } // -v
	  	  } // -a
	         } // -A
	        } // -c
	       } // -t
              } // -jqe
 	     } // -jq 
           }	
          }
         }
        }
      } // xaf
    }
  // as class load  parameters from file
  par inputpars("./mcphas.j",do_verbose);
  
  inimcdis ini("mcdisp.par",spinfile,prefix,do_jqfile,inputpars.cs.abc);
  if(ini.nofcomponents!=inputpars.cs.nofcomponents){fprintf(stderr,"Error mcdisp: number of components read from mcdisp.par (%i) and mcphas.j (%i) not equal\n",ini.nofcomponents,inputpars.cs.nofcomponents);exit(EXIT_FAILURE);}
  if(do_Erefine&&calc_rixs){fprintf(stderr,"Error mcdisp: Option -r not possible in combination with option -x -xa -xaf\n");exit(EXIT_FAILURE);}
  if(do_jqfile&&do_readtrs){fprintf(stderr,"Error mcdisp: Option -t and -jq are cannot be used at the same time\n");exit(EXIT_FAILURE);}
  if(ini.nofatoms!=inputpars.cs.nofatoms){fprintf(stderr,"Error mcdisp: number of atoms in crystal unit cell read from mcdisp.par (%i) and mcphas.j (%i) not equal\n",ini.nofatoms,inputpars.cs.nofatoms);exit(EXIT_FAILURE);}
  strcpy(prefix,"./results/_");strcpy(prefix+11,ini.prefix);  inputpars.save_sipfs(prefix); 
  strcpy(prefix+11+strlen(ini.prefix),"mcdisp.j");            inputpars.save(prefix,0);

int do_phonon=1;
//calculate dispersion and save to files
dispcalc(ini,inputpars,calc_rixs,do_phonon,calc_beyond,do_Erefine,do_jqfile,do_createtrs,do_readtrs,do_verbose,do_ignore_non_hermitian_matrix_error,maxlevels,minE,maxE,ninit,pinit,epsilon,filemode);
  
 printf("#RESULTS saved in directory ./results/  - files:\n");
  if(do_jqfile){
   printf("#  %smcdisp.jq  - Fourier Transfor J(Q) of the Interaction Parmeters\n",ini.prefix);
  }else{
  if(calc_rixs){printf("#  %smcdisp.qex  - T,H,qvector vs energies and resonant inelastic X-ray (RIXS) intensities\n",ini.prefix);}
  else{ printf("#  %smcdisp.qei  - T,H,qvector vs energies and neutron intensities\n",ini.prefix);
   printf("#  %smcdisp.dsigma.tot  - T,H,qvector vs total intensity (sum of all modes)\n",ini.prefix);
   printf("#  %smcdisp.dsigma      - (option -r) T,H,qvector,E vs intensity obtained from dyn susz\n",ini.prefix);
      }
   printf("#  %smcdisp.qom  - T,H,qvector vs all mode energies in one line\n",ini.prefix);
   printf("#  %smcdisp.qee,qsd,qod,qep,qem,qes,qel  - T,H,qvector,E vs extended eigenvectors (more components to plot observables.)\n",ini.prefix);
   printf("#  %smcdisp.trs  - single ion transitions used\n",ini.prefix);
   }
   printf("#  _%smcdisp.par - input parameters read from mcdisp.par\n",ini.prefix);
   printf("#  _%smcdisp.mf  - input parameters read from %s\n",ini.prefix,spinfile);
   printf("#  _%smcdisp.j   - input parameters read from mcphas.j\n",ini.prefix);
   printf("#  ...         - and a copy of the single ion parameter files used.\n\n");
   double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
   std::cout << "#! Finished in cputime=" << cpu_duration << " seconds [CPU Clock] " << std::endl;
   std::cout << "#! nofhkls=" << ini.nofhkls << " different q vectors calculated " << std::endl;

#ifdef _THREADS
std::cout << "#! nofthreads= " << NUM_THREADS << " threads were used in parallel processing " << std::endl;
#else
std::cout << "# mcdisp was compiled without parallel processing option " << std::endl;
#endif

   fprintf(stderr,"#************************************************************\n");
   fprintf(stderr,"#                    End of Program mcdisp\n");
   fprintf(stderr,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(stderr,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(stderr,"#************************************************************\n");

 return(0);
 
}





