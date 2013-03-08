/***********************************************************************
 *
 * mcdisp - program to calculate the dispersion of magnetic excitations
 *
 * reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751
 *            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400
 ***********************************************************************/

#define SMALL 1e-6    // deviation from single ion gap delta to take energy into account as not being equal to
                      // delta and therefore being included into output 
		      // transitions of single ions less then SMALL have in Mijkl wn/kT instead of wn-wn'
		      // !!! must match SMALL in jjjpar.cpp !!!!
#define SMALLINT 1e-4 // small intensity treshhold in barn/f.u. - only peaks larger than this threshold will
                      // be included in sta_int calculation
#define ANTIPEAK_CUTOFF 100 // if antipeaks are fitted,  sta contribution never increases above this cutoff
#define SMALLEDIF 1e-5 // small difference in calculation of transition energy
                       // used to give error if recalculation of mcdisp.trs
		       // energies gives different results than file


#include <mcdisp.h>
#include "../../version"
#include "myev.c"

#ifdef _THREADS
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
#define NUM_THREADS 2

// ----------------------------------------------------------------------------------- //
// Declares a struct to store all the information needed for each disp_calc iteration
// ----------------------------------------------------------------------------------- //
typedef struct{
   ComplexMatrix **chi, **chibey;
   ComplexMatrix **Echargedensity;mfcf  **qee_real, **qee_imag;
   ComplexMatrix **Espindensity;mfcf  **qsd_real, **qsd_imag;
   ComplexMatrix **Eorbmomdensity;mfcf  **qod_real, **qod_imag;
   ComplexMatrix **Ep;mfcf  **qep_real, **qep_imag;
   ComplexMatrix **Em;mfcf  **qem_real, **qem_imag;
   ComplexMatrix **ES;mfcf  **qes_real, **qes_imag;
   ComplexMatrix **EL;mfcf  **qel_real, **qel_imag;
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
   int dimA, level, do_verbose;
   double En, intensity, intensitybey, QQ;
   double epsilon; int iE;
   intcalcapr_input(int _dimA, int _tid, int _level, int _doverb, double _En)
   { 
      thread_id = _tid; dimA = _dimA; level = _level; do_verbose = _doverb; En = _En;
   }
};
// ----------------------------------------------------------------------------------- //
// Declares these variables global, so all threads can see them
// ----------------------------------------------------------------------------------- //
intcalcapr_thread_data thrdat;
MUTEX_TYPE mutex_loop;
MUTEX_TYPE mutex_index;
EVENT_TYPE checkfinish;

// ----------------------------------------------------------------------------------- //
// Routine to calculate the matrix product and sum for threaded calculation of jsss
// ----------------------------------------------------------------------------------- //
#endif // if _THREADS

#include "mcdisp_intcalc.c"
#include "mcdisp_observables.c"

#ifdef _THREADS
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
#ifdef _THREADS
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
         s=J.in(i1,j1,k1); 

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
	 ss=J.in(i,j,k);
//          if (do_verbose==1) {printf("#s=%i %i %i  s'=%i %i %i\n",i,j,k,i1,j1,k1);}
          // sum up 

//         mdl - Changed 110710 - To speed up computation by calculating exp(-2i.Pi.Q.d) real and imag parts separately, 
//                                and put into J.mati(s,ss) directly without using intermediate jsss matrix.
           complex<double> **jsss = J.mati(s,ss).M;

//         ComplexMatrix jsss(1,ini.nofcomponents*md.baseindex_max(i1,j1,k1),1,ini.nofcomponents*md.baseindex_max(i,j,k));
//         jsss=0;
        
	  sl=(*inputpars.jjj[ll]).sublattice[l]; // the whole loop has also to be done 
                                                 // for all the other transitions of sublattice sl

           
          // therefore calculate offset of the set of transitions
          for(tl=1;tl<=md.noft(i1,j1,k1,ll);++tl){ jsi = ini.nofcomponents*(md.baseindex(i1,j1,k1,ll,tl)-1);
	  for(tll=1;tll<=md.noft(i,j,k,sl);++tll){ jsj = ini.nofcomponents*(md.baseindex(i,j,k,sl,tll)-1);
	  
         
	     for(m=1;m<=ini.nofcomponents;++m){for(n=1;n<=ini.nofcomponents;++n){ //this should also be ok for nofcomponents > 3 !!! (components 1-3 denote the magnetic moment)
         jjval = (*inputpars.jjj[ll]).jij[l](m,n);  
         jsss[jsi+m][jsj+n] += complex<double>(jjval*REexpqd, jjval*IMexpqd);
                                              }                                 } // but orbitons should be treated correctly by extending 3 to n !!
	                                         }}        

          ++nofneighbours; // count neighbours summed up
	 }}}
   }
#ifdef _THREADS
   myinput->dimA=nofneighbours;
#endif
}
#ifdef _THREADS
#undef inputpars
#undef ini
#undef md
#undef J
#undef q
#endif

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


// some function to efficiently output standard deviations
void staout(FILE*fout,double & sta,double & sta_int,double & sta_without_antipeaks,double & sta_int_without_antipeaks,double & sta_without_weights,double & sta_int_without_weights,double & sta_without_antipeaks_weights,double & sta_int_without_antipeaks_weights)
    {
    fprintf(fout,"#definitions:\n#\n");
    fprintf(fout,"#sta                               = sum_i |weight(i)|*[Eexp(i) - nearestEcalc]^[2*sign(weight(i))]\n");
    fprintf(fout,"#sta_int                           = sum_i |weight(i)|*[Eexp(i) - nearestEcalc_with_Int>%gb/srf.u.]^[2*sign(weight(i))]\n",SMALLINT);
    fprintf(fout,"#sta_without_antipeaks             = sum_i_with_weight(i)>0  weight(i)*[Eexp(i) - nearestEcalc]^2\n");
    fprintf(fout,"#sta_int_without_antipeaks         = sum_i_with_weight(i)>0  weight(i)*[Eexp(i) - nearestEcalc_with_Int>%gb/srf.u.]^2\n",SMALLINT);
    fprintf(fout,"#sta_without_weights               = sum_i [Eexp(i) - nearestEcalc]^[2*sign(weight(i))]\n");
    fprintf(fout,"#sta_int_without_weights           = sum_i [Eexp(i) - nearestEcalc_with_Int>%gb/srf.u.]^[2*sign(weight(i))]\n",SMALLINT);
    fprintf(fout,"#sta_without_antipeaks_weights     = sum_i_with_weight(i)>0 [Eexp(i) - nearestEcalc]^2\n");
    fprintf(fout,"#sta_int_without_antipeaks_weights = sum_i_with_weight(i)>0 [Eexp(i) - nearestEcalc_with_Int>%gb/srf.u.]^2\n#\n",SMALLINT);
    fprintf(fout,"#!sta= %8.6g \n",sta);
    fprintf(fout,"#!sta_int= %8.6g \n",sta_int);
    fprintf(fout,"#!sta_without_antipeaks= %8.6g \n",sta_without_antipeaks);
    fprintf(fout,"#!sta_int_without_antipeaks= %8.6g \n",sta_int_without_antipeaks);
    fprintf(fout,"#!sta_without_weights= %8.6g \n",sta_without_weights);
    fprintf(fout,"#!sta_int_without_weights= %8.6g \n",sta_int_without_weights);
    fprintf(fout,"#!sta_without_antipeaks_weights= %8.6g \n",sta_without_antipeaks_weights);
    fprintf(fout,"#!sta_int_without_antipeaks_weights= %8.6g \n",sta_int_without_antipeaks_weights);
    }



// procedure to calculate the dispersion
void dispcalc(inimcdis & ini,par & inputpars,int calc_fast, int do_gobeyond,int do_Erefine,int do_jqfile,int do_createtrs,int do_readtrs, int do_verbose,int maxlevels,double minE,double maxE,double ninit,double pinit,double epsilon, const char * filemode)
{ int i,j,k,l,ll,s,ss,i1,i2,j1,j2,k1,k2,l1,l2,t1,t2,b,bb,m,n,tn;
  FILE * fin;
  FILE * fout;
  FILE * foutqei;
  FILE * foutqep;
  FILE * foutqee;
  FILE * foutqem;
  FILE * foutqes;
  FILE * foutqel;
  FILE * foutqsd;
  FILE * foutqod;
  FILE * fout1;
  FILE * foutds;
  FILE * foutdstot;
  FILE * foutds1;
  FILE * jqfile;
  float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
  double E;
  double sta=0,sta_int=0,sta_without_antipeaks=0,sta_int_without_antipeaks=0;
  double sta_without_weights=0,sta_int_without_weights=0,sta_without_antipeaks_weights=0,sta_int_without_antipeaks_weights=0;
  double jqsta=-1.0e10;  double jqsta_int=0;
  double jq0=0;
  Vector hkl(1,3),q(1,3),qold(1,3),qijk(1,3);                 
  Vector mf(1,ini.nofcomponents);
  int jmin;
  IntVector noftransitions(1,inputpars.nofatoms); // vector to remember how many transitions are on each atom
  //int offset[inputpars.nofatoms+1]; // vector to remember where higher  transitions are stored
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

  double chargedensity_gamman; 
  Vector chargedensity_gamma(1,CHARGEDENS_EV_DIM);ComplexVector chargedensity_coeff1(1,CHARGEDENS_EV_DIM);
  // extended transition matrix Mij,transformation matrix Uij
  ComplexMatrix chargedensity_Mijkl(1,CHARGEDENS_EV_DIM,1,CHARGEDENS_EV_DIM),chargedensity_Uijkl(1,CHARGEDENS_EV_DIM,1,CHARGEDENS_EV_DIM);

  double spindensity_gamman; 
  Vector spindensity_gamma(1,SPINDENS_EV_DIM);ComplexVector spindensity_coeff1(1,SPINDENS_EV_DIM);
  // extended transition matrix Mij,transformation matrix Uij
  ComplexMatrix spindensity_Mijkl(1,SPINDENS_EV_DIM,1,SPINDENS_EV_DIM),spindensity_Uijkl(1,SPINDENS_EV_DIM,1,SPINDENS_EV_DIM);

  double orbmomdensity_gamman; 
  Vector orbmomdensity_gamma(1,ORBMOMDENS_EV_DIM);ComplexVector orbmomdensity_coeff1(1,ORBMOMDENS_EV_DIM);
  // extended transition matrix Mij,transformation matrix Uij
  ComplexMatrix orbmomdensity_Mijkl(1,ORBMOMDENS_EV_DIM,1,ORBMOMDENS_EV_DIM),orbmomdensity_Uijkl(1,ORBMOMDENS_EV_DIM,1,ORBMOMDENS_EV_DIM);

  double magmom_gamman; 
  Vector magmom_gamma(1,MAGMOM_EV_DIM);ComplexVector magmom_coeff1(1,MAGMOM_EV_DIM);
  // extended transition matrix Mij,transformation matrix Uij
  ComplexMatrix magmom_Mijkl(1,MAGMOM_EV_DIM,1,MAGMOM_EV_DIM),magmom_Uijkl(1,MAGMOM_EV_DIM,1,MAGMOM_EV_DIM);

  double spin_gamman; 
  Vector spin_gamma(1,SPIN_EV_DIM);ComplexVector spin_coeff1(1,SPIN_EV_DIM);
  // extended transition matrix Mij,transformation matrix Uij
  ComplexMatrix spin_Mijkl(1,SPIN_EV_DIM,1,SPIN_EV_DIM),spin_Uijkl(1,SPIN_EV_DIM,1,SPIN_EV_DIM);

  double orbmom_gamman; 
  Vector orbmom_gamma(1,ORBMOM_EV_DIM);ComplexVector orbmom_coeff1(1,ORBMOM_EV_DIM);
  // extended transition matrix Mij,transformation matrix Uij
  ComplexMatrix orbmom_Mijkl(1,ORBMOM_EV_DIM,1,ORBMOM_EV_DIM),orbmom_Uijkl(1,ORBMOM_EV_DIM,1,ORBMOM_EV_DIM);

  double phonon_gamman; 
  Vector phonon_gamma(1,PHONON_EV_DIM);ComplexVector phonon_coeff1(1,PHONON_EV_DIM);
  // extended transition matrix Mij,transformation matrix Uij
  ComplexMatrix phonon_Mijkl(1,PHONON_EV_DIM,1,PHONON_EV_DIM),phonon_Uijkl(1,PHONON_EV_DIM,1,PHONON_EV_DIM);

  //calculate single ion properties of every atom in magnetic unit cell
  mdcf md(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),inputpars.nofatoms,ini.nofcomponents);

 if(calc_fast==1){minE=1e-10;} // in case of fast calculation only consider positive deltas !!
  
 //initialize output of transitions
 if (do_readtrs==0)
 {printf("saving ./results/mcdisp.trs\n");
  fout = fopen_errchk ("./results/mcdisp.trs","w");
   fprintf(fout, "#output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcdisp.trs-->\n");
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(fout,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"#(*)The unpolarized powder average neutron cross section sigma for each transition \n");
   fprintf(fout,"#   is calculated neglecting the formfactor, the Debye Wallerfactor, factor k'/k as follows:\n");
   fprintf(fout,"#-------------------------------------------------------------- \n");
   fprintf(fout,"#            Transition intensities in barn/sr.                |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"# =                    - E /T     -----                        |\n");
   fprintf(fout,"# |                   e   i       \\                   2        |\n");
   fprintf(fout,"# |     = const --------------     >    |<i,r|M |k,s>|         |\n");
   fprintf(fout,"# |             ----    - E /T    /            T               |\n");
   fprintf(fout,"# =             >    n e   i      -----                        |\n");
   fprintf(fout,"#  E -> E       ----  i            r,s                         |\n");
   fprintf(fout,"#   i    k       i                                             |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                            with                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                -----                         |\n");
   fprintf(fout,"#                      2     2   \\                        2    |\n");
   fprintf(fout,"#        |<i,r|M |k,s>|   = ---   >     |<i,r|M -<M >|k,s>|    |\n");
   fprintf(fout,"#               T            3   /             u   u           |\n");
   fprintf(fout,"#                                -----                         |\n");
   fprintf(fout,"#                             u = x,y,z                        |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                             and                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                   1       2                  |\n");
   fprintf(fout,"#                  const  =      ( --- r   )                   |\n");
   fprintf(fout,"#                                   2   0                      |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                      -12                     |\n");
   fprintf(fout,"#                  r     = -0.53908* 10    cm                  |\n");
   fprintf(fout,"#                   0                                          |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                 M   =  L  + 2 S  = g  J                      |\n");
   fprintf(fout,"#                                     J                        |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#--------------------------------------------------------------|\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                       1.Sum rule :                           |\n");
   fprintf(fout,"#                                                 - E /T       |\n");
   fprintf(fout,"#                                            n   e   i         |\n");
   fprintf(fout,"#  ----  =       2                            i                |\n");
   fprintf(fout,"#  >     |     =--- *g *g *const * J(J+1) * ----------------   |\n");
   fprintf(fout,"#  ----  =       3    J  J                  ----     - E /T    |\n");
   fprintf(fout,"#   k     E -> E                            >    n  e   i      |\n");
   fprintf(fout,"#          i    k                           ----  i            |\n");
   fprintf(fout,"#                                            i                 |\n");
   fprintf(fout,"#--------------------------------------------------------------|\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                       2. sum rule :                          |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#                                                              |\n");
   fprintf(fout,"#            ----  =            2                              |\n");
   fprintf(fout,"#            >     |         = --- * const*g *g *J(J+1)        |\n");
   fprintf(fout,"#            ----  =            3           J  J               |\n");
   fprintf(fout,"#             k,i   E -> E                                     |\n");
   fprintf(fout,"#                    i    k                                    |\n");
   fprintf(fout,"#-------------------------------------------------------------- \n");
   fprintf(fout,"#! ninit= %g (max number of initial states)\n",ninit);
   fprintf(fout,"#! pinit= %g (minimum population number of initial states)\n",pinit);
   fprintf(fout,"#! T= %g K Ha=%g Hb=%g Hc=%g T\n",ini.T,ini.Hext(1),ini.Hext(2),ini.Hext(3));
   fprintf(fout,"#*********************************************************************\n");
          fprintf (fout, "#i j k ionnr transnr energy |gamma_s|  sigma [barn/sr](*)\n");
  for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  for(l=1;l<=inputpars.nofatoms;++l){
   fprintf(stdout,"trying du1calc for ion %i in crystallographic unit cell %i %i %i:\n",l,i,j,k);
      for(ll=1;ll<=ini.nofcomponents;++ll)
       {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} //mf ... mean field vector of atom s in first 
                                                              //crystallographic unit of magnetic unit cell
   
   md.est_ini(i,j,k,l,(*inputpars.jjj[l]).eigenstates(mf,ini.Hext,ini.T)); 
   (*inputpars.jjj[l]).transitionnumber=1;
   fprintf(stdout,"transition number %i: ",(*inputpars.jjj[l]).transitionnumber);
   d=1e10;u1(1)=complex <double> (ninit,pinit);i1=(*inputpars.jjj[l]).du1calc(ini.T,mf,ini.Hext,u1,d,md.est(i,j,k,l));
   Mijkl = u1^u1;
      // here Mijkl is a nxn matrix n ... numberofcomponents
   noftransitions(l)=0;
    while (minE>d||d>maxE) //only consider transition if it is in interval minE/maxE
     {//first and following  transitions out of energy range ... do not consider them
     fprintf(stdout," .... transition not stored because out of interval [minE,maxE]=[%g,%g]meV\n",minE,maxE);
     ++(*inputpars.jjj[l]).transitionnumber;
     fprintf(stdout,"transition number %i: ",(*inputpars.jjj[l]).transitionnumber);
     d=1e10;u1(1)=complex <double> (ninit,pinit);(*inputpars.jjj[l]).du1calc(ini.T,mf,ini.Hext,u1,d,md.est(i,j,k,l));
     Mijkl=u1^u1;
     if((*inputpars.jjj[l]).transitionnumber>i1){fprintf(stderr,"ERROR mcdisp.par: no transition found within energy in range [minE,maxE]=[%g,%g] found\n (within first crystallographic unit of magnetic unit cell)\n please increase energy range in option -maxE and -minE\n",minE,maxE);
                            exit(EXIT_FAILURE);}
     }
     // calculate powder neutron intensities 
// old: removed 25.2.2013 and replaced by dm1calc derived intensities
//     if ((*inputpars.jjj[l]).gJ!=0)
//     { for(k1=1;k1<=3;++k1){intensityp+=(*inputpars.jjj[l]).gJ*(*inputpars.jjj[l]).gJ*real(Mijkl(k1,k1));}}
//     else
//     { for(k1=1;k1<=3;++k1){intensityp+=4*real(Mijkl(2*k1-1,2*k1-1))+real(Mijkl(2*k1,2*k1))+2*real(Mijkl(2*k1-1,2*k1))+2*real(Mijkl(2*k1,2*k1-1));}}
     double intensityp=0, intensitym=0;int k1; ComplexVector dm1(1,3);
     d=1e10;dm1(1)=complex <double> (ninit,pinit);
     if((*inputpars.jjj[l]).dm1calc(ini.T,mf,ini.Hext,dm1,d,md.est(i,j,k,l))) // if dm1calc is implemented for this ion
     {intensityp+=Norm2(dm1); // Norm2 ... sum of modulus squared
     intensityp*=0.048434541067;intensitym=intensityp;// prefactor for intensity in barn/sr is 2/3*0.53908*0.53908/4= 0.048434541067
     if (d>SMALL){if(d/ini.T/KB<20){intensitym=-intensityp/(1-exp(d/ini.T/KB));intensityp/=(1-exp(-d/ini.T/KB));}else{intensitym=0;}}
                                  else{intensityp=intensityp*ini.T*KB;intensitym=intensityp;}
     }
     else
     {intensityp=-1;intensitym=-1;}     


                

   jmin=(*inputpars.jjj[l]).transitionnumber;
  if (do_verbose==1){fprintf(stdout,"Matrix M(s=%i %i %i)\n",i,j,k);
                  myPrintComplexMatrix(stdout,Mijkl);fprintf(stdout,"...diagonalising\n");
                    } 
     // diagonalizeMs to get unitary transformation matrix Us and eigenvalues gamma
     //myEigenSystemHermitean (Mijkl,gamma,Uijkl,sort=1,maxiter);
     gamma=0;gamma(ini.nofcomponents)=real(Trace(Mijkl));
     if(minE<d&&d<maxE)
    { fprintf(fout,"%i %i %i  %i     %i     %g  %g  %g\n",i,j,k,l,jmin,myround(d),myround(gamma(ini.nofcomponents)),myround(intensityp));
     ++noftransitions(l);}
    if(d>=0&&minE<-d&&-d<maxE) // do not print negative energy transition if d<0 (d<0 means transiton to the same level)
    { fprintf(fout,"%i %i %i  %i     %i     %g  %g  %g\n",i,j,k,l,jmin,myround(-d),myround(gamma(ini.nofcomponents)),myround(intensitym));
    ++noftransitions(l);}

   for(j1=jmin+1;j1<=i1;++j1) 
   // for every transition add new "atom" to list ... changed to "setnumberoftransitions" for
   // ion l in cryst unit cell ijk
   {
      if (noftransitions(l)>=maxlevels){if (do_verbose) fprintf(stdout,"Maximum number of transitions for ion %i reached\n",l);
                       break;} //check if number of transitions  bigger than maximal number 
                                //(given by -max option in command line)
      
        (*inputpars.jjj[l]).transitionnumber=-j1; // try calculation for transition  j
      fprintf(stdout,"transition number %i: ",(*inputpars.jjj[l]).transitionnumber);
      d=maxE;u1(1)=complex <double> (ninit,pinit);(*inputpars.jjj[l]).du1calc(ini.T,mf,ini.Hext,u1,d,md.est(i,j,k,l));
      Mijkl = u1^u1;
   //printf("noftransitions read by mcdisp: %i",i1);
      
      if (minE<d&&d<maxE) //only consider transition if it is in interval emin/emax
     { 
      // calculate powder neutron intensities
     intensityp=0; intensitym=0;
    // old: removed 25.2.2013 and replaced by dm1calc derived intensities
    //if ((*inputpars.jjj[l]).gJ!=0)
     //{ for(k1=1;k1<=3;++k1){intensityp+=(*inputpars.jjj[l]).gJ*(*inputpars.jjj[l]).gJ*real(Mijkl(k1,k1));}}
     //else
    //{ for(k1=1;k1<=3;++k1){intensityp+=4*real(Mijkl(2*k1-1,2*k1-1))+real(Mijkl(2*k1,2*k1))+2*real(Mijkl(2*k1-1,2*k1))+2*real(Mijkl(2*k1,2*k1-1));}}
    //intensityp*=0.048434541067;intensitym=intensityp;// prefactor for intensity in barn/sr is 2/3*0.53908*0.53908/4= 0.048434541067
    // if (d>SMALL){if(d/ini.T/KB<20){intensitym=-intensityp/(1-exp(d/ini.T/KB));intensityp/=(1-exp(-d/ini.T/KB));}else{intensitym=0;}}
    //                              else{intensityp=intensityp*ini.T*KB;intensitym=intensityp;}
     d=maxE; dm1(1)=complex <double> (ninit,pinit);
     if((*inputpars.jjj[l]).dm1calc(ini.T,mf,ini.Hext,dm1,d,md.est(i,j,k,l))) // if dm1calc is implemented for this ion
     {intensityp+=Norm2(dm1); // Norm2 ... sum of modulus squared
     intensityp*=0.048434541067;intensitym=intensityp;// prefactor for intensity in barn/sr is 2/3*0.53908*0.53908/4= 0.048434541067
     if (d>SMALL){if(d/ini.T/KB<20){intensitym=-intensityp/(1-exp(d/ini.T/KB));intensityp/=(1-exp(-d/ini.T/KB));}else{intensitym=0;}}
                                  else{intensityp=intensityp*ini.T*KB;intensitym=intensityp;}
     }
     else
     {intensityp=-1;intensitym=-1;}     
         
      if (do_verbose==1){fprintf(stdout,"Matrix M(s=%i %i %i)\n",i,j,k);
                  myPrintComplexMatrix(stdout,Mijkl);
                 //fprintf(stdout,"Frobenius Norm(M)=sqrt(sum_i,j |Mij|^2)=%g\n",NormFro(Mijkl));
                  fprintf(stdout,"...diagonalising\n");
                    }  // diagonalizeMs to get unitary transformation matrix Us and eigenvalues gamma
     //if(NormFro(Mijkl)>SMALL*1e-5) {myEigenSystemHermitean (Mijkl,gamma,Uijkl,sort=1,maxiter); } \\ diagonalisation only possible for nonzero matrices !!!
     //             else {gamma=0;gamma(ini.nofcomponents)=Trace(Mijkl);}
      gamma=0;gamma(ini.nofcomponents)=real(Trace(Mijkl));
     if(minE<d&&d<maxE)
     {fprintf(fout,"%i %i %i  %i     %i     %g  %g  %g\n",i,j,k,l,j1,myround(d),myround(gamma(ini.nofcomponents)),myround(intensityp));
      ++noftransitions(l);}
     if(d>=0&&minE<-d&&-d<maxE)// do not print negative energy transition if d<0 (d<0 means transiton to the same level)
     {fprintf(fout,"%i %i %i  %i     %i     %g  %g  %g\n",i,j,k,l,j1,myround(-d),myround(gamma(ini.nofcomponents)),myround(intensitym));
      ++noftransitions(l);}
     }else
     {fprintf(stdout," .... transition not stored because  out of interval [minE,maxE]=[%g,%g]meV\n",minE,maxE);
     }

   (*inputpars.jjj[l]).transitionnumber=jmin; // put back transition number for 1st transition

   }  
  }}}}
   fclose(fout);
 } // do_readtrs==0
  if (do_createtrs==1){fprintf(stdout,"single ion transition file ./results/mcdisp.trs created - please comment transitions which should not enter the calculation and restart with option -t\n");exit(0);}

  printf("\n#reading ./results/mcdisp.trs\n\n");
// read transitions to be considered from file
 for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  fin = fopen_errchk ("./results/mcdisp.trs","rb");
  noftransitions=0;
 int nparread=0;double Tr,Har,Hbr,Hcr;
 char instr[MAXNOFCHARINLINE];
 while(fgets(instr,MAXNOFCHARINLINE,fin)!=NULL&&nparread<6)
 {nparread+=1-extract(instr,"ninit",ninit);
  nparread+=1-extract(instr,"pinit",pinit);
  nparread+=1-extract(instr,"T",Tr);
  nparread+=1-extract(instr,"Ha",Har);
  nparread+=1-extract(instr,"Hb",Hbr);
  nparread+=1-extract(instr,"Hc",Hcr);
 }

 if (Tr!=ini.T||Har!=ini.Hext(1)||Hbr!=ini.Hext(2)||Hcr!=ini.Hext(3)||nparread!=6){fprintf(stderr,"ERROR: reading mcdisp.trs one of the parameters not set or not in line with mcdisp.ini: ninit pinit T Ha Hb Hc ! \n");exit(EXIT_FAILURE);}
  while (feof(fin)==0)
  {if ((i1=inputline(fin,nn))>=5)
   {if(i==(int)nn[1]&&j==(int)nn[2]&&k==(int)nn[3])
    {l=(int)nn[4];++noftransitions(l);}
  }} 
     fclose(fin);
  //just read dimensions of matrices in md
    md.set_noftransitions(i,j,k,noftransitions);
  // for later use:
    for(l=1;l<=inputpars.nofatoms;++l){ //save eigenstates of ions (if possible)
      for(ll=1;ll<=ini.nofcomponents;++ll)
       {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} //mf ... mean field vector of atom s in first 
                                                              //crystallographic unit of magnetic unit cell
       if(do_readtrs!=0)md.est_ini(i,j,k,l,(*inputpars.jjj[l]).eigenstates(mf,ini.Hext,ini.T)); // initialize ests if not already done above
       }

    md.U(i,j,k)=0; // initialize transformation matrix U
    md.M(i,j,k)=0; // initialize matrix M
    md.sqrt_gamma(i,j,k)=0; // and sqrt(gamma^s) matrix sqrt_gamma
 }}}

// determine the dimension of the dynamical matrix Ass' s,s'=1....dimA
int dimA=0;
for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 dimA+=md.baseindex_max(i1,j1,k1);
 }}}

// matrix E^s_alpha' used to store the coefficients for extending the eigenvector (see manual)
ComplexMatrix Echargedensity(1,dimA,1,CHARGEDENS_EV_DIM);Echargedensity=0;
ComplexMatrix Espindensity(1,dimA,1,SPINDENS_EV_DIM);Espindensity=0;
ComplexMatrix Eorbmomdensity(1,dimA,1,ORBMOMDENS_EV_DIM);Eorbmomdensity=0;
ComplexMatrix Ephonon(1,dimA,1,PHONON_EV_DIM);Ephonon=0;
ComplexMatrix Emagmom(1,dimA,1,MAGMOM_EV_DIM);Emagmom=0;
ComplexMatrix Espin(1,dimA,1,SPIN_EV_DIM);Espin=0;
ComplexMatrix Eorbmom(1,dimA,1,ORBMOM_EV_DIM);Eorbmom=0;
  
 for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  for(l=1;l<=inputpars.nofatoms;++l){
  fin = fopen_errchk ("./results/mcdisp.trs","rb");
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

      fprintf(stdout,"#transition %i of ion %i of cryst. unit cell at pos  %i %i %i in mag unit cell:\n",tn,l,i,j,k);
      if(nn[6]<SMALL){fprintf(stdout,"#-");}else{fprintf(stdout,"#+");}
      
        j1=(*inputpars.jjj[l]).transitionnumber; // try calculation for transition  j
        (*inputpars.jjj[l]).transitionnumber=-tn; // try calculation for transition  tn with printout
      d=1e10;u1(1)=complex <double> (ninit,pinit);(*inputpars.jjj[l]).du1calc(ini.T,mf,ini.Hext,u1,d,md.est(i,j,k,l));
        Mijkl = u1^u1;gamman=Norm2(u1);u1/=sqrt(gamman);
       if(fabs((fabs(d)-fabs(nn[6]))/(fabs(nn[6])+1.0))>SMALLEDIF)
        {fprintf(stderr,"ERROR mcdisp: reading mcdisp.trs with transition energy delta %g meV different from internal calculation %g meV\n",nn[6],d);	 
         exit(EXIT_FAILURE);}
           if (do_verbose==1){
                  fprintf(stdout,"#Matrix M(s=%i %i %i %i)\n",i,j,k,l);
                   myPrintComplexMatrix(stdout,Mijkl);
                              }
      
     // diagonalizeMs to get unitary transformation matrix Us
      myEigenSystemHermitean (Mijkl,gamma,Uijkl,sort=1,maxiter);
       if (fabs(gamman-gamma(ini.nofcomponents))>SMALL){fprintf(stderr,"ERROR eigenvalue of single ion matrix M inconsistent: analytic value gamma= %g numerical diagonalisation of M gives gamma= %g\n",gamman,gamma(ini.nofcomponents));
                           exit(EXIT_FAILURE);}
// take highest eigenvector to be the same phase as u1
for(int ii=Uijkl.Rlo(); ii<=Uijkl.Rhi(); ii++){if (fabs(abs(u1(ii))-abs(Uijkl(ii,ini.nofcomponents)))>SMALL)
                                                {fprintf(stderr,"ERROR eigenvector of single ion matrix M inconsistent\n");
                                                 myPrintComplexVector(stderr,u1);u1=Uijkl.Column(ini.nofcomponents);myPrintComplexVector(stderr,u1);exit(EXIT_FAILURE);}
                                               Uijkl(ii,ini.nofcomponents)=u1(ii);}
         // treat correctly case for neutron energy loss
	 if (nn[6]<0) // if transition energy is less than zero do a conjugation of the matrix
	 {  for(int ii=Uijkl.Rlo(); ii<=Uijkl.Rhi(); ii++)       for(int jj=Uijkl.Clo(); jj<=Uijkl.Chi(); jj++)       Uijkl[ii][jj]=conj(Uijkl[ii][jj]);
	 }
       j1=md.baseindex(i,j,k,l,jmin); 
       md.delta(i,j,k)(j1)=nn[6]; // set delta


//----------------------------------OBSERVABLES -------------------------------------------------
if (do_verbose==1){ fprintf(stdout,"# ... recalculate now M(s=%i %i %i %i) with eigenvector dimension for observable=%i\n",i,j,k,l,CHARGEDENS_EV_DIM);}
//-----------------------------------------------------------------------------------
if(ini.calculate_chargedensity_oscillation){(*inputpars.jjj[l]).transitionnumber=-tn;d=1e10;chargedensity_coeff1(1)=complex <double> (ninit,pinit);
   if((*inputpars.jjj[l]).dchargedensity_coeff1(ini.T,mf,ini.Hext,chargedensity_coeff1,d,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,CHARGEDENS_EV_DIM,chargedensity_coeff1,inputpars,chargedensity_Mijkl,md,
             chargedensity_gamma,chargedensity_gamman,chargedensity_Uijkl,maxiter,nn,ini, gamma,Echargedensity);}
if(ini.calculate_spindensity_oscillation){(*inputpars.jjj[l]).transitionnumber=-tn;d=1e10;spindensity_coeff1(1)=complex <double> (ninit,pinit);
   if((*inputpars.jjj[l]).dspindensity_coeff1(ini.T,mf,ini.Hext,spindensity_coeff1,d,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,SPINDENS_EV_DIM,spindensity_coeff1,inputpars,spindensity_Mijkl,md,
             spindensity_gamma,spindensity_gamman,spindensity_Uijkl,maxiter,nn,ini, gamma,Espindensity);}
if(ini.calculate_orbmomdensity_oscillation){(*inputpars.jjj[l]).transitionnumber=-tn;d=1e10;orbmomdensity_coeff1(1)=complex <double> (ninit,pinit);
   if((*inputpars.jjj[l]).dorbmomdensity_coeff1(ini.T,mf,ini.Hext,orbmomdensity_coeff1,d,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,ORBMOMDENS_EV_DIM,orbmomdensity_coeff1,inputpars,orbmomdensity_Mijkl,md,
             orbmomdensity_gamma,orbmomdensity_gamman,orbmomdensity_Uijkl,maxiter,nn,ini, gamma,Eorbmomdensity);}
if(ini.calculate_phonon_oscillation){(*inputpars.jjj[l]).transitionnumber=-tn;d=1e10;phonon_coeff1(1)=complex <double> (ninit,pinit);
   if((*inputpars.jjj[l]).dp1calc(ini.T,mf,ini.Hext,phonon_coeff1,d,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,PHONON_EV_DIM,phonon_coeff1,inputpars,phonon_Mijkl,md,
             phonon_gamma,phonon_gamman,phonon_Uijkl,maxiter,nn,ini, gamma,Ephonon);}
if(ini.calculate_magmoment_oscillation){(*inputpars.jjj[l]).transitionnumber=-tn;d=1e10;magmom_coeff1(1)=complex <double> (ninit,pinit);
   if((*inputpars.jjj[l]).dm1calc(ini.T,mf,ini.Hext,magmom_coeff1,d,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,MAGMOM_EV_DIM,magmom_coeff1,inputpars,magmom_Mijkl,md,
             magmom_gamma,magmom_gamman,magmom_Uijkl,maxiter,nn,ini, gamma,Emagmom);}
if(ini.calculate_spinmoment_oscillation){(*inputpars.jjj[l]).transitionnumber=-tn;d=1e10;spin_coeff1(1)=complex <double> (ninit,pinit);
   if((*inputpars.jjj[l]).dS1calc(ini.T,mf,ini.Hext,spin_coeff1,d,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,SPIN_EV_DIM,spin_coeff1,inputpars,spin_Mijkl,md,
             spin_gamma,spin_gamman,spin_Uijkl,maxiter,nn,ini, gamma,Espin);}
if(ini.calculate_orbmoment_oscillation){(*inputpars.jjj[l]).transitionnumber=-tn;d=1e10;orbmom_coeff1(1)=complex <double> (ninit,pinit);
   if((*inputpars.jjj[l]).dL1calc(ini.T,mf,ini.Hext,orbmom_coeff1,d,md.est(i,j,k,l))!=0)
       fillE(jmin,i,j,k,l,ORBMOM_EV_DIM,orbmom_coeff1,inputpars,orbmom_Mijkl,md,
             orbmom_gamma,orbmom_gamman,orbmom_Uijkl,maxiter,nn,ini, gamma,Eorbmom);}
//----------------------------------------------------------------------------------------------

         if (gamma(ini.nofcomponents)>=0&&fabs(gamma(ini.nofcomponents-1))<SMALL) 
                           // mind in manual the 1st dimension alpha=1 corresponds
			   // to the nth dimension here, because myEigensystmHermitean
			   // sorts the eigenvalues according to ascending order !!!
                           {if (nn[6]>SMALL)
			    {md.sqrt_gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=sqrt(gamma(ini.nofcomponents));// gamma(ini.nofcomponents)=sqr(gamma^s)
                            }
			    else if (nn[6]<-SMALL)
                            {md.sqrt_gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=imaginary*sqrt(gamma(ini.nofcomponents));// gamma(ini.nofcomponents)=sqr(gamma^s)
                            }
 			    else
			    { //quasielastic line needs gamma=SMALL .... because Mijkl and therefore gamma have been set to 
			      // wn/kT instead of wn-wn'=SMALL*wn/kT (in jjjpar.cpp -mdcalc routines)
			      //set fix delta but keep sign
			          if (nn[6]>0){md.delta(i,j,k)(j1)=SMALL;
  			     md.sqrt_gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=sqrt(SMALL*gamma(ini.nofcomponents));
                                              }
				  else        {md.delta(i,j,k)(j1)=-SMALL;
                             md.sqrt_gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=imaginary*sqrt(SMALL*gamma(ini.nofcomponents));
			                      }
			    }
			   }else 
                           {fprintf(stderr,"ERROR eigenvalue of single ion matrix <0: ev1=%g ev2=%g ev3=%g ... evn=%g\n",gamma(1),gamma(2),gamma(3),gamma(ini.nofcomponents));
                            exit(EXIT_FAILURE);}
        (* inputpars.jjj[l]).transitionnumber=j1; // put back transition number for 1st transition
        for(m=1;m<=ini.nofcomponents;++m){for(n=1;n<=ini.nofcomponents;++n){
        md.U(i,j,k)(ini.nofcomponents*(j1-1)+m,ini.nofcomponents*(j1-1)+n)=Uijkl(m,n);
        md.M(i,j,k)(ini.nofcomponents*(j1-1)+m,ini.nofcomponents*(j1-1)+n)=Mijkl(m,n);
        }}    
if (do_verbose==1){
                  fprintf(stdout,"#Matrix M(s=%i %i %i)\n",i,j,k);
                  myPrintComplexMatrix(stdout,Mijkl); 
                  fprintf(stdout,"#Eigenvalues:\n");
                  myPrintVector(stdout,gamma); 
                  fprintf(stdout,"#Matrix U(s=%i%i%i)\n",i,j,k);
                  myPrintComplexMatrix(stdout,Uijkl); 
                 }

    }}}
    fclose(fin);

  }}}}

//initialize output files
  errno = 0;
if (do_jqfile==0)
{ printf("#saving mcdisp.qom and mcdisp.qei\n");
  fout = fopen_errchk ("./results/mcdisp.qom",filemode);
  foutqei = fopen_errchk ("./results/mcdisp.qei",filemode);
  foutdstot = fopen_errchk ("./results/mcdisp.dsigma.tot",filemode);
  printf("#saving mcdisp.dsigma.tot\n");

   if (do_Erefine==1){
          errno = 0;
  foutds = fopen_errchk ("./results/mcdisp.dsigma",filemode);
  printf("#saving mcdisp.dsigma\n");
                      }
}
 
// initialize file with jq matrix
if (do_jqfile==1)
{  printf("#saving mcdisp.jq\n");
 jqfile = fopen_errchk ("./results/mcdisp.jq",filemode);
 writeheader(inputpars,jqfile); printf("#saving mcdisp.jq\n");
   fprintf(jqfile,"#!<--mcphas.mcdisp.dsigma.jq-->\n");
   fprintf (jqfile, "#Fourier Transform of 2 Ion Interaction - sta is calculated by comparing the larges eigenvalue\n# to that of the first q vector of the calculation");
   fputs (asctime(loctime),jqfile);
  if (do_verbose==1){   fprintf (jqfile, "#q=(hkl)\n #spin s() - spin s'()\n #3x3 matrix jss'(q) real im .... [meV]\n");}
  else {fprintf(jqfile,"#h  vs  k  vs  l  vs largest eigenvalue of J(hkl) matrix vs components of corresponding eigenvector re im re im re im re im\n");}
}

//MAIN LOOP - do calculation of excitation energy for every Q vector     
int counter;qijk=0;double qincr=0;
for (hkl(1)=ini.qmin(1);hkl(1)<=ini.qmax(1);hkl(1)+=ini.deltaq(1)){
for (hkl(2)=ini.qmin(2);hkl(2)<=ini.qmax(2);hkl(2)+=ini.deltaq(2)){
for (hkl(3)=ini.qmin(3);hkl(3)<=ini.qmax(3);hkl(3)+=ini.deltaq(3)){
 // transform hkl to primitive lattice
 if (ini.hkllist==1){counter=(int)hkl(1);
		     hkl(1)=ini.hkls[counter][1];
		     hkl(2)=ini.hkls[counter][2];
		     hkl(3)=ini.hkls[counter][3];
 		    }
 q=inputpars.r.Transpose()*hkl;

fprintf(stdout,"#q=(%g,%g,%g)\n",hkl(1),hkl(2),hkl(3));
 if(do_verbose==1){fprintf(stdout,"#Setting up J(q) matrix .... \n");}
 // calculate J(q)
 jq J(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),md);
 jq Jl(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),md);

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

#ifdef _THREADS
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
   #endif
   intcalcapr_input *tin[NUM_THREADS];      thrdat.q = q; thrdat.thread_id = -1;
   thrdat.ini = new inimcdis*[NUM_THREADS]; thrdat.inputpars = new par*[NUM_THREADS];
   thrdat.md  = new mdcf*[NUM_THREADS];     thrdat.J         = new jq*[NUM_THREADS];
   for (int ithread=0; ithread<NUM_THREADS; ithread++) 
   {
      tin[ithread] = new intcalcapr_input(dimA,ithread,1,do_verbose,0.); thrdat.J[ithread] = new jq(J); tin[ithread]->dimA=0; 
      thrdat.md[ithread] = new mdcf(md); thrdat.ini[ithread] = new inimcdis(ini); thrdat.inputpars[ithread] = new par(inputpars);
   }
   int thrcount=0, ithread=0, num_threads_started=-1;
#endif
 // calculate Js,ss(Q) summing up contributions from the l=1-paranz parameters
   
   for(ll=1;ll<=inputpars.nofatoms;++ll)
   { //sum up l.th neighbour interaction of crystallographic atom ll
     // 1. transform dn(l) to primitive lattice and round it to integer value
  #ifndef _THREADS
      jsss_mult(ll,nofneighbours,q,inputpars,ini,J,md);
  #else
      thrcount++;
      tin[ithread]->level=ll;
      #if defined  (__linux__) || defined (__APPLE__)
      rc = pthread_create(&threads[ithread], &attr, jsss_mult, (void *) tin[ithread]);
      if(rc) { printf("Error return code %i from thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
      #else
      threads[ithread] = CreateThread(NULL, 0, jsss_mult, (void *) tin[ithread], 0, &tid[ithread]);
      if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %i from thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
      #endif
      num_threads_started = ithread+1;
      if(thrcount%NUM_THREADS==0)
      {
         #if defined  (__linux__) || defined (__APPLE__)
         for(int th=0; th<NUM_THREADS; th++)
            rc = pthread_join(threads[th], &status);
         #else
         WaitForMultipleObjects(NUM_THREADS,threads,TRUE,INFINITE);
         #endif
         ithread=0;
      }
      else ithread++;
#endif
   }
#ifdef _THREADS
    #if defined  (__linux__) || defined (__APPLE__)
    for(int th=0; th<ithread; th++)
       rc = pthread_join(threads[th], &status);
    #else
    WaitForMultipleObjects(ithread,threads,TRUE,INFINITE);
    #endif
    for(int th=0; th<NUM_THREADS; th++) 
    {
       nofneighbours += tin[th]->dimA;
       for(int i1=1;i1<=ini.mf.na();++i1) for(int j1=1;j1<=ini.mf.nb();++j1) for(int k1=1;k1<=ini.mf.nc();++k1)
          for(int i2=1;i2<=ini.mf.na();++i2) for(int j2=1;j2<=ini.mf.nb();++j2) for(int k2=1;k2<=ini.mf.nc();++k2)
             J.mat(i1,j1,k1,i2,j2,k2)+=(*thrdat.J[th]).mat(i1,j1,k1,i2,j2,k2); 
    }
    for (ithread=0; ithread<NUM_THREADS; ithread++) {
       delete thrdat.J[ithread]; delete tin[ithread]; }
    delete[] thrdat.J; //delete tin;
#endif


if (do_jqfile==1){if (do_verbose==1){fprintf (jqfile, "#q=(%g, %g, %g) ",hkl(1),hkl(2),hkl(3));
                                     fprintf(jqfile,"nofneighbours= %li\n",nofneighbours);
                                    }
                  else
                  {fprintf (jqfile, "%g  %g  %g ",hkl(1),hkl(2),hkl(3));}
                  }


if(do_verbose==1){fprintf(stdout,"#Transform J(q) matrix  with U...\n");}

// transform J(s,ss) (with md.U) and multiply 
// (with eigenvalues sqrt(gamma) [here md.sqrt_gamma]
// of matrix Malphabeta ) matrix J to L [here Jl.mati]  ... compare manual
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 s=ini.mf.in(i1,j1,k1);
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  ss=ini.mf.in(i2,j2,k2);
  Jl.mati(s,ss)=md.sqrt_gamma(i1,j1,k1)*md.U(i1,j1,k1).Hermitean()*J.mati(s,ss)*md.U(i2,j2,k2)*md.sqrt_gamma(i2,j2,k2).Conjugate();
//if (do_verbose==1){
//                  fprintf(stdout,"#J(s=%i%i%i,s''=%i%i%i)=\n",i1,j1,k1,i2,j2,k2);
//                  myPrintComplexMatrix(stdout,J.mati(s,ss)); 
//                  fprintf(stdout,"#sqr(gamma_s=%i%i%i)=\n",i1,j1,k1);
//                  myPrintComplexMatrix(stdout,md.sqrt_gamma(i1,j1,k1));
//                  fprintf(stdout,"#U(s=%i%i%i)=\n",i1,j1,k1);
//                  myPrintComplexMatrix(stdout,md.U(i1,j1,k1));
//                  fprintf(stdout,"#sqr(gamma_s=%i%i%i)=\n",i2,j2,k2);
//                  myPrintComplexMatrix(stdout,md.sqrt_gamma(i2,j2,k2));
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
//   s=(((i1-1)*ini.mf.nb()+(j1-1))*ini.mf.nc()+(k1-1))*inputpars.nofatoms;
//    for(l1=1;l1<=inputpars.nofatoms;++l1){
//      Ac(s+l1,s+l1)=md.delta(i1,j1,k1)(l1);}

    for(l1=1;l1<=inputpars.nofatoms;++l1){
     for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
      s=index_s(i1,j1,k1,l1,t1,md,ini);
      b=md.baseindex(i1,j1,k1,l1,t1);
      if(calc_fast==1){ // this is fast algorithms dynamical matrix assignment
                       Lambda(s,s)=+0.5/md.delta(i1,j1,k1)(b);
                       Ac(s,s)=+0.5*md.delta(i1,j1,k1)(b);
                      }
      else            { // this is standard DMD as described in the review rotter et al JPcondMat 2012
                       if(md.delta(i1,j1,k1)(b)<0){Lambda(s,s)=-1;}else{Lambda(s,s)=+1;}
                       Ac(s,s)=md.delta(i1,j1,k1)(b)*Lambda(s,s);
                      }
//      if(do_verbose==1){fprintf(stdout,"#i=%i j=%i k=%i atomnr=%i trans=%i ... s=%i ",i1,j1,k1,l1,t1,s);
//                        fprintf(stdout,"#lambda(%i,%i)xdelta(%i)=%g + i %g\n",s,s,s,real(Ac(s,s)),imag(Ac(s,s)));
//                       }
      }}

   for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
//   ss=(((i2-1)*ini.mf.nb()+(j2-1))*ini.mf.nc()+(k2-1))*inputpars.nofatoms;

    for(l1=1;l1<=inputpars.nofatoms;++l1){ 
     for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
    for(l2=1;l2<=inputpars.nofatoms;++l2){ 
     for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2){
      s=index_s(i1,j1,k1,l1,t1,md,ini);
      ss=index_s(i2,j2,k2,l2,t2,md,ini);
      b=md.baseindex(i1,j1,k1,l1,t1);
      bb=md.baseindex(i2,j2,k2,l2,t2);
     Ac(s,ss)-=Jl.mati(Jl.in(i1,j1,k1),Jl.in(i2,j2,k2))(ini.nofcomponents*(b-1)+ini.nofcomponents,ini.nofcomponents*(bb-1)+ini.nofcomponents);
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
if (do_jqfile==1){
       if (do_verbose==1)
       {//fprintf (jqfile, "#spin (%i*r1 %i*r2 %i*r3) - spin (%i*r1 %i*r2 %i*r3)\n",i1,j1,k1,i2,j2,k2);
         myPrintComplexMatrix(jqfile,J_Q); 
       }

	// diagonalize JQ to get eigenvalues (biggest corresponds to Tn) !!!
         Vector Tn(1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms);
         ComplexMatrix eigenvectors(1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms,1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms);
         myEigenSystemHermitean (J_Q,Tn,eigenvectors,sort=1,maxiter);
         i2=ini.nofcomponents*ini.mf.n()*inputpars.nofatoms;
       if(do_verbose==1)
       {fprintf(jqfile,"#eigenvalues(highest corresponds to Tn, predicted magstructure)\n");
         myPrintVector(jqfile,Tn); 
        fprintf(jqfile,"#eigenvectors(moment direction):\n");
         myPrintComplexMatrix(jqfile,eigenvectors); 
       }
       else
       {
        fprintf(jqfile," %g ",Tn(i2));
        for (i1=1;i1<=i2;++i1)
        {fprintf(jqfile," %6.3g ",real(eigenvectors(i1,i2)));
         fprintf(jqfile," %6.3g ",imag(eigenvectors(i1,i2)));
        }
        fprintf(jqfile,"\n");
       }
       if (jqsta<-0.9e10){jq0=Tn(i2);jqsta=-1e9;}
       else           {if(Tn(i2)>jq0)
                         {if(jqsta<0){jqsta=0;}
                          jqsta+=(Tn(i2)-jq0)*(Tn(i2)-jq0);}
                       else{if((jqsta<0)&(Tn(i2)-jq0>jqsta)){jqsta=Tn(i2)-jq0;}}
                      }
       if (ini.hkllist==1)
	             {double test;
                      for (j1=1;j1<=ini.hkls[counter][0]-3;++j1)
	              {test=fabs(Tn(i2+1-j1)-ini.hkls[counter][j1+3]);
                      jqsta_int+=test*test;}
	             }
 }
 else
 {// no jqfile but excitations to be calculated
 if(do_verbose==1){fprintf(stdout,"#diagonalizing %ix%i matrix A...\n",dimA,dimA);
              //             myPrintComplexMatrix(stdout,Ac); 
                   }
   // diagonalize Ac to get energies  and eigenvectors !!!
   Vector En(1,dimA);
   Vector ints(1,dimA);
   Vector intsbey(1,dimA);
//   myEigenValuesHermitean (Ac,En,sort,maxiter);
  
//   myEigenSystemHermitean (Ac,En,Tau,sort,maxiter);
//    myPrintVector(stdout,En);
   int eigrval = myEigenSystemHermiteanGeneral (Lambda,Ac,En,Tau,sort=0,maxiter);
   if(eigrval>0) {  // 0==sucess. +1==Lambda not hermitian. +2==Ac not hermitian, +3==Lambda AND Ac not hermitian
      if(eigrval>1) {
         fprintf(stderr,"# Dynamical Matrix Ac not hermitian. Check exchange parameter file mcphas.j is consistent\n");
         fprintf(stderr,"#   Note that each interaction between pairs of ions must match - e.g. Interaction between\n");
         fprintf(stderr,"#   Atom 1 with Neighbour 2 (which is atom 2) must equal interaction between Atom 2 with Neighbour 1\n");
         fprintf(stderr,"#   Press q to quit, or any other key to ignore this error.\n"); if(getchar()=='q') exit(1); 
      }
      if(eigrval%2==1) fprintf(stderr,"# Warning: Trace matrix Lambda is not diagonal\n");
   }
   En=1.0/En;if(calc_fast==1){for(i1=1;i1<=dimA;++i1)En(i1)=sqrt(En(i1));}
   sortE(En,Tau);
        //   Tau=Tau.Conjugate();
  	// conjugate inserted 31.10.05, because when calculating simple AF - I noticed
	// that the eigensystemhgermitean returns eigenvectors as column vectors, but
	// the components need to be complex conjugated
         // conjugate removed again MR 11.4.2011 because now done correctly in myev.c
 // check normalisation of eigenvectors --------------------
   ComplexMatrix test(1,dimA,1,dimA);
   test=Tau.Conjugate().Transpose()*Ac*Tau;
   ComplexMatrix unit(1,dimA,1,dimA);unit=1;
   if( NormFro(unit-test)>SMALL){
 if(do_verbose==1){myPrintComplexMatrix(stdout,test); }
 fprintf(stderr,"Error: eigenvectors t not correctly normalised\n"); 
 fprintf(stderr,"   Press q to quit, or any other key to ignore this error.\n"); if(getchar()=='q') exit(1); }
//-------------------------------------------------------
 if(do_verbose==1){// fprintf(stdout,"#eigenvectors (matrix Tau):\n");
                   // myPrintComplexMatrix(stdout,Tau); 
                    fprintf(stdout,"#saving the following eigenvalues (meV) to mcdisp.qom:\n");}
   for (i=1;i<=dimA;++i){
               if(do_verbose==1){fprintf(stdout, " %4.4g",En(i));}
                         }


   // calculate and printout intensities [the energies have already
   // been printed out above, so any refinement of energies during intcalc
   // is not included in the output file]
#ifndef _THREADS  
  double QQ; 
             mfcf qee_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,CHARGEDENS_EV_DIM);
             mfcf qee_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,CHARGEDENS_EV_DIM);
             mfcf qsd_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPINDENS_EV_DIM);
             mfcf qsd_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPINDENS_EV_DIM);
             mfcf qod_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOMDENS_EV_DIM);
             mfcf qod_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOMDENS_EV_DIM);
             mfcf qep_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,PHONON_EV_DIM);
             mfcf qep_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,PHONON_EV_DIM);
             mfcf qem_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,MAGMOM_EV_DIM);
             mfcf qem_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,MAGMOM_EV_DIM);
             mfcf qes_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPIN_EV_DIM);
             mfcf qes_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPIN_EV_DIM);
             mfcf qel_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOM_EV_DIM);
             mfcf qel_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOM_EV_DIM);
#endif
  double diffint=0,diffintbey=0;
  if(do_verbose==1){fprintf(stdout,"\n#calculating  intensities approximately ...\n");}
                  fprintf (fout, " > ");
  diffint=0;diffintbey=0;
                  intcalc_ini(ini,inputpars,md,do_verbose,do_gobeyond,hkl,ninit,pinit);
   if(qincr==0) { //write header for output files
                   fprintf(fout,"#!<--mcphas.mcdisp.qom-->\n");
                   writeheader(inputpars,fout);
                   fprintf (fout, "#dispersion \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energies[meV] > intensities (dipapprox vs full calc) [barn/sr/f.u.]   f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
          
                   fprintf(foutqei,"#!<--mcphas.mcdisp.qei-->\n");
                   writeheader(inputpars,foutqei);
                   fprintf (foutqei, "#dispersion displayytext=E(meV)\n#displaylines=false \n#Ha[T] Hb[T] Hc[T] T[K] h k l Q[A^-1] energy[meV] Imag_dip [barn/sr/f.u.] Imag [barn/sr/f.u.]  f.u.=crystallogrpaphic unit cell (r1xr2xr3)  vs qincrement[1/A] (for plotting)}\n");

  //------------observables-----------------------------------
if(ini.calculate_chargedensity_oscillation)foutqee=evfileinit(filemode,"./results/mcdisp.qee",inputpars,"qee",CHARGEDENS_EV_DIM);
if(ini.calculate_spindensity_oscillation)foutqsd=evfileinit(filemode,"./results/mcdisp.qsd",inputpars,"qsd",SPINDENS_EV_DIM);
if(ini.calculate_orbmomdensity_oscillation)foutqod=evfileinit(filemode,"./results/mcdisp.qod",inputpars,"qod",ORBMOMDENS_EV_DIM);
if(ini.calculate_phonon_oscillation)foutqep=evfileinit(filemode,"./results/mcdisp.qep",inputpars,"qep",PHONON_EV_DIM);
if(ini.calculate_magmoment_oscillation)foutqem=evfileinit(filemode,"./results/mcdisp.qem",inputpars,"qem",MAGMOM_EV_DIM);
if(ini.calculate_spinmoment_oscillation)foutqes=evfileinit(filemode,"./results/mcdisp.qes",inputpars,"qes",SPIN_EV_DIM);
if(ini.calculate_orbmoment_oscillation)foutqel=evfileinit(filemode,"./results/mcdisp.qel",inputpars,"qel",ORBMOM_EV_DIM);
  //-----------------------------------------------------------
                        fprintf(foutdstot,"#!<--mcphas.mcdisp.dsigma.tot-->\n");
                        writeheader(inputpars,foutdstot);
                        fprintf (foutdstot, "#!Total Scattering Cross Section in energy range [emin=%g ; emax=%g]\n#Ha[T] Hb[T] Hc[T] T[K] h k l  dsigma_mag_dip/dOmeg dsigma_mag/dOmeg[barn/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3) vs qincrement[1/A] (for plotting)}",ini.emin,ini.emax);

                        if (do_Erefine==1){
                         fprintf(foutds,"#!<--mcphas.mcdisp.dsigma-->\n");
                         writeheader(inputpars,foutds);
                         fprintf (foutds, "#Scattering Cross Section \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energy[meV] dsigma/dOmegadE' [barn/mev/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)\n");
                         fprintf (foutdstot, "for DMD algorithm  vs summing dsigma for diff energies");
                                          }  
                         fprintf (foutdstot, "\n");                       
                }

   int dim=3;
   if (ini.hkllist==1){dim=(int)((ini.hkls[counter][0]-3)/4);
                       if(ini.hklfile_start_index[0]>0)for(int is=1;is<=ini.hklfile_start_index[0];++is)if(ini.hklfile_start_index[is]==counter)
                       {fprintf(foutqei,"#!hklfile_number=%i\n",is);
                        fprintf(fout,"#!hklfile_number=%i\n",is);
                        if(ini.calculate_chargedensity_oscillation)fprintf(foutqee,"#!hklfile_number=%i\n",is);
                        if(ini.calculate_spindensity_oscillation)fprintf(foutqsd,"#!hklfile_number=%i\n",is);
                        if(ini.calculate_orbmomdensity_oscillation)fprintf(foutqod,"#!hklfile_number=%i\n",is);
                        if(ini.calculate_phonon_oscillation)fprintf(foutqep,"#!hklfile_number=%i\n",is);
                        if(ini.calculate_magmoment_oscillation)fprintf(foutqem,"#!hklfile_number=%i\n",is);
                        if(ini.calculate_spinmoment_oscillation)fprintf(foutqes,"#!hklfile_number=%i\n",is);
                        if(ini.calculate_orbmoment_oscillation)fprintf(foutqel,"#!hklfile_number=%i\n",is);
                        fprintf(foutdstot,"#!hklfile_number=%i\n",is);
                        if (do_Erefine==1){fprintf(foutds,"#!hklfile_number=%i\n",is);}
                       }
                      }
   fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
   for (i=1;i<=dimA;++i){
	       fprintf (fout, " %4.4g ",myround(En(i)));
               if(do_verbose==1){fprintf(stdout, " %4.4g",En(i));}
                         }


                  Vector dd(1,dim),dd_int(1,dim);  dd+=100000.0;dd_int+=100000.0;
                  Vector dd1(1,dim),dd1_int(1,dim);  dd1+=100000.0;dd1_int+=100000.0;
                  Vector dd_without_antipeaks(1,dim),dd_int_without_antipeaks(1,dim);  dd_without_antipeaks+=100000.0;dd_int_without_antipeaks+=100000.0;
                  Vector dd_without_weights(1,dim),dd_int_without_weights(1,dim);  dd_without_weights+=100000.0;dd_int_without_weights+=100000.0;
                  Vector dd_without_antipeaks_weights(1,dim),dd_int_without_antipeaks_weights(1,dim);  dd_without_antipeaks_weights+=100000.0;dd_int_without_antipeaks_weights+=100000.0;

 #ifndef _THREADS
                     ComplexMatrix chi(1,3*dimA,1,3*dimA);
                     ComplexMatrix chibey(1,3*dimA,1,3*dimA);
                     Matrix pol(1,3,1,3);
                     
#else
                  // Populates the thread data structure
                  thrdat.qee_real = new mfcf*[NUM_THREADS];          thrdat.qee_imag = new mfcf*[NUM_THREADS];
                  thrdat.qsd_real = new mfcf*[NUM_THREADS];          thrdat.qsd_imag = new mfcf*[NUM_THREADS];
                  thrdat.qod_real = new mfcf*[NUM_THREADS];          thrdat.qod_imag = new mfcf*[NUM_THREADS];
                  thrdat.qep_real = new mfcf*[NUM_THREADS];          thrdat.qep_imag = new mfcf*[NUM_THREADS];
                  thrdat.qem_real = new mfcf*[NUM_THREADS];          thrdat.qem_imag = new mfcf*[NUM_THREADS];
                  thrdat.qes_real = new mfcf*[NUM_THREADS];          thrdat.qes_imag = new mfcf*[NUM_THREADS];
                  thrdat.qel_real = new mfcf*[NUM_THREADS];          thrdat.qel_imag = new mfcf*[NUM_THREADS];
                  thrdat.chi      = new ComplexMatrix*[NUM_THREADS]; thrdat.chibey   = new ComplexMatrix*[NUM_THREADS];
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
                     tin[ithread] = new intcalcapr_input(dimA,ithread,1,do_verbose,En);
                     thrdat.chi[ithread] = new ComplexMatrix(1,3*dimA,1,3*dimA);
                     thrdat.chibey[ithread] = new ComplexMatrix(1,3*dimA,1,3*dimA);
                     thrdat.pol[ithread] = new Matrix(1,3,1,3);
                     

                     thrdat.qee_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,CHARGEDENS_EV_DIM);
                     thrdat.qee_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,CHARGEDENS_EV_DIM);
                     thrdat.Echargedensity[ithread] = new ComplexMatrix(1,dimA,1,CHARGEDENS_EV_DIM); *thrdat.Echargedensity[ithread]=Echargedensity;
                     thrdat.qsd_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPINDENS_EV_DIM);
                     thrdat.qsd_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPINDENS_EV_DIM);
                     thrdat.Espindensity[ithread] = new ComplexMatrix(1,dimA,1,CHARGEDENS_EV_DIM); *thrdat.Espindensity[ithread]=Espindensity;
                     thrdat.qod_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOMDENS_EV_DIM);
                     thrdat.qod_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOMDENS_EV_DIM);
                     thrdat.Eorbmomdensity[ithread] = new ComplexMatrix(1,dimA,1,CHARGEDENS_EV_DIM); *thrdat.Eorbmomdensity[ithread]=Eorbmomdensity;
                     thrdat.qep_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,PHONON_EV_DIM);
                     thrdat.qep_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,PHONON_EV_DIM);
                     thrdat.Ephonon[ithread] = new ComplexMatrix(1,dimA,1,CHARGEDENS_EV_DIM); *thrdat.Ephonon[ithread]=Ephonon;
                     thrdat.qem_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,MAGMOM_EV_DIM);
                     thrdat.qem_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,MAGMOM_EV_DIM);
                     thrdat.Emagmom[ithread] = new ComplexMatrix(1,dimA,1,CHARGEDENS_EV_DIM); *thrdat.Emagmom[ithread]=Emagmom;
                     thrdat.qes_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPIN_EV_DIM);
                     thrdat.qes_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,SPIN_EV_DIM);
                     thrdat.Espin[ithread] = new ComplexMatrix(1,dimA,1,CHARGEDENS_EV_DIM); *thrdat.Espin[ithread]=Espin;
                     thrdat.qel_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOM_EV_DIM);
                     thrdat.qel_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ORBMOM_EV_DIM);
                     thrdat.Eorbmom[ithread] = new ComplexMatrix(1,dimA,1,CHARGEDENS_EV_DIM); *thrdat.Eorbmom[ithread]=Eorbmom;

                     thrdat.Tau[ithread] = new ComplexMatrix(1,dimA,1,dimA); *thrdat.Tau[ithread]=Tau;
                  }
                  ithread=0; num_threads_started=-1; int oldi=-1; double QQ;// Vector vQQ(1,dimA); removed MR 14.1.2013
#endif
                      Vector abc(1,6); abc(1)=inputpars.a; abc(2)=inputpars.b; abc(3)=inputpars.c;
                                       abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;
                      qold=qijk;hkl2ijk(qijk,hkl, abc);QQ=Norm(qijk);qincr+=Norm(qijk-qold); 

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
                     for(ithread=0; ithread<NUM_THREADS; ithread++)
                     {
                        i=oldi+ithread; if(i>dimA) break;
                        if(do_gobeyond==0) intsbey(i)=-1.1; else intsbey(i)=+1.1;
                      if (En(i)<=ini.emax&&En(i)>=ini.emin) // only do intensity calculation if within energy range
                      {
                        tin[ithread]->En=En(i); tin[ithread]->intensitybey=intsbey(i); tin[ithread]->level = i;
                        #if defined  (__linux__) || defined (__APPLE__)
                        rc = pthread_create(&threads[ithread], &attr, intcalc_approx, (void *) tin[ithread]);
                        if(rc) { printf("Error return code %i from thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
                        #else
                        threads[ithread] = CreateThread(NULL, 0, intcalc_approx, (void *) tin[ithread], 0, &tid[ithread]);
                        if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %i from thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
                        #endif
                        num_threads_started = ithread+1;
                      }
                      else {ints(i)=-1;intsbey(i)=-1;}
                     }
                     #if defined  (__linux__) || defined (__APPLE__)
                     for(int th=0; th<num_threads_started; th++)
                        rc = pthread_join(threads[th], &status);
                     #else
                     WaitForMultipleObjects(num_threads_started,threads,TRUE,INFINITE);
                     #endif
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
                     for(ithread=0; ithread<NUM_THREADS; ithread++)
                     {
                         i=oldi+ithread; if(i>dimA) break;
                        if (En(i)<=ini.emax&&En(i)>=ini.emin) // only do intensity calculation if within energy range
                      {ints(tin[ithread]->level) = tin[ithread]->intensity; 
                       intsbey(tin[ithread]->level) = tin[ithread]->intensitybey;
                      }       
#else
                     if(do_gobeyond==0){intsbey(i)=-1.1;}else{intsbey(i)=+1.1;}
                     if (En(i)<=ini.emax&&En(i)>=ini.emin) // only do intensity calculation if within energy range
                     {
                     ints(i)=intcalc_approx(chi,chibey,pol,intsbey(i),
                                            qee_real,qee_imag,Echargedensity,
                                            qsd_real,qsd_imag,Espindensity,
                                            qod_real,qod_imag,Eorbmomdensity,
                                            qep_real,qep_imag,Ephonon,
                                            qem_real,qem_imag,Emagmom,
                                            qes_real,qes_imag,Espin,
                                            qel_real,qel_imag,Eorbmom,
                                            dimA,Tau,i,En(i),ini,inputpars,hkl,md,do_verbose,QQ);
                     }
                     else
                     {ints(i)=-1;intsbey(i)=-1;
                     }
#endif
                     if (ini.hkllist==1)
	             {double test; // add to sta distance to nearest measured peak squared
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
                       if ((test=fabs(En(i)-ini.hkls[counter][4*j1]))<dd1_int(j1)&&ints(i)>SMALLINT){dd1_int(j1)=test;
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
	             }

                     if(intsbey(i)<0)intsbey(i)=-1;
                     //printout rectangular function to .mdcisp.qom
	             fprintf (fout, " %4.4g %4.4g",myround(ints(i)),myround(intsbey(i)));
                     fprintf (foutqei, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g  %4.4g  %4.4g %4.4g\n",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)),
                                         myround(QQ),myround(En(i)),myround(1e-8,ints(i)),myround(1e-8,intsbey(i)),qincr);
                     // printout eigenvectors only if evaluated during intensity calculation...
                  if(ints(i)>-1){

if(ini.calculate_chargedensity_oscillation)print_ev(foutqee,i,ini,hkl,QQ,En,ints,intsbey,qee_real,qee_imag);
if(ini.calculate_spindensity_oscillation)print_ev(foutqsd,i,ini,hkl,QQ,En,ints,intsbey,qsd_real,qsd_imag);
if(ini.calculate_orbmomdensity_oscillation)print_ev(foutqod,i,ini,hkl,QQ,En,ints,intsbey,qod_real,qod_imag);
if(ini.calculate_phonon_oscillation)print_ev(foutqep,i,ini,hkl,QQ,En,ints,intsbey,qep_real,qep_imag);
if(ini.calculate_magmoment_oscillation)print_ev(foutqem,i,ini,hkl,QQ,En,ints,intsbey,qem_real,qem_imag);
if(ini.calculate_spinmoment_oscillation)print_ev(foutqes,i,ini,hkl,QQ,En,ints,intsbey,qes_real,qes_imag);
if(ini.calculate_orbmoment_oscillation)print_ev(foutqel,i,ini,hkl,QQ,En,ints,intsbey,qel_real,qel_imag);
                          }
                 if(do_verbose==1){fprintf(stdout, "#IdipFF= %4.4g Ibeyonddip=%4.4g\n",ints(i),intsbey(i));}
                     if(En(i)>=ini.emin&&En(i)<=ini.emax){diffint+=ints(i);diffintbey+=intsbey(i);}
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
                  for (ithread=0; ithread<NUM_THREADS; ithread++) 
                  {
                     delete thrdat.chi[ithread]; delete thrdat.chibey[ithread]; 
                     delete thrdat.pol[ithread]; 
                     delete thrdat.qee_real[ithread]; delete thrdat.qee_imag[ithread];delete thrdat.Echargedensity[ithread]; 
                     delete thrdat.qsd_real[ithread]; delete thrdat.qsd_imag[ithread];delete thrdat.Espindensity[ithread]; 
                     delete thrdat.qod_real[ithread]; delete thrdat.qod_imag[ithread];delete thrdat.Eorbmomdensity[ithread]; 
                     delete thrdat.qep_real[ithread]; delete thrdat.qep_imag[ithread];delete thrdat.Ephonon[ithread]; 
                     delete thrdat.qem_real[ithread]; delete thrdat.qem_imag[ithread];delete thrdat.Emagmom[ithread]; 
                     delete thrdat.qes_real[ithread]; delete thrdat.qes_imag[ithread];delete thrdat.Espin[ithread]; 
                     delete thrdat.qel_real[ithread]; delete thrdat.qel_imag[ithread];delete thrdat.Eorbmom[ithread]; 
                     delete thrdat.Tau[ithread]; delete tin[ithread]; 
                  }
                  delete[] thrdat.Echargedensity;  
                  delete[] thrdat.Espindensity;  
                  delete[] thrdat.Eorbmomdensity;  
                  delete[] thrdat.Ephonon;  
                  delete[] thrdat.Emagmom;  
                  delete[] thrdat.Espin;  
                  delete[] thrdat.Eorbmom;  
                  delete[] thrdat.Tau;
                  delete[] thrdat.chi; delete[] thrdat.chibey;
                  delete[] thrdat.pol; 
                  delete[] thrdat.qee_real; delete[] thrdat.qee_imag; 
                  delete[] thrdat.qsd_real; delete[] thrdat.qsd_imag; 
                  delete[] thrdat.qod_real; delete[] thrdat.qod_imag; 
                  delete[] thrdat.qep_real; delete[] thrdat.qep_imag; 
                  delete[] thrdat.qem_real; delete[] thrdat.qem_imag; 
                  delete[] thrdat.qes_real; delete[] thrdat.qes_imag; 
                  delete[] thrdat.qel_real; delete[] thrdat.qel_imag; 
#endif
    fprintf (foutdstot, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %4.4g %4.4g",ini.Hext(1),ini.Hext(2),ini.Hext(3),ini.T,hkl(1),hkl(2),hkl(3),diffint,diffintbey,qincr);
    sta+=dd*dd;sta_int+=dd_int*dd_int;
    sta_without_antipeaks+=dd_without_antipeaks*dd_without_antipeaks;
    sta_int_without_antipeaks+=dd_int_without_antipeaks*dd_int_without_antipeaks;
    sta_without_weights+=dd_without_weights*dd_without_weights;
    sta_int_without_weights+=dd_int_without_weights*dd_int_without_weights;
    sta_without_antipeaks_weights+=dd_without_antipeaks_weights*dd_without_antipeaks_weights;
    sta_int_without_antipeaks_weights+=dd_int_without_antipeaks_weights*dd_int_without_antipeaks_weights;

              //initialize output file for display
            errno = 0;
            fout1 = fopen_errchk ("./results/.mcdisp.qom","w");
            fprintf (fout1, "#{%s ",MCDISPVERSION);
            curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout1);
            fprintf (fout1, "#displayytext=I(barns/meV/sr/f.u.)\n");
            fprintf (fout1, "#displayxtext=E(meV)\n");
            fprintf (fout1, "#displaytitle=(%4.4f %4.4f %4.4f) blue: DMD_Dipapprox red: DMD_exact green: Minv_Dipapprox\n",hkl(1),hkl(2),hkl(3));
            fprintf (fout1,"#Ha[T] Hb[T] Hc[T] T[K] h k l  energies[meV] intensities(dip approx for FF) [barn/meV/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
		     if (do_Erefine==0) epsilon=(Max(En)-Min(En))/100;
		     if (epsilon<=0) epsilon=0.1;
                  for (i=1;i<=dimA;++i)
		    { 
		     if (ints(i)>SMALLINT)  // draw triangles to show calculated intensity
		      {for (E=0;E<=ints(i)/epsilon;E+=ints(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	                fprintf (fout1, " %4.4g %4.4g %4.4g 0\n",myround(En(i)-epsilon+E*epsilon*epsilon/ints(i)),myround(E),myround(En(i)));
		       }
		       for (E=ints(i)/epsilon;E>=0;E-=ints(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	                fprintf (fout1, " %4.4g %4.4g %4.4g 0\n",myround(En(i)+epsilon-E*epsilon*epsilon/ints(i)),myround(E),myround(En(i)));
		       }
                       fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	               fprintf (fout1, " %4.4g 0 %4.4g 0\n",myround(En(i)+epsilon),myround(En(i)));
		      }
		     if (intsbey(i)>SMALLINT)  // draw triangles to show calculated intensity
		      {for (E=0;E<=intsbey(i)/epsilon;E+=intsbey(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	                fprintf (fout1, " %4.4g 0 %4.4g %4.4g \n",myround(En(i)),myround(En(i)-epsilon+E*epsilon*epsilon/intsbey(i)),E);
		       }
		       for (E=intsbey(i)/epsilon;E>=0;E-=intsbey(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	                fprintf (fout1, " %4.4g 0 %4.4g %4.4g\n",myround(En(i)),myround(En(i)+epsilon-E*epsilon*epsilon/intsbey(i)),E);
		       }
                       fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	               fprintf (fout1, " %4.4g 0 %4.4g 0 \n",myround(En(i)),myround(En(i)+epsilon));
		      }
		    }
	  fclose(fout1);   
	    
		    
                if(do_verbose==1){fprintf(stdout, "\n");}
		    
   // do refinement of energies by output of scattering cross section vs enrgy transfer if required
  if (do_Erefine==1){double totint=0;
                if(do_verbose==1){fprintf(stdout, "#refining calculation with exact calculation of energy dependence of scattering cross section\n");}
          errno = 0;
          foutds1 = fopen_errchk ("./results/.mcdisp.dsigma","w");
          fprintf (foutds1, "#{%s ",MCDISPVERSION);
          curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutds1);
          fprintf (foutds1, "#Scattering Cross Section \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energy[meV] dsigma/dOmegadE'[barn/mev/sr/f.u.] (dipolar approx for FF) f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
#ifdef _THREADS
          thrdat.J = new jq*[NUM_THREADS]; thrdat.hkl = hkl; thrdat.q = q; thrdat.thread_id = -1;
          for (int ithread=0; ithread<NUM_THREADS; ithread++) 
          {
             tin[ithread] = new intcalcapr_input(dimA,ithread,1,do_verbose,0.); tin[ithread]->epsilon=epsilon; thrdat.J[ithread] = new jq(J);
          }
          int ithread=0; double oldE=0.;
          Vector vIntensity(1,(int)((ini.emax-ini.emin)/(epsilon/2)+1)); int iE=1;
#endif
		     double intensity;
#ifdef _THREADS
                     tin[ithread]->En=ini.emin; tin[ithread]->iE=iE;
                     #if defined  (__linux__) || defined (__APPLE__)
                     rc = pthread_create(&threads[ithread], &attr, intcalc, (void *) tin[ithread]); rc = pthread_join(threads[ithread], &status); 
                     if(rc) { printf("Error return code %i from joining thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
                     #else
                     threads[ithread] = CreateThread(NULL, 0, intcalc, (void *) tin[ithread], 0, &tid[ithread]);
                     if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %i from thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
                     if(WaitForSingleObject(threads[ithread],INFINITE)==0xFFFFFFFF) { printf("Error in waiting for thread %i to end\n",ithread+1); exit(EXIT_FAILURE); }
                     #endif
                     intensity=tin[ithread]->intensity;
#else
		     intensity=intcalc(dimA,ini.emin,ini,inputpars,J,q,hkl,md,do_verbose,epsilon);   
#endif
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	             fprintf (foutds1, " %4.4g %4.4g \n",ini.emin,myround(intensity));
#ifdef _THREADS
                     tin[ithread]->En=ini.emax; tin[ithread]->iE=iE;
                     #if defined  (__linux__) || defined (__APPLE__)
                     rc = pthread_create(&threads[ithread], &attr, intcalc, (void *) tin[ithread]); rc = pthread_join(threads[ithread], &status); 
                     if(rc) { printf("Error return code %i from thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
                     #else
                     threads[ithread] = CreateThread(NULL, 0, intcalc, (void *) tin[ithread], 0, &tid[ithread]);
                     if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %i from thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
                     if(WaitForSingleObject(threads[ithread],INFINITE)==0xFFFFFFFF) { printf("Error in waiting for thread %i to end\n",ithread+1); exit(EXIT_FAILURE); }
                     #endif
                     intensity=tin[ithread]->intensity;
#else
		     intensity=intcalc(dimA,ini.emax,ini,inputpars,J,q,hkl,md,do_verbose,epsilon);   
#endif
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	             fprintf (foutds1, " %4.4g %4.4g \n",ini.emax,myround(intensity));
	  fclose(foutds1);
#ifdef _THREADS
	  for(E=ini.emin;E<=ini.emax;E+=(epsilon/2)*NUM_THREADS)
#else
	  for(E=ini.emin;E<=ini.emax;E+=epsilon/2)
#endif
	   {
#ifndef _THREADS
		     intensity=intcalc(dimA,E,ini,inputpars,J,q,hkl,md,do_verbose,epsilon);   
		     totint+=intensity*epsilon/2;
#else
                     oldE=E;
                     for(ithread=0; ithread<NUM_THREADS; ithread++)
                     {
                        E=oldE+(epsilon/2)*ithread; if(E>ini.emax) break;
                        tin[ithread]->En=E; tin[ithread]->iE=iE++;
                        #if defined  (__linux__) || defined (__APPLE__)
                        rc = pthread_create(&threads[ithread], &attr, intcalc, (void *) tin[ithread]);
                        if(rc) { printf("Error return code %i from thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
                        #else
                        threads[ithread] = CreateThread(NULL, 0, intcalc, (void *) tin[ithread], 0, &tid[ithread]);
                        if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %i from thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
                        #endif
                        num_threads_started = ithread+1;
                     }
                     #if defined  (__linux__) || defined (__APPLE__)
                     for(int th=0; th<num_threads_started; th++)
                        rc = pthread_join(threads[th], &status);
                     #else
                     WaitForMultipleObjects(num_threads_started,threads,TRUE,INFINITE);
                     #endif
                     for(ithread=0; ithread<NUM_THREADS; ithread++) vIntensity(tin[ithread]->iE) = tin[ithread]->intensity;
                     E=oldE;
	   }
	  iE=1; for(E=ini.emin;E<=ini.emax;E+=epsilon/2)
	   {
                     intensity = vIntensity(iE++); totint+=intensity*epsilon/2;
#endif

          foutds1 = fopen_errchk ("./results/.mcdisp.dsigma","a");
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	             fprintf (foutds1, " %4.4g %4.4g \n",myround(E),myround(intensity));
          fclose(foutds1);	   
                     fprintf (foutds, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)));
	             fprintf (foutds, " %4.4g %4.4g \n",myround(E),myround(intensity));
	   }

#ifdef _THREADS
          for (ithread=0; ithread<NUM_THREADS; ithread++) 
          {
             delete thrdat.J[ithread]; delete tin[ithread];
          }
          delete[] thrdat.J; //delete tin;
#endif
	             fprintf (foutdstot, " %4.4g ",totint);
                     }  

   fprintf (foutdstot, "\n");              
   fprintf (fout, "\n");
   }
#ifdef _THREADS
   for (ithread=0; ithread<NUM_THREADS; ithread++) 
   {
      delete thrdat.md[ithread]; delete thrdat.ini[ithread]; delete thrdat.inputpars[ithread]; 
   }
   delete[] thrdat.inputpars; delete[] thrdat.md; delete[] thrdat.ini;
#endif

   if (ini.hkllist==1){hkl(1)=(double)counter;}
               
}}} // next hkl
    if (do_jqfile==1) 
     {fprintf(jqfile,"#it follows the standard deviation sta defined as:\n");
      fprintf(jqfile,"#the sum of squared differences between the highest eigenvalue\n");
      fprintf(jqfile,"#of a q vector and that of the first q-vector in the list in mcdisp.par.\n");
      fprintf(jqfile,"#only those eigenvalues are taken into account in the sum, which are larger\n");
      fprintf(jqfile,"#than that of the first q-vector in the list in mcdisp.par - this is usefule\n");
      fprintf(jqfile,"#for obtaining an exchange interaction with maximum at the first q-vector\n");
      fprintf(jqfile,"#in the list in mcdisp.par\n");
      fprintf(jqfile,"# ... if the first q vector has the largest eigenvalue, then sta is negative and contains the\n");
      fprintf(jqfile,"#distance to the closest eigenvalue\n");
      fprintf(jqfile,"#!sta=%g\n",jqsta);
      fprintf(jqfile,"#another standard deviation is given below: calculated as squared sum of differences between\n");
      fprintf(jqfile,"#the highest eigenvalue of J(Q) and energies in column 4 of mcdisp.par, if column 5 6 etc \n");
      fprintf(jqfile,"#in mcdisp.par contain values, then these are compared to the other eigenvalues of J(Q)\n");
      fprintf(jqfile,"#!sta4=%g\n",jqsta_int);fclose(jqfile);}
    else
     {

      staout(fout,sta,sta_int,sta_without_antipeaks,sta_int_without_antipeaks,sta_without_weights,sta_int_without_weights,sta_without_antipeaks_weights,sta_int_without_antipeaks_weights);
      staout(foutqei,sta,sta_int,sta_without_antipeaks,sta_int_without_antipeaks,sta_without_weights,sta_int_without_weights,sta_without_antipeaks_weights,sta_int_without_antipeaks_weights);
      staout(stdout,sta,sta_int,sta_without_antipeaks,sta_int_without_antipeaks,sta_without_weights,sta_int_without_weights,sta_without_antipeaks_weights,sta_int_without_antipeaks_weights);

    fclose(foutqei);

                        if(ini.calculate_chargedensity_oscillation)fclose(foutqee);
                        if(ini.calculate_spindensity_oscillation)fclose(foutqsd);
                        if(ini.calculate_orbmomdensity_oscillation)fclose(foutqod);
                        if(ini.calculate_phonon_oscillation)fclose(foutqep);
                        if(ini.calculate_magmoment_oscillation)fclose(foutqem);
                        if(ini.calculate_spinmoment_oscillation)fclose(foutqes);
                        if(ini.calculate_orbmoment_oscillation)fclose(foutqel);
                        
    fclose(fout);
                 if (do_Erefine==1){fclose(foutds);}
    fclose(foutdstot);
     }
}

//*************************************************************************************************
// main program
int main (int argc, char **argv)
{int i,do_Erefine=0,do_jqfile=0,do_verbose=0,maxlevels=10000000,do_createtrs=0,do_readtrs=0,calc_beyond=1,calc_fast=0;
 const char * spinfile="mcdisp.mf"; //default spin-configuration-input file
 const char * filemode="w";
 double epsilon; //imaginary part of omega to avoid divergence
 double minE=-100000.0,maxE=+100000.0,pinit=0,ninit=1e10;
 fprintf(stderr,"***********************************************************************\n");
 fprintf(stderr,"*\n");
 fprintf(stderr,"* mcdisp - program to calculate the dispersion of magnetic excitations\n");
 fprintf(stderr,"*\n");
 fprintf(stderr,"* reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
 fprintf(stderr,"*            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
 fprintf(stderr,"***********************************************************************\n\n");


// check command line and initialize parameters ini
for (i=1;i<=argc-1;++i){
   if(strcmp(argv[i],"-r")==0) {do_Erefine=1; if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -r needs argument epsilon\n");exit(EXIT_FAILURE);}
		                                                epsilon=strtod(argv[i+1],NULL);++i;
							        fprintf(stdout,"#epsilon= %g\n",epsilon);
				     }		
       else {if(strcmp(argv[i],"-d")==0) {calc_beyond=0;}
        else {if(strcmp(argv[i],"-f")==0) {calc_fast=0;} // calc_fast=1: MR 3.1.2013 attempt to speed up DMD by halfing matrix size, however this does not work - is there a possibility ?
         else {if(strcmp(argv[i],"-jq")==0) {do_jqfile=1;minE=SMALL;maxlevels=1;}
          else {if(strcmp(argv[i],"-t")==0) do_readtrs=1;       
           else {if(strcmp(argv[i],"-c")==0) do_createtrs=1;       
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
           	   else{spinfile=argv[i];}
		     }
		    }
		   }
		  }          
		 }
		}
	      }
	     }
	    }
           }	
          }
         }
    }
inimcdis ini("mcdisp.par",spinfile);

if (argc > 10) {ini.errexit();}
  // as class load  parameters from file
  par inputpars("./mcphas.j");

  if(ini.nofcomponents!=inputpars.nofcomponents){fprintf(stderr,"Error mcdisp: number of components read from mcdisp.par (%i) and mcphas.j (%i) not equal\n",ini.nofcomponents,inputpars.nofcomponents);exit(1);}
  if(ini.nofatoms!=inputpars.nofatoms){fprintf(stderr,"Error mcdisp: number of atoms in crystal unit cell read from mcdisp.par (%i) and mcphas.j (%i) not equal\n",ini.nofatoms,inputpars.nofatoms);exit(1);}
  inputpars.save("./results/_mcdisp.j");
  inputpars.save_sipfs("./results/_"); 


//calculate dispersion and save to files
dispcalc(ini,inputpars,calc_fast,calc_beyond,do_Erefine,do_jqfile,do_createtrs,do_readtrs,do_verbose,maxlevels,minE,maxE,ninit,pinit,epsilon,filemode);
  
 printf("RESULTS saved in directory ./results/  - files:\n");
   printf("  mcdisp.qei  - T,H,qvector vs energies and neutron intensities\n");
   printf("  mcdisp.qom  - T,H,qvector vs all mode energies in one line (and neutron intensities)\n");
   printf("  mcdisp.qee,qsd,qod,qep,qem,qes,qel  - T,H,qvector,E vs extended eigenvectors (more components to plot observables.)\n");
   printf("  mcdisp.dsigma.tot  - T,H,qvector vs total intensity (sum of all modes)\n");
   printf("  mcdisp.dsigma      - (option -r) T,H,qvector,E vs intensity obtained from dyn susz\n");
   printf("  mcdisp.trs  - single ion transitions used\n");
   printf("  _mcdisp.par - input parameters read from mcdisp.par\n");
   printf("  _mcdisp.mf  - input parameters read from mcdisp.mf\n");
   printf("  _mcdisp.j   - input parameters read from mcphas.j\n");
   printf("  ...         - and a copy of the single ion parameter files used.\n\n");
   fprintf(stderr,"************************************************************\n");
   fprintf(stderr,"                    End of Program mcdisp\n");
   fprintf(stderr," reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(stderr,"            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(stderr,"************************************************************\n");

 return(0);
 
}





