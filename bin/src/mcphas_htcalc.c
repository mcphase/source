// routines for mcphas for calculation of magnetic phases
// htcalc.c


#ifdef _THREADS
#if defined  (__linux__) || defined (__APPLE__)
#include <pthread.h>
#define MUTEX_LOCK     pthread_mutex_lock
#define MUTEX_UNLOCK   pthread_mutex_unlock
#define MUTEX_TYPE     pthread_mutex_t
#define MUTEX_INIT(m)  pthread_mutex_init (&m, NULL)
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
#define MUTEX_INIT(m)  InitializeCriticalSection (&m)
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
#define NUM_THREADS ini.nofthreads

// ----------------------------------------------------------------------------------- //
// Declares a struct to store all the information needed for each htcalc iteration
// ----------------------------------------------------------------------------------- //
typedef struct{
   Vector Happ;
   double T;
   qvectors * testqs;  
   testspincf * testspins;
   inipar * ini;
   physproperties * physprops; 
   double femin;
   spincf spsmin;
   int thread_id;
} htcalc_thread_data;
class htcalc_input { public:
   int j; 
   int thread_id;
   par *inputpars;
   htcalc_input(int  _j, int & _tid, par *pars_in) 
   { 
      thread_id = _tid; j = _j; inputpars = new par(*pars_in);
   }
   ~htcalc_input(){
    delete inputpars;}
};
// ----------------------------------------------------------------------------------- //
// Declares these variables global, so all threads can see them
// ----------------------------------------------------------------------------------- //
htcalc_thread_data thrdat;
htcalc_input *tin[256];  // Max number of threads - hard coded because global variable.
MUTEX_TYPE mutex_loop;
MUTEX_TYPE mutex_tests;
MUTEX_TYPE mutex_min;
MUTEX_TYPE mutex_ini_nofcalls;
MUTEX_TYPE mutex_ini_nofmaxspinchangeDIV;
MUTEX_TYPE mutex_ini_nofmaxloopDIV;
MUTEX_TYPE mutex_ini_successrate;
MUTEX_TYPE mutex_ini_calcmf_duration;
MUTEX_TYPE mutex_ini_calcsps_duration;

EVENT_TYPE checkfinish;
THRLC_TYPE threadSpecificKey;

#endif // def _THREADS



#ifdef _THREADS
#define ini (*thrdat.ini)
#define physprops (*thrdat.physprops)
#define inputpars (*myinput->inputpars)
#define testqs (*thrdat.testqs)
#define testspins (*thrdat.testspins)
#define T thrdat.T
#define femin thrdat.femin
#if defined  (__linux__) || defined (__APPLE__)
void *htcalc_iteration(void *input)
#else
DWORD WINAPI htcalc_iteration(void *input)
#endif
#else
int htcalc_iteration(int & j, double &femin, spincf &spsmin, Vector Happ, double T,inipar & ini, par &inputpars, qvectors &testqs, testspincf &testspins, physproperties &physprops)
#endif
{
 #ifdef _THREADS
 htcalc_input *myinput; myinput = (htcalc_input *) input; int j = myinput->j, thread_id = myinput->thread_id; Vector Happ(1,HEXT_DIMENSION); Happ = thrdat.Happ;
 THRLC_SET(threadSpecificKey, myinput); int tlsfemin=0;  // Thread local variable to judge whether to print output
 #else
 int thread_id=1;
 #endif 
 { // <--- this bracket is necessary to embrace ComplexMatrix  definitions
  // necessary in this routine and make them live within this bracket. at closing the bracket
  // destructor is called correctly. On Apple without this 
  // bracket the compiler will exit intcalc_Erefine (pthread_exit) without freeing memory ...
 
 int i,ii,iii,tryrandom,nr,rr,ri,is,r;
 double fe,fered,Eelastic,U,sc;
 double u,lnz; // free- and magnetic energy per ion [meV]
 Vector momentq0(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms),phi(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
 Vector nettom(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms),q(1,3);
 Vector mmom(1,inputpars.cs.nofcomponents);
 Vector h1(1,inputpars.cs.nofcomponents),Happzero(1,Happ.Hi()),hkl(1,3);
 Happzero=0;
 char text[MAXNOFCHARINLINE];
 char outfilename[MAXNOFCHARINLINE];
 spincf  sps(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents),sps1(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 FILE * felog; // logfile for q dependence of fe
 FILE * fout;
int s1=1,s2=2;
  for (tryrandom=0;(tryrandom<=ini.nofrndtries)&&j!=0;++tryrandom)
   {if (j>0){sps=(*testspins.configurations[j]);// take test-spinconfiguration
             #ifndef _THREADS
	     if (tryrandom==0&&verbose==1) { printf ( "str%i(%ix%ix%i)< "  ,j,sps.na(),sps.nb(),sps.nc()); fflush(stdout); }
             #else
	     if (tryrandom==0&&verbose==1) { printf ( "str%i(%ix%ix%i)<[%i] "  ,j,sps.na(),sps.nb(),sps.nc(),thread_id+1); fflush(stdout); }
             #endif 
            
            while(sps.na()<ini.minnr1)sps.extend(s2,s1,s1); 
            while(sps.nb()<ini.minnr2)sps.extend(s1,s2,s1); 
            while(sps.nc()<ini.minnr3)sps.extend(s1,s1,s2); 
            }
    else     // take q vector and choose phase and mom dir randomly
            {int mj=-j;q=testqs.q(mj);  
	     if (tryrandom==0)
	     {nettom=testqs.nettom(mj);momentq0=testqs.momentq0(mj);phi=testqs.phi(mj);
	     }
	     else
	     {for(i=1;i<=inputpars.cs.nofatoms;++i)
	      {for(ii=1;ii<=inputpars.cs.nofcomponents;++ii)
	        {iii=inputpars.cs.nofcomponents*(i-1)+ii;h1=0;h1(ii)=10*MU_B;
                 (*inputpars.jjj[i]).Icalc(mmom,T,h1,Happzero,lnz,u,(*inputpars.jjj[i]).Icalc_parstorage);
		 nettom(iii)=mmom(ii)*rnd(1);
	         momentq0(iii)=rnd(1)*mmom(ii);
	         phi(iii)=rnd(1)*3.1415;
		}
	      }
	     }int ta=testqs.na(mj),tb=testqs.nb(mj),tc=testqs.nc(mj);
	     sps.spinfromq(ta,tb,tc,
	                   q,nettom,momentq0,phi);
             while(sps.na()<ini.minnr1)sps.extend(s2,s1,s1); 
             while(sps.nb()<ini.minnr2)sps.extend(s1,s2,s1); 
             while(sps.nc()<ini.minnr3)sps.extend(s1,s1,s2); 
                // for Monte Carlo we have minimum number of spins
             hkl=inputpars.rez.Transpose()*q;  
             #ifndef _THREADS
   	     if (tryrandom==0&&verbose==1) { printf ( "(%g %g %g)(%ix%ix%i)< ",hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc()); fflush(stdout); }
             #else
   	     if (tryrandom==0&&verbose==1) { printf ( "(%g %g %g)(%ix%ix%i)<[%i] ",hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc(),thread_id+1); fflush(stdout); }
             #endif 
	    }	 
    if (tryrandom>0){
                      int g1=sps.n()*inputpars.cs.nofatoms;nr=rndint(g1);
	             for (i=1;i<=nr;++i) // randomize nr spins
                      {g1=sps.n();rr=rndint(g1);
		       ri=inputpars.cs.nofcomponents*arc4random_uniform(inputpars.cs.nofatoms);
	               for(ii=1;ii<=inputpars.cs.nofcomponents;++ii)
		       {sps.mi(rr)(ri+ii)*=(2*rnd(1.0)-1) ;}
		       } // randomize spin rr
                      
                    }
 
      //!!!calculate free energy - this is the heart of this loop !!!!
      mfcf mf(sps.na(),sps.nb(),sps.nc(),inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
      fe=fecalc(U,Eelastic,r,sc,Happ ,T,ini,inputpars,sps,mf);
          if (fe>=2*FEMIN_INI && verbose==1) {
	       if(j>0) printf ( ">for_str_%i(%ix%ix%i) "  ,j,sps.na(),sps.nb(),sps.nc());
               else    printf ( ">for(%g %g %g)(%ix%ix%i) ",hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc()); 
                              fflush(stdout);}
      // test spinconfiguration  and remember it                                    
      if (fe<femin)
            {if(ini.nofMCsteps>0)
              {
             #ifndef _THREADS
	       femin=fe; spsmin=sps;	   
               //printout fe
	        if (verbose==1) printf("fe=%gmeV, str %i(%i)",fe,physprops.j,j);
             #else
             MUTEX_LOCK (&mutex_min); if(fe<femin) { femin=fe; thrdat.spsmin=sps; } MUTEX_UNLOCK (&mutex_min); tlsfemin=1;
             #endif 
               }
             else
             {  
               // first - reduce the spinconfiguration if possible
               sps1=sps;if(1==sps1.reduce()){ // if reduction is successful, try if the energy is less or equal for reduced spoinconfigurations
                   mfcf mf1(sps1.na(),sps1.nb(),sps1.nc(),inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
               if ((fered=fecalc(U,Eelastic,r,sc,Happ ,T,ini,inputpars,sps1,mf1))<=fe*(1.0000000000001)){mf=mf1;
                                 if (verbose==1){fprintf(stdout,">[%i](%ix%ix%i)r%i->(%ix%ix%i)fe=%f->%fmeV ",thread_id+1,sps.na(),sps.nb(),sps.nc(),tryrandom,sps1.na(),sps1.nb(),sps1.nc(),fe,fered); fflush(stdout);}
                                                                                     sps=sps1;fe=fered;}
                                                                                                  else {
                                 if (verbose==1){fprintf(stdout,">[%i](%ix%ix%i)r%ife=%.15gmeV<(%ix%ix%i)fered=%.15g ",thread_id+1,sps.na(),sps.nb(),sps.nc(),tryrandom,fe,sps1.na(),sps1.nb(),sps1.nc(),fered);fflush(stdout);}
                                                                                                       }
                                                }
                                       else {
                                 if (verbose==1){fprintf(stdout,">[%i](%ix%ix%i)r%ife=%fmeV ",thread_id+1,sps.na(),sps.nb(),sps.nc(),tryrandom,fe);fflush(stdout);}
                                            }
                    spincf magmom(sps.na(),sps.nb(),sps.nc(),inputpars.cs.nofatoms,3);
                   int i1,j1,k1,l1,m1;Vector mom(1,3),d1(1,inputpars.cs.nofcomponents);
                   for (l1=1;l1<=inputpars.cs.nofatoms;++l1){
                    // go through magnetic unit cell and sum up the contribution of every atom
                  for(i1=1;i1<=sps.na();++i1){for(j1=1;j1<=sps.nb();++j1){for(k1=1;k1<=sps.nc();++k1){
                   for(m1=1;m1<=inputpars.cs.nofcomponents;++m1){d1[m1]=mf.mf(i1,j1,k1)[inputpars.cs.nofcomponents*(l1-1)+m1];}
                   (*inputpars.jjj[l1]).mcalc(mom,T,d1,Happ,(*inputpars.jjj[l1]).Icalc_parstorage);
                   for(m1=1;m1<=3;++m1){magmom.m(i1,j1,k1)(3*(l1-1)+m1)=mom(m1);}
                    }}}} 
                  
                 // display spinstructure
                if (verbose==1)
                {float * x;x=new float[inputpars.cs.nofatoms+1];float *y;y=new float[inputpars.cs.nofatoms+1];float*z;z=new float[inputpars.cs.nofatoms+1];
		 
		 for (is=1;is<=inputpars.cs.nofatoms;++is)
		   {
                    x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
                     snprintf(text,MAXNOFCHARINLINE,"fe=%g<femin=%g:T=%gK, |H|=%gT,Ha=%gT, Hb=%gT, Hc=%gT,  %i spins",fe,femin,T,Norm(Happ),Happ(1),Happ(2),Happ(3),sps.n());
                    strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dab.eps");
                    fout = fopen_errchk (outfilename, "w");
                     sps.eps3d(fout,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,4,magmom);
                    fclose (fout);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dac.eps");
                    fout = fopen_errchk (outfilename, "w");
                     sps.eps3d(fout,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,5,magmom);
                    fclose (fout);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dbc.eps");
                    fout = fopen_errchk (outfilename, "w");
                     sps.eps3d(fout,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,6,magmom);
                    fclose (fout);
		   
                    strcpy(outfilename+11+strlen(ini.prefix),"spins.eps");
                    fout = fopen_errchk (outfilename, "w");
                     magmom.eps(fout,text);
                    fclose (fout);
		delete[]x;delete []y; delete []z;
	        }
               // printf("C");fflush(stdout);delete magmom;   
                           // see if spinconfiguration is already stored
             #ifndef _THREADS
	     if (0==checkspincf(j,sps,testqs,nettom,momentq0,phi,testspins,physprops,ini))//0 means error in checkspincf/addspincf
	        {if(isfull==0){fprintf(stderr,"Warning !FT! htcalc: table of spinconfigurations full - cannot add a new configuration, which has been found.");
                 isfull=1;}else{fprintf(stderr,"!FT!");}
                physprops.j=j; 
                }
	     femin=fe; spsmin=sps;	   
            //printout fe
	    if (verbose==1) printf("fe=%gmeV, str %i(%i)",fe,physprops.j,j);
             #else
             MUTEX_LOCK(&mutex_tests); 
             int checksret = checkspincf(j,sps,testqs,nettom,momentq0,phi,testspins,physprops,ini); //0 means error in checkspincf/addspincf
             MUTEX_UNLOCK(&mutex_tests); 
	     if (checksret==0) {physprops.j=j;
                 if(isfull==0){fprintf(stderr,"%iWarning !FT! htcalc: table of spinconfigurations full - cannot add a new configuration, which has been found.",thread_id);
                 isfull=1;}else{fprintf(stderr,"%i!FT!",thread_id);}}
             MUTEX_LOCK (&mutex_min); if(fe<femin) { femin=fe; thrdat.spsmin=sps; } MUTEX_UNLOCK (&mutex_min); tlsfemin=1;
             #endif 
	     }}
            //delete mf;
             //printout fe
            #ifdef _THREADS
	    if (tryrandom==ini.nofrndtries)if(verbose==1) {
               if(tlsfemin) printf("[%i]femin=%gmeV str %i(%i)-",thread_id+1,femin,physprops.j,j); 
	       if(j>0) printf ( ">[%i]str %i(%ix%ix%i)done "  ,thread_id+1,j,sps.na(),sps.nb(),sps.nc());
               else    printf ( ">[%i](%g %g %g)(%ix%ix%i)done ",thread_id+1,hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc()); 
                                                          }
            #endif
            if (tryrandom==ini.nofrndtries)if(verbose==1){printf("\n");}
 
	    
  // log fe if required
   if (ini.logfevsQ==1) {
                 ComplexVector a(1,3*inputpars.cs.nofatoms),b(1,3*inputpars.cs.nofatoms);
                 ComplexVector b1(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
                 float inmax=0;int qh,qk,ql,l,nk=0;
                 ComplexVector * mq;  
                 int na=sps.na(),nb=sps.nb(),nc=sps.nc();
                 mq = new ComplexVector [sps.in(na,nb,nc)+2];
                 for(l=0;l<=sps.in(na,nb,nc)+1;++l){mq[l]=ComplexVector(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);}
                 Vector sq2(1,3*inputpars.cs.nofatoms),qs(1,3),qt(1,3);float in;qs(1)=1000;
                 sps.FT(mq); //Fourier trafo of spincf
		 // get the main propagation vector by looking for the
		 // biggest Fourier component of the magnetic moment arrangement 
                 for(qh=0;qh<sps.na();++qh){for(qk=0;qk<sps.nb();++qk){for(ql=0;ql<sps.nc();++ql)
                  {// get magnetic moment from momentum fouriercomponent into b 
		   b=0;na=sps.na()-qh,nb=sps.nb()-qk,nc=sps.nc()-ql;
		   b1 = mq[sps.in(na,nb,nc)];
                   for(l=1;l<=inputpars.cs.nofatoms;++l)
		   {int m1,m1max=3; if ((*inputpars.jjj[l]).gJ==0){m1max=6;}
		    for (m1=1;m1<=m1max;++m1)
		     {if((*inputpars.jjj[l]).gJ==0)
		      {if(m1==2||m1==4||m1==6){b(3*(l-1)+(m1+1)/2)+=b1(inputpars.cs.nofcomponents*(l-1)+m1);}
		       else                   {b(3*(l-1)+(m1+1)/2)+=2.0*b1(inputpars.cs.nofcomponents*(l-1)+m1);}
		      }
		      else
		      {b(3*(l-1)+m1)=b1(inputpars.cs.nofcomponents*(l-1)+m1)*(*inputpars.jjj[l]).gJ;
		      }
		     }    
		    }
		   a = b.Conjugate();
		   b1 = mq[sps.in(qh,qk,ql)];
		   b=0;
                   for(l=1;l<=inputpars.cs.nofatoms;++l)
		   {int m1,m1max=3; if ((*inputpars.jjj[l]).gJ==0){m1max=6;}
		    for (m1=1;m1<=m1max;++m1)
		     {if((*inputpars.jjj[l]).gJ==0)
		      {if(m1==2||m1==4||m1==6){b(3*(l-1)+(m1+1)/2)+=b1(inputpars.cs.nofcomponents*(l-1)+m1);}
		       else                   {b(3*(l-1)+(m1+1)/2)+=2.0*b1(inputpars.cs.nofcomponents*(l-1)+m1);}
		      }
		      else
		      {b(3*(l-1)+m1)=b1(inputpars.cs.nofcomponents*(l-1)+m1)*(*inputpars.jjj[l]).gJ;
		      }
		     }    
		    }                   
		   // inner product
                   sq2=Abs(b+a)/(double)sps.n()/(double)inputpars.cs.nofatoms;
                   Vector q(1,3);
		   q(1)=1.0*qh/sps.na();
	           q(2)=1.0*qk/sps.nb();
                   q(3)=1.0*ql/sps.nc();
                   qt=inputpars.rez.Transpose()*q;
		   in=Norm(sq2)*Norm(sq2);
	           if (in>inmax-0.001)
                    {if(in<inmax+0.001){++nk;}else{nk=1;}
                     inmax=in;qs=q;}
                   }}}

// inserted 26.11.2015 to get rec vector R such that R+q is smallest
     // try different Q vectors corresponding to q !!
    int i1,j1,k1;
    double QQmin=1e10,QQ;
// inserted 10.5.10 to make compatible with nonortholattices
     Matrix abc_in_ijk(1,3,1,3),p(1,3,1,3),pstar(1,3,1,3);
        get_abc_in_ijk(abc_in_ijk,inputpars.cs.abc);
     p=abc_in_ijk*inputpars.cs.r; // p is the primitive crystal unit cell in ijk coordinates
     pstar=2*PI*p.Inverse().Transpose();
     Vector nmin(1,3),nmax(1,3),hkl(1,3),hkls(1,3),Q(1,3),qeuklid(1,3);
     nlimits_calc(nmin, nmax, ini.maxQ, pstar);
     // problem: we want to find all lattice vectors Rn=ni*ai which are within a
     // sphere of radius r from the origin (ai = column vectors of matrix a)
     // this routine returns the maximum and minimum values of ni i=1,2,3
     // by probing the corners of a cube
              for (i1=(int)nmin(1);i1<=nmax(1);++i1){
              for (j1=(int)nmin(2);j1<=nmax(2);++j1){
              for (k1=(int)nmin(3);k1<=nmax(3);++k1){
       Q(1)=qs(1)+i1;Q(2)=qs(2)+j1;Q(3)=qs(3)+k1;
        //project back to big lattice
       hkl=inputpars.rez.Transpose()*Q;

      // qeuklid is Q in ijk coordinate system !
      hkl2ijk(qeuklid,hkl,inputpars.cs.abc);//qeuklid=ri;//qeuklid(1)=ri(1);qeuklid(2)=ri(2);qeuklid(3)=ri(3);
      QQ=Norm(qeuklid);if(QQ<QQmin){QQmin=QQ;hkls=hkl;}
 }}}
//------------------- end if insert 26.11.2015 -->> output is hkls with smallest |Q|

                   strcpy(outfilename,"./results/");strcpy(outfilename+10,ini.prefix);
                   strcpy(outfilename+10+strlen(ini.prefix),"mcphas.log");
                   felog=fopen_errchk(outfilename,"a");
                   if (verbose==1||fe>FEMIN_INI){fprintf(felog,"#! nofatoms=%i nofcomponents=%i fe=%10.6g\n#",inputpars.cs.nofatoms,inputpars.cs.nofcomponents,fe);}
      #ifndef _THREADS
                   fprintf(felog,"%10.6g %10.6g %10.6g %3i %10.6g %3i %3i %3i %3i %i %10.6g 1\n",hkls(1),hkls(2),hkls(3),nk,fe,j,sps.na(),sps.nb(),sps.nc(),r,sc);
      #else
                   fprintf(felog,"%10.6g %10.6g %10.6g %3i %10.6g %3i %3i %3i %3i %i %10.6g %i\n",hkls(1),hkls(2),hkls(3),nk,fe,j,sps.na(),sps.nb(),sps.nc(),r,sc,thread_id);
      #endif	     
                   if (verbose==1&&fe<2*FEMIN_INI){sps.print(felog);}
	           fclose(felog);
                  delete []mq;
                 }
      }
} // <--- this bracket is necessary to embrace ComplexMatrix  definitions
  // necessary in this routine and make them live within this bracket. at closing the bracket
  // destructor is called correctly. On Apple without this 
  // bracket the compiler will exit intcalc_Erefine (pthread_exit) without freeing memory ...
 
      #ifndef _THREADS
      return 1;
      #else
      MUTEX_LOCK(&mutex_loop);
      thrdat.thread_id = thread_id;
      EVENT_SIG(checkfinish);
      MUTEX_UNLOCK(&mutex_loop);
      #undef ini
      #undef physprops
      #undef inputpars
      #undef testqs
      #undef testspins
      #undef H
      #undef T
      #undef femin
      #if defined  (__linux__) || defined (__APPLE__)
      pthread_exit(NULL);
      #else
      return 0;
      #endif	     
      #endif // def _THREADS
}

int  htcalc (Vector Happ,double T,inipar & ini,par & inputpars,qvectors & testqs,
             testspincf & testspins, physproperties & physprops,int tracetest)
{/* calculates magnetic structure at a given HT- point  
  on input: 
    T	Temperature[K]
    Happ        Vector of Applied Magnetic Field [T] in ijk coordinates
                which if demag=1 corresponds to external applied field
                and if demag=0 corresponds to internal applied field (=field in the sample)
                which is desired for the calculation
    inputpars	Input parameters (exchange constants etc...)
    testqs	Set of propagation vectors to be tested 
    testspins	Set of Spinconfigurations to be tested
    tracetest   if nonzero indicate which structure of structure table should be calculated and traced
  on return:
    physprops	physical properties at (HT) point (i.e. magnetic structure
		neutron intensities, thermal expansion ...)	
 // returns 0 if successfull
 //  --> if no spinconfiguration has been found at ht point
 // returns 1 if recalculation of fe yields too different value
 // returns 2 if prevailing problem is maxnofmfloops reached
 // returns 3 if prevailing problem is maxspinchange is reached  
 */

 int i,j,k,is;
 int start_nofmaxspinchangeDIV=ini.nofmaxspinchangeDIV;
 int start_nofmaxloopDIV=ini.nofmaxloopDIV;
 Vector momentq0(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms),phi(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
 Vector nettom(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms),q(1,3);
 Vector h1(1,inputpars.cs.nofcomponents),hkl(1,3);
               
 double femin=FEMIN_INI;char text[MAXNOFCHARINLINE];char outfilename[MAXNOFCHARINLINE];
 spincf  sps(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents),sps1(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 spincf  spsmin(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 FILE * felog; // logfile for q dependence of fe
 FILE * fout;

if (T<=0.01){fprintf(stderr," ERROR htcalc - temperature too low - please check mcphas.ini !");exit(EXIT_FAILURE);}

 srand(time(0)); // initialize random number generator
 checkini(ini); // check if user pressed a button
 if (ini.logfevsQ==1) {strcpy(outfilename,"./results/");strcpy(outfilename+10,ini.prefix);
                       strcpy(outfilename+10+strlen(ini.prefix),"mcphas.log");
                       felog=fopen_errchk(outfilename,"a");
               fprintf(felog,"#Logging of h k l multiplicity fe[meV] spinconf_nr n1xn2xn3 nof_mf_loops spinchange threadid at T=%g Hi=%g Hj=%g Hk=%g\n",T,Happ(1),Happ(2),Happ(3));
               fclose(felog);
	      }
 if (verbose==1)
 { strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
   strcpy(outfilename+11+strlen(ini.prefix),"fe_status.dat");
   fout= fopen_errchk (outfilename,"w");
   #ifndef _THREADS
   fprintf(fout,"#displayxtext=time(s)\n");
   fprintf(fout,"#displaytitle=2:log(iterations) 3:log(sta)|E(meV) 4:log(spinchange)|U(meV) 5:stepratio 6:successrate 7:freeenergy(%%)\n");
   fprintf(fout,"#time(s) log(iteration) log(sta) log(spinchange+1e-10) stepratio  successrate=(nof stabilised structures)/(nof initial spinconfigs) freenergy(meV)\n");
   fprintf(fout,"#or for Monte Carlo option \n");
   fprintf(fout,"#time(s) log(iteration) E(meV)  U(meV)\n");
   fprintf(fout,"%i 0 0 0 0 0 0\n",(int)time(0));
   fprintf(fout,"%i 1 1 1 1 1 1\n",(int)time(0)+1);
   #else
   fprintf(fout,"#displayxtext=time(s)\n");
   fprintf(fout,"#displaytitle=2:log(iterations) 3:log(sta)|E(meV) 4:log(spinchange)|U(meV) 5:stepratio|threadid 6:successrate 7:freeenergy 8:threadID(%%)\n");
   fprintf(fout,"#time(s) log(iteration) log(sta) log(spinchange+1e-10) stepratio  successrate=(nof stabilised structures)/(nof initial spinconfigs) freenergy(meV)  thread_id \n");
   fprintf(fout,"#or for Monte Carlo option \n");
   fprintf(fout,"#time(s) log(iteration) E(meV)  U(meV) thread_id \n");
   fprintf(fout,"%i 0 0 0 0 0 0 0\n",(int)time(0));
   fprintf(fout,"%i 1 1 1 1 1 1 1\n",(int)time(0)+1);
   #endif
   fclose(fout);
   Vector P(1,3);Vector M(1,3);P=0;M=0;      
   printf("\n starting  "); ini.print_usrdefcols(stdout,physprops.x,physprops.y,T,Happ,inputpars.cs.abc,M,P,true);printf("\n");
   printf("with %i spinconfigurations read from mcphas.tst and table \nand\n %i spinconfigurations created from hkl's\n\n",testspins.n,testqs.nofqs());
   printf("Notation: < >            ...begin / end of mean field loop\n");
   printf("          ->             ...reduction of stabilised structure possible into ...\n");
   printf("          =              ...no reduction of stabilised structure found\n");
   printf("          [n]            ...thread number n\n");
   printf("          rn             ...random try n\n");
   printf("          fe(n)          ...free energy for random try n\n");
   printf("          (hkl)          ...Miller indizes (for abc unit cell)\n");
   printf("          (n1 x n2 x n3) ...supercell of primitive unit cell \n");   
   printf("          str s(r)       ...structure nr. s, initial values from number r\n"); 
 }

// Here we choose how to begin the scan of different initial configuarions
    // 1. randomly
    //j=-testqs.nofqs()+(int)rint(rnd(testspins.n+testqs.nofqs())); 
    //begin with j a random number, j<0 means test spinconfigurations 
    //constructed from q vector set testqs, j>0 means test spinconfigurations from
    //set testspins
    // 2. starting with the table loaded from mcphas.tst  into testspins
    j=tracetest-1;  //uncomment this for debugging purposes
    // 3. with the hkl - supercells generated from hmin hmax kmin kmax lmin lmax in mcphas.ini
    //j = -testqs.nofqs()-1;
#ifdef _THREADS
// ----------------------------------------------------------------------------------- //
// Populates the thread data structure
// ----------------------------------------------------------------------------------- //
   thrdat.Happ = Happ;
   thrdat.T = T;
   thrdat.ini=&ini;
   thrdat.testqs = &testqs; 
   thrdat.testspins = &testspins;
   thrdat.physprops = &physprops;
   thrdat.femin = femin;
   thrdat.spsmin = spsmin; 
   thrdat.thread_id = -1;
//   htcalc_input *tin[NUM_THREADS];
 /*  static int washere=0;
   if(washere==0){washere=1;if(NUM_THREADS>256){fprintf(stderr,"Error mcphas: too many threads required - change hardcode limit 256 in mcphas_htcalc.c line 69 and recompile\n");exit(EXIT_FAILURE);}
                  for (int ithread=0; ithread<NUM_THREADS; ithread++) 
                    tin[ithread] = new htcalc_input(0,ithread,&inputpars);
                  }
*/ // moved to mcphas.c
 MUTEX_INIT(mutex_loop);
 MUTEX_INIT(mutex_tests);
 MUTEX_INIT(mutex_min);
 MUTEX_INIT(mutex_ini_nofcalls);
 MUTEX_INIT(mutex_ini_nofmaxspinchangeDIV);
 MUTEX_INIT(mutex_ini_nofmaxloopDIV);
 MUTEX_INIT(mutex_ini_successrate);
 MUTEX_INIT(mutex_ini_calcmf_duration);
 MUTEX_INIT(mutex_ini_calcsps_duration);
 EVENT_INIT(checkfinish);
 THRLC_INIT(threadSpecificKey);
 #if defined  (__linux__) || defined (__APPLE__)
 pthread_t threads[NUM_THREADS]; int rc; void *status;
 pthread_attr_t attr;
 pthread_attr_init(&attr);
 pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
 #else
 HANDLE threads[NUM_THREADS];
 DWORD tid[NUM_THREADS], dwError;
 #endif
 bool all_threads_started = false; int ithread=0;
#endif
 for (k= -testqs.nofqs();(tracetest!=0 ? j<tracetest : k<=testspins.n );++k)
 {++j; if (j>testspins.n) j=-testqs.nofqs();
#ifndef _THREADS
       htcalc_iteration(j, femin, spsmin, H, T,ini, inputpars, testqs, testspins, physprops);
#else
        (*tin[ithread]).j = j;
       #if defined  (__linux__) || defined (__APPLE__)
       rc = pthread_create(&threads[ithread], &attr, htcalc_iteration, (void *) tin[ithread]);
       if(rc) 
       { // rc not zero - error means the system does not have the resources or permission to create thread
          if(rc) { printf("Warning return code %i when creating thread %i - trying again ...\n",rc,ithread+1);  }
          rc = pthread_join(threads[ithread], &status); 
          if(rc) { printf("Warning return code %i when joining thread %i\n",rc,ithread+1);  }
          rc = pthread_create(&threads[ithread], &attr, htcalc_iteration, (void *) tin[ithread]);
          if(rc) { printf("Error return code %i from creating thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
       }
       #else
       threads[ithread] = CreateThread(NULL, 0, htcalc_iteration, (void *) tin[ithread], 0, &tid[ithread]);
       if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %lu from thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
       #endif
        ithread++;
       if(ithread%NUM_THREADS==0 || all_threads_started)
       {  all_threads_started = true;
          #if defined  (__linux__) || defined (__APPLE__)
          pthread_mutex_lock (&mutex_loop); 
          while(thrdat.thread_id==-1) pthread_cond_wait(&checkfinish, &mutex_loop);
          ithread = thrdat.thread_id;
          thrdat.thread_id=-1; 
          pthread_mutex_unlock (&mutex_loop); 
          rc = pthread_join(threads[ithread], &status); 
          if(rc) { printf("Error return code %i from joining thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
     
          #else
          WaitForSingleObject(checkfinish,INFINITE);
          ithread = thrdat.thread_id;
          CloseHandle(threads[ithread]);
          thrdat.thread_id=-1; 
          ResetEvent(checkfinish);
          #endif
       }
#endif
    }
#ifdef _THREADS
// Wait for all threads to finish, before moving on to calculate physical properties!
  for(int th=0; th<(all_threads_started?NUM_THREADS:ithread); th++)
  if(th!=ithread){ // thread ithread does not exist or has already been closed above    
     #if defined  (__linux__) || defined (__APPLE__)
     rc = pthread_join(threads[th], &status); 
     if(rc) { printf("Error return code %i from joining thread %i\n",rc,th+1); exit(EXIT_FAILURE); }
     #else
     if(WaitForSingleObject(threads[th],INFINITE)==0xFFFFFFFF) { printf("Error in waiting for thread %i to end\n",th+1); exit(EXIT_FAILURE); }
     CloseHandle(threads[th]);
     #endif
                
                }
  femin = thrdat.femin;

 #if defined  (__linux__) || defined (__APPLE__)
 pthread_attr_destroy(&attr);
 pthread_mutex_destroy(&mutex_loop);
 pthread_mutex_destroy(&mutex_tests);
 pthread_mutex_destroy(&mutex_min);
 pthread_mutex_destroy(&mutex_ini_nofcalls);
 pthread_mutex_destroy(&mutex_ini_nofmaxspinchangeDIV);
 pthread_mutex_destroy(&mutex_ini_nofmaxloopDIV);
 pthread_mutex_destroy(&mutex_ini_successrate);
 pthread_mutex_destroy(&mutex_ini_calcmf_duration);
 pthread_mutex_destroy(&mutex_ini_calcsps_duration);

 #endif
 THRLC_FREE(threadSpecificKey);
#endif

if (femin>=FEMIN_INI) // did we find a stable structure ??
 {if(ini.nofmaxspinchangeDIV-start_nofmaxspinchangeDIV<ini.nofmaxloopDIV-start_nofmaxloopDIV)
          return 2; else return 3;
 }
else // if yes ... then
 {if(verbose==1){printf("... calculating physical properties ");}
 if (physprops.j>0){ // take spinconfiguration ----
                     sps=(*testspins.configurations[physprops.j]);
                       if (sps.wasstable==0)
                       {// go through qvectors and spinfconfigurations and see if periodicity matches
                        for (i=1;i<=testqs.nofqs();++i)
                         {if (testqs.na(i)==sps.na()&&testqs.nb(i)==sps.nb()&&testqs.nc(i)==sps.nc())
                             {sps.wasstable=-i;break;}
                         }
		        if (sps.wasstable==0)
                         {for (i=1;i<=testspins.n;++i)
                          {if ((*testspins.configurations[i]).na()==sps.na()&&
			       (*testspins.configurations[i]).nb()==sps.nb()&&
			       (*testspins.configurations[i]).nc()==sps.nc())
                             {sps.wasstable=i;break;}
                          }
                         }
			if (sps.wasstable==0){fprintf(stderr,"internal ERROR htcalc - calculating periodicity not possible");exit(EXIT_FAILURE);}
			//---mark it as stable with periodicity key---
			(*testspins.configurations[physprops.j]).wasstable=sps.wasstable;    
                       }
	      }
 if(verbose==1){printf("<");}
/*    else     // ---- or take q vector 
            // removed because not necessary MR 15.12.15
            { sps.spinfromq(testqs.na(-physprops.j),testqs.nb(-physprops.j),
	              testqs.nc(-physprops.j),testqs.q(-physprops.j),
		      testqs.nettom(-physprops.j),testqs.momentq0(-physprops.j),
		      testqs.phi(-physprops.j));
	      }
*/
     #ifndef _THREADS
     sps=spsmin;//take spinconfiguration which gave minimum free energy as starting value
     #else
     sps=thrdat.spsmin;//take spinconfiguration which gave minimum free energy as starting value
     #endif
   //MR 120221 removed spinconf invert in case nettoI is negative
  // now really calculate the physical properties
      mfcf mf(sps.na(),sps.nb(),sps.nc(),inputpars.cs.nofatoms,inputpars.cs.nofcomponents);int r;double sc;
      physprops.fe=fecalc(physprops.u,physprops.Eelastic,r,sc,Happ ,T,ini,inputpars,sps,mf,&physprops); 

      spincf magmom(sps.na(),sps.nb(),sps.nc(),inputpars.cs.nofatoms,3);
                   int i1,j1,k1,l1,m1;Vector mom(1,3),d1(1,inputpars.cs.nofcomponents);
                   for (l1=1;l1<=inputpars.cs.nofatoms;++l1){
                    // go through magnetic unit cell and sum up the contribution of every atom
                  for(i1=1;i1<=sps.na();++i1){for(j1=1;j1<=sps.nb();++j1){for(k1=1;k1<=sps.nc();++k1){
                  for(m1=1;m1<=inputpars.cs.nofcomponents;++m1){d1[m1]=mf.mf(i1,j1,k1)[inputpars.cs.nofcomponents*(l1-1)+m1];}                  
                   (*inputpars.jjj[l1]).mcalc(mom,T,d1,Happ,(*inputpars.jjj[l1]).Icalc_parstorage);
                    for(m1=1;m1<=3;++m1){magmom.m(i1,j1,k1)(3*(l1-1)+m1)=mom(m1);}
                    }}}}
             // display spinstructure
                if (verbose==1)
                {//printf("");
		 float * x;x=new float[inputpars.cs.nofatoms+1];float *y;y=new float[inputpars.cs.nofatoms+1];float*z;z=new float[inputpars.cs.nofatoms+1];
		 for (is=1;is<=inputpars.cs.nofatoms;++is)
		   {x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
                     snprintf(text,MAXNOFCHARINLINE,"recalculated: fe=%g,femin=%g:T=%gK,|H|=%gT,Ha=%gT, Hb=%gT, Hc=%gT, %i spins",physprops.fe,femin,T,Norm(Happ),Happ(1),Happ(2),Happ(3),sps.n());
                    strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dab.eps");
                     fout = fopen_errchk (outfilename, "w");
                     sps.eps3d(fout,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,4,magmom);
                    fclose (fout);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dac.eps");
                    fout = fopen_errchk (outfilename, "w");
                     sps.eps3d(fout,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,5,magmom);
                    fclose (fout);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dbc.eps");
                    fout = fopen_errchk (outfilename, "w");
                     sps.eps3d(fout,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,6,magmom);
                    fclose (fout);
		    strcpy(outfilename+11+strlen(ini.prefix),"spins.eps");
                    fout = fopen_errchk (outfilename, "w");
                     magmom.eps(fout,text);
                    fclose (fout);
                delete[]x;delete []y; delete []z;
		}
  //printf("G");fflush(stdout);delete magmom;
   if(verbose==1){printf(">");}
 //check if fecalculation gives again correct result
   if (physprops.fe>femin+(0.00001*fabs(femin))&&ini.maxnofmfloops>2&&ini.nofMCsteps==0){int eq=0;
   #ifndef _THREADS
   if(spsmin==sps){eq=1;};//take spinconfiguration which gave minimum free energy as starting value
     #else
   if(thrdat.spsmin==sps){eq=1;};//take spinconfiguration which gave minimum free energy as starting value
     #endif
   
   if(verbose){fprintf(stderr,"Warning htcalc.c: at T=%g K /  H= %g Tfemin=%4.9g was calc.(conf no %i),\n but recalculation  gives fe= %4.9gmeV -> no structure saved\n",
                            T,Norm(Happ),femin,physprops.j,physprops.fe);
   fprintf(stderr,"recalculation converged after %i loops and initial and final spin structures are ",r);
   if(eq==1){fprintf(stderr,"equal\n");}else{fprintf(stderr,"not equal\n");}}
if (ini.logfevsQ==1) {strcpy(outfilename,"./results/");strcpy(outfilename+10,ini.prefix);
                       strcpy(outfilename+10+strlen(ini.prefix),"mcphas.log");
                       felog=fopen_errchk(outfilename,"a");
               fprintf(felog,"#Warning htcalc.c: at T=%g K /  H= %g Tfemin=%4.9g was calc.(conf no %i),\n# but recalculation  gives fe= %4.9gmeV -> no structure saved\n",
                T,Norm(Happ),femin,physprops.j,physprops.fe);fprintf(felog,"#recalculation converged after %i loops and initial and final spin structures are ",r);
   if(eq==1){fprintf(felog,"equal\n");}else{fprintf(stderr,"not equal\n#initial values as converged from femin=%4.9g meV calculation:\n",femin);
#ifndef _THREADS
     spsmin.print(felog);//take spinconfiguration which gave minimum free energy as starting value
     #else
     thrdat.spsmin.print(felog);//take spinconfiguration which gave minimum free energy as starting value
     #endif
     }
               fclose(felog);
	      }
                             physprops.sps.epsilon=0;physprops.Eelastic=0;
                             physprops.m=0;
//                printf("I\n");fflush(stdout);delete mf;
                              return 1;
                             }
 //if(verbose==1){printf(".\n");}
if(ini.nofMCsteps==0) physpropclc(Happ,T,sps,mf,physprops,ini,inputpars);
else physprops.sps=sps;

//    printf("H"); fflush(stdout); delete mf;
 }
 


return 0; // ok we are done with this (HT) point- return ok
// #if defined __linux__ && defined _THREADS
// pthread_exit(NULL);
// #endif
}





/*****************************************************************************/
// this sub checks if a spinconfiguration has already been added to
// table testspins and adds it if necessary
int checkspincf(int & j,spincf & sps1,qvectors & testqs,Vector & nettom,
		     Vector & momentq0, Vector & phi, 
                     testspincf & testspins,physproperties & physprops,inipar & ini)
{ int i;
  spincf sps(1,1,1,sps1.nofatoms,sps1.nofcomponents);
  sps=sps1;sps.reduce();// reduce inserted MR 20120907

// compare spinconfigurations stabilized by 
// index j with existing spinconfigurations in testspins
  spincf spq(1,1,1,sps.nofatoms,sps.nofcomponents);

// compare new configuration to all stored configurations 
//check all spinconfigurations

 for (i=testspins.ninitial;i>=-testqs.nofqs();--i)
 {
  if (i>0) 
   {if (sps==(*testspins.configurations[i])) 
	 {
	 physprops.j=i;return 1;} //ok
   }
  if (i<0)
  { int mi=-i;int ta=testqs.na(mi),tb=testqs.nb(mi),tc=testqs.nc(mi);
    spq.spinfromq(ta,tb,tc,testqs.q(mi),
                   testqs.nettom(mi),testqs.momentq0(mi),testqs.phi(mi));
    if (spq==sps){
    physprops.j=i;return 1;} //ok
   
  }
 } 

   //  check initial config: take just used nettom,momentq0,phi for comparison
   if (j<0)
   {int mj=-j;int ta=testqs.na(mj),tb=testqs.nb(mj),tc=testqs.nc(mj);
    spq.spinfromq(ta,tb,tc,testqs.q(mj),
                  nettom,momentq0,phi);
    if (spq==sps) {physprops.j=j;testqs.nettom(mj)=nettom;
                  testqs.momentq0(mj)=momentq0;testqs.phi(mj)=phi;return 1;} //ok
   } 

// check newly added configuration
for (i=testspins.ninitial+1;i<=testspins.n;++i)
 {if (sps==(*testspins.configurations[i])) 
	 {
	 physprops.j=i;return 1;} //ok
   }
// if it gets here, the spins sps configuration has not been found
// -. add configuration to testspins 
//- first make sure that wasstable is 0 [might be nonzero from initialisation](MR 13.1.2015)
sps.wasstable=0;
return (physprops.j=testspins.addspincf(sps));  //ok=1
}




