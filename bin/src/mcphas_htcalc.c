// routines for mcphas for calculation of magnetic phases
// htcalc.c

#define FEMIN_INI     1e6

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
   Vector H;
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
   htcalc_input(int _j, int _tid, par *pars_in) 
   { 
      thread_id = _tid; j = _j; inputpars = new par(*pars_in);
   }
   ~htcalc_input(){delete inputpars;}
};
// ----------------------------------------------------------------------------------- //
// Declares these variables global, so all threads can see them
// ----------------------------------------------------------------------------------- //
htcalc_thread_data thrdat;
htcalc_input *tin[256];  // Max number of threads - hard coded because global variable.
MUTEX_TYPE mutex_loop;
MUTEX_TYPE mutex_tests;
MUTEX_TYPE mutex_min;
EVENT_TYPE checkfinish;
THRLC_TYPE threadSpecificKey;

#endif // def _THREADS

void checkini(testspincf & testspins,qvectors & testqs,inipar & ini)
{struct stat filestatus;
 static time_t last_modify_time;
 static int washere=0;
 int loaderr;
  errno = 0;

  if (stat(ini.savfilename,&filestatus)!=0)
    {fprintf (stderr, "Error checking mcphas.ini: Couldn't read status of file %s: %s\n",
              ini.savfilename, strerror (errno));exit (EXIT_FAILURE);
     }

  if(washere==0){washere=1;last_modify_time=filestatus.st_mtime;}
  
   
    if (filestatus.st_mtime!=last_modify_time) //check if file has been modified
    {again:
     last_modify_time=filestatus.st_mtime;
     fprintf(stdout,"mcphas.ini has been modified - reading new mcphas.ini\n");
      sleep(1);
      loaderr=ini.load();
      if(ini.exit_mcphas==1)
        {testspins.save(filemode);  //exit normally
         testqs.save(filemode);
      printf("RESULTS saved in directory ./results/  - files:\n");
   printf("  mcphas.fum  - total magnetic moment, energy at different T,H\n");
   printf("  mcphas.sps  - stable configurations at different T,H\n");
   printf("  mcphas.mf   - mean fields at different T,H\n");
   printf("  mcphas.hkl  - strong magnetic satellites, neutron diffraction intensity\n");
   printf("  mcphas*.hkl - strong magnetic satellites, Fourier Comp.of moment in * dir\n");
   printf("  mcphas*.j*  - JJ correlation functions (for exchange magnetostriction)\n");
   printf("  mcphas.xyt  - phasediagram (stable conf.nr, angular and multipolar moments)\n");
   printf("  mcphas.qvc  - ...corresponding table of all qvector generated test configs\n");
   printf("  mcphas.phs  - ...corresponding table of all test configurations (except qvecs)\n");
   printf("  _mcphas.*   - parameters read from input parameter files (.tst,.ini,.j)\n");
   printf("  ...         - and a copy of the single ion parameter files used.\n\n");
   fprintf(stderr,"**********************************************\n");
   fprintf(stderr,"          End of Program mcphas\n");
   fprintf(stderr," reference: M. Rotter JMMM 272-276 (2004) 481\n");
   fprintf(stderr,"**********************************************\n");
          exit(0);
	}

      while(ini.pause_mcphas==1||loaderr==1) // wait until pause button is released and no loaderror occurs
       {fprintf(stdout,"Pausing ...\n");
        while(filestatus.st_mtime==last_modify_time)  //wait until filestatus changes again
          {sleep(1);
            if (stat(ini.savfilename,&filestatus)!=0)
               {fprintf (stderr, "Error checking file mcphas.ini: Couldn't read status of file %s: %s\n",
                ini.savfilename, strerror (errno));exit (EXIT_FAILURE);
               }
          }
	goto again;  
       }
       

    }
}

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
int htcalc_iteration(int j, double &femin, spincf &spsmin, Vector H, double T,inipar & ini, par &inputpars, qvectors &testqs, testspincf &testspins, physproperties &physprops)
#endif
{
 fflush(stderr); fflush(stdout);
 #ifdef _THREADS
 htcalc_input *myinput; myinput = (htcalc_input *) input; int j = myinput->j, thread_id = myinput->thread_id; Vector H(1,3); H = thrdat.H;
 THRLC_SET(threadSpecificKey, myinput); int tlsfemin=0;  // Thread local variable to judge whether to print output
 #else
 int thread_id=1;
 #endif 
 int i,ii,iii,tryrandom,nr,rr,ri,is,r;
 double fe,fered,Eel,U,sc;
 double u,lnz; // free- and magnetic energy per ion [meV]
 Vector momentq0(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms),phi(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
 Vector nettom(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms),q(1,3);
 Vector mmom(1,inputpars.cs.nofcomponents);
 Vector h1(1,inputpars.cs.nofcomponents),h1ext(1,3),hkl(1,3);
 h1ext=0;
 char text[MAXNOFCHARINLINE];
 char outfilename[MAXNOFCHARINLINE];
 spincf  sps(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents),sps1(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 mfcf * mf;
 mfcf * mf1;
 spincf * magmom;
 FILE * felog; // logfile for q dependence of fe
 FILE * fin_coq;

  for (tryrandom=0;tryrandom<=ini.nofrndtries&&j!=0;++tryrandom)
   {if (j>0){sps=(*testspins.configurations[j]);// take test-spinconfiguration
             #ifndef _THREADS
	     if (tryrandom==0&&verbose==1) { printf ( "str%i(%ix%ix%i)< "  ,j,sps.na(),sps.nb(),sps.nc()); fflush(stdout); }
             #else
	     if (tryrandom==0&&verbose==1) { printf ( "str%i(%ix%ix%i)<[%i] "  ,j,sps.na(),sps.nb(),sps.nc(),thread_id); fflush(stdout); }
             #endif 
            }
    else     // take q vector and choose phase and mom dir randomly
            {q=testqs.q(-j);  
	     if (tryrandom==0)
	     {nettom=testqs.nettom(-j);momentq0=testqs.momentq0(-j);phi=testqs.phi(-j);
	     }
	     else
	     {for(i=1;i<=inputpars.cs.nofatoms;++i)
	      {for(ii=1;ii<=inputpars.cs.nofcomponents;++ii)
	        {iii=inputpars.cs.nofcomponents*(i-1)+ii;h1=0;h1(ii)=10*MU_B;
                 (*inputpars.jjj[i]).Icalc(mmom,T,h1,h1ext,lnz,u,(*inputpars.jjj[i]).Icalc_parstorage);
		 nettom(iii)=mmom(ii)*rnd(1);
	         momentq0(iii)=rnd(1)*mmom(ii);
	         phi(iii)=rnd(1)*3.1415;
		}
	      }
	     }
	     sps.spinfromq(testqs.na(-j),testqs.nb(-j),testqs.nc(-j),
	                   q,nettom,momentq0,phi);
             hkl=inputpars.rez.Transpose()*q;  
             #ifndef _THREADS
   	     if (tryrandom==0&&verbose==1) { printf ( "(%g %g %g)(%ix%ix%i)< ",hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc()); fflush(stdout); }
             #else
   	     if (tryrandom==0&&verbose==1) { printf ( "(%g %g %g)(%ix%ix%i)<[%i] ",hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc(),thread_id); fflush(stdout); }
             #endif 
	    }	 
    if (tryrandom>0){nr=(int)(rint(rnd(1.0)*(sps.n()*inputpars.cs.nofatoms-1)))+1;
	             for (i=1;i<=nr;++i) //MonteCarlo randomize nr spins
                      {rr=(int)rint(rnd(1.0)*(sps.n()-1))+1;
		       ri=inputpars.cs.nofcomponents*(int)rint(rnd(1.0)*(inputpars.cs.nofatoms-1));
	               for(ii=1;ii<=inputpars.cs.nofcomponents;++ii)
		       {sps.mi(rr)(ri+ii)*=(2*rnd(1.0)-1) ;}
		       } // randomize spin rr
                    }
 
      //!!!calculate free energy - this is the heart of this loop !!!!
      mf=new mfcf(sps.na(),sps.nb(),sps.nc(),inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
      fe=fecalc(U,Eel,r,sc,H ,T,ini,inputpars,sps,(*mf),testspins,testqs);
          if (fe>=2*FEMIN_INI && verbose==1) {
	       if(j>0) printf ( ">for_str_%i(%ix%ix%i) "  ,j,sps.na(),sps.nb(),sps.nc());
               else    printf ( ">for(%g %g %g)(%ix%ix%i) ",hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc()); 
                                             }
          
      // test spinconfiguration  and remember it                                    
      if (fe<femin)
            {               // first - reduce the spinconfiguration if possible
               sps1=sps;if(1==sps1.reduce()){ // if reduction is successful, try if the energy is less or equal for reduced spoinconfigurations
                   mf1=new mfcf(sps1.na(),sps1.nb(),sps1.nc(),inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
               if ((fered=fecalc(U,Eel,r,sc,H ,T,ini,inputpars,sps1,(*mf1),testspins,testqs))<=fe*(1.0000000000001)){(*mf)=(*mf1);
                                 if (verbose==1){fprintf(stdout,">[%i](%ix%ix%i)r%i->(%ix%ix%i)fe=%f->%fmeV ",thread_id,sps.na(),sps.nb(),sps.nc(),tryrandom,sps1.na(),sps1.nb(),sps1.nc(),fe,fered); fflush(stdout);}
                                                                                     sps=sps1;fe=fered;}
                                                                                                  else {
                                 if (verbose==1){fprintf(stdout,">[%i](%ix%ix%i)r%ife=%.15gmeV<(%ix%ix%i)fered=%.15g ",thread_id,sps.na(),sps.nb(),sps.nc(),tryrandom,fe,sps1.na(),sps1.nb(),sps1.nc(),fered);fflush(stdout);}
                                                                                                       }
                               delete mf1;  }
                                       else {
                                 if (verbose==1){fprintf(stdout,">[%i](%ix%ix%i)r%ife=%fmeV ",thread_id,sps.na(),sps.nb(),sps.nc(),tryrandom,fe);fflush(stdout);}
                                            }
                   magmom=new spincf(sps.na(),sps.nb(),sps.nc(),inputpars.cs.nofatoms,3);
                   int i1,j1,k1,l1,m1;Vector mom(1,3),d1(1,inputpars.cs.nofcomponents);
                   for (l1=1;l1<=inputpars.cs.nofatoms;++l1){
                    // go through magnetic unit cell and sum up the contribution of every atom
                  for(i1=1;i1<=sps.na();++i1){for(j1=1;j1<=sps.nb();++j1){for(k1=1;k1<=sps.nc();++k1){
                   for(m1=1;m1<=inputpars.cs.nofcomponents;++m1){d1[m1]=(*mf).mf(i1,j1,k1)[inputpars.cs.nofcomponents*(l1-1)+m1];}
                   (*inputpars.jjj[l1]).mcalc(mom,T,d1,H,(*inputpars.jjj[l1]).Icalc_parstorage);
                   for(m1=1;m1<=3;++m1){(*magmom).m(i1,j1,k1)(3*(l1-1)+m1)=mom(m1);}
                    }}}} 
                  
                 // display spinstructure
                if (verbose==1)
                {float * x;x=new float[inputpars.cs.nofatoms+1];float *y;y=new float[inputpars.cs.nofatoms+1];float*z;z=new float[inputpars.cs.nofatoms+1];
		 
		 for (is=1;is<=inputpars.cs.nofatoms;++is)
		   {
                    x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
                     snprintf(text,MAXNOFCHARINLINE,"fe=%g<femin=%g:T=%gK, |H|=%gT,Ha=%gT, Hb=%gT, Hc=%gT,  %i spins",fe,femin,T,Norm(H),H(1),H(2),H(3),sps.n());
                    strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dab.eps");
                    fin_coq = fopen_errchk (outfilename, "w");
                     sps.eps3d(fin_coq,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,4,(*magmom));
                    fclose (fin_coq);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dac.eps");
                    fin_coq = fopen_errchk (outfilename, "w");
                     sps.eps3d(fin_coq,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,5,(*magmom));
                    fclose (fin_coq);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dbc.eps");
                    fin_coq = fopen_errchk (outfilename, "w");
                     sps.eps3d(fin_coq,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,6,(*magmom));
                    fclose (fin_coq);
		   
                    strcpy(outfilename+11+strlen(ini.prefix),"spins.eps");
                    fin_coq = fopen_errchk (outfilename, "w");
                     (*magmom).eps(fin_coq,text);
                    fclose (fin_coq);
		delete[]x;delete []y; delete []z;
	        }
                delete magmom;   
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
	     }
            delete mf;
            //printout fe
            #ifdef _THREADS
	    if (tryrandom==ini.nofrndtries && verbose==1) {
               if(tlsfemin) printf("[%i]femin=%gmeV str %i(%i)-",thread_id,fe,physprops.j,j); 
	       if(j>0) printf ( ">[%i]str %i(%ix%ix%i)done "  ,thread_id,j,sps.na(),sps.nb(),sps.nc());
               else    printf ( ">[%i](%g %g %g)(%ix%ix%i)done ",thread_id,hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc()); 
                                                          }
            #endif
            if (tryrandom==ini.nofrndtries&&verbose==1){printf("\n");}
 
	    
  // log fe if required
   if (ini.logfevsQ==1) {
                 ComplexVector a(1,3*inputpars.cs.nofatoms),b(1,3*inputpars.cs.nofatoms);
                 ComplexVector b1(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
                 float inmax=0;int qh,qk,ql,l,nk=0;
                 ComplexVector * mq;  
                 mq = new ComplexVector [sps.in(sps.na(),sps.nb(),sps.nc())+2];for(l=0;l<=sps.in(sps.na(),sps.nb(),sps.nc())+1;++l){mq[l]=ComplexVector(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);}
                 Vector sq2(1,3*inputpars.cs.nofatoms),qs(1,3),qt(1,3);float in;qs(1)=1000;
                 sps.FT(mq); //Fourier trafo of spincf
		 // get the main propagation vector by looking for the
		 // biggest Fourier component of the magnetic moment arrangement 
                 for(qh=0;qh<sps.na();++qh){for(qk=0;qk<sps.nb();++qk){for(ql=0;ql<sps.nc();++ql)
                  {// get magnetic moment from momentum fouriercomponent into b 
		   b=0;
		   b1 = mq[sps.in(sps.na()-qh,sps.nb()-qk,sps.nc()-ql)];
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

int  htcalc (Vector Habc,double T,inipar & ini,par & inputpars,qvectors & testqs,
             testspincf & testspins, physproperties & physprops)
{/* calculates magnetic structure at a given HT- point  
  on input: 
    T	Temperature[K]
    Habc	Vector of External Magnetic Field [T] (components along crystal axes abc)
    inputpars	Input parameters (exchange constants etc...)
    testqs	Set of propagation vectors to be tested 
    testspins	Set of Spinconfigurations to be tested
  on return:
    physprops	physical properties at (HT) point (i.e. magnetic structure
		neutron intensities, thermal expansion ...)	
 // returns 0 if successfull
  // returns 2 if no spinconfiguration has been found at ht point
 */
 int i,j,k,is;
 Vector momentq0(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms),phi(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
 Vector nettom(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms),q(1,3);
 Vector h1(1,inputpars.cs.nofcomponents),hkl(1,3);
 Vector H(1,3); // magnetic field in ijk coordinate system
 Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1; // trick to get Habc as components along a,b,c
                  abc(4)=inputpars.cs.alpha(); abc(5)=inputpars.cs.beta(); abc(6)=inputpars.cs.gamma();
 dadbdc2ijk(H,Habc,abc); // transform Habc to ijk coordinates ... this is H
                
 double femin=FEMIN_INI;char text[MAXNOFCHARINLINE];char outfilename[MAXNOFCHARINLINE];
 spincf  sps(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents),sps1(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 spincf  spsmin(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 mfcf * mf;
 spincf * magmom;
 FILE * felog; // logfile for q dependence of fe
 FILE * fin_coq;

if (T<=0.01){fprintf(stderr," ERROR htcalc - temperature too low - please check mcphas.ini !");exit(EXIT_FAILURE);}

 srand(time(0)); // initialize random number generator
 checkini(testspins,testqs,ini); // check if user pressed a button
 if (ini.logfevsQ==1) {strcpy(outfilename,"./results/");strcpy(outfilename+10,ini.prefix);
                       strcpy(outfilename+10+strlen(ini.prefix),"mcphas.log");
                       felog=fopen_errchk(outfilename,"a");
               fprintf(felog,"#Logging of h k l multiplicity fe[meV] spinconf_nr n1xn2xn3 nof_mf_loops spinchange threadid at T=%g Ha=%g Hb=%g Hc=%g\n",T,Habc(1),Habc(2),Habc(3));
               fclose(felog);
	      }
 if (verbose==1)
 { strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
   strcpy(outfilename+11+strlen(ini.prefix),"fe_status.dat");
   fin_coq= fopen_errchk (outfilename,"w");
   #ifndef _THREADS
   fprintf(fin_coq,"#displayxtext=time(s)\n");
   fprintf(fin_coq,"#displaytitle=2:log(iterations) 3:log(sta) 4:log(spinchange) 5:stepratio 6:successrate 7:freeenergy(%%)\n");
   fprintf(fin_coq,"#time(s) log(iteration) log(sta) log(spinchange+1e-10) stepratio  successrate=(nof stabilised structures)/(nof initial spinconfigs) freenergy(meV)\n");
   fprintf(fin_coq,"%i 0 0 0 0 0 0\n",(int)time(0));
   fprintf(fin_coq,"%i 1 1 1 1 1 1\n",(int)time(0)+1);
   #else
   fprintf(fin_coq,"#displayxtext=time(s)\n");
   fprintf(fin_coq,"#displaytitle=2:log(iterations) 3:log(sta) 4:log(spinchange) 5:stepratio 6:successrate 7:freeenergy 8:threadID(%%)\n");
   fprintf(fin_coq,"#time(s) log(iteration) log(sta) log(spinchange+1e-10) stepratio  successrate=(nof stabilised structures)/(nof initial spinconfigs) freenergy(meV)  thread_id \n");
   fprintf(fin_coq,"%i 0 0 0 0 0 0 0\n",(int)time(0));
   fprintf(fin_coq,"%i 1 1 1 1 1 1 1\n",(int)time(0)+1);
   #endif
   fclose(fin_coq);	      
   printf("\n starting T=%g H=%g Ha=%g Hb=%g Hc=%g ",T,Norm(H),Habc(1),Habc(2),Habc(3));
   if (inputpars.cs.alpha()!=90||inputpars.cs.beta()!=90||inputpars.cs.gamma()!=90){printf("Hi=%g Hj=%g Hk=%g",H(1),H(2),H(3));}
   printf("\n");
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
    j=0;  //uncomment this for debugging purposes
    // 3. with the hkl - supercells generated from hmin hmax kmin kmax lmin lmax in mcphas.ini
    //j = -testqs.nofqs()-1;

#ifdef _THREADS
// ----------------------------------------------------------------------------------- //
// Populates the thread data structure
// ----------------------------------------------------------------------------------- //
   thrdat.H = H;
   thrdat.T = T;
   thrdat.ini=&ini;
   thrdat.testqs = &testqs; 
   thrdat.testspins = &testspins;
   thrdat.physprops = &physprops;
   thrdat.femin = femin;
   thrdat.spsmin = spsmin; 
   thrdat.thread_id = -1;
//   htcalc_input *tin[NUM_THREADS];
   static int washere=0;

   if(washere==0){washere=1;
                  for (int ithread=0; ithread<NUM_THREADS; ithread++) 
                    tin[ithread] = new htcalc_input(0,ithread,&inputpars);
                  }
   
 MUTEX_INIT(mutex_loop);
 MUTEX_INIT(mutex_tests);
 MUTEX_INIT(mutex_min);
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
 for (k= -testqs.nofqs();k<=testspins.n;++k)
 {++j; if (j>testspins.n) j=-testqs.nofqs();
#ifndef _THREADS
       htcalc_iteration(j, femin, spsmin, H, T,ini, inputpars, testqs, testspins, physprops);
#else
        (*tin[ithread]).j = j;
       #if defined  (__linux__) || defined (__APPLE__)
       rc = pthread_create(&threads[ithread], &attr, htcalc_iteration, (void *) tin[ithread]);
       if(rc) 
       {
          pthread_join(threads[ithread], &status); 
          rc = pthread_create(&threads[ithread], &attr, htcalc_iteration, (void *) tin[ithread]);
          if(rc) { printf("Error return code %i from thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
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

          #else
          WaitForSingleObject(checkfinish,INFINITE);
          ithread = thrdat.thread_id;
          thrdat.thread_id=-1; 
          ResetEvent(checkfinish);
          #endif
       }
#endif
    }
#ifdef _THREADS
// Wait for all threads to finish, before moving on to calculate physical properties!
  for(int th=0; th<(all_threads_started?NUM_THREADS:ithread); th++)
  {
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
 #endif
 THRLC_FREE(threadSpecificKey);
#endif

if (femin>=FEMIN_INI) // did we find a stable structure ??
 {fprintf(stderr,"Warning propcalc: femin positive ... no stable structure found at  T= %g K / Ha= %g Hb= %g Hc= %g  T\n",
                 physprops.T,physprops.H(1),physprops.H(2),physprops.H(3));return 2;}
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
      mf=new mfcf(sps.na(),sps.nb(),sps.nc(),inputpars.cs.nofatoms,inputpars.cs.nofcomponents);int r;double sc;
      physprops.fe=fecalc(physprops.u,physprops.Eel,r,sc,H ,T,ini,inputpars,sps,(*mf),testspins,testqs); 

      magmom=new spincf(sps.na(),sps.nb(),sps.nc(),inputpars.cs.nofatoms,3);
                   int i1,j1,k1,l1,m1;Vector mom(1,3),d1(1,inputpars.cs.nofcomponents);
                   for (l1=1;l1<=inputpars.cs.nofatoms;++l1){
                    // go through magnetic unit cell and sum up the contribution of every atom
                  for(i1=1;i1<=sps.na();++i1){for(j1=1;j1<=sps.nb();++j1){for(k1=1;k1<=sps.nc();++k1){
                  for(m1=1;m1<=inputpars.cs.nofcomponents;++m1){d1[m1]=(*mf).mf(i1,j1,k1)[inputpars.cs.nofcomponents*(l1-1)+m1];}                  
                   (*inputpars.jjj[l1]).mcalc(mom,T,d1,H,(*inputpars.jjj[l1]).Icalc_parstorage);
                    for(m1=1;m1<=3;++m1){(*magmom).m(i1,j1,k1)(3*(l1-1)+m1)=mom(m1);}
                    }}}}
             // display spinstructure
                if (verbose==1)
                {//printf("");
		 float * x;x=new float[inputpars.cs.nofatoms+1];float *y;y=new float[inputpars.cs.nofatoms+1];float*z;z=new float[inputpars.cs.nofatoms+1];
		 for (is=1;is<=inputpars.cs.nofatoms;++is)
		   {x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
                     snprintf(text,MAXNOFCHARINLINE,"recalculated: fe=%g,femin=%g:T=%gK,|H|=%gT,Ha=%gT, Hb=%gT, Hc=%gT, %i spins",physprops.fe,femin,T,Norm(H),physprops.H(1),physprops.H(2),physprops.H(3),sps.n());
                    strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dab.eps");
                     fin_coq = fopen_errchk (outfilename, "w");
                     sps.eps3d(fin_coq,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,4,(*magmom));
                    fclose (fin_coq);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dac.eps");
                    fin_coq = fopen_errchk (outfilename, "w");
                     sps.eps3d(fin_coq,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,5,(*magmom));
                    fclose (fin_coq);
                    strcpy(outfilename+11+strlen(ini.prefix),"spins3dbc.eps");
                    fin_coq = fopen_errchk (outfilename, "w");
                     sps.eps3d(fin_coq,text,inputpars.cs.abc,inputpars.cs.r,x,y,z,6,(*magmom));
                    fclose (fin_coq);
		    strcpy(outfilename+11+strlen(ini.prefix),"spins.eps");
                    fin_coq = fopen_errchk (outfilename, "w");
                     (*magmom).eps(fin_coq,text);
                    fclose (fin_coq);
                delete[]x;delete []y; delete []z;
		}
  delete magmom;if(verbose==1){printf(">");}
 //check if fecalculation gives again correct result
   if (physprops.fe>femin+(0.00001*fabs(femin))&&ini.maxnofmfloops>2){int eq=0;
   #ifndef _THREADS
   if(spsmin==sps){eq=1;};//take spinconfiguration which gave minimum free energy as starting value
     #else
   if(thrdat.spsmin==sps){eq=1;};//take spinconfiguration which gave minimum free energy as starting value
     #endif
   
   fprintf(stderr,"Warning htcalc.c: at T=%g K /  H= %g Tfemin=%4.9g was calc.(conf no %i),\n but recalculation  gives fe= %4.9gmeV -> no structure saved\n",
                            T,Norm(H),femin,physprops.j,physprops.fe);
   fprintf(stderr,"recalculation converged after %i loops and initial and final spin structures are ",r);
   if(eq==1){fprintf(stderr,"equal\n");}else{fprintf(stderr,"not equal\n");}
if (ini.logfevsQ==1) {strcpy(outfilename,"./results/");strcpy(outfilename+10,ini.prefix);
                       strcpy(outfilename+10+strlen(ini.prefix),"mcphas.log");
                       felog=fopen_errchk(outfilename,"a");
               fprintf(felog,"#Warning htcalc.c: at T=%g K /  H= %g Tfemin=%4.9g was calc.(conf no %i),\n# but recalculation  gives fe= %4.9gmeV -> no structure saved\n",
                T,Norm(H),femin,physprops.j,physprops.fe);fprintf(felog,"#recalculation converged after %i loops and initial and final spin structures are ",r);
   if(eq==1){fprintf(felog,"equal\n");}else{fprintf(stderr,"not equal\n#initial values as converged from femin=%4.9g meV calculation:\n",femin);
#ifndef _THREADS
     spsmin.print(felog);//take spinconfiguration which gave minimum free energy as starting value
     #else
     thrdat.spsmin.print(felog);//take spinconfiguration which gave minimum free energy as starting value
     #endif
     }
               fclose(felog);
	      }
                             physprops.sps.epsilon=0;physprops.Eel=0;
                             physprops.m=0;delete mf;return 2;
                             }
 //if(verbose==1){printf(".\n");}
 physpropclc(H,T,sps,(*mf),physprops,ini,inputpars);
      delete mf;
 }

return 0; // ok we are done with this (HT) point- return ok
 #if defined __linux__ && defined _THREADS
 pthread_exit(NULL);
 #endif
}





/*****************************************************************************/
// this sub checks if a spinconfiguration has already been added to
// table testspins and adds it if necessary
int checkspincf(int j,spincf & sps1,qvectors & testqs,Vector & nettom,
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
  { spq.spinfromq(testqs.na(-i),testqs.nb(-i),testqs.nc(-i),testqs.q(-i),
                   testqs.nettom(-i),testqs.momentq0(-i),testqs.phi(-i));
    if (spq==sps){
    physprops.j=i;return 1;} //ok
   
  }
 } 

   //  check initial config: take just used nettom,momentq0,phi for comparison
   if (j<0)
   {spq.spinfromq(testqs.na(-j),testqs.nb(-j),testqs.nc(-j),testqs.q(-j),
                  nettom,momentq0,phi);
    if (spq==sps) {physprops.j=j;testqs.nettom(-j)=nettom;
                  testqs.momentq0(-j)=momentq0;testqs.phi(-j)=phi;return 1;} //ok
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




