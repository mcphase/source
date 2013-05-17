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
#define NUM_THREADS 4

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
htcalc_input *tin[NUM_THREADS];
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
 THRLC_SET(threadSpecificKey, myinput); double tlsfemin=1e10;  // Thread local variable to judge whether to print output
 #endif 
 int i,ii,iii,tryrandom,nr,rr,ri,is;
 double fe,fered;
 double u,lnz; // free- and magnetic energy per ion [meV]
 Vector momentq0(1,inputpars.nofcomponents*inputpars.nofatoms),phi(1,inputpars.nofcomponents*inputpars.nofatoms);
 Vector nettom(1,inputpars.nofcomponents*inputpars.nofatoms),q(1,3);
 Vector mmom(1,inputpars.nofcomponents);
 Vector h1(1,inputpars.nofcomponents),h1ext(1,3),hkl(1,3);
 h1ext=0;
 char text[10000];
 spincf  sps(1,1,1,inputpars.nofatoms,inputpars.nofcomponents),sps1(1,1,1,inputpars.nofatoms,inputpars.nofcomponents);
 mfcf * mf;
 mfcf * mf1;
 spincf * magmom;
 FILE * felog; // logfile for q dependence of fe
 FILE * fin_coq;

  for (tryrandom=0;tryrandom<=ini.nofrndtries&&j!=0;++tryrandom)
   {if (j>0){sps=(*testspins.configurations[j]);// take test-spinconfiguration
             #ifndef _THREADS
	     if (tryrandom==0&&verbose==1) { printf ( "conf. no %i (%ix%ix%i spins)"  ,j,sps.na(),sps.nb(),sps.nc()); fflush(stdout); }
             #endif 
            }
    else     // take q vector and choose phase and mom dir randomly
            {q=testqs.q(-j);  
	     if (tryrandom==0)
	     {nettom=testqs.nettom(-j);momentq0=testqs.momentq0(-j);phi=testqs.phi(-j);
	     }
	     else
	     {for(i=1;i<=inputpars.nofatoms;++i)
	      {for(ii=1;ii<=inputpars.nofcomponents;++ii)
	        {iii=inputpars.nofcomponents*(i-1)+ii;h1=0;h1(ii)=10*MU_B;
                 (*inputpars.jjj[i]).Icalc(mmom,T,h1,h1ext,lnz,u,(*inputpars.jjj[i]).Icalc_parstorage);
		 nettom(iii)=mmom(ii)*rnd(1);
	         momentq0(iii)=rnd(1);
	         phi(iii)=rnd(1)*3.1415;
		}
	      }
	     }
	     sps.spinfromq(testqs.na(-j),testqs.nb(-j),testqs.nc(-j),
	                   q,nettom,momentq0,phi);
             hkl=inputpars.rez.Transpose()*q;  
             #ifndef _THREADS
   	     if (tryrandom==0&&verbose==1) { printf ( "(hkl)=(%g %g %g)..(%ix%ix%i primitive unit cells) ",hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc()); fflush(stdout); }
             #endif 
	    }	 
    if (tryrandom>0){nr=(int)(rint(rnd(1.0)*(sps.n()*inputpars.nofatoms-1)))+1;
	             for (i=1;i<=nr;++i) //MonteCarlo randomize nr spins
                      {rr=(int)rint(rnd(1.0)*(sps.n()-1))+1;
		       ri=inputpars.nofcomponents*(int)rint(rnd(1.0)*(inputpars.nofatoms-1));
	               for(ii=1;ii<=inputpars.nofcomponents;++ii)
		       {sps.mi(rr)(ri+ii)=rnd(1.0) ;}
		       } // randomize spin rr
                    }
 
      //!!!calculate free energy - this is the heart of this loop !!!!
      mf=new mfcf(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,inputpars.nofcomponents);
      fe=fecalc(H ,T,ini,inputpars,sps,(*mf),u,testspins,testqs);
          
      // test spinconfiguration  and remember it                                    
      if (fe<femin)
            {               // first - reduce the spinconfiguration if possible
               if (verbose==1){fprintf(stdout,"fe(tryrandom=%i)= %f meV\n",tryrandom,fe); fflush(stdout);}
	       sps1=sps;sps1.reduce(); 
                   mf1=new mfcf(sps1.na(),sps1.nb(),sps1.nc(),inputpars.nofatoms,inputpars.nofcomponents);
               if ((fered=fecalc(H ,T,ini,inputpars,sps1,(*mf1),u,testspins,testqs))<=fe+1e-141){(*mf)=(*mf1);sps=sps1;}
                   magmom=new spincf(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,3);
                   int i1,j1,k1,l1,m1;Vector mom(1,3),d1(1,inputpars.nofcomponents);
                   for (l1=1;l1<=inputpars.nofatoms;++l1){
                    // go through magnetic unit cell and sum up the contribution of every atom
                  for(i1=1;i1<=sps.na();++i1){for(j1=1;j1<=sps.nb();++j1){for(k1=1;k1<=sps.nc();++k1){
                   for(m1=1;m1<=inputpars.nofcomponents;++m1){d1[m1]=(*mf).mf(i1,j1,k1)[inputpars.nofcomponents*(l1-1)+m1];}
                   (*inputpars.jjj[l1]).mcalc(mom,T,d1,H,(*inputpars.jjj[l1]).Icalc_parstorage);
                   for(m1=1;m1<=3;++m1){(*magmom).m(i1,j1,k1)(3*(l1-1)+m1)=mom(m1);}
                    }}}} 
                   delete mf1; 
                 // display spinstructure
                if (verbose==1)
                {Vector abc(1,6); abc(1)=inputpars.a; abc(2)=inputpars.b; abc(3)=inputpars.c;
                  abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;
                 float * x;x=new float[inputpars.nofatoms+1];float *y;y=new float[inputpars.nofatoms+1];float*z;z=new float[inputpars.nofatoms+1];
		 
		 for (is=1;is<=inputpars.nofatoms;++is)
		   {
                    x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
                     sprintf(text,"fe=%g,fered=%g<femin=%g:T=%gK, |H|=%gT,Ha=%gT, Hb=%gT, Hc=%gT,  %i spins",fe,fered,femin,T,Norm(H),H(1),H(2),H(3),sps.n());
                    fin_coq = fopen_errchk ("./results/.spins3dab.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,4,inputpars.gJ,(*magmom));
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dac.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,5,inputpars.gJ,(*magmom));
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dbc.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,6,inputpars.gJ,(*magmom));
                    fclose (fin_coq);
		   
                    fin_coq = fopen_errchk ("./results/.spins.eps", "w");
                     (*magmom).eps(fin_coq,text);
                    fclose (fin_coq);
		delete[]x;delete []y; delete []z;
	        }
                delete magmom;   
                           // see if spinconfiguration is already stored
             #ifndef _THREADS
	     if (0==checkspincf(j,sps,testqs,nettom,momentq0,phi,testspins,physprops,ini))//0 means error in checkspincf/addspincf
	        {fprintf(stderr,"Error htcalc: too many spinconfigurations created");
                 testspins.save(filemode);testqs.save(filemode); delete mf;
		 return 1;}
	     femin=fe; spsmin=sps;	   
            //printout fe
	    if (verbose==1) printf("fe=%gmeV, struc no %i in struct-table (initial values from struct %i)",fe,physprops.j,j);
             #else
             MUTEX_LOCK(&mutex_tests); 
             int checksret = checkspincf(j,sps,testqs,nettom,momentq0,phi,testspins,(*thrdat.physprops),ini); //0 means error in checkspincf/addspincf
             MUTEX_UNLOCK(&mutex_tests); 
	     if (checksret==0) {fprintf(stderr,"Error htcalc: too many spinconfigurations created");
                 testspins.save(filemode);testqs.save(filemode); delete mf;
		 goto ret;}
             MUTEX_LOCK (&mutex_min); if(fe<femin) { femin=fe; thrdat.spsmin=sps; } MUTEX_UNLOCK (&mutex_min); tlsfemin=femin;
             #endif
	     }
            delete mf;
            //printout fe
            #ifdef _THREADS
	    if (tryrandom==ini.nofrndtries && verbose==1) {
	       if(j>0) printf ( "conf. no %i (%ix%ix%i spins)"  ,j,sps.na(),sps.nb(),sps.nc());
               else    printf ( "(hkl)=(%g %g %g)..(%ix%ix%i primitive unit cells) ",hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc()); 
               if(tlsfemin==femin) printf("fe=%gmeV, struc no %i in struct-table (initial values from struct %i)",fe,(*thrdat.physprops).j,j); }
            #endif
            if (tryrandom==ini.nofrndtries&&verbose==1){printf("\n");}
 
	    
  // log fe if required
   if (ini.logfevsQ==1) {
                 ComplexVector a(1,3*inputpars.nofatoms),b(1,3*inputpars.nofatoms);
                 ComplexVector b1(1,inputpars.nofcomponents*inputpars.nofatoms);
                 float inmax=0;int qh,qk,ql,l;
                 ComplexVector * mq;  
                 mq = new ComplexVector [sps.in(sps.na(),sps.nb(),sps.nc())+2];for(l=0;l<=sps.in(sps.na(),sps.nb(),sps.nc())+1;++l){mq[l]=ComplexVector(1,inputpars.nofcomponents*inputpars.nofatoms);}
                 Vector sq2(1,3*inputpars.nofatoms),qs(1,3),qt(1,3);float in;qs(1)=1000;
                 sps.FT(mq); //Fourier trafo of spincf
		 // get the main propagation vector by looking for the
		 // biggest Fourier component of the magnetic moment arrangement 
                 for(qh=0;qh<=sps.na()/2;++qh){for(qk=0;qk<=sps.nb()/2;++qk){for(ql=0;ql<=sps.nc()/2;++ql)
                  {// get magnetic moment from momentum fouriercomponent into b 
		   b=0;
		   b1 = mq[sps.in(sps.na()-qh,sps.nb()-qk,sps.nc()-ql)];
                   for(l=1;l<=inputpars.nofatoms;++l)
		   {int m1,m1max=3; if ((*inputpars.jjj[l]).gJ==0){m1max=6;}
		    for (m1=1;m1<=m1max;++m1)
		     {if((*inputpars.jjj[l]).gJ==0)
		      {if(m1==2||m1==4||m1==6){b(3*(l-1)+(m1+1)/2)+=b1(inputpars.nofcomponents*(l-1)+m1);}
		       else                   {b(3*(l-1)+(m1+1)/2)+=2.0*b1(inputpars.nofcomponents*(l-1)+m1);}
		      }
		      else
		      {b(3*(l-1)+m1)=b1(inputpars.nofcomponents*(l-1)+m1)*(*inputpars.jjj[l]).gJ;
		      }
		     }    
		    }
		   a = b.Conjugate();
		   b1 = mq[sps.in(qh,qk,ql)];
		   b=0;
                   for(l=1;l<=inputpars.nofatoms;++l)
		   {int m1,m1max=3; if ((*inputpars.jjj[l]).gJ==0){m1max=6;}
		    for (m1=1;m1<=m1max;++m1)
		     {if((*inputpars.jjj[l]).gJ==0)
		      {if(m1==2||m1==4||m1==6){b(3*(l-1)+(m1+1)/2)+=b1(inputpars.nofcomponents*(l-1)+m1);}
		       else                   {b(3*(l-1)+(m1+1)/2)+=2.0*b1(inputpars.nofcomponents*(l-1)+m1);}
		      }
		      else
		      {b(3*(l-1)+m1)=b1(inputpars.nofcomponents*(l-1)+m1)*(*inputpars.jjj[l]).gJ;
		      }
		     }    
		    }                   
		   // inner product
                   sq2=Abs(b+a)/(double)sps.n()/(double)inputpars.nofatoms;
                   Vector q(1,3);
		   q(1)=1.0*qh/sps.na();
	           q(2)=1.0*qk/sps.nb();
                   q(3)=1.0*ql/sps.nc();
                   qt=inputpars.rez.Transpose()*q;
		   in=Norm(sq2)*Norm(sq2);
	           if ((in>inmax-0.01&&Norm(qt)<Norm(qs))
		      ||in>inmax)
                    {inmax=in;qs=qt;}
                   }}}
                   felog=fopen_errchk("./results/mcphas.log","a");
                   if (verbose==1||fe>10000){fprintf(felog,"#");}
                   fprintf(felog,"%10.6g %10.6g %10.6g %10.6g %3i %3i %3i %3i \n",qs(1),qs(2),qs(3),fe,j,sps.na(),sps.nb(),sps.nc());
                   if (verbose==1&&fe<20000){sps.print(felog);}
	           fclose(felog);
                  delete []mq;
                 }

      }
      #ifndef _THREADS
      return 1;
      #else
      ret:;
      MUTEX_LOCK(&mutex_loop);
      thrdat.thread_id = thread_id;
      EVENT_SIG(checkfinish);
      MUTEX_UNLOCK(&mutex_loop);
      #undef ini
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
 // returns 1 if too maxnofspinconfigurations is exceeded 
 // returns 2 if no spinconfiguration has been found at ht point
 */
 int i,j,k,is;
 Vector momentq0(1,inputpars.nofcomponents*inputpars.nofatoms),phi(1,inputpars.nofcomponents*inputpars.nofatoms);
 Vector nettom(1,inputpars.nofcomponents*inputpars.nofatoms),q(1,3);
 Vector h1(1,inputpars.nofcomponents),hkl(1,3);
 Vector H(1,3); // magnetic field in ijk coordinate system
 Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1; // trick to get Habc as components along a,b,c
                  abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;
 dadbdc2ijk(H,Habc,abc); // transform Habc to ijk coordinates ... this is H
                  abc(1)=inputpars.a; abc(2)=inputpars.b; abc(3)=inputpars.c;
 double femin=10000;char text[1000];
 spincf  sps(1,1,1,inputpars.nofatoms,inputpars.nofcomponents),sps1(1,1,1,inputpars.nofatoms,inputpars.nofcomponents);
 spincf  spsmin(1,1,1,inputpars.nofatoms,inputpars.nofcomponents);
 mfcf * mf;
 spincf * magmom;
 FILE * felog; // logfile for q dependence of fe
 FILE * fin_coq;

if (T<=0.01){fprintf(stderr," ERROR htcalc - temperature too low - please check mcphas.ini !");exit(EXIT_FAILURE);}

 srand(time(0)); // initialize random number generator
 checkini(testspins,testqs,ini); // check if user pressed a button
 if (ini.logfevsQ==1) {felog=fopen_errchk("./results/mcphas.log","a");
               fprintf(felog,"#Logging of h k l fe[meV] spinconf_nr n1xn2xn3 T=%g Ha=%g Hb=%g Hc=%g\n",T,Habc(1),Habc(2),Habc(3));
               fclose(felog);
	      }
 if (verbose==1)
 { fin_coq= fopen_errchk ("./results/.fe_status.dat","w");
   #ifndef _THREADS
   fprintf(fin_coq,"#displayxtext=time(s)\n#displaytitle=2:log(iterations) 3:log(sta) 4:spinchange 5:stepratio 6:successrate(%%)\n#time(s) log(iteration) log(sta) spinchange stepratio  successrate=(nof stabilised structures)/(nof initial spinconfigs)\n");
   #else
   fprintf(fin_coq,"#displayxtext=time(s)\n#displaytitle=2:log(iterations) 3:log(sta) 4:spinchange 5:stepratio 6:successrate 7:threadID(%%)\n#time(s) log(iteration) log(sta) spinchange stepratio  successrate=(nof stabilised structures)/(nof initial spinconfigs)  thread_id \n");
   #endif
   fclose(fin_coq);	      
   printf("\n starting T=%g Ha=%g Hb=%g Hc=%g with \n %i spinconfigurations read from mcphas.tst and table \nand\n %i spinconfigurations created from hkl's\n\n",T,Habc(1),Habc(2),Habc(3),testspins.n,testqs.nofqs());
 }

 j=-testqs.nofqs()+(int)rint(rnd(testspins.n+testqs.nofqs())); 
    //begin with j a random number, j<0 means test spinconfigurations 
    //constructed from q vector set testqs, j>0 means test spinconfigurations from
    //set testspins
    //j=0;  //uncomment this for debugging purposes
    j = -testqs.nofqs()-1;

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
                  for (int ithread=0; ithread<NUM_THREADS; ithread++) tin[ithread] = new htcalc_input(0,ithread,&inputpars);}
   
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
       if(rc) { printf("Error return code %i from thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
       #else
       threads[ithread] = CreateThread(NULL, 0, htcalc_iteration, (void *) tin[ithread], 0, &tid[ithread]);
       if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %i from thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
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

if (femin>=10000) // did we find a stable structure ??
 {fprintf(stderr,"Warning propcalc: femin positive ... no stable structure found at  T= %g K / Ha= %g Hb= %g Hc= %g  T\n",
                 physprops.T,physprops.H(1),physprops.H(2),physprops.H(3));return 2;}
else // if yes ... then
 {// calculate physical properties
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
    else     // ---- or take q vector
            { sps.spinfromq(testqs.na(-physprops.j),testqs.nb(-physprops.j),
	              testqs.nc(-physprops.j),testqs.q(-physprops.j),
		      testqs.nettom(-physprops.j),testqs.momentq0(-physprops.j),
		      testqs.phi(-physprops.j));
	      }

     #ifndef _THREADS
     sps=spsmin;//take spinconfiguration which gave minimum free energy as starting value
     #else
     sps=thrdat.spsmin;//take spinconfiguration which gave minimum free energy as starting value
     #endif
   //MR 120221 removed spinconf invert in case nettoI is negative
  // now really calculate the physical properties
      mf=new mfcf(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,inputpars.nofcomponents);
      physprops.fe=fecalc(H ,T,ini,inputpars,sps,(*mf),physprops.u,testspins,testqs); 
      magmom=new spincf(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,3);
                   int i1,j1,k1,l1,m1;Vector mom(1,3),d1(1,inputpars.nofcomponents);
                   for (l1=1;l1<=inputpars.nofatoms;++l1){
                    // go through magnetic unit cell and sum up the contribution of every atom
                  for(i1=1;i1<=sps.na();++i1){for(j1=1;j1<=sps.nb();++j1){for(k1=1;k1<=sps.nc();++k1){
                  for(m1=1;m1<=inputpars.nofcomponents;++m1){d1[m1]=(*mf).mf(i1,j1,k1)[inputpars.nofcomponents*(l1-1)+m1];}                  
                   (*inputpars.jjj[l1]).mcalc(mom,T,d1,H,(*inputpars.jjj[l1]).Icalc_parstorage);
                    for(m1=1;m1<=3;++m1){(*magmom).m(i1,j1,k1)(3*(l1-1)+m1)=mom(m1);}
                    }}}}
             // display spinstructure
                if (verbose==1)
                {
		 float * x;x=new float[inputpars.nofatoms+1];float *y;y=new float[inputpars.nofatoms+1];float*z;z=new float[inputpars.nofatoms+1];
		 for (is=1;is<=inputpars.nofatoms;++is)
		   {x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
                     sprintf(text,"recalculated: fe=%g,femin=%g:T=%gK,|H|=%gT,Ha=%gT, Hb=%gT, Hc=%gT, %i spins",physprops.fe,femin,T,Norm(H),physprops.H(1),physprops.H(2),physprops.H(3),sps.n());
                    fin_coq = fopen_errchk ("./results/.spins3dab.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,4,inputpars.gJ,(*magmom));
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dac.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,5,inputpars.gJ,(*magmom));
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dbc.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,6,inputpars.gJ,(*magmom));
                    fclose (fin_coq);
		 fin_coq = fopen_errchk ("./results/.spins.eps", "w");
                     (*magmom).eps(fin_coq,text);
                    fclose (fin_coq);
                delete[]x;delete []y; delete []z;
		}
  delete magmom;
 //check if fecalculation gives again correct result
   if (physprops.fe>femin+(0.00001*fabs(femin))){fprintf(stderr,"Warning htcalc.c: at T=%g K /  H= %g Tfemin=%4.9g was calc.(conf no %i),\n but recalculation  gives fe= %4.9gmeV -> no structure saved\n",
                            T,Norm(H),femin,physprops.j,physprops.fe);delete mf;return 2;}
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
return (physprops.j=testspins.addspincf(sps));  //ok=1
}




