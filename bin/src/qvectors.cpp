#include "qvectors.hpp"
#include "myev.c"
 // *************************************************************************
 // ************************ qvectors *********************************
 // *************************************************************************
 // methods for class qvectors 

//save qvectors on file
void qvectors::save(const char * filemode)
{ FILE * fout;
 int k;
 spincf sps(1,1,1,nofatoms,3);
    
 fout=fopen_errchk (savfilename,filemode);
 fprintf (fout, "#{(h k l) of test-Qvectors (nr1 nr2 nr3) unitcellsize (no) structure number / momentum arrangement <J(i)>}\n");
 for (k=1;k<=nofqs();++k)
 {Vector hkl(1,3);
  hkl=rez.Transpose()*q(k);
  fprintf (fout, "%6.4g  %6.4g %6.4g      %4i %4i %4i %i\n",   
                  hkl(1),hkl(2),hkl(3),
		  na(k),nb(k),nc(k),-k);
    if (verbose){ sps.spinfromq(na(k),nb(k),nc(k),q(k),nettom(k),momentq0(k),phi(k));
                  sps.print(fout);
		}  
 }      
 fclose(fout);
 printf("file %s saved\n",savfilename);
}

// subs to calculate 3dim indices of qvector component
int qvectors::ia(int j){div_t result; result = div(j-1,hchkn[2][0]*hchkn[3][0]);
              return result.quot+1;}

int qvectors::ib(int j){div_t result; result = div(j-1,hchkn[2][0]*hchkn[3][0]);
              result=div(result.rem,hchkn[3][0]);
              return result.quot+1;}

int qvectors::ic(int j){div_t result; result = div(j-1,hchkn[2][0]*hchkn[3][0]);
              result=div(result.rem,hchkn[3][0]);
              return result.rem+1;}
	      

int qvectors::maxnofqs () //returns maximum nofqvectors
{return (int)hchkn[1][0]*hchkn[2][0]*hchkn[3][0];
}    

int qvectors::nofqs () //returns nofqvectors
{return nofq;
}    

Vector & qvectors::q(int i) // returns (hkl) of qvector i
{return (*q0[i]);
}

Vector & qvectors::nettom(int i) // returns pointer to nettom[i]
{return (*nm[i]);}

Vector & qvectors::momentq0(int i) // returns pointer to momentq0[i]
{ return (*mq0[i]);}

Vector & qvectors::phi(int i) // returns pointer to phi[i]
{return (*ph[i]);}

int qvectors::na (int i) // returns period for i.th qvector
{return (int)(*n[i])(1);
}

int qvectors::nb (int i) // returns period for i.th qvector
{return (int)(*n[i])(2);
}

int qvectors::nc (int i) // returns period for i.th qvector
{return (int)(*n[i])(3);
}

bool qvectors::is_in_1stBZ(Vector & hkl, Vector & abc,Matrix & rezijk)
 {int n=3;Vector qijk(1,3),R(1,3); 
  // rezijk is the primitive reciprocal lattice in ijk coordinate system
  hkl2ijk(qijk, hkl, abc); // qijk is Vector hkl in ijk coordinate system
  
   for(int n1=-n;n1<=n;++n1)
    for(int n2=-n;n2<=n;++n2)
     for(int n3=-n;n3<=n;++n3)
      {R=(double)n1*rezijk.Column(1)+(double)n2*rezijk.Column(2)+(double)n3*rezijk.Column(3);
//myPrintVector(stdout,qijk);myPrintVector(stdout,R);
//fprintf(stdout,"---\n");
       if(2*(qijk*R)>(R*R)+0.00001)return false;
      }
  return true;
 }
//constructor - generate set of qvectors
/* input 
       	filename	file to save q vector values on
        na              number of atoms in cryst. unit cell
	nmm             nofcomponents of moment vector
        v               verbose switch
*/
qvectors::qvectors (inipar & ini,par & inputpars,
                    Vector & mmax,const char * savfile, int v)
{ savfilename= new char [strlen(savfile)+1];
  strcpy(savfilename,savfile);
  verbose=v;
  float tryperiode, tryz, tryq;
  float delta,min,max;
  int i,j,maxdim[4],k,i1,i2,i3;
  hkl=Vector(1,3); 
  rez=Matrix(1,3,1,3);r=Matrix(1,3,1,3);
  rez=inputpars.rez;r=rez.Inverse();
  Matrix rt(1,3,1,3),rijk(1,3,1,3),rezijk(1,3,1,3);
 // calculate rezprim (primitive reciprocal basis in euclidean ijk coordinates)
  dadbdc2ijk(rijk,r,inputpars.cs.abc); // rijk is the primitive unit cell in ijk coordinates
  rezijk=2*PI*rijk.Inverse().Transpose();
  
  
  nofatoms=inputpars.cs.nofatoms;
  nofcomponents=inputpars.cs.nofcomponents;

if (verbose){printf ("\n");
printf ("  initialize qvector range: - min and maximum of components are calculated by\n");
printf ("          transforming points of the quader described by mcphas.ini\n"); 
printf ("	  (hklmin hklmax) to primitve recirpocal lattice coordinates and looking\n"); 
printf ("	  for the minimum/maximum of coordinates\n");
printf ("          - the minimal distance of qvectors with respect to reciprocal primitive\n"); 
printf ("	  lattice vector k is calculated the following way: deltahkl of nonprimitive\n"); 
printf ("	  lattice forms quader, the cornerpoints of the quader deltahkl are transformed\n");
printf ("	  to the primitive reciprocal lattice and maximum distant corner is used to\n");
printf ("	  determine the minimal distance of qvectors: NOTE -  consequently the\n");
printf ("          qvector density is lower than specified in the nonprimitive lattice !\n");
printf ("\n");
printf ("nonprimitive reciprocal lattice vectors:\n");
printf ("   REZ1   REZ2   REZ3\n");
printf ("    1      0      0\n");
printf ("    0      1      0\n");
printf ("    0      0      1\n");
printf ("range and spacings with respect to nonprimitive reciprocal lattice:\n");
printf ("vector	min     max    delta \n");
 for (k=1;k<=3;++k)
 {printf ("REZ%i     %g      %g     %g\n",k,ini.qmin(k),ini.qmax(k),ini.deltaq(k));
  if (ini.deltaq(k)==0){fprintf(stderr,"Error qvector.cpp: hkl-stepwidth(%i) is zero\n",k);exit(EXIT_FAILURE);}
 }
printf ("\n");
printf ("primitive reciprocal lattice vectors:\n");
printf ("  rez1    rez2    rez3\n");
  rt=rez.Transpose();myPrintMatrix(stdout,rt); 
printf ("\n");
printf ("range and spacings with respect to primitive reciprocal lattice:\n");
printf ("vector	min     max    delta \n");
}

// k=1,2,3 means create h,k,l net with respect to primitive reciprocal lattice  
 for (k=1;k<=3;++k)
 {
    // determine also minimal distance of qvectors with respect to reciprocal primitive lattice vector k
  // this is done the same way as below for the q-range(deltaq of nonprimitive lattice forms quader) consequently the
  // qvector density is lower than specified in the nonprimitive lattice
  min=1000; max=-1000; //initialize
  for (i1=0;i1<=1;++i1){hkl(1)=ini.deltaq(1)*i1; // loop cornerpoints
  for (i2=0;i2<=1;++i2){hkl(2)=ini.deltaq(2)*i2;
  for (i3=0;i3<=1;++i3){hkl(3)=ini.deltaq(3)*i3;
    if((tryq=(r.Transpose()*hkl)(k))<min){min=tryq;} //determine min/max 
    if (tryq                        >max){max=tryq;}
  }}}
  delta=max-min;
  //initialize range - min and maximum of components are calculated by
  // transforming points of the quader determined by ini.qmin and ini.qmax
  // to primitve recirpocal lattice coordinates and looking for the minimum/maximum
  // of coordinates
  min=1000; max=-1000; //initialize
  for (i1=0;i1<=(ini.qmax(1)-ini.qmin(1))/ini.deltaq(1)+1;++i1)
   {hkl(1)=ini.qmin(1)+ini.deltaq(1)*i1; // loop points in stepwidth 
  for (i2=0;i2<=(ini.qmax(2)-ini.qmin(2))/ini.deltaq(2)+1;++i2)
   {hkl(2)=ini.qmin(2)+ini.deltaq(2)*i2; // given by ini.deltaq
  for (i3=0;i3<=(ini.qmax(3)-ini.qmin(3))/ini.deltaq(3)+1;++i3)
   {hkl(3)=ini.qmin(3)+ini.deltaq(3)*i3;
    if((tryq=(r.Transpose()*hkl)(k))<min){min=tryq;} //determine min/max 
    if (tryq                        >max){max=tryq;}
  }}}

  if (verbose){printf ("rez%i     %g      %g     %g\nvalues: ",k,min,max,delta);}

  // generate memory for arrays
  maxdim[k]=(int)rint((max-min)/delta)*3+2;
  hchk[k]  = new float[maxdim[k]];
   if (hchk[k] == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  hchkz[k]  = new int[maxdim[k]];
   if (hchkz[k] == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  hchkn[k]  = new int[maxdim[k]];
   if (hchkn[k] == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  // generate hkl(k) of q-vectors
  hchk[k][0]=min;hchk[k][1]=max;
  hchkn[k][0]=0; //initialize hchkn
  
    for (tryperiode=1;tryperiode<ini.maxqperiod;tryperiode+=1.0)
    {for (tryz=integer(tryperiode*min/1.0);
          tryz<=integer(tryperiode*max/1.0)+1;++tryz)
     {tryq=1.0*tryz/tryperiode; 
      if (min<=tryq&&tryq<=max)
       {for (i=1;i<=abs(hchkn[k][0])&&hchk[k][i]<tryq;++i){};
          if (hchk[k][i-1]+delta<hchk[k][i] &&
	            (hchkn[k][0]==0||hchk[k][i-1]<tryq-delta/1.0) &&
	            (i>abs(hchkn[k][0])||tryq+delta/1.0<hchk[k][i]))
	      { // insert this  tryq into the table
	       if (hchkn[k][0]<maxdim[k]-1) 
	        {++hchkn[k][0];for(j=hchkn[k][0]+1;j>i;--j)
	            {hchk[k][j]=hchk[k][j-1];hchkz[k][j]=hchkz[k][j-1];hchkn[k][j]=hchkn[k][j-1];
		     };
	           hchk[k][i]=tryq;hchkz[k][i]=(int)tryz;hchkn[k][i]=(int)tryperiode;
		}
		 else
		 {fprintf(stderr,"dimension of q vector array to small");exit(EXIT_FAILURE);
		  }
	       }  
       }
     }
   }
 if (hchkn[k][0]==0){hchkn[k][0]=1;hchk[k][1]=0;hchkz[k][1]=0;hchkn[k][1]=1;}
 if (verbose){for(i=1;i<=hchkn[k][0];++i){printf(" %g ",hchk[k][i]);}printf("\n");}
 } //next k

 int l,m;
 Vector ddd(1,nofcomponents*nofatoms);
 Vector dd(1,3); // here memory could be saved by calculating nofq in advance
 // and reserving for the following fields only what is needed (nofq)
 q0 = new Vector * [maxnofqs()+2];        
 n = new Vector * [maxnofqs()+2];        
 nm = new Vector * [maxnofqs()+2];        
 mq0 = new Vector * [maxnofqs()+2];        
 ph = new Vector * [maxnofqs()+2];        
 // see what qvectors we can use (whether they lie in the specified region)
 // and store them
 nofq=0;if (verbose){printf("#Generating hkl in 1st Brillouin Zone\n#nr1 nr2 nr3 (h k l) \t (hprim kprim lprim)\n");}
 for (k=1;k<=maxnofqs();++k)
 {dd(1)=hchk[1][ia(k)];dd(2)=hchk[2][ib(k)];dd(3)=hchk[3][ic(k)];
//  printf("iaibic(%i %i %i)\n",ia(k),ib(k),ic(k));
  hkl=rez.Transpose()*dd;
  if (hchkn[1][ia(k)]*hchkn[2][ib(k)]*hchkn[3][ic(k)]<=ini.maxnofspins&& //not too many spins
      ini.qmin(1)-0.00001<=hkl(1)&&hkl(1)<=ini.qmax(1)+0.00001&&
      ini.qmin(2)-0.00001<=hkl(2)&&hkl(2)<=ini.qmax(2)+0.00001&&
      ini.qmin(3)-0.00001<=hkl(3)&&hkl(3)<=ini.qmax(3)+0.00001) //yes they are in the region-> increment nofq and store 
   { // !! here we could check if the qvector is in the 1st BZ and only
     // use it if it is ...!!
    if(is_in_1stBZ(hkl,inputpars.cs.abc,rezijk)){
    ++nofq;
    q0[nofq]=new Vector (1,3); (*q0[nofq])=dd;
    if (verbose){printf("# %2i  %2i  %2i (%g %g %g) \t (%g %g %g)\n",hchkn[1][ia(k)],hchkn[2][ib(k)],hchkn[3][ic(k)],hkl(1),hkl(2),hkl(3),dd(1),dd(2),dd(3));}
    
   dd(1)=hchkn[1][ia(k)];dd(2)=hchkn[2][ib(k)];dd(3)=hchkn[3][ic(k)];
    n[nofq]=new Vector (1,3); (*n[nofq])=dd;
   
   for(l=1;l<=nofatoms;++l)for(m=1;m<=nofcomponents;++m)
    {ddd(m+nofcomponents*(l-1))=mmax(m+nofcomponents*(l-1))*rnd(1);}  
   nm[nofq]=new Vector (1,nofcomponents*nofatoms); (*nm[nofq])=ddd;
   
   for(l=1;l<=nofatoms;++l)for(m=1;m<=nofcomponents;++m)
    {ddd(m+nofcomponents*(l-1))=rnd(1);}
   mq0[nofq]=new Vector (1,nofcomponents*nofatoms); (*mq0[nofq])=ddd;
   
   for(l=1;l<=nofatoms;++l)for(m=1;m<=nofcomponents;++m)
    {ddd(m+nofcomponents*(l-1))=rnd(1);}
    ph[nofq]=new Vector (1,nofcomponents*nofatoms); (*ph[nofq])=ddd;
       }
   }
  }

}



 //kopier-konstruktor
qvectors::qvectors (const qvectors & q)
{ savfilename= new char [strlen(q.savfilename)+1];
  strcpy(savfilename,q.savfilename);
  int k;
  hkl=Vector(1,3);
  rez=q.rez;r=q.r;nofq=q.nofq;
  nofatoms=q.nofatoms;
  nofcomponents=q.nofcomponents;
  verbose=q.verbose;
  
 for (k=1;k<=3;++k)
 {
  hchk[k] = new float[sizeof(*q.hchk[k])];
   if (hchk[k] == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
   memcpy (hchk[k], q.hchk[k], sizeof (*hchk[k]));
  hchkz[k] = new int[sizeof(*q.hchkz[k])];
   if (hchkz[k] == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
   memcpy (hchkz[k], q.hchkz[k], sizeof (*hchkz[k]));
  hchkn[k] = new int[sizeof(*q.hchkn[k])];
   if (hchkn[k] == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
   memcpy (hchkn[k], q.hchkn[k], sizeof (*hchkn[k])); 

  q0 = new Vector * [sizeof(*q.q0)];        
  n = new Vector * [sizeof(*q.n)];        
  nm = new Vector * [sizeof(*q.nm)];        
  mq0 = new Vector * [sizeof(*q.mq0)];        
  ph = new Vector * [sizeof(*q.ph)];        
 for (k=1;k<=nofqs();++k)
  {q0[k]=new Vector (1,3); (*q0[k])=(*q.q0[k]);
   n[k]=new Vector (1,3); (*n[k])=(*q.n[k]);
   nm[k]=new Vector (1,nofcomponents*nofatoms); (*nm[k])=(*q.nm[k]);
   mq0[k]=new Vector (1,nofcomponents*nofatoms); (*mq0[k])=(*q.mq0[k]);
   ph[k]=new Vector (1,nofcomponents*nofatoms); (*ph[k])=(*q.ph[k]);
  }

 }
}

//destruktor
qvectors::~qvectors ()
{//printf("hello destruktor qvectors\n");  
 int k;
 for (k=1;k<=nofqs();++k)
  {delete q0[k];delete n[k];delete nm[k];delete mq0[k];delete ph[k];
  }
 for (k=1;k<=3;++k)
 {delete []hchk[k]; delete []hchkz[k];delete []hchkn[k];}
  delete []q0;delete []n;
  delete []nm;delete []mq0;delete []ph;
  delete []savfilename;
//printf("hello destruktor qvectors\n");  
 
}
