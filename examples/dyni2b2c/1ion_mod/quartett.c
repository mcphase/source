  
#include <cmath>
#include <cstring>
#include <vector.h>
#include <complex>
#include <martin.h>
#include <mcdisp.h>

#define K_B  0.08617343183 

class ionpars  //class for loading ion parameters
{private: 
 public:
   double a,b,c,d,e,delta;
   ComplexMatrix Ja; ComplexMatrix Jb; ComplexMatrix Jc; 
   ionpars(const char *file);
   ~ionpars();
   ionpars(const ionpars & p);
};

ionpars::ionpars (const ionpars & p) //copy constructor
 { Ja=p.Ja; Jb=p.Jb; Jc=p.Jc; a=p.a;b=p.b;c=p.c;d=p.d;e=p.e;delta=p.delta;}

ionpars::~ionpars(){} //destructor

ionpars::ionpars(const char * file) //constructor
{ FILE * fin;
  char instr[MAXNOFCHARINLINE]; 
  errno = 0;
  printf("# reading file %s\n",file);
  fin = fopen_errchk (file, "r+");
  while (fgets(instr,MAXNOFCHARINLINE,fin)!=NULL)
  {extract(instr,"a",a);
   extract(instr,"b",b);extract(instr,"c",c);
   extract(instr,"d",d);extract(instr,"e",e);
   extract(instr,"delta",delta);
   }
  fclose (fin);
  printf("# initializing  single ion module\n");
   Ja = ComplexMatrix(1,4,1,4); Ja=0;
   Ja(1,2)=b;Ja(2,1)=b;
   Ja(1,4)=c;Ja(2,3)=c;
   Ja(3,2)=c;Ja(4,1)=c;
   Ja(3,4)=e;Ja(4,3)=e;

   Jb = ComplexMatrix(1,4,1,4); Jb=0;complex<double> im(0,1);
   Jb(1,2)=-im*b;
   Jb(2,1)=b*im;
   Jb(1,4)=c*im;Jb(2,3)=-c*im;
   Jb(3,2)=c*im;Jb(4,1)=-c*im;
   Jb(3,4)=-e*im;Jb(4,3)=e*im;

   Jc = ComplexMatrix(1,4,1,4); Jc=0;
   Jc(1,1)=a;Jc(2,2)=-a;
   Jc(3,3)=-d;Jc(4,4)=d;
printf("Co2p: a=%4.6g b=%4.6g c=%4.6g d=%4.6g e=%4.6g delta=%4.6g\n",a,b,c,d,e,delta);
   
   a=-a;b=-b;c=-c;d=-d;e=-e;   
}

// example c file for dynamically loadable module of program mcphas
// routine Icalc for quasi-quartett in Co2p
  static ionpars Co2p("./Dy.sipf");  // get 1ion parameters from file (once at loading of quartett.so)
#ifdef __MINGW32__
extern "C" __declspec(dllexport) void Icalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      double * lnZ,double * U,ComplexMatrix & est)
#else
extern "C" void Icalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      double * lnZ,double * U,ComplexMatrix & est)
#endif
{//MIND: ABC not used here !!!!

if(J.Hi()!=3||gjmbHxc.Hi()!=3)
   {fprintf(stderr,"Error loadable module quartett.so: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}
Vector gjmbH(1,gjmbHxc.Hi());
gjmbH=gjmbHxc+(*g_J)*MU_B*Hext(1,3);
// check dimensions of vector

//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
// 			   int sort, int maxiter)
   // setup hamiltonian
   Matrix Ham(1,4,1,4); 
   Ham(1,1)=gjmbH[3]*Co2p.a-Co2p.delta/2; Ham(1,2)=gjmbH[2]*Co2p.b;Ham(1,3)=0;Ham(1,4)=-gjmbH[2]*Co2p.c;
   Ham(2,1)=gjmbH[1]*Co2p.b;Ham(2,2)=-gjmbH[3]*Co2p.a-Co2p.delta/2;Ham(2,3)=gjmbH[2]*Co2p.c;Ham(2,4)=0;
   Ham(3,1)=0;Ham(3,2)=gjmbH[1]*Co2p.c;Ham(3,3)=-gjmbH[3]*Co2p.d+Co2p.delta/2;Ham(3,4)=gjmbH[2]*Co2p.e;
   Ham(4,1)=gjmbH[1]*Co2p.c;Ham(4,2)=0;Ham(4,3)=gjmbH[1]*Co2p.e;Ham(4,4)=gjmbH[3]*Co2p.d+Co2p.delta/2;
/*   int i1,j1; //printout matrix
   for (i1=1;i1<=4;++i1){
    for (j1=1;j1<=4;++j1) printf ("%4.6g ",Ham(i1,j1));
    printf ("\n");
    }*/
    
   // diagonalize
   Vector En(1,4);Matrix zr(1,4,1,4);Matrix zi(1,4,1,4);
   int sort=0;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);
   // calculate Z and wn (occupation probability)
     Vector wn(1,4);
     double x,y,Z;int i;
     x=Min(En);
     for (i=1;i<=4;++i)
     {if ((y=(En(i)-x)/K_B/(*T))<700) wn[i]=exp(-y); 
      else wn[i]=0.0;}
     Z=Sum(wn);wn/=Z;  
     (*lnZ)=log(Z)-x/K_B/(*T);
   // calculate U
     (*U)=En*wn;
   // calculate ma,mb,mc
     ComplexMatrix z(1,4,1,4);
     ComplexMatrix za(1,4,1,4);
     ComplexMatrix zb(1,4,1,4);
     ComplexMatrix zc(1,4,1,4);
  
     z=ComplexMatrix(zr,zi);
     
     za=Co2p.Ja*z;
     zb=Co2p.Jb*z;
     zc=Co2p.Jc*z;

   J[1]=0;J[2]=0;J[3]=0;
    ComplexVector ddd;
    for (i=1;i<=4;++i)
    {
     J[1]+=wn(i)*real(z.Column(i)*za.Column(i));
     J[2]+=wn(i)*real(z.Column(i)*zb.Column(i));
     J[3]+=wn(i)*real(z.Column(i)*zc.Column(i));
    }     
  
return;
}
/**************************************************************************/
/**************************************************************************/
#ifdef __MINGW32__
extern "C" __declspec(dllexport) void mcalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      ComplexMatrix & est)
#else
extern "C" void mcalc(Vector & J,double * T, Vector & gjmbHxc,Vector & Hext,double * g_J, Vector & ABC,char ** sipffile,
                      ComplexMatrix & est)
#endif
{double lnZ,U;
 Icalc(J,T,gjmbHxc,Hext,g_J,ABC,sipffile,&lnZ,&U,est);
 J*=(*g_J);
}
/**************************************************************************/
// for mcdisp this routine is needed
#ifdef __MINGW32__
extern "C" __declspec(dllexport) int du1calc(int & tn,double & T, Vector & Hxc,Vector & Hext,double * g_J,Vector & MODPAR, char ** sipffile,
                       ComplexVector & u1,float & delta,int & n, int & nd,ComplexMatrix & est)
#else
extern "C" int du1calc(int & tn,double & T, Vector & Hxc,Vector & Hext,double * g_J,Vector & MODPAR, char ** sipffile,
                       ComplexVector & u1,float & delta,int & n, int & nd,ComplexMatrix & est)
#endif
{ 
  /*on input
    tn          transition-number - meaningless for kramers doublet, because there is only one transition
    MODPAR         A,M,Ci...saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    g_J		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splitting of kramers doublet [meV]
    u1(i)	<-|Ji-<Ji>|+> sqrt(tanh(delta/2kT))
*/
  double Z,lnz,u;
  static int pr;
  static Vector gjmbHin(1,3);
gjmbHin=Hxc+(*g_J)*MU_B*Hext;
  static Vector Jin(1,3);
  static Vector J(1,3);
  // clalculate thermal expectation values (needed for quasielastic scattering)
  Jin=0;if(T>0){ Icalc(Jin,&T,Hxc,Hext,g_J,MODPAR,sipffile,&lnz,&u,est);}
                else {T=-T;}
 pr=0;
  if (tn<0) {pr=1;tn*=-1;}
if(J.Hi()!=3||Hxc.Hi()!=3)
   {fprintf(stderr,"Error loadable module quartett.so: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}
Vector gjmbH(1,Hxc.Hi());
gjmbH=Hxc+(*g_J)*MU_B*Hext;

   // setup hamiltonian
   Matrix Ham(1,4,1,4); 
   Ham(1,1)=gjmbH[3]*Co2p.a-Co2p.delta/2; Ham(1,2)=gjmbH[2]*Co2p.b;Ham(1,3)=0;Ham(1,4)=-gjmbH[2]*Co2p.c;
   Ham(2,1)=gjmbH[1]*Co2p.b;Ham(2,2)=-gjmbH[3]*Co2p.a-Co2p.delta/2;Ham(2,3)=gjmbH[2]*Co2p.c;Ham(2,4)=0;
   Ham(3,1)=0;Ham(3,2)=gjmbH[1]*Co2p.c;Ham(3,3)=-gjmbH[3]*Co2p.d+Co2p.delta/2;Ham(3,4)=gjmbH[2]*Co2p.e;
   Ham(4,1)=gjmbH[1]*Co2p.c;Ham(4,2)=0;Ham(4,3)=gjmbH[1]*Co2p.e;Ham(4,4)=gjmbH[3]*Co2p.d+Co2p.delta/2;
    
   // diagonalize
   Vector En(1,4);Matrix zr(1,4,1,4);Matrix zi(1,4,1,4);
   int sort=0;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);
   // calculate Z and wn (occupation probability)
     Vector wn(1,4);
     double x,y;int i,j,k;
     x=Min(En);
     for (i=1;i<=4;++i)
     {if ((y=(En(i)-x)/K_B/T)<700) wn[i]=exp(-y); 
      else wn[i]=0.0;}
     Z=Sum(wn);wn/=Z;  
   // calculate ma,mb,mc
     ComplexMatrix z(1,4,1,4);
     ComplexMatrix za(1,4,1,4);
     ComplexMatrix zb(1,4,1,4);
     ComplexMatrix zc(1,4,1,4);
  
     z=ComplexMatrix(zr,zi);
     za=Co2p.Ja*z;
     zb=Co2p.Jb*z;
     zc=Co2p.Jc*z;

int dj=4;  
// calculate u1 and delta for transition number tn
// 1. get i and j from tn
k=0;for(i=1;i<=dj;++i){for(j=i;j<=dj;++j)
{++k;if(k==tn)break;}if(k==tn)break;}
n=i;nd=j;
//printf("tnij=%i %i %i\n",tn,i,j);
// 2. set delta
delta=En(j)-En(i);

if (delta<-0.000001){fprintf(stderr,"ERROR module quartett.so - du1calc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(j==i)delta=-SMALL_QUASIELASTIC_ENERGY; //if transition within the same level: take negative delta !!- this is needed in routine intcalc


// 3. set mat

if(i==j){//take into account thermal expectation values <Jl>
          u1(1)=(za.Column(j)*z.Column(i))-Jin(1);
          u1(2)=(zb.Column(j)*z.Column(i))-Jin(2);
          u1(3)=(zc.Column(j)*z.Column(i))-Jin(3);
     
        }
 else    {u1(1)=(za.Column(j)*z.Column(i));
          u1(2)=(zb.Column(j)*z.Column(i));
          u1(3)=(zc.Column(j)*z.Column(i));
     }
           // ... in complex vector scalar product a*b is defined 
          //as: a.conj(b) !!! (see cvector.cc)

if (delta>SMALL_QUASIELASTIC_ENERGY)
   {u1*=sqrt(wn(i)-wn(j)); // occupation factor
    if(pr==1)
      {printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
       printf(" |<%i|Ja|%i>|^2=%4.4g |<%i|Jb|%i>|^2=%4.4g |<%i|Jc|%i>|^2=%4.4g",i,j,abs(u1(1))*abs(u1(1)),i,j,abs(u1(2))*abs(u1(2)),i,j,abs(u1(3))*abs(u1(3)));
       printf(" n%i-n%i=%4.4g\n",i,j,wn(i)-wn(j));
      }
   }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
      if (pr==1)
      {printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
       printf(" |<%i|Ja-<Ja>|%i>|^2=%4.4g |<%i|Jb-<Jb>|%i>|^2=%4.4g |<%i|Jc-<Jc>|%i>|^2=%4.4g",i,j,abs(u1(1))*abs(u1(1)),i,j,abs(u1(2))*abs(u1(2)),i,j,abs(u1(3))*abs(u1(3)));
       printf(" n%i=%4.4g\n",i,wn(i));
      }
    u1*=sqrt(wn(i)/K_B/T);
   }
     
return 10; // return number of all transitions
}

/**************************************************************************/
// for mcdisp this routine is needed
#ifdef __MINGW32__
extern "C" __declspec(dllexport) int dm1(int & tn,double & T, Vector & Hxc,Vector & Hext,double * g_J,Vector & MODPAR, char ** sipffile,
                       ComplexVector & m1,float & maxE,ComplexMatrix & est)
#else
extern "C" int dm1(int & tn,double & T, Vector & Hxc,Vector & Hext,double * g_J,Vector & MODPAR, char ** sipffile,
                       ComplexVector & m1,float & maxE,ComplexMatrix & est)
#endif
{int nnt,n,nd;
nnt=du1calc(tn,T,Hxc,Hext,g_J,MODPAR,sipffile,m1,maxE,n,nd,est);m1*=(*g_J);return nnt; 
}