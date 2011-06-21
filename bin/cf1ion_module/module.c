//\begin{verbatim}
// example c file for dynamically loadable module of program
// mcphas ... it must not be c++, but pure c compiled with gcc and linked 
// with ld  !! The calculation has been compared to the internal (doublet)
// routine of mcphas
#include <cstdio>
#include <cmath>
#include <complex>
#include "vector.h"
#include <stdlib.h>

#define MU_B 0.05788
#define K_B  0.0862
#define SMALL 1e-10

#define UNUSED_PARAMETER(a) (void)a

// use cfield as a module loaded at runtime

extern "C" void cfield0_mcphas(double ** hcfr,double ** hcfi, double ** Jxr,double ** Jxi,  double ** Jyr, double ** Jyi, double ** Jzr, double ** Jzi,
                              double ** mo20r, double ** mo20i,
                              double ** mo22r, double ** mo22i,
                              double ** mo40r, double ** mo40i,
                              double ** mo42r, double ** mo42i,
                              double ** mo44r, double ** mo44i,
                              double ** mo60r, double ** mo60i,
                              double ** mo62r, double ** mo62i,
                              double ** mo64r, double ** mo64i,
                              double ** mo66r, double ** mo66i,
			      int * dimj);

class module_ionpars  //class for loading ion parameters
{private: 
 public:
   double J; // momentum quantum number
   Matrix Ja; Matrix Jb; Matrix Jc; Matrix Hcf;
   ComplexMatrix Jaa;
   ComplexMatrix Jbb;
   ComplexMatrix Jcc;
 
   Matrix O20;ComplexMatrix OO20;
   Matrix O22;ComplexMatrix OO22;
   Matrix O40;ComplexMatrix OO40;
   Matrix O42;ComplexMatrix OO42;
   Matrix O44;ComplexMatrix OO44;
   Matrix O60;ComplexMatrix OO60;
   Matrix O62;ComplexMatrix OO62;
   Matrix O64;ComplexMatrix OO64;
   Matrix O66;ComplexMatrix OO66;

   module_ionpars(const char * text);
   ~module_ionpars();
   module_ionpars(const module_ionpars & p);
};

module_ionpars::module_ionpars (const module_ionpars & p) //copy constructor
 {Ja=p.Ja; Jb=p.Jb; Jc=p.Jc;Hcf=p.Hcf;
  O20=p.O20;OO20=p.OO20;
  O22=p.O22;OO22=p.OO22;
  O40=p.O40;OO40=p.OO40;
  O42=p.O42;OO42=p.OO42;
  O44=p.O44;OO44=p.OO44;
  O60=p.O60;OO60=p.OO60;
  O62=p.O62;OO62=p.OO62;
  O64=p.O64;OO64=p.OO64;
  O66=p.O66;OO66=p.OO66;
  
 }

module_ionpars::~module_ionpars(){} //destructor

module_ionpars::module_ionpars(const char * text) //constructor
{ double ** hcfr,**hcfi,**Jxr,**Jxi,**Jyr,**Jyi,**Jzr,**Jzi;
  int dimj;complex<double> im(0,1);
  int i,j,dj=30; //30 ... maximum number of 2j+1
  printf("%s",text);
  
    Jxr=new double*[dj+1];Jxi=new double*[dj+1];
    Jyr=new double*[dj+1];Jyi=new double*[dj+1];
    Jzr=new double*[dj+1];Jzi=new double*[dj+1];
    hcfr=new double*[dj+1];hcfi=new double*[dj+1];

    for (i=1;i<=dj;++i)
     {Jxr[i]=new double [dj];Jxi[i]=new double [dj];
      Jyr[i]=new double [dj];Jyi[i]=new double [dj];
      Jzr[i]=new double [dj];Jzi[i]=new double [dj];
      hcfr[i]=new double [dj];hcfi[i]=new double [dj];
     }

  double ** mo20r,**mo20i;mo20r=new double*[dj+1];mo20i=new double*[dj+1];
    for (i=1;i<=dj;++i){mo20r[i]=new double [dj];mo20i[i]=new double [dj];}
  double ** mo22r,**mo22i;mo22r=new double*[dj+1];mo22i=new double*[dj+1];
    for (i=1;i<=dj;++i){mo22r[i]=new double [dj];mo22i[i]=new double [dj];}
  double ** mo40r,**mo40i;mo40r=new double*[dj+1];mo40i=new double*[dj+1];
    for (i=1;i<=dj;++i){mo40r[i]=new double [dj];mo40i[i]=new double [dj];}
  double ** mo42r,**mo42i;mo42r=new double*[dj+1];mo42i=new double*[dj+1];
    for (i=1;i<=dj;++i){mo42r[i]=new double [dj];mo42i[i]=new double [dj];}
  double ** mo44r,**mo44i;mo44r=new double*[dj+1];mo44i=new double*[dj+1];
    for (i=1;i<=dj;++i){mo44r[i]=new double [dj];mo44i[i]=new double [dj];}
  double ** mo60r,**mo60i;mo60r=new double*[dj+1];mo60i=new double*[dj+1];
    for (i=1;i<=dj;++i){mo60r[i]=new double [dj];mo60i[i]=new double [dj];}
  double ** mo62r,**mo62i;mo62r=new double*[dj+1];mo62i=new double*[dj+1];
    for (i=1;i<=dj;++i){mo62r[i]=new double [dj];mo62i[i]=new double [dj];}
  double ** mo64r,**mo64i;mo64r=new double*[dj+1];mo64i=new double*[dj+1];
    for (i=1;i<=dj;++i){mo64r[i]=new double [dj];mo64i[i]=new double [dj];}
  double ** mo66r,**mo66i;mo66r=new double*[dj+1];mo66i=new double*[dj+1];
    for (i=1;i<=dj;++i){mo66r[i]=new double [dj];mo66i[i]=new double [dj];}



  cfield0_mcphas(hcfr,hcfi, Jxr,Jxi,  Jyr, Jyi, Jzr, Jzi,
  mo20r,mo20i,
  mo22r,mo22i,
  mo40r,mo40i,
  mo42r,mo42i,
  mo44r,mo44i,
  mo60r,mo60i,
  mo62r,mo62i,
  mo64r,mo64i,
  mo66r,mo66i,
  &dimj);

   J=((double)dimj-1)/2;

  printf("#J=%g\n",J);

   Ja = Matrix(1,dimj,1,dimj); 
   Jaa = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jaa(i,j)=im*(Jxi[i])[j]+(Jxr[i])[j];
    if(i<j){Ja(i,j)=(Jxi[j])[i];}else{Ja(i,j)=(Jxr[i])[j];}
   }


   Jb = Matrix(1,dimj,1,dimj); 
   Jbb = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jbb(i,j)=im*(Jyi[i])[j]+(Jyr[i])[j];
    if(i<j){Jb(i,j)=(Jyi[j])[i];}else{Jb(i,j)=(Jyr[i])[j];}
   }


   Jc = Matrix(1,dimj,1,dimj); 
   Jcc = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jcc(i,j)=im*(Jzi[i])[j]+(Jzr[i])[j];
    if(i<j){Jc(i,j)=(Jzi[j])[i];}else{Jc(i,j)=(Jzr[i])[j];}
   }
   
   O20 = Matrix(1,dimj,1,dimj); 
   OO20 = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {OO20(i,j)=im*(mo20i[i])[j]+(mo20r[i])[j];
    if(i<j){O20(i,j)=(mo20i[j])[i];}else{O20(i,j)=(mo20r[i])[j];}
   }}

   O22 = Matrix(1,dimj,1,dimj); 
   OO22 = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {OO22(i,j)=im*(mo22i[i])[j]+(mo22r[i])[j];
    if(i<j){O22(i,j)=(mo22i[j])[i];}else{O22(i,j)=(mo22r[i])[j];}
   }}

   O40 = Matrix(1,dimj,1,dimj); 
   OO40 = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {OO40(i,j)=im*(mo40i[i])[j]+(mo40r[i])[j];
    if(i<j){O40(i,j)=(mo40i[j])[i];}else{O40(i,j)=(mo40r[i])[j];}
   }}
   
   O42 = Matrix(1,dimj,1,dimj); 
   OO42 = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {OO42(i,j)=im*(mo42i[i])[j]+(mo42r[i])[j];
    if(i<j){O42(i,j)=(mo42i[j])[i];}else{O42(i,j)=(mo42r[i])[j];}
   }}
   
   O44 = Matrix(1,dimj,1,dimj); 
   OO44 = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {OO44(i,j)=im*(mo44i[i])[j]+(mo44r[i])[j];
    if(i<j){O44(i,j)=(mo44i[j])[i];}else{O44(i,j)=(mo44r[i])[j];}
   }}
   
   O60 = Matrix(1,dimj,1,dimj); 
   OO60 = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {OO60(i,j)=im*(mo60i[i])[j]+(mo60r[i])[j];
    if(i<j){O60(i,j)=(mo60i[j])[i];}else{O60(i,j)=(mo60r[i])[j];}
   }}
   
   O62 = Matrix(1,dimj,1,dimj); 
   OO62 = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {OO62(i,j)=im*(mo62i[i])[j]+(mo62r[i])[j];
    if(i<j){O62(i,j)=(mo62i[j])[i];}else{O62(i,j)=(mo62r[i])[j];}
   }}
   
   O64 = Matrix(1,dimj,1,dimj); 
   OO64 = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {OO64(i,j)=im*(mo64i[i])[j]+(mo64r[i])[j];
    if(i<j){O64(i,j)=(mo64i[j])[i];}else{O64(i,j)=(mo64r[i])[j];}
   }}
   
   O66 = Matrix(1,dimj,1,dimj); 
   OO66 = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {OO66(i,j)=im*(mo66i[i])[j]+(mo66r[i])[j];
    if(i<j){O66(i,j)=(mo66i[j])[i];}else{O66(i,j)=(mo66r[i])[j];}
   }}

   
   Hcf= Matrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {if(i<j){Hcf(i,j)=(hcfi[j])[i];}else{Hcf(i,j)=(hcfr[i])[j];}
   }
   
    for (i=1;i<=dj;++i)
     {delete[]Jxr[i];delete[]Jxi[i];
      delete[]Jyr[i];delete[]Jyi[i];
      delete[]Jzr[i];delete[]Jzi[i];
      delete[]hcfr[i];delete[]hcfi[i];
   delete[]mo20r[i];delete[]mo20i[i];
   delete[]mo22r[i];delete[]mo22i[i];
   delete[]mo40r[i];delete[]mo40i[i];
   delete[]mo42r[i];delete[]mo42i[i];
   delete[]mo44r[i];delete[]mo44i[i];
   delete[]mo60r[i];delete[]mo60i[i];
   delete[]mo62r[i];delete[]mo62i[i];
   delete[]mo64r[i];delete[]mo64i[i];
   delete[]mo66r[i];delete[]mo66i[i];}
     
     delete[]Jxr;delete[]Jxi;
     delete[]Jyr;delete[]Jyi;
     delete[]Jzr;delete[]Jzi;
     delete[]hcfr;delete[]hcfi;
   delete []mo20r;delete []mo20i;
   delete []mo22r;delete []mo22i;
   delete []mo40r;delete []mo40i;
   delete []mo42r;delete []mo42i;
   delete []mo44r;delete []mo44i;
   delete []mo60r;delete []mo60i;
   delete []mo62r;delete []mo62i;
   delete []mo64r;delete []mo64i;
   delete []mo66r;delete []mo66i;

//ATTENTION FOR NDCU2 the AXES xyz are parallel to cab
Matrix dummy(1,dimj,1,dimj);
dummy=Jb;Jb=Jc;Jc=Ja;Ja=dummy;
ComplexMatrix dummyc(1,dimj,1,dimj);
dummyc=Jbb;Jbb=Jcc;Jcc=Jaa;Jaa=dummyc;


}

module_ionpars iops("#ATTENTION in module cfield.so the AXES xyz are parallel to cab\n#The higher order interactions are described by the  PKQ Operators defined in cfield:\n#O20(c) .... Jd\n#O22(c) .... Je\n#O40(c) .... Jf\n#O42(c) .... Jg\n#O44(c) .... Jh\n#O60(c) .... Ji\n#O62(c) .... Jj\n#O64(c) .... Jk\n#O66(c) .... Jl\n");  // get 1ion parameters - operator matrices
#ifdef __declspec
extern "C" __declspec(dllexport) void mcalc(Vector & J,double & T, Vector & gjmbH,double * g_J, Vector & ABC,char ** sipffile,
                      double & lnZ,double & U)
#else
extern "C" void mcalc(Vector & J,double & T, Vector & gjmbH,double * g_J, Vector & ABC,char ** sipffile,
                      double & lnZ,double & U)
#endif  

{//ABC not used !!!
    /*on input
    T		temperature[K]
    gJmbH	vector of effective field [meV]
    gJ          Lande factor
    ABC         single ion parameter values (A, B, C corresponding to <+|Ja|->,<-|Jb|->,<+|Jc|->/i
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/
    UNUSED_PARAMETER(g_J);
    UNUSED_PARAMETER(ABC);
    UNUSED_PARAMETER(sipffile);

// check dimensions of vector
if(J.Hi()>12||gjmbH.Hi()>12)
   {fprintf(stderr,"Error loadable module cfield.so: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

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
   int dj;
   dj=iops.Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
    
   Ham=iops.Hcf-gjmbH(1)*iops.Ja-gjmbH(2)*iops.Jb-gjmbH(3)*iops.Jc;
   if(J.Hi()>=4){Ham-=gjmbH(4)*iops.O20;
   if(J.Hi()>=5){Ham-=gjmbH(5)*iops.O22;
   if(J.Hi()>=6){Ham-=gjmbH(6)*iops.O40;
   if(J.Hi()>=7){Ham-=gjmbH(7)*iops.O42;
   if(J.Hi()>=8){Ham-=gjmbH(8)*iops.O44;
   if(J.Hi()>=9){Ham-=gjmbH(9)*iops.O60;
   if(J.Hi()>=10){Ham-=gjmbH(10)*iops.O62;
   if(J.Hi()>=11){Ham-=gjmbH(11)*iops.O64;
   if(J.Hi()>=12){Ham-=gjmbH(12)*iops.O66;
   }}}}}}}}}
/*   int i1,j1; //printout matrix
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) printf ("%4.6g ",iops.O20(i1,j1));
    printf ("\n");
    }*/
      
    
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=0;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);
   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);
     double x,y;int i;
     x=Min(En);
     for (i=1;i<=dj;++i)
     {if ((y=(En(i)-x)/K_B/T)<700) wn[i]=exp(-y); 
      else wn[i]=0.0;
    //  printf("%g\n",En(i));
      }
     double Zs;
     Zs=Sum(wn);wn/=Zs;  
     lnZ=log(Zs)-x/K_B/T;
   // calculate U
     U=En*wn;
   // calculate Ja,Jb,Jc
     ComplexMatrix z(1,dj,1,dj);
     ComplexMatrix za(1,dj,1,dj);
     ComplexMatrix zb(1,dj,1,dj);
     ComplexMatrix zc(1,dj,1,dj);
     z=ComplexMatrix(zr,zi);
     
     za=iops.Jaa*z;
     zb=iops.Jbb*z;
     zc=iops.Jcc*z;

    
     J=0;
//    ComplexVector ddd;
    for (i=1;i<=dj;++i)
    {
     J[1]+=wn(i)*real(z.Column(i)*za.Column(i));
     J[2]+=wn(i)*real(z.Column(i)*zb.Column(i));
     J[3]+=wn(i)*real(z.Column(i)*zc.Column(i));
    }
     
    if(J.Hi()>=4){ComplexMatrix zo20(1,dj,1,dj);zo20=iops.OO20*z;for (i=1;i<=dj;++i) J[4]+=wn(i)*real(z.Column(i)*zo20.Column(i));
    if(J.Hi()>=5){ComplexMatrix zo22(1,dj,1,dj);zo22=iops.OO22*z;for (i=1;i<=dj;++i) J[5]+=wn(i)*real(z.Column(i)*zo22.Column(i));
    if(J.Hi()>=6){ComplexMatrix zo40(1,dj,1,dj);zo40=iops.OO40*z;for (i=1;i<=dj;++i) J[6]+=wn(i)*real(z.Column(i)*zo40.Column(i));
    if(J.Hi()>=7){ComplexMatrix zo42(1,dj,1,dj);zo42=iops.OO42*z;for (i=1;i<=dj;++i) J[7]+=wn(i)*real(z.Column(i)*zo42.Column(i));
    if(J.Hi()>=8){ComplexMatrix zo44(1,dj,1,dj);zo44=iops.OO44*z;for (i=1;i<=dj;++i) J[8]+=wn(i)*real(z.Column(i)*zo44.Column(i));
    if(J.Hi()>=9){ComplexMatrix zo60(1,dj,1,dj);zo60=iops.OO60*z;for (i=1;i<=dj;++i) J[9]+=wn(i)*real(z.Column(i)*zo60.Column(i));
    if(J.Hi()>=10){ComplexMatrix zo62(1,dj,1,dj);zo62=iops.OO62*z;for (i=1;i<=dj;++i) J[10]+=wn(i)*real(z.Column(i)*zo62.Column(i));
    if(J.Hi()>=11){ComplexMatrix zo64(1,dj,1,dj);zo64=iops.OO64*z;for (i=1;i<=dj;++i) J[11]+=wn(i)*real(z.Column(i)*zo64.Column(i));
    if(J.Hi()>=12){ComplexMatrix zo66(1,dj,1,dj);zo66=iops.OO66*z;for (i=1;i<=dj;++i) J[12]+=wn(i)*real(z.Column(i)*zo66.Column(i));
       }}}}}}}}}
       
       
         
  
return;
}
/**************************************************************************/

/**************************************************************************/
// for mcdisp this routine is needed
#ifdef __declspec
extern "C" __declspec(dllexport) int du1calc(int & tn,double & T,Vector & gjmbH,double * gJ,Vector & ABC, char ** sipffile,
                       ComplexVector & u1,float & delta)
#else
extern "C" int du1calc(int & tn,double & T,Vector & gjmbH,double * gJ,Vector & ABC, char ** sipffile,
                       ComplexVector & u1,float & delta)
#endif
{//ABC not used !!!
    /*on input
    tn          transitionnumber
    T		temperature[K]
    gJmbH	vector of effective field [meV]
    gJ          Lande factor
  on output    
    mat         transition element matrix
    delta       energy of transition
*/
// check dimensions of vector
if(gjmbH.Hi()>12)
   {fprintf(stderr,"Error loadable module cfield.so: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

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
   static Vector J(1,gjmbH.Hi());
   double lnz,u;
   J=0;
   if (T>0){ mcalc(J,T,gjmbH,gJ,ABC,sipffile,lnz,u);} else {T=-T;}
   // setup hamiltonian
   int dj;
   dj=iops.Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
    
   Ham=iops.Hcf-gjmbH(1)*iops.Ja-gjmbH(2)*iops.Jb-gjmbH(3)*iops.Jc;
   if(gjmbH.Hi()>=4){Ham-=gjmbH(4)*iops.O20;
   if(gjmbH.Hi()>=5){Ham-=gjmbH(5)*iops.O22;
   if(gjmbH.Hi()>=6){Ham-=gjmbH(6)*iops.O40;
   if(gjmbH.Hi()>=7){Ham-=gjmbH(7)*iops.O42;
   if(gjmbH.Hi()>=8){Ham-=gjmbH(8)*iops.O44;
   if(gjmbH.Hi()>=9){Ham-=gjmbH(9)*iops.O60;
   if(gjmbH.Hi()>=10){Ham-=gjmbH(10)*iops.O62;
   if(gjmbH.Hi()>=11){Ham-=gjmbH(11)*iops.O64;
   if(gjmbH.Hi()>=12){Ham-=gjmbH(12)*iops.O66;
   }}}}}}}}}

/*   int i1,j1; //printout matrix
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) printf ("%4.6g ",iops.O20(i1,j1));
    printf ("\n");
    }*/
      
int pr=1;
if(tn<0) {tn=-tn;pr=0;}
    
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=1;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);
   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);double Z;
     double x,y;int i,j=0,k,l;
     x=Min(En);
     for (i=1;i<=dj;++i)
     {if ((y=(En(i)-x)/K_B/T)<700) wn[i]=exp(-y); 
      else wn[i]=0.0;
//      printf("%g\n",En(i));
      }
     Z=Sum(wn);wn/=Z;  
     Z*=exp(-x/K_B/T);
   // calculate Ja,Jb,Jc
     ComplexMatrix z(1,dj,1,dj);
     ComplexMatrix * zp[gjmbH.Hi()+1];
     for(l=1;l<=gjmbH.Hi();++l)
      {zp[l]= new ComplexMatrix(1,dj,1,dj);}
     z=ComplexMatrix(zr,zi);

     (*zp[1])=iops.Jaa*z;
     (*zp[2])=iops.Jbb*z;
     (*zp[3])=iops.Jcc*z;
     if(gjmbH.Hi()>=4){(*zp[4])=iops.OO20*z;}
     if(gjmbH.Hi()>=5){(*zp[5])=iops.OO22*z;}
     if(gjmbH.Hi()>=6){(*zp[6])=iops.OO40*z;}
     if(gjmbH.Hi()>=7){(*zp[7])=iops.OO42*z;}
     if(gjmbH.Hi()>=8){(*zp[8])=iops.OO44*z;}
     if(gjmbH.Hi()>=9){(*zp[9])=iops.OO60*z;}
     if(gjmbH.Hi()>=10){(*zp[10])=iops.OO62*z;}
     if(gjmbH.Hi()>=11){(*zp[11])=iops.OO64*z;}
     if(gjmbH.Hi()>=12){(*zp[12])=iops.OO66*z;}

// calculate u1 and delta for transition number tn
// 1. get i and j from tn
k=0;
for(i=1;i<=dj;++i){for(j=i;j<=dj;++j)
{++k;if(k==tn)break;
}if(k==tn)break;}

// 2. set delta
delta=En(j)-En(i);

if (delta<-0.000001){fprintf(stderr,"ERROR module cfield.so - dmcalc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(j==i)delta=-SMALL; //if transition within the same level: take negative delta !!- this is needed in routine intcalc


// 3. set mat
for(l=1;l<=gjmbH.Hi();++l)
{if(i==j){//take into account thermal expectation values <Jl>
          u1(l)=((*zp[l]).Column(j)*z.Column(i))-J(l);}
 else    {u1(l)=(*zp[l]).Column(j)*z.Column(i);}}
           // ... in complex vector scalar product a*b is defined as: a.conj(b) !!! (see cvector.cc)


if (delta/K_B/T>0.000001)
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

//clean up memory
     for(l=1;l<=gjmbH.Hi();++l)
      {delete zp[l];}


     
return (int)((iops.J+1)*(2*iops.J+1)); // return number of all transitions
//return (int)(2*iops.J); // only exc from groundstate are counted
}

// this is called directly after loading it into memory from dlopen
void _init(void)
{  fprintf(stdout,"cfield.so: is loaded\n");}

// called just before removing from memory
void _fini(void)
{  fprintf(stdout,"cfield.so: is removed\n");}
