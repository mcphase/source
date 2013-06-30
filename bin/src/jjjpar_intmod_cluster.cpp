void jjjpar::cluster_ini_Imat() // to be called on initializing the cluster module
{// vector of matrix pointers which index the I1,I2,...In matrices of the cluster
 ComplexMatrix * Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+3*(*clusterpars).nofatoms+1];
 // initialize these matrices
 dim=1; Vector Hxc(1,(*clusterpars).nofcomponents);Vector Hext(1,3);
 dnn= new int [(*clusterpars).nofatoms+1];
 // determine dimension of H matrix
 for (int n=1;n<=(*clusterpars).nofatoms;++n)
 {dnn[n]=(*(*clusterpars).jjj[n]).opmat(1,Hxc,Hext).Rhi();
  dim*=dnn[n];
 }
 printf("# cluster module - number of matrices: %i  dimension=%i  \n",(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+3*(*clusterpars).nofatoms,dim);
 // initialize matrices
 Iaa[0]=new ComplexMatrix(1,dim,1,dim);(*Iaa[0])=1;
 for(int a=1;a<=(*clusterpars).nofatoms;++a)
 for(int i=1;i<=(*clusterpars).nofcomponents;++i)
  { Iaa[(a-1)*(*clusterpars).nofcomponents+i]=new ComplexMatrix(1,dim,1,dim);
    (*Iaa[(a-1)*(*clusterpars).nofcomponents+i])=0;
    Matrix Jai((*(*clusterpars).jjj[a]).opmat(i,Hxc,Hext));
    // myPrintMatrix(stdout,Jai);printf("\n");

  //1. determine dimensions
  int da=dnn[a];
  int dx=1,dz=1;
  for(int m=1  ;m<a;++m)dx*=dnn[m];
  for(int m=a+1;m<=(*clusterpars).nofatoms;
                    ++m)dz*=dnn[m];
  int mx=0,ma=1;
  if(a>1)mx=1;
  for(int m=a;  m<=(*clusterpars).nofatoms;++m)mx*=dnn[m];
  for(int m=a+1;m<=(*clusterpars).nofatoms;++m)ma*=dnn[m];

  //2. insertion loop
  for(int ix=0;ix<=dx-1;++ix)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=da-1;++ai)
  for(int bi=0;bi<=da-1;++bi)
  {int r=ix*mx+ai*ma+iz+1;
   int s=ix*mx+bi*ma+iz+1;
 //printf("hello %i %i %i %i\n",r,s,ai,bi);
  // here we should fill the matrices with values corresponding to the
  // I1 I2 of the individual atoms of the cluster
   if(r<s)        {(*Iaa[(a-1)*(*clusterpars).nofcomponents+i])(r,s)+=complex<double>(Jai(bi+1,ai+1),-Jai(ai+1,bi+1));}
   else if (r==s) {(*Iaa[(a-1)*(*clusterpars).nofcomponents+i])(r,s)+=complex<double>(Jai(ai+1,bi+1),0.0);}   
              else{(*Iaa[(a-1)*(*clusterpars).nofcomponents+i])(r,s)+=complex<double>(Jai(ai+1,bi+1),Jai(bi+1,ai+1));}   
  }
 }

      
  // her we should make possible the perlparsing the cluster sipf file
  // in order to flexibly change the Ia operator sequence ...
  double *numbers[1]; numbers[0]=&gJ;
  char numnam[]="gJ\0";
  char * numbernames[2];numbernames[0]=numnam;numbernames[1]=NULL;

  char * strings[1];char * stringnames[1];
  stringnames[0]=NULL;

 char * operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+3*(*clusterpars).nofatoms+2];
 // char opnam[]="one\0";   
  operatornames[0]=new char[5];strcpy(operatornames[0],"one\0");
 for(int a=1;a<=(*clusterpars).nofatoms;++a)
 for(int i=1;i<=(*clusterpars).nofcomponents;++i)
 {operatornames[(a-1)*(*clusterpars).nofcomponents+i]=new char[10];
  sprintf(operatornames[(a-1)*(*clusterpars).nofcomponents+i],"I%i_%i",a,i);
 }     
 // now initialize the interaction operators
 for(int i=1;i<=nofcomponents;++i)
 {operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+i]=new char[5];
  sprintf(operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+i],"I%i",i);
  Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+i]=new ComplexMatrix(1,dim,1,dim);                     
 }
 // now initialize the cluster total magnetic moment operators
 for(int i=1;i<=3;++i) 
 {operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+i]=new char[5];
  sprintf(operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+i],"M%i",i);
  Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+ i]=new ComplexMatrix(1,dim,1,dim);                     
 }

 // now initialize the cluster individual magnetic moment operators 
 //(for neutron form factor in dipole approx - function MQ dMQ1)
 for(int a=1;a<=(*clusterpars).nofatoms;++a)for(int n=1;n<=3;++n)
 {int index=(a-1)*3+n;
  operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+index]=new char[5];
  sprintf(operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+index],"M%i_%i",a,n);
  Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+index]=new ComplexMatrix(1,dim,1,dim);                     
 } 
  operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+3*(*clusterpars).nofatoms+1]=NULL;

 if(perlparse(sipffilename,numbers,numbernames,strings,stringnames,Iaa,operatornames)==false)
   {printf("Error perl parsing sipf file %s\n",sipffilename);exit(EXIT_FAILURE);}
for(int i=0;i<=(*clusterpars).nofatoms*(*clusterpars).nofcomponents;++i)
  {delete Iaa[i];delete operatornames[i];}

// here we fill the interaction operator matrices with values
Ia= new Matrix * [nofcomponents+1];
for(int n = 1;n<=nofcomponents;++n){Ia[n]=new Matrix(1,dim,1,dim);
for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){
   if(i<j){(*Ia[n])(i,j)=imag((*Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+n])(j,i));}
      else{(*Ia[n])(i,j)=real((*Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+n])(i,j));}
}
//(*Ia[i])=(*Iaa[i]);
delete           Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+n];
delete operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+n]; 
}

// here we should initialize set the operator matrices cluster_M  ... total magnetic moment
// operator storage and individual magnetic moments: this should then be used
// in the mcalc, dm1calc, MQ and dMQ1 functions of the cluster module ...
// to be done !!!

cluster_M= new Matrix * [3+3*(*clusterpars).nofatoms+1];
// first do the total moment
for(int n = 1;n<=3;++n){cluster_M[n]=new Matrix(1,dim,1,dim);
for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){
   if(i<j){(*cluster_M[n])(i,j)=imag((*Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+n])(j,i));}
      else{(*cluster_M[n])(i,j)=real((*Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+n])(i,j));}
}
delete           Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+n];
delete operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+n]; 
}
// now do the individual atom magnetic moments
 for(int a=1;a<=(*clusterpars).nofatoms;++a)for(int n=1;n<=3;++n)
{int index=(a-1)*3+n; // a .... atom index  n ... xyz components of magnetic moment
 int index_M=3*a+n;
 cluster_M[index_M]=new Matrix(1,dim,1,dim);
 for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){
   if(i<j){(*cluster_M[index_M])(i,j)=imag((*Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+index])(j,i));}
      else{(*cluster_M[index_M])(i,j)=real((*Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+index])(i,j));}
}
delete           Iaa[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+index];
delete operatornames[(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+index]; 
}

// for(int i=0;i<=(*clusterpars).nofatoms*(*clusterpars).nofcomponents+nofcomponents+3+3*(*clusterpars).nofatoms;++i)
//  {delete Iaa[i];delete operatornames[i];}
printf("#module cluster initialized\n");
}
//------------------------------------------------------------------------------------------------
//routine Icalc for cluster
//------------------------------------------------------------------------------------------------
void jjjpar::cluster_Icalc_mcalc_Micalc (int code,Vector & Jret,double & T, Vector &  Hxc,Vector & Hext, double & lnZ, double & U)
{ /*on input
   code         defining, what should be calculated
         1....   Ia (interaction operators for Icalc)
         2....    m (total magnetic moment for mcalc), or
         3....    Mi (individual magnetic moments (for MQ)
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/


Vector En(1,dim);
Matrix zr(1,dim,1,dim);
Matrix zc(1,dim,1,dim);
cluster_calcH_and_diagonalize(En,zr,zc,Hxc,Hext);

// calculate Z and wn (occupation probability)
     Vector wn(1,dim);double Zs;
     double x,y;
     x=Min(En);
    if(T>0)
     {
     for (int i=1;i<=dim;++i)
     {if ((y=(En(i)-x)/KB/T)<700) wn[i]=exp(-y);
      else wn[i]=0.0;
//      printf("%4.4g\n",En(i));
      }
     Zs=Sum(wn);wn/=Zs;
     lnZ=log(Zs)-x/KB/T;
     }
     else
     { printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                         printf ("Number   Energy     Excitation Energy\n");
     for (int i=1;i<=dim;++i) printf ("%i    %4.4g meV   %4.4g meV\n",i,En(i),En(i)-x);
     char instr[MAXNOFCHARINLINE];
     for (int i=1;i<=dim;++i)
      {printf("eigenstate %i: %4.4g meV %4.4g meV  - please enter probability w(%i):",i,En(i),En(i)-x,i);
       fgets(instr, MAXNOFCHARINLINE, stdin);

       wn(i)=strtod(instr,NULL);
      }
       Zs=Sum(wn);wn/=Zs;

       lnZ=log(Zs);
                         printf ("\n\nNumber   Energy     Excitation Energy   Probability\n");
     for (int i=1;i<=dim;++i) printf ("%i    %4.4g meV   %4.4g meV %4.4g  \n",i,En(i),En(i)-x,wn(i));
     }

   // calculate U
     U=En*wn;

 switch(code)

{case 1:for(int a=1;a<=Hxc.Hi();++a)
        {// calculate expecation Value of J
         //printf("Matrix of Operator");
         // myPrintMatrix(stdout,Ja);
         // determine expectation value
         Jret(a)=0;
         for(int i=1;i<=dim&&wn[i]>0.00001;++i)
         {Jret(a)+=wn[i]*aMb_real((*Ia[a]),zr,zc,i,i); 
          // if(fabs(aMb_imag((*Ia[a]),zr,zc,i,i))>SMALL){fprintf(stderr,"ERROR module cluster - Icalc: expectation value imaginary\n");exit(EXIT_FAILURE);}
         }
        }break;
 case 2:for(int n=1;n<=3;++n)
        {// calculate expecation Value of m
         Jret(n)=0;
         for(int i=1;i<=dim&&wn[i]>0.00001;++i)
         {Jret(n)+=wn[i]*aMb_real((*cluster_M[n]),zr,zc,i,i); 
         }
        }break;
 case 3:for(int a=1;a<=(*clusterpars).nofatoms;++a)for(int n=1;n<=3;++n)
        {int index_M=a*3+n; // a .... atom index  n ... xyz components of magnetic moment
         int index=(a-1)*3+n; // calculate expecation Value of m
         Jret(index)=0;
         for(int i=1;i<=dim&&wn[i]>0.00001;++i)
         {Jret(index)+=wn[i]*aMb_real((*cluster_M[index_M]),zr,zc,i,i); 
         }
        }break;
}

//  printf ("Ha=%g Hb=%g Hc=%g Ja=%g Jb=%g Jc=%g \n", 
//     gjmbH[1]/MU_B/gjJ, gjmbH[2]/MU_B/gjJ, gjmbH[3]/MU_B/gjJ, J[1], J[2], J[3]);
}
//------------------------------------------------------------------------------------------------
//routine ducalc for cluster
//------------------------------------------------------------------------------------------------

int jjjpar::cluster_dm(int code,int & tn,double & T,ComplexVector & u1,float & delta,ComplexMatrix & ests)
{ 
  /*on input
   code         defining, what should be calculated
         1....   Ia (interaction operators for du1calc)
         2....    m (total magnetic moment for dm1calc), or
         3....    Mi (individual magnetic moments (for dMQ1)
    transitionnumber ... number of transition to be computed - meaningless for kramers doublet, because there is only 1 transition
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splitting of kramers doublet [meV]
    u1(i)	<-|(Ji-<Ji>)|+> sqrt(tanh(delta/2kT))
*/
  int pr=0,subtractexpvalue=1;if(T<0){subtractexpvalue=0;T=-T;}
  if (tn<0) {pr=1;tn*=-1;}
   double ninit=u1[1].real();
   double pinit=u1[1].imag();
 
Vector En(1,dim);
Vector wn(1,dim);
Matrix zr(1,dim,1,dim);
Matrix zc(1,dim,1,dim);
 for(int i=1;i<=dim;++i){En(i)=real(ests(0,i));wn(i)=imag(ests(0,i));}

   for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){zr(i,j)=real(ests(i,j));zc(i,j)=imag(ests(i,j));}

// calculate mat and delta for transition number tn
// 1. get i and j from tn
int k=0,ii=1,jj=1;
for(ii=1;ii<=dim;++ii){for(jj=ii;jj<=dim;++jj){++k;if(k==tn)break;}if(k==tn)break;}
// 2. set delta
delta=En(jj)-En(ii);
if (delta<-0.000001){fprintf(stderr,"ERROR module cluster - du1calc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(jj==ii)delta=-SMALL_QUASIELASTIC_ENERGY; //if transition within the same level: take negative delta !!- this is needed in routine intcalc

int maxn;
 switch(code)
{case 1: maxn=u1.Hi();break;
 case 2: maxn=3; break;
 case 3: maxn=(*clusterpars).nofatoms*3;break;
}

Vector Jret(1,maxn);Jret=0;
ComplexVector iJj(1,maxn);
for(int a=1;a<=maxn;++a)
{// transition matrix element
 switch(code)
 {case 1:  iJj(a)=complex<double>(aMb_real((*Ia[a]),zr,zc,ii,jj),aMb_imag((*Ia[a]),zr,zc,ii,jj));break;
  case 2:  iJj(a)=complex<double>(aMb_real((*cluster_M[a]),zr,zc,ii,jj),aMb_imag((*cluster_M[a]),zr,zc,ii,jj)); break;
  case 3:  iJj(a)=complex<double>(aMb_real((*cluster_M[a+3]),zr,zc,ii,jj),aMb_imag((*cluster_M[a+3]),zr,zc,ii,jj));break;
 }

 // determine expectation value
 Jret(a)=0;
 if (subtractexpvalue==1&&ii==jj)
 { for(int i=1;i<=dim&&wn[i]>0.00001;++i)
  { switch(code)
   {case 1: Jret(a)+=wn[i]*aMb_real((*Ia[a]),zr,zc,i,i);break;
    case 2: Jret(a)+=wn[i]*aMb_real((*cluster_M[a]),zr,zc,i,i); break;
    case 3: Jret(a)+=wn[i]*aMb_real((*cluster_M[a+3]),zr,zc,i,i);break;
   }
  }
 }
}

// 3. set u1
for(int l=1;l<=maxn;++l)
{if(ii==jj){//take into account thermal expectation values <Jret>
          u1(l)=(iJj(l)-Jret(l));}
 else    {u1(l)=iJj(l);}
}

//printf("code %i: %g %g %g\n",code,Jret(1),Jret(2),Jret(3));
//myPrintMatrix(stdout,(*cluster_M[1]));
//myPrintMatrix(stdout,(*cluster_M[2]));
//myPrintMatrix(stdout,(*cluster_M[3]));



if (delta>SMALL_QUASIELASTIC_ENERGY)
   { if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",ii,jj,delta);
      printf(" |<%i|Ja|%i>|^2=%4.4g |<%i|Jb|%i>|^2=%4.4g |<%i|Jc|%i>|^2=%4.4g",ii,jj,abs(u1(1))*abs(u1(1)),ii,jj,abs(u1(2))*abs(u1(2)),ii,jj,abs(u1(3))*abs(u1(3)));
      printf(" n%i-n%i=%4.4g\n",ii,jj,wn(ii)-wn(jj));}
    u1*=sqrt(wn(ii)-wn(jj)); // occupation factor
     }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
     if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",ii,jj,delta);
      printf(" |<%i|Ja-<Ja>|%i>|^2=%4.4g |<%i|Jb-<Jb>|%i>|^2=%4.4g |<%i|Jc-<Jc>|%i>|^2=%4.4g",ii,jj,abs(u1(1))*abs(u1(1)),ii,jj,abs(u1(2))*abs(u1(2)),ii,jj,abs(u1(3))*abs(u1(3)));
      printf(" n%i=%4.4g\n",ii,wn(ii));}
    u1*=sqrt(wn(ii)/KB/T);
   }


if (pr==1) printf ("delta=%4.6g meV\n",delta);

// return number of all transitions
   double n=ninit;
   if (n>dim)n=dim;
   double zsum=0,zii,x;
   int noft=0;
   if(T>0)for(int i=1;(i<=n)&((((x=(real(est(0,i))-real(est(0,1)))/KB/T)<200)? zii=exp(-x):zii=0)>=(pinit*zsum));++i)
   {noft+=dim-i+1;zsum+=zii;}
   
   if(T<0)for(int i=1;i<=n;++i)
   {noft+=dim-i+1;}
//printf("nt=%i",noft);
 return noft;



}

//------------------------------------------------------------------------------------------------
//routine estates for cluster
//------------------------------------------------------------------------------------------------
void jjjpar::cluster_est(ComplexMatrix * eigenstates,Vector &Hxc,Vector &Hext,double & T)
{ /*on input
    gJmbH	vector of effective field [meV]
      on output
    Matrix containing the eigenvalues and eigenfunctions of the crystalfield problem
    eigenvalues ares stored as real part of row zero
    boltzmann population numbers are stored as imaginary part of row zero
*/
 fprintf(stderr,"# calculating eigenstates of cluster ...");
(*eigenstates) = ComplexMatrix(0,dim,1,dim);
 Vector En(1,dim);
 Matrix zr(1,dim,1,dim);
 Matrix zc(1,dim,1,dim);
 cluster_calcH_and_diagonalize(En,zr,zc,Hxc,Hext);
     Vector wn(1,dim);double Zs;
     double x,y;
     x=Min(En);
     for (int i=1;i<=dim;++i)
     {if ((y=(En(i)-x)/KB/T)<700) wn[i]=exp(-y);
      else wn[i]=0.0;
      }
     Zs=Sum(wn);wn/=Zs;
 for(int i=1;i<=dim;++i)(*eigenstates)(0,i)=complex<double>(En(i),wn(i));

 for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j)
 {(*eigenstates)(i,j)=complex<double>(zr(i,j),zc(i,j));
 }
 fprintf(stderr," ... done\n");
}

void jjjpar::cluster_calcH_and_diagonalize(Vector & En,Matrix & zr, Matrix & zc,Vector & Hxc,Vector & Hext)
{Matrix H(1,dim,1,dim);Vector ZeroHxc(1,1);ZeroHxc=0;
 H=0; // initialize to zero

// fill H matrix with sum over Hi of individual spins
for (int i=1;i<=(*clusterpars).nofatoms;++i)
{Matrix Hi((*(*clusterpars).jjj[i]).opmat(0,ZeroHxc,Hext)); // here we need ZeroHxc because
                                                            // exchange operators are set by user (Ia)
                                                            // and not in the single ion module
  //1. determine dimensions
  int di=dnn[i];
  int dx=1,dz=1;
  for(int m=1  ;m<i;++m)dx*=dnn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;
                    ++m)dz*=dnn[m];
  int mx=0,mi=1;
  if(i>1)mx=1;
  for(int m=i;  m<=(*clusterpars).nofatoms;++m)mx*=dnn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;++m)mi*=dnn[m];

  //2. insertion loop
  for(int ix=0;ix<=dx-1;++ix)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=di-1;++ai)
  for(int bi=0;bi<=di-1;++bi)
  {int r=ix*mx+ai*mi+iz+1;
   int s=ix*mx+bi*mi+iz+1;
//   printf("hello %i %i %i %i\n",r,s,ai,bi);
   H(r,s)+=Hi(ai+1,bi+1);
  }
}

// myPrintMatrix(stdout,H);
//printf("now enter interactions\n");
// put interactions in H matrix
for (int n=1;n<=(*clusterpars).nofatoms;++n)
for (int nn=1;nn<=(*(*clusterpars).jjj[n]).paranz;++nn)
{// determine interaction operator of atom n with neighbour nn
 // get i and j from n and sublattice[nn]
 int i=n;int j=(*(*clusterpars).jjj[n]).sublattice[nn];
 int ilj=1;
 if(i>j){i=j;j=n;ilj=0;} // if n>nn exchange i and j ... so that we always have i < j
 Matrix SinS(-0.5*(*(*clusterpars).jjj[n]).jij[nn](1,1)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(1,ZeroHxc,Hext),
                          (*(*clusterpars).jjj[j]).opmat(1,ZeroHxc,Hext)
                         )
             );

 for(int a=1;a<=(*(*clusterpars).jjj[i]).nofcomponents;++a)
 for(int b=1;b<=(*(*clusterpars).jjj[j]).nofcomponents;++b)
 {
  if(a==1&&b==1) break;
  if(ilj){
  SinS+=-0.5*(*(*clusterpars).jjj[n]).jij[nn](a,b)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(a,ZeroHxc,Hext),
                          (*(*clusterpars).jjj[j]).opmat(b,ZeroHxc,Hext)
                    );
          } else { // if n>nn then take transpose of exchange parameter jij !! (J(ij)=JT(ji))
  SinS+=-0.5*(*(*clusterpars).jjj[n]).jij[nn](b,a)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(a,ZeroHxc,Hext),
                          (*(*clusterpars).jjj[j]).opmat(b,ZeroHxc,Hext)
                    );
          }
 }
//printf("hello n=%i nn=%i i=%i j=%i ...\n",n,nn,i,j);  myPrintMatrix(stdout,SinS);
 // insert SinS into H
 //1. determine dimensions
  int di=dnn[i];
  int dj=dnn[j];
 int dx=1,dy=1,dz=1;
 for(int m=1  ;m<i;++m)dx*=dnn[m];
 for(int m=i+1;m<j;++m)dy*=dnn[m];
 for(int m=j+1;m<=(*clusterpars).nofatoms;
                   ++m)dz*=dnn[m];

 int mx=0,my=0,mi=1,mj=1;
 if(i>1)mx=1;
 if(i+1<=(*clusterpars).nofatoms&&i+1<j)my=1;
 
 for(int m=i;  m<=(*clusterpars).nofatoms;++m)mx*=dnn[m];
 for(int m=i+1;m<=(*clusterpars).nofatoms;++m)mi*=dnn[m];
 for(int m=j  ;m<=(*clusterpars).nofatoms;++m)my*=dnn[m];
 for(int m=j+1;m<=(*clusterpars).nofatoms;++m)mj*=dnn[m];

 //2. insertion loop
  for(int ix=0;ix<=dx-1;++ix)
  for(int iy=0;iy<=dy-1;++iy)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=di-1;++ai)
  for(int aj=0;aj<=dj-1;++aj)
  for(int bi=0;bi<=di-1;++bi)
  for(int bj=0;bj<=dj-1;++bj)
  {int r=ix*mx+ai*mi+iy*my+aj*mj+iz+1;
   int s=ix*mx+bi*mi+iy*my+bj*mj+iz+1;
   int k=ai*dj+aj+1;// ??
   int l=bi*dj+bj+1;// ??
   H(r,s)+=SinS(k,l);
//printf("hello %i %i %i %i %g %g\n",r,s,k,l,H(r,s),SinS(k,l));
  } 
}

// insert exchange field
for(int i =1;i<=Hxc.Hi();++i)H-=Hxc(i)*(*Ia[i]);

// diagonalize H
int sort=1;int maxiter=1000000;
// myPrintMatrix(stdout,H);
//printf("now diagonalise\n");
EigenSystemHermitean (H,En,zr,zc,sort,maxiter);
//printf("Eigenvector real\n");
// myPrintMatrix(stdout,zr);
//printf("Eigenvector imag\n"); 
// myPrintMatrix(stdout,zc);

}

/**************************************************************************/
//                   OBSERVABLES
/**************************************************************************/

void jjjpar::cluster_Micalc (Vector &mom,ComplexMatrix & ests)
{
Vector En(1,dim);
Vector wn(1,dim);
Matrix zr(1,dim,1,dim);
Matrix zc(1,dim,1,dim);
 for(int i=1;i<=dim;++i){En(i)=real(ests(0,i));wn(i)=imag(ests(0,i));}

   for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){zr(i,j)=real(ests(i,j));zc(i,j)=imag(ests(i,j));}

for(int a=1;a<=(*clusterpars).nofatoms;++a)for(int n=1;n<=3;++n)
        {int index_M=a*3+n; // a .... atom index  n ... xyz components of magnetic moment
         int index=(a-1)*3+n; // calculate expecation Value of m
         mom(index)=0;
         for(int i=1;i<=dim&&wn[i]>0.00001;++i)
         {mom(index)+=wn[i]*aMb_real((*cluster_M[index_M]),zr,zc,i,i); 
         }
        }
}