//------------------------------------------------------------------------------------------------
//routine mcalc for cluster
//------------------------------------------------------------------------------------------------
void jjjpar::cluster_mcalc (Vector & Jret,double & T, Vector & gjmbH, double & lnZ, double & U)
{ /*on input
    ABC(1...3)  A,M,Ci....saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    gJ		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/
int dim=1; int *d0= new int [(*clusterpars).nofatoms+1];
           int *dn= new int [(*clusterpars).nofatoms+1];
// determine dimension of H matrix
for (int n=1;n<=(*clusterpars).nofatoms;++n)
{d0[n]=dim;
 dn[n]=(*(*clusterpars).jjj[n]).opmat(0,gjmbH).Rhi();
 dim*=dn[n];
 if(gjmbH.Hi()!=(*(*clusterpars).jjj[n]).nofcomponents)
 {fprintf(stderr,"Error module cluster - mcalc: in current version dimensions of meanfield must match number of components of each atom in cluster, individual coupling between atoms in different clusters are not implemented yet\n");exit(EXIT_FAILURE);}
}

 Matrix H(1,dim,1,dim);
 H=0; // initialize to zero

// fill H matrix with sum over Hi of individual spins
for (int i=1;i<=(*clusterpars).nofatoms;++i)
{Matrix Hi((*(*clusterpars).jjj[i]).opmat(0,gjmbH));
  //1. determine dimensions
  int di=dn[i];
  int dx=1,dz=1;
  for(int m=1  ;m<i;++m)dx*=dn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;
                    ++m)dz*=dn[m];
  int mx=0,mi=1;
  if(i>1)mx=1;
  for(int m=i;  m<=(*clusterpars).nofatoms;++m)mx*=dn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;++m)mi*=dn[m];

  //2. insertion loop
  for(int ix=0;ix<=dx-1;++ix)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=di-1;++ai)
  for(int bi=0;bi<=di-1;++bi)
  {int r=ix*mx+ai*mi+iz+1;
   int s=ix*mx+bi*mi+iz+1;
 //  printf("hello %i %i %i %i\n",r,s,ai,bi);
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
 if(i>j){i=j;j=n;}
 Matrix SinS(-0.5*(*(*clusterpars).jjj[n]).jij[nn](1,1)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(1,gjmbH),
                          (*(*clusterpars).jjj[j]).opmat(1,gjmbH)
                         )
             );

 for(int a=1;a<=(*(*clusterpars).jjj[i]).nofcomponents;++a)
 for(int b=1;b<=(*(*clusterpars).jjj[j]).nofcomponents;++b)
 {
  if(a==1&&b==1) break;
  SinS+=-0.5*(*(*clusterpars).jjj[n]).jij[nn](a,b)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(a,gjmbH),
                          (*(*clusterpars).jjj[j]).opmat(b,gjmbH)
                    );
 }
 
 // insert SinS into H
 //1. determine dimensions
  int di=dn[i];
  int dj=dn[j];
 int dx=1,dy=1,dz=1;
 for(int m=1  ;m<i;++m)dx*=dn[m];
 for(int m=i+1;m<j;++m)dy*=dn[m];
 for(int m=j+1;m<=(*clusterpars).nofatoms;
                   ++m)dz*=dn[m];

 int mx=0,my=0,mi=1,mj=1;
 if(i>1)mx=1;
 if(i+1<=(*clusterpars).nofatoms&&i+1<j)my=1;
 
 for(int m=i;  m<=(*clusterpars).nofatoms;++m)mx*=dn[m];
 for(int m=i+1;m<=(*clusterpars).nofatoms;++m)mi*=dn[m];
 for(int m=j  ;m<=(*clusterpars).nofatoms;++m)my*=dn[m];
 for(int m=j+1;m<=(*clusterpars).nofatoms;++m)mj*=dn[m];

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
   int k=ai*di+aj+1;
   int l=bi*di+bj+1;
   H(r,s)+=SinS(k,l);
//printf("hello %i %i %i %i %g %g\n",r,s,k,l,H(r,s),SinS(k,l));
  } 
}


// diagonalize H
Vector En(1,dim);
Matrix zr(1,dim,1,dim);
Matrix zc(1,dim,1,dim);
int sort=1;int maxiter=1000000;
// myPrintMatrix(stdout,H);
//printf("now diagonalise\n");
EigenSystemHermitean (H,En,zr,zc,sort,maxiter);
// myPrintMatrix(stdout,zr);
// myPrintMatrix(stdout,zc);
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

Matrix Ja(1,dim,1,dim);
for(int a=1;a<=gjmbH.Hi();++a)
{// calculate expecation Value of J
 // determine Matrix Ja= sum_i Jai
 Ja=0;// reset Ja
 for(int i=1;i<=(*clusterpars).nofatoms;++i)
 {Matrix Jai((*(*clusterpars).jjj[i]).opmat(a,gjmbH));

// myPrintMatrix(stdout,Jai);

  //1. determine dimensions
  int di=dn[i];
  int dx=1,dz=1;
  for(int m=1  ;m<i;++m)dx*=dn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;
                    ++m)dz*=dn[m];
  int mx=0,mi=1;
  if(i>1)mx=1;
  for(int m=i;  m<=(*clusterpars).nofatoms;++m)mx*=dn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;++m)mi*=dn[m];

  //2. insertion loop
  for(int ix=0;ix<=dx-1;++ix)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=di-1;++ai)
  for(int bi=0;bi<=di-1;++bi)
  {int r=ix*mx+ai*mi+iz+1;
   int s=ix*mx+bi*mi+iz+1;
 //  printf("hello %i %i %i %i\n",r,s,ai,bi);
   Ja(r,s)+=Jai(ai+1,bi+1);
  }
 }
//printf("hello %i %i ");
// myPrintMatrix(stdout,Ja);
 // determine expectation value
 Jret[a]=0;
 for(int i=1;i<=dim&&wn[i]>0.00001;++i)
 {Jret[a]+=wn[i]*aMb_real(Ja,zr,zc,i,i); }
}

//  printf ("Ha=%g Hb=%g Hc=%g Ja=%g Jb=%g Jc=%g \n", 
//     gjmbH[1]/MU_B/gjJ, gjmbH[2]/MU_B/gjJ, gjmbH[3]/MU_B/gjJ, J[1], J[2], J[3]);
delete []d0;delete []dn;
}

int jjjpar::cluster_dm(int & tn,double & T,Vector & gjmbH,ComplexVector & u1,float & delta)
{ 
  /*on input
    transitionnumber ... number of transition to be computed - meaningless for kramers doublet, because there is only 1 transition
    ABC[i]	saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    gJ		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splitting of kramers doublet [meV]
    u1(i)	<-|(Ji-<Ji>)|+> sqrt(tanh(delta/2kT))
*/
  int pr=1,subtractexpvalue=1;if(T<0){subtractexpvalue=0;T=-T;}
  if (tn<0) {pr=0;tn*=-1;}
int dim=1; int *d0= new int [(*clusterpars).nofatoms+1];
           int *dn= new int [(*clusterpars).nofatoms+1];
// determine dimension of H matrix
for (int n=1;n<=(*clusterpars).nofatoms;++n)
{d0[n]=dim;
 dn[n]=(*(*clusterpars).jjj[n]).opmat(0,gjmbH).Rhi();
 dim*=dn[n];
 if(gjmbH.Hi()!=(*(*clusterpars).jjj[n]).nofcomponents)
 {fprintf(stderr,"Error module cluster - du1calc: in current version dimensions of meanfield must match number of components of each atom in cluster, individual coupling between atoms in different clusters are not implemented yet\n");exit(EXIT_FAILURE);}
}

 Matrix H(1,dim,1,dim);
 H=0; // initialize to zero
// fill H matrix with sum over Hi of individual spins
for (int i=1;i<=(*clusterpars).nofatoms;++i)
{Matrix Hi((*(*clusterpars).jjj[i]).opmat(0,gjmbH));
  //1. determine dimensions
  int di=dn[i];
  int dx=1,dz=1;
  for(int m=1  ;m<i;++m)dx*=dn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;
                    ++m)dz*=dn[m];
  int mx=0,mi=1;
  if(i>1)mx=1;
  for(int m=i;  m<=(*clusterpars).nofatoms;++m)mx*=dn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;++m)mi*=dn[m];

  //2. insertion loop
  for(int ix=0;ix<=dx-1;++ix)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=di-1;++ai)
  for(int bi=0;bi<=di-1;++bi)
  {int r=ix*mx+ai*mi+iz+1;
   int s=ix*mx+bi*mi+iz+1;
   H(r,s)+=Hi(ai+1,bi+1);
  }
}

// put interactions in H matrix
for (int n=1;n<=(*clusterpars).nofatoms;++n)
for (int nn=1;nn<=(*(*clusterpars).jjj[n]).paranz;++nn)
{// determine interaction operator of atom n with neighbour nn
 // get i and j from n and sublattice[nn]
 int i=n;int j=(*(*clusterpars).jjj[n]).sublattice[nn];
 if(i>j){i=j;j=n;}
 Matrix SinS(-0.5*(*(*clusterpars).jjj[n]).jij[nn](1,1)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(1,gjmbH),
                          (*(*clusterpars).jjj[j]).opmat(1,gjmbH)
                         )
             );

 for(int a=1;a<=(*(*clusterpars).jjj[i]).nofcomponents;++a)
 for(int b=1;b<=(*(*clusterpars).jjj[j]).nofcomponents;++b)
 {
  if(a==1&&b==1) break;
  SinS+=-0.5*(*(*clusterpars).jjj[n]).jij[nn](a,b)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(a,gjmbH),
                          (*(*clusterpars).jjj[j]).opmat(b,gjmbH)
                    );
 }

 //1. determine dimensions
  int di=dn[i];
  int dj=dn[j];
 int dx=1,dy=1,dz=1;
 for(int m=1  ;m<i;++m)dx*=dn[m];
 for(int m=i+1;m<j;++m)dy*=dn[m];
 for(int m=j+1;m<=(*clusterpars).nofatoms;
                   ++m)dz*=dn[m];

 int mx=0,my=0,mi=1,mj=1;
 if(i>1)mx=1;
 if(i+1<=(*clusterpars).nofatoms&&i+1<j)my=1;

 for(int m=i;  m<=(*clusterpars).nofatoms;++m)mx*=dn[m];
 for(int m=i+1;m<=(*clusterpars).nofatoms;++m)mi*=dn[m];
 for(int m=j  ;m<=(*clusterpars).nofatoms;++m)my*=dn[m];
 for(int m=j+1;m<=(*clusterpars).nofatoms;++m)mj*=dn[m];

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
   int k=ai*di+aj+1;
   int l=bi*di+bj+1;
   H(r,s)+=SinS(k,l);
  }

}


// diagonalize H
Vector En(1,dim);
Matrix zr(1,dim,1,dim);
Matrix zc(1,dim,1,dim);
int sort=1;int maxiter=1000000;
EigenSystemHermitean (H,En,zr,zc,sort,maxiter);

// calculate mat and delta for transition number tn
// 1. get i and j from tn
int k=0,ii,jj;
for(ii=1;ii<=dim;++ii){for(jj=ii;jj<=dim;++jj)
{++k;if(k==tn)break;
}if(k==tn)break;}

// 2. set delta
delta=En(jj)-En(ii);

if (delta<-0.000001){fprintf(stderr,"ERROR module cluster - du1calc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(jj==ii)delta=-SMALL; //if transition within the same level: take negative delta !!- this is needed in routine intcalc

// calculate Z and wn (occupation probability)
     Vector wn(1,dim);double Zs;
     double x,y;
     x=Min(En);
     for (int i=1;i<=dim;++i)
     {if ((y=(En(i)-x)/KB/T)<700) wn[i]=exp(-y);
      else wn[i]=0.0;
//      printf("%4.4g\n",En(i));
      }
     Zs=Sum(wn);wn/=Zs;

Vector Jret(1,gjmbH.Hi());Jret=0;
ComplexVector iJj(1,gjmbH.Hi());

Matrix Ja(1,dim,1,dim);
for(int a=1;a<=gjmbH.Hi();++a)
{// calculate expecation Value of J
 // determine Matrix Ja= sum_i Jai
 Ja=0;// reset Ja
 for(int i=1;i<=(*clusterpars).nofatoms;++i)
 {Matrix Jai((*(*clusterpars).jjj[i]).opmat(a,gjmbH));
  //1. determine dimensions
  int di=dn[i];
  int dx=1,dz=1;
  for(int m=1  ;m<i;++m)dx*=dn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;
                    ++m)dz*=dn[m];
  int mx=0,mi=1;
  if(i>1)mx=1;
  for(int m=i;  m<=(*clusterpars).nofatoms;++m)mx*=dn[m];
  for(int m=i+1;m<=(*clusterpars).nofatoms;++m)mi*=dn[m];

  //2. insertion loop
  for(int ix=0;ix<=dx-1;++ix)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=di-1;++ai)
  for(int bi=0;bi<=di-1;++bi)
  {int r=ix*mx+ai*mi+iz+1;
   int s=ix*mx+bi*mi+iz+1;
   Ja(r,s)+=Jai(ai+1,bi+1);
  }
 }
 // transition matrix element
 iJj(a)=complex<double>(aMb_real(Ja,zr,zc,ii,jj),aMb_imag(Ja,zr,zc,ii,jj));

 // determine expectation value
 Jret(a)=0;
if (subtractexpvalue==1)
{ for(int i=1;i<=dim&&wn[i]>0.00001;++i)
 {Jret(a)+=wn[i]*aMb_real(Ja,zr,zc,i,i);
 }
}
}

delete []d0;delete []dn;


// 3. set u1
for(int l=1;l<=gjmbH.Hi();++l)
{if(ii==jj){//take into account thermal expectation values <Jret>
          u1(l)=(iJj(l)-Jret(l));}
 else    {u1(l)=iJj(l);}
}

if (delta>SMALL)
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
 return (int)(dim*(dim+1)/2);



}
/**************************************************************************/

