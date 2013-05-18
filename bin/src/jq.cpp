 // *************************************************************************
 // ************************ jq *************************************
 // *************************************************************************
// methods for class jq
#include <cerrno>
#include <cstdio>
#include <cmath>
#include <martin.h>
#include <vector.h>
#include <complex>
#include "jq.hpp"

#define MAXNOFSPINS  200

// returns jq of spin [i=na,j=nb,k=nc] 
ComplexMatrix & jq::mat(int na, int nb, int nc,int ma,int mb,int mc)
{ return (*jj[iin(in(na,nb,nc),in(ma,mb,mc))]);
}
// the same but for spin number "i,j"
ComplexMatrix & jq::mati(int i,int j)
{ return (*jj[iin(i,j)]);
}
// get index ijk=iv(1-3)  of spinconfiguration number in 
int * jq::ijk(int in)
{div_t result; result=div(in,mxb*mxc); 
 iv[1]= result.quot;
 result=div(result.rem,mxc);
 iv[2]= result.quot;
 iv[3]= result.rem;
 return iv;}

// the inverse: get number of spin from indizes i,j,k
int jq::in(int i, int j, int k)
{int ret;ret=(i*mxb+j)*mxc+k;
 if(ret>mx-1){fprintf(stderr,"ERROR mcdisp - internal indexing error\n");exit(EXIT_FAILURE);}
 return ret;}

int jq::iin(int i, int j)
{return (i*mx+j);}


// return number of spins
int jq::n()
{return (nofa*nofb*nofc);
}
int jq::na()
{return nofa;
}
int jq::nb()
{return nofb;
}
int jq::nc()
{return nofc;
}


/**************************************************************************/

//constructors
jq::jq (int n1,int n2,int n3,int nofcomp,mdcf & m)
{  nofa=n1;nofb=n2;nofc=n3;
   nofatoms=m.nofatoms;nofcomponents=nofcomp;

   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
   mx=mxa*mxb*mxc+1;
   
//dimension arrays
  jj = new  ComplexMatrix * [mx*mx+1];
//  (1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (jj == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  int mbm1,mbm2;
 int i1,j1,k1,i2,j2,k2;
 for (i1=1;i1<=nofa;++i1){
 for (j1=1;j1<=nofb;++j1){
 for (k1=1;k1<=nofc;++k1){mbm1=m.baseindex_max(i1,j1,k1);if (mbm1<1){mbm1=1;}
 for (i2=1;i2<=nofa;++i2){
 for (j2=1;j2<=nofb;++j2){
 for (k2=1;k2<=nofc;++k2){mbm2=m.baseindex_max(i2,j2,k2);if(mbm2<1){mbm2=1;} // MR 2011.08.08 for zero transitions in subsystem
 jj[iin(in(i1,j1,k1),in(i2,j2,k2))]= new ComplexMatrix(1,nofcomponents*mbm1,1,nofcomponents*mbm2);
 }}}}}} 
}

//kopier-konstruktor
jq::jq (const jq & p)
{ int ind;
  nofa=p.nofa;nofb=p.nofb;nofc=p.nofc;
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  mxa=p.mxa; mxb=p.mxb; mxc=p.mxc;
  mx=p.mx;
  
//dimension arrays
  jj = new ComplexMatrix* [mx*mx+1];
  if (jj == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  //(1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);

 for (int i1=1;i1<=nofa;++i1){ for (int j1=1;j1<=nofb;++j1){ for (int k1=1;k1<=nofc;++k1){
 for (int i2=1;i2<=nofa;++i2){ for (int j2=1;j2<=nofb;++j2){ for (int k2=1;k2<=nofc;++k2){
    ind = iin(in(i1,j1,k1),in(i2,j2,k2)); 
    jj[ind]= new ComplexMatrix(1,(*p.jj[ind]).Rhi(),1,(*p.jj[ind]).Chi()); 
    (*jj[ind])=(*p.jj[ind]);
 }}}}}} 

}

//destruktor
jq::~jq ()
{
 int i1,j1,k1,i2,j2,k2;
 for (i1=1;i1<=nofa;++i1){
 for (j1=1;j1<=nofb;++j1){
 for (k1=1;k1<=nofc;++k1){
 for (i2=1;i2<=nofa;++i2){
 for (j2=1;j2<=nofb;++j2){
 for (k2=1;k2<=nofc;++k2){
 delete jj[iin(in(i1,j1,k1),in(i2,j2,k2))];
 }}}}}} 
  delete []jj;
}


