 // *************************************************************************
 // ****** md cf - single ion transition matrix elements storage ************
 // *************************************************************************
// methods for class mdcf 
#include <cerrno>
#include <cstdio>
#include <cmath>
#include <martin.h>
#include <vector.h>
#include <complex>
#include "mdcf.hpp"

#define MAXNOFSPINS  200

// returns sum of all components of intvector
int sum(IntVector & v)
 {int i,sum=0;
 for (i=v.Lo();i<=v.Hi();++i)
 {sum+=v(i);}
  return sum;
}

// returns md of cf [i=na,j=nb,k=nc] 
ComplexMatrix & mdcf::U(int na, int nb, int nc) const
{ return (*s[in(na,nb,nc)]);
}
ComplexMatrix & mdcf::V(int na, int nb, int nc) const
{ return (*sb[in(na,nb,nc)]);
}
ComplexVector & mdcf::dMQs(int na, int nb, int nc) const
{ return (*dmqs[in(na,nb,nc)]);
}
ComplexVector & mdcf::dMQ_dips(int na, int nb, int nc) const
{ return (*dmq_dips[in(na,nb,nc)]);
}
ComplexMatrix & mdcf::M(int na, int nb, int nc)
{ return (*m[in(na,nb,nc)]);
}

ComplexMatrix & mdcf::sqrt_gamma(int na, int nb, int nc) const
{ return (*l[in(na,nb,nc)]);
}
ComplexVector & mdcf::sqrt_Gamma(int na, int nb, int nc) const
{ return (*lb[in(na,nb,nc)]);
}
ComplexVector & mdcf::sqrt_Gamma_dip(int na, int nb, int nc) const
{ return (*lb_dip[in(na,nb,nc)]);
}
Vector & mdcf::delta(int na, int nb, int nc)
{ return (*d[in(na,nb,nc)]);
}
// the same but for cf number "i"
ComplexMatrix & mdcf::Ui(int i)
{ return (*s[i]);
}
ComplexMatrix & mdcf::Vi(int i)
{ return (*sb[i]);
}
ComplexMatrix & mdcf::Mi(int i)
{ return (*m[i]);
}

ComplexMatrix & mdcf::sqrt_gammai(int i)
{ return (*l[i]);
}
ComplexVector & mdcf::sqrt_Gammai(int i)
{ return (*lb[i]);
}
ComplexVector & mdcf::sqrt_Gamma_dipi(int i)
{ return (*lb_dip[i]);
}
Vector  & mdcf::deltai(int i)
{ return (*d[i]);
}
// get index ijk=iv(1-3)  of cf configuration number in 
int * mdcf::ijk(int in)
{div_t result; result=div(in,mxb*mxc); 
 iv[1]= result.quot;
 result=div(result.rem,mxc);
 iv[2]= result.quot;
 iv[3]= result.rem;
 return iv;}

// the inverse: get number of cf from indizes i,j,k
int mdcf::in(int i, int j, int k) const
{return ((i*mxb+j)*mxc+k);}

// get number of cf from indizes i,j,k,l
int mdcf::ind(int i, int j, int k, int l)
{int indd=(((i*mxb+j)*mxc+k)*nofatoms+l);
 if(indd<0||indd>mxa*mxb*mxc*(nofatoms+1)) {fprintf(stderr,"mdcf indexing error");exit(EXIT_FAILURE);}
 return indd;}

// return number of cfs
int mdcf::n()
{return (nofa*nofb*nofc);
}
int mdcf::na()
{return nofa;
}
int mdcf::nb()
{return nofb;
}
int mdcf::nc()
{return nofc;
}

/**************************************************************************/

//constructors
mdcf::mdcf (int n1,int n2,int n3,int n,int nc)
{  int i;
   nofa=n1;nofb=n2;nofc=n3;
   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
   nofatoms=n;nofcomponents=nc;

//dimension arrays
  s = new ComplexMatrix * [mxa*mxb*mxc+1];//(1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (s == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  m = new ComplexMatrix * [mxa*mxb*mxc+1];
  if (m == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  l = new ComplexMatrix * [mxa*mxb*mxc+1];
  if (l == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  sb = new ComplexMatrix * [mxa*mxb*mxc+1];
  if (sb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  dmqs = new ComplexVector * [mxa*mxb*mxc+1];
  if (dmqs == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  dmq_dips = new ComplexVector * [mxa*mxb*mxc+1];
  if (dmq_dips == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  lb = new ComplexVector * [mxa*mxb*mxc+1];
  if (lb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  lb_dip = new ComplexVector * [mxa*mxb*mxc+1];
  if (lb_dip == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  d = new Vector * [mxa*mxb*mxc+1]; 
  if (d == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  nt= new IntVector * [mxa*mxb*mxc+1];
  if (nt == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  for(i=0;i<=mxa*mxb*mxc;++i){nt[i]=new IntVector(1,nofatoms);}

  eigenstates= new ComplexMatrix * [mxa*mxb*mxc*(nofatoms+1)+1];   
  if (eigenstates == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  
  Ug=0; gU=0; bUg=0; bgU=0;
}

ComplexMatrix & mdcf::est(int i, int j, int k, int l)
{return (*eigenstates[ind(i,j,k,l)]);}


void mdcf::est_ini(int i, int j, int k, int l,ComplexMatrix & M) // initialize est
{eigenstates[ind(i,j,k,l)]=new ComplexMatrix(M.Rlo(),M.Rhi(),M.Clo(),M.Chi());
 (*eigenstates[ind(i,j,k,l)])=M;
}

// has to be called before mdcf object can be used for calculation
void mdcf::set_noftransitions(int i, int j, int k, IntVector & notr)
{      (*nt[in(i,j,k)])=notr;
       int sumnt=sum((*nt[in(i,j,k)]));
      if (sumnt<1){sumnt=1;} // MR 2011.08.09. in case there is no transition for this subsystem then initialize all matrices/Vector to 1...
     s[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sumnt,1,nofcomponents*sumnt);
     m[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sumnt,1,nofcomponents*sumnt);
     l[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sumnt,1,nofcomponents*sumnt);
     sb[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sumnt,1,sumnt);// second index only integer nofcomponents needed, so runs from 1-sumnt MR 14.9.2011
     dmqs[in(i,j,k)]= new ComplexVector(1,3*sumnt);
     dmq_dips[in(i,j,k)]= new ComplexVector(1,3*sumnt);
     lb[in(i,j,k)]= new ComplexVector(1,3*sumnt);
     lb_dip[in(i,j,k)]= new ComplexVector(1,3*sumnt);
     d[in(i,j,k)]= new Vector(1,sumnt);
      
}

int mdcf::baseindex(int i, int j, int k, int l, int tn) const
{// the baseindex is used to number the rows an columns of the
 // matrices associated with the crystallographic unit number ijk
 // it starts at one and combines indices l (atom number) and t (number of
 // transition for this atom) into a single index bi, its maximum value is
 // baseindex_max
  int bi=0;
  int i1;
for(i1=1;i1<=l;++i1)
{bi+=(*nt[in(i,j,k)])(i1);}
bi-=(*nt[in(i,j,k)])(l);
bi+=tn;
return bi;
}

int mdcf::baseindex_max(int i, int j, int k) const
{return sum((*nt[in(i,j,k)]));}

int mdcf::noft(int i, int j, int k,int l) const
{return (*nt[in(i,j,k)])(l);}

//destruktor
mdcf::~mdcf ()
{int i,j,k;
 for (i=1;i<=nofa;++i){
 for (j=1;j<=nofb;++j){
 for (k=1;k<=nofc;++k){
 int id = in(i,j,k);
 delete s[id];
 delete m[id];
 delete l[id];
 delete sb[id];
 delete dmqs[id];
 delete dmq_dips[id];
 delete lb[id];
 delete lb_dip[id];
 delete d[id];
 // For caching values in calculation of transform of chi''
 if(Ug[id]!=0) delete Ug[id]; if(gU[id]!=0) delete gU[id]; if(bUg[id]!=0) delete bUg[id]; if(bgU[id]!=0) delete bgU[id];
 }}}
 delete []s;delete []m;delete []d;delete []l;delete []nt;
 delete []sb;
 delete []dmqs;
 delete []dmq_dips;
delete []lb;
delete []lb_dip;
 if(Ug!=0) { delete []Ug; Ug=0; } if(bUg!=0) { delete []bUg; bUg=0; }
 if(gU!=0) { delete []gU; gU=0; } if(bgU!=0) { delete []bgU; bgU=0; }
}



//kopier-konstruktor
mdcf::mdcf (const mdcf & p)
{ int i,j,k;
  nofa=p.nofa;nofb=p.nofb;nofc=p.nofc;
  mxa=p.mxa; mxb=p.mxb; mxc=p.mxc;
  nofatoms=p.nofatoms;nofcomponents=p.nofcomponents;
//dimension arrays
  s = new ComplexMatrix*[mxa*mxb*mxc+1]; //(1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (s == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  m = new ComplexMatrix*[mxa*mxb*mxc+1];
  if (m == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  l = new ComplexMatrix*[mxa*mxb*mxc+1];
  if (l == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  sb = new ComplexMatrix *[mxa*mxb*mxc+1];//(1,nofcomponents*nofatoms,1,nofatoms);
  if (sb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  dmqs = new ComplexVector *[mxa*mxb*mxc+1];
  if (dmqs == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  dmq_dips = new ComplexVector *[mxa*mxb*mxc+1];
  if (dmq_dips == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  lb = new ComplexVector *[mxa*mxb*mxc+1];
  if (lb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  lb_dip = new ComplexVector *[mxa*mxb*mxc+1];
  if (lb_dip == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  d = new Vector*[mxa*mxb*mxc+1];//(1,nofatoms);
  if (d == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  nt = new IntVector*[mxa*mxb*mxc+1];//(1,nofatoms);
  if (nt == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  eigenstates= new ComplexMatrix * [mxa*mxb*mxc*nofatoms+1];   
  if (eigenstates == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 


 for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {
      nt[in(i,j,k)]=p.nt[in(i,j,k)];

      s[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      m[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      l[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      sb[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,sum((*nt[in(i,j,k)])));// second index only integer nofcomponents needed, so runs from 1-sumnt MR 14.9.2011
      dmqs[in(i,j,k)]= new ComplexVector(1,3*sum((*nt[in(i,j,k)])));
      dmq_dips[in(i,j,k)]= new ComplexVector(1,3*sum((*nt[in(i,j,k)])));
      lb[in(i,j,k)]= new ComplexVector(1,3*sum((*nt[in(i,j,k)])));
      lb_dip[in(i,j,k)]= new ComplexVector(1,3*sum((*nt[in(i,j,k)])));
      d[in(i,j,k)]= new Vector(1,sum((*nt[in(i,j,k)])),1,sum((*nt[in(i,j,k)])));

      *d[in(i,j,k)]=*p.d[in(i,j,k)];
      *s[in(i,j,k)]=*p.s[in(i,j,k)];
      *m[in(i,j,k)]=*p.m[in(i,j,k)];
      *l[in(i,j,k)]=*p.l[in(i,j,k)];
      *sb[in(i,j,k)]=*p.sb[in(i,j,k)];
      *dmqs[in(i,j,k)]=*p.dmqs[in(i,j,k)];
      *dmq_dips[in(i,j,k)]=*p.dmq_dips[in(i,j,k)];
      *lb[in(i,j,k)]=*p.lb[in(i,j,k)];
      *lb_dip[in(i,j,k)]=*p.lb_dip[in(i,j,k)];
     } 
    }
  }           

  Ug=0; gU=0; bUg=0; bgU=0;
}
