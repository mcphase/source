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

//ComplexMatrix & mdcf::V(int na, int nb, int nc) const
//{ return (*sb[in(na,nb,nc)]);}

ComplexVector & mdcf::dPs(int na, int nb, int nc) const
{ return (*dps[in(na,nb,nc)]);
}
ComplexVector & mdcf::dMQs(int na, int nb, int nc) const
{ return (*dmqs[in(na,nb,nc)]);
}
ComplexVector & mdcf::dMQ_dips(int na, int nb, int nc) const
{ return (*dmq_dips[in(na,nb,nc)]);
}
ComplexMatrix & mdcf::M(int na, int nb, int nc)
{ if(!mr){errexit();}
return (*m[in(na,nb,nc)]);
}

ComplexVector & mdcf::sqrt_gamma(int na, int nb, int nc) const
{ return (*l[in(na,nb,nc)]);
}
ComplexVector & mdcf::sqrt_GammaP(int na, int nb, int nc) const
{ return (*Pb[in(na,nb,nc)]);
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
//ComplexMatrix & mdcf::Vi(int i)
//{ return (*sb[i]);}

ComplexMatrix & mdcf::Mi(int i)
{ if(!mr){ errexit();}
return (*m[i]);
}

ComplexVector & mdcf::sqrt_gammai(int i)
{ return (*l[i]);
}
ComplexVector & mdcf::sqrt_GammaPi(int i)
{ return (*Pb[i]);
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
int mdcf::inM(int i, int j, int k, int l)
{//int indd=(((i*mxb+j)*mxc+k)*nofatoms+l);
 //if(indd<0||indd>mxa*mxb*mxc*(nofatoms+1)) {fprintf(stderr,"mdcf indexing error");exit(EXIT_FAILURE);}
 int indd=((((i-1)*nofb+j-1)*nofc+k-1)*nofatoms+l-1);
 return indd;}

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


ComplexMatrix & mdcf::est(int i, int j, int k, int l)
{return (*eigenstates[ind(i,j,k,l)]);}

ComplexMatrix ** mdcf::chi0pointer(int i, int j, int k, int l)
{return chi0s[ind(i,j,k,l)];}


void mdcf::est_ini(int i, int j, int k, int l,ComplexMatrix & M) // initialize est
{eigenstates[ind(i,j,k,l)]=new ComplexMatrix(M.Rlo(),M.Rhi(),M.Clo(),M.Chi());
 (*eigenstates[ind(i,j,k,l)])=M;
}

// has to be called before mdcf object can be used for calculation
void mdcf::set_noftransitions(int i, int j, int k, IntVector & notr,int mqd)
{      (*nt[in(i,j,k)])=notr;mqdim=mqd;
 if(storage){
       int sumnt=sum((*nt[in(i,j,k)]));
      if (sumnt<1){sumnt=1;} // MR 2011.08.09. in case there is no transition for this subsystem then initialize all matrices/Vector to 1...
     s[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sumnt,1,sumnt);if(s[in(i,j,k)]==NULL)errexit();
  if(mr){m[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sumnt,1,nofcomponents*sumnt);if(m[in(i,j,k)]==NULL)errexit();}
     l[in(i,j,k)]= new ComplexVector(1,sumnt);if(l[in(i,j,k)]==NULL)errexit();
  //   sb[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sumnt,1,sumnt);if(sb[in(i,j,k)]==NULL)errexit();// second index only integer nofcomponents needed, so runs from 1-sumnt MR 14.9.2011
     dps[in(i,j,k)]= new ComplexVector(1,1*sumnt);if(dps[in(i,j,k)]==NULL)errexit();
     dmqs[in(i,j,k)]= new ComplexVector(1,mqdim*sumnt);if(dmqs[in(i,j,k)]==NULL)errexit();
     dmq_dips[in(i,j,k)]= new ComplexVector(1,mqdim*sumnt);if(dmq_dips[in(i,j,k)]==NULL)errexit();
     Pb[in(i,j,k)]= new ComplexVector(1,sumnt);if(Pb[in(i,j,k)]==NULL)errexit();
     lb[in(i,j,k)]= new ComplexVector(1,sumnt);if(lb[in(i,j,k)]==NULL)errexit();
     lb_dip[in(i,j,k)]= new ComplexVector(1,sumnt);if(lb_dip[in(i,j,k)]==NULL)errexit();
     d[in(i,j,k)]= new Vector(1,sumnt);   if(d[in(i,j,k)]==NULL)errexit();  
     } 
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

void mdcf::errexit()
{ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

//constructors
mdcf::mdcf (int n1,int n2,int n3,int n,int nc,int stps,int do_Erefine)
{  int i; mr=do_Erefine;
   nofa=n1;nofb=n2;nofc=n3;
   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
   nofatoms=n;nofcomponents=nc;
   nofEstps=stps;
   storage=1;
//dimension arrays
  s = new ComplexMatrix * [mxa*mxb*mxc+1];if(s==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)s[i]=NULL;
  m = new ComplexMatrix * [mxa*mxb*mxc+1];if(m==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)m[i]=NULL;
  l = new ComplexVector * [mxa*mxb*mxc+1];if(l==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)l[i]=NULL;
 // sb = new ComplexMatrix * [mxa*mxb*mxc+1];if(sb==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)sb[i]=NULL;
  dps = new ComplexVector * [mxa*mxb*mxc+1];if(dps==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)dps[i]=NULL;
  dmqs = new ComplexVector * [mxa*mxb*mxc+1];if(dmqs==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)dmqs[i]=NULL;
  dmq_dips = new ComplexVector * [mxa*mxb*mxc+1];if(dmq_dips==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)dmq_dips[i]=NULL;
  Pb = new ComplexVector * [mxa*mxb*mxc+1];if(Pb==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)Pb[i]=NULL;
  lb = new ComplexVector * [mxa*mxb*mxc+1];if(lb==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)lb[i]=NULL;
  lb_dip = new ComplexVector * [mxa*mxb*mxc+1];if(lb_dip==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)lb_dip[i]=NULL;
  d = new Vector * [mxa*mxb*mxc+1]; if(d==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)d[i]=NULL;
  nt= new IntVector * [mxa*mxb*mxc+1];if(nt==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i){nt[i]=new IntVector(1,nofatoms);}
  eigenstates= new ComplexMatrix * [mxa*mxb*mxc*(nofatoms+1)+1]; if(eigenstates==NULL)errexit();
                           for(i=0;i<=mxa*mxb*mxc*(nofatoms+1);++i)eigenstates[i]=NULL; 
  if(nofEstps>0){
  chi0s= new ComplexMatrix ** [mxa*mxb*mxc*(nofatoms+1)+1]; if(chi0s==NULL)errexit();
  for(i=0;i<=mxa*mxb*mxc*(nofatoms+1);++i){chi0s[i]=new ComplexMatrix*[nofEstps+1];  
    for(int j=0;j<nofEstps;++j)chi0s[i][j]=new ComplexMatrix(1,nofcomponents,1,nofcomponents);
                                           }
                } 
  Ug=0; gU=0; bUg=0; bgU=0; PUg=0; PgU=0;
}

//kopier-konstruktor
mdcf::mdcf (const mdcf & p,int store)
{ int i,j,k;mr=p.mr;
  nofa=p.nofa;nofb=p.nofb;nofc=p.nofc;
  mxa=p.mxa; mxb=p.mxb; mxc=p.mxc;
  mqdim=p.mqdim;
  nofEstps=p.nofEstps;
  nofatoms=p.nofatoms;nofcomponents=p.nofcomponents;
  storage=store;
  nt= new IntVector * [mxa*mxb*mxc+1];if(nt==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i){nt[i]=new IntVector(1,nofatoms);}
  for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k)
     {int id=in(i,j,k); 
      (*nt[id])=(*p.nt[id]);}}}
  if(storage){
//dimension arrays
  s = new ComplexMatrix * [mxa*mxb*mxc+1];if(s==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)s[i]=NULL;
  m = new ComplexMatrix * [mxa*mxb*mxc+1];if(m==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)m[i]=NULL;
  l = new ComplexVector * [mxa*mxb*mxc+1];if(l==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)l[i]=NULL;
 // sb = new ComplexMatrix * [mxa*mxb*mxc+1];if(sb==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)sb[i]=NULL;
  dps = new ComplexVector * [mxa*mxb*mxc+1];if(dps==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)dps[i]=NULL;
  dmqs = new ComplexVector * [mxa*mxb*mxc+1];if(dmqs==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)dmqs[i]=NULL;
  dmq_dips = new ComplexVector * [mxa*mxb*mxc+1];if(dmq_dips==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)dmq_dips[i]=NULL;
  Pb = new ComplexVector * [mxa*mxb*mxc+1];if(Pb==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)Pb[i]=NULL;
  lb = new ComplexVector * [mxa*mxb*mxc+1];if(lb==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)lb[i]=NULL;
  lb_dip = new ComplexVector * [mxa*mxb*mxc+1];if(lb_dip==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)lb_dip[i]=NULL;
  d = new Vector * [mxa*mxb*mxc+1]; if(d==NULL)errexit();for(i=0;i<=mxa*mxb*mxc;++i)d[i]=NULL;    
  eigenstates= new ComplexMatrix * [mxa*mxb*mxc*(nofatoms+1)+1]; if(eigenstates==NULL)errexit();
                           for(i=0;i<=mxa*mxb*mxc*(nofatoms+1);++i){eigenstates[i]=NULL;       
                                                                    if(p.eigenstates[i]!=NULL)eigenstates[i]=new ComplexMatrix((*p.eigenstates[i]));
                                                                    }

  if(nofEstps>0){
  chi0s= new ComplexMatrix ** [mxa*mxb*mxc*(nofatoms+1)+1]; if(chi0s==NULL)errexit();
  for(i=0;i<=mxa*mxb*mxc*(nofatoms+1);++i){chi0s[i]=new ComplexMatrix*[nofEstps+1];  
    for(int j=0;j<nofEstps;++j){if(p.chi0s[i][j]!=NULL)
                               {chi0s[i][j]=new ComplexMatrix((*p.chi0s[i][j]));if(chi0s[i][j]==NULL)errexit();}
                                }
                                           }
                }
 for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k)
     {int id=in(i,j,k);       
      s[id]= new ComplexMatrix(1,nofcomponents*sum((*nt[id])),1,sum((*nt[id])));*s[id]=*p.s[id];
      if(mr){m[id]= new ComplexMatrix(1,nofcomponents*sum((*nt[id])),1,nofcomponents*sum((*nt[id])));*m[id]=*p.m[id];}
      l[id]= new ComplexVector(1,sum((*nt[id])));*l[id]=*p.l[id];
    //  sb[id]= new ComplexMatrix(1,nofcomponents*sum((*nt[id])),1,sum((*nt[id])));*sb[id]=*p.sb[id];// second index only integer nofcomponents needed, so runs from 1-sumnt MR 14.9.2011
      dps[id]= new ComplexVector(1,1*sum((*nt[id])));*dps[id]=*p.dps[id];
      dmqs[id]= new ComplexVector(1,mqdim*sum((*nt[id])));*dmqs[id]=*p.dmqs[id];
      dmq_dips[id]= new ComplexVector(1,mqdim*sum((*nt[id])));*dmq_dips[id]=*p.dmq_dips[id];
      Pb[id]= new ComplexVector(1,sum((*nt[id])));*Pb[id]=*p.Pb[id];
      lb[id]= new ComplexVector(1,sum((*nt[id])));*lb[id]=*p.lb[id];
      lb_dip[id]= new ComplexVector(1,sum((*nt[id]))); *lb_dip[id]=*p.lb_dip[id];
      d[id]= new Vector(1,sum((*nt[id])));*d[id]=*p.d[id];     
     } 
    }
  }           
 }
  Ug=0; gU=0; bUg=0; bgU=0;PUg=0; PgU=0;
} 
//destruktor
mdcf::~mdcf ()
{int i,j,k;
 for(i=0;i<=mxa*mxb*mxc;++i){delete nt[i];}
 delete []nt;
 if(storage)
 { for(i=0;i<=mxa*mxb*mxc;++i){
 if(s[i]!=NULL)delete s[i];
 if(mr){if(m[i]!=NULL)delete m[i];}
 if(l[i]!=NULL)delete l[i];
 // if(sb[i]!=NULL)delete sb[i];
 if(dps[i]!=NULL)delete dps[i];
 if(dmqs[i]!=NULL)delete dmqs[i];
 if(dmq_dips[i]!=NULL)delete dmq_dips[i];
 if(Pb[i]!=NULL)delete Pb[i];
 if(lb[i]!=NULL)delete lb[i];
 if(lb_dip[i]!=NULL)delete lb_dip[i];
 if(d[i]!=NULL)delete d[i];
    }
 delete []s;
 delete []m;
 delete []l;
 // delete []sb;
 delete []dps;
 delete []dmqs;
 delete []dmq_dips;
 delete []Pb;
 delete []lb;
 delete []lb_dip;
 delete []d;
 
 for(i=0;i<=mxa*mxb*mxc*(nofatoms+1);++i)if(eigenstates[i]!=NULL)delete eigenstates[i];
 delete []eigenstates;

  if(nofEstps>0){
  for(i=0;i<=mxa*mxb*mxc*(nofatoms+1);++i){  
                                           for(int j=0;j<nofEstps;++j)if(chi0s[i][j]!=NULL)delete chi0s[i][j];
                                           if(chi0s[i]!=NULL)delete []chi0s[i];
                                           }
    delete []chi0s;
                }

 for (i=1;i<=nofa;++i){ for (j=1;j<=nofb;++j){ for (k=1;k<=nofc;++k){
 int id = in(i,j,k);
 // For caching values in calculation of transform of chi''
 if(Ug!=0){if(Ug[id]!=0) delete Ug[id];}
 if(gU!=0){if(gU[id]!=0) delete gU[id]; }
 if(bUg!=0){if(bUg[id]!=0) delete bUg[id]; }
 if(bgU!=0){if(bgU[id]!=0) delete bgU[id];}
 if(PUg!=0){if(PUg[id]!=0) delete PUg[id]; }
 if(PgU!=0){if(PgU[id]!=0) delete PgU[id];}
 }}}
 if(Ug!=0) { delete []Ug; Ug=0; } 
 if(gU!=0) { delete []gU; gU=0; } 
 if(bUg!=0) { delete []bUg; bUg=0; }
 if(bgU!=0) { delete []bgU; bgU=0; }
 if(PUg!=0) { delete []PUg; PUg=0; }
 if(PgU!=0) { delete []PgU; PgU=0; }
 }
}
