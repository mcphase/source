 // *************************************************************************
 // ************************ spincf *************************************
 // *************************************************************************
// methods for class spincf 
#include <martin.h>
#include <complex>
#include "spincf.hpp"
#include "density.hpp"
#include "jjjpar.hpp"

#define MAXNOFSPINS  200
#define SMALLSPINDIFF 0.03
// output functions
#include "spincf_out.cpp"

int spincf::spequal(Vector a,Vector b)
{// in this function we look if two spins are equal or not (used in order to compare or reduce
 // spinconfigurations
 //if(a*b/Norm(a)/Norm(b)<COSINE){return false;} // this is by the inner product

 //this is by looking if any of the components of the 2 spins have different signs
 int i;
 for (i=a.Lo();i<=a.Hi();++i)
 {if (fabs(a(i))>SMALLSPINDIFF||fabs(b(i))>SMALLSPINDIFF)
  {if(signum(a(i))!=signum(b(i))){return false;}
  }
 }
return true;
}

// returns momentum of ion [i=na,j=nb,k=nc] 
Vector & spincf::m(int na, int nb, int nc)
{ return mom[in(na,nb,nc)];  // return <J(ijk)>
}

// the same but for spin number "i"
Vector & spincf::mi(int i)
{ return mom[i];
}

Vector spincf::moment(int i,int j,int k,int l) // returns moment of atom l (1,nofcomponents)
{Vector xyz(1,nofcomponents);
 int m;
 for(m=1;m<=nofcomponents;++m){xyz(m)=mom[in(i,j,k)](nofcomponents*(l-1)+m);}
 return xyz;
}

// (Slow) Fourier Transform of momentumconfiguration <I> 
void  spincf::FT(ComplexVector * mq)
{int i,j,k,l;
 int qh,qk,ql;
    complex<double> piq(0,2*3.1415926535),g;
  for(qh=0;qh<na();++qh){for(qk=0;qk<nb();++qk){for(ql=0;ql<nc();++ql)
     {mq[in(qh,qk,ql)]=0;
     for(i=0;i<na();++i){for(j=0;j<nb();++j){for(k=0;k<nc();++k)
      { g=exp(piq*((double)qh*i/na()+(double)qk*j/nb()+(double)ql*k/nc()));
	for (l=1;l<=nofcomponents*nofatoms;++l) {mq[in(qh,qk,ql)](l)+=g*m(i+1,j+1,k+1)(l);}
      }}}
//      fprintf(stdout,"%g %g %g %g %g \n",qh,qk,ql,imag(mq[in(qh,qk,ql)](1)),real(mq[in(qh,qk,ql)](1)));
     }}}
 for (qk=1;qk<=nb();++qk)
 {for (ql=1;ql<=nc();++ql) mq[in(na(),qk,ql)]=mq[in(0,qk,ql)];}
 for (qh=1;qh<=na();++qh)
 {for (ql=1;ql<=nc();++ql) mq[in(qh,nb(),ql)]=mq[in(qh,0,ql)];}
 for (qh=1;qh<=na();++qh)
 {for (qk=1;qk<=nb();++qk) mq[in(qh,qk,nc())]=mq[in(qh,qk,0)];}
 mq[in(na(),nb(),nc())]=mq[in(0,0,0)];
 return;
}

// get index ijk=iv(1-3)  of spinconfiguration number in 
int * spincf::ijk(int in)
{div_t result; result=div(in,mxb*mxc); 
 iv[1]= result.quot;
 result=div(result.rem,mxc);
 iv[2]= result.quot;
 iv[3]= result.rem;
 return iv;}

// the inverse: get number of spin from indizes i,j,k
int spincf::in(int i, int j, int k)
{while (i>nofa)i-=nofa;
 while (j>nofb)j-=nofb;
 while (k>nofc)k-=nofc;
 return ((i*mxb+j)*mxc+k);
}

// this subtracts n2 if n1>n2
int spincf::mod(int n1,int n2)
{if (n1>n2) return n1-n2;
 else return n1;
}

// return number of spins
int spincf::n()
{return (nofa*nofb*nofc);
}
int spincf::na()
{return nofa;
}
int spincf::nb()
{return nofb;
}
int spincf::nc()
{return nofc;
}


Vector spincf::totalJ() // returns total moment <J>
{Vector ret(1,nofcomponents);
 ret=0;
 int i,j,k,l,m; 
 for (i=1;i<=nofa;++i)
 { for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {for (m=1;m<=nofcomponents;++m)
      {for (l=1;l<=nofatoms;++l)
       { ret(m)+=mom[in(i,j,k)](nofcomponents*(l-1)+m);
       }      
      }    
     }
    }
 }     
 ret/=n()*nofatoms;
 return ret;
}


// invert all spins (AND higher order moments)
void spincf::invert()
{
 int i,j,k; 
 for (i=1;i<=nofa;++i)
 { for (j=1;j<=nofb;++j)
   {for (k=1;k<=nofc;++k)
    {mom[in(i,j,k)]= -mom[in(i,j,k)];}
   }
 }  
}

// reduce spinconfiguration if it is periodic
// returns 1 if successful, 0 if no reduction is possible
int spincf::reduce()
{int perprob,nn,pp,j,k,success=0;
 div_t result;

//reduce nofa
pp=(int)integer(nofa/2);
for (perprob=1;perprob<=pp;++perprob)
  {result=div(nofa,perprob);
   if (result.rem==0)
    {
     for (j=1;j<=nofb;++j)
      {for (k=1;k<=nofc;++k)
       {for (nn=1;nn<=nofa-perprob;++nn)
        {if (!spequal(mom[in(nn,j,k)],mom[in(nn+perprob,j,k)])){goto nexta;}
	 if(nn==nofa-perprob&&j==nofb&&k==nofc)
        	{this->nofa=perprob;
	// here we could reduce the memory of this confuration
	// however for the moment we leave it for it does not disturb
                 success=1;
	         goto reduceda;
	        }  
	}
       }
      } 	
    }
nexta:;
  }
reduceda:


//reduce nofb
pp=(int)integer(nofb/2);
for (perprob=1;perprob<=pp;++perprob)
  {result=div(nofb,perprob);
   if (result.rem==0)
    {
     for (j=1;j<=nofa;++j)
      {for (k=1;k<=nofc;++k)
       {for (nn=1;nn<=nofb-perprob;++nn)
        {if(!spequal(mom[in(j,nn,k)],mom[in(j,nn+perprob,k)])){goto nextb;}
	 if(nn==nofb-perprob&&j==nofa&&k==nofc)
        	{this->nofb=perprob;
	// here we could reduce the memory of this confuration
	// however for the moment we leave it for it does not disturb
                 success=1;
	         goto reducedb;
	        }  
	}
       } 	
      }
    }
nextb:;
  }    
reducedb:
//reduce nofc
pp=(int)integer(nofc/2);
for (perprob=1;perprob<=pp;++perprob)
  {result=div(nofc,perprob);
   if (result.rem==0)
    {
     for (j=1;j<=nofa;++j)
      {for (k=1;k<=nofb;++k)
       {for (nn=1;nn<=nofc-perprob;++nn)
        {if(!spequal(mom[in(j,k,nn)],mom[in(j,k,nn+perprob)])){goto nextc;}
	 if(nn==nofc-perprob&&j==nofa&&k==nofb)
        	{this->nofc=perprob;
	// here we could reduce the memory of this confuration
	// however for the moment we leave it for it does not disturb
                 success=1;
	         return success;
	        }  
	}
       }
      } 	
    }
nextc:;
  }

return success;
}


  
// create spinconfig from qvector, momentq0
// n1 n2 n3 .. periodicity of supercell
// nettom .... saturation moment (positive)
// qvector ... wave vector in units of reciprocal lattice
// momentq0 .. ferromagnetic component (between 0 and 1)
// phi ....... phase (for each component)
void spincf::spinfromq (int n1,int n2, int n3,Vector & qvector,Vector & nettom, 
       Vector & momentq0, Vector & phi)
{ int rra,rrb,rrc,qv;
  qv=1; // evtl koennte hier noch ein vorzeichen uebergeben werden (sodass nettom positiv)
  Vector rr(1,3);
  delete []mom;
  epsilon=0;
  nofa=n1; nofb=n2; nofc=n3;
  mxa=n1+1; mxb=n2+1; mxc=n3+1;
  int l;
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];for(l=0;l<=mxa*mxb*mxc;++l){mom[l]=Vector(1,nofcomponents*nofatoms);}
  if (mom == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 

// fill spinconfiguration with values
   for (rra=1;rra<=nofa;++rra)
     { for (rrb=1;rrb<=nofb;++rrb)
        {for (rrc=1;rrc<=nofc;++rrc)
             	{rr(1)=rra-1;rr(2)=rrb-1;rr(3)=rrc-1;
		 for(l=1;l<=nofcomponents*nofatoms;++l)
	          {
//		  mom[in(rra,rrb,rrc)](l)=copysign(nettom(l),qv*(0.01+
//		  copysign(1.0,momentq0(l)+cos(phi(l)+qvector*2.0*3.141529*rr)))
//				                    )                    ;	
		  mom[in(rra,rrb,rrc)](l)=nettom(l)*qv*(0.01+
		  momentq0(l)+cos(phi(l)+qvector*2.0*3.141529*rr))
				                                        ;	
                  }
	         }
	}
     }
        
}




// load spinconfiguration from file returns 1 on success and 0 on failure
int spincf::load(FILE * fin_coq)	
{ int i,j,k,l,nn1,nn2;  
  char * s;
  float na[MAXNOFSPINS];
  float nb[MAXNOFSPINS];
  float nc[MAXNOFSPINS];
  na[0]=MAXNOFSPINS;
  nb[0]=MAXNOFSPINS;
  nc[0]=MAXNOFSPINS;
  long int pos;
  char instr[MAXNOFCHARINLINE];
  // input comment lines 
  instr[0]='#';
   while (instr[strspn(instr," \t")]=='#') 
  {pos=ftell(fin_coq);if (pos==-1) return 0; //end of file reached
   s=fgets(instr,MAXNOFCHARINLINE,fin_coq);if (s==NULL) return 0;
  }
  j=fseek(fin_coq,pos,SEEK_SET); if (j!=0) return 0;
nn1=0;
j=0;
while ((i=inputline (fin_coq, na))!=0)  
   {for (l=1;l<=nofcomponents*nofatoms-2;++l){inputline (fin_coq, nb);}
    nn1=i;i=inputline (fin_coq, nc);++j;}
// j (nofspins in b direction)determined
if (nn1==0) return 0; // no block to read any more 
nn2=j;

//how much is nofc ??
for (k=0;j==nn2;++k)
  {for (j=0;(i=inputline(fin_coq,na))==nn1;++j)
   {for (l=1;l<=nofcomponents*nofatoms-2;++l){inputline (fin_coq, nb);}
    i=inputline (fin_coq, nc);}
  }

nofa=nn1;
nofb=nn2;    
nofc=k;

//  printf("n1=%i n2=%i n3=%i\n",nofa,nofb,nofc);
  
  j=fseek(fin_coq,pos,SEEK_SET); if (j!=0) return 0;

  delete []mom;
  mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
  
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];for(l=0;l<=mxa*mxb*mxc;++l){mom[l]=Vector(1,nofcomponents*nofatoms);}
  if (mom == NULL)
    {fprintf (stderr, "Out of memory\n");
     exit (EXIT_FAILURE);} 

 for (k=1;k<=nofc;++k)
 {for (j=1;j<=nofb;++j)
  {for (l=1;l<=nofcomponents*nofatoms;++l)
   {inputline (fin_coq, na);
    for (i=1;i<=nofa;++i)
       {mom[in(i,j,k)](l)=na[i];}
    }
   }
  if(inputline(fin_coq, na))
  {fprintf(stderr,"load spincf error - no  empty line between [001] planes");exit(EXIT_FAILURE);}  
  }       
  if(inputline(fin_coq, na))
  {fprintf(stderr,"load spincf error - no  2 empty lines after spinconfiguration");exit(EXIT_FAILURE);}  
         
 return 1;
}

// take vector dd and calculate nearest atom in spinconfiguration
double spincf::nndist(float * x, float * y, float * z,Vector & abc,Matrix & p,Vector &dd)
{int i,j,k,l;double d,mindist=1e10;Vector ddl(1,3);Vector dd0(1,3);
 Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,abc);

 for(i=-1;i<=1;++i)for(j=-1;j<=1;++j)for(k=-1;k<=1;++k)for(l=1;l<=nofatoms;++l)
 {dd0(1)=x[l];
  dd0(2)=y[l];
  dd0(3)=z[l];
  ddl=abc_in_ijk*dd0;
  ddl+=p.Column(1)*((double)(i-1)/nofa)+p.Column(2)*((double)(j-1)/nofb)+p.Column(3)*((double)(k-1)/nofc);
  d=Norm(dd-ddl);if(d>0.0001&&d<mindist)mindist=d;
 }
return mindist;
}

Vector spincf::pos(int i, int j, int k, int l,Vector & abc,Matrix & r,float * x,float *y,float*z)
{//returns position 
 // of atom l at lattice site (i j k) (Angstrom) as vector components in Euclidean ijk coordinate system
 // where i j k are defined by  j||b, k||(a x b) and i normal to k and j
Vector dd(1,3),dd0(1,3);
if(l>nofatoms){fprintf(stderr,"Error spincf.cpp::pos - index l > nofatoms\n");exit(1);}
  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);
  Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,abc);
         dd0(1)=x[l];
         dd0(2)=y[l];
         dd0(3)=z[l];
         dd=abc_in_ijk*dd0;
         dd+=p.Column(1)*((double)(i-1)/nofa)+p.Column(2)*((double)(j-1)/nofb)+p.Column(3)*((double)(k-1)/nofc);
return dd;
}
Vector spincf::pos(int i, int j, int k, int l,cryststruct & cs)
{if(cs.nofatoms!=nofatoms){fprintf(stderr,"Error spincf.cpp::pos nofatoms=%i not equal to cs.nofatoms=%i\n",nofatoms,cs.nofatoms);exit(1);}
return pos(i,j,k,l,cs.abc,cs.r,cs.x,cs.y,cs.z);

}
    
Vector spincf::pos_dabc(int i, int j, int k, int l,cryststruct & cs)
                      //returns position of atom l at lattice site (i j k) 
                      // as vector components  refering to lattice vectors abc
{Vector dd(1,3),dd0(1,3);
 if(cs.nofatoms!=nofatoms){fprintf(stderr,"Error spincf.cpp::pos_dabc nofatoms not equal to cs.nofatoms\n");exit(1);}
        dd0=pos(i,j,k,l,cs.abc,cs.r,cs.x,cs.y,cs.z);
         Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,cs.abc);
         dd=abc_in_ijk.Inverse()*dd0;
 return dd;
}    
Vector spincf::pos_dr123 (int i, int j, int k, int l,cryststruct & cs)
                      //returns position of atom l at lattice site (i j k) as
                      // vector components refering to primitive lattice vectors r1 r2 r3
{Vector dd(1,3),dd0(1,3);
  Matrix p(1,3,1,3);
if(l>nofatoms){fprintf(stderr,"Error spincf.cpp::pos - index l > nofatoms\n");exit(1);}
if(cs.nofatoms!=nofatoms){fprintf(stderr,"Error spincf.cpp::pos_dr123 nofatoms not equal to cs.nofatoms\n");exit(1);}

  calc_prim_mag_unitcell(p,cs.abc,cs.r);
  Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,cs.abc);
         dd0(1)=cs.x[l];
         dd0(2)=cs.y[l];
         dd0(3)=cs.z[l];
         dd=abc_in_ijk*dd0;
         dd+=p.Column(1)*((double)(i-1)/nofa)+p.Column(2)*((double)(j-1)/nofb)+p.Column(3)*((double)(k-1)/nofc);
         dd0=p.Inverse()*dd;dd0(1)*=nofa;dd0(2)*=nofb;dd0(3)*=nofc;
return dd0;
}
/**************************************************************************/

//zuweisung
spincf & spincf::operator= (const spincf & op2)
{if(this!=&op2){ // this is to avoid problems when copying to itself !!
 int i,j,k;
 nofa=op2.nofa; nofb=op2.nofb; nofc=op2.nofc;
 mxa=op2.mxa; mxb=op2.mxb; mxc=op2.mxc;
 nofatoms=op2.nofatoms;
 wasstable=op2.wasstable;
 nofcomponents=op2.nofcomponents;
  epsilon=Vector(1,6);
  epsilon=op2.epsilon;
  delete []mom;
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];for(k=0;k<=mxa*mxb*mxc;++k){mom[k]=Vector(1,nofcomponents*nofatoms);}
  if (mom == NULL)
    {fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);}
  for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {mom[in(i,j,k)]=op2.mom[in(i,j,k)];} 
    }
  }  
}
  return *this;
}

//addition
void spincf::add (spincf & op1, const spincf & op2)
{int i,j,k;
 if (op1.nofa!=op2.nofa||op1.nofb!=op2.nofb||op1.nofc!=op2.nofc||op1.nofatoms!=op2.nofatoms||op1.nofcomponents!=op2.nofcomponents)
 {fprintf (stderr,"Error in adding spincfonfigurations - not equal dimension\n");
 exit (EXIT_FAILURE);}
 for (i=1;i<=op1.nofa;++i)
  {for (j=1;j<=op1.nofb;++j)
    {for (k=1;k<=op1.nofc;++k)
     {op1.mom[op1.in(i,j,k)]+=op2.mom[op1.in(i,j,k)];} 
    }
  }           
}

spincf spincf::operator + (const spincf & op2)
{
    spincf copy = *this;
    add(copy, op2);
    return copy;
}

spincf & spincf::operator += (const spincf & op2) {
     add(*this, op2); 
     return (*this);
}

// multiplication of spinconfiguration with constant
void spincf::multiply (spincf& op1, const double factor)
{int i,j,k;
  for (i=1;i<=op1.nofa;++i)
  {for (j=1;j<=op1.nofb;++j)
    {for (k=1;k<=op1.nofc;++k)
     {op1.mom[op1.in(i,j,k)] *= factor;
      } 
    }
  }  
}

spincf spincf::operator * (const double factor)
{
    spincf copy = *this;
    multiply(copy, factor);
    return copy;
}

spincf & spincf::operator *= (const double factor) {
     multiply(*this, factor); 
     return (*this);
}



//vergleichsoperator
int spincf::operator==(spincf & op1)
{//untersuchen der periodizitaet
int nna,nnb,nnc,chka,chkb,chkc,ma,mb,mc;
if (nofatoms!=op1.nofatoms||nofa!=op1.nofa||nofb!=op1.nofb||nofc!=op1.nofc) return 0;
// vergleichen der spins unter zuhilfename von spequal
for (chka=1;chka<=nofa;++chka)
 {for (chkb=1;chkb<=nofb;++chkb)
    {for (chkc=1;chkc<=nofc;++chkc)
    
     {for (nna=1;nna<=nofa;++nna)
         {for (nnb=1;nnb<=nofb;++nnb)
	     {for (nnc=1;nnc<=nofc;++nnc)
              {ma=mod(nna+chka,nofa);
	       mb=mod(nnb+chkb,nofb);
	       mc=mod(nnc+chkc,nofc);
               if (!spequal(mom[in(nna,nnb,nnc)],op1.mom[op1.in(ma,mb,mc)]))
                  {goto nextchk;}
	       if (nna==nofa&&nnb==nofb&&nnc==nofc) return 1;
               }
	      }
	  }     	  
nextchk: ;
     }
    }
  }    
   
return 0; //not equal
}

//constructors
//from n1xn2xn3 unit cells with na number in the cryst. basis and nc number of components of the spin of each atom
spincf::spincf (int n1,int n2,int n3,int na,int nc)
{ wasstable=0;
  int l;epsilon=Vector(1,6);epsilon=0;
  nofa=n1;nofb=n2;nofc=n3;
   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
  nofatoms=na;
  nofcomponents=nc;
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];
  for(l=0;l<=mxa*mxb*mxc;++l){mom[l]=Vector(1,nofcomponents*nofatoms);mom[l]=0;}
  if (mom == NULL)
    { fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);} 
}

//kopier-konstruktor
spincf::spincf (const spincf & p)
{ int i,j,k;
  nofa=p.nofa;nofb=p.nofb;nofc=p.nofc;
  mxa=p.mxa; mxb=p.mxb; mxc=p.mxc;
  wasstable=p.wasstable;
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  epsilon=Vector(1,6);
  epsilon=p.epsilon;
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];for(k=0;k<=mxa*mxb*mxc;++k){mom[k]=Vector(1,nofcomponents*nofatoms);}
  if (mom == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    } 
 for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {mom[in(i,j,k)]=p.mom[in(i,j,k)];} 
    }
  }           

}


//destruktor
spincf::~spincf ()
{//printf("hello destruktor spincf\n");  
  delete []mom;
//printf("hello destruktor spincf\n");  
}


