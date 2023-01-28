 // *************************************************************************
 // ************************ meanfieldcf *************************************
 // *************************************************************************
// methods for class mfcf 
#include <cerrno>
#include <cstdio>
#include <cmath>
#include <martin.h>
#include <vector.h>
#include <complex>
#include "mfcf.hpp"

#define MAXNOFSPINS  200

// returns mf of spin [i=na,j=nb,k=nc] 
Vector & mfcf::mf(int na, int nb, int nc)
{ return mfi[in(na,nb,nc)];
}

// the same but for spin number "i"
Vector & mfcf::mi(int i)
{ return mfi[i];
}

// get index ijk=iv(1-3)  of mfonfiguration number in 
int * mfcf::ijk(int in)
{div_t result; result=div(in,mxb*mxc); 
 iv[1]= result.quot;
 result=div(result.rem,mxc);
 iv[2]= result.quot;
 iv[3]= result.rem;
 return iv;}

// the inverse: get number of spin from indizes i,j,k
int mfcf::in(int i, int j, int k)
{return ((i*mxb+j)*mxc+k);}

// this subtracts n2 if n1>n2
int mfcf::mod(int n1,int n2)
{if (n1>n2) return n1-n2;
 else return n1;
}

// return number of spins
int mfcf::n() const
{return (nofa*nofb*nofc);
}
int mfcf::na() const
{return nofa;
}
int mfcf::nb() const
{return nofb;
}
int mfcf::nc() const
{return nofc;
}


// invert all mean fields
void mfcf::invert()
{
 int i,j,k; 
 for (i=1;i<=nofa;++i)
 { for (j=1;j<=nofb;++j)
   {for (k=1;k<=nofc;++k)
    {mfi[in(i,j,k)]= -mfi[in(i,j,k)];}
   }
 }  
}

// load mfconfiguration from file
int mfcf::load(FILE * fin_coq)	
{ int i,j,k,l,nn1,nn2;  
  float na[MAXNOFSPINS];
  float nb[MAXNOFSPINS];
  float nc[MAXNOFSPINS];
  na[0]=MAXNOFSPINS;
  nb[0]=MAXNOFSPINS;
  nc[0]=MAXNOFSPINS;
  long int pos;
  
  pos=ftell(fin_coq);if (pos==-1) return 0;

nn1=0;j=0;
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


  j=fseek(fin_coq,pos,SEEK_SET); if (j!=0) return 0;

  delete []mfi;
  mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
  
//dimension arrays
  mfi = new Vector[mxa*mxb*mxc+1];for(i=0;i<=mxa*mxb*mxc;++i){mfi[i]=Vector(1,nofcomponents*nofatoms);}
  if (mfi == NULL)
    {fprintf (stderr, "Out of memory\n");
     exit (EXIT_FAILURE);} 

 for (k=1;k<=nofc;++k)
 {for (j=1;j<=nofb;++j)
  {for (l=1;l<=nofcomponents*nofatoms;++l)
   {inputline (fin_coq, na);
    for (i=1;i<=nofa;++i)
       {mfi[in(i,j,k)](l)=na[i];}
    }
   }
  if(inputline(fin_coq, na))
  {fprintf(stderr,"load mfcf error - no  empty line between [001] planes");exit(EXIT_FAILURE);}  
  }       
  if(inputline(fin_coq, na))
  {fprintf(stderr,"load mfcf error - no  empty line between at end of configuration");exit(EXIT_FAILURE);}  
         
 return 1;
}


//-----------------------------------------------------------------------
//  numeric output of mfconfiguration to file
void mfcf::print(FILE * fout) //print mfconfiguration to stream
{int i,j,k,l;
 for (k=1;k<=nofc;++k)
 {for (j=1;j<=nofb;++j)
  {for (l=1;l<=nofcomponents*nofatoms;++l)
   {for (i=1;i<=nofa;++i)
      {fprintf(fout," %4.4f",myround(1e-5,mfi[in(i,j,k)](l)));
       }
    fprintf(fout,"\n");
    }
   }
 fprintf(fout,"\n"); //new line to separate ab planes
 }
// fprintf(fout,"\n"); //new line to end configuration - removed aug 07
               
}

/**************************************************************************/

//zuweisung
mfcf & mfcf::operator= (const mfcf & op2)
{int i,j,k;
 nofa=op2.nofa; nofb=op2.nofb; nofc=op2.nofc;
 mxa=op2.mxa; mxb=op2.mxb; mxc=op2.mxc;
 nofatoms=op2.nofatoms;
 nofcomponents=op2.nofcomponents;
 epsmf=Vector(1,6);epsmf=op2.epsmf;
 wasstable=op2.wasstable;
  delete []mfi;
//dimension arrays
  mfi = new Vector[mxa*mxb*mxc+1];for(i=0;i<=mxa*mxb*mxc;++i){mfi[i]=Vector(1,nofcomponents*nofatoms);}
  if (mfi == NULL)
    {fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);}
  for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {mfi[in(i,j,k)]=op2.mfi[in(i,j,k)];} 
    }
  }           
  return *this;
}

//zuweisung of the same meanfield vector to all atoms
mfcf & mfcf::operator= (const Vector & vec)
{int i,j,k;
  for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {mfi[in(i,j,k)]=vec;} 
    }
  }         
  return *this;
}

//set all meanfields zero
void mfcf::clear()
{int i,j,k;
  for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {mfi[in(i,j,k)]=0;} 
    }
  }           
}


//constructors
mfcf::mfcf (int n1,int n2,int n3,int na,int nc)
{ wasstable=0;
  nofa=n1;nofb=n2;nofc=n3;
   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
  nofatoms=na;
  nofcomponents=nc;
  epsmf=Vector(1,6);
  int i;
//dimension arrays
  mfi = new Vector[mxa*mxb*mxc+1];for(i=0;i<=mxa*mxb*mxc;++i){mfi[i]=Vector(1,nofcomponents*nofatoms);}
  if (mfi == NULL)
    { fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);} 
}

void mfcf::resetnofc(int n)
{Vector * mfn;
 //dimension arrays
  mfn = new Vector[mxa*mxb*mxc+1];
  if (mfn == NULL)
    { fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);} 
  for(int i=0;i<=mxa*mxb*mxc;++i){mfn[i]=Vector(1,n*nofatoms);mfn[i]=0;
   for(int k=1;k<=nofatoms;++k){
    for(int l=1;l<=n&&l<=nofcomponents;++l){
    mfn[i]((k-1)*n+l)=mfi[i]((k-1)*nofcomponents+l);
   }}
  }
  delete []mfi;
  mfi=mfn;  
  nofcomponents=n;
}

//kopier-konstruktor
mfcf::mfcf (const mfcf & p)
{ int i,j,k;
  nofa=p.nofa;nofb=p.nofb;nofc=p.nofc;
  mxa=p.mxa; mxb=p.mxb; mxc=p.mxc;
  wasstable=p.wasstable;
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  epsmf=Vector(1,6);epsmf=p.epsmf;
//dimension arrays
  mfi = new Vector[mxa*mxb*mxc+1];for(i=0;i<=mxa*mxb*mxc;++i){mfi[i]=Vector(1,nofcomponents*nofatoms);}
  if (mfi == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    } 
 for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {mfi[in(i,j,k)]=p.mfi[in(i,j,k)];} 
    }
  }           

}


//destruktor
mfcf::~mfcf ()
{
  delete []mfi;
}


