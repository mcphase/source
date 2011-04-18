// classe jq to store fourier transform of exchange
#include<cstdlib>
#include<cerrno>
#include<vector.h>
#include<mdcf.hpp>
#include<complex>

#ifndef JQ_H
#define JQ_H

class jq
{
  private:
 // number of spins  
   int nofa,nofb,nofc;
   int mxa,mxb,mxc,mx;
   ComplexMatrix ** jj;
   int iv[4];
  // iindex
    int iin(int i,int j);
       
  public:
  // index
    int in(int i,int j, int k); 
    int nofatoms,nofcomponents;
       
    ComplexMatrix & mat(int i,int j,int k,int i1,int j1,int k1); // returns pointer to matrix (ijk)(i1,j1,k1) 
    ComplexMatrix & mati(int in,int jn); // returns pointer to matrix ij
    int * ijk(int in);  // returns mf indizes (ijk)(in): in=0,...,n(=na*nb*nc)
    
    int n(); // returns total number of spins
    int na(); // returns number of spins
    int nb(); // returns number of spins
    int nc(); // returns number of spins
  

    jq (int n1,int n2,int n3,mdcf & m);	//konstruktor 
    jq (const jq & spins);	// kopier-konstruktor
    ~jq ();		//destruktor

};

#endif
