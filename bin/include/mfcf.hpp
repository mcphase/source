// klasse testmfcf to store mean field configuration

#include<cstdlib>
#include<cerrno>
#include<vector.h>
#include<complex>

#ifndef MFCF_H
#define MFCF_H

class mfcf 
{
  private:
 // number of spins  
   int nofa,nofb,nofc;
 // this subtracts n2 if n1>n2
   int mod(int n1,int n2);
   int mxa,mxb,mxc;
 // frame of display
   Vector * mfi; // mean fields = gjmbHeff [meV]
   int iv[4];
          
  public:
   int nofatoms;
   int nofcomponents;
    Vector epsmf;
   void resetnofc(int n);
 // array of spins 
   int in(int i,int j, int k); 
    int wasstable; // index to remember if it was stable
   
    Vector & mf(int i,int j,int k); // returns pointer to mean field (ijk) 
    Vector & mi(int in); // returns pointer to mean field i
    int * ijk(int in);  // returns mf indizes (ijk)(in): in=0,...,n(=na*nb*nc)
    
    int n() const; // returns total number of spins
    int na() const; // returns number of spins
    int nb() const; // returns number of spins
    int nc() const; // returns number of spins
    void invert();// inverts all spins

    void print(FILE * fout);
    void print_human_readable(FILE * fout);
    int  load(FILE * fin_coq);	
     
    void clear(); // set all meanfields to zero    
   mfcf & operator= (const mfcf & op2); // zuweisung
   //zuweisung of the same meanfield vector to all atoms
   mfcf & operator= (const Vector & vec);
   
   mfcf (int n1=1,int n2=1,int n3=1, int nofatoms=1,int nofcomponents=3);	//konstruktor mit initialisierung (wenn noetig)
   mfcf (const mfcf & spins);	// kopier-konstruktor
   
~mfcf ();		//destruktor

};

#endif
