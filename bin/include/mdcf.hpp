// class mdcf to store single ion transition matrix elements etc. - for use in mcdisp

#include<cstdlib>
#include<cerrno>
#include<vector.h>
#include<complex>

#ifndef MDCF_H
#define MDCF_H

class mdcf 
{
  private:
 // mode -r switch to determine if Matrix m needs to be initalised or not 
 // (in case -r option is not set memory can be saved by not using memory for this matrices)
   int mr;
 // number of spins  
   int nofa,nofb,nofc;
   int mxa,mxb,mxc;
   int mqdim; // dimension of transition operator (unpol neutrons M(Q): dim=3, rixs: some tensor for e e' components of polarisation i.e. dim =9)
   ComplexMatrix ** m; //matrix to store M
   ComplexMatrix ** s; //matrix to store U
   ComplexVector ** l; //matrix to store eigenvalues sqrt_gamma
   // ComplexMatrix ** sb; //matrix to store V
   ComplexVector ** dps; // big vector to store all dmq1 
   ComplexVector ** dmqs; // big vector to store all dmq1 
   ComplexVector ** dmq_dips; // big vector to store all dmq_dip1 
   ComplexVector ** lb; //matrix to store eigenvalues sqrt_gamma
   ComplexVector ** Pb; //matrix to store eigenvalues sqrt_gammaP
   ComplexVector ** lb_dip; //matrix to store eigenvalues sqrt_gamma
   Vector ** d;
   IntVector ** nt; // vector to store number of transitions for each atom
   ComplexMatrix ** eigenstates; // matrix to store the eigenstates of ions
   ComplexMatrix *** chi0s; // matrix to store the chi0 of ions
   int iv[4];
   void errexit();

  public:
 // array of spins 
    int nofatoms,nofcomponents;
    int in(int i,int j, int k) const; // indexing functions
    int ind(int i,int j, int k,int l); 
    int inM(int i,int j, int k,int l); 
    int nofEstps;
    int ncel;
    int storage; // indicates if storage except nt should be created !
    ComplexMatrix **Ug,**gU, **bUg,**bgU,**PUg,**PgU; // Cache for U*sqrt(gamma) and sqrt(gamma)*U values, and beyond equiv.
    
    ComplexMatrix & M(int i,int j,int k); // returns pointer to  matrix M(ijk) 
    ComplexMatrix & Mi(int in); // returns pointer to matrix M(i)
    ComplexMatrix & U(int i,int j,int k) const; // returns pointer to eigenvector matrix (ijk) 
    ComplexMatrix & Ui(int in); // returns pointer to eigenvector matrix i
    ComplexVector & sqrt_gamma(int i,int j,int k) const; // returns pointer to eigenvaluematrix (ijk) 
    ComplexVector & sqrt_gammai(int in); // returns pointer to eigenvalue matrix i
  //  ComplexMatrix & V(int i,int j,int k) const; // returns pointer to eigenvector matrix (ijk) 
  //  ComplexMatrix & Vi(int in); // returns pointer to eigenvector matrix i
    ComplexVector & dPs(int i,int j,int k) const; // returns pointer to vector P (ijk) 
    ComplexVector & dMQs(int i,int j,int k) const; // returns pointer to vector MQ (ijk) 
    ComplexVector & dMQ_dips(int i,int j,int k) const; // returns pointer to vector MQ_dip (ijk) 
    ComplexVector & sqrt_GammaP(int i,int j,int k) const; // returns pointer to eigenvaluevector (ijk) 
    ComplexVector & sqrt_GammaPi(int in); // returns pointer to eigenvalue vector i
    ComplexVector & sqrt_Gamma(int i,int j,int k) const; // returns pointer to eigenvaluevector (ijk) 
    ComplexVector & sqrt_Gammai(int in); // returns pointer to eigenvalue vector i
    ComplexVector & sqrt_Gamma_dip(int i,int j,int k) const; // returns pointer to eigenvaluevector (ijk) 
    ComplexVector & sqrt_Gamma_dipi(int in); // returns pointer to eigenvalue vector i
    ComplexMatrix & est(int i, int j, int k, int l); // returns pointer to eigenstate matrix for atom ijkl
    ComplexMatrix ** chi0pointer(int i, int j, int k, int l); // returns pointer to chi0matrices for atom ijkl
    void est_ini(int i, int j, int k, int l,ComplexMatrix & M); // initialize est

    Vector & delta(int i,int j,int k); // returns pointer to matrix (ijk) 
    Vector & deltai(int in); // returns pointer to mean field i
    int * ijk(int in);  // returns mf indizes (ijk)(in): in=0,...,n(=na*nb*nc)
    
    int n(); // returns total number of spins
    int na(); // returns number of spins
    int nb(); // returns number of spins
    int nc(); // returns number of spins
    int baseindex (int i, int j, int k, int l, int tn) const; // returns base index for atom l and transition tn
                                                       // which is needed for  setting up matrix U,M, sqrt_gamma and Vector D
    int baseindex_max(int i, int j, int k) const;
    int noft(int i, int j, int k, int l) const; // returns number of transitions of ion l in cryst unit ijk
   
mdcf (int n1,int n2,int n3,int n,int nc,int nstps,int do_Erefine);	//konstruktor

// initialisierung 
    void set_noftransitions (int i, int j, int k, IntVector & notr,int mqd);

     mdcf (const mdcf & spins,int store);	// kopier-konstruktor
   
~mdcf ();		//destruktor

};

#endif
