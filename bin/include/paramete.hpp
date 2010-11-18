// klasse parameters zum einlesen der parameter
#ifndef PARAMETE
#define PARAMETE


#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<ctime>
#include<martin.h>
#include<vector.h>

class parameters
{
  private:
  
  float *rij;
  Vector  *jijerr;
  char *rems[40];
   
  public:
  int nofbcplanes;
  int paranz;
  Vector *jj,*jij,*dn; // denotes coupling between bc planes  
   
  //lattice
  float a,b,c,alpha,beta,gamma;
  // switches
  int diagonalexchange;
  // lande factor and momentum factor
  float gJ,J;
  float A, M, C;
  float m0, alphaa, alphab, alphac, tferro, ex;
  float m0xa, m0xb, m0xc, sta, staerr;   
  float tpara, dmag, tneel, hstar;

   int save(FILE * file); // save parameters to file (coq format)
   int savjjj(FILE * file); // save exchange parameters to file (jjj format)
   int savcf(FILE * file); // save cf parameters to file (cf format)
   
   parameters (char *file);	//konstruktor
   parameters (const parameters & pars);	// kopier-konstruktor
   ~parameters ();		//destruktor
};

#endif
