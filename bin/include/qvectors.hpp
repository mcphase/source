// class for the generation and storage of the set
// of test q vectors used in mcphase

#ifndef QVECTORS
#define QVECTORS


#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<cmath>
#include<ctime>
#include<martin.h>
#include<spincf.hpp>
#include<vector.h>
#include<inipar.hpp>

class qvectors
{
  private:
  
    char * savfilename;   
    float *hchk[4]; // variable to store hkl
    int *hchkz[4];
    int *hchkn[4];
    Vector hkl; //variable to return vector
    int nofq; // number of qvectors
    Matrix rez,r;
    
    int ia(int j); //index functions
    int ib(int j);
    int ic(int j);
    Vector **q0; // array of pointers to the q-vectors
    Vector **n; // array of pointers to the periodicity vectors (na,nb,nc)... in case of memory
	        // problems this could be subsituted by integer fields
    Vector **nm; //array of pointers to the vectors  storing nettom
    Vector **mq0; //variable to store momentq0
    Vector **ph; //variable to store phi
    int maxnofqs(); // returns maximal number of qvectors
        
  public:
    int nofqs (); //returns nofqvectors
    Vector & q(int i); // returns (hkl) (i)
    int na (int i); // returns period for i.th qvector
    int nb (int i); // returns period for i.th qvector
    int nc (int i); // returns period for i.th qvector
    int nofatoms; //number of atoms in primitive cryst unit cell
    int nofcomponents; //number of moments in spin vector
    int verbose;  /*    switch: if 1 the generation of q vectors is commented explicitely on stdout
	                and qvectors::save(const char * filemode) puts a large set of information into mcphas.qom (including
			magnetic structures corresponding to qvector)*/

    Vector & nettom(int i); // nettomoment
    Vector & momentq0(int i); // moment amplitude
    Vector & phi(int i); // phase
       
    //save table of all qvectors on file
    void save(const char * filemode);
    
    //constructor - generate set of qvectors
    /* input 
	ini		set of initial parameters of mcphase - such as Vectors defining a cube 
	                in qspace which is to be covered,Vector defining a increment in qspace
	filename	file to save q vector values on
        rz              reciproval lattice
	A,B,C           maximum values of the moments in a b c direction
	verbose         switch: if 1 the generation of q vectors is commented explicitely on stdout
	                and qvectors::save(char * filemode) puts a large set of information into mcphas.qom (including
			magnetic structures corresponding to qvector)
    */
   
 qvectors (inipar & ini,Matrix & rz,
              Vector & mmax,const char * filename,int nofatoms,int nofcomponents,int v);	//konstruktor

    qvectors (const qvectors & qs);	// kopier-konstruktor

    
~qvectors ();		//destruktor

};


#endif
