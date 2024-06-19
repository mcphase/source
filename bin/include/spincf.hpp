//class for the generation and storage of a spinconfiguration
// used in mcphase


#include<cstdlib>
#include<cerrno>
#include<vector.h>
#include<complex>
#include<cstdio>
#include<graphic_parameters.hpp>
#include<cryststruct.hpp>

#ifndef SPINCF_H
#define SPINCF_H

class spincf 
{
  // OUTPUT to FILE ... in spincf_out.cpp --------------------------------
public:
    Vector epsilon; // strain
    void print(FILE * fout);
    void print(FILE * fout, int nofcomp);
    void print_commented(FILE * fout,const char * string,int min, int max, int maxnofpars, double & absvallimit);
//print list of atoms + positions + moments
    void eps(FILE * fout);
    void eps(FILE * fout,const char * text);
    void eps3d(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z,int orientation,spincf & magmom);
    void fst(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z,spincf & magmom);

              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase

    void jvx_cd(FILE * fout,char * text,cryststruct & cs,
              graphic_parameters & gp,double phase,spincf & savev_real,spincf & savev_imag,
              Vector & hkl,double & T, Vector &  gjmbHxc,Vector & Hext,spincf & magmom,spincf & magmomev_real, spincf & magmomev_imag);
    void jvx_cd(FILE * fout,char * text,cryststruct & cs,
              graphic_parameters & gp, double phase, spincf  savev_real, spincf  savev_imag,
              Vector & hkl, double & T, Vector &  gjmbHxc, Vector & Hext, cryststruct & cs4, spincf  magmom, spincf magmomev_real, spincf  magmomev_imag);
    void jvx_cd(FILE * fout,char * text,cryststruct & cs,
              graphic_parameters & gp,double phase,spincf & savev_real,spincf & savev_imag,
              Vector & hkl,double & T, Vector &  gjmbHxc,Vector & Hext,cryststruct & cs4,spincf & magmom,spincf & magmomev_real, spincf & magmomev_imag,spincf & pev_real, spincf & pev_imag);

    void cd(FILE * fout,cryststruct & cs,graphic_parameters & gp,
                spincf & savev_real,spincf & savev_imag,double phase,Vector & hkl,double & T,Vector &  gjmbHxc,Vector & Hext);

    void fstprim(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, spincf & magmom);
    void calc_prim_mag_unitcell(Matrix & p,Vector & abc, Matrix & r);
private:
 // frame of display
   void epsarrow(FILE * fout,Vector a,Vector b);
   Vector xy(Vector xyz,int orientation,Vector min,Vector max,float bbwidth,float bbheight);
   void calc_minmax(Vector & min,Vector & max,Vector & ijkmin,Vector & ijkmax,Matrix & p,Vector & abc);
   void calc_minmax_scale(Vector & min,Vector & max,Vector & ijkmin,Vector & ijkmax,Matrix & p,Vector & abc,double scale_view_1,double scale_view_2,double scale_view_3);
    void calc_prim_mag_unitcell_old(Matrix & p,Vector & abc, Matrix & r);
// ----------------------------------------------------------


  private:
 // number of spins  
   int nofa,nofb,nofc;
 // this subtracts n2 if n1>n2
   int mod(int n1,int n2);
   int mxa,mxb,mxc;
   Vector * mom; // momentums <J>
   int iv[4];
   int spequal(Vector a,Vector b);// routine to compare spins

   static void multiply(spincf& op1, const double factor);
   static void add(spincf& op1, const spincf & op2);
     
   // take vector dd and calculate distance nearest atom in spinconfiguration
   double nndist(float * x, float * y, float * z,Vector & abc,Matrix & p,Vector &dd);
    Vector pos(int i, int j, int k, int l,Vector & abc,Matrix & r,float * x,float *y,float*z);
                      //returns position of atom l at lattice site (i j k) (Angstrom)
                      // as vector components in Euclidean ijk coordinate system
                      // defined by  j||b, k||(a x b) and i normal to k and j
  
 public:
    Vector moment(int i,int j,int k,int l); // returns moment of atom l (1,nofcomponents)

    Vector pos(int i, int j, int k, int l,cryststruct & cs);
                      //returns position of atom l at lattice site (i j k) (Angstrom)
                      // as vector components in Euclidean ijk coordinate system
                      // defined by  j||b, k||(a x b) and i normal to k and j

    Vector pos_dabc(int i, int j, int k, int l,cryststruct & cs);
                      //returns position of atom l at lattice site (i j k) 
                      // as vector components  refering to lattice vectors abc
    Vector pos_dr123 (int i, int j, int k, int l,cryststruct & cs);
                      //returns position of atom l at lattice site (i j k) as
                      // vector components refering to primitive lattice vectors r1 r2 r3


    int  load(FILE * fin_coq);	// load spincf from file returns 1 on success and 0 on failure
 // array of spins 
   int in(int i,int j, int k); 
    int wasstable; // index to remember if it was stable: if a sinconfiguration is set stable, its periodicity key is stored in wasstable
   
    int nofatoms;
    int nofcomponents;
    Vector & m(int i,int j,int k); // returns pointer to spin (ijk) 
    Vector & mi(int in); // returns pointer to spin i
    void  FT(ComplexVector * mq); // returns Fourier transform mq of spins (for use see htcalc.c)
    int * ijk(int in);  // returns spin indizes (ijk)(in): in=0,...,n(=na*nb*nc)
    
    int n(); // returns total number of primitive crystal basis in supercell
    int na(); // returns number of primitive crystal basis in supercell along a
    int nb(); // returns number of  primitive crystal basis in supercell along b
    int nc(); // returns number of  primitive crystal basis in supercell along c
 
    Vector totalJ (); // returns nettomoment <J>
    void invert();// inverts all spins (AND higher order moments)
    int reduce();// reduces spinconfiguration, if reduction is possible returns 1, otherwise 0
    void spinfromq (int n1,int n2, int n3,Vector & qvector, Vector & nettom,Vector & momentq0, Vector & phi);


    spincf operator + (const spincf & op2); // addition    
    spincf & operator += (const spincf & op2); // addition    
    spincf operator * (const double factor); // multiplication with constant
    spincf & operator *= (const double factor); 
    spincf & operator= (const spincf & op2); // zuweisung
    int operator== (spincf & op2); // vergleich

   
spincf (int n1=1,int n2=1,int n3=1,int nofatoms=1,int nofcomponents=3);
                                   //konstruktor mit initialisierung (wenn noetig)
   spincf (const spincf & spins);	// kopier-konstruktor
   
~spincf ();		//destruktor

};

#endif
