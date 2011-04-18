/* states.hpp
 * 
 * Header file for states.cpp
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#ifndef STATES_H
#define STATES_H

#include<vector>
#include<iostream>
#include<sstream>

// --------------------------------------------------------------------------------------------------------------- //
// Defines an enumeration to contain both the spectroscopic notation and integer value of the orbital quantum number
// --------------------------------------------------------------------------------------------------------------- //
enum orbital {  S=0,		//	 0 = sharp
		P,		//       1 = principal
		D,		//	 2 = diffuse
		F,		//	 3 = fine
		Fp=-3,		//	 3		// For states where S,L,v,U not enough quantum numbers
		G=4,		//	 4
		Gp=-4,		//	 4		// For states where S,L,v,U not enough quantum numbers
		H=5,		//	 5
		Hp=-5,		//	 5		// For states where S,L,v,U not enough quantum numbers
		I=6,		//	 6
		Ip=-6,		//	 6		// For states where S,L,v,U not enough quantum numbers
		K=7,		//	 7
		Kp=-7,		//	 7		// For states where S,L,v,U not enough quantum numbers
		L=8,		//	 8
		Lp=-8,		//	 8		// For states where S,L,v,U not enough quantum numbers
		M=9, 		//	 9
		N,		//	10
		O,		//	11
		Q  };		//	12

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in states.cpp, which are not constructors/members of any class below
// --------------------------------------------------------------------------------------------------------------- //
std::string strtoupper(std::string instring);           // Converts the characters of a string to upper case
std::string Lout(const orbital & l);                    // Converts the enumerated value of l into a string
orbital Lin(std::string & l);                           // Converts a string of the spectroscopy notation to l
// --------------------------------------------------------------------------------------------------------------- //
// qR7 racah_vtow(int S2, int v);                       // Converts a seniority number v to W quantum numbers.
// int racah_wtov(int S2, qR7 W);                       // Converts a (W1 W2 W3) quantum number to seniority num.
// --------------------------------------------------------------------------------------------------------------- //

// --------------------------------------------------------------------------------------------------------------- //
// Defines a class to contain and manipulate the three numbers that label the irreducible representation of R7, 
//    after Racah. These are the W quantum numbers that help label states of f-configurations.
// --------------------------------------------------------------------------------------------------------------- //
class qR7 
{
   private:

   public:
      int w1;
      int w2;
      int w3;  
      
      // Constructors //
      qR7(int w1_=0, int w2_=0, int w3_=0) : w1(w1_), w2(w2_), w3(w3_) {};	// Constructor 	(Default)
      qR7(int w) : w1(w), w2(w), w3(w) {};					//		(constant)

      // Other member functions //
      int set(int w1_, int w2_, int w3_);					// Set value
      bool isequal (const char * Wstr);						// (isequal to string)

      // Overloaded operators //
      bool operator == (const qR7 & wp) const;					// (isequal)
      bool operator != (const qR7 & wp) const;					// (~isequal)
      
      // Friend function to print out value of W
      friend std::ostream & operator << (std::ostream & o, const qR7 & W)
         { o << "(" << W.w1 << W.w2 << W.w3 << ")"; return o; }
};


// --------------------------------------------------------------------------------------------------------------- //
// Defines a class to contain and manipulate the two numbers that label the irreducible representation of G7, 
//    after Racah. These are the U quantum numbers that help label states of f-configurations.
// --------------------------------------------------------------------------------------------------------------- //
class qG2 
{
   private:

   public:
      int u1;
      int u2;
      
      // Constructors //
      qG2(int u1_=0, int u2_=0) : u1(u1_), u2(u2_) {};				// Constructor	(Default)
      qG2(int u) : u1(u), u2(u) {};						//		(constant)

      // Other member functions //
      int set(int u1_, int u2_);						// Set value
      bool isequal (const char * Ustr);						// (isequal to string)

      // Overloaded operators //
      bool operator == (const qG2 & up) const;					// (isequal)
      bool operator != (const qG2 & up) const;					// (isequal)

      // Friend function to print out value of U
      friend std::ostream & operator << (std::ostream & o, const qG2 & U)
         { o << "(" << U.u1 << U.u2 << ")"; return o; }
};

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in states.cpp, which are not constructors/members of any class below
// --------------------------------------------------------------------------------------------------------------- //
qR7 racah_vtow(int S2, int v);                                  // Converts a seniority number v to W quantum numbers.
int racah_wtov(int S2, qR7 W);                                  // Converts a (W1 W2 W3) quantum number to seniority num.
// --------------------------------------------------------------------------------------------------------------- //

// --------------------------------------------------------------------------------------------------------------- //
// Defines a class containing all the quantum numbers which characterise an f-electron state.
//   ref: Neilson and Koster, 1963. B.R. Judd, 1961 [see ic1ion.hpp]
// --------------------------------------------------------------------------------------------------------------- //
class fstates_t {
   private:

   public:

   // -----------------  Member variables --------------- //

      int S2;		// 2*S   twice spin quantum number
      orbital L;	//   L   orbital quantum number
      int v;		//   v   seniority quantum number
      qR7 W;		//   W   quantum number
      qG2 U;		//   U   quantum number
      std::string id;	//       identification string, after Neilson and Koster 1963.
      bool t;		// tau   quantum number to distinguish states of same values of all other q. num.
      int J2;		// 2*J   twice total angular momentum quantum number
      int mJ2;          // 2*mJ  twice z-component of total ang. mom. quantum number

   // -----------------  Member functions --------------- //
    //fstates_t set(const char *id_);				// Sets the values of the quantum numbers depeding on label
    //qR7 racah_vtow(int S2, int v);				// Converts the seniority (v) to Racah's label for R7 (W)
    //int racah_wtov(int S2, qR7 W);				// Converts W back to v

      // Constructor functions //
      fstates_t() : S2(1), L(S), v(1), id("nlvd") {};
      fstates_t(int S2_, orbital L_, int v_, qG2 U_, std::string id_) : S2(S2_), L(L_), v(v_), U(U_), id(id_), t(0)
      { 
         W = racah_vtow(S2_,v_); 
      }
      fstates_t(int S2_, orbital L_, int v_, qG2 U_, std::string id_, int J2_) 
         : S2(S2_), L(L_), v(v_), U(U_), id(id_), J2(J2_) { W = racah_vtow(S2_,v_); }
      fstates_t(int S2_, orbital L_, int v_, qG2 U_, std::string id_, int J2_, int mJ2_) 
         : S2(S2_), L(L_), v(v_), U(U_), id(id_), J2(J2_), mJ2(mJ2_) { W = racah_vtow(S2_,v_); }
      fstates_t(int S2_, orbital L_, int v_, qG2 U_, std::string id_, bool _t) 
         : S2(S2_), L(L_), v(v_), U(U_), id(id_), t(_t) { W = racah_vtow(S2_,v_); }
      fstates_t(int S2_, orbital L_, int v_, std::string id_) : S2(S2_), L(L_), v(v_), id(id_) {}; // For d-electrons
      fstates_t(int S2_, orbital L_, int v_, std::string id_, int J2_)                             // For d-electrons
         : S2(S2_), L(L_), v(v_), id(id_), J2(J2_) {};
      fstates_t(int S2_, orbital L_, int v_, std::string id_, int J2_, int mJ2_)                   // For d-electrons
         : S2(S2_), L(L_), v(v_), id(id_), J2(J2_), mJ2(mJ2_) {};
      fstates_t(int S2_, orbital L_, std::string id_) : S2(S2_), L(L_), id(id_) {};                // For p-electrons
      fstates_t(int S2_, orbital L_, std::string id_, int J2_) : S2(S2_),L(L_),id(id_),J2(J2_) {}; // For p-electrons
      fstates_t(int S2_, orbital L_, std::string id_, int J2_, int mJ2_)                           // For p-electrons
         : S2(S2_), L(L_), id(id_), J2(J2_), mJ2(mJ2_) {};


    //~fstates_t() {};                                          // Destructor
};

// --------------------------------------------------------------------------------------------------------------- //
// Defines a class to hold all the states of a particular f-configuration
// --------------------------------------------------------------------------------------------------------------- //
class fconf
{
   private:

   public:
      std::vector<fstates_t> states;                            // Defines a standard vector to contain all the states

      // Constructors //
      fconf(orbital l=F);                                       // Default looks up |vLS> (d-) or |vULS> (f-electrons)
      fconf(int n, orbital l=F);                                // n = number of equivalent electrons
      fconf(int n, bool mJflag, orbital l=F);                   // Construct matrix in |aLSJ> or |aLSmJ> basis

};

#endif
