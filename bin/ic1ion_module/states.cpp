/* states.cpp
 *
 * Defines the states of a particular f^n or d^n configurations on a group theoretical basis, after Racah.
 *
 * Class:       fstates_t                                       // A class to hold the quantum numbers of a f-electron LS state
 *              fstates                                         // A class to hold a vector of LS states for an f^n configuration
 * Functions:   std::string strtoupper(std::string instring);   // Converts the characters of a string to upper case
 *              std::string Lout(const orbital & l);            // Converts the enumerated value of l into a string
 *              orbital Lin(std::string & l);                   // Converts a string of the spectroscopy notation to l
 *              qR7 racah_vtow(int S2, int v);                  // Converts a seniority number v to W quantum numbers.
 *              int racah_wtov(int S2, qR7 W);                  // Converts a (W1 W2 W3) quantum number to seniority num.
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 *
 * Refs: Nielson and Koster, Spectroscopic Coefficients for the $p^n$, $d^n$, and $f^n$ Configurations, MIT Press, NY (1963)
 *
 */

#include "ic1ion.hpp"
#include <cctype>
#include <cstdio>

// --------------------------------------------------------------------------------------------------------------- //
// Converts a C++ string to upper case
// --------------------------------------------------------------------------------------------------------------- //
std::string strtoupper(std::string instring)
{
   std::string outstring;
   int i,strlength = instring.length();
   for(i=0; i<strlength; i++)
      outstring[i] = std::toupper(instring[i]);
   return outstring;
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to output the correct labels for the orbital enumeration
// --------------------------------------------------------------------------------------------------------------- //
std::string Lout(const orbital & l)
{
   std::stringstream o;
   switch(l)
   {
      case  S: o << "S";  break;
      case  P: o << "P";  break;
      case  D: o << "D";  break;
      case  F: o << "F";  break;
      case Fp: o << "Fp"; break;
      case  G: o << "G";  break;
      case Gp: o << "Gp"; break;
      case  H: o << "H";  break;
      case Hp: o << "Hp"; break;
      case  I: o << "I";  break;
      case Ip: o << "Ip"; break;
      case  K: o << "K";  break;
      case Kp: o << "Kp"; break;
      case  L: o << "L";  break;
      case Lp: o << "Lp"; break;
      case  M: o << "M";  break;
      case  N: o << "N";  break;
      case  O: o << "O";  break;
      case  Q: o << "Q";  break;
      default: o << ""; 
   }
   return o.str();
}

orbital Lin(std::string & l)
{
   orbital ls; int lint;
   std::istringstream iss(l); 

        if(l.compare("S")==0  || l.compare("s")==0)  ls = S;
   else if(l.compare("P")==0  || l.compare("p")==0)  ls = P;
   else if(l.compare("D")==0  || l.compare("d")==0)  ls = D;
   else if(l.compare("F")==0  || l.compare("f")==0)  ls = F;
   else if(l.compare("Fp")==0 || l.compare("fp")==0) ls = Fp;
   else if(l.compare("G")==0  || l.compare("g")==0)  ls = G;
   else if(l.compare("Gp")==0 || l.compare("gp")==0) ls = Gp;
   else if(l.compare("H")==0  || l.compare("h")==0)  ls = H;
   else if(l.compare("Hp")==0 || l.compare("hp")==0) ls = Hp;
   else if(l.compare("I")==0  || l.compare("i")==0)  ls = I;
   else if(l.compare("Ip")==0 || l.compare("ip")==0) ls = Ip;
   else if(l.compare("K")==0  || l.compare("k")==0)  ls = K;
   else if(l.compare("Kp")==0 || l.compare("kp")==0) ls = Kp;
   else if(l.compare("L")==0  || l.compare("l")==0)  ls = L;
   else if(l.compare("Lp")==0 || l.compare("lp")==0) ls = Lp;
   else if(l.compare("M")==0  || l.compare("m")==0)  ls = M;
   else if(l.compare("N")==0  || l.compare("n")==0)  ls = N;
   else if(l.compare("O")==0  || l.compare("o")==0)  ls = O;
   else if(l.compare("Q")==0  || l.compare("q")==0)  ls = Q;
   else                                  // If string lin is not a letter as expected, assume it is an integer
   {
     iss >> lint; ls = (orbital)lint;
   }
   return ls;
}

// --------------------------------------------------------------------------------------------------------------- //
// Member functions and overloaded operators for the qR7 class, for the labels W of the group R7
// --------------------------------------------------------------------------------------------------------------- //
int qR7::set(int w1_, int w2_, int w3_)
{
   w1 = w1_; w2 = w2_; w3 = w3_; return 0;
}

bool qR7::operator == (const qR7 & wp) const
{
   if (w1==wp.w1 && w2==wp.w2 && w3==wp.w3) { return true; }
   else { return false; }
}

bool qR7::operator != (const qR7 & wp) const
{
   if (w1!=wp.w1 || w2!=wp.w2 || w3!=wp.w3) { return true; }
   else { return false; }
}

bool qR7::isequal (const char * Wstr)
{
   char ws=0;
   int  w1n=0,w2n=1,w3n=2;

   ws = Wstr[0]; w1n = atoi(&ws);
   ws = Wstr[1]; w2n = atoi(&ws);
   ws = Wstr[2]; w3n = atoi(&ws);

   if (w1==w1n && w2==w2n && w3==w3n) { return true; }
   else { return false; }
}

// --------------------------------------------------------------------------------------------------------------- //
// Member functions and overloaded operators for the qG2 class, for the labels U of the group G2
// --------------------------------------------------------------------------------------------------------------- //
int qG2::set(int u1_, int u2_)
{
   u1 = u1_; u2 = u2_; return 0;
}

bool qG2::operator == (const qG2 & up) const
{
   if (u1==up.u1 && u2==up.u2) { return true; }
   else { return false; }
}

bool qG2::operator != (const qG2 & up) const
{
   if (u1!=up.u1 || u2!=up.u2) { return true; }
   else { return false; }
}

bool qG2::isequal (const char * Ustr)
{
   char us=0;
   int  u1n=0,u2n=1;

   us = Ustr[0]; u1n = atoi(&us);
   us = Ustr[1]; u2n = atoi(&us);

   if (u1==u1n && u2==u2n) { return true; }
   else { return false; }
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to convert from the seniority quantum number, v, to Racah's [W1 W2 W3] quantum number which labels the
//    irreducible represenation of group R7. Reference: Racah IV, eqn. 16
// --------------------------------------------------------------------------------------------------------------- //
qR7 racah_vtow(int S2, int v)
{
   int a2 = v-S2;				// a = v/2 - S
   int b  = 2*3+1-v;
   if (S2<b) b = S2;				// b = min(2S, 2l+1-v) - l=3 for f-electrons
   int b2 = 2*b;
   qR7 W;					// Initialise W to zeros

   if (a2>=0) 
   {
      if (a2>1) W.w1 = 2;
      if (a2>3) W.w2 = 2;			// w_1,...,w_a = 2
      if (a2>5) W.w3 = 2;
      if (b>0)
      {
         if (a2<1 && (a2+b2)>1) W.w1 = 1;
         if (a2<3 && (a2+b2)>3) W.w2 = 1;	// w_{a+1},...,w_{a+b} = 1
         if (a2<5 && (a2+b2)>5) W.w3 = 1;
      }
   }
   return W;
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to convert from the Racah's [W1 W2 W3] quantum number back to the seniority number v
// --------------------------------------------------------------------------------------------------------------- //
int racah_wtov(int S2, qR7 W)
{
   int a=0;	// w_1,...,w_a = 2		// a = v/2 - S
   int b=0;	// w_{a+1},...,w_{a+b} = 1	// b = min(2S, 2l+1-v)
   int v;

   switch(W.w1)
   {
      case 2: a = 1; break;
      case 1: b = 1; break;
   }
   switch(W.w2)
   {
      case 2: a = 2; break;
      case 1: b = 2; break;
   }
   switch(W.w3)
   {
      case 2: a = 3; break;
      case 1: b = 3; break;
   }

   if (a>0) { v = 2*a + S2; }			// i.e. v>2S so b=2S, so v = 2(a+S)
   else { v = 2*3+1 - b; }			// i.e. v<2S so a=0, b=2l+1-v so v = 2l+1-b
  
   return v;
}

// --------------------------------------------------------------------------------------------------------------- //
// Member functions and overloaded operators for the fconf class, which represents a particular l^n configuration
// --------------------------------------------------------------------------------------------------------------- //
fconf::fconf(orbital l)
{
   if(l==P) 
      states.push_back(fstates_t(1,P,"2P"));
   else if(l==D) 
      states.push_back(fstates_t(1,D,1,"2D"));
   else if(l==F) 
   {
      qG2 U(1,0);
      states.push_back(fstates_t(1,F,1,U,"2F")); 
   }
   else
      std::cerr << "fconf::fconf() - error, only the case of l=1, l=2 and l=3, p-, d- and f-electrons implemented.\n";
}

fconf::fconf(int n, orbital l)
{
   qG2 U; 
   
   switch(l) {

   case S:
   if (n!=1)		// Checks to see if number of s-electrons valid
   {
      std::cerr << "fconf::fconf() - Invalid value of n = number of s-electrons (must be 1)\n";
   }
   else { states.push_back(fstates_t(1,S,"2S")); }
   break; // case S:

   case P:
   if (n<1 || n>5)	// Checks to see if number of p-electrons valid
   {
      std::cerr << "fconf::fconf() - Invalid value of n = number of p-electrons (must be 1<= n <=5)\n";
   }
   else {
      if (n>3) n=6-n;	// Checks to see if we are in the second half of the series
      switch(n)
      {
         case 1:
            states.push_back(fstates_t(1,P,"2P"));                        // p1    2P
            break;
         case 2:
	    states.reserve(3);
            states.push_back(fstates_t(2,P,"3P"));                        // p2    3P
            states.push_back(fstates_t(0,S,"1S"));                        // p2    1S
            states.push_back(fstates_t(0,D,"1D"));                        // p2    1D
            break;
         case 3:
	    states.reserve(3);
            states.push_back(fstates_t(3,S,"4S"));                        // p3    4S
            states.push_back(fstates_t(1,P,"2P"));                        // p3    2P
            states.push_back(fstates_t(1,D,"2D"));                        // p3    2D
            break;
      } // switch(n)
   }    // else (n>0 && n<11)

   break; // case P:
   
   case D:

   if (n<1 || n>9)	// Checks to see if number of d-electrons valid
   {
      std::cerr << "fconf::fconf() - Invalid value of n = number of d-electrons (must be 1<= n <=9)\n";
   }
   else {
      if (n>5) n=10-n;	// Checks to see if we are in the second half of the series
      switch(n)
      {
         case 1: 
            states.push_back(fstates_t(1,D,1,"2D"));                      // d1    2D   1 (10)
            break;
         case 2:         
	    states.reserve(5);
            states.push_back(fstates_t(2,P,2,"3P"));                      // d2    3P   2 (11)
            states.push_back(fstates_t(2,F,2,"3F"));                      //       3F   2 (11)
            states.push_back(fstates_t(0,S,0,"1S"));                      //       1S   0 (00)
            states.push_back(fstates_t(0,D,2,"1D"));                      //       1D   2 (20)
            states.push_back(fstates_t(0,G,2,"1G"));                      //       1G   2 (20)
            break;
         case 3:
	    states.reserve(8);
            states.push_back(fstates_t(3,P,3,"4P"));                      // d3    4P   3 (11)
            states.push_back(fstates_t(3,F,3,"4F"));                      //       4F   3 (11)
            states.push_back(fstates_t(1,P,3,"2P"));                      //       2P   3 (21)
            states.push_back(fstates_t(1,D,1,"2D1"));                     //       2D1  1 (10)
            states.push_back(fstates_t(1,D,3,"2D2"));                     //       2D2  3 (21)
            states.push_back(fstates_t(1,F,3,"2F"));                      //       2F   3 (21)
            states.push_back(fstates_t(1,G,3,"2G"));                      //       2G   3 (21)
            states.push_back(fstates_t(1,H,3,"2H"));                      //       2H   3 (21)
            break;
         case 4:
	    states.reserve(16);
            states.push_back(fstates_t(4,D,4,"5D"));                      // d4    5D   4 (10)
            states.push_back(fstates_t(2,P,2,"3P1"));                     //       3P1  2 (11)
            states.push_back(fstates_t(2,P,4,"3P2"));                     //       3P2  4 (21)
            states.push_back(fstates_t(2,D,4,"3D"));                      //       3D   4 (21)
            states.push_back(fstates_t(2,F,2,"3F1"));                     //       3F1  2 (11)
            states.push_back(fstates_t(2,F,4,"3F2"));                     //       3F2  4 (21)
            states.push_back(fstates_t(2,G,4,"3G"));                      //       3G   4 (21)
            states.push_back(fstates_t(2,H,4,"3H"));                      //       3H   4 (21)
            states.push_back(fstates_t(0,S,0,"1S1"));                     //       1S1  0 (00)
            states.push_back(fstates_t(0,S,4,"1S2"));                     //       1S2  4 (22)
            states.push_back(fstates_t(0,D,2,"1D1"));                     //       1D1  2 (20)
            states.push_back(fstates_t(0,D,4,"1D2"));                     //       1D2  4 (22)
            states.push_back(fstates_t(0,F,4,"1F"));                      //       1F   4 (22)
            states.push_back(fstates_t(0,G,2,"1G1"));                     //       1G1  2 (20)
            states.push_back(fstates_t(0,G,4,"1G2"));                     //       1G2  4 (22)
            states.push_back(fstates_t(0,I,4,"1I"));                      //       1I   4 (22)
            break;
         case 5:
	    states.reserve(16);
            states.push_back(fstates_t(5,S,5,"6S"));                      // d5    6S   5 (00)
            states.push_back(fstates_t(3,P,3,"4P"));                      //       4P   3 (11)
            states.push_back(fstates_t(3,D,5,"4D"));                      //       4D   5 (20)
            states.push_back(fstates_t(3,F,3,"4F"));                      //       4F   3 (11)
            states.push_back(fstates_t(3,G,5,"4G"));                      //       4G   5 (20)
            states.push_back(fstates_t(1,S,5,"2S"));                      //       2S   5 (22)
            states.push_back(fstates_t(1,P,3,"2P"));                      //       2P   3 (21)
            states.push_back(fstates_t(1,D,1,"2D1"));                     //       2D1  1 (10)
            states.push_back(fstates_t(1,D,3,"2D2"));                     //       2D2  3 (21)
            states.push_back(fstates_t(1,D,5,"2D3"));                     //       2D3  5 (22)
            states.push_back(fstates_t(1,F,3,"2F1"));                     //       2F1  3 (21)
            states.push_back(fstates_t(1,F,5,"2F2"));                     //       2F2  5 (22)
            states.push_back(fstates_t(1,G,3,"2G1"));                     //       2G1  3 (21)
            states.push_back(fstates_t(1,G,5,"2G2"));                     //       2G2  5 (22)
            states.push_back(fstates_t(1,H,3,"2H"));                      //       2H   3 (21)
            states.push_back(fstates_t(1,I,5,"2I"));                      //       2I   5 (22)
            break;
      } // switch(n)
   }    // else (n>0 && n<11)

   break; // case D:
   
   case F:

   if (n<1 || n>13)	// Checks to see if number of f-electrons valid
   {
      std::cerr << "fconf::fconf(int n) - Invalid value of n = number of f-electrons (must be 1<= n <=13)\n";
   }
   else {

      if (n>7) n=14-n;	// Checks to see if we are in the second half of the series

      switch(n)
      {
         case 1:
            U.set(1,0); states.push_back(fstates_t(1,F,1,U,"2F"));        // f1    2F	1 100 10
            break;
         case 2:
	    states.reserve(7);
            U.set(1,1); states.push_back(fstates_t(2,P,2,U,"3P"));        // f2    3P	2 110 11
            U.set(1,0); states.push_back(fstates_t(2,F,2,U,"3F"));        //       3F	2 110 10
            U.set(1,1); states.push_back(fstates_t(2,H,2,U,"3H"));        //       3H	2 110 11
            U.set(0,0); states.push_back(fstates_t(0,S,0,U,"1S"));        //       1S	0 000 00
            U.set(2,0); states.push_back(fstates_t(0,D,2,U,"1D"));        //       1D	2 200 20
            U.set(2,0); states.push_back(fstates_t(0,G,2,U,"1G"));        //       1G	2 200 20
            U.set(2,0); states.push_back(fstates_t(0,I,2,U,"1I"));        //       1I	2 200 20
            break;
         case 3:
	    states.reserve(17);
            U.set(0,0); states.push_back(fstates_t(3,S,3,U,"4S"));        // f3    4S	3 111 00
            U.set(2,0); states.push_back(fstates_t(3,D,3,U,"4D"));        //       4D	3 111 20
            U.set(1,0); states.push_back(fstates_t(3,F,3,U,"4F"));        //       4F	3 111 10
            U.set(2,0); states.push_back(fstates_t(3,G,3,U,"4G"));        //       4G	3 111 20
            U.set(2,0); states.push_back(fstates_t(3,I,3,U,"4I"));        //       4I	3 111 20
            U.set(1,1); states.push_back(fstates_t(1,P,3,U,"2P"));        //       2P	3 210 11
            U.set(2,0); states.push_back(fstates_t(1,D,3,U,"2D1"));       //       2D1	3 210 20
            U.set(2,1); states.push_back(fstates_t(1,D,3,U,"2D2"));       //       2D2	3 210 21
            U.set(1,0); states.push_back(fstates_t(1,F,1,U,"2F1"));       //       2F1	1 100 10
            U.set(2,1); states.push_back(fstates_t(1,F,3,U,"2F2"));       //       2F2	3 210 21
            U.set(2,0); states.push_back(fstates_t(1,G,3,U,"2G1"));       //       2G1	3 210 20
            U.set(2,1); states.push_back(fstates_t(1,G,3,U,"2G2"));       //       2G2	3 210 21
            U.set(1,1); states.push_back(fstates_t(1,H,3,U,"2H1"));       //       2H1	3 210 11
            U.set(2,1); states.push_back(fstates_t(1,H,3,U,"2H2"));       //       2H2	3 210 21
            U.set(2,0); states.push_back(fstates_t(1,I,3,U,"2I"));        //       2I	3 210 20
            U.set(2,1); states.push_back(fstates_t(1,K,3,U,"2K"));        //       2K	3 210 21
            U.set(2,1); states.push_back(fstates_t(1,L,3,U,"2L"));        //       2L	3 210 21
            break;
         case 4:
	    states.reserve(47);
            U.set(0,0); states.push_back(fstates_t(4,S,4,U,"5S"));        // f4    5S	4 111 00
            U.set(2,0); states.push_back(fstates_t(4,D,4,U,"5D"));        //       5D	4 111 20
            U.set(1,0); states.push_back(fstates_t(4,F,4,U,"5F"));        //       5F	4 111 10
            U.set(2,0); states.push_back(fstates_t(4,G,4,U,"5G"));        //       5G	4 111 20
            U.set(2,0); states.push_back(fstates_t(4,I,4,U,"5I"));        //       5I	4 111 20
            U.set(1,1); states.push_back(fstates_t(2,P,2,U,"3P1"));       //       3P1	2 110 11
            U.set(1,1); states.push_back(fstates_t(2,P,4,U,"3P2"));       //       3P2	4 211 11
            U.set(3,0); states.push_back(fstates_t(2,P,4,U,"3P3"));       //       3P3	4 211 30
            U.set(2,0); states.push_back(fstates_t(2,D,4,U,"3D1"));       //       3D1	4 211 20
            U.set(2,1); states.push_back(fstates_t(2,D,4,U,"3D2"));       //       3D2	4 211 21
            U.set(1,0); states.push_back(fstates_t(2,F,2,U,"3F1"));       //       3F1	2 110 10
            U.set(1,0); states.push_back(fstates_t(2,F,4,U,"3F2"));       //       3F2	4 211 10
            U.set(2,1); states.push_back(fstates_t(2,F,4,U,"3F3"));       //       3F3	4 211 21
            U.set(3,0); states.push_back(fstates_t(2,F,4,U,"3F4"));       //       3F4	4 211 30
            U.set(2,0); states.push_back(fstates_t(2,G,4,U,"3G1"));       //       3G1	4 211 20
            U.set(2,1); states.push_back(fstates_t(2,G,4,U,"3G2"));       //       3G2	4 211 21
            U.set(3,0); states.push_back(fstates_t(2,G,4,U,"3G3"));       //       3G3	4 211 30
            U.set(1,1); states.push_back(fstates_t(2,H,2,U,"3H1"));       //       3H1	2 110 11
            U.set(1,1); states.push_back(fstates_t(2,H,4,U,"3H2"));       //       3H2	4 211 11
            U.set(2,1); states.push_back(fstates_t(2,H,4,U,"3H3"));       //       3H3	4 211 21
            U.set(3,0); states.push_back(fstates_t(2,H,4,U,"3H4"));       //       3H4	4 211 30
            U.set(2,0); states.push_back(fstates_t(2,I,4,U,"3I1"));       //       3I1	4 211 20
            U.set(3,0); states.push_back(fstates_t(2,I,4,U,"3I2"));       //       3I2	4 211 30
            U.set(2,1); states.push_back(fstates_t(2,K,4,U,"3K1"));       //       3K1	4 211 21
            U.set(3,0); states.push_back(fstates_t(2,K,4,U,"3K2"));       //       3K2	4 211 30
            U.set(2,1); states.push_back(fstates_t(2,L,4,U,"3L"));        //       3L	4 211 21
            U.set(3,0); states.push_back(fstates_t(2,M,4,U,"3M"));        //       3M	4 211 30
            U.set(0,0); states.push_back(fstates_t(0,S,0,U,"1S1"));       //       1S1	0 000 00
            U.set(2,2); states.push_back(fstates_t(0,S,4,U,"1S2"));       //       1S2	4 220 22
            U.set(2,0); states.push_back(fstates_t(0,D,2,U,"1D1"));       //       1D1	2 200 20
            U.set(2,0); states.push_back(fstates_t(0,D,4,U,"1D2"));       //       1D2	4 220 20
            U.set(2,1); states.push_back(fstates_t(0,D,4,U,"1D3"));       //       1D3	4 220 21
            U.set(2,2); states.push_back(fstates_t(0,D,4,U,"1D4"));       //       1D4	4 220 22
            U.set(2,1); states.push_back(fstates_t(0,F,4,U,"1F"));        //       1F	4 220 21
            U.set(2,0); states.push_back(fstates_t(0,G,2,U,"1G1"));       //       1G1	2 200 20
            U.set(2,0); states.push_back(fstates_t(0,G,4,U,"1G2"));       //       1G2	4 220 20
            U.set(2,1); states.push_back(fstates_t(0,G,4,U,"1G3"));       //       1G3	4 220 21
            U.set(2,2); states.push_back(fstates_t(0,G,4,U,"1G4"));       //       1G4	4 220 22
            U.set(2,1); states.push_back(fstates_t(0,H,4,U,"1H1"));       //       1H1	4 220 21
            U.set(2,2); states.push_back(fstates_t(0,H,4,U,"1H2"));       //       1H2	4 220 22
            U.set(2,0); states.push_back(fstates_t(0,I,2,U,"1I1"));       //       1I1	2 200 20
            U.set(2,0); states.push_back(fstates_t(0,I,4,U,"1I2"));       //       1I2	4 220 20
            U.set(2,2); states.push_back(fstates_t(0,I,4,U,"1I3"));       //       1I3	4 220 22
            U.set(2,1); states.push_back(fstates_t(0,K,4,U,"1K"));        //       1K	4 220 21
            U.set(2,1); states.push_back(fstates_t(0,L,4,U,"1L1"));       //       1L1	4 220 21
            U.set(2,2); states.push_back(fstates_t(0,L,4,U,"1L2"));       //       1L2	4 220 22
            U.set(2,2); states.push_back(fstates_t(0,N,4,U,"1N"));        //       1N	4 220 22
   	    break;
         case 5:
	    states.reserve(73);
            U.set(1,1); states.push_back(fstates_t(5,P,5,U,"6P"));        // f5    6P	5 110 11
            U.set(1,0); states.push_back(fstates_t(5,F,5,U,"6F"));        //       6F	5 110 10
            U.set(1,1); states.push_back(fstates_t(5,H,5,U,"6H"));        //       6H	5 110 11
            U.set(0,0); states.push_back(fstates_t(3,S,3,U,"4S"));        //       4S	3 111 00
            U.set(1,1); states.push_back(fstates_t(3,P,5,U,"4P1"));       //       4P1	5 211 11
            U.set(3,0); states.push_back(fstates_t(3,P,5,U,"4P2"));       //       4P2	5 211 30
            U.set(2,0); states.push_back(fstates_t(3,D,3,U,"4D1"));       //       4D1	3 111 20
            U.set(2,0); states.push_back(fstates_t(3,D,5,U,"4D2"));       //       4D2	5 211 20
            U.set(2,1); states.push_back(fstates_t(3,D,5,U,"4D3"));       //       4D3	5 211 21
            U.set(1,0); states.push_back(fstates_t(3,F,3,U,"4F1"));       //       4F1	3 111 10
            U.set(1,0); states.push_back(fstates_t(3,F,5,U,"4F2"));       //       4F2	5 211 10
            U.set(2,1); states.push_back(fstates_t(3,F,5,U,"4F3"));       //       4F3	5 211 21
            U.set(3,0); states.push_back(fstates_t(3,F,5,U,"4F4"));       //       4F4	5 211 30
            U.set(2,0); states.push_back(fstates_t(3,G,3,U,"4G1"));       //       4G1	3 111 20
            U.set(2,0); states.push_back(fstates_t(3,G,5,U,"4G2"));       //       4G2	5 211 20
            U.set(2,1); states.push_back(fstates_t(3,G,5,U,"4G3"));       //       4G3	5 211 21
            U.set(3,0); states.push_back(fstates_t(3,G,5,U,"4G4"));       //       4G4	5 211 30
            U.set(1,1); states.push_back(fstates_t(3,H,5,U,"4H1"));       //       4H1	5 211 11
            U.set(2,1); states.push_back(fstates_t(3,H,5,U,"4H2"));       //       4H2	5 211 21
            U.set(3,0); states.push_back(fstates_t(3,H,5,U,"4H3"));       //       4H3	5 211 30
            U.set(2,0); states.push_back(fstates_t(3,I,3,U,"4I1"));       //       4I1	3 111 20
            U.set(2,0); states.push_back(fstates_t(3,I,5,U,"4I2"));       //       4I2	5 211 20
            U.set(3,0); states.push_back(fstates_t(3,I,5,U,"4I3"));       //       4I3	5 211 30
            U.set(2,1); states.push_back(fstates_t(3,K,5,U,"4K1"));       //       4K1	5 211 21
            U.set(3,0); states.push_back(fstates_t(3,K,5,U,"4K2"));       //       4K2	5 211 30
            U.set(2,1); states.push_back(fstates_t(3,L,5,U,"4L"));        //       4L	5 211 21
            U.set(3,0); states.push_back(fstates_t(3,M,5,U,"4M"));        //       4M	5 211 30
            U.set(1,1); states.push_back(fstates_t(1,P,3,U,"2P1"));       //       2P1	3 210 11
            U.set(1,1); states.push_back(fstates_t(1,P,5,U,"2P2"));       //       2P2	5 211 11
            U.set(3,0); states.push_back(fstates_t(1,P,5,U,"2P3"));       //       2P3	5 211 30
            U.set(3,1); states.push_back(fstates_t(1,P,5,U,"2P4"));       //       2P4	5 211 31
            U.set(2,0); states.push_back(fstates_t(1,D,3,U,"2D1"));       //       2D1	3 210 20
            U.set(2,1); states.push_back(fstates_t(1,D,3,U,"2D2"));       //       2D2	3 210 21
            U.set(2,0); states.push_back(fstates_t(1,D,5,U,"2D3"));       //       2D3	5 211 20
            U.set(2,1); states.push_back(fstates_t(1,D,5,U,"2D4"));       //       2D4 	5 211 21
            U.set(3,1); states.push_back(fstates_t(1,D,5,U,"2D5"));       //       2D5	5 211 31
            U.set(1,0); states.push_back(fstates_t(1,F,1,U,"2F1"));       //       2F1	1 100 10
            U.set(2,1); states.push_back(fstates_t(1,F,3,U,"2F2"));       //       2F2	3 210 21
            U.set(1,0); states.push_back(fstates_t(1,F,5,U,"2F3"));       //       2F3	5 221 10
            U.set(2,1); states.push_back(fstates_t(1,F,5,U,"2F4"));       //       2F4	5 221 21
            U.set(3,0); states.push_back(fstates_t(1,F,5,U,"2F5"));       //       2F5	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,F,5,U,"2F6"));       //       2F6	5 221 31A
            U.set(3,1); states.push_back(fstates_t(1,Fp,5,U,"2F7",1));    //       2F7	5 221 31B
            U.set(2,0); states.push_back(fstates_t(1,G,3,U,"2G1"));       //       2G1	3 210 20
            U.set(2,1); states.push_back(fstates_t(1,G,3,U,"2G2"));       //       2G2	3 210 21
            U.set(2,0); states.push_back(fstates_t(1,G,5,U,"2G3"));       //       2G3	5 221 20
            U.set(2,1); states.push_back(fstates_t(1,G,5,U,"2G4"));       //       2G4 	5 221 21
            U.set(3,0); states.push_back(fstates_t(1,G,5,U,"2G5"));       //       2G5	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,G,5,U,"2G6"));       //       2G6	5 221 31
            U.set(1,1); states.push_back(fstates_t(1,H,3,U,"2H1"));       //       2H1	3 210 11
            U.set(2,1); states.push_back(fstates_t(1,H,3,U,"2H2"));       //       2H2	3 210 21
            U.set(1,1); states.push_back(fstates_t(1,H,5,U,"2H3"));       //       2H3	5 221 11
            U.set(2,1); states.push_back(fstates_t(1,H,5,U,"2H4"));       //       2H4 	5 221 21
            U.set(3,0); states.push_back(fstates_t(1,H,5,U,"2H5"));       //       2H5	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,H,5,U,"2H6"));       //       2H6	5 221 31A
            U.set(3,1); states.push_back(fstates_t(1,Hp,5,U,"2H7",1));    //       2H7	5 221 31B
            U.set(2,0); states.push_back(fstates_t(1,I,3,U,"2I1"));       //       2I1	3 210 20
            U.set(2,0); states.push_back(fstates_t(1,I,5,U,"2I2"));       //       2I2	5 221 20
            U.set(3,0); states.push_back(fstates_t(1,I,5,U,"2I3"));       //       2I3	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,I,5,U,"2I4"));       //       2I4	5 221 31A
            U.set(3,1); states.push_back(fstates_t(1,Ip,5,U,"2I5",1));    //       2I5	5 221 31B
            U.set(2,1); states.push_back(fstates_t(1,K,3,U,"2K1"));       //       2K1	3 210 21
            U.set(2,1); states.push_back(fstates_t(1,K,5,U,"2K2"));       //       2K2	5 221 21
            U.set(3,0); states.push_back(fstates_t(1,K,5,U,"2K3"));       //       2K3	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,K,5,U,"2K4"));       //       2K4	5 221 31A
            U.set(3,1); states.push_back(fstates_t(1,Kp,5,U,"2K5",1));    //       2K5	5 221 31B
            U.set(2,1); states.push_back(fstates_t(1,L,3,U,"2L1"));       //       2L1	3 210 21
            U.set(2,1); states.push_back(fstates_t(1,L,5,U,"2L2"));       //       2L2	5 221 21
            U.set(3,1); states.push_back(fstates_t(1,L,5,U,"2L3"));       //       2L3	5 221 31
            U.set(3,0); states.push_back(fstates_t(1,M,5,U,"2M1"));       //       2M1	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,M,5,U,"2M2"));       //       2M2	5 221 31
            U.set(3,1); states.push_back(fstates_t(1,N,5,U,"2N"));        //       2N	5 221 31
            U.set(3,1); states.push_back(fstates_t(1,O,5,U,"2O"));        //       2O	5 221 31
   	    break;
         case 6:
	    states.reserve(119);
            U.set(1,0); states.push_back(fstates_t(6,F,6,U,"7F"));        // f6    7F	6 100 10
            U.set(0,0); states.push_back(fstates_t(4,S,4,U,"5S"));        //       5S	4 111 00
            U.set(1,1); states.push_back(fstates_t(4,P,6,U,"5P"));        //       5P	6 210 11
            U.set(2,0); states.push_back(fstates_t(4,D,4,U,"5D1"));       //       5D1	4 111 20
            U.set(2,0); states.push_back(fstates_t(4,D,6,U,"5D2"));       //       5D2	6 210 20
            U.set(2,1); states.push_back(fstates_t(4,D,6,U,"5D3"));       //       5D3	6 210 21
            U.set(1,0); states.push_back(fstates_t(4,F,4,U,"5F1"));       //       5F1	4 111 10
            U.set(2,1); states.push_back(fstates_t(4,F,6,U,"5F2"));       //       5F2	6 210 21
            U.set(2,0); states.push_back(fstates_t(4,G,4,U,"5G1"));       //       5G1	4 111 20
            U.set(2,0); states.push_back(fstates_t(4,G,6,U,"5G2"));       //       5G2	6 210 20
            U.set(2,1); states.push_back(fstates_t(4,G,6,U,"5G3"));       //       5G3	6 210 21
            U.set(1,1); states.push_back(fstates_t(4,H,6,U,"5H1"));       //       5H1	6 210 11
            U.set(2,1); states.push_back(fstates_t(4,H,6,U,"5H2"));       //       5H2	6 210 21
            U.set(2,0); states.push_back(fstates_t(4,I,4,U,"5I1"));       //       5I1	4 111 20
            U.set(2,0); states.push_back(fstates_t(4,I,6,U,"5I2"));       //       5I2	6 210 20
            U.set(2,1); states.push_back(fstates_t(4,K,6,U,"5K"));        //       5K	6 210 21
            U.set(2,1); states.push_back(fstates_t(4,L,6,U,"5L"));        //       5L	6 210 21
            U.set(1,1); states.push_back(fstates_t(2,P,2,U,"3P1"));       //       3P1	2 110 11
            U.set(1,1); states.push_back(fstates_t(2,P,4,U,"3P2"));       //       3P2	4 211 11
            U.set(3,0); states.push_back(fstates_t(2,P,4,U,"3P3"));       //       3P3	4 211 30
            U.set(1,1); states.push_back(fstates_t(2,P,6,U,"3P4"));       //       3P4	6 221 11
            U.set(3,0); states.push_back(fstates_t(2,P,6,U,"3P5"));       //       3P5	6 221 30
            U.set(3,1); states.push_back(fstates_t(2,P,6,U,"3P6"));       //       3P6	6 221 31
            U.set(2,0); states.push_back(fstates_t(2,D,4,U,"3D1"));       //       3D1	4 211 20
            U.set(2,1); states.push_back(fstates_t(2,D,4,U,"3D2"));       //       3D2	4 211 21
            U.set(2,0); states.push_back(fstates_t(2,D,6,U,"3D3"));       //       3D3	6 221 20
            U.set(2,1); states.push_back(fstates_t(2,D,6,U,"3D4"));       //       3D4	6 221 21
            U.set(3,1); states.push_back(fstates_t(2,D,6,U,"3D5"));       //       3D5	6 221 31
            U.set(1,0); states.push_back(fstates_t(2,F,2,U,"3F1"));       //       3F1	2 110 10
            U.set(1,0); states.push_back(fstates_t(2,F,4,U,"3F2"));       //       3F1	4 211 10
            U.set(2,1); states.push_back(fstates_t(2,F,4,U,"3F3"));       //       3F3	4 211 21
            U.set(3,0); states.push_back(fstates_t(2,F,4,U,"3F4"));       //       3F4	4 211 30
            U.set(1,0); states.push_back(fstates_t(2,F,6,U,"3F5"));       //       3F5	6 221 10
            U.set(2,1); states.push_back(fstates_t(2,F,6,U,"3F6"));       //       3F6	6 221 21
            U.set(3,0); states.push_back(fstates_t(2,F,6,U,"3F7"));       //       3F7	6 221 30
            U.set(3,1); states.push_back(fstates_t(2,F,6,U,"3F8"));       //       3F8	6 221 31A
            U.set(3,1); states.push_back(fstates_t(2,Fp,6,U,"3F9",1));    //       3F9	6 221 31B
            U.set(2,0); states.push_back(fstates_t(2,G,4,U,"3G1"));       //       3G1	4 211 20
            U.set(2,1); states.push_back(fstates_t(2,G,4,U,"3G2"));       //       3G2	4 211 21
            U.set(3,0); states.push_back(fstates_t(2,G,4,U,"3G3"));       //       3G3	4 211 30
            U.set(2,0); states.push_back(fstates_t(2,G,6,U,"3G4"));       //       3G4	6 221 20
            U.set(2,1); states.push_back(fstates_t(2,G,6,U,"3G5"));       //       3G5	6 221 21
            U.set(3,0); states.push_back(fstates_t(2,G,6,U,"3G6"));       //       3G6	6 221 30
            U.set(3,1); states.push_back(fstates_t(2,G,6,U,"3G7"));       //       3G7	6 221 31
            U.set(1,1); states.push_back(fstates_t(2,H,2,U,"3H1"));       //       3H1	2 110 11
            U.set(1,1); states.push_back(fstates_t(2,H,4,U,"3H2"));       //       3H2	4 211 11
            U.set(2,1); states.push_back(fstates_t(2,H,4,U,"3H3"));       //       3H3	4 211 21
            U.set(3,0); states.push_back(fstates_t(2,H,4,U,"3H4"));       //       3H4	4 211 30
            U.set(1,1); states.push_back(fstates_t(2,H,6,U,"3H5"));       //       3H5	6 221 11
            U.set(2,1); states.push_back(fstates_t(2,H,6,U,"3H6"));       //       3H6	6 221 21
            U.set(3,0); states.push_back(fstates_t(2,H,6,U,"3H7"));       //       3H7	6 221 30
            U.set(3,1); states.push_back(fstates_t(2,H,6,U,"3H8"));       //       3H8	6 221 31A
            U.set(3,1); states.push_back(fstates_t(2,Hp,6,U,"3H9",1));    //       3H9	6 221 31B
            U.set(2,0); states.push_back(fstates_t(2,I,4,U,"3I1"));       //       3I1	4 211 20
            U.set(3,0); states.push_back(fstates_t(2,I,4,U,"3I2"));       //       3I2	4 211 30
            U.set(2,0); states.push_back(fstates_t(2,I,6,U,"3I3"));       //       3I3	6 221 20
            U.set(3,0); states.push_back(fstates_t(2,I,6,U,"3I4"));       //       3I4	6 221 30
            U.set(3,1); states.push_back(fstates_t(2,I,6,U,"3I5"));       //       3I5	6 221 31A
            U.set(3,1); states.push_back(fstates_t(2,Ip,6,U,"3I6",1));    //       3I6	6 221 31B
            U.set(2,1); states.push_back(fstates_t(2,K,4,U,"3K1"));       //       3K1	4 211 21
            U.set(3,0); states.push_back(fstates_t(2,K,4,U,"3K2"));       //       3K2	4 211 30
            U.set(2,1); states.push_back(fstates_t(2,K,6,U,"3K3"));       //       3K3	6 221 21
            U.set(3,0); states.push_back(fstates_t(2,K,6,U,"3K4"));       //       3K4	6 221 30
            U.set(3,1); states.push_back(fstates_t(2,K,6,U,"3K5"));       //       3K5	6 221 31A
            U.set(3,1); states.push_back(fstates_t(2,Kp,6,U,"3K6",1));    //       3K6	6 221 31B
            U.set(2,1); states.push_back(fstates_t(2,L,4,U,"3L1"));       //       3L1	4 211 21
            U.set(2,1); states.push_back(fstates_t(2,L,6,U,"3L2"));       //       3L2	6 221 21
            U.set(3,1); states.push_back(fstates_t(2,L,6,U,"3L3"));       //       3L3	6 221 31
            U.set(3,0); states.push_back(fstates_t(2,M,4,U,"3M1"));       //       3M1	4 211 30
            U.set(3,0); states.push_back(fstates_t(2,M,6,U,"3M2"));       //       3M2	6 221 30
            U.set(3,1); states.push_back(fstates_t(2,M,6,U,"3M3"));       //       3M3	6 221 31
            U.set(3,1); states.push_back(fstates_t(2,N,6,U,"3N"));        //       3N	6 221 31
            U.set(3,1); states.push_back(fstates_t(2,O,6,U,"3O"));        //       3O	6 221 31
            U.set(0,0); states.push_back(fstates_t(0,S,0,U,"1S1"));       //       1S1	0 000 00
            U.set(2,2); states.push_back(fstates_t(0,S,4,U,"1S2"));       //       1S2	4 220 22
            U.set(0,0); states.push_back(fstates_t(0,S,6,U,"1S3"));       //       1S3	6 222 00
            U.set(4,0); states.push_back(fstates_t(0,S,6,U,"1S4"));       //       1S4	6 222 40
            U.set(3,0); states.push_back(fstates_t(0,P,6,U,"1P"));        //       1P	6 222 30
            U.set(2,0); states.push_back(fstates_t(0,D,2,U,"1D1"));       //       1D1	2 200 20
            U.set(2,0); states.push_back(fstates_t(0,D,4,U,"1D2"));       //       1D2	4 220 20
            U.set(2,1); states.push_back(fstates_t(0,D,4,U,"1D3"));       //       1D3	4 220 21
            U.set(2,2); states.push_back(fstates_t(0,D,4,U,"1D4"));       //       1D4	4 220 22
            U.set(2,0); states.push_back(fstates_t(0,D,6,U,"1D5"));       //       1D5	6 222 20
            U.set(4,0); states.push_back(fstates_t(0,D,6,U,"1D6"));       //       1D6	6 222 40
            U.set(2,1); states.push_back(fstates_t(0,F,4,U,"1F1"));       //       1F1	4 220 21
            U.set(1,0); states.push_back(fstates_t(0,F,6,U,"1F2"));       //       1F2	6 222 10
            U.set(3,0); states.push_back(fstates_t(0,F,6,U,"1F3"));       //       1F3	6 222 30
            U.set(4,0); states.push_back(fstates_t(0,F,6,U,"1F4"));       //       1F4	6 222 40
            U.set(2,0); states.push_back(fstates_t(0,G,2,U,"1G1"));       //       1G1	2 200 20
            U.set(2,0); states.push_back(fstates_t(0,G,4,U,"1G2"));       //       1G2	4 220 20
            U.set(2,1); states.push_back(fstates_t(0,G,4,U,"1G3"));       //       1G3	4 220 21
            U.set(2,2); states.push_back(fstates_t(0,G,4,U,"1G4"));       //       1G4	4 220 22
            U.set(2,0); states.push_back(fstates_t(0,G,6,U,"1G5"));       //       1G5	6 222 20
            U.set(3,0); states.push_back(fstates_t(0,G,6,U,"1G6"));       //       1G6	6 222 30
            U.set(4,0); states.push_back(fstates_t(0,G,6,U,"1G7"));       //       1G7	6 222 40A
            U.set(4,0); states.push_back(fstates_t(0,Gp,6,U,"1G8",1));    //       1G8	6 222 40B
            U.set(2,1); states.push_back(fstates_t(0,H,4,U,"1H1"));       //       1H1	4 220 21
            U.set(2,2); states.push_back(fstates_t(0,H,4,U,"1H2"));       //       1H2	4 220 22
            U.set(3,0); states.push_back(fstates_t(0,H,6,U,"1H3"));       //       1H3	6 222 30
            U.set(4,0); states.push_back(fstates_t(0,H,6,U,"1H4"));       //       1H4	6 222 40
            U.set(2,0); states.push_back(fstates_t(0,I,2,U,"1I1"));       //       1I1	2 200 20
            U.set(2,0); states.push_back(fstates_t(0,I,4,U,"1I2"));       //       1I2	4 220 20
            U.set(2,2); states.push_back(fstates_t(0,I,4,U,"1I3"));       //       1I3	4 220 22
            U.set(2,0); states.push_back(fstates_t(0,I,6,U,"1I4"));       //       1I4	6 222 20
            U.set(3,0); states.push_back(fstates_t(0,I,6,U,"1I5"));       //       1I5	6 222 30
            U.set(4,0); states.push_back(fstates_t(0,I,6,U,"1I6"));       //       1I6	6 222 40A
            U.set(4,0); states.push_back(fstates_t(0,Ip,6,U,"1I7",1));    //       1I7	6 222 40B
            U.set(2,1); states.push_back(fstates_t(0,K,4,U,"1K1"));       //       1K1	4 220 21
            U.set(3,0); states.push_back(fstates_t(0,K,6,U,"1K2"));       //       1K2	6 222 30
            U.set(4,0); states.push_back(fstates_t(0,K,6,U,"1K3"));       //       1K3	6 222 40
            U.set(2,1); states.push_back(fstates_t(0,L,4,U,"1L1"));       //       1L1	4 220 21
            U.set(2,2); states.push_back(fstates_t(0,L,4,U,"1L2"));       //       1L2	4 220 22
            U.set(4,0); states.push_back(fstates_t(0,L,6,U,"1L3"));       //       1L3	6 222 40A
            U.set(4,0); states.push_back(fstates_t(0,Lp,6,U,"1L4",1));    //       1L4	6 222 40B
            U.set(3,0); states.push_back(fstates_t(0,M,6,U,"1M1"));       //       1M1	6 222 30
            U.set(4,0); states.push_back(fstates_t(0,M,6,U,"1M2"));       //       1M2	6 222 40
            U.set(2,2); states.push_back(fstates_t(0,N,4,U,"1N1"));       //       1N1	4 220 22
            U.set(4,0); states.push_back(fstates_t(0,N,6,U,"1N2"));       //       1N2	6 222 40
            U.set(4,0); states.push_back(fstates_t(0,Q,6,U,"1Q"));        //       1Q	6 222 40
   	    break;
         case 7:
	    states.reserve(119);
            U.set(0,0); states.push_back(fstates_t(7,S,7,U,"8S"));        // f7    8S	7 000 00
            U.set(1,1); states.push_back(fstates_t(5,P,5,U,"6P"));        //       6P	5 110 11
            U.set(2,0); states.push_back(fstates_t(5,D,7,U,"6D"));        //       6D	7 200 20
            U.set(1,0); states.push_back(fstates_t(5,F,5,U,"6F"));        //       6F	5 110 10
            U.set(2,0); states.push_back(fstates_t(5,G,7,U,"6G"));        //       6G	7 200 20
            U.set(1,1); states.push_back(fstates_t(5,H,5,U,"6H"));        //       6H	5 110 11
            U.set(2,0); states.push_back(fstates_t(5,I,7,U,"6I"));        //       6I	7 200 20
            U.set(0,0); states.push_back(fstates_t(3,S,3,U,"4S1"));       //       4S1	3 111 00
            U.set(2,2); states.push_back(fstates_t(3,S,7,U,"4S2"));       //       4S2	7 220 22
            U.set(1,1); states.push_back(fstates_t(3,P,5,U,"4P1"));       //       4P1	5 211 11
            U.set(3,0); states.push_back(fstates_t(3,P,5,U,"4P2"));       //       4P2	5 211 30
            U.set(2,0); states.push_back(fstates_t(3,D,3,U,"4D1"));       //       4D1	3 111 20
            U.set(2,0); states.push_back(fstates_t(3,D,5,U,"4D2"));       //       4D2	5 211 20
            U.set(2,1); states.push_back(fstates_t(3,D,5,U,"4D3"));       //       4D3	5 211 21
            U.set(2,0); states.push_back(fstates_t(3,D,7,U,"4D4"));       //       4D4	7 220 20
            U.set(2,1); states.push_back(fstates_t(3,D,7,U,"4D5"));       //       4D5	7 220 21
            U.set(2,2); states.push_back(fstates_t(3,D,7,U,"4D6"));       //       4D6	7 220 22
            U.set(1,0); states.push_back(fstates_t(3,F,3,U,"4F1"));       //       4F1	3 111 10
            U.set(1,0); states.push_back(fstates_t(3,F,5,U,"4F2"));       //       4F2	5 211 10
            U.set(2,1); states.push_back(fstates_t(3,F,5,U,"4F3"));       //       4F3	5 211 21
            U.set(3,0); states.push_back(fstates_t(3,F,5,U,"4F4"));       //       4F4	5 211 30
            U.set(2,1); states.push_back(fstates_t(3,F,7,U,"4F5"));       //       4F5	7 220 21
            U.set(2,0); states.push_back(fstates_t(3,G,3,U,"4G1"));       //       4G1	3 111 20
            U.set(2,0); states.push_back(fstates_t(3,G,5,U,"4G2"));       //       4G2	5 211 20
            U.set(2,1); states.push_back(fstates_t(3,G,5,U,"4G3"));       //       4G3	5 211 21
            U.set(3,0); states.push_back(fstates_t(3,G,5,U,"4G4"));       //       4G4	5 211 30
            U.set(2,0); states.push_back(fstates_t(3,G,7,U,"4G5"));       //       4G5	7 220 20
            U.set(2,1); states.push_back(fstates_t(3,G,7,U,"4G6"));       //       4G6	7 220 21
            U.set(2,2); states.push_back(fstates_t(3,G,7,U,"4G7"));       //       4G7	7 220 22
            U.set(1,1); states.push_back(fstates_t(3,H,5,U,"4H1"));       //       4H1	5 211 11
            U.set(2,1); states.push_back(fstates_t(3,H,5,U,"4H2"));       //       4H2	5 211 21
            U.set(3,0); states.push_back(fstates_t(3,H,5,U,"4H3"));       //       4H3	5 211 30
            U.set(2,1); states.push_back(fstates_t(3,H,7,U,"4H4"));       //       4H4	7 220 21
            U.set(2,2); states.push_back(fstates_t(3,H,7,U,"4H5"));       //       4H5	7 220 22
            U.set(2,0); states.push_back(fstates_t(3,I,3,U,"4I1"));       //       4I1	3 111 20
            U.set(2,0); states.push_back(fstates_t(3,I,5,U,"4I2"));       //       4I2	5 211 20
            U.set(3,0); states.push_back(fstates_t(3,I,5,U,"4I3"));       //       4I3	5 211 30
            U.set(2,0); states.push_back(fstates_t(3,I,7,U,"4I4"));       //       4I4	7 220 20
            U.set(2,2); states.push_back(fstates_t(3,I,7,U,"4I5"));       //       4I5	7 220 22
            U.set(2,1); states.push_back(fstates_t(3,K,5,U,"4K1"));       //       4K1	5 211 21
            U.set(3,0); states.push_back(fstates_t(3,K,5,U,"4K2"));       //       4K2	5 211 30
            U.set(2,1); states.push_back(fstates_t(3,K,7,U,"4K3"));       //       4K3	7 220 21
            U.set(2,1); states.push_back(fstates_t(3,L,5,U,"4L1"));       //       4L1	5 211 21
            U.set(2,1); states.push_back(fstates_t(3,L,7,U,"4L2"));       //       4L2	7 220 21
            U.set(2,2); states.push_back(fstates_t(3,L,7,U,"4L3"));       //       4L3	7 220 22
            U.set(3,0); states.push_back(fstates_t(3,M,5,U,"4M"));        //       4M	5 211 30
            U.set(2,2); states.push_back(fstates_t(3,N,7,U,"4N"));        //       4N	7 220 22
            U.set(0,0); states.push_back(fstates_t(1,S,7,U,"2S1"));       //       2S1	7 222 00
            U.set(4,0); states.push_back(fstates_t(1,S,7,U,"2S2"));       //       2S2	7 222 40
            U.set(1,1); states.push_back(fstates_t(1,P,3,U,"2P1"));       //       2P1	3 210 11
            U.set(1,1); states.push_back(fstates_t(1,P,5,U,"2P2"));       //       2P2	5 221 11
            U.set(3,0); states.push_back(fstates_t(1,P,5,U,"2P3"));       //       2P3	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,P,5,U,"2P4"));       //       2P4	5 221 31
            U.set(3,0); states.push_back(fstates_t(1,P,7,U,"2P5"));       //       2P5	7 222 30
            U.set(2,0); states.push_back(fstates_t(1,D,3,U,"2D1"));       //       2D1	3 210 20
            U.set(2,1); states.push_back(fstates_t(1,D,3,U,"2D2"));       //       2D2	3 210 21
            U.set(2,0); states.push_back(fstates_t(1,D,5,U,"2D3"));       //       2D3	5 221 20
            U.set(2,1); states.push_back(fstates_t(1,D,5,U,"2D4"));       //       2D4	5 221 21
            U.set(3,1); states.push_back(fstates_t(1,D,5,U,"2D5"));       //       2D5	5 221 31
            U.set(2,0); states.push_back(fstates_t(1,D,7,U,"2D6"));       //       2D6	7 222 20
            U.set(4,0); states.push_back(fstates_t(1,D,7,U,"2D7"));       //       2D7	7 222 40
            U.set(1,0); states.push_back(fstates_t(1,F,1,U,"2F1"));       //       2F1	1 100 10
            U.set(2,1); states.push_back(fstates_t(1,F,3,U,"2F2"));       //       2F2	3 210 21
            U.set(1,0); states.push_back(fstates_t(1,F,5,U,"2F3"));       //       2F3	5 221 10
            U.set(2,1); states.push_back(fstates_t(1,F,5,U,"2F4"));       //       2F4	5 221 21
            U.set(3,0); states.push_back(fstates_t(1,F,5,U,"2F5"));       //       2F5	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,F,5,U,"2F6"));       //       2F6	5 221 31A
            U.set(3,1); states.push_back(fstates_t(1,Fp,5,U,"2F7",1));    //       2F7	5 221 31B
            U.set(1,0); states.push_back(fstates_t(1,F,7,U,"2F8"));       //       2F8	7 222 10
            U.set(3,0); states.push_back(fstates_t(1,F,7,U,"2F9"));       //       2F9	7 222 30
            U.set(4,0); states.push_back(fstates_t(1,F,7,U,"2F10"));      //       2F10	7 222 40
            U.set(2,0); states.push_back(fstates_t(1,G,3,U,"2G1"));       //       2G1	3 210 20
            U.set(2,1); states.push_back(fstates_t(1,G,3,U,"2G2"));       //       2G2	3 210 21
            U.set(2,0); states.push_back(fstates_t(1,G,5,U,"2G3"));       //       2G3	5 221 20
            U.set(2,1); states.push_back(fstates_t(1,G,5,U,"2G4"));       //       2G4	5 221 21
            U.set(3,0); states.push_back(fstates_t(1,G,5,U,"2G5"));       //       2G5	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,G,5,U,"2G6"));       //       2G6	5 221 31
            U.set(2,0); states.push_back(fstates_t(1,G,7,U,"2G7"));       //       2G7	7 222 20
            U.set(3,0); states.push_back(fstates_t(1,G,7,U,"2G8"));       //       2G8	7 222 30
            U.set(4,0); states.push_back(fstates_t(1,G,7,U,"2G9"));       //       2G9	7 222 40A
            U.set(4,0); states.push_back(fstates_t(1,Gp,7,U,"2G10",1));   //       2G10	7 222 40B
            U.set(1,1); states.push_back(fstates_t(1,H,3,U,"2H1"));       //       2H1	3 210 11
            U.set(2,1); states.push_back(fstates_t(1,H,3,U,"2H2"));       //       2H2	3 210 21
            U.set(1,1); states.push_back(fstates_t(1,H,5,U,"2H3"));       //       2H3	5 221 11
            U.set(2,1); states.push_back(fstates_t(1,H,5,U,"2H4"));       //       2H4	5 221 21
            U.set(3,0); states.push_back(fstates_t(1,H,5,U,"2H5"));       //       2H5	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,H,5,U,"2H6"));       //       2H6	5 221 31A
            U.set(3,1); states.push_back(fstates_t(1,Hp,5,U,"2H7",1));    //       2H7	7 222 31B
            U.set(3,0); states.push_back(fstates_t(1,H,7,U,"2H8"));       //       2H8	7 222 30
            U.set(4,0); states.push_back(fstates_t(1,H,7,U,"2H9"));       //       2H9	7 222 40
            U.set(2,0); states.push_back(fstates_t(1,I,3,U,"2I1"));       //       2I1	3 210 20
            U.set(2,0); states.push_back(fstates_t(1,I,5,U,"2I2"));       //       2I2	5 221 20
            U.set(3,0); states.push_back(fstates_t(1,I,5,U,"2I3"));       //       2I3	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,I,5,U,"2I4"));       //       2I4	5 221 31A
            U.set(3,1); states.push_back(fstates_t(1,Ip,5,U,"2I5",1));    //       2I5	5 221 31B
            U.set(2,0); states.push_back(fstates_t(1,I,7,U,"2I6"));       //       2I6	7 222 20
            U.set(3,0); states.push_back(fstates_t(1,I,7,U,"2I7"));       //       2I7	7 222 30
            U.set(4,0); states.push_back(fstates_t(1,I,7,U,"2I8"));       //       2I8	7 222 40A
            U.set(4,0); states.push_back(fstates_t(1,Ip,7,U,"2I9",1));    //       2I9	7 222 40B
            U.set(2,1); states.push_back(fstates_t(1,K,3,U,"2K1"));       //       2K1	3 210 21
            U.set(2,1); states.push_back(fstates_t(1,K,5,U,"2K2"));       //       2K2	5 221 21
            U.set(3,0); states.push_back(fstates_t(1,K,5,U,"2K3"));       //       2K3	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,K,5,U,"2K4"));       //       2K4	5 221 31A
            U.set(3,1); states.push_back(fstates_t(1,Kp,5,U,"2K5",1));    //       2K5	5 221 31B
            U.set(3,0); states.push_back(fstates_t(1,K,7,U,"2K6"));       //       2K6	7 222 30
            U.set(4,0); states.push_back(fstates_t(1,K,7,U,"2K7"));       //       2K7	7 222 40
            U.set(2,1); states.push_back(fstates_t(1,L,3,U,"2L1"));       //       2L1	3 210 21
            U.set(2,1); states.push_back(fstates_t(1,L,5,U,"2L2"));       //       2L2	5 221 21
            U.set(3,1); states.push_back(fstates_t(1,L,5,U,"2L3"));       //       2L3	5 221 31
            U.set(4,0); states.push_back(fstates_t(1,L,7,U,"2L4"));       //       2L4	7 222 40A
            U.set(4,0); states.push_back(fstates_t(1,Lp,7,U,"2L5",1));    //       2L5	7 222 40B
            U.set(3,0); states.push_back(fstates_t(1,M,5,U,"2M1"));       //       2M1	5 221 30
            U.set(3,1); states.push_back(fstates_t(1,M,5,U,"2M2"));       //       2M2	5 221 31
            U.set(3,0); states.push_back(fstates_t(1,M,7,U,"2M3"));       //       2M3	7 222 30
            U.set(4,0); states.push_back(fstates_t(1,M,7,U,"2M4"));       //       2M4	7 222 40
            U.set(3,1); states.push_back(fstates_t(1,N,5,U,"2N1"));       //       2N1	5 221 31
            U.set(4,0); states.push_back(fstates_t(1,N,7,U,"2N2"));       //       2N2	7 222 40
            U.set(3,1); states.push_back(fstates_t(1,O,5,U,"2O"));        //       2O	5 221 31
            U.set(4,0); states.push_back(fstates_t(1,Q,7,U,"2Q"));        //       2Q	7 222 40
   	    break;
      } // switch(n)
   }    // else (n>0 && n<15)
   break; // case D:
   default:
      std::cerr << "fconf::fconf() - error, only the case of l=0,1,2, and 3, s-, p-, d- and f-electrons implemented.\n";
   }    // switch(l)
}

fconf::fconf(int n, bool mJflag, orbital l)
{
   if(l!=S && l!=P && l!=D && l!=F) { 
      std::cerr << "fconf::fconf() - error, only the case of l=0,1,2, and 3, s-, p-, d- and f-electrons implemented.\n"; return; }
     
   fconf confLS(n,l);
   int num_states = (int)confLS.states.size();
   int i,j,J2min,J2max,mj;
   char Jlabel[12]; std::string id;

   for(i=0; i<num_states; i++)
   {
      J2min = abs(abs(2*confLS.states[i].L) - confLS.states[i].S2);
      J2max =     abs(2*confLS.states[i].L) + confLS.states[i].S2;
      for(j=J2min; j<=J2max; j+=2)
      {
         if(mJflag)
         {
	    for(mj=-j; mj<=j; mj+=2)
            {
               id.assign(confLS.states[i].id); id.append("_");
               if(j%2==0) sprintf(Jlabel,"%i",j/2); else sprintf(Jlabel,"%i/2",j); id.append(Jlabel);
               if(mj%2==0) sprintf(Jlabel,",mJ=%i",mj/2); else sprintf(Jlabel,",mJ=%i/2",mj); id.append(Jlabel);
               if(l==D) states.push_back(fstates_t(confLS.states[i].S2,confLS.states[i].L,confLS.states[i].v,id,j,mj));
               else
               states.push_back(fstates_t(confLS.states[i].S2,
                                          confLS.states[i].L,
                                          confLS.states[i].v,
                                          confLS.states[i].U,
                                          id,
                                          j,mj));
            }
         }
	 else
	 {
            id.assign(confLS.states[i].id); id.append("_");
            if(j%2==0) sprintf(Jlabel,"%i",j/2); else sprintf(Jlabel,"%i/2",j); id.append(Jlabel);
            if(l==D) states.push_back(fstates_t(confLS.states[i].S2,confLS.states[i].L,confLS.states[i].v,id,j));
            else
            states.push_back(fstates_t(confLS.states[i].S2,
                                       confLS.states[i].L,
                                       confLS.states[i].v,
                                       confLS.states[i].U,
                                       id,
                                       j));
         }
      }
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile: g++ states.cpp; ./a.out 1
// --------------------------------------------------------------------------------------------------------------- //
/*int main(int argc, char *argv[])
{
   int n = atoi(argv[1]);
   fconf conf(n);
   bool WisW;

   // NB. As I'm switching this project entirely to C++, so I'm moving from <cstdio> to <iostream>

   //fprintf(stdout,"First state has 2S=%i, L=%i, v=%i, U=[%i %i], W=[%i %i %i]\n",
   //        conf.states[0].S2,conf.states[0].L,conf.states[0].v,conf.states[0].U.u1,conf.states[0].U.u2,
   //        conf.states[0].W.w1,conf.states[0].W.w2,conf.states[0].W.w3);
   std::cout << "First state has 2S=" << conf.states[0].S2 << ", L=" << conf.states[0].L
             << ", U=[" << conf.states[0].U.u1 << " " << conf.states[0].U.u2
	     << "], W=[" << conf.states[0].W.w1 << " " << conf.states[0].W.w2 << " " << conf.states[0].W.w3 << "]\n";

   WisW = conf.states[0].W.isequal("100");
   //fprintf(stdout,"W=[%i %i %i] == [1 0 0] is %s\n",conf.states[0].W.w1,conf.states[0].W.w2,conf.states[0].W.w3,
   //        WisW ? "true" : "false");
   std::cout << "W=[" << conf.states[0].W.w1 << " " << conf.states[0].W.w2 << " " << conf.states[0].W.w3
             << "] == [1 0 0] is " << (WisW ? "true" : "false") << "\n";

   fconf confJ(n,false);
   std::cout << "2*J=" << confJ.states[0].J2 << "\n";

   fconf confmJ(n,true);
   std::cout << "2*mJ=" << confJ.states[0].mJ2 << "\n";

   return 0;
}*/
