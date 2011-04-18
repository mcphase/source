/* njsyms.cpp
 *
 * Calculates values of the n-j symbols, by general explicit algebraic formulae
 *
 * Functions:
 *
 * long int factorial(int i);                                                   // Calculates the factorial operation !
 * double racahW(int a, int b, int c, int d, int e, int f);                     // Calculates Racah's W(abcd;ef) symbol
 * double wigner(int a, int b, int c, int d, int e, int f);                     // Calculates Wigner coeff. (abcd|ef)
 * double wigner(int a, int b, int c, int d, int e, int f, int g, int h);       // Calculates Wigner coeff. (abcd|abef)
 * double sixj(int a, int b, int c, int d, int e, int f);                       // Calculates the 6j symbol {abc;def}
 * double threej(int a, int b, int c, int d, int e, int f);                     // Calculates the 3j symbol (abd;def)
 * double ninej(int a, int b, int c, int d, int e, int f, int g, int h, int i); // Calculates the 9j symbol {abd;def;ghi}
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

// --------------------------------------------------------------------------------------------------------------- //
// NB. Note that all the functions in this file uses twice the value of any arguments - in order to represent
//     half-integral arguments as integers.
// --------------------------------------------------------------------------------------------------------------- //

#include<cmath>
#include<cstdlib>

// --------------------------------------------------------------------------------------------------------------- //
// Function to calculate the factorial operation
// --------------------------------------------------------------------------------------------------------------- //
//long int factorial(int i)
double factorial(int i)
{
   int l;
 //long int j=1;
   double j=1.;
   if(i==0) return 1;
   else if(i<0) return 0;
   else
      for(l=2; l<=i; l++)
         j*=l;
   return j;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the value of Racah's W symbols W(abcd;ef) which are related to the 6-j symbols
// --------------------------------------------------------------------------------------------------------------- //
double racahW(int a, int b, int c, int d, int e, int f)
{
   int z,minz,maxz;
   int F[] = {(a+b-e)/2, (c+d-e)/2, (a+c-f)/2, (b+d-f)/2, (e+f-a-d)/2, (e+f-b-c)/2};
   double sumz = 0.;
   // The formula for the Racah W function is:       where the triangle functions are:
   //                                           
   //                                        1/2                                        (a+b-c)!(a-b+c)!(-a+b+c)!
   // W(abcd;ef) = [T(abe)T(cde)T(acf)T(bdf)]                                  T(abc) = -------------------------
   //                                                                                       (a + b + c + 1)!
   //               ---                     (a + b + c + d + 1 - z)!
   //             * >    ------------------------------------------------------------------
   //               ---z z!(a+b-e-z)!(c+d-e-z)!(a+c-f-z)!(b+d-f-z)!(e+f-a-d+z)!(e+f-b-c+z)!
   //
   // The W functions are only defined for integer or half integer arguments that satisfy the selection rule that 
   // the four triads (abe), (cde), (acf), and (bdf) has an integral sum and satisfy the triangular inequality 
   // |a-b| <= c <= |a+b|. Otherwise W is zero.

   // Selection rules: The triads (a b e), (c d e), (a c f), (b d f)
   // 1. All satisfy the triangular inequality: |a-b| <= c <= a+b 
#define TRI(A,B,C) ( (C > (A+B)) || (C < abs(A-B)) )
   if(TRI(a,b,e) || TRI(c,d,e) || TRI(a,c,f) || TRI(b,d,f)) return 0.;
   // 2. Elements of each triad sum to an integer  [a,b,c,d,e,f are actually twice their nominal values]
   else if( ((a+b+e)%2!=0) || ((c+d+e)%2!=0) || ((a+c+f)%2!=0) || ((b+d+f)%2!=0) ) return 0.;
   
   minz = 0; for(z=4; z<6; z++) if(-F[z]>minz) minz = -F[z];        //z = max([0 -F(5:6)]):min([F(1:4) a+b+c+d+1]);
   maxz = a+b+c+d+1; for(z=0; z<4; z++) if(F[z]<maxz) maxz = F[z];
   for(z=minz; z<=maxz; z++)
      sumz += pow(-1.,z) * factorial((a+b+c+d)/2+1-z) / ( factorial(z)*factorial(F[0]-z)*factorial(F[1]-z) 
                          *factorial(F[2]-z)*factorial(F[3]-z)*factorial(F[4]+z)*factorial(F[5]+z) );
#define T(A,B,C) ( factorial((A+B-C)/2)*factorial((A-B+C)/2)*factorial((-A+B+C)/2) ) / factorial((A+B+C)/2+1)
   sumz *= sqrt(T(a,b,e)*T(c,d,e)*T(a,c,f)*T(b,d,f));
   return sumz;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the wigner coefficient (a b c d | a b e f) == (a  b  c  d | e  f) == (j1 j2 m1 m2 | j -m)
// --------------------------------------------------------------------------------------------------------------- //
double wigner(int a, int b, int c, int d, int e, int f)
{
   int t,mint,maxt;
   int F[] = {(e-b+c)/2, (e-a-d)/2, (a+b-e)/2, (a-c)/2, (b+d)/2};
   double sum_t = 0.;
   // The analytical expression for the Wigner coefficient is given by: (cf. Racah II)
   //
   // (j j m m |j j jm)      ___________________
   //   1 2 1 2  1 2    =  \/ (2j+1) T(j1,j2,j) 
   //                                          ______________________________________________
   //                        ---     z       \/ (j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(j+m)!(j-m)! 
   //                      * >   (-1)   ----------------------------------------------------------
   //                        ---        z!(j1+j2-j-z)!(j1-m1-z)!(j2+m2-z)!(j-j2+m1+z)!(j-j1-m2+z)!
   //                         z

   // Selection rules:
   // 1. m1+m2-m = 0
   if((c+d)!=f) return 0.;
   // 2. |j1-j2| <= j <= j1+j2
   else if(TRI(a,b,e)) return 0.;
   // 3. -j1 <= m1 <= j1, -j2 <= m2 <= j2, -j <= -m <= j
   else if((c>a||c<-a) || (d>b||d<-b) || (-f>e||-f<-e)) return 0.;
   // 4. Integer perimeter rule: j1 + j2 + J is integer
   else if((a+b+e)%2!=0) return 0.;

   mint = 0;    for(t=0; t<2; t++) if(-F[t]>mint) mint = -F[t];
   maxt = F[2]; for(t=3; t<5; t++) if( F[t]<maxt) maxt =  F[t];
   for(t=mint; t<=maxt; t++)
      sum_t += pow(-1.,(double)t) / ( factorial(t)*factorial(F[0]+t)*factorial(F[1]+t)*factorial(F[2]-t)
                                     *factorial(F[3]-t)*factorial(F[4]-t) );
   sum_t *= sqrt( factorial((a+c)/2)*factorial((a-c)/2)*factorial((b+d)/2)*factorial((b-d)/2)
                 *factorial((e+f)/2)*factorial((e-f)/2) );
   return sqrt((e+1.)*T(a,b,e)) * sum_t; 
}
double wigner(int a, int b, int c, int d, int e, int f, int g, int h)
{  if(e==a && f==b) return wigner(a,b,c,d,g,h); else return 0.; }

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the 6-j symbol.
// --------------------------------------------------------------------------------------------------------------- //
double sixj(int a, int b, int c, int d, int e, int f)
{
   int z,minz,maxz;
   int F[] = {(a+b+c)/2, (a+e+f)/2, (d+b+f)/2, (d+e+c)/2, (a+b+d+e)/2, (b+c+e+f)/2, (a+c+d+f)/2};
   double sumz = 0.;
   // The 6-j symbols can be computed by the Racah formula:
   //
   //                                       1/2
   // { a b c } = [T(abc)T(aef)T(dbf)T(dec)]
   // { d e f }
   //                                                     z
   //              ---                                (-1)  * (z+1)!
   //            * >    ----------------------------------------------------------------------------
   //              ---z (z-a-b-c)!(z-a-e-f)!(z-d-b-f)!(z-d-e-c)!(a+b+d+e-z)!(b+c+e+f-z)!(a+c+d+f-z)!

   // Selection rules: The triads (a b c), (a e f), (d b f), (d e c)
   // 1. All satisfy the triangular inequality: |a-b| <= c <= a+b 
   if(TRI(a,b,c) || TRI(a,e,f) || TRI(d,b,f) || TRI(d,e,c)) return 0.;
   // 2. Elements of each triad sum to an integer  [a,b,c,d,e,f are actually twice their nominal values]
   else if( ((a+b+c)%2!=0) || ((a+e+f)%2!=0) || ((d+b+f)%2!=0) || ((d+e+c)%2!=0) ) return 0.;

   minz = F[0]; for(z=1; z<4; z++) if(F[z]>minz) minz = F[z];   // for z = max(F(1:4)):min(F(5:7))
   maxz = F[4]; for(z=4; z<7; z++) if(F[z]<maxz) maxz = F[z];
   for(z=minz; z<=maxz; z++)
      sumz += ( pow(-1.,z) * factorial(z+1) ) / ( factorial(z-F[0])*factorial(z-F[1])*factorial(z-F[2])*factorial(z-F[3])
                                                 *factorial(F[4]-z)*factorial(F[5]-z)*factorial(F[6]-z) );
   return sqrt(T(a,b,c)*T(a,e,f)*T(d,b,f)*T(d,e,c)) * sumz;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the 3-j symbol.
// --------------------------------------------------------------------------------------------------------------- //
double threej(int a, int b, int c, int d, int e, int f)
{
   int t,mint,maxt;
   int F[] = {(c-b+d)/2, (c-a-e)/2, (a+b-c)/2, (a-d)/2, (b+e)/2};
   double sum_t = 0.;
   // The analytical expression for the Wigner 3j symbol is given by the Racah formula
   // ( j1 j2 J )        (j1-j2-M)  
   // ( m1 m2 M ) = (-1)^          sqrt(T(j1,j2,J)) sqrt((j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(J+M)!(J-M)!)
   //                                               t
   //               * SUM _____________________(-1)^________________________________
   //                  t  t!(J-j2+t+m1)!(J-j1+t-m2)!(j1+j2-J-t)!(j1-t-m1)!(j2-t+m2)!

   // Selection rules:
   // 1. d+e+f = 0
   if((d+e+f)!=0) return 0.;
   // 2. |a-b| <= c <= a+b
   else if(TRI(a,b,c)) return 0.;
   // 3. -a <= d <= a, -b <= e <= b, -c <= -f <= c
   else if((d>a||d<-a) || (e>b||e<-b) || (-f>c||-f<-c)) return 0.;
   // 4. Integer perimeter rule: j1 + j2 + J is integer
   else if((a+b+c)%2!=0) return 0.;

   mint = 0;    for(t=0; t<2; t++) if(-F[t]>mint) mint = -F[t];
   maxt = F[2]; for(t=2; t<5; t++) if( F[t]<maxt) maxt =  F[t];
   for(t=mint; t<=maxt; t++)
      sum_t += pow(-1.,t) / ( factorial(t)*factorial(F[0]+t)*factorial(F[1]+t)*factorial(F[2]-t)*factorial(F[3]-t)*factorial(F[4]-t) );
   return pow(-1.,(a-b-f)/2.) * sqrt(T(a,b,c) * factorial((a+d)/2)*factorial((a-d)/2)*factorial((b+e)/2)*factorial((b-e)/2) 
                                             * factorial((c+f)/2)*factorial((c-f)/2) ) * sum_t;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the 9j symbols
// --------------------------------------------------------------------------------------------------------------- //
double ninej(int j1, int j2, int J12, int j3, int j4, int J34, int J13, int J24, int J)
{
// The equation for the 9-j symbol is:
//
// { j1  j2  J12 }   ---     2g
// { j3  j4  J34 } = >   (-1)    (2g+1)  { j1  j2 J12 } { j3 j4 J34 } { J13 J24 J  }
// { J13 J24  J  }   ---                 { J34  J  g  } { j2 g  J24 } {  g  j1  j3 }
//                    g
//
// where the two-row curly brackets indicate the 6-j symbols. The 6-j symbols are equal
// to a permutation of their columns, so:
//
// { j3 j4 J34 } = { J34 j3 j4 }  and  { J13 J24 J  } = { J24 J  J13 }
// { j2 g  J24 }   { J24 j2 g  }       {  g  j1  j3 }   { j1  j3  g  }
//
// Now, the 6-j symbols { j1 j2 j3 } have triads (j1 j2 j3), (j1 J2 J3), (J1 j2 J3) and (J1 J2 j3) which satisfy the triangular
//                      { J1 J2 J3 }   inequalities. In particular, for those involving g [e.g. (j1 J2 J3)], |j1+J2| <= J3 <= j1+J2
//
// Thus for the 9j symbols we must have:   |J-j1|  <= g <=  J+j1
//                                        |J34-j2| <= g <= J34+j2
//                                        |J24-j3| <= g <= J24+j3

   int g,min_g,max_g;
   int G[] = {abs(J-j1), abs(J34-j2), abs(J24-j3), J+j1, J34+j2, J24+j3};
   double out = 0.;

   // Finds the allowed values of g
   min_g = G[0]; for(g=1; g<3; g++) if(G[g]<min_g) min_g = G[g];
   max_g = G[3]; for(g=4; g<6; g++) if(G[g]>max_g) max_g = G[g];
   for(g=min_g; g<=max_g; g++)  // g==2g, as all integers here represents twice their values (to accomodate half integral values).
      out += pow(-1.,g) * (g+1) * sixj(j1,j2,J12,J34,J,g) * sixj(j3,j4,J34,j2,g,J24) * sixj(J13,J24,J,g,j1,j3);

   return out;

}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile: g++ njsyms.cpp && ./a.out
// --------------------------------------------------------------------------------------------------------------- //
/*#include<iostream>
int main(int argc, char *argv[])
{
   int a,b,c,d,e,f,g,h,i;

   if(argc<7)
   {
      a = 3; b = 4; c = 1; d = 2; e = 1; f = 4;  // for racahW and sixj
   }
   else
   {
      a = atoi(argv[1])*2; b = atoi(argv[2])*2; c = atoi(argv[3])*2;
      d = atoi(argv[4])*2; e = atoi(argv[5])*2; f = atoi(argv[6])*2;
   }
   
   // Hmmm... problems with sizes of ints with factorials > 12 - had to modify factorial() to use long ints - again probs at >20!
   std::cout << "3! = " << factorial(3) << "; 0! = " << factorial(0) << "; 1! = " << factorial(1) << "\n";
   std::cout << "11! = " << factorial(11); std::cout << "; 12! = " << factorial(12) << "; 13! = " << factorial(13) << "\n";
   std::cout << "14! = " << factorial(14); std::cout << "; 15! = " << factorial(15) << "; 16! = " << factorial(16) << "\n";
   std::cout << "20! = " << factorial(20); std::cout << "; 21! = " << factorial(21) << "; 22! = " << factorial(22) << "\n";

   std::cout << "{" << a << "/2 " << b << "/2 " << c << "/2}\n" << "{" << d << "/2 " << e << "/2 " << f << "/2} = ";
   std::cout << sixj(a,b,c,d,e,f) << "\n\n";

   std::cout << "W(" << a << "/2 " << b << "/2 " << c << "/2 " << d << "; " << e << "/2 " << f << "/2) = " << racahW(a,b,c,d,e,f) << "\n";

   std::cout << "(" << a << "/2 " << b << "/2 " << c << "/2 " << d << "/2 | " << e << "/2 " << f << "/2) = " << wigner(a,b,c,d,e,f) << "\n\n";

   std::cout << "{" << a << "/2 " << b << "/2 " << c << "/2}\n" << "{" << e << "/2 " << d << "/2 " << c << "/2}\n" << "{" << f << "/2 " << f << "/2 " << " 0} = ";
   std::cout << ninej(a,b,c,e,d,c,f,f,0) << "\n\n";

   if(argc<7)
   {
      a = 1; b = 2; c = 1; d = 1; e = 0; f = -1;   // for threej
   }
   std::cout << "(" << a << "/2 " << b << "/2 " << c << "/2)\n" << "(" << d << "/2 " << e << "/2 " << f << "/2) = ";
   std::cout << threej(a,b,c,d,e,f) << "\n\n";

   return 0;
}*/
