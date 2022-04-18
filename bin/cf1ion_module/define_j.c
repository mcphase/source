#ifndef DEFINED_J
#define DEFINED_J
 
/*-----------------------------------------------------------------------------
 
                             I  N  C  L  U  D  E
 
                                  DEFINE_J C
 
-------------------------------------------------------------------------------
 
Aufgabe               : Drehimpulsoperatoren definieren
 
-----------------------------------------------------------------------------*/

#define J   ((DOUBLE)dimj - 1)/2   /* dimj = 2J+1 */
#define J2  J*(J+1)
#define J4  J2*J2
#define J6  J4*J2
 
#define mj  ((DOUBLE)m - J - 1)     /* m =2J+1,...,1  => mj =J,...,-J */
#define nj  ((DOUBLE)n - J - 1)     /* n =2J+1,...,1  => nj =J,...,-J */
 
#define mj2  mj*mj
#define mj3  mj*mj2
#define mj4  mj2*mj2
#define mj6  mj4*mj2
 
#define nj2  nj*nj
#define nj3  nj*nj2
#define nj4  nj2*nj2
#define nj6  nj4*nj2
 
#define JP(mj)   sqrt( J2 - (mj)*((mj)+1) )
#define JM(mj)   sqrt( J2 - (mj)*((mj)-1) )

#define JPMMP(mj)  JP(mj-1)*JM(mj)+JM(mj+1)*JP(mj)
 
#define JP1(mj)  JP(mj)
#define JM1(mj)  JM(mj)
 
#define JP2(mj)  JP(mj+1)*JP1(mj)
#define JM2(mj)  JM(mj-1)*JM1(mj)
 
#define JP3(mj)  JP(mj+2)*JP2(mj)
#define JM3(mj)  JM(mj-2)*JM2(mj)
 
#define JP4(mj)  JP(mj+3)*JP3(mj)
#define JM4(mj)  JM(mj-3)*JM3(mj)
 
#define JP5(mj)  JP(mj+4)*JP4(mj)
#define JM5(mj)  JM(mj-4)*JM4(mj)
 
#define JP6(mj)  JP(mj+5)*JP5(mj)
#define JM6(mj)  JM(mj-5)*JM5(mj)
 
#define S(a)    ( (a)>=0 ? 1 : -1 )
#define FAC(i)   fac(  ABS( (INT)(i) )  )
#define JZ(n,mj) sqrt(  FAC( J-S(n)*mj ) * FAC( J+S(n)*mj+ABS(n) )  )
#define JN(n,mj) sqrt(  FAC( J+S(n)*mj ) * FAC( J-S(n)*mj-ABS(n) )  )
#define JPM(n,mj)   ( JZ(n,mj)/JN(n,mj) )
 
#define D(nj,mj) (  ( (nj)==(mj) && ABSD(mj)<= J && ABSD(nj)<=J )? 1.0 : 0.0  )
 
 
#define P2(s)  ( (s)*(s) )
#define P3(s)  ( (s)*P2(s) )
#define P4(s)  ( (s)*P3(s) )
#define P5(s)  ( (s)*P4(s) )
#define P6(s)  ( (s)*P5(s) )
 
/*       m       ...    m                    v
     |    1,1            1,2J+1     |     |   1     |
     |    .              .          |     |   .     |
     |    .              .          |     |   .     |
     |    .              .          |     |   .     |
     |                              |     |  v      |
     |   m       ...    m           |   , |   2j+1  |
          2J+1,1         2J+1,2J+1
 
 
     dem entspricht formal :
 
 
         m       ...    m                   v
     |    -J,-J          -J,J      |     |   -J    |
     |    .              .         |     |   .     |
     |    .              .         |     |   .     |
     |    .              .         |     |   .     |
     |                             |     |  v      |
     |   m       ...    m          |   , |    J    |
          J,-J           J,J
*/
/*------------------------------------------------------------------------------
ENDEMODUL    D E F I N E _ J    C
------------------------------------------------------------------------------*/
#endif
