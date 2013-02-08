/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                                   spline c
 
 
 
-----------------------------------------------------------------------------*/
 
/*----------------------------------------------------------------------------
Includedateien holen
-----------------------------------------------------------------------------*/
#include <stdio.h>          /* damit FILE definiert wird               */
#include <stdlib.h>
#include <math.h>           /* damit sqrt in define_j.c definiert wird */
#define pi (4.0*atan(1.0))  /* atan() braucht <math.h>                 */
#include "types.c"          /* benutze Datentypen laden                */
 
#define ZERO 0.
#define ONE  1.
#define SQRT sqrt
#define PRINTF printf
#define IF if
#define FOR for
#define ELSE else
#define GOTO goto
#define RETURN return
#define CONTINUE ;}
#define STOP exit(1);
#define END
 
 
#define IFI(I,A,B,C) {IF((I)<0)GOTO A;IF((I)==0)GOTO B;IF((I)>0)GOTO C;}
#define IFD(I,A,B,C) {IF((I)<0.)GOTO A;IF((I)==0.)GOTO B;IF((I)>0.)GOTO C;}
#define DO(A,I,IA,IE) FOR( (I)=IA;(I)<=IE;++(I) ){
#define GOTO5(A,B,C,D,M) {IF((M)==1)GOTO A;IF((M)==2)GOTO B;IF((M)==3)GOTO C;IF((M)==4)GOTO D;}
#define GOTO4(A,B,C,M)   {IF((M)==1)GOTO A;IF((M)==2)GOTO B;IF((M)==3)GOTO C;}
#define GOTO3(A,B,M)     {IF((M)==1)GOTO A;IF((M)==2)GOTO B;}
#define SIGN(a1,a2) ( (a2)>=0. ? ABSD(a1) : (-ABSD(a1)) )
 
extern INT      *sort();         /* definiert in DIAHERMX.C */
 
void TB04A(INT N, DOUBLE *XX, DOUBLE *FF, DOUBLE *DD, DOUBLE *AA);

/*----------------------------------------------------------------------------
                               spline()
-----------------------------------------------------------------------------*/
DOUBLE spline(n,xx,ff,x,macheps)
INT n;
DOUBLE *xx,*ff,x,macheps;
{
   #define x(i) (*(xx+i))
   #define f(i) (*(ff+i))
 
   DOUBLE *dd,*aa,TG01B(),value;
   DOUBLE x1,x2,x3,f1,f2,f3,a,b,c,det;
   INT    *nummer,i;
 
   if( n==0 ) return(0.0);
   if( n==1 ) return(f(1));
   if( n==2 ){
       x1 = x(1); f1 = f(1); x2 = x(2); f2 = f(2);
       if( x1== x2 )  return(f1);
       if( x<MIN(x1,x2) || x>MAX(x1,x2) ) return(0.0);
       a = (f2-f1)/(x2-x1);
       b = f1-a*x1;
       return(a*x+b);
   }
   if( n==3 ){
       x1 = x(1); f1 = f(1); x2 = x(2); f2 = f(2); x3 = x(3); f3 = f(3);
       if(x1==x2||x1==x3||x2==x3) return(0.0);
       if( x<MIN(MIN(x1,x2),x3) || x>MAX(MAX(x1,x2),x3) ) return(0.0);
       det=x1*x1*(x2-x3)-x2*x2*(x1-x3)+x3*x3*(x1-x2);
 a = (x2      -x3      )*f1-(x1      -x3      )*f2+(x1      -x2      )*f3;
 b =-(x2*x2   -x3*x3   )*f1+(x1*x1   -x3*x3   )*f2-(x1*x1   -x2*x2   )*f3;
 c =+(x2*x2*x3-x3*x3*x2)*f1-(x1*x1*x3-x3*x3*x1)*f2+(x1*x1*x2-x2*x2*x1)*f3;
       return((a*x*x+b*x+c)/det);
   }
 
 
   dd = DOUBLE_ALLOC(  n);
   aa = DOUBLE_ALLOC(3*n);
   /* ordne nach x(n) >...> x(1) */
   nummer = INT_ALLOC(n);
   nummer = sort(xx,nummer,n);
   for(i=1;i<=n;++i){
       *(aa+i) = *(xx+*(nummer+i));
       *(dd+i) = *(ff+*(nummer+i));
   }
   for(i=1;i<=n;++i){
       *(xx+i) = *(aa+*(nummer+i));
       *(ff+i) = *(dd+*(nummer+i));
   }
 
   TB04A(n,xx,ff,dd,aa);
   value=TG01B(-1,n,xx,ff,dd,x,macheps);
   free_(aa);
   free_(dd);
   free_(nummer);
 
   return(value);
}
/*----------------------------------------------------------------------------
                               TB04A()
-----------------------------------------------------------------------------*/
/* C/     ADD NAME=TB04A   HSL                     SINGLE           */
/* C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.          */
/* C######ALIAS TB04A                                               */
/*    SUBROUTINE TB04A (N, X, F, D, A)                              */
/*    DIMENSION X(N),F(N),D(N),A(*)                                 */
/*    DATA NP/6/                                                    */
/*C F(I) ARE THE FUNCTION VALUES AT THE POINTS X(I) FOR I=1,N AND   */
/*C THE SPLINE DERIVATIVES D(I) ARE FOUND.  THE DIMENSION OF A MUST */
/*C NOT BE LESS THAN 3*N. PERIPHERAL NP MUST BE AN OUTPUT MEDIUM.   */
 
void TB04A(N,XX,FF,DD,AA)
INT N;
DOUBLE *XX,*FF,*DD,*AA;
{
      #define X(I) (*(XX+I))
      #define F(I) (*(FF+I))
      #define D(I) (*(DD+I))
      #define A(I) (*(AA+I))
 
      INT I,J,K;
      DOUBLE H1,H2,P;
 
      DO(5,I,2,N)
      IFD(X(I)-X(I-1),C1,C1,C5);
C1:   PRINTF(" RETURN FROM TB04A IN SPLINE.C : \n");
      PRINTF(" THE ARRAY X IS OUT OF ORDER ! \n");
      PRINTF(" LOOK TO THE %d-ELEMENT. \n",I);
/*    WRITE(NP,3)I                                                */
/*    FORMAT(28H RETURN FROM TB04A BECAUSE X,I3,13H OUT OF ORDER) */
      A(1)=1.;
      RETURN;
C5:   CONTINUE
      DO(30,I,1,N)
      J=2;
      IFI(I-1,C6,C10,C6);
C6:   J=N-1;
      IF(I==N) GOTO C10;
      H1=1./(X(I)-X(I-1));
      H2=1./(X(I+1)-X(I));
      A(3*I-2)=H1;
      A(3*I-1)=2.*(H1+H2);
      A(3*I)=H2;
      D(I)=3.0*(F(I+1)*H2*H2+F(I)*(H1*H1-H2*H2)-F(I-1)*H1*H1);
      GOTO C30;
C10:  H1=1./(X(J)-X(J-1));
      H2=1./(X(J+1)-X(J));
      A(3*I-2)=H1*H1;
      A(3*I-1)=H1*H1-H2*H2;
      A(3*I)=-H2*H2;
      D(I)=2.*(F(J)*(H2*H2*H2+H1*H1*H1)-F(J+1)*H2*H2*H2-F(J-1)*H1*H1*H1);
C30:  CONTINUE
      P=A(4)/A(1);
      A(5)=A(5)-P*A(2);
      A(6)=A(6)-P*A(3);
      D(2)=D(2)-P*D(1);
      DO(50,I,3,N)
      K=3*I-4;
      P=A(K+2)/A(K);
      A(K+3)=A(K+3)-P*A(K+1);
      D(I)=D(I)-P*D(I-1);
      IF(I != N-1) GOTO C50;
      P=A(K+5)/A(K);
      A(K+5)=A(K+6)-P*A(K+1);
      A(K+6)=A(K+7);
      D(N)=D(N)-P*D(N-2);
C50:  CONTINUE
      D(N)=D(N)/A(3*N-1);
      DO(60,I,3,N)
      J=N+2-I;
      D(J)=(D(J)-A(3*J)*D(J+1))/A(3*J-1);
/*C60:*/  CONTINUE
      D(1)=(D(1)-D(2)*A(2)-D(3)*A(3))/A(1);
      A(1)=0.;
      RETURN;
      END
}
/*----------------------------------------------------------------------------
                               TG01B()
-----------------------------------------------------------------------------*/
/* C/     ADD NAME=TG01B   HSL                     SINGLE         */
/* C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.        */
/* C######ALIAS TG01B                                             */
/* REAL FUNCTION TG01B(IX,N,U,S,D,X)                              */
DOUBLE TG01B(IX,N,UU,SSS,DD,X,EPS)
INT IX,N;
DOUBLE *UU,*SSS,*DD,X,EPS;
{
      #define U(I) (*(UU+I))
      #define S(I) (*(SSS+I))
      #define D(I) (*(DD+I))
 
      INT /* I,*/IFLG=0,J;
      DOUBLE H,Q1,Q2,SS,B,A,Z;
/*
 C
 C**********************************************************************
 C
 C      TG01B -  FUNCTION ROUTINE TO EVALUATE A CUBIC SPLINE GIVEN SPLINE
 C     VALUES AND FIRST DERIVATIVE VALUES AT THE GIVEN KNOTS.
 C
 C     THE SPLINE VALUE IS DEFINED AS ZERO OUTSIDE THE KNOT RANGE,WHICH
 C     IS EXTENDED BY A ROUNDING ERROR FOR THE PURPOSE.
 C
 C                  F = TG01B(IX,N,U,S,D,X)
 C
 C       IX    ALLOWS CALLER TO TAKE ADVANTAGE OF SPLINE PARAMETERS SET
 C             ON A PREVIOUS CALL IN CASES WHEN X POINT FOLLOWS PREVIOUS
 C             X POINT. IF IX < 0 THE WHOLE RANGE IS SEARCHED FOR KNOT
 C             INTERVAL; IF IX > 0 IT IS ASSUMED THAT X IS GREATER THAN
 C             THE X OF THE PREVIOUS CALL AND SEARCH STARTED FROM THERE.
 C       N     THE NUMBER OF KNOTS.
 C       U     THE KNOTS.
 C       S     THE SPLINE VALUES.
 C       D     THE FIRST DERIVATIVE VALUES OF THE SPLINE AT THE KNOTS.
 C       X     THE POINT AT WHICH THE SPLINE VALUE IS REQUIRED.
 C       F     THE VALUE OF THE SPLINE AT THE POINT X.
 C
 C                                      MODIFIED JULY 1970
 C
 C**********************************************************************
 C
 C    ALLOWABLE ROUNDING ERROR ON POINTS AT EXTREAMS OF KNOT RANGE
 C    IS EPS*MAX(|U(1)|,|U(N)|).
      DIMENSION U(1),S(1),D(1)
      SAVE
      DATA IFLG/0/
      EPS=FD05A(1)*2.0
*/
 
/*      TEST WETHER POINT IN RANGE. */
      IF(X <  U(1)) GOTO C990;
      IF(X >  U(N)) GOTO C991;
 
/*      JUMP IF KNOT INTERVAL REQUIRES RANDOM SEARCH. */
      IF(IX <  0 || IFLG == 0) GOTO C12;
/*      JUMP IF KNOT INTERVAL SAME AS LAST TIME.      */
      IF(X <= U(J+1)) GOTO C8;
/*      LOOP TILL INTERVAL FOUND.                     */
  C1: J=J+1;
 C11: IF(X >  U(J+1)) GOTO C1;
      GOTO C7;
 
/*      ESTIMATE KNOT INTERVAL BY ASSUMING EQUALLY SPACED KNOTS.  */
 C12: J=ABSD(X-U(1))/(U(N)-U(1))*(N-1)+1;
/*      ENSURE CASE X=U(N) GIVES J=N-1.  */
      J=MIN(J,N-1);
/*      INDICATE THAT KNOT INTERVAL INSIDE RANGE HAS BEEN USED. */
      IFLG=1;
/*      SEARCH FOR KNOT INTERVAL CONTAINING X.  */
      IF(X >= U(J)) GOTO C11;
  C2: J=J-1;
      IF(X <= U(J)) GOTO C2;
 
/*      CALCULATE SPLINE PARAMETERS FOR JTH INTERVAL.  */
  C7: H=U(J+1)-U(J);
      Q1=H*D(J);
      Q2=H*D(J+1);
      SS=S(J+1)-S(J);
      B=3.0*SS-2.0*Q1-Q2;
      A=Q1+Q2-2.0*SS;
 
/*      CALCULATE SPLINE VALUE.    */
  C8: Z=(X-U(J))/H;
      return(((A*Z+B)*Z+Q1)*Z+S(J));
/*      TEST IF X WITHIN ROUNDING ERROR OF U(1).     */
C990: IF(X <= U(1)-EPS*MAX(ABSD(U(1)),ABSD(U(N)))) GOTO C99;
      J=1;
      GOTO C7;
/*      TEST IF X WITHIN ROUNDING ERROR OF U(N).     */
C991: IF(X >= U(N)+EPS*MAX(ABSD(U(1)),ABSD(U(N)))) GOTO C99;
      J=N-1;
      GOTO C7;
 C99: IFLG=0;
/*      FUNCTION VALUE SET TO ZERO FOR POINTS OUTSIDE THE RANGE. */
      return(0.0);
      END
}
/*------------------------------------------------------------------------------
ENDEMODUL    S P L I N E   C
------------------------------------------------------------------------------*/
