/*     ADD NAME=MB11A   HSL                     SINGLE
C@PROCESS DIRECTIVE('IBMD')
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS MB11A  */
/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                                   MB11A C
 
 
 
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
#define RETURN return 0
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
 

INT MB11D(DOUBLE *A, INT ZA, INT SA, INT ZB, INT SB, INT ZC, INT SC, INT IA, DOUBLE WKK, INT MKK, INT NKK);
INT MB11E(DOUBLE *A, INT ZA, INT SA, INT ZB, INT SB, INT IA, DOUBLE *W, DOUBLE WKK, INT MKK, INT NKK);
INT MB11F(DOUBLE *A, INT ZA, INT SA, INT ZB, INT SB, INT ZC, INT SC, INT IA, DOUBLE WKK, DOUBLE AKK, INT MKK, INT NKK);

/*----------------------------------------------------------------------------
                               MB11A()
-----------------------------------------------------------------------------*/
 
/*    SUBROUTINE MB11A (M,N,A,IA,W)  */
      INT        MB11A (M,N,WW,PA,IA,PW)
      INT    M,N,IA,PA,PW;
      DOUBLE *WW;
{
      #define A(I,J) (*(WW+PA+(I-1)+(J-1)*IA))
      #define W(I)   (*(WW+PW+(I-1)         ))
      DOUBLE           BSQ,RMAX,SIGMA,SUM;
      DOUBLE           WKK,AKK;
      CHAR *TEXT1,*TEXT2;
      INT    NRW,NCW,N1,N2,N3,N4,N5,N6,I,J/*,K*/,IR,KK,MMK,KP;
printf("am in mb11a()\n");
      IR=0;
/*    COMMON/MB11B /LP
      DIMENSION A(IA,1),W(1)
      EXTERNAL MB11C
      DATA ONE/1.0E0/,ZERO/0.0E0/
C     PARTITION THE WORKING SPACE ARRAY W
C     THE FIRST PARTITION HOLDS THE FIRST COMPONENTS OF THE VECTORS OF
C     THE ELEMENTARY TRANSFORMATIONS */
      NRW=M;
/*    THE SECOND PARTITION RECORDS ROW INTERCHANGES */
      NCW=M+M;
/*    THE THIRD PARTITION RECORDS COLUMN INTERCHANGES
C     SET THE INITIAL RECORDS OF ROW AND COLUMN INTERCHANGES */
      DO(1,I,1,M)
      N1=NRW+I;
      W(N1)=(DOUBLE)0.5+(DOUBLE)I;
/* 1*/CONTINUE
      DO(2,I,1,N)
      N1=NCW+I;
      W(N1)=(DOUBLE)0.5+(DOUBLE)I;
/* 2*/CONTINUE
/*    'KK' COUNTS THE SEPARATE ELEMENTARY TRANSFORMATIONS  */
      KK=1;
/*    FIND LARGEST ROW AND MAKE ROW INTERCHANGES */
   C3:RMAX=ZERO;
      DO(4,I,KK,M)
      SUM=ZERO;
      DO(5,J,KK,N)
      SUM=SUM+A(I,J)*A(I,J);
/* 5*/CONTINUE
      IFD(RMAX-SUM,C6,C4,C4);
   C6:RMAX=SUM;
      IR=I;
   C4:CONTINUE
      IF(RMAX == ZERO) GOTO C81;
      IFI(IR-KK,C7,C7,C8);
   C8:N3=NRW+KK;
      SUM=W(N3);
      N4=NRW+IR;
      W(N3)=W(N4);
      W(N4)=SUM;
      DO(9,J,1,N)
      SUM=A(KK,J);
      A(KK,J)=A(IR,J);
      A(IR,J)=SUM;
/* 9*/CONTINUE
/*    FIND LARGEST ELEMENT OF PIVOTAL ROW, AND MAKE COLUMN INTERCHANGES*/
   C7:RMAX=ZERO;
      SUM=ZERO;
      DO(10,J,KK,N)
      SUM=SUM+A(KK,J)*A(KK,J);
      IFD(RMAX- ABSD(A(KK,J)),C11,C10,C10);
  C11:RMAX= ABSD(A(KK,J));
      IR=J;
  C10:CONTINUE
      IFI(IR-KK,C12,C12,C13);
  C13:N5=NCW+KK;
      RMAX=W(N5);
      N6=NCW+IR;
      W(N5)=W(N6);
      W(N6)=RMAX;
      DO(14,I,1,M)
      RMAX=A(I,KK);
      A(I,KK)=A(I,IR);
      A(I,IR)=RMAX;
/*14*/CONTINUE
/*    REPLACE THE PIVOTAL ROW BY THE VECTOR OF THE TRANSFORMATION*/
  C12:SIGMA= SQRT(SUM);
      BSQ= SQRT(SUM+SIGMA* ABSD(A(KK,KK)));
      W(KK)= SIGN(SIGMA+ ABSD(A(KK,KK)),A(KK,KK))/BSQ;
      A(KK,KK)=- SIGN(SIGMA,A(KK,KK));
      KP=KK+1;
      IFI(KP-N,C15,C15,C16);
  C15:DO(17,J,KP,N)
      A(KK,J)=A(KK,J)/BSQ;
/*17*/CONTINUE
/*    APPLY THE TRANSFORMATION TO THE REMAINING ROWS OF A */
      IFI(KP-M,C18,C18,C16)
  C18:WKK=W(KK);
/*    MB11D (&A(KK+1,KK+1),&A(KK+1,KK),&A(KK,KK+1),IA,WKK,M-KK,N-KK); */
      MB11D (&A(1,1),KK+1,KK+1,KK+1,KK,KK,KK+1, IA,WKK,M-KK,N-KK);
      KK=KP;
      GOTO C3;
/*    AT THIS STAGE THE REDUCTION OF A IS COMPLETE
C     NOW WE BUILD UP THE GENERALIZED INVERSE
C     APPLY THE FIRST ELEMENTARY TRANSFORMATION  */
  C16:KK=M;
      KP=M+1;
      SUM=-W(M)/A(M,M);
      IFI(N-M,C33,C33,C34);
  C34:DO(35,J,KP,N)
      A(M,J)=SUM*A(M,J);
/*35*/CONTINUE
  C33:A(M,M)=ONE/A(M,M)+SUM*W(M);
/*    NOW APPLY THE OTHER (M-1) TRANSFORMATIONS */
  C36:KP=KK;
      KK=KP-1;
      IFI(KK,C37,C37,C38);
/*    FIRST TRANSFORM THE LAST (M-KK) ROWS */
  C38:WKK=W(KK);
/*    MB11E(&A(KK+1,KK+1),&A(KK,KK+1),IA,&W(KK+1),WKK,M-KK,N-KK); */
      MB11E(&A(1,1),KK+1,KK+1 ,   KK,KK+1 ,IA,&W(KK+1),WKK,M-KK,N-KK);
/*    THEN CALCULATE THE NEW ROW IN POSITION KK */
      AKK=ONE/A(KK,KK);
/*    MB11F(&A(KK+1,KK+1),&A(KK,KK+1),&A(KK+1,KK),IA,WKK,AKK,M-KK,N-KK);*/
      MB11F(&A(1,1),KK+1,KK+1,KK,KK+1 ,   KK+1,KK ,IA,WKK,AKK,M-KK,N-KK);
/*    AND REVISE THE COLUMN IN POSITION KK  */
      SUM=ONE-WKK*WKK;
      DO(44,I,KP,M)
      SUM=SUM-A(I,KK)*W(I);
      A(I,KK)=W(I);
/*44*/CONTINUE
      A(KK,KK)=SUM/A(KK,KK);
      GOTO C36;
/*    RESTORE THE ROW INTERCHANGES  */
  C37:DO(45,I,1,M)
  C46:N1=NRW+I;
      IR=(INT)W(N1);
      IFI(I-IR,C47,C45,C45);
  C47:SUM=W(N1);
      N2=NRW+IR;
      W(N1)=W(N2);
      W(N2)=SUM;
      DO(48,J,1,N)
      SUM=A(I,J);
      A(I,J)=A(IR,J);
      A(IR,J)=SUM;
/*48*/CONTINUE
      GOTO C46;
  C45:CONTINUE
/*    RESTORE THE COLUMN INTERCHANGES */
      DO(49,J,1,N)
  C50:N1=NCW+J;
      IR=(INT)W(N1);
      IFI(J-IR,C51,C49,C49);
  C51:SUM=W(N1);
      N2=NCW+IR;
      W(N1)=W(N2);
      W(N2)=SUM;
      DO(52,I,1,M)
      SUM=A(I,J);
      A(I,J)=A(I,IR);
      A(I,IR)=SUM;
/*52*/CONTINUE
      GOTO C50;
  C49:CONTINUE
/*80*/RETURN;
  C81:MMK=M-KK;
 /*   WRITE(LP,82) MMK
   82 FORMAT(1H0,22H *** MB11A  ERROR *** ,I3,8H REDUCED,
     1 22H ROWS FOUND TO BE ZERO) */
      TEXT1=" *** MB11A  ERROR *** ";
      TEXT2=" REDUCED, ROWS FOUND TO BE ZERO) ";
      PRINTF("%s %d %s\n",TEXT1,MMK,TEXT2);
      STOP
      END
}
/*----------------------------------------------------------------------------
                               MB11C
-----------------------------------------------------------------------------*/
/*    BLOCK DATA MB11C
      COMMON/MB11B /LP
      DATA LP/6/
      END            */
/*----------------------------------------------------------------------------
                               MB11D
-----------------------------------------------------------------------------*/
/*    C@PROCESS DIRECTIVE('IBMD')  */
/*    SUBROUTINE MB11D (A,B,C,IA,WKK,MKK,NKK) */
      INT        MB11D (A ,ZA,SA,ZB,SB,ZC,SC, IA,WKK,MKK,NKK)
      DOUBLE *A ,WKK;
      INT    IA,MKK,NKK,ZA,SA,ZB,SB,ZC,SC;
{
/*    REAL             A(IA,NKK),B(*),C(IA,NKK),SUM,WKK */
      DOUBLE                                    SUM;
      INT I,J;
      #define AA(I,J)  (*(A +(ZA-1+I-1)+(SA-1+J-1)*IA))
      #define CC(I,J)  (*(A +(ZC-1+I-1)+(SC-1+J-1)*IA))
      #define BB(I)    (*(A +(ZB-1+I-1)+(SB-1    )*IA))
 
printf("am in mb11d()\n");
      DO(19,I,1,MKK)
      SUM=WKK*BB(I);
/*    CIBMD PREFER SCALAR */
      DO(20,J,1,NKK)
      SUM=SUM+CC(1,J)*AA(I,J);
/*20*/CONTINUE
      SUM=-SUM;
      BB(I)=BB(I)+SUM*WKK;
/*    CIBMD PREFER SCALAR  */
      DO(21,J,1,NKK)
      AA(I,J)=AA(I,J)+SUM*CC(1,J);
/*21*/CONTINUE
/*19*/CONTINUE
      RETURN;
      END
}
/*----------------------------------------------------------------------------
                               MB11E()
-----------------------------------------------------------------------------*/
/*    C@PROCESS DIRECTIVE('IBMD')  */
/*    SUBROUTINE MB11E ( A, B,IA, W,WKK,MKK,NKK)  */
      INT        MB11E (A ,ZA,SA,ZB,SB,IA,W ,WKK,MKK,NKK)
      DOUBLE *A ,*W,WKK;
      INT    IA,MKK,NKK,ZA,SA,ZB,SB;
{
      #define AAA(I,J) (*(A +(ZA-1+I-1)+(SA-1+J-1)*IA))
      #define BBB(I,J) (*(A +(ZB-1+I-1)+(SB-1+J-1)*IA))
      #define WWW(I)   (*(W +(I-1)         ))
/*    REAL             A(IA,NKK),B(IA,NKK),W(*),WKK,SUM  */
      DOUBLE                                        SUM;
      INT I,J;
/*    REAL             ZERO
      DATA ZERO/0.0E0/         */
printf("am in mb11e()\n");
      DO(39,I,1,MKK)
      SUM=ZERO;
/*    CIBMD PREFER SCALAR   */
      DO(40,J,1,NKK)
      SUM=SUM+BBB(1,J)*AAA(I,J);
/*40*/CONTINUE
      SUM=-SUM;
/*    CIBMD PREFER SCALAR  */
      DO(41,J,1,NKK)
      AAA(I,J)=AAA(I,J)+SUM*BBB(1,J);
/*41*/CONTINUE
      WWW(I)=SUM*WKK;
/*38*/CONTINUE
      RETURN;
      END
}
/*----------------------------------------------------------------------------
                               MB11F()
-----------------------------------------------------------------------------*/
/*    @PROCESS DIRECTIVE('IBMD')   */
/*    SUBROUTINE MB11F ( A, B, C,IA,WKK,AKK,MKK,NKK)  */
      INT        MB11F (A ,ZA,SA,ZB,SB,ZC,SC,IA,WKK,AKK,MKK,NKK)
      DOUBLE *A ,WKK,AKK;
      INT    ZA,ZB,ZC,SA,SB,SC;
      INT    IA,MKK,NKK;
{
      #define AAA(I,J) (*(A +(ZA-1+I-1)+(SA-1+J-1)*IA))
      #define BBB(I,J) (*(A +(ZB-1+I-1)+(SB-1+J-1)*IA))
      #define CCC(I)   (*(A +(ZC-1+I-1)+(SC-1    )*IA))
 /*   REAL             A(IA,*),B(IA,*),C(*),SUM,WKK,AKK  */
      DOUBLE                                SUM;
      INT I,J;
printf("am in mb11f()\n");
      DO(42,J,1,NKK)
      SUM=WKK*BBB(1,J);
      DO(43,I,1,MKK)
      SUM=SUM+CCC(I)*AAA(I,J);
/*43*/CONTINUE
      SUM=-SUM;
      BBB(1,J)=SUM*AKK;
/*42*/CONTINUE
      RETURN;
      END
}
/*------------------------------------------------------------------------------
ENDEMODUL    M B 1 1 A   C
------------------------------------------------------------------------------*/
