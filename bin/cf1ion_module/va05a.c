/*     ADD NAME=VA05A   HSL                     SINGLE
*######DATE 14 MAR 1990     COPYRIGHT UKAEA, HARWELL.
C######ALIAS VA05A
C###### CALLS   MB11
*/
 
/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                                   VA05A C
 
 
Minimiering einer N-dimensionalen Funktion
 
-----------------------------------------------------------------------------*/
 
/*----------------------------------------------------------------------------
Includedateien holen
-----------------------------------------------------------------------------*/
#include <stdio.h>          /* damit FILE definiert wird               */
#include <stdlib.h>
#include <math.h>           /* damit sqrt in define_j.c definiert wird */
#define pi (4.0*atan(1.0))  /* atan() braucht <math.h>                 */
#include "types.c"          /* benutze Datentypen laden                */
 
#define SQRT sqrt
#define PRINTF printf
#define IF if
#define FOR for
#define ELSE else
#define GOTO goto
#define RETURN return 0
#define CONTINUE ;}
#define END
 
#define IFI(I,A,B,C) {IF((I)<0)GOTO A;IF((I)==0)GOTO B;IF((I)>0)GOTO C;}
#define IFD(I,A,B,C) {IF((I)<0.)GOTO A;IF((I)==0.)GOTO B;IF((I)>0.)GOTO C;}
#define DO(A,I,IA,IE) FOR( (I)=IA;(I)<=IE;++(I) ){
#define GOTO5(A,B,C,D,M) {IF((M)==1)GOTO A;IF((M)==2)GOTO B;IF((M)==3)GOTO C;IF((M)==4)GOTO D;}
#define GOTO4(A,B,C,M)   {IF((M)==1)GOTO A;IF((M)==2)GOTO B;IF((M)==3)GOTO C;}
#define GOTO3(A,B,M)     {IF((M)==1)GOTO A;IF((M)==2)GOTO B;}
 
#define AMAX1(a,b) ( ((a)>=(b  ))?  (a) : (b) )
#define AMIN1(a,b) ( ((a)<=(b  ))?  (a) : (b) )
/*----------------------------------------------------------------------------
Extern definierte Funktionen
-----------------------------------------------------------------------------*/
 
extern INT CALFUN();   /*definiert in minima.c*/
extern INT MB11A();    /*definiert in mb11a.c */
 
/*    SUBROUTINE VA05A (M,N,F,   X,DSTEP,DMAX,ACC,MAXFUN,IPRINT, W) */
      INT        VA05A (INT *PM,INT *PN,DOUBLE *FF,DOUBLE *XX,DOUBLE *PDSTEP,DOUBLE *PDMAX,DOUBLE *PACC,INT *PMAXFUN,INT *PIPRINT,DOUBLE *WW)
{
      #define X(I) (*(XX+(I)))
      #define F(I) (*(FF+(I)))
      #define W(I) (*(WW+(I)))
 
      DOUBLE DTEST,FMIN,DD,DSS,DM,PARM,DPAR,TINC;
 
      DOUBLE FSQ,PAR,PPAR,DS,DN,SP,DMULT,PRED,AP,AD,PTM;
      DOUBLE ANMULT,DW,FNP,SPP,SS,PJ,ST,DSTEP,DMAX,ACC;
      INT    MAXC,MPN,NT,NTEST,IS,NWI,NWX,NWF,NWC,NWD,NWW,NWV,NWT/*,NWU*/;
      INT    IC,IPC,KS,NTPAR,KK,I,J,K,M,N,MAXFUN,IPRINT;
      CHAR   *TEXT1,*TEXT2,*TEXT3/*,*TEXT4*/,*TEXTX,*TEXTF,*TEXTM;
      M =*PM; N=*PN; MAXFUN=*PMAXFUN; IPRINT=*PIPRINT;
      DSTEP=*PDSTEP; DMAX=*PDMAX; ACC=*PACC;
      TEXTX=" X(%5d) = %17.8f\n";
      TEXTF=" F(%5d) = %17.8f\n";
      TEXTM="     THE SUM OF SQUARES IS %17.8f\n";
 
      NTPAR=1; FNP = 0; PTM = 0; 
 
/*    NOTE THAT THE INSTRUCTION CALLING SUBROUTINE 'MB11A',
C      ON LINE NUMBER '138' ,IS NOT STANDARD FORTRAN
C     DIMENSION F(1),X(1),W(1)
      SET VARIOUS PARAMETERS     */
      MAXC=0;
/*    'MAXC' COUNTS THE NUMBER OF CALLS OF CALFUN   */
      MPN=M+N;
      NT=N+2;
      NTEST=0;
/*    'NT' AND 'NTEST' CAUSE AN ERROR RETURN IF F(X) DOES NOT DECREASE*/
      DTEST=(DOUBLE)(N+N)-(DOUBLE)0.5;
/*    'DTEST' IS USED IN A TEST TO MAINTAIN LINEAR INDEPENDENCE
C     PARTITION THE WORKING SPACE ARRAY W
C     THE FIRST PARTITION HOLDS THE JACOBIAN APPROXIMATION  */
      NWI=M*N;
/*    THE NEXT PARTITION HOLDS THE GENERALIZED INVERSE */
      NWX=NWI+MPN*N;
/*    THE NEXT PARTITION HOLDS THE BEST VECTOR X */
      NWF=NWX+N;
/*    THE NEXT PARTITION HOLDS THE BEST VECTOR F */
      NWC=NWF+M;
/*    THE NEXT PARTITION HOLDS THE COUNTS OF THE INDEPENDENT DIRECTIONS*/
      NWD=NWC+N;
/*    THE NEXT PARTITION HOLDS THE INDEPENDENT DIRECTIONS */
      NWW=NWD+N*N;
/*    THE REMAINDER OF W IS USED FOR SCRATCH VECTORS  */
      NWV=NWW+N;
      NWT=NWV+M;
/*    NWU=NWT+N; */
      FMIN=-1.;
/*    USUALLY 'FMIN' IS THE LEAST CALCULATED VALUE OF F(X) */
      DD=0.;
/*    USUALLY 'DD' IS THE SQUARE OF THE CURRENT STEP LENGTH */
      DSS=DSTEP*DSTEP;
      DM=DMAX*DMAX;
      PARM=SQRT(ACC)/DMAX;
/*    'PARM' IS THE LEAST VALUE OF THE MARQUARDT PARAMETER */
      DPAR=10.*DM;
      PAR = PPAR = PARM; PPAR *= PAR;
/*    'DPAR' AND 'NTPAR' ARE USED TO REGULATE THE MARQUARDT PARAMETER*/
      IS=4;
/*    'IS' CONTROLS A GOTO  STATEMENT FOLLOWING A CALL OF CALFUN */
      IC=0;
      TINC=1.;
/*    'TINC' IS USED IN THE CRITERION TO INCREASE THE STEP LENGTH
C     START A NEW PAGE FOR PRINTING  */
      IFI(IPRINT,C1,C3,C1);
/*  1 WRITE(6,2)
    2 FORMAT (1H1,4X,'THE FOLLOWING OUTPUT IS PROVIDED BY SUBROUTINE',
     1' VA05A'//) */
   C1:TEXT1="THE FOLLOWING OUTPUT IS PROVIDED BY SUBROUTINE VA05A";
      PRINTF("1    %s\n",TEXT1);
      IPC=0;
      GOTO C3;
/*    TEST WHETHER THERE HAVE BEEN MAXFUN CALLS OF CALFUN */
   C4:IFI(MAXFUN-MAXC,C5,C5,C3);
   C5:IFI(IPRINT,C139,C140,C139);
 C140:IPRINT=2;
      GOTO C19;
/*139 WRITE(6,6)MAXC
    6 FORMAT (///5X,'ERROR RETURN FROM VA05A BECAUSE THERE HAVE BEEN',
     1I5,' CALLS OF CALFUN') */
 C139:TEXT1="ERROR RETURN FROM VA05A BECAUSE THERE HAVE BEEN";
      TEXT2="CALLS OF CALFUN";
      PRINTF("%s%d%s\n",TEXT1,MAXC,TEXT2);
      GOTO C7;
/*    CALL THE SUBROUTINE CALFUN  */
   C3:MAXC=MAXC+1;
      CALFUN (M,N,FF,XX);
/*    CALCULATE THE SUM OF SQUARES */
      FSQ=0.;
      DO(8,I,1,M)
      FSQ=FSQ+F(I)*F(I);
/* 8*/CONTINUE
/*    TEST FOR ERROR RETURN BECAUSE F(X) DOES NOT DECREASE */
      GOTO5(C9,C10,C9,C10,IS);
   C9:IFD(FSQ-FMIN,C11,C12,C12);
  C12:IFD(DD-DSS,C13,C13,C10);
  C13:NTEST=NTEST-1;
      IFI(NTEST,C14,C14,C10);
  C14:IFI(IPRINT,C15,C17,C15);
  C17:IPRINT=1;
      GOTO C19;
/* 15 WRITE(6,16)
   16 FORMAT (///5X,'ERROR RETURN FROM VA05A BECAUSE F(X) NO LONGER',
     1' DECREASES'//5X,'THIS MAY BE DUE TO THE VALUES OF DSTEP',
     2' AND ACC, OR TO LOSS OF RANK IN THE JACOBIAN MATRIX')
C     PROVIDE PRINTING OF FINAL SOLUTION IF REQUESTED  */
C15:  TEXT1="      ERROR  RETURN FROM  VA05A  BECAUSE  F(X) NO LONGER";
      TEXT2="      DECREASES. THIS MAY BE DUE  TO THE VALUES OF DSTEP";
      TEXT3="      AND ACC, OR TO LOSS OF RANK IN THE JACOBIAN MATRIX.";
      PRINTF("\n%s\n%s\n%s\n",TEXT1,TEXT2,TEXT3);
   C7:IFI(IPRINT,C18,C19,C18);
/* 18 WRITE(6,20)MAXC
   20 FORMAT (///5X,'THE FINAL SOLUTION CALCULATED BY VA05A REQUIRED',
     1I5,' CALLS OF CALFUN, AND IS') */
  C18:
      TEXT1="      THE FINAL SOLUTION CALCULATED BY VA05A REQUIRED";
      TEXT2="      CALLS OF CALFUN, AND IS";
      PRINTF("%s %d %s\n",TEXT1,MAXC,TEXT2);
/*    WRITE(6,21)(I,W(NWX+I),I=1,N)
   21 FORMAT (//4X,'I',7X,'X(I)',10X,'I',7X,'X(I)',10X,'I',7X,'X(I)',
     110X,'I',7X,'X(I)',10X,'I',7X,'X(I)'//5(I5,E17.8))*/
      FOR(I=1;I<=N;++I)PRINTF(TEXTX,I,W(NWX+I));
 
/*    WRITE(6,22)(I,W(NWF+I),I=1,M)
   22 FORMAT (//4X,'I',7X,'F(I)',10X,'I',7X,'F(I)',10X,'I',7X,'F(I)',
     110X,'I',7X,'F(I)',10X,'I',7X,'F(I)'//5(I5,E17.8))*/
      FOR(I=1;I<=M;++I)PRINTF(TEXTF,I,W(NWF+I));
 
/*    WRITE(6,23)FMIN
   23 FORMAT (/5X,'THE SUM OF SQUARES IS',E17.8) */
      PRINTF(TEXTM,FMIN);
 
/*    RESTORE THE BEST VALUES OF X AND F */
  C19:DO(135,I,1,N)
      X(I)=W(NWX+I);
/*135*/CONTINUE
      DO(136,I,1,M)
      F(I)=W(NWF+I);
/*136*/CONTINUE
      RETURN;
  C11:NTEST=NT;
/*    PROVIDE ORDINARY PRINTING IF REQUESTED */
  C10:IFI(ABS(IPRINT)-1,C39,C38,C40);
/* 38 WRITE(6,41)MAXC
   41 FORMAT (///5X,'AT THE',I5,' TH CALL OF CALFUN WE HAVE') */
  C38:TEXT1="     AT THE";
      TEXT2="      TH CALL OF CALFUN WE HAVE";
      PRINTF("%s%5d%s\n",TEXT1,MAXC,TEXT2);
 
/* 42 WRITE(6,21)(I,X(I),I=1,N)   */
  C42:FOR(I=1;I<=N;++I)PRINTF(TEXTX,I,X(I));
/*    WRITE(6,23)FSQ      */
      PRINTF(TEXTM,FSQ );
 
      IFI(IPRINT,C39,C39,C142);
/*142 WRITE(6,22) (I,F(I),I=1,M)  */
 C142:FOR(I=1;I<=M;++I)PRINTF(TEXTF,I,F(I));
 
      GOTO  C39;
  C40:IPC=IPC-1;
      IFI(IPC,C43,C43,C39);
/* 43 WRITE(6,44)MAXC
   44 FORMAT (///5X,'THE BEST ESTIMATE AFTER',I5,' CALLS OF CALFUN IS')*/
 C43: TEXT1="     THE BEST ESTIMATE AFTER";
      TEXT2="     CALLS OF CALFUN IS";
      PRINTF("%s%5d%s\n",TEXT1,MAXC,TEXT2);
      IPC=ABS(IPRINT);
      IFD(FSQ-FMIN,C42,C45,C45);
  C45:IFD(FMIN,C42,C46,C46);
/* 46 WRITE(6,21)(I,W(NWX+I),I=1,N)  */
  C46:FOR(I=1;I<=N;++I)PRINTF(TEXTX,I,W(NWX+I));
/*    WRITE(6,23)FMIN      */
      PRINTF(TEXTM,FMIN);
      IFI(IPRINT,C39,C39,C143);
/*143 WRITE(6,22) (I,W(NWF+I),I=1,M)  */
 C143:FOR(I=1;I<=M;++I)PRINTF(TEXTF,I,W(NWF+I));
  C39:GOTO5(C49,C47,C47,C48,IS);
/*    STORE THE INITIAL VECTORS X AND F */
  C48:IFI(IC,C50,C50,C51);
  C50:DO(52,I,1,N)
      W(NWX+I)=X(I);
/*52*/CONTINUE
      GOTO C54;
/*    CALCULATE THE INITIAL JACOBIAN APPROXIMATION */
  C51:K=IC;
      DO(55,I,1,M)
      W(K)=(F(I)-W(NWF+I))/DSTEP;
      K=K+N;
/*55*/CONTINUE
/*    TEST WHETHER THE MOST RECENT X IS BEST  */
      IFD(FMIN-FSQ,C56,C56,C57);
  C56:X(IC)=W(NWX+IC);
      GOTO C58;
  C57:W(NWX+IC)=X(IC);
  C54:DO(53,I,1,M)
      W(NWF+I)=F(I);
/*53*/CONTINUE
      FMIN=FSQ;
/*    SET X FOR THE NEXT CALL OF CALFUN */
  C58:IC=IC+1;
      IFI(IC-N,C59,C59,C60);
  C59:X(IC)=W(NWX+IC)+DSTEP;
      GOTO C3;
/*    SET THE DIRECTION MATRIX */
  C60:K=NWD;
      DO(61,I,1,N)
      DO(62,J,1,N)
      K=K+1;
      W(K)=0.;
/*62*/CONTINUE
      W(K+I-N)=1.;
      W(NWC+I)=1.+(DOUBLE)(N-I);
/*61*/CONTINUE
/*    SET THE MARQUARDT PARAMETER TO ITS LEAST VALUE */
  C24:PAR=PARM;
/*    COPY THE JACOBIAN AND APPEND THE MARQUARDT MATRIX */
  C25:PPAR=PAR*PAR;
      NTPAR=0;
  C63:KK=0;
      K=NWI+NWI;
      DO(26,I,1,N)
      DO(141,J,1,M)
      KK=KK+1;
      W(KK+NWI)=W(KK);
/*141*/CONTINUE
      DO(27,J,1,N)
      K=K+1;
      W(K)=0.;
/*27*/CONTINUE
      W(K+I-N)=PAR;
/*26*/CONTINUE
/*    CALCULATE THE GENERALIZED INVERSE OF J */
      MB11A (N,MPN,WW,NWI+1,N, NWW+1 );
/*    NOTE THAT THE THIRD AND FIFTH ENTRIES OF THIS ARGUMENT LIST
C     STAND FOR ONE-DIMENSIONAL ARRAYS.
C     START THE ITERATION BY TESTING FMIN  */
  C64:IFI(FMIN-ACC,C7,C7,C65);
/*    NEXT PREDICT THE DESCENT AND MARQUARDT MINIMA  */
  C65:DS=0.;
      DN=0.;
      SP=0.;
      DO(66,I,1,N)
      X(I)=0.;
      F(I)=0.;
      K=I;
      DO(67,J,1,M)
      X(I)=X(I)-W(K)*W(NWF+J);
      F(I)=F(I)-W(NWI+K)*W(NWF+J);
      K=K+N;
/*67*/CONTINUE
      DS=DS+X(I)*X(I);
      DN=DN+F(I)*F(I);
      SP=SP+X(I)*F(I);
/*66*/CONTINUE
/*    PREDICT THE REDUCTION IN F(X) DUE TO THE MARQUARDT STEP
C     AND ALSO PREDICT THE LENGTH OF THE STEEPEST DESCENT STEP */
      PRED=SP+SP;
      DMULT=0.;
      K=0;
      DO(68,I,1,M)
      AP=0.;
      AD=0.;
      DO(69,J,1,N)
      K=K+1;
      AP=AP+W(K)*F(J);
      AD=AD+W(K)*X(J);
/*69*/CONTINUE
      PRED=PRED-AP*AP;
      DMULT=DMULT+AD*AD;
/*68*/CONTINUE
/*    TEST FOR CONVERGENCE   */
      IFD(DN-DM,C28,C28,C29);
  C28:AP=SQRT(DN);
      IFD(PRED+2.*PPAR*AP*(DMAX-AP)-ACC,C7,C7,C70);
  C29:IFD(PRED+PPAR*(DM-DN)-ACC,C7,C7,C70);
/*    TEST WHETHER TO APPLY THE FULL MARQUARDT CORRECTION */
  C70:DMULT=DS/DMULT;
      DS=DS*DMULT*DMULT;
  C71:IS=2;
      IFD(DN-DD,C72,C72,C73);
/*    TEST THAT THE MARQUARDT PARAMETER HAS ITS LEAST VALUE */
  C72:IFD(PAR-PARM,C30,C30,C24);
  C30:DD=AMAX1(DN,DSS);
      DS=0.25*DN;
      TINC=1.;
      IFD(DN-DSS,C74,C132,C132);
  C74:IS=3;
      GOTO C103;
/*    TEST WHETHER TO INCREASE THE MARQUARDT PARAMETER */
  C73:IFD(DN-DPAR,C31,C31,C32);
  C31:NTPAR=0;
      GOTO C33;
  C32:IFI(NTPAR,C34,C34,C35);
  C34:NTPAR=1;
      PTM=DN;
      GOTO C33;
  C35:NTPAR=NTPAR+1;
      PTM=AMIN1(PTM,DN);
      IFI(NTPAR-NT,C33,C36,C36);
/*    SET THE LARGER VALUE OF THE MARQUARDT PARAMETER */
/*C36:PAR=PAR*(PTM/DM)**0.25;*/
  C36:PAR=PAR*SQRT(SQRT((PTM/DM)));
      IFD(6.*DD-DM,C137,C25,C25);
 C137:AP=SQRT(PRED/DN);
      IFD(AP-PAR,C25,C25,C138);
/*138:PAR=AMIN1(AP,PAR*(DM/(6.*DD))**0.25);*/
 C138:PAR=AMIN1(AP,PAR* SQRT(SQRT(DD/(6.*DD))) );
      GOTO C25;
/*    TEST WHETHER TO USE THE STEEPEST DESCENT DIRECTION   */
  C33:IFD(DS-DD,C75,C76,C76);
/*    TEST WHETHER THE INITIAL VALUE OF DD HAS BEEN SET */
  C76:IFD(DD,C77,C77,C78);
  C77:DD=AMIN1(DM,DS);
      IFD(DD-DSS,C79,C78,C78);
  C79:DD=DSS;
      GOTO C71;
/*    SET THE MULTIPLIER OF THE STEEPEST DESCENT DIRECTION */
  C78:ANMULT=0.;
      DMULT=DMULT*SQRT(DD/DS);
      GOTO C80;
/*    INTERPOLATE BETWEEN THE STEEPEST DESCENT AND MARQUARDT DIRECTIONS*/
  C75:SP=SP*DMULT;
/*    ANMULT=(DD-DS)/((SP-DS)+SQRT((SP-DD)**2+(DN-DD)*(DD-DS)))*/
      ANMULT=(DD-DS)/((SP-DS)+SQRT((SP-DD)*(SP-DD)+(DN-DD)*(DD-DS)));
      DMULT=DMULT*(1.-ANMULT);
/*    CALCULATE THE CORRECTION TO X, AND ITS ANGLE WITH THE FIRST
C     DIRECTION */
  C80:DN=0.;
      SP=0.;
      DO(81,I,1,N);
      F(I)=DMULT*X(I)+ANMULT*F(I);
      DN=DN+F(I)*F(I);
      SP=SP+F(I)*W(NWD+I);
/*81*/CONTINUE
      DS=0.25*DN;
/*    TEST WHETHER AN EXTRA STEP IS NEEDED FOR INDEPENDENCE */
      IFD(W(NWC+1)-DTEST,C132,C132,C82);
  C82:IFD(SP*SP-DS,C83,C132,C132);
/*    TAKE THE EXTRA STEP AND UPDATE THE DIRECTION MATRIX */
  C83:DO(84,I,1,N)
      X(I)=W(NWX+I)+DSTEP*W(NWD+I);
      W(NWC+I)=W(NWC+I+1)+1.;
/*84*/CONTINUE
      W(NWD)=1.;
      IF(N <= 1)GOTO C4;
      DO(85,I,1,N)
      K=NWD+I;
      SP=W(K);
      DO(86,J,2,N)
      W(K)=W(K+N);
      K=K+N;
/*86*/CONTINUE
      W(K)=SP;
/*85*/CONTINUE;
      GOTO C4;
/*    EXPRESS THE NEW DIRECTION IN TERMS OF THOSE OF THE DIRECTION
C     MATRIX, AND UPDATE THE COUNTS IN W(NWC+1) ETC. */
 C132:IF(N >= 2)GOTO C153;
      IS=1;
      GOTO C152;
 C153:SP=0.;
      K=NWD;
      DW=0.0;
      DO(87,I,1,N)
      X(I)=DW;
      DW=0.;
      DO(88,J,1,N)
      K=K+1;
      DW=DW+F(J)*W(K);
/*88*/CONTINUE
      GOTO3(C89,C90,IS);
  C90:W(NWC+I)=W(NWC+I)+1.;
      SP=SP+DW*DW;
      IFD(SP-DS,C87,C87,C91);
  C91:IS=1;
      KK=I;
      X(1)=DW;
      GOTO C92;
  C89:X(I)=DW;
  C92:W(NWC+I)=W(NWC+I+1)+1.;
  C87:CONTINUE
      W(NWD)=1.;
/*    REORDER THE DIRECTIONS SO THAT KK IS FIRST   */
      IFI(KK-1,C93,C93,C94);
  C94:KS=NWC+KK*N;
      DO(95,I,1,N)
      K=KS+I;
      SP=W(K);
      DO(96,J,2,KK)
      W(K)=W(K-N);
      K=K-N;
/*96*/CONTINUE
      W(K)=SP;
/*95*/CONTINUE
/*    GENERATE THE NEW ORTHOGONAL DIRECTION MATRIX */
  C93:DO(97,I,1,N)
      W(NWW+I)=0.;
/*97*/CONTINUE
      SP=X(1)*X(1);
      K=NWD;
      DO(98,I,2,N)
      DS=SQRT(SP*(SP+X(I)*X(I)));
      DW=SP/DS;
      DS=X(I)/DS;
      SP=SP+X(I)*X(I);
      DO(99,J,1,N)
      K=K+1;
      W(NWW+J)=W(NWW+J)+X(I-1)*W(K);
      W(K)=DW*W(K+N)-DS*W(NWW+J);
/*99*/CONTINUE
/*98*/CONTINUE
      SP=1./SQRT(DN);
      DO(100,I,1,N)
      K=K+1;
      W(K)=SP*F(I);
/*100*/CONTINUE
/*    PREDICT THE NEW RIGHT HAND SIDES   */
 C152:FNP=0.;
      K=0;
      DO(101,I,1,M)
      W(NWW+I)=W(NWF+I);
      DO(102,J,1,N)
      K=K+1;
      W(NWW+I)=W(NWW+I)+W(K)*F(J);
/*102*/CONTINUE
/*    FNP=FNP+W(NWW+I)**2 */
      FNP=FNP+W(NWW+I)*W(NWW+I);
/*101*/CONTINUE
/*    CALCULATE THE NEXT VECTOR X, AND THEN CALL CALFUN  */
 C103:DO(104,I,1,N)
      X(I)=W(NWX+I)+F(I);
/*104*/CONTINUE
      GOTO C4;
/*    UPDATE THE STEP SIZE  */
  C49:DMULT=0.9*FMIN+0.1*FNP-FSQ;
      IFD(DMULT,C105,C108,C108);
 C105:DD=AMAX1(DSS,0.25*DD);
      TINC=1.;
      IFD(FSQ-FMIN,C106,C107,C107);
/*    TRY THE TEST TO DECIDE WHETHER TO INCREASE THE STEP LENGTH */
 C108:SP=0.;
      SS=0.;
      DO(109,I,1,M)
      SP=SP+ABSD(F(I)*(F(I)-W(NWW+I)));
/*    SS=SS+(F(I)-W(NWW+I))**2   */
      SS=SS+(F(I)-W(NWW+I))*W(NWW+I);
/*109*/CONTINUE
      PJ=1.+DMULT/(SP+SQRT(SP*SP+DMULT*SS));
      SP=AMIN1(4.,AMIN1(TINC,PJ));
      TINC=PJ/SP;
      DD=AMIN1(DM,SP*DD);
      GOTO C106;
/*    IF F(X) IMPROVES STORE THE NEW VALUE OF X */
  C47:IFD(FSQ-FMIN,C106,C110,C110);
 C106:FMIN=FSQ;
      DO(111,I,1,N)
      SP=X(I);
      X(I)=W(NWX+I);
      W(NWX+I)=SP;
/*111*/CONTINUE
      DO(112,I,1,M)
      SP=F(I);
      F(I)=W(NWF+I);
      W(NWF+I)=SP;
/*112*/CONTINUE
 C110:GOTO4(C107,C107,C113,IS);
 C113:IS=2;
      IFD(FMIN-ACC,C7,C7,C83);
/*    CALCULATE THE CHANGES IN X AND IN F */
 C107:DS=0.;
      DO(114,I,1,N)
      X(I)=X(I)-W(NWX+I);
      DS=DS+X(I)*X(I);
/*114*/CONTINUE
      DO(115,I,1,M)
      F(I)=F(I)-W(NWF+I);
/*115*/CONTINUE
/*    CALCULATE THE GENERALIZED INVERSE TIMES THE CHANGE IN X */
      K=NWI;
      SS=0.;
      DO(116,I,1,MPN)
      SP=0.;
      DO(117,J,1,N)
      K=K+1;
      SP=SP+W(K)*X(J);
/*117*/CONTINUE
      W(NWV+I)=SP;
      SS=SS+SP*SP;
/*116*/CONTINUE
/*    CALCULATE J TIMES THE CHANGE IN F
C     ALSO APPLY PROJECTION TO THE GENERALIZED INVERSE */
      DO(118,I,1,N)
      ST=0.;
      K=NWI+I;
      DO(119,J,1,MPN)
      ST=ST+W(K)*W(J+NWV);
      K=K+N;
/*119*/CONTINUE
      ST=ST/SS;
      K=NWI+I;
      DO(120,J,1,MPN)
      W(K)=W(K)-ST*W(J+NWV);
      K=K+N;
/*120*/CONTINUE
      ST=PPAR*X(I);
      K=I;
      DO(121,J,1,M)
      ST=ST+W(K)*F(J);
      K=K+N;
/*121*/CONTINUE
      W(NWW+I)=ST;
/*118*/CONTINUE
/*    REVISE J AND CALCULATE ROW VECTOR FOR CORRECTION TO INVERSE */
      IC=0;
      K=0;
      KK=NWI;
      SP=0.;
      SPP=0.;
      DO(122,I,1,M)
      SS=F(I);
      ST=F(I);
      DO(123,J,1,N)
      IC=IC+1;
      KK=KK+1;
      SS=SS-W(IC)*X(J);
      ST=ST-W(KK)*W(NWW+J);
/*123*/CONTINUE
      SS=SS/DS;
      W(NWV+I)=ST;
      SP=SP+F(I)*ST;
      SPP=SPP+ST*ST;
      DO(124,J,1,N)
      K=K+1;
      W(K)=W(K)+SS*X(J);
/*124*/CONTINUE
/*122*/CONTINUE
      DO(125,I,1,N)
      ST=PAR*X(I);
      DO(126,J,1,N)
      KK=KK+1;
      ST=ST-W(KK)*W(NWW+J);
/*126*/CONTINUE
      W(NWT+I)=ST;
      SP=SP+PAR*X(I)*ST;
      SPP=SPP+ST*ST;
/*125*/CONTINUE
/*    TEST THAT THE SCALAR PRODUCT IS SUFFICIENTLY ACCURATE */
      IFD(0.01*SPP-ABSD(SP-SPP),C63,C63,C127);
/*    CALCULATE THE NEW GENERALIZED INVERSE   */
 C127:DO(128,I,1,N)
      K=NWI+I;
      ST=X(I);
      DO(129,J,1,M)
      ST=ST-W(K)*F(J);
      K=K+N;
/*129*/CONTINUE
      SS=0.;
      DO(130,J,1,N)
      SS=SS+W(K)*X(J);
      K=K+N;
/*130*/CONTINUE
      ST=(ST-PAR*SS)/SP;
      K=NWI+I;
      DO(131,J,1,MPN)
      W(K)=W(K)+ST*W(NWV+J);
      K=K+N;
/*131*/CONTINUE
/*128*/CONTINUE
      GOTO C64;
      END
}
/*------------------------------------------------------------------------------
ENDEMODUL    V A 0 5 A   C
------------------------------------------------------------------------------*/
