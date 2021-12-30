      PROGRAM RCN2
C     PROGRAM RCN2(IN2,TAPE10=IN2,OUT2,TAPE9=OUT2,
C    1  TAPE2N,TAPE2=TAPE2N,OUT2ING,TAPE11=OUT2ING,TAPE31,TAPE32)
C
C     THIS PROGRAM COPYRIGHT 1980 BY LOS ALAMOS NATIONAL LABORATORY
C          AND THE U.S. GOVERNMENT.  THE PROGRAM MAY BE FREELY
C          DISTRIBUTED GRATIS ANYWHERE IN THE WORLD, BUT MAY NOT BE SOLD.
C
C
C          THIS PROGRAM READS RADIAL WAVEFUNCTIONS FROM TAPE (OR DISK)
C               UNIT NUMBER 2 PREPARED BY PROGRAM RCN OR HF, AND (USING
C               INPUT OPTION G5INP) CALCULATES CONFIGURATION-INTERACTION
C               PARAMETERS (RK) AND RADIAL ELECTRIC MULTIPOLE INTEGRALS,
C               SCALES THE ENERGY-LEVEL PARAMETER VALUES (FK, GK, ZETA,
C               AND RK), AND PUNCHES COMPLETE INPUT CARD DECKS FOR
C               PROGRAM RCG MOD 5, 6, 7, 8, OR 9.
C          OPTION TO INCLUDE OVERLAP INTEGRALS OF SPECTATOR ELECTRONS
C               IN COMPUTED CONFIGURATION-INTERACTION AND ELECTRIC-
C               MULTIPOLE INTEGRALS ADDED AT THE ZEEMAN LABORATORY,
C               DECEMBER 1987.
C
C          THIS VERSION, EASILY ADAPTABLE TO THE VAX 8200,
C               MADE AT LUND UNIVERSITY MARCH 1988,
C                    MODIFIED AT LANL AUGUST-DECEMBER 1988.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KBGKS=41,KCASE=50,KM=1801,KO=20)
      PARAMETER (KCO=100,KP=KCO+1)
      DIMENSION PA(KM),PB(KM),PC(KM),PD(KM),CALSQ(5)
C
      COMMON NLZ(4),SUBR0,SUBR,NCFG(2),MMIN,NCDES,NKK(4),L(4),
     1  M,IEXCH,IZ1,IDUM0,OVERLP,X(KM),DUMMY(KM,7)
C
      COMMON/C1/FACT(5),IFACT(5),IDUMMY2,DUMMY2(20)
      COMMON/C3/ICON,ISLI,IDIP,IC,ICN(2),NCONFN(2),NC0,ICASE,
     1  NSCONF(3,2),DEL(KCO),IDIPP,IOVFACT,OVFACT
      COMMON/CP/P(KM,4),WNL(2),RM1(KM),RMKP1(KM)
      COMMON/LCM1/NT2(KO),LT2(KO),WT2(KO),NCSPVS,IDUM2,VPAR(50),
     1  NCONF,IZ,Z,ION,IREL,HXID,C,MESH,IDB,EXF,CORRF,KUT,NPAR
      COMMON/LCM1A/NLBCD
      COMMON/LCM1B/NNLZ(KO),NKKK(KO),EE(KO),EE8,R(KM),RU(KM),
     1  PNL(KM,KO),NPNL
      COMMON/GENOSC/IGEN,NBIGKS,RJ,BIGKS(KBGKS),IBK
      COMMON/TIMING/TINIT,ITINIT
      CHARACTER DH*80,OPTION*28,NLH(4)*3,HA*3,HB*3,HC*3,HD*3,HA1*3,
     1  HB1*3,HE*3,NLBCD(KO)*8,NLHT(8)*3,NLB(8,KCO)*3,KST(KCASE,2)*3
      CHARACTER*6 CONF(3),CST(KCASE,6),CONF1(3),CONF2(3),
     1  CF1(KP),CF2(KP),CF3(KP)
      COMMON/CDH/DH,OPTION,NLH,HA1,HB1,HE,NLHT,NLB,KST,
     1  CONF,CST,CONF1,CONF2,CF1,CF2,CF3
      EQUIVALENCE (P,PA),(P(1,2),PB),(P(1,3),PC),(P(1,4),PD),
     1  (L(1),LA),(L(2),LB),(L(3),LC),(L(4),LD),
     2  (NLH(1),HA),(NLH(2),HB),(NLH(3),HC),(NLH(4),HD)  
      CHARACTER*20 INF
C
      DATA (CALSQ(I),I=1,5)/3HDIP,3HOVE,3HSLI,3HZET,3HG5I/
C
      ITINIT=0                   ! INITIALIZE TIMING FOR IBM-PC
C
C          OPEN FILES
C               (SET III=0 TO NAME FILES INTERACTIVELY, FOR VAX)
C               (SET III=1 TO USE STANDARD FILE NAMES, FOR VAX)
C               (SET III>1 AND UNCOMMENT PROGRAM CARD AT
C                     START OF PROGRAM FOR CRAY AND CDC COMPUTERS)
C
      III=1
      IF (III.EQ.1) GO TO 50
      IF (III.GT.1) GO TO 100
C
      WRITE(6,*) '  INPUT FILE NAME  '
      READ(5,11) INF
   11 FORMAT(A)
      OPEN(10,FILE=INF,STATUS='OLD')
      WRITE(6,*) '  OUTPUT FILE NAME  '
      READ(5,11) INF
      OPEN(9,FILE=INF,STATUS='UNKNOWN')
      WRITE(6,*) '  GIVE FILE NAME FOR WAVEFUNC. FILE '
      READ(5,11) INF
      OPEN(2,FILE=INF,STATUS='OLD',FORM='UNFORMATTED')
      WRITE(6,*) '  GIVE FILE NAME FOR RCG INPUT FILE '
      READ(5,11) INF
      OPEN(11,FILE=INF,STATUS='UNKNOWN')
      GO TO 55
C
C          USE THESE OPEN STATEMENTS FOR STANDARD FILE NAMES
C               (SET III=1 ON STATEMENT ABOVE)
C
   50 OPEN(10,FILE='IN2',STATUS='OLD')
      OPEN(9,FILE='OUT2',STATUS='UNKNOWN')
      OPEN(2,FILE='TAPE2N',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(11,FILE='ING11',STATUS='UNKNOWN')
   55 OPEN(31,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(32,STATUS='SCRATCH',FORM='UNFORMATTED')
C
  100 REWIND 2
      IC=2
      ICN(1)=2
      ICN(2)=2
      NCONF=99
      NCONFN(1)=99
      NCONFN(2)=99
      SUBR0=0.D0
      P(1,1)=-99.D0
C
  200 READ (10,19) DH
   19 FORMAT (A)
C     DECODE (3,20,DH) SUBR
      READ (DH,20) SUBR
   20 FORMAT (A3)
      IF (SUBR.EQ.CALSQ(5)) GO TO 800
C     DECODE (45,21,DH) SUBR,IZ1,NCFG(1),NCFG(2),HA,HB,HC,HD,MMIN
      READ (DH,21) SUBR,IZ1,NCFG(1),NCFG(2),HA,HB,HC,HD,MMIN
   21 FORMAT (A3,4X,I3,2I5, 4(2X,A3), I5)
      IF (IZ1.LT.0) WRITE (11,22)
   22 FORMAT (5H   -1,74X,1H )
      IF (IZ1.LT.0) THEN 
        CLOSE (31,STATUS='DELETE')
        CLOSE (32,STATUS='DELETE')
        STOP '(NORMAL EXIT)'
      ENDIF
      DO 300 I=1,5
      IF(SUBR.EQ.CALSQ(I)) GO TO (500,550,600,700,800) I
  300 CONTINUE
      GO TO 200
  500 CALL DIP                                                          DIPOLE
      GO TO 200
  550 CALL OVER                                                         OVERLAP
      GO TO 200
  600 CALL SLI2                                                         SLATER I
      GO TO 200
  700 CALL ZETA2                                                        ZETA
      GO TO 200
  800 CALL G5INP                                                        G5INP
      WRITE (9,80)
   80 FORMAT (//15H0FINISHED G5INP/)
      GO TO 200
      END
      SUBROUTINE OVER
C
C          EVALUATE OVERLAP INTEGRALS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KBGKS=41,KCASE=50,KM=1801,KO=20)
      PARAMETER (KCO=100,KP=KCO+1)
      DIMENSION PA(KM),PB(KM),PC(KM),PD(KM)
C
      COMMON NLZ(4),SUBR0,SUBR,NCFG(2),MMIN,NCDES,NKK(4),L(4),
     1  M,IEXCH,IZ1,IDUM0,OVERLP,X(KM),ABSOVR,FRAC
C
      COMMON/C3/ICON,ISLI,IDIP,IC,ICN(2),NCONFN(2),NC0,ICASE,
     1  NSCONF(3,2),DEL(KCO),IDIPP,IOVFACT,OVFACT
      COMMON/CP/P(KM,4),WNL(2),RM1(KM),RMKP1(KM)
      COMMON/LCM1/NT2(KO),LT2(KO),WT2(KO),NCSPVS,IDUM2,VPAR(50),
     1  NCONF,IZ,Z,ION,IREL,HXID,C,MESH,IDB,EXF,CORRF,KUT,NPAR
      COMMON/LCM1A/NLBCD
      COMMON/LCM1B/NNLZ(KO),NKKK(KO),EE(KO),EE8,R(KM),RU(KM),
     1  PNL(KM,KO),NPNL
      CHARACTER DH*80,OPTION*28,NLH(4)*3,HA*3,HB*3,HC*3,HD*3,HA1*3,
     1  HB1*3,HE*3,NLBCD(KO)*8,NLHT(8)*3,NLB(8,KCO)*3,KST(KCASE,2)*3
      CHARACTER*6 CONF(3),CST(KCASE,6),CONF1(3),CONF2(3),
     1  CF1(KP),CF2(KP),CF3(KP)
      COMMON/CDH/DH,OPTION,NLH,HA1,HB1,HE,NLHT,NLB,KST,
     1  CONF,CST,CONF1,CONF2,CF1,CF2,CF3
      EQUIVALENCE (P,PA),(P(1,2),PB),(P(1,3),PC),(P(1,4),PD),
     1  (L(1),LA),(L(2),LB),(L(3),LC),(L(4),LD),
     2  (NLH(1),HA),(NLH(2),HB),(NLH(3),HC),(NLH(4),HD)
C
      K=1
      IF (SUBR) 100,105,100
  100 IF (SUBR.EQ.SUBR0) GO TO 150
  101 WRITE (9,10)
   10 FORMAT (1H1)
  105 SUBR0=SUBR
      IF (ISLI.GT.0) GO TO 200
      WRITE (9,11)
   11 FORMAT (1H0//15H  Z  KUT       ,9X,5HCONF1,17X,5HCONF2,15X,
     1  16HOVERLAP INTEGRAL,11X,8HFRACTION/17H --- ---         ,
     2  22H-------------------   ,23H-------------------    ,
     3  24H------------------------,7X,8H--------/)
      GO TO 200
C
  150 IF (NCFG(1).NE.NCONFN(K).OR.HA.NE.HA1) GO TO 200
      IF (IOVFACT.GT.2) GO TO 200
      N=2
      GO TO 251
C
C          READ WAVEFUNCTIONS FROM TAPE IC AND SELECT THOSE DESIRED
C
  200 N=1
      IZ0=IZ1
      HA1=HA
      IF (ISLI.GT.0) GO TO 205
      WRITE (9,20)
   20 FORMAT (1H )
  205 IC=ICN(N)
      IF (IC.GT.30) K=IC-30
      IF (NCFG(N).EQ.NCONF) GO TO 220
      IF (NCFG(N)-NCONFN(K)) 210,210,211
C     IF (NCFG(N)-NCONFN(K)) 210,206,211
C 206 BACKSPACE IC
C     GO TO 211
  210 REWIND IC
  211 READ (IC) NT2,LT2,WT2,NCSPVS,VPAR,CONF,NCONF,IZ,Z,
     1  ION,IREL,HXID,MESH,C,IDB,EXF,CORRF,KUT,NPAR,
     2  NLBCD,NNLZ,NKKK,EE,EE8,NPNL,R,RU, ((PNL(I,M),I=1,KM),M=1,NPNL)
      IF (NCONF.NE.NCFG(N)) GO TO 211
      NCONFN(K)=NCONF
C
  220 IF (N.GT.1) GO TO 223
  221 DO 222 I=1,3
  222 CONF1(I)=CONF(I)
      GO TO 225
  223 DO 224 I=1,3
  224 CONF2(I)=CONF(I)
  225 DO 226 M=1,NCSPVS
      IF (NLBCD(M)(1:3).EQ.NLH(N)) GO TO 230
  226 CONTINUE
      GO TO 1000
C
  230 NKK(N)=NKKK(M)
      NLZ(N)=NNLZ(M)
      WNL(N)=WT2(M)
      L(N)=NNLZ(M)-100*(NNLZ(M)/100)
      J=M+NPNL-NCSPVS
      IF (J.LE.0) GO TO 1000
      IF (NLH(1).EQ.NLH(2).AND.NNLZ(M).GE.9900) GO TO 260
      DO 231 I=1,MESH
  231 P(I,N)=PNL(I,J)
C
  250 N=N+1
  251 IF (N-2) 200,205,300
C          SET OVERLP=1.0 FOR CONTINUUM ELECTRONS
  260 OVERLP=1.D0
      ABSOVR=1.D0
      GO TO 400
C
C          BEGIN INTEGRATION LOOP
C
  300 NF=-1
      OVERLP=0.D0
      ABSOVR=0.D0
      IMAX=MIN0(NKK(1),NKK(2))
      DO 399 I=1,IMAX
      IF (NF) 340,350,360
  340 ERAS=1.5D0
      IF (I-1) 342,341,342
  341 ERAS=1.D0
      GO TO 344
  342 IF (I-IMAX) 344,343,344
  343 ERAS=0.5D0
  344 ERAS=PA(I)*PB(I)*ERAS
      OVERLP=(OVERLP+ERAS)*0.5D0
      ABSOVR=(ABSOVR+DABS(ERAS))*0.5D0
      NF=0
      NOCYC=40
      IF (I.GE.IDB) NOCYC=IMAX-I
      GO TO 399
C
  350 ERAS=PA(I)*PB(I)
      ERAS=ERAS+ERAS
      OVERLP=OVERLP+ERAS
      ABSOVR=ABSOVR+DABS(ERAS)
      NF=1
      NOCYC=NOCYC-2
      IF (NOCYC) 399,351,399
  351 NF=-1
      GO TO 399
C
  360 ERAS=PA(I)*PB(I)
      OVERLP=OVERLP+ERAS
      ABSOVR=ABSOVR+DABS(ERAS)
      NF=0
C
  399 CONTINUE
C
C          COMPLETE CALCULATION OF INTEGRAL
C
  400 FRAC=OVERLP/ABSOVR
      OVERLP=OVERLP*(R(IMAX)-R(IMAX-1))/0.75D0
      IF (ISLI.GT.0) GO TO 450
  410 WRITE (9,41) IZ,KUT,CONF1,CONF2,HA,HB,OVERLP,FRAC
   41 FORMAT (2I4,6X, 2(3X,3A6,1X), 5H    (,A3,1H/,A3,2H)=,F14.9,F15.6)
  450 RETURN
C
 1000 WRITE (9,1001) NLH(N),CONF
 1001 FORMAT (//9H ORBITAL ,A3, 23H NOT FOUND FOR CONFIG  ,3A6)
      GO TO 450
C
      END
      SUBROUTINE ZETA2
C
C          CALCULATE SPIN-ORBIT CONFIGURATION-INTERACTION INTEGRALS
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KBGKS=41,KCASE=50,KM=1801,KO=20)
      PARAMETER (KCO=100,KP=KCO+1)
      DIMENSION PA(KM),PB(KM),PC(KM),PD(KM)
C
      COMMON NLZ(4),SUBR0,SUBR,NCFG(2),MMIN,NCDES,NKK(4),L(4),
     1  M,IEXCH,IZ1,IDUM0,OVERLP,X(KM),DRUDR(KM),DERU(KM),
     2  V(KM),ZEI(KM),ZETA,ZETAK
C
      COMMON/C1/FACT(5),IFACT(5),IDUMMY2,DUMMY2(20)
      COMMON/C3/ICON,ISLI,IDIP,IC,ICN(2),NCONFN(2),NC0,ICASE,
     1  NSCONF(3,2),DEL(KCO),IDIPP,IOVFACT,OVFACT
      COMMON/CP/P(KM,4),WNL(2),RM1(KM),RMKP1(KM)
      COMMON/LCM1/NT2(KO),LT2(KO),WT2(KO),NCSPVS,IDUM2,VPAR(50),
     1  NCONF,IZ,Z,ION,IREL,HXID,C,MESH,IDB,EXF,CORRF,KUT,NPAR
      COMMON/LCM1A/NLBCD
      COMMON/LCM1B/NNLZ(KO),NKKK(KO),EE(KO),EE8,R(KM),RU(KM),
     1  PNL(KM,KO),NPNL
      CHARACTER DH*80,OPTION*28,NLH(4)*3,HA*3,HB*3,HC*3,HD*3,HA1*3,
     1  HB1*3,HE*3,NLBCD(KO)*8,NLHT(8)*3,NLB(8,KCO)*3,KST(KCASE,2)*3
      CHARACTER*6 CONF(3),CST(KCASE,6),CONF1(3),CONF2(3),
     1  CF1(KP),CF2(KP),CF3(KP)
      COMMON/CDH/DH,OPTION,NLH,HA1,HB1,HE,NLHT,NLB,KST,
     1  CONF,CST,CONF1,CONF2,CF1,CF2,CF3
      EQUIVALENCE (P,PA),(P(1,2),PB),(P(1,3),PC),(P(1,4),PD),
     1  (L(1),LA),(L(2),LB),(L(3),LC),(L(4),LD),
     2  (NLH(1),HA),(NLH(2),HB),(NLH(3),HC),(NLH(4),HD)
C
  111 IF (SUBR-SUBR0) 114,113,114
  113 IF (NCFG(2)-NCONF) 114,150,114
  114 SUBR0=SUBR
      WRITE (9,10)
   10 FORMAT (11H1RCN MOD 36///)
C
C
C          READ WAVEFUNCTIONS FROM TAPE IC AND SELECT THOSE DESIRED
C
  150 N=1
      K=1
  205 IC=ICN(N)
      IF (IC.GT.30) K=IC-30
      IF (NCFG(N).EQ.NCONF) GO TO 220
      IF (NCFG(N)-NCONFN(K)) 210,210,211
C     IF (NCFG(N)-NCONFN(K)) 210,206,211
C 206 BACKSPACE IC
C     GO TO 211
  210 REWIND IC
  211 READ (IC) NT2,LT2,WT2,NCSPVS,VPAR,CONF,NCONF,IZ,Z,
     1  ION,IREL,HXID,MESH,C,IDB,EXF,CORRF,KUT,NPAR,
     2  NLBCD,NNLZ,NKKK,EE,EE8,NPNL,R,RU, ((PNL(I,M),I=1,KM),M=1,NPNL)
      IF (NCONF.NE.NCFG(N)) GO TO 211
      NCONFN(K)=NCONF
C
  220 DO 221 I=1,3
      CONF1(I)=CONF2(I)
  221 CONF2(I)=CONF(I)
C
      DO 222 M=1,NCSPVS
      IF (NLBCD(M)(1:3).EQ.NLH(N)) GO TO 225
  222 CONTINUE
      GO TO 1000
  225 NKK(N)=NKKK(M)
      NLZ(N)=NNLZ(M)
      L(N)=NNLZ(M)-100*(NNLZ(M)/100)
C
      N=N+1
      JJ=M+NPNL-NCSPVS
      IF (JJ.LE.0) GO TO 1000
      IF (N.GT.2) GO TO 350
  230 DO 231 I=1,MESH
  231 PA(I)=PNL(I,JJ)
      GO TO 205
C
  350 WRITE (9,35) CONF1,CONF2,IZ,KUT
   35 FORMAT ( 1H ,3A6,5X,3A6,6X,2HZ=,I3,6X,4HKUT=,I2///
     1  6X,32HNL   NL  ZETA(RYD)    ZETA(CM-1) )
C
C          CALCULATE V AND DERIVATIVE OF RU
C
  400 DR3=3.D0*(R(2)-R(1))
      KKK=MIN0(NKK(1),NKK(2))
      NBLOCK=KKK/40
      I=1
      DO 420 NB=1,NBLOCK
      DO 410 J=1,10
      I=I+4
      DRUDR(I-3)=(-5.5D0*RU(I-3)+9.D0*RU(I-2)-4.5D0*RU(I-1)+RU(I))/DR3
      DRUDR(I-2)=(-RU(I-3)-1.5D0*RU(I-2)+3.D0*RU(I-1)-0.5D0*RU(I))/DR3
      DRUDR(I-1)=(0.5D0*RU(I-3)-3.D0*RU(I-2)+1.5D0*RU(I-1)+RU(I))/DR3
  410 DRUDR(I)=(-RU(I-3)+4.5D0*RU(I-2)-9.D0*RU(I-1)+5.5D0*RU(I))/DR3
      IF (I.GT.IDB) GO TO 420
      DR3=DR3+DR3
  420 CONTINUE
C
      TWOZ=Z+Z
      DO 450 I=2,KKK
      V(I)=RU(I)/R(I)
  450 DERU(I)=(DRUDR(I)-V(I))/(R(I)**2)
C
      DO 520 I=2,KKK
      PSQ=PA(I)*PNL(I,JJ)
  520 ZEI(I)=PSQ*DERU(I)
C
C          INTEGRATE
C
      ZEI(1)=0.D0
      XIF=0.00125D0*C*5.D0/288.D0
      J=1
      S1Z=0.D0
  530 MM=8
      IF (I.GT.IDB) GO TO 531
      XIF=XIF+XIF
  531 CONTINUE
      S2Z=0.D0
  532 S2Z=S2Z+19.D0*(ZEI(J)+ZEI(J+5))+75.D0*(ZEI(J+1)+ZEI(J+4))
     1       +50.D0*(ZEI(J+2)+ZEI(J+3))
  533 J=J+5
      MM=MM-1
      IF (MM) 534,534,532
  534 S1Z=S1Z+XIF*S2Z
      IF (J-KKK) 530,535,535
  535 ZETA=2.66246D-5*S1Z
      ZETAK=ZETA*109737.31D0
C
      WRITE (9,70) HA,HB,ZETA,ZETAK
   70 FORMAT (5X,A3,2X,A3,F11.5,F14.3////)
C
C
  900 RETURN
C
 1000 WRITE (9,1001) NLH(N),CONF
 1001 FORMAT (//9H ORBITAL ,A3, 23H NOT FOUND FOR CONFIG  ,3A6)
       GO TO 900
      END
      SUBROUTINE DIP
C
C          CALCULATE ELECTRIC DIPOLE OR QUADRUPOLE RADIAL INTEGRALS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KBGKS=41,KCASE=50,KM=1801,KO=20)
      PARAMETER (KCO=100,KP=KCO+1)
      DIMENSION PA(KM),PB(KM),PC(KM),PD(KM)
C
      COMMON NLZ(4),SUBR0,SUBR,NCFG(2),MMIN,NCDES,NKK(4),L(4),
     1  M,IEXCH,IZ1,IDUM0,OVERLP,X(KM),DIPOLE,ABSDIP,FRAC,
     2  REDDIP
C
      COMMON/C2/EX(KP),STN(KP),WT(KP),GAMMA(KP),Q(KP),
     1  Q2(KP),PDIPSC(KP),RRSC(KP,12),NCD,NCPD,NFIRST,NCONT1,IPWR
      COMMON/C3/ICON,ISLI,IDIP,IC,ICN(2),NCONFN(2),NC0,ICASE,
     1  NSCONF(3,2),DEL(KCO),IDIPP,IOVFACT,OVFACT
      COMMON/CP/P(KM,4),WNL(2),RM1(KM),RMKP1(KM)
      COMMON/LCM1/NT2(KO),LT2(KO),WT2(KO),NCSPVS,IDUM2,VPAR(50),
     1  NCONF,IZ,Z,ION,IREL,HXID,C,MESH,IDB,EXF,CORRF,KUT,NPAR
      COMMON/LCM1A/NLBCD
      COMMON/LCM1B/NNLZ(KO),NKKK(KO),EE(KO),EE8,R(KM),RU(KM),
     1  PNL(KM,KO),NPNL
      COMMON/LCM2/RDIPS(KBGKS,KCASE),ST(KCASE,5)
      COMMON/GENOSC/IGEN,NBIGKS,RJ,BIGKS(KBGKS),IBK
      CHARACTER DH*80,OPTION*28,NLH(4)*3,HA*3,HB*3,HC*3,HD*3,HA1*3,
     1  HB1*3,HE*3,NLBCD(KO)*8,NLHT(8)*3,NLB(8,KCO)*3,KST(KCASE,2)*3
      CHARACTER*6 CONF(3),CST(KCASE,6),CONF1(3),CONF2(3),
     1  CF1(KP),CF2(KP),CF3(KP)
      COMMON/CDH/DH,OPTION,NLH,HA1,HB1,HE,NLHT,NLB,KST,
     1  CONF,CST,CONF1,CONF2,CF1,CF2,CF3
      EQUIVALENCE (P,PA),(P(1,2),PB),(P(1,3),PC),(P(1,4),PD),
     1  (L(1),LA),(L(2),LB),(L(3),LC),(L(4),LD),
     2  (NLH(1),HA),(NLH(2),HB),(NLH(3),HC),(NLH(4),HD)
C
      ITEST=1
  100 IF (SUBR.EQ.SUBR0) GO TO 150
  101 SUBR0=SUBR
      IF (IGEN.LE.8) WRITE (9,10)
   10 FORMAT (1H1//19H  Z  KUT           ,7X,5HCONF1,15X,5HCONF2,8X,
     1  10HR INTEGRAL,6X,22HREDUCED MUPOLE ELEMENT,7X,4HFRAC,4X,
     2  8H OVFACT ,12H  OVER   FR1/
     3  19H --- --- ---- ---- ,20H------------------- ,
     4  20H------------------- ,10H----------,3X,14H--------------,
     5  14H--------------,3X,6H------,2X, 9H---------,12H ----- -----/)
      GO TO 200
C
  150 IF (NCFG(1).NE.NCFG1.OR.HA.NE.HA1) GO TO 200
      N=2
      GO TO 251
C
C          READ WAVEFUNCTIONS FROM TAPE IC AND SELECT THOSE DESIRED
C
  200 N=1
      IC=ICN(1)
      NCFG1=NCFG(1)
      HA1=HA
      K=1
      WRITE (9,20)
   20 FORMAT (1H )
  205 IF (IC.GT.30) K=IC-30
      IF (NCFG(N).EQ.NCONF) GO TO 220
      IF (NCFG(N)-NCONFN(K)) 210,210,211
C     IF (NCFG(N)-NCONFN(K)) 210,206,211
C 206 BACKSPACE IC
C     GO TO 211
  210 REWIND IC
  211 READ (IC) NT2,LT2,WT2,NCSPVS,VPAR,CONF,NCONF,IZ,Z,
     1  ION,IREL,HXID,MESH,C,IDB,EXF,CORRF,KUT,NPAR,
     2  NLBCD,NNLZ,NKKK,EE,EE8,NPNL,R,RU, ((PNL(I,M),I=1,KM),M=1,NPNL)
      IF (NCONF.NE.NCFG(N)) GO TO 211
      NCONFN(K)=NCONF
C
  220 IF (N.GT.1) GO TO 223
  221 DO 222 I=1,3
  222 CONF1(I)=CONF(I)
      GO TO 225
  223 DO 224 I=1,3
  224 CONF2(I)=CONF(I)
  225 IF (MMIN.LT.0) GO TO 250
      DO 226 M=1,NCSPVS
      IF (NLBCD(M)(1:3).EQ.NLH(N)) GO TO 230
  226 CONTINUE
      GO TO 1000
C
  230 NKK(N)=NKKK(M)
      NLZ(N)=NNLZ(M)
      WNL(N)=WT2(M)
      L(N)=NNLZ(M)-100*(NNLZ(M)/100)
      J=M+NPNL-NCSPVS
      IF (J.LE.0) GO TO 1000
      DO 231 I=1,MESH
  231 P(I,N)=PNL(I,J)
C
  250 N=N+1
  251 IC=ICN(2)
      IF (N-2) 200,205,300
C
C          BEGIN INTEGRATION LOOP
C
  300 PWRT=2*IPWR+1
      PWRT=DSQRT(PWRT)
      DO 450 IBK=1,NBIGKS
      IF (MMIN.LT.0) THEN
  500   HA='   '
        HB='   '
        HA1='   '
        DIPOLE=0.D0
        REDDIP=0.D0
        FRAC=0.D0
        SIGSQ=0.D0
        OVERLP=0.D0
        FRAC1=0.D0
        GO TO 410
      ENDIF
      DIPOLE=0.D0
      ABSDIP=0.D0
      RPWR=0.D0
      NF=-1
      OVERLP=0.D0
      ABSOVR=0.D0
      IMAX=MIN0(NKK(1),NKK(2))
      DO 399 I=1,IMAX
      RIP=R(I)
      IF (I.EQ.1) GO TO 320
      RPWR=RIP**IPWR
      IF (IGEN.LT.5) GO TO 320
C          CALCULATE SPHERICAL BESSEL FUNCTION
      BIGK=BIGKS(IBK)
      BIGK2=BIGK**2
      CALL SBESS(RIP,ITEST,IPWR,BIGK2,SZ,SZP)
      RPWR=PWRT*SZ/(RIP*BIGK2)
C
  320 CONTINUE
      IF (NF) 340,350,360
  340 ERAS=1.5D0
      IF (I-1) 342,341,342
  341 ERAS=1.D0
      GO TO 344
  342 IF (I-IMAX) 344,343,344
  343 ERAS=0.5D0
  344 ERAS1=PA(I)*PB(I)*ERAS
      ERAS=RPWR*ERAS1
      DIPOLE=(DIPOLE+ERAS)*0.5D0
      ABSDIP=(ABSDIP+DABS(ERAS))*0.5D0
      OVERLP=0.5D0*(OVERLP+ERAS1)
      ABSOVR=0.5D0*(ABSOVR+DABS(ERAS1))
      NF=0
      NOCYC=40
      IF (I.GE.IDB) NOCYC=IMAX-I
      GO TO 399
C
  350 ERAS1=2.D0*PA(I)*PB(I)
      ERAS=RPWR*ERAS1
      DIPOLE=DIPOLE+ERAS
      ABSDIP=ABSDIP+DABS(ERAS)
      OVERLP=OVERLP+ERAS1
      ABSOVR=ABSOVR+DABS(ERAS1)
      NF=1
      NOCYC=NOCYC-2
      IF (NOCYC) 399,351,399
  351 NF=-1
      GO TO 399
C
  360 ERAS1=PA(I)*PB(I)
      ERAS=RPWR*ERAS1
      DIPOLE=DIPOLE+ERAS
      ABSDIP=ABSDIP+DABS(ERAS)
      OVERLP=OVERLP+ERAS1
      ABSOVR=ABSOVR+DABS(ERAS1)
      NF=0
C
  399 CONTINUE
C
C          COMPLETE CALCULATION OF INTEGRAL
C
  400 FRAC=DIPOLE/ABSDIP
      FRAC1=OVERLP/ABSOVR
      FR1=1.333333333D0*(R(IMAX)-R(IMAX-1))
      DIPOLE=DIPOLE*FR1
      OVERLP=OVERLP*FR1
      IF (IPWR.EQ.0) DIPOLE=DIPOLE-OVERLP/BIGK
      DIPOLE=DIPOLE*OVFACT
      A=LA
      B=LB
      C=IPWR
          ZERO=0.D0
          TEMPS=S3J(A,B,C,ZERO,ZERO)
      REDDIP=((-1.D0)**LA)*DSQRT((2.D0*A+1.D0)*(2.D0*B+1.D0))
     1  *TEMPS*DIPOLE
      PDIPSC(NCPD)=REDDIP
      LGR=MAX0(LA,LB)
      FLGR=LGR
      SIGSQ=DIPOLE*DIPOLE/(4.D0*FLGR*FLGR-1.D0)
  410 IF (IGEN.LT.5) GO TO 420
      IF ((IGEN.LE.5.AND.REDDIP.NE.0.D0).OR.(IGEN.LE.8.AND.IBK.EQ.1))
     1  GO TO 420
      GO TO 430
  420 WRITE (9,41) IZ,KUT,CONF1,CONF2,
     1  DIPOLE,  HA,RJ,IPWR,HB, REDDIP , FRAC, OVFACT,OVERLP,FRAC1
   41 FORMAT (1X,2I3,11X,    2(1X,3A6,1X), F11.6, 3X, 1H(,A3,
     1  2H//,A1,I1,2H//,A3, 2H)=, F10.5, 3H A0, F9.4, F11.6,2F6.3)
  430 IF (NCFG(2).GE.NCONT1+NC0) REDDIP=0.D0
      IOVER=100.4D0*OVERLP
      IF1=100.4D0*FRAC1
      IF (IGEN.GE.5) GO TO 440
      IF (IDIP.NE.8) GO TO 432
      FRAC=1.D0
      GO TO 435
  432 IF (IDIP.LT.7) GO TO 435
      I=NCFG(2)-NC0
      REDDIP=DEL(I)*REDDIP
      IF (IDIP.LT.9) FRAC=DEL(I)
      WRITE (9,42) DEL(I),REDDIP
   42 FORMAT (30X,9H SQRTDEL=,F10.4,38X,F10.5)
  435 IF (FRAC.GE.0.D0.AND.FRAC.LT.9.99D0) GO TO 437
      WRITE (11,43) CONF1,CONF2,REDDIP,HA,IPWR,HB,FRAC,HXID,IOVER,IF1
   43 FORMAT (3A6,2X,3A6,F12.5,1H(,A3,3H//R,I1,2H//,A3,1H),F6.3,A2,2I4)
      GO TO 450
  437 WRITE (11,44) CONF1,CONF2,REDDIP,HA,IPWR,HB,FRAC,HXID,IOVER,IF1
   44 FORMAT (3A6,2X,3A6,F12.5,1H(,A3,3H//R,I1,2H//,A3,1H),F6.4,A2,2I4)
      GO TO 450
C          STORE RESULTS IN BESSEL-FUNCTION CASE
  440 IF (IBK.GT.1) GO TO 449
      ICASE=ICASE+1
      IF (ICASE.GT.KCASE) STOP 44
      DO 447 ICF=1,3
      CST(ICASE,ICF)=CONF1(ICF)
  447 CST(ICASE,ICF+3)=CONF2(ICF)
      ST(ICASE,1)=IPWR
      ST(ICASE,2)=FRAC
      ST(ICASE,3)=HXID
      ST(ICASE,4)=IOVER
      ST(ICASE,5)=IF1
      KST(ICASE,1)=HA
      KST(ICASE,2)=HB
C
  449 RDIPS(IBK,ICASE)=REDDIP
  450 CONTINUE
      RETURN
C
 1000 WRITE (9,1001) NLH(N),CONF
 1001 FORMAT (//9H ORBITAL ,A3, 23H NOT FOUND FOR CONFIG  ,3A6)
      RETURN
C
      END
C
      SUBROUTINE SLI2
C
C          CALCULATE CONFIGURATION-INTERACTION COULOMB INTEGRALS RK
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KBGKS=41,KCASE=50,KM=1801,KO=20)
      PARAMETER (KCO=100,KP=KCO+1)
      DIMENSION PA(KM),PB(KM),PC(KM),PD(KM)
C
      COMMON NLZ(4),SUBR0,SUBR,NCFG(2),MMIN,NCDES,NKK(4),L(4),
     1  M,IEXCH,IZ1,IDUM0,OVERLP,X(KM),
     2  XI(KM),XJ(KM),XA(KM),XB(KM),F1(KM),F2(KM),RKP1(KM)
      COMMON/C1/FACT(5),IFACT(5),NRK,RKAYSC(20)
      COMMON/C2/EX(KP),STN(KP),WT(KP),GAMMA(KP),Q(KP),
     1  Q2(KP),PDIPSC(KP),RRSC(KP,12),NCD,NCPD,NFIRST,NCONT1,IPWR
      COMMON/C3/ICON,ISLI,IDIP,IC,ICN(2),NCONFN(2),NC0,ICASE,
     1  NSCONF(3,2),DEL(KCO),IDIPP,IOVFACT,OVFACT
      COMMON/C4/IRHO,ISIG,IRHOP,ISIGP,N12,IDUM4,T0
      COMMON/CP/P(KM,4),WNL(2),RM1(KM),RMKP1(KM)
      COMMON/LCM1/NT2(KO),LT2(KO),WT2(KO),NCSPVS,IDUM2,VPAR(50),
     1  NCONF,IZ,Z,ION,IREL,HXID,C,MESH,IDB,EXF,CORRF,KUT,NPAR
      COMMON/LCM1A/NLBCD
      COMMON/LCM1B/NNLZ(KO),NKKK(KO),EE(KO),EE8,R(KM),RU(KM),
     1  PNL(KM,KO),NPNL
      COMMON/LCM3/PE(KM),WA1,WB1
      CHARACTER DH*80,OPTION*28,NLH(4)*3,HA*3,HB*3,HC*3,HD*3,HA1*3,
     1  HB1*3,HE*3,NLBCD(KO)*8,NLHT(8)*3,NLB(8,KCO)*3,KST(KCASE,2)*3
      CHARACTER*6 CONF(3),CST(KCASE,6),CONF1(3),CONF2(3),
     1  CF1(KP),CF2(KP),CF3(KP)
      COMMON/CDH/DH,OPTION,NLH,HA1,HB1,HE,NLHT,NLB,KST,
     1  CONF,CST,CONF1,CONF2,CF1,CF2,CF3
      EQUIVALENCE (P,PA),(P(1,2),PB),(P(1,3),PC),(P(1,4),PD),
     1  (L(1),LA),(L(2),LB),(L(3),LC),(L(4),LD),
     2  (NLH(1),HA),(NLH(2),HB),(NLH(3),HC),(NLH(4),HD)
C
      DIMENSION K1(10,2),RRYD(10,2),RKAY(10,2),FRAC(10,2),EAV(2),EPS(2)
      DIMENSION RR1(10,2,3),RK1(10,2,3),KK(4),WW(4),NCFGO(2),NN(4)
      EQUIVALENCE (KK(1),KA),(KK(2),KB),(KK(3),KC),(KK(4),KD),
     1            (WW(1),WA),(WW(2),WB),(WW(3),WC),(WW(4),WD)
C
C
      IF (P(1,1).EQ.(-99.D0)) NCFGO(1)=-7
      IF (P(1,1).EQ.(-99.D0)) NCSO=-7
      IF (NCFG(1).NE.NCFGO(1).OR.NCFG(2).NE.NCFGO(2)) ITEST=0
      NCS=1000*NCFG(1)+NCFG(2)
      IF (ISLI.EQ.0) GO TO 110
      WRITE (9,95)
      GO TO 120
  110 WRITE (9,10)
   10 FORMAT (11H1RCN MOD 36//)
C          DETERMINE WHETHER INTERACTION IS CLASS 1B,7, OR 8
C
  120 ICLASS=0
      LA1=LA
      LB1=LB
C
      DO 125 I=1,4
C 125 DECODE (3,12,NLH(I)) L(I)
  125 READ (NLH(I),12) L(I)
   12 FORMAT (2X,A1)
      IF (HA.EQ.HC) GO TO 150
      IF (HB.EQ.HD) GO TO 130
      IF (HB.EQ.HC.OR.HA.EQ.HD) GO TO 140
      GO TO 150
  130 IF (LA.NE.LC.OR.HB(1:2).EQ.'99') GO TO 150
      HE=HB
      HB=HA
      HA=HE
      HE=HD
      HD=HC
      HC=HE
      ICLASS=8
      IF (HA.EQ.HB.OR.HC.EQ.HD) ICLASS=1
      GO TO 148
  140 IF (LA.NE.LD.OR.LB.NE.LC) GO TO 150
      ICLASS=7
      IF (HA.EQ.HD) GO TO 145
      HE=HB
      HB=HA
      HA=HE
      GO TO 148
  145 HE=HD
      HD=HC
      HC=HE
  148 N12=1
C
C          READ WAVEFUNCTIONS FROM TAPE IC AND SELECT THOSE DESIRED
C
  150 IF (NCFG(2).GT.NCFG(1)+1.AND.ICLASO.EQ.7) N12=1
      N=N12
      NPOS=0
      IF (N12.EQ.2) NPOS=NPOS1
      LA=LA1
      LB=LB1
      WA=WA1
      WB=WB1
      EEE=0.D0
      IF (N12.EQ.2) EEE=EE1
      CALL SECONDS(T1)
      IREAD=0
  151 IF (NCFG(N12)-NCONF) 160,220,170
  160 REWIND IC
  170 READ (IC) NT2,LT2,WT2,NCSPVS,VPAR,CONF,NCONF,IZ,Z,
     1  ION,IREL,HXID,MESH,C,IDB,EXF,CORRF,KUT,NPAR,
     2  NLBCD,NNLZ,NKKK,EE,EE8,NPNL,R,RU, ((PNL(I,M),I=1,KM),M=1,NPNL)
      EAV(N)=VPAR(1)
      EPS(N)=EE(NCSPVS)
      IREAD=IREAD+1
      IF (NCONF.NE.NCFG(N)) GO TO 170
      K=1
      IF (IC.GT.30) K=IC-30
      NCONFN(K)=NCONF
C
  220 DO 221 I=1,3
      IF (N12.EQ.1) CONF1(I)=CONF2(I)
  221 CONF2(I)=CONF(I)
C
      DO 229 J=1,2
      K=2*N-2+J
      DO 222 M=1,NCSPVS
      IF (NLBCD(M)(1:3).EQ.NLH(K)) GO TO 225
  222 CONTINUE
      GO TO 1000
  225 NKK(K)=NKKK(M)
      NLZ(K)=NNLZ(M)
      L(K)=NNLZ(M)-100*(NNLZ(M)/100)
      WW(K)=WT2(M)
      NN(K)=NT2(M)
      KK(K)=M
      IP=MESH
      JJ=M+NPNL-NCSPVS
      IF (JJ.LE.0) GO TO 1000
      DO 226 I=1,IP
  226 P(I,K)=PNL(I,JJ)
  229 CONTINUE
      M=KK(2*N)
      IF (NT2(M).GE.99) NPOS=NPOS+1
      IF (N.EQ.1) NPOS1=NPOS
      EEE=EEE+0.5D0*EE(M)
      IF (ISLI.LE.0) WRITE (9,22) EE(M),EEE
   22 FORMAT (8H EE,EEE=,2F10.5)
C
      N=N+1
      WA1=WA
      WB1=WB
      IF (N.GT.2) GO TO 250
      EE1=EEE
  230 IF (NCFG(2)-NCONF) 160,231,170
  231 HC=HA
      HD=HB
      GO TO 220
C
C          DETERMINE WHETHER INTERACTION IS CLASS 1A,6, OR 11
C
  250 IF (ICLASS.NE.0) GO TO 300
      IF (HA.NE.HC) GO TO 300
      IF (LB.NE.LD) GO TO 300
      IF ((WB+WD).EQ.2.D0) ICLASS=11
      IF (ISIG.NE.ISIGP) ICLASS=6
      IF (HA.EQ.HB.OR.HC.EQ.HD) ICLASS=1
C
C          SET UP K TABLES FOR DIRECT (D) AND EXCHANGE (E) INTEGRALS
C
  300 KDMIN=MAX0(IABS(LA-LC),IABS(LB-LD))
      KEMIN=MAX0(IABS(LA-LD),IABS(LB-LC))
      KDMAX=MIN0(LA+LC, LB+LD)
      KEMAX=MIN0(LA+LD, LB+LC)
      KMIN=MIN0(KDMIN,KEMIN)
      KMAX=MAX0(KDMAX,KEMAX)
      NOKD=(KDMAX-KDMIN)/2+1
      NOKE=(KEMAX-KEMIN)/2+1
      NOKMAX=MAX0(NOKD,NOKE)
      CALL SECONDS(T2)
C
C          CALCULATE INTEGRATION LIMITS
C
      I1MAX=MAX0( MIN0(NKK(1),NKK(3)), MIN0(NKK(1),NKK(4)) )
      I2MAX=MAX0( MIN0(NKK(2),NKK(4)), MIN0(NKK(2),NKK(3)) )
      IMAX=MAX0(I1MAX,I2MAX)
C                                                                       SLI10650
C          CALC PNL-PRODUCTS FOR DIRECT OR EXCHANGE CASE                SLI10660
C                                                                       SLI10670
  400 DO 700 N1=1,2                                                     SLI10680
      IF (N1-1) 420,420,412                                             SLI10690
  412 DO 415 I=1,IMAX                                                   SLI10710
      F1(I)=PA(I)*PD(I)                                                 SLI10720
  415 F2(I)=PB(I)*PC(I)                                                 SLI10730
      K=KEMIN-2                                                         SLI10740
      GO TO 500                                                         SLI10750
  420 DO 425 I=1,IMAX                                                   SLI10760
      F1(I)=PA(I)*PC(I)                                                 SLI10770
  425 F2(I)=PB(I)*PD(I)                                                 SLI10780
      K=KDMIN-2                                                         SLI10790
C                                                                       SLI10810
C          CALC FOR EACH VALUE OF K                                     SLI10820
C                                                                       SLI10830
  500 DO 700 N=1,NOKMAX                                                 SLI10840
      K1(N,N1)=0
  520 K=IABS(K+2)                                                       SLI10890
      K1(N,N1)=K                                                        SLI10900
C                                                                       SLI10910
      IF (MMIN.GE.0) GO TO 522
      A1=0.D0
      B1=1.D0
      ICLASS=0
      GO TO 650
  522 RKP1(1)=0.001D0                                                   SLI10920
      RM1(1)=1.D0/R(1)
      DO 530 I=2,IMAX                                                   SLI10930
      RM1(I)=1.D0/R(I)
      IF (N-1) 525,525,529                                              SLI10940
  525 RKP1(I)=R(I)                                                      SLI10950
      IF (K) 530,530,526                                                SLI10960
  526 DO 527 M=1,K                                                      SLI10970
  527 RKP1(I)=RKP1(I)*R(I)                                              SLI10980
      GO TO 530                                                         SLI10990
  529 RKP1(I)=RKP1(I)*(R(I)**2)                                         SLI11000
  530 RMKP1(I)=1.D0/RKP1(I)                                             SLI11010
C                                                                       SLI11020
C          CALC R1 INTEGRAL                                             SLI11030
C                                                                       SLI11040
      XI(1)=0.D0                                                        SLI11050
      XJ(1)=0.D0                                                        SLI11060
      XA(1)=0.D0                                                        SLI11070
      XB(1)=0.D0                                                        SLI11080
      A2=0.D0                                                           SLI11090
      B2=0.D0                                                           SLI11100
      HO12=R(2)/12.D0                                                   SLI11110
      NI=40                                                             SLI11120
      DO 550 I=3,I1MAX,2                                                SLI11130
      A0=A2                                                             SLI11310
      B0=B2                                                             SLI11320
      ERAS=8.D0*F1(I-1)                                                 SLI11140
      A1=ERAS*RKP1(I-1)*RM1(I-1)                                        SLI11150
      B1=ERAS*RMKP1(I-1)                                                SLI11160
      A2=F1(I)*RKP1(I)*RM1(I)                                           SLI11170
      B2=F1(I)*RMKP1(I)                                                 SLI11180
      ERAS=HO12*(5.D0*A0+A1-A2)                                         SLI11190
      XI(I-1)=XI(I-2)+ERAS                                              SLI11200
      XA(I-1)=XA(I-2)+DABS(ERAS)                                        SLI11210
      ERAS=HO12*(5.D0*B0+B1-B2)                                         SLI11220
      XJ(I-1)=XJ(I-2)+ERAS                                              SLI11230
      XB(I-1)=XB(I-2)+DABS(ERAS)                                        SLI11240
      ERAS=HO12*(-A0+A1+5.D0*A2)                                        SLI11250
      XI(I)=XI(I-1)+ERAS                                                SLI11260
      XA(I)=XA(I-1)+DABS(ERAS)                                          SLI11270
      ERAS=HO12*(-B0+B1+5.D0*B2)                                        SLI11280
      XJ(I)=XJ(I-1)+ERAS                                                SLI11290
      XB(I)=XB(I-1)+DABS(ERAS)                                          SLI11300
      NI=NI-2                                                           SLI11330
      IF (NI) 540,540,550                                               SLI11340
  540 HO12=HO12+HO12                                                    SLI11350
      NI=40                                                             SLI11360
      IF (I.GE.IDB) NI=10000
  550 CONTINUE                                                          SLI11370
C                                                                       SLI11380
      DO 560 I=2,I1MAX                                                  SLI11390
      ERAS=RKP1(I)*RM1(I)                                               SLI11400
      XI(I)=XI(I)*RMKP1(I)+(XJ(I1MAX)-XJ(I))*ERAS                       SLI11410
  560 XA(I)=XA(I)*RMKP1(I)+(XB(I1MAX)-XB(I))*ERAS                       SLI11420
      IF (I2MAX-I1MAX) 600,600,569                                      SLI11430
  569 A0=XI(I1MAX)*RKP1(I1MAX)                                          SLI11440
      B0=XA(I1MAX)*RKP1(I1MAX)                                          SLI11450
      DO 570 I=I1MAX,I2MAX                                              SLI11460
      XI(I)=A0*RMKP1(I)                                                 SLI11470
  570 XA(I)=B0*RMKP1(I)                                                 SLI11480
C                                                                       SLI11490
C          CALC R2 INTEGRAL                                             SLI11500
C                                                                       SLI11510
  600 A0=0.D0                                                           SLI11520
      B0=0.D0                                                           SLI11530
      A1=0.D0                                                           SLI11540
      B1=0.D0                                                           SLI11550
      HO3=R(2)/1.5D0                                                    SLI11560
      NI=40                                                             SLI11570
      DO 620 I=3,I2MAX,2                                                SLI11580
      A2=F2(I)*XI(I)                                                    SLI11590
      B2=DABS(F2(I))*XA(I)                                              SLI11600
      ERAS=4.D0*F2(I-1)                                                 SLI11610
      A1=A1+HO3*(A0+ERAS*XI(I-1)+A2)                                    SLI11620
      B1=B1+HO3*(B0+DABS(ERAS)*XA(I-1)+B2)                              SLI11630
      XJ(I)=A1
      A0=A2                                                             SLI11640
      B0=B2                                                             SLI11650
      NI=NI-2                                                           SLI11660
      IF (NI) 610,610,620                                               SLI11670
  610 HO3=HO3+HO3                                                       SLI11680
      NI=40                                                             SLI11690
      IF (I.GE.IDB) NI=10000
  620 CONTINUE                                                          SLI11700
      IF (NPOS.NE.2) GO TO 650
      AMIN=1.0D30
      AMAX=-1.0D30
      JN=0
      JX=0
      I=I2MAX+2
      IX=I2MAX/4
      DO 630 J=1,IX
      I=I-2
      IF (XJ(I).LE.AMAX) GO TO 622
      AMAX=XJ(I)
      JX=JX+1
  622 IF (XJ(I).GE.AMIN) GO TO 624
      AMIN=XJ(I)
      JN=JN+1
  624 IF (XJ(I).LT.AMAX.AND.XJ(I).GT.AMIN.AND.JN.GT.3.AND.JX.GT.3)
     1  GO TO 640
  630 CONTINUE
  640 A1=0.5D0*(AMIN+AMAX)
      IF (ISLI.GT.0) GO TO 650
      WRITE (9,8764) AMIN,AMAX,A1
 8764 FORMAT ( 3F15.6)
C                                                                       
  650 FRAC(N,N1)=A1/B1                                                  
      RRYD(N,N1)=A1                                                     
  700 RKAY(N,N1)=A1*109737.31D0                                         
C
C          OUTPUT
C
      MMINP1=MAX0(MMIN+1,1)
      DO 890 N=1,NOKMAX
      DO 890 N1=1,2
      RR1(N,N1,MMINP1)=RRYD(N,N1)
  890 RK1(N,N1,MMINP1)=RKAY(N,N1)
      IF (ISLI.GT.0) GO TO 935
      IF (MMIN.GT.0) GO TO 930
  900 WRITE (9,90) CONF1,NCFG(1),IZ,ION,EXF,IREL,CONF2,NCFG(2),KUT,
     1  CORRF,ICLASS
   90 FORMAT (/23H0 SLATER INTEGRALS FOR   ,3A6,6X,6HNCONF=,I2,
     1      14X,2HZ=,I3,3X,4HION=,I2,5X,4HEXF=,F5.3,11X,5HIREL=,I2/
     2    23X,3A6,6X,6HNCONF=,I2,22X,4HKUT=,I2,3X,6HCORRF=,F5.3,
     3  10X,6HCLASS ,I2)
  905 IF (NCFG(1)-NCFG(2)) 920,910,920
  910 WRITE (9,91) HA,HB,HA,HB
   91 FORMAT (/1H0,8X,1HK,14X,3HFK(,A3,1H,,A3,1H),21X,4HFRAC,
     1 9X,1HK,14X,3HGK(,A3,1H,,A3,1H),21X,4HFRAC)
      GO TO 930
  920 WRITE (9,92) HA,HB,HC,HD,HA,HB,HC,HD
   92 FORMAT (
     1  /1H0,8X,1HK,10X,4HRDK(,A3,1H,,A3,1H,,A3,1H,,A3,1H),16X,4HFRAC,
     2       9X,1HK,10X,4HREK(,A3,1H,,A3,1H,,A3,1H,,A3,1H),16X,4HFRAC)
  930 WRITE (9,93)
   93 FORMAT (1H ,8X,43H-     -------------------------------------,
     1 8H   -----,9X,43H-     -------------------------------------,
     2 8H   -----//)
  935 IF (ISLI.GT.1) GO TO 945
      IF (ISLI.EQ.1.AND.ICLASS.NE.0.AND.MMIN.GE.0.AND.MMIN.LT.3)
     1  GO TO 945
      DO 940 N=1,NOKMAX
  940 WRITE (9,94) (K1(N,N1),RRYD(N,N1),RKAY(N,N1),
     1      FRAC(N,N1), N1=1,2)
   94 FORMAT (I10,F16.8,7H RYD  =,F14.4,5H CM-1,F8.3,
     1        I10,F16.8,7H RYD  =,F14.4,5H CM-1,F8.3)
      WRITE (9,95)
   95 FORMAT (1H0)
C
  945 AA=0.D0
      A=0.D0
      IF (ICLASS.EQ.0) GO TO 980
      MMIN=MMIN+1
      GO TO (950,960,965,970), MMIN
  950 DO 951 I=1,IP
      PE(I)=PC(I)
  951 PC(I)=PA(I)
      GO TO 400
  960 DO 962 I=1,IP
      PC(I)=PE(I)
      PE(I)=PA(I)
  962 PA(I)=PC(I)
      GO TO 400
  965 DO 966 N=1,NOKMAX
      DO 966 N1=1,2
      RRYD(N,N1)=(RRYD(N,N1)+RR1(N,N1,2))/2.D0
  966 RKAY(N,N1)=(RKAY(N,N1)+RK1(N,N1,2))/2.D0
      IF (ISLI.GT.0) GO TO 935
      GO TO 930
  970 HE=HA
      HA=HB
      HB=HD
      FLA=LA
      FLB=LB
      HB=HA
      HA=HE
      DO 972 I=1,IP
  972 PA(I)=PE(I)
C
      B=(2.D0*FLB+1.D0)/(4.D0*FLB+1.D0)
      IF (ITEST.NE.0) GO TO 975
      ITEST=7
      AA=0.D0
      A=AA
  975 FK=KEMIN-2
      DO 979 N=1,NOKE
      FK=FK+2.D0
      IF (ICLASS.EQ.1) GO TO 976
      D=WW(1)
      C=D*RRYD(N,2)
      IF (ICLASS.EQ.11) GO TO 1978
      GO TO 978
  976 N2=2
      N3=3
      IF (HA.EQ.HB) GO TO 977
      N2=3
      N3=2
  977 D=WW(N2)-1.D0
      C=0.5D0*RR1(N,2,N3)
      IF (N.GT.1) C=C+B*RR1(N,1,N2)
      C=C*D
  978 A=A+0.5D0*C*S3J0SQ(FLA,FK,FLB)
 1978 CONTINUE
      IF (ISLI.GT.0) GO TO 979
      WRITE (9,96) D,C,FLA,FK,FLB,A,EEE,OVERLP
   96 FORMAT (/8F10.5)
  979 CONTINUE
      RR1(1,1,1)=A/D
      AA=109.73731D0*AA/D
C
  980 M=0
      DO 981 N1=1,2
      NCFGO(N1)=NCFG(N1)
      N2=N1
      IF (ICLASS.EQ.7) N2=3-N1
      DO 981 N=1,NOKE
      M=M+1
  981 RKAY(M,1)=109.73731D0*RR1(N,N2,1)
      N3=(2-N2)*NOKE+1
      M=1
      N=2*NOKD
      IF (HA.EQ.HB.OR.HC.EQ.HD) N=NOKD
      N1=N-M+1
      N2=5
      IF (ISLI.LT.3) WRITE (9,97) CONF1(1),CONF1(2),
     1      CONF1(3),CONF2(2),CONF2(3),N1, (RKAY(I,1),N2, I=M,N),OVFACT
   97 FORMAT (/1H0,2A6,A4,1H-,A6,A4,4X,I2,5(F11.4,1X,I1)/7(F11.4,1X,I1))
      IF (NCS.NE.NCSO) NO=0
C
C         THIS WRITE STATEMENT INCLUDED TO CURE AN APPARENT COMPILER BUG
            IF (ISLI.GT.5) WRITE (9,2356) NCSO,NCS,NO,M,N
 2356 FORMAT (' NCSO,NCS,NO,M,N=',5I5)
C
      DO 982 I=M,N
      J=NO+I
      RKAYSC(J)=FACT(5)*RKAY(I,1)*OVFACT
      IF (NCD.NE.NFIRST) GO TO 982
      RRSC(NCPD,J)=RKAYSC(J)/109.73731D0
  982 CONTINUE
      N=N+NO
      NP1=N+1
      IF (IDIPP.EQ.5.AND.EPS(2).GT.0.D0.AND.DABS(EAV(1)-EAV(2)).
     1  GT.0.01D0) NP1=1
      IF (IDIP.NE.6.OR.(NN(2).LT.99.OR.NN(4).LT.99)) GO TO 983
      NP1=1
      IF (EPS(1).EQ.EPS(2).AND.L(2).EQ.L(4)) NP1=NO+1
  983 DO 985 I=NP1,20
  985 RKAYSC(I)=0.D0
      NRK=N
  990 NO=N
      NCSO=NCS
C
      CALL SECONDS(T3)
      T4=T2-T1
      T5=T3-T2
      T3=T3-T0
      WRITE (9,9988) NCFG,IREAD,T4,T5,T3
 9988 FORMAT (3I5,2F10.3,F10.3,4H SEC)
      ICLASO=ICLASS
      IF (ICLASS.EQ.0.OR.ICLASS.EQ.6.OR.ICLASS.EQ.11) RETURN
      NCONF=1000
      HA='   '
      RETURN
C
 1000 WRITE (9,1001) NLH(K),CONF
 1001 FORMAT (//9H ORBITAL ,A3, 23H NOT FOUND FOR CONFIG  ,3A6)
      HA='   '
      GO TO 990
C
      END
      SUBROUTINE G5INP
C
C          PUNCH INPUT FOR PROGRAM RCG MOD 11 (DIFFERENT FROM MOD 5-9,12)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KBGKS=41,KCASE=50,KM=1801,KO=20)
      PARAMETER (KCO=100,KP=KCO+1)
      DIMENSION NORB(KP),N(25,KP),L(25,KP),IWW(8),
     1  W(25,KP),NN(8,KCO),WW(8,KCO),LSYMB(22),
     2  LB(8,KCO),LL(8,KCO),WWO(8,KCO),LBO(8,KCO),ETOT(KCO),EE8NC(KCO)
C
      COMMON NLZ(4),SUBR0,SUBR,NCFG(2),MMIN,NCDES,
     1  JDUM(8),IDUM(4),OVERLP,
     2  IX(8),FK(8,10),ZE(8),FGK(8,8,10),KPAR(50),IVPAR(50)
      COMMON/C1/FACT(5),IFACT(5),NRK,RKAYSC(20)
      COMMON/C2/EE(KP),STN(KP),WT(KP),GAMMA(KP),Q(KP),
     1  Q2(KP),PDIPSC(KP),RRSC(KP,12),NCD,NCPD,NFIRST,NCONT1,IPWR
      COMMON/C3/ICON,ISLI,IDIP,IC,ICN(2),NCONFN(2),NC0,ICASE,
     1  NSCONF(3,2),DEL(KCO),IDIPP,IOVFACT,OVFACT
      COMMON/C4/IRHO,ISIG,IRHOP,ISIGP,N12,IDUM4,T0
      COMMON/CP/P(KM,4),WNL(2),RM1(KM),RMKP1(KM)
      COMMON/LCM1/NT2(KO),LT2(KO),WT2(KO),NCSPVS,IDUM2,VPAR(50),
     1  NCONF,IZ,Z,ION,IREL,HXID,C,MESH,IDB,EXF,CORRF,KUT,NPAR
      COMMON/LCM1A/NLBCD
      COMMON/LCM1B/NNLZ(KO),NKKK(KO),E(KO),EE8,R(KM),RU(KM),
     1  PNL(KM,KO),NPNL
      COMMON/LCM2/RDIPS(KBGKS,KCASE),ST(KCASE,5)
      COMMON/GENOSC/IGEN,NBIGKS,RJ,BIGKS(KBGKS),IBK
      CHARACTER DH*80,OPTION*28,NLH(4)*3,HA*3,HB*3,HC*3,HD*3,HA1*3,
     1  HB1*3,HE*3,NLBCD(KO)*8,NLHT(8)*3,NLB(8,KCO)*3,KST(KCASE,2)*3
      CHARACTER*6 CONF(3),CST(KCASE,6),CONF1(3),CONF2(3),
     1  CF1(KP),CF2(KP),CF3(KP)
      CHARACTER HA1T*3
      COMMON/CDH/DH,OPTION,NLH,HA1,HB1,HE,NLHT,NLB,KST,
     1  CONF,CST,CONF1,CONF2,CF1,CF2,CF3
      EQUIVALENCE (NLH(1),HA),(NLH(2),HB),(NLH(3),HC),(NLH(4),HD),
     1  (PNL,WWO),(PNL(1,4),LBO),(RDIPS,ETOT),(RDIPS(1,10),EE8NC)
C
      DIMENSION EAV(KP)
      DATA (LSYMB(I),I=1,22)/1HS,1HP,1HD,1HF,1HG,1HH,1HI,1HK,1HL,1HM,
     1  1HN,1HO,1HQ,1HR,1HT,1HU,1HV,1HW,1HX,1HY,1HZ,1HA/
      DATA SLIN,DIPN/3HSLI,3HDIP/
      DATA I2S,BCDR,BCDJ/3H 2S,1HR,1HJ/
C
C
C     KCO=100
C     KP=KCO+1
      CALL SECONDS(T0)
C 100 DECODE (80,10,DH) NCK,IOVFACT,NOCSET,NSCONF(3,1),NSCONF(3,2),
C    1  EAV11,IABG,
C    2  OPTION,IQUAD,IFACT,DMIN,IPRINT,IENGYD,ISPECC,ICON,ISLI,IDIP,ALF
  100 READ(DH,10) NCK,IOVFACT,NOCSET,NSCONF(3,1),NSCONF(3,2),EAV11,IABG,
     1  OPTION,IQUAD,IFACT,DMIN,IPRINT,IENGYD,ISPECC,ICON,ISLI,IDIP,ALF
   10 FORMAT (5X,A2,I1,I2,I1,I2,F7.4,I1,A28,I1,5I2,2A5,2A1,3I1,F5.0)
C
C               CHECK FOR GENERALISED OSCILLATOR STRENGTH OPTION
C                    AND COMPUTE MOMENTUM TRANSFER VALUES
C
C     DECODE (18,11,OPTION) NBIGKS,NENRGS,IGEN,IRNQ,IRXQ,IRND,IRXD
      READ (OPTION,11) NBIGKS,NENRGS,IGEN,IRNQ,IRXQ,IRND,IRXD
   11 FORMAT (1X,I3,I2,5X,5I1)
      IF (IGEN.LT.5.AND.IQUAD.GT.0) IRNQ=2
      IRXQ=MAX0(IRXQ,IRNQ)
      IRND=MAX0(IRND,1)
      IRXD=MAX0(IRXD,IRND)
      RJ=BCDR
      IF (IGEN.LT.5) NBIGKS=1
      IF (IGEN.LT.5) GO TO 80
      IF (NBIGKS.LE.10) NBIGKS=KBGKS
      IF (NENRGS.LE.0) NENRGS=21
      RJ=BCDJ
      OPTION(1:8)='        '
      XMIN=1.01D0
      IF (ICON.EQ.1) XMIN=1.1D0
      IF (ICON.GT.2) XMIN=ICON
      IF (ICON.EQ.9) XMIN=10.D0
      ICON=0
      XMAX=100.D0
      IF (IDIP.EQ.1) XMAX=20.D0
      IF (IDIP.GT.1) XMAX=100*IDIP
      IF (IDIP.EQ.9) XMAX=1000.D0
      IDIP=0
   80 NC0=0
      IDIPP=0
      IF (IDIP.NE.5) GO TO 82
      IDIPP=5
      IDIP=9
   82 NCMAX=KCO
      NOCSET=-IABS(NOCSET)
      IABG=IABS(IABG)
      IQUAD=IABS(IQUAD)
      DO 101 I=1,5
      IF (IFACT(I).EQ.0.AND.I.EQ.2) IFACT(I)=95
      IF (IFACT(I).EQ.0) IFACT(I)=85
      ERAS=IFACT(I)
      IF (ERAS.LT.1.5D0) ERAS=0.1D0
      FACT(I)=0.01D0*ERAS
  101 IF (IFACT(I).EQ.99) FACT(I)=1.D0
      IF (ICON.EQ.0) ICON=8
      NSCF31=IABS(NSCONF(3,1))
      NS32=NSCONF(3,2)
      IF (NSCONF(3,2).EQ.0.AND.IDIP.EQ.8) NSCONF(3,2)=-8
      NSCF32=IABS(NSCONF(3,2))
      IDIP0=IDIP
      NOSUBC=0
      NC3=0
      REWIND 2
      IF (NC0.EQ.0) GO TO 102
      DO 105 NC=1,NC0
  105 READ (2)
C
  102 DO 90 I=1,KP
      STN(I)=0.D0
      WT(I)=0.D0
      GAMMA(I)=0.D0
      Q(I)=0.D0
      Q2(I)=0.D0
      PDIPSC(I)=0.D0
      EAV(I)=0.D0
      DO 90 J=1,12
   90 RRSC(I,J)=0.D0
      NCONT1=1000
      NCONT2=1000
      IDELE=1000
      EMAX=0.D0
C
C          DETERMINE NUMBER OF CONFIGS OF EACH PARITY
C
      REWIND 31
      REWIND 32
      IC=31
      ICN(1)=31
      ICN(2)=31
      NC=0
      NMAX=0
      LMAX=0
      NORBX=0
      IPCT=0
      VPAR1=-99999999.D0
      NC2=0
  110 NC=NC+1
      IF (NC+NC0.GT.NCMAX) GO TO 130
      READ (2,END=112) NT2,LT2,WT2,NCSPVS,VPAR,CONF,NCONF,IZ,Z,
     1  ION,IREL,HXID,MESH,C,IDB,EXF,CORRF,KUT,NPAR,
     2  NLBCD,NNLZ,NKKK,E,EE8,NPNL,R,RU, ((PNL(I,M),I=1,KM),M=1,NPNL),
     3  IW6
      GO TO 115
  112 NCMAX=NC0+NC-1
      GO TO 130
  115 CF1(NC)=CONF(1)
      CF2(NC)=CONF(2)
      CF3(NC)=CONF(3)
      ETOT(NC)=109.73731D0*VPAR(1)
      EE(NC)=E(NCSPVS)
      EE8NC(NC)=EE8
      NORB(NC)=NCSPVS
      IF (IC.NE.32.OR.NS32.NE.-7.OR.IPCT.NE.0) GO TO 118
      IF (NC.GT.NC2+7.AND.VPAR(1).LT.VPAR1) IPCT=NC-NC2-1
  118 VPAR1=VPAR(1)
      IF (IDIP.LT.9.OR.NT2(NCSPVS).LT.99) GO TO 120
      IF (NCONT1.EQ.1000) NCONT1=NC
      IF (EE(NC).EQ.EE(NCONT1)) NCONT2=NC
      IERAS=NCONT2-NCONT1
      IF (IERAS.GT.0) IDELE=1000/IERAS
      IERAS=NC-NCONT2-1
      IF (IERAS.GT.0) IDELE=MIN0(IDELE,3000/IERAS)
  120 IF (NC.EQ.1) IONT=ION
      IF (NC.EQ.1) IZT=IZ
      IF (IONT.EQ.ION.AND.IZT.EQ.IZ) GO TO 132
  130 NC3=NC-1
      IF (NC3.EQ.0) GO TO 900
      BACKSPACE 2
      GO TO 160
  132 IF (NCSPVS.GT.NORBX) NORBX=NCSPVS
      PRTY=0.D0
      DO 135 I=1,NCSPVS
      N(I,NC)=NT2(I)
      L(I,NC)=LT2(I)
      W(I,NC)=WT2(I)
      IF (I.EQ.NCSPVS.AND.W(I,NC).LE.1.D0) GO TO 134
      NMAX=MAX0(NMAX,N(I,NC))
      LMAX=MAX0(LMAX,L(I,NC))
  134 ERAS=L(I,NC)
  135 PRTY=PRTY+W(I,NC)*ERAS
      TWO=2.D0
      PRTY=MOD(PRTY,TWO)
      IF (NC.EQ.1) PRTY1=PRTY
      IF (PRTY.EQ.PRTY1) GO TO 150
      IF (IC.EQ.32) GO TO 130
      PRTY1=PRTY
      NC2=NC-1
      IC=32
      ICN(2)=32
  150 WRITE (IC) NT2,LT2,WT2,NCSPVS,VPAR,CONF,NCONF,IZ,Z,
     1  ION,IREL,HXID,MESH,C,IDB,EXF,CORRF,KUT,NPAR,
     2  NLBCD,NNLZ,NKKK,E,EE8,NPNL,R,RU, ((PNL(I,M),I=1,KM),M=1,NPNL)
      GO TO 110
C
  160 IF (NC3.EQ.1.AND.CF1(1).EQ.'DUMMY ') GO TO 700
      WRITE (9,13)
   13 FORMAT (34H1RCN2---INPUT FOR RCG MOD 10 OR 11///)
      IF (NC2.EQ.0) NC2=NC3
      NSCONF(2,1)=NC2
      NSCONF(2,2)=NC3-NC2
      IF (NSCF31.EQ.0) NSCONF(3,1)=NC2
      IF (NSCF32.EQ.0) NSCONF(3,2)=NSCONF(2,2)
      DO 162 I=1,8
      IWW(I)=0
      DO 162 NC=1,KCO
      LBO(I,NC)=LB(I,NC)
      LB(I,NC)=LSYMB(1)
      LL(I,NC)=0
      NN(I,NC)=0
      WWO(I,NC)=WW(I,NC)
  162 WW(I,NC)=0.D0
C
C          SET UP PARTIALLY-FILLED OR INTERACTING SHELLS
C
C     FIRST, FIND ORBITALS NOT ALWAYS FULL, BUT WITH MORE
C        THAN ONE ELECTRON IN SOME CONFIGURATION,
C          OR FIND SINGLY-OCCUPIED ORBITAL WHICH IS NOT THE OUTERMOST ORBITAL
      J1=0
      IF (NMAX.EQ.0) GO TO 170
      NMIN=1
      DO 169 NNN=NMIN,NMAX
      LMAXP1=MIN0(NNN,LMAX+1)
      DO 169 LP1=1,LMAXP1
      LLL=LP1-1
      WMAX=1.D0
      A=4*LLL+2
      WMIN=A
      DO 166 NC=1,NC3
      NCSPVS=NORB(NC)
      DO 165 I=1,NCSPVS
      IF (N(I,NC).NE.NNN.OR.L(I,NC).NE.LLL) GO TO 165
      WMIN=DMIN1(WMIN,W(I,NC))
      WMAX=DMAX1(WMAX,W(I,NC))
      IF (W(I,NC).EQ.1.D0.AND.I.NE.NCSPVS) WMAX=2.D0
      GO TO 166
  165 CONTINUE
      WMIN=0.D0
  166 CONTINUE
      IF (WMIN.EQ.A.OR.WMAX.EQ.1.D0) GO TO 169
      J1=J1+1
      DO 168 NC=1,NC3
      NN(J1,NC)=NNN
  168 LL(J1,NC)=LLL
  169 CONTINUE
C
C     SECOND, FIND REMAINING (SINGLY-OCCUPIED) ORBITALS
C        (FOR ONE PARITY AT A TIME)
  170 NSUBCO=NOSUBC
      NOSUBC=0
      K=1
      NCN=1
      NCX=NC2
      IF (ICON.NE.1) ICON=MAX0(ICON,J1+1)
  171 J3=J1
      DO 180 NC=NCN,NCX
      NC1=MAX0(NC-1,NCN)
      DO 172 J=1,8
  172 NLB(J,NC)='   '
      NCSPVS=NORB(NC)
      DO 179 I=1,NCSPVS
      IF (J3.EQ.0.AND.I.EQ.NCSPVS.AND.NC.EQ.NCX) GO TO 177
      IF(J3.EQ.0) GO TO 174
      DO 173 J=1,J3
      LLLL=LL(J,NC1)
      LLL2=L(I,NC)
      NNNN=NN(J,NC1)
      NNN2=N(I,NC)
      IF (N(I,NC).EQ.NN(J,NC1).AND.L(I,NC).EQ.LL(J,NC1)) GO TO 178      *
  173 CONTINUE
  174 IF (J3.EQ.J1.OR.NC.EQ.NCN) GO TO 176
      J=J3
      IF (L(I,NC).EQ.LL(J,NC1).AND.I.EQ.NCSPVS) GO TO 178
  176 A=4*L(I,NC)+2
      IF (W(I,NC).EQ.A) GO TO 179
  177 J3=J3+1
      IF (ICON.GT.1) J3=MIN0(J3,ICON)
      DO 1177 NC1=NC,NCX
 1177 LL(J3,NC1)=L(I,NC)
      J=J3
  178 WW(J,NC)=W(I,NC)
      NN(J,NC)=N(I,NC)
      IF (NN(J,NCN).EQ.0) NN(J,NCN)=N(I,NC)                             *
      JJ=LL(J,NC)
      LB(J,NC)=LSYMB(JJ+1)
C     ENCODE (3,16,NLB(J,NC)) NN(J,NC),LB(J,NC)
      WRITE (NLB(J,NC),16) NN(J,NC),LB(J,NC)
   16 FORMAT (I2,A1)
  179 CONTINUE
  180 CONTINUE
      NOSUBC=MAX0(NOSUBC,J3)
      IF (NOSUBC.LE.8) GO TO 181
      WRITE (9,17)
   17 FORMAT (/'*****MORE THAN 8 SUBSHELLS!')
      GO TO 700
  181 IF (NC2.EQ.NC3.OR.K.EQ.2) GO TO 182
      K=2
      NCN=NC2+1
      NCX=NC3
      GO TO 171
C
C          IF CONFIGURATIONS SAME AS FOR PREVIOUS CASE, SET ICTC=1
C
  182 ICTC=1
      IF (NOSUBC.NE.NSUBCO) ICTC=0
      IF (ICTC.EQ.0) GO TO 185
      DO 184 I=1,NOSUBC
      DO 183 NC=1,NC3
      IF (WW(I,NC).NE.WWO(I,NC)) ICTC=0
      IF (WW(I,NC).EQ.0.D0) GO TO 183
      IF (LB(I,NC).NE.LBO(I,NC)) ICTC=0
  183 CONTINUE
  184 CONTINUE
      IF (ICTC.EQ.0) GO TO 185
      BACKSPACE 11
      GO TO 200
C          PUNCH CONTROL CARD
  185 NSCONF(1,1)=NOSUBC
      NSCONF(1,2)=NOSUBC
      IF (NC2.EQ.NC3) NSCONF(1,2)=0
      IF (IW6.GE.0) IW6=6
      IF (IW6.LT.0) IW6=-6
      IF (IDIP.EQ.6) IW6=8*IW6/ABS(IW6)
      WRITE (9,18) NCK,NOCSET,NSCONF,IABG,OPTION,IQUAD,DMIN,IPRINT,
     1  IENGYD,ISPECC,IW6,ICON,ISLI,IDIP,IPCT,ICTC,IREL
   18 FORMAT (4X,1H1,A2,I3,I1,2I2,I1,2I2,I1,A28,I1,1X,9H1000.0000,
     1  2A5,2A1,6I2,I8)
      IF (ICON.NE.1) WRITE(11,18) NCK,NOCSET,NSCONF,IABG,OPTION,IQUAD,
     1  DMIN,IPRINT,IENGYD,ISPECC,IW6,IPCT,ICTC
      IF (NOCSET.LT.0) NOCSET=NOCSET-1
C
C          PUNCH CONFIGURATION CARDS
C
  191 K=1
      NCN=1
      NCX=NC2
  192 IF (NCN.EQ.NCX) GO TO 196
      DO 195 I=1,NOSUBC
      DO 194 NC=NCN,NCX
      IF (WW(I,NC).EQ.0.D0) GO TO 194
      DO 193 M=NCN,NCX
      IF (WW(I,M).NE.0.D0) GO TO 193
      LL(I,M)=LL(I,NC)
      LB(I,M)=LB(I,NC)
  193 CONTINUE
      GO TO 195
  194 CONTINUE
  195 CONTINUE
  196 DO 198 NC=NCN,NCX
      DO 197 I=1,NOSUBC
  197 IWW(I)=WW(I,NC)
      WRITE  (9,19) (LB(I,NC),IWW(I), I=1,8),CF1(NC),CF2(NC),CF3(NC),
     1  ETOT(NC),EE8NC(NC)
   19 FORMAT (1X,8(A1,I2,4X),5X,3A6,5X,4HEAV=,F14.3,3H KK,
     1  5X,5HEKIN=,F9.4,3H RY)
      IF (ICON.NE.1) WRITE (11,20) (LB(I,NC),IWW(I), I=1,8),
     1  CF1(NC),CF2(NC),CF3(NC),ETOT(NC),EE8NC(NC)
   20 FORMAT (8(A1,I2,2X),3A6,F13.3,F9.4)
  198 CONTINUE
      IF (NC2.EQ.NC3.OR.K.EQ.2) GO TO 200
      K=2
      NCN=NC2+1
      NCX=NC3
      GO TO 192
C
C
C          REARRANGE, SCALE, AND PUNCH SINGLE-CONF PARAMETERS
C
  200 REWIND 31
      REWIND 32
      IFLAG=0                                                           
      KK=1
      NFIRST=1
      NLAST=NC2
      IC=31
      DO 495 NC=1,NC3
        ICN(1)=IC
        ICN(2)=IC
        READ (IC) NT2,LT2,WT2,NCSPVS,VPAR,CONF,NCONF,IZ,Z,
     1    ION,IREL,HXID,MESH,C,IDB,EXF,CORRF,KUT,NPAR,
     2    NLBCD,NNLZ,NKKK,E,EE8
        I1=NCSPVS+1
        DO 201 I=I1,25
  201   L(I,NC)=77
        NCFG2=NCONF
        NCONF=99
        NCONFN(KK)=99
        IF (NC.EQ.1) E1=VPAR(1)
        KPAR(1)=0
        VPAR(1)=EAV11+109.73731D0*(VPAR(1)-E1)
        EAV(NC)=VPAR(1)
C            DETERMINE ENERGY WIDTHS FOR PSEUDO-DISCRETE CONFIGURATIONS
        DEL(NC)=1.D0
        IF (NC.GE.NCONT1) DEL(NC)=0.D0
        IF (IDIP.EQ.8.AND.NT2(NCSPVS).GE.99) DEL(NC)=0.D0
        IF (NC.EQ.NFIRST.AND.NT2(NCSPVS).LT.99) GO TO 202
        IF (IDIP.NE.7.OR.W(NCSPVS,NC).NE.1.D0) GO TO 206
        IF (L(NCSPVS,NC).NE.L(NCSPVS,NC-1))
     1    EMAX=EAV(NC)-109.737D0*E(NCSPVS)
  202   ZC=ION+1
        IF (ZC.LE.0.D0) ZC=1.D0
        FN=99.D0
        IF (E(NCSPVS).LT.0.D0) FN=ZC/DSQRT(-E(NCSPVS))
        DEL0=109.73731D0
        IF (FN.EQ.99.D0) GO TO 203
        DEL0=2.D0*DEL0*(ZC**2)/(FN**3)
        IF (NC.EQ.NFIRST) GO TO 1202
        IF (L(NCSPVS,NC).NE.L(NCSPVS,NC-1)) GO TO 206
        IF (N(NCSPVS,NC).GT.N(NCSPVS,NC-1)+1) GO TO 203
 1202   EMAX=EAV(NC)+0.5D0*DEL0
        IF (NC.EQ.NFIRST.OR.NC.EQ.NC2.OR.NC.EQ.NC3) GO TO 206
        NCS1=NCSPVS
        VPAR1=VPAR(1)
        READ (IC) NT2,LT2,WT2,NCSPVS,VPAR(1)
        BACKSPACE IC
        EAVP1=EAV11+109.73731D0*(VPAR(1)-E1)
        NCSPVS=NCS1
        VPAR(1)=VPAR1
        IF (EAVP1.LE.EAV(NC)) GO TO 206
        X=(EAVP1-EAV(NC))/(EAVP1-EAV(NC-1))
        EMAX=MIN(EMAX,X*EAV(NC)+(1.D0-X)*EAVP1)
        GO TO 206
  203   EMIN=EMAX
        EMAX=2.D0*EAV(NC)-EMIN
        IF (NC.EQ.NC2.OR.NC.EQ.NC3) GO TO 205
        NCS1=NCSPVS
        VPAR1=VPAR(1)
        READ (IC) NT2,LT2,WT2,NCSPVS,VPAR(1)
        BACKSPACE IC
        IB=LT2(NCSPVS)
        EAVP1=EAV11+109.73731D0*(VPAR(1)-E1)
        NCS2=NCSPVS
        NCSPVS=NCS1
        VPAR(1)=VPAR1
        IF (NC.EQ.NFIRST.OR.EAV(NC).LT.EAV(NC-1)) GO TO 1203
        IF (NCSPVS.NE.NCS2) GO TO 205
        IF (EE(NC).LT.EE(NC-1)) GO TO 1203
        IF (NC.GT.NFIRST.AND.IB.NE.L(NCSPVS,NC)) GO TO 205
        IF (NC.GT.NFIRST.AND.IB.EQ.L(NCSPVS,NC).AND.
     1    IB.EQ.L(NCSPVS,NC-1)) GO TO 204
 1203   EMAX=0.5D0*(EAV(NC)+EAVP1)
        EMIN=EMAX+EAV(NC)-EAVP1
        GO TO 205
  204   IF (EAVP1.LE.EAV(NC)) GO TO 205
        X=(EAVP1-EAV(NC))/(EAVP1-MIN(EAV(NC-1),EMIN))
        HALF=0.5D0
        X=DMAX1(X,HALF)
        EMAX=DMIN1(EMAX,X*EAV(NC)+(1.D0-X)*EAVP1)
  205   IF (EMAX.LE.EMIN) GO TO 206
        DEL(NC)=DSQRT((EMAX-EMIN)/DEL0)
C
  206   DO 210 I=1,8
        IX(I)=0
        NLHT(I)='   '
        ZE(I)=-1.D0
          DO 210 K=1,10
          FK(I,K)=-1.D0
          DO 210 J=1,8
  210   FGK(I,J,K)=-1.D0
C
        NPAR=1
        JX=0
        DO 220 I=1,NCSPVS
        A=4*L(I,NC)+1
        IF (W(I,NC).LT.1.D0.OR.W(I,NC).GT.A) GO TO 220
        JX=JX+1
        NLHT(JX)=NLBCD(I)(1:3)
        L(JX,NC)=L(I,NC)
        W(JX,NC)=W(I,NC)
        IF (W(I,NC).EQ.1.D0.OR.W(I,NC).EQ.A) GO TO 220
        NK=L(I,NC)
        DO 216 K=1,NK
        NPAR=NPAR+1
  216   FK(JX,K)=VPAR(NPAR)
  220   CONTINUE
C
        IF (JX.EQ.0) GO TO 291
        DO 230 I=1,JX
        IF (L(I,NC).EQ.0) GO TO 230
        NPAR=NPAR+1
        ZE(I)=VPAR(NPAR)
  230   CONTINUE
C                  Calc Fk(Li,Lj)
        IF (JX.EQ.1) GO TO 251
        JM1=JX-1
        NLPAR=NPAR                  ! NUMBER OF LEGAL PARAMETERS
        ILLEG=0                     ! SET 'ILLEGAL' FLAG
        IF (IABG.EQ.2.OR.IABG.EQ.4) ILLEG=1
        DO 240 I=1,JM1
        J1=I+1
        DO 240 J=J1,JX
          NK=MIN0(L(I,NC),L(J,NC))
          IF (NK.EQ.0) GO TO 240
C            DOUBLE MAX NUMBER OF FK - FOR INCLUDING ILLEGAL K'S
          NK=2*NK
          DO 235 K=1,NK
            MODK2=MOD(K,2)
            IF (ILLEG.EQ.0.AND.MODK2.NE.0) GO TO 235
            NPAR=NPAR+1
            IF (MODK2.NE.0) THEN
              IF (ILLEG.NE.0) FGK(I,J,K)=0.D0
            ELSE
              NLPAR=NLPAR+1
              FGK(I,J,K)=VPAR(NLPAR)
            ENDIF
  235     CONTINUE
  240   CONTINUE
C                   Calc Gk(Li,Lj)
        DO 250 I=1,JM1
          J1=I+1
          DO 250 J=J1,JX
          NK=MIN0(L(I,NC),L(J,NC))+1
C             DOUBLE MAX NUMBER OF GK - FOR INCLUDING ILLEGAL K'S
          NK=2*NK-1
          DO 245 K=1,NK
            MODK2=MOD(K,2)
            IF (ILLEG.EQ.0.AND.MODK2.EQ.0) GO TO 245
            NPAR=NPAR+1
            IF (MODK2.EQ.0) THEN
              IF (ILLEG.NE.0) FGK(J,I,K)=0
            ELSE
              NLPAR=NLPAR+1
              FGK(J,I,K)=VPAR(NLPAR)
            ENDIF
  245     CONTINUE
  250   CONTINUE
C              SCALE F(Li,Li)
  251   NPAR=1
        DO 260 I=1,NOSUBC
          DO 256 J=1,JX
            IF (NLHT(J).NE.NLB(I,NC)) GO TO 256
            IX(I)=J
            DO 252 K=1,5
              IF (FK(J,K).LT.0.D0) GO TO 253
              NPAR=NPAR+1
              KPAR(NPAR)=1
  252       VPAR(NPAR)=FK(J,K)*FACT(1)
  253       IF (FK(J,1).LT.0.D0) GO TO 260
            IF (IABG.EQ.0) GO TO 260
            NK=LL(I,NC)
            IF (NK.NE.2) GO TO 254
            NK=NK+1
            IF (WW(I,NC).GT.2.AND.WW(I,NC).LT.8) NK=NK+2
C              FORM ALPHA,BETA,T
  254       DO 255 K=1,NK
              NPAR=NPAR+1
              KPAR(NPAR)=1
              VPAR(NPAR)=0.0
  255       IF (K.EQ.1) VPAR(NPAR)=ALF
            GO TO 260
  256     CONTINUE
  260   CONTINUE
C             SCALE ZETA
        JJ=0                                                                *
        DO 270 I=1,NOSUBC
          II=IX(I)
          IF (II.EQ.0) GO TO 270
          IF (II.LE.JJ) IFLAG=1                                               *
          JJ=II                                                               *
          IF (ZE(II).LT.0.0) GO TO 270
          NPAR=NPAR+1
          KPAR(NPAR)=2
          VPAR(NPAR)=ZE(II)*FACT(2)
  270   CONTINUE
C
        IF (JX.EQ.1) GO TO 291
        JM1=NOSUBC-1
        DO 280 I=1,JM1
          IF (IX(I).EQ.0) GO TO 280
          J1=I+1
          DO 279 J=J1,NOSUBC
            II=IX(I)
            JJ=IX(J)
            IF (JJ.EQ.0) GO TO 279
            IF (II.LT.JJ) GO TO 272
            JJJ=II
            II=JJ
            JJ=JJJ
C              SCALE Fk(Li,Lj)
  272       DO 275 K=1,10
              IF (FGK(II,JJ,K).LT.0.D0) GO TO 275
              NPAR=NPAR+1
              KPAR(NPAR)=3
              VPAR(NPAR)=FGK(II,JJ,K)*FACT(3)
  275       CONTINUE
  279     CONTINUE
  280   CONTINUE
C
        DO 290 I=1,JM1
          IF (IX(I).EQ.0) GO TO 290
          J1=I+1
          DO 289 J=J1,NOSUBC
            II=IX(I)
            JJ=IX(J)
            IF (JJ.EQ.0) GO TO 289
            IF (II.LT.JJ) GO TO 282
            JJJ=II
            II=JJ
            JJ=JJJ
C             SCALE Gk(Li,Lj)
  282       DO 285 K=1,10
              IF (FGK(JJ,II,K).LT.0.D0) GO TO 285
              NPAR=NPAR+1
              KPAR(NPAR)=4
              VPAR(NPAR)=FGK(JJ,II,K)*FACT(4)
  285       CONTINUE
  289     CONTINUE
  290   CONTINUE
C
  291   IF (NPAR.GE.5) GO TO 294
        N1=NPAR+1
        DO 292 I=N1,5
          KPAR(I)=0
          IVPAR(I)=0.D0
  292   VPAR(I)=0.D0
  294   WRITE (9,28) CONF, (VPAR(I), I=1,5), HXID,FACT
   28   FORMAT (//1H0,3A6,6X,5F12.5,5X,A2,5X,5HFACT=,5F6.2)
        IF (NC.GE.NCONT1) VPAR(1)=-9500+IDELE*(NC-NCONT1)
        IF (NC.GT.NCONT2) VPAR(1)=-7500+IDELE*(NC-NCONT2-1)
        DO 295 I=1,NPAR
  295   IVPAR(I)=10000.D0*VPAR(I)
        IF (VPAR(1).GT.99999.99D0) GO TO 296
        IF (VPAR(1).LT.-9999.99D0) GO TO 296
        WRITE (11,29) CONF,NPAR, (IVPAR(I),KPAR(I), I=1,5), HXID,
     1    (IFACT(I), I=1,4)
   29   FORMAT (3A6,I2,5(I9,I1),A2,4I2)
        GO TO 298
  296   WRITE (11,27) CONF,NPAR, (VPAR(I),KPAR(I), I=1,5), HXID,
     1    (IFACT(I), I=1,4)
   27   FORMAT (3A6,I2,F9.1,I1,4(F9.2,I1),A2,4I2)
  298   IF (NPAR.LE.5) GO TO 300
        WRITE ( 9,30) (VPAR(I), I=6,NPAR)
   30   FORMAT (7F12.5)
        WRITE (11,31) (IVPAR(I),KPAR(I), I=6,NPAR)
   31   FORMAT (7(I9,I1))
C
C            CALC CONFIGURATION INTERACTION PARAMETERS
C
  300   IF (NC.LT.NLAST) GO TO 495
        IF (NFIRST.EQ.NLAST) GO TO 350
        IF (ISLI.GT.4) GO TO 350
        SUBR0=0.D0
        P(1,1)=-99.D0
        SUBR=SLIN
        NLM1=NLAST-1
        NERAS=NSCONF(3,KK)
        IF (NERAS.GE.0.OR.NERAS.LE.-8) GO TO 301
        NLM1=MIN0(NLM1,NFIRST-1+MAX0(-NERAS,IPCT))
        NCD1=0
        HA1='   '
        HB1='   '
  301   DO 329 NCD=NFIRST,NLM1
        NCFG(1)=NCD+NC0
        NCPDF=NCD+1
        NCPDL=NLAST
        IF (NCD.NE.NFIRST.AND.NERAS.EQ.-9) NCPDL=NCPDF
        IRHO1=-7
        ISIG1=-7
        DO 329 NCPD=NCPDF,NCPDL
        NCFG(2)=NCPD+NC0
        NRK=0
        DO 328 IRHO=1,NOSUBC
        DO 328 ISIG=IRHO,NOSUBC
        DO 328 IRHOP=1,NOSUBC
        DO 327 ISIGP=IRHOP,NOSUBC
        MMIN=0
        IF (IRHO.NE.IRHOP.OR.ISIG.NE.ISIGP) GO TO 303
        WWW=0.D0
        DO 302 I=ISIG,NOSUBC
  302   WWW=WWW+WW(I,NCD)+WW(I,NCPD)
        IF (WWW.GT.2.D0) GO TO 327
  303   DO 304 I=1,NOSUBC
        WW1=WW(I,NCD)
        WW2=WW(I,NCPD)
        IF (I.EQ.IRHO) WW1=WW1-1.D0
        IF (I.EQ.ISIG) WW1=WW1-1.D0
        IF (I.EQ.IRHOP) WW2=WW2-1.D0
        IF (I.EQ.ISIGP) WW2=WW2-1.D0
        IF (WW1.NE.WW2) GO TO 327
        IF (WW1.LT.0.D0) GO TO 327
        IF (I.LE.ISIG.OR.I.LE.ISIGP.OR.WW1.LE.0.D0) GO TO 304
        IF (LL(I,NCD).NE.LL(I,NCPD)) GO TO 327
        IF (NN(I,NCD).NE.NN(I,NCPD)) MMIN=-1
        IF (NN(I,NCD)+NN(I,NCPD).EQ.198.AND.EE(NCD).NE.EE(NCPD)) 
     1    MMIN=-1
  304   CONTINUE
        A=4*LL(IRHO,NCD)+2
        IF (WW(IRHO,NCD).EQ.A.AND.WW(IRHO,NCPD).EQ.A) GO TO 327
        A=4*LL(ISIG,NCD)+2
        IF (WW(ISIG,NCD).EQ.A.AND.WW(ISIG,NCPD).EQ.A) GO TO 327
        KDMIN=MAX0(IABS(LL(IRHO,NCD)-LL(IRHOP,NCPD)),IABS(LL(ISIG,
     1    NCD)-LL(ISIGP,NCPD)))
        KDMAX=MIN0(LL(IRHO,NCD)+LL(IRHOP,NCPD),LL(ISIG,NCD)
     1    +LL(ISIGP,NCPD))
        IF(KDMIN.GT.KDMAX) GO TO 327
        OVFACT=1.D0
        IF (MMIN.LT.0) GO TO 320
        IF (IOVFACT.LE.0.OR.IOVFACT.EQ.2) GO TO 320
C            CALCULATE OVERLAP FACTOR FOR SPECTATOR ELECTRONS
        HA1T=HA1
        DO 310 I=1,NOSUBC
        WW1=WW(I,NCD)
        IF (I.EQ.IRHO) WW1=WW1-1.D0
        IF (I.EQ.ISIG) WW1=WW1-1.D0
        IF (WW1.LE.0.D0) GO TO 310
        HA=NLB(I,NCD)
        HB=HA
        CALL OVER
        OVFACT=OVFACT*(OVERLP**WW1)
        IF (IOVFACT.GT.3) WRITE (9,1031) NCFG,I,HA,WW1,OVERLP,OVFACT
 1031   FORMAT (3I3,2X,A3,F5.1,2F10.6)
  310   CONTINUE
        HA1=HA1T
        NCD1=-7
  320   CONTINUE
        HA=NLB(IRHO,NCD)
        HB=NLB(ISIG,NCD)
        HC=NLB(IRHOP,NCPD)
        HD=NLB(ISIGP,NCPD)
        N12=2
        IF (NCD.NE.NCD1) N12=1
        IF (HA.NE.HA1.OR.HB.NE.HB1) N12=1
        IF (IRHO.NE.IRHO1.OR.ISIG.NE.ISIG1) N12=1
        CALL SLI2
        NCD1=NCD
        HA1=HA
        HB1=HB
        IRHO1=IRHO
        ISIG1=ISIG
        NCFG2=NCFG(2)
  327   CONTINUE
  328   CONTINUE
        IF (NRK.EQ.0) GO TO 329
        IF (IDIP.EQ.7.OR.IDIP.EQ.8) GO TO 324
        IF (NCD.LT.NCONT1) GO TO 326
  324   DO 325 I=1,NRK
  325   RKAYSC(I)=DEL(NCD)*DEL(NCPD)*RKAYSC(I)
  326   WRITE (9,32) DEL(NCD),DEL(NCPD),(RKAYSC(I),I=1,NRK)
   32   FORMAT (9H SQRTDEL=,2F10.4,5X,5(F11.4,2X)/7(F11.4,2X))
        N2=5
        WRITE (11,33) CONF1(2),CONF1(3),CONF2(2),CONF2(3),NRK,
     1    (RKAYSC(I),N2,I=1,5),HXID,(IFACT(I), I=1,3),IFACT(5)
   33   FORMAT (A6,A3,1H-,A6,A2,I2,5(F9.4,I1),A2,4I2)
        IF (NRK.LE.5) GO TO 329
        WRITE (11,34) (RKAYSC(I),N2, I=6,NRK)
   34   FORMAT (7(F9.4,I1))
  329   CONTINUE
C
  350   IF (IQUAD.LT.3.AND.IQUAD.NE.KK) GO TO 390
        N1=NFIRST
        N2=NLAST
        N3=N1
        N4=N2
        IRMN1=IRNQ+1
        IRMX1=IRXQ+1
        ICN(1)=IC
        ICN(2)=IC
        GO TO 402
  390   IF (KK.EQ.2) GO TO 400
        IF (IW6.LT.0) WRITE (6,38) NFIRST,NLAST,KK
   38   FORMAT (' FINISH G5INP FOR CONFIGURATIONS',I4,'  TO',I4,
     1    ',   PARITY=',I3)
        IF (NC2.EQ.NC3) GO TO 496
        KK=2
        NFIRST=NC2+1
        NLAST=NC3
        IC=32
        WRITE (9,39)
   39   FORMAT (1H1)
        GO TO 495
C
C          CALC QUADRUPOLE OR DIPOLE REDUCED MATRIX ELEMENTS
C
  400   N1=1
        N2=NC2
        N3=NC2+1
        N4=NC3
        IRMN1=IRND+1
        IRMX1=IRXD+1
        ICN(1)=31
        ICN(2)=32
  402   SUBR0=0.D0
        SUBR=DIPN
        IPWR=IRMN1-1
        IF (IDIP.EQ.1) GO TO 460
        J1=ICN(1)-30
        J2=ICN(2)-30
        IF (IGEN.LT.5) GO TO 3402
        ELO=EAV(1)
        EHI=ELO
        DO 1402 I=1,NC3
        EA=EAV(I)
        IF (EA.LT.ELO) ELO=EA
        IF (EA.GT.EHI) EHI=EA
 1402   CONTINUE
        IF (EHI.LE.ELO) EHI=ELO+VPAR(2)
        ENMAX=(EHI-ELO)/109.73731D0
        ENMAX=ENMAX*1.1D0*XMAX
        IF (ENMAX.LT.0.001D0) STOP 140
        EKMAX=2.D0*DSQRT(ENMAX)
        EKMIN=1.0D-3
        EKMAXL=DLOG(EKMAX)
        EKMINL=DLOG(EKMIN)
        DELEKL=(EKMAXL-EKMINL)/DFLOAT(NBIGKS-1)
        WRITE (11,42) NBIGKS,EKMINL,EKMAXL,DELEKL,NENRGS,XMIN,XMAX
   42   FORMAT (I5,3F15.10,I10,2F10.2)
        DO 2402 I=1,NBIGKS
        EKL=EKMINL+DELEKL*(I-1)
 2402   BIGKS(I)=DEXP(EKL)
        ICASE=0
        WRITE (9,43) ENMAX,NBIGKS
   43   FORMAT (/7H ENMAX=,F14.7,8H NBIGKS=,I4,//9H K-VALUES/)
        WRITE (9,44) (BIGKS(I),I=1,NBIGKS)
   44   FORMAT (8F14.7)
 3402   CONTINUE
        DO 450 IPWR1=IRMN1,IRMX1,2
        IPWR=IPWR1-1
        DO 450 NCD=N1,N2
        NCFG(1)=NCD+NC0
        DO 450 NCPD=N3,N4
        NCFG(2)=NCPD+NC0
        DO 439 IRHO=1,NOSUBC
        DO 438 IRHOP=1,NOSUBC
        IF (IABS(LL(IRHO,NCD)-LL(IRHOP,NCPD)).GT.IPWR) GO TO 438
        LSUM=LL(IRHO,NCD)+LL(IRHOP,NCPD)
        IF (LSUM.LT.IPWR) GO TO 438
        IF (MOD(LSUM+IPWR,2).NE.0) GO TO 438
        IF (IRHO.NE.IRHOP) GO TO 410
        IF (NCD.NE.NCPD) GO TO 404
        DO 403 I=IRHO,NOSUBC
        ERAS=4*LL(I,NCD)+2
        WWW=WW(I,NCD)
        IF ((WWW.EQ.0.D0.OR.WWW.EQ.ERAS).AND.I.EQ.IRHO) GO TO 438
        IF (I.EQ.IRHO) GO TO 403
        IF (WWW.EQ.0.D0.OR.WWW.EQ.ERAS) GO TO 403
        IF (LL(I,NCD).EQ.0) GO TO 403
        GO TO 438
  403   CONTINUE
        GO TO 410
  404   WWW=0.D0
        DO 405 I=IRHO,NOSUBC
  405   WWW=WWW+WW(I,NCD)+WW(I,NCPD)
        IF (WWW.NE.2.D0) GO TO 438
  410   DO 420 I=1,NOSUBC
        WW1=WW(I,NCD)
        WW2=WW(I,NCPD)
        IF (I.NE.IRHO.AND.I.NE.IRHOP) GO TO 412
        IF (I.EQ.IRHO) WW1=WW1-1.D0
        IF (I.EQ.IRHOP) WW2=WW2-1.D0
        IF (WW1.LT.0.D0) GO TO 438
        GO TO 414
  412   IF (WW1+WW2.EQ.0.D0) GO TO 420
        IF (NLB(I,NCD).NE.NLB(I,NCPD)) GO TO 438
  414   IF (WW1.NE.WW2) GO TO 438
  420   CONTINUE
        OVFACT=1.D0
        IF (IOVFACT.LT.2) GO TO 430
C            CALCULATE OVERLAP FACTOR FOR SPECTATOR ELECTRONS
        HA1T=HA1
        DO 425 I=1,NOSUBC
        WW1=WW(I,NCD)
        IF (I.EQ.IRHO) WW1=WW1-1.D0
        IF (WW1.LE.0.D0) GO TO 425
        HA=NLB(I,NCD)
        HB=HA
        CALL OVER
        OVFACT=OVFACT*(OVERLP**WW1)
        IF (IOVFACT.GT.3) WRITE (9,1031) NCFG,I,HA,WW1,OVERLP,OVFACT
  425   CONTINUE
        HA1=HA1T
  430   CONTINUE
        HA=NLB(IRHO,NCD)
        HB=NLB(IRHOP,NCPD)
        MMIN=0
        CALL DIP
        GO TO 450
  438   CONTINUE
  439   CONTINUE
        OVFACT=1.D0
        MMIN=-7
        CALL DIP
  450   CONTINUE
        IF (IGEN.LT.5) GO TO 460
        DO 455 IBK=1,NBIGKS
        DO 454 IC=1,ICASE
        DO 452 ICC=1,3
        CONF1(ICC)=CST(IC,ICC)
  452   CONF2(ICC)=CST(IC,ICC+3)
        IPWR=ST(IC,1)
        FRAC=ST(IC,2)
        HXID=ST(IC,3)
        IOVER=ST(IC,4)
        IF1=ST(IC,5)
        HA=KST(IC,1)
        HB=KST(IC,2)
        REDDIP=RDIPS(IBK,IC)
        WRITE (11,46) CONF1,CONF2,REDDIP,HA,RJ,IPWR,HB,FRAC,HXID,IOVER,
     1    IF1
   46   FORMAT(3A6,2X,3A6,F12.5,1H(,A3,2H//,A1,I1,2H//,A3,1H),F6.3,A2,
     1    2I4)
  454   CONTINUE
  455   CONTINUE
        IF (IW6.LT.0) WRITE (6,47)
   47   FORMAT(' FINISH CALCULATION OF BESSEL-FUNCTION MATRIX ELEMENTS')
  460   IF (MOD(IPWR,2).EQ.0) GO TO 390
C
  495 CONTINUE
      IF (IW6.LT.0.AND.KK.EQ.2) WRITE (6,38) NFIRST,NLAST,KK
C
  496 IF (IFLAG.EQ.0) GO TO 500                                          
      DO 497 I=1,5                                                       
  497 WRITE (9,40)
   40 FORMAT (11H **********)
      WRITE (9,41)
   41 FORMAT (34H0*WARNING --- ORBITALS REARRANGED*)
C
C
C
  500 CALL SECONDS(T)
      T=T-T0
      WRITE (9,50) T
   50 FORMAT (//20X,10H-99999999.,10X,2HT=,F7.3,4H SEC)
      IF (NOCSET.LT.0) WRITE (11,51)
   51 FORMAT (20X,10H-55555555.,49X,1H )
      WRITE (11,52)
   52 FORMAT (20X,10H-99999999.,49X,1H )
C
C          FOR CONTINUUM RUNS, CALC FANO GAMMA AND Q
C
C 600 IF (NC2.EQ.NC3.OR.NC2.GT.1) GO TO 700
C     IF (NCX.LE.NCN+2) GO TO 700
C     IF (EE(NCX).LT.0.D0) GO TO 700
C     IF (IDIP.EQ.8) GO TO 700
C     IN=NCN+1
C     DO 610 I=IN,NCX
C     J=NOSUBC
C 606 IF (WW(J,I)-1.D0) 608,610,700
C 608 J=J-1
C     GO TO 606
C 610 CONTINUE
C     WRITE (9,61)
C  61 FORMAT (////15X,2HEE,7X,2HN*,8X,3HWT=,6X,
C    1  10HGAMMA/RKSQ,3X,4HQ*RK,2X,12HPI*(QRKSQ)/2,2X,6HPDIPSC,7X,
C    2  13HRKSC(RYD)----/27X,17HRT(0.5*NST3/ZST2) /)
C     PI=3.141592654D0
C     ZST=ION+1
C     DO 660 I=NCN,NCX
C     IF (I.EQ.NCN) GO TO 630
C     IF (EE(I).GE.0.D0) GO TO 620
C     STN(I)=ZST/DSQRT(-EE(I))
C     WT(I)=DSQRT(0.5D0*STN(I))*STN(I)/ZST
C     PDIPSC(I)=WT(I)*PDIPSC(I)
C     DO 615 J=1,12
C 615 RRSC(I,J)=WT(I)*RRSC(I,J)
C 620 RMAX=0.D0
C     JX=0
C     DO 625 J=1,12
C     IF (RRSC(I,J).NE.0.D0) JX=J
C     IF (DABS(RRSC(I,J)).GT.DABS(RMAX)) RMAX=RRSC(I,J)
C 625 CONTINUE
C     JX=MAX0(JX,1)
C     JX=MIN0(JX,4)
C     GAMMA(I)=2.D0*PI*(RMAX**2)
C     IF (PDIPSC(I).NE.0.D0.AND.RMAX.NE.0.D0) GO TO 628
C     Q(I)=999.999D0
C     Q2(I)=999.999D0
C     GO TO 630
C 628 Q(I)=PDIPSC(NCN)/(PI*RMAX*PDIPSC(I))
C     Q2(I)=0.5D0*PI*(Q(I)**2)
C 630 WRITE (9,63) CF2(I),EE(I),STN(I),WT(I),GAMMA(I),Q(I),Q2(I),
C    1  PDIPSC(I), (RRSC(I,J), J=1,JX)
C  63 FORMAT (1X,A6,F12.6,F9.3,F11.5,F11.5,F11.3,F11.3,F10.5,3X,4F12.6)
C 660 CONTINUE
C
C
  700 NC0=NC0+NC3
      IF (NC0.EQ.NCMAX) GO TO 900
      GO TO 102
C
  900 RETURN
      END
      FUNCTION S3J0SQ(FJ1,FJ2,FJ3)                                      S3J00050
C                                                                       03J000 0
C          CALCULATE SQUARE OF 3-J SYMBOL WITH ZERO MAGNETIC QUANTUM NOS.
C                                                                       S3J00040
C                                                                       S3J00060
      IMPLICIT REAL*8 (A-H,O-Z)
      FJ=FJ1+FJ2+FJ3                                                    S3J00070
      A=FJ-FJ1-FJ1                                                      S3J00080
      B=FJ-FJ2-FJ2                                                      S3J00090
      C=FJ-FJ3-FJ3                                                      S3J00100
      S3J0SQ=FCTRL(A)*FCTRL(B)*FCTRL(C)/FCTRL(FJ+1.D0)                  S3J00110
      A=FCTRL(A*0.5D0)*FCTRL(B*0.5D0)*FCTRL(C*0.5D0)/FCTRL(FJ*0.5D0)    S3J00120
      S3J0SQ=S3J0SQ/(A*A)                                               S3J00130
      RETURN                                                            S3J00140
      END                                                               S3J00150
      FUNCTION S3J (FJ1, FJ2, FJ3, FM1, FM2)                            C3150001
C                                                                       C3150002
C          CALCULATE 3-J SYMBOL
C
      IMPLICIT REAL*8 (A-H,O-Z)
      S3J=0.D0                                                          C3150003
      FM3=FM1+FM2                                                       C3150004
      A=FJ2+FJ3+FM1                                                     C3150005
      B=FJ1-FM1                                                         C3150006
      C=-FJ1+FJ2+FJ3                                                    C3150007
      D=FJ3+FM3                                                         C3150008
      E=FJ1-FJ2-FM3                                                     C3150009
      F=FJ1-FJ2+FJ3                                                     C3150010
      G=FJ1+FJ2-FJ3                                                     C3150011
      H=FJ1+FJ2+FJ3+1.D0                                                C3150012
      E2=FCTRL(B)*FCTRL(FJ1+FM1)*FCTRL(FJ2-FM2)*FCTRL(FJ2+FM2)          C3150013
      E1=(FCTRL(C)*FCTRL(F)/FCTRL(H))*FCTRL(G)*FCTRL(D)*FCTRL(FJ3-FM3)  C3150014
      IF (E1.LE.0.D0) GO TO 153
      E1=DSQRT(E2/E1)                                                   C3150015
      IF (E1) 153,153,150                                               C3150016
  150 ZERO=0.D0
      I1=DMAX1(ZERO,-E)                                                 C3   017
      I2=DMIN1(A,C,D)                                                   C3150018
      IF (I2-I1) 153,151,151                                            C3150019
  151 DO 152 I=I1,I2                                                    C3150020
      FI=I                                                              C3150021
      E2=FCTRL(FI)*FCTRL(C-FI)*FCTRL(D-FI)/FCTRL(A-FI)                  C3150022
  152 S3J=S3J+(((-1.D0)** MOD (I,2))/E2)*FCTRL(B+FI)/FCTRL(E+FI)        C3150023
      IERAS=FJ1+FM2+FM3
      S3J=S3J*((-1.D0)**MOD(IERAS,2))/E1                                C3150024
  153 RETURN                                                            C3150025
      END                                                               C3150026
      FUNCTION FCTRL(A)                                                 FCTR0030
C                                                                       FCTR0020
C          CALCULATE FACTORIAL OF A
C                                                                       FCTR0040
      IMPLICIT REAL*8 (A-H,O-Z)
      FCTRL=1.D0                                                        FCTR0050
      IF (DABS(A)-0.1D0) 10, 10, 11                                     FCTR0060
   10 A=0.D0                                                            FCTR0070
   11 IF (A) 12, 15, 13                                                 FCTR0080
   12 FCTRL=0.D0                                                        FCTR0090
      GO TO 15                                                          FCTR0100
   13 IMAX=A+0.1D0                                                      FCTR0110
      DO 14 I=1,IMAX                                                    FCTR0120
   14 FCTRL=FCTRL*DFLOAT(I)                                             FCTR0130
   15 RETURN                                                            FCTR0140
      END                                                               FCTR0150
      SUBROUTINE SBESS(RB,ITEST,LS,EKSQ,SZ,SZP)
C
C     ROUTINE FOR THE EVALUATION OF SPHERICAL BESSEL FUNCTIONS
C           USING RECURRENCE RELATIONS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /BLK2/ IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6,LSWT,LDEL,LMIN
      COMMON/INFORM/ IREAD,IWRITE,IPUNCH
      DIMENSION SB(200), S(200)
C
      IF (ITEST.EQ.3) EKSQ=-EKSQ
      AK=DSQRT(EKSQ)
      R=AK*RB
      R_1=1.D0/R
      R_2=R_1*R_1
      IF (LS-R) 10,20,20
   10 GO TO (30,60,90), ITEST
   20 GO TO (180,60,90), ITEST
   30 SB(1)=DSIN(R)
      SB(2)=DSIN(R)*R_1-DCOS(R)
      IF (LS-1) 40,50,140
   40 SZ=SB(1)
      SZP=AK*DCOS(R)
      GO TO 210
   50 SZ=SB(2)
      SZP=AK*(DCOS(R)*R_1-DSIN(R)*R_2+DSIN(R))
      GO TO 210
   60 SB(1)=-DCOS(R)
      SB(2)=-DCOS(R)*R_1-DSIN(R)
      IF (LS-1) 70,80,140
   70 SZ=SB(1)
      SZP=AK*DSIN(R)
      GO TO 210
   80 SZ=SB(2)
      SZP=AK*(DSIN(R)*R_1+DCOS(R)*R_2-DCOS(R))
      GO TO 210
   90 SB(1)=1.D0
      SB(2)=-(1.D0+R_1)
      IF (LS-1) 100,110,120
  100 SZ=DEXP(-R)
      SZP=-AK*DEXP(-R)
      GO TO 210
  110 SZ=-(1.D0+R_1)*DEXP(-R)
      SZP=(1.D0+R_1+R_2)*DEXP(-R)
      GO TO 210
  120 LS1=LS+1
      ERAS=-1.D0
      DO 130 L=2,LS1
      ERAS=-ERAS
  130 SB(L+1)=(2*L-1)*ERAS*SB(L)*R_1-SB(L-1)
      L=LS1
      SZ=SB(L+1)*DEXP(-R)
      SZP=AK*((L+1)*(-1)**L*SB(L+1)*R_1+SB(L))*DEXP(-R)
      GO TO 210
  140 LS1=LS+1
      DO 150 L=2,LS1
        SB(L+1)=(2*L-1)*SB(L)*R_1-SB(L-1)
  150 CONTINUE
      IF (ITEST-2) 200,160,200
  160 IF (LS-LS/2*2) 170,200,170
  170 SZ=-SZ
      SZP=-SZP
      GO TO 200
  180 M=LS+12
      IF (R.LT.0.05D0) M=LS+6
      S(M+1)=0.D0
      S(M)=1.0D-20
      MM1=M-1
      DO 190 I=1,MM1
      J=M+1-I
      S(J-1)=(2*J-1)*R_1*S(J)-S(J+1)
  190 CONTINUE
      SB(M)=DSIN(R)/S(1)
      SB(LS)=S(LS)*SB(M)
      SB(LS+1)=S(LS+1)*SB(M)
      SB(LS+2)=S(LS+2)*SB(M)
  200 SZ=SB(LS+1)
      SZP=AK*(SB(LS+1)*R_1+(LS*SB(LS)-(LS+1)*SB(LS+2))/(2*LS+1))
  210 IF (IBUG5-1) 240,220,240
  220 WRITE(IWRITE,230)SZ,SZP
  230 FORMAT (//11H SZ,SZP  = ,(2E16.8)//)
  240 CONTINUE
      RETURN
      END
      SUBROUTINE SECONDS(T)
C
C------------------------------------
C          USE THIS CODE FOR A VAX OR MACINTOSH
      REAL*8 T
      REAL*4 T1
      T1=SECNDS(0.0)
      T=DBLE(T1)
C------------------------------------
C          USE THIS CODE FOR A CRAY
C     T=SECOND(T)
C------------------------------------
C          USE THIS CODE FOR A SUN
C     REAL ETIME,T1
C     REAL*8 T
C     DIMENSION T1(2)
C     T=ETIME(T1)
C-------------------------------------
C          USE THIS CODE FOR IBM RISC
C     REAL*8 T
C     IT=MCLOCK()
C     T=IT
C     T=T/100.0
C-------------------------------------
C          USE THIS CODE FOR IBM-PC LAHEY FORTRAN COMPILER
C     REAL*8 T,TINIT,T1
C     COMMON/TIMING/TINIT,ITINIT
C     CHARACTER TIMESTR*11
C     INTEGER HH,MM,SS,HS
C     CALL TIME(TIMESTR)
C     READ (TIMESTR,10) HH,MM,SS,HS
C  10 FORMAT (3(I2,1X),I2)
C     T1=DFLOAT(3600*HH+60*MM+SS)+HS*1.D-2
C     IF (ITINIT.EQ.0) THEN
C       TINIT=T1
C       ITINIT=1
C     ENDIF
C     T=T1-TINIT
C-------------------------------------
      RETURN
      END
