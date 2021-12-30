      PROGRAM RCE20K   ! Version designed to work with Gerrit van het Hof's
C                      ! TROUT program (TROUTK version).
C     PROGRAM RCE20K(INPUT=OUTGINE=TAPE10=INPUT,OUTPUT=TAPE9=OUTE,
C    1  TAPE2=TAPE2E,TAPE11=PARVALS,TAPE8=LEVELS1,
C    2  TAPE7=LEVELS2,TAPE12=LEVELS3)
C
C     ---------------------------------------------------------------
C          THIS PROGRAM COPYRIGHT BY LOS ALAMOS NATIONAL LABORATORY
C               AND THE U.S. GOVERNMENT.  THE PROGRAM MAY BE FREELY
C               DISTRIBUTED GRATIS ANYWHERE IN THE WORLD, BUT
C               MAY NOT BE SOLD.
C     ---------------------------------------------------------------
C
C     RCE MOD 19  NOV, 1971  (C-SET INPUT VIA TAPE FROM RCG MODS 5 TO 9)
C               MOD 19 IS THE SAME AS MOD 18 EXCEPT WITH LSQM=+1 REMOVED
C               AND WITH RECODING TO REDUCE ARRAY SIZES
C               MOD 20 IS DATED MARCH 1990
C
C
C     PROGRAM TO VARY AS MANY AS IB PARAMETERS TO MAKE
C          LEAST-SQUARES FIT OF ATOMIC ENERGY LEVELS FOR UP TO
C          IC CONFIGURATIONS INVOLVING NO MORE THAN IX J-MATRICES,
C          EACH UP TO IA-BY-IA IN SIZE; SUM OF ALL MATRIX SIZES
C          (TOTAL NUMBER OF ENERGY LEVELS) NO GREATER THAN IV.
C
C          NOTE---REQUIRED STORAGE FOR D(IB,IB+1) SPILLS OVER
C               INTO V(IV,IB).   DO NOT REORDER D AND V IN MEMORY.
C               IV MUST BE .GE. IB.
C
C          THE VERSION FOR A VAX 8200 MADE AT LUND UNIVERSITY,
C               MARCH 1988
C
C               THE VERSION FOR IBM-PC 386 OR 486 MADE AT THE
C          INSTITUTE FOR SPECTROSCOPY RUSSIAN ACADEMY OF SCIENSES
C                     BY A.E.KRAMIDA, NOVEMBER 1993
C          MODIFICATIONS MADE: 
C               - BINARY TAPE 2 TAKEN FROM RCG IS IN COMPRESSED
C                 FORM; 
C               - FOR ITERATIONS, MATRIXES ARE NOT SWAPPED TO DISK,
C                 BUT TO MEMORY BUFFERS IN COMPRESSED FORM (SEE
C                 SUBROUTINES WR1 AND RD1). BUFFERS SIZES ARE MAXB1 
C                 AND MAXB2.
C               - FOR EIGENVECTOR REORDERING AFTER ITERATIONS, A NEW
C                 ALGORITHM IS IMPLEMENTED (SEE SUBROUTINE ORDER)
C               - ALTERNATE INPUT FILE IN SO-CALLED VECTOR FORMAT 
C                 CAN BE USED INSTEAD OF THE OLD STANDARD 'OUTGINE' FILE.
C                 NAME OF THE FILE IS 'RCEINP'; IT IS CREATED BY THIS
C                 RCE PROGRAM IN THE FIRST RUN, AND THEN MAY BE EDITED.
C
C          THIS VERSION 20K IS MADE BY A.KRAMIDA AT NIST, USA,
C              SEPTEMBER 1994 - JANUARY 1995
C
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D(JB,1),V(IV,JB),B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB)
      EQUIVALENCE (D,DELX)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW(IA,IX),EV(IA,3),T(IA,IX),EVAL(IA,IX),W(IA,IX),
     1  FJ(IX),KW(IA,IX),NF(IA,IX),NC(IA,IX),N3(2),S
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,2,IX)*8
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
      COMMON/CMLC1/C(IA,IA)
      COMMON/CMLC2/A(IA,IA)
C     COMMON/CMLC3/Y(IA,IA)
C     COMMON/CMLC4/EVEC(IA,IA)
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC(MAXB1),BUFA(MAXB2),BUFY(MAXBY),
     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2
      PARAMETER (MAXC=10)
      INTEGER*2 NW,NW0
      INTEGER*2 IEVV,IEVV0
      INTEGER*1 II1,II2
      COMMON/VECT/ IEVV(IA,MAXC,IX),IEVV0(IA,MAXC,IX),
     1  NW(IA,MAXC,IX),NW0(IA,MAXC,IX)
      CHARACTER C1*8,Y1*12
      DIMENSION C1(MAXC),Y1(MAXC)
      COMMON/TIMING/TINIT,ITINIT
      DIMENSION KCPL(2)
C
      ITINIT=0                      ! INITIALIZE IBM-PC TIMING
C
C     CALL ERRSET(208,256,-1,1,0,0)
C     CALL ERRSET(212,256 -1,1,0,0)
C
C
C          READ PRINT, REORDER, AND END-OF-RUN CONTROL NUMBERS
C
      OPEN(10,FILE='OUTGINE',STATUS='OLD')
      OPEN(9,FILE='OUTE',STATUS='UNKNOWN')
      OPEN(2,FILE='TAPE2E',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(11,FILE='PARVALS',STATUS='UNKNOWN')
      OPEN(8,FILE='LEVELS1',STATUS='UNKNOWN')
      OPEN(7,FILE='LEVELS2',STATUS='UNKNOWN')
      OPEN(12,FILE='LEVELS3',STATUS='UNKNOWN')
      REWIND 2
      REWIND 8
      OPEN(21,FILE='RCEOUT',STATUS='UNKNOWN')
      IDENOM=0
      NOCSET=-999
      K3=0
      IW=9
      GO TO 94
C
   95 CALL SORT8TO7
      CALL CLOCK(TIME)
      DELT=TIME-TIME0
      WRITE (IW,7) DELT,TIME
      STOP '(NORMAL EXIT)'
C
   94 READ (10,8,END=95) NOCYPR,NOCYCR,NOCYCE,IW6,INV1,IORD,ISAMPL
    8 FORMAT (8I5)
      IF (IW6.GT.0) IW=6           ! DIRECT OUTPUT TO THE SCREEN
C ******************************************************* A.KRAMIDA v
      IF (NOCSET.EQ.-999) THEN
        INVECT=INV1
        IF (INVECT.EQ.0) THEN
C
C         PREPARE NEW VECTOR-INPUT FILE FOR SUBSEQUENT FITTING
C
          OPEN(20,FILE='RCEINP',STATUS='UNKNOWN')
          REWIND(20)
        ELSE
C
C         OPEN THE EXISTING VECTOR-INPUT FILE FOR THE FITTING
C
          OPEN(20,FILE='RCEINP',STATUS='OLD')
        ENDIF
      ENDIF
C ******************************************************* A.KRAMIDA ^
      IF (NOCYPR.EQ.-1.AND.NOCYCR.EQ.0.AND.NOCYCE.EQ.0) GO TO 94
      IF (I216.NE.2) I216=16
      K3=K3+1
C
C     LOCATE THE DESIRED CSET ON UNIT 2
C
      IDN=IDENOM
      IF (INVECT.EQ.0) THEN  ! A.KRAMIDA
        READ (10,10) NCDES,IDENOM,IPRNV,IPRNA,WMAX,IPRNX,IPRNSQ,S ! R.COWAN
        WRITE (20,10) NCDES,IDENOM,IPRNV,IPRNA,WMAX,IPRNX,IPRNSQ,S ! A.KRAMIDA
      ELSE                   ! A.KRAMIDA
        READ (20,10) NCDES,IDENOM,IPRNV,IPRNA,WMAX,IPRNX,IPRNSQ,S !A.KRAMIDA
      ENDIF
      WRITE (21,10) NCDES,IDENOM,IPRNV,IPRNA,WMAX,IPRNX,IPRNSQ,S ! A.KRAMIDA
   10 FORMAT(4I5,F10.4,2I5,F10.4)
      IF (WMAX.EQ.0.D0) WMAX=1.0D10
      IF (S.EQ.0.D0) S=1.D0
      IDN1=0
      IF (NCDES.LT.0) THEN
        GO TO 95
      ELSE
        IF (NCDES.GT.0) THEN
  949     IF (NOCSET.GT.NCDES) THEN
            REWIND 2
          ELSE
            IF (NOCSET.EQ.NCDES) NCDES=0
          ENDIF
C
          IF (NCDES.NE.0) THEN
   96       READ (2) NOCSET,MMAX,LMAX,NMAX,
     *        ((NPAR(L,M),M=L,NMAX),L=1,NMAX)
            IF (NOCSET.LT.NCDES) THEN
              LMAX = LMAX-NMAX+1
              DO 961 M=1,MMAX
                READ (2) N
                READ (2) N       ! ADDITIONAL READ BECAUSE OF USE OF PK1
                DO 961 L=1,LMAX
  961         READ (2) F
              GO TO 96
            ELSE
              IF (NOCSET.GT.NCDES) GO TO 949
            ENDIF
          ENDIF
        ELSE
          IF (IDENOM.LE.0) THEN
            IF (IDENOM.EQ.0) IDENOM=IDN
            IDN1=1
          ENDIF
        ENDIF
      ENDIF
      IF (IDENOM.NE.0.AND.IDN1.EQ.0.AND.INVECT.EQ.0)
     1  READ (10,14) (DENOM(L), L=1,LMAX)
      N3(K3)=NMAX
C
      DO 1098 M=1,MMAX
        DO 1098 K=1,IA
 1098 KW(K,M)=0
      CALL CLOCK(TIME0)
      TIME=0.D0
      WRITE (IW,7) TIME,TIME0
    7 FORMAT (/6H0TIME=,F9.4,19H MIN,   TOTAL TIME=,F9.4,
     1  F14.0,16H LEVELS FIT WITH,F5.0,11H PARAMETERS)
      IF (IW6.LT.0) WRITE (6,7) TIME,TIME0
C
C     REWIND 4
C     REWIND 15
      DO 100 I=1,IC
  100 CONFIG(I,K3)=' '
      IF (INVECT.EQ.0) THEN                        ! A.KRAMIDA
        READ (10,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)   ! R.COWAN
        WRITE (20,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)   ! A.KRAMIDA
      ELSE
        READ (20,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)   ! A.KRAMIDA
      ENDIF
      WRITE (21,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)   ! A.KRAMIDA
   12 FORMAT (7(A8,2X))
      WRITE (11,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)
C     I=2*LMAX
      WRITE (*,9710) K3
 9710 FORMAT (14H Conf. set No ,I1)            ! A.KRAMIDA
      IF (INVECT.EQ.0) THEN                    ! A.KRAMIDA
        READ (10,13) (XIDEN(L), L=1,LMAX)      ! R.COWAN
      ENDIF
   13 FORMAT (7A10)
C9714 CONTINUE
  971 TNO=0.D0
      IMSUM=0
      IMMAX=0
      JLMAX=0
      KKLM=0
C ****************************************** A.KRAMIDA v
C               INITIALIZE BUFFER INDEX OFFSETS
      IAOFFS=0
      ICOFFS=0
      ICOFS1=0
      IAOFS1=0
      IA2=1
C ****************************************** A.KRAMIDA ^
C
      DO 105 M=1,MMAX
        IF (NCDES.NE.0) THEN
C
C      TRANSFER DATA FROM TAPE 2 TO DISC 4.  USE A FOR TRANSFORMATION
C      MATRICES HERE AND IN 'ORDER'.
C
C         READ (2) FJ(M),IM,((A(L,LP),L=1,IM),LP=1,IM),     ! R.COWAN
C    1      (NC(L,M),L=1,IM), ((YIDEN(L,N,M),L=1,IM),N=1,2) ! R.COWAN
C         IMSQ=IM*IM                                        ! R.COWAN
C  *********************************************************  A.KRAMIDA v
          READ (2) KCPL(1),KCPL(2),FJ(M),IM,
     1      (NC(L,M),L=1,IM), ((YIDEN(L,N,M),L=1,IM),N=1,2)
          IMSQ=IM*IM
          CALL UNPK2(2,IM,A,IA,BUFA,IMSQ,IAOFS1,MAXB2)
*********************************************************  A.KRAMIDA ^
          IMAX(M) = IM
          IMSUM=IMSUM+IM
          IMMAX=MAX0(IMMAX,IM)
          IMSQ1=(IMSQ+IM)/2
          N=1
          K=0
C         IFS=4
          LMAXP1=LMAX+1
          DO 3145 L1=1,LMAXP1
            K=K+1
            IF ((K.GT.1).OR.(N.GT.NMAX)) THEN
C             READ (2) ((C(I,J),I=J,IM), J=1,IM)           ! R.COWAN
C  **************************************************
              CALL UNPK2(2,IM,C,IA,BUFC,IMSQ,ICOFS1,MAXB1) ! A.KRAMIDA
C  **************************************************
              IF (L1.EQ.LMAXP1) THEN
                DO 3126 L=1,IM
 3126           GKK(L,M)=C(L,L)
C               IFS=15
              ENDIF
            ELSE
              DO 3120 I=1,IM
                DO 3120 J=1,I
                  C(I,J)=0.D0
                  IF (J.EQ.I.AND.NC(I,M).EQ.N) C(I,J)=1.D0
 3120         CONTINUE
            ENDIF
C           CALL WR(IFS,IM,C,IA,Y,IMSQ1)
            CALL WR1(M,L1,IM,C,IA,BUFC,IBC,IMSQ1,
     *               ICOFFS,MAXB1,MMAX,LMAXP1)
            ICOFS1=ICOFFS
            IF (K.GE.NPAR(N,N)) THEN
              N=N+1
              K=0
            ENDIF
 3145     CONTINUE
C         CALL WR(5,IM,A,IA,Y,IMSQ)
          CALL WR1(M,IA2,IM,A,IA,BUFA,IBA,IMSQ,
     *             IAOFFS,MAXB2,MMAX,IA2)
          IAOFS1=IAOFFS
        ENDIF
C
 1000   IM=IMAX(M)
        KKLM=MAX0(KKLM,IM)
        IF (INVECT.EQ.0) THEN       ! A.KRAMIDA
          KX=0
 3190     KN=KX+1
          KX=KX+7
          IF (KX.GT.IM) KX=IM
          READ (10,1014) (W(K,M), K=KN,KX), TEMP
 1014     FORMAT (8F10.5)
          IF (TEMP.LE.0.D0.AND.KX.LT.IM) GO TO 3190
          READ (10,15) (NF(K,M), K=1,KX)
   15     FORMAT (7I10)
          IF (KX.LT.IM) THEN
            KN=KX+1
            DO 3192 K=KN,IM
              W(K,M)=TEMP
 3192       NF(K,M)=-K
          ENDIF
C ****************************************** A.KRAMIDA v
        ELSE
          IF (M.EQ.1)
     1      READ(20,19)F0,FM,F1,DELF,NOFMAX,GCOE,DXMAX,CRIT1,CRIT,ICRIT2
C
C           READ THE EIGENVECTORS FROM THE VECTOR-INPUT FILE
C
          MAXI=MIN0(MAXC,IM)
          DO 3300 K=1,IM
            DO 3194 I=1,MAXI
              NW(K,I,M)=0 ! -128
 3194       IEVV(K,I,M)=0
            READ (20,3193) T(K,M),ERAS(1),W(K,M),J1,J2,
     1        (IEVV(K,I,M),C1(I),Y1(I), I=1,MAXI)
 3193       FORMAT (F13.6,A1,F10.3,I4,1X,I1,10(I4,1X,A6,1X,A12))
C  16   FORMAT (F13.6,A1,F10.3,I4,A2,10(I4,1H ,A6,1X,A12))
            DO I=1,MAXI
              IEVV(K,I,M)=IEVV(K,I,M)*100
            ENDDO
            NF(K,M)=K
            IF (ERAS(1)(1:1).EQ.'*') THEN
              NF(K,M)=-K
            ELSE
              W(K,M)=T(K,M)
            ENDIF
            DO 3250 I1=1,MAXI
              IF (IEVV(K,I1,M).EQ.0) GO TO 3250
              I2=0
              DO 3195 I=1,IM
                IF (YIDEN(I,1,M).NE.Y1(I1)(1:8)) GO TO 3195
                IF (CONFIG(NC(I,M),K3)(1:6).NE.C1(I1)(1:6)) GO TO 3195
                I2=I ! -128
                GO TO 3200
 3195         CONTINUE
 3200         NW(K,I1,M)=I2
 3250       CONTINUE
 3300     CONTINUE
        ENDIF
C ****************************************** A.KRAMIDA ^
C
        IFIXL=0
        DO 105 K=1,IM
          NFKM=NF(K,M)
          IF (NFKM.LT.0) IFIXL=IFIXL+1
          I=IABS(NFKM)
          T(K,M)=W(I,M)
          TW(K,M)=0.D0
          IF (NFKM.GE.0) THEN
            IF (NFKM.GT.0) NF(K,M)=0
            TW(K,M)=T(K,M)
            TNO=TNO+1.D0
          ENDIF
  105 CONTINUE
      WRITE (IW,16) NMAX,IC,LMAX,IB,MMAX,IX,IMMAX,IA,IMSUM,IV,ICOFFS,
     1   MAXB1,IAOFFS,MAXB2
   16 FORMAT (1H ,22X,'SIZE',2X,'DIMENSION'/
     1  ' NO. OF CONFIGURATIONS=',I3,I9/
     2  ' NO. OF PARAMETERS=',I7,I9/
     3  ' NO. OF J-VALUES=',I9,I9/
     4  ' LARGEST J-MATRIX=',I8,I9/
     5  ' TOTAL NO. OF LEVELS=',I5,I9/
     7  ' C-BUFFER SIZE=',I11,I9/
     8  ' A-BUFFER SIZE=',I11,I9)
 1104 IF (NMAX.GT.IC.OR.LMAX.GT.IB.OR.MMAX.GT.IX)
     1  STOP '(DIMENSION ERROR: IC, IB OR IX TOO SMALL)'
      IF (IMMAX.GT.IA.OR.IMSUM.GT.IV) 
     1  STOP '(DIMENSION ERROR: IA OR IV TOO SMALL)'
C
C       NF(K,M) - ARRAY OF LEVEL FLAGS FOR MOMENTUM J=M:
C         NF(K,M)=0 FOR KNOWN LEVELS
C         NF(K,M)=-K FOR UNKNOWN LEVELS
C       T(K,M)  - ARRAY OF INPUT LEVEL ENERGIES
C       W(K,M) - ARRAY OF CALCULATED LEVELS
C       TW(K,M) - DIFFERENCES T-W FOR KNOWN LEVELS, 0 FOR UNKNOWN
C       LF(L) - ARRAY OF PARAMETER FLAGS
C
      IF (INVECT.EQ.0) THEN             ! A.KRAMIDA
        READ (10,15) (LF(L), L=1,LMAX)
        READ (10,14) (XMAX(L), L=1,LMAX)
   14   FORMAT (7F10.5)
        READ (10,14) (X(L), L=1,LMAX)
      ELSE
C ****************************************** A.KRAMIDA v
        READ (20,3411) L
 3411   FORMAT (A1)
        IF (IDENOM.EQ.0) THEN
          DO 3400 L=1,LMAX
 3400     READ (20, 3410) XIDEN(L),LF(L),X(L),XMAX(L)
 3410     FORMAT (A10,I5,2F14.4)
        ELSE
          DO 3420 L=1,LMAX
 3420     READ (20, 3430) XIDEN(L),LF(L),X(L),XMAX(L),
     1      IDENOM
 3430     FORMAT (A10,3F14.4)
        ENDIF
      ENDIF
      I=0
      DO 3450 L=1,LMAX
        IF (XIDEN(L)(1:3).EQ.'EAV') I=I+1
        IF (XIDEN(L)(1:1).EQ.'1') I=0
        IF (I.NE.0) THEN
          WRITE (NX(L),3444) I
 3444     FORMAT (I2)
        ELSE
          WRITE (NX(L),3445) '  '
 3445     FORMAT (A2)
        ENDIF
 3450 CONTINUE
      IF (LMAX.EQ.1.AND.X(1).EQ.0.D0) X(1)=1.D-6  ! TO AVOID 1/0 IN EIGEN
C ****************************************** A.KRAMIDA ^
      NOCYC=0
      LSQM=0
      JLMAX=LMAX
      SUMXSQ=0.D0
      IFIXED=0
      DO 116 L=1,LMAX
        LF0(L)=LF(L)
        XMAX0(L)=XMAX(L)
        SDX(L) = 0.D0
        X0(L)=0.D0
        DELXC(L)=0.D0
        DELX(L)=0.D0
        IF (XMAX(L).LE.0) THEN
          XMAX(L)=1.0D10
        ELSE
          XM=XMAX(L)            ! CONVERT REAL*4 TO REAL*8
          XMAX(L)=DMAX1(XM,X(L))
        ENDIF
        ISUM=0
        IF (LF(L).LT.-50) THEN
          XMIN(L)=-XMAX(L)
          IF (LF(L).GT.-99) THEN
            LF(L)=LF(L)+50
            LFL=LF(L)           ! CONVERT INTEREG*2 TO INTEGER*4
            JLMAX=MIN0(JLMAX,LFL)
            IFIXED=IFIXED+1
            ISUM=1
          ELSE
            IF (LF(L).EQ.-99) THEN
              LF(L)=0
              IFIXED=IFIXED+1
              ISUM=1
            ELSE
              GO TO 116
            ENDIF
          ENDIF
        ENDIF
        IF (ISUM.EQ.0) THEN
          XMIN(L)=DMAX1(1.0D-6, X(L)-XMAX(L)+X(L) )
          IF (LF(L).GE.100) GO TO 116
          IF (LF(L).LT.0) THEN
            LFL=LF(L)           ! CONVERT INTEREG*2 TO INTEGER*4
            JLMAX=MIN0(JLMAX,LFL)
          ELSE
            LF(L)=0
          ENDIF
          ISUM=1
          IFIXED=IFIXED+1
        ENDIF
        IF (ISUM.NE.0) SUMXSQ=SUMXSQ+X(L)*X(L)
  116 CONTINUE
C**********************************************    A.KRAMIDA v
C           IF NO PARAMETERS ARE VARIED,
C           OR NO LEVELS ARE TO FIT, OR
C           LEVELS ARE TO BE FITTED WITH TOO MANY FREE PARAMETERS,
C           OR IT SEEMS THAT NOTHING WAS CHANGED IN OUTGINE FILE
C           AFTER RUNNING RCG,
C           THEN EXIT AFTER THE 0-TH CYCLE
C        IFIXED = NUMBER OF VARIED PARAMETERS
C        IFIXL  = NUMBER OF UNKNOWN LEVELS (NOT "EXPERIMENTAL")
      IFIRST=0
      IF (IFIXED.GT.IMSUM-IFIXL) IFIRST=1
      IF (IFIXED.EQ.LMAX.AND.IFIXL.EQ.0.AND.
     1  (NOCYCE .EQ. 5 .OR. NOCYCE .EQ. 0) ) IFIRST=1
      IF (IFIXED.EQ.0.OR.IFIXL.EQ.IMSUM) IFIRST=1
      IF (IFIRST.NE.0) NOCYCE=0
C**********************************************    A.KRAMIDA ^
      MIK=0
      IF (SUMXSQ.EQ.0.D0) THEN
C       REWIND I216
C       DO 1115 I=1,KKLM
C         DO 1114 K=1,KKLM
C1114     Y(I,K)=0.D0
C1115   Y(I,I)=1.D0
C       DO 1116 M=1,MMAX
C         IM=IMAX(M)
C         IMSQ=IM*IM
C         CALL WR(I216,IM,Y,IA,C,IMSQ)
C1116   CONTINUE
        DO 1116 M=1,MMAX
          IM=IMAX(M)
          DO K=1,IM
            DO I=1,IM
              MIK=MIK+1
              IF (MIK.LE.MAXBY) THEN
                IF (I.NE.K) THEN
                  BUFY(MIK)=0.D0
                ELSE
                  BUFY(MIK)=1.D0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
 1116   CONTINUE
      ELSE
        DO 1118 M=1,MMAX
          IM=IMAX(M)
          DO 1118 K=1,IM
          DO 1118 I=1,IM
            MIK=MIK+1
 1118   CONTINUE
      ENDIF
      WRITE (IW,1117) MIK,MAXBY
 1117 FORMAT (' Y-BUFFER SIZE=',I11,I9)
      IF (MIK.GT.MAXBY) THEN
        WRITE (*,1119) MIK,MAXBY
 1119   FORMAT (15H Y buffer size=,I10, 12H, dimension=,I10)
        STOP '(DIMENSION OF Y-BUFFER EXCEEDED.)'
      ENDIF
C
C
C          READ ITERATION CONTROL CARD AND FIND STARTING POINT
C
  120 IF (INVECT.EQ.0) THEN    ! A.KRAMIDA
        READ (10,19) F0,FM,F1,DELF,NOFMAX,GCOE,DXMAX,CRIT1,CRIT,ICRIT2 !R.COWAN
   19   FORMAT (3F10.1,F5.2,I5,F10.3,F10.2,2F5.2,I2)            !R.COWAN
C *************************************************** A.KRAMIDA v
        WRITE (20,19) F0,FM,F1,DELF,NOFMAX,GCOE,DXMAX,CRIT1,CRIT,ICRIT2
      ENDIF
      WRITE (21,19) F0,FM,F1,DELF,NOFMAX,GCOE,DXMAX,CRIT1,CRIT,ICRIT2
C *************************************************** A.KRAMIDA ^
 2120 F=F0
      CALL CLOCK(TIME)
      DELT=TIME-TIME0
      WRITE (IW,7) DELT,TIME
C
C        SET UP FOR CALCN WITH LINKED PARAMS
C
      IF (JLMAX.LT.0.AND.SUMXSQ.GT.0) THEN
        F=1.D0
        LSQM=-1
        JLMAX=-JLMAX
        DO L=1,LMAX
          IABSLF=IABS(LF(L))-100
          IF ((IABSLF.LT.0.AND.LF(L).GE.0).OR.
     *        (IABSLF.GE.0.AND.XMAX(L).LE.0.9D10)) THEN
            JLMAX=JLMAX+1
            LF(L)=-JLMAX
          ENDIF
        ENDDO
        PNO=JLMAX
      ENDIF
C
      FF=-1.D0
      G=1.D0
      IF (F.GE.F1) THEN
        FF=0.D0
        IF (F.GE.0.9999D0) THEN
          NOFMAX=1000
          FF=1.D0
          IF (LSQM.GE.0) THEN
            DO L=1,LMAX
              IF (IABS(LF(L)).GE.100) THEN
                IF (XMAX(L).LE.0.9D10) LF(L)=0
              ELSE
                XMAX(L)=1.0D10
                IF (XMIN(L).GE.0) THEN
                  XMIN(L)=1.0D-6
                ELSE
                  XMIN(L)=-1.0D10
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
C
  400 IF (LSQM.GE.0) THEN
        JLMAX=0
        PNO=0.D0
        DO L=1,LMAX
          IF (IABS(LF(L)).LT.100) THEN
            JLMAX=JLMAX+1
            LF(L)=JLMAX
            PNO=PNO+1.D0
          ENDIF
        ENDDO
      ENDIF
C
C     REWIND 3
C
      NOF=0
C     N5=-1
C ****************************************v  A.KRAMIDA
      DO 135 M=1,MMAX
      DO 135 I=1,MAXC
      DO 135 K=1,IMMAX
        IEVV0(K,I,M)=IEVV(K,I,M)
  135 NW0(K,I,M)=NW(K,I,M)
C ****************************************^  A.KRAMIDA
      IF (SUMXSQ.NE.0.OR.NOCYCE.EQ.0) THEN
        SUMXSQ=0.D0
        DETB=0.D0
C       N5=-2
C       IF (NOCYCE.LE.2) N5=1-NOCYCE
C       IF (NOCYCE.EQ.0) N5=10
        ABDXMX=0.D0
        CALL CLOCK(TIME)
        DELT=TIME-TIME0
        WRITE (IW,7) DELT,TIME
        GO TO 180
      ENDIF
C
C
C                    ITERATION CONTROL
C
C
  140 NOCYC=NOCYC+1
      NOF=NOF+1
      CALL CLOCK(TIME)
      DELT=TIME-TIME0
      WRITE (IW,7) DELT,TIME
      IF (NOF.LE.NOFMAX) GO TO 170
      NOF=1
      IF (FF.GT.0) GO TO 170
      IF (FF.LT.0) THEN
        F=F*FM
        IF (F.LT.F1) GO TO 170
        F=F1-DELF
        FF=0.D0
      ENDIF
      IF (F.LT.0.9999D0) THEN
        F=DMIN1((F+DELF),1.0D0)
        IF (F.GT.0.9999D0) THEN
          F=1.D0
          NOFMAX=NOFMAX+4
        ENDIF
        GO TO 170
      ENDIF
      NOFMAX=1000
      FF=1.D0
      DO 162 L=1,LMAX
        IF (IABS(LF(L)).GE.100) THEN
          IF (XMAX(L).LE.0.9D10) THEN
            JLMAX=JLMAX+1
            IF (LSQM.LT.0) THEN
              LF(L)=-JLMAX
            ELSE
              LF(L)=JLMAX
            ENDIF
  158       PNO=PNO+1.D0
          ENDIF
        ELSE
          XMAX(L)=1.0D10
          IF (XMIN(L).LT.0) THEN
            XMIN(L)=-1.0D10
          ELSE
            XMIN(L)=1.0D-6
          ENDIF
        ENDIF
  162 CONTINUE
C
C
C                    CALCULATE NEW PARAMETERS AND EIGENVALUES
C
  170 CALL CALCV                ! READS Y ARRAY FROM Y-BUFFER 
                                ! AND C ARRAY FROM C-BUFFER
      DO 175 M=1,MMAX
      DO 175 I=1,MAXC
      DO 175 K=1,IMMAX
        IEVV0(K,I,M)=IEVV(K,I,M)
  175 NW0(K,I,M)=NW(K,I,M)
      CALL CALCX
C
  180 NSQ=0
      CALL EIGEN              ! WRITES Y ARRAY TO Y-BUFFER 
                              ! READS C ARRAY FROM C-BUFFER
C
      DO L=1,LMAX
        SDX(L)=STDEV*SDX(L)
      ENDDO
C
C
C
C                    OUTPUT
C
      CALL CLOCK(TIME)
      DELT=TIME-TIME0
      WRITE (IW,7) DELT,TIME,TNO,PNO
C
      LEAVE=0                             ! SET END-OF-CYCLE CONDITION
      IF (NOCYC.GE.NOCYCE) THEN
        LEAVE=1
      ELSE
        IF (NOCYC.NE.0.AND.F.GT.0.9999D0) THEN
          ABDXMX=0.D0
          DO L=1,LMAX
            Z=DABS(DELX(L))
            IF (Z.GT.ABDXMX) ABDXMX=Z
          ENDDO
          IF (ABDXMX.LE.3.D-3) LEAVE=1
        ENDIF
      ENDIF
C
C     IF (IPRNX.EQ.0.OR.LEAVE) THEN
C            WRITE NEW AND OLD PARAMETERS TO 'OUTE'
      WRITE (IW,71) ELIDEN, (CONFIG(L,K3),L=1,18),NOCYC,F,GEFF,CRIT,
     1  DETB,AVDEV,STDEV  !, SD1,SD2
   71 FORMAT (1H , 16A8/9X,3A8/ 8H CYC NO=,I2, 4H  F=,F5.3,
     1  4H  G=,F5.3, 7H  CRIT=,F4.2, 8H  DET B=,1PE10.3,
     2  8H  AVDEV=,0PF8.3, 8H  STDEV=,F8.3/)
C     WRITE (IW,18)
C  18 FORMAT (1H0)
C     WRITE (IW,20)
   20 FORMAT (1H )
      WRITE (IW,72) NOCSET
   72 FORMAT (44H SLATER PRMTR    LF      NEW X         DEL X,
     1  39H        DEL XC      OLD X           SDX,
     2  34H         RSS W      DIFF     CSET=,I3)
C     DO L=1,LMAX
C
C
C       WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
C    1    DELXC(L), X0(L),   SDX(L)
C     ENDDO
      WRITE (IW,20)
C
      DO L=1,LMAX
C
C
        IF (L.LE.MMAX) THEN
          DIFF=ANORM(L)-SEVSQ(L)
          WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
     1      DELXC(L), X0(L),   SDX(L), SEVSQ(L), DIFF
   73     FORMAT (1H ,A10,2X,I6, 4F13.4, F13.3, F15.3, F9.3)
        ELSE
          WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
     1      DELXC(L), X0(L),   SDX(L)
        ENDIF
      ENDDO
      IF (IDENOM.NE.0.AND.LEAVE.NE.0) THEN
C
C       CONVERT SUPERSCRIPTED PARAMETERS TO SUBSCRIPTED
C
        WRITE (IW,71) ELIDEN, (CONFIG(L,K3),L=1,15),NOCYC,F,GEFF,
     1    CRIT,DETB,AVDEV,STDEV  !, SD1,SD2
        WRITE (IW,80)
   80   FORMAT (13H  SUBSCRIPTED )
        WRITE (IW,72) NOCSET
        DO L=1,LMAX
          IF (DENOM(L).NE.0.D0) THEN
            X(L) = X(L)/DENOM(L)
            DELX(L) = DELX(L)/DENOM(L)
            DELXC(L) = DELXC(L)/DENOM(L)
            X0(L) = X0(L)/DENOM(L)
            SDX(L) = SDX(L)/DENOM(L)
          ENDIF
C
C
          IF (L.LE.MMAX) THEN
            DIFF=ANORM(L)-SEVSQ(L)
            WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
     1        DELXC(L), X0(L),   SDX(L), SEVSQ(L), DIFF
          ELSE
            WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
     1        DELXC(L), X0(L),   SDX(L)
          ENDIF
        ENDDO
      ENDIF
      IF (IPRNV.EQ.0.AND.NOCYCE.EQ.0) CALL CALCV  ! FOR V-PRINTING
C
      WRITE (11,39) NOCYC,AVDEV    ! WRITE TO 'PARVALS'
   39 FORMAT (32H PARAMETER VALUES,  CYCLE NUMBER,I3,9H,  AVDEV=,F8.3)
      WRITE (11,40) (X(L), L=1,LMAX)
C
C *************************************************** A.KRAMIDA v
      WRITE (*,9711) NOCYC,AVDEV   ! WRITE TO THE SCREEN
9711  FORMAT (14H Iteration No ,I1,21H finished with AVDEV=,0PF8.3)
   40 FORMAT (7F10.4)
C
      IF (LEAVE.EQ.0) GO TO 140
C
C     IF (LEAVE.NE.0.OR.IPRNSQ.NE.2) THEN
C
C       REORDER EIGENVECTORS BY CALCULATED ENERGY
C
        DO M=1,MMAX
          IM=IMAX(M)
          DO K=1,IM
            KW(K,M)=0
          ENDDO
          DO K1=1,IM
            GRE=1.0D30
            DO 7255 N1=1,IM
              K2=K1-1
              DO K=1,K2
                IF (KW(K,M).EQ.N1) GO TO 7255
C                   ! SKIP ALREADY MAPPED LEVELS
              ENDDO
              IF (W(N1,M).LE.GRE) THEN
C                   ! FIND THE LEAST FROM THE REST
                GRE=W(N1,M)
                KW3=N1
              ENDIF
 7255       CONTINUE
            KW(K1,M)=KW3          ! MAP THE K1-TH LEVEL ONTO KW3
          ENDDO
        ENDDO
C
C               WRITE TO 'OUTE' AND 'LEVELS1' AND 'RCEOUT'/'RCEINP'
C
        DO NOUT=1,2
          DO NSQ1=0,1
            IF (NSQ1.EQ.1) THEN
C                   WRITE TO 'LEVELS1', FOR BOTH COUPLING SCHEMES
              WRITE (8,77) ELIDEN, (CONFIG(I,K3),I=1,NMAX)
   77         FORMAT (A8,10(2X,A8)/(8X,10(2X,A8)))
            ENDIF
            IF (IPRNSQ.NE.2) THEN
C
C               PRINT EIGENVECTOR COMPONENTS FOR BOTH COUPLING SCHEMES
C
              IF (NSQ1.EQ.0) THEN
                WRITE (IW,74) NORUN,NOCYC
   74           FORMAT (/1H ,A2,I3,
     1            '  T         W        T-W    GKK    G    CONF',
     2            '     EIGENVECTOR COMPONENTS')
              ELSE
                IF (IPRNSQ.EQ.0) THEN
                  WRITE (IW,75) NORUN,NOCYC
   75             FORMAT (/1H ,A2,I3,
     1              '  T         W        T-W    GKK    G    CONF',
     2              '     EIGENVECTOR COMPONENTS SQUARED')
                ENDIF
              ENDIF
            ENDIF
C           MIK=0
            NYOFS=0
            DO M=1,MMAX
              IM=IMAX(M)
              IMSQ=IM*IM
              ITRANS=0
              IF (NOUT.EQ.2.AND.KCPL(1).EQ.1) ITRANS=1
              IF (NOUT.EQ.1.AND.KCPL(1).EQ.2) ITRANS=1
              IF (KCPL(2).EQ.1.AND.KCPL(1).NE.2) THEN
                IF (NOUT.EQ.1) THEN
                  ITRANS=2
                ELSE
                  ITRANS=0
                ENDIF
              ENDIF
              IF (ITRANS.NE.0.AND.NSQ1.EQ.0) THEN
c             IF (NOUT.EQ.2.AND.NSQ1.EQ.0) THEN
                IMSQ1=(IMSQ+IM)/2
C               CALL RD(I216,IM,Y,IA,EVEC,IMSQ)  ! R.COWAN
C               CALL RD(15,IM,C,IA,EVEC,IMSQ1)   ! R.COWAN
C               CALL RD(15,IM,A,IA,EVEC,IMSQ)    ! R.COWAN
C               DO K=1,IM
C                 DO I=1,IM
C                   MIK=MIK+1
C                   Y(I,K)=BUFY(MIK)
C                 ENDDO
C               ENDDO
                CALL RD1(M,LMAXP1,IM,C,IA,BUFC,IBC,IMSQ1,   ! A.KRAMIDA
     *            ICOFFS,MAXB1,MMAX,LMAXP1)
                CALL RD1(M,IA2,IM,A,IA,BUFA,IBA,IMSQ,       ! A.KRAMIDA
     *            IAOFFS,MAXB2,MMAX,IA2)
                DO K=1,ITRANS
                  CALL ORDOUT2
                ENDDO
              ENDIF
C
C               Find pure-coupling-term labels for levels
              IF (CRIT1.GT.0.AND.NSQ1.EQ.1) CALL ORDOUT1
C
              DO K=1,IM
                J=IABS(KW(K,M))
                IF (W(J,M).LE.WMAX.OR.TW(J,M).NE.0.D0) THEN
C                   Find configuration name for each level 
C                         (by max. contribution)
                  DO KK=1,IC
                    VC(KK)=0.D0
                  ENDDO
                  JOFS=NYOFS+(J-1)*IM
                  DO I=1,IM
                    NCC=NC(I,M)
                    VC(NCC)=VC(NCC)+BUFY(JOFS+I)**2
C                   IF (NSQ1.EQ.0) THEN
C                     VC(NCC)=VC(NCC)+EVEC(I,J)**2
C                   ELSE
C                     VC(NCC)=VC(NCC)+EVEC(I,J)
C                   ENDIF
                  ENDDO
                  GRE=0.D0
                  DO KK=1,IC
                    IF (VC(KK).GT.GRE) THEN
                      GRE=VC(KK)
                      NJM=KK
                    ENDIF
                  ENDDO
C
C                   Print the level to all output files needed
                  ERAS(1)=' '
                  IF (NF(J,M).LT.0.OR.IFIRST.NE.0) ERAS(1)='*'
C                 IF (IPRNSQ.EQ.0.OR.(NSQ1.EQ.0.AND.IPRNSQ.NE.2))
C    1              WRITE (IW,76) T(J,M),ERAS(1), W(J,M), TW(J,M),
C    2                GKK(J,M),GFAC(J,M), CONFIG(NJM,K3),
C    3                (EVEC(I,J), I=1,IM)
C  76               FORMAT (1H ,F10.3,A1,F10.3,F9.3,F7.3,F6.3,1X,
C    1                A6,8F10.3:/(49X,8F10.3:))
C                 IF (LEAVE.NE.0)
                  CALL OUT8(J,NJM,NOUT,NSQ1,LEAVE,IPRNSQ,KCPL(NOUT))
                ENDIF
              ENDDO
              NYOFS=NYOFS+IMSQ
            ENDDO
          ENDDO
          WRITE (IW,78) NOUT
   78     FORMAT (29H Finished output for coupling,I2)
        ENDDO
C     ENDIF
      IF (LEAVE.EQ.0) GO TO 140
C
C          WRITE PARAMETER VALUES FOR RCG INPUT
C
  900 WRITE (11,91)
   91 FORMAT (31H PARAMETER VALUES FOR RCG INPUT)
      L2=0
      DO NAZ=1,NMAX
        L1=L2+1
        L2=L1+NPAR(NAZ,NAZ)-1
        WRITE (11,92) ELIDEN,CONFIG(NAZ,K3), (X(L),L=L1,L2)
   92   FORMAT (A6,A8,6X,F10.4,4(F9.4,1X),10HLS99999999:/
     1    (7(F9.4,1X),10HLS99999999))
      ENDDO
      NXM1=NMAX-1
      IF (NXM1.EQ.0) GO TO 969
      DO N1=1,NXM1
        N1P1=N1+1
        DO N2=N1P1,NMAX
          NP12=NPAR(N1,N2)
          IF (NP12.NE.0) THEN
            L1=L2+1
            L2=L1+NP12-1
            WRITE (11,93) CONFIG(N1,K3),CONFIG(N2,K3), (X(L),L=L1,L2)
   93       FORMAT (A8,1H-,A8,3X,5(F9.4,1X),10HLS99999999:/
     1        (7(F9.4,1X),10HLS99999999))
          ENDIF
        ENDDO
      ENDDO
C *************************************************** A.KRAMIDA v
C
C          WRITE PARAMETER VALUES FOR RCE INPUT
C
  969 IF (INVECT.EQ.0) WRITE (20,3490)
      WRITE (21,3490)
 3490 FORMAT (51H PARAMETER FLAG      VALUE     MAX.VALUE      DENOM)
      LCI=0
      IF (IDENOM.NE.0) THEN
        DO 3500 L=1,LMAX
 3500   X(L)=X(L)*DENOM(L)
      ENDIF
      DO L=1,LMAX
        IF (XIDEN(L)(1:1).EQ.'1') LCI=1
        C1(1)=' '
        C1(2)=' '
        IF (LCI.NE.0) THEN
          READ (XIDEN(L),3505) II1,II2
 3505     FORMAT(2A1)
          IF (II1.EQ.42) THEN
            II1=62
          ELSE
            IF (II1.LE.122.AND.II1.GE.97) THEN
              II1=II1-61
            ELSE
              IF (II1.LE.90.AND.II1.GE.65) THEN
                II1=II1-55
              ELSE
                IF (II1.LE.57.AND.II1.GE.49) II1=II1-48
              ENDIF
            ENDIF
          ENDIF
          IF (II2.EQ.42) THEN
            II2=62
          ELSE
            IF (II2.LE.122.AND.II2.GE.97) THEN
              II2=II2-61
            ELSE
              IF (II2.LE.90.AND.II2.GE.65) THEN
                II2=II2-55
              ELSE
                IF (II2.LE.57.AND.II2.GE.49) II2=II2-48
              ENDIF
            ENDIF
          ENDIF
          C1(1)=CONFIG(II1,K3)
          C1(2)(1:1)='-'
          C1(2)(2:7)=CONFIG(II2,K3)
        ENDIF
        IF (IDENOM.NE.0) THEN
          IF (INVECT.EQ.0) WRITE (20, 3520) XIDEN(L),
     1      LF0(L),X(L),XMAX0(L),DENOM(L),C1(1),C1(2)
          WRITE (21, 3520) XIDEN(L),LF0(L),X(L),
     1      XMAX0(L),DENOM(L),C1(1),C1(2)
 3520     FORMAT (A10,I5,3F14.4,A6,A8)
        ELSE
          IF (INVECT.EQ.0) WRITE (20, 3540) XIDEN(L),
     1      LF0(L),X(L),XMAX0(L),C1(1),C1(2)
          WRITE (21, 3540) XIDEN(L),LF0(L),X(L),
     1      XMAX0(L),C1(1),C1(2)
 3540     FORMAT (A10,I5,2F14.4,14X,A6,A8)
        ENDIF
      ENDDO
C *************************************************** A.KRAMIDA ^
C
      WRITE (11,93)
C
C          SET FOR NEXT VALUE OF LSQM OR GO TO NEXT RUN
C
C 970 IF (LSQM.LT.0) THEN             ! A.KRAMIDA: ICRIT2 IS DISABLED
C       LSQM=LSQM+1                   ! BECAUSE THE EIGENVECTORS SAMPLE
C       IF (LSQM.LE.ICRIT2) GO TO 400 ! IS IN 1ST COUPLING SCHEME BUT
C     ENDIF                           ! CURRENT EIGENVECTORS HAVE BEEN
C                                     ! REARRANGED TO 2ND COUPLING SCHEME
C *************************************************** A.KRAMIDA v
      IF (INVECT.NE.0) THEN
 9709   READ (10,2110) ERAS(1)
 2110   FORMAT (A8)
        IF (ERAS(1).NE.'   -1   ') GO TO 9709
      ENDIF
C *************************************************** A.KRAMIDA ^
      GO TO 94
C
      END
C
      SUBROUTINE UNPK2(LD,IM,A,IA,B,IMSQ,IOFFS,MAXLEN)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*12 B,BB
      INTEGER*2 BB1,I,J
      REAL*8 R8
      DIMENSION A(IA,IA),B(MAXLEN),BB1(6)
      EQUIVALENCE (R8,BB1(1)),(BB,BB1(1)),(I,BB1(5)),(J,BB1(6))
C
      READ (LD) IOFFS2,(B(K),K=IOFFS+1,IOFFS+IOFFS2)
C
      DO 10 J=1,IM
        DO 10 I=1,IM
   10 A(I,J)=0.D0
C
      IM2=IM*IM
      IF (IMSQ.LT.IM2) THEN
        DO K=IOFFS+1,IOFFS+IOFFS2
          BB=B(K)
          A(I,J)=R8
          A(J,I)=R8
        ENDDO
      ELSE
        DO K=IOFFS+1,IOFFS+IOFFS2
          BB=B(K)
          A(I,J)=R8
        ENDDO
      ENDIF
      IOFFS=IOFFS+IOFFS2
  300 RETURN
      END
C
C
      SUBROUTINE OUT8(J,NJM,NOUT,NSQ1,LEAVE,IPRNSQ,KCPL)
C        OUTPUT LEVELS LIST
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D(JB,1),V(IV,JB),B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB)
      EQUIVALENCE (D,DELX)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW(IA,IX),EV(IA,3),T(IA,IX),EVAL(IA,IX),W(IA,IX),
     1  FJ(IX),KW(IA,IX),NF(IA,IX),NC(IA,IX),N3(2),S
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,2,IX)*8
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
      COMMON/CMLC1/C(IA,IA)
      COMMON/CMLC2/A(IA,IA)
C     COMMON/CMLC3/Y(IA,IA)
C     COMMON/CMLC4/EVEC(IA,IA)
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C  *******************************************************v A.KRAMIDA
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC(MAXB1),BUFA(MAXB2),BUFY(MAXBY),
     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2
      PARAMETER (MAXC=10)
      INTEGER*2 NW,NW0
      INTEGER*2 IEVV,IEVV0
      INTEGER*1 II1,II2
      COMMON/VECT/ IEVV(IA,MAXC,IX),IEVV0(IA,MAXC,IX),
     1  NW(IA,MAXC,IX),NW0(IA,MAXC,IX)
C  *******************************************************^ A.KRAMIDA
C
      CHARACTER C1*8,Y1*12,Y11*3,AJ1*2,AJ2*2
      DIMENSION IP1(MAXC),IP2(MAXC),C1(MAXC),Y1(MAXC)
C *********************************************** R.COWAN v
C     DO 120 L=1,3
C       GRE=0.D0
C       DO 110 I=1,IM
C         IF (L.GT.1.AND.I.EQ.I1(1)) GO TO 110
C         IF (L.GT.2.AND.I.EQ.I1(2)) GO TO 110
C         IF (EVEC(I,J).LT.GRE) GO TO 110
C         GRE=EVEC(I,J)
C         I1(L)=I
C 110   CONTINUE
C 120 CONTINUE
C *********************************************** R.COWAN ^
C *********************************************** A.KRAMIDA v
C     IF (LEAVE.EQ.0.AND.NOUT.EQ.2) RETURN
C     IF (LEAVE.EQ.0.AND.IPRNSQ.NE.0.AND.NSQ1.EQ.1) RETURN
C     IF (LEAVE.EQ.0.AND.IPRNSQ.GT.1) RETURN
C
      IF (LEAVE.EQ.0) RETURN
C        IPRNSQ: 0 - PRINT AMPLITUDES AND SQUARES TO 'OUTE' : NSQ1=0,1
C                1 - PRINT ONLY AMPLITUDES TO 'OUTE'        : NSQ1=0
C                2 - NO EIGENVECTORS PRINT TO 'OUTE'
C                3-  PRINT TO 'OUTE' IN FORMAT READABLE BY TROUT
C     IF (NOUT.EQ.2.AND.NSQ1.EQ.0) RETURN
      L=MIN0(IM,MAXC,5)
      LM=0
      LM1=0
      LM2=0
      DO 128 I2=1,L
C       I=I1(I2)
C       IP1(I2)=100.D0*EVEC(I,J)
        I=NW(J,I2,M)
C       I=I+128
        IF (I.EQ.0) GO TO 129
        IP1(I2)=IEVV(J,I2,M)
C       P1=EVEC(I,J)
C       IP2(I2)=1.D2*P1+0.5D0
        IP2(I2)=DFLOAT(IP1(I2))*(IP1(I2)*1.D-6)+0.5D0
        IP1(I2)=IP1(I2)*1.D-2+0.5D0
        IF (IP1(I2).EQ.0) GO TO 129
        IF (I2.LE.3.AND.IP2(I2).GT.0) LM=I2
        IF (IP2(I2).NE.0) LM2=I2
        LM1=I2
        LC=NC(I,M)
        C1(I2)=CONFIG(LC,K3)
        Y1(I2)=YIDEN(I,NOUT,M)
        IF (KCPL.EQ.2.AND.NSQ1.EQ.1) THEN

C             UNPACK JJ-NOTATION ONLY FOR PRINTING TO 'LEVELS1..3'
C
          READ (Y1(I2),1101) Y11,J1,J2
 1101     FORMAT (A3,I2,1X,I2)
C          YJ1=DFLOAT(J1)/2
C          YJ2=DFLOAT(J2)/2
          IF (MOD(J1,2).EQ.0) THEN
            J1=J1/2
            AJ1='  '
          ELSE
            AJ1='/2'
          ENDIF
          IF (MOD(J2,2).EQ.0) THEN
            J2=J2/2
            AJ2='  '
          ELSE
            AJ2='/2'
          ENDIF
          WRITE (Y1(I2),1102) Y11,J1,AJ1,J2,AJ2
 1102     FORMAT (A3,I2,A2,1H),I2,A2)
C         WRITE (Y1(I2),1102) Y11,YJ1,YJ2
C1102     FORMAT (A3,F4.1,1H),F4.1)
        ENDIF
  128 CONTINUE
  129 TS=T(J,M)
      WS=W(J,M)
      TWS=TW(J,M)
C     IF (LEAVE.NE.0.AND.NSQ1.EQ.1) THEN
      IF (NSQ1.EQ.1) THEN
C           WRITE TO 'LEVELS1'
        IF (S.LE.1.D0) THEN
          WRITE (8,13) TS,ERAS(1),WS,TWS,GFAC(J,M),CONFIG(NJM,K3),
     1      FJ(M), (IP2(I),C1(I),Y1(I), I=1,LM)
   13     FORMAT (F13.6,A1,F10.3,F11.6,F7.3,1X,A6,A8,4X,
     1      3(I4,2H% ,A6,1X,A12))
        ELSE
          TS1=S*TS
          WS1=S*WS
          TWS1=S*TWS
          WRITE (8,14) TS1,ERAS(1),WS1,TWS1,GFAC(J,M),CONFIG(NJM,K3),
     1      FJ(M), (IP2(I),C1(I),Y1(I), I=1,LM)
   14     FORMAT (F13.3,A1,F10.0,F11.3,F7.3,1X,A6,A8,4X,
     1      3(I4,2H% ,A6,1X,A12))
        ENDIF
      ENDIF
C
      IF (NOUT.EQ.1.AND.NSQ1.EQ.0) THEN
C       IF (LEAVE.NE.0) THEN
C            WRITE TO 'RCEOUT'; 
C            WRITE TO 'RCEINP' IF INVECT.EQ.0
          READ (FJ(M),15) FJ1
   15     FORMAT (4X,F4.1)
          J1=2.D0*FJ1
          IF (MOD(J1,2).EQ.0) THEN
            J1=J1/2
            AJ1='  '
          ELSE
            AJ1='/2'
          ENDIF
          IM1=MIN0(LM1,5)
          IF (INVECT.EQ.0) WRITE (20,16)
     1      TS,ERAS(1),WS,J1,AJ1,(IP1(I),C1(I),Y1(I), I=1,IM1)
   16     FORMAT (F13.6,A1,F10.3,I4,A2,5(I4,1H ,A6,1X,A12:))
          WRITE (21,16)TS,ERAS(1),WS,J1,AJ1,(IP1(I),C1(I),Y1(I),I=1,IM1)
C       ENDIF
      ENDIF
      IF (NSQ1.EQ.0.AND.IPRNSQ.NE.2) THEN
        IF (IPRNSQ.EQ.3) THEN
          IF (J.EQ.1) THEN 
            WRITE (IW,91)
   91       FORMAT (1H )
            WRITE (IW,90) FJ(M), (YIDEN(I,NOUT,M), I=1,IM)
   90       FORMAT(1H ,A10,41X,8A10:/(50X,8A10:))
          ENDIF
          JOFS=NYOFS+(J-1)*IM
          WRITE (IW,76) TS,ERAS(1), WS, TWS, GKK(J,M),
     1      GFAC(J,M), CONFIG(NJM,K3), (BUFY(I+JOFS), I=1,IM)
   76     FORMAT (1H ,F10.3,A1,F10.3,F9.3,F7.3,F6.3,1X,A6,8F10.3:/
     1      (49X,8F10.3:))
        ELSE
          IF (J.EQ.1) WRITE (IW,20) FJ(M)
   20     FORMAT (A8)
          WRITE (IW,17) TS,ERAS(1),WS,TWS,GKK(J,M),GFAC(J,M),
     1      CONFIG(NJM,K3),(IP1(I),C1(I),Y1(I), I=1,LM1)
   17     FORMAT (F10.3,A1,F10.3,F9.3,F7.3,F6.3,1X,A6,
     1      5(I4,1H ,A6,1X,A12:),
     2      (50X,5(I4,1H ,A6,1X,A12:)))
        ENDIF
      ENDIF
C
C     IF (NSQ1.EQ.1.AND.IPRNSQ.EQ.0.AND.NOUT.EQ.1) THEN
      IF (NSQ1.EQ.1.AND.IPRNSQ.EQ.0) THEN
        IF (J.EQ.1) WRITE (IW,20) FJ(M)
        WRITE (IW,19) TS,ERAS(1),WS,TWS,GKK(J,M),GFAC(J,M),
     1    CONFIG(NJM,K3),(IP2(I),C1(I),Y1(I), I=1,LM2)
   19   FORMAT (F10.3,A1,F10.3,F9.3,F7.3,F6.3,1X,A6,
     1   5(I4,1H%,A6,1X,A12:),
     2   (/50X,5(I4,1H%,A6,1X,A12:)))
      ENDIF
C *********************************************** A.KRAMIDA ^
  900 RETURN
      END
C
      SUBROUTINE SORT8TO7
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW(IA,IX),EV(IA,3),T(IA,IX),EVAL(IA,IX),W(IA,IX),
     1  FJ(IX),KW(IA,IX),NF(IA,IX),NC(IA,IX),N3(2),S
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,
     1  YIDEN(IA,2,IX)*8
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,YIDEN
      DIMENSION F(IV),JNK(4),JXK(4)
      CHARACTER A(IV)*136,FCH(IV)*6
C
      JNK(3)=0
      JXK(3)=0
      K=1
      REWIND 8
      REWIND 7
C     WRITE (IW,9)
C   9 FORMAT (1H1)
      I=0
      JN=1
      JX=0
      IEND=0
  100 I=I+1
      READ (8,10,END=500) A(I)
   10 FORMAT (A136)
      IF (A(I)(52:53).EQ.'J=') THEN
        JX=I
        GO TO 100
      ENDIF
      IF (JX.LE.0) THEN
        JN=I+1
        GO TO 100
      ENDIF
      GO TO 120
  500 IEND=1
  120 IF (JX.NE.JN) THEN
        DO 130 J=JN,JX
  130   READ (A(J),13) F(J)
   13   FORMAT (14X,F10.3)
        N=JX-JN+1
        CALL SORT(N,F(JN),A(JN))
      ENDIF
      DO 150 J=1,JX
  150 WRITE (7,10) A(J)
      JNM1=JN-1
C     DO 160 J=1,JNM1
C 160 WRITE (IW,16) A(J)
C  16 FORMAT (1X,A136)
C     DO 170 J=JN,JX
C 170 WRITE (IW,10) A(J)
      JNK(K)=JN
      JXK(K)=JX
      IF (IEND.LE.0) THEN
        K=K+1
        A(1)=A(I)
        I=1
        JN=2
        JX=0
        GO TO 100
      ENDIF
C
C
  600 REWIND 7
      REWIND 12
      DO 650 K=1,4
        IF (JXK(K).EQ.0) GO TO 900
        K33=1+(K-1)/2
        NMAX=N3(K33)
        JN=JNK(K)
        JX=JXK(K)
        JNM1=JN-1
        DO 610 J=1,JX
  610     READ (7,10) A(J)
          DO 620 J=1,JNM1
  620     WRITE (12,10) A(J)
C         DO 625 J=JN,JX
C 625     READ (A(J),62) FCH(J)
C  62     FORMAT (43X,A6)
          DO 640 N=1,NMAX
            ERAS1=CONFIG(N,K33)
            DO 630 J=JN,JX
              IF (A(J)(44:51).EQ.ERAS1) WRITE (12,10) A(J)
  630       CONTINUE
            WRITE (12,10)
  640     CONTINUE
  650   CONTINUE
  900 STOP 'NORMAL EXIT'
      END
C
      SUBROUTINE SORT(N,F,A)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
      DIMENSION F(IV)
      CHARACTER*136 A(IV),A1
C
      J=1
   60 IF (F(J).GT.F(J+1)) THEN
        M=J
        ERAS=F(M+1)
        A1=A(M+1)
   70   F(M+1)=F(M)
        A(M+1)=A(M)
        M=M-1
        IF (M.EQ.0) GO TO 90
        IF (F(M).GT.ERAS) GO TO 70
   90   F(M+1)=ERAS
        A(M+1)=A1
      ENDIF
   80 J=J+1
      IF (J.LT.N) GO TO 60
      RETURN
      END
C
      SUBROUTINE CLOCK(T)
C
C          TIMING ROUTINE (RESULTS FOR INFORMATION ONLY)
C               FOR OTHER MACHINES, SET T=0.0 OR ADD
C               APPROPRIATE CODE.
C
      REAL*8 T
      REAL ETIME,T1,T2
      DIMENSION T2(2)
C
C------------------------------------
C          USE THIS CODE FOR A CRAY
C     T=SECOND(T)
C------------------------------------
C          USE THIS CODE FOR A VAX OR MACINTOSH
C     T1=SECNDS(0.0)
C     T=DBLE(T1)
C------------------------------------
C          USE THIS CODE FOR A SUN
C     T=ETIME(T2)
C-------------------------------------
C          USE THIS CODE FOR IBM RISC
C     IT=MCLOCK()
C     T=IT
C     T=T/100.0
C-------------------------------------
C          USE THIS CODE FOR IBM-PC LAHEY FORTRAN COMPILER
      REAL*8 TINIT,T3
      COMMON/TIMING/TINIT,ITINIT
      INTEGER HH,MM,SS,HS
      CHARACTER TIMESTR*11
      CALL TIME(TIMESTR)
      READ (TIMESTR,10) HH,MM,SS,HS
   10 FORMAT (3(I2,1X),I2)
      T3=DFLOAT(3600*HH+60*MM+SS)+HS*1.D-2
      IF (ITINIT.EQ.0) THEN
        TINIT=T3
        ITINIT=1
      ENDIF
      T=T3-TINIT
C-------------------------------------
      T=T/60.D0
      RETURN
      END
C
      SUBROUTINE CALCV
C
C          CALCULATES V MATRIX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D(JB,1),V(IV,JB),B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB)
      EQUIVALENCE (D,DELX)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW(IA,IX),EV(IA,3),T(IA,IX),EVAL(IA,IX),W(IA,IX),
     1  FJ(IX),KW(IA,IX),NF(IA,IX),NC(IA,IX),N3(2),S
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,2,IX)*8
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
      COMMON/CMLC1/C(IA,IA)
      COMMON/CMLC2/A(IA,IA)
C     COMMON/CMLC3/Y(IA,IA)
C     COMMON/CMLC4/EVEC(IA,IA)
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC(MAXB1),BUFA(MAXB2),BUFY(MAXBY),
     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2
      INTEGER*2 LF1
      DIMENSION LF1(IB)
C
C
C          CALCULATE V MATRIX
C
      IM=0
      FF = 1.D0/F
C     REWIND I216                    ! R.COWAN
C     REWIND 4                      ! R.COWAN
      KMX=0
      LMAXP1=LMAX+1
      MIK=0
      NYOFS=0
      DO M=1,MMAX
        IM=IMAX(M)
        IMSQ=IM*IM
        IMSQ1=(IMSQ+IM)/2
C       CALL RD(I216,IM,Y,IA,EVEC,IMSQ)            ! R.COWAN
C       DO 100 K=1,IM                             ! A.KRAMIDA
C       DO 100 I=1,IM
C         MIK=MIK+1
C 100   Y(I,K)=BUFY(MIK)
        JL=0
        DO L=1,LMAX
          IF (IABS(LF(L)).LT.100.OR.NOCYCE.EQ.0) THEN
            JL=JL+1
            IF (JL.GT.JB) STOP 'JB DIMENSION TOO SMALL'
C           CALL RD(4,IM,C,IA,EVEC,IMSQ1)           ! R.COWAN
            CALL RD1(M,L,IM,C,IA,BUFC,IBC,IMSQ1,    ! A.KRAMIDA
     *               ICOFFS,MAXB1,MMAX,LMAXP1)
            KM=KMX
            DO K=1,IM
              IF (NF(K,M).GE.0) THEN
                KM=KM+1
                V(KM,JL)=0.D0
                SUMIJ=0.D0
                KOFS=NYOFS+(K-1)*IM
                DO I=1,IM
                  SUMJ=0.D0
                  DO J=1,IM
                    TEMP=C(I,J)
                    IF (TEMP.NE.0.D0) THEN
                      TEMP=TEMP*BUFY(KOFS+J)       !Y(J,K)
                      IF (I.EQ.J) TEMP=TEMP*FF
                    ENDIF
                    SUMJ=SUMJ+TEMP
                  ENDDO
                  SUMIJ=SUMIJ+SUMJ*BUFY(KOFS+I)    !Y(I,K)
                ENDDO
                V(KM,JL)=SUMIJ*F
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        KMX=KM
        NYOFS=NYOFS+IMSQ
      ENDDO
C
C     IF (IPRNV.EQ.0.AND.N5.NE.0) THEN
      IF (IPRNV.EQ.0) THEN
        WRITE (IW,20)
   20   FORMAT (/34H      T        T-W         V(KM,L))
        KM=0
        DO M=1,MMAX
          IM=IMAX(M)
          WRITE (IW,21) FJ(M)
   21     FORMAT (A10)
          DO K=1,IM
            IF (NF(K,M).GE.0) THEN
              KM=KM+1
              LKM=0
              JL=0
              DO L=1,LMAX
                IF (IABS(LF(L)).LT.100) THEN
                  JL=JL+1
                  LF1(L)=JL
                  IF (DABS(V(KM,JL)).GE.1.D-4) THEN
                    LKM=LKM+1
                    LDEG(LKM)=L
                  ENDIF
                ENDIF
              ENDDO
              WRITE (IW,22) T(K,M),TW(K,M), ( V( KM,LF1( LDEG(L) ) ), 
     1          NX( LDEG(L) ), XIDEN( LDEG(L) ) , L=1,LKM)
   22         FORMAT (2F10.3,3(2X,F9.4,1H ,A2,A8)/
     1          1(20X,3(2X,F9.4,1H ,A2,A8)) )
            ENDIF
C
          ENDDO
        ENDDO
      ENDIF
C
C          IF LSQM NON-ZERO, CALCULATE V*
C
  230 IF (LSQM.NE.0) THEN
C
  260   DO KM=1,KMX
          DO J=1,JLMAX
            VC(J)=0.D0
          ENDDO
          JL=0
          DO L=1,LMAX
            IF (LF(L).GT.-100.AND.LF(L).LT.0) THEN
              JL=JL+1
              J=-LF(L)
              VC(J)=VC(J)+V(KM,JL)*X(L)
            ENDIF
          ENDDO
          DO J=1,JLMAX
            V(KM,J)=VC(J)
          ENDDO
        ENDDO
        IF (IPRNV.EQ.0) THEN
          WRITE (IW,27)
   27     FORMAT (/35H      T        T-W         V*(KM,L))
          KM=0
          DO M=1,MMAX
            IM=IMAX(M)
            WRITE (IW,21) FJ(M)
            DO K=1,IM
              IF (NF(K,M).GE.0) THEN
                KM=KM+1
                LKM=0
                DO L=1,LMAX
                  J=IABS(LF(L))
                  IF (J.LT.100.AND.DABS(V(KM,J)).GE.1.D-4) THEN
                    LKM=LKM+1
                    LDEG(LKM)=L
                  ENDIF
                ENDDO
                WRITE (IW,22) T(K,M),TW(K,M), (V(KM,IABS(LF(LDEG(L)))),
     1             NX(LDEG(L)), XIDEN(LDEG(L)), L=1,LKM)
              ENDIF
C
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C
C          CALCULATE B AND D MATRICES
C
      DO L=1,JLMAX
        DL1=0.D0
        DO J=L,JLMAX
          BLJ=0.D0
          KM=0
          DO M=1,MMAX
            IM=IMAX(M)
            DO K=1,IM
              IF (NF(K,M).GE.0) THEN
                KM=KM+1
                IF (J.EQ.L) THEN
                  DL1=DL1+V(KM,L)*TW(K,M)
                ENDIF
                BLJ=BLJ+V(KM,L)*V(KM,J)
              ENDIF
            ENDDO
          ENDDO
          B(L,J)=BLJ
          B(J,L)=BLJ
        ENDDO
        D(L,1)=DL1
        IF (B(L,L).EQ.0.D0) B(L,L)=1.D0
      ENDDO
C
C     IF (IPRNV.EQ.0.AND.N5.NE.0) THEN
C       WRITE (IW,32)
C  32   FORMAT (///30H           D            B(L,J))
C       WRITE (IW,21)
C       DO 325 L=1,JLMAX
C         IF (LSQM.GE.0) THEN
C           WRITE (IW,33) D(L,1), (B(L,J), J=1,JLMAX)
C         ELSE
C           WRITE (IW,34) D(L,1), (B(L,J), J=1,JLMAX)
C  34       FORMAT (1H ,F15.4,1PE14.3,8E11.3/(E32.3,8E11.3))
C       ENDIF
C 325   CONTINUE
C  33   FORMAT (1H , F15.4, F12.3, 10F9.3/(F30.3, 10F9.3))
C     ENDIF
      RETURN
      END
C
      SUBROUTINE CALCX
C
C                 LEAST-SQUARE CALCN OF NEW PARAMETERS (MOD 19)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D(JB,1),V(IV,JB),B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB),
      DIMENSION DELX(IB)
      EQUIVALENCE (D,DELX)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW(IA,IX),EV(IA,3),T(IA,IX),EVAL(IA,IX),W(IA,IX),
     1  FJ(IX),KW(IA,IX),NF(IA,IX),NC(IA,IX),N3(2),S
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,2,IX)*8
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
      COMMON/CMLC1/C(IA,IA)
      COMMON/CMLC2/A(IA,IA)
C     COMMON/CMLC3/Y(IA,IA)
C     COMMON/CMLC4/EVEC(IA,IA)
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC(MAXB1),BUFA(MAXB2),BUFY(MAXBY),
     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2
C
C          CALCULATE DELX AND NEW X
C
  340 DO 343 L=1,JLMAX
        DO 343 J=L,JLMAX
          IF (L.EQ.J) THEN
            D(L,L+1)=1.D0
          ELSE
            D(L,J+1)=0.D0
            D(J,L+1)=0.D0
          ENDIF
  343 CONTINUE
      JLMAX1=JLMAX+1
  350 CALL LSS(JLMAX,JLMAX1,JB,B,D,ERAS,DETB)
      SUMDX0=SUMDX
      SUMDX=0.D0
      DO 355 L=1,LMAX
        IF (IABS(LF(L)).GE.100) THEN
          DELXC(L)=0.D0
          SDX(L)=0.D0
        ELSE
          IF (LSQM.LT.0) THEN
            J=-LF(L)
            IF (DELX(J).LE.-1.D0.AND.XMIN(L).GT.0.D0) DELX(J)=-0.99D0
            DELXC(L)=X(L)*DELX(J)
            DLL=D(J,J+1)
            IF (DLL.LT.0) DLL=0
            SDX(L)=X(L)*DSQRT(DLL)
          ELSE
            J=LF(L)
            DELXC(L)=DELX(J)
            DLL=D(J,J+1)
            IF (DLL.LT.0) DLL=0
            SDX(L)=DSQRT(DLL)
          ENDIF
        ENDIF
  355 SUMDX=SUMDX+DELXC(L)
      GEFF=1.D0
      IF (NOF.GT.3) THEN
        ERS1=0.5D0*GCOE*SUMDX/SUMDX0
        GEFF=1.D0+SIGN(DMIN1(0.5D0,DABS(ERS1)),ERS1)
      ENDIF
  357 DXMAX1=DXMAX
      IF (F.LT.0.9999D0) THEN
        DO 360 L=1,LMAX
  360   GEFF=DMIN1(GEFF, DXMAX1/(ABS(DELXC(L))+0.001E0))
      ELSE
        IF (F.GE.0.7D0) THEN
          DXMAX1=2.D0*DXMAX
          DO 364 L=1,LMAX
  364     GEFF=DMIN1(GEFF, DXMAX1/(ABS(DELXC(L))+0.001E0))
        ENDIF
      ENDIF
  365 DO 366 L=1,LMAX
        X0(L)=X(L)
        DELX(L)=DELXC(L)*GEFF
  366 X(L)=X0(L)+DELX(L)
  370 DO 375 L=1,LMAX
        IF (X(L).GT.XMAX(L)) THEN
          X(L)=XMAX(L)
          DELX(L)=X(L)-X0(L)
        ELSE
          IF (X(L).LE.XMIN(L)) THEN
            X(L)=XMIN(L)
            DELX(L)=X(L)-X0(L)
          ENDIF
        ENDIF
  375 CONTINUE
C     DELX15=0.D0
C     DO 380 L=1,LMAX
C       IF (LF(L).EQ.-15) DELX15=DELX(L)
C 380 CONTINUE
C     DO 390 L=1,LMAX
C       IF (LF(L).EQ.115) THEN
C         DELX(L)=DELX15
C         X(L)=X0(L)+DELX15
C       ENDIF
C 390 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE EIGEN
C
C                CALCN OF NEW EIGENVALUES AND VECTORS
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D(JB,1),V(IV,JB),B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB),A1(IA,IA)
      EQUIVALENCE (D,DELX,A1)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW(IA,IX),EV(IA,3),T(IA,IX),EVAL(IA,IX),W(IA,IX),
     1  FJ(IX),KW(IA,IX),NF(IA,IX),NC(IA,IX),N3(2),S
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,2,IX)*8
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
      COMMON/CMLC1/C(IA,IA)
      COMMON/CMLC2/A(IA,IA)
C     COMMON/CMLC3/Y(IA,IA)
C     COMMON/CMLC4/EVEC(IA,IA)
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC(MAXB1),BUFA(MAXB2),BUFY(MAXBY),
     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2
C
      R=0.D0
      TNO=0.D0
C ======================================== R.COWAN
C     REWIND I216
C     REWIND 4
C     REWIND 15
C ======================================== R.COWAN
      LMAXP1=LMAX+1
C     MIK=0
      NYOFS=0                  ! NYOFS IS OFFSET IN THE Y-BUFFER
      DO 520 M=1,MMAX
        IM=IMAX(M)
        IMSQ=IM*IM
        IMSQ1=(IMSQ+IM)/2
C
C            CALCULATE A MATRIX
C
        IF (IPRNA.EQ.0) THEN
          DO L=1,IA
            DO I=1,IM
              A1(I,L)=0.D0
            ENDDO
          ENDDO
        ENDIF
        DO I=1,IM
          DO J=1,I
            A(I,J)=0.D0
          ENDDO
        ENDDO
        DO 402 L=1,LMAX
C ======================================== R.COWAN
C         CALL RD(4,IM,C,IA,EVEC,IMSQ1)
C ======================================== R.COWAN
C *********************************************** A.KRAMIDA v
          CALL RD1(M,L,IM,C,IA,BUFC,IBC,IMSQ1,
     *           ICOFFS,MAXB1,MMAX,LMAXP1)
C *********************************************** A.KRAMIDA ^
C *********************************************** A.KRAMIDA  v
          IF (IPRNA.EQ.0) THEN
            DO I=1,IM
C                SET UP A1(I,J) TO HOLD NON-ZERO ANGULAR COEFFICIENTS
C                FOR SLATER INTEGRALS (IN THE DIAGONAL OF THE MATRIX
C                TO BE DIAGONALIZED)
              IF (DABS(C(I,I)).GT.1.D-5) THEN
                IF (A1(I,1).EQ.0.D0) A1(I,1)=DFLOAT(L)
                I1=A1(I,1)
                A1(I,2)=DFLOAT(L)
                I2=3+L-I1
                IF (I2.GT.IA) I2=IA
                A1(I,I2)=C(I,I)
              ENDIF
            ENDDO
          ENDIF
C *********************************************** A.KRAMIDA  ^
          EVKM=X(L)
          DO I=1,IM
            DO J=1,I
              TEMP=C(I,J)
              IF (TEMP.NE.0.D0) A(I,J)=A(I,J)+TEMP*EVKM
            ENDDO
          ENDDO
  402   CONTINUE
        EVKM=0.D0
        DO I=1,IM
          DO J=1,I-1
            TEMP=A(I,J)
            IF (F.NE.1.D0) THEN
              TEMP=TEMP*F
              A(I,J)=TEMP
              A(J,I)=TEMP
            ENDIF
            TEMP=TEMP*TEMP
            EVKM=EVKM+TEMP+TEMP
          ENDDO
          TEMP=A(I,I)
          EVKM=EVKM+TEMP*TEMP
        ENDDO
        ANORM(M)=DSQRT(EVKM)
C
C            IF IPRNA IS ZERO, WRITE TAPE 9  ('OUTE')
C
C       IF (IPRNA.EQ.0.AND.N5.NE.0) THEN
C         IF (M.EQ.1) THEN
C           WRITE (IW,40) NORUN,NOCYC
C  40       FORMAT (1H1,'RUN NO=',A6,'     CYC NO=',I3,'        A(I,J)')
C         ENDIF
C         WRITE (IW,41)
C  41     FORMAT (1H0)
C         DO 408 I=1,IM
C 408     WRITE (IW,42) (A(I,J), J=1,IM)
C  42     FORMAT (1H ,F10.3,11F11.3/(4X,11F11.3))
C       ENDIF
C *********************************************** A.KRAMIDA  v
        IF (IPRNA.EQ.0) THEN
          IF (M.EQ.1) THEN
            WRITE (IW,40)
   40       FORMAT (/46H      T         ANG.COEFS FOR RADIAL INTEGRALS)
          ENDIF
          WRITE (IW,41) FJ(M)
   41     FORMAT (A8)
          DO 420 I=1,IM
            IF (NF(I,M).LT.0) GO TO 420
            I1=A1(I,1)
            I2=A1(I,2)
            IF (I1.LE.0) GO TO 420
            WRITE (IW,4003)
 4003       FORMAT (1H )
            IF (I2.GT.IA) I2=IA
            WRITE (IW,4004) T(I,M),(A1(I,J),XIDEN(I1+J-3), J=3,3+I2-I1)
 4004       FORMAT (F10.3,3(2X,F9.4,1H ,A8)/
     1        1(10X,3(2X,F9.4,1H ,A8)))
  420     CONTINUE
C
        ENDIF
C *********************************************** A.KRAMIDA  ^
C
C            CALCULATE EIGENVALUES AND VECTORS
C
C 500   CALL TRED2(IA,IM,A,EVAL(1,M),EV,EVEC)
C       CALL TQL2(IA,IM,EVAL(1,M),EV,EVEC,MM)
  500   CALL TRED2(IA,IM,A,EVAL(1,M),EV,BUFY(NYOFS+1))
        CALL TQL2(IA,IM,EVAL(1,M),EV,BUFY(NYOFS+1),MM)
        SEVSQ(M)=0.D0
        DO 506 K=1,IM
          EVKM=EVAL(K,M)
  506   SEVSQ(M)=SEVSQ(M)+EVKM*EVKM
        SEVSQ(M)=DSQRT(SEVSQ(M))
        IF (SEVSQ(M).EQ.0.D0.OR.ANORM(M).NE.0.D0.AND.
     *  DABS(SEVSQ(M)/ANORM(M)-1.D0).GT.0.000002D0) THEN
          WRITE (IW,50) M
   50     FORMAT ('0SEVSQ AND ANORM DISAGREE FOR M=',I2)
        ENDIF
C
C            READ G-VALUE MATRIX INTO C ARRAY
C
C ======================================== R.COWAN
C       CALL RD(15,IM,C,IA,A,IMSQ1)
C       READ (15)
C ======================================== R.COWAN
C
C *********************************************** A.KRAMIDA
        CALL RD1(M,LMAXP1,IM,C,IA,BUFC,IBC,IMSQ1,
     *           ICOFFS,MAXB1,MMAX,LMAXP1)
C *********************************************** A.KRAMIDA
C
        N=-10
        CALL ORDER
        NYOFS=NYOFS+IMSQ         ! INCREASE THE Y-BUFFER OFFSET
C
C       CALL WR(I216,IM,Y,IA,C,IMSQ)       ! R.COWAN
C       DO 515 K=1,IM
C       DO 515 I=1,IM
C         MIK=MIK+1
C 515   BUFY(MIK)=Y(I,K)
        CALL CLOCK(TIME)
        DELT=TIME-TIME0
        IF (IW6.LT.0) WRITE (6,52) NOCYC,FJ(M),DELT
   52   FORMAT (' CYCLE NO.=',I2, '    FINISHED DIAGONALIZATION FOR',
     1    A8,'  AT TIME=',F8.4,' MIN')
  520 CONTINUE
      IPRNA=1                   ! A.KRAMIDA: A-MATRIX DOES NOT CHANGE.
C                               ! IT IS ENOUGH TO HAVE IT PRINTED ONCE.
C
      IF (TNO.NE.0.D0) THEN
        AVDEV=DSQRT(R/TNO)
      ELSE
        AVDEV=0.D0
        NOCYCE=0
      ENDIF
      IF (TNO.GT.PNO) THEN
        STDEV=DSQRT(R/(TNO-PNO))
      ELSE
        STDEV=0.D0
      ENDIF
      RETURN
C
      END
C
      SUBROUTINE ORDER
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D(JB,1),V(IV,JB),B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB)
      EQUIVALENCE (D,DELX)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW(IA,IX),EV(IA,3),T(IA,IX),EVAL(IA,IX),W(IA,IX),
     1  FJ(IX),KW(IA,IX),NF(IA,IX),NC(IA,IX),N3(2),S
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,2,IX)*8
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
      COMMON/CMLC1/C(IA,IA)
      COMMON/CMLC2/A(IA,IA)
C     COMMON/CMLC3/Y(IA,IA)
C     COMMON/CMLC4/EVEC(IA,IA)
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C ********************************************************v  A.KRAMIDA
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC(MAXB1),BUFA(MAXB2),BUFY(MAXBY),
     1   IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2
C       DETERMINE EIGENVECTORS ORDERING FROM FIRST MAXC COMPONENTS
C       WITH LARGEST AMPLITUDES
      PARAMETER (MAXC=10)
      INTEGER*2 NW,NW0
      INTEGER*2 IEVV,IEVV0,IEVV1,NW1
      INTEGER*1 II1,II2
      COMMON/VECT/ IEVV(IA,MAXC,IX),IEVV0(IA,MAXC,IX),
     1   NW(IA,MAXC,IX),NW0(IA,MAXC,IX)
C ********************************************************^  A.KRAMIDA
C
      REAL*4 GFACT
      INTEGER*1 ICOMP
      DIMENSION GFACT(IA),IDM(IA),ICOMP(IA,IA),IEVV1(IA,MAXC),
     1  NW1(IA,MAXC)
C
      DO K=1,IM
        GFK=0.D0
        KOFS=NYOFS+(K-1)*IM
        DO J=1,IM
C         EJK=EVEC(J,K)
          EJK=BUFY(J+KOFS)
          DO I=1,IM
            GRE=C(I,J)
C           IF (GRE.NE.0.D0) GFK=GFK+EVEC(I,K)*GRE*EJK
            IF (GRE.NE.0.D0) GFK=GFK+BUFY(I+KOFS)*GRE*EJK
          ENDDO
        ENDDO
        GFAC(K,M)=GFK
        GFACT(K)=GFK
      ENDDO
C
C ========================================= R.COWAN v
C     IF (N.GT.0) GO TO 670   ! OUTPUT PHASE IS MOVED TO ORDOUT PROCEDURES
C
C          SET TO ORDER EIGENVECTORS BY LARGEST COMPONENT
C
C     DO 603 I=1,IM
C 603 KW(I,M)=0
C 604 DO 614 J=1,IM
C       GRE=0.D0
C       DO 606 I=1,IM
C         EV(I,1)=EVEC(I,J)**2
C         IF (EV(I,1).GT.GRE) THEN
C           GRE=EV(I,1)
C           K=I
C         ENDIF
C 606   CONTINUE
C       DO 608 I=1,IM
C         IF (I.NE.K.AND.EV(I,1)/GRE.GE.CRIT) GO TO 614
C 608   CONTINUE
C       IF (GRE.GT.0.5D0) THEN
C         KW(K,M)=J
C       ELSE
C         KWKM=KW(K,M)
C         IF (KWKM.LE.0) THEN
C           IF (KWKM.LT.0) THEN
C             KWKM=-100
C           ELSE
C             KWKM=-J
C           ENDIF
C           KW(K,M)=KWKM
C         ENDIF
C       ENDIF
C 614 CONTINUE
C     DO 622 K=1,IM
C       IF (KW(K,M).LE.-100) THEN
C         KW(K,M)=0
C       ELSE
C         KW(K,M)=IABS(KW(K,M))
C       ENDIF
C 622 CONTINUE
C
C          FOR EIGENVECTORS WITH NO DOMINANT COMPONENT,
C              SET TO ORDER BY MAGNITUDE OF EVAL (N NEG) OR GFAC (N POS)
C
C     N1I=1
C 631 DO 633 N1=N1I,IM
C       IF (NC(N1,M).EQ.NC(N1I,M)) N1F=IM
C 633 CONTINUE
C 631 N1F=IM
C 634 ICNTR=0
C     GRE=-1.0D30
C     DO 641 N1=N1I,N1F
C       IF (KW(N1,M).EQ.0) THEN
C         ICNTR=ICNTR+1
C         IF (N.GT.0) THEN
C           IF (GKK(N1,M).GE.GRE) THEN
C             GRE=GKK(N1,M)
C             K=N1
C           ENDIF
C         ELSE
C           IF (T(N1,M).GE.GRE) THEN
C             GRE=T(N1,M)
C             K=N1
C           ENDIF
C         ENDIF
C       ENDIF
C 641 CONTINUE
C     IF (ICNTR.NE.0) THEN
C       GRE=-1.0D30
C       DO 649 N1=N1I,N1F
C         DO 643 I=N1I,N1F
C           IF (IABS(KW(I,M)).EQ.N1) GO TO 649
C 643     CONTINUE
C         IF (N.LE.0) THEN
C           IF (EVAL(N1,M).GT.GRE) THEN
C             GRE=EVAL(N1,M)
C             J=N1
C           ENDIF
C         ELSE
C           IF (GFAC(N1,M).GT.GRE) THEN
C             GRE=GFAC(N1,M)
C             J=N1
C           ENDIF
C         ENDIF
C 649   CONTINUE
C       KW(K,M)=-J
C       GO TO 634
C     ENDIF
C 650 N1I=N1F+1
C     IF (N1I.LE.IM) GO TO 631
C ========================================= R.COWAN  ^
C
C ***************************************** A.KRAMIDA v
C       MAP EIGENVECTORS INTO IEVV(I,J,M) AND NW(I,J,M), J=1..MAXC.
C       K=NW(I,J,M) KEEPS THE REFERENCE TO EVEC(I,K) FOR CURRENT M:
C          IEVV(I,J,M)=EVEC(I,K)
C       IEVV(I,J,M) IS SORTED BY AMPLITUDE: IEVV(I,J,M)>=IEVV(I,J+1,M)
C
      N1F=MIN0(IM,MAXC)
      DO N2=1,IM
        DO K=1,IM
          KW(K,M)=-K
        ENDDO
        N2OFS=NYOFS+(N2-1)*IM
        DO K=1,N1F
          GRE=0
          K1=0
          DO N1=1,IM
            IF (KW(N1,M).LE.0) THEN
C             EV1=DABS(EVEC(N1,N2))
              EV1=DABS(BUFY(N1+N2OFS))
              IF (EV1.GT.GRE) THEN
                K1=N1
                GRE=EV1
              ENDIF
            ENDIF
          ENDDO
          NW(N2,K,M)=K1 ! -128 FOR 1-BYTE INTEGER STORAGE
          IF (K1.GT.0) THEN
            KW(K1,M)=K1
C           IEVV(N2,K,M)=100.D0*EVEC(K1,N2)+0.5
            IEVV(N2,K,M)=1.D4*BUFY(K1+N2OFS)+0.5D0
          ELSE
            IEVV(N2,K,M)=0
          ENDIF
        ENDDO
        DO I=N1F+1,MAXC
          IEVV(N2,I,M)=0
          NW(N2,I,M)=0 ! -128 FOR 1-BYTE INTEGER STORAGE
        ENDDO
      ENDDO
      N1F=MIN0(IM,MAXC)
      IF ((NOCYC.EQ.0.AND.INVECT.EQ.0).OR.NOCYC.LT.IORD) THEN
C           NO REORDERING!
        DO K=1,IM
          KW(K,M)=0
        ENDDO
  632   ICNTR=0
        GRE=-1.0D30
        DO N1=1,IM
          IF (KW(N1,M).EQ.0) THEN
            ICNTR=1
            IF (T(N1,M).GE.GRE) THEN
              GRE=T(N1,M)
              K=N1
            ENDIF
          ENDIF
        ENDDO
        IF (ICNTR.EQ.0) GO TO 636
        GRE=-1.0D30
        DO 635 N1=1,IM
          DO I=1,IM
            IF (IABS(KW(I,M)).EQ.N1) GO TO 635
          ENDDO
          IF (EVAL(N1,M).GT.GRE) THEN
            GRE=EVAL(N1,M)
            J=N1
          ENDIF
  635   CONTINUE
        KW(K,M)=J
        GO TO 632
  636   CONTINUE
      ELSE
C       IDENTIFY THE NEW EIGENVECTORS WITH THE OLD ONES
        DO J=1,IM
          IDM(J)=0
          DO I=1,IM
            ICOMP(I,J)=1
          ENDDO
          KW(J,M)=0
        ENDDO
        I=0
C       DO 660 I=1,IM
  641   IF (I.EQ.IM) GO TO 660
        I=I+1
  642   IDM(I)=1000000000
        JMIN=0
C
C         FIND IEVV0(J,...,M) CLOSEST TO IEVV(I,...,M)
        IFAIL=0
        DO 652 JS=1,2
        DO 652 J=1,IM
          IF (ICOMP(I,J).EQ.0) GO TO 652  ! SKIP COMPARISOM
C           CALC DISTANCE OF VECTOR IEVV0(I,...,M) FROM IEVV(J,...,M)
C           IN IM-DIMENSIONAL SPACE
C
          IDD=0
          NFOUND=0
          DO K=1,N1F
            N1=NW(I,K,M)
            IA0=IEVV(I,K,M)
            N2=0
            DO K1=1,N1F
              NK1=NW0(J,K1,M)
              IF (NK1.EQ.0) GOTO 646   ! -128 FOR INTEGER*1 NW,NW0
              IF (NK1.EQ.N1) THEN
                N2=K1
                NFOUND=NFOUND+1
                GOTO 646
              ENDIF
            ENDDO
  646       IA1=0
            IF (N2.NE.0) IA1=IEVV0(J,N2,M)
            IF (JS.EQ.2) IA1=-IA1     ! CHECK FOR INVERTED SIGN
            IA10=IA1-IA0
            IDD=IDD+IA10*IA10
          ENDDO
          IF (NFOUND.LT.N1F) THEN
C                ADD COMPONENTS OF NW0(J) MISSING IN NW(I)
            DO 650 K=1,N1F
              N1=NW0(J,K,M)
              IA0=IEVV0(J,K,M)
              N2=0
              DO K1=1,N1F
                NK1=NW(I,K1,M)
                IF (NK1.EQ.0) GOTO 649   ! -128 FOR INTEGER*1 NW,NW0
                IF (NK1.EQ.N1) GOTO 650  ! SKIP COMPONENT PROCESSED EARLIER
              ENDDO
  649         IDD=IDD+IA0*IA0
  650       CONTINUE
          ENDIF
          IF (IDD.LT.IDM(I)) THEN
            JMIN=J
            IDM(I)=IDD
C        DISTANCE OF 0.10 IS WELL ENOUGH TO IDENTIFY THE EIGENVECTORS.
C        IT CORRESPONDS TO VALUE OF PROJECTION (SCALAR MULT.) OF 0.95
            IF (IDM(I).LT.1000000) GO TO 653
          ENDIF
  652   CONTINUE
        IF (JMIN.EQ.0) IFAIL=1    ! COMPARISON LIST IS EMPTY FOR I
  653   IF (IFAIL.EQ.0) THEN
          JK=IABS(KW(JMIN,M))
          IF (JK.NE.0) THEN
C           IF (ISORT.EQ.0) THEN
C            COMPARE DISTANCES IDM(I) AND IDM(JK)
              ID1=IDM(JK)
              IF (IDM(I).GT.ID1) THEN
C                 EXCLUDE JMIN FROM COMPARE LIST FOR CURRENT I
C                 AND REDO THE CYCLE FOR CURRENT I
                ICOMP(I,JMIN)=0
              ELSE
C                 EXCLUDE JMIN FROM COMPARE LIST FOR I=JK
C                 AND RETURN BACK TO I=JK.
C                 NOTE THAT JK CAN NOT BE GREATER THAN I
                ICOMP(JK,JMIN)=0
                DO J=1,IM
                  KWJM=KW(J,M)
                  IF (KWJM.NE.0.AND.KWJM.GE.JK) KW(J,M)=0
                ENDDO
                I=JK
              ENDIF
              GO TO 642
C           ELSE
C             KW(JMIN,M)=-JK
C           ENDIF
          ENDIF
        ENDIF
        IF (IFAIL.NE.0) THEN
C
C             IDENTIFICATION FAILED
C
          WRITE (*,*) 'Failed to identify eigenvectors in RCEINP.'
          STOP '(Fix some of parameters or change MAXC and recompile.)'
        ENDIF
        IF (KW(JMIN,M).GE.0) KW(JMIN,M)=I
        GO TO 641
  660   CONTINUE
C       IF (ISORT.NE.0) THEN
C            FOR EIGENVECTORS NOT MAPPED AT THIS POINT,
C            SET TO ORDER BY MAGNITUDE OF EVAL (N NEG) OR GFAC (N POS)
C  
C         DO I=1,IM
C           IF (KW(I,M).LT.0) KW(I,M)=0
C         ENDDO
C 662     ICNTR=0
C         GRE=-1.0D30
C         DO N1=1,IM
C             FIND THE GREATEST EIGENVALUE (NUMBER K) NOT MAPPED YET
C           IF (KW(N1,M).EQ.0) THEN
C             ICNTR=ICNTR+1
C             IF (N.GT.0) THEN
C               IF (GKK(N1,M).GE.GRE) THEN
C                 GRE=GKK(N1,M)
C                 K=N1
C               ENDIF
C             ELSE 
C               IF (T(N1,M).GE.GRE) THEN
C                 GRE=T(N1,M)
C                 K=N1
C               ENDIF
C             ENDIF
C           ENDIF
C         ENDDO
C         IF (ICNTR.NE.0) THEN
C           GRE=-1.0D30
C           DO 665 N1=1,IM
C             DO I=1,IM
C               IF (KW(I,M).EQ.N1) GO TO 665
C             ENDDO
C             IF (N.LE.0) THEN
C               IF (EVAL(N1,M).GT.GRE) THEN
C                 GRE=EVAL(N1,M)
C                 J=N1
C               ENDIF
C             ELSE
C               IF (GFAC(N1,M).GT.GRE) THEN
C                 GRE=GFAC(N1,M)
C                 J=N1
C               ENDIF
C             ENDIF
C 665       CONTINUE
C           KW(K,M)=J
C           GO TO 662
C         ENDIF
C       ENDIF
      ENDIF
C ***************************************** A.KRAMIDA ^
C
C          REARRANGE EIGENVALUES, VECTORS, AND GFACTORS
C
      N1F=MIN0(IM,MAXC)               ! A.KRAMIDA
      DO K=1,IM
        J=IABS(KW(K,M))
        W(K,M)=EVAL(J,M)
        GFAC(K,M)=GFACT(J)
C       KOFS=NYOFS+(K-1)*IM
C       JOFS=NYOFS+(J-1)*IM
C       DO I=1,IM
C         Y(I,K)=EVEC(I,J)           ! Y(I,K) IS CHANGED
C         BUFY(I+KOFS)=EVEC(I,J)
C         BUFY(I+KOFS)=EV(I,1)
C       ENDDO
C ***************************************** A.KRAMIDA v
        DO N1=1,N1F
          IEVV1(K,N1)=IEVV(J,N1,M)
          NW1(K,N1)=NW(J,N1,M)
        ENDDO
C ***************************************** A.KRAMIDA ^
        IF (NF(K,M).EQ.0) THEN
          TW(K,M)=T(K,M)-W(K,M)
          GRE=TW(K,M)
          R=R+GRE*GRE
          TNO=TNO+1.D0
        ENDIF
      ENDDO
C           REORDER EIGENVECTORS DIRECTLY IN THE Y-BUFFER
      DO K=1,IM
        J=IABS(KW(K,M))
        IF (J.NE.K) THEN
          KOFS=NYOFS+(K-1)*IM
          JOFS=NYOFS+(J-1)*IM
          DO I=1,IM
            GRE=BUFY(I+KOFS)
            BUFY(I+KOFS)=BUFY(I+JOFS)
            BUFY(I+JOFS)=GRE
          ENDDO
          ICTRL=0
          DO I=K,IM
            IF (KW(I,M).EQ.K) THEN
              ICTRL=1
              KW(I,M)=J
              GO TO 655
            ENDIF
          ENDDO
655       CONTINUE
          IF (ICTRL.EQ.0) THEN
            STOP '**** EGENVECTORS REORDERING FAILED. ****'
          ENDIF
        ENDIF
      ENDDO
      DO N1=1,N1F
        DO K=1,IM
          IEVV(K,N1,M)=IEVV1(K,N1)
          NW(K,N1,M)=NW1(K,N1)
        ENDDO
      ENDDO
      IF (NOCYC.EQ.IORD.OR.(NOCYC.GT.IORD.AND.ISAMPL.EQ.1)) THEN
        DO N1=1,N1F
          DO K=1,IM
            IEVV0(K,N1,M)=IEVV(K,N1,M)
            NW0(K,N1,M)=NW(K,N1,M)
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
C
      SUBROUTINE ORDOUT1
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D(JB,1),V(IV,JB),B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB)
      EQUIVALENCE (D,DELX)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW(IA,IX),EV(IA,3),T(IA,IX),EVAL(IA,IX),W(IA,IX),
     1  FJ(IX),KW(IA,IX),NF(IA,IX),NC(IA,IX),N3(2),S
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,2,IX)*8
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
      COMMON/CMLC1/C(IA,IA)
      COMMON/CMLC2/A(IA,IA)
C     COMMON/CMLC3/Y(IA,IA)
C     COMMON/CMLC4/EVEC(IA,IA)
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C ********************************************************v  A.KRAMIDA
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC(MAXB1),BUFA(MAXB2),BUFY(MAXBY),
     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2
C       DETERMINE EIGENVECTORS ORDERING FROM FIRST MAXC COMPONENTS
C       WITH LARGEST AMPLITUDES
      PARAMETER (MAXC=10)
      INTEGER*2 NW,NW0,KW1
      INTEGER*2 IEVV,IEVV0
      INTEGER*1 II1,II2,IEX
      COMMON/VECT/ IEVV(IA,MAXC,IX),IEVV0(IA,MAXC,IX),
     1  NW(IA,MAXC,IX),NW0(IA,MAXC,IX)
C ********************************************************^  A.KRAMIDA
      INTEGER*1 ICOMP
      DIMENSION IDM(IA),ICOMP(IA,IA),KW1(IA)
C
      N1F=MIN0(IM,MAXC)
      DO J=1,IM
        NW0(J,1,M)=J
        IEVV0(J,1,M)=10000
        IDM(J)=0
        KW1(J)=0
        DO NM1=2,N1F
          NW0(J,NM1,M)=0
          IEVV0(J,NM1,M)=0
        ENDDO
        DO I=1,IM
          ICOMP(I,J)=0
        ENDDO
      ENDDO
C       IDENTIFY THE EIGENVECTORS WITH THE SAMPLE ONES
      I=0
  641 IF (I.EQ.IM) GO TO 660
      I=I+1
  642 IDM(I)=1000000000
      JMIN=0
C
C       FIND IEVV0(J,...,M) CLOSEST TO IEVV(I,...,M)
      IFAIL=0
      DO 652 JS=1,2
      DO 652 J=1,IM
        IF (ICOMP(I,J).GT.0) GO TO 652  ! SKIP COMPARISOM
C         CALC DISTANCE OF VECTOR IEVV0(I,...,M) FROM IEVV(J,...,M)
C         IN IM-DIMENSIONAL SPACE
C
        IDD=0
        NFOUND=0
C       NK1=NW0(J,1,M)
        IA11=IEVV0(J,1,M)
        IAM=IEVV(I,1,M)
        DO K=1,N1F
          N1=NW(I,K,M)
          IA0=IEVV(I,K,M)
          IA1=0
          IF (N1.EQ.J) THEN
            NFOUND=1
            IF (DFLOAT(IA0*IA0)/DFLOAT(IAM*IAM).LT.CRIT1) GO TO 652
            IF (JS.EQ.2) THEN 
              IA1=-IA11     ! CHECK FOR INVERTED SIGN
            ELSE
              IA1=IA11
            ENDIF
          ENDIF
          IA10=IA1-IA0
          IDD=IDD+IA10*IA10
        ENDDO
        IF (NFOUND.EQ.0) GO TO 652   !  IDD=IDD+IA11*IA11
        IF (IDD.LT.IDM(I)) THEN
          JMIN=J
          IDM(I)=IDD
C
C           IDD is a square of (distance*10000).
C           Distance of 0.76 (IDD=58000000) is well enough to identify 
C           the eigenvectors. It corresponds to value of projection 
C           (scalar mult.) of 0.71, and to 50% composition of the component.
C           Projection = 1.0 - (distance**2)/2
C
          IF (IDM(I).LT.58000000) GO TO 653
        ENDIF
  652 CONTINUE
      IF (JMIN.EQ.0) IFAIL=1    ! COMPARISON LIST IS EMPTY FOR I
  653 IF (IFAIL.EQ.0) THEN
        JK=IABS(KW1(JMIN))
        IF (JK.NE.0) THEN
C           COMPARE DISTANCES IDM(I) AND IDM(JK)
          ID1=IDM(JK)
          IF (IDM(I).GT.ID1) THEN
C               EXCLUDE JMIN FROM COMPARE LIST FOR CURRENT I
C               AND REDO THE CYCLE FOR CURRENT I
            ICOMP(I,JMIN)=1
          ELSE
C               EXCLUDE JMIN FROM COMPARE LIST FOR I=JK
C               AND RETURN BACK TO I=JK.
C               NOTE THAT JK CAN NOT BE GREATER THAN I
            DO J=1,IM
              IF (KW1(J).GE.JK) KW1(J)=0
            ENDDO
            I=JK
            ICOMP(JK,JMIN)=1
          ENDIF
          GO TO 642
        ENDIF
      ENDIF
 654  IF (IFAIL.NE.0) THEN
C
C             IDENTIFICATION FAILED
C
          WRITE (IW,670) FJ(M)
 670      FORMAT 
     1      (34H Failed to label eigenvectors for ,A8)
          RETURN
      ENDIF
C       IF WE GET HERE, IDENTIFICATION I -> JMIN IS FOUND
      IF (KW1(JMIN).GE.0) KW1(JMIN)=I
      GO TO 641
 660  CONTINUE
      DO J=1,IM            ! LOOP THRU TERMS
        K=KW1(J)           ! K IS INDEX OF EIGENVECTOR LABELED BY TERM J
        DO N1=1,N1F    ! FIND THE JTH TERM AMONG ALL COMPONENTS
          NN1=NW(K,N1,M)
          IF (NN1.EQ.0) GO TO 680
          IF (NN1.EQ.J.AND.N1.NE.1) THEN
            IA1=IEVV(K,N1,M)
            NK1=NW(K,N1,M)
            DO N2=N1-1,1,-1
              IEVV(K,N2+1,M)=IEVV(K,N2,M)
              NW(K,N2+1,M)=NW(K,N2,M)
            ENDDO
            IEVV(K,1,M)=IA1
            NW(K,1,M)=NK1
            GO TO 680
          ENDIF
        ENDDO
 680    CONTINUE
      ENDDO
      RETURN
      END

      SUBROUTINE ORDOUT2
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY 
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCEPAR.FOR'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D(JB,1),V(IV,JB),B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB)
      EQUIVALENCE (D,DELX)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW(IA,IX),EV(IA,3),T(IA,IX),EVAL(IA,IX),W(IA,IX),
     1  FJ(IX),KW(IA,IX),NF(IA,IX),NC(IA,IX),N3(2),S
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,2,IX)*8
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
      COMMON/CMLC1/C(IA,IA)
      COMMON/CMLC2/A(IA,IA)
C     COMMON/CMLC3/Y(IA,IA)
C     COMMON/CMLC4/EVEC(IA,IA)
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C ********************************************************v  A.KRAMIDA
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC(MAXB1),BUFA(MAXB2),BUFY(MAXBY),
     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2
C       DETERMINE EIGENVECTORS ORDERING FROM FIRST MAXC COMPONENTS
C       WITH LARGEST AMPLITUDES
      PARAMETER (MAXC=10)
      INTEGER*2 NW,NW0
      INTEGER*2 IEVV,IEVV0
      INTEGER*1 II1,II2
      COMMON/VECT/ IEVV(IA,MAXC,IX),IEVV0(IA,MAXC,IX),
     1  NW(IA,MAXC,IX),NW0(IA,MAXC,IX)
C ********************************************************^  A.KRAMIDA
C
      REAL*4 GFACT
      INTEGER*2 KE
      DIMENSION GFACT(IA),KE(IA)
C
C       REARRANGE EVEC(I,J),IEVV(I,K,J),NW(I,K,J) FOR THE SECOND
C       COUPLING SCHEME, USING TRANSFORMATION MATRIX A
C
      DO K=1,IM
        GKKM=0.D0
        KOFS=NYOFS+(K-1)*IM
        DO J=1,IM
          EJK=0.D0
          AJK=A(J,K)
          DO I=1,IM
            GKKM=GKKM+A(I,K)*C(I,J)*AJK
C           EJK=EJK+A(I,J)*Y(I,K)
          ENDDO
          NF1=MIN0(MAXC,IM)
          DO I=1,NF1
            I1=NW(K,I,M)
            IF (I1.GT.0) EJK=EJK+A(I1,J)*DFLOAT(IEVV(K,I,M))*1.D-4
          ENDDO
C         EVEC(J,K)=EJK
          BUFY(J+KOFS)=EJK
        ENDDO
        GKK(K,M)=GKKM
      ENDDO
C ***************************************** A.KRAMIDA v
C       MAP EIGENVECTORS INTO IEVV(I,J,M) AND NW(I,J,M), J=1..MAXC.
C       K=NW(I,J,M) KEEPS THE REFERENCE TO EVEC(I,K) FOR CURRENT M:
C          IEVV(I,J,M)=EVEC(I,K);
C       IEVV(I,J,M) IS SORTED BY AMPLITUDE: IEVV(I,J,M)>=IEVV(I,J+1,M);
C       EVEC(I,K) is located inside linear array BUFY using offset pointer.
C
      N1F=MIN0(IM,MAXC)
      DO N2=1,IM
        DO K=1,IM
          KE(K)=-K
        ENDDO
        N2OFS=NYOFS+(N2-1)*IM
        DO K=1,N1F
          GRE=0
          K1=0
          DO N1=1,IM
            IF (KE(N1).LE.0) THEN
C             EV1=DABS(EVEC(N1,N2))
              EV1=DABS(BUFY(N1+N2OFS))
              IF (EV1.GT.GRE) THEN
                K1=N1
                GRE=EV1
              ENDIF
            ENDIF
          ENDDO
          NW(N2,K,M)=K1       ! -128 for 1-byte type
          IF (K1.GT.0) THEN
            KE(K1)=K1
C           IEVV(N2,K,M)=100.D0*EVEC(K1,N2)+0.5D0
            IEVV(N2,K,M)=1.D4*BUFY(K1+N2OFS)+0.5D0
          ELSE
            IEVV(N2,K,M)=0
          ENDIF
        ENDDO
        DO I=N1F+1,MAXC
          IEVV(N2,I,M)=0
          NW(N2,I,M)=0 ! -128 FOR INTEGER*1 NW
        ENDDO
      ENDDO
  680 RETURN
      END
C
C
C     SUBROUTINE WR(IFS,IM,A,IA,B,IMSQ)
C
C     IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSION A(IA,IA),B(IMSQ)
C
C     IM2=IM*IM
C     K=0
C     DO 100 J=1,IM
C       I1=J
C       IF (IMSQ.GE.IM2) I1=1
C       DO 100 I=I1,IM
C         K=K+1
C 100 B(K)=A(I,J)
C     WRITE (IFS) (B(I),I=1,IMSQ)
C     RETURN
C     END
C
C     SUBROUTINE RD(IFS,IM,A,IA,B,IMSQ)
C
C     IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSION A(IA,IA),B(IMSQ)
C
C     READ (IFS) (B(I),I=1,IMSQ)
C     IM2=IM*IM
C     K=0
C     DO 100 J=1,IM
C       I1=J
C       IF (IMSQ.GE.IM2) I1=1
C       DO 100 I=I1,IM
C         K=K+1
C 100 A(I,J)=B(K)
C     IF (IMSQ.LT.IM2) THEN
C       DO 200 J=1,IM
C       DO 200 I=J,IM
C 200     A(J,I)=A(I,J)
C     ENDIF
C 300 RETURN
C     END
C
      SUBROUTINE WR1(IND1,IND2,IM,A,IA,B,IBIND,IMSQ,
     *               IOFFS,MAXLEN,MAXI1,MAXI2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*12 B,BB
      INTEGER*2 BB1,I,J,I1
      REAL*8 R8
      DIMENSION A(IA,IA),B(MAXLEN),BB1(6),IBIND(MAXI1,MAXI2)
      EQUIVALENCE (R8,BB1(1)),(BB,BB1(1)),(I,BB1(5)),(J,BB1(6))
C
      IM2=IM*IM
      IBIND(IND1,IND2)=IOFFS+1
      IFILL=0
      DO 100 J=1,IM
        I1=J
        IF (IMSQ.GE.IM2) I1=1
        DO 100 I=I1,IM
          R8=A(I,J)
          IF (DABS(R8).GT.1.D-30) THEN
            IF (IOFFS.GE.MAXLEN) THEN
              WRITE (*,10) MAXLEN
   10         FORMAT(11H More than ,I9,28H non-zero elements in matrix)
              STOP '(Dimension error)'
            ENDIF
            IOFFS=IOFFS+1
            B(IOFFS)=BB
            IFILL=1
          ENDIF
  100 CONTINUE
      IF (IFILL.EQ.0) THEN
        IOFFS=IOFFS+1
        I=0
        J=0
        R8=0.D0
        B(IOFFS)=BB
      ENDIF
      RETURN
      END
C
      SUBROUTINE RD1(IND1,IND2,IM,A,IA,B,IBIND,IMSQ,
     *               MOFFS,MAXLEN,MAXI1,MAXI2)
C
C     PARAMETER (IB=200)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*12 B,BB
      INTEGER*2 BB1,I,J
      REAL*8 R8
      DIMENSION A(IA,IA),B(MAXLEN),BB1(6),IBIND(MAXI1,MAXI2)
      EQUIVALENCE (R8,BB1(1)),(BB,BB1(1)),(I,BB1(5)),(J,BB1(6))
C
      IOFFS1=IBIND(IND1,IND2)
      IF (IND2.LT.MAXI2) THEN
        IOFFS2=IBIND(IND1,IND2+1)
      ELSE
        IF (IND1.LT.MAXI1) THEN
          IOFFS2=IBIND(IND1+1,1)
        ELSE
          IOFFS2=MOFFS+1
        ENDIF
      ENDIF
      DO 10 J=1,IM
        DO 10 I=1,IM
   10 A(I,J)=0.D0
C
      IM2=IM*IM
      IF (IMSQ.LT.IM2) THEN
        DO K=IOFFS1,IOFFS2-1
          BB=B(K)
          IF (I.EQ.0) RETURN
          A(I,J)=R8
          A(J,I)=R8
        ENDDO
      ELSE
        DO K=IOFFS1,IOFFS2-1
          BB=B(K)
          IF (I.EQ.0) RETURN
          A(I,J)=R8
        ENDDO
      ENDIF
  300 RETURN
      END
C
C
      SUBROUTINE LSS(N,M,I,A,B,ERAS,DET)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 ERAS
      DIMENSION A(I,N), B(I,M)    !,D(N)
C
      NN=N
      IF (NN.LE.0.OR.NN.GT.I) THEN
        WRITE (IW,99) N
   99   FORMAT ('0LSS CALLED WITH N.LE.0 OR N.GT.I.  N= ',I4)
        DET=0.0D0
      ELSE
        MM=M
        SN=1.0D0
        DO 9 J=1,NN
          L=J-1
          IF (J.NE.NN) THEN
            T=DABS(A(J,J))
            M1=J
            M2=J+1
            DO 1 K=M2,NN
              X=DABS(A(K,J))
              IF (X.GT.T) THEN
                T=X
                M1=K
              ENDIF
1           CONTINUE
            IF (M1.NE.J) THEN
              DO 2 K=1,NN
                T=A(J,K)
                A(J,K)=A(M1,K)
2             A(M1,K)=T
              SN=-SN
              IF (MM.GT.0) THEN
                DO 3 K=1,MM
                  T=B(J,K)
                  B(J,K)=B(M1,K)
3               B(M1,K)=T
              ENDIF
            ENDIF
4           IF (A(J,J).EQ.0.D0) GO TO 13
            AJJ1=1.D0/A(J,J)
            DO 6 K=M2,NN
              S1=0.0D0
              S2=0.0D0
              IF (L.NE.0) S1=DOTPRO(L,A(J,1),I,A(1,K),1)
              A(J,K)=(A(J,K)-S1)*AJJ1
              S2=DOTPRO(J,A(K,1),I,A(1,M2),1)
6           A(K,M2)=A(K,M2)-S2
          ENDIF
          IF (MM.GT.0) THEN
            IF (A(J,J).EQ.0.D0) GO TO 13
            AJJ1=1.D0/A(J,J)
            DO 8 K=1,MM
              S1=0.0D0
              IF (L.NE.0) S1=DOTPRO(L,A(J,1),I,B(1,K),1)
8           B(J,K)=(B(J,K)-S1)*AJJ1
          ENDIF
9       CONTINUE
        DET=A(1,1)*SN
        IF (DET.EQ.0.D0) GO TO 13
        IF (N.EQ.1) GO TO 15
        DO 10 J=2,NN
10      DET=DET*A(J,J)
        IF (DET.EQ.0.D0) GO TO 13
        IF (MM.EQ.0) GO TO 15
        M3=NN-1
        DO 12 J=1,MM
          DO 11 L=1,M3
            M1=NN-L
            S1=0.0D0
            M2=M1+1
            K=NN-M2+1
            S1=DOTPRO(K,A(M1,M2),I,B(M2,J),1)
11        B(M1,J)=B(M1,J)-S1
12      CONTINUE
        GO TO 15
   13   WRITE (IW,98) N
   98   FORMAT ('0LSS SINGULAR SYSTEM. INPUT DESTROYED.  N= ',I4)
        DET=0.0D0
      ENDIF
15    RETURN
      END
C
      FUNCTION DOTPRO(N,X,IX,Y,IY)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(1),Y(1)
C
      DOTPRO=0.0D0
      IF (N.LE.0) RETURN
      SUM=0.0D0
      J=1-IX
      K=1-IY
      DO 100 I=1,N
        J=J+IX
        K=K+IY
  100 SUM=SUM+X(J)*Y(K)
      DOTPRO=SUM
      RETURN
      END
C
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                    TRED2  3
C     DIMENSION A(NM,N),D(N),E(N),Z(NM,N)                               TRED2 56
C***BEGIN PROLOGUE  TRED2                                               TRED2  4
C***DATE WRITTEN   760101   (YYMMDD)                                    TRED2  5
C***REVISION DATE  830518   (YYMMDD)                                    TRED2  6
C***CATEGORY NO.  D4C1B1                                                TRED2  7
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK                          TRED2  8
C***AUTHOR  SMITH, B. T., ET AL.                                        TRED2  9
C***PURPOSE  REDUCE REAL SYMMETRIC MATRIX TO SYMMETRIC TRIDIAGONAL      TRED2 10
C            MATRIX USING AND ACCUMULATING ORTHOGONAL TRANSFORMATION    TRED2 11
C***DESCRIPTION                                                         TRED2 12
C                                                                       TRED2 13
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,    TRED2 14
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   TRED2 15
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   TRED2 16
C                                                                       TRED2 17
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A              TRED2 18
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING               TRED2 19
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            TRED2 20
C                                                                       TRED2 21
C     ON INPUT                                                          TRED2 22
C                                                                       TRED2 23
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TRED2 24
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TRED2 25
C          DIMENSION STATEMENT.                                         TRED2 26
C                                                                       TRED2 27
C        N IS THE ORDER OF THE MATRIX.                                  TRED2 28
C                                                                       TRED2 29
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          TRED2 30
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               TRED2 31
C                                                                       TRED2 32
C     ON OUTPUT                                                         TRED2 33
C                                                                       TRED2 34
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    TRED2 35
C                                                                       TRED2 36
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         TRED2 37
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      TRED2 38
C                                                                       TRED2 39
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                TRED2 40
C          PRODUCED IN THE REDUCTION.                                   TRED2 41
C                                                                       TRED2 42
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.            TRED2 43
C                                                                       TRED2 44
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        TRED2 45
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         TRED2 46
C     ------------------------------------------------------------------TRED2 47
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW, TRED2 48
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-    TRED2 49
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,    TRED2 50
C                 1976.                                                 TRED2 51
C***ROUTINES CALLED  (NONE)                                             TRED2 52
C***END PROLOGUE  TRED2                                                 TRED2 53
C                                                                       TRED2 54
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,K,L,N,II,NM,JP1                                       TRED2 55
      DIMENSION A(NM,N),D(N),E(N),Z(N,N)                                TRED2 56
C     REAL A(NM,N),D(N),E(N),Z(NM,N)                                    TRED2 56
C     REAL F,G,H,HH,SCALE                                               TRED2 57
C                                                                       TRED2 58
C***FIRST EXECUTABLE STATEMENT  TRED2                                   TRED2 59
      DO 100 I = 1, N                                                   TRED2 60
C                                                                       TRED2 61
         DO 100 J = 1, I                                                TRED2 62
            Z(I,J) = A(I,J)                                             TRED2 63
  100 CONTINUE                                                          TRED2 64
C                                                                       TRED2 65
      IF (N .EQ. 1) GO TO 320                                           TRED2 66
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........               TRED2 67
      DO 300 II = 2, N                                                  TRED2 68
         I = N + 2 - II                                                 TRED2 69
         L = I - 1                                                      TRED2 70
         H = 0.D0                                                       TRED2 71
         SCALE = 0.D0                                                   TRED2 72
         IF (L .LT. 2) GO TO 130                                        TRED2 73
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       TRED2 74
         DO 120 K = 1, L                                                TRED2 75
  120    SCALE = SCALE + DABS(Z(I,K))                                   TRED2 7
C                                                                       TRED2 77
         IF (SCALE .NE. 0.D0) GO TO 140                                 TRED2 78
  130    E(I) = Z(I,L)                                                  TRED2 79
         GO TO 290                                                      TRED2 80
C                                                                       TRED2 81
  140    DO 150 K = 1, L                                                TRED2 82
            Z(I,K) = Z(I,K) / SCALE                                     TRED2 83
            H = H + Z(I,K) * Z(I,K)                                     TRED2 84
  150    CONTINUE                                                       TRED2 85
C                                                                       TRED2 86
         F = Z(I,L)                                                     TRED2 87
         G = -SIGN(DSQRT(H),F)                                          TRED2 88
         E(I) = SCALE * G                                               TRED2 89
         H = H - F * G                                                  TRED2 90
         Z(I,L) = F - G                                                 TRED2 91
         F = 0.D0                                                       TRED2 92
C                                                                       TRED2 93
         DO 240 J = 1, L                                                TRED2 94
            Z(J,I) = Z(I,J) / H                                         TRED2 95
            G = 0.D0                                                    TRED2 96
C     .......... FORM ELEMENT OF A*U ..........                         TRED2 97
            DO 180 K = 1, J                                             TRED2 98
  180       G = G + Z(J,K) * Z(I,K)                                     TRED2 99
C                                                                       TRED2100
            JP1 = J + 1                                                 TRED2101
            IF (L .LT. JP1) GO TO 220                                   TRED2102
C                                                                       TRED2103
            DO 200 K = JP1, L                                           TRED2104
  200       G = G + Z(K,J) * Z(I,K)                                     TRED2105
C     .......... FORM ELEMENT OF P ..........                           TRED2106
  220       E(J) = G / H                                                TRED2107
            F = F + E(J) * Z(I,J)                                       TRED2108
  240    CONTINUE                                                       TRED2109
C                                                                       TRED2110
         HH = F / (H + H)                                               TRED2111
C     .......... FORM REDUCED A ..........                              TRED2112
         DO 260 J = 1, L                                                TRED2113
            F = Z(I,J)                                                  TRED2114
            G = E(J) - HH * F                                           TRED2115
            E(J) = G                                                    TRED2116
C                                                                       TRED2117
            DO 260 K = 1, J                                             TRED2118
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)                  TRED2119
  260    CONTINUE                                                       TRED2120
C                                                                       TRED2121
  290    D(I) = H                                                       TRED2122
  300 CONTINUE                                                          TRED2123
C                                                                       TRED2124
  320 D(1) = 0.D0                                                       TRED2125
      E(1) = 0.D0                                                       TRED2126
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........     TRED2127
      DO 500 I = 1, N                                                   TRED2128
         L = I - 1                                                      TRED2129
         IF (D(I) .EQ. 0.D0) GO TO 380                                  TRED2130
C                                                                       TRED2131
         DO 360 J = 1, L                                                TRED2132
            G = 0.D0                                                    TRED2133
C                                                                       TRED2134
            DO 340 K = 1, L                                             TRED2135
  340       G = G + Z(I,K) * Z(K,J)                                     TRED2136
C                                                                       TRED2137
            DO 360 K = 1, L                                             TRED2138
               Z(K,J) = Z(K,J) - G * Z(K,I)                             TRED2139
  360    CONTINUE                                                       TRED2140
C                                                                       TRED2141
  380    D(I) = Z(I,I)                                                  TRED2142
         Z(I,I) = 1.D0                                                  TRED2143
         IF (L .LT. 1) GO TO 500                                        TRED2144
C                                                                       TRED2145
         DO 400 J = 1, L                                                TRED2146
            Z(I,J) = 0.D0                                               TRED2147
            Z(J,I) = 0.D0                                               TRED2148
  400    CONTINUE                                                       TRED2149
C                                                                       TRED2150
  500 CONTINUE                                                          TRED2151
C                                                                       TRED2152
      RETURN                                                            TRED2153
      END                                                               TRED2154
C
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  TQL2   3
C***BEGIN PROLOGUE  TQL2                                                TQL2   4
C***DATE WRITTEN   760101   (YYMMDD)                                    TQL2   5
C***REVISION DATE  830518   (YYMMDD)                                    TQL2   6
C***CATEGORY NO.  D4A5,D4C2A                                            TQL2   7
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK                          TQL2   8
C***AUTHOR  SMITH, B. T., ET AL.                                        TQL2   9
C***PURPOSE  COMPUTE EIGENVALUES AND EIGENVECTORS OF SYMMETRIC          TQL2  10
C            TRIDIAGONAL MATRIX.                                        TQL2  11
C***DESCRIPTION                                                         TQL2  12
C                                                                       TQL2  13
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     TQL2  14
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     TQL2  15
C     WILKINSON.                                                        TQL2  16
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   TQL2  17
C                                                                       TQL2  18
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            TQL2  19
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               TQL2  20
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              TQL2  21
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  TQL2  22
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  TQL2  23
C                                                                       TQL2  24
C     ON INPUT                                                          TQL2  25
C                                                                       TQL2  26
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TQL2  27
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TQL2  28
C          DIMENSION STATEMENT.                                         TQL2  29
C                                                                       TQL2  30
C        N IS THE ORDER OF THE MATRIX.                                  TQL2  31
C                                                                       TQL2  32
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.          TQL2  33
C                                                                       TQL2  34
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        TQL2  35
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.               TQL2  36
C                                                                       TQL2  37
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           TQL2  38
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      TQL2  39
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        TQL2  40
C          THE IDENTITY MATRIX.                                         TQL2  41
C                                                                       TQL2  42
C      ON OUTPUT                                                        TQL2  43
C                                                                       TQL2  44
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          TQL2  45
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          TQL2  46
C          UNORDERED FOR INDICES 1,2,...,IERR-1.                        TQL2  47
C                                                                       TQL2  48
C        E HAS BEEN DESTROYED.                                          TQL2  49
C                                                                       TQL2  50
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           TQL2  51
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     TQL2  52
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       TQL2  53
C          EIGENVALUES.                                                 TQL2  54
C                                                                       TQL2  55
C        IERR IS SET TO                                                 TQL2  56
C          ZERO       FOR NORMAL RETURN,                                TQL2  57
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               TQL2  58
C                     DETERMINED AFTER 30 ITERATIONS.                   TQL2  59
C                                                                       TQL2  60
C     CALLS PYTHAG(A,B) FOR SQRT(A**2 + B**2).                          TQL2  61
C                                                                       TQL2  62
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        TQL2  63
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         TQL2  64
C     ------------------------------------------------------------------TQL2  65
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW, TQL2  66
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-    TQL2  67
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,    TQL2  68
C                 1976.                                                 TQL2  69
C***ROUTINES CALLED  PYTHAG                                             TQL2  70
C***END PROLOGUE  TQL2                                                  TQL2  71
C                                                                       TQL2  72
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR                          TQL2  73
      DIMENSION D(N),E(N),Z(N,N)                                        TQL2  74
C     REAL D(N),E(N),Z(NM,N)                                            TQL2  74
C     REAL B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2                             TQL2  75
C     REAL PYTHAG                                                       TQL2  76
C                                                                       TQL2  77
C***FIRST EXECUTABLE STATEMENT  TQL2                                    TQL2  78
      IERR = 0                                                          TQL2  79
      IF (N .EQ. 1) GO TO 1001                                          TQL2  80
C                                                                       TQL2  81
      ONE=1.D0
      DO 100 I = 2, N                                                   TQL2  82
  100 E(I-1) = E(I)                                                     TQL2  83
C                                                                       TQL2  84
      F = 0.D0                                                          TQL2  85
      B = 0.D0                                                          TQL2  86
      E(N) = 0.D0                                                       TQL2  87
C                                                                       TQL2  88
      DO 240 L = 1, N                                                   TQL2  89
         J = 0                                                          TQL2  90
         H = DABS(D(L)) + DABS(E(L))                                    TQL2  91
         IF (B .LT. H) B = H                                            TQL2  92
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........         TQL2  93
         DO 110 M = L, N                                                TQL2  94
            IF (B + DABS(E(M)) .EQ. B) GO TO 120                        TQL2  95
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               TQL2  96
C                THROUGH THE BOTTOM OF THE LOOP ..........              TQL2  97
  110    CONTINUE                                                       TQL2  98
C                                                                       TQL2  99
  120    IF (M .EQ. L) GO TO 220                                        TQL2 100
  130    IF (J .EQ. 30) GO TO 1000                                      TQL2 101
         J = J + 1                                                      TQL2 102
C     .......... FORM SHIFT ..........                                  TQL2 103
         L1 = L + 1                                                     TQL2 104
         L2 = L1 + 1                                                    TQL2 105
         G = D(L)                                                       TQL2 106
         P = (D(L1) - G) / (2.D0 * E(L))                                TQL2 107
C        R = PYTHAG(P,1.D0)                                             TQL2 108
           PP=DMAX1(DABS(P),ONE)                                        PYTHAG14
           QQ=DMIN1(DABS(P),ONE)                                        PYTHAG15
           IF (QQ.EQ.0.D0) GO TO 20                                     PYTHAG16
   10      CONTINUE                                                     PYTHAG17
           RR=(QQ/PP)**2                                                PYTHAG18
           TT=4.D0+RR                                                   PYTHAG19
           IF (TT.EQ.4.D0) GO TO 20                                     PYTHAG20
           SS=RR/TT                                                     PYTHAG21
           PP=PP+2.D0*PP*SS                                             PYTHAG22
           QQ=QQ*SS                                                     PYTHAG23
           GO TO 10                                                     PYTHAG24
   20    R=PP                                                           PYTHAG25
         D(L) = E(L) / (P + SIGN(R,P))                                  TQL2 109
         D(L1) = E(L) * (P + SIGN(R,P))                                 TQL2 110
         DL1 = D(L1)                                                    TQL2 111
         H = G - D(L)                                                   TQL2 112
         IF (L2 .GT. N) GO TO 145                                       TQL2 113
C                                                                       TQL2 114
         DO 140 I = L2, N                                               TQL2 115
  140    D(I) = D(I) - H                                                TQL2 116
C                                                                       TQL2 117
  145    F = F + H                                                      TQL2 118
C     .......... QL TRANSFORMATION ..........                           TQL2 119
         P = D(M)                                                       TQL2 120
         C = 1.D0                                                       TQL2 121
         C2 = C                                                         TQL2 122
         EL1 = E(L1)                                                    TQL2 123
         S = 0.D0                                                       TQL2 124
         MML = M - L                                                    TQL2 125
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........             TQL2 126
         DO 200 II = 1, MML                                             TQL2 127
            C3 = C2                                                     TQL2 128
            C2 = C                                                      TQL2 129
            S2 = S                                                      TQL2 130
            I = M - II                                                  TQL2 131
            G = C * E(I)                                                TQL2 132
            H = C * P                                                   TQL2 133
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150                      TQL2 134
            C = E(I) / P                                                TQL2 135
            R = DSQRT(C*C+1.D0)                                         TQL2 136
            E(I+1) = S * P * R                                          TQL2 137
            S = C / R                                                   TQL2 138
            C = 1.D0 / R                                                TQL2 139
            GO TO 160                                                   TQL2 140
  150       C = P / E(I)                                                TQL2 141
            R = DSQRT(C*C+1.D0)                                         TQL2 142
            E(I+1) = S * E(I) * R                                       TQL2 143
            S = 1.D0 / R                                                TQL2 144
            C = C * S                                                   TQL2 145
  160       P = C * D(I) - S * G                                        TQL2 146
            D(I+1) = H + S * (C * G + S * D(I))                         TQL2 147
C     .......... FORM VECTOR ..........                                 TQL2 148
            DO 180 K = 1, N                                             TQL2 149
               H = Z(K,I+1)                                             TQL2 150
               Z(K,I+1) = S * Z(K,I) + C * H                            TQL2 151
               Z(K,I) = C * Z(K,I) - S * H                              TQL2 152
  180       CONTINUE                                                    TQL2 153
C                                                                       TQL2 154
  200    CONTINUE                                                       TQL2 155
C                                                                       TQL2 156
         P = -S * S2 * C3 * EL1 * E(L) / DL1                            TQL2 157
         E(L) = S * P                                                   TQL2 158
         D(L) = C * P                                                   TQL2 159
         IF (B + DABS(E(L)) .GT. B) GO TO 130                           TQL2 160
  220    D(L) = D(L) + F                                                TQL2 161
  240 CONTINUE                                                          TQL2 162
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........          TQL2 163
      DO 300 II = 2, N                                                  TQL2 164
         I = II - 1                                                     TQL2 165
         K = I                                                          TQL2 166
         P = D(I)                                                       TQL2 167
C                                                                       TQL2 168
         DO 260 J = II, N                                               TQL2 169
            IF (D(J) .GE. P) GO TO 260                                  TQL2 170
            K = J                                                       TQL2 171
            P = D(J)                                                    TQL2 172
  260    CONTINUE                                                       TQL2 173
C                                                                       TQL2 174
         IF (K .EQ. I) GO TO 300                                        TQL2 175
         D(K) = D(I)                                                    TQL2 176
         D(I) = P                                                       TQL2 177
C                                                                       TQL2 178
         DO 280 J = 1, N                                                TQL2 179
            P = Z(J,I)                                                  TQL2 180
            Z(J,I) = Z(J,K)                                             TQL2 181
            Z(J,K) = P                                                  TQL2 182
  280    CONTINUE                                                       TQL2 183
C                                                                       TQL2 184
  300 CONTINUE                                                          TQL2 185
C                                                                       TQL2 186
      GO TO 1001                                                        TQL2 187
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      TQL2 188
C                EIGENVALUE AFTER 30 ITERATIONS ..........              TQL2 189
 1000 IERR = L                                                          TQL2 190
 1001 RETURN                                                            TQL2 191
      END                                                               TQL2 192
