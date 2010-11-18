DECLARE SUB cfhamilton (j!, gJ!, b#(), hx!, hy!, hz!, pxy!, pxz!, pyz!, hr#(), hc#())
DECLARE SUB cfpar (b#(), q!, xNN#, yNN#, zNN#, site%)
DECLARE SUB diagonalize (hr#(), hc#(), D%, en#(), cr#(), cc#())
DECLARE SUB htribk (NM%, N%, ar#(), ai#(), tau#(), m%, zR#(), ZI#())
DECLARE SUB htridi (NM%, N%, ar#(), ai#(), D#(), e#(), E2#(), tau#())
DECLARE SUB imtql2 (NM%, N%, D#(), e#(), z#(), Ierr%)
DECLARE SUB susz (T!, g!, j!, cr#(), cc#(), en#(), chi#())
DECLARE SUB gsmoment (j!, gJ!, cr#(), cc#(), ewr#)
DECLARE SUB transprob (j, cr#(), cc#(), ent!, i!, px!, py!, pz!)
DECLARE SUB detber (ar#(), ac#(), rang%, det#)
DECLARE SUB prtewuev (en#(), cr#(), cc#(), j!)
DECLARE SUB termerww (hr#(), hc#(), D%, T, en#(), cr#(), cc#(), tew#)
DECLARE FUNCTION delta (x)
DECLARE FUNCTION multr# (ar#, ac#, br#, bc#)
DECLARE FUNCTION divc# (ar#, ac#, br#, bc#)
DECLARE FUNCTION divr# (ar#, ac#, br#, bc#)
DECLARE FUNCTION multc# (ar#, ac#, br#, bc#)
DECLARE SUB svlmadd (x#, y#, z#, r#)
  DEFDBL A-D, R-S, X-Z
  DEF FNDELTA (x) = 1 - SGN(x) * SGN(x)
  DEF FNV20 (z, r) = .25 * SQR(5 / pi#) * (3 * z * z - r * r) / r / r
  DEF FNV22 (x, y, z) = 1 / 4 * SQR(15 / pi#) * (x * x - y * y) / (x * x + y * y + z * z)
  DEF FNV40 (z, r) = 3 / 16 * SQR(1 / pi#) * (35 * z * z * z * z - 30 * z * z * r * r + 3 * r * r * r * r) / r / r / r / r
  DEF FNV42 (x, y, z) = 3 / 8 * SQR(5 / pi#) * (7 * z * z - (x * x + y * y + z * z)) * (x * x - y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV43 (x, y, z) = 3 / 8 * SQR(70 / pi#) * z * x * (x * x - 3 * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV44 (x, y, z) = 3 / 16 * SQR(35 / pi#) * (x * x * x * x - 6 * x * x * y * y + y * y * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV60 (z, r) = 1 / 32 * SQR(13 / pi#) * (231 * z * z * z * z * z * z - 315 * z * z * z * z * r * r + 105 * z * z * r * r * r * r - 5 * r * r * r * r * r * r) / r / r / r / r / r / r
  DEF FNV62 (x, y, z) = 1 / 64 * SQR(2730 / pi#) * (16 * z * z * z * z - 16 * (x * x + y * y) * z * z + (x * x + y * y) * (x * x + y * y)) * (x * x - y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV63 (x, y, z) = 1 / 32 * SQR(2730 / pi#) * z * x * (11 * z * z - 3 * (x * x + y * y + z * z)) * (x * x - 3 * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV64 (x, y, z) = 21 / 32 * SQR(13 / 7 / pi#) * (11 * z * z - x * x - y * y - z * z) * (x * x * x * x - 6 * x * x * y * y + y * y * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV66 (x, y, z) = 231 / 64 * SQR(26 / 231 / pi#) * (x * x * x * x * x * x - 15 * x * x * x * x * y * y + 15 * x * x * y * y * y * y - y * y * y * y * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV66S (x, y, z) = 231 / 32 * SQR(26 / 231 / pi#) * y * x * (3 * x * x - y * y) * (x * x - 3 * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)

DIM hr#(20, 20), hc#(30, 30), cr#(20, 20), cc#(20, 20), energy(10)
DIM b#(6, 6), en#(100), nul#(6, 6), chi(3, 3)

PRINT "POINTC POINTC POINTC POINTC POINTC POINTC POINTC POINTC POINTC POINTC"
IF LTRIM$(COMMAND$) = "" GOTO 333
1 DATA "*************************************************************************"
 DATA "Program to calculate Crystalfield Parameters from Point Charges"
 DATA "Usage: pointc 0.2 Ce3+ 4 1 5.3"
 DATA " ... meaning calculate Blms for one pointcharge of 0.2|e| in distance"
 DATA " x=4 A y=1 A z=5.3 A from a Ce3+ ion, results written to pointc.out"
  DATA "*************************************************************************"
 pi# = 3.141592654#

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
A$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
q = VAL(RTRIM$(LEFT$(A$, INSTR(A$, " ")))): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - INSTR(A$, " ")))
ion$ = LCASE$(LTRIM$(RTRIM$(LEFT$(A$, INSTR(A$, " "))))): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - INSTR(A$, " ")))
xcharge# = VAL(RTRIM$(LEFT$(A$, INSTR(A$, " ")))): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - INSTR(A$, " ")))
ycharge# = VAL(RTRIM$(LEFT$(A$, INSTR(A$, " ")))): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - INSTR(A$, " ")))
zcharge# = VAL(A$)
OPEN "o", 1, "pointc.out"
PRINT "Calculating: point charge "; q; "|e| at position ("; xcharge#; "A /"; ycharge#; " A / "; zcharge#; "A ) of "; ion$
PRINT #1, "{point charge "; q; "|e| at position ("; xcharge#; "A /"; ycharge#; " A / "; zcharge#; "A ) of "; ion$
REM site% legt das ion welches berechnet wird fest

IF ion$ = "ce3+" THEN site% = 1: j = 5 / 2: gJ = 6 / 7: GOTO 2 ' drehimpulsquantenzahl J,landefaktor gj ce3+
IF ion$ = "nd3+" THEN site% = 2: j = 9 / 2: gJ = 8 / 11: GOTO 2' drehimpulsquantenzahl J,landefaktor gj nd3+
IF ion$ = "er3+" THEN site% = 3: j = 15 / 2: gJ = 6 / 5: GOTO 2' drehimpulsquantenzahl J,landefaktor gj er3+
IF ion$ = "tm3+" THEN site% = 4: j = 6: gJ = 7 / 6: GOTO 2'drehimpulsquantenzahl J,landefaktor gj tm3+
IF ion$ = "pr3+" THEN site% = 5: j = 4: gJ = 4 / 5: GOTO 2' drehimpulsquantenzahl J,landefaktor gj pr3+
IF ion$ = "tb3+" THEN site% = 6: j = 6: gJ = 3 / 2: GOTO 2' drehimpulsquantenzahl J,landefaktor gj pr3+
IF ion$ = "u3+" THEN site% = 7: j = 4: gJ = 4 / 5: GOTO 2' drehimpulsquantenzahl J,landefaktor gj pr3+
IF ion$ = "ho3+" THEN site% = 8: j = 8: gJ = 5 / 4: GOTO 2' drehimpulsquantenzahl J,landefaktor gj pr3+
IF ion$ = "dy3+" THEN site% = 9: j = 15 / 2: gJ = 4 / 3: GOTO 2' drehimpulsquantenzahl J,landefaktor gj pr3+
PRINT "Error program Pointc: ion "; ion$; " not implemented": END


2 CALL cfpar(b#(), q, xcharge#, ycharge#, zcharge#, site%)
PRINT #1, j; "{momentum J}"
  PRINT #1, gJ; "{landefactor gJ}"
'RINT #1, "No{   vs Blm[meV]   }"
PRINT #1, b#(2, 0); "1  {B20[meV]}"
PRINT #1, b#(2, 2); "2  {B22[meV]}"
PRINT #1, b#(4, 0); "3  {B40[meV]}"
PRINT #1, b#(4, 2); "4  {B42[meV]}"
PRINT #1, b#(4, 3); "5  {B43[meV]}"
PRINT #1, b#(4, 4); "6  {B44[meV]}"
PRINT #1, b#(6, 0); "7  {B60[meV]}"
PRINT #1, b#(6, 2); "8  {B62[meV]}"
PRINT #1, b#(6, 3); "9  {B63[meV]}"
PRINT #1, b#(6, 4); "10 {B64[meV]}"
PRINT #1, b#(6, 6); "11 {B66[meV]}"
CLOSE 1
END

REM ****************check und kontrolle des parametersatzes******************************
'CALL cfpar(b#(), q(2, 0), q(2, 2), q(4, 0), q(4, 2), q(4, 3), q(4, 4), q(6, 0), q(6, 2), q(6, 3), q(6, 4), q(6, 6), site%)
'CALL cfhamilton(j, gJ, b#(), hx, hy, hz, -dpxy, -dpxz, -dpyz, hr#(), hc#())
'CALL diagonalize(hr#(), hc#(), 2 * j + 1, en#(), cr#(), cc#())
'CALL prtewuev(en#(), cr#(), cc#(), j)
'T = 4: CALL susz(T, gJ, j, cr#(), cc#(), en#(), chi#())
END
333 FOR i = 1 TO 6: READ A$: PRINT A$: NEXT i: END

  
    A$ = INKEY$
    IF A$ <> "" THEN
    IF ASC(A$) = 27 THEN END
    END IF

END

DEFSNG A-D, R-S, X-Z
SUB cfhamilton (j, gJ, b#(), hx, hy, hz, pxy, pxz, pyz, hr#(), hc#())

 REM ********************************************************************
 REM diese routine dient zur berechnung des kristallfeldhamiltonoperators
 REM eingabe:
 REM J............drehimpulsquantenzahl
 REM gJ...........Landefaktor(nur fuer H<>0 von bedeutung
 REM b#(l,m)......kristallfeldpar.in mev(implementier:l gerade,m>0,m<>1,5)
 REM hx,hy,hz.....magnetfeld in tesla
 REM pxy,pxz,pyz..quadrupolfeld in meV
 REM augabe:
 REM hr#(1..2J+1,1..2J+1),hc#(1..2J+1,1..2J+1)........real- und imaginaerteil
 REM der kristallfeldmatrix <m|Hcf+Hze|n>=hr#(m+J+1,n+J+1)+ihc#(m+J+1,n+J+1)
 REM ********************************************************************

 DIM o#(6, 6)
 mb# = 5.788378E-02: REM Bohrmagneton in meV/tesla
REM PRINT "m'= "; : FOR x = -J TO J: PRINT USING "###.##"; x; : NEXT x: PRINT
REM PRINT "m= ";
 FOR x = -j TO j
 REM PRINT : PRINT USING "##.#"; x;
  FOR y = -j TO j

   REM addieren der einzelnen beitrÑge v(m,n).<J,x:o(m,n):J,y>
   hr#(x + j + 1, y + j + 1) = 0: hc#(x + j + 1, y + j + 1) = 0

   GOSUB 5000


REM  IF ABS(hr#(x + J + 1, y + J + 1)) > .001 THEN PRINT USING "###.##"; hr#(x + J + 1, y + J + 1); :     ELSE PRINT "  0   ";
  NEXT y
 REM PRINT : PRINT "    ";
REM  FOR y = -J TO J
REM   IF ABS(hc#(x + J + 1, y + J + 1)) > .001 THEN PRINT " i"; : PRINT USING "##.#"; hc#(x + J + 1, y + J + 1); :     ELSE PRINT "      ";
REM  NEXT y
  NEXT x
REM PRINT
 GOTO 281

5000 REM Berechnung von <Jx:o(m,n):Jy>=o#(m,n)
5010 REM fÅr geringe Symmetrie sind eventuell zusÑtzliche
5020 REM stevenson-Operatoren hier hinzuzufÅgen und am Programmanfang
5030 REM die entsprechenden b#(m,n) zu ergÑnzen
5040 REM der gewîhnliche Ausdruck fÅr einen Stevenson-op. wird hier
5050 REM so verÑndert,dass direkt statt des operators jz entweder die
5060 REM zahl x od y bzw. statt j+ und j- entsprechende zahlen in der
5070 REM formel vorkommen.
5080 REM
 REM ********************************************************************
 REM berechnung der den potenzen von j+ und j- entsprechenden zahlen
 REM *****     Formeln fÅr die Stevenson - Operatoren              ******
 REM                       Summation
 REM ********************************************************************
IF x = y THEN
 o#(2, 0) = (3 * x * x - j * (j + 1))
 o#(4, 0) = (35 * x * x * x * x - 30 * j * (j + 1) * x * x + 25 * x * x - 6 * j * (j + 1) + 3 * j * j * (j + 1) * (j + 1))
 o#(6, 0) = (231 * x * x * x * x * x * x - 315 * j * (j + 1) * x * x * x * x + 735 * x * x * x * x + 105 * j * j * (j + 1) * (j + 1) * x * x - 525 * j * (j + 1) * x * x + 294 * x * x - 5 * j * j * j * (j + 1) * (j + 1) * (j + 1) + 40 * j * j * (j +  _
1) * (j + 1) - 60 * j * (j + 1))
 FOR m = 2 TO 6 STEP 2
   hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(m, 0) * o#(m, 0)
    REM da in einem array nur positive argumente zugelassen sind, entspricht
    REM hcf(m,m')..HCF#(X+J+1,Y+J+1) mit x,y=m,m'!!!!!
 NEXT m
 REM Zeemannterm
  hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) - gJ * mb# * hz * x
END IF

IF x = y + 1 THEN
 JP1# = SQR((j - y) * (j + y + 1))
 Jx# = JP1# / 2: REM reell
 Jy# = -JP1# / 2: REM imaginaer
   REM Zeemannterm
     hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) - gJ * mb# * hx * Jx#
     hc#(x + j + 1, y + j + 1) = hc#(x + j + 1, y + j + 1) - gJ * mb# * hy * Jy#
   REM pxy,pxz,pyz
 opxzr# = (x * Jx# + Jx# * y) / 2
 opyzc# = (x * Jy# + Jy# * y) / 2
     hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + pxz * opxzr#
     hc#(x + j + 1, y + j + 1) = hc#(x + j + 1, y + j + 1) + pyz * opyzc#
END IF

IF x = y - 1 THEN
 JM1# = SQR((j + y) * (j - y + 1))
 Jx# = JM1# / 2: REM reell
 Jy# = JM1# / 2: REM imaginaer
   REM Zeemannterm
     hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) - gJ * mb# * hx * Jx#
     hc#(x + j + 1, y + j + 1) = hc#(x + j + 1, y + j + 1) - gJ * mb# * hy * Jy#
   REM pxy,pxz,pyz
 opxzr# = (x * Jx# + Jx# * y) / 2
 opyzc# = (x * Jy# + Jy# * y) / 2
     hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + pxz * opxzr#
     hc#(x + j + 1, y + j + 1) = hc#(x + j + 1, y + j + 1) + pyz * opyzc#
END IF

IF x = y + 2 THEN
JP2# = SQR((j - y) * (j + y + 1) * (j - y - 1) * (j + y + 2))
   REM pxy
 opxyc# = -JP2# / 4
 hc#(x + j + 1, y + j + 1) = hc#(x + j + 1, y + j + 1) + pxy * opxyc#
 o#(2, 2) = JP2# / 2
 o#(4, 2) = ((7 * x * x - j * (j + 1) - 5) * JP2# + JP2# * (7 * y * y - j * (j + 1) - 5)) / 4
 o#(6, 2) = ((33 * x * x * x * x - (18 * j * (j + 1) + 123) * x * x + j * j * (j + 1) * (j + 1) + 10 * j * (j + 1) + 102) * JP2# + JP2# * (33 * y * y * y * y - (18 * j * (j + 1) + 123) * y * y + j * j * (j + 1) * (j + 1) + 10 * j * (j + 1) + 102)) / _
 4
   FOR m = 2 TO 6 STEP 2
    hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(m, 2) * o#(m, 2)
    REM da in einem array nur positive argumente zugelassen sind, entspricht
    REM hcf(m,m')..HCF#(X+J+1,Y+J+1) mit x,y=m,m'!!!!!
   NEXT m
  END IF

IF x = y - 2 THEN
JM2# = SQR((j + y) * (j - y + 1) * (j + y - 1) * (j - y + 2))
   REM pxy
 opxyc# = JM2# / 4
 hc#(x + j + 1, y + j + 1) = hc#(x + j + 1, y + j + 1) + pxy * opxyc#
 o#(2, 2) = JM2# / 2
 o#(4, 2) = ((7 * x * x - j * (j + 1) - 5) * JM2# + JM2# * (7 * y * y - j * (j + 1) - 5)) / 4
 o#(6, 2) = ((33 * x * x * x * x - (18 * j * (j + 1) + 123) * x * x + j * j * (j + 1) * (j + 1) + 10 * j * (j + 1) + 102) * JM2# + JM2# * (33 * y * y * y * y - (18 * j * (j + 1) + 123) * y * y + j * j * (j + 1) * (j + 1) + 10 * j * (j + 1) + 102)) / _
 4
   FOR m = 2 TO 6 STEP 2
    hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(m, 2) * o#(m, 2)
    REM da in einem array nur positive argumente zugelassen sind, entspricht
    REM hcf(m,m')..HCF#(X+J+1,Y+J+1) mit x,y=m,m'!!!!!
   NEXT m
END IF

IF x = y + 3 THEN
JP3# = SQR((j - y) * (j - y - 1) * (j - y - 2) * (j + y + 1) * (j + y + 2) * (j + y + 3))
o#(4, 3) = (x * JP3# + JP3# * y) / 4
o#(6, 3) = ((11 * x * x * x - 3 * j * (j + 1) * x - 59 * x) * JP3# + JP3# * (11 * y * y * y - 3 * j * (j + 1) * y - 59 * y)) / 4
 hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(4, 3) * o#(4, 3)
 hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(6, 3) * o#(6, 3)
END IF

IF x = y - 3 THEN
JM3# = SQR((j + y) * (j + y - 1) * (j + y - 2) * (j - y + 1) * (j - y + 2) * (j - y + 3))
o#(4, 3) = (x * JM3# + JM3# * y) / 4
o#(6, 3) = ((11 * x * x * x - 3 * j * (j + 1) * x - 59 * x) * JM3# + JM3# * (11 * y * y * y - 3 * j * (j + 1) * y - 59 * y)) / 4
 hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(4, 3) * o#(4, 3)
 hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(6, 3) * o#(6, 3)
END IF

IF x = y + 4 THEN
JP4# = SQR((j - y) * (j - y - 1) * (j - y - 2) * (j - y - 3) * (j + y + 1) * (j + y + 2) * (j + y + 3) * (j + y + 4))
 o#(4, 4) = JP4# / 2
 o#(6, 4) = ((11 * x * x - j * (j + 1) - 38) * JP4# + JP4# * (11 * y * y - j * (j + 1) - 38)) / 4
  hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(4, 4) * o#(4, 4)
  hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(6, 4) * o#(6, 4)
END IF

IF x = y - 4 THEN
JM4# = SQR((j + y) * (j + y - 1) * (j + y - 2) * (j + y - 3) * (j - y + 1) * (j - y + 2) * (j - y + 3) * (j - y + 4))
 o#(4, 4) = JM4# / 2
 o#(6, 4) = ((11 * x * x - j * (j + 1) - 38) * JM4# + JM4# * (11 * y * y - j * (j + 1) - 38)) / 4
  hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(4, 4) * o#(4, 4)
  hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(6, 4) * o#(6, 4)
END IF

IF x = y + 6 THEN
JP6# = SQR((j - y) * (j - y - 1) * (j - y - 2) * (j - y - 3) * (j - y - 4) * (j - y - 5) * (j + y + 1) * (j + y + 2) * (j + y + 3) * (j + y + 4) * (j + y + 5) * (j + y + 6))
 o#(6, 6) = JP6# / 2
  hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(6, 6) * o#(6, 6)
END IF

IF x = y - 6 THEN
JM6# = SQR((j + y) * (j + y - 1) * (j + y - 2) * (j + y - 3) * (j + y - 4) * (j + y - 5) * (j - y + 1) * (j - y + 2) * (j - y + 3) * (j - y + 4) * (j - y + 5) * (j - y + 6))
 o#(6, 6) = JM6# / 2
  hr#(x + j + 1, y + j + 1) = hr#(x + j + 1, y + j + 1) + b#(6, 6) * o#(6, 6)
END IF

5810 RETURN

281 END SUB

DEFDBL A-D, R-S, X-Z
SUB cfpar (b#(), q, xNN, yNN, zNN, site%)
REM************************************************************************
REM sub zur berechnung der kristallfeldparameter b#()
REM eingabe
REM q     ...........parameter - reduced charges [|e|]
'   xNN,yNN,zNN ..... position of charge [A]
REM ausgabe
REM b#(l,m).........kristallfeldparameter
REM************************************************************************
REM******berechnung der b#() aus den reduced charges q..*******
pi# = 3.141592654#

REM clear b#(),svlm
FOR l% = 0 TO 6: FOR m% = 0 TO l%: b#(l%, m%) = 0: NEXT m%: NEXT l%
sv20 = 0: sv22 = 0: sv40 = 0: sv42 = 0: sv43 = 0: sv44 = 0: sv60 = 0: sv62 = 0: sv63 = 0: sv64 = 0: sv66 = 0

REM berechnung der Entwicklungskoeffizienten svlm der Zlm
xR = 0: yR = 0: zR = 0' .... position of rare earth

 GOSUB 66
REM umrechnung in Koeffizienten von Olm
     REM <r^n>-werte in bohrradius^n

IF site% = 1 THEN r2 = 1.2: r4 = 3.455: r6 = 21.226: REM ce3+
IF site% = 2 THEN r2 = 1.001: r4 = 2.401: r6 = 12.396: REM nd3+
IF site% = 3 THEN r2 = .666: r4 = 1.126: r6 = 3.987: REM er3+
IF site% = 4 THEN r2 = .6372: r4 = 1.041: r6 = 3.5913: REM tm3+
IF site% = 5 THEN r2 = 1.086: r4 = 2.822: r6 = 15.726: REM Pr3+
IF site% = 6 THEN r2 = .822: r4 = 1.651: r6 = 6.852: REM Tb3+ values from FREEMAN/DESCLAUX JMMM79
IF site% = 7 THEN r2 = 2.028: r4 = 8.4: r6 = 57.818: REM u3+
IF site% = 8 THEN r2 = .745: r4 = 1.379: r6 = 5.379: REM ho3+ Freeman desclaux JMMM 12(1979) 11
IF site% = 9 THEN r2 = .726: r4 = 1.322: r6 = 5.102: REM dy3+


PRINT #1, "<r^2>="; r2; " a0^2  <r^4>="; r4; " a0^4  <r^6>="; r6; " a0^6    a0=0.5292 Angstroem}"

e = 4.80325E-10: REM elementarladung
     REM stevensfaktoren
IF site% = 1 THEN alpha = -2 / (5 * 7): beta = 2 / 9 / 5 / 7: gamma = 0: REM ce3+
IF site% = 2 THEN alpha = -7 / (9 * 121): beta = -8 * 17 / 27 / 11 / 121 / 13: gamma = -5 * 17 * 19 / 27 / 7 / 121 / 11 / 13 / 13: REM Nd3+
IF site% = 3 THEN alpha = 4 / (9 * 25 * 7): beta = 2 / 9 / 35 / 11 / 13: gamma = 8 / 27 / 7 / 121 / 13 / 13: REM er3+
IF site% = 4 THEN alpha = 1 / (9 * 11): beta = 8 / 81 / 5 / 121: gamma = -5 / 81 / 7 / 121 / 13: REM tm3+
IF site% = 5 THEN alpha = -4 * 13 / (9 * 25 * 11): beta = -4 / (9 * 5 * 121): gamma = 16 * 17 / 81 / 5 / 7 / 121 / 13: REM pr3+
IF site% = 6 THEN alpha = -1 / 99: beta = 2 / (27 * 5 * 121): gamma = -1 / 81 / 7 / 121 / 13: REM Tb3+ values from HUTCHINGS
IF site% = 7 THEN alpha = -4 * 13 / (9 * 25 * 11): beta = -4 / (9 * 5 * 121): gamma = 16 * 17 / 81 / 5 / 7 / 121 / 13: REM up3+
IF site% = 8 THEN alpha = -1 / 2 / 9 / 25: beta = -1 / 2 / 3 / 5 / 7 / 11 / 13: gamma = -5 / 27 / 7 / 121 / 13 / 13 'ho3+
IF site% = 9 THEN alpha = -2 / 9 / 5 / 7: beta = -8 / 27 / 5 / 7 / 11 / 13: gamma = 4 / 27 / 7 / 121 / 13 / 13 'dy3+


PRINT #1, alpha; "{alpha Stevens factor}"
PRINT #1, beta; "{beta Stevens factor}"
PRINT #1, gamma; "{gamma Stevens factor}"

     REM einheit von r in <r^n> ist Bohrradius^n = a0^n in angstroem^n
     a0 = .5292#: REM (Angstroem)
   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in mJ = 10^4
   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in THz =1.509166084e22
   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in meV =0.624146e23
   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in K =0.72429024e24
     umr = 6.24146E+22
     ehv2 = a0 ^ 2 * umr
     ehv4 = a0 ^ 4 * umr
     ehv6 = a0 ^ 6 * umr

  b#(2, 0) = -e ^ 2 * r2 * alpha * sv20 * SQR(5 / 16 / pi#) * ehv2
  b#(2, 2) = -e ^ 2 * r2 * alpha * sv22 * SQR(15 / 16 / pi#) * ehv2

  b#(4, 0) = -e ^ 2 * r4 * beta * sv40 * 3 / 16 * SQR(1 / pi#) * ehv4
  b#(4, 2) = -e ^ 2 * r4 * beta * sv42 * 3 / 8 * SQR(5 / pi#) * ehv4
  b#(4, 3) = -e ^ 2 * r4 * beta * sv43 * 3 / 8 * SQR(70 / pi#) * ehv4
  b#(4, 4) = -e ^ 2 * r4 * beta * sv44 * 3 / 16 * SQR(35 / pi#) * ehv4

  b#(6, 0) = -e ^ 2 * r6 * gamma * sv60 / 32 * SQR(13 / pi#) * ehv6
  b#(6, 2) = -e ^ 2 * r6 * gamma * sv62 / 64 * SQR(2730 / pi#) * ehv6
  b#(6, 3) = -e ^ 2 * r6 * gamma * sv63 / 32 * SQR(2730 / pi#) * ehv6
  b#(6, 4) = -e ^ 2 * r6 * gamma * sv64 * 21 / 32 * SQR(13 / 7 / pi#) * ehv6
  b#(6, 6) = -e ^ 2 * r6 * gamma * sv66 * 231 / 64 * SQR(26 / 231 / pi#) * ehv6
 
 
GOTO 77

66 x = xNN - xR: y = yNN - yR: z = zNN - zR
r = SQR(x * x + y * y + z * z)


REM svlmadd
sv20 = sv20 + FNV20(z, r) / r ^ 3 * 4 * pi# / 5 * q     '= gamma20 in cgs
sv22 = sv22 + FNV22(x, y, z) / r ^ 3 * 4 * pi# / 5 * q
sv40 = sv40 + FNV40(z, r) / r ^ 5 * 4 * pi# / 9 * q
sv42 = sv42 + FNV42(x, y, z) / r ^ 5 * 4 * pi# / 9 * q
sv43 = sv43 + FNV43(x, y, z) / r ^ 5 * 4 * pi# / 9 * q
sv44 = sv44 + FNV44(x, y, z) / r ^ 5 * 4 * pi# / 9 * q
sv60 = sv60 + FNV60(z, r) / r ^ 7 * 4 * pi# / 13 * q
sv62 = sv62 + FNV62(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
sv63 = sv63 + FNV63(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
sv64 = sv64 + FNV64(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
sv66 = sv66 + FNV66(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
RETURN

77
END SUB

DEFSNG A-D, R-S, X-Z
FUNCTION delta (x)
delta = 1 - SGN(x) * SGN(x)
END FUNCTION

DEFINT I-N
DEFDBL A, D
SUB detber (ar(), ac(), rang%, det#)
 REM ********************************************************************
 REM subroutine zur berechnung der determinante der matrix ar()+i.ac()
 REM eingabe
 REM ar(1..rang%,1..rang%)+i.ac(1..rang%,1..rang%)...matrix
 REM rang%.........rang der matrix
 REM ausgabe
 REM det#..........determinante der matrix
 REM ********************************************************************

lim4# = .1#: REM untere grenze fuer diagonalelement bei gauss elimination
               REM wenn fuer ein bestimmtes K unweigerlich |a(K,K)|<lim4#
               REM auch wenn man dieses element durch
               REM zeilentauschen zu vergroessern sucht,dann wird gauss
               REM elimination abgebrochen, und der wert dieses diagonal-
               REM elements in det# eingetragen
lim3# = .00000000000001#: REM ist weiters |a(K,K)| sogar <lim3#,wird det#=0 gesetzt

 det# = 1
 FOR K = 1 TO rang%
    IF SQR(ar(K, K) * ar(K, K) + ac(K, K) * ac(K, K)) >= lim4# GOTO 7030
    FOR l = K + 1 TO rang%
      IF SQR(ar(l, K) * ar(l, K) + ac(l, K) * ac(l, K)) > lim4# GOTO 7024
    NEXT l
    REM dieses K.diagonalelement kann auch durch tauschen einer zeile nicht
    REM >lim4# gemacht werden, daher:
    det# = SQR(ar(K, K) * ar(K, K) + ac(K, K) * ac(K, K)) * lim4#
    IF SQR(ar(K, K) * ar(K, K) + ac(K, K) * ac(K, K)) < lim3# THEN det# = 0
    ar(K, K) = 0: ac(K, K) = 0
    GOTO 7200
7024 REM tauschen der K. mit der L. Zeile
    FOR i = K TO rang%
      SAVEA# = ar(K, i): ar(K, i) = ar(l, i): ar(l, i) = SAVEA#
      SAVEA# = ac(K, i): ac(K, i) = ac(l, i): ac(l, i) = SAVEA#
    NEXT i
     det# = -det#

7030 FOR l = K + 1 TO rang%
      IF SQR(ar(l, K) * ar(l, K) + ac(l, K) * ac(l, K)) = 0 GOTO 7100
      FOR i = K + 1 TO rang%
       dr = multr#(ar(K, i), ac(K, i), ar(l, K), ac(l, K))
       dc = multc#(ar(K, i), ac(K, i), ar(l, K), ac(l, K))
       ar(l, i) = ar(l, i) - divr#(dr, dc, ar(K, K), ac(K, K))
       ac(l, i) = ac(l, i) - divc#(dr, dc, ar(K, K), ac(K, K))
      NEXT i
      ar(l, K) = 0: ac(l, K) = 0
7100 NEXT l
7110 NEXT K
 FOR K = 1 TO rang%
  det# = det# * SGN(ar(K, K)) * (ar(K, K) * ar(K, K) + ac(K, K) * ac(K, K))
 NEXT K
 det# = SGN(det#) * SQR(ABS(det#))
7200 END SUB

DEFINT D
DEFDBL C, E-H, P, R-T, Z
SUB diagonalize (hr#(), hc#(), D%, en(), cr#(), cc#())
DIM e(30), E2(30), tau(2, 30)
CALL htridi(D%, D%, hr#(), hc#(), en(), e(), E2(), tau())
FOR j = 1 TO D%
FOR K = 1 TO D%
cr#(j, K) = 0
IF j = K THEN cr#(j, j) = 1
NEXT: NEXT
 CALL imtql2(D%, D%, en(), e(), cr#(), Ierr)
IF Ierr <> 0 THEN PRINT "Error in diagonalise routine": STOP
 CALL htribk(D%, D%, hr#(), hc#(), tau(), D%, cr#(), cc#())
END SUB

DEFSNG A, C-N, P, R-T, Z
FUNCTION divc# (ar#, ac#, br#, bc#)
 divc# = (ac# * br# - ar# * bc#) / (br# * br# + bc# * bc#)
END FUNCTION

FUNCTION divr# (ar#, ac#, br#, bc#)
 divr# = (ar# * br# + ac# * bc#) / (br# * br# + bc# * bc#)
END FUNCTION

DEFDBL A-D, R-S, X-Z
SUB gsmoment (j, gJ, cr#(), cc#(), ewr#)
 DIM nul#(6, 6), hr#(30, 30), hc#(30, 30)
 REM fit momente im grundzustand
 mb# = 5.788378E-02: REM Bohrmagneton in meV/tesla
 nul#(2, 0) = 0: hz = 1 / mb#
 CALL cfhamilton(j, gJ, nul#(), hx, hy, hz, -dpxy, -dpxz, -dpyz, hr#(), hc#())
 nul#(2, 0) = 0: hz = 0
  REM berechnung des erwartungswertes von hr#()+i.hc#() im zustand I
  ewr# = 0: ewc# = 0
  FOR m = 1 TO 2 * j + 1: FOR N = 1 TO 2 * j + 1
   zR# = multr#(cr#(1, m), -cc#(1, m), hr#(m, N), hc#(m, N))
   zc# = multc#(cr#(1, m), -cc#(1, m), hr#(m, N), hc#(m, N))
   ewr# = ewr# + multr#(zR#, zc#, cr#(1, N), cc#(1, N))
   ewc# = ewc# + multc#(zR#, zc#, cr#(1, N), cc#(1, N))
  NEXT N: NEXT m
'PRINT "my in z-Richtung(grundzustand)"; ewr#
END SUB

DEFINT I-N
DEFSNG B, X-Y
DEFDBL E-H, P, T
SUB htribk (NM, N, ar(), ai(), tau(), m, zR(), ZI())
'     -------------------------------------------------------------------
'     This subroutine is a translation of a complex analogue of the
'     algol procedure TRED1, Num. Math. 11, 181-195 (1968) by
'     Martin, Reinsch, and Wilkinson.
'     Handbook for Auto. Comp., Vol.II - Linear Algebra, 212-226 (1971).
'
'     This subroutine forms the eigenvectors of a complex Hermitean
'     matrix by back transforming those of the corresponding real
'     symmetric tridiagonal matrix determined by HTRIDI.
'
'     On input:
'
'          NM must be set to the row dimension of the two-dimensional
'          array parameters as declared in the calling program
'          dimension statement.
'
'          N is the order of the matrix.
'
'          AR and AI contain information about the unitary trans-
'          formations used in the reduction of HTRIDI in their full
'          lower triangles except for the diagonal of AR.
'
'          TAU contains further information about the transformation.
'
'          M is the number of eigenvectors to be back transformed.
'
'          ZR contains the eigenvectors to be back transformed
'          in its first M columns.
'
'     On output:
'
'          ZR and ZI contain the real and imaginary parts,
'          respectively, of the transformed eigenvectors in
'          their first M columns.
'
'     Note that the last component of each returned vector is
'     real and that vector euclidean norms are preserved.
'
'     Qustions and comments should be directed to B.S. Garbow,
'     Applied Mathematics Division, Argonne National Laboratory.
'
'     ------------------------------------------------------------------
'
'     :::::::: Transform the eigenvectors of the real symmetric
'              tridiagonal matrix to those of the Hermitean
'              tridiagonal matrix. ::::::::::::::::::::::::::::
'
      FOR K = 1 TO N
      FOR j = 1 TO m
      ZI(j, K) = -zR(j, K) * tau(2, K)
      zR(j, K) = zR(j, K) * tau(1, K)
50    NEXT: NEXT
      IF N = 1 GOTO 200
'     :::::::::::: RECOVER AND APPLY THE HOUSEHOLDER MATRICES :::::::
      FOR i = 2 TO N
      l = i - 1
      H = ai(i, i)
      IF H = 0! GOTO 140
      FOR j = 1 TO m
      S = 0!
      si = 0!
      FOR K = 1 TO l
      S = S + ar(i, K) * zR(j, K) - ai(i, K) * ZI(j, K)
      si = si + ar(i, K) * ZI(j, K) + ai(i, K) * zR(j, K)
      NEXT
      S = S / H
      si = si / H
      FOR K = 1 TO l
      zR(j, K) = zR(j, K) - S * ar(i, K) - si * ai(i, K)
      ZI(j, K) = ZI(j, K) - si * ar(i, K) + S * ai(i, K)
    NEXT
    NEXT
140 NEXT
'     :::::::::::: LAST LINE OF HTRIBK :::::::::::::::::::::::::::
200 END SUB

DEFDBL M
SUB htridi (NM, N, ar(), ai(), D(), e(), E2(), tau())
'
'
'-----------------------------------------------------------------------
'
'     This subroutine is a translation of a complex analogue of
'     the algol procedure TRED1, Num. Math. 11, 181-195 (1968)
'     by Martin, Reinsch, and Wilkinson.
'     Handbook for Auto. Comp., Vol.II-Linear Algebra, 212-226 (1971).
'
'     This subroutine reduces a complex Hermitean matrix
'     to a real symmetri' tridiagonal matrix using
'     unitary similarity transformations.
'
'     On input:
'
'          NM must be set to the row dimension of two-dimensional
'          array parameters as declared in the calling program
'          dimension statements.
'
'          N is the order of the matrix.
'
'          AR and AI contain the real and imaginary parts,
'          respectively, of the complex Hermitean input matrix.
'          Only the lower triangle of the matrix need be supplied.
'
'     On output:
'
'          AR and AI contain information about the unitary trans-
'          formations used in the reduction in their full lower
'          triangles. Their strict upper triangles and the diagonal
'          of AR (and AI ?) are unaltered.
'
'          D contains the diagonal elements of the tridiagonal matrix.
'
'          E contains the subdiagonal elements of the tridiagonal
'          matrix in its last N-1 positions. E(1) is set to zero.
'
'          E2 contains the squares of the corresponding elements of E.
'          E2 may coincide with E if the squares are not needed.
'
'          TAU contains further information about the transformations
'
'     Arithmeti' is real......
'     Questions and comments should be directed to B.S. Garbow,
'     Applied Mathematics Division, Argonne National Laboratory.
'     ------------------------------------------------------------------
'
      tau(1, N) = 1!
      tau(2, N) = 0!
      FOR i = 1 TO N
100   D(i) = ar(i, i)
      NEXT
      FOR ii = 1 TO N
      i = N + 1 - ii
      l = i - 1
      H = 0!
      scale = 0!
      IF l < 1 GOTO 130
      FOR K = 1 TO l
120   scale = scale + ABS(ar(i, K)) + ABS(ai(i, K))
      NEXT K
      IF scale <> 0! GOTO 141
      tau(1, l) = 1!
      tau(2, l) = 0!
130   e(i) = 0!
      E2(i) = 0!
      GOTO 290
141   FOR K = 1 TO l
      ar(i, K) = ar(i, K) / scale
      ai(i, K) = ai(i, K) / scale
      H = H + ar(i, K) * ar(i, K) + ai(i, K) * ai(i, K)
150   NEXT K
      E2(i) = scale * scale * H
      g = SQR(H)
      e(i) = scale * g
      F = SQR(ar(i, l) * ar(i, l) + ai(i, l) * ai(i, l))
'     ::::::::: FORM NEXT DIAGONAL ELEMENT OF MATRIX T ::::::::::
      IF F = 0! GOTO 160
      tau(1, l) = (ai(i, l) * tau(2, i) - ar(i, l) * tau(1, i)) / F
      si = (ar(i, l) * tau(2, i) + ai(i, l) * tau(1, i)) / F
      H = H + F * g
      g = 1! + g / F
      ar(i, l) = g * ar(i, l)
      ai(i, l) = g * ai(i, l)
      IF l = 1 GOTO 270
      GOTO 170
160   tau(1, l) = -tau(1, i)
      si = tau(2, i)
      ar(i, l) = g
170   F = 0!
      FOR j = 1 TO l
      g = 0!
      GI = 0!
'     :::::::::::: FORM ELEMENT OF A*U ::::::::::::::::::::
      FOR K = 1 TO j
      g = g + ar(j, K) * ar(i, K) + ai(j, K) * ai(i, K)
      GI = GI - ar(j, K) * ai(i, K) + ai(j, K) * ar(i, K)
180 NEXT K
      JP1 = j + 1
      IF l < JP1 GOTO 220
      FOR K = JP1 TO l
      g = g + ar(K, j) * ar(i, K) - ai(K, j) * ai(i, K)
      GI = GI - ar(K, j) * ai(i, K) - ai(K, j) * ar(i, K)
     NEXT K
'     ::::::::::::: FORM ELEMENT OF P ::::::::::::::::::::::
220   e(j) = g / H
      tau(2, j) = GI / H
      F = F + e(j) * ar(i, j) - tau(2, j) * ai(i, j)
240  NEXT j
      HH = F / (H + H)
'     :::::::::::::: FORM REDUCED A ::::::::::
      FOR j = 1 TO l
      F = ar(i, j)
      g = e(j) - HH * F
      e(j) = g
      FI = -ai(i, j)
      GI = tau(2, j) - HH * FI
      tau(2, j) = -GI
      FOR K = 1 TO j
      ar(j, K) = ar(j, K) - F * e(K) - g * ar(i, K) + FI * tau(2, K) + GI * ai(i, K)
      ai(j, K) = ai(j, K) - F * tau(2, K) - g * ai(i, K) - FI * e(K) - GI * ar(i, K)
260  NEXT K: NEXT j
270   FOR K = 1 TO l
      ar(i, K) = scale * ar(i, K)
      ai(i, K) = scale * ai(i, K)
280 NEXT K
      tau(2, l) = -si
290   HH = D(i)
      D(i) = ar(i, i)
      ar(i, i) = HH
      ai(i, i) = scale * scale * H
300   NEXT ii
'     ::::::::::::: LAST LINE OF HTRIDI ::::::::::::::::::::::::::::::

END SUB

DEFSNG A, H, T
SUB imtql2 (NM, N, D(), e(), z(), Ierr)
'
'     --------------------------------------------------------------
'
'     This subroutine is a translation of the algol procedure IMTQL2,
'     Num. Math. 12, 377-383 (1968) by Martin and Wilkinson, as
'     modified in Num. Math. 15, 450 (1970) by Dubrulle, Handbook
'     for Auto. Comp. Vol.II - Linear Algebra, 241-248 (1971).
'
'     This subroutine finds the eigenvalues and the eigenvectors of
'     a symmetri' tridiagonal matrix by the "implicit QL method."
'     The eigenvectors of a full symmetri' matrix can also be found
'     if TRED2 has been used to reduced the full matrix to
'     tridiagonal form.
'
'     On input:
'
'          NM must be set to the row dimension of two-dimensional
'          array parameters as declared in the calling program
'          dimension statement.
'
'          N is the order of the matrix.
'
'          D contains the diagonal elements of the input matrix.
'
'          E contains the subdiagonal elements of the input matrix
'          in its last N-1 positions. E(1) is set to arbitrary.
'
'          Z contains the transformation matrix produced in the
'          reduction by TRED2, if performed. If the eigenvectors
'          of the tridiagonal matrix are desired, Z must contain
'          the identity matrix.                   **************
'          ********************
'     On output:
'
'          D contains the eigenvalues in ascending order. If an
'          error exit is made, the eigenvalues are correct but
'          unordered for indices 1,2,...,IERR-1.
'
'          E has been destroyed.
'
'          Z contains orthonormal eigenvectors of the symmetric
'          tridiagonal (or full) matrix. If an error exit is made
'          Z contains the eigenvectors associated with the stored
' eigenvalues.
'
'          IERR is set to
'               0    for normal return
'               J    if the J'th eigenvalue has not been
'                         determined after 30 iterations.
'
'     Questions and comments should be directed to B.S. Gabow,
'     Applied Mathematics Division, Argonne National Laboratory.
'
'     ------------------------------------------------------------------
'
'     :::: MACHEP is (or should have been) a machine dependent
'          parameter specifying the relative precision of floating
'          point aritmetic.(MACHEP=4*8**(-13) for single precision
'          arithmeticon B6700).  ::::::::::::::::::::::::::::::::::
'
       MACHEP = 1E-11
      Ierr = 0
      IF N = 1 GOTO 1001
      FOR i = 2 TO N
      e(i - 1) = e(i)
      NEXT
      e(N) = 0!
      FOR l = 1 TO N
      j = 0
'     ::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::
105   FOR m = l TO N
      IF m = N GOTO 121
      IF ABS(e(m)) <= MACHEP * (ABS(D(m)) + ABS(D(m + 1))) GOTO 121
110  NEXT
121   P = D(l)
      IF m = l GOTO 241
      IF j = 30 GOTO 1000
      j = j + 1
': : : : : : : : : : : : : : : : : FORM SHIFT: : : : : : : : : : : : : : : : : : :
      g = (D(l + 1) - P) / (2! * e(l))
      r = SQR(g * g + 1!)
      Signum = 1: IF g < 0 THEN Signum = -1
      g = D(m) - P + e(l) / (g + Signum * ABS(r))
      S = 1!
      c = 1!
      P = 0!
      MML = m - l
      FOR ii = 1 TO MML
      i = m - ii
      F = S * e(i)
      b = c * e(i)
      IF ABS(F) < ABS(g) GOTO 151
      c = g / F
      r = SQR(c * c + 1!)
      e(i + 1) = F * r
      S = 1! / r
      c = c * S
      GOTO 161
151   S = F / g
      r = SQR(S * S + 1!)
      e(i + 1) = g * r
      c = 1! / r
      S = S * c
161   g = D(i + 1) - P
      r = (D(i) - g) * S + 2! * c * b
      P = S * r
      D(i + 1) = g + P
      g = c * r - b
'     ::::::::::::::::::  FORM VECTOR :::::::::::::::
      FOR K = 1 TO N
      F = z(i + 1, K)
      z(i + 1, K) = S * z(i, K) + c * F
      z(i, K) = c * z(i, K) - S * F
181   NEXT
201  NEXT
      D(l) = D(l) - P
      e(l) = g
      e(m) = 0!
      GOTO 105
241 NEXT
'     :::::::ORDER EIGENVALUES AND EIGENVECTORS :::::::::::
      FOR ii = 2 TO N
      i = ii - 1
      K = i
      P = D(i)
      FOR j = ii TO N
      IF D(j) >= P GOTO 261
      K = j
      P = D(j)
261  NEXT
      IF K = i GOTO 301
      D(K) = D(i)
      D(i) = P
      FOR j = 1 TO N
      P = z(i, j)
      z(i, j) = z(K, j)
      z(K, j) = P
     NEXT
301  NEXT
      GOTO 1001
'     :::::::::::: SET ERROR--NO CONVERGENCE TO AN
'                   EIGRNVALUE AFTER 30 ITERATIONS :::::::::::
1000  Ierr = l
'     ::::::::::::::: LAST LINE OF ITMQL2 ::::::::::::::::::::::::::

1001 END SUB

DEFSNG C-G, I-N, P, R-S, Z
FUNCTION multc# (ar#, ac#, br#, bc#)
 multc# = ac# * br# + ar# * bc#
END FUNCTION

FUNCTION multr# (ar#, ac#, br#, bc#)
 multr# = ar# * br# - ac# * bc#
END FUNCTION

SUB prtewuev (en#(), cr#(), cc#(), j)
REM**********************************************************************
REM ausgabe der eigenwerte und eigenvektoren
REM eingabe:
REM EN#(1..2J+1).....................eigenwerte
REM cr#(I,1..2J+1)+i.cc#(I,1..2J+1)..eigenvektoren
REM J................................drehimpulsquantenzahl
REM**********************************************************************
DIM un(4), unit$(4)
un(3) = 11.6045036#: un(4) = 1: un(2) = .24179696#: un(1) = 806.5479#
unit$(3) = " (K)": unit$(4) = " (meV)": unit$(2) = " (THz)": unit$(1) = " (/m)"
PRINT "Eigenwerte und Eigenvektoren:"
FOR N = 1 TO 4
 PRINT : PRINT "Eclc";
 FOR i = 1 TO 2 * j + 1: PRINT USING "###.###"; en#(i) * un(N); : NEXT i
 PRINT unit$(N): PRINT "Eclc";
 FOR i = 1 TO 2 * j + 1: PRINT USING "###.###"; (en#(i) - en#(1)) * un(N); : NEXT i
 PRINT unit$(N)
NEXT N
INPUT A$
PRINT "m="
 FOR i = 1 TO 2 * j + 1: PRINT USING "##.#"; i - j - 1;
    FOR l = 1 TO 2 * j + 1
     IF ABS(cr#(l, i)) > .0001 THEN PRINT USING "###.###"; cr#(l, i); :                                                                        ELSE PRINT "  0    ";
    NEXT l: PRINT : PRINT "    ";
    FOR l = 1 TO 2 * j + 1
     IF ABS(cc#(l, i)) > .0001 THEN PRINT " +i"; : PRINT USING "#.##"; cc#(l, i); :                                                                        ELSE PRINT "       ";
    NEXT l: PRINT
 NEXT i

END SUB

DEFDBL A-D, R-S, X-Z
SUB susz (T, g, j, cr#(), cc#(), en#(), chi#())
REM********************************************************************
REM berechnung der suszeptibilitaet chi(1..3,1..3)
REM********************************************************************
 DIM jir#(30, 30, 3), jic#(30, 30, 3)
 kb# = .086173528#: REM boltzmann-konstante in meV/K
 mb# = 5.788378E-02: REM Bohrmagneton in meV/tesla

    REM zustandssumme und ausgleichen von numerischen ungenauigkeiten
    REM in EN#()
    z# = 0: jeff = j
    FOR H = 2 * j + 1 TO 1 STEP -1
     z# = z# + EXP(-en#(H) / kb# / T)
    IF en#(H) - en#(1) > 100 * kb# * T THEN jeff = (H - 1) / 2
    IF ABS(en#(H) - en#(H - 1)) < .001 THEN en#(H) = en#(H - 1)
    NEXT H
jeff = j

 REM********************************************************************
 REM* berechnung <H:JI:L>=JIr#+i.jrc#(H,L,I)
 REM********************************************************************
FOR H = 1 TO 2 * jeff + 1: FOR l = 1 TO 2 * jeff + 1
  jir#(H, l, 1) = 0: jic#(H, l, 1) = 0
  jir#(H, l, 2) = 0: jic#(H, l, 2) = 0
  jir#(H, l, 3) = 0: jic#(H, l, 3) = 0
  FOR x = -j TO j: FOR y = -j TO j
   JP1 = 0: IF x = y + 1 THEN JP1 = SQR((j - y) * (j + y + 1))
   JM1 = 0: IF x = y - 1 THEN JM1 = SQR((j + y) * (j - y + 1))
r# = multr#(cr#(l, y + j + 1), cc#(l, y + j + 1), cr#(H, x + j + 1), -cc#(H, x + j + 1))
c# = multc#(cr#(l, y + j + 1), cc#(l, y + j + 1), cr#(H, x + j + 1), -cc#(H, x + j + 1))
jir#(H, l, 1) = jir#(H, l, 1) + r# / 2 * (JP1 + JM1)
jic#(H, l, 1) = jic#(H, l, 1) + c# / 2 * (JP1 + JM1)
jir#(H, l, 2) = jir#(H, l, 2) + c# / 2 * (JP1 - JM1)
jic#(H, l, 2) = jic#(H, l, 2) - r# / 2 * (JP1 - JM1)
jir#(H, l, 3) = jir#(H, l, 3) + r# * x * delta(x - y)
jic#(H, l, 3) = jic#(H, l, 3) + c# * x * delta(x - y)
  NEXT y: NEXT x
NEXT l: NEXT H


510 REM********************************************************************
    REM* berechnung und ausgabe der suszep-werte
520 REM********************************************************************
FOR m = 1 TO 3
FOR N = 1 TO 3

    REM summation
    chi#(m, N) = 0
730 FOR H = 1 TO 2 * jeff + 1: FOR l = 1 TO 2 * jeff + 1
r# = multr#(jir#(H, l, m), jic#(H, l, m), jir#(l, H, N), jic#(l, H, N))
740 IF en#(H) = en#(l) THEN chi#(m, N) = chi#(m, N) + r# * EXP(-en#(H) / kb# / T) / z# / kb# / T
750 IF en#(H) <> en#(l) THEN chi#(m, N) = chi#(m, N) + r# * (EXP(-en#(H) / kb# / T) - EXP(-en#(l) / kb# / T)) / z# / (en#(l) - en#(H))
760 NEXT l: NEXT H

765 JIm# = 0: jin# = 0
770 FOR H = 1 TO 2 * jeff + 1
 JIm# = JIm# + jir#(H, H, m) * EXP(-en#(H) / kb# / T) / z#
 jin# = jin# + jir#(H, H, N) * EXP(-en#(H) / kb# / T) / z#
    NEXT H

790 chi#(m, N) = chi#(m, N) - JIm# * jin# / kb# / T
NEXT N: NEXT m
     PRINT "single-ion-suszeptibilitaet chi(1..3,1..3) bei "; T; "kelvin(in mb/Tesla):"
   FOR i = 1 TO 3
   FOR K = 1 TO 3
   chi#(i, K) = chi#(i, K) * g * g * mb#
    PRINT USING "####.###########"; chi#(i, K);
   NEXT K: PRINT : NEXT i


END SUB

DEFSNG A-D
SUB svlmadd (x, y, z, r)
SHARED sv20, sv22, sv40, sv42, sv44, sv43, sv60, sv62, sv63, sv64, sv66, pi#, q2, q4, q6

  REM drehung des koordinatensystems um die z-achse um fi#
  REM fi# = 0 / 180 * 3.1415
  REM g# = x: f# = y: x = COS(fi#) * f# - SIN(fi#) * g#: y = SIN(fi#) * f# + COS(fi#) * g#

sv20 = sv20 + FNV20(z, r) / r / r / r * 4 * pi# / 5 * q2
sv22 = sv22 + FNV22(x, y, z) / r / r / r * 4 * pi# / 5 * q2
sv40 = sv40 + FNV40(z, r) / r / r / r / r / r * 4 * pi# / 9 * q4
sv42 = sv42 + FNV42(x, y, z) / r / r / r / r / r * 4 * pi# / 9 * q4
sv43 = sv43 + FNV43(x, y, z) / r / r / r / r / r * 4 * pi# / 9 * q4
sv44 = sv44 + FNV44(x, y, z) / r / r / r / r / r * 4 * pi# / 9 * q4
sv60 = sv60 + FNV60(z, r) / r / r / r / r / r / r / r * 4 * pi# / 13 * q6
sv62 = sv62 + FNV62(x, y, z) / r / r / r / r / r / r / r * 4 * pi# / 13 * q6
sv63 = sv63 + FNV63(x, y, z) / r / r / r / r / r / r / r * 4 * pi# / 13 * q6
sv64 = sv64 + FNV64(x, y, z) / r / r / r / r / r / r / r * 4 * pi# / 13 * q6
sv66 = sv66 + FNV66(x, y, z) / r / r / r / r / r / r / r * 4 * pi# / 13 * q6
END SUB

DEFINT I-N
DEFSNG R-S, X-Z
SUB termerww (hr#(), hc#(), D%, T, en#(), cr#(), cc#(), tew#)
 REM ********************************************************************
 REM diese sub berechnet den thermischen erwartungswert tew# des hermiteschen
 REM operators hr#()+i.hc#() bei der temperatur T
 REM eingabe
 REM hr#(1..d%,1..d%)+i.hc#(1..d%.1..d%)...operatorkomponenten
 REM d%....................................dimension des zustandraums
 REM T.....................................temperatur in kelvin
 REM EN#(1..d%)............................energieeigenwerte in meV
 REM cr#(I,1..d%)+i.cc#(I,1..d%)...........eigenzustandskomponenten (I=1..d%)
 REM ausgabe:
 REM tew#..................................thermischer erwartungswert
 REM ********************************************************************
 kb# = .086173528#: REM boltzmann-konstante in meV/K

 REM berechnung der zustandssumme
 z# = 0
 FOR i = 1 TO D%
  z# = z# + EXP(-en#(i) / kb# / T)
 NEXT i

 REM berechnung des thermischen erwartungswertes
 tew# = 0
 FOR i = 1 TO D%
  REM berechnung des erwartungswertes von hr#()+i.hc#() im zustand I
  ewr# = 0: ewc# = 0
  FOR m = 1 TO D%: FOR N = 1 TO D%
   zR# = multr#(cr#(i, m), -cc#(i, m), hr#(m, N), hc#(m, N))
   zc# = multc#(cr#(i, m), -cc#(i, m), hr#(m, N), hc#(m, N))
   ewr# = ewr# + multr#(zR#, zc#, cr#(i, N), cc#(i, N))
   ewc# = ewc# + multc#(zR#, zc#, cr#(i, N), cc#(i, N))
  NEXT N: NEXT m
  IF ABS(ewc#) > 1E-10 THEN PRINT "fehler termerww: zustand"; i; "erwartungswert nicht reell!"
  tew# = tew# + EXP(-en#(i) / kb# / T) / z# * ewr#
 NEXT i
END SUB

DEFSNG I-N
DEFDBL A-D, R-S, X-Z
SUB transprob (j, cr#(), cc#(), ent, i, px, py, pz)
REM berechnung der uebergangswahrscheinlichkeiten px,py,pz
REM vom Grundzustand in den zustand i
px = 0: py = 0: pz = 0: FOR K = 1 TO ent
 REM transitiponprobability px,py,pz  k-->i
 jzr# = 0: jzc# = 0: jxr# = 0: jxc# = 0: jyr# = 0: jyc# = 0
 FOR x = -j TO j: FOR y = -j TO j
  JP1# = 0: IF x = y + 1 THEN JP1# = SQR((j - y) * (j + y + 1))
  JM1# = 0: IF x = y - 1 THEN JM1# = SQR((j + y) * (j - y + 1))
  jzr# = jzr# + multr(cr#(K, x + j + 1), -cc#(K, x + j + 1), cr#(i, y + j + 1), cc#(i, y + j + 1)) * x * delta(x - y)
  jzc# = jzc# + multc(cr#(K, x + j + 1), -cc#(K, x + j + 1), cr#(i, y + j + 1), cc#(i, y + j + 1)) * x * delta(x - y)
  jxr# = jxr# + multr(cr#(K, x + j + 1), -cc#(K, x + j + 1), cr#(i, y + j + 1), cc#(i, y + j + 1)) * (JP1# + JM1#) / 2
  jxc# = jxc# + multc(cr#(K, x + j + 1), -cc#(K, x + j + 1), cr#(i, y + j + 1), cc#(i, y + j + 1)) * (JP1# + JM1#) / 2
  jyr# = jyr# + multc(cr#(K, x + j + 1), -cc#(K, x + j + 1), cr#(i, y + j + 1), cc#(i, y + j + 1)) * (-JP1# + JM1#) / 2: REM imaginaer
  jyc# = jyc# + multr(cr#(K, x + j + 1), -cc#(K, x + j + 1), cr#(i, y + j + 1), cc#(i, y + j + 1)) * (-JP1# + JM1#) / 2: REM imaginaer
 NEXT: NEXT
 px = px + jxr# * jxr# + jxc# * jxc#
 py = py + jyr# * jyr# + jyc# * jyc#
 pz = pz + jzr# * jzr# + jzc# * jzc#
NEXT K
'PRINT "tr.prob.state"; I; "="; px, py, pz

END SUB

