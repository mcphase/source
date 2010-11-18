DECLARE SUB neutint (lambda!, thetamax!, ovalltemp!, lorenz!, r1!(), r2!(), r3!(), n!, x!(), y!(), z!(), mx!(), my!(), mz!(), sl!(), slcode%(), ff!(), ffcode%(), occ!(), m!, h!(), k!(), l!(), D!(), theta!(), imag!(), ikern!(), sf!(), lpg!())
DECLARE SUB printeln (filename$, unitcell$, lambda!, ovalltemp!, lorenz!, r1!(), r2!(), r3!(), n!, x!(), y!(), z!(), mx!(), my!(), mz!(), sl!(), slcode%(), ff!(), ffcode%(), occ!(), m!, h!(), k!(), l!(), ikern!(), imag!(), D!(), theta!(), sf!(),  _
lpg!())
DECLARE SUB rezcalc (r1!(), r2!(), r3!(), rez1!(), rez2!(), rez3!())
DECLARE SUB elneutrons ()
DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB getvar (a$, var$, var!)
DECLARE SUB inner (a!, b!(), c!())
DECLARE SUB vec (a!(), b!(), c!())
DECLARE SUB splitstring (ala$, d1#(), col1%, d2$(), col2%)
DIM sp(120)
REM this programm was checked and compared with powls programm for phase AF1
REM on 27.7.1995- agreement was found except multiplicity of 101 reflection
REM which seems to show wrong intensities (correspondence
REM with kockelmann confirmed this error in powls program)
REM subsequently enhanced to read output of spins.c (mcphase package)
IF COMMAND$ = "" OR COMMAND$ = "-h" OR COMMAND$ = "-help" THEN
PRINT "program eln - use as: eln filename"
PRINT "             calculates neutron diffraction pattern"
PRINT "             for format of input file filename see mcphase manual"
END
END IF



programpath$ = ".\"
'programpath$ = "c:\martin\gddip_pr\gdsachot\gd\gdag2\"


'scattering lengths [10^-12cm=10^-14m]
2 DATA 51
'       ^----number of elements
DATA "Ag",0.5922
DATA "Am",0.83
DATA "As",0.658
DATA "Au",0.763
DATA "B",0.530
DATA "Ba",0.525
DATA "Bi",0.8533
DATA "Br",0.679
DATA "C",0.66484
DATA "Cd",0.5
DATA "Ce",0.484
'sl Ce nach marshall lovesey 0.816 barn(i.e. 10^-12cm)
DATA "Cm",0.95
DATA "Cs",0.542
DATA "Cu",0.7718
DATA "Dy",1.69
DATA "Er",0.803
DATA "Eu",0.53
 ' sl Eu energieabh. !!
DATA "Fe",0.9452
DATA "Ge",0.8193
DATA "Gd",1.1
 'sl gd nach marshall lovesey 0.95 barn [energieabh!!!]
 ' Gd scattering length: lt chattopadhyay et al Sol Stat. Comm. 100 (1996) 117-122
 '                        ist slGd=11.3 fm bei lambda=0.5 A
 '                            slGd=6.5 fm bei lambda=1.8A
DATA "Hg",1.266
DATA "Ho",0.808
DATA "I",0.528
DATA "In",0.4065
 'slIn=0.4065-0.00539i   -- is complex
DATA "Kr",0.785
DATA "Mo",0.695
DATA "Mn",-0.375
DATA "Nb",0.7054
DATA "Nd",0.769
DATA "Ni",1.03
DATA "Pa",0.91
DATA "Pb",0.9401
DATA "Pd",0.591
DATA "Pr",0.458
DATA "Ra",1.0
DATA "Rb",0.708
DATA "Rh",0.593
DATA "Ru",0.721
DATA "Sb",0.564
DATA "Se",0.797
DATA "Si",0.41491
DATA "Sn",0.6228
DATA "Sr",0.702
DATA "Tb",0.738
'sl tb nach marshall lovesey 0.738 barn
DATA "Tc",0.68
DATA "Te",0.580
DATA "Th",0.984
DATA "U",0.842
DATA "Xe",0.489
DATA "Y",0.775
DATA "Zr",0.716

'magnetic formfactors (international tables)
3 DATA 8
'       ^----number of elements
'DATA "symbol",gJ(Landeefactor),
'     ff(i,1...7)= <j0(kr)>-terms A,a,B,b,C,c,D
'     ff(i,8..14)= <j2(kr)>-terms A,a,B,b,C,c,D

DATA "Ce3+", 0.85714
'gJ=6/7
DATA 0.2953, 17.685, 0.2923, 6.733, 0.4313, 5.383, -0.0194
DATA 0.9809, 18.063, 1.8413, 7.769, 0.9905, 2.845, 0.012

'magnetische formfaktoren nach international tables (eigentl nach Mcewen)
'rem mag formfaktor fÅr Ce2+ ion
'     ff(i,1...7)= <j0(kr)>-terms A,a,B,b,C,c,D
'ffCe(1) = .2953: ffCe(2) = 17.685: ffCe(3) = .2923: ffCe(4) = 6.733
'ffCe(5) = .4313: ffCe(6) = 5.383: ffCe(7) = -.0194
'     ff(i,8..14)= <j2(kr)>-terms A,a,B,b,C,c,D
'ffCe(8) = .9809: ffCe(9) = 18.063: ffCe(10) = 1.8413: ffCe(11) = 7.769
'ffCe(12) = .9905: ffCe(13) = 2.845: ffCe(14) = .012

DATA "Pr3+", 0.8
'gJ=0.8

DATA 0.0504, 24.9989, 0.2572, 12.0377, 0.7142, 5.0039, -0.0219
DATA 0.8734, 18.9876, 1.5594, 6.0872,  0.8142, 2.4150, 0.0111

DATA "Dy3+",1.33333
DATA 0.1157, 15.073, 0.327, 6.799, 0.5821, 3.02, -0.0249
DATA 0.2523, 18.517, 1.0914, 6.736, 0.9345, 2.208, 0.025
'magnetische formfaktoren nach international tables
'     ff(i,1...7)= <j0(kr)>-terms A,a,B,b,C,c,D
'ffDy(1) = .1157: ffDy(2) = 15.073: ffDy(3) = .327: ffDy(4) = 6.799
'ffDy(5) = .5821: ffDy(6) = 3.02: ffDy(7) = -.0249
'     ff(i,8..14)= <j2(kr)>-terms A,a,B,b,C,c,D
'ffDy(8) = .2523: ffDy(9) = 18.517: ffDy(10) = 1.0914: ffDy(11) = 6.736
'ffDy(12) = .9345: ffDy(13) = 2.208: ffDy(14) = .025


DATA "Gd3+",2
DATA 0.0186, 25.387, 0.2895, 11.142, 0.7135, 3.752, -.0217
DATA 0.3347, 18.476, 1.2465, 6.877, 0.9537, 2.318, 0.0217
'magnetische formfaktoren nach international tables
'     ff(i,1...7)= <j0(kr)>-terms A,a,B,b,C,c,D
'ffGd(1) = .0186: ffGd(2) = 25.387: ffGd(3) = .2895: ffGd(4) = 11.142
'ffGd(5) = .7135: ffGd(6) = 3.752: ffGd(7) = -.0217
'     ff(i,8..14)= <j2(kr)>-terms A,a,B,b,C,c,D
'ffGd(8) = .3347: ffGd(9) = 18.476: ffGd(10) = 1.2465: ffGd(11) = 6.877
'ffGd(12) = .9537: ffGd(13) = 2.318: ffGd(14) = .0217

DATA "Mn",2
DATA 0.2438, 24.9629, 0.1472, 15.6728, 0.6189, 6.5403, -.0105
DATA 2.6681, 16.0601, 1.7561, 5.6396, 0.3675, 2.0488, 0.0017

DATA "Nd3+",0.72727272
DATA 0.054, 25.029, 0.3101, 12.102, 0.6575, 4.722, -0.0216
DATA 0.6751, 18.342, 1.6272, 7.26, 0.9644, 2.602, 0.015
'magnetische formfaktoren nach international tables
'     ff(i,1...7)= <j0(kr)>-terms A,a,B,b,C,c,D
'ffNd(1) = .054: ffNd(2) = 25.029: ffNd(3) = .3101: ffNd(4) = 12.102
'ffNd(5) = .6575: ffNd(6) = 4.722: ffNd(7) = -.0216
'     ff(i,8..14)= <j2(kr)>-terms A,a,B,b,C,c,D
'ffNd(8) = .6751: ffNd(9) = 18.342: ffNd(10) = 1.6272: ffNd(11) = 7.26
'ffNd(12) = .9644: ffNd(13) = 2.602: ffNd(14) = .015


DATA "Tb3+",1.5
DATA 0.0177, 25.51,0.2921, 10.577, 0.7133, 3.512, -0.0231
DATA 0.2892, 18.497, 1.1678, 6.797, 0.9437, 2.257, 0.0232
'magnetische formfaktoren nach international tables
'     ff(i,1...7)= <j0(kr)>-terms A,a,B,b,C,c,D
'fftb(1) = .0177: fftb(2) = 25.51: fftb(3) = .2921: fftb(4) = 10.577
'fftb(5) = .7133: fftb(6) = 3.512: fftb(7) = -.0231
'     ff(i,8..14)= <j2(kr)>-terms A,a,B,b,C,c,D
'fftb(8) = .2892: fftb(9) = 18.497: fftb(10) = 1.1678: fftb(11) = 6.797
'fftb(12) = .9437: fftb(13) = 2.257: fftb(14) = .0232
DATA "Tm3+",1.16666667
DATA 0.0581, 15.0922, 0.2787,7.8015,0.6854,2.7931,-0.0224
DATA 0.176, 18.5417, 0.9105, 6.5787,0.897,2.0622,0.0294


CALL elneutrons
END

SUB elneutrons

DIM h(4120), k(4120), l(4120), D(4120), imag(4120), ikern(4120)
DIM sf(4120), lpg(4120), theta(4120), ffGd(14)
DIM r1(3), r2(3), r3(3)
DIM rez1(3), rez2(3), rez3(3)
SHARED programpath$
REM beschreibung der elemtarzelle input from file eln.cel
pi = 3.141592654#

'---------------------------------------------------------------
'read scattering lengths
RESTORE 2
READ slnumber%
DIM sl(slnumber%), sl$(slnumber%)
FOR i% = 1 TO slnumber%: READ sl$(i%), sl(i%): NEXT i%

'read form factors
RESTORE 3
READ ffnumber%
DIM ff(ffnumber%, 14), ff$(ffnumber%)
FOR i% = 1 TO ffnumber%: READ ff$(i%), ff(i%, 0)
   FOR ii% = 1 TO 14: READ ff(i%, ii%): NEXT ii%
NEXT i%
'-------------------------------------------------------------

'READ FROM FILE -----------------------------
'enn = 327.2: REM energy in meV
' lambda = .58   '1 / SQR(enn / 81.81): REM wavelength [A]
''thetamax = 12
' a = 6.62: b = 6.62: c = 6.62:   'gitterkonstante [A] atRT (lindbaum rt)
'nna = sp(2): nnb = sp(4): nnc = sp(6)
'r1(1) = a * nna 'primitive unit cell vector (read from file)
'r2(2) = b * nnb
'r3(3) = c * nnc
'---------------------------------------------
IF LTRIM$(COMMAND$) = "" THEN file$ = programpath$ + "eln.cel" ELSE file$ = COMMAND$


OPEN "i", 1, file$
PRINT "ELN: reading file "; file$
'check for commas
WHILE EOF(1) = 0: a$ = INPUT$(1, #1): IF a$ = "," THEN PRINT "ERROR ELN: input file contains commas ','": END
 WEND
SEEK 1, 1

'input section 1 *******************************************************
a$ = "#"
WHILE (LEFT$(LTRIM$(a$), 1) = "#" AND INSTR(a$, "%SECTION 2%") = 0)
aa = SEEK(1): INPUT #1, a$
 CALL getvar(a$, "lambda", lambda)
 CALL getvar(a$, "thetamax", thetamax)
 CALL getvar(a$, "nat", nat)
 CALL getvar(a$, "ovalltemp", ovalltemp)
 CALL getvar(a$, "lorentz", lorenz)
WEND
SEEK 1, aa

IF lorenz = 0 THEN lorenz = 1
IF lambda = 0 THEN PRINT "ERROR ELN: no wavelength lambda given or line does not start with # in section 1": END
IF thetamax = 0 THEN PRINT "ERROR ELN: no thetamax given or line does not start with # in section 1": END
PRINT "             section 1 - lambda="; lambda; " A thetamax="; thetamax; "deg"
PRINT "                         ovalltemp="; ovalltemp; " A^2 lorentz-type="; lorenz

'input section 2 *********************************************************
DIM d1#(40), d2$(4)
a$ = "#"
WHILE (LEFT$(LTRIM$(a$), 1) = "#" AND a = 0 AND b = 0 AND c = 0)
aa = SEEK(1): INPUT #1, a$
 CALL getvar(a$, "nat", nat)
 CALL getvar(a$, " a", a)
 CALL getvar(a$, " b", b)
 CALL getvar(a$, " c", c)
WEND
SEEK 1, aa
PRINT "             section 2 - nat="; nat

DIM x1(nat), y1(nat), z1(nat), nam$(nat)
IF nat <> 0 THEN
FOR i% = 1 TO nat
31 CALL inputline(1, d1#(), col1%, d2$(), col2%)
         IF col1% = -1 THEN
          INPUT #1, a$: IF INSTR(a$, "%SECTION 3%") > 0 THEN PRINT "ERROR ELN: Section 3 started before all nat="; nat; "atoms of crystallographic unit cell were listed !": END
          GOTO 31
         END IF
 IF col1% < 6 THEN PRINT "ERROR ELN: Section 2 - Nonmagnetic Atoms: too few positional parameters for atom "; i%; "!": END
 x1(i%) = d1#(4)
 y1(i%) = d1#(5)
 z1(i%) = d1#(6)
 nam$(i%) = d2$(1)
 PRINT "                         "; nam$(i%); " at "; x1(i%); " * r1 +"; y1(i%); " * r2 +"; z1(i%); " * r3"
NEXT i%
END IF
'input section 3 *********************************************************

a$ = "#"
WHILE (LEFT$(LTRIM$(a$), 1) = "#" AND nna * nnb * nnc = 0)
 INPUT #1, a$
 CALL getvar(a$, " a", a)
 CALL getvar(a$, " b", b)
 CALL getvar(a$, " c", c)
 CALL getvar(a$, "r1x", r1(1))
 CALL getvar(a$, "r1y", r1(2))
 CALL getvar(a$, "r1z", r1(3))
 CALL getvar(a$, "r2x", r2(1))
 CALL getvar(a$, "r2y", r2(2))
 CALL getvar(a$, "r2z", r2(3))
 CALL getvar(a$, "r3x", r3(1))
 CALL getvar(a$, "r3y", r3(2))
 CALL getvar(a$, "r3z", r3(3))
 CALL getvar(a$, "nr1", nna)
 CALL getvar(a$, "nr2", nnb)
 CALL getvar(a$, "nr3", nnc)
 CALL getvar(a$, "nat", natmagnetic)
WEND
IF a = 0 THEN PRINT "ERROR ELN: no lattice constant a given in section 3 or line does not start with #: "; a$: END
IF b = 0 THEN PRINT "ERROR ELN: no lattice constant b given in section 3 or line does not start with #: "; a$: END
IF c = 0 THEN PRINT "ERROR ELN: no lattice constant c given in section 3 or line does not start with #: "; a$: END
PRINT "             section 3 - a="; a; " A  b="; b; " A c="; c; "A"
unitcell$ = "# a=" + STR$(a) + " A  b=" + STR$(b) + " A c=" + STR$(c) + " A  alpha=90  beta=90 gamma=90"
PRINT USING "                            / ##.###a \     / ##.###a\     / ##.###a \"; r1(1); r2(1); r3(1)
PRINT USING "                         r1=| ##.###b |  r2=| ##.###b|  r3=| ##.###b |"; r1(2); r2(2); r3(2)
PRINT USING "                            \ ##.###c /     \ ##.###c/     \ ##.###c /"; r1(3); r2(3); r3(3)
r1(1) = a * r1(1) * nna
r2(1) = a * r2(1) * nnb
r3(1) = a * r3(1) * nnc
r1(2) = b * r1(2) * nna
r2(2) = b * r2(2) * nnb
r3(2) = b * r3(2) * nnc
r1(3) = c * r1(3) * nna
r2(3) = c * r2(3) * nnb
r3(3) = c * r3(3) * nnc

'input section 4 *********************************************************

IF nna = 0 THEN PRINT "ERROR ELN: nr1 not given or line does not start with # in section 4": END
IF nnb = 0 THEN PRINT "ERROR ELN: nr2 not given or line does not start with # in section 4": END
IF nnc = 0 THEN PRINT "ERROR ELN: nr3 not given or line does not start with # in section 4": END
PRINT "             section 4 - nr1="; nna; "  nr2="; nnb; "  nr3="; nnc
PRINT "                         nat="; natmagnetic; " magnetic atoms"

'----------------------------------------------------------------

n = nna * nnb * nnc * nat + natmagnetic: 'atoms in der magnetic unit cell
DIM x(n), y(n), z(n), mx(n), my(n), mz(n), slcode%(n), ffcode%(n), occ(n)
PRINT "                         reading magnetic atoms and moments ..."

FOR i% = 1 TO natmagnetic
41 CALL inputline(1, d1#(), col1%, d2$(), col2%)
         IF col1% = -1 THEN INPUT #1, a$: GOTO 41
 IF col1% < 6 THEN PRINT "ERROR ELN: Section 4 - Nonmagnetic Atoms: too few positional parameters for atom "; i%; "!": END
 x(i%) = d1#(4) / nna
 y(i%) = d1#(5) / nnb
 z(i%) = d1#(6) / nnc
 
 nam$ = d2$(1)
' PRINT "                         "; nam$; " at "; x1(i%); " * r1 +"; y1(i%); " * r2 +"; z1(i%); " * r3 ... m=gJ*(";
' PRINT d1#(7); ","; d1#(8); ","; d1#(9); ") m_Bohr"
 'determine occupancy
 occ(i%) = VAL(nam$): IF occ(i%) = 0 THEN occ(i%) = 1

 'delete the leading number from nam$
 ii% = 1: WHILE (ASC(MID$(nam$, ii%, 1)) < 58)
          ii% = ii% + 1
          WEND
 nam$ = MID$(nam$, ii%)


 'determine slcode%
 FOR ii% = 1 TO slnumber%: IF LEFT$(LTRIM$(nam$), 2) = sl$(ii%) THEN slcode%(i%) = ii%
 NEXT ii%:
 IF slcode%(i%) = 0 THEN
  PRINT "ERROR ELN: scattering length of element "; nam$; " not found in internal tables"
  PRINT "Implemented are Scattering lengths of "; : FOR ii% = 1 TO slnumber%: PRINT sl$(ii%); " "; : NEXT
  END
 END IF
 'determine ffcode%
 FOR ii% = 1 TO ffnumber%: IF nam$ = ff$(ii%) THEN ffcode%(i%) = ii%
 NEXT ii%:
  IF ffcode%(i%) = 0 THEN
   PRINT "error form factor of element "; nam$; " not found in internal tables"
   PRINT "Implemented are form factors of "; : FOR ii% = 1 TO ffnumber%: PRINT ff$(ii%); " "; : NEXT
   END
  END IF
 mx(i%) = d1#(7) * ff(ffcode%(i%), 0)
 my(i%) = d1#(8) * ff(ffcode%(i%), 0)
 mz(i%) = d1#(9) * ff(ffcode%(i%), 0)
NEXT i%

CLOSE 1
PRINT "     ... calculating ... "
'now insert also nonmagnetic elements into the unit cell
ncryst = natmagnetic
FOR na = 1 TO nna
FOR nb = 1 TO nnb
FOR nc = 1 TO nnc
IF nat <> 0 THEN
 FOR i% = 1 TO nat
 ncryst = ncryst + 1
 'determine occupancy
 occ(ncryst) = VAL(nam$(i%)): IF occ(ncryst) = 0 THEN occ(ncryst) = 1

 'delete the leading number from nam$
 ii% = 1: WHILE (ASC(MID$(nam$(i%), ii%, 1)) < 58)
          ii% = ii% + 1
          WEND
 nam$ = MID$(nam$(i%), ii%)

 'determine slcode%
 FOR ii% = 1 TO slnumber%: IF nam$ = sl$(ii%) THEN slcode%(ncryst) = ii%
 NEXT ii%: IF slcode%(ncryst) = 0 THEN PRINT "error scattering length of element "; nam$; "not found in internal tables": END
 x(ncryst) = (na + x1(i%) - 1) / nna
 y(ncryst) = (nb + y1(i%) - 1) / nnb
 z(ncryst) = (nc + z1(i%) - 1) / nnc
 NEXT i%
END IF
NEXT nc: NEXT nb: NEXT na

'________________________________________________________________

CALL neutint(lambda, thetamax, ovalltemp, lorenz, r1(), r2(), r3(), n, x(), y(), z(), mx(), my(), mz(), sl(), slcode%(), ff(), ffcode%(), occ(), m, h(), k(), l(), D(), theta(), imag(), ikern(), sf(), lpg())
filename$ = programpath$ + "eln.out"

CALL rezcalc(r1(), r2(), r3(), rez1(), rez2(), rez3())

REM transformieren der millerindizes auf kristallographische einheitszelle
FOR i = 1 TO m:
hh = h(i): kk = k(i): ll = l(i)
 h(i) = (hh * rez1(1) + kk * rez2(1) + ll * rez3(1)) / 2 / pi * a
 k(i) = (hh * rez1(2) + kk * rez2(2) + ll * rez3(2)) / 2 / pi * b
 l(i) = (hh * rez1(3) + kk * rez2(3) + ll * rez3(3)) / 2 / pi * c
NEXT



CALL printeln(filename$, unitcell$, lambda, ovalltemp, lorenz, r1(), r2(), r3(), n, x(), y(), z(), mx(), my(), mz(), sl(), slcode%(), ff(), ffcode%(), occ(), m, h(), k(), l(), ikern(), imag(), D(), theta(), sf(), lpg())


END SUB

SUB getvar (a$, var$, var)
' gets variable var$ out of string a$
po% = INSTR(a$, var$)
eq% = INSTR(po% + 1, a$, "=")
IF eq% - po% > LEN(var$) THEN
   betw$ = LTRIM$(MID$(a$, po% + LEN(var$), eq% - po% - LEN(var$)))
   IF betw$ <> "" THEN eq% = 0
END IF


IF po% = 0 OR eq% = 0 THEN
' var = 0
ELSE
 var = VAL(MID$(a$, eq% + 1))
END IF

END SUB

SUB inner (a, b(), c())
'inner product a=(b.c)
a = 0
FOR i% = 1 TO 3
a = a + b(i%) * c(i%)
NEXT i%

END SUB

SUB inputline (n, d1#(), col1%, d2$(), col2%)

'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
'comments in [] are stored in d2$
'if a comment is started somewhere in this line by "[" and not finished
'then col1% is set -1 and the filepointer is set to the beginning of the
'line (with seek)
'# means comment line and col1% is set to -1

a$ = INKEY$: IF a$ <> "" THEN IF ASC(a$) = 27 THEN END


aa = SEEK(n)
LINE INPUT #n, ala$
IF LEFT$(LTRIM$(ala$), 1) = "#" THEN
  col1% = -1: SEEK n, aa 'a comment line ... no data read
ELSE
WHILE INSTR(ala$, CHR$(9)) > 0  'abandon tabs
 i% = INSTR(ala$, CHR$(9))
 ala$ = LEFT$(ala$, i% - 1) + " " + MID$(ala$, i% + 1)
WEND
'treat {} comments in input line
klauf% = INSTR(ala$, "{")
klzu% = INSTR(ala$, "}")
col2% = 0
WHILE klauf% < klzu% AND klauf% > 0   'take out closed bracketed expressions
 col2% = col2% + 1
 d2$(col2%) = MID$(ala$, klauf% + 1, klzu% - klauf% - 1)
 ala$ = LEFT$(ala$, klauf% - 1) + " " + MID$(ala$, klzu% + 1)
 klauf% = INSTR(ala$, "{")
 klzu% = INSTR(ala$, "}")
WEND
IF klauf% > 0 THEN
 col1% = -1: SEEK n, aa 'a comment bracket is not closed ... no data read
ELSE
 'treat [] comments in input line
 klauf% = INSTR(ala$, "[")
 klzu% = INSTR(ala$, "]")
 WHILE klauf% < klzu% AND klauf% > 0   'take out closed bracketed expressions
  col2% = col2% + 1
  d2$(col2%) = MID$(ala$, klauf% + 1, klzu% - klauf% - 1)
  ala$ = LEFT$(ala$, klauf% - 1) + " " + MID$(ala$, klzu% + 1)
  klauf% = INSTR(ala$, "[")
  klzu% = INSTR(ala$, "]")
 WEND
 

 IF klauf% > 0 THEN
  col1% = -1: SEEK n, aa 'a comment bracket is not closed ... no data read
 ELSE
  col1% = 0
  WHILE LEN(ala$) > 0
     col1% = col1% + 1
     ala$ = LTRIM$(ala$) + " "
     d1#(col1%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
     ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " ")))
  WEND
 END IF
END IF
END IF
END SUB

SUB neutint (lambda, thetamax, ovalltemp, lorenz, r1(), r2(), r3(), n, x(), y(), z(), mx(), my(), mz(), sl(), slcode%(), ff(), ffcode%(), occ(), m, h(), k(), l(), D(), theta(), imag(), ikern(), sf(), lpg())
REM****************************************************************************
REM this routine calculates the intensity of elastic neutrons
REM for a given magnetic unit cell (crystal axis orthogonal)
REM the magnetic scattering is treated in the dipole approximation
REM input:
 REM lambda                             wavelength[A]
 REM ovalltemp                          overall temperature factor [A]
 REM r1(1..3),r2(),r3()                 vectors of primitive unit cell[A]
 REM n                                  number of atoms per unit cell
 REM x(1...n),y(1...n),z(1...n)         atomic positional parameters
                                        '(with respect to primitive lattice)
 REM mx(1...n),my(1...n),mz(1...n)      atomic magnetic moment [mb]
                                       ' (with respect to ortholatice xyz)
 REM sl(1...c)                          nuclear scattering length[10^-12cm]
 REM slcode%(1...n)                     scattering length codefor atom 1...n
 REM ff(1...c,0)                          Landefactors gJ
 REM ffcode%(1..n)                      formfactor code for atom 1...n
 REM ff(1...c,1...14)                          4f magnetic form factor
 REM                                according to  Volume C of the new edition
'                                  of the international Tables
'                                  for Crystallography published by the
'                                   International Union of Crystallography
'     ff(i,1...7)= <j0(kr)>-terms A,a,B,b,C,c,D
'     ff(i,8..14)= <j2(kr)>-terms A,a,B,b,C,c,D
'     <jl(kr)> is defined as = integral[0,inf] U^2(r) jl(kr) 4 pi r^2 dr
' where U(r) is the Radial wave function for the unpaired electrons in the atom

REM output
 REM m                                  number of calculated reflections
 REM d(1...m)                           d spacing
 REM theta(1...m)                       scattering angle theta
 REM imag(1...m)                        magnetic intensity
 REM ikern(1...m)                       nuclear intensity
 REM sf(1...m)                          nuclear structurfactor |sf|
 REM lpg(1...m)                         lorentzfactor
REM****experimental parameters*************************************************
scale = 1 / n / n: REM scalingfactor of intensities
REM***************************************************************************

pi = 3.141592654#: D(0) = 10000
m = 0: 'reset m

'calculate reciprocal lattice vectors from r1,r2,r3
DIM rez1(3), rez2(3), rez3(3), Qvec(3)
CALL rezcalc(r1(), r2(), r3(), rez1(), rez2(), rez3())

qmax = 4 * pi * SIN(thetamax / 180 * pi) / lambda
CALL inner(rr, r1(), r1()): hmax = INT(qmax / 2 / pi * SQR(rr)) + 1
CALL inner(rr, r2(), r2()): kmax = INT(qmax / 2 / pi * SQR(rr)) + 1
CALL inner(rr, r3(), r3()): lmax = INT(qmax / 2 / pi * SQR(rr)) + 1

FOR ahi = 0 TO hmax
FOR aki = 0 TO kmax
FOR ali = 0 TO lmax

FOR sh = -1 TO 1 STEP 2
FOR sk = -1 TO 1 STEP 2
FOR sl = -1 TO 1 STEP 2
IF ahi = 0 THEN sh = 1: hi = 0 ELSE hi = sh * ahi
IF aki = 0 THEN sk = 1: ki = 0 ELSE ki = sk * aki
IF ali = 0 THEN sl = 1: li = 0 ELSE li = sl * ali


IF hi = 0 AND li = 0 AND ki = 0 THEN ltrue% = 1: ktrue% = 1: htrue% = 1: GOTO 30

'calculate d spacing and intensity for h,k,l triple (d,imag,ikern)************
'd = 1 / SQR(hi * hi / a / a + ki * ki / b / b + li * li / c / c): 'dspacing
FOR i% = 1 TO 3
 Qvec(i%) = hi * rez1(i%) + ki * rez2(i%) + li * rez3(i%)
NEXT i%
Q = SQR(Qvec(1) * Qvec(1) + Qvec(2) * Qvec(2) + Qvec(3) * Qvec(3)): 'dspacing
D = 2 * pi / Q
s = 1 / 2 / D: sintheta = lambda * s
IF sintheta > SIN(thetamax / 180 * pi) GOTO 30
 theta = 180 / pi * ATN(sintheta / SQR(1 - sintheta * sintheta))
'nuclear(|nsfr+i nsfc|^2) and magnetic structure factor(msf) calculation
 nsfr = 0: nsfc = 0: msfrx = 0: msfcx = 0
 msfry = 0: msfcy = 0: msfrz = 0: msfcz = 0
 FOR i% = 1 TO n
  qr = hi * x(i%) + ki * y(i%) + li * z(i%)
  'nuclear structure factor nsfr,nsfc
  nsfr = nsfr + COS(-2 * pi * qr) * occ(i%) * sl(slcode%(i%))
  nsfc = nsfc + SIN(-2 * pi * qr) * occ(i%) * sl(slcode%(i%))

 'magnetic structure factors
  j0 = ff(ffcode%(i%), 1) * EXP(-ff(ffcode%(i%), 2) * s * s) + ff(ffcode%(i%), 3) * EXP(-ff(ffcode%(i%), 4) * s * s)
  j0 = j0 + ff(ffcode%(i%), 5) * EXP(-ff(ffcode%(i%), 6) * s * s) + ff(ffcode%(i%), 7)
  j2 = ff(ffcode%(i%), 8) * s * s * EXP(-ff(ffcode%(i%), 9) * s * s) + ff(ffcode%(i%), 10) * s * s * EXP(-ff(ffcode%(i%), 11) * s * s)
  j2 = j2 + ff(ffcode%(i%), 12) * s * s * EXP(-ff(ffcode%(i%), 13) * s * s) + s * s * ff(ffcode%(i%), 14)

  IF ff(ffcode%(i%), 0) <> 0 THEN
  FQ = occ(i%) * (j0 + j2 * (2 / ff(ffcode%(i%), 0) - 1)): REM occupance * formfactor F(Q)
  msfrx = msfrx + COS(-2 * pi * qr) * FQ / 2 * mx(i%)
  msfcx = msfcx + SIN(-2 * pi * qr) * FQ / 2 * mx(i%)
  msfry = msfry + COS(-2 * pi * qr) * FQ / 2 * my(i%)
  msfcy = msfcy + SIN(-2 * pi * qr) * FQ / 2 * my(i%)
  msfrz = msfrz + COS(-2 * pi * qr) * FQ / 2 * mz(i%)
  msfcz = msfcz + SIN(-2 * pi * qr) * FQ / 2 * mz(i%)
END IF
NEXT i%

'magnetic structure factors + polarisation factor===>msf
msf = msfrx * msfrx + msfcx * msfcx + msfry * msfry + msfcy * msfcy + msfrz * msfrz + msfcz * msfcz
'msf = msf - 2 * hi / a * d * ki / b * d * (msfrx * msfry + msfcx * msfcy)
'msf = msf - 2 * hi / a * d * li / c * d * (msfrx * msfrz + msfcx * msfcz)
'msf = msf - 2 * li / c * d * ki / b * d * (msfry * msfrz + msfcy * msfcz)
msf = msf - 2 * Qvec(1) * Qvec(2) / Q / Q * (msfrx * msfry + msfcx * msfcy)
msf = msf - 2 * Qvec(1) * Qvec(3) / Q / Q * (msfrx * msfrz + msfcx * msfcz)
msf = msf - 2 * Qvec(2) * Qvec(3) / Q / Q * (msfry * msfrz + msfcy * msfcz)

'msf = msf - hi / a * d * hi / a * d * (msfrx * msfrx + msfcx * msfcx)
'msf = msf - ki / b * d * ki / b * d * (msfry * msfry + msfcy * msfcy)
'msf = msf - li / c * d * li / c * d * (msfrz * msfrz + msfcz * msfcz)
msf = msf - Qvec(1) * Qvec(1) / Q / Q * (msfrx * msfrx + msfcx * msfcx)
msf = msf - Qvec(2) * Qvec(2) / Q / Q * (msfry * msfry + msfcy * msfcy)
msf = msf - Qvec(3) * Qvec(3) / Q / Q * (msfrz * msfrz + msfcz * msfcz)


'lorentzfactor*************************************************************
sin2theta = 2 * sintheta * SQR(1 - sintheta * sintheta)
IF lorenz = 1 THEN lorentzf = 1 / sin2theta / sin2theta  'powder flat sample
IF lorenz = 2 THEN lorentzf = 1 / sin2theta / sintheta   'powder cyl. sample
IF lorenz = 3 THEN lorentzf = 1 / sin2theta              'single crystal
IF lorenz = 4 THEN lorentzf = D * D * D     'TOF powder cyl sample... log scaled d-pattern
IF lorenz = 5 THEN lorentzf = D * D * D * D 'TOF powder cyl sample... d-pattern


'overall temperature factor*************************************************
ovallt = EXP(-2 * ovalltemp * (sintheta * sintheta / lambda / lambda))
'***************************************************************************
    'A)nuclear intenisty

sf = SQR(nsfr * nsfr + nsfc * nsfc)
ikern = sf * sf * lorentzf * scale * ovallt

    'B)magnetic intensity
imag = msf * 3.65 / 4 / pi * lorentzf * scale * ovallt

'sort according to descending d spacing
IF imag + ikern > .01 THEN

m = m + 1: IF m > 4120 THEN PRINT "ERROR eln: out of memory - too many reflections - chose smaller thetamax": END
 msort = m
20 IF D(msort - 1) <= D THEN
    D(msort) = D(msort - 1)
    theta(msort) = theta(msort - 1)
    h(msort) = h(msort - 1)
    k(msort) = k(msort - 1)
    l(msort) = l(msort - 1)
    imag(msort) = imag(msort - 1)
    ikern(msort) = ikern(msort - 1)
    sf(msort) = sf(msort - 1)
    lpg(msort) = lpg(msort - 1)
   msort = msort - 1: GOTO 20
END IF
h(msort) = hi: k(msort) = ki: l(msort) = li
D(msort) = D: theta(msort) = theta
imag(msort) = imag: ikern(msort) = ikern
sf(msort) = sf: lpg(msort) = lorentzf
END IF
30 NEXT sl: NEXT sk: NEXT sh

 
 NEXT ali
 NEXT aki
32  PRINT USING "### "; ahi;
 NEXT ahi

33 END SUB

SUB printeln (filename$, unitcell$, lambda, ovalltemp, lorenz, r1(), r2(), r3(), n, x(), y(), z(), mx(), my(), mz(), sl(), slcode%(), ff(), ffcode%(), occ(), m, h(), k(), l(), ikern(), imag(), D(), theta(), sf(), lpg())
REM ausgabe auf file filename$
OPEN "o", 1, filename$
PRINT #1, "#{"; filename$; " "; TIME$; DATE$; "   unit cell:"
PRINT #1, unitcell$
PRINT #1, "#r1x="; r1(1); " A r2x="; r2(1); " A r3x="; r3(1)
PRINT #1, "#r1y="; r1(2); " A r2y="; r2(2); " A r3y="; r3(2)
PRINT #1, "#r1z="; r1(3); " A r2z="; r2(3); " A r3z="; r3(3)
PRINT #1, "#wavelength="; lambda; "A   number of atoms:"; n
PRINT #1, "#overall temperature factor: exp(-2*"; ovalltemp; "A*(sin(theta)/lambda)^2)"
IF lorenz = 1 THEN l$ = "1 / sin^2(2theta)   neutron powder flat sample"
IF lorenz = 2 THEN l$ = "1 / sin(2theta) / sin(theta)    neutron powder cyl. sample"
IF lorenz = 3 THEN l$ = "1 / sin(2theta)     neutron single crystal"
IF lorenz = 4 THEN l$ = "d^3  neutron TOF powder cyl sample... log scaled d-pattern"
IF lorenz = 5 THEN l$ = "d^4  neutron TOF powder cyl sample... d-pattern"
PRINT #1, "#Lorentz Factor: "; l$
PRINT #1, "#"
PRINT #1, "#    x(r1)    y(r2)  z(r3)  mx(mb) my(mb) mz(mb) occupancy  sl     gJ   ff"
FOR i = 1 TO n
PRINT #1, "# ";
PRINT #1, USING "###.###"; x(i);
PRINT #1, USING " ###.###"; y(i);
PRINT #1, USING " ###.###"; z(i);
PRINT #1, USING " ###.###"; mx(i);
PRINT #1, USING " ###.###"; my(i);
PRINT #1, USING " ###.###"; mz(i);
PRINT #1, USING " ###.###"; occ(i);
PRINT #1, USING " ###.###"; sl(slcode%(i));
FOR j = 0 TO 14
PRINT #1, USING " ###.###"; ff(ffcode%(i), j);
NEXT j
PRINT #1, "#"
NEXT
PRINT #1, "# h       k       l      d[A]    |Q|[A^-1]    2theta    Inuc(2th) Imag(2th) Itot(2t) |sf|   lpg }"
h(0) = 0: k(0) = 0: l(0) = 0: D(0) = 100: theta(0) = 0: ikern(0) = 0: imag(0) = 0: sf(0) = 0: lpg(0) = 0
FOR i = 0 TO m
PRINT #1, USING "###.###"; h(i);
PRINT #1, USING " ###.###"; k(i);
PRINT #1, USING " ###.###"; l(i);
PRINT #1, USING " ###.#####"; D(i);
PRINT #1, USING " ###.#####"; 2 * 3.14152 / D(i);
PRINT #1, USING " ####.###"; 2 * theta(i);
PRINT #1, USING " #####.###"; ikern(i);
PRINT #1, USING " #####.###"; imag(i);
PRINT #1, USING " #####.###"; imag(i) + ikern(i);
PRINT #1, USING " #####.###"; sf(i);
PRINT #1, USING " #####.###"; lpg(i)
NEXT
CLOSE 1
END SUB

SUB rezcalc (r1(), r2(), r3(), rez1(), rez2(), rez3())
' calculate reciprocal lattice rezi from real lattice ri
DIM dd(3)
pi = 3.141592654#
CALL vec(rez1(), r2(), r3()): CALL vec(dd(), r2(), r3())
                              CALL inner(vol, r1(), dd())
FOR i% = 1 TO 3: rez1(i%) = rez1(i%) * 2 * pi / vol: NEXT i%

CALL vec(rez2(), r1(), r3()): CALL vec(dd(), r1(), r3())
                              CALL inner(vol, r2(), dd())
FOR i% = 1 TO 3: rez2(i%) = rez2(i%) * 2 * pi / vol: NEXT i%

CALL vec(rez3(), r1(), r2()): CALL vec(dd(), r1(), r2())
                              CALL inner(vol, r3(), dd())
FOR i% = 1 TO 3: rez3(i%) = rez3(i%) * 2 * pi / vol: NEXT i%

END SUB

SUB splitstring (ala$, d1#(), col1%, d2$(), col2%)
'split ala$ into numbers and strings

WHILE INSTR(ala$, CHR$(9)) > 0  'abandon tabs
 i% = INSTR(ala$, CHR$(9))
 ala$ = LEFT$(ala$, i% - 1) + " " + MID$(ala$, i% + 1)
WEND
'treat comments in input line
klauf% = INSTR(ala$, "{")
klzu% = INSTR(ala$, "}")
col2% = 0
WHILE klauf% < klzu% AND klauf% > 0   'take out closed bracketed expressions
 col2% = col2% + 1
 d2$(col2%) = MID$(ala$, klauf% + 1, klzu% - klauf% - 1)
 ala$ = LEFT$(ala$, klauf% - 1) + " " + MID$(ala$, klzu% + 1)
 klauf% = INSTR(ala$, "{")
 klzu% = INSTR(ala$, "}")
WEND


IF klauf% > 0 THEN
 col1% = -1: SEEK n, aa 'a comment bracket is not closed ... no data read
ELSE
 col1% = 0
 WHILE LEN(ala$) > 0
    col1% = col1% + 1
    ala$ = LTRIM$(ala$) + " "
    d1#(col1%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
    ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " ")))
 WEND
END IF


END SUB

SUB vec (a(), b(), c())
'vector product a=bxc
a(1) = b(2) * c(3) - b(3) * c(2)
a(2) = b(3) * c(1) - b(1) * c(3)
a(3) = b(1) * c(2) - b(2) * c(1)

END SUB

