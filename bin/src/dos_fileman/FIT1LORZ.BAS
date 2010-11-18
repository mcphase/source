DECLARE SUB calcsta (gr%, firstvalue!, par#(), sta#, ii%, jj%, xmin!, xmax!, ymin!, ymax!)
DECLARE SUB simanneal (firstvalue!, par#(), sta#, ii%, jj%, nmax!, speed!)
DECLARE SUB getpar (x!, x$, a$)
DECLARE SUB calcstartvalues (firstvalue!, par#(), parstp!(), parmax!(), parmin!(), ii%, jj%, xmin!, xmax!, ymin!, ymax!)
DECLARE SUB pointdraw3d (x!, xmin!, xmax!, y!, ymin!, ymax!, z!, zmin!, zmax!, size!)
DECLARE SUB minmax (firstvalue!, ii%, jj%, ymin!, ymax!, xmin!, xmax!, xmid!)
DECLARE SUB calcfunction (par#(), y#, x#)
DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB iteratefiles (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "FIT1LORZ FIT1LORZ FIT1LORZ FIT1LORZ FIT1LORZ FIT1LORZ FIT1LORZ FIT1LORZ FIT1LORZ FIT1LORZ"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program FIT1LORZ - use it like FIT1LORZ *.* 23"
DATA "    (23 means calc.1-lineLORentZian fit of 3rd columns with respect to second column) - "
DATA " ---> the result is written in file *.*,    /G..... view graphics when fitting  "
DATA " the LORentZianfit  is calculated by least squares method              "
DATA " the result is written to the last column of file *.*         "
DATA " format of file:            valid parameters are nmax=1e4 ..noofiterations"
DATA " { header: this is the                    speed=1  ..speed of Tvariation
DATA "   file header delimited                  /b ... fit also const bkgrd                               "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                      /c... do not wait for return-continue "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "*************************************************************************"
DIM text$(300), x#(30, 100), d1#(30), d2$(10), par#(30), av#(30)

'aaaaa analyse command string (get filename$ - ii% - jj% and interval) aaaaaa
a$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))
i% = INSTR(a$, " ")
filename$ = RTRIM$(LEFT$(a$, i% - 1))
a$ = LEFT$(RTRIM$(LTRIM$(MID$(a$, i%))), 2)

jj% = ASC(RIGHT$(a$, 1)) MOD 48: ii% = ASC(MID$(a$, LEN(a$) - 1, 1)) MOD 48  'get column numbers
 
CALL getpar(nmax, "NMAX", COMMAND$)
IF nmax = 0 THEN nmax = 10000
CALL getpar(speed, "SPEED", COMMAND$)
IF speed = 0 THEN speed = 1


999 CALL iteratefiles(filename$)

'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'oooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooooooooooooooooooooo
'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)
firstvalue = SEEK(1)
IF EOF(1) <> 0 THEN PRINT "no data points in file "; filename$; " - exiting": END
CALL inputline(1, d1#(), col%, d2$(), col2%)
'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
'make fit
 CALL simanneal(firstvalue, par#(), sta#, ii%, jj%, nmax, speed)'this is not good
CLOSE 1
  'CALL simplex(firstvalue, par#(), sta#, ii%, jj%)
IF firstvalue = -1 THEN PRINT "error fitting " + filename$: GOTO 998

'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

'  write to output file  ddddddddddddddddddddddddddd
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "FIT1LORZ.lin"
PRINT #2, "#{";
PRINT #2, "fitp:";
FOR i% = 1 TO par#(0): PRINT #2, USING "##.###^^^^"; par#(i%); : NEXT i%
PRINT #2, USING "sta=##.###^^^^"; sta#
FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, "#"; DATE$; " "; TIME$; " column"; jj%; "FIT1LORZ calc. with respect to column"; ii%;
PRINT #2, "#fit put into column"; col% + 1;
PRINT #2, "calculation done by program FIT1LORZ.bas";
PRINT #2, "}"

i = 0: colmax% = 0: WHILE EOF(1) = 0
CALL inputline(1, d1#(), col%, d2$(), col2%)
IF col% > colmax% THEN colmax% = col%

i = i + 1: FOR coll% = 1 TO col%: PRINT #2, d1#(coll%);
av#(coll%) = ((i - 1) * av#(coll%) + d1#(coll%)) / i
NEXT

CALL calcfunction(par#(), y#, d1#(ii%)): PRINT #2, y#;
FOR coll% = 1 TO col2%: PRINT #2, " {" + d2$(coll%) + "}"; : NEXT: PRINT #2,
WEND
'ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

CLOSE 1, 2

'loglogloglogloglogloglogloglogloglogloglogloglogloglogloglogloglogloglog
'write fit to logfile
OPEN "a", 1, "fit1LORZ.log"
IF LOF(1) < 2 THEN PRINT #1, "#{par#(1) vs par#(2) vs ... vs sta vs av(1)(average of columns 1) vs av(2) vs ... vs filename (logfile of fit1LORZ.log)}"
FOR i% = 1 TO par#(0): PRINT #1, par#(i%); : NEXT i%
PRINT #1, sta#;
FOR i% = 1 TO colmax%: PRINT #1, av#(i%); : av#(i%) = 0: NEXT i%
PRINT #1, " {" + filename$ + "}"
CLOSE 1
'loglogloglogloglogloglogloglogloglogloglogloglogloglogloglogloglogloglog

SHELL "copy FIT1LORZ.lin " + filename$
SHELL "del FIT1LORZ.lin"
PRINT
PRINT "END FIT1LORZ in file "; filename$;
PRINT " column"; jj%; "FIT1LORZ calc. with respect to column"; ii%;
PRINT " fit put into column"; col% + 1

998 IF INSTR(COMMAND$, "*") <> 0 GOTO 999

END
333 FOR i = 1 TO 14: READ a$: PRINT a$: NEXT i: END

SUB calcfunction (par#(), y#, x#)
' calculate the function y#=y#(x#,par#(1),par#(2),...)

'one Lorentzian
par#(0) = 3   '3 parameters:
              '1...center 2...fwhm 3...height
              '(if /B switch is set: 4...constant background)

IF par#(2) = 0 THEN
 y# = 0
ELSE
 dx# = 2 * (x# - par#(1)) / par#(2)
  y# = par#(3) / (1 + dx# * dx#)
END IF

IF INSTR(COMMAND$, "/B") > 0 THEN par#(0) = 4: y# = y# + par#(4)

END SUB

SUB calcsta (gr%, firstvalue, par#(), sta#, ii%, jj%, xmin, xmax, ymin, ymax)
STATIC timesav
DIM d1#(30), d2$(10)
'sub to calculate standard deviation sta#
REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
SEEK #1, firstvalue
 IF INSTR(COMMAND$, "/G") > 0 AND TIMER - timesav > .5 THEN gr% = 1: timesav = TIMER
sta# = 0
'input data and calculate linear regression  iiiiiiiiiiiiiiiiiiiii
WHILE EOF(1) = 0
CALL inputline(1, d1#(), col1%, d2$(), col2%)

CALL calcfunction(par#(), yy#, d1#(ii%))
d# = yy# - d1#(jj%)
sta# = sta# + d# * d#
x = d1#(ii%): y = d1#(jj%)
IF gr% = 1 THEN CALL pointdraw3d(x, xmin, xmax, y, ymin, ymax, z, 1, 1, 11.05)
WEND

IF gr% = 1 THEN
  FOR i% = 0 TO 100
   xx# = xmin + i% * (xmax - xmin) / 100
   CALL calcfunction(par#(), yy#, xx#)
   x = xx#: y = yy#
   IF i% = 0 THEN
    CALL pointdraw3d(x, xmin, xmax, y, ymin, ymax, 0, 1, 1, 1.03)
   ELSE
    CALL pointdraw3d(x, xmin, xmax, y, ymin, ymax, 0, 1, 1, -11.23)
   END IF
  NEXT
  LOCATE 1, 1: PRINT "F"
END IF

END SUB

SUB calcstartvalues (firstvalue, par#(), parstp(), parmax(), parmin(), ii%, jj%, xmin, xmax, ymin, ymax)
DIM d1#(30), d2$(10)
REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'if less than 5 values are there - firstvalue returns -1

 SEEK #1, firstvalue

 xmin = 1E+20
 xmax = -1E+20
 ymin = 1E+20
 ymax = -1E+20
 n = 0: bkg = 0
 WHILE EOF(1) = 0
 CALL inputline(1, d1#(), col1%, d2$(), col2%)

 IF d1#(ii%) < xmin THEN xmin = d1#(ii%)
 IF d1#(ii%) > xmax THEN xmax = d1#(ii%)
 IF d1#(jj%) < ymin THEN ymin = d1#(jj%)
 IF d1#(jj%) > ymax THEN ymax = d1#(jj%): xmid = d1#(ii%)
 IF n = 0 OR bkg = 0 THEN weight = 1 ELSE weight = EXP(-(d1#(jj%) - bkg) * (d1#(jj%) - bkg) / bkg / bkg)

 IF n + weight <> 0 THEN bkg = (bkg * n + weight * d1#(jj%)) / (n + weight)
 IF bkg = 0 THEN bkg = 1
 n = n + weight
 WEND
IF xmax - xmin < 1E-10 THEN PRINT "error - all points have equal x": END

FOR trytwice% = 1 TO 3

222 SEEK #1, firstvalue
sumxy = 0
sumxxy = 0
sumy = 0
sumyy = 0
n = 0
area = 0
'input data and calculate linear regression  iiiiiiiiiiiiiiiiiiiii
WHILE EOF(1) = 0
CALL inputline(1, d1#(), col1%, d2$(), col2%)

 IF INSTR(COMMAND$, "/B") > 0 THEN
  d1#(jj%) = d1#(jj%) - bkg
 END IF
 IF trytwice% = 1 OR ABS(d1#(ii%) - par#(1)) < fwhm * 1 THEN
  sumy = sumy + d1#(jj%)
  sumyy = sumyy + d1#(jj%) * d1#(jj%)
  sumxy = sumxy + d1#(ii%) * d1#(jj%)
  sumxxy = sumxxy + d1#(ii%) * d1#(ii%) * d1#(jj%)
  IF n > 0 THEN dx = ABS(d1#(ii%) - xold): area = area + d1#(jj%) * dx
  n = n + 1
  xold = d1#(ii%)
 END IF
WEND
IF trytwice% = 1 AND n < 8 THEN firstvalue = -1: GOTO 6699

IF INSTR(COMMAND$, "/B") > 0 THEN
 IF n < 2 THEN bkg = bkg * .9: fwhm = xmax - xmin: GOTO 222
 IF bkg <> 0 AND n <> 0 THEN
  IF ABS(sumy / bkg / n) < .01 THEN bkg = bkg * .9: fwhm = xmax - xmin: GOTO 222
 END IF
END IF

'2 ... fwhm
IF sumy <> 0 THEN sigma = SQR(ABS((sumxxy / sumy - sumxy * sumxy / sumy / sumy)))
fwhm = SQR(8 * LOG(2)) * sigma
IF fwhm > 0 THEN par#(2) = fwhm ELSE fwhm = .1
IF n > 0 THEN parstp(2) = sigma / SQR(2 * n) ELSE parstp(2) = fwhm / 10
parmax(2) = xmax - xmin
parmin(2) = 0

IF sigma = 0 AND INSTR(COMMAND$, "/B") > 0 THEN bkg = bkg * .9: fwhm = xmax - xmin: GOTO 222


'1...center
'par#(1) = RND(1) * (xmax - xmin) + xmin
'par#(1) = xmid
IF sumy <> 0 THEN par#(1) = sumxy / sumy ELSE par#(1) = sumxy
parmax(1) = xmax
parmin(1) = xmin
IF n > 0 THEN parstp(1) = fwhm / SQR(8 * LOG(2) * n) ELSE parstp(1) = fwhm / 10

'3 ... height
par#(3) = area / fwhm * 2 * SQR(LOG(2) / 3.1415)
parstp(3) = .1 * par#(3)
parmax(3) = 1.5 * par#(3)
parmin(3) = .5 * par#(3)


'4 bkground
IF INSTR(COMMAND$, "/B") > 0 THEN
 par#(4) = bkg
 parstp(4) = .1 * ABS(sumy / n)
 parmax(4) = bkg + SQR(sumyy) / n
 parmin(4) = bkg - SQR(sumyy) / n
END IF
IF area < 0 THEN bkg = 0
NEXT trytwice%

6699 END SUB

SUB getpar (x, x$, a$)
'gets parameter x (name: x$) out of a$ (if present) - case match)
x = 0
in% = INSTR(a$, x$)
IF in% > 0 THEN
 in% = in% + LEN(x$)
 b$ = LTRIM$(MID$(a$, in%))
 IF LEFT$(b$, 1) = "=" THEN b$ = MID$(b$, 2)
 x = VAL(b$)
END IF


END SUB

SUB headerinput (text$(), j, n)
'*********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'*********************************************************************
   j = 0
   IF EOF(1) <> 0 THEN PRINT "no data points in file "; filename$; " - exiting": END
   LINE INPUT #n, a$
   IF LEFT$(LTRIM$(a$), 1) = "#" THEN i = 0: GOTO 2
   i = INSTR(a$, "{")    'look for "{"
   IF i = 0 THEN SEEK #1, 1: GOTO 1
2  j = j + 1: text$(j) = RIGHT$(a$, LEN(a$) - i)
   k = INSTR(text$(j), "}"): IF k > 0 THEN i = k: GOTO 3'look for "}" in first line
   k = j + 1
   FOR j = k TO 300
   IF EOF(1) <> 0 THEN PRINT "no data points in file "; filename$; " - exiting": END
   f = SEEK(1): LINE INPUT #n, text$(j)
   IF LEFT$(LTRIM$(text$(j)), 1) <> "#" AND i = 0 THEN SEEK #1, f: j = j - 1: GOTO 1
   k = INSTR(text$(j), "}"): IF k > 0 THEN i = k: GOTO 3   'look for "}"
   NEXT j: PRINT "text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

1 END SUB

SUB inputline (n, d1#(), col1%, d2$(), col2%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
'comments in {} are stored in d2$
'if a comment is started somewhere in this line by "{" and not finished
'then col1% is set -1 and the filepointer is set to the beginning of the
'line (with seek)

a$ = INKEY$: IF a$ <> "" THEN IF ASC(a$) = 27 THEN END


aa = SEEK(n)
LINE INPUT #n, ala$
IF LEFT$(LTRIM$(ala$), 1) = "#" THEN col1% = 0: col2% = 1: d2$(1) = ala$: GOTO 111

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

111 END SUB

SUB iteratefiles (file1$)
STATIC washere%, path$

  IF INSTR(COMMAND$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere% = 0 THEN
      washere% = 1
      path$ = ""
      IF INSTR(file1$, "\") > 0 THEN
       i% = INSTR(file1$, "\")
       WHILE INSTR(i% + 1, file1$, "\") <> 0: i% = INSTR(i% + 1, file1$, "\"): WEND
       path$ = LEFT$(file1$, i%)
      END IF
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > fact.dir"
      OPEN "i", 9, "fACT.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del fact.dir": END
   INPUT #9, file1$
   IF LCASE$(file1$) = "fact.dir" THEN
    IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del fact.dir": END
    INPUT #9, file1$
   END IF
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del fact.dir": END
   file1$ = path$ + file1$
  END IF


END SUB

SUB pointdraw3d (x, xmin, xmax, y, ymin, ymax, z, zmin, zmax, size)
STATIC mins(), maxs(), i, s()
REM*************************************************************************
REM*************************************************************************
REM diese routine zeichnet einen punkt der groesse size
REM in einer 3dgrafik und mahlt die skala, wenn sie das erste mal aufgerufen
REM wird bzw. wenn die skala sich aendert werden auch alle anderen punkte
REM neu skaliert und gemahlt
REM nach dem aufruf erscheint in der linke oberen ecke des bildschrims
REM eine 0, wird diese ueberschrieben, so wird bei einem neuerlichen
REM aufruf der routine ein neuer plot angefangen und die alten punkte
REM vergessen
REM
REM eingabe
REM x,y,z               punktkoordinaten
REM x,y,z,min,xyzmax       bereich (ist xyzmin>xyzmax, so erfolgt automatische
REM                                     bereichswahl,
REM                                     ist zmin=zmax so folgt 2d darstellung)
REM size                punktgroesse (>0), (<0..zeichnet linie von vorhergehendem punkt)
      'z.B 3.14 means color 3 and circle diameter 14

REM some fix parameters of the routine XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
alpha = 45 / 180 * 3.1415: REM winkel des schraegrisses
zsmax = 210: REM laenge der z-achse
x$ = " ": y$ = " ": z$ = " " 'achsenbeschriftung
height = 360: widt = 620: REM hoehe und breite des bildes
h0 = 40: w0 = 0: REM linke obere ecke des bildes
'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

REM*************************************************************************
REM*************************************************************************
DIM min(3), max(3): min(1) = xmin: min(2) = ymin: min(3) = zmin 'rename xyzmin
                    max(1) = xmax: max(2) = ymax: max(3) = zmax 'rename xyzmax

'ccccccccc calculate maximum length of x,y (z) axis on the screen cccccccccc
IF min(3) = max(3) THEN zsmax = 20 / COS(alpha)  'wenn 2dim plot dann zsmax fast 0
xsmax = widt - zsmax * SIN(alpha)
ysmax = height - 10 - zsmax * COS(alpha)
'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

'iiiii this is only done once on first call(initialize plotroutine) iiiiiiiiiiiiii
IF i = 0 THEN
 DIM s(3, 1000), mins(3), maxs(3)' dimension array(static) for saving data points
 SCREEN 12, 6:  WIDTH 80, 60: CLS
 VIEW SCREEN (w0, h0 + 15)-(w0 + widt + 11, h0 + height + 15), 0, 3
                       'setup output screen
END IF 'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

'look if pointdraw should make new plot or add point to old one
IF SCREEN(1, 1, 0) <> 48 THEN makenewgraph% = 1 ELSE makenewgraph% = 0
IF makenewgraph% = 1 THEN i = 0

i = i + 1: 'a new point is to be drawn, give it number i
s(1, i) = x: s(2, i) = y: s(3, i) = z: s(0, i) = size   'store the new point

REM abfrage automatic scale for x y and z axis  sssssssssssssssssssssssssssss
FOR axis% = 1 TO 3
IF min(axis%) > max(axis%) THEN
 min(axis%) = s(axis%, 1): max(axis%) = s(axis%, 1)
 FOR j = 1 TO i 'calculate extremums
  IF s(axis%, j) > max(axis%) THEN max(axis%) = s(axis%, j)
  IF s(axis%, j) < min(axis%) THEN min(axis%) = s(axis%, j)
 NEXT: IF min(axis%) = max(axis%) THEN max(axis%) = max(axis%) + 10 ^ 20
 delta = (max(axis%) - min(axis%)) / 9: exponent% = INT(LOG(delta) / LOG(10))
 delta = (1 + INT(delta / 10 ^ exponent%)) * 10 ^ exponent%'round delta to sensible number
 IF min(axis%) <> 0 THEN exponentmin% = INT(LOG(ABS(min(axis%))) / LOG(10)) ELSE exponentmin% = -100
 IF max(axis%) <> 0 THEN exponentmax% = INT(LOG(ABS(max(axis%))) / LOG(10)) ELSE exponentmax% = -100
 expon% = exponentmin%: IF exponentmin% < exponentmax% THEN expon% = exponentmax%
 IF expon% < eponent% + 5 THEN min(axis%) = INT(min(axis%) / delta) * delta
 max(axis%) = min(axis%) + 10 * delta
END IF
NEXT axis% 'Sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

redrawscale% = makenewgraph% 'redrawscale% tells if the scales has to be redrawn-------
FOR axis% = 1 TO 3
IF min(axis%) <> mins(axis%) OR max(axis%) <> maxs(axis%) THEN redrawscale% = 1
  mins(axis%) = min(axis%): maxs(axis%) = max(axis%)
NEXT axis% '----------------------------------------------------------------

IF redrawscale% = 1 THEN
 CLS : COLOR 15: VIEW PRINT: LOCATE 1, 1: PRINT "0";  'skala machen
             
 ' labels and scale for x,y axis xxxxxyyyyyyyyyyyyyyyxxxxxxxxxxxxyyyyyyyyyy
 LOCATE CINT((h0 + height + 8) * .125), 1: PRINT USING "##.##^^^^"; min(1);
  FOR j = .2 TO 1.1 STEP .2
   LOCATE CINT((h0 + height + 8) * .125), CINT(.12 * w0 + .12 * j * xsmax): PRINT USING "##.##^^^^"; min(1) + j * (max(1) - min(1));
  NEXT: PRINT x$;
 LOCATE CINT((h0 + height - 8) * .125), 1 + CINT(.12 * w0): PRINT USING "##.####^^^^"; min(2);
 LOCATE CINT((h0 + height) * .125) - CINT(.125 * ysmax), 1 + CINT(.12 * w0): PRINT USING "##.####^^^^"; max(2); : PRINT y$;

 LINE (5 + w0, h0 + height - 5)-(w0 + xsmax, h0 + height - 5), 3
 LINE (5 + w0, h0 + height - 5 - ysmax)-(w0 + 5, h0 + height - 5), 3
 FOR j = 0 TO 1.01 STEP .1
 LINE (w0 + xsmax * j + 5, h0 + height - 5)-(w0 + xsmax * j + 5, h0 + height - 2), 3
 LINE (w0 + 2, h0 + height - 5 - j * ysmax)-(w0 + 5, h0 + height - 5 - j * ysmax), 3
 NEXT 'xxxxxxxxxyyyyyyyyyyyyyxxxxxxxxxxxxxxxxyyyyyyyyyyyyyyyxxxxxxxxxxyyyyy

 IF min(3) <> max(3) THEN  'scale z-axis zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
  LOCATE CINT((h0 + height) * .132) + 1 - CINT(.16 * zsmax * COS(alpha)), CINT(.12 * (w0 + zsmax * SIN(alpha)))
  PRINT USING "##.#^^^^"; max(3); : PRINT z$;
  LOCATE CINT((h0 + height) * .132) + 1 - CINT(.1 * zsmax * COS(alpha)), 1 + CINT(.12 * w0): 'CINT(.07 * zsmax * SIN(alpha)) - 7
  PRINT USING "##.#^^^^"; min(3) + .5 * (max(3) - min(3))
  LINE (w0 + 5, h0 + height - 5)-(w0 + 5 + SIN(alpha) * zsmax, h0 + height - 5 - COS(alpha) * zsmax), 3
  FOR j = 0 TO 1.01 STEP .1
  LINE (w0 + 2 + j * zsmax * SIN(alpha), h0 + height - 5 - COS(alpha) * j * zsmax)-(5 + j * zsmax * SIN(alpha), height - 5 - COS(alpha) * j * zsmax), 3
  NEXT
 END IF 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

  FOR j = 1 TO i: GOSUB 3199: NEXT j 'redraw old points if there are
END IF

j = i: GOSUB 3199  'draw new point i if there are already more than one
LOCATE CINT(height * .123) + 1, 1
GOTO 3200   'goto end sub

'$$$$$$$$$$$$$$$$$$$$$$$$$  GOSUB drawpoint $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
3199  '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
REM zeichnen des punktes nummer j drawpoint $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
     'calculate position of point xp,yp on screen
IF min(3) <> max(3) THEN    'make 3d-plot
 xp = 5 + (s(1, j) - mins(1)) / (maxs(1) - mins(1)) * xsmax + (s(3, j) - mins(3)) / (maxs(3) - mins(3)) * SIN(alpha) * zsmax
 yp = height - 5 - (s(2, j) - mins(2)) / (maxs(2) - mins(2)) * ysmax - (s(3, j) - mins(3)) / (maxs(3) - mins(3)) * COS(alpha) * zsmax
ELSE   'make 2d plot
 xp = 5 + (s(1, j) - mins(1)) / (maxs(1) - mins(1)) * xsmax
 IF maxs(2) = mins(2) THEN maxs(2) = mins(2) + 1
 yp = height - 5 - (s(2, j) - mins(2)) / (maxs(2) - mins(2)) * ysmax
END IF

IF s(0, j) >= 0 THEN      'draw circle
CIRCLE (w0 + xp, h0 + yp), (s(0, j) - INT(s(0, j))) * 100, INT(s(0, j))
ELSE
 IF j > 1 THEN            'draw line
  IF min(3) <> max(3) THEN
   xpm = 5 + (s(1, j - 1) - mins(1)) / (maxs(1) - mins(1)) * xsmax + (s(3, j - 1) - mins(3)) / (maxs(3) - mins(3)) * SIN(alpha) * zsmax
   ypm = height - 5 - (s(2, j - 1) - mins(2)) / (maxs(2) - mins(2)) * ysmax - (s(3, j - 1) - mins(3)) / (maxs(3) - mins(3)) * COS(alpha) * zsmax
  ELSE
   xpm = 5 + (s(1, j - 1) - mins(1)) / (maxs(1) - mins(1)) * xsmax
   ypm = height - 5 - (s(2, j - 1) - mins(2)) / (maxs(2) - mins(2)) * ysmax
  END IF
  LINE (w0 + xp, h0 + yp)-(w0 + xpm, h0 + ypm), INT(s(0, j))
 END IF
END IF
RETURN '$$$$$$$$$$$$$$$$$$ end sub drawpoint $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

3200 END SUB

SUB simanneal (firstvalue, par#(), sta#, ii%, jj%, nmax, speed)

DIM parstp(30), parmax(30), parmin(30), parsav(30)




'***********************************************************************
 ' simulated annealing algorithm

'startwerte

  CALL calcstartvalues(firstvalue, par#(), parstp(), parmax(), parmin(), ii%, jj%, xmin, xmax, ymin, ymax)
    IF firstvalue = -1 GOTO 11
  gr% = 0: CALL calcsta(gr%, firstvalue, par#(), sta#, ii%, jj%, xmin, xmax, ymin, ymax)
  stattemp = sta# * 1: stps = 1

FOR n = 1 TO nmax
  stasave = sta#
  FOR i = 1 TO par#(0): parsav(i) = par#(i): NEXT i
  FOR i = 1 TO par#(0): par#(i) = par#(i) + (RND(1) - .5) * parstp(i)
   IF par#(i) > parmax(i) THEN par#(i) = parmax(i)
   IF par#(i) < parmin(i) THEN par#(i) = parmin(i)
  NEXT i
   rndm! = RND(1)
   gr% = 0: CALL calcsta(gr%, firstvalue, par#(), sta#, ii%, jj%, xmin, xmax, ymin, ymax)
   
  IF sta# > stasave THEN
   IF rndm! > EXP(-(sta# - stasave) / stattemp) THEN
     'simulated annealing -recover old paramters
     FOR i = 1 TO par#(0): par#(i) = parsav(i): NEXT i
     sta# = stasave
     FOR i = 1 TO par#(0)
      parstp(i) = parstp(i) * .95: stps = stps * .95 ^ speed
     NEXT i
   ELSE
    stattemp = stattemp * .95 ^ speed
   END IF
  ELSE
   stattemp = stattemp * .9 ^ speed
  END IF
  stav = (stav * (n - 1) + sta#) / n
IF (n - 1) / 10 = INT((n - 1) / 10) THEN LOCATE 1, 1: PRINT USING "sta=##.####^^^^"; sta#;
IF stattemp < stav * .001 GOTO 11
NEXT n
 
  gr% = 1: CALL calcsta(gr%, firstvalue, par#(), sta#, ii%, jj%, xmin, xmax, ymin, ymax)
  IF INSTR(COMMAND$, "/C") = 0 THEN INPUT ala

11 END SUB

