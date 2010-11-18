DECLARE SUB analizecommand (file1$)
DECLARE SUB calcsta (x#(), y2#(), ix%, iy%, nofpoints%, simp!(), i!, paranz%)
DECLARE SUB splint (xx#, xl#, xh#, yl#, yh#, y2l#, y2h#, yy#)
DECLARE SUB fitspline (x#(), y2#(), ix%, iy%, nofpoints%)
DECLARE SUB y2calc (x#(), y2#(), ix%, iy%, n%, yp1#, ypn#, sta!)
DECLARE SUB inputline (n!, xm#(), col%)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "SPLINE SPLINE SPLINE SPLINE SPLINE SPLINE SPLINE SPLINE SPLINE SPLINE"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program spline - use it like 'SPLINE *.* 2 spl=3'                      "
DATA "    (means take file *.* second column as xaxis 3rd col as yaxis and calc."
DATA "    spline interpolation curvature, add a column and store it there) -    "
DATA " ---> the result is written in file *.*, Note the file must be            "
DATA " sorted before starting this programm according to ascending or descending"
DATA " values of the xaxis -column                                              "
DATA " format of file                                                           "
DATA "                                                                          "
DATA " { header: this is the                                                    "
DATA "   file header delimited                                                  "
DATA "   by brackets after this header there follow 2 or more data columns }    "
DATA " 11 3.14235 65367                                                         "
DATA "  .    .     .                                                            "
DATA "  .    .     .                                                            "
DATA "  .    .     .    .  .   .                                                "
DATA "  .    .     .                                                            "
DATA " 32 2412.34 324.2                                                         "
DATA "                                                                          "
DATA "************************************************************************* "
DIM text$(100), xm#(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(a$, INSTR(a$, " "))): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
ix% = ASC(LEFT$(a$, 1)) MOD 48: iy% = ASC(MID$(a$, INSTR(a$, "spl=") + 4, 1)) MOD 48
                                'get x and y columns
999 CALL analizecommand(filename$)

IF RIGHT$(filename$, 4) = ".rcp" THEN
   PRINT "you should never ever change data in *.rcp files": PLAY "dgdgdg"
222 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 222
END IF
'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

REM input data columns iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
CALL inputline(1, xm#(), coll%)
nofpoints% = 1: col% = coll%: DIM x#(3000, 1)
x#(nofpoints%, 1) = xm#(ix%): x#(nofpoints%, 0) = xm#(iy%)

12 WHILE EOF(1) = 0
 CALL inputline(1, xm#(), col%)
IF xm#(ix%) = x#(nofpoints%, 1) GOTO 12  'ignore points with same x values
IF coll% <> col% GOTO 110 'see if we need to stop inputting
nofpoints% = nofpoints% + 1: IF nofpoints% > 3000 THEN PRINT "error spline: too many datapoints(>1000)": END
x#(nofpoints%, 1) = xm#(ix%): x#(nofpoints%, 0) = xm#(iy%)
WEND
110 CLOSE
'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
DIM y2#(nofpoints%, 0)


'calculate spline ccccccccccccccccccccccccccccccccccccccccccccccccccc
CALL fitspline(x#(), y2#(), 1, 0, nofpoints%)
'cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 'write result to file wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
'open input file and input file header (just to get the other columns right)
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)
' open output file and write fileheader
OPEN "o", 2, "spline.spl"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " column "; ix%; " and "; iy%; " were taken as the x and y axis for calculating"
PRINT #2, "a cubic spline - for the corresponding curvature values a "; coll% + 1; ". column was created by program spline.bas}"
FOR i% = 1 TO nofpoints%:
112 CALL inputline(1, xm#(), col%): IF xm#(ix%) <> x#(i%, 1) GOTO 112
                                      'ignore points with same x values
    FOR k% = 1 TO coll%: PRINT #2, xm#(k%); : NEXT: PRINT #2, y2#(i%, 0)
NEXT i%
CLOSE 'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

SHELL "copy spline.spl " + filename$
SHELL "del spline.spl"

PRINT
PRINT "END SPLINE in file "; filename$; " column "; ix%; "-"; iy%; " were taken as the x-y axis"
PRINT "for calculating a cubic spline - for the corresponding curvature values a ";
PRINT USING "##."; coll% + 1; : PRINT " column was created": PRINT
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END
333 FOR i = 1 TO 20: READ a$: PRINT a$: NEXT i: END

SUB analizecommand (file1$)
STATIC washere%

  IF INSTR(COMMAND$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere% = 0 THEN
      washere% = 1
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
  END IF


END SUB

SUB calcsta (x#(), y2#(), ix%, iy%, nofpoints%, simp(), i, paranz%)
'*******************************************
' this sub is only needed to keep track of variables in the
' simplex algorithm - it calls directly the y2calc to calculate
' standard deviation
'*******************************************
sl# = simp(i, 1)
sh# = simp(i, 2)
CALL y2calc(x#(), y2#(), ix%, iy%, nofpoints%, sl#, sh#, sta)
simp(i, 3) = sta

END SUB

SUB fitspline (x#(), y2#(), ix%, iy%, nofpoints%)
'************************************************************************
' this sub fits the slope of a splinefunction at upper and lower
' boundary of the dataset {x#(,ix%) vs. x#(,iy%)} such
' that the standard deviation of the splinefunction in each interval
' from the linear interpolation is minimized, simplex algorithm is
' used
'************************************************************************
'calculate slope of linear interpolation at end of dataset as initial
' values
DIM s(3)
s(1) = (x#(2, iy%) - x#(1, iy%)) / (x#(2, ix%) - x#(1, ix%)) 'slopelow
s(2) = (x#(nofpoints%, iy%) - x#(nofpoints% - 1, iy%)) / (x#(nofpoints%, ix%) - x#(nofpoints% - 1, ix%))
                                'slopehigh
'----------------------------------------------------------------------

REM****************SIMPLEX ALGORITHMUS**************************************
paranz% = 2
pap1% = paranz% + 1

DIM h(pap1%), l(pap1%), nxt(1, pap1%), ctr(pap1%)
DIM mean(pap1%), er(pap1%), maxerr(pap1%), p(pap1%)
DIM q(pap1%), stp(pap1%)
DIM simp(pap1%, pap1%)
alfa = 1!: REM { reflection coeff.>0}
beta = .5: REM { contraction coeff.0to1}
gamma = 2!: REM { expansion coeff.>1}
maxiter = 100 'maximum number od iterations

REM anfangswerte
FOR i% = 1 TO paranz%: simp(1, i%) = s(i%): NEXT

REM anfangsschrittweiten
FOR i% = 1 TO paranz%: stp(i%) = s(i%) / 10: IF s(i%) = 0 THEN stp(i%) = .01
NEXT

REM abbruchfehler
FOR i% = 1 TO paranz%: maxerr(i%) = stp(i%) / 100: NEXT
maxerr(pap1%) = .001: REM sta-maxerr

CALL calcsta(x#(), y2#(), ix%, iy%, nofpoints%, simp(), 1, paranz%)
   FOR i = 1 TO paranz%: REM {compute offset of the vertexes of the start simplex}
      p(i) = stp(i) * (SQR(pap1%) + 1) / (paranz% * SQR(2))
      q(i) = stp(i) * (SQR(pap1%) - 1) / (paranz% * SQR(2))
   NEXT i
   FOR i = paranz% TO pap1%: REM     {all vertexes of the starting simplex}
           FOR j = 1 TO paranz%: simp(i, j) = simp(1, j) + q(j): NEXT j
           simp(i, i - 1) = simp(1, i - 1) + p(i - 1)
           CALL calcsta(x#(), y2#(), ix%, iy%, nofpoints%, simp(), i, paranz%)
   NEXT i
   FOR i = 1 TO pap1%: l(i) = 1: h(i) = 1: NEXT i
     FOR j = 1 TO pap1%: FOR i = 1 TO pap1%
     IF simp(i, j) < simp(l(j), j) THEN l(j) = i
     IF simp(i, j) > simp(h(j), j) THEN h(j) = i
     NEXT i: NEXT j
   niter = 0
10 REM Repeat            { Hauptschleife }
   done = 1: REM 1..true
   niter = niter + 1
   FOR i = 1 TO pap1%: ctr(i) = 0!: NEXT i
   FOR i = 1 TO pap1%: REM     {compute centroid}
   IF i <> h(pap1%) THEN FOR j = 1 TO pap1%: ctr(j) = ctr(j) + simp(i, j): NEXT j
   NEXT i
   FOR i = 1 TO pap1%: ctr(i) = ctr(i) / paranz%
   nxt(1, i) = (1! + alfa) * ctr(i) - alfa * simp(h(pap1%), i)
   NEXT i
           CALL calcsta(x#(), y2#(), ix%, iy%, nofpoints%, nxt(), 1, paranz%)
IF nxt(1, pap1%) >= simp(l(pap1%), pap1%) GOTO 20
FOR i = 1 TO pap1%: simp(h(pap1%), i) = nxt(1, i): NEXT i
FOR i = 1 TO paranz%: nxt(1, i) = gamma * simp(h(pap1%), i) + (1! - gamma) * ctr(i): NEXT i
           CALL calcsta(x#(), y2#(), ix%, iy%, nofpoints%, nxt(), 1, paranz%)
IF nxt(1, pap1%) <= simp(l(pap1%), pap1%) THEN FOR i = 1 TO pap1%: simp(h(pap1%), i) = nxt(1, i): NEXT i: REM  {expansion accepted}
GOTO 30
20 IF nxt(1, pap1%) <= simp(h(pap1%), pap1%) THEN FOR i = 1 TO pap1%: simp(h(pap1%), i) = nxt(1, i): NEXT i: GOTO 30
FOR i = 1 TO paranz%: nxt(1, i) = beta * simp(h(pap1%), i) + (1! - beta) * ctr(i): NEXT i
           CALL calcsta(x#(), y2#(), ix%, iy%, nofpoints%, nxt(), 1, paranz%)
IF nxt(1, pap1%) <= simp(h(pap1%), pap1%) THEN FOR i = 1 TO pap1%: simp(h(pap1%), i) = nxt(1, i): NEXT i: GOTO 30: REM  {contraction accepted}
FOR i = 1 TO pap1%: FOR j = 1 TO paranz%
 simp(i, j) = (simp(i, j) + simp(l(pap1%), j)) * beta
 NEXT j
           CALL calcsta(x#(), y2#(), ix%, iy%, nofpoints%, simp(), i, paranz%)
NEXT i
30 FOR j = 1 TO pap1%: FOR i = 1 TO pap1%
    IF simp(i, j) < simp(l(j), j) THEN l(j) = i
    IF simp(i, j) > simp(h(j), j) THEN h(j) = i
   NEXT i: NEXT j
   FOR j = 1 TO pap1%
     IF ABS(simp(h(j), j)) < 1E-20 THEN simp(h(j), j) = 1E-20
     er(j) = (simp(h(j), j) - simp(l(j), j)) / simp(h(j), j)
     IF er(j) > maxerr(j) THEN done = 0
   NEXT j
IF done <> 1 AND niter < maxiter GOTO 10: REM Ende der Hauptschleife}

   FOR i = 1 TO pap1%: mean(i) = 0!
        FOR j = 1 TO pap1%: mean(i) = mean(i) + simp(j, i): NEXT j
        mean(i) = mean(i) / pap1%: REM mean(i)....angepasste parameter
      'PRINT mean(i)
   NEXT i
REM ***********************END SIMPLEX**************************************
   
    ' calculate sta with fitted parameters
FOR i = 1 TO pap1%: simp(1, i) = mean(i): NEXT
CALL calcsta(x#(), y2#(), ix%, iy%, nofpoints%, simp(), 1, paranz%)
END SUB

SUB headerinput (text$(), j, n)
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'**********************************************************************

1 INPUT #n, a$
   i = INSTR(a$, "{"): IF i > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(a$, LEN(a$) - i)
  j = 1: i = INSTR(text$(j), "}"): IF i > 0 GOTO 3  'look for "}" in first line
                                           
   FOR j = 2 TO 100
   INPUT #n, text$(j)
   i = INSTR(text$(j), "}"): IF i > 0 GOTO 3      'look for "}"
   NEXT j: PRINT "text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

END SUB

SUB inputline (n, xm#(), col%)
'**********************************************************************
'diese sub liest eine zeile von file #n, bestimmt die
' anzahl der datenspalten col% und sichert die daten in xm#(1,..,col%)
'**********************************************************************
 INPUT #n, ala$        'input data point line as string and split into numbers
          col% = 0: WHILE LEN(ala$) > 0: ala$ = LTRIM$(ala$) + " ": col% = col% + 1
 xm#(col%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
          ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " "))): WEND
END SUB

SUB splint (xx#, xl#, xh#, yl#, yh#, y2l#, y2h#, yy#)
'************************************************************************
' this sub calculates a cubic splinefunction at position xx# in the interval
' [xl#,xh#] with given values [yl#,yh#] and curvature [y2l#,y2h#]
' at the endpoints in the interval
' the result is calculated as yy#
'************************************************************************
h# = xh# - xl#
a# = (xh# - xx#) / h#
B# = (xx# - xl#) / h#
yy# = a# * yl# + B# * yh# + ((a# * a# * a# - a#) * y2l# + (B# * B# * B# - B#) * y2h#) * h# * h# / 6!

END SUB

SUB y2calc (x#(), y2#(), ix%, iy%, n%, yp1#, ypn#, sta)
'*************************************************************************
' this sub calculates the curvature y2#(,iy%) of a cubic spline function
' for a dataset of n% points x#(,ix%) vs. x#(,iy%) with
' slope yp1# (ypn#) at the lower(upper) end
' the standard deviation sta of the interval midpoints interpolated
' linear and using cubic spline is calculated
'*************************************************************************
DIM u#(n%)
'-------------------------------------------------------------------------
y2#(1, iy%) = -.5    'the lower boundary condition is set to have specified
                     ' first derivative
u#(1) = 3! / (x#(2, ix%) - x#(1, ix%)) * ((x#(2, iy%) - x#(1, iy%)) / (x#(2, ix%) - x#(1, ix%)) - yp1#)
'--------------------------------------------------------------------------

'lllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
' decomposition loop of the tridiagonal algorithm, y2# and u# are used
' for temporary storageof the decomposed factors
FOR i% = 2 TO n% - 1
 sig# = (x#(i%, ix%) - x#(i% - 1, ix%)) / (x#(i% + 1, ix%) - x#(i% - 1, ix%))
 p# = sig# * y2#(i% - 1, iy%) + 2!
 y2#(i%, iy%) = (sig# - 1!) / p#
 u#(i%) = (x#(i% + 1, iy%) - x#(i%, iy%)) / (x#(i% + 1, ix%) - x#(i%, ix%)) - (x#(i%, iy%) - x#(i% - 1, iy%)) / (x#(i%, ix%) - x#(i% - 1, ix%))
 u#(i%) = (6! * u#(i%) / (x#(i% + 1, ix%) - x#(i% - 1, ix%)) - sig# * u#(i% - 1)) / p#
NEXT i%
'lllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll

'-------------------------------------------------------------------------
' the upper boundary condition is set to have specified first derivative
qn# = .5
un# = 3! / (x#(n%, ix%) - x#(n% - 1, ix%)) * (ypn# - (x#(n%, iy%) - x#(n% - 1, iy%)) / (x#(n%, ix%) - x#(n% - 1, ix%)))
y2#(n%, iy%) = (un# - qn# * u#(n% - 1)) / (qn# * y2#(n% - 1, iy%) + 1!)
'-------------------------------------------------------------------------

'lllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
' this is the backsubstitution loop of the tridiagonal algorithm
' + calculation of the standard deviation sta
sta = 0: FOR k% = n% - 1 TO 1 STEP -1
y2#(k%, iy%) = y2#(k%, iy%) * y2#(k% + 1, iy%) + u#(k%)

xx# = (x#(k% + 1, ix%) + x#(k%, ix%)) / 2
CALL splint(xx#, x#(k%, ix%), x#(k% + 1, ix%), x#(k%, iy%), x#(k% + 1, iy%), y2#(k%, iy%), y2#(k% + 1, iy%), yy#)
dlta = yy# - (x#(k%, iy%) + (x#(k% + 1, iy%) - x#(k%, iy%)) * (xx# - x#(k%, ix%)) / (x#(k% + 1, ix%) - x#(k%, ix%)))
sta = sta + dlta * dlta
         NEXT k%
'lllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll


END SUB

