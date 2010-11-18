DECLARE SUB analizecommand (file1$)
DECLARE SUB newton (yy#, xx#, xsave#(), nn%, i%, ix%)
DECLARE SUB splint (xx#, xl#, xh#, yl#, yh#, y2l#, y2h#, yy#)
DECLARE SUB inputline (N!, xm#(), col%)
DECLARE SUB headerinput (text$(), j!, N!)
PRINT "CONVOLUT CONVOLUT CONVOLUT CONVOLUT CONVOLUT CONVOLUT CONVOLUT"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "program CONVOLUT-use it like CONVOLUT *.* 23 stp=0.014 mode=gauss fwhm=0.2"
DATA "    (means take file *.* second column as xaxis and calculate            "
DATA "    CONVOLUTation for y axis=3rd column  with stepwidth 0.014) -    "
DATA " ---> the result is written in file *.cvt                         "
DATA "   use the mode= switch to specify the type of convolution function     "
DATA "   possible modes are: gauss, fullprof2t, fullprofq"
DATA "   pars: for gauss fwhm=, for fullprof* u= v= w=, fir fullprofq lamdba=[A]"
DATA "    [u,v,w used to calculate fwhm from 2theta: for D9-ILL u=16.08575,
DATA "     v=-6.2195, w=0.83903]
DATA "   used (gauss) and then specify the function parameters of f(x)    "
DATA "   convolution is done according to the formula:             "
DATA "    conv(x')=sum_{x} col3(x)*f(x'-x)"
DATA " "
DATA " format of file                                                          "
DATA ""                                       
DATA " { header: this is the                                                  "
DATA "   file header delimited                                                "
DATA "   by brackets after this header there follow 2 or more data columns }  "
DATA " 11 3.14235 65367                                                       "
DATA "  .    .     .                                                          "
DATA "  .    .     .    .  .   .                                              "
DATA "  .    .     .                                                          "
DATA " 32 2412.34 324.2                                                       "
DATA "*************************************************************************"
DIM text$(600), x#(1000), y#(1000), in#(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(a$, INSTR(a$, " "))): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
ix% = ASC(LEFT$(a$, 1)) MOD 48: iy% = ASC(MID$(a$, 2, 1)) MOD 48
stp = VAL(MID$(a$, INSTR(a$, "stp=") + 4))  'get stepwidth
mode$ = LTRIM$(MID$(a$, INSTR(a$, "mode=") + 5))
IF LEFT$(mode$, 5) = "gauss" THEN fwhm = VAL(MID$(a$, INSTR(a$, "fwhm=") + 5))
IF LEFT$(mode$, 8) = "fullprof" THEN
 u = VAL(MID$(a$, INSTR(a$, "u=") + 2))
 v = VAL(MID$(a$, INSTR(a$, "v=") + 2))
 w = VAL(MID$(a$, INSTR(a$, "w=") + 2))
END IF
IF LEFT$(mode$, 9) = "fullprofq" THEN
 lambda = VAL(MID$(a$, INSTR(a$, "lambda=") + 7))
END IF
999 CALL analizecommand(filename$)
'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'open input file and input file header
SHELL "copy " + filename$ + " CONVOLUT.cvt"
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)
CALL inputline(1, in#(), col%)
xmin# = in#(ix%): xmax# = in#(ix%)
WHILE EOF(1) = 0
IF xmax# < in#(ix%) THEN xmax# = in#(ix%)
IF xmin# > in#(ix%) THEN xmin# = in#(ix%)
CALL inputline(1, in#(), col%)
WEND: CLOSE 1


IF INSTR(filename$, ".") > 0 THEN
  cvtfile$ = LEFT$(filename$, INSTR(filename$, ".")) + "cvt"
 ELSE
  cvtfile$ = filename$ + ".cvt"
END IF
' open output file and write fileheader
OPEN "o", 2, cvtfile$
PRINT #2, "#{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, "#"; DATE$; " "; TIME$; " column "; ix%; "  was taken as the x axis for ";
PRINT #2, "CONVOLUTing the y axis with stepwidth "; stp;
PRINT #2, iy%; " was taken as the yaxis - convolution was done with a  "; mode$; "-function"
PRINT #2, "#the command was: "; COMMAND$;
PRINT #2, "- used program: CONVOLUT.bas}"

REM write output data columns IN BLOCKS OF 1000 iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
xmin# = xmin# - 1.1 * (xmax# - xmin#)
xmax# = xmax# + 1.1 * (xmax# - xmin#)
x#(0) = xmin#: sign% = SGN(xmax# - xmin#)
22 h% = 0: WHILE x#(h%) * sign% < xmax# * sign% AND h% < 1000'setup block of xdata
             h% = h% + 1: x#(h%) = x#(h% - 1) + stp * sign%: y#(h%) = 0
           WEND: x#(0) = x#(h%)


OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)
'CALL inputline(1, in#(), col%)
WHILE EOF(1) = 0
CALL inputline(1, in#(), col%)
FOR i% = 1 TO h%    'sum up data block

IF LEFT$(mode$, 5) = "gauss" THEN    'GAUSS CONVOLUTION
 IF ABS((x#(i%) - in#(ix%)) / fwhm) < 1E+19 THEN
   dd = (x#(i%) - in#(ix%)) * (x#(i%) - in#(ix%)) / (fwhm * fwhm / 4 / LOG(.5))
   y#(i%) = y#(i%) + in#(iy%) * EXP(dd) / fwhm
 END IF
END IF

IF LEFT$(mode$, 10) = "fullprof2t" THEN    'fullprof convolution (2theta)
IF ABS((x#(i%) - in#(ix%)) / .1) < 1E+19 THEN
   zt = in#(ix%) * 3.14159 / 180 'zwei theta in rad
'   u = 16.08575
'   v = -6.2195
'   w = .83903
  IF u * TAN(zt / 2) * TAN(zt / 2) + v * TAN(zt / 2) + w > 0 THEN
    fwhm = SQR(u * TAN(zt / 2) * TAN(zt / 2) + v * TAN(zt / 2) + w)
    dd = (x#(i%) - in#(ix%)) * (x#(i%) - in#(ix%)) / (fwhm * fwhm / 4 / LOG(.5))
    y#(i%) = y#(i%) + in#(iy%) * EXP(dd) / fwhm
  ELSE
   PRINT "warning: at x="; in#(ix%); "fwhm given by parameters uvw gets negative ... continuing"
  END IF
 END IF

END IF

IF LEFT$(mode$, 9) = "fullprofq" THEN    'fullprof convolution (q)
IF ABS((x#(i%) - in#(ix%)) / .1) < 1E+19 THEN
   d = 2 * 3.1415 / in#(ix%)
 '  lambda = .58
   tantheta = lambda / SQR(4 * d * d - lambda * lambda)
   'zt = in#(ix%) * 3.14159 / 180 'zwei theta in rad
 '  u = 26 '16.08575
 '  v = 10 ' -6.2195
 '  w = .83903
 IF u * tantheta * tantheta + v * tantheta + w > 0 THEN
   fwhm = SQR(u * tantheta * tantheta + v * tantheta + w)
   fwhm = 3.1415 / 180 * fwhm'transform to rad
   fwhm = fwhm * 2 * 3.1415 / lambda * SQR(1 - lambda * lambda / 2 / 2 / d / d)
                'transform to q
   dd = (x#(i%) - in#(ix%)) * (x#(i%) - in#(ix%)) / (fwhm * fwhm / 4 / LOG(.5))
   y#(i%) = y#(i%) + in#(iy%) * EXP(dd) / fwhm
 ELSE
  PRINT "warning: at x="; in#(ix%); "fwhm given by parameters uvw gets negative ... continuing"
 END IF
 

 END IF

END IF

NEXT i%
WEND: CLOSE 1
FOR i% = 1 TO h%
IF LEFT$(mode$, 10) = "fullprof2t" THEN    'fullprof convolution (2theta)
zt = x#(i%)
bkpos = 2
klammer = zt / bkpos - 1
bk = 1466.6 * klammer ^ 0
bk = bk - 198.7 * klammer
bk = bk + 16.942 * klammer ^ 2
bk = bk + .53499 * klammer ^ 3
bk = bk - .12783 * klammer ^ 4
bk = bk + .0041838 * klammer ^ 5
'y#(i%) = y#(i%) + bk
END IF
nn$ = STR$(x#(i%)): IF INSTR(nn$, "D") > 0 THEN MID$(nn$, INSTR(nn$, "D"), 1) = "E"
PRINT #2, " " + nn$;
nn$ = STR$(y#(i%)): IF INSTR(nn$, "D") > 0 THEN MID$(nn$, INSTR(nn$, "D"), 1) = "E"
PRINT #2, " " + nn$

NEXT i%
IF h% >= 999 GOTO 22

CLOSE 2 'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

SHELL "del CONVOLUT.cvt"

PRINT
PRINT "END CONVOLUT in file "; filename$; " column "; ix%; "  was taken as the x axis for";
PRINT "convolution ofaxis  "; iy%; " axis with stepwidth "; stp
PRINT
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR i = 1 TO 21: READ a$: PRINT a$: NEXT i: END

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

SUB headerinput (text$(), j, N)
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'**********************************************************************

1 INPUT #N, a$
   i = INSTR(a$, "{"): IF i > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(a$, LEN(a$) - i)
  j = 1: i = INSTR(text$(j), "}"): IF i > 0 GOTO 3  'look for "}" in first line
                                           
   FOR j = 2 TO 600
   INPUT #N, text$(j)
   i = INSTR(text$(j), "}"): IF i > 0 GOTO 3      'look for "}"
   NEXT j: PRINT "text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

END SUB

SUB inputline (N, xm#(), col%)
'**********************************************************************
'diese sub liest eine zeile von file #n, bestimmt die
' anzahl der datenspalten col% und sichert die daten in xm#(1,..,col%)
'**********************************************************************
 INPUT #N, ala$        'input data point line as string and split into numbers
          col% = 0: WHILE LEN(ala$) > 0: ala$ = LTRIM$(ala$) + " ": col% = col% + 1
 xm#(col%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
          ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " "))): WEND
END SUB

SUB newton (yy#, xx#, xsave#(), nn%, i%, ix%)
'this sub calculates the newton CONVOLUTation yy# at
'xx# using the data  points stored in xsave#
DIM a%(2 * nn%), fr#(nn%, nn% - 1)
FOR S% = 1 TO 2 * nn%: a%(S%) = S%: NEXT S%

'look for nn% nearest data points to xx#
44 dd# = 0: FOR S% = 1 TO 2 * nn% - 1
ddold# = dd#: dd# = ABS(xx# - xsave#(a%(S%), ix%))
IF ddold# > dd# THEN asave% = a%(S%): a%(S%) = a%(S% - 1): a%(S% - 1) = asave%: GOTO 44
NEXT S%

'calculaste divided differences
FOR S% = 1 TO nn%: fr#(S%, 0) = xsave#(a%(S%), i%): NEXT S%

FOR r% = 1 TO nn% - 1
FOR S% = 1 TO nn% - r%
fr#(S%, r%) = (fr#(S% + 1, r% - 1) - fr#(S%, r% - 1)) / (xsave#(a%(S% + r%), ix%) - xsave#(a%(S%), ix%))
NEXT: NEXT

'sum up newton CONVOLUTation formula
yy# = 0: FOR r% = nn% - 1 TO 1 STEP -1
yy# = (yy# + fr#(1, r%)) * (xx# - xsave#(a%(r%), ix%))
NEXT r%: yy# = yy# + fr#(1, 0)

REM yy# = xl#(i%) + (xx# - xl#(ix%)) / (xh#(ix%) - xl#(ix%)) * (xh#(i%) - xl#(i%))

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

