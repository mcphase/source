DECLARE SUB analizecommand (file1$)
DECLARE SUB newton (yy#, xx#, xsave#(), nn%, i%, ix%)
DECLARE SUB splint (xx#, xl#, xh#, yl#, yh#, y2l#, y2h#, yy#)
DECLARE SUB inputline (N!, xm#(), col%)
DECLARE SUB headerinput (text$(), j!, N!)
PRINT "INTERPOL INTERPOL INTERPOL INTERPOL INTERPOL INTERPOL INTERPOL"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program interpol - use it like INTERPOL *.* 2 stp=0.014 spl=3"
DATA "    (means take file *.* second column as xaxis and calculate            "
DATA "    linear interpolation for all other column with stepwidth 0.014) -    "
DATA " ---> the result is written in file *.*                                  "
DATA "   to do a spline interpolation use programm spline to calculate         "
DATA "   the curvature column and the add the option spl=3 (means y column     "
DATA "   is column 3 curvature column is the last column)                      "
DATA "   newton=5 instead means use newtons interpolation formula with
DATA "   5 nearest points
DATA " format of file                                                          "
DATA ""                                                                   
DATA " { header: this is the                                                  "
DATA "   file header delimited                                                "
DATA "   by brackets after this header there follow 2 or more data columns }  "
DATA " 11 3.14235 65367                                                       "
DATA "  .    .     .                                                          "
DATA "  .    .     .                                                          "
DATA "  .    .     .    .  .   .                                              "
DATA "  .    .     .                                                          "
DATA " 32 2412.34 324.2                                                       "
DATA ""
DATA "*************************************************************************"
DIM text$(300), xl#(30), xh#(30), xsave#(30, 30), in#(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(a$, INSTR(a$, " "))): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
ix% = ASC(LEFT$(a$, 1)) MOD 48
i = INSTR(a$, "spl="): IF i > 0 THEN iy% = ASC(MID$(a$, i + 4, 1)) MOD 48
                                        'get spline columns
nn% = 2: N = INSTR(a$, "newton="): IF N > 0 THEN nn% = VAL(MID$(a$, N + 7, 10))
                                        'see if newton interpolation is desired
a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, "stp=") - 3))
stp = VAL(a$)  'get stepwidth

999 CALL analizecommand(filename$)

IF RIGHT$(filename$, 4) = ".rcp" OR RIGHT$(filename$, 4) = ".mrc" THEN
   PRINT "you should never ever change data in *.rcp/mrc files": PLAY "dgdgdg"
222 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 222
END IF
'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "interpol.int"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " column "; ix%; "  was taken as the x axis for linear";
PRINT #2, "interpolating the other axis with stepwidth "; stp;
IF iy% <> 0 THEN PRINT #2, iy%; " and "; coll%; "col were taken as the yaxis and curvature for interpolating a cubic spline";
IF nn% > 2 THEN PRINT #2, nn%; " neighbouring next data points were taken for calculating a newton interpolation";
PRINT #2, "using program interpol.bas}"

FOR S% = 1 TO nn%: CALL inputline(1, in#(), coll%)
FOR i% = 1 TO coll%: xsave#(S%, i%) = in#(i%): NEXT i%
NEXT S%: S% = nn%: H% = 1
FOR i% = 1 TO coll%: xl#(i%) = xsave#(1, i%): NEXT i%
REM input data columns iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

WHILE H% <> S%
22 IF EOF(1) = 0 THEN
 insave# = in#(ix%): CALL inputline(1, in#(), col%): IF in#(ix%) = insave# GOTO 22
 S% = S% + 1: IF S% > 2 * nn% THEN S% = 1
END IF
H% = H% + 1: IF H% > 2 * nn% THEN H% = 1
FOR i% = 1 TO col%: xsave#(S%, i%) = in#(i%): xh#(i%) = xsave#(H%, i%): NEXT i%

IF SGN(xh#(ix%) - xl#(ix%)) <> 0 THEN stp = ABS(stp) * SGN(xh#(ix%) - xl#(ix%))
estp = stp: IF ABS(stp) > ABS(xh#(ix%) - xl#(ix%)) THEN estp = 0
 FOR xx# = xl#(ix%) TO xh#(ix%) - estp STEP stp
  FOR i% = 1 TO coll%
   IF i% <> ix% THEN
    IF i% = iy% THEN   'spline interpolation
     CALL splint(xx#, xl#(ix%), xh#(ix%), xl#(iy%), xh#(iy%), xl#(coll%), xh#(coll%), yy#)
    ELSE               'newton/linear interpolation
     CALL newton(yy#, xx#, xsave#(), nn%, i%, ix%)
     REM yy# = xl#(i%) + (xx# - xl#(ix%)) / (xh#(ix%) - xl#(ix%)) * (xh#(i%) - xl#(i%))
    END IF
   
    IF iy% = 0 THEN 'rem print last column only if spl is not used
      PRINT #2, yy#;
    ELSE
     IF i% < coll% THEN PRINT #2, yy#;
    END IF

   ELSE
    PRINT #2, xx#;
   END IF
  NEXT i%: PRINT #2,
 NEXT xx#


FOR i% = 1 TO coll%: xl#(i%) = xh#(i%): NEXT i%
WEND

'print last point
FOR i% = 1 TO coll% - 1: PRINT #2, xh#(i%); : NEXT
IF iy% = 0 THEN PRINT #2, xh#(coll%);
PRINT #2, : CLOSE 'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

SHELL "copy interpol.int " + filename$
SHELL "del interpol.int"

PRINT
PRINT "END INTERPOL in file "; filename$; " column "; ix%; "  was taken as the x axis for linear ";
PRINT "interpolating the other axis with stepwidth "; stp;
IF iy% <> 0 THEN PRINT ", the "; iy%; " and "; coll%; "col were taken as the yaxis and curvature for interpolating a cubic spline"
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
                                           
   FOR j = 2 TO 300
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
'this sub calculates the newton interpolation yy# at
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

'sum up newton interpolation formula
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
H# = xh# - xl#
a# = (xh# - xx#) / H#
B# = (xx# - xl#) / H#
yy# = a# * yl# + B# * yh# + ((a# * a# * a# - a#) * y2l# + (B# * B# * B# - B#) * y2h#) * H# * H# / 6!

END SUB

