DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB iterate (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "INT INT INT INT INT INT INT INT INT"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program INT - use it like INT *.* 23                   "
DATA "    (means integrate function y(x) with xcolumn 2 and y column 3 "
DATA "     from file *.*) -                                                 "
DATA " format of file                                                          "
DATA ""
DATA " { header: this is the                                                  "
DATA "   file header delimited                                                "
DATA "   by brackets after this header there follow 3 or more data columns }  "
DATA " 11 3.14235 65367                                                       "
DATA "  .    .     .                                                          "
DATA "  .    .     .    .  .   .                                              "
DATA "  .    .     .                                                          "
DATA " 32 2412.34 324.2                                                       "
DATA ""
DATA "*************************************************************************"
DIM text$(300), xm#(300), x$(300)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
A$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(A$, INSTR(A$, " "))): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - INSTR(A$, " ")))
x% = ASC(LEFT$(A$, 1)) MOD 48 'get columns
y% = ASC(MID$(A$, 2, 1)) MOD 48

999 CALL iterate(filename$)

IF RIGHT$(filename$, 4) = ".rcp" OR RIGHT$(filename$, 4) = ".mrc" THEN
   PRINT "you should never ever change data in *.rcp/mrc files": PLAY "dgdgdg"
222 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 222
END IF
'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "INT.del"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; "column"; y%;
PRINT #2, "has been integrated with respect to column "; x%; " using program INT.bas}"

firstpoint% = 1

 REM input data columns
WHILE EOF(1) = 0
 CALL inputline(1, xm#(), col%, x$(), cols%)
 IF col% = -1 THEN  'take care of inserted comments
  CALL headerinput(text$(), j, 1)
  PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
  PRINT #2, DATE$; " "; TIME$; " column"; y%;
  PRINT #2, "has been integrated with respect to column "; y%; " using program INT.bas}"
  firstpoint% = 1
 END IF

IF firstpoint% = 1 THEN
 integral# = 0
 firstpoint% = 0
 dx = 0
ELSE
 dx = xm#(x%) - xold#
 integral# = integral# + (xm#(y%) + yold#) / 2 * dx
END IF
 
 xold# = xm#(x%)
 yold# = xm#(y%)
 xm#(y%) = integral#

 'write result to file
  FOR coll% = 1 TO col%: PRINT #2, xm#(coll%); : NEXT
  FOR coll% = 1 TO cols%: PRINT #2, "{"; x$(coll%); "}"; : NEXT
  PRINT #2,

WEND

 CLOSE 1, 2
SHELL "copy INT.del " + filename$
SHELL "del INT.del"

PRINT
PRINT "END INT in file "; filename$; "column "; y%; "has been integrated with respect to column "; x%; " deleted"
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR i = 1 TO 16: READ A$: PRINT A$: NEXT i: END

SUB headerinput (text$(), j, n)
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'**********************************************************************

1 LINE INPUT #n, A$
   i = INSTR(A$, "{"): IF i > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(A$, LEN(A$) - i)
  j = 1: i = INSTR(text$(j), "}"): IF i > 0 GOTO 3  'look for "}" in first line
                                          
   FOR j = 2 TO 300
   LINE INPUT #n, text$(j)
   i = INSTR(text$(j), "}"): IF i > 0 GOTO 3      'look for "}"
   NEXT j: PRINT "text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

END SUB

SUB inputline (n, d1#(), col1%, d2$(), col2%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
'comments in {} are stored in d2$
'if a comment is started somewhere in this line by "{" and not finished
'then col1% is set -1 and the filepointer is set to the beginning of the
'line (with seek)

A$ = INKEY$: IF A$ <> "" THEN IF ASC(A$) = 27 THEN END
IF SCREEN(24, 10) <> 45 THEN PRINT "-";  ELSE LOCATE 24, 1: PRINT SPACE$(15); : LOCATE 24, 1

aa = SEEK(n)
LINE INPUT #n, ala$
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

SUB iterate (file1$)
STATIC washere%

  IF INSTR(COMMAND$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere% = 0 THEN
      washere% = 1
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > int.dir"
      OPEN "i", 9, "int.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del int.dir": END
   INPUT #9, file1$
   IF LCASE$(file1$) = "int.dir" THEN
    IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del int.dir": END
    INPUT #9, file1$
   END IF
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del int.dir": END
  END IF


END SUB

