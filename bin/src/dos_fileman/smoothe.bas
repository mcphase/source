DECLARE SUB analizecommand (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "SMOOTHE SMOOTHE SMOOTHE SMOOTHE SMOOTHE SMOOTHE SMOOTHE SMOOTHE SMOOTHE SMOOTHE SMOOTHE SMOOTHE SMOOTHE"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program SMOOTHE - use it like SMOOTHE *.* 23 stp=4.46 "
DATA "    (23 means SMOOTHE 3rd column with respect to second column) - "
DATA " ---> the result is written in file *.*,                                 "
DATA " smoothing is done in col2-intervals of 2.34 (+1 datapoint extra)        "
DATA " te SMOOTHing  is done by sampling the xaxis(2nd column) in      "
DATA " steps of 2.34(or point by point) and caculating the linear regression   "
DATA " at the average <x> of such an xaxis interval, if a data point exceeds "
DATA " the standard deviation of this linear regression then its value is changed"
DATA " to the standard deviation but only within a small interval around <x>   "
DATA " format of file:                                                          "
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "*************************************************************************"
DIM text$(300), x#(30, 200), xm#(30)

'aaaaa analyse command string (get filename$ - ii% - jj% and interval) aaaaaa
A$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))
I = INSTR(LCASE$(A$), "stp=")      'check if option stp is used
interval = ABS(VAL(RIGHT$(A$, LEN(A$) - I - 3))): A$ = RTRIM$(LEFT$(A$, I - 1))

filename$ = RTRIM$(LEFT$(A$, LEN(A$) - 2))
jj% = ASC(RIGHT$(A$, 1)) MOD 48: ii% = ASC(MID$(A$, LEN(A$) - 1, 1)) MOD 48  'get column numbers

999 CALL analizecommand(filename$)

IF RIGHT$(filename$, 4) = ".mrc" OR RIGHT$(filename$, 4) = ".rcp" THEN
   PRINT "you should never ever change data in *.rcp/mrc files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 22
END IF
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'oooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooooooooooooooooooooo
'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "SMOOTHE.SMO"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " column"; jj%; "SMOOTHED with respect to column"; ii%; "stp="; interval
PRINT #2, "by program SMOOTHE.bas}"
'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
I = 0

'input column y until abs(x#(ii%,i)-x#(ii%,1))>interval iiiiiiiiiiii
4 WHILE ABS(x#(ii%, I) - x#(ii%, 1)) <= interval
     IF EOF(1) <> 0 GOTO 5
INPUT #1, ala$        'input data point line as string and split into numbers
          col% = 0: WHILE LEN(ala$) > 0: ala$ = LTRIM$(ala$) + " ": col% = col% + 1
x#(col%, I + 1) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
          ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " "))): WEND
 I = I + 1: IF I > 199 THEN PRINT "ERROR - stp too big - program SMOOTHE terminated !!!!": END
  WEND
'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    
     IF I < 4 GOTO 6   'continue if there are too few points in the interval

' SMOOTHE and write to output file  ddddddddddddddddddddddddddddddddddddd
FOR coll% = 1 TO col%
  xm# = 0: FOR k% = 1 TO I: xm# = xm# + x#(ii%, k%): NEXT: xm# = xm# / I
NEXT coll%

' SMOOTHEerentiate column jj% with respect to col ii%
ym# = 0: FOR k% = 1 TO I: ym# = ym# + x#(jj%, k%): NEXT: ym# = ym# / I
xym# = 0: FOR k% = 1 TO I: xym# = xym# + x#(ii%, k%) * x#(jj%, k%): NEXT: xym# = xym# / I
xxm# = 0: FOR k% = 1 TO I: xxm# = xxm# + x#(ii%, k%) * x#(ii%, k%): NEXT: xxm# = xxm# / I

' here calculate the slope of the column jj% with respect to column ii%
slope# = (xym# - ym# * xm#) / (xxm# - xm# * xm#)
offset# = ym# - slope# * xm#
sta# = 0: FOR k% = 1 TO I
dy# = x#(jj%, k%) - (offset# + slope# * x#(ii%, k%)): sta# = sta# + dy# * dy#
         NEXT k%: sta# = SQR(sta# / (I - 2))
'smoothe
FOR k% = 1 TO I: dy# = x#(jj%, k%) - (offset# + slope# * x#(ii%, k%))
IF ABS(x#(ii%, k%) - xm#) < interval / 4 AND ABS(dy#) > sta# THEN
     x#(jj%, k%) = x#(jj%, k%) - dy# / 2
END IF
NEXT k%

6 'write result to file
FOR coll% = 1 TO col%: PRINT #2, x#(coll%, 1); : NEXT: PRINT #2,
'ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

'delete first point of interval
FOR k% = 1 TO I - 1: FOR coll% = 1 TO col%
x#(coll%, k%) = x#(coll%, k% + 1)
NEXT coll%: NEXT k%
I = I - 1
GOTO 4 '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

5 FOR k% = 1 TO I: FOR coll% = 1 TO col%: PRINT #2, x#(coll%, k%); : NEXT: PRINT #2, : NEXT
CLOSE 1, 2
SHELL "copy SMOOTHE.SMO " + filename$
SHELL "del SMOOTHE.SMO"
PRINT
PRINT "END SMOOTHE in file "; filename$; " column "; jj%; "has been SMOOTHEd with respect to";
PRINT "column "; ii%;
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END
333 FOR I = 1 TO 18: READ A$: PRINT A$: NEXT I: END

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

SUB headerinput (text$(), j, n)
'*********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'*********************************************************************

1 INPUT #n, A$
   I = INSTR(A$, "{"): IF I > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(A$, LEN(A$) - I)
  j = 1: I = INSTR(text$(j), "}"): IF I > 0 GOTO 3  'look for "}" in first line
                                            
   FOR j = 2 TO 300
   INPUT #n, text$(j)
   I = INSTR(text$(j), "}"): IF I > 0 GOTO 3      'look for "}"
   NEXT j: PRINT "text in data file too long": END
3 text$(j) = LEFT$(text$(j), I - 1)

END SUB

