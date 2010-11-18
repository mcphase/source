DECLARE SUB analizecommand (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "DIFF DIFF DIFF DIFF DIFF DIFF DIFF DIFF DIFF DIFF DIFF DIFF DIFF"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program diff - use it like diff *.* 23"                     
DATA "    (23 means differentiate 3rd columns with respect to second column) - "
DATA " ---> the result is written in file *.*,                                 "
DATA " the differentiation is done point by point unless you enter             "
DATA " the programm like diff *.* 23 stp=2.34, which means the differentiation"
DATA " is done in col2-intervals of 2.34 (+1 datapoint extra)                  "
DATA " te differentiation is done by sampling the xaxis(2nd column) in         "
DATA " steps of 2.34(or point by point) and caculating the linear regression   "
DATA " at the average <x> of such an xaxis interval for all other columns !!   "
DATA " to the yaxis (3rd column in our example) the slope of Its linear regress."
DATA " is written, to all other axis (z1,z2..) the value of the regression     "
DATA " at <x> is written (which corresponds to <z1>,<z2> ...in the interval)   "
DATA " NOTE: all columns are changed by this operation -                       "
DATA " format of file:                                                          "
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "*************************************************************************"
DIM text$(300), x#(30, 100), xm#(30)

'aaaaa analyse command string (get filename$ - ii% - jj% and interval) aaaaaa
A$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))
i = INSTR(LCASE$(A$), "stp=")      'check if option stp is used
IF i > 0 THEN interval = ABS(VAL(RIGHT$(A$, LEN(A$) - i - 3))): A$ = RTRIM$(LEFT$(A$, i - 1))

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
OPEN "o", 2, "diff.dif"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " column"; jj%; "differentiated with respect to column"; ii%;
PRINT #2, "by program diff.bas}"
'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
i = 0: x#(ii%, i) = 0: x#(ii%, 1) = 0

'input column y until abs(x#(ii%,i)-x#(ii%,1))>interval iiiiiiiiiiii
4 WHILE ABS(x#(ii%, i) - x#(ii%, 1)) <= interval
     IF EOF(1) <> 0 GOTO 5
INPUT #1, ala$        'input data point line as string and split into numbers
          col% = 0: WHILE LEN(ala$) > 0: ala$ = LTRIM$(ala$) + " ": col% = col% + 1
x#(col%, i + 1) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
          ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " "))): WEND
 i = i + 1: IF i > 99 THEN PRINT "ERROR - stp too big - program diff terminated !!!!": END
  WEND
'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

' differentiate and write to output file  ddddddddddddddddddddddddddd
' for the untouched columns just take the average of the intervaldddd
FOR coll% = 1 TO col%
 IF coll% <> jj% THEN
  xm#(coll%) = 0: FOR k% = 1 TO i: xm#(coll%) = xm#(coll%) + x#(coll%, k%): NEXT: xm#(coll%) = xm#(coll%) / i
 END IF
NEXT coll%

' differentiate column jj% with respect to col ii%
ym# = 0: FOR k% = 1 TO i: ym# = ym# + x#(jj%, k%): NEXT: ym# = ym# / i
xym# = 0: FOR k% = 1 TO i: xym# = xym# + x#(ii%, k%) * x#(jj%, k%): NEXT: xym# = xym# / i
xxm# = 0: FOR k% = 1 TO i: xxm# = xxm# + x#(ii%, k%) * x#(ii%, k%): NEXT: xxm# = xxm# / i

' here calculate the slope of the column jj% with respect to column ii%
xm#(jj%) = (xym# - ym# * xm#(ii%)) / (xxm# - xm#(ii%) * xm#(ii%))

'write result to file
IF washere% = 0 THEN
washere% = 1: FOR coll% = 1 TO col%: IF coll% <> jj% THEN PRINT #2, x#(coll%, 1);  ELSE PRINT #2, xm#(coll%);
NEXT: PRINT #2,
END IF
FOR coll% = 1 TO col%:
nn$ = STR$(xm#(coll%)): MID$(nn$, INSTR(nn$, "D"), 1) = "E"
PRINT #2, " " + nn$; : NEXT: PRINT #2,
'ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

'delete first point of interval
FOR k% = 1 TO i - 1: FOR coll% = 1 TO col%
x#(coll%, k%) = x#(coll%, k% + 1)
NEXT coll%: NEXT k%
i = i - 1
GOTO 4 '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

5 FOR coll% = 1 TO col%: IF coll% <> jj% THEN PRINT #2, x#(coll%, i);  ELSE PRINT #2, xm#(coll%);
NEXT: PRINT #2,
CLOSE 1, 2
SHELL "copy diff.dif " + filename$
SHELL "del diff.dif"
PRINT
PRINT "END DIFF in file "; filename$; " column "; jj%; "has been differentiated with respect to";
PRINT "column "; ii%;

IF INSTR(COMMAND$, "*") <> 0 GOTO 999

END
333 FOR i = 1 TO 22: READ A$: PRINT A$: NEXT i: END

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
   i = INSTR(A$, "{"): IF i > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(A$, LEN(A$) - i)
  j = 1: i = INSTR(text$(j), "}"): IF i > 0 GOTO 3  'look for "}" in first line
                                            
   FOR j = 2 TO 300
   INPUT #n, text$(j)
   i = INSTR(text$(j), "}"): IF i > 0 GOTO 3      'look for "}"
   NEXT j: PRINT "text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

END SUB

