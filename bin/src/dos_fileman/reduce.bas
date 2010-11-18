DECLARE SUB analizecommand (file1$)
DECLARE SUB inputline (n!, d1#(), col1%)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "REDUCE REDUCE REDUCE REDUCE REDUCE REDUCE REDUCE REDUCE REDUCE REDUCE"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*********************************************************************"
DATA "this program is designed to reduce data by deleting close data points"
DATA "use as:
DATA " REDUCE *.* 2 dmin=0.2 [options] .... takes 2nd axis and looks if there
DATA "			       are points closer than 0.2, if yes:"
DATA "			       option  
DATA "			                only the first point of them is kept
DATA "			       /av      median of points is calculated and kept
DATA " format of file:                                                          "
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "*************************************************************************"
DIM text$(300), x#(30, 100), xm#(30), in#(30)

'aaaaa analyse command string (get filename$ - ii% and dmin) aaaaaa
A$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))
I = INSTR(LCASE$(A$), "dmin=")      'check dmin
dmin = ABS(VAL(RIGHT$(A$, LEN(A$) - I - 4)))

filename$ = RTRIM$(LEFT$(A$, INSTR(A$, " "))): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - INSTR(A$, " ")))
ii% = ASC(LEFT$(A$, 1)) MOD 48               'get column number

999 CALL analizecommand(filename$)

IF RIGHT$(filename$, 4) = ".rcp" OR RIGHT$(filename$, 4) = ".mrc" THEN
   PRINT "you should never ever change data in *.rcp/mrc files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ALA$: IF LCASE$(ALA$) = "n" THEN END
   IF LCASE$(ALA$) <> "y" GOTO 22
END IF
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'oooooooooooooooooo OPEN FILE oooooooooooooooooooooooooooooooooooooooooooooo
'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "REDUCE.red"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " column"; ii%; "REDUCED to minimal step "; dmin; " in col "; ii%;
PRINT #2, "by program REDUCE.bas}"
'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
I = 0
' x#(ii%, 0) = 1E+30: x#(ii%, 1) = 1E+30

'input column y until abs(x#(ii%,i)-x#(ii%,1))>dmin iiiiiiiiiiii
4 WHILE ABS(x#(ii%, I) - x#(ii%, 1)) <= dmin
     IF EOF(1) <> 0 GOTO 5
 I = I + 1: IF I > 99 THEN PRINT "ERROR - dmin too big - program REDUCE terminated !!!!": END
 CALL inputline(1, in#(), col%): FOR coll% = 1 TO col%: x#(coll%, I) = in#(coll%): NEXT coll%
  WEND
'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

' REDUCE and write to output file  ddddddddddddddddddddddddddd
FOR coll% = 1 TO col%
 'if no options just take a point in about the middle of the interval
  xm#(coll%) = x#(coll%, INT(I / 2) + 1)
 
 'if option /av just take the average of the interval dmin
  IF INSTR(A$, "/av") > 0 THEN
   xm#(coll%) = 0: FOR k% = 1 TO I - 1: xm#(coll%) = xm#(coll%) + x#(coll%, k%): NEXT: xm#(coll%) = xm#(coll%) / (I - 1)
  END IF
NEXT coll%

'write result to file
IF ABS(x#(ii%, 1) - x#(ii%, 0)) > dmin THEN
 FOR coll% = 1 TO col%: PRINT #2, x#(coll%, 1); : NEXT: PRINT #2,
 IF ABS(x#(ii%, 1) - xm#(ii%)) > dmin THEN FOR coll% = 1 TO col%: PRINT #2, xm#(coll%); : NEXT: PRINT #2,
ELSE
 FOR coll% = 1 TO col%: PRINT #2, xm#(coll%); : NEXT: PRINT #2,
END IF
'ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
FOR coll% = 1 TO col%: x#(coll%, 0) = xm#(coll%)
                       x#(coll%, 1) = x#(coll%, I): NEXT coll%: I = 1
GOTO 4 '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

5 FOR coll% = 1 TO col%: PRINT #2, x#(coll%, I); : NEXT: PRINT #2,
CLOSE 1, 2
SHELL "copy REDUCE.red " + filename$
SHELL "del REDUCE.red"
PRINT
PRINT "END REDUCE in file "; filename$;
PRINT " column"; ii%; "REDUCED to minimal step "; dmin; " in col "; ii%;
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END
333 FOR I = 1 TO 16: READ A$: PRINT A$: NEXT I: END

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

SUB inputline (n, d1#(), col1%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
A$ = INKEY$: IF A$ <> "" THEN IF ASC(A$) = 27 THEN END

INPUT #n, ALA$
col1% = 0
 WHILE LEN(ALA$) > 0
    col1% = col1% + 1
    ALA$ = LTRIM$(ALA$) + " "
    d1#(col1%) = VAL(LEFT$(ALA$, INSTR(ALA$, " ")))
    ALA$ = LTRIM$(RIGHT$(ALA$, LEN(ALA$) - INSTR(ALA$, " ")))
 WEND

END SUB

